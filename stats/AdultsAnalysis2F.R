# LIBRARY IMPORTS ==================================================================================
library(lme4)
library(lmerTest)
library(emmeans)
library(brms)
library(coda)
library(nortest)
library(tidyverse)
library(broom)
library(ggeffects)
library(eyetrackingR)
library(RColorBrewer)

source("Routines.R")
source("StatTools.R")
source("geom_flat_violin.R")

# GATHER DATA ======================================================================================
d <- LT_data.gather("adults_2f")
# Unload snow packages so that parallel works for brms
detach("package:doSNOW")
detach("package:snow")
# Get behavioural and LT data, excluding outliers in terms of
# number of blocks before learning (graphically, from boxplot)
behaviour <- d[[2]] %>%
  mutate(Condition = relevel(Condition, ref = "No Label")) %>%
  subset((Condition == "Label" & NBlocks < 8) |
           (Condition == "No Label" & NBlocks < 5))
rote_learning_check <- behaviour %>%
  subset(Phase == "Test") %>%
  group_by(Participant, Stimulus) %>%
  summarise(a = ACC) %>%
  subset(a == 0) # Participants with mistakes at test (including missed item)
# No rote learners
LT.clean <- d[[4]] %>%
  mutate(Condition = relevel(Condition, ref = "No Label")) %>%
  subset((Condition == "Label" & NBlocks < 8) |
           (Condition == "No Label" & NBlocks < 5)) %>%
  make_eyetrackingr_data(participant_column = "Participant",
                         trial_column = "TrialId",
                         time_column = "TimeStamp",
                         trackloss_column = "TrackLoss",
                         aoi_columns = c("Head", "Tail"),
                         treat_non_aoi_looks_as_missing = T)

# TRAINING LT ANALYSIS: PROP TAIL LOOKING BY PARTICIPANT ===========================================
save_path <- "../results/adults_2f/PropTail/TrialAverage_"
# DATA PREPARATION
prop_tail.per_fstlst <- LT.clean %>%
  drop_na(FstLst) %>%
  make_time_window_data(aois="Tail",
                        predictor_columns=c("Condition",
                                            "FstLst",
                                            "ACC",
                                            "Stimulus",
                                            "StimLabel"))

# Testing ArcSin ~ FstLst*Condition
run_model <- F # Running the models takes around 20 minutes on a 4.40GHz 12-core
if(run_model){
  ## Run lmer (Sampling Theory Based)
  t <- proc.time()
  prop_tail.per_fstlst.lmer.model <- lmer(ArcSin ~ FstLst*Condition +
                                            (1 + FstLst | Participant) +
                                            (1 | Stimulus) +
                                            (1 | StimLabel),
                                          data = prop_tail.per_fstlst)
  prop_tail.per_fstlst.lmer.anova <- anova(prop_tail.per_fstlst.lmer.model, type = 1)
  ## Run brms (Bayesian)
  ### Set priors for models other than intercept-only
  priors.prop_tail.per_fstlst <- list(set_prior("uniform(0,1.6)",
                                                class = "Intercept"),
                                      c(set_prior("uniform(0,1.6)",
                                                  class = "Intercept"),
                                        set_prior("normal(0,.5)", class = "b")))
  ### Set all nested formulas for model comparisons
  formulas.prop_tail.per_fstlst <- list(ArcSin ~ 1 +
                                          (1 | Participant) +
                                          (1 | Stimulus) +
                                          (1 | StimLabel),
                                        ArcSin ~ 1 + FstLst +
                                          (1 + FstLst | Participant) +
                                          (1 | Stimulus) +
                                          (1 | StimLabel),
                                        ArcSin ~ 1 + FstLst + Condition +
                                          (1 + FstLst | Participant) +
                                          (1 | Stimulus) +
                                          (1 | StimLabel),
                                        ArcSin ~ 1 + FstLst + Condition + FstLst:Condition +
                                          (1 + FstLst | Participant) +
                                          (1 | Stimulus) +
                                          (1 | StimLabel))
  ### Get brms results
  brms.results <- bayes_factor.brm_fixef(formulas.prop_tail.per_fstlst,
                                         prop_tail.per_fstlst,
                                         priors.prop_tail.per_fstlst,
                                         controls = list(adapt_delta = 0.999999999999999,
                                                         max_treedepth = 15))
  prop_tail.per_fstlst.brms.models <- brms.results[[1]]
  prop_tail.per_fstlst.brms.bayes_factors <- brms.results[[2]]
  prop_tail.time <- proc.time() - t
  ## Save all the results
  saveRDS(prop_tail.per_fstlst.lmer.model, paste0(save_path, "FstLst_lmerModel.rds"))
  saveRDS(prop_tail.per_fstlst.lmer.anova, paste0(save_path, "FstLst_lmerAnova.rds"))
  lapply(seq_along(prop_tail.per_fstlst.brms.models),
         function(i){
           saveRDS(prop_tail.per_fstlst.brms.models[[i]],
                   paste0(save_path, "FstLst_brmsModel", i, ".rds"))
         })
  saveRDS(prop_tail.per_fstlst.brms.bayes_factors, paste0(save_path, "FstLst_brmsBF.rds"))
}else{
  prop_tail.per_fstlst.lmer.model <- readRDS(paste0(save_path, "FstLst_lmerModel.rds"))
  prop_tail.per_fstlst.lmer.anova <- readRDS(paste0(save_path, "FstLst_lmerAnova.rds"))
  prop_tail.per_fstlst.brms.models <- lapply(1:4,
                                             function(i){
                                               readRDS(paste0(save_path,
                                                              "FstLst_brmsModel", i, ".rds"))
                                             })
  prop_tail.per_fstlst.brms.bayes_factors <- readRDS(paste0(save_path, "FstLst_brmsBF.rds"))
}

# PLOTTING
generate_plots <- F
## Plot jitter + mean&se + lines
if(generate_plots){
  ## Get brm predicted values (using three levels of HPDI to better appreciate data shape)
  prop_tail.raw_predictions <- last(prop_tail.per_fstlst.brms.models) %>%
    predict(summary = F,
            transform = function(x){sin(x)^2}) %>%
    t() %>%
    as_tibble() %>%
    mutate(RowNames = 1:nrow(.))
  prop_tail.predicted <- prop_tail.per_fstlst %>%
    mutate(RowNames = 1:nrow(.)) %>%
    select(FstLst, Condition, RowNames) %>%
    inner_join(prop_tail.raw_predictions) %>%
    select(-RowNames) %>%
    gather(key = Sample, value = Predicted, -c(FstLst, Condition))
  prop_tail.predicted.hpdi.97 <- prop_tail.predicted %>%
    select(-Sample) %>%
    split(list(.$FstLst, .$Condition)) %>%
    lapply(function(df){
      hpdi <- as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.97)
      df.summary <- df %>%
        group_by(FstLst, Condition) %>%
        summarise(Mean = mean(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows()
  prop_tail.predicted.hpdi.89 <- prop_tail.predicted %>%
    select(-Sample) %>%
    split(list(.$FstLst, .$Condition)) %>%
    lapply(function(df){
      hpdi <- as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.89)
      df.summary <- df %>%
        group_by(FstLst, Condition) %>%
        summarise(Mean = mean(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows()
  prop_tail.predicted.hpdi.67 <- prop_tail.predicted %>%
    select(-Sample) %>%
    split(list(.$FstLst, .$Condition)) %>%
    lapply(function(df){
      hpdi <- as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.67)
      df.summary <- df %>%
        group_by(FstLst, Condition) %>%
        summarise(Mean = mean(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows()
  ## Plot raincloud + predicted mean&HPDIs per FstLst
  prop_tail.per_fstlst.plot <- ggplot(prop_tail.per_fstlst,
                                      aes(x = Condition, y = Prop,
                                          colour = Condition,
                                          fill = Condition)) +
    theme_bw() + ylab("Looking to Tail (Prop)") +
    theme(legend.position = "top",
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    coord_flip() + facet_grid(.~FstLst) +
    geom_flat_violin(position = position_nudge(x = .2), colour = "black", alpha = .5) +
    geom_point(position = position_jitter(width = .15),
               size = 1, alpha = .6,
               show.legend = F) +
    geom_boxplot(width = .1, alpha = .3, outlier.shape = NA, colour = "black",
                 show.legend = F) +
    # geom_pointrange(data = prop_tail.predicted.hpdi.67,
    #                 aes(x = Condition, y = Mean, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = 1.5, size = 1.5,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    # geom_pointrange(data = prop_tail.predicted.hpdi.89,
    #                 aes(x = Condition, y = Mean, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = .5, size = 1,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    # geom_pointrange(data = prop_tail.predicted.hpdi.97,
    #                 aes(x = Condition, y = Mean, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = .5, size = .5,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")
  ## Save plot
  ggsave(paste0(save_path, "FstLst_data.pdf"),
         prop_tail.per_fstlst.plot,
         width = 5.5, height = 3, dpi = 600)
}

# TRAINING LT ANALYSIS: PROP TAIL LOOKING BY TRIAL PART ============================================
save_path <- "../results/adults_2f/PropTail/TrialParts_"
# DATA PREPARATION
trial_parts.per_fstlst <- LT.clean %>%
  drop_na(FstLst) %>%
  make_time_window_data(aois="Tail",
                        predictor_columns=c("Condition",
                                            "FstLst",
                                            "CurrentObject",
                                            "Stimulus",
                                            "StimLabel"))

# Testing ArcSin ~ FstLst*CurrentObject(TrialParts)*Condition
run_model <- F # Running the models takes around 4h30 on a 4.40GHz 12-core
if(run_model){
  ## Run lmer (Sampling Theory Based)
  t <- proc.time()
  trial_parts.per_fstlst.lmer.model <- lmer(ArcSin ~ FstLst*CurrentObject*Condition +
                                                     (1 + FstLst*CurrentObject | Participant) +
                                                     (1 | Stimulus) +
                                                     (1 | StimLabel),
                                                   data = trial_parts.per_fstlst)
  trial_parts.per_fstlst.lmer.anova <- anova(trial_parts.per_fstlst.lmer.model, type = 1)
  ## Run brms (Bayesian)
  ### Set priors for models other than intercept-only
  priors.trial_parts.per_fstlst <- list(set_prior("uniform(0,1.6)",
                                                  class = "Intercept"),
                                        c(set_prior("uniform(0,1.6)",
                                                    class = "Intercept"),
                                          set_prior("normal(0,.5)", class = "b")))
  ### Set all nested formulas for model comparisons
  formulas.trial_parts.per_fstlst <- list(ArcSin ~ 1 +
                                            (1 | Participant) +
                                            (1 | Stimulus) +
                                            (1 | StimLabel),
                                          ArcSin ~ 1 + FstLst +
                                            (1 + FstLst | Participant) +
                                            (1 | Stimulus) +
                                            (1 | StimLabel),
                                          ArcSin ~ 1 + FstLst + CurrentObject +
                                            (1 + FstLst + CurrentObject | Participant) +
                                            (1 | Stimulus) +
                                            (1 | StimLabel),
                                          ArcSin ~ 1 + FstLst + CurrentObject + Condition +
                                            (1 + FstLst + CurrentObject | Participant) +
                                            (1 | Stimulus) +
                                            (1 | StimLabel),
                                          ArcSin ~ 1 + FstLst + CurrentObject + Condition +
                                            FstLst:CurrentObject +
                                            (1 + FstLst*CurrentObject | Participant) +
                                            (1 | Stimulus) +
                                            (1 | StimLabel),
                                          ArcSin ~ 1 + FstLst + CurrentObject + Condition +
                                            FstLst:CurrentObject + FstLst:Condition +
                                            (1 + FstLst*CurrentObject | Participant) +
                                            (1 | Stimulus) +
                                            (1 | StimLabel),
                                          ArcSin ~ 1 + FstLst + CurrentObject + Condition +
                                            FstLst:CurrentObject + FstLst:Condition +
                                            CurrentObject:Condition +
                                            (1 + FstLst*CurrentObject | Participant) +
                                            (1 | Stimulus) +
                                            (1 | StimLabel),
                                          ArcSin ~ 1 + FstLst + CurrentObject + Condition +
                                            FstLst:CurrentObject + FstLst:Condition +
                                            CurrentObject:Condition +
                                            FstLst:CurrentObject:Condition +
                                            (1 + FstLst*CurrentObject | Participant) +
                                            (1 | Stimulus) +
                                            (1 | StimLabel))
  ### Get brms results
  brms.results <- bayes_factor.brm_fixef(formulas.trial_parts.per_fstlst,
                                         trial_parts.per_fstlst,
                                         priors.trial_parts.per_fstlst,
                                         controls = list(adapt_delta = .999999999999999,
                                                         max_treedepth = 15))
  trial_parts.per_fstlst.brms.models <- brms.results[[1]]
  trial_parts.per_fstlst.brms.bayes_factors <- brms.results[[2]]
  trial_parts.time <- proc.time() - t
  ## Save all the results
  saveRDS(trial_parts.per_fstlst.lmer.model, paste0(save_path, "FstLst_lmerModel.rds"))
  saveRDS(trial_parts.per_fstlst.lmer.anova, paste0(save_path, "FstLst_lmerAnova.rds"))
  lapply(seq_along(trial_parts.per_fstlst.brms.models),
         function(i){
           saveRDS(trial_parts.per_fstlst.brms.models[[i]],
                   paste0(save_path, "FstLst_brmsModel", i, ".rds"))
         })
  saveRDS(trial_parts.per_fstlst.brms.bayes_factors, paste0(save_path, "FstLst_brmsBF.rds"))
}else{
  ## Read all the results
  trial_parts.per_fstlst.lmer.model <- readRDS(paste0(save_path, "FstLst_lmerModel.rds"))
  trial_parts.per_fstlst.lmer.anova <- readRDS(paste0(save_path, "FstLst_lmerAnova.rds"))
  trial_parts.per_fstlst.brms.models <- lapply(1:8,
                                             function(i){
                                               readRDS(paste0(save_path,
                                                              "FstLst_brmsModel", i, ".rds"))
                                             })
  trial_parts.per_fstlst.brms.bayes_factors <- readRDS(paste0(save_path, "FstLst_brmsBF.rds"))
}

generate_plots <- F
if(generate_plots){
  ## Get brm predicted values (using three levels of HPDI to better appreciate data shape)
  trial_parts.raw_predictions <- last(trial_parts.per_fstlst.brms.models) %>%
    predict(summary = F,
            transform = function(x){sin(x)^2}) %>%
    t() %>%
    as_tibble() %>%
    mutate(RowNames = 1:nrow(.))
  trial_parts.predicted <- trial_parts.per_fstlst %>%
    mutate(RowNames = 1:nrow(.)) %>%
    select(CurrentObject, FstLst, Condition, RowNames) %>%
    inner_join(trial_parts.raw_predictions) %>%
    select(-RowNames) %>%
    gather(key = Sample, value = Predicted, -c(CurrentObject, FstLst, Condition))
  trial_parts.predicted.hpdi.97 <- trial_parts.predicted %>%
    select(-Sample) %>%
    split(list(.$CurrentObject, .$FstLst, .$Condition)) %>%
    lapply(function(df){
      hpdi <- as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.97)
      df.summary <- df %>%
        group_by(CurrentObject, FstLst, Condition) %>%
        summarise(Mean = mean(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows()
  trial_parts.predicted.hpdi.89 <- trial_parts.predicted %>%
    select(-Sample) %>%
    split(list(.$CurrentObject, .$FstLst, .$Condition)) %>%
    lapply(function(df){
      hpdi <- as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.89)
      df.summary <- df %>%
        group_by(CurrentObject, FstLst, Condition) %>%
        summarise(Mean = mean(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows()
  trial_parts.predicted.hpdi.67 <- trial_parts.predicted %>%
    select(-Sample) %>%
    split(list(.$CurrentObject, .$FstLst, .$Condition)) %>%
    lapply(function(df){
      hpdi <- as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.67)
      df.summary <- df %>%
        group_by(CurrentObject, FstLst, Condition) %>%
        summarise(Mean = mean(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows()
  ## Plot raincloud + predicted mean&HPDIs per FstLst
  trial_parts.per_fstlst.plot <- ggplot(trial_parts.per_fstlst,
                                      aes(x = Condition, y = Prop,
                                          colour = Condition,
                                          fill = Condition)) +
    theme_bw() + ylab("Looking to Tail (Prop)") +
    theme(legend.position = "top",
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    coord_flip() + facet_grid(CurrentObject~FstLst) +
    geom_flat_violin(position = position_nudge(x = .2), colour = "black", alpha = .5) +
    geom_point(position = position_jitter(width = .15),
               size = 1, alpha = .6,
               show.legend = F) +
    geom_boxplot(width = .1, alpha = .3, outlier.shape = NA, colour = "black",
                 show.legend = F) +
    # geom_pointrange(data = trial_parts.predicted.hpdi.67,
    #                 aes(x = Condition, y = Mean, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = 1.5, size = 1.5,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    # geom_pointrange(data = trial_parts.predicted.hpdi.89,
    #                 aes(x = Condition, y = Mean, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = .5, size = 1,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    # geom_pointrange(data = trial_parts.predicted.hpdi.97,
    #                 aes(x = Condition, y = Mean, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = .5, size = .5,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")
  ## Save plot
  ggsave(paste0(save_path, "FstLst_data.pdf"),
         trial_parts.per_fstlst.plot,
         width = 5.5, height = 7)
}

# TRAINING LT ANALYSIS: TAIL LOOKING TIME COURSE ===================================================
save_path <- "../results/adults_2f/PropTail/TimeCourse_"
# DATA PREPARATION
prop_tail.time_course.per_fstlst <- LT.clean %>%
  drop_na(FstLst) %>%
  subset_by_window(window_start_time = -1000, rezero = F) %>%
  make_time_sequence_data(time_bin_size = 100,
                          aois = c("Tail"),
                          predictor_columns=c("Condition",
                                              "FstLst"),
                          summarize_by = "Participant")

# BOOTSTRAPPED CLUSTER-BASED PERMUTATION ANALYSIS
run_model <- F # Running the model takes around 2 minutes on a 4.40GHz 12-core
if(run_model){
  t <- proc.time()
  ## Determine threshold based on alpha = .05 two-tailed
  num_sub = length(unique((prop_tail.time_course.per_fstlst$Participant)))
  threshold_t = qt(p = 1 - .05/2,
                   df = num_sub-1)
  ## Determine clusters
  prop_tail.time_cluster.per_fstlst <- prop_tail.time_course.per_fstlst %>%
    split(.$FstLst) %>%
    lapply(make_time_cluster_data,
           predictor_column = "Condition",
           treatment_level = "No Label",
           aoi = "Tail",
           test = "t.test",
           threshold = threshold_t)
  ## Run the analysis
  prop_tail.time_cluster.per_fstlst.analysis <- prop_tail.time_cluster.per_fstlst %>%
    lapply(analyze_time_clusters,
           within_subj = F,
           parallel = T)
  prop_tail.bcbp.time <- proc.time() - t
  ## Save results
  saveRDS(prop_tail.time_cluster.per_fstlst, paste0(save_path, "FstLst_bcbpClusters.rds"))
  saveRDS(prop_tail.time_cluster.per_fstlst.analysis, paste0(save_path, "FstLst_bcbpAnalysis.rds"))
}else{
  ## Read the results
  prop_tail.time_cluster.per_fstlst <- readRDS(paste0(save_path, "FstLst_bcbpClusters.rds"))
  prop_tail.time_cluster.analysis <- readRDS(paste0(save_path, "FstLst_bcbpAnalysis.rds"))
}
# PLOTTING
## Plot prop_tail time-course for first block and last block
generate_plots <- F
if(generate_plots){
  intercept <- tibble(Part = c(rep("First Block", 2), rep("Last Block", 2)),
                      x_int = c(0, 2000, 0, 2000))
  prop_tail.time_course.per_fstlst.plot <- ggplot(prop_tail.time_course.per_fstlst,
                                                  aes(x = Time, y=Prop,
                                                      colour=Condition,
                                                      fill=Condition)) +
    xlab('Time in Trial') + ylab("Looking to Tail (Prop)") + theme_bw() +
    theme(legend.pos = "top",
          axis.text.x = element_text(angle=45, vjust=1, hjust = 1)) +
    facet_grid(.~FstLst) + ylim(0,1) +
    scale_x_continuous(breaks = c(-1000, 0, 1000, 2000, 3000)) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    geom_vline(data = intercept, aes(xintercept = x_int), linetype = "62", alpha = .5) +
    stat_summary(fun.y='mean', geom='line', linetype = '61') +
    stat_summary(fun.data="mean_se", geom='ribbon', alpha= .33, colour=NA) +
    geom_hline(yintercept = .5)
  ggsave(paste0(save_path, "FstLst_data.pdf"),
         plot = prop_tail.time_course.per_fstlst.plot,
         width = 5.5, height = 3, dpi = 600)
}

# TRAINING LT ANALYSIS: AOI SWITCHES ===============================================================
save_path <- "../results/adults_2f/FamSwitches/FstLst_"
# Prepare dataset
fam_switches.fstlst <- LT.clean %>%
  drop_na(Tail, FstLst) %>%
  group_by(Participant, TrialId) %>%
  summarise(Switches = sum(Tail != lag(Tail), na.rm = T), # Count switches per trial per participant
            # Keep columns for analysis (only one value per trial per participant)
            FstLst = first(FstLst),
            Condition = first(Condition)) %>%
  ungroup()
# Testing Switches ~ Condition*FstLst
run_model <- F # Running the models takes around XXX minutes on a 4.40GHz 12-core
if(run_model){
  t <- proc.time()
  ## Run (g)lmer
  fam_switches.per_fstlst.glmer.model <- glmer(Switches ~ FstLst*Condition +
                                                 (1 + FstLst | Participant),
                                               data = fam_switches.fstlst,
                                               family = poisson())
  # Current p-values from summary may not be the best. Do something else?
  ## Run brms
  ### Set priors for models other than intercept-only
  priors.fam_switches.per_fstlst <- list(NULL,
                                         set_prior("normal(0,.5)", class = "b"))
  ### Set all nested formulas for model comparisons
  formulas.fam_switches.per_fstlst <- list(Switches ~ 1 +
                                             (1 | Participant),
                                           Switches ~ FstLst +
                                             (1 + FstLst | Participant),
                                           Switches ~ FstLst + Condition +
                                             (1 + FstLst | Participant),
                                           Switches ~ FstLst + Condition +
                                             FstLst:Condition +
                                             (1 + FstLst | Participant))
  ### Get brms results
  brms.results <- bayes_factor.brm_fixef(formulas.fam_switches.per_fstlst,
                                         fam_switches.fstlst,
                                         priors.fam_switches.per_fstlst,
                                         family = poisson())
                                         # iter = 4000,
                                         # controls = list(adapt_delta = .95))

  fam_switches.per_fstlst.brms.models <- brms.results[[1]]
  fam_switches.per_fstlst.brms.bayes_factors <- brms.results[[2]]
  fam_switches.time <- proc.time() - t
  ## Save all the results
  saveRDS(fam_switches.per_fstlst.glmer.model, paste0(save_path, "glmerModel.rds"))
  lapply(seq_along(fam_switches.per_fstlst.brms.models),
         function(i){
           saveRDS(fam_switches.per_fstlst.brms.models[[i]],
                   paste0(save_path, "brmsModel", i, ".rds"))
         })
  saveRDS(fam_switches.per_fstlst.brms.bayes_factors, paste0(save_path, "brmsBF.rds"))
}else{
  ## Read all the results
  fam_switches.per_fstlst.glmer.model <- readRDS(paste0(save_path, "glmerModel.rds"))
  fam_switches.per_fstlst.brms.models <- lapply(1:4,
                                                function(i){
                                                  readRDS(paste0(save_path,
                                                                 "brmsModel", i, ".rds"))
                                                })
  fam_switches.per_fstlst.brms.bayes_factors <- readRDS(paste0(save_path, "brmsBF.rds"))
}
# BEHAVIOURAL ANALYSIS: BLOCKS PER PARTICIPANTS ====================================================
save_path <- "../results/adults_2f/nBlocks/"
# Get number of blocks to learning per participant
blocks_per_part <- behaviour %>%
  select(c(Participant, Condition, NBlocks)) %>%
  unique()
# Test NBlocks ~ Condition
run_model <- F
if(run_model){
  t <- proc.time()
  ## STB stats
  blocks_per_part.normality <- ad.test(blocks_per_part$NBlocks)
  blocks_per_part.wilcox <- wilcox.test(NBlocks ~ Condition,
                                        data = blocks_per_part)
  ## brm
  ### Set priors for models other than intercept-only
  priors.blocks_per_part <- list(NULL,
                                 set_prior("normal(0,.5)", class = "b"))
  ### Set all nested formulas for model comparisons
  formulas.blocks_per_part <- list(NBlocks ~ 1,
                                   NBlocks ~ 1 + Condition)
  ### Get brms results
  brms.results <- bayes_factor.brm_fixef(formulas.blocks_per_part,
                                         blocks_per_part,
                                         priors.blocks_per_part,
                                         family = poisson())
  blocks_per_part.brms.models <- brms.results[[1]]
  blocks_per_part.brms.bayes_factors <- brms.results[[2]]
  blocks_per_part.time <- proc.time() - t
  ## Save all the results
  saveRDS(blocks_per_part.normality, paste0(save_path, "normality.rds"))
  saveRDS(blocks_per_part.wilcox, paste0(save_path, "wilcox.rds"))
  lapply(seq_along(blocks_per_part.brms.models),
         function(i){
           saveRDS(blocks_per_part.brms.models[[i]],
                   paste0(save_path, "brmsModel", i, ".rds"))
         })
  saveRDS(blocks_per_part.brms.bayes_factors, paste0(save_path, "brmsBF.rds"))
}else{
  ## Read all the results
  blocks_per_part.normality <- readRDS(paste0(save_path, "normality.rds"))
  blocks_per_part.wilcox <- readRDS(paste0(save_path, "wilcox.rds"))
  blocks_per_part.brms.models <- lapply(1:2,
                                        function(i){
                                          readRDS(paste0(save_path, "brmsModel", i, ".rds"))
                                        })
  blocks_per_part.brms.bayes_factors <- readRDS(paste0(save_path, "brmsBF.rds"))
}

generate_plots <- F
if(generate_plots){
  ## Get brm predicted values (using three levels of HPDI to better appreciate data shape)
  blocks_per_part.raw_predictions <- last(blocks_per_part.brms.models) %>%
    predict(summary = F) %>%
    t() %>%
    as_tibble() %>%
    mutate(RowNames = 1:nrow(.))
  blocks_per_part.predicted <- blocks_per_part %>%
    mutate(RowNames = 1:nrow(.)) %>%
    select(Condition, RowNames) %>%
    inner_join(blocks_per_part.raw_predictions) %>%
    select(-RowNames) %>%
    gather(key = Sample, value = Predicted, -c(Condition))
  blocks_per_part.predicted.hpdi.97 <- blocks_per_part.predicted %>%
    select(-Sample) %>%
    split(list(.$Condition)) %>%
    lapply(function(df){
      hpdi <- as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.97)
      df.summary <- df %>%
        group_by(Condition) %>%
        summarise(Mean = mean(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows()
  blocks_per_part.predicted.hpdi.89 <- blocks_per_part.predicted %>%
    select(-Sample) %>%
    split(list(.$Condition)) %>%
    lapply(function(df){
      hpdi <- as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.89)
      df.summary <- df %>%
        group_by(Condition) %>%
        summarise(Mean = mean(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows()
  blocks_per_part.predicted.hpdi.67 <- blocks_per_part.predicted %>%
    select(-Sample) %>%
    split(list(.$Condition)) %>%
    lapply(function(df){
      hpdi <- as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.67)
      df.summary <- df %>%
        group_by(Condition) %>%
        summarise(Mean = mean(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows()
  ## Plot raincloud + predicted mean&HPDIs per FstLst
  blocks_per_part.plot <- ggplot(blocks_per_part,
                                        aes(x = Condition, y = NBlocks,
                                            colour = Condition,
                                            fill = Condition)) +
    theme_bw() + ylab("Number of Blocks to Learning") +
    theme(legend.position = "top",
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    coord_flip() +
    geom_flat_violin(position = position_nudge(x = .2), colour = "black", alpha = .5) +
    geom_point(position = position_jitter(width = .15, height = .05),
               size = 1, alpha = .6,
               show.legend = F) +
    geom_boxplot(width = .1, alpha = .3, outlier.shape = NA, colour = "black",
                 show.legend = F) +
    # geom_pointrange(data = blocks_per_part.predicted.hpdi.67,
    #                 aes(x = Condition, y = Mean, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = 1.5, size = 1.5,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    # geom_pointrange(data = blocks_per_part.predicted.hpdi.89,
    #                 aes(x = Condition, y = Mean, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = .5, size = 1,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    # geom_pointrange(data = blocks_per_part.predicted.hpdi.97,
    #                 aes(x = Condition, y = Mean, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = .5, size = .5,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")
  ## Save plot
  ggsave(paste0(save_path, "NBlocks_data.pdf"),
         blocks_per_part.plot,
         width = 4, height = 3, dpi = 600)
}
# BEHAVIOURAL ANALYSIS: ACCURACY ~ CONDITION*RT) ===================================================
save_path <- "../results/adults_2f/ACC/"
# Get datasets for training and test
behaviour.training <- behaviour %>%
  subset(Phase == "Familiarisation") %>%
  mutate(BlockZero = Block - 1) # Useful for models
behaviour.test <- behaviour %>%
  subset(Phase == "Test")
# Test ACC ~ Condition * RT
run_model <- T
if(run_model){
  ## Run binomial glmer
  ### During training
  ACC_by_cond_by_RT_by_block.training.glmer <- glmer(ACC ~ Condition*zLogRT*BlockZero +
                                                       (1 + zLogRT + BlockZero | Participant) +
                                                       (1 | Stimulus) +
                                                       (1 | StimLabel),
                                                     family = binomial,
                                                     control = glmerControl(optimizer = "bobyqa"),
                                                     data = behaviour.training)
  ### At test
  ACC_by_cond_by_RT_by_block.test.glmer <- glmer(ACC ~ Condition*zLogRT +
                                                   (1 + zLogRT | Participant),
                                                 family = binomial,
                                                 data = behaviour.test)
  ## Save results
  saveRDS(ACC_by_cond_by_RT_by_block.training.glmer, paste0(save_path, "Training.rds"))
  saveRDS(ACC_by_cond_by_RT_by_block.test.glmer, paste0(save_path, "Test.rds"))
}else{
  ACC_by_cond_by_RT_by_block.training.glmer <- readRDS(paste0(save_path, "Training.rds"))
  ACC_by_cond_by_RT_by_block.test.glmer <- readRDS(paste0(save_path, "Test.rds"))
}

generate_plots <- T
if(generate_plots){
  # Prepare and plot data
  ## During training
  ACC_by_diag_by_RT.training <- behaviour.training %>%
    drop_na(FstLst) %>%
    group_by(Participant, Condition, FstLst) %>%
    summarise(Accuracy = sum(ACC == 1)/n())
  ACC_by_diag_by_RT.training.plot <- ggplot(ACC_by_diag_by_RT.training,
                                            aes(x = Condition,
                                                y = Accuracy,
                                                colour = Condition,
                                                fill = Condition)) +
    ylim(0,1) + theme_bw() +
    theme(legend.pos = "top",
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    coord_flip() + facet_grid(.~FstLst) +
    #scale_fill_discrete(labels = c("Label", "No Label")) +
    geom_flat_violin(position = position_nudge(x = .2), colour = "black", alpha = .5) +
    geom_point(position = position_jitter(width = .15),
               size = 1, alpha = .6,
               show.legend = F) +
    geom_boxplot(width = .1, alpha = .3, outlier.shape = NA, colour = "black",
                 show.legend = F) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")
  ggsave(paste0(save_path, "RTbyFstLst_data.pdf"), plot = ACC_by_diag_by_RT.training.plot,
         width = 5.5, height = 3)
  ## At test
  # ACC_by_diag_by_RT.test <- behaviour.test %>%
  #   group_by(Participant, Condition) %>%
  #   summarise(Accuracy = sum(ACC == 1)/n())
  # ACC_by_diag_by_RT.test.plot <- ggplot(ACC_by_diag_by_RT.test,
  #                                       aes(x = Condition,
  #                                           y = Accuracy,
  #                                           fill = Condition)) +
  #   ylim(0,1) + theme_apa(legend.pos = "bottomright") +
  #   scale_x_discrete(labels = c("Label", "No Label")) +
  #   geom_violin() +
  #   geom_boxplot(alpha = 0, outlier.alpha = 1,
  #                width = .15, position = position_dodge(.9))
  # ggsave("../results/adults_2f/ACCbyRT_test.pdf", plot = ACC_by_diag_by_RT.test.plot,
  #        width = 3.5, height = 2.7)
}
