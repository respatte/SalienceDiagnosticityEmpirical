# LIBRARY IMPORTS ==================================================================================
library(lme4)
library(lmerTest)
library(emmeans)
library(brms)
library(coda)
library(nortest)
library(tidyverse)
library(broom)
library(eyetrackingR)
library(RColorBrewer)

source("Routines.R")
source("StatTools.R")
source("geom_flat_violin.R")

# GATHER DATA ======================================================================================
d <- LT_data.gather("adults_3f")
# Unload snow packages so that parallel works for brms
detach("package:doSNOW")
detach("package:snow")
# Get behavioural and LT data, excluding outliers in terms of
# number of blocks before learning (graphically, from boxplot)
behaviour <- d[[2]] %>%
  subset((Condition == "Label" & NBlocks < 21) |
           (Condition == "No Label" & NBlocks < 10))
LT.clean <- d[[4]] %>%
  subset((Condition == "Label" & NBlocks < 21) |
           (Condition == "No Label" & NBlocks < 10)) %>%
  make_eyetrackingr_data(participant_column = "Participant",
                         trial_column = "TrialId",
                         time_column = "TimeStamp",
                         trackloss_column = "TrackLoss",
                         aoi_columns = c("Head","Tail", "Feet"),
                         treat_non_aoi_looks_as_missing = T) %>%
  subset_by_window(window_start_time = -1000, rezero = F)

# TRAINING LT ANALYSIS: PROP TAIL LOOKING BY PARTICIPANT ===========================================
save_path <- "../results/adults_3f/PropTail/TrialAverage_"
# DATA PREPARATION
prop_tail.per_fstlst <- LT.clean %>%
  drop_na(FstLst) %>%
  make_time_window_data(aois=c("Tail","Feet","Head"),
                        predictor_columns=c("Condition",
                                            "FstLst",
                                            "Stimulus",
                                            "StimLabel"))

# Testing ArcSin ~ FstLst*AOI*Condition
run_model <- T # Running the models takes around XX minutes on a 4.40GHz 12-core
if(run_model){
  ## Run lmer (Sampling Theory Based)
  t <- proc.time()
  prop_tail.per_fstlst.lmer.model <- lmer(ArcSin ~ FstLst*AOI*Condition +
                                            (1 + AOI*FstLst | Participant) +
                                            (1 + AOI | Stimulus) +
                                            (1 + AOI | StimLabel),
                                          data = prop_tail.per_fstlst)
  prop_tail.per_fstlst.lmer.anova <- anova(prop_tail.per_fstlst.lmer.model, type = 1)
  ## Run brms (Bayesian)
  ### Set priors for models other than intercept-only
  priors.prop_tail.per_fstlst <- list(set_prior("normal(.8,.5)",
                                                class = "Intercept"),
                                      c(set_prior("normal(.8,.5)",
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
                                        ArcSin ~ 1 + FstLst + AOI +
                                          (1 + FstLst + AOI | Participant) +
                                          (1 + AOI | Stimulus) +
                                          (1 + AOI | StimLabel),
                                        ArcSin ~ 1 + FstLst + AOI + Condition +
                                          (1 + FstLst + AOI | Participant) +
                                          (1 + AOI | Stimulus) +
                                          (1 + AOI | StimLabel),
                                        ArcSin ~ 1 + FstLst + Condition +
                                          FstLst:AOI +
                                          (1 + FstLst*AOI | Participant) +
                                          (1 + AOI | Stimulus) +
                                          (1 + AOI | StimLabel),
                                        ArcSin ~ 1 + FstLst + Condition +
                                          FstLst:AOI + FstLst:Condition +
                                          (1 + FstLst*AOI | Participant) +
                                          (1 + AOI | Stimulus) +
                                          (1 + AOI | StimLabel),
                                        ArcSin ~ 1 + FstLst + Condition +
                                          FstLst:AOI + FstLst:Condition + AOI:Condition +
                                          (1 + FstLst*AOI | Participant) +
                                          (1 + AOI | Stimulus) +
                                          (1 + AOI | StimLabel),
                                        ArcSin ~ 1 + FstLst + Condition +
                                          FstLst:AOI + FstLst:Condition + AOI:Condition +
                                          FstLst:AOI:Condition +
                                          (1 + FstLst*AOI | Participant) +
                                          (1 + AOI | Stimulus) +
                                          (1 + AOI | StimLabel))
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
    coord_flip() + facet_grid(AOI~FstLst) +
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
         width = 5.5, height = 5.5, dpi = 600)
}

# TRAINING LT ANALYSIS: TAIL LOOKING TIME COURSE ===================================================
save_path <- "../results/adults_3f/PropTail/TimeCourse_"
# DATA PREPARATION
prop_tail.time_course.per_fstlst <- LT.clean %>%
  drop_na(FstLst) %>%
  subset_by_window(window_start_time = -1000, rezero = F) %>%
  make_time_sequence_data(time_bin_size = 100,
                          aois = c("Tail", "Feet", "Head"),
                          predictor_columns=c("Condition",
                                              "FstLst"),
                          summarize_by = "Participant")

# BOOTSTRAPPED CLUSTER-BASED PERMUTATION ANALYSIS
run_model <- T # Running the model takes around 2 minutes on a 4.40GHz 12-core
if(run_model){
  t <- proc.time()
  ## Determine threshold based on alpha = .05 two-tailed
  num_sub = length(unique((prop_tail.time_course.per_fstlst$Participant)))
  threshold_t = qt(p = 1 - .05/2,
                   df = num_sub-1)
  ## Determine clusters
  prop_tail.time_cluster.per_fstlst <- prop_tail.time_course.per_fstlst %>%
    split(list(.$FstLst, .$AOI)) %>%
    lapply(function(df){
      return(make_time_cluster_data(df,
                                    predictor_column = "Condition",
                                    treatment_level = "No Label",
                                    aoi = first(df$AOI),
                                    test = "t.test",
                                    threshold = threshold_t))
    })
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
generate_plots <- T
if(generate_plots){
  intercept <- tibble(Part = c(rep("First Block", 2), rep("Last Block", 2)),
                      x_int = c(0, 2000, 0, 2000))
  fstlst_by_aoi <- names(prop_tail.time_cluster.per_fstlst.analysis)
  prop_tail.time_course.plot.clusters <-  fstlst_by_aoi %>%
    lapply(function(name){
      tmp <- strsplit(name, "\\.")[[1]]
      fstlst <- tmp[[1]]
      aoi <- tmp[[2]]
      df <- prop_tail.time_cluster.per_fstlst.analysis[[name]]$clusters %>%
        mutate(FstLst = fstlst,
               AOI = aoi) %>%
        subset(Probability < .05)
    }) %>%
    bind_rows()
  prop_tail.time_course.per_fstlst.plot <- ggplot(prop_tail.time_course.per_fstlst,
                                                  aes(x = Time, y=Prop,
                                                      colour=Condition,
                                                      fill=Condition)) +
    xlab('Time in Trial') + ylab("Looking to Tail (Prop)") + theme_bw() +
    theme(legend.pos = "top",
          axis.text.x = element_text(angle=45, vjust=1, hjust = 1)) +
    facet_grid(AOI~FstLst) + ylim(0,1) +
    scale_x_continuous(breaks = c(-1000, 0, 1000, 2000, 3000)) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    geom_vline(data = intercept, aes(xintercept = x_int), linetype = "62", alpha = .5) +
    stat_summary(fun.y='mean', geom='line', linetype = '61') +
    stat_summary(fun.data="mean_se", geom='ribbon', alpha= .33, colour=NA) +
    geom_rect(data = prop_tail.time_course.plot.clusters,
              inherit.aes = F,
              aes(xmin = StartTime, xmax = EndTime,
                  ymin = 0, ymax = 1),
              alpha = 0.5,
              fill = brewer.pal(3, "Dark2")[[3]]) +
    geom_hline(yintercept = 1/3)
  ggsave(paste0(save_path, "FstLst_data.pdf"),
         plot = prop_tail.time_course.per_fstlst.plot,
         width = 5.5, height = 5.5, dpi = 600)
}

# LOOKING TIME ANALYSIS: TIME COURSE ===============================================================
# DATA PREPARATION
LT.time_course_aois.first_last <- LT.clean %>%
  make_time_sequence_data(time_bin_size = 50,
                          aois = c("Head","Tail","Feet"),
                          predictor_columns=c("Condition",
                                              "FstLst",
                                              "ACC",
                                              "Stimulus",
                                              "StimLabel",
                                              "Diagnostic")) %>%
  drop_na(FstLst)
# GROWTH CURVE ANALYSIS
run_model = F # Running the model takes around 12h30 on a [check office CPU specs]
if(run_model){
  ## Run model
  # Analysing proportions => main effect of Condition or Part nonsensical,
  # we can only expect differences between AOIs, and between AOIs on different levels
  t <- proc.time()
  LT.time_course_aois.GCA <- lmer(ArcSin ~ (AOI + Condition:AOI + AOI:FstLst +
                                              AOI:Condition:FstLst)*
                                    (ot1 + ot2 + ot3 + ot4 + ot5 + ot6 + ot7) +
                                    (1 + AOI + FstLst + ot1 + ot2 + ot3 +
                                       ot4 + ot5 + ot6 + ot7 | Participant) +
                                    (1 + AOI | Stimulus) +
                                    (1 + AOI | StimLabel),
                                  data = LT.time_course_aois.first_last, REML = F,
                                  control = lmerControl(optCtrl = list(maxfun = 100000)))
  ## Run ANOVA for the model effects
  LT.time_course_aois.GCA.anova <- anova(LT.time_course_aois.GCA, type = 1)
  gca.time <- proc.time() - t
  ## Save model and ANOVA
  saveRDS(LT.time_course_aois.GCA, file = "../results/adults_3f/GCA.rds")
  saveRDS(LT.time_course_aois.GCA.anova, file = "../results/adults_3f/GCA_anova.rds")
}else{
  LT.time_course_aois.GCA <- readRDS("../results/adults_3f/GCA.rds")
  LT.time_course_aois.GCA.anova <- readRDS("../results/adults_3f/GCA_anova.rds")
}
# BOOTSTRAPPED CLUSTER-BASED PERMUTATION ANALYSIS
# needs fixing to test each AOI separately (or together?)
run_model <- F # Running the model takes around 24h30 on a [check office CPU specs]
if(run_model){
  t <- proc.time()
  ## Determine clusters
  LT.time_cluster_aois.first_last <- LT.time_course_aois.first_last %>%
    split(list(.$FstLst, .$AOI)) %>%
    lapply(make_time_cluster_data,
           predictor_column = "Condition",
           treatment_level = "NoLabel",
           test = "lmer",
           threshold = 1.5,
           formula = ArcSin ~ Condition +
             (1 | Participant) +
             (1 | Stimulus))
  ## Run the analysis
  LT.time_cluster_aois.first_last.analysis <- LT.time_cluster_aois.first_last %>%
    lapply(analyze_time_clusters, within_subj = T, parallel = T)
  bcbp.time <- proc.time() - t
  ## Save results
  saveRDS(LT.time_cluster_aois.first_last,
          "../results/adults_3f/BCBP_clusters.rds")
  saveRDS(LT.time_cluster_aois.first_last.analysis,
          "../results/adults_3f/BCBP_analysis.rds")
}else{
  ## Read the results
  LT.time_cluster_aois.first_last <- readRDS("../results/adults_3f/BCBP_clusters.rds")
  LT.time_cluster_aois.first_last.analysis <- readRDS("../results/adults_3f/BCBP_analysis.rds")
}
# PLOTTING
# Plotting eye-tracking data and GCA predictions for all AOIs, for first block, last block, and test
intercept <- tibble(FstLst = c(rep("First Block", 2), rep("Last Block", 2)),
                    x_int = c(0, 2000, 0, 2000))
LT.clean.time_course.first_last.plot <- ggplot(LT.time_course_aois.first_last,
                                               aes(x = Time, y=Prop,
                                                   colour=Condition,
                                                   fill=Condition)) +
  xlab('Time in Trial') + ylab("Looking to AOI (Prop)") + theme_apa(legend.pos = "top") +
  scale_colour_discrete(labels = c("Label", "No Label")) +
  scale_fill_discrete(labels = c("Label", "No Label")) +
  facet_grid(AOI~FstLst, scales = "free_x") + ylim(0,1) +
  scale_x_continuous(breaks = c(-1000, 0, 1000, 2000, 3000)) +
  geom_vline(data = intercept, aes(xintercept = x_int), linetype = "62", alpha = .5) +
  stat_summary(fun.y='mean', geom='line', linetype = '61') +
  stat_summary(fun.data=mean_se, geom='ribbon', alpha= .33, colour=NA) +
  geom_hline(yintercept = 1/3)
ggsave("../results/adults_3f/LookingTimeCourseFirstLast.pdf",
       plot = LT.clean.time_course.first_last.plot,
       width = 7, height = 5.7)

# LOOKING TIME ANALYSIS: PROP AOI LOOKING BY PARTICIPANT BY PART ===================================
# DATA PREPARATION
## Prepare dataset with all blocks
LT.prop_aois.per_block <- make_time_window_data(LT.clean,
                                                aois=c("Tail","Feet","Head"),
                                                predictor_columns=c("Condition",
                                                                    "Block",
                                                                    "FstLst",
                                                                    "ACC",
                                                                    "Stimulus",
                                                                    "StimLabel"))
## Comparing first block against last block
LT.prop_aois.first_last <- LT.prop_aois.per_block %>%
  drop_na(FstLst)
# MIXED-EFFECTS MODELS FOR PROP ~ CONDITION*PART*AOI
run_model = T # Running the model takes around 30 seconds on a [check office CPU specs]
if(run_model){
  ## Run and save the model
  #- No main effect of Part since looking at proportions,
  #- so overall they should all equate to one when collapsing
  t <- proc.time()
  LT.prop_aois.first_last.lmer <- lmer(ArcSin ~ FstLst*AOI*Condition - (FstLst*Condition) +
                                         (1 + AOI + FstLst:AOI | Participant) +
                                         (1 + AOI | Stimulus) +
                                         (1 + AOI | StimLabel),
                                       data = LT.prop_aois.first_last)
  saveRDS(LT.prop_aois.first_last.lmer, "../results/adults_3f/PropAOI.rds")
  ## Run and save the ANOVA for model effects
  LT.prop_aois.first_last.lmer.anova <- anova(LT.prop_aois.first_last.lmer, type = 1)
  saveRDS(LT.prop_aois.first_last.lmer.anova, "../results/adults_3f/PropAOI_anova.rds")
  prop_aois.time <- proc.time() - t
}else{
  LT.prop_aois.first_last.lmer <- readRDS("../results/adults_3f/PropAOI.rds")
  LT.prop_aois.first_last.lmer.anova <- readRDS("../results/adults_3f/PropAOI_anova.rds")
}
# PLOTTING
## Plot jitter + mean&se + lines
LT.prop_aois.first_last.plot <- ggplot(LT.prop_aois.first_last,
                                       aes(x = FstLst, y = Prop,
                                           colour = Condition,
                                           fill = Condition)) +
  facet_wrap(~AOI) + theme_apa(legend.pos = "top") + ylab("Looking to AOI (Prop)") +
  scale_colour_discrete(labels = c("Label", "No Label")) +
  scale_fill_discrete(labels = c("Label", "No Label")) +
  geom_point(position = position_jitterdodge(dodge.width = .8,
                                             jitter.width = .2),
             alpha = .25) +
  geom_errorbar(stat = "summary",
                width = .2, colour = "black",
                position = position_dodge(.1)) +
  geom_line(aes(x = FstLst, y = Prop, group = Condition),
            stat = "summary", fun.y = "mean",
            colour = "black",
            position = position_dodge(.1)) +
  geom_point(stat = "summary", fun.y = "mean",
             shape = 18, size = 3,
             position = position_dodge(.1))
ggsave("../results/adults_3f/AOILookingFirstLast.pdf",
       LT.prop_aois.first_last.plot,
       width = 7, height = 5.4)

# BEHAVIOURAL ANALYSIS: PARTICIPANTS AND BLOCKS ====================================================
# Get how many participants for each block, and make a bar plot
behaviour.parts_per_block <- behaviour %>%
  subset(Block > 0) %>%
  group_by(Block, Condition) %>%
  summarise(N_Participants = n_distinct(Participant))
behaviour.parts_per_block.plot <- ggplot(behaviour.parts_per_block,
                                         aes(x = Block, y = N_Participants, fill = Condition)) +
  theme_apa(legend.pos = "topright") + scale_fill_discrete(labels = c("Label", "No Label")) +
  scale_x_continuous(breaks = 1:10) + ylab("Participants") +
  geom_col(position = "dodge")
ggsave("../results/adults_3f/ParticipantsPerBlock.pdf", plot = behaviour.parts_per_block.plot,
       width = 3.5, height = 2.7)
# Get number of blocks to learning per participant, plot a violin
behaviour.blocks_per_part <- behaviour %>%
  select(c(Participant, Condition, NBlocks)) %>%
  unique()
behaviour.blocks_per_part.plot <- ggplot(behaviour.blocks_per_part,
                                         aes(x = Condition, y = NBlocks, fill = Condition)) +
  theme_apa(legend.pos = "topright") + ylab("Blocks") +
  scale_fill_discrete(labels = c("Label", "No Label")) +
  geom_violin() +
  geom_boxplot(alpha = 0, outlier.alpha = 1, width = .15)
ggsave("../results/adults_3f/BlocksPerParticipant.pdf", plot = behaviour.blocks_per_part.plot,
       width = 3.5, height = 2.9)
# Stats for the number of blocks per participant
behaviour.blocks_per_part.normality <- ad.test(behaviour.blocks_per_part$NBlocks)
behaviour.blocks_per_part.freq_test <- wilcox.test(NBlocks ~ Condition,
                                                   data = behaviour.blocks_per_part)

# BEHAVIOURAL ANALYSIS: ACCURACY ~ CONDITION*DIAG*RT) ==============================================
# Get datasets for training and test
behaviour.training <- behaviour %>%
  subset(Phase == "Familiarisation")
behaviour.test <- behaviour %>%
  subset(Phase == "Test")
# Run binomial glmer
## During training
ACC_by_diag_by_RT.training.glmer <- glmer(ACC ~ Condition*Diagnostic*zLogRT +
                                      (1 + Diagnostic + zLogRT | Participant) +
                                      (1 + zLogRT | Stimulus) +
                                      (1 + zLogRT | StimLabel),
                                    family = binomial,
                                    control = glmerControl(optimizer = "bobyqa"),
                                    data = behaviour.training)
## At test !!! NOT CONVERGING !!! ==> All participants at ceiling
##ACC_by_diag_by_RT.test.glmer <- glm(ACC ~ Condition*Diagnostic +
##                                  (1 | Participant),
##                                family = binomial,
##                                data = behaviour.test)
# Prepare and plot data
## During training
ACC_by_diag_by_RT.training <- behaviour.training %>%
  group_by(Participant, Diagnostic, Condition) %>%
  summarise(Accuracy = sum(ACC == 1)/n())
ACC_by_diag_by_RT.training.plot <- ggplot(ACC_by_diag_by_RT.training,
                                    aes(x = Condition,
                                        y = Accuracy,
                                        fill = Diagnostic)) +
  ylim(0,1) + theme_apa(legend.pos = "bottomright") +
  scale_x_discrete(labels = c("Label", "No Label")) +
  geom_violin() +
  geom_boxplot(alpha = 0, outlier.alpha = 1,
               width = .15, position = position_dodge(.9))
ggsave("../results/adults_3f/ACCbyRTbyDiag_training.pdf", plot = ACC_by_diag_by_RT.training.plot,
       width = 3.5, height = 2.7)
## At test
ACC_by_diag_by_RT.test <- behaviour.test %>%
  group_by(Participant, Diagnostic, Condition) %>%
  summarise(Accuracy = sum(ACC == 1)/n())
ACC_by_diag_by_RT.test.plot <- ggplot(ACC_by_diag_by_RT.test,
                                aes(x = Condition,
                                    y = Accuracy,
                                    fill = Diagnostic)) +
  ylim(0,1) + theme_apa(legend.pos = "bottomright") +
  scale_x_discrete(labels = c("Label", "No Label")) +
  geom_violin() +
  geom_boxplot(alpha = 0, outlier.alpha = 1,
               width = .15, position = position_dodge(.9))
ggsave("../results/adults_3f/ACCbyRTbyDiag_test.pdf", plot = ACC_by_diag_by_RT.test.plot,
       width = 3.5, height = 2.7)
