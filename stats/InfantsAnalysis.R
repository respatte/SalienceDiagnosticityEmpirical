# LIBRARY IMPORTS ==================================================================================
library(eyetrackingR)
library(lme4)
library(lmerTest)
library(emmeans)
library(brms)
library(coda)
library(modeest)
library(tidyverse)
library(RColorBrewer)

source("Routines.R")
source("StatTools.R")
source("geom_flat_violin.R")

# GATHER DATA ======================================================================================
  # Load data and run general checks
  d <- LT_data.gather("infants")
  # Unload snow packages so that parallel works for brms
  detach("package:doSNOW")
  detach("package:snow")
  # Remove participants from label condition with wrong setup
  d[[4]] <- d[[4]] %>%
    subset(!(Condition == "Label" & as.numeric(as.character(Participant)) < 72))
  # Check for counterbalancing, gender balancing, and age balancing
  pres_seq <- d[[4]] %>%
    group_by(PresentationSequence, Participant) %>%
    summarise(T = sum(TrackLoss)/n())
  # Keep only one participant when multiple infants saw the same presentation sequence.
  # Current choice: improve gender balance in total and between conditions.
  # remove: P03, P53, P55
  # keep:   P51, P06, P08
  d[[4]] <- d[[4]] %>%
    subset(!(Participant %in% c("03","53","55"))) %>%
    droplevels()
  gender <- d[[4]] %>%
    group_by(Gender, Condition) %>%
    summarise(N = n_distinct(Participant))
  age <- d[[4]] %>%
    group_by(Participant, Condition, Gender) %>%
    summarise(Age = first(Age))
  age_gender_condition <- age %>%
    group_by(Condition, Gender) %>%
    summarise(mu = mean(Age),
              l = min(Age),
              u = max(Age))
  # Checking gaze-data offset (checking on Familiarisation only for simplicity)
  ## Saving heatmaps per participant pre-correction
  generate_plots <- F
  if(generate_plots){
    head_fam <- tibble(AOI_type = factor(c("Reg", "Flip")),
                       xmin = c(1031,1920-1031-450), xmax = c(1031+450,1920-1031),
                       ymin = c(197,197), ymax = c(197+450,197+450))
    tail_fam <- tibble(AOI_type = factor(c("Reg", "Flip")),
                       xmin = c(390,1920-390-450), xmax = c(390+450,1920-390),
                       ymin = c(299,299), ymax = c(299+450,299+450))
    LT.gaze_offset.data.pre <- d[[4]] %>%
      subset(Phase == "Familiarisation")
    for(participant in levels(LT.gaze_offset.data.pre$Participant)){
      LT.gaze_offset.plot.pre <- LT.gaze_offset.data.pre %>%
        subset(Participant == participant) %>%
        ggplot(aes(x=CursorX, y=CursorY)) +
        xlim(c(0, 1920)) + scale_y_reverse(limits = c(1080, 0)) +
        facet_grid(AOI_type~.) +
        geom_bin2d(binwidth = c(20,20)) +
        geom_rect(data = head_fam,
                  inherit.aes = F,
                  aes(xmin = xmin, xmax = xmax,
                      ymin = ymin, ymax = ymax),
                  fill = NA, colour = "red") +
        geom_rect(data = tail_fam,
                  inherit.aes = F,
                  aes(xmin = xmin, xmax = xmax,
                      ymin = ymin, ymax = ymax),
                  fill = NA, colour = "red")
      ggsave(paste0("../results/infants/cleaning/GazeOffsets/PreCorrection/",participant,".png"),
             plot = LT.gaze_offset.plot.pre,
             width = 3.5, height = 2.9)
    }
  }
  ## Apply linear correction to gaze data by participant when necessary
  LT.gaze_offset.data.correction <- d[[4]] %>%
    mutate(CursorX = case_when(Participant == "01" ~ CursorX - 60,
                               Participant == "07" ~ CursorX - 80,
                               Participant == "08" ~ CursorX + 80,
                               Participant == "26" ~ CursorX - 80,
                               Participant == "31" ~ CursorX + 30,
                               Participant == "42" ~ CursorX - 30,
                               Participant == "51" ~ CursorX + 60,
                               Participant == "56" ~ CursorX + 30,
                               Participant == "60" ~ CursorX - 100,
                               Participant == "61" ~ CursorX - 80,
                               Participant == "72" ~ CursorX + 50,
                               Participant == "74" ~ CursorX + 150,
                               Participant == "88" ~ CursorX -30,
                               Participant == "100" ~ CursorX -50,
                               T ~ as.double(CursorX)),
           CursorY = case_when(Participant == "07" ~ CursorY + 100,
                               Participant == "08" ~ CursorY - 50,
                               Participant == "14" ~ CursorY - 80,
                               Participant == "15" ~ CursorY - 150,
                               Participant == "31" ~ CursorY - 100,
                               Participant == "36" ~ CursorY - 50,
                               Participant == "37" ~ CursorY - 50,
                               Participant == "42" ~ CursorY - 50,
                               Participant == "47" ~ CursorY - 50,
                               Participant == "52" ~ CursorY - 50,
                               Participant == "59" ~ CursorY - 150,
                               Participant == "60" ~ CursorY - 150,
                               Participant == "61" ~ CursorY - 100,
                               Participant == "64" ~ CursorY - 100,
                               Participant == "65" ~ CursorY - 50,
                               Participant == "66" ~ CursorY - 50,
                               Participant == "67" ~ CursorY - 100,
                               Participant == "72" ~ CursorY -30,
                               Participant == "74" ~ CursorY -30,
                               Participant == "90" ~ CursorY -50,
                               Participant == "92" ~ CursorY -30,
                               Participant == "100" ~ CursorY -50,
                               T ~ as.double(CursorY))) %>%
    select(-c(Head, Tail, NewTail, OldTail, NewHead, OldHead, Centre, Target, Distractor))
    # Removing AOIs as they need updating
  ## Saving heatmaps per participant post-correction
  generate_plots <- F
  if(generate_plots){
    head_fam <- tibble(AOI_type = factor(c("Reg", "Flip")),
                       xmin = c(1031,1920-1031-450), xmax = c(1031+450,1920-1031),
                       ymin = c(197,197), ymax = c(197+450,197+450))
    tail_fam <- tibble(AOI_type = factor(c("Reg", "Flip")),
                       xmin = c(390,1920-390-450), xmax = c(390+450,1920-390),
                       ymin = c(299,299), ymax = c(299+450,299+450))
    LT.gaze_offset.data.post <- LT.gaze_offset.data.correction %>%
      subset(Phase == "Familiarisation") %>%
      droplevels()
    for(participant in levels(LT.gaze_offset.data.post$Participant)){
      LT.gaze_offset.plot.post <- LT.gaze_offset.data.post %>%
        subset(Participant == participant) %>%
        ggplot(aes(x=CursorX, y=CursorY)) +
        xlim(c(0, 1920)) + scale_y_reverse(limits = c(1080, 0)) +
        facet_grid(AOI_type~.) +
        geom_bin2d(binwidth = c(20,20)) +
        geom_rect(data = head_fam,
                  inherit.aes = F,
                  aes(xmin = xmin, xmax = xmax,
                      ymin = ymin, ymax = ymax),
                  fill = NA, colour = "red") +
        geom_rect(data = tail_fam,
                  inherit.aes = F,
                  aes(xmin = xmin, xmax = xmax,
                      ymin = ymin, ymax = ymax),
                  fill = NA, colour = "red")
      ggsave(paste0("../results/infants/cleaning/GazeOffsets/PostCorrection/",participant,".png"),
             plot = LT.gaze_offset.plot.post,
             width = 3.5, height = 2.9)
    }
  }
  ## Updating AOI hits
  infants.AOIs <- grep("infants", names(AOIs))
  LT.gaze_offset.data.corrected <- LT.gaze_offset.data.correction
  for(AOI in names(AOIs[infants.AOIs])){
    AOI.name <- sub("infants\\.", "", AOI)
    LT.gaze_offset.data.corrected <- LT.gaze_offset.data.corrected %>%
      left_join(AOIs[[AOI]]) %>%
      mutate(!!AOI.name := CursorX>Left & CursorX<Right & CursorY>Top & CursorY<Bottom) %>%
      select(-c(Left, Right, Top, Bottom))
  }
  # Creating general datasets for analysis (separating phases, general window sub-setting)
  ## Familiarisation
  LT.fam <- LT.gaze_offset.data.corrected %>%
    subset(Phase == "Familiarisation") %>%
    mutate_at("PrePost", parse_factor,
              levels = c("Pre Label Onset", "Post Label Onset"),
              include_na = F) %>%
    group_by(Participant) %>%
    mutate(First = min(TrialNum, na.rm = T),
           Last = max(TrialNum, na.rm = T),
           FstLst = case_when(TrialNum <= First + 2 ~ "First Trials",
                              TrialNum >= Last - 2 ~ "Last Trials")) %>%
    # Useful to compare beginning-end of experiment per infant
    select(-c(First, Last)) %>%
    make_eyetrackingr_data(participant_column = "Participant",
                           trial_column = "TrialId",
                           time_column = "TimeStamp",
                           trackloss_column = "TrackLoss",
                           aoi_columns = c("Head","Tail"),
                           treat_non_aoi_looks_as_missing = T) %>%
    subset_by_window(window_start_time = 1500,        # Start after stimulus moving in,
                     window_end_col = "TrialEnd") %>% # and end 3000ms after LabelOnset
    mutate(LabelOnset = LabelOnset - 1500) # Update LabelOnset after window sub-setting
  ## Contrast tests
  LT.test.ctr <- LT.gaze_offset.data.corrected %>%
    subset(Phase == "Test - Contrast") %>%
    mutate(NewFeature = ifelse(ContrastType == "Relative" | (is.na(NewTail) & is.na(NewHead)),
                               NA, (NewTail | NewHead) %in% T),
           OldFeature = ifelse(ContrastType == "Relative" | (is.na(OldTail) & is.na(NewTail)),
                               NA, (OldTail | OldHead) %in% T)) %>%
    make_eyetrackingr_data(participant_column = "Participant",
                           trial_column = "TrialId",
                           time_column = "TimeStamp",
                           trackloss_column = "TrackLoss",
                           aoi_columns = c("NewHead","OldHead","NewTail","OldTail",
                                           "NewFeature", "OldFeature"),
                           # Ignore looks twoards the `centre' AOI`
                           treat_non_aoi_looks_as_missing = T)
  ## Word learning tests
  LT.test.wl <- LT.gaze_offset.data.corrected %>%
    subset(Phase == "Test - Word Learning") %>%
    droplevels() %>%
    mutate_at("PrePost", parse_factor,
              levels = c("Pre Label Onset", "Post Label Onset"),
              include_na = F) %>%
    make_eyetrackingr_data(participant_column = "Participant",
                           trial_column = "TrialId",
                           time_column = "TimeStamp",
                           trackloss_column = "TrackLoss",
                           aoi_columns = c("Target","Distractor"),
                           treat_non_aoi_looks_as_missing = T)
# FAMILIARISATION ANALYSIS: PROP TAIL LOOKING BY TRIAL/BLOCK/FSTLST ================================
save_path <- "../results/infants/PropTail/TrialAverage_"
# Prepare dataset
prop_tail.fstlst <- LT.fam %>%
  drop_na(FstLst) %>%
  subset_by_window(window_start_col = "LabelOnset") %>%
  make_time_window_data(aois=c("Tail"),
                        predictor_columns=c("Condition",
                                            "FstLst",
                                            "Stimulus",
                                            "CategoryName")) %>%
  drop_na(ArcSin)
# Testing Prop ~ FstLst*Condition
run_model <- F # Running the models takes around 4 minutes on a 4.40GHz 12-core
if(run_model){
  t <- proc.time()
  ## Run lmer
  prop_tail.per_fstlst.lmer.model <- lmer(ArcSin ~ FstLst*Condition +
                                            (1 + FstLst | Participant) +
                                            (1 | Stimulus),
                                          data = prop_tail.fstlst)
  prop_tail.per_fstlst.lmer.anova <- anova(prop_tail.per_fstlst.lmer.model, type = 1)
  ## Run brms
  ### Set priors for models other than intercept-only
  priors.prop_tail.per_fstlst <- list(set_prior("uniform(0,1.6)",
                                                class = "Intercept"),
                                      c(set_prior("uniform(0,1.6)",
                                                  class = "Intercept"),
                                        set_prior("normal(0,.5)", class = "b")))
  ### Set all nested formulas for model comparisons
  formulas.prop_tail.per_fstlst <- list(ArcSin ~ 1 +
                                          (1 | Participant) +
                                          (1 | Stimulus),
                                        ArcSin ~ FstLst +
                                          (1 + FstLst | Participant) +
                                          (1 | Stimulus),
                                        ArcSin ~ FstLst + Condition +
                                          (1 + FstLst | Participant) +
                                          (1 | Stimulus),
                                        ArcSin ~ FstLst + Condition +
                                          FstLst:Condition +
                                          (1 + FstLst | Participant) +
                                          (1 | Stimulus))
  ### Get brms results
  brms.results <- bayes_factor.brm_fixef(formulas.prop_tail.per_fstlst,
                                         prop_tail.fstlst,
                                         priors.prop_tail.per_fstlst)
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
  ## Read all the results
  prop_tail.per_fstlst.lmer.model <- readRDS(paste0(save_path, "FstLst_lmerModel.rds"))
  prop_tail.per_fstlst.lmer.anova <- readRDS(paste0(save_path, "FstLst_lmerAnova.rds"))
  prop_tail.per_fstlst.brms.models <- lapply(1:4,
                                             function(i){
                                               readRDS(paste0(save_path,
                                                              "FstLst_brmsModel", i, ".rds"))
                                             })
  prop_tail.per_fstlst.brms.bayes_factors <- readRDS(paste0(save_path, "FstLst_brmsBF.rds"))
}

# Plot jitter + mean&se + lines
generate_plots <- F
if(generate_plots){
  ## Get brm predicted values (using three levels of HPDI to better appreciate data shape)
  prop_tail.raw_predictions <- last(prop_tail.per_fstlst.brms.models) %>%
    predict(summary = F,
            transform = function(x){sin(x)^2}) %>%
    t() %>%
    as_tibble() %>%
    mutate(RowNames = 1:nrow(.))
  prop_tail.predicted <- prop_tail.fstlst %>%
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
        summarise(Mode = mlv(df$Predicted)) %>%
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
        summarise(Mode = mlv(df$Predicted)) %>%
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
        summarise(Mode = mlv(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows()
  ## Plot raincloud + predicted mean&HPDIs per FstLst
  prop_tail.per_fstlst.plot <- ggplot(prop_tail.fstlst,
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
    #                 aes(x = Condition, y = Mode, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = 1.5, size = 1.5,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    # geom_pointrange(data = prop_tail.predicted.hpdi.89,
    #                 aes(x = Condition, y = Mode, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = .5, size = 1,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    # geom_pointrange(data = prop_tail.predicted.hpdi.97,
    #                 aes(x = Condition, y = Mode, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = .5, size = .5,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")
  ## Save plot
  ggsave(paste0(save_path, "FstLst_data.pdf"),
         prop_tail.per_fstlst.plot,
         width = 5.5, height = 3)
}

# FAMILIARISATION ANALYSIS: PROP TAIL LOOKING TIME COURSE BY FSTLST  ===============================
save_path <- "../results/infants/PropTail/TimeCourse_"
# Data preparation
prop_tail.time_course <- LT.fam %>%
  drop_na(FstLst) %>%
  subset_by_window(window_start_col = "LabelOnset") %>%
  make_time_sequence_data(time_bin_size = 100,
                          aois = "Tail",
                          predictor_columns=c("Condition",
                                              "FstLst"),
                          summarize_by = "Participant")
# BOOTSTRAPPED CLUSTER-BASED PERMUTATION ANALYSIS
run_model <- F # Running the model takes around 50 seconds on a 4.40GHz 12-core
if(run_model){
  t <- proc.time()
  ## Determine threshold based on alpha = .05 two-tailed
  num_sub = length(unique((prop_tail.time_course$Participant)))
  threshold_t = qt(p = 1 - .05/2,
                   df = num_sub-1)
  ## Determine clusters
  prop_tail.time_cluster <- prop_tail.time_course %>%
    split(.$FstLst) %>%
    lapply(make_time_cluster_data,
           predictor_column = "Condition",
           treatment_level = "No Label",
           aoi = "Tail",
           test = "t.test",
           threshold = threshold_t)
  ## Run analysis
  prop_tail.time_cluster.analysis <- prop_tail.time_cluster %>%
    lapply(analyze_time_clusters,
           within_subj = F,
           parallel = T)
  prop_tail.time_course.time <- proc.time() - t
  ## Save clusters and analysis
  saveRDS(prop_tail.time_cluster, paste0(save_path, "FstLst_bcbpClusters.rds"))
  saveRDS(prop_tail.time_cluster.analysis, paste0(save_path, "FstLst_bcbpAnalysis.rds"))
}else{
  ## Read the results
  prop_tail.time_cluster <- readRDS(paste0(save_path, "FstLst_bcbpClusters.rds"))
  prop_tail.time_cluster.analysis <- readRDS(paste0(save_path, "FstLst_bcbpAnalysis.rds"))
}

# PLOT
generate_plots <- F
if(generate_plots){
  prop_tail.time_course.plot.clusters <- prop_tail.time_cluster.analysis %>%
    {lapply(seq_along(.),
            function(i){
              df <- prop_tail.time_cluster.analysis[[i]]$clusters %>%
                mutate(FstLst = if(i == 1){"First Trials"}else{"Last Trials"}) %>%
                subset(Probability < .05)
            })} %>%
    bind_rows()
  prop_tail.time_course.plot.per_fstlst <- ggplot(prop_tail.time_course,
                                                  aes(x = Time, y=Prop,
                                                      colour=Condition,
                                                      fill=Condition)) +
    xlab('Time in Trial') + ylab("Looking to Tail (Prop)") + theme_bw() +
    facet_grid(.~FstLst) + ylim(0,1) +
    theme(legend.position = "top",
          axis.text.x = element_text(angle=45, vjust=1, hjust = 1)) +
    stat_summary(fun.y='mean', geom='line', linetype = '61') +
    stat_summary(fun.data=mean_se, geom='ribbon', alpha= .25, colour=NA) +
    # NO SIGNIFICANT CLUSTERS TO PRINT IN prop_tail.time_course.plot.clusters
    # geom_rect(data = prop_tail.time_course.plot.clusters,
    #           inherit.aes = F,
    #           aes(xmin = StartTime, xmax = EndTime,
    #               ymin = 0, ymax = 1),
    #           alpha = 0.5,
    #           fill = brewer.pal(3, "Dark2")[[3]]) +
    geom_hline(yintercept = .5) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")
  ggsave(paste0(save_path, "FstLst_data.pdf"),
         plot = prop_tail.time_course.plot.per_fstlst,
         width = 5.5, height = 2.5)
}

# FAMILIARISATION ANALYSIS: PROP AOI LOOKING PRE/POST LABEL ONSET ==================================
save_path <- "../results/infants/PrePost/TrialAverage_"
# Prepare dataset
prop_tail.pre_post.fstlst <- LT.fam %>%
  drop_na(PrePost, FstLst) %>%
  make_time_window_data(aois=c("Tail"),
                        predictor_columns=c("Condition",
                                            "FstLst",
                                            "PrePost",
                                            "Stimulus",
                                            "CategoryName")) %>%
  drop_na(ArcSin)
# Testing Prop ~ FstLst*PrePost*Condition
run_model <- F # Running the models takes around 8 minutes on a 4.40GHz 12-core
if(run_model){
  t <- proc.time()
  ## Run lmer
  pre_post.per_fstlst.lmer.model <- lmer(ArcSin ~ FstLst*PrePost*Condition +
                                              (1 + FstLst | Participant) +
                                              (1 | Stimulus),
                                            data = prop_tail.pre_post.fstlst)
  pre_post.per_fstlst.lmer.anova <- anova(pre_post.per_fstlst.lmer.model, type = 1)
  ## Run brms
  ### Set priors for models other than intercept-only
  priors.pre_post.per_fstlst <- list(set_prior("uniform(0,1.6)",
                                               class = "Intercept"),
                                     c(set_prior("uniform(0,1.6)",
                                                 class = "Intercept"),
                                       set_prior("normal(0,.5)", class = "b")))
  ### Set all nested formulas for model comparisons
  formulas.pre_post.per_fstlst <- list(ArcSin ~ 1 +
                                         (1 | Participant) +
                                         (1 | Stimulus),
                                       ArcSin ~ 1 + FstLst +
                                         (1 + FstLst | Participant) +
                                         (1 | Stimulus),
                                       ArcSin ~ 1 + FstLst + PrePost +
                                         (1 + FstLst | Participant) +
                                         (1 | Stimulus),
                                       ArcSin ~ 1 + FstLst + PrePost + Condition +
                                         (1 + FstLst | Participant) +
                                         (1 | Stimulus),
                                       ArcSin ~ 1 + FstLst + PrePost + Condition +
                                         FstLst:PrePost +
                                         (1 + FstLst | Participant) +
                                         (1 | Stimulus),
                                       ArcSin ~ 1 + FstLst + PrePost + Condition +
                                         FstLst:PrePost + FstLst:Condition +
                                         (1 + FstLst | Participant) +
                                         (1 | Stimulus),
                                       ArcSin ~ 1 + FstLst + PrePost + Condition +
                                         FstLst:PrePost + FstLst:Condition +
                                         PrePost:Condition +
                                         (1 + FstLst | Participant) +
                                         (1 | Stimulus),
                                       ArcSin ~ 1 + FstLst + PrePost + Condition +
                                         FstLst:PrePost + FstLst:Condition +
                                         PrePost:Condition +
                                         FstLst:PrePost:Condition +
                                         (1 + FstLst | Participant) +
                                         (1 | Stimulus))
  ### Get brms results
  brms.results <- bayes_factor.brm_fixef(formulas.pre_post.per_fstlst,
                                         prop_tail.pre_post.fstlst,
                                         priors.pre_post.per_fstlst)
  pre_post.per_fstlst.brms.models <- brms.results[[1]]
  pre_post.per_fstlst.brms.bayes_factors <- brms.results[[2]]
  pre_post.time <- proc.time() - t
  ## Save all the results
  saveRDS(pre_post.per_fstlst.lmer.model, paste0(save_path, "FstLst_lmerModel.rds"))
  saveRDS(pre_post.per_fstlst.lmer.anova, paste0(save_path, "FstLst_lmerAnova.rds"))
  lapply(seq_along(pre_post.per_fstlst.brms.models),
         function(i){
           saveRDS(pre_post.per_fstlst.brms.models[[i]],
                   paste0(save_path, "FstLst_brmsModel", i, ".rds"))
         })
  saveRDS(pre_post.per_fstlst.brms.bayes_factors, paste0(save_path, "FstLst_brmsBF.rds"))
}else{
  ## Read all the results
  pre_post.per_fstlst.lmer.model <- readRDS(paste0(save_path, "FstLst_lmerModel.rds"))
  pre_post.per_fstlst.lmer.anova <- readRDS(paste0(save_path, "FstLst_lmerAnova.rds"))
  pre_post.per_fstlst.brms.models <- lapply(1:8,
                                             function(i){
                                               readRDS(paste0(save_path,
                                                              "FstLst_brmsModel", i, ".rds"))
                                             })
  pre_post.per_fstlst.brms.BF <- readRDS(paste0(save_path, "FstLst_brmsBF.rds"))
}

# Plot jitter + mean&se + lines
generate_plots <- F
if(generate_plots){
  ## Get brm predicted values (using three levels of HPDI to better appreciate data shape)
  pre_post.raw_predictions <- last(pre_post.per_fstlst.brms.models) %>%
    predict(summary = F,
            transform = function(x){sin(x)^2}) %>%
    t() %>%
    as_tibble() %>%
    mutate(RowNames = 1:nrow(.))
  pre_post.predicted <- prop_tail.pre_post.fstlst %>%
    mutate(RowNames = 1:nrow(.)) %>%
    select(FstLst, PrePost, Condition, RowNames) %>%
    inner_join(pre_post.raw_predictions) %>%
    select(-RowNames) %>%
    gather(key = Sample, value = Predicted, -c(FstLst, PrePost, Condition))
  pre_post.predicted.hpdi.97 <- pre_post.predicted %>%
    select(-Sample) %>%
    split(list(.$FstLst, .$PrePost, .$Condition)) %>%
    lapply(function(df){
      hpdi <- if(length(df$Predicted)>1){
        as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.97)
      }else{
        matrix(ncol = 2, dimnames = list(NULL, c("lower", "upper"))) # When no cases
      }
      df.summary <- df %>%
        group_by(FstLst, PrePost, Condition) %>%
        summarise(Mode = mlv(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows()
  pre_post.predicted.hpdi.89 <- pre_post.predicted %>%
    select(-Sample) %>%
    split(list(.$FstLst, .$PrePost, .$Condition)) %>%
    lapply(function(df){
      hpdi <- if(length(df$Predicted)>1){
        as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.89)
      }else{
        matrix(ncol = 2, dimnames = list(NULL, c("lower", "upper"))) # When no cases
      }
      df.summary <- df %>%
        group_by(FstLst, PrePost, Condition) %>%
        summarise(Mode = mlv(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows()
  pre_post.predicted.hpdi.67 <- pre_post.predicted %>%
    select(-Sample) %>%
    split(list(.$FstLst, .$PrePost, .$Condition)) %>%
    lapply(function(df){
      hpdi <- if(length(df$Predicted)>1){
        as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.67)
      }else{
        matrix(ncol = 2, dimnames = list(NULL, c("lower", "upper"))) # When no cases
      }
      df.summary <- df %>%
        group_by(FstLst, PrePost, Condition) %>%
        summarise(Mode = mlv(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows()
  ## Plot raincloud + predicted mean&sd per FstLst
  pre_post.per_fstlst.plot <- ggplot(prop_tail.pre_post.fstlst,
                                    aes(x = Condition, y = Prop,
                                        colour = Condition,
                                        fill = Condition)) +
    theme_bw() + ylab("Looking to Tail (Prop)") +
    coord_flip() + facet_grid(FstLst~PrePost) + guides(fill = F, colour = F) +
    geom_flat_violin(position = position_nudge(x = .2),
                     colour = "black", alpha = .5, width = .7) +
    geom_point(position = position_jitter(width = .15),
               size = 1, alpha = .6) +
    geom_boxplot(width = .1, alpha = .3, outlier.shape = NA, colour = "black") +
    geom_pointrange(data = pre_post.predicted.hpdi.67,
                    aes(x = Condition, y = Mode, ymin = lb, ymax = ub),
                    colour = brewer.pal(3, "Dark2")[[3]],
                    fatten = 1.5, size = 1.5,
                    position = position_nudge(x = -.23)) +
    geom_pointrange(data = pre_post.predicted.hpdi.89,
                    aes(x = Condition, y = Mode, ymin = lb, ymax = ub),
                    colour = brewer.pal(3, "Dark2")[[3]],
                    fatten = .5, size = 1,
                    position = position_nudge(x = -.23)) +
    geom_pointrange(data = pre_post.predicted.hpdi.97,
                    aes(x = Condition, y = Mode, ymin = lb, ymax = ub),
                    colour = brewer.pal(3, "Dark2")[[3]],
                    fatten = .5, size = .5,
                    position = position_nudge(x = -.23)) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")
  ## Save plot
  ggsave(paste0(save_path, "FstLst_data.pdf"),
         pre_post.per_fstlst.plot,
         width = 7, height = 5)
}

# CONTRAST TEST ANALYSIS: PROP NEW FEATURE BY HEAD/TAIL CONTRAST TEST ==============================
save_path <- "../results/infants/OldNew/TrialAverage_"
# Prepare dataset
new_old <- LT.test.ctr %>%
  subset(ContrastType %in% c("Tail", "Head")) %>%
  make_time_window_data(aois = "NewFeature",
                        predictor_columns = c("Condition",
                                              "ContrastType")) %>%
  drop_na(ArcSin) %>%
  mutate(ChanceArcsin = ArcSin - asin(sqrt(.5)))
## Check for amount of data available
new_old.participants <- new_old %>%
  group_by(Participant, Condition) %>%
  summarise(nTrials = n_distinct(TrialId)) %>%
  group_by(Condition, nTrials) %>%
  summarise(Participants = n_distinct(Participant))
# Testing Prop ~ ContrastType*Condition
run_model <- F # Running the models takes around 5 minutes on a 4.40GHz 12-core
if(run_model){
  t <- proc.time()
  ## Run lmer
  new_old.lmer.model <- lmer(ChanceArcsin ~ ContrastType*Condition +
                                  (1 | Participant),
                                data = new_old)
  new_old.lmer.anova <- anova(new_old.lmer.model, type = 1)
  new_old.lmer.emmeans <- emmeans(new_old.lmer.model, ~ ContrastType | Condition,
                                  options = list(infer = c(T, T),
                                                 null = 0,
                                                 adjust = "bonferroni"))
  ## Run brm
  ### Set priors for models other than intercept-only
  priors.new_old <- list(set_prior("uniform(-.8,.8)",
                                   class = "Intercept"),
                         c(set_prior("uniform(-.8,.8)",
                                     class = "Intercept"),
                           set_prior("normal(0,.5)", class = "b")))
  ### Set all nested formulas for model comparisons
  formulas.new_old <- list(ChanceArcsin ~ 0 +
                             (1 | Participant),
                           ChanceArcsin ~ 1 +
                             (1 | Participant),
                           ChanceArcsin ~ ContrastType +
                             (1 | Participant),
                           ChanceArcsin ~ ContrastType + Condition +
                             (1 | Participant),
                           ChanceArcsin ~ ContrastType + Condition +
                             ContrastType:Condition +
                             (1 | Participant))
  ### Get brms results
  brms.results <- bayes_factor.brm_fixef(formulas.new_old,
                                         new_old,
                                         priors.new_old,
                                         no_intercept = 1)
  new_old.brms.models <- brms.results[[1]]
  new_old.brms.bayes_factors <- brms.results[[2]]
  new_old.brms.emmeans.model <- brm(ChanceArcsin ~ ContrastType + Condition +
                                      ContrastType:Condition +
                                      (1 | Participant),
                                    data = new_old,
                                    prior = last(priors.new_old),
                                    family = gaussian(),
                                    chains = 4, cores = 4, iter = 2000,
                                    sample_prior = "yes")
  new_old.brms.emmeans <- emmeans(new_old.brms.emmeans.model, ~ ContrastType | Condition,
                                  options = list(level = .89))
  new_old.brms.emmeans.bayes_factor <- hypothesis(new_old.brms.emmeans.model,
                                                  c("Intercept > 0",                 # No Label Head
                                                    "Intercept + ContrastTypeTail > 0", # NL Tail
                                                    "Intercept + ConditionLabel > 0",   # Label Head
                                                    paste("Intercept +",                # Label Tail
                                                          "ConditionLabel +",
                                                          "ContrastTypeTail +",
                                                          "ContrastTypeTail:ConditionLabel",
                                                          "> 0")),
                                                  alpha = .11)
  new_old.time <- proc.time() - t
  ## Save all the results
  saveRDS(new_old.lmer.model, paste0(save_path, "lmerModel.rds"))
  saveRDS(new_old.lmer.anova, paste0(save_path, "lmerAnova.rds"))
  saveRDS(new_old.lmer.emmeans, paste0(save_path, "lmerEMmeans.rds"))
  lapply(seq_along(new_old.brms.models),
         function(i){
           saveRDS(new_old.brms.models[[i]],
                   paste0(save_path, "brmsModel", i, ".rds"))
         })
  saveRDS(new_old.brms.bayes_factors, paste0(save_path, "brmsBF.rds"))
  saveRDS(new_old.brms.emmeans, paste0(save_path, "brmsEMmeans.rds"))
  saveRDS(new_old.brms.emmeans.bayes_factor, paste0(save_path, "brmsEMmeansBF.rds"))
}else{
  ## Read all the results
  new_old.lmer.model <- readRDS(paste0(save_path, "lmerModel.rds"))
  new_old.lmer.anova <- readRDS(paste0(save_path, "lmerAnova.rds"))
  new_old.lmer.emmeans <- readRDS(paste0(save_path, "lmerEMmeans.rds"))
  new_old.brms.models <- lapply(1:5,
                                function(i){
                                  readRDS(paste0(save_path,
                                                 "brmsModel", i, ".rds"))
                                })
  new_old.brms.bayes_factors <- readRDS(paste0(save_path, "brmsBF.rds"))
  new_old.brms.emmeans <- readRDS(paste0(save_path, "brmsEMmeans.rds"))
  new_old.brms.emmeans.bayes_factors <- readRDS(paste0(save_path, "brmsEMmeansBF.rds"))
}
# Plot jitter + mean&se
generate_plots <- F
if(generate_plots){
  ## Get brm predicted values
  new_old.raw_predictions <- last(new_old.brms.models) %>%
    predict(summary = F,
            transform = function(x){sin(x + asin(sqrt(.5)))^2}) %>%
    t() %>%
    as_tibble() %>%
    mutate(RowNames = 1:nrow(.))
  new_old.predicted <- new_old %>%
    mutate(RowNames = 1:nrow(.)) %>%
    select(ContrastType, Condition, RowNames) %>%
    inner_join(new_old.raw_predictions) %>%
    select(-RowNames) %>%
    gather(key = Sample, value = Predicted, -c(ContrastType, Condition))
  new_old.predicted.hpdi.97 <- new_old.predicted %>%
    select(-Sample) %>%
    split(list(.$ContrastType, .$Condition)) %>%
    lapply(function(df){
      hpdi <- if(length(df$Predicted)>1){
        as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.97)
      }else{
        matrix(ncol = 2, dimnames = list(NULL, c("lower", "upper"))) # When no cases
      }
      df.summary <- df %>%
        group_by(ContrastType, Condition) %>%
        summarise(Mode = mlv(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows() %>%
    drop_na()
  new_old.predicted.hpdi.89 <- new_old.predicted %>%
    select(-Sample) %>%
    split(list(.$ContrastType, .$Condition)) %>%
    lapply(function(df){
      hpdi <- if(length(df$Predicted)>1){
        as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.89)
      }else{
        matrix(ncol = 2, dimnames = list(NULL, c("lower", "upper"))) # When no cases
      }
      df.summary <- df %>%
        group_by(ContrastType, Condition) %>%
        summarise(Mode = mlv(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows() %>%
    drop_na()
  new_old.predicted.hpdi.67 <- new_old.predicted %>%
    select(-Sample) %>%
    split(list(.$ContrastType, .$Condition)) %>%
    lapply(function(df){
      hpdi <- if(length(df$Predicted)>1){
        as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.67)
      }else{
        matrix(ncol = 2, dimnames = list(NULL, c("lower", "upper"))) # When no cases
      }
      df.summary <- df %>%
        group_by(ContrastType, Condition) %>%
        summarise(Mode = mlv(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows() %>%
    drop_na()
  ## Plot raincloud + predicted mean&sd per FstLst
  new_old.plot <- ggplot(new_old,
                         aes(x = Condition, y = Prop,
                             colour = Condition,
                             fill = Condition)) +
    theme_bw() + ylab("Looking to New Feature (Prop)") +
    geom_hline(yintercept = .5, colour = "black", linetype = 2) +
    theme(legend.position = "top",
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    coord_flip() + facet_grid(.~ContrastType) +
    geom_flat_violin(position = position_nudge(x = .2),
                     colour = "black", alpha = .5, width = .7) +
    geom_point(position = position_jitter(width = .15),
               size = 1, alpha = .6,
               show.legend = F) +
    geom_boxplot(width = .1, alpha = .3, outlier.shape = NA, colour = "black",
                 show.legend = F) +
    # geom_pointrange(data = new_old.predicted.hpdi.67,
    #                 aes(x = Condition, y = Mode, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = 1.5, size = 1.5,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    # geom_pointrange(data = new_old.predicted.hpdi.89,
    #                 aes(x = Condition, y = Mode, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = .5, size = 1,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    # geom_pointrange(data = new_old.predicted.hpdi.97,
    #                 aes(x = Condition, y = Mode, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = .5, size = .5,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")
  ## Save plot
  ggsave(paste0(save_path, "data.pdf"),
         new_old.plot,
         width = 5.5, height = 3)
}

# CONTRAST TEST ANALYSIS: PROP NEW FEATURE LOOKING TIME COURSE BY FSTLST  ==========================
save_path <- "../results/infants/OldNew/TimeCourse_"
# Data preparation
LT.new_old.time_course <- LT.test.ctr %>%
  subset(ContrastType %in% c("Tail", "Head")) %>%
  make_time_sequence_data(time_bin_size = 50,
                          aois = "NewFeature",
                          predictor_columns=c("Condition",
                                              "ContrastType"),
                          summarize_by = "Participant")
# Run bootstrapped cluster-based permutation analysis
#### NOT ENOuGH DATA

# PLOT
generate_plots <- F
if(generate_plots){
  LT.new_old.time_course.plot <- ggplot(LT.new_old.time_course,
                                        aes(x = Time, y=Prop,
                                            colour=Condition,
                                            fill=Condition)) +
    xlab('Time in Trial') + ylab("Looking to Tail (Prop)") +
    facet_grid(ContrastType~.) +
    theme(legend.position = "top") + ylim(0,1) +
    stat_summary(fun.y='mean', geom='line', linetype = '61') +
    stat_summary(fun.data=mean_se, geom='ribbon', alpha= .25, colour=NA) +
    geom_hline(yintercept = .5)
  ggsave(paste0(save_path, "data.pdf"),
         plot = LT.new_old.time_course.plot,
         width = 3.5, height = 5)
}

# WORD LEARNING TEST ANALYSIS: PROP TARGET FOR LABEL CONDITION =====================================
save_path <- "../results/infants/WordLearning/TrialAverage_"
# Prepare dataset
prop_target <- LT.test.wl %>%
  subset(Condition == "Label") %>%
  subset_by_window(window_end_col = "TrialEnd", rezero = F) %>%
  drop_na(PrePost) %>%
  make_time_window_data(aois = "Target",
                        predictor_columns = c("PrePost",
                                              "CategoryName")) %>%
  drop_na(ArcSin) %>%
  mutate(ChanceArcsin = ArcSin - asin(sqrt(.5))) # Value centered on chance looking, useful for test
## Check for amount of data available, and word-learning score
prop_target.participants <- prop_target %>%
  group_by(Participant) %>%
  summarise(nTrials = n_distinct(TrialId),
            Score = sum(SamplesInAOI/SamplesTotal)/nTrials,
            ScoreAboveChance = Score > .5,
            nCorrect = sum((SamplesInAOI/SamplesTotal) > .5),
            Perfect = nTrials == nCorrect)
# Testing in general
run_model <- F # Running the models takes around 2 minutes on a 4.40GHz 12-core
if(run_model){
  t <- proc.time()
  ## Run lmer
  prop_target.lmer.model <- lmer(ChanceArcsin ~ PrePost + (1 + PrePost | Participant),
                                    data = prop_target)
  prop_target.lmer.anova <- anova(prop_target.lmer.model)
  ## Run brms
  ### Set priors for models other than intercept-only
  priors.prop_target <- list(set_prior("uniform(-.8,.8)",
                                       class = "Intercept"),
                             c(set_prior("uniform(-.8,.8)",
                                         class = "Intercept"),
                               set_prior("normal(0,.5)", class = "b")))
  ### Set all nested formulas for model comparisons
  formulas.prop_target <- list(ChanceArcsin ~ 1 + (1 | Participant),
                               ChanceArcsin ~ PrePost + (1 + PrePost | Participant))
  ### Get brms results
  prop_target.brms.results <- bayes_factor.brm_fixef(formulas.prop_target,
                                                     prop_target,
                                                     priors.prop_target,
                                                     control = list(adapt_delta = .9999,
                                                                    max_treedepth = 20))
  prop_target.brms.models <- prop_target.brms.results[[1]]
  prop_target.brms.bayes_factor <- prop_target.brms.results[[2]]
  prop_target.time <- proc.time() - t
  ## Save all the results
  saveRDS(prop_target.lmer.model, paste0(save_path, "lmerModel.rds"))
  saveRDS(prop_target.lmer.anova, paste0(save_path, "lmerAnova.rds"))
  lapply(seq_along(prop_target.brms.models),
         function(i){
           saveRDS(prop_target.brms.models[[i]],
                   paste0(save_path, "brmsModel", i, ".rds"))
         })
  saveRDS(prop_target.brms.bayes_factor, paste0(save_path, "brmsBF.rds"))
}else{
  ## Read all the results
  prop_target.lmer.model <- readRDS(paste0(save_path, "lmerModel.rds"))
  prop_target.lmer.anova <- readRDS(paste0(save_path, "lmerAnova.rds"))
  prop_target.brms.models <- lapply(1:2,
                                    function(i){
                                      readRDS(paste0(save_path,
                                                     "brmsModel", i, ".rds"))
                                    })
  prop_target.brms.bayes_factor <- readRDS(paste0(save_path, "brmsBF.rds"))
}

# Plot jitter + mean&se
generate_plots <- F
if(generate_plots){
  ## Get brm predicted values
  prop_target.raw_predictions <- last(prop_target.brms.models) %>%
    predict(summary = F,
            transform = function(x){sin(x + asin(sqrt(.5)))^2}) %>%
    t() %>%
    as_tibble() %>%
    mutate(RowNames = 1:nrow(.))
  prop_target.predicted <- prop_target %>%
    mutate(RowNames = 1:nrow(.)) %>%
    select(PrePost, RowNames) %>%
    inner_join(prop_target.raw_predictions) %>%
    select(-RowNames) %>%
    gather(key = Sample, value = Predicted, -c(PrePost))
  prop_target.predicted.hpdi.97 <- prop_target.predicted %>%
    select(-Sample) %>%
    split(list(.$PrePost)) %>%
    lapply(function(df){
      hpdi <- if(length(df$Predicted)>1){
        as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.97)
      }else{
        matrix(ncol = 2, dimnames = list(NULL, c("lower", "upper"))) # When no cases
      }
      df.summary <- df %>%
        group_by(PrePost) %>%
        summarise(Mode = mlv(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows()
  prop_target.predicted.hpdi.89 <- prop_target.predicted %>%
    select(-Sample) %>%
    split(list(.$PrePost)) %>%
    lapply(function(df){
      hpdi <- if(length(df$Predicted)>1){
        as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.89)
      }else{
        matrix(ncol = 2, dimnames = list(NULL, c("lower", "upper"))) # When no cases
      }
      df.summary <- df %>%
        group_by(PrePost) %>%
        summarise(Mode = mlv(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows()
  prop_target.predicted.hpdi.67 <- prop_target.predicted %>%
    select(-Sample) %>%
    split(list(.$PrePost)) %>%
    lapply(function(df){
      hpdi <- if(length(df$Predicted)>1){
        as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.67)
      }else{
        matrix(ncol = 2, dimnames = list(NULL, c("lower", "upper"))) # When no cases
      }
      df.summary <- df %>%
        group_by(PrePost) %>%
        summarise(Mode = mlv(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows()
  ## Plot raincloud + predicted mean&sd per FstLst
  prop_target.plot <- ggplot(prop_target,
                         aes(x = PrePost, y = Prop,
                             colour = PrePost,
                             fill = PrePost)) +
    theme_bw() + ylab("Looking to Target (Prop)") +
    geom_hline(yintercept = .5, colour = "black", linetype = 2) +
    theme(legend.position = "top",
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    scale_x_discrete(limits = c("Post Label Onset", "Pre Label Onset")) + coord_flip() +
    geom_flat_violin(position = position_nudge(x = .2),
                     colour = "black", alpha = .5, width = .7) +
    geom_point(position = position_jitter(width = .15),
               size = 1, alpha = .6,
               show.legend = F) +
    geom_boxplot(width = .1, alpha = .3, outlier.shape = NA, colour = "black",
                 show.legend = F) +
    geom_pointrange(data = prop_target.predicted.hpdi.67,
                    aes(x = PrePost, y = Mode, ymin = lb, ymax = ub),
                    colour = brewer.pal(3, "Dark2")[[3]],
                    fatten = 1.5, size = 1.5,
                    position = position_nudge(x = -.23),
                    show.legend = F) +
    geom_pointrange(data = prop_target.predicted.hpdi.89,
                    aes(x = PrePost, y = Mode, ymin = lb, ymax = ub),
                    colour = brewer.pal(3, "Dark2")[[3]],
                    fatten = .5, size = 1,
                    position = position_nudge(x = -.23),
                    show.legend = F) +
    geom_pointrange(data = prop_target.predicted.hpdi.97,
                    aes(x = PrePost, y = Mode, ymin = lb, ymax = ub),
                    colour = brewer.pal(3, "Dark2")[[3]],
                    fatten = .5, size = .5,
                    position = position_nudge(x = -.23),
                    show.legend = F) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")
  ## Save plot
  ggsave(paste0(save_path, "data.pdf"),
         prop_target.plot,
         width = 4, height = 3.5)
}

# WORD LEARNING TEST ANALYSIS: PROP TARGET TIME COURSE FOR LABEL CONDITION  ========================
save_path <- "../results/infants/WordLearning/TimeCourse_"
# Data preparation
## For overall analysis
prop_target.time_course <- LT.test.wl %>%
  subset_by_window(window_start_col = "PhraseOnset", remove = F) %>%
  mutate(TrialEnd = TrialEnd - PhraseOnset) %>%
  subset_by_window(window_end_col = "TrialEnd", rezero = F) %>%
  mutate(Chance = F) %>%
  make_time_sequence_data(time_bin_size = 100,
                          aois = "Target",
                          predictor_columns=c("Chance"),
                          summarize_by = c("Participant"))
prop_target.time_course.chance <- prop_target.time_course %>%
  mutate(Chance = T,
         Participant = paste0("Chance", Participant),
         Prop = .5)
prop_target.time_course.chance_test <- rbind(prop_target.time_course,
                                             prop_target.time_course.chance) %>%
  mutate(Chance = as.character(Chance)) %>%
  mutate_at("Chance", parse_factor, levels = NULL)
## For analysis by label
prop_target.time_course.by_label <- LT.test.wl %>%
  subset_by_window(window_start_col = "LabelOnset", remove = F) %>%
  mutate(TrialEnd = TrialEnd - LabelOnset) %>%
  subset_by_window(window_end_col = "TrialEnd", rezero = F) %>%
  mutate(Chance = F) %>%
  make_time_sequence_data(time_bin_size = 100,
                          aois = "Target",
                          predictor_columns=c("Chance", "Stimulus"),
                          summarize_by = c("Participant"))
prop_target.time_course.by_label.chance <- prop_target.time_course.by_label %>%
  mutate(Chance = T,
         Participant = paste0("Chance", Participant),
         Prop = .5)
prop_target.time_course.by_label.chance_test <- rbind(prop_target.time_course.by_label,
                                                      prop_target.time_course.by_label.chance) %>%
  mutate(Chance = as.character(Chance)) %>%
  mutate_at("Chance", parse_factor, levels = NULL)
# BOOTSTRAPPED CLUSTER-BASED PERMUTATION ANALYSIS
run_model <- F
if(run_model){
  t <- proc.time()
  ## Determine threshold based on alpha = .05 two-tailed
  num_sub = length(unique((prop_target.time_course.chance_test$Participant)))
  threshold_t = qt(p = 1 - .05/2,
                   df = num_sub-1)
  ## Determine clusters overall and by label
  prop_target.time_cluster <- prop_target.time_course.chance_test %>%
    make_time_cluster_data(predictor_column = "Chance",
                           treatment_level = "TRUE",
                           aoi = "Target",
                           test = "t.test",
                           threshold = threshold_t)
  prop_target.time_cluster.by_label <- prop_target.time_course.by_label.chance_test %>%
    split(.$Stimulus) %>%
    lapply(make_time_cluster_data,
           predictor_column = "Chance",
           treatment_level = "TRUE",
           aoi = "Target",
           test = "t.test",
           threshold = threshold_t)
  ## Run analysis overall and by label
  prop_target.time_cluster.analysis <- prop_target.time_cluster %>%
    analyze_time_clusters(within_subj = F,
                          parallel = T)
  prop_target.time_cluster.by_label.analysis <- prop_target.time_cluster.by_label %>%
    lapply(analyze_time_clusters,
           within_subj = F,
           parallel = T)
  prop_target.time_course.time <- proc.time() - t
  ## Save clusters and analysis
  saveRDS(prop_target.time_cluster, paste0(save_path, "bcbpClusters.rds"))
  saveRDS(prop_target.time_cluster.analysis, paste0(save_path, "bcbpAnalysis.rds"))
  saveRDS(prop_target.time_cluster.by_label, paste0(save_path, "bcbpClusters_byLabel.rds"))
  saveRDS(prop_target.time_cluster.by_label.analysis, paste0(save_path, "bcbpAnalysis_byLabel.rds"))
}else{
  ## Read the results
  prop_target.time_cluster <- readRDS(paste0(save_path, "bcbpClusters.rds"))
  prop_target.time_cluster.analysis <- readRDS(paste0(save_path, "bcbpAnalysis.rds"))
  prop_target.time_cluster.by_label <- readRDS(paste0(save_path, "bcbpClusters_byLabel.rds"))
  prop_target.time_cluster.by_label.analysis <- readRDS(paste0(save_path,
                                                               "bcbpAnalysis_byLabel.rds"))
}

# PLOT
generate_plots <- F
if(generate_plots){
  # Plot overall
  prop_target.time_course.plot.clusters <- prop_target.time_cluster.analysis$clusters %>%
    subset(Probability < .05)
  prop_target.time_course.plot <- ggplot(prop_target.time_course,
                                         aes(x = Time, y=Prop,
                                             colour = Chance,
                                             fill = Chance)) +
    xlab('Time in Trial') + ylab("Looking to Tail (Prop)") + theme_bw() +
    theme(legend.position = "none") + ylim(0,1) +
    stat_summary(fun.y='mean', geom='line', linetype = '61') +
    stat_summary(fun.data=mean_se, geom='ribbon', alpha= .25, colour=NA) +
    geom_rect(data = prop_target.time_course.plot.clusters,
              inherit.aes = F,
              aes(xmin = StartTime, xmax = EndTime,
                  ymin = 0, ymax = 1),
              alpha = 0.5,
              fill = brewer.pal(3, "Dark2")[[3]]) +
    geom_hline(yintercept = .5)
  ggsave(paste0(save_path, "data.pdf"),
         plot = prop_target.time_course.plot,
         width = 5.5, height = 2.5)
  # Plot by label
  prop_target.time_course.plot.clusters.by_label <- prop_target.time_cluster.by_label.analysis %>%
  {lapply(seq_along(.),
          function(i){
            df <- prop_target.time_cluster.by_label.analysis[[i]]$clusters %>%
              mutate(Stimulus = if(i == 1){"WLG"}else{"WLS"})
          })} %>%
    bind_rows() %>%
    subset(Probability < .05)
  prop_target.time_course.plot <- ggplot(prop_target.time_course.by_label,
                                         aes(x = Time, y=Prop,
                                             colour=Stimulus,
                                             fill=Stimulus)) +
    xlab('Time in Trial') + ylab("Looking to Tail (Prop)") +
    theme(legend.position = "top") + ylim(0,1) + facet_grid(.~Stimulus) +
    stat_summary(fun.y='mean', geom='line', linetype = '61') +
    stat_summary(fun.data=mean_se, geom='ribbon', alpha= .25, colour=NA) +
    geom_rect(data = prop_target.time_course.plot.clusters.by_label,
              inherit.aes = F,
              aes(xmin = StartTime, xmax = EndTime,
                  ymin = 0, ymax = 1),
              alpha = 0.5,
              fill = brewer.pal(3, "Dark2")[[3]]) +
    geom_hline(yintercept = .5)
  ggsave(paste0(save_path, "byLabel_data.pdf"),
         plot = prop_target.time_course.plot,
         width = 5.5, height = 2.5)
}

# FAMILIARISATION: NUMBER OF SWITCHES ==============================================================
save_path <- "../results/infants/FamSwitches/FstLst_"
# Prepare dataset
fam_switches.fstlst <- LT.fam %>%
  drop_na(Tail, FstLst) %>%
  group_by(Participant, TrialId) %>%
  summarise(Switches = sum(Tail != lag(Tail), na.rm = T), # Count switches per trial per participant
            # Keep columns for analysis (only one value per trial per participant)
            FstLst = first(FstLst),
            Condition = first(Condition)) %>%
  ungroup()
# Testing Switches ~ Condition*FstLst
run_model <- F # Running the models takes around 5 minutes on a 4.40GHz 12-core
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
                                         family = poisson(),
                                         iter = 4000,
                                         controls = list(adapt_delta = .95))

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

# Plotting boxplots
generate_plots <- F
if(generate_plots){
  ## Get brm predicted values
  fam_switches.raw_predictions <- last(fam_switches.per_fstlst.brms.models) %>%
    predict(summary = F) %>%
    t() %>%
    as_tibble() %>%
    mutate(RowNames = 1:nrow(.))
  fam_switches.predicted <- fam_switches.fstlst %>%
    mutate(RowNames = 1:nrow(.)) %>%
    select(FstLst, Condition, RowNames) %>%
    inner_join(fam_switches.raw_predictions) %>%
    select(-RowNames) %>%
    gather(key = Sample, value = Predicted, -c(FstLst, Condition))
  fam_switches.predicted.hpdi.97 <- fam_switches.predicted %>%
    select(-Sample) %>%
    split(list(.$FstLst, .$Condition)) %>%
    lapply(function(df){
      hpdi <- as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.97)
      df.summary <- df %>%
        group_by(FstLst, Condition) %>%
        summarise(Mode = mlv(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows()
  fam_switches.predicted.hpdi.89 <- fam_switches.predicted %>%
    select(-Sample) %>%
    split(list(.$FstLst, .$Condition)) %>%
    lapply(function(df){
      hpdi <- as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.89)
      df.summary <- df %>%
        group_by(FstLst, Condition) %>%
        summarise(Mode = mlv(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows()
  fam_switches.predicted.hpdi.67 <- fam_switches.predicted %>%
    select(-Sample) %>%
    split(list(.$FstLst, .$Condition)) %>%
    lapply(function(df){
      hpdi <- as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.67)
      df.summary <- df %>%
        group_by(FstLst, Condition) %>%
        summarise(Mode = mlv(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = hpdi[1,"upper"])
      return(df.summary)
    }) %>%
    bind_rows()
  ## Plot raincloud + predicted mean&sd per FstLst
  fam_switches.per_fstlst.plot <- ggplot(fam_switches.fstlst,
                                         aes(x = Condition, y = Switches,
                                             colour = Condition,
                                             fill = Condition)) +
    theme_bw() + ylab("Number of switches between AOIs") +
    theme(legend.position = "top",
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    coord_flip() + facet_grid(.~FstLst) +
    geom_flat_violin(position = position_nudge(x = .2),
                     colour = "black", alpha = .5, width = .7) +
    geom_point(position = position_jitter(width = 0.15, height = 0.1),
               size = 1, alpha = .6,
               show.legend = F) +
    geom_boxplot(width = .1, alpha = .3, outlier.shape = NA, colour = "black",
                 show.legend = F) +
    # geom_pointrange(data = fam_switches.predicted.hpdi.67,
    #                 aes(x = Condition, y = Mode, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = 1.5, size = 1.5,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    # geom_pointrange(data = fam_switches.predicted.hpdi.89,
    #                 aes(x = Condition, y = Mode, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = .5, size = 1,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    # geom_pointrange(data = fam_switches.predicted.hpdi.97,
    #                 aes(x = Condition, y = Mode, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = .5, size = .5,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")
  ## Save plot
  ggsave(paste0(save_path, "data.pdf"),
         fam_switches.per_fstlst.plot,
         width = 5.5, height = 3)
}

# FAMILIARISATION: FIRST LOOK ======================================================================
save_path <- "../results/infants/FirstLook/"
# Prepare datasets
first_look <- LT.fam %>%
  drop_na(Tail) %>%
  group_by(Participant, TrialId, Tail) %>%
  summarise(FirstAOILook = first(TimeStamp), # First time to each AOI
            AOI = parse_factor(ifelse(first(Tail), "Tail", "Head"),
                               levels = c("Head", "Tail")), # Explicit AOI flag
            # Keep columns for analysis (only one value per trial per participant)
            Condition = first(Condition),
            FstLst = first(FstLst)) %>%
  select(-Tail) %>%
  ungroup()
first_aoi.fstlst <- first_look %>%
  drop_na(FstLst) %>%
  group_by(Participant, TrialId) %>%
  arrange(FirstAOILook) %>%
  summarise(AOI = first(AOI),
            # Keep columns for analysis (only one value per trial per participant)
            Condition = first(Condition),
            FstLst = first(FstLst)) %>%
  ungroup()
first_tail.fstlst <- first_look %>%
  drop_na(FstLst) %>%
  subset(AOI == "Tail", select = -AOI) %>% # Discards trials with no Tail look: is it okay?
  mutate(logFirstAOILook = log(FirstAOILook))

# Testing (First)AOI ~ Condition*FstLst
run_model <- F # Running the models takes around 5 minutes on a 4.40GHz 12-core
if(run_model){
  t <- proc.time()
  ## Run (g)lmer
  first_aoi.per_fstlst.glmer.model <- glmer(AOI ~ FstLst*Condition +
                                              (1 + FstLst | Participant),
                                            data = first_aoi.fstlst,
                                            family = binomial())
  # Current p-values from summary may not be the best. Do something else?
  ## Run brms
  ### Set priors for models other than intercept-only
  priors.first_aoi.per_fstlst <- list(NULL,
                                     set_prior("normal(0,.5)", class = "b"))
  ### Set all nested formulas for model comparisons
  formulas.first_aoi.per_fstlst <- list(AOI ~ 1 +
                                          (1 | Participant),
                                        AOI ~ FstLst +
                                          (1 + FstLst | Participant),
                                        AOI ~ FstLst + Condition +
                                          (1 + FstLst | Participant),
                                        AOI ~ FstLst + Condition +
                                          FstLst:Condition +
                                          (1 + FstLst | Participant))
  ### Get brms results
  brms.results <- bayes_factor.brm_fixef(formulas.first_aoi.per_fstlst,
                                         first_aoi.fstlst,
                                         priors.first_aoi.per_fstlst,
                                         family = bernoulli())
  first_aoi.per_fstlst.brms.models <- brms.results[[1]]
  first_aoi.per_fstlst.brms.bayes_factors <- brms.results[[2]]
  first_aoi.time <- proc.time() - t
  ## Save all the results
  saveRDS(first_aoi.per_fstlst.glmer.model, paste0(save_path, "FirstAOI_glmerModel.rds"))
  lapply(seq_along(first_aoi.per_fstlst.brms.models),
         function(i){
           saveRDS(first_aoi.per_fstlst.brms.models[[i]],
                   paste0(save_path, "FirstAOI_brmsModel", i, ".rds"))
         })
  saveRDS(first_aoi.per_fstlst.brms.bayes_factors, paste0(save_path, "FirstAOI_brmsBF.rds"))
}else{
  ## Read all the results
  first_aoi.per_fstlst.glmer.model <- readRDS(paste0(save_path, "FirstAOI_glmerModel.rds"))
  first_aoi.per_fstlst.brms.models <- lapply(1:4,
                                             function(i){
                                               readRDS(paste0(save_path,
                                                              "FirstAOI_brmsModel", i, ".rds"))
                                             })
  first_aoi.per_fstlst.brms.bayes_factors <- readRDS(paste0(save_path, "FirstAOI_brmsBF.rds"))
}
# Testing FirstAOI(Tail)Look ~ Condition*FstLst
run_model <- F # Running the models takes around 4 minutes on a 4.40GHz 12-core
if(run_model){
  t <- proc.time()
  ## Run lmer
  first_tail.per_fstlst.lmer.model <- lmer(logFirstAOILook ~ FstLst*Condition +
                                               (1 + FstLst | Participant),
                                             data = first_tail.fstlst)
  first_tail.per_fstlst.lmer.anova <- anova(first_tail.per_fstlst.lmer.model, type = 1)
  # Current p-values from summary may not be the best. Do something else?
  ## Run brms
  ### Set priors for models other than intercept-only
  priors.first_tail.per_fstlst <- list(NULL,
                                       set_prior("normal(0,.5)", class = "b"))
  ### Set all nested formulas for model comparisons
  formulas.first_tail.per_fstlst <- list(logFirstAOILook ~ 1 +
                                           (1 | Participant),
                                         logFirstAOILook ~ FstLst +
                                           (1 + FstLst | Participant),
                                         logFirstAOILook ~ FstLst + Condition +
                                           (1 + FstLst | Participant),
                                         logFirstAOILook ~ FstLst + Condition +
                                           FstLst:Condition +
                                           (1 + FstLst | Participant))
  ### Get brms results
  brms.results <- bayes_factor.brm_fixef(formulas.first_tail.per_fstlst,
                                         first_tail.fstlst,
                                         priors.first_tail.per_fstlst)
  first_tail.per_fstlst.brms.models <- brms.results[[1]]
  first_tail.per_fstlst.brms.bayes_factors <- brms.results[[2]]
  first_tail.time <- proc.time() - t
  ## Save all the results
  saveRDS(first_tail.per_fstlst.lmer.model, paste0(save_path, "FirstTail_lmerModel.rds"))
  saveRDS(first_tail.per_fstlst.lmer.anova, paste0(save_path, "FirstTail_lmerAnova.rds"))
  lapply(seq_along(first_tail.per_fstlst.brms.models),
         function(i){
           saveRDS(first_tail.per_fstlst.brms.models[[i]],
                   paste0(save_path, "FirstTail_brmsModel", i, ".rds"))
         })
  saveRDS(first_tail.per_fstlst.brms.bayes_factors, paste0(save_path, "FirstTail_brmsBF.rds"))
}else{
  ## Read all the results
  first_tail.per_fstlst.lmer.model <- readRDS(paste0(save_path, "FirstTail_lmerModel.rds"))
  first_tail.per_fstlst.lmer.anova <- readRDS(paste0(save_path, "FirstTail_lmerAnova.rds"))
  first_tail.per_fstlst.brms.models <- lapply(1:4,
                                             function(i){
                                               readRDS(paste0(save_path,
                                                              "FirstTail_brmsModel", i, ".rds"))
                                             })
  first_tail.per_fstlst.brms.bayes_factors <- readRDS(paste0(save_path, "FirstTail_brmsBF.rds"))
}

# Plotting
generate_plots <- F
if(generate_plots){
  ## First AOI
  ### Get data for plot
  first_aoi.fstlst.to_plot <- first_aoi.fstlst %>%
    drop_na(FstLst) %>%
    group_by(FstLst, AOI, Participant, Condition) %>%
    summarise(N = n()) %>%
    ungroup()
  ### Get brm predicted values ## CURRENTLY NOT WORKING PROPERLY
  # first_aoi.raw_predictions <- last(first_aoi.per_fstlst.brms.models) %>%
  #   predict(summary = F) %>%
  #   t() %>%
  #   as_tibble() %>%
  #   mutate(RowNames = 1:nrow(.))
  # first_aoi.predicted <- first_aoi.fstlst %>%
  #   mutate(RowNames = 1:nrow(.)) %>%
  #   select(FstLst, Condition, RowNames) %>%
  #   inner_join(first_aoi.raw_predictions) %>%
  #   select(-RowNames) %>%
  #   gather(key = Sample, value = AOI, -c(FstLst, Condition)) %>%
  #   mutate(AOI = ifelse(AOI == 0, "Head", "Tail"),
  #          AOI = parse_factor(AOI, levels = NULL))
  # first_aoi.predicted.hpdi.97 <- first_aoi.predicted %>%
  #   group_by(FstLst, Condition, Sample, AOI) %>%
  #   summarise(N = n()) %>%
  #   ungroup() %>%
  #   select(-Sample) %>%
  #   split(list(.$FstLst, .$Condition, .$AOI)) %>%
  #   lapply(function(df){
  #     hpdi <- as.mcmc(df$N) %>% HPDinterval(prob = 0.97)
  #     df.summary <- df %>%
  #       group_by(FstLst, Condition, AOI) %>%
  #       summarise(Mode = first(mlv(df$N))) %>%
  #       mutate(lb = hpdi[1,"lower"],
  #              ub = hpdi[1,"upper"])
  #     return(df.summary)
  #   }) %>%
  #   bind_rows()
  # first_aoi.predicted.hpdi.89 <- first_aoi.predicted %>%
  #   group_by(FstLst, Condition, Sample, AOI) %>%
  #   summarise(N = n()) %>%
  #   ungroup() %>%
  #   select(-Sample) %>%
  #   split(list(.$FstLst, .$Condition, .$AOI)) %>%
  #   lapply(function(df){
  #     hpdi <- as.mcmc(df$N) %>% HPDinterval(prob = 0.89)
  #     df.summary <- df %>%
  #       group_by(FstLst, Condition, AOI) %>%
  #       summarise(Mode = first(mlv(df$N))) %>%
  #       mutate(lb = hpdi[1,"lower"],
  #              ub = hpdi[1,"upper"])
  #     return(df.summary)
  #   }) %>%
  #   bind_rows()
  # first_aoi.predicted.hpdi.67 <- first_aoi.predicted %>%
  #   group_by(FstLst, Condition, Sample, AOI) %>%
  #   summarise(N = n()) %>%
  #   ungroup() %>%
  #   select(-Sample) %>%
  #   split(list(.$FstLst, .$Condition, .$AOI)) %>%
  #   lapply(function(df){
  #     hpdi <- as.mcmc(df$N) %>% HPDinterval(prob = 0.67)
  #     df.summary <- df %>%
  #       group_by(FstLst, Condition, AOI) %>%
  #       summarise(Mode = first(mlv(df$N))) %>%
  #       mutate(lb = hpdi[1,"lower"],
  #              ub = hpdi[1,"upper"])
  #     return(df.summary)
  #   }) %>%
  #   bind_rows()
  ### Plot raincloud + predicted mean&sd per FstLst
  first_aoi.per_fstlst.plot <- ggplot(first_aoi.fstlst.to_plot,
                                      aes(x = Condition, y = N,
                                          colour = Condition,
                                          fill = Condition)) +
    theme_bw() + ylab("Number of first look to AOI") +
    theme(legend.position = "top",
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    scale_y_continuous(breaks = c(0, 1, 2, 3)) +
    coord_flip() + facet_grid(AOI~FstLst) +
    geom_flat_violin(position = position_nudge(x = .2),
                     colour = "black", alpha = .5, width = .7) +
    geom_point(position = position_jitter(width = 0.15, height = 0.05),
               size = 1, alpha = .6,
               show.legend = F) +
    geom_boxplot(width = .1, alpha = .3, outlier.shape = NA, colour = "black",
                 show.legend = F) +
    # geom_pointrange(data = first_aoi.predicted.hpdi.67,
    #                 aes(x = Condition,
    #                     y = Mode, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = 1.5, size = 1.5,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    # geom_pointrange(data = first_aoi.predicted.hpdi.89,
    #                 aes(x = Condition,
    #                     y = Mode, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = .5, size = 1,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    # geom_pointrange(data = first_aoi.predicted.hpdi.97,
    #                 aes(x = Condition,
    #                     y = Mode, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = .5, size = .5,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")
  ### Save plot
  ggsave(paste0(save_path, "FirstAOI_data.pdf"),
         first_aoi.per_fstlst.plot,
         width = 5.5, height = 5)
  ## Time to first tail look (boxplot)
  ### Get brm predicted values
  first_tail.raw_predictions <- last(first_tail.per_fstlst.brms.models) %>%
    predict(summary = F,
            transform = exp) %>%
    t() %>%
    as_tibble() %>%
    mutate(RowNames = 1:nrow(.))
  first_tail.predicted <- first_tail.fstlst %>%
    mutate(RowNames = 1:nrow(.)) %>%
    select(FstLst, Condition, RowNames) %>%
    inner_join(first_tail.raw_predictions) %>%
    select(-RowNames) %>%
    gather(key = Sample, value = Predicted, -c(FstLst, Condition))
  first_tail.predicted.hpdi.97 <- first_tail.predicted %>%
    select(-Sample) %>%
    split(list(.$FstLst, .$Condition)) %>%
    lapply(function(df){
      hpdi <- as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.97)
      df.summary <- df %>%
        group_by(FstLst, Condition) %>%
        summarise(Mode = mlv(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = min(hpdi[1,"upper"], 5000)) # Correct for abnormaly large estimations
      return(df.summary)
    }) %>%
    bind_rows()
  first_tail.predicted.hpdi.89 <- first_tail.predicted %>%
    select(-Sample) %>%
    split(list(.$FstLst, .$Condition)) %>%
    lapply(function(df){
      hpdi <- as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.89)
      df.summary <- df %>%
        group_by(FstLst, Condition) %>%
        summarise(Mode = mlv(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = min(hpdi[1,"upper"], 5000)) # Correct for abnormaly large estimations
      return(df.summary)
    }) %>%
    bind_rows()
  first_tail.predicted.hpdi.67 <- first_tail.predicted %>%
    select(-Sample) %>%
    split(list(.$FstLst, .$Condition)) %>%
    lapply(function(df){
      hpdi <- as.mcmc(df$Predicted) %>% HPDinterval(prob = 0.67)
      df.summary <- df %>%
        group_by(FstLst, Condition) %>%
        summarise(Mode = mlv(df$Predicted)) %>%
        mutate(lb = hpdi[1,"lower"],
               ub = min(hpdi[1,"upper"], 5000)) # Correct for abnormaly large estimations
      return(df.summary)
    }) %>%
    bind_rows()
  ### Plot raincloud + predicted mean&sd per FstLst
  first_tail.per_fstlst.plot <- ggplot(first_tail.fstlst,
                                         aes(x = Condition, y = FirstAOILook,
                                             colour = Condition,
                                             fill = Condition)) +
    theme_bw() + ylab("Time to first look to Tail") + ylim(0, 5000) +
    theme(legend.position = "top",
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    coord_flip() + facet_grid(.~FstLst) +
    geom_flat_violin(position = position_nudge(x = .2),
                     colour = "black", alpha = .5, width = .7) +
    geom_point(position = position_jitter(width = 0.15, height = 0),
               size = 1, alpha = .6,
               show.legend = F) +
    geom_boxplot(width = .1, alpha = .3, outlier.shape = NA, colour = "black",
                 show.legend = F) +
    # geom_pointrange(data = first_tail.predicted.hpdi.67,
    #                 aes(x = Condition, y = Mode, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = 1.5, size = 1.5,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    # geom_pointrange(data = first_tail.predicted.hpdi.89,
    #                 aes(x = Condition, y = Mode, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = .5, size = 1,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    # geom_pointrange(data = first_tail.predicted.hpdi.97,
    #                 aes(x = Condition, y = Mode, ymin = lb, ymax = ub),
    #                 colour = brewer.pal(3, "Dark2")[[3]],
    #                 fatten = .5, size = .5,
    #                 position = position_nudge(x = -.23),
    #                 show.legend = F) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2")
  ### Save plot
  ggsave(paste0(save_path, "FirstTail_data.pdf"),
         first_tail.per_fstlst.plot,
         width = 5.5, height = 3)
}
