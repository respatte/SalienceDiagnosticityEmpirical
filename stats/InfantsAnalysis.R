# LIBRARY IMPORTS ==================================================================================
library(eyetrackingR)
library(lme4)
library(lmerTest)
library(brms)
library(tidyverse)

source("Routines.R")
source("StatTools.R")

# GATHER DATA ======================================================================================
# Load data and run general checks
d <- LT_data.gather("infants")
# Unload snow packages so that parallel works for brms
detach("package:doSNOW")
detach("package:snow")
# Check for counterbalancing, gender balancing, and age balancing
pres_seq <- d[[4]] %>%
  group_by(PresentationSequence, Participant) %>%
  summarise(T = sum(TrackLoss)/n())
# Keep only one participant when multiple infants saw the same presentation sequence.
# Current choice: improve gender balance in total and between conditions.
# remove: P03, P53, P55
# keep:   P51, P06, P08
d[[4]] <- d[[4]] %>%
  subset(!(Participant %in% c("P03","P53","P55"))) %>%
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
  mutate(CursorX = case_when(Participant == "P01" ~ CursorX - 60,
                             Participant == "P07" ~ CursorX - 80,
                             Participant == "P08" ~ CursorX + 80,
                             Participant == "P26" ~ CursorX - 80,
                             Participant == "P31" ~ CursorX + 30,
                             Participant == "P42" ~ CursorX - 30,
                             Participant == "P51" ~ CursorX + 60,
                             Participant == "P56" ~ CursorX + 30,
                             Participant == "P60" ~ CursorX - 100,
                             Participant == "P61" ~ CursorX - 80,
                             T ~ as.double(CursorX)),
         CursorY = case_when(Participant == "P07" ~ CursorY + 100,
                             Participant == "P08" ~ CursorY - 50,
                             Participant == "P14" ~ CursorY - 80,
                             Participant == "P15" ~ CursorY - 150,
                             Participant == "P31" ~ CursorY - 100,
                             Participant == "P36" ~ CursorY - 50,
                             Participant == "P37" ~ CursorY - 50,
                             Participant == "P42" ~ CursorY - 50,
                             Participant == "P47" ~ CursorY - 50,
                             Participant == "P52" ~ CursorY - 50,
                             Participant == "P59" ~ CursorY - 150,
                             Participant == "P60" ~ CursorY - 150,
                             Participant == "P61" ~ CursorY - 100,
                             Participant == "P64" ~ CursorY - 100,
                             Participant == "P65" ~ CursorY - 50,
                             Participant == "P66" ~ CursorY - 50,
                             Participant == "P67" ~ CursorY - 100,
                             T ~ as.double(CursorY))) %>%
  select(-c(Head, Tail, NewTail, OldTail, NewHead, OldHead, Centre, Target, Distractor))
  # Removing AOIs as they need updating
## Saving heatmaps per participant pre-correction
generate_plots <- F
if(generate_plots){
  head_fam <- tibble(AOI_type = factor(c("Reg", "Flip")),
                     xmin = c(1031,1920-1031-450), xmax = c(1031+450,1920-1031),
                     ymin = c(197,197), ymax = c(197+450,197+450))
  tail_fam <- tibble(AOI_type = factor(c("Reg", "Flip")),
                     xmin = c(390,1920-390-450), xmax = c(390+450,1920-390),
                     ymin = c(299,299), ymax = c(299+450,299+450))
  LT.gaze_offset.data.post <- LT.gaze_offset.data.correction %>%
    subset(CursorX != d[[4]]$CursorX | CursorY != d[[4]]$CursorY) %>%
    # Only generate graphs with corrected gaze values
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
                                            "CategoryName"))
# Testing Prop ~ FstLst*Condition
run_model <- T
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
                                         prop_tail.per_fstlst,
                                         priors.prop_tail.per_fstlst)
  prop_tail.per_fstlst.brms.models <- brms.results[[1]]
  prop_tail.per_fstlst.brms.bayes_factors <- brms.results[[2]]
  prop_tail.time <- proc.time() - t
  ## Save all the results
  saveRDS(prop_tail.per_fstlst.lmer.model, paste0(save_path, "FstLst_lmerModel.rds"))
  saveRDS(prop_tail.per_fstlst.lmer.anova, paste0(save_path, "FstLst_lmerAnova.rds"))
  saveRDS(prop_tail.per_fstlst.brms.models, paste0(save_path, "FstLst_brmsModels.rds"))
  saveRDS(prop_tail.per_fstlst.brms.bayes_factors, paste0(save_path, "FstLst_brmsBF.rds"))
}else{
  ## Read all the results
  prop_tail.per_fstlst.lmer.model <- readRDS(paste0(save_path, "FstLst_lmerModel.rds"))
  prop_tail.per_fstlst.lmer.anova <- readRDS(paste0(save_path, "FstLst_lmerAnova.rds"))
  prop_tail.per_fstlst.brms.models <- readRDS(paste0(save_path, "FstLst_brmsModels.rds"))
  prop_tail.per_fstlst.brms.bayes_factors <- readRDS(paste0(save_path, "FstLst_brmsBF.rds"))
}

# Plot jitter + mean&se + lines
generate_plots <- F
if(generate_plots){
  ## Plot per FstLst
  LT.prop_tail.per_part.plot <- ggplot(LT.prop_tail.fstlst,
                                       aes(x = FstLst, y = Prop,
                                           colour = Condition,
                                           fill = Condition)) +
    theme(legend.pos = "top") + ylab("Looking to Tail (Prop)") +
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
  ggsave(paste0(save_path, "FstLst_data.pdf"),
         LT.prop_tail.per_part.plot,
         width = 7, height = 5.4)
}

# FAMILIARISATION ANALYSIS: PROP TAIL LOOKING TIME COURSE BY FSTLST  ===============================
save_path <- "../results/infants/PropTail/TimeCourse_"
# Data preparation
LT.time_course_tail <- LT.fam %>%
  drop_na(FstLst) %>%
  subset_by_window(window_start_col = "LabelOnset") %>%
  make_time_sequence_data(time_bin_size = 50,
                          aois = "Tail",
                          predictor_columns=c("Condition",
                                              "FstLst"),
                          summarize_by = "Participant")
# BOOTSTRAPPED CLUSTER-BASED PERMUTATION ANALYSIS
run_model <- F
if(run_model){
  t <- proc.time()
  ## Determine threshold based on alpha = .05 two-tailed
  num_sub = length(unique((LT.time_course_tail$Participant)))
  threshold_t = qt(p = 1 - .05/2,
                   df = num_sub-1)
  ## Determine clusters
  LT.time_cluster_tail <- LT.time_course_tail %>%
    split(.$FstLst) %>%
    lapply(make_time_cluster_data,
           predictor_column = "Condition",
           treatment_level = "No Label",
           aoi = "Tail",
           test = "t.test",
           threshold = threshold_t)
  ## Run analysis
  LT.time_cluster_tail.analysis <- LT.time_cluster_tail %>%
    lapply(analyze_time_clusters,
           within_subj = F,
           parallel = T)
  bcbp.time <- proc.time() - t
  ## Save clusters and analysis
  saveRDS(LT.time_cluster_tail, paste0(save_path, "FstLst_bcbpClusters.rds"))
  saveRDS(LT.time_cluster_tail.analysis, paste0(save_path, "FstLst_bcbpAnalysis.rds"))
}else{
  ## Read the results
  LT.time_cluster_tail <- readRDS(paste0(save_path, "FstLst_bcbpClusters.rds"))
  LT.time_cluster_tail.analysis <- readRDS(paste0(save_path, "FstLst_bcbpAnalysis.rds"))
}

# PLOT
generate_plots <- F
if(generate_plots){
  LT.fam.time_course.plot.blocks <- ggplot(LT.time_course_tail,
                                           aes(x = Time, y=Prop,
                                               colour=Condition,
                                               fill=Condition)) +
    xlab('Time in Trial') + ylab("Looking to Tail (Prop)") +
    facet_grid(FstLst~.) +
    theme(legend.position = "top") + ylim(0,1) +
    stat_summary(fun.y='mean', geom='line', linetype = '61') +
    stat_summary(fun.data=mean_se, geom='ribbon', alpha= .25, colour=NA) +
    geom_hline(yintercept = .5)
  ggsave(paste0(save_path, "FstLst_data.pdf"),
         plot = LT.fam.time_course.plot.blocks,
         width = 3.5, height = 5)
}

# FAMILIARISATION ANALYSIS: PROP AOI LOOKING PRE/POST LABEL ONSET ==================================
save_path <- "../results/infants/PrePost/TrialAverage_"
# Prepare dataset
prop_tail.pre_post.fstLst <- LT.fam %>%
  drop_na(PrePost, FstLst) %>%
  make_time_window_data(aois=c("Tail"),
                        predictor_columns=c("Condition",
                                            "FstLst",
                                            "PrePost",
                                            "Stimulus",
                                            "CategoryName")) %>%
# Testing Prop ~ FstLst*Condition
run_model <- T
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
                                         (1 | Stimulus) +
                                         (1 | StimLabel),
                                       ArcSin ~ 1 + FstLst +
                                         (1 + FstLst | Participant) +
                                         (1 | Stimulus) +
                                         (1 | StimLabel),
                                       ArcSin ~ 1 + FstLst + PrePost +
                                         (1 + FstLst | Participant) +
                                         (1 | Stimulus) +
                                         (1 | StimLabel),
                                       ArcSin ~ 1 + FstLst + PrePost + Condition +
                                         (1 + FstLst | Participant) +
                                         (1 | Stimulus) +
                                         (1 | StimLabel),
                                       ArcSin ~ 1 + FstLst + PrePost + Condition +
                                         FstLst:PrePost +
                                         (1 + FstLst | Participant) +
                                         (1 | Stimulus) +
                                         (1 | StimLabel),
                                       ArcSin ~ 1 + FstLst + PrePost + Condition +
                                         FstLst:PrePost + FstLst:Condition +
                                         (1 + FstLst | Participant) +
                                         (1 | Stimulus) +
                                         (1 | StimLabel),
                                       ArcSin ~ 1 + FstLst + PrePost + Condition +
                                         FstLst:PrePost + FstLst:Condition +
                                         PrePost:Condition +
                                         (1 + FstLst | Participant) +
                                         (1 | Stimulus) +
                                         (1 | StimLabel),
                                       ArcSin ~ 1 + FstLst + PrePost + Condition +
                                         FstLst:PrePost + FstLst:Condition +
                                         PrePost:Condition +
                                         FstLst:PrePost:Condition +
                                         (1 + FstLst | Participant) +
                                         (1 | Stimulus) +
                                         (1 | StimLabel))
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
  saveRDS(pre_post.per_fstlst.brms.models, paste0(save_path, "FstLst_brmsModels.rds"))
  saveRDS(pre_post.per_fstlst.brms.bayes_factors, paste0(save_path, "FstLst_brmsBF.rds"))
}else{
  ## Read all the results
  pre_post.per_fstlst.lmer.model <- readRDS(paste0(save_path, "FstLst_lmerModel.rds"))
  pre_post.per_fstlst.lmer.anova <- readRDS(paste0(save_path, "FstLst_lmerAnova.rds"))
  pre_post.per_fstlst.brms.models <- readRDS(paste0(save_path, "FstLst_brmsModels.rds"))
  pre_post.per_fstlst.brms.BF <- readRDS(paste0(save_path, "FstLst_brmsBF.rds"))
}

# Plot jitter + mean&se + lines
generate_plots <- F
if(generate_plots){
  ## Plot per FstLst
  LT.pre_post.per_fstlst.plot <- ggplot(LT.pre_post.fstlst,
                                        aes(x = PrePost, y = Prop,
                                            colour = Condition,
                                            fill = Condition)) +
    theme(legend.pos = "top") + ylab("Looking to Tail (Prop)") + facet_grid(.~FstLst) +
    geom_point(position = position_jitterdodge(dodge.width = .8,
                                               jitter.width = .2),
               alpha = .25) +
    geom_errorbar(stat = "summary",
                  width = .2, colour = "black",
                  position = position_dodge(.1)) +
    geom_line(aes(x = PrePost, y = Prop, group = Condition),
              stat = "summary", fun.y = "mean",
              colour = "black",
              position = position_dodge(.1)) +
    geom_point(stat = "summary", fun.y = "mean",
               shape = 18, size = 3,
               position = position_dodge(.1))
  ggsave(paste0(save_path, "FstLst_data.pdf"),
         LT.pre_post.per_fstlst.plot,
         width = 7, height = 3)
}

# CONTRAST TEST ANALYSIS: PROP NEW FEATURE BY HEAD/TAIL CONTRAST TEST ==============================
save_path <- "../results/infants/OldNew/TrialAverage_"
# Prepare dataset
new_old <- LT.test.ctr %>%
  subset(ContrastType %in% c("Tail", "Head")) %>%
  make_time_window_data(aois = "NewFeature",
                        predictor_columns = c("Condition",
                                              "ContrastType")) %>%
  mutate(ChanceArcsin = ArcSin - asin(sqrt(.5)))
## Check for amount of data available
participants.new_old <- new_old %>%
  group_by(Participant, Condition) %>%
  summarise(nTrials = n_distinct(TrialId)) %>%
  group_by(Condition, nTrials) %>%
  summarise(Participants = n_distinct(Participant))
# Testing Prop ~ ContrastType*Condition
run_model <- T
if(run_model){
  t <- proc.time()
  ## Run lmer
  new_old.lmer.model <- lmer(ChanceArcsin ~ ContrastType*Condition +
                                  (1 | Participant),
                                data = new_old)
  new_old.lmer.anova <- anova(new_old.lmer.model, type = 1)
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
  new_old.time <- proc.time() - t
  ## Save all the results
  saveRDS(new_old.lmer.model, paste0(save_path, "lmerModel.rds"))
  saveRDS(new_old.lmer.anova, paste0(save_path, "lmerAnova.rds"))
  saveRDS(new_old.brms.models, paste0(save_path, "brmsModels.rds"))
  saveRDS(new_old.brms.bayes_factors, paste0(save_path, "brmsBF.rds"))
}else{
  ## Read all the results
  new_old.lmer.model <- readRDS(paste0(save_path, "lmerModel.rds"))
  new_old.lmer.anova <- readRDS(paste0(save_path, "lmerAnova.rds"))
  new_old.brms.models <- readRDS(paste0(save_path, "brmsModels.rds"))
  new_old.brms.bayes_factors <- readRDS(paste0(save_path, "brmsBF.rds"))
}
# Plot jitter + mean&se
generate_plots <- F
if(generate_plots){
  LT.new_old.plot.data <- ggplot(LT.new_old,
                            aes(x = ContrastType, y = Prop,
                                colour = Condition,
                                fill = Condition)) +
    theme(legend.pos = "top") + ylab("Looking to New Feature (Prop)") +
    geom_point(position = position_jitterdodge(dodge.width = .8,
                                               jitter.width = .2),
               alpha = .25) +
    geom_errorbar(stat = "summary",
                  width = .2, colour = "black",
                  position = position_dodge(.1)) +
    geom_point(stat = "summary", fun.y = "mean",
               shape = 18, size = 3,
               position = position_dodge(.1)) +
    geom_hline(yintercept = 0.5)
  ggsave(paste0(save_path, "data.pdf"),
         LT.new_old.plot.data,
         width = 7, height = 5.4)
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
generate_plots <- T
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
  subset_by_window(window_start_col = "LabelOnset",
                   window_end_col = "TrialEnd") %>%
  make_time_window_data(aois = "Target",
                        predictor_columns = "CategoryName") %>%
  mutate(ChanceArcsin = ArcSin - asin(sqrt(.5))) # Value centered on chance looking, useful for test
## Check for amount of data available
participants <- prop_target %>%
  group_by(Participant) %>%
  summarise(nTrials = n_distinct(TrialId))
# Testing in general
run_model <- T
if(run_model){
  t <- proc.time()
  ## Run lmer
  prop_target.lmer.model <- lmer(ChanceArcsin ~ 1 + (1 | Participant),
                                    data = prop_target)
  prop_target.lmer.null <- lmer(ChanceArcsin ~ 0 + (1 | Participant),
                                   data = prop_target)
  prop_target.lmer.anova <- anova(prop_target.lmer.null,
                                     prop_target.lmer.model)
  ## Run brms
  prop_target.brms.model <- brm(ChanceArcsin ~ 1 + (1 | Participant),
                                   data = prop_target,
                                   chains = 4, cores = 4, iter = 2000,
                                   prior = set_prior("uniform(-.8,.8)",
                                                     class = "Intercept"),
                                   control = list(adapt_delta = .999,
                                                  max_treedepth = 20),
                                   save_all_pars = T)
  prop_target.brms.null <- brm(ChanceArcsin ~ 0 + (1 | Participant),
                                   data = prop_target,
                                   chains = 4, cores = 4, iter = 2000,
                                   control = list(adapt_delta = .999,
                                                  max_treedepth = 20),
                                  save_all_pars = T)
  prop_target.brms.models <- list(prop_target.brms.null, prop_target.brms.model)
  prop_target.brms.bayes_factor <- bayes_factor(prop_target.brms.model,
                                                prop_target.brms.null)
  prop_target.time <- proc.time() - t
  ## Save all the results
  saveRDS(prop_target.lmer.model, paste0(save_path, "lmerModel.rds"))
  saveRDS(prop_target.lmer.anova, paste0(save_path, "lmerAnova.rds"))
  saveRDS(prop_target.brms.models, paste0(save_path, "brmsModels.rds"))
  saveRDS(prop_target.brms.bayes_factor, paste0(save_path, "brmsBF.rds"))
}else{
  ## Read all the results
  prop_target.lmer.model <- readRDS(paste0(save_path, "lmerModel.rds"))
  prop_target.lmer.anova <- readRDS(paste0(save_path, "lmerAnova.rds"))
  prop_target.brms.models <- readRDS(paste0(save_path, "brmsModels.rds"))
  prop_target.brms.bayes_factor <- readRDS(paste0(save_path, "brmsBF.rds"))
}

# Plot jitter + mean&se
generate_plots <- F
if(generate_plots){
  LT.prop_target.plot.data <- ggplot(LT.prop_target,
                                 aes(x = AOI, y = Prop)) +
    theme(legend.pos = "top") + ylab("Looking to New Feature (Prop)") +
    geom_jitter(alpha = .25) +
    geom_errorbar(stat = "summary",
                  width = .2, colour = "black") +
    geom_point(stat = "summary", fun.y = "mean",
               shape = 18, size = 3) +
    geom_hline(yintercept = 0.5)
  ggsave(paste0(save_path, "data.pdf"),
         LT.prop_target.plot.data,
         width = 7, height = 5.4)
}

# WORD LEARNING TEST ANALYSIS: PROP TARGET TIME COURSE FOR LABEL CONDITION  ========================
save_path <- "../results/infants/WordLearning/TimeCourse_"
# Data preparation
LT.prop_target.time_course <- LT.test.wl %>%
  subset_by_window(window_start_col = "LabelOnset",
                   window_end_col = "TrialEnd") %>%
  mutate(Chance = F) %>%
  make_time_sequence_data(time_bin_size = 50,
                          aois = "Target",
                          predictor_columns=c("Chance"),
                          summarize_by = "Participant")
LT.prop_target.time_course.chance <- LT.prop_target.time_course %>%
  mutate(Chance = T,
         Participant = paste0("Chance", Participant),
         Prop = .5)
LT.prop_target.time_course.chance_test <- rbind(LT.prop_target.time_course,
                                                LT.prop_target.time_course.chance) %>%
  mutate_at("Chance", parse_factor, levels = NULL)
# BOOTSTRAPPED CLUSTER-BASED PERMUTATION ANALYSIS
run_model <- T
if(run_model){
  t <- proc.time()
  ## Determine threshold based on alpha = .05 two-tailed
  num_sub = length(unique((LT.prop_target.time_course.chance_test$Participant)))
  threshold_t = qt(p = 1 - .05/2,
                   df = num_sub-1)
  ## Determine clusters
  LT.prop_target.time_cluster <- LT.prop_target.time_course.chance_test %>%
    make_time_cluster_data(predictor_column = "Chance",
                           aoi = "Target",
                           test = "t.test",
                           threshold = threshold_t)
  ## Run analysis
  LT.prop_target.time_cluster.analysis <- LT.prop_target.time_cluster %>%
    analyze_time_clusters(within_subj = F,
                          parallel = T)
  bcbp.time <- proc.time() - t
  ## Save clusters and analysis
  saveRDS(LT.prop_target.time_cluster, paste0(save_path, "bcbpClusters.rds"))
  saveRDS(LT.prop_target.time_cluster.analysis, paste0(save_path, "bcbpAnalysis.rds"))
}else{
  ## Read the results
  LT.prop_target.time_cluster <- readRDS(paste0(save_path, "bcbpClusters.rds"))
  LT.prop_target.time_cluster.analysis <- readRDS(paste0(save_path, "bcbpAnalysis.rds"))
}

# PLOT
generate_plots <- T
if(generate_plots){
  LT.prop_target.time_course.plot <- ggplot(LT.prop_target.time_course,
                                            aes(x = Time, y=Prop)) +
    xlab('Time in Trial') + ylab("Looking to Tail (Prop)") +
    theme(legend.position = "top") + ylim(0,1) +
    stat_summary(fun.y='mean', geom='line', linetype = '61') +
    stat_summary(fun.data=mean_se, geom='ribbon', alpha= .25, colour=NA) +
    geom_hline(yintercept = .5)
  ggsave(paste0(save_path, "data.pdf"),
         plot = LT.prop_target.time_course.plot,
         width = 3.5, height = 2.5)
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
            Condition = first(Condition))
# Testing Switches ~ Condition*FstLst
run_model <- T
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
  saveRDS(fam_switches.per_fstlst.brms.models, paste0(save_path, "brmsModels.rds"))
  saveRDS(fam_switches.per_fstlst.brms.bayes_factors, paste0(save_path, "brmsBF.rds"))
}else{
  ## Read all the results
  fam_switches.per_fstlst.glmer.model <- readRDS(paste0(save_path, "glmerModel.rds"))
  fam_switches.per_fstlst.brms.models <- readRDS(paste0(save_path, "brmsModels.rds"))
  fam_switches.per_fstlst.brms.bayes_factors <- readRDS(paste0(save_path, "brmsBF.rds"))
}

# Plotting boxplots
generate_plots <- F
if(generate_plots){
  fam_switches.per_fstlst.plot <- LT.fam_switches %>%
    drop_na(FstLst) %>%
    ggplot(aes(y = Switches,
               x = FstLst,
               fill = Condition)) +
    theme(legend.position = "top",
          axis.title.x = element_blank()) +
    ylab("Switches between AOIs") +
    geom_boxplot()
  ggsave(paste0(save_path, "data.pdf"),
         fam_switches.per_fstlst.plot,
         width = 3.5, height = 3.5)
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
  select(-Tail)
first_aoi <- first_look %>%
  group_by(Participant, TrialId) %>%
  arrange(FirstAOILook) %>%
  summarise(AOI = first(AOI),
            # Keep columns for analysis (only one value per trial per participant)
            Condition = first(Condition),
            FstLst = first(FstLst))
first_tail <- first_look %>%
  subset(AOI == "Tail", select = -AOI) %>% # Discards trials with no Tail look: is it okay?
  mutate(logFirstAOILook = log(FirstAOILook))

# Testing (First)AOI ~ Condition*FstLst
run_model <- T
if(run_model){
  t <- proc.time()
  ## Run (g)lmer
  first_aoi.per_fstlst.glmer.model <- glmer(AOI ~ FstLst*Condition +
                                              (1 + FstLst | Participant),
                                            data = first_aoi,
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
                                         first_aoi,
                                         priors.first_aoi.per_fstlst,
                                         family = bernoulli())
  first_aoi.per_fstlst.brms.models <- brms.results[[1]]
  first_aoi.per_fstlst.brms.bayes_factors <- brms.results[[2]]
  first_aoi.time <- proc.time() - t
  ## Save all the results
  saveRDS(first_aoi.per_fstlst.glmer.model, paste0(save_path, "FirstAOI_glmerModel.rds"))
  saveRDS(first_aoi.per_fstlst.brms.models, paste0(save_path, "FirstAOI_brmsModels.rds"))
  saveRDS(first_aoi.per_fstlst.brms.bayes_factors, paste0(save_path, "FirstAOI_brmsBF.rds"))
}else{
  ## Read all the results
  first_aoi.per_fstlst.glmer.model <- readRDS(paste0(save_path, "FirstAOI_glmerModel.rds"))
  first_aoi.per_fstlst.brms.models <- readRDS(paste0(save_path, "FirstAOI_brmsModels.rds"))
  first_aoi.per_fstlst.brms.bayes_factors <- readRDS(paste0(save_path, "FirstAOI_brmsBF.rds"))
}
# Testing FirstAOI(Tail)Look ~ Condition*FstLst
run_model <- T
if(run_model){
  t <- proc.time()
  ## Run lmer
  first_tail.per_fstlst.lmer.model <- lmer(logFirstAOILook ~ FstLst*Condition +
                                               (1 + FstLst | Participant),
                                             data = first_tail)
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
                                         first_tail,
                                         priors.first_tail.per_fstlst)
  first_tail.per_fstlst.brms.models <- brms.results[[1]]
  first_tail.per_fstlst.brms.bayes_factors <- brms.results[[2]]
  first_tail.time <- proc.time() - t
  ## Save all the results
  saveRDS(first_tail.per_fstlst.lmer.model, paste0(save_path, "FirstTail_lmerModel.rds"))
  saveRDS(first_tail.per_fstlst.lmer.anova, paste0(save_path, "FirstTail_lmerAnova.rds"))
  saveRDS(first_tail.per_fstlst.brms.models, paste0(save_path, "FirstTail_brmsModels.rds"))
  saveRDS(first_tail.per_fstlst.brms.bayes_factors, paste0(save_path, "FirstTail_brmsBF.rds"))
}else{
  ## Read all the results
  first_tail.per_fstlst.lmer.model <- readRDS(paste0(save_path, "FirstTail_lmerModel.rds"))
  first_tail.per_fstlst.lmer.anova <- readRDS(paste0(save_path, "FirstTail_lmerAnova.rds"))
  first_tail.per_fstlst.brms.models <- readRDS(paste0(save_path, "FirstTail_brmsModels.rds"))
  first_tail.per_fstlst.brms.bayes_factors <- readRDS(paste0(save_path, "FirstTail_brmsBF.rds"))
}

# Plotting
generate_plots <- F
if(generate_plots){
  ## First AOI (boxplot)
  first_aoi.per_fstlst.plot <- LT.first_aoi %>%
    drop_na(FstLst) %>%
    group_by(Participant, FstLst, AOI) %>%
    summarise(N = n(),
              Condition = first(Condition)) %>%
    ggplot(aes(y = N,
               x = AOI,
               fill = Condition)) +
    theme(legend.position = "top",
          axis.title.x = element_blank()) +
    ylab("Number of first look to AOI") +
    facet_grid(.~FstLst) +
    geom_boxplot()
  ggsave(paste0(save_path, "FirstAOI_data.pdf"),
         first_aoi.per_fstlst.plot,
         width = 7, height = 3.5)
  ## Time to first tail look (boxplot)
  first_tail.per_fstlst.plot <- LT.first_tail %>%
    drop_na(FstLst) %>%
    ggplot(aes(y = FirstAOILook,
               x = FstLst,
               fill = Condition)) +
    theme(legend.position = "top",
          axis.title.x = element_blank()) +
    ylab("Time to first Tail look") +
    geom_boxplot()
  ggsave(paste0(save_path, "FirstTail_data.pdf"),
         first_tail.per_fstlst.plot,
         width = 3.5, height = 3.5)
}
