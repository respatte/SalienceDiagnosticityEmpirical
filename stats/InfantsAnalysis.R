# LIBRARY IMPORTS ==================================================================================
library(eyetrackingR)
library(lme4)
library(lmerTest)
library(brms)
library(tidyverse)

source("Routines.R")

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
LT.prop_tail <- LT.fam %>%
  subset_by_window(window_start_col = "LabelOnset") %>%
  make_time_window_data(aois=c("Tail"),
                        predictor_columns=c("Condition",
                                            "TrialId",
                                            "TrialNum",
                                            "FamPart",
                                            "FstLst",
                                            "Stimulus",
                                            "CategoryName"))
# Testing Prop ~ Trial*Condition
run_model <- F
if(run_model){
  ## Run lmer (Sampling Theory Based)
  LT.prop_tail.per_trial.lmer.model <- lmer(ArcSin ~ TrialNum*Condition +
                                              (1 | Participant) +
                                              (1 | Stimulus),
                                            data = LT.prop_tail)
  LT.prop_tail.per_trial.lmer.anova <- anova(LT.prop_tail.per_trial.lmer.model, type = 1)
  ## Run brms (Bayesian)
  prior.prop_tail.per_trial <- c(set_prior("uniform(0,1.6)",
                                           class = "Intercept"),
                                 set_prior("normal(0,.5)", class = "b"))
  LT.prop_tail.per_trial.brms.model <- brm(ArcSin ~ TrialNum*Condition +
                                             (1 | Participant) +
                                             (1 | Stimulus),
                                           prior = prior.prop_tail.per_trial,
                                           data = LT.prop_tail,
                                           chains = 4, cores = 4)
  ## Save all results
  saveRDS(LT.prop_tail.per_trial.lmer.model, paste0(save_path, "Trials_lmerModel.rds"))
  saveRDS(LT.prop_tail.per_trial.lmer.anova, paste0(save_path, "Trials_lmerAnova.rds"))
  saveRDS(LT.prop_tail.per_trial.brms.model, paste0(save_path, "Trials_brmsModel.rds"))
}else{
  ## Read all the results
  LT.prop_tail.per_trial.lmer.model <- readRDS(paste0(save_path, "Trials_lmerModel.rds"))
  LT.prop_tail.per_trial.lmer.anova <- readRDS(paste0(save_path, "Trials_lmerAnova.rds"))
  LT.prop_tail.per_trial.brms.model <- readRDS(paste0(save_path, "Trials_brmsModel.rds"))
}
# Testing Prop ~ Part*Condition
run_model <- F
if(run_model){
  ## Run lmer
  LT.prop_tail.per_part.lmer.model <- lmer(ArcSin ~ FamPart*Condition +
                                             (1 + FamPart | Participant) +
                                             (1 | Stimulus),
                                           data = LT.prop_tail)
  LT.prop_tail.per_part.lmer.anova <- anova(LT.prop_tail.per_part.lmer.model, type = 1)
  ## Run brms
  prior.prop_tail.per_part <- c(set_prior("uniform(0,1.6)",
                                          class = "Intercept"),
                                set_prior("normal(0,.5)", class = "b"))
  LT.prop_tail.per_part.brms.model <- brm(ArcSin ~ FamPart*Condition +
                                            (1 + FamPart | Participant) +
                                            (1 | Stimulus),
                                          data = LT.prop_tail,
                                          prior = prior.prop_tail.per_part,
                                          chains = 4, cores = 4)
  ## Save all the results
  saveRDS(LT.prop_tail.per_part.lmer.model, "../results/infants/PropTail/Part_lmerModel.rds")
  saveRDS(LT.prop_tail.per_part.lmer.anova, "../results/infants/PropTail/Part_lmerAnova.rds")
  saveRDS(LT.prop_tail.per_part.brms.model, "../results/infants/PropTail/Part_brmsModel.rds")
}else{
  ## Read all the results
  LT.prop_tail.per_part.lmer.model <- readRDS("../results/infants/PropTail/Part_lmerModel.rds")
  LT.prop_tail.per_part.lmer.anova <- readRDS("../results/infants/PropTail/Part_lmerAnova.rds")
  LT.prop_tail.per_part.brms.model <- readRDS("../results/infants/PropTail/Part_brmsModel.rds")
}
# Testing Prop ~ FstLst*Condition
run_model <- F
if(run_model){
  ## Select data
  LT.prop_tail.fstlst <- LT.prop_tail %>%
    drop_na(FstLst)
  ## Run lmer
  LT.prop_tail.per_fstlst.lmer.model <- lmer(ArcSin ~ FstLst*Condition +
                                               (1 + FstLst | Participant) +
                                               (1 | Stimulus),
                                             data = LT.prop_tail.fstlst)
  LT.prop_tail.per_fstlst.lmer.anova <- anova(LT.prop_tail.per_fstlst.lmer.model, type = 1)
  ## Run brms
  prior.prop_tail.per_fstlst <- c(set_prior("uniform(0,1.6)",
                                           class = "Intercept"),
                                 set_prior("normal(0,.5)", class = "b"))
  LT.prop_tail.per_fstlst.brms.model.3 <- brm(ArcSin ~ FstLst*Condition +
                                                (1 + FstLst | Participant) +
                                                (1 | Stimulus),
                                              data = LT.prop_tail.fstlst,
                                              prior = prior.prop_tail.per_fstlst,
                                              chains = 4, cores = 4,
                                              save_all_pars = T)
  LT.prop_tail.per_fstlst.brms.model.2 <- brm(ArcSin ~ FstLst + Condition +
                                                (1 + FstLst | Participant) +
                                                (1 | Stimulus),
                                              data = LT.prop_tail.fstlst,
                                              prior = prior.prop_tail.per_fstlst,
                                              chains = 4, cores = 4,
                                              save_all_pars = T)
  LT.prop_tail.per_fstlst.brms.model.1 <- brm(ArcSin ~ FstLst +
                                                (1 + FstLst | Participant) +
                                                (1 | Stimulus),
                                              data = LT.prop_tail.fstlst,
                                              prior = prior.prop_tail.per_fstlst,
                                              chains = 4, cores = 4,
                                              save_all_pars = T)
  LT.prop_tail.per_fstlst.brms.model.0 <- brm(ArcSin ~ 1 +
                                                (1 | Participant) +
                                                (1 | Stimulus),
                                              data = LT.prop_tail.fstlst,
                                              prior = set_prior("uniform(0,1.6)",
                                                                class = "Intercept"),
                                              chains = 4, cores = 4,
                                              save_all_pars = T)
  LT.prop_tail.per_fstlst.brms.bf.3_2 <- bayes_factor(LT.prop_tail.per_fstlst.brms.model.3,
                                                      LT.prop_tail.per_fstlst.brms.model.2)
  LT.prop_tail.per_fstlst.brms.bf.2_1 <- bayes_factor(LT.prop_tail.per_fstlst.brms.model.2,
                                                      LT.prop_tail.per_fstlst.brms.model.1)
  LT.prop_tail.per_fstlst.brms.bf.1_0 <- bayes_factor(LT.prop_tail.per_fstlst.brms.model.1,
                                                      LT.prop_tail.per_fstlst.brms.model.0)
  LT.prop_tail.per_fstlst.brms.bayes_factors <- list(LT.prop_tail.per_fstlst.brms.bf.1_0,
                                                     LT.prop_tail.per_fstlst.brms.bf.2_1,
                                                     LT.prop_tail.per_fstlst.brms.bf.3_2)
  ## Save all the results
  saveRDS(LT.prop_tail.per_fstlst.lmer.model, "../results/infants/PropTail/FstLst_lmerModel.rds")
  saveRDS(LT.prop_tail.per_fstlst.lmer.anova, "../results/infants/PropTail/FstLst_lmerAnova.rds")
  saveRDS(LT.prop_tail.per_fstlst.brms.model.3, "../results/infants/PropTail/FstLst_brmsModel.rds")
  saveRDS(LT.prop_tail.per_fstlst.brms.bayes_factors,
          "../results/infants/PropTail/FstLst_brmsBF.rds")
}else{
  ## Read all the results
  LT.prop_tail.per_fstlst.lmer.model <- readRDS("../results/infants/PropTail/FstLst_lmerModel.rds")
  LT.prop_tail.per_fstlst.lmer.anova <- readRDS("../results/infants/PropTail/FstLst_lmerAnova.rds")
  LT.prop_tail.per_fstlst.brms.model.3 <- readRDS("../results/infants/PropTail/FstLst_brmsModel.rds")
  LT.prop_tail.per_fstlst.brms.bayes_factors <- readRDS("../results/infants/PropTail/FstLst_brmsBF.rds")
}

# Plot jitter + mean&se + lines
generate_plots <- F
if(generate_plots){
  ## Plot per trial
  LT.prop_tail.per_trial.plot <- ggplot(LT.prop_tail,
                                        aes(x = TrialId, y = Prop,
                                            colour = Condition,
                                            fill = Condition)) +
    theme(legend.pos = "top") + ylab("Looking to Tail (Prop)") +
    geom_point(position = position_jitterdodge(dodge.width = .8,
                                               jitter.width = .2),
               alpha = .25) +
    geom_errorbar(stat = "summary",
                  width = .2, colour = "black",
                  position = position_dodge(.1)) +
    geom_line(aes(x = TrialId, y = Prop, group = Condition),
              stat = "summary", fun.y = "mean",
              colour = "black",
              position = position_dodge(.1)) +
    geom_point(stat = "summary", fun.y = "mean",
               shape = 18, size = 3,
               position = position_dodge(.1))
  ggsave("../results/infants/PropTail/TrialAverage_Trials.pdf",
         LT.prop_tail.per_trial.plot,
         width = 7, height = 5.4)
  ## Plot per part
  LT.prop_tail.per_part.plot <- ggplot(LT.prop_tail,
                                       aes(x = FamPart, y = Prop,
                                           colour = Condition,
                                           fill = Condition)) +
    theme(legend.pos = "top") + ylab("Looking to Tail (Prop)") +
    geom_point(position = position_jitterdodge(dodge.width = .8,
                                               jitter.width = .2),
               alpha = .25) +
    geom_errorbar(stat = "summary",
                  width = .2, colour = "black",
                  position = position_dodge(.1)) +
    geom_line(aes(x = FamPart, y = Prop, group = Condition),
              stat = "summary", fun.y = "mean",
              colour = "black",
              position = position_dodge(.1)) +
    geom_point(stat = "summary", fun.y = "mean",
               shape = 18, size = 3,
               position = position_dodge(.1))
  ggsave("../results/infants/PropTail/TrialAverage_Parts.pdf",
         LT.prop_tail.per_part.plot,
         width = 7, height = 5.4)
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
  ggsave("../results/infants/PropTail/TrialAverage_FstLst.pdf",
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
# GROWTH CURVE ANALYSIS
## TODO?
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
LT.pre_post <- LT.fam %>%
  make_time_window_data(aois=c("Tail"),
                        predictor_columns=c("Condition",
                                            "TrialId",
                                            "TrialNum",
                                            "FamPart",
                                            "FstLst",
                                            "PrePost",
                                            "Stimulus",
                                            "CategoryName")) %>%
  drop_na(PrePost)
# Testing Prop ~ Trial*PrePost*Condition
run_model <- F
if(run_model){
  ## Run lmer (Sampling Theory Based)
  LT.pre_post.per_trial.lmer.model <- lmer(ArcSin ~ TrialNum*PrePost*Condition +
                                             (1 | Participant) +
                                             (1 | Stimulus),
                                           data = LT.pre_post)
  LT.pre_post.per_trial.lmer.anova <- anova(LT.prop_tail.per_trial.lmer.model, type = 1)
  ## Run brms (Bayesian)
  prior.pre_post.per_trial <- c(set_prior("uniform(0,1.6)",
                                          class = "Intercept"),
                                set_prior("normal(0,.5)", class = "b"))
  LT.pre_post.per_trial.brms.model <- brm(ArcSin ~ TrialNum*PrePost*Condition +
                                            (1 | Participant) +
                                            (1 | Stimulus),
                                          data = LT.pre_post,
                                          prior = prior.pre_post.per_trial,
                                          chains = 4, cores = 4)
  ## Save all results
  saveRDS(LT.pre_post.per_trial.lmer.model, paste0(save_path, "Trials_lmerModel.rds"))
  saveRDS(LT.pre_post.per_trial.lmer.anova, paste0(save_path, "Trials_lmerAnova.rds"))
  saveRDS(LT.pre_post.per_trial.brms.model, paste0(save_path, "Trials_brmsModel.rds"))
}else{
  ## Read all the results
  LT.pre_post.per_trial.lmer.model <- readRDS(paste0(save_path, "Trials_lmerModel.rds"))
  LT.pre_post.per_trial.lmer.anova <- readRDS(paste0(save_path, "Trials_lmerAnova.rds"))
  LT.pre_post.per_trial.brms.model <- readRDS(paste0(save_path, "Trials_brmsModel.rds"))
}
# Testing Prop ~ Part*PrePost*Condition
run_model <- F
if(run_model){
  ## Run lmer
  LT.pre_post.per_part.lmer.model <- lmer(ArcSin ~ FamPart*PrePost*Condition +
                                            (1 + FamPart | Participant) +
                                            (1 | Stimulus),
                                          data = LT.pre_post)
  LT.pre_post.per_part.lmer.anova <- anova(LT.pre_post.per_part.lmer.model, type = 1)
  ## Run brms
  prior.pre_post.per_part <- c(set_prior("uniform(0,1.6)",
                                         class = "Intercept"),
                               set_prior("normal(0,.5)", class = "b"))
  LT.pre_post.per_part.brms.model <- brm(ArcSin ~ FamPart*PrePost*Condition +
                                           (1 + FamPart | Participant) +
                                           (1 | Stimulus),
                                         data = LT.pre_post,
                                         prior = prior.pre_post.per_part,
                                         chains = 4, cores = 4)
  ## Save all the results
  saveRDS(LT.pre_post.per_part.lmer.model, paste0(save_path, "Parts_lmerModel.rds"))
  saveRDS(LT.pre_post.per_part.lmer.anova, paste0(save_path, "Parts_lmerAnova.rds"))
  saveRDS(LT.pre_post.per_part.brms.model, paste0(save_path, "Parts_brmsModel.rds"))
}else{
  ## Read all the results
  LT.pre_post.per_part.lmer.model <- readRDS(paste0(save_path, "Parts_lmerModel.rds"))
  LT.pre_post.per_part.lmer.anova <- readRDS(paste0(save_path, "Parts_lmerAnova.rds"))
  LT.pre_post.per_part.brms.model <- readRDS(paste0(save_path, "Parts_brmsModel.rds"))
}
# Testing Prop ~ FstLst*Condition
run_model <- F
if(run_model){
  ## Select data
  LT.pre_post.fstlst <- LT.pre_post %>%
    drop_na(FstLst)
  ## Run lmer
  LT.pre_post.per_fstlst.lmer.model <- lmer(ArcSin ~ FstLst*PrePost*Condition +
                                              (1 + FstLst | Participant) +
                                              (1 | Stimulus),
                                            data = LT.pre_post.fstlst)
  LT.pre_post.per_fstlst.lmer.anova <- anova(LT.pre_post.per_fstlst.lmer.model, type = 1)
  ## Run brms
  prior.pre_post.per_fstlst <- c(set_prior("uniform(0,1.6)",
                                           class = "Intercept"),
                                 set_prior("normal(0,.5)", class = "b"))
  LT.pre_post.per_fstlst.brms.model <- brm(ArcSin ~ FstLst*PrePost*Condition +
                                             (1 + FstLst | Participant) +
                                             (1 | Stimulus),
                                           data = LT.pre_post.fstlst,
                                           prior = prior.pre_post.per_fstlst,
                                           chains = 4, cores = 4)
  ## Save all the results
  saveRDS(LT.pre_post.per_fstlst.lmer.model, paste0(save_path, "FstLst_lmerModel.rds"))
  saveRDS(LT.pre_post.per_fstlst.lmer.anova, paste0(save_path, "FstLst_lmerAnova.rds"))
  saveRDS(LT.pre_post.per_fstlst.brms.model, paste0(save_path, "FstLst_brmsModel.rds"))
}else{
  ## Read all the results
  LT.pre_post.per_fstlst.lmer.model <- readRDS(paste0(save_path, "FstLst_lmerModel.rds"))
  LT.pre_post.per_fstlst.lmer.anova <- readRDS(paste0(save_path, "FstLst_lmerAnova.rds"))
  LT.pre_post.per_fstlst.brms.model <- readRDS(paste0(save_path, "FstLst_brmsModel.rds"))
}

# Plot jitter + mean&se + lines
generate_plots <- F
if(generate_plots){
  ## Plot per part
  LT.pre_post.per_part.plot <- ggplot(LT.pre_post,
                                      aes(x = PrePost, y = Prop,
                                          colour = Condition,
                                          fill = Condition)) +
    theme(legend.pos = "top") + ylab("Looking to Tail (Prop)") + facet_grid(.~FamPart) +
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
  ggsave(paste0(save_path, "Parts_data.pdf"),
         LT.pre_post.per_part.plot,
         width = 7, height = 3)
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
LT.new_old <- LT.test.ctr %>%
  subset(ContrastType %in% c("Tail", "Head")) %>%
  make_time_window_data(aois = "NewFeature",
                        predictor_columns = c("Condition",
                                              "ContrastType")) %>%
  mutate(ChanceArcsin = ArcSin - asin(sqrt(.5)))
## Check for amount of data available
participants.new_old <- LT.new_old %>%
  group_by(Participant, Condition) %>%
  summarise(nTrials = n_distinct(TrialId)) %>%
  group_by(Condition, nTrials) %>%
  summarise(Participants = n_distinct(Participant))
# Testing Prop ~ ContrastType*Condition
run_model <- F
if(run_model){
  ## Run lmer
  LT.new_old.lmer.model <- lmer(ChanceArcsin ~ ContrastType*Condition +
                                  (1 | Participant),
                                data = LT.new_old)
  LT.new_old.lmer.anova <- anova(LT.new_old.lmer.model, type = 1)
  ## Run brm
  prior.new_old <- c(set_prior("uniform(-0.8,0.8)",
                               class = "Intercept"),
                     set_prior("normal(0,.5)", class = "b"))
  LT.new_old.brms.model.4 <- brm(ChanceArcsin ~ ContrastType*Condition +
                                   (1 | Participant),
                                 data = LT.new_old,
                                 prior = prior.new_old,
                                 chains = 4, cores = 4,
                                 save_all_pars = T)
  LT.new_old.brms.model.3 <- brm(ChanceArcsin ~ ContrastType + Condition +
                                   (1 | Participant),
                                 data = LT.new_old,
                                 prior = prior.new_old,
                                 chains = 4, cores = 4,
                                 save_all_pars = T)
  LT.new_old.brms.model.2 <- brm(ChanceArcsin ~ ContrastType +
                                   (1 | Participant),
                                 data = LT.new_old,
                                 prior = prior.new_old,
                                 chains = 4, cores = 4,
                                 save_all_pars = T)
  LT.new_old.brms.model.1 <- brm(ChanceArcsin ~ 1 +
                                   (1 | Participant),
                                 data = LT.new_old,
                                 prior = set_prior("uniform(-.8,.8)",
                                                   class = "Intercept"),
                                 chains = 4, cores = 4,
                                 save_all_pars = T)
  LT.new_old.brms.model.0 <- brm(ChanceArcsin ~ 0 +
                                   (1 | Participant),
                                 data = LT.new_old,
                                 chains = 4, cores = 4,
                                 save_all_pars = T)
  LT.new_old.brms.bf.4_3 <- bayes_factor(LT.new_old.brms.model.4,
                                         LT.new_old.brms.model.3)
  LT.new_old.brms.bf.3_2 <- bayes_factor(LT.new_old.brms.model.3,
                                         LT.new_old.brms.model.2)
  LT.new_old.brms.bf.2_1 <- bayes_factor(LT.new_old.brms.model.2,
                                         LT.new_old.brms.model.1)
  LT.new_old.brms.bf.1_0 <- bayes_factor(LT.new_old.brms.model.1,
                                         LT.new_old.brms.model.0)
  LT.new_old.brms.bayes_factors <- list(LT.new_old.brms.bf.1_0,
                                        LT.new_old.brms.bf.2_1,
                                        LT.new_old.brms.bf.3_2,
                                        LT.new_old.brms.bf.4_3)
  ## Save all the results
  saveRDS(LT.new_old.lmer.model, paste0(save_path, "lmerModel.rds"))
  saveRDS(LT.new_old.lmer.anova, paste0(save_path, "lmerAnova.rds"))
  saveRDS(LT.new_old.brms.model.4, paste0(save_path, "brmsModel.rds"))
  saveRDS(LT.new_old.brms.bayes_factors, paste0(save_path, "brmsBF.rds"))
}else{
  ## Read all the results
  LT.new_old.lmer.model <- readRDS(paste0(save_path, "lmerModel.rds"))
  LT.new_old.lmer.anova <- readRDS(paste0(save_path, "lmerAnova.rds"))
  LT.new_old.brms.model.4 <- readRDS(paste0(save_path, "brmsModel.rds"))
  LT.new_old.brms.bayes_factors <- readRDS(paste0(save_path, "brmsBF.rds"))
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
# Run growth curve analysis
#### NOT ENOUGH DATA
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
LT.prop_target <- LT.test.wl %>%
  subset(Condition == "Label") %>%
  subset_by_window(window_start_col = "LabelOnset",
                   window_end_col = "TrialEnd") %>%
  make_time_window_data(aois = "Target",
                        predictor_columns = "CategoryName") %>%
  mutate(ChanceArcsin = ArcSin - asin(sqrt(.5))) # Value centered on chance looking, useful for test
## Check for amount of data available
participants <- LT.prop_target %>%
  group_by(Participant) %>%
  summarise(nTrials = n_distinct(TrialId))
# Testing in general
run_model <- F
if(run_model){
  ## Run lmer
  LT.prop_target.lmer.model <- lmer(ChanceArcsin ~ 1 + (1 | Participant),
                                    data = LT.prop_target)
  LT.prop_target.lmer.null <- lmer(ChanceArcsin ~ 0 + (1 | Participant),
                                   data = LT.prop_target)
  LT.prop_target.lmer.anova <- anova(LT.prop_target.lmer.null,
                                     LT.prop_target.lmer.model)
  ## Run brms
  LT.prop_target.brms.model <- brm(ChanceArcsin ~ 1 + (1 | Participant),
                                   data = LT.prop_target,
                                   chains = 4, cores = 4, iter = 2000,
                                   prior = set_prior("uniform(-.8,.8)",
                                                     class = "Intercept"),
                                   control = list(adapt_delta = .999,
                                                  max_treedepth = 20),
                                   save_all_pars = T)
  LT.prop_target.brms.null <- brm(ChanceArcsin ~ 0 + (1 | Participant),
                                   data = LT.prop_target,
                                   chains = 4, cores = 4, iter = 2000,
                                   control = list(adapt_delta = .999,
                                                  max_treedepth = 20),
                                  save_all_pars = T)
  LT.prop_target.brms.bayes_factor <- bayes_factor(LT.prop_target.brms.model,
                                                   LT.prop_target.brms.null)
  ## Save all the results
  saveRDS(LT.prop_target.lmer.model, paste0(save_path, "lmerModel.rds"))
  saveRDS(LT.prop_target.lmer.anova, paste0(save_path, "lmerAnova.rds"))
  saveRDS(LT.prop_target.brms.model, paste0(save_path, "brmsModel.rds"))
  saveRDS(LT.prop_target.brms.bayes_factor, paste0(save_path, "brmsBF.rds"))
}else{
  ## Read all the results
  LT.prop_target.lmer.model <- readRDS(paste0(save_path, "lmerModel.rds"))
  LT.prop_target.lmer.anova <- readRDS(paste0(save_path, "lmerAnova.rds"))
  LT.prop_target.brms.model <- readRDS(paste0(save_path, "brmsModel.rds"))
  LT.prop_target.brms.bayes_factor <- readRDS(paste0(save_path, "brmsBF.rds"))
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
# GROWTH CURVE ANALYSIS
## TODO?
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
LT.fam_switches <- LT.fam %>%
  drop_na(Tail) %>%
  group_by(Participant, TrialId) %>%
  summarise(Switches = sum(Tail != lag(Tail), na.rm = T), # Count switches per trial per participant
            # Keep columns for analysis (only one value per trial per participant)
            FstLst = first(FstLst),
            Condition = first(Condition))
# Testing Switches ~ Condition*FstLst
run_model <- T
if(run_model){
  fam_switches.per_fstlst.glmer.model <- glmer(Switches ~ FstLst*Condition +
                                                 (1 + FstLst | Participant),
                                               data = LT.fam_switches,
                                               family = poisson())
  # Current p-values from summary may not be the best. Do something else?
  ## Run brms
  prior.fam_switches.per_fstlst <- c(set_prior("normal(0,.5)", class = "b"))
  fam_switches.per_fstlst.brms.model.3 <- brm(Switches ~ FstLst*Condition +
                                                (1 + FstLst | Participant),
                                              data = LT.fam_switches,
                                              prior = prior.fam_switches.per_fstlst,
                                              family = poisson(),
                                              chains = 4, cores = 4, iter = 4000,
                                              control = list(adapt_delta = .95),
                                              save_all_pars = T)
  fam_switches.per_fstlst.brms.model.2 <- brm(Switches ~ FstLst + Condition +
                                                (1 + FstLst | Participant),
                                              data = LT.fam_switches,
                                              prior = prior.fam_switches.per_fstlst,
                                              family = poisson(),
                                              chains = 4, cores = 4, iter = 4000,
                                              control = list(adapt_delta = .95),
                                              save_all_pars = T)
  fam_switches.per_fstlst.brms.model.1 <- brm(Switches ~ FstLst +
                                                (1 + FstLst | Participant),
                                              data = LT.fam_switches,
                                              prior = prior.fam_switches.per_fstlst,
                                              family = poisson(),
                                              chains = 4, cores = 4, iter = 4000,
                                              control = list(adapt_delta = .95),
                                              save_all_pars = T)
  fam_switches.per_fstlst.brms.model.0 <- brm(Switches ~ 1 +
                                                (1 | Participant),
                                              data = LT.fam_switches,
                                              family = poisson(),
                                              chains = 4, cores = 4, iter = 4000,
                                              control = list(adapt_delta = .95),
                                              save_all_pars = T)
  fam_switches.per_fstlst.brms.bf.3_2 <- bayes_factor(fam_switches.per_fstlst.brms.model.3,
                                                      fam_switches.per_fstlst.brms.model.2)
  fam_switches.per_fstlst.brms.bf.2_1 <- bayes_factor(fam_switches.per_fstlst.brms.model.2,
                                                      fam_switches.per_fstlst.brms.model.1)
  fam_switches.per_fstlst.brms.bf.1_0 <- bayes_factor(fam_switches.per_fstlst.brms.model.1,
                                                      fam_switches.per_fstlst.brms.model.0)
  fam_switches.per_fstlst.brms.bayes_factors <- list(fam_switches.per_fstlst.brms.bf.1_0,
                                                     fam_switches.per_fstlst.brms.bf.2_1,
                                                     fam_switches.per_fstlst.brms.bf.3_2)
  ## Save all the results
  saveRDS(fam_switches.per_fstlst.glmer.model, paste0(save_path, "glmerModel.rds"))
  saveRDS(fam_switches.per_fstlst.brms.model.3, paste0(save_path, "brmsModel.rds"))
  saveRDS(fam_switches.per_fstlst.brms.bayes_factors, paste0(save_path, "brmsBF.rds"))
}else{
  ## Read all the results
  fam_switches.per_fstlst.glmer.model <- readRDS(paste0(save_path, "glmerModel.rds"))
  fam_switches.per_fstlst.brms.model.3 <- readRDS(paste0(save_path, "brmsModel.rds"))
  fam_switches.per_fstlst.brms.bayes_factors <- readRDS(paste0(save_path, "brmsBF.rds"))
}

# FAMILIARISATION: FIRST LOOK ======================================================================
save_path <- "../results/infants/FirstLook/"
# Prepare datasets
LT.first_look <- LT.fam %>%
  drop_na(Tail) %>%
  group_by(Participant, TrialId, Tail) %>%
  summarise(FirstAOILook = first(TimeStamp), # First time to each AOI
            AOI = parse_factor(ifelse(first(Tail), "Tail", "Head"),
                               levels = c("Head", "Tail")), # Explicit AOI flag
            # Keep columns for analysis (only one value per trial per participant)
            Condition = first(Condition),
            FstLst = first(FstLst)) %>%
  select(-Tail)
LT.first_aoi <- LT.first_look %>%
  group_by(Participant, TrialId) %>%
  arrange(FirstAOILook) %>%
  summarise(AOI = first(AOI),
            # Keep columns for analysis (only one value per trial per participant)
            Condition = first(Condition),
            FstLst = first(FstLst))
LT.first_tail <- LT.first_look %>%
  subset(AOI == "Tail", select = -AOI) # Discards trials with no Tail look: is it okay?

# Testing (First)AOI ~ Condition*FstLst
run_model <- F
if(run_model){
  first_aoi.per_fstlst.glmer.model <- glmer(AOI ~ FstLst*Condition +
                                              (1 + FstLst | Participant),
                                            data = LT.first_aoi,
                                            family = binomial())
  # Current p-values from summary may not be the best. Do something else?
  ## Run brms
  prior.first_aoi.per_fstlst <- c(set_prior("normal(0,.5)", class = "b"))
  first_aoi.per_fstlst.brms.model.3 <- brm(AOI ~ FstLst*Condition +
                                             (1 + FstLst | Participant),
                                           data = LT.first_aoi,
                                           prior = prior.first_aoi.per_fstlst,
                                           family = bernoulli(),
                                           chains = 4, cores = 4,
                                           save_all_pars = T)
  first_aoi.per_fstlst.brms.model.2 <- brm(AOI ~ FstLst + Condition +
                                             (1 + FstLst | Participant),
                                           data = LT.first_aoi,
                                           prior = prior.first_aoi.per_fstlst,
                                           family = bernoulli(),
                                           chains = 4, cores = 4,
                                           save_all_pars = T)
  first_aoi.per_fstlst.brms.model.1 <- brm(AOI ~ FstLst +
                                             (1 + FstLst | Participant),
                                           data = LT.first_aoi,
                                           prior = prior.first_aoi.per_fstlst,
                                           family = bernoulli(),
                                           chains = 4, cores = 4,
                                           save_all_pars = T)
  first_aoi.per_fstlst.brms.model.0 <- brm(AOI ~ 1 +
                                             (1 | Participant),
                                           data = LT.first_aoi,
                                           family = bernoulli(),
                                           chains = 4, cores = 4,
                                           save_all_pars = T)
  first_aoi.per_fstlst.brms.bf.3_2 <- bayes_factor(first_aoi.per_fstlst.brms.model.3,
                                                   first_aoi.per_fstlst.brms.model.2)
  first_aoi.per_fstlst.brms.bf.2_1 <- bayes_factor(first_aoi.per_fstlst.brms.model.2,
                                                   first_aoi.per_fstlst.brms.model.1)
  first_aoi.per_fstlst.brms.bf.1_0 <- bayes_factor(first_aoi.per_fstlst.brms.model.1,
                                                   first_aoi.per_fstlst.brms.model.0)
  first_aoi.per_fstlst.brms.bayes_factors <- list(first_aoi.per_fstlst.brms.bf.1_0,
                                                  first_aoi.per_fstlst.brms.bf.2_1,
                                                  first_aoi.per_fstlst.brms.bf.3_2)
  ## Save all the results
  saveRDS(first_aoi.per_fstlst.glmer.model, paste0(save_path, "FirstAOI_glmerModel.rds"))
  saveRDS(first_aoi.per_fstlst.brms.model.3, paste0(save_path, "FirstAOI_brmsModel.rds"))
  saveRDS(first_aoi.per_fstlst.brms.bayes_factors, paste0(save_path, "FirstAOI_brmsBF.rds"))
}else{
  ## Read all the results
  first_aoi.per_fstlst.glmer.model <- readRDS(paste0(save_path, "FirstAOI_glmerModel.rds"))
  first_aoi.per_fstlst.brms.model.3 <- readRDS(paste0(save_path, "FirstAOI_brmsModel.rds"))
  first_aoi.per_fstlst.brms.bayes_factors <- readRDS(paste0(save_path, "FirstAOI_brmsBF.rds"))
}
