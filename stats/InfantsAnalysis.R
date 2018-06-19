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
  subset(!(Participant %in% c("P03","P53","P55")))
gender <- d[[4]] %>%
  group_by(Gender, Condition) %>%
  summarise(N = n_distinct(Participant))
age <- d[[4]] %>%
  group_by(Participant, Condition) %>%
  summarise(Age = first(Age))
# Creating general datasets for analysis (separating phases, general window sub-setting)
## Familiarisation
LT.fam <- d[[4]] %>%
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
  subset_by_window(window_start_time = 1500,      # Start after stimulus moving in,
                   window_end_col = "TrialEnd") %>% # and end 4000ms after LabelOnset
  mutate(LabelOnset = LabelOnset - 1500) # Update LabelOnset after window sub-setting
## Contrast tests
LT.test.ctr <- d[[4]] %>%
  subset(Phase == "Test - Contrast") %>%
  make_eyetrackingr_data(participant_column = "Participant",
                         trial_column = "TrialId",
                         time_column = "TimeStamp",
                         trackloss_column = "TrackLoss",
                         aoi_columns = c("NewHead","OldHead","NewTail","OldTail","Centre"),
                         treat_non_aoi_looks_as_missing = T)
## Word learning tests
LT.test.wl <- d[[4]] %>%
  subset(Phase == "Test - Word Learning") %>%
  make_eyetrackingr_data(participant_column = "Participant",
                         trial_column = "TrialId",
                         time_column = "TimeStamp",
                         trackloss_column = "TrackLoss",
                         aoi_columns = c("Target","Distractor"),
                         treat_non_aoi_looks_as_missing = T)

# LOOKING TIME ANALYSIS: PROP AOI LOOKING BY PARTICIPANT BY TRIAL/BLOCK ============================
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
run_model <- T
if(run_model){
  ## Run lmer (Sampling Theory Based)
  LT.prop_tail.per_trial.lmer.model <- lmer(ArcSin ~ TrialNum*Condition +
                                              (1 | Participant) +
                                              (1 | Stimulus),
                                            data = LT.prop_tail)
  LT.prop_tail.per_trial.lmer.anova <- anova(LT.prop_tail.per_trial.lmer.model, type = 1)
  ## Run brms (Bayesian)
  LT.prop_tail.per_trial.brms.model <- brm(ArcSin ~ TrialNum*Condition +
                                             (1 | Participant) +
                                             (1 | Stimulus),
                                           data = LT.prop_tail,
                                           chains = 4, cores = 4)
  ## Save all results
  saveRDS(LT.prop_tail.per_trial.lmer.model, "../results/infants/Trial_lmerModel.rds")
  saveRDS(LT.prop_tail.per_trial.lmer.anova, "../results/infants/Trial_lmerAnova.rds")
  saveRDS(LT.prop_tail.per_trial.brms.model, "../results/infants/Trial_brmsModel.rds")
}else{
  ## Read all the results
  LT.prop_tail.per_trial.lmer.model <- readRDS("../results/infants/Trial_lmerModel.rds")
  LT.prop_tail.per_trial.lmer.anova <- readRDS("../results/infants/Trial_lmerAnova.rds")
  LT.prop_tail.per_trial.brms.model <- readRDS("../results/infants/Trial_brmsModel.rds")
}
# Testing Prop ~ Part*Condition
run_model <- T
if(run_model){
  ## Run lmer
  LT.prop_tail.per_part.lmer.model <- lmer(ArcSin ~ FamPart*Condition +
                                             (1 + FamPart | Participant) +
                                             (1 | Stimulus),
                                           data = LT.prop_tail)
  LT.prop_tail.per_part.lmer.anova <- anova(LT.prop_tail.per_part.lmer.model, type = 1)
  ## Run brms
  LT.prop_tail.per_part.brms.model <- brm(ArcSin ~ FamPart*Condition +
                                            (1 + FamPart | Participant) +
                                            (1 | Stimulus),
                                          data = LT.prop_tail,
                                          chains = 4, cores = 4)
  ## Save all the results
  saveRDS(LT.prop_tail.per_part.lmer.model, "../results/infants/Part_lmerModel.rds")
  saveRDS(LT.prop_tail.per_part.lmer.anova, "../results/infants/Part_lmerAnova.rds")
  saveRDS(LT.prop_tail.per_part.brms.model, "../results/infants/Part_brmsModel.rds")
}else{
  ## Read all the results
  LT.prop_tail.per_part.lmer.model <- readRDS("../results/infants/Part_lmerModel.rds")
  LT.prop_tail.per_part.lmer.anova <- readRDS("../results/infants/Part_lmerAnova.rds")
  LT.prop_tail.per_part.brms.model <- readRDS("../results/infants/Part_brmsModel.rds")
}
# Testing Prop ~ FstLst*Condition
run_model <- T
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
  LT.prop_tail.per_fstlst.brms.model <- brm(ArcSin ~ FstLst*Condition +
                                              (1 + FstLst | Participant) +
                                              (1 | Stimulus),
                                            data = LT.prop_tail.fstlst,
                                            chains = 4, cores = 4)
  ## Save all the results
  saveRDS(LT.prop_tail.per_fstlst.lmer.model, "../results/infants/FstLst_lmerModel.rds")
  saveRDS(LT.prop_tail.per_fstlst.lmer.anova, "../results/infants/FstLst_lmerAnova.rds")
  saveRDS(LT.prop_tail.per_fstlst.brms.model, "../results/infants/FstLst_brmsModel.rds")
}else{
  ## Read all the results
  LT.prop_tail.per_fstlst.lmer.model <- readRDS("../results/infants/FstLst_lmerModel.rds")
  LT.prop_tail.per_fstlst.lmer.anova <- readRDS("../results/infants/FstLst_lmerAnova.rds")
  LT.prop_tail.per_fstlst.brms.model <- readRDS("../results/infants/FstLst_brmsModel.rds")
}
#### NO SIGNIFICANT EFFECTS

# Plot jitter + mean&se + lines
## Plot per trial
LT.prop_tail.per_trial.plot <- ggplot(LT.prop_tail,
                                      aes(x = TrialId, y = Prop,
                                          colour = Condition,
                                          fill = Condition)) +
  theme(legend.pos = "top") + ylab("Looking to Tail (Prop)") +
  scale_colour_discrete(labels = c("Label", "No Label")) +
  scale_fill_discrete(labels = c("Label", "No Label")) +
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
ggsave("../results/infants/AOILookingPerTrial.pdf",
       LT.prop_tail.per_trial.plot,
       width = 7, height = 5.4)
## Plot per part
LT.prop_tail.per_part.plot <- ggplot(LT.prop_tail,
                                     aes(x = FamPart, y = Prop,
                                         colour = Condition,
                                         fill = Condition)) +
  theme(legend.pos = "top") + ylab("Looking to Tail (Prop)") +
  scale_colour_discrete(labels = c("Label", "No Label")) +
  scale_fill_discrete(labels = c("Label", "No Label")) +
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
ggsave("../results/infants/AOILookingPerParts.pdf",
       LT.prop_tail.per_part.plot,
       width = 7, height = 5.4)
## Plot per FstLst
LT.prop_tail.per_part.plot <- ggplot(LT.prop_tail.fstlst,
                                     aes(x = FstLst, y = Prop,
                                         colour = Condition,
                                         fill = Condition)) +
  theme(legend.pos = "top") + ylab("Looking to Tail (Prop)") +
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
ggsave("../results/infants/AOILookingPerFstLst.pdf",
       LT.prop_tail.per_part.plot,
       width = 7, height = 5.4)

# LOOKING TIME ANALYSIS: PROP AOI LOOKING PRE/POST LABEL ONSET =====================================
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
run_model <- T
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
                                set_prior("normal(0,.5)",
                                          coef = "TrialNum"),
                                set_prior("normal(0,.5)",
                                          coef = "PrePostPostLabelOnset"),
                                set_prior("normal(0,.5)",
                                          coef = "ConditionLabel"),
                                set_prior("normal(0,.5)",
                                          coef = "TrialNum:PrePostPostLabelOnset"),
                                set_prior("normal(0,.5)",
                                          coef = "TrialNum:ConditionLabel"),
                                set_prior("normal(0,.5)",
                                          coef = "PrePostPostLabelOnset:ConditionLabel"),
                                set_prior("normal(0,.5)",
                                          coef = "TrialNum:PrePostPostLabelOnset:ConditionLabel"))
  LT.pre_post.per_trial.brms.model <- brm(ArcSin ~ TrialNum*PrePost*Condition +
                                            (1 | Participant) +
                                            (1 | Stimulus),
                                          data = LT.pre_post,
                                          prior = prior.pre_post.per_trial,
                                          chains = 4, cores = 4)
  ## Save all results
  saveRDS(LT.pre_post.per_trial.lmer.model, "../results/infants/PrePost_Trial_lmerModel.rds")
  saveRDS(LT.pre_post.per_trial.lmer.anova, "../results/infants/PrePost_Trial_lmerAnova.rds")
  saveRDS(LT.pre_post.per_trial.brms.model, "../results/infants/PrePost_Trial_brmsModel.rds")
}else{
  ## Read all the results
  LT.pre_post.per_trial.lmer.model <- readRDS("../results/infants/PrePost_Trial_lmerModel.rds")
  LT.pre_post.per_trial.lmer.anova <- readRDS("../results/infants/PrePost_Trial_lmerAnova.rds")
  LT.pre_post.per_trial.brms.model <- readRDS("../results/infants/PrePost_Trial_brmsModel.rds")
}
# Testing Prop ~ Part*PrePost*Condition
run_model <- T
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
                               set_prior("normal(0,.5)",
                                         coef = "FamPart"),
                               set_prior("normal(0,.5)",
                                         coef = "PrePostPostLabelOnset"),
                               set_prior("normal(0,.5)",
                                         coef = "ConditionLabel"),
                               set_prior("normal(0,.5)",
                                         coef = "FamPart:PrePostPostLabelOnset"),
                               set_prior("normal(0,.5)",
                                         coef = "FamPart:ConditionLabel"),
                               set_prior("normal(0,.5)",
                                         coef = "PrePostPostLabelOnset:ConditionLabel"),
                               set_prior("normal(0,.5)",
                                         coef = "FamPart:PrePostPostLabelOnset:ConditionLabel"))
  LT.pre_post.per_part.brms.model <- brm(ArcSin ~ FamPart*PrePost*Condition +
                                           (1 + FamPart | Participant) +
                                           (1 | Stimulus),
                                         data = LT.pre_post,
                                         prior = prior.pre_post.per_part,
                                         chains = 4, cores = 4)
  ## Save all the results
  saveRDS(LT.pre_post.per_part.lmer.model, "../results/infants/PrePost_Part_lmerModel.rds")
  saveRDS(LT.pre_post.per_part.lmer.anova, "../results/infants/PrePost_Part_lmerAnova.rds")
  saveRDS(LT.pre_post.per_part.brms.model, "../results/infants/PrePost_Part_brmsModel.rds")
}else{
  ## Read all the results
  LT.pre_post.per_part.lmer.model <- readRDS("../results/infants/PrePost_Part_lmerModel.rds")
  LT.pre_post.per_part.lmer.anova <- readRDS("../results/infants/PrePost_Part_lmerAnova.rds")
  LT.pre_post.per_part.brms.model <- readRDS("../results/infants/PrePost_Part_brmsModel.rds")
}
# Testing Prop ~ FstLst*Condition
run_model <- T
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
                                 set_prior("normal(0,.5)",
                                           coef = "FstLstLastTrials"),
                                 set_prior("normal(0,.5)",
                                           coef = "PrePostPostLabelOnset"),
                                 set_prior("normal(0,.5)",
                                           coef = "ConditionLabel"),
                                 set_prior("normal(0,.5)",
                                           coef = "FstLstLastTrials:PrePostPostLabelOnset"),
                                 set_prior("normal(0,.5)",
                                           coef = "FstLstLastTrials:ConditionLabel"),
                                 set_prior("normal(0,.5)",
                                           coef = "PrePostPostLabelOnset:ConditionLabel"),
                                 set_prior("normal(0,.5)",
                                           coef = "FstLstLastTrials:PrePostPostLabelOnset:ConditionLabel"))
  LT.pre_post.per_fstlst.brms.model <- brm(ArcSin ~ FstLst*PrePost*Condition +
                                             (1 + FstLst | Participant) +
                                             (1 | Stimulus),
                                           data = LT.pre_post.fstlst,
                                           prior = prior.pre_post.per_fstlst,
                                           chains = 4, cores = 4)
  ## Save all the results
  saveRDS(LT.pre_post.per_fstlst.lmer.model, "../results/infants/PrePost_FstLst_lmerModel.rds")
  saveRDS(LT.pre_post.per_fstlst.lmer.anova, "../results/infants/PrePost_FstLst_lmerAnova.rds")
  saveRDS(LT.pre_post.per_fstlst.brms.model, "../results/infants/PrePost_FstLst_brmsModel.rds")
}else{
  ## Read all the results
  LT.pre_post.per_fstlst.lmer.model <- readRDS("../results/infants/PrePost_FstLst_lmerModel.rds")
  LT.pre_post.per_fstlst.lmer.anova <- readRDS("../results/infants/PrePost_FstLst_lmerAnova.rds")
  LT.pre_post.per_fstlst.brms.model <- readRDS("../results/infants/PrePost_FstLst_brmsModel.rds")
}

# Plot jitter + mean&se + lines
## Plot per part
LT.pre_post.per_part.plot <- ggplot(LT.pre_post,
                                    aes(x = PrePost, y = Prop,
                                        colour = Condition,
                                        fill = Condition)) +
  theme(legend.pos = "top") + ylab("Looking to Tail (Prop)") + facet_grid(.~FamPart) +
  scale_colour_discrete(labels = c("Label", "No Label")) +
  scale_fill_discrete(labels = c("Label", "No Label")) +
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
ggsave("../results/infants/AOILookingPrePostPerParts.pdf",
       LT.pre_post.per_part.plot,
       width = 7, height = 3)
## Plot per FstLst
LT.pre_post.per_fstlst.plot <- ggplot(LT.pre_post.fstlst,
                                      aes(x = PrePost, y = Prop,
                                          colour = Condition,
                                          fill = Condition)) +
  theme(legend.pos = "top") + ylab("Looking to Tail (Prop)") + facet_grid(.~FstLst) +
  scale_colour_discrete(labels = c("Label", "No Label")) +
  scale_fill_discrete(labels = c("Label", "No Label")) +
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
ggsave("../results/infants/AOILookingPrePostPerFstLst.pdf",
       LT.pre_post.per_fstlst.plot,
       width = 7, height = 3)

# LOOKING TIME ANALYSIS: TIME COURSE ===============================================================
# DATA PREPARATION
LT.time_course_tail <- LT.fam %>%
  drop_na(FstLst) %>%
  subset_by_window(window_start_col = "LabelOnset", remove = F) %>%
  make_time_sequence_data(time_bin_size = 50,
                          aois = "Tail",
                          predictor_columns=c("Condition",
                                              "Stimulus",
                                              "FstLst"))
# GROWTH CURVE ANALYSIS
## TODO
# BOOTSTRAPPED CLUSTER-BASED PERMUTATION ANALYSIS
run_model <- T
if(run_model){
  t <- proc.time()
  ## Determine clusters
  LT.time_cluster_tail <- LT.time_course_tail %>%
    split(.$FstLst) %>%
    lapply(make_time_cluster_data,
           predictor_column = "Condition",
           treatment_level = "No Label",
           aoi = "Tail",
           test = "lmer",
           threshold = 1,
           formula = ArcSin ~ Condition +
             (1 | Participant) +
             (1 | Stimulus))
  ## Run analysis
  LT.time_cluster_tail.analysis <- LT.time_cluster_tail %>%
    lapply(analyze_time_clusters,
           formula = ArcSin ~ Condition +
             (1 | Participant) +
             (1 | Stimulus),
           within_subj = T,
           parallel = T)
  # #### TESTING TO REPORT ON eyetrackingR GITHUB
  #   ## Determine clusters
  #   LT.time_cluster_tail <- make_time_cluster_data(LT.time_course_tail,
  #                                                  predictor_column = "ConditionNo Label",
  #                                                  aoi = "Tail",
  #                                                  test = "lmer",
  #                                                  threshold = 1.5,
  #                                                  formula = ArcSin ~ Condition +
  #                                                    (1 | Participant) +
  #                                                    (1 | Stimulus))
  #   ## Run analysis
  #   LT.time_cluster_tail.analysis <- analyze_time_clusters(LT.time_cluster_tail,
  #                                                          formula = ArcSin ~ Condition +
  #                                                            (1 | Participant) +
  #                                                            (1 | Stimulus),
  #                                                          within_subj = T,
  #                                                          parallel = T,
  #                                                          samples = 200)
  bcbp.time <- proc.time() - t
  ## Save clusters and analysis
  saveRDS(LT.time_cluster_tail,
          "../results/infants/BCBP_clusters.rds")
  saveRDS(LT.time_cluster_tail.analysis,
          "../results/infants/BCBP_analysis.rds")
}else{
  ## Read the results
  LT.time_cluster_tail <- readRDS("../results/infants/BCBP_clusters.rds")
  LT.time_cluster_tail.analysis <- readRDS("../results/infants/BCBP_analysis.rds")
}
# PLOT
intercept <- tibble(Part = 0:2,
                    x_int = rep(1500,3)) # Label onset ish (second half trials includes "the")
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
#geom_vline(data = intercept, aes(xintercept = x_int), linetype = "62", alpha = .5)
ggsave("../results/infants/LookingTimeCourseFirstLast.pdf",
       plot = LT.fam.time_course.plot.blocks,
       width = 3.5, height = 5)
