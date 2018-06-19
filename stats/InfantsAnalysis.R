# LIBRARY IMPORTS ==================================================================================
library(eyetrackingR)
library(lme4)
library(lmerTest)
library(brms)
library(jtools)
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
  group_by(Participant) %>%
  mutate(First = min(TrialNum, na.rm = T),
         Last = max(TrialNum, na.rm = T),
         FstLst = case_when(TrialNum <= First + 2 ~ "First Trials",
                            TrialNum >= Last - 2 ~ "Last Trials")) %>%
         # Useful to compare beginning-end of experiment per infant
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
# Prepare dataset, include only familiarisation
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
  saveRDS(LT.prop_tail.per_trial.lmer.model, "../results/infants/TrialByCondition_lmerModel.rds")
  saveRDS(LT.prop_tail.per_trial.lmer.anova, "../results/infants/TrialByCondition_lmerAnova.rds")
  saveRDS(LT.prop_tail.per_trial.brms.model, "../results/infants/TrialByCondition_brmsModel.rds")
}else{
  ## Read all the results
  LT.prop_tail.per_trial.lmer.model <- readRDS("../results/infants/TrialByCondition_lmerModel.rds")
  LT.prop_tail.per_trial.lmer.anova <- readRDS("../results/infants/TrialByCondition_lmerAnova.rds")
  LT.prop_tail.per_trial.brms.model <- readRDS("../results/infants/TrialByCondition_brmsModel.rds")
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
  saveRDS(LT.prop_tail.per_part.lmer.model, "../results/infants/PartByCondition_lmerModel.rds")
  saveRDS(LT.prop_tail.per_part.lmer.anova, "../results/infants/PartByCondition_lmerAnova.rds")
  saveRDS(LT.prop_tail.per_part.brms.model, "../results/infants/PartByCondition_brmsModel.rds")
}else{
  ## Read all the results
  LT.prop_tail.per_part.lmer.model <- readRDS("../results/infants/PartByCondition_lmerModel.rds")
  LT.prop_tail.per_part.lmer.anova <- readRDS("../results/infants/PartByCondition_lmerAnova.rds")
  LT.prop_tail.per_part.brms.model <- readRDS("../results/infants/PartByCondition_brmsModel.rds")
}
# Testing Prop ~ FstLst*Condition
run_model <- T
if(run_model){
  ## Run lmer
  LT.prop_tail.per_part.lmer.model <- lmer(ArcSin ~ FstLst*Condition +
                                             (1 + FstLst | Participant) +
                                             (1 | Stimuli),
                                           data = LT.prop_tail)
  LT.prop_tail.per_part.lmer.anova <- anova(LT.prop_tail.per_part.lmer.model, type = 1)
  ## Run brms
  LT.prop_tail.per_part.brms.model <- brm(ArcSin ~ FstLst*Condition +
                                            (1 + FstLst | Participant) +
                                            (1 | Stimulus),
                                          data = LT.prop_tail,
                                          chains = 4, cores = 4)
  ## Save all the results
  saveRDS(LT.prop_tail.per_part.lmer.model, "../results/infants/PartByCondition_lmerModel.rds")
  saveRDS(LT.prop_tail.per_part.lmer.anova, "../results/infants/PartByCondition_lmerAnova.rds")
  saveRDS(LT.prop_tail.per_part.brms.model, "../results/infants/PartByCondition_brmsModel.rds")
}else{
  ## Read all the results
  LT.prop_tail.per_part.lmer.model <- readRDS("../results/infants/PartByCondition_lmerModel.rds")
  LT.prop_tail.per_part.lmer.anova <- readRDS("../results/infants/PartByCondition_lmerAnova.rds")
  LT.prop_tail.per_part.brms.model <- readRDS("../results/infants/PartByCondition_brmsModel.rds")
}
#### NO SIGNIFICANT EFFECTS
# Plot jitter + lmer mean&se + lines
## Plot per trial
LT.prop_tail.per_trial.plot <- ggplot(LT.prop_tail,
                                     aes(x = TrialId, y = Prop,
                                         colour = Condition,
                                         fill = Condition)) +
  theme_apa(legend.pos = "top") + ylab("Looking to Tail (Prop)") +
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
            colour = "black") +
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
  theme_apa(legend.pos = "top") + ylab("Looking to Tail (Prop)") +
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
            colour = "black") +
  geom_point(stat = "summary", fun.y = "mean",
             shape = 18, size = 3,
             position = position_dodge(.1))
ggsave("../results/infants/AOILookingPerParts.pdf",
       LT.prop_tail.per_part.plot,
       width = 7, height = 5.4)

# LOOKING TIME ANALYSIS: TIME COURSE ===============================================================
# DATA PREPARATION
LT.time_course_tail <- LT.fam %>%
  drop_na(FstLst) %>%
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
