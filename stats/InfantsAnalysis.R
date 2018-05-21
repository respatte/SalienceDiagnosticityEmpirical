# LIBRARY IMPORTS ==================================================================================
library(eyetrackingR)
library(lme4)
library(tidyverse)

source("Routines.R")

# GATHER DATA ======================================================================================
d <- LT_data.gather("infants")
LT.fam.clean <- d[[4]] %>%
  subset(Phase = "Familiarisation") %>%
  make_eyetrackingr_data(participant_column = "Participant",
                         trial_column = "TrialId",
                         time_column = "TimeStamp",
                         trackloss_column = "TrackLoss",
                         aoi_columns = c("Head","Tail"),
                         treat_non_aoi_looks_as_missing = T) %>%
  subset_by_window(window_start_col = "LabelOnset", window_end_col = "TrialEnd")
# Check for available PresentationSequence
pres_seq <- LT.fam.clean %>%
  group_by(PresentationSequence, Participant) %>%
  summarise(T = sum(TrackLoss)/n())
gender <- LT.fam.clean %>%
  group_by(Gender, Condition) %>%
  summarise(N = n_distinct(Participant))
age <- LT.fam.clean %>%
  group_by(Participant, Condition) %>%
  summarise(Age = first(Age))

# LOOKING TIME ANALYSIS: PROP AOI LOOKING BY PARTICIPANT BY TRIAL/BLOCK ============================
# Prepare dataset, include only familiarisation
LT.prop_tail <- make_time_window_data(LT.fam.clean,
                                                aois=c("Tail"),
                                                predictor_columns=c("Condition",
                                                                    "TrialId",
                                                                    "Stimulus",
                                                                    "CategoryName")) %>%
  mutate_at("TrialId", as.numeric) %>%
  subset(TrialId < 25) %>%
  mutate(Part = (TrialId-1) %/% 8,
         Trial = (TrialId-1) %/% 2)
## LMER for Prop ~ Condition*Part*AOI
# No main effect of Part since looking at proportions,
# so overall they should all equate to one when collapsing
LT.prop_tail.per_trial.lmer.0 <- lmer(ArcSin ~ Trial*Condition +
                                        (1 | Participant) +
                                        (1 | Stimulus),
                                        data = LT.prop_tail)
LT.prop_tail.per_trial.lmer.1 <- update(LT.prop_tail.per_trial.lmer.0,
                                        . ~ . - Trial:Condition)
LT.prop_tail.per_trial.lmer.2 <- update(LT.prop_tail.per_trial.lmer.1,
                                        . ~ . - Condition)
LT.prop_tail.per_trial.lmer.3 <- update(LT.prop_tail.per_trial.lmer.2,
                                        . ~ . - Trial)
LT.prop_tail.per_trial.lmer.comp <- anova(LT.prop_tail.per_trial.lmer.3,
                                          LT.prop_tail.per_trial.lmer.2,
                                          LT.prop_tail.per_trial.lmer.1,
                                          LT.prop_tail.per_trial.lmer.0)
#### NO SIGNIFICANT EFFECT
# Plot jitter + lmer mean&se + lines
LT.prop_tail.per_part.plot <- ggplot(LT.prop_tail,
                                       aes(x = Part, y = Prop,
                                           colour = Condition,
                                           fill = Condition)) +
  theme_apa(legend.pos = "top") + ylab("Looking to AOI (Prop)") +
  scale_colour_discrete(labels = c("Label", "No Label")) +
  scale_fill_discrete(labels = c("Label", "No Label")) +
  geom_point(position = position_jitterdodge(dodge.width = .8,
                                             jitter.width = .2),
             alpha = .25) +
  geom_errorbar(stat = "summary",
                width = .2, colour = "black",
                position = position_dodge(.1)) +
  geom_line(aes(x = Part, y = Prop, group = Condition),
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
LT.time_course_aois <- LT.fam.clean %>%
  subset(Phase == "Familiarisation") %>%
  make_time_sequence_data(time_bin_size = 50,
                          aois = "Tail",
                          predictor_columns=c("Condition",
                                              "Stimulus"))
# GROWTH CURVE ANALYSIS
## TODO
# BOOTSTRAPPED CLUSTER-BASED PERMUTATION ANALYSIS
LT.time_cluster_aois <- make_time_cluster_data(LT.time_course_aois, 
                                               predictor_column = "Condition",
                                               treatment_level = "No Label",
                                               aoi = "Tail",
                                               test = "lmer",
                                               threshold = 1.5,
                                               formula = ArcSin ~ TrialId*Condition +
                                                 (1 | Participant) +
                                                 (1 | Stimulus))
LT.time_cluster_aois.analysis <- analyze_time_clusters(LT.time_cluster_aois,
                                                       within_subj = F)
# PLOT
intercept <- tibble(Part = 0:2,
                    x_int = rep(1500,3)) # Label onset ish (second half trials includes "the")
LT.fam.clean.time_course.plot.blocks <- ggplot(LT.time_course_aois,
                                           aes(x = Time, y=Prop,
                                               colour=Condition,
                                               fill=Condition)) +
  xlab('Time in Trial') + ylab("Looking to AOI (Prop)") +
  facet_grid(.~Part) +
  theme(legend.position = "top") + ylim(0,1) +
  stat_summary(fun.y='mean', geom='line', linetype = '61') +
  stat_summary(fun.data=mean_se, geom='ribbon', alpha= .25, colour=NA) +
  geom_hline(yintercept = .5) +
  geom_vline(data = intercept, aes(xintercept = x_int), linetype = "62", alpha = .5)
ggsave("../results/infants/LookingTimeCoursePerPart.pdf",
       plot = LT.fam.clean.time_course.plot.blocks,
       width = 7, height = 2.5)
