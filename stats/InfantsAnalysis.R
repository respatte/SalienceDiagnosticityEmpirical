library(eyetrackingR)
library(Matrix)
library(lme4)
library(tidyverse)

source("Routines.R")

# GATHER DATA ======================================================================================
d <- LT_data.gather("infants")
LT.clean <- d[[4]] %>%
  make_eyetrackingr_data(participant_column = "Participant",
                         trial_column = "TrialId",
                         time_column = "TimeStamp",
                         trackloss_column = "TrackLoss",
                         aoi_columns = c("Head","Tail"),
                         treat_non_aoi_looks_as_missing = T)

# LOOKING TIME ANALYSIS: PROP AOI LOOKING BY PARTICIPANT BY TRIAL/BLOCK ============================
# Prepare dataset, include only familiarisation
LT.prop_tail_per_trial <- make_time_window_data(LT.clean,
                                                aois=c("Tail"),
                                                predictor_columns=c("Condition",
                                                                    "TrialId")) %>%
  mutate_at("TrialId", as.numeric) %>%
  subset(TrialId < 25) %>%
  mutate(Block = ((TrialId-1) %/% 8) + 1)
# Plot points and smoother accross trials
LT.prop_tail_per_trial.plot <- ggplot(LT.prop_tail_per_trial,
                                      aes(x = TrialId, y = Prop,
                                          colour = Condition)) +
  geom_jitter() +
  geom_smooth()
# Plot violin + boxplot for each block
LT.prop_tail_per_block.plot <- ggplot(LT.prop_tail_per_trial,
                                      aes(x = as.factor(Block), y = Prop,
                                          fill = Condition)) +
  geom_violin() +
  geom_boxplot(alpha = 0, outlier.alpha = 1, width = .15,
               position = position_dodge(.9))
