library(eyetrackingR)
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
                         treat_non_aoi_looks_as_missing = T) %>%
  subset_by_window(window_start_time = 1500, window_end_time = 7000)

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

# LOOKING TIME ANALYSIS: TIME COURSE ===============================================================
# Preparing data for analysis and plot
LT.time_course_aois <- LT.clean %>%
  make_time_sequence_data(time_bin_size = 50,
                          aois = c("Head","Tail"),
                          predictor_columns=c("Condition")) %>%
  mutate_at("TrialId", as.numeric) %>%
  subset(TrialId < 25) %>%
  mutate(Part = ((TrialId-1) %/% 8))
# Plotting first/middle/last 8 trials
intercept <- tibble(Part = 0:2,
                    x_int = rep(2250,3)) # Label onset ish (second half trials includes "the")
LT.clean.time_course.plot.blocks <- ggplot(LT.time_course_aois,
                                           aes(x = Time, y=Prop,
                                               colour=Condition,
                                               fill=Condition)) +
  xlab('Time in Trial') + ylab("Looking to AOI (Prop)") +
  facet_grid(AOI~Part) + theme(legend.position = "top") + ylim(0,1) +
  stat_summary(fun.y='mean', geom='line', linetype = '61') +
  stat_summary(fun.data=mean_se, geom='ribbon', alpha= .25, colour=NA) +
  geom_hline(yintercept = .5) +
  geom_vline(data = intercept, aes(xintercept = x_int), linetype = "62", alpha = .5)
ggsave("../results/infants/LookingTimeCoursePerPart.pdf",
       plot = LT.clean.time_course.plot.blocks,
       width = 7, height = 4)
