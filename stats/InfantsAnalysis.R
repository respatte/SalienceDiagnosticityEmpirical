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
  mutate(Part = (TrialId-1) %/% 8)
# Plot jitter + lmer mean&se + lines
LT.prop_aois.first_last.plot <- ggplot(LT.prop_aois.first_last,
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
       LT.prop_aois.first_last.plot,
       width = 7, height = 5.4)th = .15,
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
