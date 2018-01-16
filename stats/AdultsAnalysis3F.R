library(lme4)
library(nortest)
library(tidyverse)
library(eyetrackingR)

source("Routines.R")

# ==================================================================================================
# GATHER DATA
# ==================================================================================================
d <- LT_data.gather("adults_3f")
behaviour <- d[[2]]
LT.clean <- d[[4]] %>%
  subset(Block <= 17) %>%
  make_eyetrackingr_data(participant_column = "Participant",
                         trial_column = "TrialId",
                         time_column = "TimeStamp",
                         trackloss_column = "TrackLoss",
                         aoi_columns = c("Head","Tail"),
                         treat_non_aoi_looks_as_missing = T) %>%
  subset_by_window(window_start_time = -1000, rezero = F)

# ==================================================================================================
# LOOKING TIME ANALYSIS: TIME COURSE
# ==================================================================================================
# Plotting eye-tracking data for all AOIs, averaged across all trials
LT.time_course_aois <- make_time_sequence_data(LT.clean, time_bin_size = 100,
                                               aois = c("Head","Tail","Feet"),
                                               predictor_columns=c("Condition",
                                                                   "Block"))
LT.clean.time_course.plot.blocks <- ggplot(LT.time_course_aois,
                                           aes(x = Time, y=Prop,
                                               colour=Condition,
                                               fill=Condition)) +
  xlab('Time in Trial') + ylab("Looking to AOI (Prop)") +
  facet_grid(Block~AOI) + theme(legend.position = "top") + ylim(0,1) +
  stat_summary(fun.y='mean', geom='line', linetype = 'F1') +
  stat_summary(fun.data=mean_se, geom='ribbon', alpha= .25, colour=NA) +
  geom_hline(yintercept = 1/3)
ggsave("../results/adults_3f/LookingTimeCoursePerBlock.pdf",
       plot = LT.clean.time_course.plot.blocks,
       width = 7, height = 30)
# Making data time-window-analysis ready
LT.clean.total_per_AOI <- make_time_window_data(LT.clean,
                                                 aois=c("Head","Tail","Feet"),
                                                 predictor_columns=c("Condition",
                                                                     "Stimulus",
                                                                     "TrialId",
                                                                     "CategoryName"))
# ==================================================================================================
# LOOKING TIME ANALYSIS: PROP AOI LOOKING BY PARTICIPANT BY BLOCK
# ==================================================================================================
LT.prop_tail_per_block <- make_time_window_data(LT.clean,
                                                aois=c("Tail","Feet","Head"),
                                                predictor_columns=c("Condition",
                                                                    "Block",
                                                                    "ACC")) %>%
  subset(Block > 0) %>%
  group_by(Participant) %>%
  mutate(N_Blocks = max(Block),
         OppBlock = Block - N_Blocks,
         NormBlock = Block/N_Blocks) %>%
  gather("BlockTransformation","Block", Block, OppBlock, NormBlock)
LT.prop_tail_per_block.plot <- ggplot(LT.prop_tail_per_block,
                                      aes(x = Block, y = ArcSin,
                                          colour = Condition,
                                          fill = Condition)) +
  facet_grid(AOI~BlockTransformation, scales = "free_x") +
  theme(aspect.ratio = 1.618/1, legend.position = "top") +
  stat_smooth(linetype = "61", level = 0.87) +
  geom_hline(yintercept = asin(sqrt(1/3)))
ggsave("../results/adults_3f/AOILookingEvolution.png",
       plot = LT.prop_tail_per_block.plot,
       width = 7, height = 13)

# ==================================================================================================
# BEHAVIOURAL ANALYSIS: PARTICIPANTS AND BLOCKS
# ==================================================================================================
# Get how many participants for each block, and make a bar plot
behaviour.parts_per_block <- behaviour %>%
  subset(Block > 0) %>%
  group_by(Block, Condition) %>%
  summarise(N_Participants = n_distinct(Participant))
behaviour.parts_per_block.plot <- ggplot(behaviour.parts_per_block,
                                         aes(x = Block, y = N_Participants, fill = Condition)) +
  geom_col(position = "dodge")
ggsave("../results/adults_3f/ParticipantsPerBlock.png",
       plot = behaviour.parts_per_block.plot)
# Get number of blocks to learning per participant, plot a violin
behaviour.blocks_per_part <- behaviour %>%
  group_by(Participant, Condition) %>%
  summarise(N_Blocks = max(Block))
behaviour.blocks_per_part.plot <- ggplot(behaviour.blocks_per_part,
                                         aes(x = Condition, y = N_Blocks, fill = Condition)) +
  geom_violin() +
  geom_boxplot(alpha = 0, outlier.alpha = 1, width = .15)
ggsave("../results/adults_3f/BlocksPerParticipant.png",
       plot = behaviour.blocks_per_part.plot)
