library(lme4)
library(nortest)
library(tidyverse)
library(eyetrackingR)

source("Routines.R")

# ==================================================================================================
# GATHER DATA
# ==================================================================================================
# TODO - AOIs plot needs updating
# AOIs.plot <- ggplot(AOIs.adults, aes(xmin = L, xmax = R, ymin = T, ymax = B)) +
#   xlim(c(0,640)) + scale_y_reverse(limits = c(480,0)) +
#   geom_rect(aes(fill = name))
# ggsave("../results/adults_2f/data_cleaning_graphs/AOIs.png",
#        plot = AOIs.plot , width = 3.2, height = 2.95)

# Data import and storage into variables
d <- LT_data.gather("adults_2f")
behaviour <- d[[2]]
LT.clean <- d[[4]] %>%
  subset(Block <= 8) %>%
  make_eyetrackingr_data(participant_column = "Participant",
                         trial_column = "TrialId",
                         time_column = "TimeStamp",
                         trackloss_column = "TrackLoss",
                         aoi_columns = c("Head","Tail"),
                         treat_non_aoi_looks_as_missing = T) %>%
  subset_by_window(window_start_time = -1000, rezero = F)

# ==================================================================================================
# LOOKING TIME ANALYSIS: PROP TAIL LOOKING BY PARTICIPANT BY BLOCK
# ==================================================================================================
# Plot proportion looking to tail for each participant, across time (trials or blocks?)
LT.prop_tail_per_block <- make_time_window_data(LT.clean,
                                                aois="Tail",
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
                                          colour = Condition)) +
  facet_wrap(~BlockTransformation, scales = "free_x") +
  theme(aspect.ratio = 1.618/1, legend.position = "top") +
  stat_smooth(linetype = "61", level = 0.87) +
  geom_hline(yintercept = asin(sqrt(.5)))
ggsave("../results/adults_2f/TailLookingEvolution.png",
       plot = LT.prop_tail_per_block.plot,
       width = 7, height = 5)
# lmer of the plot above
LT.prop_tail_per_block.lmer.0 <- lmer(ArcSin ~ Block*Condition*ACC +
                                        (1 + Block + ACC | Participant),
                                      data = subset(LT.prop_tail_per_block,
                                                    BlockTransformation == "NormBlock"))
LT.prop_tail_per_block.lmer.1 <- update(LT.prop_tail_per_block.lmer.0,
                                        . ~ . - Block:Condition:ACC)
LT.prop_tail_per_block.lmer.2 <- update(LT.prop_tail_per_block.lmer.1,
                                        . ~ . - Condition:ACC)
LT.prop_tail_per_block.lmer.3 <- update(LT.prop_tail_per_block.lmer.2,
                                        . ~ . - Block:ACC)
LT.prop_tail_per_block.lmer.4 <- update(LT.prop_tail_per_block.lmer.3,
                                        . ~ . - Block:Condition)
LT.prop_tail_per_block.lmer.5 <- update(LT.prop_tail_per_block.lmer.4,
                                        . ~ . - Condition)
LT.prop_tail_per_block.lmer.6 <- update(LT.prop_tail_per_block.lmer.5,
                                        . ~ . - Block)
LT.prop_tail_per_block.lmer.comp <- anova(LT.prop_tail_per_block.lmer.6,
                                          LT.prop_tail_per_block.lmer.5,
                                          LT.prop_tail_per_block.lmer.4,
                                          LT.prop_tail_per_block.lmer.3,
                                          LT.prop_tail_per_block.lmer.2,
                                          LT.prop_tail_per_block.lmer.1,
                                          LT.prop_tail_per_block.lmer.0)
# ==================================================================================================
# LOOKING TIME ANALYSIS: TIME COURSE
# ==================================================================================================
LT.time_course_tail <- make_time_sequence_data(LT.clean, time_bin_size = 100,
                                                aois="Tail",
                                                predictor_columns=c("Condition",
                                                                    "Block"))
LT.clean.time_course.plot.blocks <- ggplot(LT.time_course_tail,
                                           aes(x = Time, y=Prop,
                                               colour=Condition,
                                               fill=Condition)) +
  xlab('Time in Trial') + ylab("Looking to Tail (Prop)") +
  facet_wrap(~Block) + theme(legend.position = "top") + ylim(0,1) +
  stat_summary(fun.y='mean', geom='line', linetype = 'F1') +
  stat_summary(fun.data=mean_se, geom='ribbon', alpha= .25, colour=NA) +
  geom_hline(yintercept = .5)
ggsave("../results/adults_2f/LookingTimeCoursePerBlock.pdf",
       plot = LT.clean.time_course.plot.blocks)
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
ggsave("../results/adults_2f/ParticipantsPerBlock.png",
       plot = behaviour.parts_per_block.plot)
# Get number of blocks to learning per participant, plot a violin
behaviour.blocks_per_part <- behaviour %>%
  group_by(Participant, Condition) %>%
  summarise(N_Blocks = max(Block))
behaviour.blocks_per_part.plot <- ggplot(behaviour.blocks_per_part,
                                        aes(x = Condition, y = N_Blocks, fill = Condition)) +
  geom_violin() +
  geom_boxplot(alpha = 0, outlier.alpha = 1, width = .15)
ggsave("../results/adults_2f/BlocksPerParticipant.png",
       plot = behaviour.blocks_per_part.plot)
