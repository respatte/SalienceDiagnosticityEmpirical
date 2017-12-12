library(lme4)
library(tidyverse)
library(eyetrackingR)

source("Routines.R")

# GATHER DATA
# Define AOIs
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
  make_eyetrackingr_data(participant_column = "Participant",
                         trial_column = "TrialId",
                         time_column = "TimeStamp",
                         trackloss_column = "TrackLoss",
                         aoi_columns = c("Head","Tail"),
                         treat_non_aoi_looks_as_missing = T)

# ==================================================================================================
# LOOKING TIME ANALYSIS: PROP TAIL LOOKING BY PARTICIPANT BY BLOCK
# ==================================================================================================
# Plot proportion looking to tail for each participant, across time (trials or blocks?)
LT.prop_tail_per_block <- make_time_window_data(LT.clean,
                                                aois="Tail",
                                                predictor_columns=c("Condition",
                                                                    "Block"),
                                                summarize_by = "Participant") %>%
  subset(Block > 0)
LT.prop_tail_per_block.plot <- ggplot(LT.prop_tail_per_block,
                                      aes(x = Block, y = Prop,
                                          colour = Participant,
                                          group = Participant)) +
  geom_line() +
  geom_hline(yintercept = .5)

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
