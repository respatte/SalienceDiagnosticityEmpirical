library(eyetrackingR)
library(Matrix)
library(lme4)
library(tidyverse)

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
                         treat_non_aoi_looks_as_missing = F)

behaviour.part_per_block <- behaviour %>%
  subset(Block > 0) %>%
  group_by(Block, Condition) %>%
  summarise(N_Participants = n_distinct(Participant))
behaviour.part_per_block.plot <- ggplot(behaviour.part_per_block,
                                        aes(x = Block, y = N_Participants, fill = Condition)) +
  geom_col(position = "dodge")