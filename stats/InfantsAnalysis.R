library(eyetrackingR)
library(Matrix)
library(lme4)
library(ggplot2)
library(dplyr)

source("Routines.R")

# GATHER DATA
# Define AOIs
AOIs.infants.head <- data.frame(AOI_type=c("Reg","Flip"),
                                L=c(1031,1920-1031-450),
                                R=c(1031+450,1920-1031),
                                T=c(197,197),
                                B=c(197+450,197+450))
AOIs.infants.tail <- data.frame(AOI_type=c("Reg","Flip"),
                                L=c(390,1920-390-450),
                                R=c(390+450,1920-390),
                                T=c(299,299),
                                B=c(299+450,299+450))
AOIs.plot <- ggplot(NULL, aes(xmin = L, xmax = R, ymin = T, ymax = B)) +
  xlim(c(0,1920)) + scale_y_reverse(limits = c(1080,0)) +
  geom_rect(data = AOIs.infants.head, aes(fill = AOI_type)) +
  geom_rect(data = AOIs.infants.tail, aes(fill = AOI_type))
ggsave("../results/infants/data_cleaning_graphs/AOIs.png",
       plot = AOIs.plot , width = 6.05, height = 2.95)

# Import raw data
raw_data.infants <- LT_data.infants.import()
# Turn raw into clean eyetrackingR data
LT.infants <- raw_data.infants %>%
  #add_aoi(., AOIs.infants.head, "CursorX", "CursorY", aoi_name = "Head") %>%
  #add_aoi(., AOIs.infants.tail, "CursorX", "CursorY", aoi_name = "Tail") %>%
  LT_data.to_eyetrackingR(AOIs.adults, participants = "infants") %>%
  make_eyetrackingr_data(participant_column = "ParticipantName",
                         trial_column = "TrialN",
                         time_column = "TimeStamp",
                         trackloss_column = "TrackLoss",
                         aoi_columns = c('Head','Tail'),
                         treat_non_aoi_looks_as_missing = F) %>%
  LT_data.trackloss_clean(res.repo = "../results/infants/data_cleaning_graphs/")
