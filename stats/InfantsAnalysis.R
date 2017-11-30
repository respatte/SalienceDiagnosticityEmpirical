library(eyetrackingR)
library(Matrix)
library(lme4)
library(ggplot2)
library(dplyr)

source("Routines.R")

# GATHER DATA
# Define AOIs
AOIs.infants.head <- data.frame(Facing=c("Reg","Flip"),L=c(20,400),R=c(220,620),T=c(110,55),B=c(330,255))
AOIs.infants.tail <- data.frame(Facing=c("Reg","Flip"),L=c(20,400),R=c(220,620),T=c(110,55),B=c(330,255))
# Import raw data
raw_data.infant <- read.csv("../results/infants.tsv", sep = "\t")
# Turn raw into behavioural data, save it to a csv file
behaviour.adults <- LT_data.to_responses(raw_data.adults)
write.csv(behaviour.adults, "../results/BeviouralData.csv")
# Turn raw into clean eyetrackingR data
LT.adults <- raw_data.adults %>%
  LT_data.to_eyetrackingR(AOIs.adults) %>%
  make_eyetrackingr_data(participant_column = "Subject",
                         trial_column = "TrialId",
                         time_column = "NormTimeStamp",
                         trackloss_column = "TrackLoss",
                         aoi_columns = c('Head','Tail'),
                         treat_non_aoi_looks_as_missing = TRUE) %>%
  LT_data.trackloss_clean()