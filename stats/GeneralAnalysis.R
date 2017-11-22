library(eyetrackingR)
library(Matrix)
library(lme4)
library(ggplot2)
library(dplyr)

source("Routines.R")

# GATHER DATA
# Define AOIs
AOIs.adults <- data.frame(name=c("Tail","Head"),L=c(20,400),R=c(220,620),T=c(110,55),B=c(330,255))

# Import data and transform into eyetrackingR data
LT.adults <- LT.data.import() %>% raw.to.ET(AOIs.adults) %>% make_eyetrackingr_data(participant_column = "Subject",
                                                                                    trial_column = "TrialId",
                                                                                    time_column = "TimeStamp",
                                                                                    trackloss_column = "TrackLoss",
                                                                                    aoi_columns = c('Head','Tail'),
                                                                                    treat_non_aoi_looks_as_missing = TRUE)

# QUALITY CHECK

