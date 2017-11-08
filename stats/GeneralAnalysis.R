library("Matrix")
library("lme4")
library("ggplot2")

source("Routines.R")

LT <- raw.to.ET(LT.data.import())
# LT <- make_eyetrackingr_data(LT, 
#                              participant_column = "Subject",
#                              trial_column = "TrialId",
#                              time_column = "timestamp",
#                              trackloss_column = "trackLoss",
#                              aoi_columns = c('Head','Tail'),
#                              treat_non_aoi_looks_as_missing = TRUE
# )
# 
# LT.clean <- clean_by_trackloss(data = LT,
#                                trial_prop_thresh = .25)
# 
# data_summary <- describe_data(LT.clean, describe_column='Head', group_columns=c("StiLabel"))
# print(data_summary)