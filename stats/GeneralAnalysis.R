library(eyetrackingR)
library(Matrix)
library(lme4)
library(ggplot2)
library(dplyr)

source("Routines.R")

# GATHER DATA
# Define AOIs
AOIs.adults <- data.frame(name=c("Tail","Head"),L=c(20,400),R=c(220,620),T=c(110,55),B=c(330,255))

# Import raw data
raw_data.adults <- LT_data.import()
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

# ANALYSIS - LOOKING TIME
# Plotting eye-tracking data for all AOIs, averaged across all trials
LT.adults.time_course <- make_time_sequence_data(LT.adults, time_bin_size = 1e-2,
                                             predictor_columns = c("Condition"),
                                             aois = c("Head","Tail"))
LT.adults.time_course.plot <- plot(LT.adults.time_course, predictor_column = "Condition") + 
  theme_light() + coord_cartesian(ylim = c(0,1))
# Analysing and plotting total looking time to each AOI
LT.adults.total_per_AOI <- make_time_window_data(LT.adults, 
                                                 aois=c("Head","Tail"),
                                                 predictor_columns=c("Condition"),
                                                 summarize_by = "Subject")
LT.adults.total_per_AOI.plot <- plot(LT.adults.total_per_AOI, predictor_columns="Condition", dv = "ArcSin")
LT.adults.total_per_AOI.t_test <- t.test(ArcSin ~ Condition, data=LT.adults.total_per_AOI)

# ANALYSIS -- BEHAVIOURAL DATA -- NUMBER OF BLOCKS TO TRAINING
# Creating sub-dataframe for NBlocks
behaviour.adults.n_blocks <- subset(behaviour.adults, !duplicated(Subject))
# Plotting and analysing NBlocks by Condition
behaviour.adults.n_blocks.plot <- ggplot(behaviour.adults.n_blocks, aes(x=Condition,y=LogNBlocks)) + geom_boxplot()
behaviour.adults.n_blocks.t_test <- t.test(LogNBlocks ~ Condition, data=behaviour.adults.n_blocks)
behaviour.adults.n_blocks.lm <- lm(LogNBlocks ~ Condition*Gender*Age,
                                   data = behaviour.adults.n_blocks)
# ANALYSIS -- BEHAVIOURAL DATA -- REACTION TIME
behaviour.adults.reaction_time.plot <- ggplot(behaviour.adults, aes(x=Condition,y=LogRT)) + geom_boxplot()
behaviour.adults.reaction_time.t_test <- t.test(LogRT ~ Condition, data = behaviour.adults)
behaviour.adults.reaction_time.anova <- aov(LogRT ~ Condition*Block, data = behaviour.adults)
behaviour.adults.reaction_time.model_comparison <- drop1(behaviour.adults.reaction_time.anova, .~, test = "F")