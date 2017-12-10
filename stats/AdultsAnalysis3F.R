library(eyetrackingR)
library(Matrix)
library(lme4)
library(tidyverse)

source("Routines.R")

# GATHER DATA
# Define AOIs
# AOIs.plot <- ggplot(AOIs.adults, aes(xmin = L, xmax = R, ymin = T, ymax = B)) +
#   xlim(c(0,640)) + scale_y_reverse(limits = c(480,0)) +
#   geom_rect(aes(fill = name)) +
#   theme(legend.position = "top")
# ggsave("../results/adults_3f/data_cleaning_graphs/AOIs.png",
#        plot = AOIs.plot , width = 3.2, height = 2.95)

LT.clean <- LT_data.gather("adults_3f")

# ANALYSIS - LOOKING TIME
# Plotting heatmap for label and no-label participants
LT.adults.heatmap <- ggplot(LT.adults.clean, aes(x=CursorX,y=CursorY)) +
  xlim(c(0,640)) + scale_y_reverse(limits = c(480,0)) +
  facet_wrap(~Condition) +
  geom_bin2d(binwidth = c(20,20))
# Plotting eye-tracking data for all AOIs, averaged across all trials
LT.adults.time_course <- make_time_sequence_data(LT.adults, time_bin_size = 1e-2,
                                                 predictor_columns = c("Condition"),
                                                 aois = c("Head","Tail","Feet"))
LT.adults.time_course.plot <- plot(LT.adults.time_course, predictor_column = "Condition") +
  theme_light()
ggsave("../results/adults_3f/LookingTimeCourseNorm.pdf", plot = LT.adults.time_course.plot)
# Analysing and plotting total looking time to each AOI
# Making data time-window-analysis ready
LT.adults.total_per_AOI <- make_time_window_data(LT.adults,
                                                 aois=c("Head","Tail","Feet"),
                                                 predictor_columns=c("Condition",
                                                                     "Stimulus",
                                                                     "TrialId",
                                                                     "CategoryName"))
# Boxplots of total looking time per AOI
LT.adults.total_per_AOI.plot <- plot(LT.adults.total_per_AOI,
                                     predictor_columns=c("Condition"),
                                     dv = "ArcSin")
ggsave("../results/adults_3f/TotalPerAOI.png", plot = LT.adults.total_per_AOI.plot)
# Simple t-test
LT.adults.total_per_AOI.t_test <- t.test(ArcSin ~ Condition, data=LT.adults.total_per_AOI)

# ANALYSIS -- BEHAVIOURAL DATA -- NUMBER OF BLOCKS TO TRAINING
# Creating sub-dataframe for NBlocks
behaviour.adults.n_blocks <- subset(behaviour.adults, !duplicated(Subject))
# Plotting and analysing NBlocks by Condition
behaviour.adults.n_blocks.plot <- ggplot(behaviour.adults.n_blocks, aes(x=Condition,y=LogNBlocks)) + geom_boxplot()
ggsave("../results/adults_3f/BlocksToLearning.png", plot = behaviour.adults.n_blocks.plot)
behaviour.adults.n_blocks.t_test <- t.test(LogNBlocks ~ Condition, data=behaviour.adults.n_blocks)
# ANALYSIS -- BEHAVIOURAL DATA -- REACTION TIME
behaviour.adults.reaction_time.plot <- ggplot(behaviour.adults, aes(x=Condition,y=LogRT)) + geom_boxplot()
ggsave("../results/adults_3f/ReactionTime.png", plot = behaviour.adults.reaction_time.plot)
behaviour.adults.reaction_time.t_test <- t.test(LogRT ~ Condition, data = behaviour.adults)
behaviour.adults.reaction_time.anova <- aov(LogRT ~ Condition*Block, data = behaviour.adults)
behaviour.adults.reaction_time.model_comparison <- drop1(behaviour.adults.reaction_time.anova, ~., test = "F")
