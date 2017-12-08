library(eyetrackingR)
library(Matrix)
library(lme4)
library(ggplot2)
library(dplyr)

source("Routines.R")

# GATHER DATA
# Define AOIs
AOIs.adults <- data.frame(name=c("Tail","Head","Feet"),
                          L=c(36,363,146),
                          R=c(260,602,483),
                          T=c(115,66,364),
                          B=c(310,295,465))
AOIs.plot <- ggplot(AOIs.adults, aes(xmin = L, xmax = R, ymin = T, ymax = B)) +
  xlim(c(0,640)) + scale_y_reverse(limits = c(480,0)) +
  geom_rect(aes(fill = name)) +
  theme(legend.position = "top")
ggsave("../results/adults_3f/data_cleaning_graphs/AOIs.png",
       plot = AOIs.plot , width = 3.2, height = 2.95)

# Import raw data
raw_data.adults <- LT_data.adults.import(participants = "adults_3f", subjects = 1:40, pinfo = F)
# Turn raw into behavioural data, save it to a csv file
behaviour.adults <- LT_data.to_responses(raw_data.adults)
write.csv(behaviour.adults, "../results/adults_3f/data/BeviouralData.csv")
# Turn raw into eyetrackingR data
LT.adults <- raw_data.adults %>%
  LT_data.to_eyetrackingR(AOIs.adults) %>%
  make_eyetrackingr_data(participant_column = "Subject",
                         trial_column = "TrialId",
                         time_column = "NormTimeStamp",
                         trackloss_column = "TrackLoss",
                         aoi_columns = c('Head','Tail','Feet'),
                         treat_non_aoi_looks_as_missing = F)
LT.adults$NAOI = !LT.adults$TrackLoss & !(LT.adults$Head | LT.adults$Tail)
# Check for trackloss ratio and NAOI (non-AOI) ratio
LT.adults.gaze_summary <- LT.adults %>%
  group_by(Subject,CurrentObject) %>%
  summarise(TrackLossRatio = sum(TrackLoss)/n(),
            NAOIRatio = sum(NAOI)/(n()-sum(TrackLoss)))
LT.adults.gaze_summary.plot.TrackLossRatio <- ggplot(LT.adults.gaze_summary,
                                                     aes(x = CurrentObject, y = TrackLossRatio)) +
  geom_violin(aes(fill = CurrentObject)) +
  geom_boxplot(alpha=0, width=.3, outlier.alpha = 1) +
  guides(fill = "none")
ggsave("../results/adults_3f/data_cleaning_graphs/TrackLossRatio.png")
LT.adults.gaze_summary.plot.NAOIRatio <- ggplot(LT.adults.gaze_summary,
                                                aes(x = CurrentObject, y = NAOIRatio)) +
  geom_violin(aes(fill = CurrentObject)) +
  geom_boxplot(alpha=0, width=.3, outlier.alpha = 1) +
  guides(fill = "none")
ggsave("../results/adults_3f/data_cleaning_graphs/NonAOIRatio.png")
# Make clean
LT.adults.clean <- LT_data.trackloss_clean(LT.adults, trial_prop_thresh = .3, incl_crit = .5,
                                           participants = "adults_3f")
LT.adults.clean$TrialId <- as.numeric(LT.adults.clean$TrialId)

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
