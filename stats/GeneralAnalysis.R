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
                                                                                    time_column = "timestamp",
                                                                                    trackloss_column = "track_loss",
                                                                                    aoi_columns = c('Head','Tail'),
                                                                                    treat_non_aoi_looks_as_missing = TRUE)

# QUALITY CHECK
# Get trackloss information
trackloss <- trackloss_analysis(LT.adults)
trackloss.subject.trial <- unique(trackloss[, c('Subject','TrialId','TracklossForTrial')])
# Plot trackloss per trial per subject
trackloss.subject.trial.p <- ggplot(trackloss.subject, aes(x=TrialId, y=TracklossForTrial)) +
  facet_wrap(~Subject, nrow = 10, scales = "free_x") + geom_point()
ggsave("../results/TracklossSubjectTrial.png", trackloss.subject.trial.p, width=9, height=15)
# Remove trials with trackloss proportion greater than 0.25
LT.adults.trackloss <- clean_by_trackloss(data = LT.adults,
                                          trial_prop_thresh = .25)
# Compute and plot proportion of valid trials per subject (number of valid trials / number of trials for subject)
LT.adults.described <- describe_data(LT.adults.trackloss, 'Block', 'Subject')
LT.adults.described$ProportionTrials <- LT.adults.described$NumTrials /
                                          (12*(LT.adults.described$Max + 1))
LT.adults.described$AboveCriteria <- factor(ifelse(LT.adults.described$ProportionTrials >= .5, 1, 0)) # Change .5 to actual inclusion criteria
print(summary(LT.adults.described))
LT.adults.described.p <- ggplot(LT.adults.described,
                                aes(x=Subject, y=ProportionTrials, colour = AboveCriteria)) +
  scale_colour_manual(values = c("red","green"), guide = F) + geom_point()
ggsave("../results/ProportionTrialPerSubject.png", LT.adults.described.p, width=10, height=3)
# Select subjects to keep
LT.adults.clean <- LT.adults.trackloss[LT.adults.trackloss$Subject %in%
                                         LT.adults.described$Subject[LT.adults.described$AboveCriteria == 1],] %>%
  droplevels()
# Check how many subjects missing per condition
print(summary(unique(LT.adults.clean[,c('Subject','Condition')])))
