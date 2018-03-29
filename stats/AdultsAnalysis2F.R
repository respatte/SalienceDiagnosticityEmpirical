# LIBRARY IMPORTS ==================================================================================
library(lme4);library(lmerTest)
library(nortest)
library(tidyverse); library(broom)
library(jtools)
library(eyetrackingR)
library(rethinking)
library(texreg)

source("Routines.R")

# GATHER DATA ======================================================================================
d <- LT_data.gather("adults_2f")
# Get behavioural and LT data, excluding outliers in terms of
# number of blocks before learning (graphically, from boxplot)
behaviour <- d[[2]] %>%
  subset((Condition == "Label" & NBlocks < 8) |
           (Condition == "NoLabel" & NBlocks < 5))
rote_learning_check <- behaviour %>%
  subset(Phase == "Test") %>%
  group_by(Participant, Stimulus) %>%
  summarise(a = ACC) %>%
  subset(a == 0) # Participants with mistakes at test (including missed item)
# No rote learners
LT.clean <- d[[4]] %>%
  subset((Condition == "Label" & NBlocks < 8) |
           (Condition == "NoLabel" & NBlocks < 5)) %>%
  make_eyetrackingr_data(participant_column = "Participant",
                         trial_column = "TrialId",
                         time_column = "TimeStamp",
                         trackloss_column = "TrackLoss",
                         aoi_columns = c("Head", "Tail"),
                         treat_non_aoi_looks_as_missing = T) %>%
  subset_by_window(window_start_time = -1000, rezero = F)

# LOOKING TIME ANALYSIS: TIME COURSE ===============================================================
# DATA PREPARATION
LT.time_course_tail.first_last <- LT.clean %>%
  make_time_sequence_data(time_bin_size = 50,
                          aois = c("Tail"),
                          predictor_columns=c("Condition",
                                              "FstLst",
                                              "ACC",
                                              "Stimulus",
                                              "StimLabel")) %>%
  drop_na(FstLst)
# GROWTH CURVE ANALYSIS
run_model = F # Running the model takes around 10 hours on a [check office CPU specs]
if(run_model){
  ## Run and save the model
  #- Analysing proportions => main effect of Condition or Part nonsensical,
  #- we can only expect differences between AOIs, and between AOIs on different levels
  LT.time_course_tail.GCA <- lmer(ArcSin ~ (Condition*FstLst)*
                                    (ot1 + ot2 + ot3 + ot4 + ot5 + ot6 + ot7) +
                                    (1 + Part +
                                       ot1 + ot2 + ot3 + ot4 + ot5 + ot6 + ot7 | Participant) +
                                    (1 | Stimulus) +
                                    (1| StimLabel),
                                  data = LT.time_course_tail.first_last, REML = F,
                                  control = lmerControl(optCtrl = list(maxfun = 100000)))
  saveRDS(LT.time_course_tail.GCA, file = "../results/adults_2f/GCA.rds")
  ## Run and save the ANOVA for the model effects
  LT.time_course_tail.GCA.anova <- anova(LT.time_course_tail.GCA, type = 2)
  saveRDS(LT.time_course_tail.GCA, file = "../results/adults_2f/GCA_anova.rds")
}else{
  LT.time_course_tail.GCA <- readRDS("../results/adults_2f/GCA.rds")
  LT.time_course_tail.GCA.anova <- readRDS("../results/adults_2f/GCA_anova.rds")
}
# PLOTTING
## Plot eye-tracking data and GCA predictions for all AOIs, for first block, last block, and test
intercept <- tibble(Part = c(rep("First Block", 2), rep("Last Block", 2)),
                    x_int = c(0, 2000, 0, 2000))
LT.clean.time_course.first_last.plot <- ggplot(LT.time_course_tail.first_last,
                                               aes(x = Time, y=Prop,
                                                   colour=Condition,
                                                   fill=Condition)) +
  xlab('Time in Trial') + ylab("Looking to AOI (Prop)") + theme_apa(legend.pos = "top") +
  scale_colour_discrete(labels = c("Label", "No Label")) +
  scale_fill_discrete(labels = c("Label", "No Label")) +
  facet_grid(.~FstLst, scales = "free_x") + ylim(0,1) +
  scale_x_continuous(breaks = c(-1000, 0, 1000, 2000, 3000)) +
  geom_vline(data = intercept, aes(xintercept = x_int), linetype = "62", alpha = .5) +
  stat_summary(fun.y='mean', geom='line', linetype = '61') +
  stat_summary(fun.data=mean_se, geom='ribbon', alpha= .33, colour=NA) +
  geom_hline(yintercept = .5)
ggsave("../results/adults_2f/LookingTimeCourseFirstLast.pdf",
       plot = LT.clean.time_course.first_last.plot,
       width = 7, height = 3)

# LOOKING TIME ANALYSIS: PROP AOI LOOKING BY PARTICIPANT BY PART ===================================
# DATA PREPARATION
## Prepare dataset with all blocks
LT.prop_tail.per_block <- make_time_window_data(LT.clean,
                                                aois="Tail",
                                                predictor_columns=c("Condition",
                                                                    "Block",
                                                                    "FstLst",
                                                                    "ACC",
                                                                    "Stimulus",
                                                                    "StimLabel"))
## Comparing first block against last block
LT.prop_tail.first_last <- LT.prop_tail.per_block %>%
  drop_na(FstLst)
# MIXED-EFFECTS MODEL FOR PROP ~ CONDITION*PART
run_model <- F # Running the model takes around 10 minutes on a [check office CPU specs]
if(run_model){
  ## Run and save the model
  LT.prop_tail.first_last.lmer <- lmer(ArcSin ~ FstLst*Condition +
                                         (1 + FstLst | Participant) +
                                         (1 | Stimulus) +
                                         (1 | StimLabel),
                                       data = LT.prop_tail.first_last)
  saveRDS(LT.prop_tail.first_last.lmer, "../results/adults_2f/PropAOI.rds")
  ## Run and save the ANOVA for model effects
  LT.prop_tail.first_last.lmer.anova <- anova(LT.prop_tail.first_last.lmer, type = 2)
  saveRDS(LT.prop_tail.first_last.lmer.anova, "../results/adults_2f/PropAOI_anova.rds")
}else{
  LT.prop_tail.first_last.lmer <- readRDS("../results/adults_2f/PropAOI.rds")
  LT.prop_tail.first_last.lmer.anova <- readRDS("../results/adults_2f/PropAOI_anova.rds")
}
# PLOTTING
## Plot jitter + mean&se + lines
LT.prop_tail.first_last.plot <- ggplot(LT.prop_tail.first_last,
                                       aes(x = FstLst, y = Prop,
                                           colour = Condition,
                                           fill = Condition)) +
  theme_apa(legend.pos = "top") + ylab("Looking to AOI (Prop)") +
  scale_colour_discrete(labels = c("Label", "No Label")) +
  scale_fill_discrete(labels = c("Label", "No Label")) +
  geom_point(size = 1,
             position = position_jitterdodge(dodge.width = .8,
                                             jitter.width = .2),
             alpha = .25) +
  geom_errorbar(stat = "summary",
                width = .2, colour = "black",
                position = position_dodge(.1)) +
  geom_line(aes(x = FstLst, y = Prop, group = Condition),
            stat = "summary", fun.y = "mean",
            colour = "black",
            position = position_dodge(.1)) +
  geom_point(stat = "summary", fun.y = "mean",
             shape = 18, size = 2,
             position = position_dodge(.1))
ggsave("../results/adults_2f/AOILookingFirstLast.pdf",
       LT.prop_tail.first_last.plot,
       width = 3.5, height = 3)

# BEHAVIOURAL ANALYSIS: PARTICIPANTS AND BLOCKS ====================================================
# Get how many participants for each block, and make a bar plot
behaviour.parts_per_block <- behaviour %>%
  subset(Block > 0) %>%
  group_by(Block, Condition) %>%
  summarise(N_Participants = n_distinct(Participant))
behaviour.parts_per_block.plot <- ggplot(behaviour.parts_per_block,
                                         aes(x = Block, y = N_Participants, fill = Condition)) +
  theme_apa(legend.pos = "topright") + scale_fill_discrete(labels = c("Label", "No Label")) +
  scale_x_continuous(breaks = 1:10) + ylab("Participants") +
  geom_col(position = "dodge")
ggsave("../results/adults_2f/ParticipantsPerBlock.pdf", plot = behaviour.parts_per_block.plot,
       width = 3.5, height = 2.7)
# Get number of blocks to learning per participant, plot a violin
behaviour.blocks_per_part <- behaviour %>%
  select(c(Participant, Condition, NBlocks)) %>%
  unique()
behaviour.blocks_per_part.plot <- ggplot(behaviour.blocks_per_part,
                                         aes(x = Condition, y = NBlocks, fill = Condition)) +
  theme_apa(legend.pos = "topright") + ylab("Blocks") +
  scale_fill_discrete(labels = c("Label", "No Label")) +
  geom_violin() +
  geom_boxplot(alpha = 0, outlier.alpha = 1, width = .15)
ggsave("../results/adults_2f/BlocksPerParticipant.pdf", plot = behaviour.blocks_per_part.plot,
       width = 3.5, height = 2.9)
# Stats for the number of blocks per participant
behaviour.blocks_per_part.normality <- ad.test(behaviour.blocks_per_part$NBlocks)
behaviour.blocks_per_part.freq_test <- wilcox.test(NBlocks ~ Condition,
                                                   data = behaviour.blocks_per_part)

# BEHAVIOURAL ANALYSIS: ACCURACY ~ CONDITION*RT) ===================================================
# Get datasets for training and test
behaviour.training <- behaviour %>%
  subset(Phase == "Familiarisation")
behaviour.test <- behaviour %>%
  subset(Phase == "Test")
# Run binomial glmer
## During training
ACC_by_diag_by_RT.training.glmer <- glmer(ACC ~ Condition*zLogRT +
                                            (1 + zLogRT | Participant) +
                                            (1 | Stimulus) +
                                            (1 | StimLabel),
                                          family = binomial,
                                          control = glmerControl(optimizer = "bobyqa"),
                                          data = behaviour.training)
## At test
ACC_by_diag_by_RT.test.glmer <- glmer(ACC ~ Condition +
                                  (1 | Participant),
                                family = binomial,
                                data = behaviour.test)
# Prepare and plot data
## During training
ACC_by_diag_by_RT.training <- behaviour.training %>%
  group_by(Participant, Condition) %>%
  summarise(Accuracy = sum(ACC == 1)/n())
ACC_by_diag_by_RT.training.plot <- ggplot(ACC_by_diag_by_RT.training,
                                          aes(x = Condition,
                                              y = Accuracy,
                                              fill = Condition)) +
  ylim(0,1) + theme_apa(legend.pos = "bottomright") +
  scale_fill_discrete(labels = c("Label", "No Label")) +
  geom_violin() +
  geom_boxplot(alpha = 0, outlier.alpha = 1,
               width = .15, position = position_dodge(.9))
ggsave("../results/adults_2f/ACCbyRTbyBlock_training.pdf", plot = ACC_by_diag_by_RT.training.plot,
       width = 3.5, height = 2.7)
## At test
ACC_by_diag_by_RT.test <- behaviour.test %>%
  group_by(Participant, Condition) %>%
  summarise(Accuracy = sum(ACC == 1)/n())
ACC_by_diag_by_RT.test.plot <- ggplot(ACC_by_diag_by_RT.test,
                                      aes(x = Condition,
                                          y = Accuracy,
                                          fill = Condition)) +
  ylim(0,1) + theme_apa(legend.pos = "bottomright") +
  scale_x_discrete(labels = c("Label", "No Label")) +
  geom_violin() +
  geom_boxplot(alpha = 0, outlier.alpha = 1,
               width = .15, position = position_dodge(.9))
ggsave("../results/adults_2f/ACCbyRT_test.pdf", plot = ACC_by_diag_by_RT.test.plot,
       width = 3.5, height = 2.7)