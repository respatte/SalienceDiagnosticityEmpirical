library(lme4)
library(nortest)
library(tidyverse); library(broom)
library(jtools)
library(eyetrackingR)

source("Routines.R")

# GATHER DATA ======================================================================================
d <- LT_data.gather("adults_3f")
# Get behavioural and LT data, excluding outliers in terms of
# number of blocks before learning (graphically, from boxplot)
behaviour <- d[[2]] %>%
  subset((Condition == "Label" & NBlocks < 21) |
           (Condition == "NoLabel" & NBlocks < 10))
LT.clean <- d[[4]] %>%
  subset((Condition == "Label" & NBlocks < 21) |
           (Condition == "NoLabel" & NBlocks < 10)) %>%
  make_eyetrackingr_data(participant_column = "Participant",
                         trial_column = "TrialId",
                         time_column = "TimeStamp",
                         trackloss_column = "TrackLoss",
                         aoi_columns = c("Head","Tail", "Feet"),
                         treat_non_aoi_looks_as_missing = T) %>%
  subset_by_window(window_start_time = -1500, rezero = F)

# LOOKING TIME ANALYSIS: TIME COURSE ===============================================================
# Preparing data for analysis and plot
LT.time_course_aois.first_last <- LT.clean %>%
  make_time_sequence_data(time_bin_size = 50,
                          aois = c("Head","Tail","Feet"),
                          predictor_columns=c("Condition",
                                              "Block",
                                              "NBlocks",
                                              "ACC",
                                              "Stimulus",
                                              "StiLabel",
                                              "Diagnostic")) %>%
  mutate(Part = case_when(Block == 0 ~ "Test",
                          Block == 1 ~ "First Block",
                          Block == NBlocks ~ "Last Block")) %>%
  drop_na(Part)
# # Growth Curve Analysis of the data
# # Analysing proportions => main effect of Condition or Part nonsensical,
# # we can only expect differences between AOIs, and between AOIs on different levels
# LT.time_course_aois.GCA <- lmer(ArcSin ~ (AOI + Condition:AOI + AOI:Part + AOI:Condition:Part)*
#                                   (ot1 + ot2 + ot3 + ot4 + ot5 + ot6 + ot7) +
#                                   (1 + AOI + Part + Condition +
#                                      ot1 + ot2 + ot3 + ot4 + ot5 + ot6 + ot7 | Stimulus) +
#                                   (1 + AOI + Part + Condition +
#                                      ot1 + ot2 + ot3 + ot4 + ot5 + ot6 + ot7 | StiLabel) +
#                                   (1 + AOI + Part +
#                                      ot1 + ot2 + ot3 + ot4 + ot5 + ot6 + ot7 | Participant),
#                                 data = LT.time_course_aois.first_last, REML = F)
# ## Check convergence
# ### Recompute gradient and Hessian with Richardson extrapolation
# devfun <- update(LT.time_course_aois.GCA, devFunOnly=TRUE)
# pars <- getME(LT.time_course_aois.GCA, "theta")
# if(require("numDeriv")){
#   cat("hess:\n"); print(hess <- hessian(devfun, unlist(pars)))
#   cat("grad:\n"); print(grad <- grad(devfun, unlist(pars)))
#   cat("scaled gradient:\n")
#   print(scgrad <- solve(chol(hess), grad))
# }
# print(LT.time_course_aois.GCA@optinfo$derivs)
# ### Restart the fit from computed values
# LT.time_course_aois.GCA.restart.1 <- update(LT.time_course_aois.GCA, start = pars)
# Plotting eye-tracking data and GCA predictions for all AOIs, for first block, last block, and test
intercept <- tibble(Part = c(rep("First Block", 2), rep("Last Block", 2)),
                    x_int = c(0, 2000, 0, 2000))
LT.clean.time_course.first_last.plot <- ggplot(LT.time_course_aois.first_last,
                                               aes(x = Time, y=Prop,
                                                   colour=Condition,
                                                   fill=Condition)) +
  xlab('Time in Trial') + ylab("Looking to AOI (Prop)") + theme_apa(legend.pos = "top") +
  scale_colour_discrete(labels = c("Label", "No Label")) +
  scale_fill_discrete(labels = c("Label", "No Label")) +
  facet_grid(AOI~Part, scales = "free_x") + ylim(0,1) +
  scale_x_continuous(breaks = c(-1000, 0, 1000, 2000, 3000)) +
  geom_vline(data = intercept, aes(xintercept = x_int), linetype = "62", alpha = .5) +
  stat_summary(fun.y='mean', geom='line', linetype = '61') +
  stat_summary(fun.data=mean_se, geom='ribbon', alpha= .33, colour=NA) +
  geom_hline(yintercept = 1/3)
ggsave("../results/adults_3f/LookingTimeCourseFirstLast.pdf",
       plot = LT.clean.time_course.first_last.plot,
       width = 7, height = 5.4)

# LOOKING TIME ANALYSIS: PROP AOI LOOKING BY PARTICIPANT BY PART ===================================
# Prepare dataset with BlockTransformation
LT.prop_aois.per_block <- make_time_window_data(LT.clean,
                                                aois=c("Tail","Feet","Head"),
                                                predictor_columns=c("Condition",
                                                                    "Block",
                                                                    "NBlocks",
                                                                    "ACC",
                                                                    "Stimulus",
                                                                    "StiLabel")) %>%
  subset(Block > 0) %>%
  group_by(Participant) %>%
  mutate(OppBlock = Block - NBlocks,
         NormBlock = Block/NBlocks) %>%
  gather("BlockTransformation","Block", Block, OppBlock, NormBlock)
# Comparing first block against last block
LT.prop_aois.first_last <- LT.prop_aois.per_block %>%
  subset(BlockTransformation == "Block" &
           (Block == 1 | Block == NBlocks)) %>%
  mutate(Part = ifelse(Block == 1, "First Block", "Last Block"))
## LMER for Prop ~ Condition*Part*AOI
# No main effect of Part since looking at proportions,
# so overall they should all equate to one when collapsing
LT.prop_aois.first_last.lmer.0 <- lmer(ArcSin ~ Part*AOI*Condition - (Part + Condition) +
                                         (1 + AOI + Part:AOI | Participant) +
                                         (1 + AOI + Part:AOI | Stimulus) +
                                         (1 + AOI + Part:AOI | StiLabel),
                                       data = LT.prop_aois.first_last)
LT.prop_aois.first_last.lmer.1 <- update(LT.prop_aois.first_last.lmer.0,
                                         . ~ . - Part:AOI:Condition) # Remove
LT.prop_aois.first_last.lmer.2 <- update(LT.prop_aois.first_last.lmer.1,
                                         . ~ . - AOI:Condition)
LT.prop_aois.first_last.lmer.3 <- update(LT.prop_aois.first_last.lmer.2,
                                         . ~ . - Part:Condition) # Remove
LT.prop_aois.first_last.lmer.4 <- update(LT.prop_aois.first_last.lmer.3,
                                         . ~ . - Part:AOI)
LT.prop_aois.first_last.lmer.5 <- update(LT.prop_aois.first_last.lmer.4,
                                         . ~ . - AOI)
LT.prop_aois.first_last.lmer.comp <- anova(LT.prop_aois.first_last.lmer.5,
                                           LT.prop_aois.first_last.lmer.4,
                                           LT.prop_aois.first_last.lmer.3,
                                           LT.prop_aois.first_last.lmer.2,
                                           LT.prop_aois.first_last.lmer.1,
                                           LT.prop_aois.first_last.lmer.0)
LT.prop_aois.first_last.lmer.final <- update(LT.prop_aois.first_last.lmer.0,
                                             . ~ . - (Part:Condition + Part:AOI:Condition))
## Plot jitter + lmer mean&se + lines
LT.prop_aois.first_last.plot <- ggplot(LT.prop_aois.first_last,
                                       aes(x = Part, y = Prop,
                                           colour = Condition,
                                           fill = Condition)) +
  facet_wrap(~AOI) + theme_apa(legend.pos = "top") + ylab("Looking to AOI (Prop)") +
  scale_colour_discrete(labels = c("Label", "No Label")) +
  scale_fill_discrete(labels = c("Label", "No Label")) +
  geom_point(position = position_jitterdodge(dodge.width = .8,
                                             jitter.width = .2),
             alpha = .25) +
  geom_errorbar(stat = "summary",
                width = .2, colour = "black",
                position = position_dodge(.1)) +
  geom_line(aes(x = Part, y = Prop, group = Condition),
            stat = "summary", fun.y = "mean",
            colour = "black") +
  geom_point(stat = "summary", fun.y = "mean",
             shape = 18, size = 3,
             position = position_dodge(.1))
ggsave("../results/adults_3f/AOILookingFirstLast.pdf",
       LT.prop_aois.first_last.plot,
       width = 7, height = 5.4)

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
ggsave("../results/adults_3f/ParticipantsPerBlock.pdf", plot = behaviour.parts_per_block.plot,
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
ggsave("../results/adults_3f/BlocksPerParticipant.pdf", plot = behaviour.blocks_per_part.plot,
       width = 3.5, height = 2.9)
# Stats for the number of blocks per participant
NBlocksLabel <- behaviour.blocks_per_part %>%
  subset(Condition == "Label", select = NBlocks) %>%
  unlist()
NBlocksNoLabel <- behaviour.blocks_per_part %>%
  subset(Condition == "NoLabel", select = NBlocks) %>%
  unlist()
behaviour.blocks_per_part.t_test <- ks.test(NBlocksLabel, NBlocksNoLabel)

# BEHAVIOURAL ANALYSIS: ACCURACY ~ CONDITION*DIAG*RT) ==============================================
# Get datasets for training and test
behaviour.training <- behaviour %>%
  subset(Phase == "Familiarisation")
behaviour.test <- behaviour %>%
  subset(Phase == "Test")
# Run binomial glmer
## During training
ACC_by_diag_by_RT.training.glmer <- glmer(ACC ~ Condition*Diagnostic*zLogRT +
                                      (1 + Diagnostic + zLogRT | Participant) +
                                      (1 + Diagnostic + zLogRT | Stimulus) +
                                      (1 + Diagnostic + zLogRT | StiLabel),
                                    family = binomial,
                                    control = glmerControl(optimizer = "bobyqa"),
                                    data = behaviour.training)
## At test !!! NOT CONVERGING !!! ==> All participants at ceiling
##ACC_by_diag_by_RT.test.glmer <- glm(ACC ~ Condition*Diagnostic +
##                                  (1 | Participant),
##                                family = binomial,
##                                data = behaviour.test)
# Prepare and plot data
## During training
ACC_by_diag_by_RT.training <- behaviour.training %>%
  group_by(Participant, Diagnostic, Condition) %>%
  summarise(Accuracy = sum(ACC)/n())
ACC_by_diag_by_RT.training.plot <- ggplot(ACC_by_diag_by_RT.training,
                                    aes(x = Condition,
                                        y = Accuracy,
                                        fill = Diagnostic)) +
  ylim(0,1) + theme_apa(legend.pos = "bottomright") +
  scale_x_discrete(labels = c("Label", "No Label")) +
  geom_violin() +
  geom_boxplot(alpha = 0, outlier.alpha = 1,
               width = .15, position = position_dodge(.9))
ggsave("../results/adults_3f/ACCbyRTbyDiag_training.pdf", plot = ACC_by_diag_by_RT.training.plot,
       width = 3.5, height = 2.7)
## At test
ACC_by_diag_by_RT.test <- behaviour.test %>%
  group_by(Participant, Diagnostic, Condition) %>%
  summarise(Accuracy = sum(ACC)/n())
ACC_by_diag_by_RT.test.plot <- ggplot(ACC_by_diag_by_RT.test,
                                aes(x = Condition,
                                    y = Accuracy,
                                    fill = Diagnostic)) +
  ylim(0,1) + theme_apa(legend.pos = "bottomright") +
  scale_x_discrete(labels = c("Label", "No Label")) +
  geom_violin() +
  geom_boxplot(alpha = 0, outlier.alpha = 1,
               width = .15, position = position_dodge(.9))
ggsave("../results/adults_3f/ACCbyRTbyDiag_test.pdf", plot = ACC_by_diag_by_RT.test.plot,
       width = 3.5, height = 2.7)

# BEHAVIOURAL ANALYSIS: RT ~ CONFIDENCE ============================================================
RT_by_conf.lmer.0 <- lmer(RT ~ Confidence*Condition +
                            (1 + Confidence | Participant) +
                            (1 + Confidence | Stimulus) +
                            (1 + Confidence | StiLabel),
                          data = behaviour.training)
RT_by_conf.lmer.1 <- update(RT_by_conf.lmer.0, . ~ . - Confidence:Condition)
RT_by_conf.lmer.2 <- update(RT_by_conf.lmer.1, . ~ . - Confidence)
RT_by_conf.lmer.3 <- update(RT_by_conf.lmer.2, . ~ . - Condition)
RT_by_conf.lmer.comp <- anova(RT_by_conf.lmer.3,
                              RT_by_conf.lmer.2,
                              RT_by_conf.lmer.1,
                              RT_by_conf.lmer.0)
# No significant effect