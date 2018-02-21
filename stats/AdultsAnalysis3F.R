library(lme4)
library(nortest)
library(tidyverse); library(broom)
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
                         aoi_columns = c("Head","Tail"),
                         treat_non_aoi_looks_as_missing = T) %>%
  subset_by_window(window_start_time = -1500, rezero = F)

# LOOKING TIME ANALYSIS: TIME COURSE ===============================================================
# Preparing data for analysis and plot
LT.time_course_aois <- LT.clean %>%
  make_time_sequence_data(time_bin_size = 50,
                          aois = c("Head","Tail","Feet"),
                          predictor_columns=c("Condition",
                                              "Block",
                                              "NBlocks",
                                              "ACC",
                                              "Stimulus",
                                              "StiLabel")) %>%
  mutate(Part = case_when(Block == 0 ~ "Test",
                          Block == 1 ~ "First Block",
                          Block == NBlocks ~ "Last Block"))
# Growth Curve Analysis of the data
LT.time_course_aois.GCA <- lmer(Prop ~ AOI:(ot1 + ot2 + ot3 + ot4 + ot5 + ot6 + ot7) +
                                  Condition:AOI:(ot1 + ot2 + ot3 + ot4 + ot5 + ot6 + ot7) +
                                  (1 | Block) +
                                  (1 | Stimulus) +
                                  (1 + ot1 + ot2 + ot3 + ot4 + ot5 + ot6 + ot7 | Participant),
                                data = LT.time_course_aois, REML = F,
                                verbose = 2)
# Plotting eye-tracking data for all AOIs, averaged across trials
## Individual plot for each block
LT.clean.time_course.plot.blocks <- ggplot(LT.time_course_aois,
                                           aes(x = Time, y=Prop,
                                               colour=Condition,
                                               fill=Condition)) +
  xlab('Time in Trial') + ylab("Looking to AOI (Prop)") +
  facet_grid(Block~AOI) + theme(legend.position = "top") + ylim(0,1) +
  stat_summary(fun.y='mean', geom='line', linetype = '61') +
  stat_summary(fun.data=mean_se, geom='ribbon', alpha= .25, colour=NA) +
  geom_hline(yintercept = 1/3)
ggsave("../results/adults_3f/LookingTimeCoursePerBlock.pdf",
       plot = LT.clean.time_course.plot.blocks,
       width = 7, height = 30)
## Plot for first block, last block, and test
LT.time_course_aois.first_last <- LT.time_course_aois %>%
  drop_na(Part)
intercept <- tibble(Part = c(rep("First Block", 2), rep("Last Block", 2)),
                    x_int = c(0, 2000, 0, 2000))
LT.clean.time_course.first_last.plot <- ggplot(LT.time_course_aois.first_last,
                                               aes(x = Time, y=Prop,
                                                   colour=Condition,
                                                   fill=Condition)) +
  xlab('Time in Trial') + ylab("Looking to AOI (Prop)") +
  facet_grid(AOI~Part, scales = "free_x") + theme(legend.position = "top") + ylim(0,1) +
  geom_vline(data = intercept, aes(xintercept = x_int), linetype = "62", alpha = .5) +
  stat_summary(fun.y='mean', geom='line', linetype = '61') +
  stat_summary(fun.data=mean_se, geom='ribbon', alpha= .25, colour=NA) +
  geom_hline(yintercept = 1/3)
ggsave("../results/adults_3f/LookingTimeCourseFirstLast.pdf",
       plot = LT.clean.time_course.first_last.plot,
       width = 7, height = 10)
## Global plot
LT.clean.time_course.plot <- ggplot(subset(LT.time_course_aois, Part != "Test"),
                                    aes(x = Time, y=Prop,
                                        colour=Condition,
                                        fill=Condition)) +
  xlab('Time in Trial') + ylab("Looking to AOI (Prop)") +
  facet_wrap(~AOI) + theme(legend.position = "top") + ylim(0,1) +
  stat_summary(fun.y='mean', geom='line', linetype = '61') +
  stat_summary(fun.data=mean_se, geom='ribbon', alpha= .25, colour=NA) +
  geom_hline(yintercept = 1/3)
ggsave("../results/adults_3f/LookingTimeCourse.pdf",
       plot = LT.clean.time_course.plot,
       width = 7, height = 10)
# LOOKING TIME ANALYSIS: PROP AOI LOOKING BY PARTICIPANT BY BLOCK ==================================
# Prepare dataset with BlockTransformation
LT.prop_aois_per_block <- make_time_window_data(LT.clean,
                                                aois=c("Tail","Feet","Head"),
                                                predictor_columns=c("Condition",
                                                                    "Block",
                                                                    "NBlocks",
                                                                    "ACC")) %>%
  subset(Block > 0) %>%
  group_by(Participant) %>%
  mutate(OppBlock = Block - NBlocks,
         NormBlock = Block/NBlocks) %>%
  gather("BlockTransformation","Block", Block, OppBlock, NormBlock)
# Growth Curve Analysis (GCA) for each AOI for each BlockTransformation
## Set orthogonal polynomials for GCA
blocks <- LT.prop_aois_per_block %>%
  {sort(as.vector(unique(.$Block)))}
orth_poly <- poly(blocks, 7) %>%
  as.tibble() %>%
  mutate(Block = blocks)
colnames(orth_poly) <- c(paste0("ot", 1:7), "Block")
LT.prop_aois_per_block <- left_join(LT.prop_aois_per_block, orth_poly)
## Run lmer model for GCA, for each AOI: create subset, run model, drop effects and test sig.
### Feet
LT.prop_feet_per_block <- LT.prop_aois_per_block %>%
  subset(AOI == "Feet")
LT.prop_feet_per_block.GCA <- lmer(ArcSin ~ Condition * (ot1 + ot2 + ot3 + ot4 + ot5 + ot6 + ot7) +
                                     (1 | TrialId) +
                                     (1 | Participant),
                                   data = LT.prop_feet_per_block, REML = FALSE)
LT.prop_feet_per_block.comp <- drop1(LT.prop_feet_per_block.GCA, ~., test = "Chi")
LT.prop_feet_per_block.GCA <- lmer(ArcSin ~ Condition:ot5 +
                                     (1 | TrialId) +
                                     (1 | Participant),
                                   data = LT.prop_feet_per_block, REML = FALSE)
### Tail
LT.prop_tail_per_block <- LT.prop_aois_per_block %>%
  subset(AOI == "Tail")
LT.prop_tail_per_block.GCA <- lmer(ArcSin ~ Condition * (ot1 + ot2 + ot3 + ot4 + ot5 + ot6 + ot7) +
                                     (1 | TrialId) +
                                     (1 | Participant),
                                   data = LT.prop_tail_per_block, REML = FALSE)
LT.prop_tail_per_block.comp <- drop1(LT.prop_tail_per_block.GCA, ~., test = "Chi")
LT.prop_tail_per_block.GCA <- lmer(ArcSin ~ Condition:(ot1 + ot2 + ot7) +
                                     (1 | TrialId) +
                                     (1 | Participant),
                                   data = LT.prop_tail_per_block, REML = FALSE)
### Head
LT.prop_head_per_block <- LT.prop_aois_per_block %>%
  subset(AOI == "Head")
LT.prop_head_per_block.GCA <- lmer(ArcSin ~ Condition * (ot1 + ot2 + ot3 + ot4 + ot5 + ot6 + ot7) +
                                     (1 | TrialId) +
                                     (1 | Participant),
                                   data = LT.prop_head_per_block, REML = FALSE)
LT.prop_head_per_block.comp <- drop1(LT.prop_head_per_block.GCA, ~., test = "Chi")
LT.prop_head_per_block.GCA <- lmer(ArcSin ~ Condition:(ot1 + ot2 + ot7) +
                                     (1 | TrialId) +
                                     (1 | Participant),
                                   data = LT.prop_head_per_block, REML = FALSE)
# Plot all AOIs and all BlockTransformations
## Get predicted values from model, add to initial data frame
LT.prop_feet_per_block$Predicted <- predict(LT.prop_feet_per_block.GCA,
                                            LT.prop_feet_per_block,
                                            re.form = NA)
LT.prop_head_per_block$Predicted <- predict(LT.prop_head_per_block.GCA,
                                            LT.prop_head_per_block,
                                            re.form = NA)
LT.prop_tail_per_block$Predicted <- predict(LT.prop_tail_per_block.GCA,
                                            LT.prop_tail_per_block,
                                            re.form = NA)
LT.prop_aois_per_block <- rbind(LT.prop_feet_per_block,
                                LT.prop_head_per_block,
                                LT.prop_tail_per_block)
## Plot raw data and predicte values from model
LT.prop_aois_per_block.plot <- ggplot(LT.prop_aois_per_block,
                                      aes(x = Block, y = ArcSin,
                                          colour = Condition,
                                          fill = Condition)) +
  facet_wrap(~AOI) + theme(aspect.ratio = 1.618/1, legend.position = "top") +
  #stat_smooth(linetype = "61", level = 0.87) +
  stat_summary(fun.y = 'mean', geom = 'line', linetype = '61') +
  stat_summary(fun.data = 'mean_se', geom = 'ribbon', alpha = .25, colour = NA) +
  stat_summary(aes(y = Predicted), fun.y = 'mean', geom = "line", size = 1.2) +
  geom_hline(yintercept = asin(sqrt(1/3)))
ggsave("../results/adults_3f/AOILookingEvolution.png",
       plot = LT.prop_aois_per_block.plot,
       width = 7, height = 5)
# Comparing first block agains last block
## Plot jitter + mean&se + lines
LT.prop_aois.first_last <- LT.prop_aois_per_block %>%
  subset(BlockTransformation == "Block" &
           (Block == 1 | Block == NBlocks)) %>%
  mutate(Part = ifelse(Block == 1, "First Block", "Last Block"))
LT.prop_aois.first_last.plot <- ggplot(LT.prop_aois.first_last,
                                       aes(x = Part, y = Prop,
                                           colour = Condition,
                                           fill = Condition)) +
  facet_wrap(~AOI) + theme(aspect.ratio = 1.618/1, legend.position = "top") +
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
       width = 9, height = 6)
## LMER for Prop ~ Condition*Part*AOI
LT.prop_aois.first_last.lmer.0 <- lmer(ArcSin ~ Part*AOI*Condition - Part +
                                         (1 + AOI | Participant),
                                       data = LT.prop_aois.first_last)
# No main effect of Part since looking at proportions,
# so overall they should all equate to one when collapsing
LT.prop_aois.first_last.lmer.1 <- update(LT.prop_aois.first_last.lmer.0,
                                         . ~ . - Part:AOI:Condition) # Remove
LT.prop_aois.first_last.lmer.2 <- update(LT.prop_aois.first_last.lmer.1,
                                         . ~ . - AOI:Condition)
LT.prop_aois.first_last.lmer.3 <- update(LT.prop_aois.first_last.lmer.2,
                                         . ~ . - Part:Condition) # Remove
LT.prop_aois.first_last.lmer.4 <- update(LT.prop_aois.first_last.lmer.3,
                                         . ~ . - Part:AOI)
LT.prop_aois.first_last.lmer.5 <- update(LT.prop_aois.first_last.lmer.4,
                                         . ~ . - Condition) # Remove
LT.prop_aois.first_last.lmer.6 <- update(LT.prop_aois.first_last.lmer.5,
                                         . ~ . - AOI)
LT.prop_aois.first_last.lmer.comp <- anova(LT.prop_aois.first_last.lmer.6,
                                           LT.prop_aois.first_last.lmer.5,
                                           LT.prop_aois.first_last.lmer.4,
                                           LT.prop_aois.first_last.lmer.3,
                                           LT.prop_aois.first_last.lmer.2,
                                           LT.prop_aois.first_last.lmer.1,
                                           LT.prop_aois.first_last.lmer.0)
LT.prop_aois.first_last.lmer.final <- update(LT.prop_aois.first_last.lmer.0,
                                             . ~ . - (Condition +
                                                        Part:Condition +
                                                        Part:AOI:Condition))

# BEHAVIOURAL ANALYSIS: PARTICIPANTS AND BLOCKS ====================================================
# Get how many participants for each block, and make a bar plot
behaviour.parts_per_block <- behaviour %>%
  subset(Block > 0) %>%
  group_by(Block, Condition) %>%
  summarise(N_Participants = n_distinct(Participant))
behaviour.parts_per_block.plot <- ggplot(behaviour.parts_per_block,
                                         aes(x = Block, y = N_Participants, fill = Condition)) +
  geom_col(position = "dodge")
ggsave("../results/adults_3f/ParticipantsPerBlock.png",
       plot = behaviour.parts_per_block.plot)
# Get number of blocks to learning per participant, plot a violin
behaviour.blocks_per_part <- behaviour %>%
  select(c(Participant, Condition, NBlocks)) %>%
  unique()
behaviour.blocks_per_part.plot <- ggplot(behaviour.blocks_per_part,
                                         aes(x = Condition, y = NBlocks, fill = Condition)) +
  geom_violin() +
  geom_boxplot(alpha = 0, outlier.alpha = 1, width = .15)
ggsave("../results/adults_3f/BlocksPerParticipant.png",
       plot = behaviour.blocks_per_part.plot)
# BEHAVIOURAL ANALYSIS: ACCURACY ~ CONDITION*(DIAG || RT) ==========================================
# Get datasets for training and test
behaviour.training <- behaviour %>%
  subset(Phase == "Familiarisation")
behaviour.test <- behaviour %>%
  subset(Phase == "Test")
# Run binomial glmer
## During training
ACC_by_diag.training.glmer <- glmer(ACC ~ Condition*Diagnostic*zLogRT +
                                      (1 + Diagnostic + zLogRT | Participant),
                                    family = binomial,
                                    control = glmerControl(optimizer = "bobyqa"),
                                    data = behaviour.training)
## At test !!! NOT CONVERGING !!! ==> All participants at ceiling
##ACC_by_diag.test.glmer <- glm(ACC ~ Condition*Diagnostic +
##                                  (1 | Participant),
##                                family = binomial,
##                                data = behaviour.test)
# Prepare and plot data
## During training
ACC_by_diag.training <- behaviour.training %>%
  group_by(Participant, Diagnostic, Condition) %>%
  summarise(Accuracy = sum(ACC)/n())
ACC_by_diag.training.plot <- ggplot(ACC_by_diag.training,
                                    aes(x = Condition,
                                        y = Accuracy,
                                        fill = Diagnostic)) +
  ylim(0,1) +
  geom_violin() +
  geom_boxplot(alpha = 0, outlier.alpha = 1,
               width = .15, position = position_dodge(.9))
## At test
ACC_by_diag.test <- behaviour.test %>%
  group_by(Participant, Diagnostic, Condition) %>%
  summarise(Accuracy = sum(ACC)/n())
ACC_by_diag.test.plot <- ggplot(ACC_by_diag.test,
                                aes(x = Condition,
                                    y = Accuracy,
                                    fill = Diagnostic)) +
  ylim(0,1) +
  geom_violin() +
  geom_boxplot(alpha = 0, outlier.alpha = 1,
               width = .15, position = position_dodge(.9))
