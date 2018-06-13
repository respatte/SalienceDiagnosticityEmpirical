library(doSNOW)
library(tidyverse)
library(eyetrackingR)

# LOOKING-TIME DATA IMPORT -- ADUTLTS
# Function importing looking time data from all adult participants,
# in the ../results/adults repository by default
LT_data.import.adults <- function(participants="adults_2f"){
  single.file.import <- function(file){
    tmp <- read_tsv(file)
    return(subset(tmp, CurrentObject %in% c("Feedback","Label","Stimulus"),
                  select = c(Subject, CursorX, CursorY,
                             TimestampSec, TimestampMicrosec, TrialId, Block,
                             CRESP, RESP, ACC, RT,
                             CurrentObject, Stimulus, StiLabel)))
  }
  res.repo <- paste0("../results/",participants,"/data/")
  # Getting participant info
  participant_info <- read_csv(paste0(res.repo,"ParticipantInformation.csv"))
  # Reading all participant files
  cl <- makeCluster(4)
  registerDoSNOW(cl)
  file.names <- list.files(path=res.repo, pattern=".gazedata")
  df <- foreach(i=1:length(file.names),
                .combine = "rbind",
                .inorder = F,
                .export = "read_tsv") %dopar%
    single.file.import(paste0(res.repo, file.names[i]))
  stopCluster(cl)
  # Transforming data
  # TODO -- Check what happens when no answer in time: RT == 0 or RT == 10000? ACC? CRESP?
  df %<>% rename(Participant = Subject, StimLabel = StiLabel) %>%
    inner_join(participant_info) %>%
    drop_na(RT, CursorX, CursorY) %>%
    subset(RT > 200) %>%
    mutate(TrialId = ifelse(Block==0, TrialId + 252, TrialId),
           Phase = ifelse(TrialId < 252, "Familiarisation", "Test"),
           TrackLoss = pmin.int(CursorX,CursorY)<0,
           TimeStamp = TimestampMicrosec*1e-3 + TimestampSec*1e3,
           AOI_type = 1) %>%
    group_by(Participant) %>%
    mutate(Condition = if(first(StimLabel) == "NoLabelFeedback"){"NoLabel"}else{"Label"},
           CategoryName = if(first(Condition) == "NoLabel"){"NoName"}else{
             if((grepl("A",as.character(first(Stimulus))) &
                 first(StimLabel) == "Saldie")|
                (grepl("B",as.character(first(Stimulus))) &
                 first(StimLabel) == "Gatoo")){
               "A_Saldie"
             }else{"A_Gatoo"}},
           NBlocks = max(Block),
           LogNBlocks = log(NBlocks),
           RT = ifelse(RT < 200, NA, RT),
           LogRT = log(RT),
           zLogRT = scale(LogRT)) %>%
    ungroup() %>%
    mutate(FstLst = case_when(Block == 1 ~ "First Block",
                              Block == NBlocks ~ "Last Block")) %>%
    mutate_at(c("Participant", "Phase", "Condition", "CategoryName", "FstLst",
                "CRESP","RESP","ACC","CurrentObject","Stimulus","StimLabel", "Gender"),
              parse_factor, levels = NULL, include_na = F) %>%
    select(-one_of("TimestampMicrosec","TimestampSec"))
  if(participants == "adults_3f"){
    df <- df %>%
      mutate(Diagnostic = case_when(grepl("[12]", Stimulus) ~ "Feet",
                                    grepl("[34]", Stimulus) ~ "Both",
                                    grepl("[56]", Stimulus) ~ "Tail"))
  }
  return(df)
}

# LOOKING-TIME DATA IMPORT -- INFANTS
# Function importing looking time data from all infant participants,
# in the ../results/infants/data/infants.tsv file by default
LT_data.import.infants <- function(res.repo="../results/infants/data/", file.name="infants.tsv",
                                   participants="infants"){
  # Get participant information (Gender, DOB, DOT)
  participant_info <- read_csv(paste0(res.repo,"ParticipantInformation.csv")) %>%
    drop_na(DOB) %>%
    mutate(Age = as.numeric(difftime(DOT, DOB, units = "days")),
           DiffTo15mo = case_when(Age >= 430 & Age <= 470 ~ Age - 15*7*52/12)) %>%
    drop_na(DiffTo15mo) %>%
    # 15months * 7days/week * 52/12weeks/month = 15mo in days (455)
    # Lower and upper limits: 455 -15 and 455 + 15 days
    select(-c(DOB,DOT))
  # Get sequence information
  sequence_info <- read.csv("../scripts/infants/SequenceInfo.csv") %>%
    mutate(CategoryName = ifelse(Name.first.fam.stim == "NoLabel",
                                 "NL",
                                 ifelse(First.fam.stim == "A",
                                        paste0("A_",substr(Name.first.fam.stim,1,1)),
                                        ifelse(Name.first.fam.stim == "Saldie",
                                               "A_G",
                                               "A_S"))),
           Condition = ifelse(Name.first.fam.stim == "NoLabel",
                              "No Label",
                              "Label")) %>%
    select(c(PresentationSequence,CategoryName,Condition))
  # Read file, do all transformations
  df <- read_tsv(paste0(res.repo,file.name)) %>%
    select(Participant = ParticipantName, PresentationSequence, MediaName,
           TimeStamp = RecordingTimestamp, StudioEventIndex, ValidityLeft, ValidityRight,
           CursorX = "GazePointX (ADCSpx)", CursorY = "GazePointY (ADCSpx)") %>%
    subset(!(grepl("AG",.$MediaName) | .$MediaName == "")) %>%
    group_by(Participant, MediaName) %>%
    mutate(TrackLoss = ValidityLeft + ValidityRight == 8,
           TrialId = max(StudioEventIndex, na.rm = T)/2,
           TrialNum = TrialId - 1) %>% # Trial kept as a numeric, starting at 0 for lmer
    ungroup() %>%
    drop_na(MediaName, TrackLoss) %>%
    select(-c(ValidityLeft, ValidityRight, StudioEventIndex)) %>%
    unique() %>%
    inner_join(participant_info) %>%
    left_join(sequence_info) %>%
    mutate(Phase = case_when(grepl("Flip|Reg", MediaName) ~ "Familiarisation",
                             grepl("WL[GS]_", MediaName) ~ "Test - Word Learning",
                             grepl("[HRT]C_", MediaName) ~ "Test - Contrast"),
           FamPart = case_when(TrialId <= 8 ~ 0,
                               TrialId <= 16 ~ 1,
                               TrialId <= 24 ~ 2),
           LabelOnset = case_when(grepl("_NL1", MediaName) ~ 1940,
                                  grepl("_NL2", MediaName) ~ 2240,
                                  grepl("_G1", MediaName) ~ 2300,
                                  grepl("_G2", MediaName) ~ 2485,
                                  grepl("_S1", MediaName) ~ 2050,
                                  grepl("_S2", MediaName) ~ 2400),
           TrialEnd = LabelOnset + 4000,
           Stimulus = ifelse(Phase == "Familiarisation",
                             sapply(strsplit(as.character(MediaName), "_"), "[", 2),
                             sapply(strsplit(as.character(MediaName), "_"), "[", 1)),
           AOI_type = case_when(grepl("Flip", MediaName) ~ "Flip",
                                grepl("Reg", MediaName) ~ "Reg",
                                grepl("HC_[AB][12]L", MediaName) ~ "NewHeadR",
                                grepl("HC_[AB][12]R", MediaName) ~ "NewHeadL",
                                grepl("TC_[AB][12]L", MediaName) ~ "NewTailR",
                                grepl("TC_[AB][12]R", MediaName) ~ "NewTailL",
                                grepl("RC_[AB]L", MediaName) ~ "NewHeadL_NewTailR",
                                grepl("RC_[AB]R", MediaName) ~ "NewHeadR_NewTailL",
                                grepl("WL[SG]_A1", MediaName) ~ 
                                  ifelse(CategoryName == "NL",
                                         "HeadsIn_NoTarget",
                                         paste0("HeadsIn_Target",
                                                ifelse(CategoryName==paste0("A_",
                                                                            substr(MediaName,3,3)),
                                                       substr(MediaName,7,7),
                                                       substr(MediaName,11,11)))),
                                grepl("WL[SG]_A2", MediaName) ~
                                  ifelse(CategoryName == "NL",
                                         "HeadsOut_NoTarget",
                                         paste0("HeadsOut_Target",
                                                ifelse(CategoryName==paste0("A_",
                                                                            substr(MediaName,3,3)),
                                                       substr(MediaName,7,7),
                                                       substr(MediaName,11,11))))))
  return(df)
}

# LOOKING-TIME DATA TO BEHAVIOUR
# Function extracting all non-LT data per participant per trial
LT_data.to_behaviour <- function(df){
  df %<>% select(-c(CursorX, CursorY,
                    CurrentObject, TimeStamp,
                    TrackLoss, AOI_type)) %>%
    unique()
  return(df)
}

# LOOKING-TIME DATA TO EYETRACKINGR
# Function adding AOIs, defining trial time-windows, and returning eyetrackingR data
LT_data.to_eyetrackingR <- function(df, participants, AOIs){
  # TODO/GOAL - Get list of AOIs (one df per AOI), all defined with AOI_type specific values
  # Transform dfname into string with deparse(substitute(dfname))
  # Add AOIs to data frame, one by one
  df <- mutate(df, NonAOI = !TrackLoss)
  for (AOI in names(AOIs)){
    AOI.name <- sub("infants\\.|adults_[23]f\\.", "", AOI)
    df %<>% left_join(AOIs[[AOI]]) %>%
      mutate(!!AOI.name := CursorX>Left & CursorX<Right & CursorY>Top & CursorY<Bottom,
             NonAOI = xor(NonAOI, CursorX>Left & CursorX<Right & CursorY>Top & CursorY<Bottom)) %>%
      select(-one_of("Left", "Right", "Top", "Bottom"))
  }
  # Set starting time of all trials to 0
  df %<>% group_by(Participant, TrialId) %>%
    mutate(TimeStamp = TimeStamp - min(TimeStamp))
  # Additional useful variable depending on participants
  if(grepl("adults_[23]f", participants)){
    df %<>% subset(CurrentObject == "Stimulus") %>%
      group_by(Participant, TrialId) %>%
      summarise(FeedbackOnset = max(TimeStamp)) %>%
      inner_join(df, .) %>%
      mutate(TimeStamp = TimeStamp - FeedbackOnset)
  }
  return(df)
}

# LOOKING-TIME DATA TRACKLOSS CLEAN
# Cleans the data by trackloss with specified thresholds, saving diagnostic plots
LT_data.trackloss_clean <- function(df, participants="adults_2f", trial_prop_thresh=.3,
                                    incl_crit=.7, verbose = F, graphs = F){
  # Get trackloss information for plotting (inlc. test)
  trackloss.subject.trial <- trackloss_analysis(df)
  if(graphs){
    # Plot trackloss per trial per subject
    trackloss.subject.trial.p <- ggplot(trackloss.subject.trial,
                                        aes(x=TrialId, y=TracklossForTrial)) +
      facet_wrap(~Participant, nrow = 10, scales = "free_x") + geom_point()
    ggsave(paste0("../results/", participants, "/cleaning/TracklossSubjectTrial.png"),
           trackloss.subject.trial.p, width=9, height=15)
  }
  # Remove trials with trackloss proportion greater than 0.25
  df.trackloss <- clean_by_trackloss(data = df,
                                     trial_prop_thresh = trial_prop_thresh)
  # Compute and plot proportion of valid trials per subject
  # (number of valid trials / number of trials for subject)
  # (looking only at familiarisation phase)
  df.trackloss$TrialId <- as.numeric(df.trackloss$TrialId)
  df.described <- describe_data(df.trackloss[which(df.trackloss$Phase == "Familiarisation"),],
                                'TrialId', 'Participant')
  if(grepl("adults_[23]f", participants)){
    df.described$ProportionTrials <- df.described$NumTrials / df.described$Max
  }else{
    df.described$ProportionTrials <- df.described$NumTrials / 24
  }
  df.described$AboveCriteria <- factor(df.described$ProportionTrials >= incl_crit)
  if(verbose){
    print(summary(df.described))
  }
  if(graphs){
    df.described.p <- ggplot(df.described,
                             aes(x=Participant, y=ProportionTrials, colour = AboveCriteria)) +
      scale_colour_manual(values = c("red","green"), guide = F) + geom_point()
    ggsave(paste0("../results/", participants, "/cleaning/ProportionTrialsPerSubject.png"),
           plot = df.described.p, width = 10, height = 3)
  }
  # Select subjects to keep
  df.trackloss <- inner_join(df.trackloss, select(df.described,
                                                  one_of("Participant",
                                                         "AboveCriteria")))
  df.clean <- df.trackloss %>%
    subset(AboveCriteria == T, select = -AboveCriteria)
  # Print how many subjects remain per condition
  print(df.clean %>% group_by(Condition) %>% summarise(n_distinct(Participant)))
  return(df.clean)
}

# LOOKING-TIME DATA GATHER
# General function for importing, transforming into eyetrackingR, and cleaning data,
# as well as giving general plots for general information on the data
LT_data.gather <- function(participants, verbose = F, graphs = F){
  # Define a list of AOI dataframes for all experiments
  AOIs <<- list()
  AOIs[["adults_2f.Head"]] <<- tibble(AOI_type=1,
                                          Left=c(400), Right=c(620),
                                          Top=c(55), Bottom=c(255))
  AOIs[["adults_2f.Tail"]] <<- tibble(AOI_type=1,
                                          Left=c(20), Right=c(220),
                                          Top=c(110), Bottom=c(330))
  AOIs[["adults_3f.Head"]] <<- tibble(AOI_type=1,
                                          Left=c(363), Right=c(602),
                                          Top=c(66), Bottom=c(295))
  AOIs[["adults_3f.Feet"]] <<- tibble(AOI_type=1,
                                          Left=c(146), Right=c(483),
                                          Top=c(364), Bottom=c(465))
  AOIs[["adults_3f.Tail"]] <<- tibble(AOI_type=1,
                                          Left=c(36), Right=c(260),
                                          Top=c(115), Bottom=c(310))
  AOIs[["infants.Head"]] <<- tibble(AOI_type=c("Reg","Flip"),
                                        Left=c(1031,1920-1031-450),
                                        Right=c(1031+450,1920-1031),
                                        Top=c(197,197),
                                        Bottom=c(197+450,197+450))
  AOIs[["infants.Tail"]] <<- tibble(AOI_type=c("Reg","Flip"),
                                        Left=c(390,1920-390-450),
                                        Right=c(390+450,1920-390),
                                        Top=c(299,299),
                                        Bottom=c(299+450,299+450))
  AOIs[["infants.NewTail"]] <<- tibble(AOI_type=c("NewTailR", "NewTailL",
                                                  "NewHeadL_NewTailR",
                                                  "NewHeadR_NewTailL"),
                                       Left=c(1920-450-1, 1, 1920-450-1, 1),
                                       Right=c(1920-1, 1+450, 1920-1, 1+450),
                                       Top=rep(299, 4),
                                       Bottom=rep(299+450, 4))
  AOIs[["infants.OldTail"]] <<- tibble(AOI_type=c("NewTailR", "NewTailL"),
                                       Left=c(1, 1920-450-1),
                                       Right=c(1+450, 1920-1),
                                       Top=rep(299, 2),
                                       Bottom=c(299+450, 2))
  AOIs[["infants.NewHead"]] <<- tibble(AOI_type=c("NewHeadR", "NewHeadL",
                                                  "NewHeadL_NewTailR",
                                                  "NewHeadR_NewTailL"),
                                       Left=c(1920-450-1, 1, 1920-450-1, 1),
                                       Right=c(1920-1, 1+450, 1920-1, 1+450),
                                       Top=rep(197, 4),
                                       Bottom=rep(197+450, 4))
  AOIs[["infants.OldHead"]] <<- tibble(AOI_type=c("NewHeadR", "NewHeadL"),
                                       Left=c(1, 1920-450-1),
                                       Right=c(1+450, 1920-1),
                                       Top=rep(197, 2),
                                       Bottom=rep(197+450, 2))
  AOIs[["infants.Centre"]] <<- tibble(AOI_type=c("NewTailR", "NewTailL",
                                                 "NewHeadR", "NewHeadL",
                                                 "NewHeadL_NewTailR",
                                                 "NewHeadR_NewTailL"),
                                      Left=rep(510, 6),
                                      Right=rep(510+900, 6),
                                      Top=rep(275, 6),
                                      Bottom=rep(275+450, 6))
  AOIs[["infants.Target"]] <<- tibble(AOI_type=c("HeadsIn_TargetL", "HeadsOut_TargetL",
                                                 "HeadsIn_TargetR", "HeadsOut_TargetR"),
                                      Left=c(1, 1, 1920/2 + 10, 1920/2 + 10),
                                      Right=c(1920/2 - 10, 1920/2 - 10, 1920-1, 1920-1),
                                      Top=rep(250, 4),
                                      Bottom=rep(250+660, 4))
  AOIs[["infants.Distractor"]] <<- tibble(AOI_type=c("HeadsIn_TargetL", "HeadsOut_TargetL",
                                                     "HeadsIn_TargetR", "HeadsOut_TargetR"),
                                          Left=c(1920/2 + 10, 1920/2 + 10, 1, 1),
                                          Right=c(1920/2 - 1, 1920/2 - 1, 1920-10, 1920-10),
                                          Top=rep(250, 4),
                                          Bottom=rep(250+660, 4))
  # Import raw data
  fun_name <- paste0("LT_data.import.",
                     sub("_[23]f", "", participants))
  raw_data <- match.fun(fun_name)(participants = participants)
  # Extract behavioural data for adults
  if(grepl("adults_[23]f", participants)){
    behaviour <- LT_data.to_behaviour(raw_data)
  }
  else{
    behaviour <- NULL
  }
  # Create raw (not clean) eyetrackingR data
  current.AOIs <- grep(participants, names(AOIs)) # First get the indexes of AOIs to use
  LT.raw_data <-raw_data %>%
    LT_data.to_eyetrackingR(participants, AOIs[current.AOIs]) %>%
    make_eyetrackingr_data(participant_column = "Participant",
                           trial_column = "TrialId",
                           time_column = "TimeStamp",
                           trackloss_column = "TrackLoss",
                           aoi_columns = sub(paste0(participants,"\\."), "",
                                             names(AOIs)[current.AOIs]),
                           treat_non_aoi_looks_as_missing = F)
  # Analyse trackloss, clean data
  LT.AOI_summary <- LT.raw_data %>%
  {if(grepl("adults_[23]f", participants)){
    group_by(., Participant, CurrentObject)
  }else{group_by(., Participant, Condition)}}%>%
    summarise(TrackLossRatio = sum(TrackLoss)/n(),
              NonAOIRatio = sum(NonAOI, na.rm = T)/(n()-sum(TrackLoss)))
  # Save diagnostic plot if necessary
  if(graphs){
    # Select aesthetics for plot depending on available variables
    if(grepl("adults_[23]f", participants)){
      LT.AOI_summary.plot.TrackLossRatio <- ggplot(LT.AOI_summary,
                                                   aes(x = CurrentObject, y = TrackLossRatio,
                                                       fill = CurrentObject))
      LT.AOI_summary.plot.NonAOIRatio <- ggplot(LT.AOI_summary,
                                                aes(x = CurrentObject, y = NonAOIRatio,
                                                    fill = CurrentObject))
    }else{
      LT.AOI_summary.plot.TrackLossRatio <- ggplot(LT.AOI_summary,
                                                   aes(x = Condition, y = TrackLossRatio,
                                                       fill = Condition))
      LT.AOI_summary.plot.NonAOIRatio <- ggplot(LT.AOI_summary,
                                                aes(x = Condition, y = NonAOIRatio,
                                                    fill = Condition))
    }
    # Finishing plots with layers
    LT.AOI_summary.plot.TrackLossRatio <- LT.AOI_summary.plot.TrackLossRatio +
      geom_violin() +
      geom_boxplot(alpha=0, width=.2, outlier.alpha = 1) +
      guides(fill = "none")
    LT.AOI_summary.plot.NonAOIRatio <- LT.AOI_summary.plot.NonAOIRatio +
      geom_violin() +
      geom_boxplot(alpha=0, width=.2, outlier.alpha = 1) +
      guides(fill = "none")
    # Saving plots
    ggsave(paste0("../results/", participants, "/cleaning/TrackLossRatio.png"),
           plot = LT.AOI_summary.plot.TrackLossRatio)
    ggsave(paste0("../results/", participants, "/cleaning/NonAOIRatio.png"),
           plot = LT.AOI_summary.plot.NonAOIRatio)
  }
  # Make clean, with inclusion criteria dependent on participants tested
  if(grepl("adults_[23]f", participants)){
    LT.clean <- LT_data.trackloss_clean(LT.raw_data, participants,
                                        trial_prop_thresh = .3, incl_crit = .7,
                                        verbose = verbose, graphs = graphs)
  }else{
    LT.clean <- LT_data.trackloss_clean(LT.raw_data, participants,
                                        trial_prop_thresh = .5, incl_crit = .5,
                                        verbose = verbose, graphs = graphs)
  }
  # Plot heatmaps of remaining participants and trials by condition by phase
  if(graphs){
    x_max = min(max(LT.clean$CursorX, na.rm = T), 1920)
    y_max = min(max(LT.clean$CursorY, na.rm = T), 1080)
    LT.heatmap <- ggplot(LT.clean, aes(x=CursorX,y=CursorY)) +
      xlim(c(0, x_max)) + scale_y_reverse(limits = c(y_max, 0)) +
      theme(aspect.ratio = y_max/x_max) +
      facet_grid(Phase~Condition) +
      geom_bin2d(binwidth = c(20,20))
    ggsave(paste0("../results/", participants, "/cleaning/Heatmaps.png"),
           plot = LT.heatmap,
           width = 8)
  }
  # Return all datasets for analysis and checks
  return(list(raw_data, behaviour, LT.raw_data, LT.clean))
}
