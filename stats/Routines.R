library(doSNOW)
library(reshape2)
library(tidyverse)

# LOOKING-TIME DATA IMPORT -- ADUTLTS
# Function importing looking time data from all adult participants, in the ../results/adults repository by default
LT_data.adults.import <- function(participants="adults_2f", pinfo = T){
  single.file.import <- function(file){
    tmp <- read.delim(file)[,-c(2:5,10:23)]
    return(droplevels(tmp[tmp$CurrentObject %in% c("Feedback","Label","Stimulus"),]))
  }
  cl <- makeCluster(4)
  registerDoSNOW(cl)
  res.repo <- paste0("../results/",participants,"/data/")
  file.names <- list.files(path=res.repo, pattern=".gazedata")
  df <- foreach(i=1:length(file.names),.combine="rbind", .inorder=F) %dopar%
    single.file.import(paste0(res.repo,file.names[i]))
  stopCluster(cl)
  df <- df %>%
    rename(Participant = Subject) %>%
    mutate(TrialId = ifelse(Block==0, TrialId + 252, TrialId),
           TrackLoss = pmin.int(CursorX,CursorY)<0,
           TimeStamp = TimestampMicrosec*1e-3 + TimestampSec*1e3,
           AOI_type = 1) %>%
    group_by(Participant) %>%
    mutate(Condition = factor(if(first(StiLabel) == "NoLabelFeedback"){"NoLabel"}else{"Label"},
                              levels = c("NoLabel","Label")),
           CategoryName = factor(if(first(Condition) == "NoLabel"){"NoName"}else{
             if((grepl("A",as.character(first(Stimulus))) &
                 first(StiLabel) == "Saldie")|
                (grepl("B",as.character(first(Stimulus))) &
                 first(StiLabel) == "Gatoo")){
               "A_Saldie"
             }else{"A_Gatoo"}},
             levels = c("NoName","A_Saldie","A_Gatoo"))) %>%
    select(-one_of("TimestampMicrosec","TimestampSec"))
  if(pinfo){
    # Adding participant information
    participant_info <- read.csv(paste0(res.repo,"ParticipantInformation.csv"))
    df <- df %>% inner_join(participant_info)
  }
  return(df)
}

# LOOKING-TIME DATA IMPORT -- ADUTLTS
# Function importing looking time data from all infant participants, in the ../results/infants.tsv file by default
LT_data.infants.import <- function(res.repo="../results/infants/data/", file.name="infants.tsv"){
  # Read file, drop empty MediaName rows, delete Attention Gatherer (AG) recordings,
  # as well as unused columns from Tobii output (update with new output)
  df <- read.csv(paste0(res.repo,file.name), sep = "\t") %>%
    rename(TimeStamp = RecordingTimestamp,
           Participant = ParticipantName) %>%
    subset(!(grepl("AG",.$MediaName) | .$MediaName == "")) %>%
    droplevels() %>%
    group_by(ParticipantName, MediaName) %>%
    mutate(TrackLoss = ValidityLeft + ValidityRight == 8,
           TrialId = max(StudioEventIndex, na.rm = T)/2) %>%
    drop_na(MediaName, TrackLoss) %>%
    subset(!duplicated(TimeStamp)) %>%
    select(-one_of("X", "ValidityLeft", "ValidityRight", "StudioEventIndex"))
  # Add participant information (Gender, DOB, DOT)
  participant_info <- read.csv(paste0(res.repo,"ParticipantInformation.csv"),
                               colClasses = c("factor","factor","Date","Date")) %>%
    drop_na(DOB) %>%
    mutate(Age = as.numeric(difftime(DOT, DOB, units = "days")),
           DiffTo15mo = Age - 15*7*52/12) %>%
    select(-one_of("DOB","DOT"))
  # 15months * 7days/week * 52/12weeks/month = 15months in days
  df <- merge(df, participant_info, by = "ParticipantName")
  # Add sequence information
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
    select(one_of("PresentationSequence","CategoryName","Condition"))
  df <- merge(df, sequence_info, by = "PresentationSequence")
  # TODO - Add Reg-Flip and others for AOIs
  df <- df %>%
    mutate(Phase = case_when(grepl("Flip|Reg", MediaName) ~ "Familiarisation",
                             grepl("WL[GS]_", MediaName) ~ "Test - Word Learning",
                             grepl("[HRT]C_", MediaName) ~ "Test - Contrast"),
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
                                grepl("WL[SG]_A1", MediaName) ~ ifelse(CategoryName == "NL",
                                                                       "HeadsIn_NoTarget",
                                                                       paste0("HeadsIn_Target",
                                                                              ifelse(CategoryName == paste0("A_",
                                                                                                            substr(MediaName,3,3)),
                                                                                     substr(MediaName,7,7),
                                                                                     substr(MediaName,11,11)))),
                                grepl("WL[SG]_A2", MediaName) ~ ifelse(CategoryName == "NL",
                                                                       "HeadsOut_NoTarget",
                                                                       paste0("HeadsOut_Target",
                                                                              ifelse(CategoryName == paste0("A_",
                                                                                                            substr(MediaName,3,3)),
                                                                                     substr(MediaName,7,7),
                                                                                     substr(MediaName,11,11))))
           ))
  return(df)
}

# LOOKING-TIME DATA TO RESPONSES
# Function extracting all non-LT data per participant per trial
LT_data.to_responses <- function(df){
  df <- df[,-c(2,3,9,12,15,16)] %>%
    unique() %>%
    group_by(Subject) %>%
    mutate(NBlocks = max(Block),
           LogNBlocks = log(NBlocks),
           RT = ifelse(RT < 200, NA, RT),
           LogRT = log(RT))
  return(df)
}

# LOOKING-TIME DATA TO EYETRACKINGR
# Function adding AOIs, defining trial time-windows, and returning eyetrackingR data
LT_data.to_eyetrackingR <- function(df, AOIs){
  # TODO/GOAL - Get list of AOIs (one df per AOI),
  # all possibly defined with AOI_type specific values
  # Transform dfname into string with deparse(substitute(dfname))
  # Add AOIs to data frame, one by one
  for (AOI.name in names(AOIs)){
    df <- df %>%
      left_join(AOIs[[AOI.name]]) %>%
      mutate(!!sub("infants\\.|adults_[23]f\\.", "", AOI.name) := CursorX>Left & CursorX<Right & CursorY>Top & CursorY<Bottom) %>%
      select(-one_of("Left", "Right", "Top", "Bottom"))
  }
  # Set starting time of all trials to 0
  df <- df %>%
    group_by(Participant, TrialId) %>%
    mutate(TimeStamp = TimeStamp - min(TimeStamp))
  return(df)
}

LT_data.trackloss_clean <- function(df, participants="adults_2f", trial_prop_thresh=.25, incl_crit=.5, verbose=F){
  res.repo <- paste0("../results/", participants, "/cleaning/")
  # Get trackloss information
  trackloss.subject.trial <- trackloss_analysis(df)
  # Plot trackloss per trial per subject
  trackloss.subject.trial.p <- ggplot(trackloss.subject.trial, aes(x=TrialId, y=TracklossForTrial)) +
    facet_wrap(~Subject, nrow = 10, scales = "free_x") + geom_point()
  ggsave(paste0(res.repo,"TracklossSubjectTrial.png"), trackloss.subject.trial.p, width=9, height=15)
  # Remove trials with trackloss proportion greater than 0.25
  df.trackloss <- clean_by_trackloss(data = df,
                                     trial_prop_thresh = trial_prop_thresh)
  # Compute and plot proportion of valid trials per subject (number of valid trials / number of trials for subject)
  df.described <- describe_data(df.trackloss, 'Block', 'Subject')
  df.described$ProportionTrials <- df.described$NumTrials /
    (12*(df.described$Max + 1))
  df.described$AboveCriteria <- factor(df.described$ProportionTrials >= incl_crit)
  if(verbose){
    print(summary(df.described))
  }
  df.described.p <- ggplot(df.described,
                           aes(x=Subject, y=ProportionTrials, colour = AboveCriteria)) +
    scale_colour_manual(values = c("red","green"), guide = F) + geom_point()
  ggsave(paste0(res.repo,"ProportionTrialPerSubject.png"), df.described.p, width=10, height=3)
  # Select subjects to keep
  df.trackloss <- merge(df.trackloss, select(df.described, one_of("Subject", "AboveCriteria")), by = "Subject")
  df.clean <- subset(df.trackloss, AboveCriteria == T) %>% droplevels()
  if(verbose){
    # Check how many subjects missing per condition
    print(summary(unique(df.clean[,c('Subject','Condition')])))
  }
  return(df.clean)
}
