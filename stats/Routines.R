library(doSNOW)
library(tidyr)
library(dplyr)
library(reshape2)

# LOOKING-TIME DATA IMPORT -- ADUTLTS
# Function importing looking time data from all adult participants, in the ../results/adults repository by default
LT_data.adults.import <- function(res.repo="../results/adults/", subjects=1:60){
  single.file.import <- function(file){
    tmp <- read.delim(file)[,-c(2:5,10:23)]
    return(droplevels(tmp[tmp$CurrentObject %in% c("Feedback","Label","Stimulus"),]))
  }
  cl <- makeCluster(4)
  registerDoSNOW(cl)
  file.names <- list.files(path=res.repo, pattern=".gazedata")
  df <- foreach(i=subjects,.combine="rbind", .inorder=F) %dopar% single.file.import(paste0(res.repo,file.names[i]))
  stopCluster(cl)
  df$TrialId <- ifelse(df$Block==0, df$TrialId + 252, df$TrialId)
  df <- df %>% group_by(Subject) %>%
    mutate(Condition = factor(if(first(StiLabel) == "NoLabelFeedback"){"NoLabel"}else{"Label"},
                              levels = c("NoLabel","Label")),
           CategoryName = factor(if(first(Condition) == "NoLabel"){"NoName"}else{
             if((grepl("A",as.character(first(Stimulus))) &
                 first(StiLabel) == "Saldie")|
                (grepl("B",as.character(first(Stimulus))) &
                 first(StiLabel) == "Gatoo")){
               "A_Saldie"
             }else{"A_Gatoo"}},
             levels = c("NoName","A_Saldie","A_Gatoo")))
  df$TrackLoss <- pmin.int(df$CursorX,df$CursorY)<0
  # Creating TimeStamp in milliseconds
  df$TimeStamp <- df$TimestampMicrosec*1e-3 + df$TimestampSec*1e3
  df <- df[,-(4:5)]
  # Adding participant information
  participant_info <- read.csv(paste0(res.repo,"ParticipantInformation.csv"))
  df <- merge(df, participant_info, by="Subject")
  return(df)
}

# LOOKING-TIME DATA IMPORT -- ADUTLTS
# Function importing looking time data from all infant participants, in the ../results/infants.tsv file by default
LT_data.infants.import <- function(res.repo="../results/infants/", file.name="infants.tsv"){
  # Read file, drop empty MediaName rows, delete Attention Gatherer (AG) recordings,
  # as well as unused columns from Tobii output (update with new output)
  df <- read.csv(paste0(res.repo,file.name), sep = "\t") %>%
    drop_na(MediaName) %>%
    subset(!(grepl("AG",.$MediaName) | .$MediaName == "")) %>%
    droplevels() %>%
    group_by(ParticipantName, MediaName) %>%
    mutate(TrackLoss = ValidityLeft + ValidityRight == 8,
           TrialN = max(StudioEventIndex, na.rm = T)) %>%
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
           TrialId = ifelse(Phase == "Familiarisation",
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
LT_data.to_eyetrackingR <- function(df, AOIs, set.trial.start = T){
  # Add AOIs to data frame, one by one
  for (AOI in levels(AOIs$name)){
    row <- AOIs[AOIs$name==AOI,]
    df[,AOI] <- df$CursorX>row$L & df$CursorX<row$R & df$CursorY>row$T & df$CursorY<row$B
  }
  if(set.trial.start){
    # Set starting time of all trials to 0
    df <- df %>% group_by(Subject, TrialId) %>% mutate(TimeStamp = TimeStamp - min(TimeStamp),
                                                       NormTimeStamp = TimeStamp/max(TimeStamp))
  }
  return(df)
}

LT_data.trackloss_clean <- function(df, trial_prop_thresh=.25, incl_crit=.5, verbose=F){
  # Get trackloss information
  trackloss <- trackloss_analysis(df)
  trackloss.subject.trial <- unique(trackloss[, c('Subject','TrialId','TracklossForTrial')])
  # Plot trackloss per trial per subject
  trackloss.subject.trial.p <- ggplot(trackloss.subject.trial, aes(x=TrialId, y=TracklossForTrial)) +
    facet_wrap(~Subject, nrow = 10, scales = "free_x") + geom_point()
  ggsave("../results/TracklossSubjectTrial.png", trackloss.subject.trial.p, width=9, height=15)
  # Remove trials with trackloss proportion greater than 0.25
  df.trackloss <- clean_by_trackloss(data = df,
                                            trial_prop_thresh = trial_prop_thresh)
  # Compute and plot proportion of valid trials per subject (number of valid trials / number of trials for subject)
  df.described <- describe_data(df.trackloss, 'Block', 'Subject')
  df.described$ProportionTrials <- df.described$NumTrials /
    (12*(df.described$Max + 1))
  df.described$AboveCriteria <- factor(ifelse(df.described$ProportionTrials >= incl_crit, 1, 0))
  if(verbose){
    print(summary(df.described))
  }
  df.described.p <- ggplot(df.described,
                                  aes(x=Subject, y=ProportionTrials, colour = AboveCriteria)) +
    scale_colour_manual(values = c("red","green"), guide = F) + geom_point()
  ggsave("../results/ProportionTrialPerSubject.png", df.described.p, width=10, height=3)
  # Select subjects to keep
  df.clean <- df.trackloss[df.trackloss$Subject %in%
                                           df.described$Subject[df.described$AboveCriteria == 1],] %>%
    droplevels()
  if(verbose){
    # Check how many subjects missing per condition
    print(summary(unique(df.clean[,c('Subject','Condition')])))
  }
  return(df.clean)
}