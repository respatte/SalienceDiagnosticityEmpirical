library(doSNOW)
library(dplyr)

# LOOKING-TIME DATA IMPORT
# Function importing looking time data from all participants, in the ../results/ repository by default
LT_data.import <- function(res.repo="../results/adults/"){
  single.file.import <- function(file){
    tmp <- read.delim(file)[,-c(2:5,10:23)]
    return(droplevels(tmp[tmp$CurrentObject %in% c("Feedback","Label","Stimulus"),]))
  }
  cl <- makeCluster(4)
  registerDoSNOW(cl)
  file.names <- list.files(path=res.repo, pattern=".gazedata")
  df <- foreach(i=1:60,.combine="rbind", .inorder=F) %dopar% single.file.import(paste0(res.repo,file.names[i]))
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
  df$TrackLoss <- ifelse(pmin.int(df$CursorX,df$CursorY)<0,T,F)
  # Creating TimeStamp in milliseconds
  df$TimeStamp <- df$TimestampMicrosec*1e-3 + df$TimestampSec*1e3
  df <- df[,-(4:5)]
  # Adding participant information
  participant_info <- read.csv(paste0(res.repo,"ParticipantInformation.csv"))
  df <- merge(df, participant_info, by="Subject")
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
           LogRT = log(RT))
  return(df)
}

# LOOKING-TIME DATA TO EYETRACKINGR
# Function adding AOIs, defining trial time-windows, and returning eyetrackingR data
LT_data.to_eyetrackingR <- function(df, AOIs){
  # Add AOIs to data frame, one by one
  for (AOI in levels(AOIs$name)){
    row <- AOIs[AOIs$name==AOI,]
    df[,AOI] <- df$CursorX>row$L & df$CursorX<row$R & df$CursorY>row$T & df$CursorY<row$B
  }
  # Set starting time of all trials to 0
  df <- df %>% group_by(Subject, TrialId) %>% mutate(TimeStamp = TimeStamp - min(TimeStamp),
                                                     NormTimeStamp = TimeStamp/max(TimeStamp))
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