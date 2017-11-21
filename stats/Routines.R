library(doSNOW)
library(dplyr)

# LOOKING-TIME DATA IMPORT
# Function importing looking time data from all participants, in the ../results/ repository by default
LT.data.import <- function(res.repo="../results/adults/"){
  single.file.import <- function(file){
    tmp <- read.delim(file)[,-c(2:5,10:23)]
    tmp$TrialId <- ifelse(tmp$Block==0, tmp$TrialId + 252, tmp$TrialId)
    tmp$Condition <- factor(ifelse("NoLabelFeedback" %in% tmp$StiLabel,"NoLabel","Label"))
    return(droplevels(tmp[tmp$CurrentObject %in% c("Feedback","Label","Stimulus"),]))
  }
  cl <- makeCluster(4)
  registerDoSNOW(cl)
  file.names <- list.files(path=res.repo, pattern=".gazedata")
  df <- foreach(i=1:60,.combine="rbind", .inorder=F) %dopar% single.file.import(paste0(res.repo,file.names[i]))
  stopCluster(cl)
  df$TrackLoss <- ifelse(pmin.int(df$CursorX,df$CursorY)<0,T,F)
  df$TimeStamp <- df$TimestampMicrosec + df$TimestampSec*1e6
  df <- df[,-(4:5)]
  return(df)
}

# RAW TO EYE-TRACKING
# Function adding AOIs, defining trial time-windows, and returning eyetrackingR data
raw.to.ET <- function(df, AOIs){
  # Add AOIs to data frame, one by one
  for (AOI in levels(AOIs$name)){
    row <- AOIs[AOIs$name==AOI,]
    df[,AOI] <- df$CursorX>row$L & df$CursorX<row$R & df$CursorY>row$T & df$CursorY<row$B
  }
  # Set starting time of all trials to 0
  df <- df %>% group_by(Subject, TrialId) %>% mutate(timestamp = timestamp - min(timestamp))
  return(df)
}
