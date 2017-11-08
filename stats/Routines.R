library(doSNOW)
library(eyetrackingR)

# LOOKING-TIME DATA IMPORT
# Function importing looking time data from all participants, in the ../results/ repository by default
LT.data.import <- function(res.repo="../results/"){
  single.file.import <- function(file){
    tmp <- read.delim(file)[,-c(2:5,10:23)]
    tmp$TrialId <- ifelse(tmp$Block==0, tmp$TrialId + 252, tmp$TrialId)
    return(droplevels(tmp[tmp$CurrentObject %in% c("Feedback","Label","Stimulus"),]))
  }
  cl <- makeCluster(12)
  registerDoSNOW(cl)
  file.names <- list.files(path=res.repo, pattern=".gazedata")
  df <- foreach(i=1:60,.combine="rbind", .inorder=F) %dopar% single.file.import(paste0(res.repo,file.names[i]))
  stopCluster(cl)
  df$track_loss <- ifelse(pmin.int(df$CursorX,df$CursorY)<0,T,F)
  df$timestamp <- df$TimestampMicrosec + df$TimestampSec*1e6
  df <- df[,-(4:5)]
  return(df)
}

# RAW TO EYE-TRACKING
# Function adding AOIs, defining trial time-windows, and returning eyetrackingR data
raw.to.ET <- function(df){
  # Define AOIs
  AOIs <- data.frame(name=c("Tail","Head"),L=c(20,400),R=c(220,620),T=c(110,55),B=c(330,255))
  # Add AOIs to data frame, one by one
  for (AOI in levels(AOIs$name)){
    row <- AOIs[AOIs$name==AOI,]
    df[,AOI] <- df$CursorX>row$L & df$CursorX<row$R & df$CursorY>row$T & df$CursorY<row$B
  }
  # Set starting time of all trials to 0
  df <- df %>% group_by(Subject, TrialId) %>% mutate(modified_timestamp = timestamp - min(timestamp))
  return(df)
}