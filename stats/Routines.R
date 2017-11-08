library(eyetrackingR)

# LOOKING-TIME DATA IMPORT
# Function importing looking time data from all participants, in the ../results/ repository by default
LT.data.import <- function(res.repo="../results/"){
  df = data.frame()
  i <- 1
  for (file.name in list.files(path=res.repo, pattern=".gazedata")){
    print(i)
    i <- i+1
    tmp <- read.delim(paste0(res.repo,file.name))[,-c(2:5,10:23)]
    tmp$TrialId <- ifelse(tmp$Block==0, tmp$TrialId + 252, tmp$TrialId)
    df <- rbind(df, droplevels(tmp[tmp$CurrentObject %in% c("Feedback","Label","Stimulus"),]))
  }
  df <- df
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