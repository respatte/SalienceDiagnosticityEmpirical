# LT DATA IMPORT
# Function importing looking time data from all participants, in the ../results/ repository by default.
# This function returns a dataframe with the raw gazedata from E-Prime.
LT.data.import <- function(res.repo="../results/"){
  df = data.frame()
  for (file.name in list.files(path=res.repo, pattern=".gazedata")){
    df <- rbind(df, read.delim(paste0(res.repo,file.name)))
  }
  return(df)
}