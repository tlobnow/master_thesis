MAIN <- fread(paste0(FILES_LOC, "summaries/", FOLDER,"_fromJSON.csv"))
if (file.exists(paste0(FILES_LOC, "summaries/", FOLDER,"_fromJSON.csv"))) {
  MAIN <- read.csv(paste0(FILES_LOC, "summaries/", FOLDER,"_fromJSON.csv"))
  MAIN$RECYCLE_num <-unlist(lapply(strsplit(MAIN$RECYCLE, "_", fixed=TRUE), function(x) return(x[2])))
} else if (file.exists(paste0(FILES_LOC, "summaries/", FOLDER,"fromSLURM_.csv"))) {
  MAIN <- read.csv(paste0(FILES_LOC, "summaries/", FOLDER,"fromSLURM_.csv"))
} else print("Please make sure correct path was provided or change to ANNOTATE = F.")