if (file.exists(paste0(FILES_LOC, "summaries_annotated/", FOLDER,"_fromJSON_ANNOTATED.csv"))) {
  MAIN <- read.csv(paste0(FILES_LOC, "summaries_annotated/", FOLDER,"_fromJSON_ANNOTATED.csv"))
} else if (file.exists(paste0(FILES_LOC, "summaries_annotated/", FOLDER,"fromSLURM_ANNOTATED.csv"))) {
  MAIN <- read.csv(paste0(FILES_LOC, "summaries_annotated/", FOLDER,"fromSLURM_ANNOTATED.csv"))
} else print("Please make sure correct path was provided or change to ANNOTATE = F.")