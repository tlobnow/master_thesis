#         .
#       ":"
#     ___:____     |"\/"|
#   ,'        `.    \  /
#   |  O        \___/  |
# ~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^

load_dtConnOrLocal <- function(PATH_DT, PATH_LOCAL, FILE) {
  if ( DATA_TAY_CONNECTION == F ) { # source from local folder
    read_xlsx(path = paste0(PATH_LOCAL, FILE), sheet = 1, na = c("", " ", "NA", "NaN"))
  } else { # source from Data-Tay
    read_xlsx(path = paste0(PATH_DT, FILE), sheet = 1, na = c("", " ", "NA", "NaN"))
  }}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

retrieveAccessionIDs <- function(DF,OUT="retrievedAccessionIDs.txt") {
  # Select UniProt Accession Numbers and retain unique values
  DF_uniq <- DF %>% select(Protein.IDs) %>% unique()
  
  # filter out contaminants
  DF_fil <- DF_uniq %>% filter(!str_detect(Protein.IDs, paste("CON__")))
  
  # separate the joined protein IDs (sometimes multiple per row, separated by ";")
  DF_sep <- unlist(lapply(strsplit(DF_fil$Protein.IDs, ";", fixed=TRUE), function(x) return(x[1:50]))) %>%
    unique()
  
  # filter out NAs
  DF_sep <- as.data.frame(DF_sep) %>% drop_na()
  
  # filter out invalid Accession IDs ("REV__")
  DF_sep <- DF_sep %>% filter(!str_detect(DF_sep, paste("REV__")))
  
  # write text file
  write.table(x = DF_sep, 
              file = OUT, 
              quote = F, sep = "\t", row.names = F, col.names = F)
}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

jsonExtract_1.1_ratio <- function(JSON, OUT) {
  # REMOVE ALL "INFINITY" STRINGS IN YOUR JSON FILES USING:
  # grep -rl Infinity . | xargs sed -i 's/Infinity/9999/g'
  json <- fromJSON(JSON)
  # EXTRACT FILE
  FILE_A  <- unlist(lapply(strsplit(json[["chains"]][["A"]], "_", fixed=TRUE), function(x) return(x[1])))
  FILE_B  <- unlist(lapply(strsplit(json[["chains"]][["B"]], "_", fixed=TRUE), function(x) return(x[1])))
  FILE    <- paste(FILE_A, FILE_B, "x1", sep = "_")
  MODEL   <- unlist(lapply(strsplit(as.data.frame(json$order)[1,], "_", fixed=TRUE), function(x) return(x[2])))
  TOL     <- as.data.frame(json$tol_values) %>% pivot_longer(names_to = "RECYCLE", values_to = "TOL", cols = 1:ncol(.))
  pLDDT   <- as.data.frame(json$plddts) %>% pivot_longer(names_to = "RECYCLE", values_to = "pLDDT", cols = 1:ncol(.))
  pTM     <- as.data.frame(json$ptms) %>% pivot_longer(names_to = "RECYCLE", values_to = "pTM", cols = 1:ncol(.))
  piTM    <- as.data.frame(json$pitms) %>% pivot_longer(names_to = "RECYCLE", values_to = "piTM", cols = 1:ncol(.))
  iScore  <- as.data.frame(json$`interface score`) %>% pivot_longer(names_to = "RECYCLE", values_to = "iScore", cols = 1:ncol(.))
  iRes    <- as.data.frame(json$`interfacial residue number`) %>% pivot_longer(names_to = "RECYCLE", values_to = "iRes", cols = 1:ncol(.))
  iCnt    <- as.data.frame(json$`interficial contact number`) %>% pivot_longer(names_to = "RECYCLE", values_to = "iCnt", cols = 1:ncol(.))
  FILE_MODEL    <- paste(FILE, MODEL, sep = "_")
  NUM_CLUSTERS  <- json[["clusters"]][[iScore$RECYCLE[1]]][["num_clusters"]]
  N_MONOMERS    <- length(json[["chains"]])
  # JOIN, REMOVE DUPLICATES, WRITE TO CSV
  EXTRACT   <- cbind(FILE, MODEL, TOL, pLDDT, pTM, piTM, iScore, iRes, iCnt, FILE_MODEL, NUM_CLUSTERS, N_MONOMERS)
  EXTRACT   <- EXTRACT[, !duplicated(colnames(EXTRACT))]
  write.table(EXTRACT, file = paste0(OUT,"_fromJSON.csv"),sep = ",", append = T, quote = F, row.names = F, col.names = F)
}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

slurmExtract <- function(SLURM, OUT) {
  SLURM        <- fread(SLURM, sep = "\t", header = F)
  SLURM        <- SLURM %>% mutate(RECYCLED = F, TOP = F, INFO = F, GRAB = F, N_MON = F, MODEL = F)
  SLURM$RECYCLED[grep("*_recycled_*",      SLURM$V1)] <- T
  SLURM$INFO[grep("*Info:*",               SLURM$V1)] <- T
  SLURM$GRAB[grep("Info: num_clusters*",   SLURM$V1)] <- T
  SLURM$N_MON[grep("* to model *",         SLURM$V1)] <- T
  SLURM$MODEL[grep("*pLDDT*",              SLURM$V1)] <- T
  
  MODEL_INFO        <- SLURM %>% filter(MODEL == T & RECYCLED == F) %>% select(V1)
  MODEL_INFO        <- separate(data = MODEL_INFO, col = V1,   sep = ",",   into = c("NAME", "TOL", "pLDDT", "pTM", "piTM", "iScore", "iRes"), convert = T)
  MODEL_INFO        <- separate(data = MODEL_INFO, col = NAME, sep = " ",   into = c("INFO", "FILE", "MODEL", "PERFORMED", "X", "CYCLE"), convert = T)
  MODEL_INFO$TOL    <- as.numeric(unlist(lapply(strsplit(MODEL_INFO$TOL,   "= ", fixed=TRUE), function(x) return(x[2]))))
  MODEL_INFO$pLDDT  <- as.numeric(unlist(lapply(strsplit(MODEL_INFO$pLDDT, "= ", fixed=TRUE), function(x) return(x[2]))))
  MODEL_INFO$pTM    <- as.numeric(unlist(lapply(strsplit(MODEL_INFO$pTM,   "= ", fixed=TRUE), function(x) return(x[2]))))
  MODEL_INFO$piTM   <- as.numeric(unlist(lapply(strsplit(MODEL_INFO$piTM,  "= ", fixed=TRUE), function(x) return(x[2]))))
  MODEL_INFO$iScore <- as.numeric(unlist(lapply(strsplit(MODEL_INFO$iScore,"= ", fixed=TRUE), function(x) return(x[2]))))
  MODEL_INFO        <- separate(data = MODEL_INFO, col = iRes, sep =  "iCnt = ", into = c("iRes", "iCnt"), convert = T)
  MODEL_INFO$iRes   <- as.numeric(unlist(lapply(strsplit(MODEL_INFO$iRes,  "= ", fixed=TRUE), function(x) return(x[2]))))
  MODEL_INFO$iCnt   <- as.numeric(MODEL_INFO$iCnt)
  MODEL_INFO$Clash_Indicator <- MODEL_INFO$iRes / MODEL_INFO$iCnt
  
  CLUSTER_INFO               <- SLURM %>% filter(GRAB == T) %>% select(V1)
  CLUSTER_INFO               <- separate(data = CLUSTER_INFO, col = V1,   sep = " = ", into = c("V1", "NUM_CLUSTERS", "CLUSTER_SIZES", "CLUSTERS"), convert = T, remove = T)
  CLUSTER_INFO$NUM_CLUSTERS  <- as.numeric(unlist(lapply(strsplit(CLUSTER_INFO$NUM_CLUSTERS, ", ", fixed=TRUE), function(x) return(x[1]))))
  CLUSTER_INFO$CLUSTER_SIZES <- unlist(lapply(strsplit(CLUSTER_INFO$CLUSTER_SIZES, ",  clusters", fixed=TRUE), function(x) return(x[1])))
  CLUSTER_INFO$CLUSTER_SIZES <- str_replace(CLUSTER_INFO$CLUSTER_SIZES, ", ", "_")
  CLUSTER_INFO$CLUSTERS      <- str_replace(CLUSTER_INFO$CLUSTERS, ", ", "/")
  
  N_MONOMERS_INFO <- SLURM %>% filter(N_MON == T) %>% select(V1)
  N_MONOMERS_INFO <- separate(data = N_MONOMERS_INFO, col = V1,   sep = c("chain"), into = c("N_MONOMERS", "TRASH"), convert = T, remove = T) %>% unique()
  N_MONOMERS_INFO$N_MONOMERS <- unlist(lapply(strsplit(N_MONOMERS_INFO$N_MONOMERS, ": ", fixed=TRUE), function(x) return(x[2]))) %>% as.numeric()
  N_MONOMERS_INFO <- N_MONOMERS_INFO %>% select(N_MONOMERS)
  
  EXTRACT       <- bind_cols(MODEL_INFO, CLUSTER_INFO)
  EXTRACT       <- EXTRACT %>% 
    mutate(FILE_MODEL = paste(FILE, MODEL, sep = "_"), N_MONOMERS = N_MONOMERS_INFO$N_MONOMERS) %>% 
    select(FILE, MODEL, TOL, pLDDT, pTM, piTM, iScore, iRes, iCnt, FILE_MODEL, NUM_CLUSTERS, N_MONOMERS)
  EXTRACT$MODEL <- unlist(lapply(strsplit(EXTRACT$MODEL, "_multimer", fixed=TRUE), function(x) return(x[1])))
  EXTRACT       <- unique(EXTRACT)
  
  write.table(EXTRACT, file = paste0(OUT ,"_fromSLURM.csv"),sep = ",", append = T, quote = F, row.names = F, col.names = F)
}


join_timepoints <- function(DF1, DF2, DF3, OUT) {
  # Join the initial dataframes derived from protein X pull down
  DF <- full_join(DF1, DF2)
  DF <- full_join(DF, DF3)
}
