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
  
  write.table(EXTRACT, file = paste0(OUT ,".csv"),sep = ",", append = T, quote = F, row.names = F, col.names = F)
}