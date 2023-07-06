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


# jsonExtract <- function(JSON, OUT, FILE) {
#   # REMOVE ALL "INFINITY" STRINGS IN YOUR JSON FILES USING:
#   # grep -rl Infinity . | xargs sed -i 's/Infinity/9999/g'
#   json <- fromJSON(JSON)
#   # EXTRACT FILE
#   #FILE_A  <- unlist(lapply(strsplit(json[["chains"]][["A"]], "_", fixed=TRUE), function(x) return(x[1])))
#   #FILE_B  <- unlist(lapply(strsplit(json[["chains"]][["B"]], "_", fixed=TRUE), function(x) return(x[1])))
#   #FILE    <- paste(FILE_A, FILE_B, "x1", sep = "_")
#   FILE    <- FILE
#   MODEL   <- unlist(lapply(strsplit(as.data.frame(json$order)[1,], "_", fixed=TRUE), function(x) return(x[2])))
#   TOL     <- as.data.frame(json$tol_values) %>% pivot_longer(names_to = "RECYCLE", values_to = "TOL", cols = 1:ncol(.))
#   pLDDT   <- as.data.frame(json$plddts) %>% pivot_longer(names_to = "RECYCLE", values_to = "pLDDT", cols = 1:ncol(.))
#   pTM     <- as.data.frame(json$ptms) %>% pivot_longer(names_to = "RECYCLE", values_to = "pTM", cols = 1:ncol(.))
#   piTM    <- as.data.frame(json$pitms) %>% pivot_longer(names_to = "RECYCLE", values_to = "piTM", cols = 1:ncol(.))
#   iScore  <- as.data.frame(json$`interface score`) %>% pivot_longer(names_to = "RECYCLE", values_to = "iScore", cols = 1:ncol(.))
#   iRes    <- as.data.frame(json$`interfacial residue number`) %>% pivot_longer(names_to = "RECYCLE", values_to = "iRes", cols = 1:ncol(.))
#   iCnt    <- as.data.frame(json$`interficial contact number`) %>% pivot_longer(names_to = "RECYCLE", values_to = "iCnt", cols = 1:ncol(.))
#   FILE_MODEL    <- paste(FILE, MODEL, sep = "_")
#   NUM_CLUSTERS  <- json[["clusters"]][[iScore$RECYCLE[1]]][["num_clusters"]]
#   N_MONOMERS    <- length(json[["chains"]])
#   # JOIN, REMOVE DUPLICATES, WRITE TO CSV
#   EXTRACT   <- cbind(FILE, MODEL, TOL, pLDDT, pTM, piTM, iScore, iRes, iCnt, FILE_MODEL, NUM_CLUSTERS, N_MONOMERS)
#   EXTRACT   <- EXTRACT[, !duplicated(colnames(EXTRACT))]
#   write.table(EXTRACT, file = paste0(OUT,"_fromJSON.csv"),sep = ",", append = T, quote = F, row.names = F, col.names = F)
# }

jsonExtract <- function(JSON, OUT, FILE) {
  json       <- fromJSON(JSON)
  MODEL      <- strsplit(json$order[[1]], "_", fixed = TRUE)[[1]][2]
  TOL        <- json$tol_values %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "TOL") 
  pLDDT      <- json$plddts %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "pLDDT") 
  pTM        <- json$ptms %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "pTM") 
  piTM       <- json$pitms %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "piTM") 
  iScore     <- json$`interface score` %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "iScore") 
  iRes       <- json$`interfacial residue number` %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "iRes") 
  iCnt       <- json$`interficial contact number` %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "iCnt") 
  FILE_MODEL <- paste(FILE, MODEL, sep = "_")
  NUM_CLUSTERS <- json$clusters[[iScore$RECYCLE[1]]]$num_clusters
  N_MONOMERS <- length(json$chains)
  DATE       <- sub(".*_(\\d+)_.*", "\\1", names(json$clusters))
  
  EXTRACT <- cbind(FILE, MODEL, TOL, pLDDT, pTM, piTM, iScore, iRes, iCnt, FILE_MODEL, NUM_CLUSTERS, N_MONOMERS)
  EXTRACT$DATE = lubridate::as_date(DATE)
  EXTRACT <- EXTRACT[, !duplicated(colnames(EXTRACT))]
  
  write.table(EXTRACT, file = paste0(OUT,".csv"),sep = ",", append = T, quote = F, row.names = F, col.names = F)
  
  EXTRACT_noRecycle <- EXTRACT %>% filter(!str_detect(RECYCLE, "_recycled_"))
  write.table(EXTRACT_noRecycle, file = paste0(OUT,"_noRecycle.csv"),sep = ",", append = T, quote = F, row.names = F, col.names = F)
  
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

ELISA_Fx <- function(Input_Directory) {
  # Initialize as an empty data frame
  All_plates_data = data.frame()
  
  # Get the list of subdirectories matching the pattern "Plate_"
  subdirs <- list.files(Input_Directory, recursive = FALSE, full.names = TRUE, pattern = "^Plate_\\d+_\\d{8}$")
  
  if (length(subdirs) > 0) {
    print("Plates exist!")
    for (input_plate_dir in subdirs) {
      
      Input_plate <- input_plate_dir
      
      #Reading Plate Treatment 
      MEASUREMENTS <- fread(paste0(Input_plate, "/MEASUREMENTS.csv"), header = F)
      CELL_LINES   <- fread(paste0(Input_plate, "/CELL_LINES.csv"), header = F)
      CONDITIONS   <- fread(paste0(Input_plate, "/CONDITIONS.csv"), header = F)
      STIM_DAYS    <- fread(paste0(Input_plate, "/STIM_DAYS.csv"), header = F)
      
      
      #Converting tables into vector for to make a single table
      MEASUREMENTS <- as.vector(as.matrix(MEASUREMENTS))
      CELL_LINES   <- as.vector(as.matrix(CELL_LINES))
      CONDITIONS   <- as.vector(as.matrix(CONDITIONS))
      STIM_DAYS    <- as.vector(as.matrix(STIM_DAYS))
      
      #Creating Table containing all plate Information
      Plate <- NULL
      Plate$MEASUREMENTS <- MEASUREMENTS
      Plate$CELL_LINES <- CELL_LINES
      Plate$CONDITIONS <- CONDITIONS
      Plate$STIM_DAYS <- STIM_DAYS
      
      rm(MEASUREMENTS, CELL_LINES, CONDITIONS, STIM_DAYS)
      
      Plate <- Plate %>% as.data.table()
      
      #Removing Empty Wells
      Plate <- Plate %>% filter(CELL_LINES != "BLANK") %>% as.data.table()
      
      # Standard Curve ---------------------------------------------------------
      Plate_Standards <- Plate %>% 
        filter(CONDITIONS == "CALIBRATION") %>% 
        group_by(CELL_LINES) %>% 
        summarise(MEASUREMENTS_mean = mean(MEASUREMENTS)) %>%  #, 
        #MEASUREMENTS_median = median(MEASUREMENTS)) %>%  
        mutate(CELL_LINES = as.numeric(CELL_LINES)) %>% 
        arrange(CELL_LINES)
      
      Fit     <- lm(CELL_LINES ~ MEASUREMENTS_mean - 1, data = Plate_Standards) #linear model of the Standard curve. -1 omits the intercept
      R       <- summary(Fit)$r.squared
      Rsquare <- signif(R, digits = 4)
      
      print(paste0("IL2-Amount = slope*Intensity"))
      print(paste0("IL2-Amount = ", Fit$coefficients[1],"*Intensity"))
      
      Plate_Standards <- Plate_Standards %>% 
        mutate(Fit_Test = (Fit$coefficients[1]*MEASUREMENTS_mean))
      
      # Plotting Standard Curve
      p <- ggplot(data = Plate_Standards) +
        geom_point(aes(x = MEASUREMENTS_mean, y = CELL_LINES), size = 5) +
        geom_line(aes(x = MEASUREMENTS_mean, y = Fit_Test), linetype = "dashed") +
        annotate('text', x = 0.15, y = 700, label = paste0("R^2 = ", Rsquare), size = 10) +
        annotate('text', 
                 x = max(Plate_Standards$MEASUREMENTS_mean) - (0.25 * max(Plate_Standards$MEASUREMENTS_mean)),
                 y = 150, label = paste0("IL-Amount = \n", signif(Fit$coefficients[1], digits = 4), " * Intensity")) +
        labs(x = "Measured Values",
             y = "IL-Concentration (pg/mL)") +
        ggtitle(label = paste0(basename(Input_plate)),
                subtitle = paste0("R^2 = ", Rsquare, "\n IL-Amount = ", signif(Fit$coefficients[1], digits = 4), " * Intensity")) +
        theme_classic() +
        theme(axis.title = element_text(size = 30),
              axis.text = element_text(size = 20))
      
      # Saving the plot
      Save_Name <- paste0(basename(Input_plate), "_Standard_Curve.pdf")
      Save_Name <- file.path(Output_Directory, Save_Name)
      ggsave(Save_Name, plot = p, height = 3 * 3, width = 5 * 4)
      
      # Further processing of the Plate object if needed
      
      
      # Fitting Data To Standarad Curve ----------------------------------------
      Plate <- Plate %>% 
        filter(CONDITIONS != "CALIBRATION") %>% 
        mutate(Plate = as.numeric(gsub("^Plate_(\\d+)_\\d{8}$", "\\1", basename(input_plate_dir))),
               Date  = as_date(str_extract(basename(input_plate_dir), "\\d{8}")),
               Dilution_Factor = case_when(#Date <= "2023-01-01" ~ DILUTION_FACTOR_2,
                 Date <= "2023-05-09" ~ DILUTION_FACTOR_5,
                 Date >  "2023-05-09" ~ DILUTION_FACTOR_10),
               MEASUREMENTS = as.numeric(MEASUREMENTS),
               IL2_concentration = (Fit$coefficients[1]*MEASUREMENTS),
               IL2_concentration_DILUTION_FACTOR = IL2_concentration*Dilution_Factor)
      
      All_plates_data <- rbind(All_plates_data, Plate)
    }
  } else {
    print("No plates found!")
  }
  return(All_plates_data)
}