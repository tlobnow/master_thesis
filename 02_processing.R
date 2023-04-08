################################################################################
################################################################################
################################################################################

### LOAD LIBRARIES
library(tidyverse)
library(data.table)
library(jsonlite)
library(janitor)
library(ggrepel)
library(UniprotR)


### LOAD EXTERNAL FUNCTIONS
source(file = ifelse(exists("https://raw.githubusercontent.com/tlobnow/master_thesis/main/functions.R"), 
                     yes =  "https://raw.githubusercontent.com/tlobnow/master_thesis/main/functions.R",
                     no  =  "~/Documents/Github/master_thesis/functions.R"))

### SET MODES
SLURM_XTRCT   = F
CHIMERA_XTRCT = F
JSON_XTRCT    = F
PROCESS_SLURM = T
PROCESS_JSON  = T
ANNOTATE      = T


### DEFINE PATHS
MAIN    = ifelse(dir.exists("/Volumes/TAYLOR-LAB/Finn/RESULTS/IP_MS/"), 
                 yes =  "/Volumes/TAYLOR-LAB/Finn/RESULTS/IP_MS/",
                 no  =  "~/Documents/Github/transferGit/IP_MS/")
FOLDER  = "MYD88"
CSV_LOC = "~/Documents/Github/master_thesis/files_csv/"
TXT_LOC = "~/Documents/Github/master_thesis/files_txt/"
OUT     = paste0(CSV_LOC, FOLDER)

################################################################################
################################################################################
################################################################################

### EXTRACT SLURM FILES
if (SLURM_XTRCT == T) {
  system(command = paste0(" [ -f ", OUT ,".csv ] && rm ", OUT, ".csv"))
  LOCATION = paste0(MAIN, FOLDER)
  LIST = list.files(LOCATION, pattern = "_x")
  for (FILE in LIST) {
    SLURM = paste0(LOCATION, "/", FILE, "/slurm.out")
    OUT   = paste0(CSV_LOC, FOLDER)
    rs <- tryCatch(slurmExtract(SLURM = SLURM, OUT = OUT), 
                   #finally = print("whoops"),
                   error=function(e) NULL)
    if (is.null(rs)){
      print(paste(" -(^o^)- ", FILE, "loaded!"))
    }
    else{
      print(paste(" /(x.x)\ ", FILE, "missing!"))
      next
    }
  } 
}

### EXTRACT CHIMERA FILES
if (CHIMERA_XTRCT == T) {
  log <- read.table(paste0(TXT_LOC, FOLDER, "_log.txt"), header = F, sep = "\t")
  
  log$info <- NA
  log$info[grep("Chain information for *", log$V1)] <- T
  
  log$buried <- NA
  log$buried[grep("* buried areas: *",     log$V1)] <- T
  
  # extract the INFO lines and save in 1-col-df
  INFO   <- log %>% filter(log$info   == T) %>% mutate(INFO   = V1) %>% select(INFO) 
  
  # extract the BURIED lines and save in 1-col-df
  BURIED <- log %>% filter(log$buried == T) %>% mutate(BURIED = V1) %>% select(BURIED)
  
  # Collect and manipulate data from the LOG DATA
  LOG <- bind_cols(INFO, BURIED) # bind cols - due to same order they now align according w/ info
  LOG$INFO <- unlist(lapply(strsplit(LOG$INFO, " ", fixed=TRUE), function(x) return(x[4])))
  LOG$INFO <- unlist(lapply(strsplit(LOG$INFO, ".", fixed=TRUE), function(x) return(x[1])))
  
  # extract names as INFO
  LOG$FILE <- unlist(lapply(strsplit(LOG$INFO, "_rlx", fixed=TRUE), function(x) return(x[1])))
  
  # extract interface information
  LOG$INTERFACES <- unlist(lapply(strsplit(LOG$BURIED, ": ", fixed=TRUE), function(x) return(x[2])))
  
  # separate interface info
  LOG <- separate(data = LOG, col = INTERFACES, sep = ", ", into = c("INTERFACE_1", "INTERFACE_2", "INTERFACE_3", "INTERFACE_4", "INTERFACE_5", "INTERFACE_6", "INTERFACE_7", "INTERFACE_8", "INTERFACE_9", "INTERFACE_10","INTERFACE_11", "INTERFACE_12", "INTERFACE_13", "INTERFACE_14", "INTERFACE_15", "INTERFACE_16", "INTERFACE_17", "INTERFACE_18", "INTERFACE_19", "INTERFACE_20", "INTERFACE_21", "INTERFACE_22", "INTERFACE_23", "INTERFACE_24", "INTERFACE_25", "INTERFACE_26", "INTERFACE_27", "INTERFACE_28", "INTERFACE_29", "INTERFACE_30"), convert = T, remove = T)
  
  # extract the number of buried interfaces
  LOG$BURIED <- as.numeric(unlist(lapply(strsplit(LOG$BURIED, " ", fixed=TRUE), function(x) return(x[1]))))
  
  # retain unique samples only
  LOG       <- unique(LOG) 
  LOG$MODEL <- unlist(lapply(strsplit(LOG$INFO, "_model_", fixed=TRUE), function(x) return(x[2])))
  LOG$MODEL <- unlist(lapply(strsplit(LOG$MODEL, "_rlx", fixed=TRUE), function(x) return(x[1])))
  
  # write csv file
  write.csv(LOG, paste0(CSV_LOC, "INTERFACES_", FOLDER,".csv"), row.names = F)
}

### EXTRACT JSON FILES
if (JSON_XTRCT == T) {
  # GENERATE FILE
  system(command = paste0(" [ -f ", OUT ,"_fromJSON.csv ] && rm ", OUT, "_fromJSON.csv"))
  LOCATION = paste0(MAIN, FOLDER)
  LIST = list.files(LOCATION, pattern = "_x1")
  for (FILE in LIST) {
    for (i in 1:5) {
      JSON = paste0(LOCATION, "/" , FILE, "/JSON/", FILE, "_ranking_model_",i,".json")
      OUT   = paste0(CSV_LOC, FOLDER)
      jsonExtract_1.1_ratio(JSON = JSON, OUT = OUT)
    } 
  }
}

### PROCESS SLURM SUMMARY FILE
if (PROCESS_SLURM == T) {
  
  SLURM_EXTRACT  <- read.csv(paste0(CSV_LOC, FOLDER,".csv"), header = F) %>% mutate(ORIGIN = FOLDER) %>% setnames(old = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"), new = c("FILE", "MODEL", "TOL", "pLDDT", "pTM", "piTM", "iScore", "iRes", "iCnt", "FILE_MODEL", "NUM_CLUSTERS", "N_MONOMERS"), skip_absent = T) 
  
  # make sure that there are no "+" in the file names // MUST BE COHERENT WITH CHIMERA XTRACT NAME THOUGH!
  #SLURM_EXTRACT$FILE <- gsub("\\+", replacement = "_", x = SLURM_EXTRACT$FILE)
  # set up the FILE_MODEL column for joining with LOG file
  SLURM_EXTRACT <- SLURM_EXTRACT %>% mutate(FILE_MODEL = paste0(FILE, "_x", N_MONOMERS, "_", MODEL))  
  SLURM_EXTRACT$MODEL <- as.numeric(unlist(lapply(strsplit(SLURM_EXTRACT$MODEL, "_", fixed=TRUE), function(x) return(x[2])))) # retain model number only (as numeric column)
  
  SE  <- SLURM_EXTRACT %>% unique()
  
  # ADD RANK COLUMN
  SE <- SE %>% mutate(FILE_N = paste0(FILE, "_" , N_MONOMERS))
  SE$RANK <- ave(-SE$pTM, SE$FILE_N, FUN = rank)
  SE <- SE %>% group_by(FILE_N) %>% 
    mutate(RANK = frank(desc(iScore), ties.method = "max"),
           inv.RANK = 6-RANK,
           VIS = case_when(iScore == 0 | piTM == 0 ~ "none")) %>% 
    distinct(FILE_MODEL, .keep_all = T)
  
  write.csv(SE, paste0(CSV_LOC, "df_", FOLDER,".csv"), row.names = F)
  
} 

### PROCESS JSON SUMMARY FILE
if (PROCESS_JSON == T) {
  # READ FILE
  JSON_EXTRACT  <- read.csv(paste0(CSV_LOC, FOLDER, "_fromJSON.csv"), header = F) %>% 
    mutate(ORIGIN = FOLDER) %>% 
    setnames(old = c("V1",    "V2",     "V3",    "V4",   "V5",   "V6",   "V7",   "V8",     "V9",  "V10",  "V11",        "V12",          "V13"), 
             new = c("FILE", "MODEL", "RECYCLE", "TOL", "pLDDT", "pTM", "piTM", "iScore", "iRes", "iCnt", "FILE_MODEL", "NUM_CLUSTERS", "N_MONOMERS"), 
             skip_absent = T) 
  # ADD RANK COLUMN
  JE <- JSON_EXTRACT %>% 
    mutate(FILE_RECYCLE = paste0(FILE_MODEL, "_", RECYCLE), RANK = NA)   %>% distinct(FILE_RECYCLE, .keep_all = T) %>% 
    group_by(FILE_MODEL) %>% mutate(RANK = frank(desc(iScore), ties.method = "max")) %>% 
    distinct(FILE_MODEL, .keep_all = T)
  # SAVE CSV FILE
  write.csv(JE, paste0(CSV_LOC, "df_", FOLDER,"_fromJSON.csv"), row.names = F)
}

if (ANNOTATE == T) {
  
  if (file.exists(paste0(CSV_LOC, "df_", FOLDER,".csv"))) {
    MAIN <- read.csv(paste0(CSV_LOC, "df_", FOLDER,".csv")) %>% remove_empty("cols") %>% filter(RANK == 1)
    
    # get accession numnbers (if applicable)
    Entry <- unlist(lapply(strsplit(MAIN$FILE, "_", fixed=TRUE), function(x) return(x[2])))
    
    # bind to Main df
    MAIN  <- cbind(MAIN, Entry)
    
    #Get Taxonomy Information, names that are not accession numbers are automatically excluded with the warning:
    # "Bad request. The resource you requested doesn't exist or There is a problem with your input."
    TaxaObj <- GetNamesTaxa(MAIN$Entry) 
    
    # join TaxaObj with MAIN to create annotated df
    if (length(TaxaObj) != 0) {
      ANNOTATED <- full_join(MAIN, TaxaObj) %>% mutate(Entry.Name = case_when(!is.na(Entry.Name) ~ Entry.Name, is.na(Entry.Name) ~ FILE))
      # remove unnecessary columns
      ANNOTATED <- ANNOTATED %>% select(-c(FILE_MODEL, ORIGIN, FILE_N, inv.RANK, VIS, Entry, Virus.hosts, Gene.Names..ordered.locus., Gene.Names..ORF.))
    } else ANNOTATED <- MAIN %>% select(-c(FILE_MODEL, ORIGIN, FILE_N, inv.RANK, VIS))
    # write Annotated df to file
    write.csv(ANNOTATED, paste0(CSV_LOC, "df_", FOLDER,"_ANNOTATED.csv"), row.names = F)
  }
  
  if (file.exists(paste0(CSV_LOC, "df_", FOLDER,"_fromJSON.csv"))) {
    MAIN <- read.csv(paste0(CSV_LOC, "df_", FOLDER,"_fromJSON.csv")) %>% remove_empty("cols") %>% distinct(FILE, .keep_all = T)
    
    # get accession numnbers (if applicable)
    Entry <- unlist(lapply(strsplit(MAIN$FILE, "_", fixed=TRUE), function(x) return(x[2])))
    
    # bind to Main df
    MAIN  <- cbind(MAIN, Entry)
    
    #Get Taxonomy Information, names that are not accession numbers are automatically excluded with the warning:
    # "Bad request. The resource you requested doesn't exist or There is a problem with your input."
    TaxaObj <- GetNamesTaxa(MAIN$Entry) 
    
    # join TaxaObj with MAIN to create annotated df
    if (length(TaxaObj) != 0) {
      ANNOTATED <- full_join(MAIN, TaxaObj) %>% mutate(Entry.Name = case_when(!is.na(Entry.Name) ~ Entry.Name, is.na(Entry.Name) ~ FILE))
      # remove unnecessary columns
      ANNOTATED <- ANNOTATED %>% select(-c(FILE_MODEL, ORIGIN, FILE_N, inv.RANK, Entry, Virus.hosts, Gene.Names..ordered.locus., Gene.Names..ORF.))
    } else ANNOTATED <- MAIN %>% select(-c(FILE_MODEL, ORIGIN, FILE_N, inv.RANK, VIS))
    # SAVE ANNOTATED CSV FILE
    write.csv(ANNOTATED, paste0(CSV_LOC, "df_", FOLDER,"_ANNOTATED.csv"), row.names = F)
  }
  
  
}


