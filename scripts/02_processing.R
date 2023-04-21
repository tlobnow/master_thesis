#         .
#       ":"
#     ___:____     |"\/"|
#   ,'        `.    \  /
#   |  O        \___/  |
# ~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^

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
SLURM_XTRCT   = T
JSON_XTRCT    = T
PROCESS_SLURM = T
PROCESS_JSON  = T
ANNOTATE      = T


### DEFINE PATHS
MAIN    = ifelse(dir.exists("/Volumes/TAYLOR-LAB/Finn/RESULTS/IP_MS_2/"), 
                 yes =  "/Volumes/TAYLOR-LAB/Finn/RESULTS/IP_MS_2/",
                 no  =  "~/Documents/Github/transferGit/")

#MAIN = "~/Documents/Github/transferGit/"
#MAIN = "/Volumes/TAYLOR-LAB/Finn/RESULTS/"
### DEFINE PATHS
#MAIN    = "~/Desktop/"
FOLDER  = "MYD88"
#FOLDER  = "IRAK4"
#FOLDER  = "IRAK1"

FILES_LOC = "~/Documents/Github/master_thesis/"
OUT     = paste0(FILES_LOC, FOLDER)

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
    OUT   = paste0(FILES_LOC, "SLURM_EXTRACTED/", FOLDER)
    rs <- tryCatch(slurmExtract(SLURM = SLURM, OUT = OUT), 
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


### EXTRACT JSON FILES
if (JSON_XTRCT == T) {
  # replace all "Infinity" strings with large number (9999)
  #system(command = paste0("grep -rl Infinity ", MAIN, FOLDER, " | xargs sed -i '' -e 's/Infinity/9999/g'"))
  # remove pre-existing csv file, append would lead to duplicate rows
  system(command = paste0(" [ -f ", OUT ,"_fromJSON.csv ] && rm ", OUT, "_fromJSON.csv"))
  
  # GENERATE FILE
  LOCATION = paste0(MAIN, FOLDER)
  #LIST = list.files(LOCATION, pattern = "_x1")
  #LIST = list.files(LOCATION, pattern = "_x")
  LIST = list.files(LOCATION)
  FILE="MYD88_A2A791_x1"
  for (FILE in LIST) {
    print(paste0("processing ", FILE))
    maxJSON=list.files(paste0(LOCATION, "/" , FILE, "/JSON/"))
    #system(command = paste0("grep -rl Infinity ", MAIN, FOLDER,"/", FILE, " | xargs sed -i '' -e 's/Infinity/9999/g'"))
    if (length(maxJSON) > 0) {
      for (i in 1:length(maxJSON)) { # 1:5) {#paste0(LOCATION, "/" , FILE, "/JSON/*.json")) {
        JSON  = paste0(LOCATION, "/" , FILE, "/JSON/", FILE, "_ranking_model_",i,".json")
        OUT   = paste0(FILES_LOC, "JSON_EXTRACTED/", FOLDER)
        jsonExtract_1.1_ratio(JSON = JSON, OUT = OUT)
      } 
    } else {
      next
      print(paste0("skipped", FILE))
    }
  }
}

### PROCESS SLURM SUMMARY FILE
if (PROCESS_SLURM == T) {
  
  SLURM_EXTRACT  <- read.csv(paste0(FILES_LOC, "SLURM_EXTRACTED/", FOLDER,".csv"), header = F) %>% 
    mutate(ORIGIN = FOLDER) %>% setnames(old = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"), new = c("FILE", "MODEL", "TOL", "pLDDT", "pTM", "piTM", "iScore", "iRes", "iCnt", "FILE_MODEL", "NUM_CLUSTERS", "N_MONOMERS"), skip_absent = T) 
  
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
  
  write.csv(SE, paste0(FILES_LOC, "SLURM_PROCESSED/", FOLDER,".csv"), row.names = F)
  
} 

### PROCESS JSON SUMMARY FILE
if (PROCESS_JSON == T) {
  # READ FILE
  JSON_EXTRACT  <- read.csv(paste0(FILES_LOC, "JSON_EXTRACTED/", FOLDER, "_fromJSON.csv"), header = F) %>% 
    mutate(ORIGIN = FOLDER) %>% 
    setnames(old = c("V1",    "V2",     "V3",    "V4",   "V5",   "V6",   "V7",   "V8",     "V9",  "V10",  "V11",        "V12",          "V13"), 
             new = c("FILE", "MODEL", "RECYCLE", "TOL", "pLDDT", "pTM", "piTM", "iScore", "iRes", "iCnt", "FILE_MODEL", "NUM_CLUSTERS", "N_MONOMERS"), 
             skip_absent = T) 
  # ADD RANK COLUMN
  JE <- JSON_EXTRACT %>% 
    mutate(FILE_RECYCLE = paste0(FILE_MODEL, "_", RECYCLE), RANK = NA) %>% 
    distinct(FILE_RECYCLE, .keep_all = T) %>% 
    group_by(FILE) %>% 
    mutate(RANK = frank(desc(iScore), ties.method = "max"))
  # SAVE CSV FILE
  #write.csv(JE, paste0(FILES_LOC, "df_", FOLDER,"_fromJSON.csv"), row.names = F)
  write.csv(JE, paste0(FILES_LOC, "JSON_PROCESSED/", FOLDER,"_fromJSON.csv"), row.names = F)
}

# USE BIG ANNOTATION FILES (CONTAIN ALL SIGNIFICANT PROTEINS THAT CAME DOWN WITH MYD88, IRAK4, IRAK1 RESPECTIVELY)
if (ANNOTATE == T) {
  if (file.exists(paste0(FILES_LOC, "JSON_PROCESSED/", FOLDER,"_fromJSON.csv")) == T) {
    system(command = paste0(" [ -f ", FILES_LOC, "ANNOTATED/", FOLDER ,"_fromJSON_ANNOTATED.csv ] && rm ", FILES_LOC, "ANNOTATED/", FOLDER ,"_fromJSON_ANNOTATED.csv"))
    MAIN  <- fread(paste0(FILES_LOC, "JSON_PROCESSED/", FOLDER,"_fromJSON.csv")) %>% remove_empty("cols") %>% filter(RANK == 1)
  } else if (file.exists(paste0(FILES_LOC, "SLURM_PROCESSED/", FOLDER,".csv")) == T) {
    system(command = paste0(" [ -f ", FILES_LOC, "ANNOTATED/", FOLDER ,"_ANNOTATED.csv ] && rm ", FILES_LOC, "ANNOTATED/", FOLDER ,"_ANNOTATED.csv"))
    MAIN  <- fread(paste0((paste0(FILES_LOC, "SLURM_PROCESSED/", FOLDER,".csv")))) %>% remove_empty("cols") %>% filter(RANK == 1)
  } else print("Please make sure either processed JSON or processed SLURMs exist in the designated folder..")
    
    Entry <- unlist(lapply(strsplit(MAIN$FILE, "_", fixed=TRUE), function(x) return(x[2])))
    MAIN  <- cbind(MAIN, Entry)
    
    if (file.exists(paste0(FILES_LOC, "TAXA_INFO/", FOLDER,"_TAXA.csv")) == T) {
      TAXA <- fread(paste0(FILES_LOC, "TAXA_INFO/", FOLDER,"_TAXA.csv"))
    } else {
      TAXA  <- GetNamesTaxa(ProteinAccList = MAIN$Entry)
    }
    
    ANNOTATED <- left_join(MAIN, TAXA) %>% select(-c(FILE_MODEL, ORIGIN, Entry, Virus.hosts, Gene.Names..ordered.locus., Gene.Names..ORF.))
    
    if (file.exists(paste0(FILES_LOC, "JSON_PROCESSED/", FOLDER,"_fromJSON.csv")) == T) {
      write.csv(ANNOTATED, paste0(FILES_LOC, "ANNOTATED/", FOLDER ,"_fromJSON_ANNOTATED.csv"))
    } else if (file.exists(paste0(FILES_LOC, "SLURM_PROCESSED/", FOLDER,".csv")) == T) {
      write.csv(ANNOTATED, paste0(FILES_LOC, "ANNOTATED/", FOLDER ,"_ANNOTATED.csv"))
    } else print("Could not annotate.")
    
} else print("Please check whether processed slurm or JSON files exist in respective folders.")




