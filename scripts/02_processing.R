#          .                                                   ####  ##        #
#       ":"                               ####              ###########        #
#     ___:____     |"\/"|               ########              #######          #
#   ,'        `.    \  /                  #####                                #
#   |  O        \___/  |                                                       #
# ~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~#

### LOAD LIBRARIES
library(tidyverse)
library(data.table)
library(jsonlite)
library(janitor)
library(ggrepel)
library(UniprotR)

### LOAD FUNCTIONS
source(file = ifelse(exists("https://raw.githubusercontent.com/tlobnow/master_thesis/main/scripts/functions.R"), 
                     yes =  "https://raw.githubusercontent.com/tlobnow/master_thesis/main/scripts/functions.R",
                     no  =  "~/Documents/Github/master_thesis/scripts/functions.R"))

# STRUCTURE OF YOUR FOLDERS SHOULD LOOK AS FOLLOWS FOR JSON EXTRACTION:
#
# MAIN
#   |____FOLDER
#           |____FILE
#                 |____JSON
#                 |     |____FILE_ranking_model_1.json
#                 |     |____FILE_ranking_model_2.json
#                 |     |____FILE_ranking_model_3.json
#                 |     |____FILE_ranking_model_4.json
#                 |     |____FILE_ranking_model_5.json
#                 |     
#                 |____SLURMS
#                 |     |____ ...
#                 |
#                 |____UNRLXD
#                 |     |____ ...
#                 |
#                 |...

### SET MODES

JSON_XTRCT    = T
PROCESS_JSON  = T
ANNOTATE      = T

# avoid SLURM extraction (unfortunately quite error-prone..)
SLURM_XTRCT   = F
PROCESS_SLURM = F


### DEFINE PATHS
MAIN    = ifelse(dir.exists("/Volumes/TAYLOR-LAB/Finn/RESULTS/IP_MS_2/"), 
                 yes =  "/Volumes/TAYLOR-LAB/Finn/RESULTS/IP_MS_2/",
                 no  =  "~/Documents/Github/transferGit/")

FOLDER  = "MYD88"
# FOLDER  = "IRAK4"
# FOLDER  = "IRAK1"
# FOLDER  = "MYD88_signif"


FILES_LOC          = "~/Documents/Github/master_thesis/"
OUT                = paste0(FILES_LOC, FOLDER)
RAW_SUMMARIES_PATH = paste0(FILES_LOC, "raw_summaries/")        # preprocessed files (raw extracted json, slurms)
SUMMARIES_PATH     = paste0(FILES_LOC, "summaries/")            # processed files
TAXA_PATH          = paste0(FILES_LOC, "taxa_lists/")
ANNOTATED_PATH     = paste0(FILES_LOC, "summaries_annotated/")  # added taxa info


################################################################################
################################################################################
################################################################################

### EXTRACT JSON FILES
if (JSON_XTRCT == T) {

  # remove pre-existing csv file, append would lead to duplicate rows
  system(command = paste0(" [ -f ", RAW_SUMMARIES_PATH, FOLDER ,"_fromJSON.csv ] && rm ", RAW_SUMMARIES_PATH, FOLDER, "_fromJSON.csv"))

  LOCATION = paste0(MAIN, FOLDER)
  LIST = list.files(LOCATION)
  
  for (FILE in LIST) {
    print(paste0("processing ", FILE))
    maxJSON=list.files(paste0(LOCATION, "/" , FILE, "/JSON/"))
    
    # replace all "Infinity" strings with large number (9999)
    system(command = paste0("grep -rl Infinity ", MAIN, FOLDER,"/", FILE, "/JSON/", " | xargs sed -i '' -e 's/Infinity/9999/g'"))
    
    if (length(maxJSON) > 0) {
      for (i in 1:length(maxJSON)) {
        JSON   = paste0(LOCATION, "/" , FILE, "/JSON/", maxJSON[i])
        OUT    = paste0(RAW_SUMMARIES_PATH, FOLDER)
        jsonExtract(JSON = JSON, OUT = OUT, FILE = FILE)
      } 
    } else {
      next
      print(paste0("skipped", FILE))
    }
  }
}

### PROCESS JSON SUMMARY FILE
if (PROCESS_JSON == T) {
  # READ FILE
  JSON_EXTRACT  <- read.csv(paste0(RAW_SUMMARIES_PATH, FOLDER, "_fromJSON.csv"), header = F) %>% 
    mutate(ORIGIN = FOLDER) %>% 
    # RENAME COLUMNS
    setnames(old = c("V1",    "V2",     "V3",    "V4",   "V5",   "V6",   "V7",   "V8",     "V9",  "V10",  "V11",        "V12",          "V13"), 
             new = c("FILE", "MODEL", "RECYCLE", "TOL", "pLDDT", "pTM", "piTM", "iScore", "iRes", "iCnt", "FILE_MODEL", "NUM_CLUSTERS", "N_MONOMERS"), 
             skip_absent = T) 
  # ADD RANK COLUMN
  JE <- JSON_EXTRACT %>% 
    mutate(FILE_RECYCLE = paste0(FILE_MODEL, "_", RECYCLE), RANK = NA) %>% 
    distinct(FILE_RECYCLE, .keep_all = T) %>% 
    group_by(FILE) %>% 
    mutate(RANK = frank(desc(iScore), ties.method = "min"))
  # SAVE CSV FILE
  #write.csv(JE, paste0(FILES_LOC, "df_", FOLDER,"_fromJSON.csv"), row.names = F)
  write.csv(JE, paste0(SUMMARIES_PATH, FOLDER,"_fromJSON.csv"), row.names = F)
}

# ### EXTRACT SLURM FILES
# if (SLURM_XTRCT == T) {
#   system(command = paste0(" [ -f ", RAW_SUMMARIES_PATH, FOLDER ,"_fromSLURM.csv ] && rm ", RAW_SUMMARIES_PATH, FOLDER ,"_fromSLURM.csv"))
#   LOCATION = paste0(MAIN, FOLDER)
#   LIST = list.files(LOCATION, pattern = "_x")
#   for (FILE in LIST) {
#     SLURM = paste0(LOCATION, "/", FILE, "/slurm.out")
#     OUT   = paste0(RAW_SUMMARIES_PATH, FOLDER)
#     rs <- tryCatch(slurmExtract(SLURM = SLURM, OUT = OUT), 
#                    error=function(e) NULL)
#     if (is.null(rs)){
#       print(paste(" -(^o^)- ", FILE, "loaded!"))
#     }
#     else{
#       print(paste(" /(x.x)\ ", FILE, "missing!"))
#       next
#     }
#   } 
# }
#
# ### PROCESS SLURM SUMMARY FILE
# if (PROCESS_SLURM == T) {
#   
#   SLURM_EXTRACT  <- read.csv(paste0(RAW_SUMMARIES_PATH, FOLDER,"_fromSLURM.csv"), header = F) %>% 
#     mutate(ORIGIN = FOLDER) %>% setnames(old = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"), new = c("FILE", "MODEL", "TOL", "pLDDT", "pTM", "piTM", "iScore", "iRes", "iCnt", "FILE_MODEL", "NUM_CLUSTERS", "N_MONOMERS"), skip_absent = T) 
#   
#   # set up the FILE_MODEL column for joining with LOG file
#   SLURM_EXTRACT <- SLURM_EXTRACT %>% mutate(FILE_MODEL = paste0(FILE, "_x", N_MONOMERS, "_", MODEL))  
#   SLURM_EXTRACT$MODEL <- as.numeric(unlist(lapply(strsplit(SLURM_EXTRACT$MODEL, "_", fixed=TRUE), function(x) return(x[2])))) # retain model number only (as numeric column)
#   
#   SE  <- SLURM_EXTRACT %>% unique()
#   
#   # ADD RANK COLUMN
#   SE <- SE %>% mutate(FILE_N = paste0(FILE, "_" , N_MONOMERS))
#   SE$RANK <- ave(-SE$pTM, SE$FILE_N, FUN = rank)
#   SE <- SE %>% group_by(FILE_N) %>% 
#     #mutate(RANK = frank(desc(iScore), ties.method = "max"),
#     #       inv.RANK = 6-RANK,
#     #       VIS = case_when(iScore == 0 | piTM == 0 ~ "none")) %>% 
#     distinct(FILE_MODEL, .keep_all = T)
#   
#   write.csv(SE, paste0(SUMMARIES_PATH, FOLDER,"_fromSLURM.csv"), row.names = F)
#   
# } 


# USE BIG ANNOTATION FILES (CONTAIN ALL SIGNIFICANT PROTEINS THAT CAME DOWN WITH MYD88, IRAK4, IRAK1 RESPECTIVELY)
if (ANNOTATE == T) {
  if (file.exists(paste0(SUMMARIES_PATH, FOLDER,"_fromJSON.csv")) == T) {
    system(command = paste0(" [ -f ", ANNOTATED_PATH, FOLDER ,"_fromJSON_ANNOTATED.csv ] && rm ", ANNOTATED_PATH, FOLDER ,"_fromJSON_ANNOTATED.csv"))
    MAIN  <- fread(paste0(SUMMARIES_PATH, FOLDER,"_fromJSON.csv")) %>% remove_empty("cols") %>% filter(RANK == 1)
    Entry <- unlist(lapply(strsplit(MAIN$FILE, "_", fixed=TRUE), function(x) return(x[2])))
    MAIN  <- cbind(MAIN, Entry)
    
    if (file.exists(paste0(TAXA_PATH, FOLDER,"_TAXA.csv")) == T) {
    #if (file.exists(paste0(TAXA_PATH, "MYD88_TAXA.csv")) == T) {
      TAXA <- fread(paste0(TAXA_PATH, FOLDER,"_TAXA.csv"))
      #TAXA <- fread(paste0(TAXA_PATH, "MYD88_TAXA.csv"))
    } else {
      TAXA  <- GetNamesTaxa(ProteinAccList = MAIN$Entry)
    }
    
    ANNOTATED <- left_join(MAIN, TAXA) %>% select(-c(FILE_MODEL, ORIGIN, Entry, Virus.hosts, Gene.Names..ordered.locus., Gene.Names..ORF.))
    write.csv(ANNOTATED, paste0(ANNOTATED_PATH, FOLDER ,"_fromJSON_ANNOTATED.csv"))
  } 
  # if (file.exists(paste0(SUMMARIES_PATH, FOLDER,"_fromSLURM.csv")) == T) {
  #   system(command = paste0(" [ -f ", ANNOTATED_PATH, FOLDER ,"_fromSLURM_ANNOTATED.csv ] && rm ", ANNOTATED_PATH, FOLDER ,"_fromSLURM_ANNOTATED.csv"))
  #   MAIN  <- fread(paste0((paste0(SUMMARIES_PATH, FOLDER,"_fromSLURM.csv")))) %>% remove_empty("cols") %>% filter(RANK == 1)
  #   Entry <- unlist(lapply(strsplit(MAIN$FILE, "_", fixed=TRUE), function(x) return(x[2])))
  #   MAIN  <- cbind(MAIN, Entry)
  #   
  #   if (file.exists(paste0(TAXA_PATH, FOLDER,"_TAXA.csv")) == T) {
  #     TAXA <- fread(paste0(TAXA_PATH, FOLDER,"_TAXA.csv"))
  #   } else {
  #     TAXA  <- GetNamesTaxa(ProteinAccList = MAIN$Entry)
  #   }
  #   ANNOTATED <- left_join(MAIN, TAXA) %>% select(-c(FILE_MODEL, ORIGIN, Entry, Virus.hosts, Gene.Names..ordered.locus., Gene.Names..ORF.))
  #   write.csv(ANNOTATED, paste0(ANNOTATED_PATH, FOLDER ,"_fromSLURM_ANNOTATED.csv"))
  # }
} else print("Please check whether processed slurm or JSON files exist in respective folders.")




