### LOAD LIBRARIES
library(tidyverse)
library(data.table)
library(jsonlite)
library(janitor)
library(ggrepel)
library(UniprotR)
library(tidyverse)
library(data.table)
library(janitor)
library(ggrepel)
library(knitr)
library(svglite)
library(readxl)
library(UniprotR)
library(plotly)


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
# MAIN    = ifelse(dir.exists("/Volumes/TAYLOR-LAB/Finn/RESULTS/DHF119/"), 
#                  yes =  "/Volumes/TAYLOR-LAB/Finn/RESULTS/DHF119/",
#                  no  =  "~/Documents/Github/transferGit/")
MAIN = "~/Documents/Github/transferGit/DHF119/"
FOLDER = "DHF119"
FILES_LOC          = "~/Documents/Github/master_thesis/"
OUT                = paste0(FILES_LOC, FOLDER)
RAW_SUMMARIES_PATH = paste0(FILES_LOC, "raw_summaries/")        # preprocessed files (raw extracted json, slurms)
SUMMARIES_PATH     = paste0(FILES_LOC, "summaries/")            # processed files
TAXA_PATH          = paste0(FILES_LOC, "taxa_lists/")
ANNOTATED_PATH     = paste0(FILES_LOC, "summaries_annotated/")  # added taxa info

ANNOTATE = T

OUT    = paste0(FILES_LOC, FOLDER)

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
        JSON  = paste0(LOCATION, "/" , FILE, "/JSON/", maxJSON[i])
        OUT   = paste0(RAW_SUMMARIES_PATH, FOLDER)
        FILE  = FILE
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


if (ANNOTATE == T) {
  if (file.exists(paste0(FILES_LOC, "summaries_annotated/", FOLDER,"_fromJSON_ANNOTATED.csv"))) {
    MAIN <- read.csv(paste0(FILES_LOC, "summaries_annotated/", FOLDER,"_fromJSON_ANNOTATED.csv"))
  } else if (file.exists(paste0(FILES_LOC, "summaries_annotated/", FOLDER,"fromSLURM_ANNOTATED.csv"))) {
    MAIN <- read.csv(paste0(FILES_LOC, "summaries_annotated/", FOLDER,"fromSLURM_ANNOTATED.csv"))
  } else print("Please make sure correct path was provided or change to ANNOTATE = F.")
} else {
  MAIN <- fread(paste0(FILES_LOC, "summaries/", FOLDER,"_fromJSON.csv"))
  if (file.exists(paste0(FILES_LOC, "summaries/", FOLDER,"_fromJSON.csv"))) {
    MAIN <- read.csv(paste0(FILES_LOC, "summaries/", FOLDER,"_fromJSON.csv"))
    MAIN$RECYCLE_num <-unlist(lapply(strsplit(MAIN$RECYCLE, "_", fixed=TRUE), function(x) return(x[2])))
  } else if (file.exists(paste0(FILES_LOC, "summaries/", FOLDER,"fromSLURM_.csv"))) {
    MAIN <- read.csv(paste0(FILES_LOC, "summaries/", FOLDER,"fromSLURM_.csv"))
  } else print("Please make sure correct path was provided or change to ANNOTATE = F.")
}

JE$RECYCLE_BOOL <- F
JE$RECYCLE_BOOL[grep("*recycled*", JE$RECYCLE)] <- T

MAIN_PLOT_all_rcycld <- JE %>% 
  ggplot(aes(iScore, piTM, col = FILE)) +
  geom_abline(col = "gray") +
  geom_point(size = 3) +
  scale_x_continuous(name = "iScore", breaks = c(0, 0.4, 0.5, 0.7, 1)) +
  scale_y_continuous(name = "piTM", breaks = c(0, 0.5, 1)) +
  theme(legend.position = "bottom") +
  scale_color_discrete(name = "Control") +
  #ggtitle(paste0('Computational screening for PPI partners of ', FOLDER)) +
  ggtitle(paste0("iScores plotted for all recycles of ", FOLDER)) +
  expand_limits(x=c(0,1), y=c(0,1)) +
  annotate("text", x = 0.4, y = -0.05, label = "medium \n confidence") +
  annotate("text", x = 0.5, y = -0.05, label = "high \n confidence") +
  annotate("text", x = 0.7, y = -0.05, label = "very high \n confidence")

#MAIN_PLOT
ggplotly(MAIN_PLOT_all_rcycld)

MAIN_PLOT_pdb <- JE %>% 
  filter(RECYCLE_BOOL == F) %>%
  ggplot(aes(iScore, piTM, col = FILE)) +
  geom_abline(col = "gray") +
  geom_point(size = 3) +
  scale_x_continuous(name = "iScore", breaks = c(0, 0.4, 0.5, 0.7, 1)) +
  scale_y_continuous(name = "piTM", breaks = c(0, 0.5, 1)) +
  theme(legend.position = "bottom") +
  scale_color_discrete(name = "Control") +
  #ggtitle(paste0('Computational screening for PPI partners of ', FOLDER)) +
  ggtitle(paste0("iScores plotted for final model files (pdb available) of ", FOLDER)) +
  expand_limits(x=c(0,1), y=c(0,1)) +
  annotate("text", x = 0.4, y = -0.05, label = "medium \n confidence") +
  annotate("text", x = 0.5, y = -0.05, label = "high \n confidence") +
  annotate("text", x = 0.7, y = -0.05, label = "very high \n confidence")

#MAIN_PLOT
ggplotly(MAIN_PLOT_pdb)

## HISTOGRAM
JE %>% #group_by(FILE) %>% filter(iScore == max(iScore)) %>% 
  ggplot(aes(iScore)) +
  geom_vline(xintercept = 0.4, col = "gray40", linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = 0.5, col = "cornflowerblue", linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = 0.7, col = "lightgreen", linetype = "dotted", linewidth = 1) +
  geom_histogram(bins = 30) +
  expand_limits(x=c(0,1)) +
  annotate("text", x = 0.4, y = -0.05, label = "medium \n confidence") +
  annotate("text", x = 0.5, y = -0.05, label = "high \n confidence") +
  annotate("text", x = 0.7, y = -0.05, label = "very high \n confidence") +
  ggtitle("Model Confidence assessed by iScore")
