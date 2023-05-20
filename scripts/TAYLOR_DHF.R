### LOAD LIBRARIES
library(tidyverse)
library(data.table)
library(jsonlite)
library(janitor)
library(ggrepel)
library(UniprotR)
library(knitr)
library(svglite)
library(readxl)
library(plotly)
library(fs)
library(stringr)


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
MAIN_FOLDER        = "DHF"
FOLDER             = "DHF91"
FOLDER             = "DHF119"

FILES_LOC          = "~/Documents/Github/master_thesis/"
MAIN               = ifelse(dir.exists(paste0("/Volumes/TAYLOR-LAB/Finn/RESULTS/", MAIN_FOLDER, "/")), 
                            yes =  paste0("/Volumes/TAYLOR-LAB/Finn/RESULTS/", MAIN_FOLDER, "/"),
                            no  =  "~/Documents/Github/transferGit/")

OUT                = paste0(FILES_LOC, FOLDER)
RAW_SUMMARIES_PATH = paste0(FILES_LOC, "raw_summaries/")        # preprocessed files (raw extracted json, slurms)
SUMMARIES_PATH     = paste0(FILES_LOC, "summaries/")            # processed files
TAXA_PATH          = paste0(FILES_LOC, "taxa_lists/")
ANNOTATED_PATH     = paste0(FILES_LOC, "summaries_annotated/")  # added taxa info


### EXTRACTION & PROCESSING
if (JSON_XTRCT   == T) {source(paste0(FILES_LOC, "scripts/JSON_XTRCT.R"))}
if (PROCESS_JSON == T) {source(paste0(FILES_LOC, "scripts/PROCESS_JSON.R"))}
if (ANNOTATE     == T) {
  source(paste0(FILES_LOC, "scripts/ANNOTATE_1.R"))
} else source(paste0(FILES_LOC, "scripts/ANNOTATE_2.R"))

### PLOTS
## MAIN PLOT GEOM_POINT
source("~/Documents/Github/master_thesis/scripts/MAIN_PLOT_iScore_piTM_point.R")
ggplotly(MAIN_PLOT_POINT)

## MAIN PLOT HISTOGRAM
source("~/Documents/Github/master_thesis/scripts/MAIN_PLOT_iScore_piTM_hist.R")
MAIN_PLOT_HIST
