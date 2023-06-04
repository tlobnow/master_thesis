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
library(ggalt)


### LOAD FUNCTIONS
source(file = ifelse(exists("https://raw.githubusercontent.com/tlobnow/master_thesis/main/scripts/functions.R"), 
                     yes =  "https://raw.githubusercontent.com/tlobnow/master_thesis/main/scripts/functions.R",
                     no  =  "~/Documents/Github/master_thesis/scripts/functions.R"))

### SET MODES
JSON_XTRCT    = T
JSON_PROCESS  = T
ANNOTATE      = T
# avoid SLURM extraction (unfortunately quite error-prone..)
SLURM_XTRCT   = F
PROCESS_SLURM = F

# STRUCTURE OF YOUR FOLDERS SHOULD LOOK AS FOLLOWS FOR JSON EXTRACTION:
#
# MAIN                              # where are the main folder collections located?
#   |____MAIN_FOLDER                # upper level folder that contains the folder collection (e.g. MYD88, DHF, ..)
#           |____FOLDER             # run folders of interest
#                   |____FILE
#                         |____JSON
#                         |     |____FILE_ranking_model_1.json
#                         |     |____FILE_ranking_model_2.json
#                         |     |____FILE_ranking_model_3.json
#                         |     |____FILE_ranking_model_4.json
#                         |     |____FILE_ranking_model_5.json
#                         |     
#                         |____SLURMS
#                         |     |____ ...
#                         |
#                         |____UNRLXD
#                         |     |____ ...
#                         |
#                         |...

### DEFINE PATHS
MAIN_FOLDER        = "DHF"
FOLDER             = "DHF_ALL"
#FOLDER             = "DHF91"
#FOLDER             = "DHF119"


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
if (JSON_PROCESS == T) {source(paste0(FILES_LOC, "scripts/JSON_PROCESS.R"))}
if (ANNOTATE     == T) {
  source(paste0(FILES_LOC, "scripts/ANNOTATE_1.R"))
} else source(paste0(FILES_LOC, "scripts/ANNOTATE_2.R"))

### PLOTS

JE <- JE %>% mutate(DHF = "missing")
#JE <- separate(data = JE, col = FILE2, into = c("DHF", "rest"), sep = "_")
JE$DHF <- unlist(lapply(strsplit(JE$FILE,  "_", fixed=TRUE), function(x) return(x[1])))

plt <- JE %>%
  #filter(FILE == "DHF119_x6") %>%
  #filter(FILE %in% c("DHF91_x6", "DHF119_x6")) %>%
  ggplot(aes(iScore, piTM)) +
  geom_abline(col = "gray") +
  geom_encircle(aes(col = FILE, fill = FILE, alpha = 0.5), show.legend = F) +
  geom_point(aes(shape = as.factor(MODEL), size = 3), alpha = 1) +
  expand_limits(x=c(0,1), y=c(0,1)) +
  facet_wrap(~DHF)
#ggplotly(plt)
plt

## MAIN PLOT GEOM_POINT
#source("~/Documents/Github/master_thesis/scripts/MAIN_PLOT_iScore_piTM_point.R")
#ggplotly(MAIN_PLOT_POINT)

## MAIN PLOT HISTOGRAM
#source("~/Documents/Github/master_thesis/scripts/MAIN_PLOT_iScore_piTM_hist.R")
#MAIN_PLOT_HIST

