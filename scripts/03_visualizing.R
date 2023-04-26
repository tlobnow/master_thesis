#         .
#       ":"
#     ___:____     |"\/"|
#   ,'        `.    \  /
#   |  O        \___/  |
# ~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^
library(tidyverse)
library(data.table)
library(janitor)
library(ggrepel)
library(knitr)
library(svglite)
library(readxl)
library(UniprotR)
library(plotly)

theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())

### DEFINE PATHS
MAIN    = ifelse(dir.exists("/Volumes/TAYLOR-LAB/Finn/RESULTS/IP_MS_2/"), 
                 yes =  "/Volumes/TAYLOR-LAB/Finn/RESULTS/IP_MS_2/",
                 no  =  "~/Documents/Github/transferGit/IP_MS/")
#MAIN    = "~/Documents/Github/transferGit/"
#MAIN = "/Volumes/TAYLOR-LAB/Finn/RESULTS/"

FOLDER  = "MYD88"
#FOLDER  = "IRAK4"
FOLDER  = "IRAK1"

### DEFINE PATHS
FILES_LOC      = "~/Documents/Github/master_thesis/"
SUMMARIES_PATH = paste0(FILES_LOC, "summaries/")
ANNOTATED_PATH = paste0(FILES_LOC, "summaries_annotated/")  # added taxa info
OUT            = paste0(FILES_LOC, FOLDER)

ANNOTATED = T

if (ANNOTATED == T) {
  
  if (file.exists(paste0(ANNOTATED_PATH, FOLDER,"_fromJSON_ANNOTATED.csv"))) {
    MAIN <- read.csv(paste0(ANNOTATED_PATH, FOLDER,"_fromJSON_ANNOTATED.csv"))
  } else if (file.exists(paste0(ANNOTATED_PATH, FOLDER,"_fromSLURM_ANNOTATED.csv"))) {
    MAIN <- read.csv(paste0(ANNOTATED_PATH, FOLDER,"_fromSLURM_ANNOTATED.csv"))
  } else print("Please make sure correct path was provided or change to ANNOTATE = F.")
  
  CONTROL_SAMPLES <- MAIN %>% filter(FILE %in% c("TNFa_MOUSE_x6",
                                                 "MYD88-DD_x10",
                                                 "MYD88_x1_IRAK4_x1",
                                                 "MYD88_MOUSE_slim3_1-155_x6",
                                                 "MYD88_IL1R_x1",
                                                 "MYD88_MOUSE_x6"))
  MAIN <- MAIN %>% filter(!FILE %in% CONTROL_SAMPLES$FILE)
  
  CONTROL_SAMPLES$N_MONOMERS[ CONTROL_SAMPLES$FILE %in% "TNFa_MOUSE_x6"] <- 3
  
  ## HISTOGRAM
  MAIN %>% ggplot(aes(iScore, fill = as.factor(N_MONOMERS))) +
    geom_vline(xintercept = 0.4, col = "gray40", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0.5, col = "cornflowerblue", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0.7, col = "lightgreen", linetype = "dotted", linewidth = 1) +
    geom_histogram(bins = 30) +
    expand_limits(x=c(0,1)) +
    annotate("text", x = 0.4, y = -0.05, label = "medium \n confidence") +#, angle = 45) +
    annotate("text", x = 0.5, y = -0.05, label = "high \n confidence") +#, angle = 45) +
    annotate("text", x = 0.7, y = -0.05, label = "very high \n confidence") +#, angle = 45)
    ggtitle("Model Confidence assessed by iScore")
  
  SIGNIF_PROTEINS <- MAIN %>% filter(iScore >= 0.4, !is.na(Gene.Names)) %>% distinct(FILE, .keep_all = T)
  MAIN <- MAIN %>% mutate(SIGNIF = case_when(MAIN$FILE %in% SIGNIF_PROTEINS$FILE  ~ T,
                                             !MAIN$FILE %in% SIGNIF_PROTEINS$FILE ~ F))
  max_SIGNIF_PROTEINS <- SIGNIF_PROTEINS %>% group_by(FILE) %>% summarise(iScore = max(iScore))
  
  MAIN_PLOT <- MAIN %>% ggplot(aes(iScore, piTM, 
                      #col = as.factor(N_MONOMERS),
                      label = Entry.Name)) +
    geom_abline(col = "gray") +
    geom_point(size = 3) +
    #geom_point(data = CONTROL_SAMPLES, col = "black", size = 5) +
    #geom_point(data = CONTROL_SAMPLES, aes(col = as.factor(N_MONOMERS)), size = 3) +
    ggrepel::geom_label_repel(show.legend = F) +
    #ggrepel::geom_label_repel(data = CONTROL_SAMPLES, aes(label = Entry.Name), max.overlaps = 2, show.legend = F, col = "black", alpha = 0.3) +
    #ggrepel::geom_label_repel(data = SIGNIF_PROTEINS[SIGNIF_PROTEINS$FILE %in% max_SIGNIF_PROTEINS$FILE], show.legend = F) +
    scale_x_continuous(name = "iScore", breaks = c(0, 0.4, 0.5, 0.7, 1)) +
    scale_y_continuous(name = "piTM", breaks = c(0, 0.5, 1)) +
    theme(legend.position = "bottom") +
    scale_color_discrete(name = "Complex Size") +
    ggtitle(paste0('Computational screening for PPI partners of ', FOLDER)) +
    #facet_wrap(~SIGNIF) +
    expand_limits(x=c(0,1), y=c(0,1)) +
    annotate("text", x = 0.4, y = -0.05, label = "medium \n confidence") +#, angle = 45) +
    annotate("text", x = 0.5, y = -0.05, label = "high \n confidence") +#, angle = 45) +
    annotate("text", x = 0.7, y = -0.05, label = "very high \n confidence")#, angle = 45)
  
  MAIN_PLOT
  
  #ggplotly(MAIN_PLOT)
  
  
  
} else {
  
  MAIN <- fread(paste0(SUMMARIES_PATH, FOLDER,"_fromJSON.csv"))
  if (file.exists(paste0(SUMMARIES_PATH, FOLDER,"_fromJSON.csv"))) {
    MAIN <- read.csv(paste0(SUMMARIES_PATH, FOLDER,"_fromJSON.csv"))
    MAIN$RECYCLE_num <-unlist(lapply(strsplit(MAIN$RECYCLE, "_", fixed=TRUE), function(x) return(x[2])))
  } else if (file.exists(paste0(SUMMARIES_PATH, FOLDER,"_fromSLURM.csv"))) {
    MAIN <- read.csv(paste0(SUMMARIES_PATH, FOLDER,"_fromSLURM.csv"))
  } else print("Please make sure correct path was provided or change to ANNOTATE = F.")
  
  
  CONTROL_SAMPLES <- MAIN %>% filter(FILE %in% c("TNFa_MOUSE_x6",
                                                 "MYD88-DD_x10",
                                                 "MYD88_x1_IRAK4_x1",
                                                 "MYD88_MOUSE_slim3_1-155_x6",
                                                 "MYD88_IL1R_x1",
                                                 "MYD88_MOUSE_x6"))
  MAIN <- MAIN %>% filter(!FILE %in% CONTROL_SAMPLES$FILE)
  
  CONTROL_SAMPLES$N_MONOMERS[ CONTROL_SAMPLES$FILE %in% "TNFa_MOUSE_x6"] <- 3
  
  ## HISTOGRAM
  MAIN %>% ggplot(aes(iScore, fill = as.factor(N_MONOMERS))) +
    geom_vline(xintercept = 0.4, col = "gray40", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0.5, col = "cornflowerblue", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = 0.7, col = "lightgreen", linetype = "dotted", linewidth = 1) +
    geom_histogram(bins = 30) +
    expand_limits(x=c(0,1)) +
    annotate("text", x = 0.4, y = -0.05, label = "medium \n confidence") +#, angle = 45) +
    annotate("text", x = 0.5, y = -0.05, label = "high \n confidence") +#, angle = 45) +
    annotate("text", x = 0.7, y = -0.05, label = "very high \n confidence") +#, angle = 45)
    ggtitle("Model Confidence assessed by iScore")
  
  SIGNIF_PROTEINS <- MAIN %>% filter(iScore >= 0.4) %>% distinct(FILE, .keep_all = T)
  
  MAIN <- MAIN %>% mutate(SIGNIF = case_when(MAIN$FILE %in% SIGNIF_PROTEINS$FILE  ~ T,
                                             !MAIN$FILE %in% SIGNIF_PROTEINS$FILE ~ F))
  max_SIGNIF_PROTEINS <- SIGNIF_PROTEINS %>% group_by(FILE) %>% summarise(iScore = max(iScore))
  
  SIGNIF_PROTEINS %>% ggplot(aes(iScore, piTM,
                                 # col = FILE,
                                 label = FILE)) +
    geom_abline(col = "gray") +
    geom_point(data = MAIN, size = 3, col = "gray40", alpha =.1) +
    geom_point(size = 3) +
    #ggrepel::geom_label_repel(data = SIGNIF_PROTEINS[SIGNIF_PROTEINS$FILE %in% max_SIGNIF_PROTEINS$FILE], show.legend = F) +
    scale_x_continuous(name = "iScore", breaks = c(0, 0.4, 0.5, 0.7, 1)) +
    scale_y_continuous(name = "piTM", breaks = c(0, 0.5, 1)) +
    theme(legend.position = "none") +
    scale_color_discrete(name = "Complex Size") +
    ggtitle(paste0('Computational screening for PPI partners of ', FOLDER)) +
    expand_limits(x=c(0,1), y=c(0,1)) +
    annotate("text", x = 0.4, y = -0.05, label = "medium \n confidence") +#, angle = 45) +
    annotate("text", x = 0.5, y = -0.05, label = "high \n confidence") +#, angle = 45) +
    annotate("text", x = 0.7, y = -0.05, label = "very high \n confidence")#, angle = 45)
}

