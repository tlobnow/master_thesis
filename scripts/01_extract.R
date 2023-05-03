#          .                                                   ####  ##        #
#       ":"                               ####              ###########        #
#     ___:____     |"\/"|               ########              #######          #
#   ,'        `.    \  /                  #####                                #
#   |  O        \___/  |                                                       #
# ~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~#

### LOAD LIBRARIES
library(tidyverse)
library(data.table)
library(readxl)
library(UniprotR)
library(stringr)


### LOAD FUNCTIONS
source(file = ifelse(exists("https://raw.githubusercontent.com/tlobnow/master_thesis/main/scripts/functions.R"), 
                     yes =  "https://raw.githubusercontent.com/tlobnow/master_thesis/main/scripts/functions.R",
                     no  =  "~/Documents/Github/master_thesis/scripts/functions.R"))

### EXTRACTING INFO FROM FENJA'S IMMUNOPRECIPITATION-MASS SPECTROMETRY (IP-MS) DATA
READ_CSV          = F
JOIN_AND_RETRIEVE = F
TIMING_SUMMARY    = F

EXTRACT_TAXA      = F
WRITE_TAXA_FILES  = F

# DEFINE PATHS
FILES_LOC             = "~/Documents/Github/master_thesis/"
FOLDER_PATH           = "/Volumes/TAYLOR-LAB/Fenja/Mass Spec/Mass Spec analysis/IL-1 proteomics/lists to take a look at/MyD88_IRAK4_IRAK1-IPs/"
FOLDER_PATH2          = "/Volumes/TAYLOR-LAB/Fenja/Mass Spec/Mass Spec analysis/IL-1 proteomics/lists to take a look at/KO_IRAK4: IRAK1 MyD88-IPs/"
DATA_TAY_CONNECTION   = file.exists("/Volumes/TAYLOR-LAB")
FOLDER_PATH_LOCAL     = paste0(FILES_LOC, "raw_csv/")
PATH_ACCESSION_LISTS  = paste0(FILES_LOC, "uniprot_accession_lists/")
TAXA_PATH             = paste0(FILES_LOC, "taxa_lists/")

if (READ_CSV == T) {
  MYD88_min15    <- fread(file = paste0(FOLDER_PATH_LOCAL, "MYD88_min15.csv"))
  MYD88_min30    <- fread(file = paste0(FOLDER_PATH_LOCAL, "MYD88_min30.csv"))
  MYD88_min60    <- fread(file = paste0(FOLDER_PATH_LOCAL, "MYD88_min60.csv"))
  IRAK4_min15    <- fread(file = paste0(FOLDER_PATH_LOCAL, "IRAK4_min15.csv"))
  IRAK4_min30    <- fread(file = paste0(FOLDER_PATH_LOCAL, "IRAK4_min30.csv"))
  IRAK4_min60    <- fread(file = paste0(FOLDER_PATH_LOCAL, "IRAK4_min60.csv"))
  IRAK1_min15    <- fread(file = paste0(FOLDER_PATH_LOCAL, "IRAK1_min15.csv"))
  IRAK1_min30    <- fread(file = paste0(FOLDER_PATH_LOCAL, "IRAK1_min30.csv"))
  IRAK1_min60    <- fread(file = paste0(FOLDER_PATH_LOCAL, "IRAK1_min60.csv"))
  KO_IRAK4_min15 <- fread(file = paste0(FOLDER_PATH_LOCAL, "KO_IRAK4_min15.csv"))
  KO_IRAK4_min30 <- fread(file = paste0(FOLDER_PATH_LOCAL, "KO_IRAK4_min30.csv"))
  KO_IRAK4_min60 <- fread(file = paste0(FOLDER_PATH_LOCAL, "KO_IRAK4_min60.csv"))
  KO_IRAK1_min15 <- fread(file = paste0(FOLDER_PATH_LOCAL, "KO_IRAK1_min15.csv"))
  KO_IRAK1_min30 <- fread(file = paste0(FOLDER_PATH_LOCAL, "KO_IRAK1_min30.csv"))
  KO_IRAK1_min60 <- fread(file = paste0(FOLDER_PATH_LOCAL, "KO_IRAK1_min60.csv")) 
}

if (JOIN_AND_RETRIEVE == T) {
  ### JOIN DATA FRAMES FOR DIFFERENT TIME POINTS
  MYD88    <- join_timepoints(DF1 = MYD88_min15,    DF2 = MYD88_min30,    DF3 = MYD88_min60)    %>% unique()
  KO_IRAK4 <- join_timepoints(DF1 = KO_IRAK4_min15, DF2 = KO_IRAK4_min30, DF3 = KO_IRAK4_min60) %>% unique()
  KO_IRAK1 <- join_timepoints(DF1 = KO_IRAK1_min15, DF2 = KO_IRAK1_min30, DF3 = KO_IRAK1_min60) %>% unique()
  IRAK4    <- join_timepoints(DF1 = IRAK4_min15,    DF2 = IRAK4_min30,    DF3 = IRAK4_min60)    %>% unique()
  IRAK1    <- join_timepoints(DF1 = IRAK1_min15,    DF2 = IRAK1_min30,    DF3 = IRAK1_min60)    %>% unique()
  
  ### RETRIEVE ACCESSION IDs
  retrieveAccessionIDs(DF  = MYD88,    OUT = paste0(PATH_ACCESSION_LISTS, "MYD88_Protein.IDs.txt"))
  retrieveAccessionIDs(DF  = KO_IRAK1, OUT = paste0(PATH_ACCESSION_LISTS, "KO_IRAK1_Protein.IDs.txt"))
  retrieveAccessionIDs(DF  = KO_IRAK4, OUT = paste0(PATH_ACCESSION_LISTS, "KO_IRAK4_Protein.IDs.txt"))
  retrieveAccessionIDs(DF  = IRAK4,    OUT = paste0(PATH_ACCESSION_LISTS, "IRAK4_Protein.IDs.txt"))
  retrieveAccessionIDs(DF  = IRAK1,    OUT = paste0(PATH_ACCESSION_LISTS, "IRAK1_Protein.IDs.txt"))
}

if (TIMING_SUMMARY == T) {
  ### RETRIEVE TIME CHECK MAIN
  COL_SELECTION <- c("Significant", "-LOG(P-value)", "Difference", "Protein.IDs", "ORIGIN", "PULLED_PROTEIN", "KO")
  DF_SUMMARY <- rbind(MYD88[,    ..COL_SELECTION], 
                      KO_IRAK1[, ..COL_SELECTION], 
                      KO_IRAK4[, ..COL_SELECTION], 
                      IRAK4[,    ..COL_SELECTION], 
                      IRAK1[,    ..COL_SELECTION])
  ### MUTATE SIGNIF CLASS COLUMN
  
  # To obtain an IL-1β dependent-myddosome interactome, Fenja was exclusively interested in  proteins enriched in the stimulated conditions (Fig. 7, right side of all scatter plots). She classified proteins as significantly enriched if they met the following criteria: 
    # a) at least two-fold higher abundance after IL-1β stimulation (> 1 on the log2 scale), 
    # b) a -log10 p- value larger than 2 (Fig. 7, y-axis)
    # c) a false-discovery rate (FDR) of ≤ 0.05
  
  # TIME_CHECK <- DF_SUMMARY %>% mutate(SIGNIF_CLASS = case_when(Significant == T & Difference == 10 ~ "NEW SIGNIFICANT", 
  #                                                              Significant == T ~ "SIGNIFICANT", 
  #                                                              Significant == F ~ "ns"))
  
  TIME_CHECK <- DF_SUMMARY %>% mutate(SIGNIF_CLASS = case_when(`-LOG(P-value)` > 2 &  Difference == 10  ~ "NEW SIGNIFICANT",
                                                               `-LOG(P-value)` > 2 &  Difference > 1 ~ "SIGNIFICANT"))
  
  ### SAVE CSV
  write.csv(TIME_CHECK, file = paste0(FILES_LOC, "summaries/", "TIME_CHECK.csv"), row.names = F, quote = F)
  
  # # split TIME_CHECK df by protein IDs to remove contaminant protein IDs + other invalid accession IDs
  # TIME_CHECK_split <- as.data.frame(str_split_fixed(TIME_CHECK$Protein.IDs, ";", n = Inf))
  # TIME_CHECK_split <- 
  #   # bind data frames
  #   cbind(TIME_CHECK, TIME_CHECK_split) %>% 
  #   # remove the initial Protein.IDs column
  #   select(-Protein.IDs) %>%
  #   # pivot the data frame to pivot all individual protein names into one new Protein.IDs column
  #   pivot_longer(names_to = "COL", values_to = "Protein.IDs", cols = 8:39) %>%
  #   # remove unnecessary column "COL"
  #   select(-"COL") %>%
  #   # filter df to remove empty Protein.IDs, rows missing p-value or Differences
  #   filter(Protein.IDs != "", !is.na(`-LOG(P-value)`), !is.na(Difference),
  #          # filter out CONTAMINANTS and other invalid accession IDs
  #          !str_detect(Protein.IDs, '^CON__*'), !str_detect(Protein.IDs, '^REV__*')) 
}

if (EXTRACT_TAXA == T) {
  
  ### READ DATA FRAMES
  MYD88_TAXA <- fread(paste0(PATH_ACCESSION_LISTS, "MYD88_Protein.IDs.txt"), header = F)
  IRAK4_TAXA <- fread(paste0(PATH_ACCESSION_LISTS, "IRAK4_Protein.IDs.txt"), header = F)
  IRAK1_TAXA <- fread(paste0(PATH_ACCESSION_LISTS, "IRAK1_Protein.IDs.txt"), header = F)
  
  ### EXTRACT TAXONMY INFORMATION AND SAVE IN LISTS
  MYD88_TAXA <- GetNamesTaxa(ProteinAccList = MYD88_TAXA$V1)
  IRAK4_TAXA <- GetNamesTaxa(ProteinAccList = IRAK4_TAXA$V1)
  IRAK1_TAXA <- GetNamesTaxa(ProteinAccList = IRAK1_TAXA$V1)
  
  if (WRITE_TAXA_FILES == T) {
    write.csv(MYD88_TAXA, file = paste0(TAXA_PATH, "MYD88_TAXA.csv"))
    write.csv(IRAK4_TAXA, file = paste0(TAXA_PATH, "IRAK4_TAXA.csv"))
    write.csv(IRAK1_TAXA, file = paste0(TAXA_PATH, "IRAK1_TAXA.csv"))
  }
}

