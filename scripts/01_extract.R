#         .
#       ":"
#     ___:____     |"\/"|
#   ,'        `.    \  /
#   |  O        \___/  |
# ~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^

### EXTRACT THE PROTEINS INTERESTING FOR COMPLEX PREDICTION

# LOAD LIBRARIES
library(tidyverse)
library(data.table)
library(readxl)
library(UniprotR)

# LOAD FUNCTIONS
source(file = ifelse(exists("https://raw.githubusercontent.com/tlobnow/master_thesis/main/scripts/functions.R"), 
                     yes =  "https://raw.githubusercontent.com/tlobnow/master_thesis/main/scripts/functions.R",
                     no  =  "~/Documents/Github/master_thesis/scripts/functions.R"))

### EXTRACTING INFO FROM FENJA'S IMMUNOPRECIPITATION-MASS SPECTROMETRY (IP-MS) DATA

# DEFINE PATHS
FILES_LOC             = "~/Documents/Github/master_thesis/"
FOLDER_PATH           = "/Volumes/TAYLOR-LAB/Fenja/Mass Spec/Mass Spec analysis/IL-1 proteomics/lists to take a look at/MyD88_IRAK4_IRAK1-IPs/"
FOLDER_PATH2          = "/Volumes/TAYLOR-LAB/Fenja/Mass Spec/Mass Spec analysis/IL-1 proteomics/lists to take a look at/KO_IRAK4: IRAK1 MyD88-IPs/"
DATA_TAY_CONNECTION   = file.exists("/Volumes/TAYLOR-LAB")
FOLDER_PATH_LOCAL     = paste0(FILES_LOC, "raw_csv/")
PATH_ACCESSION_LISTS  = paste0(FILES_LOC, "uniprot_accession_lists/")
TAXA_PATH             = paste0(FILES_LOC, "taxa_lists/")


# READ PROTEINS THAT WERE SIGNIFICANTLY ENRICHED UPON STIMULATION (MISSING IN UNSTIMULATED)
MYD88_min15_signif    <- fread(file = paste0(FOLDER_PATH_LOCAL, "MYD88_min15_signif.csv"))
MYD88_min30_signif    <- fread(file = paste0(FOLDER_PATH_LOCAL, "MYD88_min30_signif.csv"))
MYD88_min60_signif    <- fread(file = paste0(FOLDER_PATH_LOCAL, "MYD88_min60_signif.csv"))
IRAK4_min15_signif    <- fread(file = paste0(FOLDER_PATH_LOCAL, "IRAK4_min15_signif.csv"))
IRAK4_min30_signif    <- fread(file = paste0(FOLDER_PATH_LOCAL, "IRAK4_min30_signif.csv"))
IRAK4_min60_signif    <- fread(file = paste0(FOLDER_PATH_LOCAL, "IRAK4_min60_signif.csv"))
IRAK1_min15_signif    <- fread(file = paste0(FOLDER_PATH_LOCAL, "IRAK1_min15_signif.csv"))
IRAK1_min30_signif    <- fread(file = paste0(FOLDER_PATH_LOCAL, "IRAK1_min30_signif.csv"))
IRAK1_min60_signif    <- fread(file = paste0(FOLDER_PATH_LOCAL, "IRAK1_min60_signif.csv"))
KO_IRAK4_min15_signif <- fread(file = paste0(FOLDER_PATH_LOCAL, "KO_IRAK4_min15_signif.csv"))
KO_IRAK4_min30_signif <- fread(file = paste0(FOLDER_PATH_LOCAL, "KO_IRAK4_min30_signif.csv"))
KO_IRAK4_min60_signif <- fread(file = paste0(FOLDER_PATH_LOCAL, "KO_IRAK4_min60_signif.csv"))
KO_IRAK1_min15_signif <- fread(file = paste0(FOLDER_PATH_LOCAL, "KO_IRAK1_min15_signif.csv"))
KO_IRAK1_min30_signif <- fread(file = paste0(FOLDER_PATH_LOCAL, "KO_IRAK1_min30_signif.csv"))
KO_IRAK1_min60_signif <- fread(file = paste0(FOLDER_PATH_LOCAL, "KO_IRAK1_min60_signif.csv"))

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


### JOIN DATA FRAMES FOR DIFFERENT TIME POINTS & RETRIEVE UNIPROT ACCESSION IDs
### MYD88
    # JOIN MYD88 SIGNIFICANT PROTEINS
    MYD88_signif <- join_timepoints(DF1 = MYD88_min15_signif, DF2 = MYD88_min30_signif, DF3 = MYD88_min60_signif)
    retrieveAccessionIDs(DF = MYD88_signif, OUT = paste0(PATH_ACCESSION_LISTS, "MYD88_Protein.IDs_signif.txt"))
    
    # JOIN MYD88 NEW SIGNIFICANT PROTEINS (undetected in unstimulated)
    MYD88_new_signif <- MYD88_signif %>% filter(Difference == 10)
    retrieveAccessionIDs(DF  = MYD88_new_signif, OUT = paste0(PATH_ACCESSION_LISTS, "MYD88_Protein.IDs_new_signif.txt"))
    
    # JOIN MYD88 ALL
    MYD88 <- join_timepoints(DF1 = MYD88_min15, DF2 = MYD88_min30, DF3 = MYD88_min60) %>% unique()
    retrieveAccessionIDs(DF  = MYD88, OUT = paste0(PATH_ACCESSION_LISTS, "MYD88_Protein.IDs.txt"))
    # many filtered out due to pulldown at multiple time points [we only want the accession IDs, so we don't care about that here]
    
### KO
    # JOIN IRAK4 KO ALL
    KO_IRAK4 <- join_timepoints(DF1 = KO_IRAK4_min15, DF2 = KO_IRAK4_min30, DF3 = KO_IRAK4_min60)
    retrieveAccessionIDs(DF  = KO_IRAK4, OUT = paste0(PATH_ACCESSION_LISTS, "KO_IRAK4_Protein.IDs.txt"))
    
    # JOIN IRAK1 KO ALL
    KO_IRAK1 <- join_timepoints(DF1 = KO_IRAK1_min15, DF2 = KO_IRAK1_min30, DF3 = KO_IRAK1_min60)
    colnames(KO_IRAK1)[colnames(KO_IRAK1)%in%"T:.Protein.IDs"]          <- "Protein.IDs"
    colnames(KO_IRAK1)[colnames(KO_IRAK1)%in%"N:.Valid.values"]         <- "Valid.values"
    colnames(KO_IRAK1)[colnames(KO_IRAK1)%in%"T:.Majority.protein.IDs"] <- "Majority.protein.IDs"
    colnames(KO_IRAK1)[colnames(KO_IRAK1)%in%"T:.Gene.names"]           <- "Gene.names"
    retrieveAccessionIDs(DF  = KO_IRAK1, OUT = paste0(PATH_ACCESSION_LISTS, "KO_IRAK1_Protein.IDs.txt"))

### IRAK4
      # JOIN IRAK4 SIGNIFICANT PROTEINS
      IRAK4_signif <- join_timepoints(DF1 = IRAK4_min15_signif, DF2 = IRAK4_min30_signif, DF3 = IRAK4_min60_signif)
      retrieveAccessionIDs(DF  = IRAK4_signif, OUT = paste0(PATH_ACCESSION_LISTS, "IRAK4_Protein.IDs_signif.txt"))
      
      # JOIN IRAK4 NEW SIGNIFICANT PROTEINS (undetected in unstimulated)
      IRAK4_new_signif <- IRAK4_signif %>% filter(Difference == 10)
      retrieveAccessionIDs(DF  = IRAK4_new_signif, OUT = paste0(PATH_ACCESSION_LISTS, "IRAK4_Protein.IDs_new_signif.txt"))
      
      # JOIN IRAK4 ALL
      IRAK4 <- join_timepoints(DF1 = IRAK4_min15, DF2 = IRAK4_min30, DF3 = IRAK4_min60)
      retrieveAccessionIDs(DF  = IRAK4, OUT = paste0(PATH_ACCESSION_LISTS, "IRAK4_Protein.IDs.txt"))

### IRAK 1
      # JOIN IRAK1 SIGNIFICANT PROTEINS
      IRAK1_signif <- join_timepoints(DF1 = IRAK1_min15_signif, DF2 = IRAK1_min30_signif, DF3 = IRAK1_min60_signif)
      colnames(IRAK1_signif)[colnames(IRAK1_signif)%in%"T:.Protein.IDs"]          <- "Protein.IDs"
      colnames(IRAK1_signif)[colnames(IRAK1_signif)%in%"N:.Valid.values"]         <- "Valid.values"
      colnames(IRAK1_signif)[colnames(IRAK1_signif)%in%"T:.Majority.protein.IDs"] <- "Majority.protein.IDs"
      colnames(IRAK1_signif)[colnames(IRAK1_signif)%in%"T:.Gene.names"]           <- "Gene.names"
      retrieveAccessionIDs(DF  = IRAK1_signif, OUT = paste0(PATH_ACCESSION_LISTS, "IRAK1_Protein.IDs_signif.txt"))
      
      # JOIN IRAK1 NEW SIGNIFICANT PROTEINS (undetected in unstimulated)
      IRAK1_new_signif <- IRAK1_signif %>% filter(Difference == 10)
      retrieveAccessionIDs(DF  = IRAK1_new_signif, OUT = paste0(PATH_ACCESSION_LISTS, "IRAK1_Protein.IDs_new_signif.txt"))
      
      # JOIN IRAK1 ALL
      IRAK1 <- join_timepoints(DF1 = IRAK1_min15, DF2 = IRAK1_min30, DF3 = IRAK1_min60)
      colnames(IRAK1)[colnames(IRAK1)%in%"T:.Protein.IDs"]          <- "Protein.IDs"
      colnames(IRAK1)[colnames(IRAK1)%in%"N:.Valid.values"]         <- "Valid.values"
      colnames(IRAK1)[colnames(IRAK1)%in%"T:.Majority.protein.IDs"] <- "Majority.protein.IDs"
      colnames(IRAK1)[colnames(IRAK1)%in%"T:.Gene.names"]           <- "Gene.names"
      retrieveAccessionIDs(DF  = IRAK1, OUT = paste0(PATH_ACCESSION_LISTS, "IRAK1_Protein.IDs.txt"))


### EXTRACT TAXONMY INFORMATION AND SAVE IN LISTS
MYD88_TAXA <- fread(paste0(PATH_ACCESSION_LISTS, "MYD88_Protein.IDs.txt"), header = F)
MYD88_TAXA <- GetNamesTaxa(ProteinAccList = MYD88_TAXA$V1)
write.csv(MYD88_TAXA, file = paste0(TAXA_PATH, "MYD88_TAXA.csv"))

IRAK4_TAXA <- fread(paste0(PATH_ACCESSION_LISTS, "IRAK4_Protein.IDs.txt"), header = F)
IRAK4_TAXA <- GetNamesTaxa(ProteinAccList = IRAK4_TAXA$V1)
write.csv(IRAK4_TAXA, file = paste0(TAXA_PATH, "IRAK4_TAXA.csv"))

IRAK1_TAXA <- fread(paste0(PATH_ACCESSION_LISTS, "IRAK1_Protein.IDs.txt"), header = F)
IRAK1_TAXA <- GetNamesTaxa(ProteinAccList = IRAK1_TAXA$V1)
write.csv(IRAK1_TAXA, file = paste0(TAXA_PATH, "IRAK1_TAXA.csv"))
