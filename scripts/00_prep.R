# READING, MANIPULATION AND PREPARATION OF IP-MS DATA

# LOAD LIBRARIES
library(tidyverse)
library(data.table)
library(readxl)

# LOAD FUNCTIONS
source("https://raw.githubusercontent.com/tlobnow/master_thesis/main/functions.R")

FOLDER_PATH  = "/Volumes/TAYLOR-LAB/Fenja/Mass Spec/Mass Spec analysis/IL-1 proteomics/lists to take a look at/MyD88_IRAK4_IRAK1-IPs/"
FOLDER_PATH2 = "/Volumes/TAYLOR-LAB/Fenja/Mass Spec/Mass Spec analysis/IL-1 proteomics/lists to take a look at/KO_IRAK4: IRAK1 MyD88-IPs/"
FOLDER_PATH_local = "~/Documents/Github/master_thesis/raw_xlsx/"
OUT_local = "~/Documents/Github/master_thesis/IP_MS_DF/"

DATA_TAY_CONNECTION <- file.exists("/Volumes/TAYLOR-LAB")

READ_RAW         = T
SIGNIF_ONLY      = T
CORRECT_COLNAMES = T
WRITE_SIGNIF_DFs = T
READ_DFs         = T


# READ IN THE RAW FILES STORED IN DATA-TAY
if (READ_RAW == T) {
  MYD88_min15    <- load_dtConnOrLocal(PATH_DT = FOLDER_PATH, PATH_LOCAL = FOLDER_PATH_local, FILE = "MyD88/21M017_MyD88_15min_matrix34_FDR005_s0=1.xlsx") %>% mutate(ORIGIN = 15, PULLED_PROTEIN = "MYD88", KO = "NONE", Significant = case_when(Significant == "+" ~ T, is.na(Significant) ~ F))
  MYD88_min30    <- load_dtConnOrLocal(PATH_DT = FOLDER_PATH, PATH_LOCAL = FOLDER_PATH_local, FILE = "MyD88/21M005_MyD88_30min_matrix74_FDR005_s0=1.xlsx") %>% mutate(ORIGIN = 30, PULLED_PROTEIN = "MYD88", KO = "NONE", Significant = case_when(Significant == "+" ~ T, is.na(Significant) ~ F))
  MYD88_min60    <- load_dtConnOrLocal(PATH_DT = FOLDER_PATH, PATH_LOCAL = FOLDER_PATH_local, FILE = "MyD88/21M036_MyD88_60min_newunstim_matrix40_FDR005_s0=1.xlsx") %>% mutate(ORIGIN = 60, PULLED_PROTEIN = "MYD88", KO = "NONE", Significant = case_when(Significant == "+" ~ T, is.na(Significant) ~ F))
  ################################################################################
  IRAK4_min15    <- load_dtConnOrLocal(PATH_DT = FOLDER_PATH, PATH_LOCAL = FOLDER_PATH_local, FILE = "IRAK4/21M017_21M014unstim_IRAK4_15min_matrix25_FDR005_s0=1.xlsx") %>% mutate(ORIGIN = 15, PULLED_PROTEIN = "IRAK4", KO = "NONE", Significant = case_when(Significant == "+" ~ T, is.na(Significant) ~ F))
  IRAK4_min30    <- load_dtConnOrLocal(PATH_DT = FOLDER_PATH, PATH_LOCAL = FOLDER_PATH_local, FILE = "IRAK4/21M017_21M014unstim_IRAK4_30min_matrix25_FDR005_s0=1.xlsx") %>% mutate(ORIGIN = 30, PULLED_PROTEIN = "IRAK4", KO = "NONE", Significant = case_when(Significant == "+" ~ T, is.na(Significant) ~ F))
  IRAK4_min60    <- load_dtConnOrLocal(PATH_DT = FOLDER_PATH, PATH_LOCAL = FOLDER_PATH_local, FILE = "IRAK4/21M036_IRAK4_60min_matrix57_FDR005_s0=1.xlsx") %>% mutate(ORIGIN = 60, PULLED_PROTEIN = "IRAK4", KO = "NONE", Significant = case_when(Significant == "+" ~ T, is.na(Significant) ~ F))
  ################################################################################
  IRAK1_min15    <- load_dtConnOrLocal(PATH_DT = FOLDER_PATH, PATH_LOCAL = FOLDER_PATH_local, FILE = "IRAK1/21M014_IRAK1_15min_matrix36_FDR005_s0=1.xlsx") %>% mutate(ORIGIN = 15, PULLED_PROTEIN = "IRAK1", KO = "NONE", Significant = case_when(Significant == "+" ~ T, is.na(Significant) ~ F))
  IRAK1_min30    <- load_dtConnOrLocal(PATH_DT = FOLDER_PATH, PATH_LOCAL = FOLDER_PATH_local, FILE = "IRAK1/21M014_IRAK1_30min_matrix36_FDR005_s0=1.xlsx") %>% mutate(ORIGIN = 30, PULLED_PROTEIN = "IRAK1", KO = "NONE", Significant = case_when(Significant == "+" ~ T, is.na(Significant) ~ F))
  IRAK1_min60    <- load_dtConnOrLocal(PATH_DT = FOLDER_PATH, PATH_LOCAL = FOLDER_PATH_local, FILE = "IRAK1/21M036_IRAK1_60min_matrix17_FDR005_s0=1.xlsx") %>% mutate(ORIGIN = 60, PULLED_PROTEIN = "IRAK1", KO = "NONE", Significant = case_when(Significant == "+" ~ T, is.na(Significant) ~ F))
  ################################################################################
  KO_IRAK4_min15 <- load_dtConnOrLocal(PATH_DT = FOLDER_PATH2, PATH_LOCAL = FOLDER_PATH_local, FILE = "IRAK4 KO/21M036_MyD88KOIRAK4_15min_matrix24_FDR005_s0=1.xlsx") %>% mutate(ORIGIN = 15, PULLED_PROTEIN = "MYD88", KO = "IRAK4", Significant = case_when(Significant == "+" ~ T, is.na(Significant) ~ F))
  KO_IRAK4_min30 <- load_dtConnOrLocal(PATH_DT = FOLDER_PATH2, PATH_LOCAL = FOLDER_PATH_local, FILE = "IRAK4 KO/21M036_MyD88KOIRAK4_30min_matrix24_FDR005_s0=1.xlsx") %>% mutate(ORIGIN = 30, PULLED_PROTEIN = "MYD88", KO = "IRAK4", Significant = case_when(Significant == "+" ~ T, is.na(Significant) ~ F))
  KO_IRAK4_min60 <- load_dtConnOrLocal(PATH_DT = FOLDER_PATH2, PATH_LOCAL = FOLDER_PATH_local, FILE = "IRAK4 KO/21M036_MyD88KOIRAK4_60min_matrix24_FDR005_s0=1.xlsx") %>% mutate(ORIGIN = 60, PULLED_PROTEIN = "MYD88", KO = "IRAK4", Significant = case_when(Significant == "+" ~ T, is.na(Significant) ~ F))
  ################################################################################
  KO_IRAK1_min15 <- load_dtConnOrLocal(PATH_DT = FOLDER_PATH2, PATH_LOCAL = FOLDER_PATH_local, FILE = "IRAK1 KO/21M036_MyD88KOIRAK1_15min_matrix18_FDR005_s0=1.xlsx") %>% mutate(ORIGIN = 15, PULLED_PROTEIN = "MYD88", KO = "IRAK1", Significant = case_when(Significant == "+" ~ T, is.na(Significant) ~ F))
  KO_IRAK1_min30 <- load_dtConnOrLocal(PATH_DT = FOLDER_PATH2, PATH_LOCAL = FOLDER_PATH_local, FILE = "IRAK1 KO/21M036_MyD88KOIRAK1_30min_matrix18_FDR005_s0=1.xlsx") %>% mutate(ORIGIN = 30, PULLED_PROTEIN = "MYD88", KO = "IRAK1", Significant = case_when(Significant == "+" ~ T, is.na(Significant) ~ F))
  KO_IRAK1_min60 <- load_dtConnOrLocal(PATH_DT = FOLDER_PATH2, PATH_LOCAL = FOLDER_PATH_local, FILE = "IRAK1 KO/21M036_MyD88KOIRAK1_60min_matrix18_FDR005_s0=1.xlsx") %>% mutate(ORIGIN = 60, PULLED_PROTEIN = "MYD88", KO = "IRAK1", Significant = case_when(Significant == "+" ~ T, is.na(Significant) ~ F))
}

# CORRECT_COLNAMES
if (CORRECT_COLNAMES == T) {
  colnames(MYD88_min15) <- gsub(" ", ".", colnames(MYD88_min15))
  colnames(MYD88_min30) <- gsub(" ", ".", colnames(MYD88_min30))
  colnames(MYD88_min60) <- gsub(" ", ".", colnames(MYD88_min60))
  colnames(IRAK4_min15) <- gsub(" ", ".", colnames(IRAK4_min15))
  colnames(IRAK4_min30) <- gsub(" ", ".", colnames(IRAK4_min30))
  colnames(IRAK4_min60) <- gsub(" ", ".", colnames(IRAK4_min60))
  colnames(IRAK1_min15) <- gsub(" ", ".", colnames(IRAK1_min15))
  colnames(IRAK1_min30) <- gsub(" ", ".", colnames(IRAK1_min30))
  colnames(IRAK1_min60) <- gsub(" ", ".", colnames(IRAK1_min60))
  colnames(KO_IRAK1_min15) <- gsub(" ", ".", colnames(KO_IRAK1_min15))
  colnames(KO_IRAK1_min30) <- gsub(" ", ".", colnames(KO_IRAK1_min30))
  colnames(KO_IRAK1_min60) <- gsub(" ", ".", colnames(KO_IRAK1_min60))
  colnames(KO_IRAK4_min15) <- gsub(" ", ".", colnames(KO_IRAK4_min15))
  colnames(KO_IRAK4_min30) <- gsub(" ", ".", colnames(KO_IRAK4_min30))
  colnames(KO_IRAK4_min60) <- gsub(" ", ".", colnames(KO_IRAK4_min60))
}

# filter for significant protein enrichment only
if (SIGNIF_ONLY == T) {
  MYD88_min15_signif <- filter(MYD88_min15, Significant == T)
  MYD88_min30_signif <- filter(MYD88_min30, Significant == T)
  MYD88_min60_signif <- filter(MYD88_min60, Significant == T)
  IRAK4_min15_signif <- filter(IRAK4_min15, Significant == T)
  IRAK4_min30_signif <- filter(IRAK4_min30, Significant == T)
  IRAK4_min60_signif <- filter(IRAK4_min60, Significant == T)
  IRAK1_min15_signif <- filter(IRAK1_min15, Significant == T)
  IRAK1_min30_signif <- filter(IRAK1_min30, Significant == T)
  IRAK1_min60_signif <- filter(IRAK1_min60, Significant == T)
  KO_IRAK4_min15_signif <- filter(KO_IRAK4_min15, Significant == T)
  KO_IRAK4_min30_signif <- filter(KO_IRAK4_min30, Significant == T)
  KO_IRAK4_min60_signif <- filter(KO_IRAK4_min60, Significant == T)
  KO_IRAK1_min15_signif <- filter(KO_IRAK1_min15, Significant == T)
  KO_IRAK1_min30_signif <- filter(KO_IRAK1_min30, Significant == T)
  KO_IRAK1_min60_signif <- filter(KO_IRAK1_min60, Significant == T)
}

if (WRITE_SIGNIF_DFs == T) {
  if (SIGNIF_ONLY == T) {
    write_csv(MYD88_min15_signif, file = paste0(OUT_local, "MYD88_min15_signif.csv"))
    write_csv(MYD88_min30_signif, file = paste0(OUT_local, "MYD88_min30_signif.csv"))
    write_csv(MYD88_min60_signif, file = paste0(OUT_local, "MYD88_min60_signif.csv"))
    write_csv(IRAK4_min15_signif, file = paste0(OUT_local, "IRAK4_min15_signif.csv"))
    write_csv(IRAK4_min30_signif, file = paste0(OUT_local, "IRAK4_min30_signif.csv"))
    write_csv(IRAK4_min60_signif, file = paste0(OUT_local, "IRAK4_min60_signif.csv"))
    write_csv(IRAK1_min15_signif, file = paste0(OUT_local, "IRAK1_min15_signif.csv"))
    write_csv(IRAK1_min30_signif, file = paste0(OUT_local, "IRAK1_min30_signif.csv"))
    write_csv(IRAK1_min60_signif, file = paste0(OUT_local, "IRAK1_min60_signif.csv"))
    write_csv(KO_IRAK4_min15_signif, file = paste0(OUT_local, "KO_IRAK4_min15_signif.csv"))
    write_csv(KO_IRAK4_min30_signif, file = paste0(OUT_local, "KO_IRAK4_min30_signif.csv"))
    write_csv(KO_IRAK4_min60_signif, file = paste0(OUT_local, "KO_IRAK4_min60_signif.csv"))
    write_csv(KO_IRAK1_min15_signif, file = paste0(OUT_local, "KO_IRAK1_min15_signif.csv"))
    write_csv(KO_IRAK1_min30_signif, file = paste0(OUT_local, "KO_IRAK1_min30_signif.csv"))
    write_csv(KO_IRAK1_min60_signif, file = paste0(OUT_local, "KO_IRAK1_min60_signif.csv"))
  } else {
    write_csv(MYD88_min15, file = paste0(OUT_local, "MYD88_min15.csv"))
    write_csv(MYD88_min30, file = paste0(OUT_local, "MYD88_min30.csv"))
    write_csv(MYD88_min60, file = paste0(OUT_local, "MYD88_min60.csv"))
    write_csv(IRAK4_min15, file = paste0(OUT_local, "IRAK4_min15.csv"))
    write_csv(IRAK4_min30, file = paste0(OUT_local, "IRAK4_min30.csv"))
    write_csv(IRAK4_min60, file = paste0(OUT_local, "IRAK4_min60.csv"))
    write_csv(IRAK1_min15, file = paste0(OUT_local, "IRAK1_min15.csv"))
    write_csv(IRAK1_min30, file = paste0(OUT_local, "IRAK1_min30.csv"))
    write_csv(IRAK1_min60, file = paste0(OUT_local, "IRAK1_min60.csv"))
    write_csv(KO_IRAK4_min15, file = paste0(OUT_local, "KO_IRAK4_min15.csv"))
    write_csv(KO_IRAK4_min30, file = paste0(OUT_local, "KO_IRAK4_min30.csv"))
    write_csv(KO_IRAK4_min60, file = paste0(OUT_local, "KO_IRAK4_min60.csv"))
    write_csv(KO_IRAK1_min15, file = paste0(OUT_local, "KO_IRAK1_min15.csv"))
    write_csv(KO_IRAK1_min30, file = paste0(OUT_local, "KO_IRAK1_min30.csv"))
    write_csv(KO_IRAK1_min60, file = paste0(OUT_local, "KO_IRAK1_min60.csv"))
  }
}

if (READ_DFs == T) {
  MYD88_min15_signif <- fread(file = paste0(OUT_local, "MYD88_min15_signif.csv"))
  MYD88_min30_signif <- fread(file = paste0(OUT_local, "MYD88_min30_signif.csv"))
  MYD88_min60_signif <- fread(file = paste0(OUT_local, "MYD88_min60_signif.csv"))
  IRAK4_min15_signif <- fread(file = paste0(OUT_local, "IRAK4_min15_signif.csv"))
  IRAK4_min30_signif <- fread(file = paste0(OUT_local, "IRAK4_min30_signif.csv"))
  IRAK4_min60_signif <- fread(file = paste0(OUT_local, "IRAK4_min60_signif.csv"))
  IRAK1_min15_signif <- fread(file = paste0(OUT_local, "IRAK1_min15_signif.csv"))
  IRAK1_min30_signif <- fread(file = paste0(OUT_local, "IRAK1_min30_signif.csv"))
  IRAK1_min60_signif <- fread(file = paste0(OUT_local, "IRAK1_min60_signif.csv"))
  KO_IRAK4_min15_signif <- fread(file = paste0(OUT_local, "KO_IRAK4_min15_signif.csv"))
  KO_IRAK4_min30_signif <- fread(file = paste0(OUT_local, "KO_IRAK4_min30_signif.csv"))
  KO_IRAK4_min60_signif <- fread(file = paste0(OUT_local, "KO_IRAK4_min60_signif.csv"))
  KO_IRAK1_min15_signif <- fread(file = paste0(OUT_local, "KO_IRAK1_min15_signif.csv"))
  KO_IRAK1_min30_signif <- fread(file = paste0(OUT_local, "KO_IRAK1_min30_signif.csv"))
  KO_IRAK1_min60_signif <- fread(file = paste0(OUT_local, "KO_IRAK1_min60_signif.csv"))
  
  MYD88_min15 <- fread(file = paste0(OUT_local, "MYD88_min15.csv"))
  MYD88_min30 <- fread(file = paste0(OUT_local, "MYD88_min30.csv"))
  MYD88_min60 <- fread(file = paste0(OUT_local, "MYD88_min60.csv"))
  IRAK4_min15 <- fread(file = paste0(OUT_local, "IRAK4_min15.csv"))
  IRAK4_min30 <- fread(file = paste0(OUT_local, "IRAK4_min30.csv"))
  IRAK4_min60 <- fread(file = paste0(OUT_local, "IRAK4_min60.csv"))
  IRAK1_min15 <- fread(file = paste0(OUT_local, "IRAK1_min15.csv"))
  IRAK1_min30 <- fread(file = paste0(OUT_local, "IRAK1_min30.csv"))
  IRAK1_min60 <- fread(file = paste0(OUT_local, "IRAK1_min60.csv"))
  KO_IRAK4_min15 <- fread(file = paste0(OUT_local, "KO_IRAK4_min15.csv"))
  KO_IRAK4_min30 <- fread(file = paste0(OUT_local, "KO_IRAK4_min30.csv"))
  KO_IRAK4_min60 <- fread(file = paste0(OUT_local, "KO_IRAK4_min60.csv"))
  KO_IRAK1_min15 <- fread(file = paste0(OUT_local, "KO_IRAK1_min15.csv"))
  KO_IRAK1_min30 <- fread(file = paste0(OUT_local, "KO_IRAK1_min30.csv"))
  KO_IRAK1_min60 <- fread(file = paste0(OUT_local, "KO_IRAK1_min60.csv")) 
}
