library(dplyr)
library(tidyr)
library(stringr)
library(fs)
library(jsonlite)
library(purrr)
library(utils)
library(data.table)

MAIN = paste0("/Volumes/TAYLOR-LAB/Finn/RESULTS/ARL8B_UN93B/ARL8B_UN93B_ALL/")
MAIN = paste0("/Users/u_lobnow/Documents/Github/transferGit/")
MAIN = paste0("/Volumes/TAYLOR-LAB/Finn/RESULTS/CHIMY_T6BM/")

FILES <- list.files(MAIN)

# if you wish to extract a single file, unhash below and provide file name
#FILES = "ARL8B_MOUSE_x1_UN93B_MOUSE_x1_TLR7_MOUSE_x1"
#FILE = "ARL8B_MOUSE_x1_UN93B_MOUSE_x1_TLR7_MOUSE_x1"


jsonExtract <- function(JSON, OUT, FILE) {
  json       <- fromJSON(JSON)
  MODEL      <- strsplit(json$order[[1]], "_", fixed = TRUE)[[1]][2]
  TOL        <- json$tol_values %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "TOL") 
  pLDDT      <- json$plddts %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "pLDDT") 
  pTM        <- json$ptms %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "pTM") 
  piTM       <- json$pitms %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "piTM") 
  iScore     <- json$`interface score` %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "iScore") 
  iRes       <- json$`interfacial residue number` %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "iRes") 
  iCnt       <- json$`interficial contact number` %>% as.data.frame() %>% pivot_longer(everything(), names_to = "RECYCLE", values_to = "iCnt") 
  FILE_MODEL <- paste(FILE, "MODEL", MODEL, sep = "_")
  NUM_CLUSTERS <- json$clusters[[iScore$RECYCLE[1]]]$num_clusters
  N_MONOMERS <- length(json$chains)
  
  EXTRACT <- cbind(FILE, MODEL, TOL, pLDDT, pTM, piTM, iScore, iRes, iCnt, FILE_MODEL, NUM_CLUSTERS, N_MONOMERS)
  EXTRACT <- EXTRACT[, !duplicated(colnames(EXTRACT))]
  
  write.table(EXTRACT, file = paste0(OUT,"_withRecycles.csv"),sep = ",", append = T, quote = F, row.names = F, col.names = F)
  
  EXTRACT_noRecycle <- EXTRACT %>% filter(!str_detect(RECYCLE, "_recycled_"))
  write.table(EXTRACT_noRecycle, file = paste0(OUT,".csv"),sep = ",", append = T, quote = F, row.names = F, col.names = F)
}

for (FILE in FILES) {
  FOLDER = paste0(MAIN, FILE, "/")
  dir_create(FOLDER,"CSV")
  csv_file <- paste0(FOLDER, "CSV/", FILE, "_withRecycles.csv")
  csv_file2 <- paste0(FOLDER, "CSV/", FILE, ".csv")
  
  if (file.exists(csv_file))  {file.remove(csv_file)}
  if (file.exists(csv_file2)) {file.remove(csv_file2)}
  
  json_folder <- ifelse(dir.exists(file.path(FOLDER, "JSON")), yes = file.path(FOLDER, "JSON"), file.path(FOLDER))
  maxJSON     <- list.files(json_folder, pattern = ".json")
   if (is_empty(maxJSON)) {
     next
   }
  
  json_files <- dir_ls(json_folder, regexp = "\\.json$", recurse = TRUE)# %>% str_subset("Infinity")
  
  for (file in json_files) {
    json <- readLines(file)
    json <- str_replace_all(json, "Infinity", "9999")
    writeLines(json, file)
  }
  
  if (length(maxJSON) > 0) {
    for (i in 1:length(maxJSON)) {
      JSON <- file.path(json_folder, maxJSON[i])
      OUT  <- paste(FOLDER, "CSV", FILE, sep = "/")
      jsonExtract(JSON = JSON, OUT = OUT, FILE = FILE)
    }
  }
  
  CSV_FILES <- c(csv_file, csv_file2)
  
  for (CSV_FILE in CSV_FILES) {
    JSON_EXTRACT <- data.table::fread(CSV_FILE, header = FALSE) %>%
      dplyr::mutate(ORIGIN = FILE) %>%
      dplyr::rename(FILE = V1, MODEL = V2, RECYCLE = V3, TOL = V4, pLDDT = V5, pTM = V6, piTM = V7, iScore = V8, iRes = V9, iCnt = V10, FILE_MODEL = V11, NUM_CLUSTERS = V12, N_MONOMERS = V13)
    
    JE <- JSON_EXTRACT %>%
      dplyr::mutate(FILE_RECYCLE = paste0(FILE_MODEL, "_RECYCLE_", RECYCLE), RANK = NA) %>%
      dplyr::distinct(FILE_RECYCLE, .keep_all = TRUE) %>%
      dplyr::group_by(FILE) %>%
      dplyr::mutate(RANK = frank(desc(iScore), ties.method = "min"))
    
    data.table::fwrite(JE, CSV_FILE, row.names = FALSE)
    
    data.table::fwrite(JE, paste0(MAIN, "BigBoySummary.csv"), row.names = FALSE, append = T)
  }
}

bigboy <- fread(paste0(MAIN, "BigBoySummary.csv")) %>% unique()
data.table::fwrite(bigboy, paste0(MAIN, "BigBoySummary.csv"), row.names = FALSE, append = T)

bigboy_noRecycle <- bigboy %>% filter(!str_detect(RECYCLE, "_recycled_"))


#JE_all <- data.table::fread(csv_file)
#JE <- data.table::fread(csv_file2)
