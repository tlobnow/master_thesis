# READ FILE
# JSON_EXTRACT  <- read.csv(paste0(RAW_SUMMARIES_PATH, FOLDER, "_fromJSON.csv"), header = F) %>% 
JSON_EXTRACT  <- read.csv(paste0(RAW_SUMMARIES_PATH, FOLDER, ".csv"), header = F) %>% 
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