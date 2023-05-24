# # remove pre-existing csv file, append would lead to duplicate rows
# system(command = paste0(" [ -f ", RAW_SUMMARIES_PATH, FOLDER ,"_fromJSON.csv ] && rm ", RAW_SUMMARIES_PATH, FOLDER, "_fromJSON.csv"))
# 
# LOCATION = paste0(MAIN, FOLDER)
# LIST = list.files(LOCATION)
# 
# for (FILE in LIST) {
#   print(paste0("processing ", FILE))
#   maxJSON=list.files(paste0(LOCATION, "/" , FILE, "/JSON/"))
#   
#   # replace all "Infinity" strings with large number (9999)
#   system(command = paste0("grep -rl Infinity ", MAIN, FOLDER,"/", FILE, "/JSON/", " | xargs sed -i '' -e 's/Infinity/9999/g'"))
#   
#   if (length(maxJSON) > 0) {
#     for (i in 1:length(maxJSON)) {
#       JSON  = paste0(LOCATION, "/" , FILE, "/JSON/", maxJSON[i])
#       OUT   = paste0(RAW_SUMMARIES_PATH, FOLDER)
#       jsonExtract(JSON = JSON, OUT = OUT, FILE = FILE)
#     } 
#   } else {
#     next
#     print(paste0("skipped", FILE))
#   }
# }

# Remove pre-existing CSV file, append would lead to duplicate rows
csv_file <- paste0(RAW_SUMMARIES_PATH, FOLDER, ".csv")
if (file.exists(csv_file)) {
  file.remove(csv_file)
}

csv_file2 <- paste0(RAW_SUMMARIES_PATH, FOLDER, "_noRecycle.csv")
if (file.exists(csv_file2)) {
  file.remove(csv_file2)
}

LOCATION <- file.path(MAIN, FOLDER)
LIST     <- list.files(LOCATION)

for (FILE in LIST) {
  print(paste("Processing", FILE))
  json_folder <- file.path(LOCATION, FILE, "JSON")
  maxJSON     <- list.files(json_folder)
  if (is_empty(maxJSON)) {
    next
  }
  
  # Replace all "Infinity" strings with large number (9999)
  # Get a list of JSON files with "Infinity" strings
  json_files <- dir_ls(json_folder, regexp = "\\.json$", recurse = TRUE) %>%
    str_subset("Infinity")
  # Replace "Infinity" with 9999 in each JSON file
  for (file in json_files) {
    json <- readLines(file)
    json <- str_replace_all(json, "Infinity", "9999")
    writeLines(json, file)
  }

  if (length(maxJSON) > 0) {
    for (i in 1:length(maxJSON)) {
      JSON <- file.path(json_folder, maxJSON[i])
      OUT  <- paste0(RAW_SUMMARIES_PATH, FOLDER)
      jsonExtract(JSON = JSON, OUT = OUT, FILE = FILE)
    } 
  } else {
    print(paste("Skipped", FILE))
  }
}