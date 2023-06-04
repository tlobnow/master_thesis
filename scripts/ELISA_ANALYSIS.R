#          .                             
#       ":"                               
#     ___:____     |"\/"|               
#   ,'        `.    \  /                
#   |  O        \___/  |               
# ~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^

#Load libraries
library(pacman)
pacman::p_load(data.table, tidyr, dplyr, ggplot2, interleave, ggpubr)

#Specify Input Folder (where plates are located)
Input_Directory <- "/Volumes/TAYLOR-LAB/Finn/ELISA/20230529_ELISA/"
Input_Directory <- "/Volumes/TAYLOR-LAB/Finn/ELISA/20230512_ELISA/"


CELL_LINES    <- "CELL_LINES.csv" # COHORT
STIM_DAYS     <- "STIMULATION_DAYS.csv" # SAMPLE_DAY
CONDITIONS    <- "CONDITIONS.csv" # STIMULATION_CONDITION
MEASUREMENTS  <- "MEASUREMENTS.csv" # VALUES_MEASURED


# each plate folder should contain the following files:
  # Cohort.csv                  -- Cell line  per well (WT, MyD88 KO, etc.)
  # Sample_Day.csv              -- Stimulation Day per well (Day_1, Day_2, Day_3)
  # Stimulation_Condition.csv   -- Treatment per well (Stimulated, Unstimulated, Calibration, NA)
  # Values_measured.csv         -- Plate measurements (machine)


#Creating Output Folder if missing
Output_Directory <- file.path(Input_Directory, "Output")
if(!file.exists(Output_Directory)){
  dir.create(Output_Directory)
}

#Specify Dilution factor of supernatant
# DILUTION_FACTOR <- 5
DILUTION_FACTOR <- 10

#Calculating Number of plates in Input Folder
for (x in 1:length(list.files(Input_Directory, pattern = "Plate_"))){
  if(file.exists(paste0(Input_Directory, "/Plate_", x))){
    Number_of_plates <- x
  }
}

ELISA_Fx <- function(plate_number){
  
  # Reading Data ------------------------------------------------------------
  #Getting Path of treatment conditions and values
  Input_plate <- paste0(Input_Directory, "/Plate_", plate_number)
  
  #Reading Plate Treatement 
  MEASUREMENTS <- fread(paste0(Input_plate, "/MEASUREMENTS.csv"), header = F)
  CELL_LINES <- fread(paste0(Input_plate, "/CELL_LINES.csv"), header = F)
  CONDITIONS <- fread(paste0(Input_plate, "/CONDITIONS.csv"), header = F)
  STIM_DAYS <- fread(paste0(Input_plate, "/STIM_DAYS.csv"), header = F)
  
  
  #Converting tables into vector for to make a single table
  MEASUREMENTS <- as.vector(as.matrix(MEASUREMENTS))
  CELL_LINES <- as.vector(as.matrix(CELL_LINES))
  CONDITIONS <- as.vector(as.matrix(CONDITIONS))
  STIM_DAYS <- as.vector(as.matrix(STIM_DAYS))
  
  #Creating Table containing all plate Information
  Plate <- NULL
  Plate$MEASUREMENTS <- MEASUREMENTS
  Plate$CELL_LINES <- CELL_LINES
  Plate$CONDITIONS <- CONDITIONS
  Plate$STIM_DAYS <- STIM_DAYS
  
  rm(
    MEASUREMENTS,
    CELL_LINES,
    CONDITIONS,
    STIM_DAYS
  )
  
  Plate <- Plate %>% as.data.table()
  
  #Removing Empty Wells
  Plate <- Plate %>% 
    filter(
      CELL_LINES != "BLANK"
    ) %>% as.data.table()
  
  
  # Standard Curve ----------------------------------------------------------
  #Creating a Standard Curve
  Plate_Standards <- Plate %>% 
    filter(
      CONDITIONS == "CALIBRATION"
    ) %>% 
    group_by(
      CELL_LINES
    ) %>% 
    summarise(
      MEASUREMENTS_mean = mean(MEASUREMENTS)
      # MEASUREMENTS_median = median(MEASUREMENTS)
    ) %>%  
    mutate(
      CELL_LINES = as.numeric(CELL_LINES)
    ) %>% 
    arrange(
      CELL_LINES
    )
  
  Fit <- lm(CELL_LINES ~ MEASUREMENTS_mean -1, data = Plate_Standards) #linear model of the Standard curve. -1 omits the intercept
  
  R <- summary(Fit)$r.squared
  
  Rsquare <- signif(R, digits = 4)
  
  rm(
    R
  )
  
  print(paste0("IL2-Amount = slope*Intensity"))
  print(paste0("IL2-Amount = ", Fit$coefficients[1],"*Intensity"))
  
  Plate_Standards <- 
    Plate_Standards %>% 
    mutate(
      Fit_Test = (Fit$coefficients[1]*MEASUREMENTS_mean)
    )
  
  ggplot(
    data = Plate_Standards,
  ) +
    geom_point(
      aes(
        x = MEASUREMENTS_mean,
        y = CELL_LINES
      ),
      size = 5
    ) +
    geom_line(
      aes(
        x = MEASUREMENTS_mean,
        y = Fit_Test
      ),
      linetype = "dashed"
    ) +
    annotate(
      'text',
      x = 0.15,
      y = 700,
      label = paste0("R^2 = ",Rsquare),
      size = 10
    ) +
    annotate(
      'text',
      x = max(Plate_Standards$MEASUREMENTS_mean) - (0.25*max(Plate_Standards$MEASUREMENTS_mean)),
      y = 150,
      label = paste0("IL-Amount = \n", signif(Fit$coefficients[1], digits = 4),"*Intensity")
    ) +
    labs(
      x = "Measured Values",
      y = "IL-Concentration (pg/mL)"
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 20)
    )
  
  Save_Name <- paste0("Plate_Number_" , plate_number, "_Standard_Curve.pdf")
  Save_Name <- file.path(Output_Directory, Save_Name)
  
  ggsave(
    Save_Name,
    plot = last_plot(),
    height = 3*3,
    width = 5*4
  )
  
  rm(
    Save_Name,
    Plate_Standards,
    Rsquare
  )
  
  
  # Fitting Data To Standarad Curve -----------------------------------------
  Plate <- Plate %>% 
    filter(
      CONDITIONS != "CALIBRATION"
    ) %>% 
    mutate(
      MEASUREMENTS = as.numeric(MEASUREMENTS),
      IL2_concentration = (Fit$coefficients[1]*MEASUREMENTS),
      IL2_concentration_DILUTION_FACTOR = IL2_concentration*DILUTION_FACTOR
    )
  
  rm(
    Fit
  )
  
  return(Plate)
}

All_plates_data = data.frame()

for (plate_number in 1:Number_of_plates){
  
  
  data_temp <- ELISA_Fx(plate_number)
  data_temp$Plate <- plate_number
  
  All_plates_data <- rbind(All_plates_data, data_temp)
  
  rm(data_temp)
}

# Graph with real pg/ml ---------------------------------------------------
# Getting mean of each column on a day to day basis
Plate_Summary_Day <- 
  All_plates_data %>% 
  group_by(
    CELL_LINES,
    CONDITIONS,
    STIM_DAYS
  ) %>% 
  summarise(
    IL2_concentration_DILUTION_FACTOR_mean = mean(IL2_concentration_DILUTION_FACTOR),
    IL2_concentration_DILUTION_FACTOR_median = median(IL2_concentration_DILUTION_FACTOR)
  ) %>% 
  as.data.table()

#Getting overall value of a single CELL_LINES
Plate_Summary_CELL_LINES <- 
  All_plates_data %>% 
  group_by(
    CELL_LINES,
    CONDITIONS
  ) %>% 
  summarise(
    IL2_concentration_DILUTION_FACTOR_mean = mean(IL2_concentration_DILUTION_FACTOR),
    IL2_concentration_DILUTION_FACTOR_median = median(IL2_concentration_DILUTION_FACTOR),
    IL2_concentration_DILUTION_FACTOR_sd = sd(IL2_concentration_DILUTION_FACTOR)
  ) %>% 
  as.data.table()

color_elisa <- c("Unstimulated" = "white",
                 "Stimulated" = "grey")

ggplot(
  data = Plate_Summary_CELL_LINES,
  aes(
    x = CELL_LINES,
    y = IL2_concentration_DILUTION_FACTOR_mean,
    fill = CONDITIONS
  )
) +
  geom_col(
    position = position_dodge(width = 1),
    color = "black"
  ) +
  geom_errorbar(
    data = Plate_Summary_CELL_LINES,
    aes(
      x = CELL_LINES,
      ymin = IL2_concentration_DILUTION_FACTOR_mean - IL2_concentration_DILUTION_FACTOR_sd,
      ymax = IL2_concentration_DILUTION_FACTOR_mean + IL2_concentration_DILUTION_FACTOR_sd
    ),
    linewidth = .75,
    position = position_dodge(width = 1),
    width = 0.5
  ) +
  geom_point(
    data = Plate_Summary_Day,
    aes(
      x = CELL_LINES,
      y = IL2_concentration_DILUTION_FACTOR_mean
    ),
    size = 2,
    position = position_dodge(width = 1)
  )+
  fill_palette(
    palette = color_elisa
  ) +
  labs(
    y = "IL-2 conc. (pg/mL)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust=1, size = 12),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 20)
  )

Plot_Save <- file.path(Output_Directory, "01_Real_IL2_pg_per_ml_conc.pdf")
ggsave(
  Plot_Save,
  plot = last_plot(),
  height = 3*3,
  width = 5*4
)


# Relative IL2 secretion --------------------------------------------------

# Getting mean of each column on a day to day basis
Plate_Summary_Day_Rel <- 
  All_plates_data %>% 
  group_by(
    CELL_LINES,
    CONDITIONS,
    STIM_DAYS
  ) %>% 
  summarise(
    MEASUREMENTS_mean = mean(MEASUREMENTS),
    MEASUREMENTS_median = median(MEASUREMENTS)
  ) %>% 
  ungroup(
    CELL_LINES,
    CONDITIONS
  ) %>% 
  group_by(
    STIM_DAYS
  ) %>% 
  mutate(
    Relative_IL2_concentration = (MEASUREMENTS_mean-min(MEASUREMENTS_mean))/(max(MEASUREMENTS_mean)-min(MEASUREMENTS_mean))
  ) %>% 
  as.data.table()

#Getting overall value of a single cohort (cell line)
Plate_Summary_CELL_LINES_Rel <- 
  Plate_Summary_Day_Rel %>% 
  group_by(
    CELL_LINES,
    CONDITIONS
  ) %>% 
  summarise(
    Relative_IL2_concentration_mean = mean(Relative_IL2_concentration),
    Relative_IL2_concentration_median = median(Relative_IL2_concentration),
    Relative_IL2_concentration_sd = sd(Relative_IL2_concentration)
  ) %>% 
  as.data.table()

ggplot(
  data = Plate_Summary_CELL_LINES_Rel,
  aes(
    x = CELL_LINES,
    y = Relative_IL2_concentration_mean,
    fill = CONDITIONS
  )
) +
  geom_col(
    position = position_dodge(width = 1),
    color = "black"
  ) +
  geom_errorbar(
    data = Plate_Summary_CELL_LINES_Rel,
    aes(
      x = CELL_LINES,
      ymin = Relative_IL2_concentration_mean - Relative_IL2_concentration_sd,
      ymax = Relative_IL2_concentration_mean + Relative_IL2_concentration_sd
    ),
    linewidth = .75,
    position = position_dodge(width = 1),
    width = 0.5
  ) +
  geom_point(
    data = Plate_Summary_Day_Rel,
    aes(
      x = CELL_LINES,
      y = Relative_IL2_concentration
    ),
    size = 2,
    position = position_jitterdodge(dodge.width = 1)
  )+
  fill_palette(
    palette = color_elisa
  ) +
  labs(
    y = "relative IL-2 conc."
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust=1, size = 12),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 20)
  )

Plot_Save <- file.path(Output_Directory, "01_Relative_IL2_conc.pdf")

ggsave(
  Plot_Save,
  plot = last_plot(),
  height = 3*3,
  width = 5*4
)















