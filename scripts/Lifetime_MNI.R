# Load libraries
# install.packages("pacman")
library(pacman)

pacman::p_load(dplyr, ggplot2, viridis, data.table, dtplyr, ggpubr, viridis, tidyr, purrr, R.utils, ggbeeswarm, RColorBrewer)

#load tables
TablePaths <- c(#"/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/chMyD88-DHF91-TRAF6-BD-GFP/20220401 6nM_cl221_chMyD88-DHF91-TRAF6-BD-GFP_002/Essential.csv.gz",
  #"/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/chMyD88-DHF91-TRAF6-BD-GFP/20220401 6nM_cl222_chMyD88-DHF91-TRAF6-BD-GFP_001/Essential.csv.gz",
  #"/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/chMyD88-DHF91-TRAF6-BD-GFP/20220401 6nM_cl222_chMyD88-DHF91-TRAF6-BD-GFP_002/Essential.csv.gz",
  #"/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/DHF58v2mid_DHF58v.2/20220401 6nM_DHF58v2mid_DHF58v.2_001/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-TRAF6-BD TRAF6/20220516 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 001/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-TRAF6-BD TRAF6/20220610 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 001/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-TRAF6-BD TRAF6/20220615 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 002/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-TRAF6-BD-3xA TRAF6/20220615 1.5nM_cl234_TRAF6_MyD88-TRAF6-BD-3xA-GFP 008/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD TRAF6/20220516 1.5nM_cl236_TRAF6_chMyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD TRAF6/20220516 1.5nM_cl236_TRAF6_chMyD88-DHF91-TRAF6-BD-GFP 002/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD TRAF6/20220516 1.5nM_cl236_TRAF6_chMyD88-DHF91-TRAF6-BD-GFP 005/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD TRAF6/20220530 1.5nM_cl236_TRAF6_MyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD TRAF6/20220615 1.5nM_cl236_TRAF6_MyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz",
  #"/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD TRAF6/20220530 1.5nM_cl237_TRAF6_MyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD-3xA TRAF6/20220615 1.5nM_cl238_TRAF6_MyD88-DHF91-TRAF6-BD-3xA-GFP 001/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-TIR-TRAF6-BD TRAF6/20220610 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-TIR-TRAF6-BD TRAF6/20220615 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-TIR-TRAF6-BD TRAF6/20220615 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 007/Essential.csv.gz",
  #"/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-TIR-TRAF6-BD TRAF6/20220516 1.5nM_cl241_TRAF6_TIR-TRAF6-BD-GFP 001/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-TIR-TRAF6-BD-3xA TRAF6/20220615 1.5nM_cl242_TRAF6_MyD88-TIR-TRAF6-BD-3xA-GFP 001/Essential.csv.gz",
  #"/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-TIR/20220516 1.5nM_cl244_TRAF6_TIR-GFP 001/Essential.csv.gz",
  file.path("~/Documents/new_pipeline/Essential.csv.gz"),
  file.path("~/Documents/new_pipeline/BDLD_10H-MyD88-TIR-TRAF6-BD-GFP TRAF6/Essential.csv.gz"),
  file.path("~/Documents/new_pipeline/BDLD_13H-MyD88-TIR-TRAF6-BD-GFP TRAF6/Essential.csv.gz"),
  file.path("~/Documents/new_pipeline/BDLD_57H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230602 6nM_cl321-BDLD57H_TRAF6_MyD88_001/Essential.csv.gz"),
  file.path("~/Documents/new_pipeline/BDLD_57H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230602 6nM_cl321-BDLD57H_TRAF6_MyD88_002/Essential.csv.gz"),
  file.path("~/Documents/new_pipeline/BDLD_62H-MyD88-TIR-TRAF6-BD-GFP TRAF6/Essential.csv.gz"),
  file.path("/Volumes/TAYLOR-LAB/Finn/new_pipeline/pending_processing/batch_1_20230519/Output/Essential.csv.gz")
)

Table <- lapply(TablePaths, 
                fread)

#fill=TRUE since not all images are two color
Table <- rbindlist(Table, 
                   fill=TRUE)

Table <- Table[COHORT != "c;028-3E10-GFP"]
Table$COHORT <- gsub("*-MyD88-TIR-TRAF6-BD-GFP TRAF6", replacement = "", x = Table$COHORT)


#get all the different cohorts in Table
unique(Table$COHORT)

unique(Table$LIGAND_DENSITY_CAT)


#start grouping and summarizing the data
Tracks<-
  Table %>% 
  filter(
    NORMALIZED_INTENSITY>=1.5,
    FRAMES_ADJUSTED<=100,
    FRAMES_SINCE_LANDING<=200
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    LIFETIME=max(TIME_ADJUSTED),
    LIFETIME_FRAMES=max(FRAMES_ADJUSTED),
    MAX_NORMALIZED_INTENSITY=max(NORMALIZED_INTENSITY)
  ) %>% 
  filter(
    LIFETIME_FRAMES>=3,
    FRAMES_ADJUSTED==0
  ) %>% 
  as.data.table()

#create data points for the lifetimes in each cell
# Cells <-
#   Tracks %>% 
#   group_by(
#     PROTEIN,
#     IMAGE,
#     COHORT,
#     CELL
#   ) %>% 
#   summarize(LIFETIME=mean(LIFETIME),
#             MAX_NORMALIZED_INTENSITY=mean(MAX_NORMALIZED_INTENSITY)
#   ) %>% 
#   as.data.table()

#Add stats for the images taken to clean up violin
Replicates<-
  Tracks %>% 
  group_by(
    PROTEIN,
    IMAGE,
    COHORT
  ) %>% 
  summarize(LIFETIME=mean(LIFETIME),
            MAX_NORMALIZED_INTENSITY=mean(MAX_NORMALIZED_INTENSITY)
  ) %>% 
  as.data.table()

#create a stats table to add mean and sem of COHORT to the violin plots
Stats<-
  Tracks %>% 
  group_by(
    COHORT,
    PROTEIN
  ) %>%
  mutate(
    SD_LIFETIME=sd(LIFETIME),
    SD_MAX_NORMALIZED_INTENSITY=sd(MAX_NORMALIZED_INTENSITY)
  ) %>% 
  summarize(
    LIFETIME=mean(LIFETIME),
    MAX_NORMALIZED_INTENSITY=mean(MAX_NORMALIZED_INTENSITY),
    SEM_LIFETIME=mean(SD_LIFETIME/sqrt(n())),
    SEM_MAX_NORMALIZED_INTENSITY=mean(SD_MAX_NORMALIZED_INTENSITY)/sqrt(n())
  ) %>% 
  as.data.table()

#set a color palette for the plots
color_violin<-c("MyD88" = "#0bb4ff", "TRAF6" = "#ffa300")

#Plot Lifetimes as Violin
#contain jitter in corresponding violin
ggplot(data=Replicates[COHORT %in% c("BDLD_27H", 
                                     "BDLD_62H",
                                     "MyD88-TIR-TRAF6-BD TRAF6",
                                     "MyD88-TRAF6-BD TRAF6")],
       aes(
         x=COHORT,
         y=LIFETIME,
         fill=PROTEIN
       )
)+
  geom_violin(
    data=
      Tracks[COHORT %in% c("BDLD_27H", 
                           "BDLD_62H",
                           "MyD88-TIR-TRAF6-BD TRAF6",
                           "MyD88-TRAF6-BD TRAF6")],
    position = 
      position_dodge(
        width = 1
      ),
    alpha=0.5
  )+
  fill_palette(
    palette = color_violin
  )+
  geom_point(
    data=Replicates[COHORT %in% c("BDLD_27H", 
                                  "BDLD_62H",
                                  "MyD88-TIR-TRAF6-BD TRAF6",
                                  "MyD88-TRAF6-BD TRAF6")],
    position = 
      position_jitterdodge(
        seed = 1, 
        dodge.width = 1),
    color="black"
  )+
  geom_errorbar(
    data = Stats[COHORT %in% c("BDLD_27H", 
                               "BDLD_62H",
                               "MyD88-TIR-TRAF6-BD TRAF6",
                               "MyD88-TRAF6-BD TRAF6")],
    aes(
      x = COHORT,
      y= LIFETIME,
      ymin = LIFETIME-SEM_LIFETIME,
      ymax = LIFETIME+SEM_LIFETIME),
    color="black",
    width=0.4,
    position= position_dodge(
      width = 1)
  )+
  stat_summary(
    data=Stats[COHORT %in% c("BDLD_27H", 
                             "BDLD_62H",
                             "MyD88-TIR-TRAF6-BD TRAF6",
                             "MyD88-TRAF6-BD TRAF6")],
    fun = "mean",
    geom = "crossbar", 
    width = 0.3,
    size = 0.4,
    position= position_dodge(
      width = 1),
    colour = "black"
  )+
  labs(
    x="Cell Lines",
    y="Lifetime of puncta (s)"
  )+
  scale_y_continuous(
    limits=c(0, 200)
  )+
  theme_bw()

#bargraph for Lifetime  
ggplot(data=Stats,
       aes(
         x=COHORT,
         y=LIFETIME,
         fill=PROTEIN
       )
)+
  fill_palette(
    palette = color_violin
  )+
  geom_col(
    position=position_dodge()
  )

#Plot max intensities
ggplot(data=Replicates,
       aes(
         x=COHORT,
         y=MAX_NORMALIZED_INTENSITY,
         fill=PROTEIN
       )
)+
  fill_palette(
    palette = color_violin
  )+
  geom_violin(
    data = Tracks,
    position = 
      position_dodge(
        width = 1
      ),
    alpha=0.5
  ) +
  geom_point(
    position = 
      position_jitterdodge(
        seed = 1, 
        dodge.width = 1),
    color="black"
  )+
  scale_y_continuous(
    limits = c(1,100),
    trans = "log10",
    breaks = c(1, 2, 4, 6, 8, 10, 20, 40, 70, 100)
  )+
  geom_errorbar(
    data = Stats,
    aes(
      x = COHORT,
      y = LIFETIME,
      ymin = MAX_NORMALIZED_INTENSITY-SEM_MAX_NORMALIZED_INTENSITY,
      ymax = MAX_NORMALIZED_INTENSITY+SEM_MAX_NORMALIZED_INTENSITY),
    width=0.4,
    position = position_dodge(
      width = 1
    ),
    color="black"
  )+
  stat_summary(
    data=Stats,
    fun = "mean",
    geom = "crossbar", 
    width = 0.3,
    size = 0.4,
    position= position_dodge(
      width = 1),
    colour = "black"
  )+
  labs(
    x="Cell line",
    y="Norm. Max Intensity per Cell (a.u.)"
  )+
  theme_bw()

###
#Plot max intensities
ggplot(data=Replicates[COHORT %in% c(#"BDLD_27H", 
                                     "BDLD_62H",
                                     "MyD88-TIR-TRAF6-BD TRAF6",
                                     "MyD88-TRAF6-BD TRAF6")],
       aes(
         x=COHORT,
         y=MAX_NORMALIZED_INTENSITY,
         fill=PROTEIN
       )
)+
  fill_palette(
    palette = color_violin
  )+
  geom_violin(
    # data = Tracks[COHORT %in% c("BDLD_27H", "MyD88-TRAF6-BD TRAF6")],
    data = Tracks[COHORT %in% c(#"BDLD_27H",
                                "BDLD_62H",
                                "MyD88-TIR-TRAF6-BD TRAF6",
                                "MyD88-TRAF6-BD TRAF6")],
    position = 
      position_dodge(
        width = 1
      ),
    alpha=0.5
  ) +
  geom_point(
    position = 
      position_jitterdodge(
        seed = 1, 
        dodge.width = 1),
    color="black"
  )+
  scale_y_continuous(
    limits = c(1,100),
    trans = "log10",
    breaks = c(1, 2, 4, 6, 8, 10, 20, 40, 70, 100)
  )+
  geom_errorbar(
    # data = Stats[COHORT %in% c("BDLD_27H", "MyD88-TRAF6-BD TRAF6")],
    data = Stats[COHORT %in% c(#"BDLD_27H",
                               "BDLD_62H",
                               "MyD88-TIR-TRAF6-BD TRAF6",
                               "MyD88-TRAF6-BD TRAF6")],
    aes(
      x = COHORT,
      y = LIFETIME,
      ymin = MAX_NORMALIZED_INTENSITY-SEM_MAX_NORMALIZED_INTENSITY,
      ymax = MAX_NORMALIZED_INTENSITY+SEM_MAX_NORMALIZED_INTENSITY),
    width=0.4,
    position = position_dodge(
      width = 1
    ),
    color="black"
  )+
  stat_summary(
    # data=Stats[COHORT %in% c("BDLD_27H", "MyD88-TRAF6-BD TRAF6")],
    data = Stats[COHORT %in% c(#"BDLD_27H",
                               "BDLD_62H",
                               "MyD88-TIR-TRAF6-BD TRAF6",
                               "MyD88-TRAF6-BD TRAF6")],
    fun = "mean",
    geom = "crossbar", 
    width = 0.3,
    size = 0.4,
    position= position_dodge(
      width = 1),
    colour = "black"
  )+
  labs(
    x="Cell line",
    y="Norm. Max Intensity per Cell (a.u.)"
  )+
  theme_bw()

###



#bargraph for Max intensities  
ggplot(data=Stats,
       aes(
         x=COHORT,
         y=MAX_NORMALIZED_INTENSITY,
         fill=PROTEIN
       )
)+
  fill_palette(
    palette = color_violin
  )+
  geom_col(
    position=position_dodge()
  )

#plot Lifetime over normalized max intensity
MyD88Int<-
  Table %>% 
  filter(
    PROTEIN=="MyD88",
    NORMALIZED_INTENSITY>=1.5,
    FRAMES_ADJUSTED<=100,
    FRAMES_SINCE_LANDING<=200
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    LIFETIME=max(TIME_ADJUSTED),
    LIFETIME_FRAMES=max(FRAMES_ADJUSTED),
    MAX_NORMALIZED_INTENSITY=max(NORMALIZED_INTENSITY)
  ) %>% 
  filter(
    LIFETIME_FRAMES>=3,
    FRAMES_ADJUSTED==0
  ) %>% 
  as.data.table()

#geom_hex
ggplot(
  data=
    MyD88Int,
  aes(
    x=LIFETIME,
    y=MAX_NORMALIZED_INTENSITY
  )
)+
  geom_hex(
    bins=15,
    alpha=0.7
  )+
  geom_smooth(
    method = "lm", 
    se = FALSE
  )+
  stat_cor(
    method="spearman",
    color="red"
  )+
  scale_y_continuous(
    limits=c(0,20)
  )+
  scale_x_continuous(
    limits=c(0, 200)
  )+
  scale_fill_viridis(
    trans = 'log10'
  )+
  labs(
    x="Lifetime of MyD88 in s",
    y="Normalized Max Intensity of MyD88 [a. u.]"
  )+
  facet_wrap(
    ~COHORT
  )+
  theme_classic()

#plot Lifetime over normalized max intensity
TRAF6Int<-
  Table %>% 
  filter(
    PROTEIN=="TRAF6",
    NORMALIZED_INTENSITY>=1.5,
    FRAMES_ADJUSTED<=100,
    FRAMES_SINCE_LANDING<=200
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    LIFETIME=max(TIME_ADJUSTED),
    LIFETIME_FRAMES=max(FRAMES_ADJUSTED),
    MAX_NORMALIZED_INTENSITY=max(NORMALIZED_INTENSITY)
  ) %>% 
  filter(
    LIFETIME_FRAMES>=3,
    FRAMES_ADJUSTED==0
  ) %>% 
  as.data.table()

#gemo_hex
ggplot(
  data=
    TRAF6Int,
  aes(
    x=LIFETIME,
    y=MAX_NORMALIZED_INTENSITY
  )
)+
  geom_hex(
    bins=15,
    alpha=0.7
  )+
  scale_y_continuous(
    limits=c(
      0,16
    )
  )+
  scale_x_continuous(
    limits=c(
      0, 200
    )
  )+
  facet_wrap(
    ~COHORT
  )+
  scale_fill_viridis(
    trans = 'log10'
  )+
  geom_smooth(
    method = "lm", 
    se = FALSE
  )+
  stat_cor(
    method="spearman",
    color="red"
  )+
  labs(
    x="Lifetime of TRAF6 in s",
    y="Normalized Max Intensity of TRAF6 [a. u.]"
  )+
  theme_classic()

ggplot(
  data=
    TRAF6Int,
  aes(
    x=LIFETIME,
    y=MAX_NORMALIZED_INTENSITY
  )
)+
  geom_point(
  )+
  scale_y_continuous(
    limits=c(
      0,16
    )
  )+
  scale_x_continuous(
    limits=c(
      0, 200
    )
  )+
  facet_wrap(
    ~COHORT
  )+
  geom_smooth(
    method = "lm", 
    se = FALSE
  )+
  stat_cor(
    method="spearman",
    color="red"
  )+
  labs(
    x="Lifetime of TRAF6 in s",
    y="Normalized Max Intensity of TRAF6 [a. u.]"
  )+
  theme_classic()

