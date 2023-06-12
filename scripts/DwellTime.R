library(pacman)

pacman::p_load(ggplot2, ggdark, data.table, dplyr, ggfx, viridis, ggridges, RColorBrewer, ggpubr, lemon, cowplot)
filter <- dplyr::filter


#load tables
TablePaths <- c(
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/chMyD88-DHF91-TRAF6-BD-GFP/20220401 6nM_cl221_chMyD88-DHF91-TRAF6-BD-GFP_002/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/chMyD88-DHF91-TRAF6-BD-GFP/20220401 6nM_cl222_chMyD88-DHF91-TRAF6-BD-GFP_001/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/chMyD88-DHF91-TRAF6-BD-GFP/20220401 6nM_cl222_chMyD88-DHF91-TRAF6-BD-GFP_002/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/DHF58v2mid_DHF58v.2/20220401 6nM_DHF58v2mid_DHF58v.2_001/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-TRAF6-BD TRAF6/20220516 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 001/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-TRAF6-BD TRAF6/20220610 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 001/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-TRAF6-BD TRAF6/20220615 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 002/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-TRAF6-BD-3xA TRAF6/20220615 1.5nM_cl234_TRAF6_MyD88-TRAF6-BD-3xA-GFP 008/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD TRAF6/20220516 1.5nM_cl236_TRAF6_chMyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD TRAF6/20220516 1.5nM_cl236_TRAF6_chMyD88-DHF91-TRAF6-BD-GFP 002/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD TRAF6/20220516 1.5nM_cl236_TRAF6_chMyD88-DHF91-TRAF6-BD-GFP 005/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD TRAF6/20220530 1.5nM_cl236_TRAF6_MyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD TRAF6/20220615 1.5nM_cl236_TRAF6_MyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD TRAF6/20220530 1.5nM_cl237_TRAF6_MyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD-3xA TRAF6/20220615 1.5nM_cl238_TRAF6_MyD88-DHF91-TRAF6-BD-3xA-GFP 001/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-TIR-TRAF6-BD TRAF6/20220610 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-TIR-TRAF6-BD TRAF6/20220615 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-TIR-TRAF6-BD TRAF6/20220615 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 007/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-TIR-TRAF6-BD TRAF6/20220516 1.5nM_cl241_TRAF6_TIR-TRAF6-BD-GFP 001/Essential.csv.gz",
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-TIR-TRAF6-BD-3xA TRAF6/20220615 1.5nM_cl242_TRAF6_MyD88-TIR-TRAF6-BD-3xA-GFP 001/Essential.csv.gz"
  # "/Volumes/TAYLOR-LAB/Mauriz/03_data_analysis/MyD88-TIR/20220516 1.5nM_cl244_TRAF6_TIR-GFP 001/Essential.csv.gz"
  # "/Volumes/TAYLOR-LAB/Finn/new_pipeline/pending_processing/batch_1_20230519/Output/Essential.csv.gz",
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
Table$COHORT <- gsub("*-MyD88-TIR-TRAF6-BD-GFP TRAF6", replacement = "", x = Table$COHORT)
Table <- Table[COHORT != "c;028-3E10-GFP"]

#Get Table with Dwell frames & Dwell Time
Dwell_table <- Table %>%
  filter(PROTEIN == "MyD88") %>% 
  arrange(UNIVERSAL_SPOT_ID) %>% 
  #threshold at which recruitment is counted
  mutate(COLOCALIZATION = COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 1 ) %>% 
  group_by( UNIVERSAL_TRACK_ID) %>%
  #new column which creates a grouping variable for continuous frames which are above threshold
  mutate(STREAK = cumsum(!COLOCALIZATION)   ) %>%
  #group continuous frames above threshold together
  group_by(COHORT,IMAGE,UNIVERSAL_TRACK_ID,STREAK) %>% 
  #number of frames that complementary protein is above threshold in a continuous stretch
  summarise( DWELL_FRAMES = sum(COLOCALIZATION)) %>% 
  #I want to look at transient recruitment events so even 1 frame is interesting for me, different cutoffs are prbly useful
  # now at least 3 frames = 12 s == stable recruitment
  filter(DWELL_FRAMES >= 3) %>% 
  mutate(DWELL_TIME = DWELL_FRAMES*4) %>% 
  as.data.table()



#create the stats variable to add mean etc. to plots
Stats<-
  Dwell_table %>% 
  group_by(
    COHORT
  ) %>%
  mutate(
    SD_DWELL_FRAMES = sd(DWELL_FRAMES),
    SD_DWELL_TIME = sd(DWELL_TIME)
  ) %>% 
  summarize(
    DWELL_FRAMES = mean(DWELL_FRAMES),
    DWELL_TIME = mean(DWELL_TIME),
    SEM_DWELL_FRAMES = mean(SD_DWELL_FRAMES/sqrt(n()))
  ) %>% 
  as.data.table()

#Add stats for the images taken to clean up violin
Replicates<-
  Dwell_table %>% 
  group_by(
    IMAGE,
    COHORT
  ) %>% 
  summarize(
    DWELL_FRAMES = mean(DWELL_FRAMES),
    DWELL_TIME = mean(DWELL_TIME),
  ) %>% 
  as.data.table()

#Set colors
#color_violin_rec<-c("+ve" = "#e60049", "-ve" = "#00bfa0")

ggplot(data=Dwell_table[COHORT %in% c("BDLD_27H",
                                      "BDLD_62H",
                                      "MyD88-TIR-TRAF6-BD TRAF6",
                                      "MyD88-TRAF6-BD TRAF6")], aes(x = DWELL_TIME,
         # y=after_stat(ndensity),
         fill = COHORT
       )
)+
  geom_histogram(
    data = Dwell_table[COHORT %in% c("BDLD_27H",
                                     #"BDLD_62H",
                                     "MyD88-TIR-TRAF6-BD TRAF6",
                                     "MyD88-TRAF6-BD TRAF6")],
    alpha = 0.5,
    color = "black",
    
    binwidth = 4,
  ) +
  scale_y_continuous(
      trans = "log2"
  ) +
  facet_rep_wrap(
    repeat.tick.labels = TRUE,
    scales = "free_y",
    ~COHORT,
    ncol = 1
  )+
  scale_x_continuous(
    limits = c(0,400),
    breaks = seq(0,360, by = 60)
  )+
  geom_vline(
    data = Stats[COHORT %in% c("BDLD_27H",
                               #"BDLD_62H",
                               "MyD88-TIR-TRAF6-BD TRAF6",
                               "MyD88-TRAF6-BD TRAF6")],
    aes(xintercept = DWELL_TIME),
    color="black", linetype="dashed", linewidth=1
  ) +
    # geom_jitter(
    # )+
    # stat_summary(
    #   data=Stats,
    #   fun = "mean",
    #   geom = "crossbar", 
    #   width = 0.3,
    #   size = 0.4,
    #   colour = "black"
    # )+
    # scale_y_log10(
    #   limits = c(1,100),
    #   breaks = scales::breaks_log(n=10, base=10)
    # )+
    # geom_errorbar(
    #   data = Stats,
    #   aes(
    #     x = COHORT,
    #     y= DWELL_FRAMES,
    #     ymin = DWELL_FRAMES-SEM_DWELL_FRAMES,
    #     ymax = DWELL_FRAMES+SEM_DWELL_FRAMES),
    #   width=0.4,
    #   color="black" 
    # )+
    # labs(
    #   x="Cell line",
    #   # y="Normalized Max Intensity (a.u.)"
  # )+
  # scale_x_continuous(limits = c(0,20)) + 
  # scale_y_continuous(
  #   trans = "log2",
  #   limits = c(0,50)) +
xlab("Dwell Time (s)") +
  ylab("count (log2)") +
#labs(
#  
#)+  
theme_cowplot() +
  theme(legend.position = "none")

ggplot(data=Dwell_table,
       aes(
         x=COHORT,
         y=DWELL_FRAMES,
         fill = COHORT
       )
)+
  geom_boxplot(
    data = Dwell_table, 
    alpha = 0.5,
  )+
  scale_y_continuous(
    limits = c(0,50)
    # )+
    # geom_vline(
    #   data = Stats,
    #   aes(xintercept = DWELL_FRAMES),
    #   color="black", linetype="dashed", linewidth=1
    # geom_jitter(
    # )+
    # stat_summary(
    #   data=Stats,
    #   fun = "mean",
    #   geom = "crossbar", 
    #   width = 0.3,
    #   size = 0.4,
    #   colour = "black"
    # )+
    # scale_y_log10(
    #   limits = c(1,100),
    #   breaks = scales::breaks_log(n=10, base=10)
    # )+
    # geom_errorbar(
    #   data = Stats,
    #   aes(
    #     x = COHORT,
    #     y= DWELL_FRAMES,
    #     ymin = DWELL_FRAMES-SEM_DWELL_FRAMES,
    #     ymax = DWELL_FRAMES+SEM_DWELL_FRAMES),
    #   width=0.4,
    #   color="black" 
    # )+
    # labs(
    #   x="Cell line",
    #   # y="Normalized Max Intensity (a.u.)"
  )+
  theme_bw()

ggplot(data=Replicates,
       aes(
         x=COHORT,
         y=DWELL_FRAMES,
         fill = COHORT
       )
)+
  geom_violin(
    data = Dwell_table, 
    alpha = 0.5,
  )+
  geom_jitter(
  )+
  stat_summary(
    data=Stats,
    fun = "mean",
    geom = "crossbar",
    width = 0.3,
    size = 0.4,
    colour = "black"
  )+
  scale_y_log10(
    limits = c(1,100),
    breaks = scales::breaks_log(n=10, base=10)
  )+
  geom_errorbar(
    data = Stats,
    aes(
      x = COHORT,
      y= DWELL_FRAMES,
      ymin = DWELL_FRAMES-SEM_DWELL_FRAMES,
      ymax = DWELL_FRAMES+SEM_DWELL_FRAMES),
    width=0.4,
    color="black"
  )+
  labs(
    x="Cell line",
    # y="Normalized Max Intensity (a.u.)"
  )+
  theme_bw()













#
Dwell_table<-
  Table %>%
  filter(
    PROTEIN == "MyD88"
  ) %>% 
  arrange(
    UNIVERSAL_SPOT_ID
  ) %>% 
  mutate(
    COLOCALIZATION = COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 1 #threshold at which recruitment is counted
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    STREAK = cumsum(!COLOCALIZATION) #new column which creates a grouping variable for continuous frames which are above threshold
  ) %>%
  group_by(
    COHORT,
    IMAGE,
    UNIVERSAL_TRACK_ID,
    STREAK #group continuous frames above threshold together
  ) %>% 
  summarise(
    DWELL_FRAMES = sum(COLOCALIZATION) #number of frames that complementary protein is above threshold in a continuous stretch
  ) %>% 
  filter(
    DWELL_FRAMES >= 3 #only count 3 continuous frames as recruitment
  ) %>% 
  as.data.table()

