library(pacman)

pacman::p_load(ggplot2, ggdark, data.table, dplyr, ggfx, viridis, ggridges, RColorBrewer, ggpubr)
filter <- dplyr::filter

data_tay <- "/Volumes/TAYLOR-LAB/"

# check if the path exists 
if (file.exists(data_tay)){
  main = file.path("/Volumes/TAYLOR-LAB")
} else { 
  # use local path
  main = file.path("~/Documents")
}

#load tables
TablePaths <- c(
  # file.path(main, "/Mauriz/03_data_analysis/chMyD88-DHF91-TRAF6-BD-GFP/20220401 6nM_cl221_chMyD88-DHF91-TRAF6-BD-GFP_002/Essential.csv.gz"),
  # file.path(main, "/Mauriz/03_data_analysis/chMyD88-DHF91-TRAF6-BD-GFP/20220401 6nM_cl222_chMyD88-DHF91-TRAF6-BD-GFP_001/Essential.csv.gz"),
  # file.path(main, "/Mauriz/03_data_analysis/chMyD88-DHF91-TRAF6-BD-GFP/20220401 6nM_cl222_chMyD88-DHF91-TRAF6-BD-GFP_002/Essential.csv.gz"),
  # file.path(main, "/Mauriz/03_data_analysis/DHF58v2mid_DHF58v.2/20220401 6nM_DHF58v2mid_DHF58v.2_001/Essential.csv.gz"),
  
  # normal synthetic
  file.path(main, "/Mauriz/03_data_analysis/MyD88-TRAF6-BD TRAF6/20220516 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 001/Essential.csv.gz"),
  file.path(main, "/Mauriz/03_data_analysis/MyD88-TRAF6-BD TRAF6/20220610 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 001/Essential.csv.gz"),
  file.path(main, "/Mauriz/03_data_analysis/MyD88-TRAF6-BD TRAF6/20220615 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 002/Essential.csv.gz"),
  # "~/Documents/Mauriz_R/MyD88-TRAF6-BD TRAF6/20220516 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 001/Essential.csv.gz"),
  # "~/Documents/Mauriz_R/MyD88-TRAF6-BD TRAF6/20220610 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 001/Essential.csv.gz"),
  # "~/Documents/Mauriz_R/MyD88-TRAF6-BD TRAF6/20220615 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 002/Essential.csv.gz")
  
  # DHF91
  # file.path(main, "/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD TRAF6/20220516 1.5nM_cl236_TRAF6_chMyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz"),
  # file.path(main, "/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD TRAF6/20220516 1.5nM_cl236_TRAF6_chMyD88-DHF91-TRAF6-BD-GFP 002/Essential.csv.gz"),
  # file.path(main, "/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD TRAF6/20220516 1.5nM_cl236_TRAF6_chMyD88-DHF91-TRAF6-BD-GFP 005/Essential.csv.gz"),
  # file.path(main, "/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD TRAF6/20220530 1.5nM_cl236_TRAF6_MyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz"),
  # file.path(main, "/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD TRAF6/20220615 1.5nM_cl236_TRAF6_MyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz"),
  # file.path(main, "/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD TRAF6/20220530 1.5nM_cl237_TRAF6_MyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz"),
  
  # negative controls 3xA
  file.path(main, "/Mauriz/03_data_analysis/MyD88-TRAF6-BD-3xA TRAF6/20220615 1.5nM_cl234_TRAF6_MyD88-TRAF6-BD-3xA-GFP 008/Essential.csv.gz"),
  file.path(main, "/Mauriz/03_data_analysis/MyD88-TIR-TRAF6-BD-3xA TRAF6/20220615 1.5nM_cl242_TRAF6_MyD88-TIR-TRAF6-BD-3xA-GFP 001/Essential.csv.gz"),
  # file.path(main, "/Mauriz/03_data_analysis/MyD88-DHF91-TRAF6-BD-3xA TRAF6/20220615 1.5nM_cl238_TRAF6_MyD88-DHF91-TRAF6-BD-3xA-GFP 001/Essential.csv.gz"),
  
  # TIR only
  file.path(main, "/Mauriz/03_data_analysis/MyD88-TIR-TRAF6-BD TRAF6/20220610 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz"),
  file.path(main, "/Mauriz/03_data_analysis/MyD88-TIR-TRAF6-BD TRAF6/20220615 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz"),
  file.path(main, "/Mauriz/03_data_analysis/MyD88-TIR-TRAF6-BD TRAF6/20220615 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 007/Essential.csv.gz"),
  file.path(main, "/Mauriz/03_data_analysis/MyD88-TIR-TRAF6-BD TRAF6/20220516 1.5nM_cl241_TRAF6_TIR-TRAF6-BD-GFP 001/Essential.csv.gz"),
  
  # BDLD - Bacterial Death-Like Domains
  file.path(main, "/Finn/new_pipeline/pending_processing/batch_1_20230519/Output/Essential.csv.gz"),
  file.path(main, "/Finn/new_pipeline/pending_processing/batch_2_20230602/Processing/06_Colocalization/BDLD_62H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230602 6nM_cl318-BDLD62H_TRAF6_MyD88_001/Essential.csv.gz")
  # "~/Documents/new_pipeline/pending_processing/batch_1_20230519/Output/Essential.csv.gz")
  )

Table <- lapply(TablePaths, 
                fread)

#fill=TRUE since not all images are two color
Table <- rbindlist(Table, 
                   fill=TRUE)

#start grouping and summarizing the data
#Find out the percentage of spots recruiting TRAF6
Colocalisation_List <-
  Table %>%
  filter(
    PROTEIN == "MyD88"
  ) %>% 
  arrange(
    UNIVERSAL_SPOT_ID
  ) %>% 
  mutate(
    COLOCALIZATION = COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 1.5 #threshold at which recruitment is counted
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
    DWELL_FRAMES >= 1 #only count 3 continuous frames as recruitment
  ) %>% 
  distinct(
    UNIVERSAL_TRACK_ID
  ) %>% 
  as.data.table()  

Track_rec <- c(unique(Colocalisation_List$UNIVERSAL_TRACK_ID))

Puncta_Cell<-
  Table %>% 
  filter(
    PROTEIN=="MyD88",
    NORMALIZED_INTENSITY >= 1.5,
    FRAMES_ADJUSTED <= 100,
    FRAMES_SINCE_LANDING<=200,
    COHORT!="chMyD88-DHF91-TRAF6-BD-GFP",
    COHORT!="DHF58v2mid_DHF58v.2",
    COHORT!="MyD88-TIR"
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    LIFETIME=max(TIME_ADJUSTED),
    LIFETIME_FRAMES=max(FRAMES_ADJUSTED),
    MAX_NORMALIZED_INTENSITY=max(NORMALIZED_INTENSITY),
    MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1=max(COMPLEMENTARY_NORMALIZED_INTENSITY_1),
    RECRUITEMENT=case_when(
      UNIVERSAL_TRACK_ID %in% Track_rec ~ "+ve",
      TRUE ~ "-ve"
    )
  ) %>% 
  filter(
    LIFETIME_FRAMES>=3,
    FRAMES_ADJUSTED==0
  ) %>%
  group_by(
    COHORT,
    IMAGE,
    CELL
  ) %>%
  summarize(
    T6p=sum(RECRUITEMENT=="+ve")/(sum(RECRUITEMENT=="+ve")+sum(RECRUITEMENT=="-ve")),
    T6n=sum(RECRUITEMENT=="-ve")/(sum(RECRUITEMENT=="+ve")+sum(RECRUITEMENT=="-ve"))
  ) %>%
  as.data.table()


Stats<-
  #create the stats variable to add mean etc. to plots
  Stats<-
  Puncta_Cell %>% 
  group_by(
    COHORT
  ) %>%
  mutate(
    SD_T6p=sd(T6p)
  ) %>% 
  summarize(
    T6p=mean(T6p),
    SEM_T6p=mean(SD_T6p/sqrt(n())),
  ) %>% 
  as.data.table()

#Add stats for the images taken to clean up violin
Replicates<-
  Puncta_Cell %>% 
  group_by(
    IMAGE,
    COHORT
  ) %>% 
  summarize(T6p=mean(T6p)
  ) %>% 
  as.data.table()  

#plot percentage of tracks recruiting TRAF6 per cell as a violin
ggplot(
  data=Puncta_Cell,
  aes(
    x=COHORT,
    y=T6p
  )
)+
  geom_violin(
    alpha=0.75,
    fill="gray",
    scale = "width"
  )+
  geom_jitter(data=Replicates,
              color="black"
  )+
  geom_errorbar(
    data = Stats,
    aes(
      x = COHORT,
      y= T6p,
      ymin = T6p-SEM_T6p,
      ymax = T6p+SEM_T6p),
    color="black",
    width=0.4
  )+
  stat_summary(
    data=Stats,
    fun = "mean",
    geom = "crossbar", 
    width = 0.3,
    linewidth = 0.4,
    position= position_dodge(
      width = 1),
    colour = "black"
  )+
  labs(
    x="Cell Lines",
    y="Percentage of Puncta recruiting TRAF6 per Cell"
  )+
  theme_classic()

#Percentage of Spots above threshold Lifetime/ Max Intensity per Cell per Cell line



#plot MyD88 intensity over TRAF6 intensity using geom_hex
MyD88Int<-
  Table %>% 
  filter(
    PROTEIN=="MyD88",
    FRAMES_ADJUSTED<=100,
    FRAMES_SINCE_LANDING<=200,
    COHORT!="chMyD88-DHF91-TRAF6-BD-GFP",
    COHORT!="DHF58v2mid_DHF58v.2",
    COHORT!="MyD88-TIR",
    NORMALIZED_INTENSITY>=1.5
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    LIFETIME=max(TIME_ADJUSTED),
    LIFETIME_FRAMES=max(FRAMES_ADJUSTED),
    MAX_NORMALIZED_INTENSITY=max(NORMALIZED_INTENSITY),
    MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1=max(COMPLEMENTARY_NORMALIZED_INTENSITY_1)
  ) %>% 
  filter(
    LIFETIME_FRAMES>=3,
    FRAMES_ADJUSTED==0
  ) %>% 
  as.data.table()

ggplot(
  data=
    MyD88Int,
  aes(
    y=MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1,
    x=MAX_NORMALIZED_INTENSITY
  )
)+
  geom_hex(
    bins=15
  )+
  scale_y_continuous(
    limits=c(
      0,20
    )
  )+
  scale_x_continuous(
    limits=c(
      0, 50
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
    x="Normalized Max Intensity of MyD88",
    y="Normalized Max Intensity of TRAF6"
  )+
  theme_classic()

#plot MyD88 intensity over TRAF6 intensity using geom_path
TRAF6col <-
  Table %>%
  filter(PROTEIN=="MyD88") %>% 
  mutate(
    COMPLEMENTARY_NORMALIZED_INTENSITY_1 = round(COMPLEMENTARY_NORMALIZED_INTENSITY_1)
  ) %>% 
  group_by(
    COHORT,
    COMPLEMENTARY_NORMALIZED_INTENSITY_1
  ) %>% 
  summarize(
    NORMALIZED_INTENSITY = mean(NORMALIZED_INTENSITY)
  ) %>%
  as.data.table()

ggplot(
  TRAF6col,
  aes(
    x=NORMALIZED_INTENSITY,
    y=COMPLEMENTARY_NORMALIZED_INTENSITY_1,
    color=COHORT
  )
)+
  geom_path(
  )+
  scale_x_continuous(
    limits=c(0, 50)
  )+
  scale_y_continuous(
    limits=c(0, 10)
  )+
  labs(
    x="Normalized Intensity of MyD88",
    y="Normalized Intensity of TRAF6"
  )+
  theme_classic()

#percentage of puncta at distinct size/ Lifetime to recruit TRAF6
PCTTRAF6<-
  Table %>% 
  filter(
    PROTEIN=="MyD88",
    NORMALIZED_INTENSITY >= 1.5,
    FRAMES_ADJUSTED <= 100,
    FRAMES_SINCE_LANDING<=200,
    COHORT!="chMyD88-DHF91-TRAF6-BD-GFP",
    COHORT!="DHF58v2mid_DHF58v.2",
    COHORT!="MyD88-TIR"
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    LIFETIME=max(TIME_ADJUSTED),
    LIFETIME_FRAMES=max(FRAMES_ADJUSTED),
    MAX_NORMALIZED_INTENSITY=max(NORMALIZED_INTENSITY),
    MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1=max(COMPLEMENTARY_NORMALIZED_INTENSITY_1),
    RECRUITEMENT=case_when(
      UNIVERSAL_TRACK_ID %in% Track_rec ~ "+ve",
      TRUE ~ "-ve"
    ),
    MAX_NORMALIZED_INTENSITY=round(MAX_NORMALIZED_INTENSITY),
    LIFETIME=round(LIFETIME)
  ) %>% 
  filter(
    LIFETIME_FRAMES>=3,
    FRAMES_ADJUSTED==0
  ) %>% 
  as.data.table()

Stats <- 
  PCTTRAF6 %>% 
  group_by(
    COHORT
  ) %>% 
  summarise(
    MAX_NORMALIZED_INTENSITY=mean(MAX_NORMALIZED_INTENSITY),
    LIFETIME=mean(LIFETIME)
  ) %>% 
  as.data.table()

#Percent of TRAF6 at different Max normalized Intensities
PCTTRAF6Int<-
  PCTTRAF6 %>% 
  group_by(
    MAX_NORMALIZED_INTENSITY,
    COHORT
  ) %>% 
  filter(n() >= 3) %>%
  summarize(
    T6p=sum(RECRUITEMENT=="+ve")/(sum(RECRUITEMENT=="+ve")+sum(RECRUITEMENT=="-ve")),
    T6n=sum(RECRUITEMENT=="-ve")/(sum(RECRUITEMENT=="+ve")+sum(RECRUITEMENT=="-ve"))
  ) %>%
  as.data.table()

ggplot(
  data=PCTTRAF6Int,
  aes(
    x=MAX_NORMALIZED_INTENSITY,
    y=T6p*100,
    fill=COHORT
  )
)+
  geom_smooth(
    level=0.90,
    color="black"
  )+
  scale_x_continuous(
    limits = c(0,50)
  )+
  scale_y_continuous(
    limits=c(0,100)
  )+
  geom_jitter(
    color="black"
  )+
  facet_wrap(
    ~COHORT
  )+
  labs(
    x="Max Normalized Intensity of Tracks",
    y="% of Tracks recruiting TRAF6"
  )+
  theme_bw()

#density distribution of Puncta to Max norm Int
ggplot(
  data=PCTTRAF6,
  aes(
    x=MAX_NORMALIZED_INTENSITY,
    y=..ndensity..,
    fill=COHORT
  )
)+
  geom_histogram(
    binwidth = 1,
    colour="black",
    alpha=0.75
  )+
  geom_vline(
    data = Stats,
    aes(xintercept = MAX_NORMALIZED_INTENSITY),
    color="black", linetype="dashed", linewidth=1
  )+
  scale_x_continuous(
    limits = c(0,50)
  )+
  facet_wrap(
    ~COHORT
  )+
  labs(
    x="Max Normalized Intensity of Tracks",
    y="Density"
  )+
  theme_bw()


#Percent of TRAF6 recruitement at different Lifetimes      
PCTTRAF6LT<-
  PCTTRAF6 %>% 
  group_by(
    LIFETIME,
    COHORT
  ) %>% 
  filter(n() >= 3) %>%
  summarize(
    T6p=sum(RECRUITEMENT=="+ve")/(sum(RECRUITEMENT=="+ve")+sum(RECRUITEMENT=="-ve")),
    T6n=sum(RECRUITEMENT=="-ve")/(sum(RECRUITEMENT=="+ve")+sum(RECRUITEMENT=="-ve"))
  ) %>%
  as.data.table()


ggplot(
  data=PCTTRAF6LT,
  aes(
    x=LIFETIME,
    y=T6p,
    fill=COHORT
  )
)+
  geom_smooth(
    level=0.90,
    color="black"
  )+
  scale_x_continuous(
    limits = c(0,400)
  )+
  geom_jitter(
    color="black"
  )+
  facet_wrap(
    ~COHORT
  )+
  labs(
    x="Lifetime of Tracks",
    y="% of Tracks recruiting TRAF6"
  )+
  theme_bw()

ggplot(
  data=PCTTRAF6,
  aes(
    x=LIFETIME,
    y=..ndensity..,
    fill=COHORT
  )
)+
  geom_histogram(
    binwidth = 5,
    colour="black",
    alpha=0.75
  )+
  geom_vline(
    data = Stats,
    aes(xintercept = LIFETIME),
    color="black", linetype="dashed", linewidth=1
  )+
  scale_x_continuous(
    limits = c(0,400)
  )+
  facet_wrap(
    ~COHORT
  )+
  labs(
    x="Lifetime of Tracks",
    y="Density"
  )+
  theme_bw()

