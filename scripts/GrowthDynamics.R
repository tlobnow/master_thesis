library(pacman)

pacman::p_load(ggplot2, ggdark, data.table, dplyr, ggfx, viridis, ggridges, RColorBrewer, ggpubr, lemon, gghalves)
filter <- dplyr::filter

#load tables
TablePaths <- 
  c(#"/Volumes/TAYLOR-LAB/Mauriz/05_data_analysis_new_parameters/MyD88-TRAF6-BD TRAF6/20220516 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 001/Essential.csv.gz",
    # "/Volumes/TAYLOR-LAB/Mauriz/05_data_analysis_new_parameters/MyD88-TRAF6-BD TRAF6/20220610 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 001/Essential.csv.gz",
    # "/Volumes/TAYLOR-LAB/Mauriz/05_data_analysis_new_parameters/MyD88-TRAF6-BD TRAF6/20220615 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 002/Essential.csv.gz",
    # "/Volumes/TAYLOR-LAB/Mauriz/05_data_analysis_new_parameters/MyD88-TIR-TRAF6-BD TRAF6/20220610 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz",
    # "/Volumes/TAYLOR-LAB/Mauriz/05_data_analysis_new_parameters/MyD88-TIR-TRAF6-BD TRAF6/20220615 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz",
    # "/Volumes/TAYLOR-LAB/Mauriz/05_data_analysis_new_parameters/MyD88-TIR-TRAF6-BD TRAF6/20220615 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 007/Essential.csv.gz",
    # "/Volumes/TAYLOR-LAB/Mauriz/05_data_analysis_new_parameters/MyD88-DHF91-TRAF6-BD TRAF6/20220615 1.5nM_cl236_TRAF6_MyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz",
    # "/Volumes/TAYLOR-LAB/Mauriz/05_data_analysis_new_parameters/MyD88-DHF91-TRAF6-BD TRAF6/20221207 4nM_cl236_TRAF6_MyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz"
    "/Volumes/TAYLOR-LAB/Finn/new_pipeline/pending_processing/batch_1_20230519/Output/Essential.csv.gz"
  )

Table <- lapply(TablePaths, fread)

#fill=TRUE since not all images are two color
Table <- rbindlist(Table, fill=TRUE)

#Get the growth for different complex sizes ie the k_on rate

Growth<- Table %>% 
  filter(
    PROTEIN == "MyD88",
    FRAMES_ADJUSTED <= 50,
    FRAMES_SINCE_LANDING <= 100) %>% 
  
  #to be able to correlate the value of the present frame to the next frame
  group_by(UNIVERSAL_TRACK_ID) %>% 
  
  #to order the spots of each UNIVERSAL_TRACK_ID by FRAME
  arrange(FRAME, .by_group = TRUE) %>% 
  mutate(ROUND_NORMALIZED_INTENSITY = round(NORMALIZED_INTENSITY),
         DELTA_NORMALIZED_INTENSITY_1 = lead(ROUND_NORMALIZED_INTENSITY, n=1) - ROUND_NORMALIZED_INTENSITY,
         DELTA_NORMALIZED_INTENSITY_2 = lead(ROUND_NORMALIZED_INTENSITY, n=2) - ROUND_NORMALIZED_INTENSITY,
         DELTA_NORMALIZED_INTENSITY_5 = lead(ROUND_NORMALIZED_INTENSITY, n=5) - ROUND_NORMALIZED_INTENSITY,
         DELTA_NORMALIZED_INTENSITY_8 = lead(ROUND_NORMALIZED_INTENSITY, n=8) - ROUND_NORMALIZED_INTENSITY) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% 
  group_by(COHORT, ROUND_NORMALIZED_INTENSITY) %>% 
  #at least 5 events to be included in summary
  filter(n() > 5 ) %>% 
  summarise(DELTA_AVG_NORMALIZED_INTENSITY_1 = mean(DELTA_NORMALIZED_INTENSITY_1, na.rm = TRUE),
            DELTA_AVG_NORMALIZED_INTENSITY_2 = mean(DELTA_NORMALIZED_INTENSITY_2, na.rm = TRUE),
            DELTA_AVG_NORMALIZED_INTENSITY_5 = mean(DELTA_NORMALIZED_INTENSITY_5, na.rm = TRUE),
            DELTA_AVG_NORMALIZED_INTENSITY_8 = mean(DELTA_NORMALIZED_INTENSITY_8, na.rm = TRUE),
            EVENTS = n()) %>% 
  mutate(
    # NORM_DELTA_AVG_NORMALIZED_INTENSITY_1 = DELTA_AVG_NORMALIZED_INTENSITY_1/ROUND_NORMALIZED_INTENSITY,
    # NORM_DELTA_AVG_NORMALIZED_INTENSITY_2 = DELTA_AVG_NORMALIZED_INTENSITY_2/ROUND_NORMALIZED_INTENSITY,
    # NORM_DELTA_AVG_NORMALIZED_INTENSITY_5 = DELTA_AVG_NORMALIZED_INTENSITY_5/ROUND_NORMALIZED_INTENSITY,
    # NORM_DELTA_AVG_NORMALIZED_INTENSITY_8 = DELTA_AVG_NORMALIZED_INTENSITY_8/ROUND_NORMALIZED_INTENSITY,
    
    NORM_DELTA_AVG_NORMALIZED_INTENSITY_1 = DELTA_AVG_NORMALIZED_INTENSITY_1/1,
    NORM_DELTA_AVG_NORMALIZED_INTENSITY_2 = DELTA_AVG_NORMALIZED_INTENSITY_2/2,
    NORM_DELTA_AVG_NORMALIZED_INTENSITY_5 = DELTA_AVG_NORMALIZED_INTENSITY_5/5,
    NORM_DELTA_AVG_NORMALIZED_INTENSITY_8 = DELTA_AVG_NORMALIZED_INTENSITY_8/8)

ggplot(
  data = Growth,
  aes(x = ROUND_NORMALIZED_INTENSITY,
      y = NORM_DELTA_AVG_NORMALIZED_INTENSITY_1)) +
  geom_path(aes(color = "1 Frame")) +
  geom_path(aes(y = NORM_DELTA_AVG_NORMALIZED_INTENSITY_2, color = "2 Frames")) +
  geom_path(aes(y = NORM_DELTA_AVG_NORMALIZED_INTENSITY_5, color = "5 Frames")) +
  facet_wrap(~COHORT) +
  # geom_path(
  #   aes(
  #     y = NORM_DELTA_AVG_NORMALIZED_INTENSITY_8),
  #   color = "orange"
  # )+
  scale_x_continuous(limits = c(0,30), breaks = scales::breaks_width(5)) +
  scale_y_continuous(limits = c(-3,3)) +
  scale_color_manual( breaks = c("1 Frame", "2 Frames", "5 Frames"), values = c("black", "red", "green")) +
  labs(x = "Normalized Intensity", y = "average change in intensity per frame")+
  theme_bw(base_size = 18) +
  theme(legend.position = "top", legend.title = element_blank())

ggplot(
  data = Growth,
  aes(
    x = ROUND_NORMALIZED_INTENSITY,
    y = NORM_DELTA_AVG_NORMALIZED_INTENSITY_1
  )
)+
  geom_path(
    color = "black"
  )+
  scale_x_continuous(
    limits = c(1,150)
  )+
  geom_path(
    aes(
      y = NORM_DELTA_AVG_NORMALIZED_INTENSITY_2),
    color = "red"
  )+
  geom_path(
    aes(
      y = NORM_DELTA_AVG_NORMALIZED_INTENSITY_5),
    color = "green"
    # )+
    # geom_path(
    #   aes(
    #     y = NORM_DELTA_AVG_NORMALIZED_INTENSITY_8),
    #   color = "orange"
  )+
  theme_bw()

ggplot(
  data = Growth,
  aes(
    x = ROUND_NORMALIZED_INTENSITY,
    y = DELTA_AVG_NORMALIZED_INTENSITY_1
  )
)+
  geom_path(
    aes(
      color = "1 Frame")
  )+
  geom_path(
    aes(
      y = DELTA_AVG_NORMALIZED_INTENSITY_2,
      color = "2 Frames")
  )+
  geom_path(
    aes(
      y = DELTA_AVG_NORMALIZED_INTENSITY_5,
      color = "5 Frames")
  )+
  # geom_path(
  #   aes(
  #     y = DELTA_AVG_NORMALIZED_INTENSITY_8),
  #   color = "orange"
  # )+
  scale_x_continuous(
    limits = c(0,60),
    breaks = scales::breaks_width(5)
  )+
  scale_y_continuous(
    limits = c(-3,6)
  )+
  scale_color_manual(
    breaks = c("1 Frame", "2 Frames", "5 Frames"),
    values = c("black", "red", "green")
  )+
  labs(
    x = "Normalized Intensity",
    y = "average change in intensity over indicated frames"
  )+
  theme_bw(
    base_size = 18
  )+
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )
