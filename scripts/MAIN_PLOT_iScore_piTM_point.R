MAIN_PLOT_POINT <- JE %>% 
  ggplot(aes(iScore, piTM, col = FILE)) +
  geom_abline(col = "gray") +
  #geom_label_repel(aes(label = FILE_MODEL, alpha = 0.3)) +
  geom_label_repel(aes(label = MODEL, alpha = 0.3)) +
  geom_point(size = 3) +
  scale_x_continuous(name = "iScore", breaks = c(0, 0.4, 0.5, 0.7, 1)) +
  scale_y_continuous(name = "piTM", breaks = c(0, 0.5, 1)) +
  theme(legend.position = "right") +
  scale_color_discrete(name = "Control") +
  ggtitle(paste0('Computational screening of ', FOLDER)) +
  expand_limits(x=c(0,1), y=c(0,1)) +
  annotate("text", x = 0.4, y = -0.05, label = "medium \n confidence") +
  annotate("text", x = 0.5, y = -0.05, label = "high \n confidence") +
  annotate("text", x = 0.7, y = -0.05, label = "very high \n confidence")

