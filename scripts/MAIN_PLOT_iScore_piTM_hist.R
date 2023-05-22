MAIN_PLOT_HIST <- JE_all %>% #group_by(FILE) %>% filter(iScore == max(iScore)) %>% 
  ggplot(aes(iScore, fill = FILE)) +
  geom_vline(xintercept = 0.4, col = "gray40", linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = 0.5, col = "cornflowerblue", linetype = "dotted", linewidth = 1) +
  geom_vline(xintercept = 0.7, col = "lightgreen", linetype = "dotted", linewidth = 1) +
  geom_histogram(bins = 30) +
  expand_limits(x=c(0,1)) +
  theme(legend.position = "none") +
  annotate("text", x = 0.4, y = -0.05, label = "medium \n confidence") +
  annotate("text", x = 0.5, y = -0.05, label = "high \n confidence") +
  annotate("text", x = 0.7, y = -0.05, label = "very high \n confidence") +
  ggtitle("Model Confidence assessed by iScore") +
  facet_wrap(~FILE)

