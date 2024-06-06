library(tidyverse)
# read sims
sims <- read_tsv("autosomal_sims", col_names = F)
# add colnames
colnames(sims) <- c("simname","gfprop","P1","P2","P3","Z","D","f4")

# order by ascending geneflow 3
sims <- sims %>%
  arrange(gfprop) %>%
  # change prop3 to percentage
  mutate(gf_percent = gfprop * 100) %>%
  # set to significant if Z >= 3
  mutate(is_significant = if_else(Z >= 3, TRUE, FALSE))

# make gf_percent factor
sims$gf_percent<- factor(sims$gf_percent, levels = unique(sims$gf_percent))
# and also make significance factor
sims$is_significant <- factor(sims$is_significant, levels=c(FALSE,TRUE))

# plot dstats
p <- ggplot() +
  geom_boxplot(outlier.shape = NA, data = sims, aes(x = gf_percent, y = D)) +
  geom_jitter(data = sims, width = .15, height = 0, aes(x = gf_percent, y = D, shape = is_significant)) +
  scale_shape_manual(values= c(1,19)) +
  #geom_boxplot() +
  # ylim(c(0,.020)) +
  labs(y = "D-statistic", x = "Percentage of gene flow C. mitis -> C. denti") +
  theme_bw() +
  theme(legend.position = "none") +
  #geom_hline(yintercept = 3, color = "red", linetype = "dashed") +
  theme(plot.background = element_blank())


# add lines with observed dstatistcs with the different p1s
observed_d.mean <- 0.00318
observed_d.max <- 0.0120

p <- p + geom_hline(yintercept = observed_d.max, color = "red", linetype = "dashed")

# save plot
ggsave("autosomal_simulations_D.pdf",p, width = 6, height = 3)