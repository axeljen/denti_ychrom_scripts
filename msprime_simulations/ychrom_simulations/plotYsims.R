library(tidyverse)

# read file
ysims <- read_tsv("msprime_ychrom_simulations")
ysims
# order by geneflow prop 3
ysims <- ysims %>%
  arrange(GeneFlowProp) %>%
  # change number of fixation events to percentage values
  mutate(percent_fixed = (Count / 1000) * 100) %>%
  # change prop to percentage
  mutate(percentage_geneflow = GeneFlowProp * 100) 

# make geneflow prop3 to factor
ysims$percentage_geneflow <- factor(ysims$percentage_geneflow, levels = unique(ysims$percentage_geneflow))

# and plot it
py <- ggplot(ysims, aes(x = percentage_geneflow, y = percent_fixed)) +
  #geom_boxplot() +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .15, height = 0, alpha = .3) +
  labs(y = "denti/mitis monophyly (%)", x = "Percentage of gene flow C. mitis -> C. denti") +
  theme_bw() +
  theme(plot.background = element_blank())

ggsave("ysims_plot.pdf",
	   py, width = 6, height = 3)