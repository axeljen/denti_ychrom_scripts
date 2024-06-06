rm(list = ls())
library(tidyverse)

# read the sims
sims <- read_csv("simple_yintro_model/simulations_merged.txt")

colnames(sims) <- c("cycle","introgressed_Y_freq",
                    "selection_coefficient","introprop",
                    "popsize","sexratio","simID")

# summarize the sims
simsum <- sims %>%
  group_by(selection_coefficient, introprop, simID) %>%
  summarize(ngen = max(cycle), maxfreq = max(introgressed_Y_freq),
            minfreq = min(introgressed_Y_freq),
            mingen = min(cycle)) %>%
  # if minfreq is 0 its lost
  mutate(outcome = ifelse(minfreq == 0, "lost", "segregating")) %>%
  # if it's one its fixed
  mutate(outcome = ifelse(maxfreq == 1, "fixed", outcome))

# count the occurrences for all different combinations of introprop and selection
catsum <- 
  simsum %>%
  group_by(introprop, selection_coefficient) %>%
  count(outcome)

# order by intro and selection
catsum <- catsum %>%
  arrange(introprop, selection_coefficient)

# make a factor of the two
catsum$introprop <- factor(catsum$introprop, levels = unique(catsum$introprop))
catsum$selection_coefficient <- factor(catsum$selection_coefficient, levels = unique(catsum$selection_coefficient))

# and a factor of the outcome
catsum$outcome <- factor(catsum$outcome, levels = c("segregating","fixed","lost"))
# make a barchart
propplot <- ggplot(catsum, aes(x = introprop, y = n,
                   fill = outcome)) +
  facet_wrap(vars(selection_coefficient), nrow = 5) +
  geom_col() +
  theme_bw() + 
  scale_fill_manual(values = c("#FFAE70","#6DB3D1","#B21531")) +
  # add labels with fixation count
  geom_text(data = catsum[which(catsum$outcome == "fixed"),],aes(label = n, y = 50)) +
  labs(y = "Simulation count", x = "Initial allele frequency") +
  theme(legend.position = "top")
# save
ggsave("slimsims_plot.pdf",
       propplot, width = 4, height = 6)