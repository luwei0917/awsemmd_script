library(tidyverse)
library(gridExtra)

# data <- read_csv("../../T0766/simulation_v2/data")
# data <- read_csv("../../T0766/simulation_iteration_1/data")
# data <- read_csv("../../T0778_single/simulation/data")
# data <- read_csv("../../T0778/simulation_iteration_1/data")
data <- read_csv("../../T0778/simulation/data")
ggplot(data) +
  aes(run, qw) +
  geom_boxplot()

data_end <- data %>% filter(step > 7000)
ggplot(data_end) +
  aes(reorder(run, qw, FUN = median), qw) +
  geom_boxplot() +
  xlab("") +
  ylab("Qw") +
  ylim(0.2, 0.8) +
  coord_flip() +
  theme_bw() +
  theme(axis.text=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) + 
  theme(axis.title.y=element_text(size=30))


p1 <- ggplot(r0) +
  aes(qw) + 
  geom_histogram()
summary(r0)
r0 <- data %>% filter(run == "run_11", step > 7000)
ggplot(r0) +
  aes(qw) + 
  geom_histogram()
summary(r0)
grid.arrange(p1,p2)
