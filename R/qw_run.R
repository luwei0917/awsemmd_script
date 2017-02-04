library(tidyverse)
library(gridExtra)

# data <- read_csv("../../T0766/simulation_v2/data")
# data <- read_csv("../../T0766/simulation_iteration_1/data")
# data <- read_csv("../../T0766_single/simulation_iteration_1/data")
data <- read_csv("../../T0766_single/simulation/data")
# data <- read_csv("../../T0778_single/simulation/data")
# data <- read_csv("../../T0778/simulation_iteration_1/data")
# data <- read_csv("../../T0778/simulation/data")

ggplot(data) +
  aes(run, qw) +
  geom_boxplot()

data_end <- data %>% filter(step > 7000)
ggplot(data_end) +
  aes(reorder(run, qw, FUN = median), qw) +
  geom_boxplot()  +
  xlab("") +
  ylab("Qw") +
  coord_flip() +
  theme_bw() +
  theme(axis.text=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) + 
  theme(axis.title.y=element_text(size=30))

# qw run
ggplot(data_end) +
  aes(reorder(run, qw, FUN = median), qw) +
  geom_boxplot() +
  coord_flip()


ggplot(data_end) +
  aes(reorder(run, qw, FUN = median), qw) +
  geom_boxplot() +
  xlab("") +
  ylab("Qw") +
  ylim(0.35, 0.55) +
  coord_flip() +
  theme_bw() +
  theme(axis.text=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) + 
  theme(axis.title.y=element_text(size=30))

# energy run
data_energy <- data %>% filter(step > 7000)
ggplot(data_energy) +
  aes(reorder(run, -energy, FUN = median), energy) +
  geom_boxplot() +
  xlab("") +
  ylab("Energy") +
  coord_flip() +
  theme_bw() +
  theme(axis.text=element_text(size=30)) +
  theme(axis.title.x=element_text(size=30)) + 
  theme(axis.title.y=element_text(size=30))

ggplot(data_last) +
  aes(qw, energy, color=run) +
  geom_point()

ggplot(data %>% filter(step > 7000)) +
  aes(energy) + 
  geom_histogram() +
  facet_wrap(~run, scales = "free") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

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


ggsave("~/Desktop/T0766_round3_run_qw.png")
ggsave("~/Desktop/T0766_round3_run_qw.png")