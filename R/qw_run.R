library(tidyverse)
library(gridExtra)
library(stringr)
# data <- read_csv("../../T0766/simulation_v2/data")
# data <- read_csv("../../T0766/simulation_iteration_1/data")
# data <- read_csv("../../T0766_single/simulation_iteration_1/data")
data <- read_csv("../../T0766_single/simulation/data")
# data <- read_csv("../../T0778_single/simulation/data")
# data <- read_csv("../../T0778/simulation_iteration_1/data")
# data <- read_csv("../../T0778/simulation/data")

setwd("/Users/weilu/Research/server/project/freeEnergy_2xov/qValue_v2/simulation/350")
setwd("/Users/weilu/Research/server/project/freeEnergy_2xov/qValue_v3/simulation/350")

setwd("/Users/weilu/Research/server/project/freeEnergy_2xov/qValue_v3/simulation/350/0")
setwd("/Users/weilu/Research/server/project/freeEnergy_2xov/2xov_folding_temperature/frag0.2/simulation")

setwd("/Users/weilu/Research/server/project/freeEnergy_2xov/2xov_folding_temperature/frag0.1/simulation")

setwd("/Users/weilu/Research/server/project/freeEnergy_2xov/2xov_folding_temperature/frag0.3/simulation")
data1 <- read_csv("data")
data2 <- read_csv("data")
data3 <- read_csv("data")
data <- as_tibble( read.table("energy.dat", header = TRUE) )
data <- data[-Shake]

data
ggsave("~/Desktop/pulling_data/frag0.01_qw_run.png")
ggsave("~/Desktop/pulling_data/frag0.01_qw_temp.png")

var = c(0.2, 0.1, 0.01)
pre = "/Users/weilu/Research/server/project/freeEnergy_2xov/2xov_folding_temperature/"
dir = str_c(pre, var, "/simulation/data")

dat <- function(var){
  pre = "/Users/weilu/Research/server/project/freeEnergy_2xov/2xov_folding_temperature/frag"
  dir = str_c(pre, var, "/simulation/data")
  data <- read_csv(dir)
  data %>% mutate(frag = var) %>% mutate(temp = 200 + step/(2*10^4))
}

set <- function(var){
  name <- str_c("frag", var)
  assign(name, value = dat(var))
}
var <- list("0.2","0.1","0.01")
var %>% 
  walk(str_c("frag",.),value = dat(.), assign)

frag0.2 <- dat(0.2)
frag0.1 <- dat(0.1)
frag0.01 <- dat(0.01)



data <- rbind(frag0.2,frag0.1,frag0.01)
ggplot(data)+
  aes(x=temp, qw) +
  geom_point()+
  geom_smooth()+
  ylim(0.2, 1) +
  facet_wrap(~frag) +
  theme_bw() +
  theme(axis.text=element_text(size=20)) +
  theme(axis.title.x=element_text(size=30)) + 
  theme(axis.title.y=element_text(size=30))

ggplot(data)+
  aes(x=temp, qw, color=run) +
  geom_point()+
  geom_smooth()+
  ylim(0.2, 1) +
  facet_wrap(~frag)+
  theme_bw() +
  theme(axis.text=element_text(size=20)) +
  theme(axis.title.x=element_text(size=30)) + 
  theme(axis.title.y=element_text(size=30))
ggsave("~/Desktop/pulling_data/frag_folding_color.png", width = 14, height =4)

plots <- data %>% 
  split(.$frag) %>% 
  map(~ggplot(., aes(temp, qw)) + geom_point()+ geom_smooth())

folder <- "~/Desktop/pulling_data/frag"
mylist = c(0.2, 0.1, 0.01)
paths <- stringr::str_c(folder, mylist, ".png")
pwalk(list(paths, plots), ggsave)



ggplot(data) +
  aes(reorder(run, qw, FUN = median), qw) +
  geom_boxplot()  +
  xlab("") +
  ylab("Qw") +
  coord_flip() +
  theme_bw()

data_end <- data %>% filter(step > 7000)
ggplot(data) +
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