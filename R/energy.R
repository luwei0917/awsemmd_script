library(tidyverse)
library(gridExtra)
library(stringr)
place <- "/Users/weilu/Research/server/project/aawsem/aawsemJan16/T0766_single/simulation_iteration_1/0"
pre <- "/Users/weilu/Research/server/project/aawsem/aawsemJan16/T0766_single/simulation_iteration_1/"
place <- str_c(pre, "13")
setwd(place)
data <- as_tibble(read.table("energy.log", header = TRUE))
data <- data[-c(3,6,12,14,15,16)]
data <- data[-1,]

test <- head(data)

d2 <- data %>% 
  gather(Chain, Chi, Rama, DSSP,P_AP, Water, Burial, Helix, Frag_Mem, VTotal, key = "part", value = "energy")

test
test %>% 
  gather(Chain, Chi, Rama, DSSP,P_AP, Water, Burial, Helix, Frag_Mem, VTotal, key = "part", value = "energy")

ggplot(d2) +
  aes(y=energy, x= part) +
  geom_boxplot()

read <- function(var) {
  pre <- "/Users/weilu/Research/server/project/aawsem/aawsemJan16/T0766_single/simulation_iteration_1/"
  place <- str_c(pre, var)
  setwd(place)
  data <- as_tibble(read.table("energy.log", header = TRUE))
  data <- data[-c(3,6,12,14,15,16)]
  data <- data[-1,]
}
data10 <- read(10)
d10 <- data10 %>% 
  gather(Chain, Chi, Rama, DSSP,P_AP, Water, Burial, Helix, Frag_Mem, VTotal, key = "part", value = "energy")


ggplot(d2) +
  aes(y=energy, x= part) +
  geom_boxplot() +
  geom_boxplot(data = d10, aes(y=energy, x= part), color = "Red", alpha = 0.1) +
  coord_flip()


# chi, chain the same
ggplot(d2) +
  aes(y=energy, x= part) +
  geom_boxplot() +
  geom_boxplot(data = d10, aes(y=energy, x= part), color = "Red", alpha = 0.1) +
  coord_flip() +
  ylim(0, 100)


ggplot(d2) +
  aes(y=energy, x= part) +
  geom_boxplot() +
  geom_boxplot(data = d10, aes(y=energy, x= part), color = "Red", alpha = 0.1) +
  coord_flip() +
  ylim(-1000, -500)


ggplot(d2) +
  aes(y=energy) +
  geom_boxplot(y = energy) +
  geom_boxplot(data = d10, aes(y=energy, x= part), color = "Red", alpha = 0.1) +
  coord_flip() +
  facet_wrap(~part)
