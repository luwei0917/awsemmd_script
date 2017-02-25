library(tidyverse)
library(gridExtra)
library(stringr)
library(xml2)
data <- read_csv("~/Downloads/Awards.csv")
data <- read_xml("~/Downloads/2016/1600011.xml")
setwd("~/Desktop/feb19/")
## Example
# we generate some random plot
require(seqLog)
## the first plot is taken from the seqLogo help ( ?seqLogo )
## I selected this example on purpose because the seqLogo function is based on the grid graphics
## and is coded in such a way that doesn't allow the use of the par() function
mFile <- system.file("Exfiles/pwm1", package="seqLogo")
m <- read.table(mFile)
pwm <- makePWM(m)
png("seqLogo1.png", width=400, height=400)
seqLogo(pwm)
dev.off()
## totally unrelated
png("plot2.png", width=400, height=400)
plot(density(2*rnorm(1000)))
dev.off()


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


ggplot(d2%>% filter(Step > 7000)) +
  aes(energy, x = "1") +
  geom_boxplot() +
  geom_boxplot(data = d10 %>% filter(Step > 7000), aes(energy, x = "1"), color = "Red", alpha = 0.1) +
  facet_wrap(~part, scales="free") + 
  xlab("Red is run 10, black is run 13")

ggsave("~/Desktop/aawsem/each_energy.png")
ggplot(d2) +
  aes(y=energy) +
  geom_boxplot() +
  geom_boxplot(data = d10, aes(y = energy), color = "Red", alpha = 0.1) +
  coord_flip() +
  facet_wrap(~part)
