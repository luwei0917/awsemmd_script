library(tidyverse)
library(gridExtra)
library(stringr)
ggsave("~/Desktop/feb14/force_dis_all.png")
ggsave("~/Desktop/feb14/force_dis.png", width = 20, height = 10)
# setwd("/Users/weilu/Research/server/project/freeEnergy_2xov/qValue_v2/simulation/350")
# setwd("/Users/weilu/Research/server/project/freeEnergy_2xov/pullingDistance/simulation")
read_table("data")
# setwd("/Users/weilu/Research/server/project/freeEnergy_2xov/pullingDistance_v3/simulation")
pre <- "/Users/weilu/Research/server/project/pulling_2xov/feb09/"
place <- str_c(pre, "")
setwd(place)
data <- read_csv("data")
data_c <- read_csv("data_cont")

data <- rbind(data, data_c)
data <- data %>% mutate(force = step * 1e-7)
data_c <- data_c %>% mutate(force = step * 1e-7)

data %>% 
  ggplot() +
  aes(dis, force) +
  geom_point(color = "grey70") +
  geom_line(data = data %>% filter(run == "run_0"), aes(dis, force), color = "blue")+
  xlab("Distance (Angstrom)") +
  ylab("Force (Kcal/mole-Angstrom)") +
  theme(axis.text=element_text(size=20)) +
  theme(axis.title.x=element_text(size=30)) + 
  theme(axis.title.y=element_text(size=30))



data <- data %>% mutate(step = step / 1000000)
ggplot(data) +
  geom_line(aes(step, qn)) +
  geom_line(aes(step, qc), color = "Red") +
  facet_wrap(~run) +
  ylab("Q") +
  xlab("Steps in millions") +
  theme(axis.text=element_text(size=20)) +
  theme(axis.title.x=element_text(size=30)) + 
  theme(axis.title.y=element_text(size=30))

ggsave("~/Desktop/feb14/step_qnqc.png", width = 25, height  = 15)
ggsave("~/Desktop/feb14/step_qnqc.png", width = 7.9, height = 8.42)


# ggtitle("Red is Qc, black is Qn")  +
#   theme(plot.title = element_text(hjust = 0.5))+
ggplot(data) +
  aes(dis, force) +
  geom_line() +
  facet_wrap(~run) +
  xlab("Distance (Angstrom)") +
  ylab("Force (Kcal/mole-Angstrom)") +
  theme(axis.text=element_text(size=20)) +
  theme(axis.title.x=element_text(size=30)) + 
  theme(axis.title.y=element_text(size=30))
  
force_line <- data %>%
  group_by(force) %>%
  summarise(average = mean(dis))
  
data %>% 
  ggplot() +
    aes(dis, force) +
    geom_point() +
    geom_smooth(method = "glm")


data %>%
  mutate(bins = cut(data$dis, breaks = 100 )) %>% 
  ggplot() +
  aes(bins, force) +
  stat_summary(fun.y = "mean", geom = "point")

data %>% 
  mutate(bins = cut(data$dis, breaks =4 )) %>% 
  ggplot() +
  aes(bins, force) +
  geom_histogram()

# 
# dis_line <- data %>% 
#   group_by(run, dis) %>% 
#   summarise(average = mean(force)) 
# ggplot(dis_line) +
#   aes(dis, average, color = run) +
#   geom_line() 


data %>% filter(step > 6000000) %>% 
  ggplot() +
  aes(dis, force) +
  geom_point(color = "black") +
  geom_point(data = data_c %>% filter(step < 8000000),
             aes(dis, force), color = "red", alpha = 0.2) +
  facet_wrap(~run)


data %>% filter(step > 6000000) %>% 
  ggplot() +
    aes(dis, force, color = run) +
    geom_line() +
    geom_point(data = data_c, aes(dis, force, color = run))


ggplot(data) +
  aes(step, force, color = run) +
  geom_line()


data <- data %>% mutate(qn = as.double(X1), qc = as.double(X2), n = 1:80000)
data <- data %>% gather(qn,qc, key="side",value ="value")


