library(tidyverse)
library(gridExtra)
library(stringr)
ggsave("~/Desktop/feb14/force_dis_all.png")
ggsave("~/Desktop/feb14/force_dis.png", width = 20, height = 10)
ggsave("~/Desktop/feb14/force_dis_all.png", width = 20, height = 10)
ggsave("~/Desktop/feb19/unfolded.png")
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
data <- data %>% mutate(force = step * 1e-7* 69.477)  # force in units of pN
data %>% filter() %>% 
  ggplot() +
  aes(dis, force) +
  geom_point(color = "grey70") +
  geom_line(data = data %>% filter(run == "run_0"), aes(dis, force), color = "blue")+
  xlab("Distance (Angstrom)") +
  ylab("Force (pN)") +
  theme(axis.text=element_text(size=20)) +
  theme(axis.title.x=element_text(size=30)) + 
  theme(axis.title.y=element_text(size=30))

data %>% filter(step %% 100000 == 0) %>% mutate(jump = c(0,diff(dis))) %>% 
  filter(jump > 25) %>% 
  ggplot() +
  geom_point(data = data  ,aes(dis, force) ,color = "grey70") +
  geom_point(aes(dis, force, size = jump) ) 

data %>% filter(step %% 100000 == 0) %>% group_by(run)%>% mutate(jump = c(0,diff(dis))) %>% 
  filter(jump > 25 & dis < 250) %>%
  ggplot() +
  geom_point(data = data  %>% filter(step %% 100000 == 0)  ,aes(dis, force) ,color = "grey70") +
  geom_point(aes(dis, force, size = jump, color = run) ) 

options(tibble.print_max = 20, tibble.print_min = 100)
data %>% filter(step %% 100000 == 0) %>% group_by(run) %>% mutate(jump = c(0,diff(dis))) %>% 
  ungroup(run)

data %>% filter(step %% 100000 == 0 & dis < 280 ) %>% group_by(run) %>% mutate(jump = c(0,diff(dis))) %>% 
  filter(jump > 20) %>% ungroup(run) %>% 
  ggplot() +
  aes(force) +
  geom_histogram(bins = 10)

data %>% filter(step %% 100000 == 0 & dis < 280 ) %>% group_by(run) %>% mutate(jump = c(0,diff(dis))) %>% 
  filter(jump > 20) %>% summarise(m = first(force)) %>% 
  ggplot() +
  aes(m) +
  geom_histogram(bins = 10)

data %>% filter(step %% 100000 == 0 & dis < 280 ) %>% group_by(run) %>% mutate(jump = c(0,diff(dis))) %>% 
  filter(jump > 25) %>% summarise(m = first(force))

data %>% filter(step %% 100000 == 0 & dis < 280 ) %>% group_by(run) %>% mutate(jump = c(0,diff(dis))) %>% 
  filter(jump > 20) %>% ungroup(run) %>%
  ggplot() +
  geom_point(data = data  %>% filter(step %% 100000 == 0)  ,aes(dis, force) ,color = "grey70") +
  geom_point(aes(dis, force, size = jump, color = run) ) 
# --------------bowie's method--------------------------
unfolded_fraction <- function(data, f) {
  a <- data %>% filter(dis < 150 & between(force, f-1, f+1)) %>% group_by(run) %>% 
    summarise(n())
  dim(a)[[1]]/20.0
}
unfolded_fraction(data, 80)

f <- seq(20,60,2)
output <- vector("double", 0)
for (i in seq_along(f)) {
  output <- c(output, unfolded_fraction(data,f[[i]]))
}
result <- tibble(x = f, y = output)

m <- nls(y ~ exp(-a/b*(exp(x*b) -1 )), result, start=list(a=1000,b=1000))
m
result <- result %>% mutate(re = predict(m))
ggplot(result)+
  aes(x,1-y)+
  geom_point(size = 5)+
  geom_line(aes(x,1-re), color = "blue") +
  xlab("Force (pN)") +
  ylab("Unfolded fraction") +
  theme(axis.text=element_text(size=20)) +
  theme(axis.title.x=element_text(size=30)) + 
  theme(axis.title.y=element_text(size=30))

data %>% group_by(force, run) %>% summarise(dis)
data %>% filter(dis < 120) %>%  
  ggplot() +
  aes(force) +
  geom_histogram()

# --------------bowie's method end--------------------------
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

data %>% 
  ggplot() +
  aes(dis, force) +
  geom_line() +
  facet_wrap(~run) + 
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


