library(tidyverse)
library(gridExtra)
library(stringr)
# mutate(run = as.integer(gsub("^[^0-9]*(\\d*)", "\\1", run))) %>% 

# test %>% separate(run, into=c("prefix","suffix")) %>% 
#   arrange(as.numeric(suffix)) %>%
#   mutate(new_test=apply( . , 1 , paste , collapse = "_" ))

# mutate(step1 = as.numeric( flatten(str_split(run, "_"))[c(FALSE, TRUE)]))
place <- "/Users/weilu/Research/server/project/freeEnergy_2xov/qValue_v4/simulation"
setwd(place)
data300 <- read_csv("300/data")
data350 <- read_csv("350/data")
data350 <- data350 %>% mutate(temp = "350")
data300 <- data300 %>% mutate(temp = "300")
data <- rbind(data350, data300)
data <- data %>% arrange(run, step)
data <- data %>%  mutate(run = as.character(run))

data %>% filter(step > 2000000) %>% 
  ggplot() +
  aes(seq_along(qw), qw) +
  geom_point(alpha = 0.2, color = "blue") +
  geom_line(aes(seq_along(qw), target))

data %>% filter(step > 2000000) %>% 
  ggplot(aes(energy, color = run)) +
  geom_density() +
  facet_wrap(~temp, ncol = 1, scales = "free")

data %>% filter(step > 2000000) %>% 
  ggplot(aes(energy)) +
  geom_density()

data350 %>% filter(step > 2000000) %>% 
  ggplot(aes(energy, color = run)) +
  geom_density()

data %>% filter(step > 2000000) %>% 
  ggplot(aes(qw, color = run)) +
  geom_density()+
  facet_wrap(~temp, ncol = 1, scales = "free")

data <- rbind(data300, data350)
data <- data %>% arrange(run, step)
data
test <- tibble(run = c("run_1", "run_11","run_2", "run_3", "run_11", "run_111", "run_4"))
test %>% arrange(run)
test %>% arrange(!desc(run)) 

mixedsort(rev(test[["run"]]))

test %>% arrange(mixedsort(rev(run)))
test %>% separate(run, into=c("prefix","suffix")) %>% 
  arrange(as.numeric(suffix)) %>%
  mutate(new_test=apply( . , 1 , paste , collapse = "_" ))

library(gtools)
test %>%  mutate(num = as.integer(gsub("^[^0-9]*(\\d*)", "\\1", run)))

test <- tibble(run = c("run_1","run_2", "run_3", "run_11", "run_111", "run_4"))
test[mixedorder(test$run),]

ggsave("~/Desktop/feb14/qvalue_qw_run.png")
