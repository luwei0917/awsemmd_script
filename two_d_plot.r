library(tidyverse)
data <- read_table("pmf-350.dat", skip =1)
data <- data[-c(1,6,7,8)]

ggplot(data)+
  aes(bin_center_1,bin_center_2,color=f)+
  geom_point()
