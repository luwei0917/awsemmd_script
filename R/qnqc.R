library(tidyverse)
library(gridExtra)

# setwd("/Users/weilu/Research/server/project/freeEnergy_2xov/qValue_v2/simulation/350")
# setwd("/Users/weilu/Research/server/project/freeEnergy_2xov/pullingDistance/simulation")
read_table("data")
setwd("/Users/weilu/Research/server/project/freeEnergy_2xov/pullingDistance_v3/simulation")
data <- read_csv("data")
data <- read_delim("data", "\t", col_names = FALSE)
data <- data %>% mutate(qn = as.double(X1), qc = as.double(X2), n = 1:80000)

data <- data %>% gather(qn,qc, key="side",value ="value")

ggplot(data) +
  aes(n,value, color = side) +
  geom_point()

ggplot(data) +
  aes(x) +
  geom_histogram()
