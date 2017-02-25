ggsave("~/Desktop/feb19/rate_force.png")

pre <- "/Users/weilu/Research/data/pulling/"
target <- str_c(pre, "p_v2")
target2 <- str_c(pre, "p2_v1")
data2 <- data %>% mutate(run = str_c(run , "_2"))


load_and_transform <- function(name, rate){
  pre <- "/Users/weilu/Research/data/pulling/"
  target <- str_c(pre, name)
  data <- read_csv(target)
  data <- data %>% mutate(force = step * rate * 69.477) 
}

dis_force <- function(data) {
  data %>% filter(step %% 10000 == 0) %>% 
    ggplot() +
    aes(dis, force) +
    geom_point(color = "grey70") +
    xlab("Distance (Angstrom)") +
    ylab("Force (pN)") +
    theme(axis.text=element_text(size=20)) +
    theme(axis.title.x=element_text(size=30)) + 
    theme(axis.title.y=element_text(size=30))
}
p2_v1 <- load_and_transform("p2_v1", 2e-7)
dis_force(p2_v1)

binf <- p2_v1 %>%  filter(step %% 100000 == 0 & dis < 280 ) %>% group_by(run) %>% mutate(jump = c(0,diff(dis))) %>% 
  filter(jump > 25) %>% summarise(m = first(force))
tmp <- hist(binf$m, breaks = 8)

unfolded_fraction <- function(data, f) {
  a <- data %>% filter(dis < 140 & between(force, f, f+5)) %>% group_by(run) %>% 
    summarise(n())
  dim(a)[[1]]
}

k <- tmp$counts / result$y / 5
f <- seq(25,55,5)
output <- vector("double", 0)
for (i in seq_along(f)) {
  output <- c(output, unfolded_fraction(p2_v1,f[[i]]))
}
result <- tibble(x = f, y = output)

data2 <- data2 %>% mutate(force = step * 1e-7* 69.477)  # force in units of pN

f_check <- d_total%>% filter(step %% 100000 == 0 & dis < 280 ) %>% group_by(run) %>% mutate(jump = c(0,diff(dis))) %>% 
  filter(jump > 25) %>% summarise(m = first(force))

b <- hist(f_check$m)
fv1 <- data %>% filter(step %% 100000 == 0 & dis < 280 ) %>% 
  group_by(run) %>% mutate(jump = c(0,diff(dis))) %>% 
  filter(jump > 25) %>% summarise(m = first(force)) 

f_total <- rbind(fv1,fv2)

k <- b$counts / result$y / 5

names(k_data) <- c("x", "y")
ggplot(k_data)+
  aes(x,y)+
  geom_point()+
  scale_y_log10()


m <- nls(y ~ k * (1- a*x)*exp(b*(1-(1- a*x)^2 )), k_data, start=list(k= 0.0036,a=0.0024,b=20))
m

test <- function(x){
  exp(80*(1- (1 - 0.0025*x)^2 ))
}
test(30)

compare <- function(k_data, k ,a, b) {
  experiment <- function(x,k,a,b){
    k* (1 - a*x)*exp(b*(1-(1- a*x)^2 ))
  }
  xx <- seq(0,80,5)
  yy <- experiment(xx,k,a,b)
  exp_data <- tibble( x = xx, y = yy)
  ggplot()+
    geom_point(data = k_data, aes(x,y))+
    geom_line(data = exp_data, aes(x,y))+
    scale_y_log10()
}

compare(k_data, 0.0043, 0.0006, 68)

compare(k_data, 0.002,0.0025,20)
compare(k_data, 0.0036,0.0024,20) +
  xlab("Force (pN)") +
  ylab("Rate (s^-1)") +
  theme(axis.text=element_text(size=20)) +
  theme(axis.title.x=element_text(size=30)) + 
  theme(axis.title.y=element_text(size=30))

experiment <- function(x,k,a,b){
  k* (1 - a*x)*exp(b*(1-(1- a*x)^2 ))
}

yy <- experiment(seq(25,55,5),k,a,b)
exp_data <- tibble( x = seq(25,55,5), y = yy)
ggplot()+
  geom_point(data = k_data, aes(x,y))+
  geom_line(data = exp_data, aes(x,y))+
  scale_y_log10()






pre <- "/Users/weilu/Research/data/pulling/"
target <- str_c(pre, "p_v2")
target2 <- str_c(pre, "p2_v1")
data <- read_csv(target)
p2_v1 <- read_csv(target2)
data2 <- data %>% mutate(run = str_c(run , "_2"))
d_total <- rbind(data, data2)

load_and_transform <- function(name, rate){
  pre <- "/Users/weilu/Research/data/pulling/"
  target <- str_c(pre, name)
  data <- read_csv(target)
  data <- data %>% mutate(force = step * rate * 69.477) 
}

dis_force <- function(data) {
  data %>% filter(step %% 10000 == 0) %>% 
    ggplot() +
    aes(dis, force) +
    geom_point(color = "grey70") +
    xlab("Distance (Angstrom)") +
    ylab("Force (pN)") +
    theme(axis.text=element_text(size=20)) +
    theme(axis.title.x=element_text(size=30)) + 
    theme(axis.title.y=element_text(size=30))
}
p2_v1 <- load_and_transform("p2_v1", 2e-7)

dis_force(p2_v1)


binf <- p2_v1 %>%  filter(step %% 100000 == 0 & dis < 280 ) %>% group_by(run) %>% mutate(jump = c(0,diff(dis))) %>% 
  filter(jump > 25) %>% summarise(m = first(force))
tmp <- hist(binf$m, breaks = 8)

unfolded_fraction <- function(data, f) {
  a <- data %>% filter(dis < 140 & between(force, f, f+5)) %>% group_by(run) %>% 
    summarise(n())
  dim(a)[[1]]
}

k <- tmp$counts / result$y / 5
f <- seq(25,55,5)
output <- vector("double", 0)
for (i in seq_along(f)) {
  output <- c(output, unfolded_fraction(p2_v1,f[[i]]))
}
result <- tibble(x = f, y = output)

data2 <- data2 %>% mutate(force = step * 1e-7* 69.477)  # force in units of pN

f_check <- d_total%>% filter(step %% 100000 == 0 & dis < 280 ) %>% group_by(run) %>% mutate(jump = c(0,diff(dis))) %>% 
  filter(jump > 25) %>% summarise(m = first(force))

b <- hist(f_check$m)
fv1 <- data %>% filter(step %% 100000 == 0 & dis < 280 ) %>% 
  group_by(run) %>% mutate(jump = c(0,diff(dis))) %>% 
  filter(jump > 25) %>% summarise(m = first(force)) 

f_total <- rbind(fv1,fv2)

k <- b$counts / result$y / 5

names(k_data) <- c("x", "y")
ggplot(k_data)+
  aes(x,y)+
  geom_point()+
  scale_y_log10()


m <- nls(y ~ k * (1- a*x)*exp(b*(1-(1- a*x)^2 )), k_data, start=list(k= 0.0036,a=0.0024,b=20))
m

test <- function(x){
  exp(80*(1- (1 - 0.0025*x)^2 ))
}
test(30)

compare <- function(k_data, k ,a, b) {
  experiment <- function(x,k,a,b){
    k* (1 - a*x)*exp(b*(1-(1- a*x)^2 ))
  }
  xx <- seq(0,80,5)
  yy <- experiment(xx,k,a,b)
  exp_data <- tibble( x = xx, y = yy)
  ggplot()+
    geom_point(data = k_data, aes(x,y))+
    geom_line(data = exp_data, aes(x,y))+
    scale_y_log10()
}

compare(k_data, 0.0043, 0.0006, 68)

compare(k_data, 0.002,0.0025,20)
compare(k_data, 0.0036,0.0024,20) +
  xlab("Force (pN)") +
  ylab("Rate (s^-1)") +
  theme(axis.text=element_text(size=20)) +
  theme(axis.title.x=element_text(size=30)) + 
  theme(axis.title.y=element_text(size=30))

experiment <- function(x,k,a,b){
  k* (1 - a*x)*exp(b*(1-(1- a*x)^2 ))
}

yy <- experiment(seq(25,55,5),k,a,b)
exp_data <- tibble( x = seq(25,55,5), y = yy)
ggplot()+
  geom_point(data = k_data, aes(x,y))+
  geom_line(data = exp_data, aes(x,y))+
  scale_y_log10()


models <- tibble(
  a1 = runif(250, -20, 40),
  a2 = runif(250, -5, 5)
)
f_check %>%  ggplot() +
  aes(m) +
  geom_histogram(bins = 8)

ggplot(sim1, aes(x, y)) + 
  geom_abline(aes(intercept = a1, slope = a2), data = models, alpha = 1/4) +
  geom_point() 
sim1 <- result
model1 <- function(a, data) {
  exp(a[1]/a[2]*(exp(-data$x*a[2]) -1 ))
}
model1(c(7, 1.5), sim1)

measure_distance <- function(mod, data) {
  diff <- data$y - model1(mod, data)
  sqrt(mean(diff ^ 2))
}
sim1_dist <- function(a1, a2) {
  measure_distance(c(a1, a2), sim1)
}

models <- tibble(
  a1 = runif(2500, 0.5, 0575),
  a2 = runif(2500, 0.65, 0.75)
)



models <- models %>% 
  mutate(dist = purrr::map2_dbl(a1, a2, sim1_dist))
models

ggplot(models, aes(a1, a2)) +
  geom_point(data = filter(models, rank(dist) <= 10), size = 4, colour = "red") +
  geom_point(aes(colour = -dist))

result %>% mutate(re = exp(0.50225/0.7272*(exp(-x*0.7272) -1 )))

ggplot(result)+
  aes(x,1-y)+
  geom_point(size = 5)+
  geom_line(x = result$x, y = model1(c(0.50225, 0.7272),result))+
  xlab("Force (pN)") +
  ylab("Unfolded fraction") +
  theme(axis.text=element_text(size=20)) +
  theme(axis.title.x=element_text(size=30)) + 
  theme(axis.title.y=element_text(size=30))
