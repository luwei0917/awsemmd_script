models <- tibble(
  a1 = runif(250, -20, 40),
  a2 = runif(250, -5, 5)
)

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
  a1 = runif(250, 0, 1),
  a2 = runif(250, 0, 1)
)



models <- models %>% 
  mutate(dist = purrr::map2_dbl(a1, a2, sim1_dist))
models

ggplot(models, aes(a1, a2)) +
  geom_point(data = filter(models, rank(dist) <= 10), size = 4, colour = "red") +
  geom_point(aes(colour = -dist))
