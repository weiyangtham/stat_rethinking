# Figure 2.7 ----

# Define grid
grid.len = 20

p_grid = seq(from = 0, to = 1, length.out = grid.len)

# Define prior
prior = rep(1, grid.len)
prior = ifelse(p_grid < 0.5, 0, 1)
prior = exp(-5*abs(p_grid - 0.5))

# compute likelihood at each value in grid
likelihood = dbinom(6, size  = 9, prob = p_grid)

# compute product of likelihood and prior
unstd.posterior = likelihood * prior

# standardize posterior to sum to 1
posterior = unstd.posterior/sum(unstd.posterior)

plot(p_grid, posterior, type = "b")

# Quadratic appoximation ----

library(rethinking)

globe.qa = map(
  alist(
    w ~ dbinom(9, p), # binomial likelihood
    p ~ dunif(0, 1) # uniform prior
  ),
  data = list(w = 6))

precis(globe.qa) # parameters of posterior IF Gaussian

# analytical calculation to compare to quad approx
w = 6
n = 9

curve(dbeta(x, w+1, n-w+1), from = 0, to = 1) # analytical
curve(dnorm(x, 0.67, 0.16), lty = 2, add = TRUE) # quadratic approximation

# Practice ----

#' I found it easiest to solve the problem by drawing out a tree, starting
#' from whether the panda's species is A or B, and the possible sequences of births

# 2H1
0.166

# 2H2
0.33

# 2H3
0.36

# 2H4.1
0.69565

# 2H4.2
0.5625


# Statistical rethinking Winter 2019 homework ----

# Q1
grid.len = 20

p_grid = seq(from = 0, to = 1, length.out = grid.len)

# Define prior
prior = rep(1, grid.len)

# compute likelihood at each value in grid
likelihood = dbinom(8, size  = 15, prob = p_grid)

# compute product of likelihood and prior
unstd.posterior = likelihood * prior

# standardize posterior to sum to 1
posterior = unstd.posterior/sum(unstd.posterior)

plot(p_grid, posterior, type = "b")

# Q2H ("Hard")

grid.len = 20

p_grid = seq(from = 0, to = 1, length.out = grid.len)

# Define prior
prior = ifelse(p_grid < 0.5, 0, 1)

# compute likelihood at each value in grid
likelihood = dbinom(8, size  = 15, prob = p_grid)

# compute product of likelihood and prior
unstd.posterior = likelihood * prior

# standardize posterior to sum to 1
posterior2 = unstd.posterior/sum(unstd.posterior)

plot(p_grid, posterior2, type = "b")

tibble(type = rep(c("prior1", "prior2"), each = grid.len),
       x = rep(p_grid, 2),
       y = c(posterior, posterior2)) %>%
  ggplot(aes(x, y, colour = type)) +
  geom_point() +
  geom_line() +
  labs(title = "Flat prior vs narrower prior")

# Q3

posterior

