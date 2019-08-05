library(tidyverse)
library(magrittr)
library(rethinking)

# Fig 3.1 ----

grid.len = 1000
p_grid = seq(from = 0, to = 1, length.out = grid.len)

# Define prior
prior = rep(1, grid.len)

# compute likelihood at each value in grid
likelihood = dbinom(6, size  = 9, prob = p_grid)

# compute product of likelihood and prior
posterior = likelihood * prior

# standardize posterior to sum to 1
posterior = posterior/sum(posterior)

sum(posterior)

samples = sample(p_grid, prob = posterior, size = 1e4, replace = TRUE)

plot(samples)

tibble(x = 1:1e4, y = samples) %>%
  ggplot(aes(x, y)) +
  geom_point(alpha = 1/10)

dens(samples)

tibble(x = 1:1e4, y = samples) %>%
  ggplot(aes(y)) +
  geom_histogram()
  # geom_density()


# posterior probability that the proportion of water is less than 0.5
# using grid approx
sum(posterior[p_grid < 0.5])

# posterior probability that the proportion of water is less than 0.5
# using samples from the posterior
sum(samples < 0.5) / 1e4

# how much posterior prob lies between 0.5 and 0.75
# using samples from the posterior
sum(samples > 0.5 & samples < 0.75) / 1e4

# Samples are useful for finding boundaries of percentiles

quantile(samples, 0.8) # lower 80 percentile
quantile(samples, c(0.1, 0.9)) # middle 80 percentile

# Figure 3.3 ----

p_grid = seq(0, 1 , length.out = 1000)
prior = rep(1, 1000)
likelihood = dbinom(3, size = 3, prob = p_grid)
posterior = likelihood * prior
posterior = posterior / sum(posterior)
samples = sample(p_grid, size = 1e4, replace = TRUE, prob = posterior)

# middle 50 percentile
PI(samples, prob = 0.95)

# highest posterior density interval
# narrowest interval containing 50% probability mass
HPDI(samples, prob = 0.95)


# Point estimates - which one? ----

# mode
p_grid[which.max(posterior)]
chainmode(samples, adj = 0.01)

# mean or median
mean(samples)
median(samples)

# Loss function ----

# Absolute loss
# Your guess is proportional to absolute value of distance from the correct answer

# Suppose we guess 0.5
sum(posterior * abs(0.5 - p_grid))

# weighted average loss for every guess
loss = sapply(p_grid, function(d) sum(posterior * abs(d - p_grid)))

p_grid[which.min(loss)] # close to posterior median
median(samples)

# Dummy data ----

# if toss globe two times, chances of each outcome: 0 water, 1 water, or 2 water

dbinom(0:2, size = 2, prob = 0.7)

# sample from binomial distribution
rbinom(1, size = 2, prob = 0.7)
rbinom(10, size = 2, prob = 0.7)

dummy_w = rbinom(1e5, size = 2, prob = 0.7)

table(dummy_w)/1e5 # similar to analytic likelihood (dbinom)

# nine tosses instead of two tosses

dummy_w = rbinom(10, size = 1000, prob = 0.7)
simplehist(dummy_w, xlab = "dummy water count")

# Model checking ----

# simulate predicted observations for a single value of p
w = rbinom(1e4, size = 9, prob = 0.6)
simplehist(w)

# Replace 0.6 with samples from the posterior
w = rbinom(1e4, size = 9, prob = samples)

simplehist(w)
