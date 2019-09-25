# Statistical Rethinking Winter 2019 HW 2
# https://github.com/rmcelreath/statrethinking_winter2019/blob/master/homework/week02.pdf



library(tidyverse)
library(magrittr)
library(rethinking)

# Question 1 ----

# repeat steps around p105 of Stat Rethinking Edition 1

data("Howell1")

d = Howell1

head(d)

# Build a model

mdl_hw1 = map(
  flist = alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*weight,
    a ~ dnorm(178,100),
    b ~ dnorm(0, 10),
    sigma ~ dunif(0, 50)
  ),
  data = d)

precis(mdl_hw1)

# Vector of weights the HW wants us to predict
w = c(45, 40, 65, 31, 53)

d_hw = data.frame(weight = w)

mu = link(mdl_hw1, data = d_hw)

# 1000 by 5 matrix
str(mu)

colnames(mu) = as.character(w)
predictions = as_tibble(mu) %>%
  gather(weight, height) %>%
  mutate(weight = as.numeric(weight))

predictions %>%
  ggplot(aes(weight, height)) +
  geom_point(alpha = 1/100)

predictions %>%
  group_by(weight) %>%
  summarise(height = mean(height)) %>%
  ggplot(aes(weight, height)) +
  geom_point()

apply(mu, 2, HPDI, prob = 0.89)

# Question 2 ----

d %<>%
  mutate(weight_log = log(weight),
         height_log = log(height),
         weight_log_c = weight_log - mean(weight_log))

head(d)

ggplot(d) +
  aes(weight_log_c, height_log) +
  geom_point()

ggplot(d) +
  aes(weight_log) +
  geom_histogram()

ggplot(d) +
  aes(height_log) +
  geom_histogram()

mdl_hw2 = map(
  flist = alist(
    height_log ~ dnorm(mu, sigma),
    mu <- a + b*weight_log_c,
    a ~ dnorm(0, 10),
    b ~ dnorm(1, 0.2),
    sigma ~ dunif(0, 50)
  ),
  data = d)

precis(mdl_hw2)

cov2cor(vcov(mdl_hw2))

# Plot posterior predictions against raw data
post = extract.samples(mdl_hw2, n = 20)

d %>%
  ggplot(aes(weight_log_c, height_log)) +
  geom_point() +
  geom_abline(data = post, aes(intercept = a, slope = b),
              alpha = 0.1)

