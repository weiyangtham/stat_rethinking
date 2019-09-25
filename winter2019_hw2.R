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



