# Chapter 6, Statistical Rethinking Edition 1

library(rethinking)
data(milk)

library(tidyverse)

detach(package:rethinking, unload = T)
library(brms)

d <-
  milk %>%
  drop_na(ends_with("_s"))
rm(milk)

d <-
  d %>%
  mutate(neocortex = neocortex.perc / 100)

inits <- list(Intercept = mean(d$kcal.per.g),
              sigma     = sd(d$kcal.per.g))

inits_list <-list(inits, inits, inits, inits)

b6.11 <-
  brm(data = d, family = gaussian,
      kcal.per.g ~ 1,
      prior = c(prior(uniform(-1000, 1000), class = Intercept),
                prior(uniform(0, 100), class = sigma)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      inits = inits_list,
      seed = 6)

inits <- list(Intercept = mean(d$kcal.per.g),
              neocortex = 0,
              sigma     = sd(d$kcal.per.g))
inits_list <-list(inits, inits, inits, inits)

b6.12 <-
  brm(data = d, family = gaussian,
      kcal.per.g ~ 1 + neocortex,
      prior = c(prior(uniform(-1000, 1000), class = Intercept),
                prior(uniform(-1000, 1000), class = b),
                prior(uniform(0, 100), class = sigma)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      inits = inits_list,
      seed = 6)

inits <- list(Intercept   = mean(d$kcal.per.g),
              `log(mass)` = 0,
              sigma       = sd(d$kcal.per.g))
inits_list <-list(inits, inits, inits, inits)

b6.13 <-
  update(b6.12,
         newdata = d,
         formula = kcal.per.g ~ 1 + log(mass),
         inits   = inits_list)

inits <- list(Intercept   = mean(d$kcal.per.g),
              neocortex   = 0,
              `log(mass)` = 0,
              sigma       = sd(d$kcal.per.g))
inits_list <-list(inits, inits, inits, inits)

b6.14 <-
  update(b6.13,
         newdata = d,
         formula = kcal.per.g ~ 1 + neocortex + log(mass),
         inits   = inits_list,
      seed = 6)

waic(b6.14)

b6.11 <- add_criterion(b6.11, "waic")

b6.11$waic

# compute and save the WAIC information for the next three models
b6.12 <- add_criterion(b6.12, "waic")
b6.13 <- add_criterion(b6.13, "waic")
b6.14 <- add_criterion(b6.14, "waic")

# compare the WAIC estimates
w <- loo_compare(b6.11, b6.12, b6.13, b6.14,
                 criterion = "waic")

print(w)

cbind(waic_diff = w[, 1] * -2,
      se        = w[, 2] * 2)
