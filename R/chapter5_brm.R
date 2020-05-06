library(rethinking)
data("WaffleDivorce")

d = WaffleDivorce
rm(WaffleDivorce)

detach(package:rethinking, unload = T)
library(brms)
library(tidyverse)
library(magrittr)

d %<>%
  mutate(MedianAgeMarriage_s = (MedianAgeMarriage - mean(MedianAgeMarriage)) /
           sd(MedianAgeMarriage),
         Marriage_s = (Marriage - mean(Marriage)) / sd(Marriage))

# Multivariate regression model
b5.3 =
  brm(
    data = d,
    family = gaussian,
    formula = Divorce ~ 1 + Marriage_s + MedianAgeMarriage_s,
    prior = c(prior(normal(10, 10), class = Intercept),
              prior(normal(0, 1), class = b),
              prior(uniform(0, 10), class = sigma)),
    iter = 2000, warmup = 500, chains = 4, cores = 4,
    seed = 5
  )

print(b5.3)

stanplot(b5.3)

fitted(b5.3) %>%
  as_tibble() %>%
  bind_cols(d) %>%

  ggplot(aes(x = Divorce, y = Estimate)) +
  geom_abline(linetype = 2, color = "grey50", size = .5) +
  geom_point(size = 1.5, color = "firebrick4", alpha = 3/4) +
  geom_linerange(aes(ymin = Q2.5, ymax = Q97.5),
                 size = 1/4, color = "firebrick4") +
  geom_linerange(aes(ymin = Estimate - Est.Error,
                     ymax = Estimate + Est.Error),
                 size = 1/2, color = "firebrick4") +
  # Note our use of the dot placeholder, here: https://magrittr.tidyverse.org/reference/pipe.html
  geom_text(data = . %>% filter(Loc %in% c("ID", "UT")),
            aes(label = Loc),
            hjust = 0, nudge_x = - 0.65) +
  labs(x = "Observed divorce",
       y = "Predicted divorce") +
  theme_bw() +
  theme(panel.grid = element_blank())

# 5.2 Masked relationship ----
library(rethinking)
data(milk)
d <- milk

rm(milk)
detach(package:rethinking, unload = T)
library(brms)
library(tidyverse)

# pairwise correlations
d %>%
  select(kcal.per.g, mass, neocortex.perc) %>%
  pairs(col = "firebrick4")

dcc <-
  d %>%
  drop_na(ends_with("_s"))

b5.5 = brm(data = dcc, family = gaussian,
  formula = kcal.per.g ~ 1 + neocortex.perc,
  prior = c(prior(normal(0, 100), class = Intercept),
                      prior(normal(0, 1), class = b),
                      prior(cauchy(0, 1), class = sigma)),
    # prior(normal(0, 100), class = Intercept),
    # prior(normal(0, 1), class = b),
    # prior(cauchy(0, 1), class = sigma),
    iter = 2000, warmup = 500, chains = 4, cores = 4,
    seed = 5)

print(b5.5, digits = 3)
fixef(b5.5)[2] * (76 - 55)


nd <- tibble(neocortex.perc = 54:80)

fitted(b5.5,
       newdata = nd,
       probs = c(.025, .975, .25, .75)) %>%
  as_tibble() %>%
  bind_cols(nd) %>%

  ggplot(aes(x = neocortex.perc, y = Estimate)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5),
              fill = "firebrick", alpha = 1/5) +
  geom_smooth(aes(ymin = Q25, ymax = Q75),
              stat = "identity",
              fill = "firebrick4", color = "firebrick4", alpha = 1/5, size = 1/2) +
  geom_point(data = dcc,
             aes(y = kcal.per.g),
             size = 2, color = "firebrick4") +
  coord_cartesian(xlim = range(dcc$neocortex.perc),
                  ylim = range(dcc$kcal.per.g)) +
  ylab("kcal.per.g") +
  theme_bw() +
  theme(panel.grid = element_blank())


dcc = dcc %>% mutate(log_mass = log(mass))

b5.6 <-
  brm(data = dcc, family = gaussian,
      kcal.per.g ~ 1 + log_mass,
      prior = c(prior(normal(0, 100), class = Intercept),
                prior(normal(0, 1), class = b),
                prior(uniform(0, 1), class = sigma)),
      iter = 2000, warmup = 500, chains = 4, cores = 4,
      control = list(adapt_delta = 0.9),
      seed = 5)

print(b5.6, digits = 3)

nd <- tibble(log_mass = seq(from = -2.5, to = 5, length.out = 30))

fitted(b5.6,
       newdata = nd,
       probs = c(.025, .975, .25, .75)) %>%
  as_tibble() %>%
  bind_cols(nd) %>%

  ggplot(aes(x = log_mass, y = Estimate)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5),
              fill = "firebrick", alpha = 1/5) +
  geom_smooth(aes(ymin = Q25, ymax = Q75),
              stat = "identity",
              fill = "firebrick4", color = "firebrick4", alpha = 1/5, size = 1/2) +
  geom_point(data = dcc,
             aes(y = kcal.per.g),
             size = 2, color = "firebrick4") +
  coord_cartesian(xlim = range(dcc$log_mass),
                  ylim = range(dcc$kcal.per.g)) +
  ylab("kcal.per.g") +
  theme_bw() +
  theme(panel.grid = element_blank())

# Fit model with both variables

b5.7 = brm(data = dcc, family = gaussian,
           formula = kcal.per.g ~ 1 + neocortex.perc + log_mass,
           prior = c(prior(normal(0, 100), class = Intercept),
                     prior(normal(0, 1), class = b),
                     prior(cauchy(0, 1), class = sigma)
                     ),
           iter = 4000, warmup = 2000, chains = 4, cores = 4,
           control = list(adapt_delta = 0.999),
           seed = 5)

nd <-
  tibble(neocortex.perc = 54:80 %>% as.double(),
         log_mass       = mean(dcc$log_mass))

b5.7 %>%
  fitted(newdata = nd,
         probs = c(.025, .975, .25, .75)) %>%
  as_tibble() %>%
  bind_cols(nd) %>%

  ggplot(aes(x = neocortex.perc, y = Estimate)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5),
              fill = "firebrick", alpha = 1/5) +
  geom_smooth(aes(ymin = Q25, ymax = Q75),
              stat = "identity",
              fill = "firebrick4", color = "firebrick4", alpha = 1/5, size = 1/2) +
  geom_point(data = dcc,
             aes(y = kcal.per.g),
             size = 2, color = "firebrick4") +
  coord_cartesian(xlim = range(dcc$neocortex.perc),
                  ylim = range(dcc$kcal.per.g)) +
  ylab("kcal.per.g") +
  theme_bw() +
  theme(panel.grid = element_blank())
f

