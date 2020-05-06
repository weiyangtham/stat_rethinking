# Chapter 7

library(rethinking)
data(rugged)
d = rugged
detach(package:rethinking, unload = T)

library(tidyverse)
library(brms)

# make the log version of criterion
d <-
  d %>%
  mutate(log_gdp = log(rgdppc_2000))

# extract countries with GDP data
dd <-
  d %>%
  filter(complete.cases(rgdppc_2000))

# split the data into countries in Africa and not in Africa
d.A1 <-
  dd %>%
  filter(cont_africa == 1)

d.A0 <-
  dd %>%
  filter(cont_africa == 0)

b7.1 <-
  brm(data = d.A1, family = gaussian,
      formula = log_gdp ~ 1 + rugged,
      prior = c(
        prior(normal(8, 100), class = Intercept),
        prior(normal(0, 1), class = b)
      ),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      seed = 7)

b7.2 <- update(b7.1, newdata = d.A0)

d.A1 %>%
  bind_cols(data = fitted(b7.1) %>% as_tibble()) %>%
  as_tibble() %>%
  ggplot(aes(rugged, log_gdp)) +
  geom_point() +
  geom_line(aes(y = Estimate))  +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.5)

d.A0 %>%
  bind_cols(data = fitted(b7.2) %>% as_tibble()) %>%
  as_tibble() %>%
  ggplot(aes(rugged, log_gdp)) +
  geom_point() +
  geom_line(aes(y = Estimate))  +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.5)


b7.3 <-
  update(b7.1,
         newdata = dd)

b7.4 <-
  update(b7.3,
         newdata = dd,
         formula = log_gdp ~ 1 + rugged + cont_africa)

b7.3 <- add_criterion(b7.3, c("loo", "waic"))
b7.4 <- add_criterion(b7.4, c("loo", "waic"))

b7.5 <-
  update(b7.4,
         formula = log_gdp ~ 1 + rugged*cont_africa)

b7.5 <- add_criterion(b7.5, c("loo", "waic"))

l <- loo_compare(b7.3, b7.4, b7.5,
                 criterion = "loo")

print(l, simplify = F)

# Plotting interaction

dd %<>% bind_cols(fitted(b7.5) %>% as_tibble()) %>%
  as_tibble()

dd %>%
  ggplot(aes(rugged, log_gdp, colour = factor(cont_africa))) +
  geom_point() +
  geom_line(aes(y = Estimate))  +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.5) +
  facet_wrap(~cont_africa) +
  ggthemes::scale_colour_pander()

# Incorporate uncertainty of interaction coefficents ----

posterior_summary(b7.5)

# Incorporate uncertainty by working with posterior samples

post <- posterior_samples(b7.5)

post %>%
  transmute(gamma_Africa    = b_rugged + `b_rugged:cont_africa`,
            gamma_notAfrica = b_rugged) %>%
  gather(key, value) %>%
  group_by(key) %>%
  summarise(mean = mean(value))

post %>%
  transmute(gamma_Africa    = b_rugged + `b_rugged:cont_africa`,
            gamma_notAfrica = b_rugged) %>%
  gather(key, value) %>%

  ggplot(aes(x = value, group = key, color = key, fill = key)) +
  geom_density(alpha = 1/4) +
  ggthemes::scale_colour_pander() +
  ggthemes::scale_fill_pander() +
  scale_x_continuous(expression(gamma), expand = c(0, 0)) +
  scale_y_continuous(NULL, breaks = NULL) +
  ggtitle("Terraine Ruggedness slopes",
          subtitle = "Blue = African nations, Green = others") +
  theme(text = element_text(family = "Times"),
        legend.position = "none")

post %>%
  mutate(gamma_Africa    = b_rugged + `b_rugged:cont_africa`,
         gamma_notAfrica = b_rugged) %>%
  mutate(diff            = gamma_Africa -gamma_notAfrica) %>%
  summarise(p = sum(diff < 0)/length(diff))

# Continuous interactions

library(rethinking)
data(tulips)
d <- tulips
str(d)

detach(package:rethinking, unload = T)

b7.6 <- brm(
  data = d, family = gaussian,
  formula = blooms ~ 1 + water + shade,
  prior = c(
    prior(normal(0, 100), class = Intercept),
    prior(normal(0, 100), class = b),
    prior(uniform(0, 100), class = sigma)
  ),
  iter = 2000, warmup = 1000, cores = 4, chains = 4,
  seed = 7)

b7.7 <- update(b7.6,
               formula = blooms ~ 1 + water + shade + water:shade)

# Try better priors, adapt_delta

b7.6 <-
  update(b7.6,
         prior = c(prior(normal(0, 100), class = Intercept),
                   prior(normal(0, 100), class = b),
                   prior(cauchy(0, 10), class = sigma)),
         control = list(adapt_delta = 0.9),
         seed = 7)

b7.7 <-
  update(b7.6,
         formula = blooms ~ 1 + water + shade + water:shade)

posterior_summary(b7.6)
posterior_summary(b7.7)

b7.6 <- add_criterion(b7.6, "waic")
b7.7 <- add_criterion(b7.7, "waic")


w <- loo_compare(b7.6, b7.7, criterion = "waic")

print(w, simplify = F)

cbind(waic_diff = w[, 1] * -2,
      se        = w[, 2] *  2)

model_weights(b7.6, b7.7, weights = "waic")














