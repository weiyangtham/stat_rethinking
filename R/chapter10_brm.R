library(tidyverse)
library(brms)
library(magrittr)

library(wesanderson)

wes_palette("Moonrise2")

library(ggthemes)
library(bayesplot)

theme_set(theme_default() +
            theme_tufte() +
            theme(plot.background = element_rect(fill  = wes_palette("Moonrise2")[3],
                                                 color = wes_palette("Moonrise2")[3])))
library(rethinking)


data("chimpanzees")
d <- chimpanzees

detach(package:rethinking, unload = T)
library(brms)
rm(chimpanzees)

b10.1 <- brm(data = d, family = binomial,
             formula = pulled_left | trials(1) ~ 1,
             prior = prior(normal(0, 10), class = Intercept),
             seed = 10)

fixef(b10.1) %>%
  round(digits = 2)

# convert from log-odds to probabilty
c(.18, .46) %>%
  inv_logit_scaled()

fixef(b10.1) %>%
  inv_logit_scaled()

b10.2 <- brm(data = d, family = binomial,
             formula = pulled_left | trials(1) ~ 1 + prosoc_left,
             prior = c(prior(normal(0, 10), class = Intercept),
                       prior(normal(0, 10), class = b)),
             seed = 10)

b10.3 <- brm(data = d, family = binomial,
             formula = pulled_left | trials(1) ~ 1 + prosoc_left + prosoc_left:condition,
             prior = c(prior(normal(0, 10), class = Intercept),
                       prior(normal(0, 10), class = b)),
             seed = 10)

b10.1 <- add_criterion(b10.1, "waic")
b10.2 <- add_criterion(b10.2, "waic")
b10.3 <- add_criterion(b10.3, "waic")

w <- loo_compare(b10.1, b10.2, b10.3, criterion = "waic")

print(w, simplify = F)

cbind(waic_diff = w[, 1] * -2,
      se        = w[, 2] *  2) %>%
  round(digits = 2)



ppa <-
  pp_average(b10.1, b10.2, b10.3,
             weights = "waic",
             method = "fitted") %>%
  as_tibble() %>%
  bind_cols(b10.3$data) %>%
  distinct(Estimate, Q2.5, Q97.5, condition, prosoc_left) %>%
  mutate(x_axis = str_c(prosoc_left, condition, sep = "/")) %>%
  mutate(x_axis = factor(x_axis, levels = c("0/0", "1/0", "0/1", "1/1"))) %>%
  rename(pulled_left = Estimate)

# the empirically-based summaries
d_plot <-
  d %>%
  group_by(actor, condition, prosoc_left) %>%
  summarise(pulled_left = mean(pulled_left)) %>%
  mutate(x_axis = str_c(prosoc_left, condition, sep = "/")) %>%
  mutate(x_axis = factor(x_axis, levels = c("0/0", "1/0", "0/1", "1/1")))

# the plot
ppa %>%
  ggplot(aes(x = x_axis)) +
  geom_smooth(aes(y = pulled_left, ymin = Q2.5, ymax = Q97.5, group = 0),
              stat = "identity",
              fill = wes_palette("Moonrise2")[2], color = "black",
              alpha = 1, size = 1/2) +
  geom_line(data = d_plot,
            aes(y = pulled_left, group = actor),
            color = wes_palette("Moonrise2")[1], size = 1/3) +
  scale_x_discrete(expand = c(.03, .03)) +
  coord_cartesian(ylim = 0:1) +
  labs(x = "prosoc_left/condition",
       y = "proportion pulled left") +
  theme(axis.ticks.x = element_blank())

# this helps us set our custom color scheme
color_scheme_set(c(wes_palette("Moonrise2")[3],
                   wes_palette("Moonrise2")[1],
                   wes_palette("Moonrise2")[2],
                   wes_palette("Moonrise2")[2],
                   wes_palette("Moonrise2")[1],
                   wes_palette("Moonrise2")[1]))

# the actual plot
mcmc_pairs(x = posterior_samples(b10.3),
           pars = c("b_Intercept", "b_prosoc_left", "b_prosoc_left:condition"),
           off_diag_args = list(size = 1/10, alpha = 1/6),
           diag_fun = "dens")

b10.4 <-
  brm(data = d, family = binomial,
      pulled_left | trials(1) ~ 0 + factor(actor) + prosoc_left + condition:prosoc_left ,
      prior(normal(0, 10), class = b),
      iter = 2500, warmup = 500, chains = 2, cores = 2,
      control = list(adapt_delta = 0.9),
      seed = 10)

print(b10.4)

post <- posterior_samples(b10.4)

post %>%
  glimpse()

post %>%
  ggplot(aes(x = b_factoractor2)) +
  geom_density(color = "transparent",
               fill = wes_palette("Moonrise2")[1]) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(x        = NULL,
       title    = "Actor 2's large and uncertain intercept",
       subtitle = "Once your log-odds are above, like, 4, it's all\npretty much a probability of 1.")

library(rethinking)
library(brms)
data(UCBadmit)

d = UCBadmit
detach(package:rethinking, unload = T)

d <- d %>% mutate(male = ifelse(applicant.gender == "male", 1, 0))


b10.6 <- brm(data = d, family = binomial,
             formula = admit | trials(applications) ~ 1 + male,
             prior = c(prior(normal(0, 10), class = Intercept),
                       prior(normal(0, 10), class = b)),
             iter = 2500, warmup = 500, cores = 2, chains = 2,
             seed = 10)

b10.7 <- brm(data = d, family = binomial,
             formula = admit | trials(applications) ~ 1,
             prior = c(prior(normal(0, 10), class = Intercept)),
             iter = 2500, warmup = 500, cores = 2, chains = 2,
             seed = 10)

b10.6 <- add_criterion(b10.6, "waic")
b10.7 <- add_criterion(b10.7, "waic")

w <- loo_compare(b10.6, b10.7, criterion = "waic")

w
print(w, simplify = F)

cbind(waic_diff = w[, 1] * -2,
      se        = w[, 2] *  2) %>%
  round(digits = 2)

model_weights(b10.6, b10.7, weights = "waic")

print(b10.6)

fixef(b10.6)[2] %>%
  exp() %>%
  round(digits = 2)

post <- posterior_samples(b10.6)

post %>%
  mutate(p_admit_male = inv_logit_scaled(b_Intercept + b_male),
         p_admit_female = inv_logit_scaled(b_Intercept),
         diff_admit     = p_admit_male - p_admit_female) %>%
  summarise(`2.5%`  = quantile(diff_admit, probs = .025),
            `50%`   = median(diff_admit),
            `97.5%` = quantile(diff_admit, probs = .975))


d <-
  d %>%
  mutate(case = factor(1:12))

p <-
  predict(b10.6) %>%
  as_tibble() %>%
  bind_cols(d)

d_text <-
  d %>%
  group_by(dept) %>%
  summarise(case  = mean(as.numeric(case)),
            admit = mean(admit / applications) + .05)

ggplot(data = d, aes(x = case, y = admit / applications)) +
  geom_pointrange(data = p,
                  aes(y    = Estimate / applications,
                      ymin = Q2.5     / applications ,
                      ymax = Q97.5    / applications),
                  color = wes_palette("Moonrise2")[1],
                  shape = 1, alpha = 1/3) +
  geom_point(color = wes_palette("Moonrise2")[2]) +
  geom_line(aes(group = dept),
            color = wes_palette("Moonrise2")[2]) +
  geom_text(data = d_text,
            aes(y = admit, label = dept),
            color = wes_palette("Moonrise2")[2],
            family = "serif") +
  coord_cartesian(ylim = 0:1) +
  labs(y     = "Proportion admitted",
       title = "Posterior validation check") +
  theme(axis.ticks.x = element_blank())

b10.8 <-
  brm(data = d, family = binomial,
      admit | trials(applications) ~ 0 + dept,
      prior(normal(0, 10), class = b),
      iter = 2500, warmup = 500, cores = 2, chains = 2,
      seed = 10)

b10.9 <-
  update(b10.8,
         newdata = d,
         formula = admit | trials(applications) ~ 0 + dept + male)

b10.8 <- add_criterion(b10.8, criterion = "waic")
b10.9 <- add_criterion(b10.9, criterion = "waic")

l_b10.6_reloo <- loo(b10.6, reloo = T)
l_b10.7_reloo <- loo(b10.7, reloo = T)
l_b10.8_reloo <- loo(b10.8, reloo = T)
l_b10.9_reloo <- loo(b10.9, reloo = T)


loo_compare(l_b10.6_reloo, l_b10.7_reloo, l_b10.8_reloo, l_b10.9_reloo)

model_weights(b10.6, b10.7, b10.8, b10.9,
              weights = "loo") %>%
  round(digits = 3)

fixef(b10.9) %>% round(digits = 2)

fixef(b10.9)[7, 1] %>% exp()

predict(b10.9) %>%
  as_tibble() %>%
  bind_cols(d) %>%

  ggplot(aes(x = case, y = admit / applications)) +
  geom_pointrange(aes(y    = Estimate / applications,
                      ymin = Q2.5     / applications ,
                      ymax = Q97.5    / applications),
                  color = wes_palette("Moonrise2")[1],
                  shape = 1, alpha = 1/3) +
  geom_point(color = wes_palette("Moonrise2")[2]) +
  geom_line(aes(group = dept),
            color = wes_palette("Moonrise2")[2]) +
  geom_text(data = d_text,
            aes(y = admit, label = dept),
            color = wes_palette("Moonrise2")[2],
            family = "serif") +
  coord_cartesian(ylim = 0:1) +
  labs(y     = "Proportion admitted",
       title = "Posterior validation check") +
  theme(axis.ticks.x = element_blank())

pairs(b10.9,
      off_diag_args = list(size = 1/10, alpha = 1/6))

# Poisson ----

library(rethinking)
data("Kline")
d <- Kline

detach(package:rethinking, unload = T)

head(d)

d %<>% mutate(log_pop = log(population),
              contact_high = contact == "high")

b10.10 <- brm(data = d, family = poisson,
              total_tools ~ 1 + log_pop + contact_high + contact_high:log_pop,
              prior = c(prior(normal(0, 100), class = Intercept),
                        prior(normal(0, 1), class = b)),
              iter = 3000, warmup = 1000, chains = 4, cores = 4,
              seed = 10)

b10.10

post <- posterior_samples(b10.10)

post %>%
  select(-lp__) %>%
  rename(b_interaction = `b_log_pop:contact_highTRUE`) %>%
  psych::lowerCor()

color_scheme_set(c(wes_palette("Moonrise2")[2],
                   wes_palette("Moonrise2")[1],
                   wes_palette("Moonrise2")[4],
                   wes_palette("Moonrise2")[2],
                   wes_palette("Moonrise2")[1],
                   wes_palette("Moonrise2")[1]))

post %>%
  select(-lp__) %>%
  rename(b_interaction = `b_log_pop:contact_highTRUE`) %>%

  mcmc_intervals(prob = .5, prob_outer = .95) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0))

post <-
  post %>%
  mutate(lambda_high = exp(b_Intercept + b_contact_highTRUE + (b_log_pop + `b_log_pop:contact_highTRUE`) * 8),
         lambda_low  = exp(b_Intercept + b_log_pop * 8)) %>%
  mutate(diff        = lambda_high - lambda_low)

post %>%
  summarise(sum = sum(diff > 0) / length(diff))

post %>%
  ggplot(aes(x = diff)) +
  geom_density(color = "transparent",
               fill = wes_palette("Moonrise2")[1]) +
  geom_vline(xintercept = 0, linetype = 2,
             color = wes_palette("Moonrise2")[2]) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(x = "lambda_high - lambda_low")

# no interaction
b10.11 <-
  update(b10.10, formula = total_tools ~ 1 + log_pop + contact_high)

# no contact rate
b10.12 <-
  update(b10.10, formula = total_tools ~ 1 + log_pop)

# no log-population
b10.13 <-
  update(b10.10, formula = total_tools ~ 1 + contact_high)

# intercept only
b10.14 <-
  update(b10.10, formula = total_tools ~ 1,
         seed = 10)

b10.10 <- add_criterion(b10.10, criterion = "waic")
b10.11 <- add_criterion(b10.11, criterion = "waic")
b10.12 <- add_criterion(b10.12, criterion = "waic")
b10.13 <- add_criterion(b10.13, criterion = "waic")
b10.14 <- add_criterion(b10.14, criterion = "waic")

w <- loo_compare(b10.10, b10.11, b10.12, b10.13, b10.14, criterion = "waic")

cbind(waic_diff = w[, 1] * -2,
      se        = w[, 2] *  2) %>%
  round(digits = 2)

model_weights(b10.10, b10.11, b10.12, b10.13, b10.14, weights = "waic") %>%
  round(digits = 2)

w %>%
  data.frame() %>%
  rownames_to_column(var = "model") %>%

  ggplot(aes(x = reorder(model, -waic),
             y    = waic,
             ymin = waic - se_waic,
             ymax = waic + se_waic,
             color = model)) +
  geom_pointrange(shape = 16, show.legend = F) +
  scale_color_manual(values = wes_palette("Moonrise2")[c(1, 2, 1, 1, 1)]) +
  coord_flip() +
  labs(x = NULL, y = NULL,
       title = "WAIC") +
  theme(axis.ticks.y    = element_blank())

nd <-
  tibble(contact_high = c(FALSE, TRUE)) %>%
  expand(contact_high,
         log_pop = seq(from = 6.5, to = 13, length.out = 50))

ppa <-
  pp_average(b10.10, b10.11, b10.12,
             weights = "loo",
             method  = "fitted",
             newdata = nd) %>%
  as_tibble() %>%
  bind_cols(nd)

ppa %>%
  ggplot(aes(x     = log_pop,
             group = contact_high)) +
  geom_smooth(aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
                  fill = contact_high, color = contact_high),
              stat = "identity",
              alpha = 1/4, size = 1/2) +
  geom_text(data = d,
            aes(y     = total_tools,
                label = total_tools,
                color = contact_high),
            size = 3.5) +
  coord_cartesian(xlim = c(7.1, 12.4),
                  ylim = c(12, 70)) +
  labs(x = "log population",
       y = "total tools",
       subtitle = "Blue is the high contact rate; black is the low.") +
  theme(legend.position = "none",
        panel.border    = element_blank())

# Center predictors
d <-
  d %>%
  mutate(log_pop_c = log_pop - mean(log_pop))

b10.10_c <-
  brm(data = d, family = poisson,
      total_tools ~ 1 + log_pop_c + contact_high + contact_high:log_pop_c,
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(normal(0, 1), class = b)),
      iter = 3000, warmup = 1000, chains = 4, cores = 4,
      seed = 10)

color_scheme_set(c(wes_palette("Moonrise2")[3],
                   wes_palette("Moonrise2")[1],
                   wes_palette("Moonrise2")[2],
                   wes_palette("Moonrise2")[2],
                   wes_palette("Moonrise2")[1],
                   wes_palette("Moonrise2")[1]))

# the actual plot
mcmc_pairs(x = posterior_samples(b10.10),
           pars = c("b_Intercept", "b_log_pop", "b_contact_highTRUE", "b_log_pop:contact_highTRUE"),
           off_diag_args = list(size = 1/10, alpha = 1/10),
           diag_fun = "dens")

mcmc_pairs(x = posterior_samples(b10.10_c),
           pars = c("b_Intercept", "b_log_pop_c",
                    "b_contact_highTRUE", "b_log_pop_c:contact_highTRUE"),
           off_diag_args = list(size = 1/10, alpha = 1/10),
           diag_fun = "dens")

set.seed(10)

num_days  <- 30
y         <- rpois(num_days, 1.5)

num_weeks <- 4
y_new     <- rpois(num_weeks, 0.5 * 7)

d <-
  tibble(y         = c(y, y_new),
         days      = c(rep(1, num_days), rep(7, num_weeks)),
         monastery = c(rep(0, num_days), rep(1, num_weeks))) %>%
  mutate(log_days  = log(days))

d %>% View()

# offset fixes the parameter value (1) for log_days

b10.15 <- brm(data = d, family = poisson,
              formula = y ~ 1 + offset(log_days) + monastery,
              prior = c(prior(normal(0, 100), class = Intercept),
                        prior(normal(0, 1), class = b)),
              iter = 2500, warmup = 500, cores = 2, chains = 2,
              seed = 10)

library(tidybayes)
posterior_samples(b10.15) %>%
  transmute(lambda_old = exp(b_Intercept),
            lambda_new = exp(b_Intercept + b_monastery)) %>%
  gather() %>%
  mutate(key = factor(key, levels = c("lambda_old", "lambda_new"))) %>%
  group_by(key) %>%
  mean_hdi(value, .width = .89) %>%
  mutate_if(is.double, round, digits = 2)



