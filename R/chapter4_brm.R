# Chapter 4 in brms

library(rethinking)
data("Howell1")
d = Howell1

rm(Howell1)
detach(package:rethinking, unload = T)
library(brms)
library(tidyverse)

d2 = d %>% filter(age >= 18)

ggplot(data = tibble(x = seq(from = 100, to = 250, by = .1)),
       aes(x = x, y = dnorm(x, mean = 178, sd = 20))) +
  geom_line() +
  ylab("density")

# Grid approximation of the posterior distribution
n = 200

d_grid =
  tibble(mu    = seq(from = 140, to = 160, length.out = n),
         sigma = seq(from = 4,   to = 9,   length.out = n)) %>%
  # we'll accomplish with `tidyr::expand()` what McElreath did with base R `expand.grid()`
  expand(mu, sigma)

head(d_grid)

grid_function = function(mu, sigma){
  dnorm(d2$height, mean = mu, sd = sigma, log = T) %>%
    sum()
}

d_grid = d_grid %>%
  mutate(log_likelihood = map2(mu, sigma, grid_function)) %>%
  unnest(cols = c(log_likelihood)) %>%
  mutate(prior_mu       = dnorm(mu,    mean = 178, sd  = 20, log = T),
         prior_sigma    = dunif(sigma, min  = 0,   max = 50, log = T)) %>%
  mutate(product        = log_likelihood + prior_mu + prior_sigma) %>%
  mutate(probability    = exp(product - max(product)))

d_grid %>%
  ggplot(aes(x = mu, y = sigma, z = probability)) +
  geom_contour() +
  labs(x = expression(mu),
       y = expression(sigma)) +
  coord_cartesian(xlim = range(d_grid$mu),
                  ylim = range(d_grid$sigma)) +
  theme(panel.grid = element_blank())

d_grid %>%
  ggplot(aes(x = mu, y = sigma)) +
  geom_raster(aes(fill = probability),
              interpolate = T) +
  scale_fill_viridis_c(option = "A") +
  labs(x = expression(mu),
       y = expression(sigma)) +
  theme(panel.grid = element_blank())

# Sampling from the posterior ---

set.seed(4)
d_grid_samples <-
  d_grid %>%
  sample_n(size = 1e4, replace = T, weight = probability)

d_grid_samples %>%
  ggplot(aes(x = mu, y = sigma)) +
  geom_point(size = .9, alpha = 1/15) +
  scale_fill_viridis_c() +
  labs(x = expression(mu[samples]),
       y = expression(sigma[samples])) +
  theme(panel.grid = element_blank())

d_grid_samples %>%
  select(mu, sigma) %>%
  gather() %>%

  ggplot(aes(x = value)) +
  geom_density(fill = "grey33", size = 0) +
  scale_y_continuous(NULL, breaks = NULL) +
  xlab(NULL) +
  theme(panel.grid = element_blank()) +
  facet_wrap(~key, scales = "free")

library(tidybayes)

d_grid_samples %>%
  select(mu, sigma) %>%
  gather() %>%
  group_by(key) %>%
  median_hdi(value)

# Fitting the model with brm ----

b4.1 =
  brm(
    data = d2,
    family = gaussian,
    height ~ 1,
    prior = c(prior(normal(178, 20), class = Intercept),
              prior(uniform(0, 50), class = sigma)),
    iter = 31000, warmup = 30000, chains = 4, cores = 4,
    seed = 4)

plot(b4.1)

# uniform takes long to "warm up"
b4.1_half_cauchy <-
  brm(data = d2, family = gaussian,
      height ~ 1,
      prior = c(prior(normal(178, 20), class = Intercept),
                prior(cauchy(0, 1), class = sigma)),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,
      seed = 4)

plot(b4.1_half_cauchy)

# Sampling from a brm fit ----

vcov(b4.1_half_cauchy)

post <- posterior_samples(b4.1_half_cauchy)

head(post)

select(post, b_Intercept:sigma) %>%
  cov()

select(post, b_Intercept:sigma) %>%
  cov() %>%
  diag()

post %>%
  select(-lp__) %>%
  gather(parameter) %>%
  group_by(parameter) %>%
  summarise(mean = mean(value),
            SD   = sd(value),
            `2.5_percentile`  = quantile(value, probs = .025),
            `97.5_percentile` = quantile(value, probs = .975)) %>%
  mutate_if(is.numeric, round, digits = 2)

posterior_summary(b4.1_half_cauchy)

post %>%
  select(-lp__) %>%
  gather(parameter) %>%
  group_by(parameter) %>%
  mean_qi(value)

# Fitting linear model ----

b4.3 =
  brm(data = d2, family = gaussian,
      height ~ 1 + weight,
      prior = c(
        prior(normal(156, 100), class = Intercept),
        prior(normal(0, 10), class = b),
        prior(cauchy(0, 1), class = sigma)),
        # prior(uniform(0, 50), class = sigma)),
      iter = 41000, warmup = 40000, chains = 4, cores = 4,
      seed = 4
      )

plot(b4.3)

posterior_summary(b4.3)[1:3, ]

posterior_samples(b4.3) %>%
  select(-lp__) %>%
  cor() %>%
  round(digits = 2)

# Centering reduces correlations

d2 <-
  d2 %>%
  mutate(weight_c = weight - mean(weight))

b4.4 <-
  brm(data = d2, family = gaussian,
      height ~ 1 + weight_c,
      prior = c(prior(normal(178, 100), class = Intercept),
                prior(normal(0, 10), class = b),
                prior(uniform(0, 50), class = sigma)),
      iter = 46000, warmup = 45000, chains = 4, cores = 4,
      seed = 4)

posterior_samples(b4.4) %>%
  select(-lp__) %>%
  cor() %>%
  round(digits = 2)

d2 %>%
  ggplot(aes(x = weight, y = height)) +
  geom_abline(intercept = fixef(b4.3)[1],
              slope     = fixef(b4.3)[2]) +
  geom_point(shape = 1, size = 2, color = "royalblue") +
  theme_bw() +
  theme(panel.grid = element_blank())

post <- posterior_samples(b4.3)
head(post)


mu_at_50 <-
  post %>%
  transmute(mu_at_50 = b_Intercept + b_weight * 50)

mu <- fitted(b4.3, summary = F)

weight_seq <- tibble(weight = seq(from = 25, to = 70, by = 1))

mu <-
  fitted(b4.3,
         summary = F,
         newdata = weight_seq) %>%
  as_tibble() %>%
  # here we name the columns after the `weight` values from which they were computed
  set_names(25:70) %>%
  mutate(iter = 1:n())

mu <-
  mu %>%
  gather(weight, height, -iter) %>%
  # We might reformat `weight` to numerals
  mutate(weight = as.numeric(weight))

head(mu)

d2 %>%
  ggplot(aes(x = weight, y = height)) +
  geom_point(data = mu %>% filter(iter < 101),
             alpha = .1)

mu_summary <-
  fitted(b4.3,
         newdata = weight_seq) %>%
  as_tibble() %>%
  # let's tack on the `weight` values from `weight_seq`
  bind_cols(weight_seq)

mu_summary

d2 %>%
  ggplot(aes(x = weight, y = height)) +
  geom_ribbon(data = mu_summary,
              aes(y = Estimate, ymin = Q2.5, ymax = Q97.5),
              alpha = 1/2) +
  geom_line(data = mu_summary,aes(y = Estimate)) +
  # geom_smooth(data = mu_summary,
  #             aes(y = Estimate, ymin = Q2.5, ymax = Q97.5),
  #             stat = "identity",
  #             fill = "grey70", color = "black", alpha = 1, size = 1/2) +
  geom_point(color = "navyblue", shape = 1, size = 1.5, alpha = 2/3) +
  coord_cartesian(xlim = range(d2$weight)) +
  theme(text = element_text(family = "Times"),
        panel.grid = element_blank())

# Incoporating uncertainty from sigma ----
# brms::predict() : rethinking::sim()

weight_seq <- tibble(weight = seq(from = 25, to = 70, by = 1))

pred_height <-
  predict(b4.3,
          newdata = weight_seq) %>%
  as_tibble() %>%
  bind_cols(weight_seq)

pred_height %>%
  slice(1:6)

d2 %>%
  ggplot(aes(x = weight)) +
  geom_ribbon(data = pred_height,
              aes(ymin = Q2.5, ymax = Q97.5),
              fill = "grey83") +
  geom_smooth(data = mu_summary,
              aes(y = Estimate, ymin = Q2.5, ymax = Q97.5),
              stat = "identity",
              fill = "grey70", color = "black", alpha = 1, size = 1/2) +
  geom_point(aes(y = height),
             color = "navyblue", shape = 1, size = 1.5, alpha = 2/3) +
  coord_cartesian(xlim = range(d2$weight),
                  ylim = range(d2$height)) +
  theme(text = element_text(family = "Times"),
        panel.grid = element_blank())
