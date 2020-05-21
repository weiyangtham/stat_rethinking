library(tidyverse)
library(rethinking)
library(brms)
library(ggthemes)


data("reedfrogs")
d <- reedfrogs

detach(package:rethinking, unload = T)

d <- d %>% mutate(tank = row_number())

b12.1 <- brm(data = d, family = binomial,
             formula = surv | trials(density) ~ 0 + factor(tank),
             prior = c(prior(normal(0, 5), class = b)),
             iter = 2000, warmup = 500, chains = 4, cores = 4,
             seed = 12)

b12.2 <- brm(data = d, family = binomial,
             formula = surv | trials(density) ~ 1 + (1 | tank),
             prior = c(prior(normal(0, 1), class = Intercept),
                       prior(cauchy(0, 1), class = sd)),
             iter = 4000, warmup = 1000, chains = 4, cores = 4,
             seed = 12)

b12.1 <- add_criterion(b12.1, "waic")
b12.2 <- add_criterion(b12.2, "waic")

w <- loo_compare(b12.1, b12.2, criterion = "waic")

print(w, simplify = F)

cbind(waic_diff = w[, 1] * -2,
      se        = w[, 2] *  2)

post <- posterior_samples(b12.2, add_chain = T)

post_mdn <-
  coef(b12.2, robust = T)$tank[, , ] %>%
  as_tibble() %>%
  bind_cols(d) %>%
  mutate(post_mdn = inv_logit_scaled(Estimate))

post_mdn

post_mdn %>%
  ggplot(aes(x = tank)) +
  geom_hline(yintercept = inv_logit_scaled(median(post$b_Intercept)), linetype = 2, size = 1/4) +
  geom_vline(xintercept = c(16.5, 32.5), size = 1/4) +
  geom_point(aes(y = propsurv), color = "orange2") +
  geom_point(aes(y = post_mdn), shape = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(breaks = c(1, 16, 32, 48)) +
  labs(title    = "Multilevel shrinkage!",
       subtitle = "The empirical proportions are in orange while the model-\nimplied proportions are the black circles. The dashed line is\nthe model-implied average survival proportion.") +
  annotate("text", x = c(8, 16 + 8, 32 + 8), y = 0,
           label = c("small tanks", "medium tanks", "large tanks")) +
  theme_fivethirtyeight() +
  theme(panel.grid = element_blank())

post <- as_tibble(post)

set.seed(12)
post %>%
  sample_n(100) %>%
  expand(nesting(iter, b_Intercept, sd_tank__Intercept),
         x = seq(from = -4, to = 5, length.out = 100)) %>%

  ggplot(aes(x = x, group = iter)) +
  geom_line(aes(y = dnorm(x, b_Intercept, sd_tank__Intercept)),
            alpha = .2, color = "orange2") +
  labs(title = "Population survival distribution",
       subtitle = "The Gaussians are on the log-odds scale.") +
  scale_y_continuous(NULL, breaks = NULL) +
  coord_cartesian(xlim = c(-3, 4)) +
  theme_fivethirtyeight() +
  theme(plot.title    = element_text(size = 13),
        plot.subtitle = element_text(size = 10))

ggplot(data = post %>% slice(1:8000),
       aes(x = rnorm(n    = 8000,
                     mean = b_Intercept,
                     sd   = sd_tank__Intercept) %>%
             inv_logit_scaled())) +
  geom_density(size = 0, fill = "orange2") +
  scale_y_continuous(NULL, breaks = NULL) +
  ggtitle("Probability of survival") +
  theme_fivethirtyeight()

# Simulate pond data ----


a <- 1.4 # grand mean
sigma <- 1.5 # std dev of grand mean
nponds <- 60

# 60 ponds with four possible sample sizes
ni <- as.integer(rep(c(5, 10, 25, 35), each = 15))

# simulate pond intercepts
set.seed(1059)
true_a = rnorm(nponds, mean = a, sd = sigma)

dsim <- tibble(pond = 1:nponds,
               ni = ni,
               true_a = true_a)

# Simulate survivors
set.seed(1059)
dsim <- dsim %>%
  mutate(si = rbinom(nponds, size = ni, prob = inv_logit_scaled(true_a)))

# No pooling
dsim <- dsim %>% mutate(p_nopool = si/ni)

# Partial pooling

b12.3 <- brm(data = dsim, family = binomial,
             formula = si | trials(ni) ~ 1 + (1 | pond),
             prior = c(prior(normal(0, 1), class = Intercept),
                       prior(cauchy(0, 1), class = sd)),
             iter = 10000, warmup = 1000, chains = 1, cores = 1,
             seed = 12)

print(b12.3)

coef(b12.3)$pond[c(1:2, 59:60), , ] %>%
  round(digits = 2)

rhat(b12.3)

neff_ratio(b12.3)

# Get effective sample
(n_iter <- (b12.3$fit@sim$iter - b12.3$fit@sim$warmup) * b12.3$fit@sim$chains)

neff_ratio(b12.3) %>%
  data.frame() %>%
  rownames_to_column() %>%
  set_names("parameter", "neff_ratio") %>%
  mutate(eff_sample = (neff_ratio * n_iter) %>% round(digits = 0)) %>%
  head()

# Calculate absolute error of estimates

p_partpool <-
  coef(b12.3)$pond[, , ] %>%
  as_tibble() %>%
  transmute(p_partpool = inv_logit_scaled(Estimate))

dsim <-
  dsim %>%
  bind_cols(p_partpool) %>%
  mutate(p_true         = inv_logit_scaled(true_a)) %>%
  mutate(nopool_error   = abs(p_nopool   - p_true),
         partpool_error = abs(p_partpool - p_true))

dfline <-
  dsim %>%
  select(ni, nopool_error:partpool_error) %>%
  gather(key, value, -ni) %>%
  group_by(key, ni) %>%
  summarise(mean_error = mean(value)) %>%
  mutate(x    = c( 1, 16, 31, 46),
         xend = c(15, 30, 45, 60))

dsim %>%
  ggplot(aes(x = pond)) +
  geom_vline(xintercept = c(15.5, 30.5, 45.4),
             color = "white", size = 2/3) +
  geom_point(aes(y = nopool_error), color = "orange2", size = 3) +
  geom_point(aes(y = partpool_error), shape = 1, size = 3) +
  geom_segment(data = dfline,
               aes(x = x, xend = xend,
                   y = mean_error, yend = mean_error),
               color = rep(c("orange2", "black"), each = 4),
               linetype = rep(1:2, each = 4)) +
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 50, 60)) +
  annotate("text", x = c(15 - 7.5, 30 - 7.5, 45 - 7.5, 60 - 7.5), y = .45,
           label = c("tiny (5)", "small (10)", "medium (25)", "large (35)")) +
  labs(y        = "absolute error",
       title    = "Estimate error by model type") +
  theme_fivethirtyeight() +
  theme(panel.grid    = element_blank(),
        plot.subtitle = element_text(size = 10))

# Multiple clusters ----

set.seed(123)
y1 <- rnorm(1e4, 10, 1)
y2 <- 10 + rnorm(1e4, 0, 1)

tibble(y = c(y1, y2), g = rep(c("a", "b"), each = 1e4)) %>%
  ggplot(aes(x = y, fill = g)) +
  geom_density(alpha = 1/2)

data("chimpanzees")
d <- chimpanzees

detach(package:rethinking, unload = T)

b12.4 <- brm(data = d, family = binomial,
             formula = pulled_left | trials(1) ~ 1 + (1 | actor) +
               prosoc_left + condition:prosoc_left,
             prior = c(prior(normal(0, 10), class = Intercept),
                       prior(normal(0, 10), class = b),
                       prior(cauchy(0, 1), class = sd)),
             # I'm using 4 cores, instead of the `cores=3` in McElreath's code
             iter = 5000, warmup = 1000, chains = 4, cores = 3,
             control = list(adapt_delta = 0.95),
             seed = 12)

post <- posterior_samples(b12.4) %>% as_tibble()

post %>%
  select(starts_with("r_actor")) %>%
  gather() %>%
  mutate(value = value + post$b_Intercept) %>%
  group_by(key) %>%
  summarise(mean = mean(value) %>% round(digits = 2))


b12.5 <- brm(data = d, family = binomial,
             formula = pulled_left | trials(1) ~ 1 + (1 | actor) + (1 | block) +
                                                 prosoc_left + condition:prosoc_left,
             prior = c(prior(normal(0, 10), class = Intercept),
                       prior(normal(0, 10), class = b),
                       prior(cauchy(0, 1), class = sd)),
             iter = 6000, warmup = 1000, cores = 3, chains = 4,
             control = list(adapt_delta = 0.99),
             seed = 12)

library(bayesplot)
color_scheme_set("orange")

post <- posterior_samples(b12.5, add_chain = T) %>% as_tibble()

post %>% select(-lp__, -iter) %>%
  mcmc_trace(facet_args = list(ncol = 4)) +
  scale_x_continuous(breaks = c(0, 2500, 5000)) +
  theme_fivethirtyeight() +
  theme(legend.position = c(.75, .06))

neff_ratio(b12.5) %>%
  mcmc_neff() +
  theme_fivethirtyeight()

print(b12.5)

ranef(b12.5)$actor[, , "Intercept"] %>%
  round(digits = 2)

ranef(b12.5)$block[, , "Intercept"] %>%
  round(digits = 2)

stanplot(b12.5, pars = c("^r_", "^b_", "^sd_")) +
  theme_fivethirtyeight() +
  theme(axis.text.y = element_text(hjust = 0))

post %>%
  select(starts_with("sd")) %>%
  gather() %>%
  ggplot(aes(value, fill = key)) +
  geom_density() +
  theme_fivethirtyeight()

b12.4 <- add_criterion(b12.4, "loo")
b12.5 <- add_criterion(b12.5, "loo")

loo_compare(b12.4, b12.5) %>%
  print(simplify = F)

