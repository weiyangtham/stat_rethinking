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

# Posterior predictions for existing clusters ----

chimp <- 2

d <- as_tibble(d)

d

nd <- tibble(prosoc_left = c(1, 0, 1, 0),
             condition = c(0, 0, 1, 1),
             actor = chimp)

fitted(b12.4, newdata = nd)

chimp_2_fitted <-
  fitted(b12.4,
         newdata = nd) %>%
  as_tibble() %>%
  mutate(condition = factor(c("0/0", "1/0", "0/1", "1/1"),
                            levels = c("0/0", "1/0", "0/1", "1/1")))

chimp_2_fitted

(
  chimp_2_d <-
    d %>%
    filter(actor == chimp) %>%
    group_by(prosoc_left, condition) %>%
    summarise(prob = mean(pulled_left)) %>%
    ungroup() %>%
    mutate(condition = str_c(prosoc_left, "/", condition)) %>%
    mutate(condition = factor(condition, levels = c("0/0", "1/0", "0/1", "1/1")))
)

chimp_2_fitted %>%
  # if you want to use `geom_line()` or `geom_ribbon()` with a factor on the x axis,
  # you need to code something like `group = 1` in `aes()`
  ggplot(aes(x = condition, y = Estimate, group = 1)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "orange1") +
  geom_line(color = "blue") +
  geom_point(data = chimp_2_d,
             aes(y = prob),
             color = "grey25") +
  ggtitle("Chimp #2",
          subtitle = "The posterior mean and 95%\nintervals are the blue line\nand orange band, respectively.\nThe empirical means are\nthe charcoal dots.") +
  coord_cartesian(ylim = c(.75, 1)) +
  theme_fivethirtyeight() +
  theme(plot.subtitle = element_text(size = 10))

post <- posterior_samples(b12.4) %>% as_tibble()

glimpse(post)

post

post %>% ggplot(aes(`r_actor[5,Intercept]`)) + geom_density()



chimp_5_fitted <- post %>%
  uncount(4) %>%
  bind_cols(uncount(nd %>% select(-actor), nrow(post))) %>%
  mutate(fitted = (b_Intercept +
                        `r_actor[5,Intercept]` +
                        b_prosoc_left * prosoc_left +
                        `b_prosoc_left:condition` * prosoc_left * condition) %>%
              inv_logit_scaled())

chimp_5_fitted <- chimp_5_fitted %>%
  mutate(condition = str_c(prosoc_left, "/", condition)) %>%
  mutate(condition = factor(condition, levels = c("0/0", "1/0", "0/1", "1/1")))

chimp_5_fitted %>%
  group_by(condition) %>%
  tidybayes::mean_qi(fitted)

chimp_5_my_fitted %>%
  ggplot(aes(x = condition, y = fitted, group = 1)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = "orange1") +
  geom_line(color = "blue") +
  geom_point(data = chimp_5_d,
             aes(y = prob),
             color = "grey25") +
  ggtitle("Chimp #5",
          subtitle = "This plot is like the last except\nwe did more by hand.") +
  coord_cartesian(ylim = 0:1) +
  theme_fivethirtyeight() +
  theme(plot.subtitle = element_text(size = 10))

# Posterior prediction for new clusters ----

post_average_actor <-
  post %>%
  # here we use the linear regression formula to get the log_odds for the 4 conditions
  transmute(`0/0` = b_Intercept,
            `1/0` = b_Intercept + b_prosoc_left,
            `0/1` = b_Intercept,
            `1/1` = b_Intercept + b_prosoc_left + `b_prosoc_left:condition`) %>%
  # with `mutate_all()` we can convert the estimates to probabilities in one fell swoop
  mutate_all(inv_logit_scaled) %>%
  # putting the data in the long format and grouping by condition (i.e., `key`)
  gather() %>%
  mutate(key = factor(key, level = c("0/0", "1/0", "0/1", "1/1"))) %>%
  group_by(key) %>%
  # here we get the summary values for the plot
  summarise(m  = mean(value),
            # note we're using 80% intervals
            ll = quantile(value, probs = .1),
            ul = quantile(value, probs = .9))

post_average_actor

p1 <-
  post_average_actor %>%
  ggplot(aes(x = key, y = m, group = 1)) +
  geom_ribbon(aes(ymin = ll, ymax = ul), fill = "orange1") +
  geom_line(color = "blue") +
  ggtitle("Average actor") +
  coord_cartesian(ylim = 0:1) +
  theme_fivethirtyeight() +
  theme(plot.title = element_text(size = 14, hjust = .5))

p1

set.seed(12.42)
ran_ef <-
  tibble(random_effect = rnorm(n = 1000, mean = 0, sd = post$sd_actor__Intercept)) %>%
  # uncount(4)
  bind_rows(., ., ., .)

fix_ef <-
  post %>%
  slice(1:1000) %>%
  transmute(`0/0` = b_Intercept,
            `1/0` = b_Intercept + b_prosoc_left,
            `0/1` = b_Intercept,
            `1/1` = b_Intercept + b_prosoc_left + `b_prosoc_left:condition`) %>%
  gather() %>%
  rename(condition    = key,
         fixed_effect = value) %>%
  mutate(condition = factor(condition, level = c("0/0", "1/0", "0/1", "1/1")))

ran_and_fix_ef <-
  bind_cols(ran_ef, fix_ef) %>%
  mutate(intercept = fixed_effect + random_effect) %>%
  mutate(prob      = inv_logit_scaled(intercept))

(
  marginal_effects <-
    ran_and_fix_ef %>%
    group_by(condition) %>%
    summarise(m  = mean(prob),
              ll = quantile(prob, probs = .1),
              ul = quantile(prob, probs = .9))
)

p2 <-
  marginal_effects %>%
  ggplot(aes(x = condition, y = m, group = 1)) +
  geom_ribbon(aes(ymin = ll, ymax = ul), fill = "orange1") +
  geom_line(color = "blue") +
  ggtitle("Marginal of actor") +
  coord_cartesian(ylim = 0:1) +
  theme_fivethirtyeight() +
  theme(plot.title = element_text(size = 14, hjust = .5))

p2

p3 <-
  ran_and_fix_ef %>%
  mutate(iter = rep(1:1000, times = 4)) %>%
  filter(iter %in% c(1:50)) %>%

  ggplot(aes(x = condition, y = prob, group = iter)) +
  theme_fivethirtyeight() +
  ggtitle("50 simulated actors") +
  coord_cartesian(ylim = 0:1) +
  geom_line(alpha = 1/2, color = "orange3") +
  theme(plot.title = element_text(size = 14, hjust = .5))

p3

# using fitted() to fit new clusters ----
nd <-
  tibble(prosoc_left = c(0, 1, 0, 1),
         condition   = c(0, 0, 1, 1))

(
  f <-
    fitted(b12.4,
           newdata = nd,
           re_formula = NA,
           probs = c(.1, .9)) %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(condition = str_c(prosoc_left, "/", condition) %>%
             factor(., levels = c("0/0", "1/0", "0/1", "1/1")))
)

p4 <-
  f %>%
  ggplot(aes(x = condition, y = Estimate, group = 1)) +
  geom_ribbon(aes(ymin = Q10, ymax = Q90), fill = "blue") +
  geom_line(color = "orange1") +
  ggtitle("Average actor") +
  coord_cartesian(ylim = 0:1) +
  theme_fivethirtyeight() +
  theme(plot.title = element_text(size = 14, hjust = .5))

p4

(
  f <-
    fitted(b12.4,
           newdata = nd,
           probs = c(.1, .9),
           allow_new_levels = T,
           sample_new_levels = "gaussian") %>%
    as_tibble() %>%
    bind_cols(nd) %>%
    mutate(condition = str_c(prosoc_left, "/", condition) %>%
             factor(., levels = c("0/0", "1/0", "0/1", "1/1")))
)

p5 <-
  f %>%
  ggplot(aes(x = condition, y = Estimate, group = 1)) +
  geom_ribbon(aes(ymin = Q10, ymax = Q90), fill = "blue") +
  geom_line(color = "orange1") +
  ggtitle("Marginal of actor") +
  coord_cartesian(ylim = 0:1) +
  theme_fivethirtyeight() +
  theme(plot.title = element_text(size = 14, hjust = .5))

p5

n_sim <- 50

(
  f <-
    fitted(b12.4,
           newdata = nd,
           # probs = c(.1, .9),
           allow_new_levels = T,
           sample_new_levels = "gaussian",
           summary = F,
           nsamples = n_sim) %>%
    as_tibble() %>%
    mutate(iter = 1:n_sim) %>%
    gather(key, value, -iter) %>%
    bind_cols(nd %>%
                transmute(condition = str_c(prosoc_left, "/", condition) %>%
                            factor(., levels = c("0/0", "1/0", "0/1", "1/1"))) %>%
                expand(condition, iter = 1:n_sim))
)

p6 <-
  f %>%

  ggplot(aes(x = condition, y = value, group = iter)) +
  geom_line(alpha = 1/2, color = "blue") +
  ggtitle("50 simulated actors") +
  coord_cartesian(ylim = 0:1) +
  theme_fivethirtyeight() +
  theme(plot.title = element_text(size = 14, hjust = .5))

p6

# Focus ----

data(Kline)
detach(package:rethinking, unload = T)

d <- Kline

b12.6 <- brm(formula = total_tools ~ 0 + Intercept + (1 | culture) + log(population), data = d,
             prior = c(prior(normal(0, 10), class = b, coef = Intercept),
                       prior(normal(0, 1), class = b),
                       prior(cauchy(0, 1), class = sd)),
             family = poisson(),
             iter = 4000, warmup = 1000, cores = 3, chains = 3,
             seed = 12)


nd <-
  tibble(population = seq(from = 1000, to = 400000, by = 5000),
         # to "simulate counterfactual societies, using the hyper-parameters" (p. 383),
         # we'll plug a new island into the `culture` variable
         culture    = "my_island")

set.seed(1038)
d_pred <- fitted(b12.6,
       newdata = nd,
       re_formula = ~(1 | culture), # want the random effects
       probs = c(.015, .055, .165, .835, .945, .985),
       allow_new_levels = T,
       sample_new_levels = "gaussian") %>%
  as_tibble() %>%
  bind_cols(nd)

d_pred

d_pred %>%
  ggplot(aes(x = log(population), y = Estimate)) +
  geom_ribbon(aes(ymin = Q1.5,  ymax = Q98.5), fill = "orange2", alpha = 1/3) +
  geom_ribbon(aes(ymin = Q5.5,  ymax = Q94.5), fill = "orange2", alpha = 1/3) +
  geom_ribbon(aes(ymin = Q16.5, ymax = Q83.5), fill = "orange2", alpha = 1/3) +
  geom_line(color = "orange4") +
  geom_text(data = d, aes(y = total_tools, label = culture),
            size = 2.33, color = "blue") +
  ggtitle("Total tools as a function of log(population)") +
  coord_cartesian(ylim = range(d$total_tools)) +
  theme_fivethirtyeight() +
  theme(plot.title = element_text(size = 12, hjust = .5))

# tidybayes::spread_draws ----

b12.4$formula

b12.4$data %>%
  glimpse()

# fit an intercepts only model
b12.7 <- brm(data = b12.4$data, family = binomial(),
             formula = pulled_left | trials(1) ~ 1 + (1 | actor),
             prior = c(prior(normal(0, 10), class = Intercept),
                       prior(cauchy(0, 1), class = sd)),
             iter = 5000, warmup = 1000, chains = 4, cores = 4,
             control = list(adapt_delta = 0.95),
             seed = 12)

print(b12.7)

b12.5$formula

b12.8 <-
  brm(data = b12.5$data, family = binomial,
      pulled_left | trials(1) ~ 1 + (1 | actor) + (1 | block),
      prior = c(prior(normal(0, 10), class = Intercept),
                prior(cauchy(0, 1), class = sd)),
      iter = 5000, warmup = 1000, chains = 4, cores = 4,
      control = list(adapt_delta = 0.95),
      seed = 12)

posterior_samples(b12.7) %>% str()

posterior_samples(b12.7) %>%
  transmute(`chimp 1's average probability of pulling left` = (b_Intercept + `r_actor[1,Intercept]`) %>% inv_logit_scaled()) %>%
  head()

p1 <-
  posterior_samples(b12.7) %>%
  transmute(`chimp 1's average probability of pulling left` = b_Intercept + `r_actor[1,Intercept]`,
            `chimp 2's average probability of pulling left` = b_Intercept + `r_actor[2,Intercept]`,
            `chimp 3's average probability of pulling left` = b_Intercept + `r_actor[3,Intercept]`,
            `chimp 4's average probability of pulling left` = b_Intercept + `r_actor[4,Intercept]`,
            `chimp 5's average probability of pulling left` = b_Intercept + `r_actor[5,Intercept]`,
            `chimp 6's average probability of pulling left` = b_Intercept + `r_actor[6,Intercept]`,
            `chimp 7's average probability of pulling left` = b_Intercept + `r_actor[7,Intercept]`) %>%
  mutate_all(inv_logit_scaled)

str(p1)

posterior_samples(b12.8) %>% glimpse()

p2 <-
  posterior_samples(b12.8) %>%
  transmute(`chimp 1's average probability of pulling left` = b_Intercept + `r_actor[1,Intercept]`,
            `chimp 2's average probability of pulling left` = b_Intercept + `r_actor[2,Intercept]`,
            `chimp 3's average probability of pulling left` = b_Intercept + `r_actor[3,Intercept]`,
            `chimp 4's average probability of pulling left` = b_Intercept + `r_actor[4,Intercept]`,
            `chimp 5's average probability of pulling left` = b_Intercept + `r_actor[5,Intercept]`,
            `chimp 6's average probability of pulling left` = b_Intercept + `r_actor[6,Intercept]`,
            `chimp 7's average probability of pulling left` = b_Intercept + `r_actor[7,Intercept]`) %>%
  mutate_all(inv_logit_scaled)

str(p2)


p3 <-
  posterior_samples(b12.8) %>%
  transmute(`chimp 1's average probability of pulling left` = b_Intercept + `r_actor[1,Intercept]` + `r_block[1,Intercept]`,
            `chimp 2's average probability of pulling left` = b_Intercept + `r_actor[2,Intercept]` + `r_block[1,Intercept]`,
            `chimp 3's average probability of pulling left` = b_Intercept + `r_actor[3,Intercept]` + `r_block[1,Intercept]`,
            `chimp 4's average probability of pulling left` = b_Intercept + `r_actor[4,Intercept]` + `r_block[1,Intercept]`,
            `chimp 5's average probability of pulling left` = b_Intercept + `r_actor[5,Intercept]` + `r_block[1,Intercept]`,
            `chimp 6's average probability of pulling left` = b_Intercept + `r_actor[6,Intercept]` + `r_block[1,Intercept]`,
            `chimp 7's average probability of pulling left` = b_Intercept + `r_actor[7,Intercept]` + `r_block[1,Intercept]`) %>%
  mutate_all(inv_logit_scaled)

str(p3)


coef(b12.7)

c1 <-
  coef(b12.7, summary = F)$actor[, , ] %>%
  as_tibble() %>%
  gather() %>%
  mutate(key   = str_c("chimp ", key, "'s average probability of pulling left"),
         value = inv_logit_scaled(value),
         # we need an iteration index for `spread()` to work properly
         iter  = rep(1:16000, times = 7)) %>%
  spread(key = key, value = value)

str(c1)

c2 <-
  coef(b12.8, summary = F)$actor[, , ] %>%
  as_tibble() %>%
  gather() %>%
  mutate(key   = str_c("chimp ", key, "'s average probability of pulling left"),
         value = inv_logit_scaled(value),
         iter  = rep(1:16000, times = 7)) %>%
  spread(key = key, value = value)

c2

coef(b12.8)

ranef(b12.8)

(nd <- b12.7$data %>% distinct(actor))

f1 <-
  fitted(b12.7,
         newdata = nd,
         summary = F,
         # within `fitted()`, this line does the same work that
         # `inv_logit_scaled()` did with the other two methods
         scale = "response") %>%
  as_tibble() %>%
  set_names(str_c("chimp ", 1:7, "'s average probability of pulling left"))

str(f1)

library(tidybayes)
b12.7 %>%
  spread_draws(b_Intercept, r_actor[actor,]) %>%
  filter(actor %in% c(1, 3))

b12.7 %>%
  spread_draws(b_Intercept, r_actor[actor,]) %>%
  ungroup() %>%
  select(.chain:.draw) %>%
  gather() %>%
  group_by(key) %>%
  summarise(min = min(value),
            max = max(value))

s1 <-
  b12.7 %>%
  spread_draws(b_Intercept, r_actor[actor,]) %>%
  mutate(p = (b_Intercept + r_actor) %>% inv_logit_scaled()) %>%
  select(.draw, actor, p) %>%
  ungroup() %>%
  mutate(actor = str_c("chimp ", actor, "'s average probability of pulling left")) %>%
  spread(value = p, key = actor)

s1
