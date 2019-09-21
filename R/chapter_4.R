library(tidyverse)
library(magrittr)
library(rethinking)

# Normal by addition ----
pos <- replicate(1000, sum(runif(16, -1, 1)^5))

hist(pos)

# Normal by multiplication ----

# small multiplicative increases are approximately additive
prod(1 + runif(12, 0, 0.1))

growth = replicate(10000, prod(1 + runif(12, 0, 0.1)))
dens(growth, norm.comp = TRUE)

big = replicate(1e4, prod(1 + runif(12, 0, 0.5)))
small = replicate(1e4, prod(1 + runif(12, 0, 0.01)))

dens(big, norm.comp = T) # big growth doesn't look normal
dens(small, norm.comp = T) # small growth looks normal

# but log of multiplicative go=rowth looks normal!
log.big = replicate(1e4, log(prod(1 + runif(12, 0, 0.5))))
dens(log.big, norm.comp = T)

# Gaussian model of height ----
data("Howell1")
d = Howell1

str(d)

d2 = d %>% filter(age >= 18)

nrow(d2)

dens(d2$height)

d2 %>%
  ggplot(aes(height)) +
  geom_histogram(aes(y = ..density..)) +
  geom_density(size = 1.2)

# Specify likelihood for height, then priors. Plot priors!

curve(dnorm(x, 178, 20), from = 100, to = 250)
curve(dunif(x, 0, 50), from = 10, to = 60)

# Priors imply a distribution (of distributions) for prior heights. Sample and plot it

sample_mu = rnorm(1e4, 178, 20)
sample_sigma = runif(1e4, 0, 10)
prior_h = rnorm(1e4, sample_mu, sample_sigma)
dens(prior_h, norm.comp = T)

# Posterior through brute force ----

mu.list = seq(from = 140, to = 160, length.out = 200)
sigma.list = seq(from = 4, to = 9, length.out = 200)
post = expand.grid(mu = mu.list, sigma = sigma.list)
post$LL = sapply(1:nrow(post), function(i) sum(dnorm(
  d2$height,
  mean = post$mu[i],
  sd = post$sigma[i],
  log = TRUE
)))

post$prod = post$LL + dnorm(post$mu, 178, 20, TRUE) +
  dunif(post$sigma, 0, 50, TRUE)
post$prob = exp(post$prod - max(post$prod))

contour_xyz(post$mu, post$sigma, post$prob)

image_xyz(post$mu, post$sigma, post$prob)

# Sampling from posterior with two parameters
sample.rows = sample(1:nrow(post), size = 1e4, replace = TRUE)
sample.mu = post$mu[sample.rows]
sample.sigma = post$sigma[sample.rows]

plot(sample.mu, sample.sigma, cex = 0.8, pch = 16,
     col = col.alpha(rangi2, 0.5))

# Sampling from posterior with two parameters: tidyverse version

n <- 200

d_grid <-
  tibble(mu    = seq(from = 140, to = 160, length.out = n),
         sigma = seq(from = 4,   to = 9,   length.out = n)) %>%
  # we'll accomplish with `tidyr::expand()` what McElreath did with base R `expand.grid()`
  expand(mu, sigma)

head(d_grid)

grid_function <- function(mu, sigma){
  dnorm(d2$height, mean = mu, sd = sigma, log = T) %>%
    sum()
}

# For each combination of (mu, sigma), draw a value from the
# prior density
d_grid <-
  d_grid %>%
  mutate(log_likelihood = map2(mu, sigma, grid_function)) %>%
  unnest() %>%
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

dens(d_grid$prior_mu)
HPDI(sample.mu)
HPDI(sample.sigma)

# Quadratic approximation ----
data(Howell1)
d = Howell1
d2 = d %>% filter(age >= 18)

flist = alist(
  height ~ dnorm(mu, sigma),
  mu ~ dnorm(178, 20),
  sigma ~ dunif(0, 50)
)

m4.1 = map(flist, data = d2)

precis(m4.1)

# Estimate with a more informative prior (smaller sd for mu's prior)
m4.2 = map(
  flist = alist(
    height ~ dnorm(mu, sigma),
    mu ~ dnorm(178, 0.1),
    sigma ~ dunif(0, 50)
  ),
  data = d2
)

precis(m4.1)
precis(m4.2)

vcov(m4.1)
diag(vcov(m4.1)) %>% sqrt()
cov2cor(vcov(m4.1))

# Sampling from a multi-dimensional posterior ----
post = extract.samples(m4.1, n = 1e4)
head(post)
precis(post)
plot(post)

# Adding a predictor ----

ggplot(d2, aes(weight, height)) + geom_point()

m4.3 = map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*weight,
    a ~ dnorm(156, 100),
    b ~ dnorm(0, 10),
    sigma ~ dunif(0, 50)
  ),
  data = d2
)

# Centering ----
# Strong correlations between parameters can make it difficult to fit
# model to data. Centering can help.

precis(m4.3)
precis(m4.3, corr = TRUE) # see correlation between a and b

d2 %<>% mutate(weight.c = weight - mean(weight))

# fit again with centered weight
m4.4 = map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*weight.c,
    a ~ dnorm(156, 100),
    b ~ dnorm(0, 10),
    sigma ~ dunif(0, 50)
  ),
  data = d2
)

precis(m4.4, corr = TRUE)

ggplot(d2, aes(weight.c, height)) + geom_point()

plot(height ~ weight, data = d2)
abline(a = coef(m4.3)["a"], b = coef(m4.3)["b"])

post = extract.samples(m4.3)

head(post)

N = 1e6
dN = d2 %>% filter(row_number() <= N)

mN = map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*weight,
    a ~ dnorm(178, 100),
    b ~ dnorm(0, 10),
    sigma ~ dunif(0, 50)
  ),
  data = dN
)

post = extract.samples(mN, n = 20)

plot(dN$weight, dN$height,
     xlim = range(d2$weight),
     ylim = range(d2$height),
     col = rangi2, xlab = "weight", ylab = "height")
mtext(concat("N = ", N))

for (i in 1:20)
  abline(a = post$a[i], b = post$b[i], col = col.alpha("black", 0.3))


mu_at_50 = post$a + post$b*50
mu_at_50
dens(mu_at_50, col = rangi2, lwd = 2, xlab = "mu|weight = 50")

HPDI(mu_at_50, prob = 0.89)

mu = link(m4.3)
str(mu)

weight.seq = seq(from = 25, to = 70, by = 1)
mu = link(m4.3, data = data.frame(weight = weight.seq))

str(mu)

plot(height ~ weight, d2, type = "n")

for (i in 1:100){
  points(weight.seq, mu[i,], pch = 16, col = col.alpha(rangi2, 0.1))
}

mu.mean = apply(mu, 2, mean)
mu.HPDI = apply(mu, 2, HPDI, prob = 0.89)

plot(height ~ weight, data = d2, col = col.alpha(rangi2, 0.5))
lines(weight.seq, mu.mean)
shade(mu.HPDI, weight.seq)

# Fitting polynomials ----
# not recommended, but illustrative

# weight-height relationship with full population looks non-linear
ggplot(d, aes(weight, height)) + geom_point()

# standardize weight
d %<>% mutate(weight.s = (weight - mean(weight))/sd(weight))

