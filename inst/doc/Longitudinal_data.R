## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  # fig.path = "img/",
  fig.align = "center",
  fig.dim = c(8, 6),
  out.width = "75%"
)

options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup--------------------------------------------------------------------
library(LaMa)

## ----data_pooling-------------------------------------------------------------
# parameters shared across individuals
mu = c(15, 60)    # state-dependent means
sigma = c(10, 40) # state-dependent standard deviations
Gamma = matrix(c(0.95, 0.05,
                 0.15, 0.85), nrow = 2, byrow = TRUE)
delta = stationary(Gamma)

set.seed(123)
K = 200 # number of individuals
N = 2   # number of states
n = 50  # observations per individual

s = x = rep(NA, n * K)
for(k in 1:K){
  sk = xk = rep(NA, n)
  sk[1] = sample(1:2, 1, prob = delta)
  xk[1] = rnorm(1, mu[sk[1]], sigma[sk[1]])
  for(t in 2:n){
    sk[t] = sample(1:2, 1, prob = Gamma[sk[t-1],])
    xk[t] = rnorm(1, mu[sk[t]], sigma[sk[t]])
  }
  s[(k-1)*n + 1:n] = sk
  x[(k-1)*n + 1:n] = xk
}

trackID = rep(1:K, each = n)

## ----nll_pool-----------------------------------------------------------------
nll_pool = function(par) {
  getAll(par, dat)
  Gamma = tpm(eta)
  delta = stationary(Gamma)
  mu = exp(log_mu); REPORT(mu)
  sigma = exp(log_sigma); REPORT(sigma)
  allprobs = matrix(1, length(x), N)
  for(j in 1:N) allprobs[,j] = dnorm(x, mu[j], sigma[j])
  -forward(delta, Gamma, allprobs, trackID)
}

## ----model_pool, warning = FALSE, cache = TRUE--------------------------------
par = list(
  eta = rep(-2, N*(N-1)),    # initial t.p.m. parameters (logit scale)
  log_mu = log(c(15, 60)),   # initial state-dependent means (log scale)
  log_sigma = log(c(10, 40)) # initial state-dependent sds (log scale)
)

dat = list(
  x = x,
  trackID = trackID,
  N = N
)

obj_pool = MakeADFun(nll_pool, par, silent = TRUE)
system.time(
  opt_pool <- nlminb(obj_pool$par, obj_pool$fn, obj_pool$gr)
)
mod_pool = report(obj_pool)

## ----pool_results-------------------------------------------------------------
mod_pool$mu     # true: 15, 60
mod_pool$sigma  # true: 10, 40
mod_pool$Gamma  # true: c(0.95, 0.05, 0.15, 0.85)

## ----data_partial-------------------------------------------------------------
K = 5   # number of individuals
N = 2   # number of states
n = 200 # observations per individual

# state-dependent parameters shared across individuals
mu = c(15, 60)
sigma = c(10, 40)

# individual-specific tpms depending on a covariate (e.g. age)
set.seed(123)
z = rnorm(K)
beta = matrix(c(-2, -2, 1, -1), nrow = N*(N-1))

# compute K individual-specific tpms
Gamma = tpm(beta, matrix(z, ncol = 1))
# stationary() on an array returns a c(K, N) matrix of stationary distributions
Delta = stationary(Gamma)

# simulate all tracks
set.seed(123)
s = x = rep(NA, n * K)
for(k in 1:K){
  sk = xk = rep(NA, n)
  sk[1] = sample(1:2, 1, prob = Delta[k,])
  xk[1] = rnorm(1, mu[sk[1]], sigma[sk[1]])
  for(t in 2:n){
    sk[t] = sample(1:2, 1, prob = Gamma[sk[t-1],,k])
    xk[t] = rnorm(1, mu[sk[t]], sigma[sk[t]])
  }
  s[(k-1)*n + 1:n] = sk
  x[(k-1)*n + 1:n] = xk
}

## ----nll_partial--------------------------------------------------------------
nll_partial = function(par) {
  getAll(par, dat)
  Gamma = tpm(beta, matrix(z, ncol = 1)) # array of K individual-specific tpms
  Delta = stationary(Gamma) # when passed an array, returns a matrix of dim c(K, N),
                            # each row being the stationary distribution for one tpm
  mu = exp(log_mu); REPORT(mu)
  sigma = exp(log_sigma); REPORT(sigma)
  REPORT(beta)
  allprobs = matrix(1, length(x), N)
  for(j in 1:N) allprobs[,j] = dnorm(x, mu[j], sigma[j])
  # Delta matrix and Gamma array supply individual-specific starts and tpms
  -forward(Delta, Gamma, allprobs, trackID)
}

## ----model_partial, warning = FALSE-------------------------------------------
trackID = rep(1:K, each = n)

par = list(
  beta = matrix(c(-2, -2, 0, 0), nrow = N*(N-1)), # regression coefficients for tpm
  log_mu = log(c(15, 60)),                         # initial state-dependent means
  log_sigma = log(c(10, 40))                       # initial state-dependent sds
)

dat = list(
  x = x,
  z = z,
  trackID = trackID,
  N = N
)

obj_partial = MakeADFun(nll_partial, par, silent = TRUE)
system.time(
  opt_partial <- nlminb(obj_partial$par, obj_partial$fn, obj_partial$gr)
)
mod_partial = report(obj_partial)

## ----partial_results----------------------------------------------------------
mod_partial$beta   # true: c(-2, -2, 1, -1)
mod_partial$mu     # true: 15, 60
mod_partial$sigma  # true: 10, 40

## ----data_re------------------------------------------------------------------
set.seed(123)
K = 50    # number of individuals
N = 2     # number of states
n = 100   # observations per individual

# Shared state-dependent parameters
mu = c(15, 60)
sigma = c(10, 40)

# Population-level intercepts and individual random effects
beta_0 = c(-2, -2)  # one per off-diagonal element
sigma_u = c(1, 1)   # one per off-diagonal element

u_true = matrix(rnorm(K * N*(N-1), 0, rep(sigma_u, each = K)),
                nrow = K, ncol = N*(N-1))

# K individual-specific tpms and their stationary distributions
Gamma_true = array(NA, c(N, N, K))
for(k in 1:K) Gamma_true[,,k] = tpm(beta_0 + u_true[k,])
Delta_true = stationary(Gamma_true)

# Simulate all tracks
s = x = rep(NA, n * K)
for(k in 1:K){
  sk = xk = rep(NA, n)
  sk[1] = sample(1:2, 1, prob = Delta_true[k,])
  xk[1] = rnorm(1, mu[sk[1]], sigma[sk[1]])
  for(t in 2:n){
    sk[t] = sample(1:2, 1, prob = Gamma_true[sk[t-1],,k])
    xk[t] = rnorm(1, mu[sk[t]], sigma[sk[t]])
  }
  s[(k-1)*n + 1:n] = sk
  x[(k-1)*n + 1:n] = xk
}

trackID = rep(1:K, each = n)

## ----jnll_re------------------------------------------------------------------
jnll_re = function(par) {
  getAll(par, dat)

  # Random effect standard deviations
  sigma_u = exp(log_sigma_u); REPORT(sigma_u); ADREPORT(sigma_u)

  # Gaussian prior on u: u[,ij] ~ N(0, sigma_u[ij])
  jnll = 0
  for(ij in 1:(N*(N-1))) {
    jnll = jnll - sum(dnorm(u[,ij], 0, sigma_u[ij], log = TRUE))
  }

  # Individual-specific tpms from population intercept + random effect
  Gamma = AD(array(0, c(N, N, K)))
  for(k in 1:K) Gamma[,,k] = tpm(beta_0 + u[k,])
  Delta = stationary(Gamma)

  # Observation model
  mu = exp(log_mu); REPORT(mu)
  sigma = exp(log_sigma); REPORT(sigma)
  allprobs = matrix(1, length(x), N)
  for(j in 1:N) allprobs[,j] = dnorm(x, mu[j], sigma[j])

  jnll - forward(Delta, Gamma, allprobs, trackID)
}

## ----model_re, warning = FALSE, cache = TRUE----------------------------------
par = list(
  beta_0 = rep(-2, N*(N-1)),      # population-level intercepts
  log_sigma_u = log(c(1, 1)),     # log-SDs of random effects
  u = matrix(0, K, N*(N-1)),      # random effects (integrated out)
  log_mu = log(c(15, 60)),        # initial state-dependent means
  log_sigma = log(c(10, 40))      # initial state-dependent sds
)

dat = list(
  x = x,
  trackID = trackID,
  K = K,
  N = N
)

obj_re = MakeADFun(jnll_re, par, random = "u", silent = TRUE)
system.time(
  opt_re <- nlminb(obj_re$par, obj_re$fn, obj_re$gr)
)
mod_re = report(obj_re)

## ----re_results---------------------------------------------------------------
mod_re$sigma_u  # true: 1, 1
mod_re$mu       # true: 15, 60
mod_re$sigma    # true: 10, 40

## ----re_sdr-------------------------------------------------------------------
sdr_re = sdreport(obj_re)
summary(sdr_re, "report") # sigma_u with SEs

