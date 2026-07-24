## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  # fig.path = "img/",
  fig.align = "center",
  fig.dim = c(8, 6),
  out.width = "85%"
)

options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup--------------------------------------------------------------------
library(LaMa)

## ----data---------------------------------------------------------------------
beta = 2 # baseline standard deviation
phi = 0.95 # AR parameter of the log-volatility process
sigma = 0.5 # variability of the log-volatility process

n = 1000
set.seed(123)
g = rep(NA, n)
g[1] = rnorm(1, 0, sigma / sqrt(1-phi^2)) # stationary distribution of AR(1) process
for(t in 2:n){
  # sampling next state based on previous state and AR(1) equation
  g[t] = rnorm(1, phi*g[t-1], sigma)
}
# sampling zero-mean observations with standard deviation given by latent process
y = rnorm(n, 0, beta * exp(g/2)) 

# share returns
color = LaMaColors(2)
oldpar = par(mar = c(5,4,3,4.5)+0.1)
plot(y, type = "l", bty = "n", ylim = c(-40,20), yaxt = "n")
# true underlying standard deviation
lines(beta*exp(g)/7 - 40, col = color[2], lwd = 2)
axis(side=2, at = seq(-20,20,by=5), labels = seq(-20,20,by=5))
axis(side=4, at = seq(0,150,by=75)/7-40, labels = seq(0,150,by=75))
mtext("standard deviation", side=4, line=3, at = -30)
par(oldpar)

## ----mllk---------------------------------------------------------------------
nll = function(par) {
  getAll(par, dat)
  phi   = plogis(logit_phi); REPORT(phi)
  sigma = exp(log_sigma);    REPORT(sigma)
  beta  = exp(log_beta);     REPORT(beta)
  b     = seq(-bm, bm, length = m + 1) # intervals for midpoint quadrature
  h     = b[2] - b[1]                  # interval width
  bstar = (b[-1] + b[-(m+1)]) / 2     # interval midpoints
  # approximating t.p.m. via midpoint quadrature
  Gamma = sapply(bstar, dnorm, mean = phi * bstar, sd = sigma) * h
  # stationary distribution of the AR(1) process: g ~ N(0, sigma^2/(1-phi^2))
  delta = h * dnorm(bstar, 0, sigma / sqrt(1 - phi^2))
  # state-dependent densities: y_t | g_t ~ N(0, (beta * exp(g_t/2))^2)
  allprobs = t(sapply(y, dnorm, mean = 0, sd = beta * exp(bstar / 2)))
  -forward(delta, Gamma, allprobs)
}

## ----model, warning=FALSE-----------------------------------------------------
TapeConfig(matmul = "atomic")

bm = 5   # truncation of state space at ±bm
m  = 100 # number of quadrature points

par = list(
  logit_phi = qlogis(0.95),
  log_sigma = log(0.3),
  log_beta  = log(1)
)
dat = list(y = y, bm = bm, m = m)

obj = MakeADFun(nll, par, silent = TRUE)
system.time(
  opt <- nlminb(obj$par, obj$fn, obj$gr)
)
mod = report(obj)

## ----results------------------------------------------------------------------
mod$phi    # true: 0.95
mod$sigma  # true: 0.5
mod$beta   # true: 2

## decoding states
b     = seq(-bm, bm, length = m + 1)
h     = b[2] - b[1]
bstar = (b[-1] + b[-(m+1)]) / 2
Gamma    = sapply(bstar, dnorm, mean = mod$phi * bstar, sd = mod$sigma) * h
delta    = h * dnorm(bstar, 0, mod$sigma / sqrt(1 - mod$phi^2))
allprobs = t(sapply(y, dnorm, mean = 0, sd = mod$beta * exp(bstar / 2)))

probs  = stateprobs(delta, Gamma, allprobs) # local/soft decoding
states = viterbi(delta, Gamma, allprobs)    # global/hard decoding

oldpar = par(mar = c(5,4,3,4.5)+0.1)
plot(y, type = "l", bty = "n", ylim = c(-50,20), yaxt = "n")
# with many states, plotting the full conditional distribution is more informative
# than just the single most-probable state
maxprobs = apply(probs, 1, max)
for(t in 1:nrow(probs)){
  colend = round((probs[t,] / (maxprobs[t] * 5)) * 100)
  colend[which(colend < 10)] = paste0("0", colend[which(colend < 10)])
  points(rep(t, m), bstar * 4 - 35, col = paste0("#FFA200", colend), pch = 20)
}
# Viterbi-decoded path as a summary line
lines(bstar[states] * 4 - 35)

axis(side = 2, at = seq(-20, 20, by = 5), labels = seq(-20, 20, by = 5))
axis(side = 4, at = seq(-5, 5, by = 5) * 4 - 35, labels = seq(-5, 5, by = 5))
mtext("g", side = 4, line = 3, at = -30)
par(oldpar)

## ----jnll---------------------------------------------------------------------
jnll = function(par) {
  getAll(par, dat)
  phi   = plogis(logit_phi); REPORT(phi)
  sigma = exp(log_sigma);    REPORT(sigma)
  beta  = exp(log_beta);     REPORT(beta)
  jnll = 0
  # stationary initial distribution: g_1 ~ N(0, sigma^2 / (1 - phi^2))
  jnll = jnll - dnorm(g[1], 0, sigma / sqrt(1 - phi^2), log = TRUE)
  # AR(1) state transitions: g_t | g_{t-1} ~ N(phi * g_{t-1}, sigma^2)
  jnll = jnll - sum(dnorm(g[-1], phi * g[-length(g)], sigma, log = TRUE))
  # observations: y_t | g_t ~ N(0, (beta * exp(g_t / 2))^2)
  jnll = jnll - sum(dnorm(y, 0, beta * exp(g / 2), log = TRUE))
  jnll
}

## ----lap_model, warning=FALSE-------------------------------------------------
dat = list(y = y)

par_lap = list(
  logit_phi = qlogis(0.95),
  log_sigma  = log(0.3),
  log_beta   = log(1),
  g          = rep(0, n)   # random effects initialised at zero
)

obj_lap = MakeADFun(jnll, par_lap, random = "g", silent = TRUE)
system.time(
  opt_lap <- nlminb(obj_lap$par, obj_lap$fn, obj_lap$gr)
)
mod_lap = report(obj_lap)

## ----lap_results--------------------------------------------------------------
mod_lap$phi    # true: 0.95
mod_lap$sigma  # true: 0.5
mod_lap$beta   # true: 2

