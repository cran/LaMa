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

## ----parameters---------------------------------------------------------------
# generator matrix Q
Q = matrix(c(-0.5, 0.5,
              1.0, -1.0), nrow = 2, byrow = TRUE)

# state-dependent (normal) distribution parameters
mu = c(5, 20)
sigma = c(2, 5)

## ----data---------------------------------------------------------------------
set.seed(123)

k = 200              # number of state switches to simulate
s = rep(NA, k)       # latent state sequence
trans_times = rep(NA, k) # cumulative times at which state switches occur

# initialise: draw first state and its dwell time
s[1] = sample(1:2, 1)
dwell = rexp(1, -Q[s[1], s[1]])  # dwell time ~ Exp(-q_ii)
trans_times[1] = dwell

# Poisson arrivals during first sojourn [0, dwell], uniform within
new_obs = sort(runif(rpois(1, dwell), 0, dwell))
obs_times = new_obs
x = rnorm(length(new_obs), mu[s[1]], sigma[s[1]])

for(t in 2:k){
  # for 2 states, the only possible next state is 3 - current
  s[t] = 3 - s[t-1]

  # dwell time in new state, and absolute transition time
  dwell = rexp(1, -Q[s[t], s[t]])
  trans_times[t] = trans_times[t-1] + dwell

  # Poisson arrivals during sojourn [t_start, t_end], uniform within
  t_start = trans_times[t-1]
  t_end   = trans_times[t]
  new_obs = sort(runif(rpois(1, dwell), t_start, t_end))

  obs_times = c(obs_times, new_obs)
  x = c(x, rnorm(length(new_obs), mu[s[t]], sigma[s[t]]))
}

## ----vis_ctHMM----------------------------------------------------------------
color = LaMaColors(2)

n = length(obs_times)
plot(obs_times[1:50], x[1:50], pch = 16, bty = "n",
     xlab = "observation times", ylab = "x", ylim = c(-5, 25))
segments(x0 = c(0, trans_times[1:48]), x1 = trans_times[1:49],
         y0 = rep(-5, 49), y1 = rep(-5, 49),
         col = color[s[1:49]], lwd = 4)
legend("topright", lwd = 2, col = color,
       legend = c("state 1", "state 2"), box.lwd = 0)

## ----nll----------------------------------------------------------------------
nll = function(par) {
  getAll(par, dat)
  mu = exp(log_mu); REPORT(mu)
  sigma = exp(log_sigma); REPORT(sigma)
  Q = generator(log_qs); REPORT(Q)  # maps unconstrained params -> valid Q
  Pi = stationary_ct(Q)             # stationary distribution of the CTMC
  Qube = tpm_ct(Q, timediff)     # exp(Q * dt) for each time difference
  allprobs = matrix(1, length(x), N)
  ind = which(!is.na(x))
  for(j in 1:N) allprobs[ind, j] = dnorm(x[ind], mu[j], sigma[j])
  -forward(Pi, Qube, allprobs)
}

## ----model, warning = FALSE---------------------------------------------------
N = 2
timediff = diff(obs_times) # time differences between consecutive observations

par = list(
  log_mu = log(c(5, 15)),    # initial state-dependent means (log scale)
  log_sigma = log(c(3, 5)),  # initial state-dependent sds (log scale)
  log_qs = log(c(1, 0.5))    # initial off-diagonal generator entries (log scale)
)

dat = list(
  x = x,
  timediff = timediff,
  N = N
)

obj = MakeADFun(nll, par, silent = TRUE)
system.time(
  opt <- nlminb(obj$par, obj$fn, obj$gr)
)
mod = report(obj)

## ----results------------------------------------------------------------------
mod$mu     # true: 5, 20
mod$sigma  # true: 2, 5
round(mod$Q, 3) # true: Q as defined above

## ----dwell_times--------------------------------------------------------------
round(1 / diag(-mod$Q), 2) # estimated mean dwell times; true: 2, 1

## ----tpm_curves---------------------------------------------------------------
dt_seq = seq(0, 10, length.out = 300)
Qube_seq = tpm_ct(mod$Q, dt_seq)
Pi_ct = stationary_ct(mod$Q)

plot(dt_seq, Qube_seq[1, 2, ], type = "l", lwd = 2, col = color[1], bty = "n",
     xlab = expression(Delta*t), ylab = "transition probability", ylim = c(0, 1))
lines(dt_seq, Qube_seq[2, 1, ], lwd = 2, col = color[2])
abline(h = Pi_ct[2], lty = 2, col = color[1])  # long-run limit of gamma_12
abline(h = Pi_ct[1], lty = 2, col = color[2])  # long-run limit of gamma_21
legend("right", lwd = 2, col = color, bty = "n",
       legend = c(expression(gamma[12](Delta*t)), expression(gamma[21](Delta*t))))

## ----state_decoding-----------------------------------------------------------
states = viterbi(mod = mod)

plot(obs_times[1:50], x[1:50], pch = 16, bty = "n",
     col = color[states[1:50]],
     xlab = "observation times", ylab = "x", ylim = c(-5, 25))
legend("topright", pch = 16, col = color,
       legend = paste("state", 1:N), box.lwd = 0)

## ----parameters2--------------------------------------------------------------
# generator matrix Q
Q = matrix(c(-0.5, 0.2, 0.3,
              1.0, -2.0, 1.0,
              0.4,  0.6, -1.0), nrow = 3, byrow = TRUE)

# state-dependent (normal) distribution parameters
mu = c(5, 15, 30)
sigma = c(2, 3, 5)

## ----data2--------------------------------------------------------------------
set.seed(123)

k = 200
s = rep(NA, k)
trans_times = rep(NA, k)

# initialise
s[1] = sample(1:3, 1)
dwell = rexp(1, -Q[s[1], s[1]])
trans_times[1] = dwell

new_obs = sort(runif(rpois(1, dwell), 0, dwell))
obs_times = new_obs
x = rnorm(length(new_obs), mu[s[1]], sigma[s[1]])

for(t in 2:k){
  # conditional transition probabilities: omega_ij = q_ij / (-q_ii)
  omega = Q[s[t-1], -s[t-1]] / -Q[s[t-1], s[t-1]]
  s[t] = sample((1:3)[-s[t-1]], 1, prob = omega)

  dwell = rexp(1, -Q[s[t], s[t]])
  trans_times[t] = trans_times[t-1] + dwell

  t_start = trans_times[t-1]
  t_end   = trans_times[t]
  new_obs = sort(runif(rpois(1, dwell), t_start, t_end))

  obs_times = c(obs_times, new_obs)
  x = c(x, rnorm(length(new_obs), mu[s[t]], sigma[s[t]]))
}

## ----model2, warning = FALSE--------------------------------------------------
N = 3
timediff = diff(obs_times)

par = list(
  log_mu = log(c(5, 10, 25)),   # initial state-dependent means
  log_sigma = log(c(2, 2, 6)),  # initial state-dependent sds
  log_qs = rep(0, N*(N-1))      # initial off-diagonal generator entries (exp(0) = 1)
)

dat = list(
  x = x,
  timediff = timediff,
  N = N
)

obj2 = MakeADFun(nll, par, silent = TRUE)
system.time(
  opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)
)
mod2 = report(obj2)

## ----results2-----------------------------------------------------------------
mod2$mu     # true: 5, 15, 30
mod2$sigma  # true: 2, 3, 5
round(mod2$Q, 3) # true: Q as defined above

## ----dwell_times2-------------------------------------------------------------
round(1 / diag(-mod2$Q), 2) # estimated mean dwell times; true: 2, 0.5, 1

## ----state_decoding2----------------------------------------------------------
color3 = LaMaColors(3)
states2 = viterbi(mod = mod2)

plot(obs_times[1:50], x[1:50], pch = 16, bty = "n",
     col = color3[states2[1:50]],
     xlab = "observation times", ylab = "x", ylim = c(-5, 40))
legend("topright", pch = 16, col = color3,
       legend = paste("state", 1:N), box.lwd = 0)

