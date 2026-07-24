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
lambda = c(2, 15)  # state-dependent Poisson rates
Q = matrix(c(-0.5, 0.5,
              2.0, -2.0), nrow = 2, byrow = TRUE)

## ----simulation---------------------------------------------------------------
set.seed(123)

k = 200              # number of state switches to simulate
s = rep(NA, k)       # latent state sequence
trans_times = rep(NA, k)  # cumulative transition times

# initialise
s[1] = sample(1:2, 1)
dwell = rexp(1, -Q[s[1], s[1]])
trans_times[1] = dwell

n_arrivals = rpois(1, lambda[s[1]] * dwell)
arrival_times = sort(runif(n_arrivals, 0, dwell))

for (t in 2:k) {
  # 2-state chain: always switches to the other state
  s[t] = 3 - s[t-1]

  dwell = rexp(1, -Q[s[t], s[t]])
  trans_times[t] = trans_times[t-1] + dwell

  t_start = trans_times[t-1]
  t_end   = trans_times[t]

  # Poisson arrivals uniformly distributed within sojourn
  n_arrivals = rpois(1, lambda[s[t]] * dwell)
  arrival_times = c(arrival_times, sort(runif(n_arrivals, t_start, t_end)))
}

## ----vis_MMPP-----------------------------------------------------------------
color = LaMaColors(2)
n = length(arrival_times)
plot(c(0, arrival_times), 0:n, type = "l", bty = "n",
     xlab = "time", ylab = "cumulative arrivals",
     ylim = c(-n * 0.05, n))
segments(x0 = c(0, trans_times[-k]), x1 = trans_times,
         y0 = rep(-n * 0.05, k), y1 = rep(-n * 0.05, k),
         col = color[s], lwd = 4)
legend("topleft", lwd = 2, col = color, legend = paste("state", 1:2), box.lwd = 0)

## ----nll----------------------------------------------------------------------
nll = function(par) {
  getAll(par, dat)
  lambda = exp(log_lambda); REPORT(lambda)
  Q = generator(log_qs); REPORT(Q)
  Pi = stationary_ct(Q); REPORT(Pi)
  Qube = tpm_ct(Q, timediff, lambda)
  allprobs = matrix(rep(lambda, length(timediff) + 1),
                    nrow = length(timediff) + 1, ncol = N, byrow = TRUE)
  allprobs[1,] = 1
  -forward(Pi, Qube, allprobs)
}

## ----model, warning=FALSE-----------------------------------------------------
N = 2
timediff = diff(arrival_times)

par = list(
  log_lambda = log(c(2, 15)),
  log_qs     = log(c(2, 0.5))
)
dat = list(
  timediff = timediff,
  N = N
)

obj = MakeADFun(nll, par, silent = TRUE)
system.time(
  opt <- nlminb(obj$par, obj$fn, obj$gr)
)
mod = report(obj)

## ----results------------------------------------------------------------------
mod$lambda           # true: 2, 15
round(mod$Q, 3)      # true: Q as defined above
round(1 / diag(-mod$Q), 2)  # mean dwell times; true: 2, 0.5
mod$Pi               # stationary distribution of the CTMC

## ----state_decoding-----------------------------------------------------------
states = viterbi(mod = mod)

plot(arrival_times[1:100], rep(0.5, 100), type = "h", bty = "n",
     ylim = c(0, 1), yaxt = "n", col = color[states[1:100]],
     xlab = "arrival times", ylab = "")
legend("top", lwd = 2, col = color, legend = paste("state", 1:N), box.lwd = 0)

## ----parameters2--------------------------------------------------------------
lambda = c(1, 5, 20)
Q = matrix(c(-0.5, 0.3, 0.2,
              0.7, -1.0, 0.3,
              1.0,  1.0, -2.0), nrow = 3, byrow = TRUE)
mu    = c(-5, 0, 5)
sigma = c(2, 1, 2)

color = LaMaColors(3)
curve(dnorm(x, mu[2], sigma[2]), xlim = c(-10, 10), bty = "n",
      lwd = 2, col = color[2], n = 200, ylab = "density", xlab = "mark")
curve(dnorm(x, mu[1], sigma[1]), add = TRUE, lwd = 2, col = color[1], n = 200)
curve(dnorm(x, mu[3], sigma[3]), add = TRUE, lwd = 2, col = color[3], n = 200)
legend("topright", lwd = 2, col = color, bty = "n", legend = paste("state", 1:3))

## ----simulation2--------------------------------------------------------------
set.seed(123)

k = 200
s = rep(NA, k)
trans_times = rep(NA, k)

s[1] = sample(1:3, 1)
dwell = rexp(1, -Q[s[1], s[1]])
trans_times[1] = dwell

n_arrivals = rpois(1, lambda[s[1]] * dwell)
new_arrivals = sort(runif(n_arrivals, 0, dwell))
arrival_times = new_arrivals
marks = rnorm(n_arrivals, mu[s[1]], sigma[s[1]])

for (t in 2:k) {
  # conditional next-state probabilities: omega_ij = q_ij / (-q_ii)
  omega = Q[s[t-1], -s[t-1]] / -Q[s[t-1], s[t-1]]
  s[t] = sample((1:3)[-s[t-1]], 1, prob = omega)

  dwell = rexp(1, -Q[s[t], s[t]])
  trans_times[t] = trans_times[t-1] + dwell

  t_start = trans_times[t-1]
  t_end   = trans_times[t]

  n_arrivals = rpois(1, lambda[s[t]] * dwell)
  new_arrivals = sort(runif(n_arrivals, t_start, t_end))
  arrival_times = c(arrival_times, new_arrivals)
  marks = c(marks, rnorm(n_arrivals, mu[s[t]], sigma[s[t]]))
}

## ----vis_MMMPP----------------------------------------------------------------
n = length(arrival_times)
plot(arrival_times[1:100], marks[1:100], pch = 16, bty = "n",
     ylim = c(-9, 9), xlab = "arrival times", ylab = "marks")
segments(x0 = c(0, trans_times[1:98]), x1 = trans_times[1:99],
         y0 = rep(-9, 99), y1 = rep(-9, 99), col = color[s[1:99]], lwd = 4)
legend("topright", lwd = 2, col = color,
       legend = paste("state", 1:3), box.lwd = 0)

## ----nll2---------------------------------------------------------------------
nll2 = function(par) {
  getAll(par, dat)
  lambda = exp(log_lambda); REPORT(lambda)
  REPORT(mu)
  sigma = exp(log_sigma); REPORT(sigma)
  Q = generator(log_qs); REPORT(Q)
  Pi = stationary_ct(Q); REPORT(Pi)
  Qube = tpm_ct(Q, timediff, lambda)
  allprobs = matrix(1, length(marks), N)
  for (j in 1:N) allprobs[, j] = dnorm(marks, mu[j], sigma[j])
  # multiply rows 2:n by lambda to fold in the arrival-rate contribution
  allprobs[-1,] = allprobs[-1,] * matrix(lambda, length(marks) - 1, N, byrow = TRUE)
  -forward(Pi, Qube, allprobs)
}

## ----model2, warning=FALSE----------------------------------------------------
N = 3
timediff = diff(arrival_times)

par = list(
  log_lambda = log(c(1, 5, 20)),
  mu         = c(-5, 0, 5),
  log_sigma  = log(c(2, 1, 2)),
  log_qs     = log(c(0.7, 1, 0.3, 1, 0.2, 0.3))
)
dat = list(
  marks    = marks,
  timediff = timediff,
  N = N
)

obj2 = MakeADFun(nll2, par, silent = TRUE)
system.time(
  opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)
)
mod2 = report(obj2)

## ----results2-----------------------------------------------------------------
mod2$lambda                  # true: 1, 5, 20
mod2$mu                      # true: -5, 0, 5
mod2$sigma                   # true: 2, 1, 2
round(mod2$Q, 3)             # true: Q as defined above
round(1 / diag(-mod2$Q), 2) # mean dwell times; true: 2, 1, 0.5
mod2$Pi                      # stationary distribution

## ----mark_density-------------------------------------------------------------
# marks arrive at rate lambda_j in state j, so the fraction of marks from
# state j is proportional to Pi[j] * lambda[j], not Pi[j] alone
w = mod2$Pi * mod2$lambda / sum(mod2$Pi * mod2$lambda)

hist(marks, prob = TRUE, breaks = 50, border = "white", main = "",
     xlab = "mark")
for (j in 1:N) {
  curve(w[j] * dnorm(x, mod2$mu[j], mod2$sigma[j]),
        add = TRUE, lwd = 2, col = color[j], n = 300)
}
curve(rowSums(sapply(1:N, function(j) w[j] * dnorm(x, mod2$mu[j], mod2$sigma[j]))),
      add = TRUE, lwd = 2, lty = 2, n = 300)
legend("topright", lwd = 2, lty = c(1, 1, 1, 2), bty = "n",
       col = c(color, "black"), legend = c(paste("state", 1:N), "marginal"))

## ----state_decoding2----------------------------------------------------------
states2 = viterbi(mod = mod2)

plot(arrival_times[1:100], marks[1:100], pch = 16, bty = "n",
     col = color[states2[1:100]],
     ylim = c(-9, 9), xlab = "arrival times", ylab = "marks")
legend("topright", pch = 16, col = color,
       legend = paste("state", 1:3), box.lwd = 0)

