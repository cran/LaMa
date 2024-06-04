## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  # fig.path = "img/",
  fig.align = "center",
  fig.dim = c(8, 6),
  out.width = "85%"
)

## ----setup--------------------------------------------------------------------
# loading the package
library("LaMa")

## ----parameters---------------------------------------------------------------
# 2-state example

# generator matrix Q:
Q = matrix(c(-0.5,0.5,1,-1), nrow = 2, byrow = TRUE)
# state 1 has a smaller rate (dwell-time in state one ~ Exp(1)), i.e. it exhibits 
# longer dwell times than state 3 with rate 3.

# parameters for the state-dependent (normal) distributions
mu = c(5, 20)
sigma = c(2, 5)

## ----data---------------------------------------------------------------------
set.seed(123)

k = 200 # number of state switches
trans_times = s = rep(NA, k) # time points where the chain transitions
s[1] = sample(1:2, 1) # initial distribuion c(0.5, 0.5)
# exponentially distributed waiting times
trans_times[1] = rexp(1, -Q[s[1],s[1]])
n_arrivals = rpois(1, trans_times[1])
obs_times = sort(runif(n_arrivals, 0, trans_times[1]))
x = rnorm(n_arrivals, mu[s[1]], sigma[s[1]])
for(t in 2:k){
  s[t] = c(1,2)[-s[t-1]] # for 2-states, always a state swith when transitioning
  # exponentially distributed waiting times
  trans_times[t] = trans_times[t-1] + rexp(1, -Q[s[t], s[t]])
  n_arrivals = rpois(1, trans_times[t]-trans_times[t-1])
  obs_times = c(obs_times, 
                sort(runif(n_arrivals, trans_times[t-1], trans_times[t])))
  x = c(x, rnorm(n_arrivals, mu[s[t]], sigma[s[t]]))
}

## ----vis_ctHMM----------------------------------------------------------------
color = c("orange", "deepskyblue")

n = length(obs_times)
plot(obs_times[1:50], x[1:50], pch = 16, bty = "n", xlab = "observation times", 
     ylab = "x", ylim = c(-5,25))
segments(x0 = c(0,trans_times[1:48]), x1 = trans_times[1:49], 
         y0 = rep(-5,50), y1 = rep(-5,50), col = color[s[1:49]], lwd = 4)
legend("topright", lwd = 2, col = color, 
       legend = c("state 1", "state 2"), box.lwd = 0)

## ----mllk---------------------------------------------------------------------
mllk = function(theta.star, timediff, x, N=2){
  mu = theta.star[1:N]
  sigma = exp(theta.star[N+1:N])
  Q = diag(N) # generator matrix
  Q[!Q] = exp(theta.star[2*N+1:(N*(N-1))])
  diag(Q) = 0
  diag(Q) = -rowSums(Q)
  delta = solve(t(Q+1), rep(1,N), tol = 1e-20) # stationary distribution of the
  # continuous-time Markov chain
  Qube = LaMa::tpm_cont(Q, timediff) # this computes exp(Q*timediff)
  allprobs = matrix(1, nrow = length(x), ncol = N)
  ind = which(!is.na(x))
  for(j in 1:N){
    allprobs[ind,j] = dnorm(x[ind], mu[j], sigma[j])
  }
  -LaMa::forward_g(delta, Qube, allprobs)
}

## ----model, warning=FALSE-----------------------------------------------------
theta.star = c(5, 15, log(3), log(5), # mu and sigma
                log(1), log(0.5)) # off-diagonals of Q

timediff = diff(obs_times)

t1 = Sys.time()
mod = nlm(mllk, theta.star, timediff=timediff, x=x, stepmax = 10)
# we often need the stepmax, as the matrix exponential can be numerically unstable
Sys.time()-t1

## ----results------------------------------------------------------------------
N = 2
# mu
round(mod$estimate[1:N],2)
# sigma
round(exp(mod$estimate[N+1:N]))
Q = diag(N) # generator matrix
Q[!Q] = exp(mod$estimate[2*N+1:(N*(N-1))])
diag(Q) = 0
diag(Q) = -rowSums(Q)
round(Q,3)

## ----parameters2--------------------------------------------------------------
# 2-state example

# generator matrix Q:
Q = matrix(c(-0.5,0.2,0.3,
             1,-2, 1,
             0.4, 0.6, -1), nrow = 3, byrow = TRUE)

# parameters for the state-dependent (normal) distributions
mu = c(5, 15, 30)
sigma = c(2, 3, 5)

## ----data2--------------------------------------------------------------------
set.seed(123)

k = 200 # number of state switches
trans_times = s = rep(NA, k) # time points where the chain transitions
s[1] = sample(1:3, 1) # uniform initial distribuion
# exponentially distributed waiting times
trans_times[1] = rexp(1, -Q[s[1],s[1]])
n_arrivals = rpois(1, trans_times[1])
obs_times = sort(runif(n_arrivals, 0, trans_times[1]))
x = rnorm(n_arrivals, mu[s[1]], sigma[s[1]])
for(t in 2:k){
  # off-diagonal elements of the s[t-1] row of Q divided by the diagonal element
  # give the probabilites of the next state
  s[t] = sample(c(1:3)[-s[t-1]], 1, prob = Q[s[t-1],-s[t-1]]/-Q[s[t-1],s[t-1]])
  # exponentially distributed waiting times
  trans_times[t] = trans_times[t-1] + rexp(1, -Q[s[t], s[t]])
  n_arrivals = rpois(1, trans_times[t]-trans_times[t-1])
  obs_times = c(obs_times, 
                sort(runif(n_arrivals, trans_times[t-1], trans_times[t])))
  x = c(x, rnorm(n_arrivals, mu[s[t]], sigma[s[t]]))
}

## ----model2, warning=FALSE----------------------------------------------------
theta.star = c(5, 10, 25, log(2), log(2), log(6), # mu and sigma
                rep(0, 6)) # off-diagonals of Q

timediff = diff(obs_times)

t1 = Sys.time()
mod2 = nlm(mllk, theta.star, timediff=timediff, x=x, N = 3, stepmax = 10)
Sys.time()-t1

## ----results2-----------------------------------------------------------------
N = 3
# mu
round(mod2$estimate[1:N],2)
# sigma
round(exp(mod2$estimate[N+1:N]),2)
Q = diag(N) # generator matrix
Q[!Q] = exp(mod2$estimate[2*N+1:(N*(N-1))])
diag(Q) = 0
diag(Q) = -rowSums(Q)
round(Q, 3)

