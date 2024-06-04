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
# state-dependent rates
lambda = c(2, 15)
# generator matrix of the underlying Markov chain
Q = matrix(c(-0.5,0.5,2,-2), nrow = 2, byrow = TRUE)

## ----simulation---------------------------------------------------------------
set.seed(123)

k = 200 # number of state switches
trans_times = s = rep(NA, k) # time points where the chain transitions
s[1] = sample(1:2, 1) # initial distribuion c(0.5, 0.5)
# exponentially distributed waiting times
trans_times[1] = rexp(1, -Q[s[1],s[1]])
# in a fixed interval, the number of arrivals is Pois(lambda * interval_length)
n_arrivals = rpois(1, lambda[s[1]]*trans_times[1]) 
# arrival times within fixed interval are uniformly distributed
arrival_times = runif(n_arrivals, 0, trans_times[1])
for(t in 2:k){
  s[t] = c(1,2)[-s[t-1]] # for 2-states, always a state swith when transitioning
  # exponentially distributed waiting times
  trans_times[t] = trans_times[t-1] + rexp(1, -Q[s[t], s[t]])
  # in a fixed interval, the number of arrivals is Pois(lambda * interval_length)
  n_arrivals = rpois(1, lambda[s[t]]*(trans_times[t]-trans_times[t-1]))
  # arrival times within fixed interval are uniformly distributed
  arrival_times = c(arrival_times, 
                    runif(n_arrivals, trans_times[t-1], trans_times[t]))
}
arrival_times = sort(arrival_times)

## ----vis_MMPP-----------------------------------------------------------------
n = length(arrival_times)
color = c("orange", "deepskyblue")
plot(arrival_times[1:100], rep(0.5,100), type = "h", bty = "n", ylim = c(0,1), 
     yaxt = "n", xlab = "arrival times", ylab = "")
segments(x0 = c(0,trans_times[1:98]), x1 = trans_times[1:99], 
         y0 = rep(0,100), y1 = rep(0,100), col = color[s[1:99]], lwd = 4)
legend("top", lwd = 2, col = color, legend = c("state 1", "state 2"), box.lwd = 0)

## ----mllk---------------------------------------------------------------------
mllk = function(theta.star, timediff, N=2){
  lambda = exp(theta.star[1:N]) # state specific rates
  Q = diag(N) # generator matrix
  Q[!Q] = exp(theta.star[N+1:(N*(N-1))])
  diag(Q) = 0
  diag(Q) = -rowSums(Q)
  Qube = LaMa::tpm_cont(Q-diag(lambda), timediff) # exp((Q-Lambda)*deltat)
  allprobs = matrix(lambda, nrow = length(timediff+1), ncol = N, byrow = T)
  allprobs[1,] = 1
  delta = solve(t(Q+1), rep(1,N), tol = 1e-20)
  -LaMa::forward_g(delta, Qube, allprobs)
}

## ----model, warning=FALSE-----------------------------------------------------
theta.star = log(c(2, 15, # lambda
                   2, 0.5)) # off-diagonals of Q

timediff = diff(arrival_times)

t1 = Sys.time()
mod = nlm(mllk, theta.star, timediff=timediff, stepmax = 10)
# we often need the stepmax, as the matrix exponential can be numerically unstable
Sys.time()-t1

## ----results------------------------------------------------------------------
exp(mod$estimate)

## ----parameters2--------------------------------------------------------------
# state-dependent rates
lambda = c(1, 5, 20)
# generator matrix of the underlying Markov chain
Q = matrix(c(-0.5,0.3,0.2,
             0.7, -1, 0.3,
             1 ,1,-2), nrow = 3, byrow = TRUE)
# parmeters for distributions of state-dependent marks
# (here normally distributed)
mu = c(-5, 0, 5)
sigma = c(2, 1, 2)

color = c("orange", "deepskyblue", "seagreen2")
curve(dnorm(x, 0, 1), xlim = c(-10,10), bty = "n", lwd = 2, col = color[2], 
      n = 200, ylab = "density", xlab = "mark")
curve(dnorm(x, -5, 2), add = TRUE, lwd = 2, col = color[1], n = 200)
curve(dnorm(x, 5, 2), add = TRUE, lwd = 2, col = color[3], n = 200)

## ----simulation2--------------------------------------------------------------
set.seed(123)
k = 200 # number of state switches
trans_times = s = rep(NA, k) # time points where the chain transitions
s[1] = sample(1:3, 1) # initial distribuion uniformly
# exponentially distributed waiting times
trans_times[1] = rexp(1, -Q[s[1],s[1]])
# in a fixed interval, the number of arrivals is Pois(lambda * interval_length)
n_arrivals = rpois(1, lambda[s[1]]*trans_times[1]) 
# arrival times within fixed interval are uniformly distributed
arrival_times = runif(n_arrivals, 0, trans_times[1])
# marks are iid in interval, given underlying state
marks = rnorm(n_arrivals, mu[s[1]], sigma[s[1]])

for(t in 2:k){
  # off-diagonal elements of the s[t-1] row of Q divided by the diagonal element
  # give the probabilites of the next state
  s[t] = sample(c(1:3)[-s[t-1]], 1, prob = Q[s[t-1],-s[t-1]]/-Q[s[t-1],s[t-1]])
  # exponentially distributed waiting times
  trans_times[t] = trans_times[t-1] + rexp(1, -Q[s[t],s[t]])
  # in a fixed interval, the number of arrivals is Pois(lambda * interval_length)
  n_arrivals = rpois(1, lambda[s[t]]*(trans_times[t]-trans_times[t-1]))
  # arrival times within fixed interval are uniformly distributed
  arrival_times = c(arrival_times, 
                    runif(n_arrivals, trans_times[t-1], trans_times[t]))
  # marks are iid in interval, given underlying state
  marks = c(marks, rnorm(n_arrivals, mu[s[t]], sigma[s[t]]))
}
arrival_times = sort(arrival_times)

## ----vis_MMMPP----------------------------------------------------------------
n = length(arrival_times)
plot(arrival_times[1:100], marks[1:100], pch = 16, bty = "n", 
     ylim = c(-9,9), xlab = "arrival times", ylab = "marks")
segments(x0 = c(0,trans_times[1:98]), x1 = trans_times[1:99], 
         y0 = rep(-9,100), y1 = rep(-9,100), col = color[s[1:99]], lwd = 4)
legend("topright", lwd = 2, col = color, 
       legend = c("state 1", "state 2", "state 3"), box.lwd = 0)

## ----mllk2--------------------------------------------------------------------
mllk = function(theta.star, y, timediff, N){
  lambda = exp(theta.star[1:N]) # state specific rates
  mu = theta.star[N+1:N]
  sigma = exp(theta.star[2*N+1:N])
  Q = diag(N) # generator matrix
  Q[!Q] = exp(theta.star[3*N+1:(N*(N-1))])
  diag(Q) = 0
  diag(Q) = -rowSums(Q)
  delta = solve(t(Q+1), rep(1,N), tol = 1e-20)
  Qube = LaMa::tpm_cont(Q-diag(lambda), timediff) # exp((Q-Lambda)*deltat)
  allprobs = matrix(1, length(y), N)
  for(j in 1:N){
    allprobs[,j] = dnorm(y, mu[j], sigma[j])
  }
  allprobs[-1,] = allprobs[-1,] * matrix(lambda, length(y)-1, N, byrow = T)
  -LaMa::forward_g(delta, Qube, allprobs)
}

## ----model2, warning=FALSE----------------------------------------------------
theta.star = c(log(c(1, 5, 20)), # lambda
                   -5, 0, 5, # mu
               log(c(2, 1, 2)), # sigma
               log(c(0.7, 1, 0.3, 1, 0.2, 0.3))) # Q
timediff = diff(arrival_times)
t1 = Sys.time()
mod2 = nlm(mllk, theta.star, y = marks, timediff=timediff, N=3, stepmax = 5)
Sys.time()-t1

## ----results2-----------------------------------------------------------------
N = 3
round(exp(mod2$estimate[1:N]),2)
# mu
round(mod2$estimate[N+1:N], 2)
# sigma
round(exp(mod2$estimate[2*N+1:N]), 2)
Q = diag(N)
Q[!Q] = exp(mod2$estimate[3*N+1:(N*(N-1))])
diag(Q) = 0
diag(Q) = -rowSums(Q)
round(Q, 3)

