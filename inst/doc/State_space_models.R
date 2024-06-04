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

## ----data---------------------------------------------------------------------
beta = 2 # baseline standard deviation
phi = 0.95 # AR parameter
sigma = 0.5 # variability of the AR process

n = 1000
set.seed(123)
g = y = rep(NA, n)
g[1] = rnorm(1, 0, sigma / sqrt(1-phi^2)) # stationary distribution of AR process
y[1] = stats::rnorm(1, 0, beta*exp(g[1]/2)) 
# conditional distribution of y_1 given underlying volatility
for(t in 2:n){
  g[t] = rnorm(1, phi*g[t-1] , sigma) # transition density
  y[t] = stats::rnorm(1, 0, beta*exp(g[t]/2)) 
  # conditional distribution of y_t given underlying volatility
}

# share returns
oldpar = par(mar = c(5,4,3,4.5)+0.1)
plot(y, type = "l", bty = "n", ylim = c(-40,20), yaxt = "n")
# true underlying standard deviation
lines(beta*exp(g)/7 - 40, col = "deepskyblue", lwd = 2)
axis(side=2, at = seq(-20,20,by=5), labels = seq(-20,20,by=5))
axis(side=4, at = seq(0,150,by=75)/7-40, labels = seq(0,150,by=75))
mtext("standard deviation", side=4, line=3, at = -30)
par(oldpar)

## ----mllk---------------------------------------------------------------------
mllk = function(theta.star, y, bm, m){
  phi = plogis(theta.star[1])
  sigma = exp(theta.star[2])
  beta = exp(theta.star[3])
  b = seq(-bm, bm, length = m+1) # intervals for midpoint quadrature
  h = b[2]-b[1] # interval width
  bstar = (b[-1] + b[-(m+1)])/2 # interval midpoints
  # approximation resulting from midpoint quadrature
  Gamma = sapply(bstar, dnorm, mean = phi*bstar, sd = sigma) * h
  Gamma = Gamma / rowSums(Gamma) # normalizing out approximation errors
  delta = h * dnorm(bstar, 0, sigma/sqrt(1-phi^2)) # stationary distribution
  # approximating state-dependent density based on midpoints
  allprobs = t(sapply(y, dnorm, mean = 0, sd = beta * exp(bstar/2)))
  # return negative for minimization
  -forward(delta, Gamma, allprobs)
}

## ----model, warning=FALSE-----------------------------------------------------
theta.star = c(qlogis(0.95), log(0.3), log(1))
bm = 5 # relevant range of underlying volatility (-5,5)
m = 50 # number of approximating states

t1 = Sys.time()
mod = stats::nlm(mllk, theta.star, y = y, bm = bm, m = m)
Sys.time()-t1

## ----results------------------------------------------------------------------
# parameter estimates
(phi = plogis(mod$estimate[1]))
(sigma = exp(mod$estimate[2]))
(beta = exp(mod$estimate[3]))

# decoding states
b = seq(-bm, bm, length = m+1) # intervals for midpoint quadrature
h = b[2]-b[1] # interval width
bstar = (b[-1] + b[-(m+1)])/2 # interval midpoints
Gamma = sapply(bstar, dnorm, mean = phi*bstar, sd = sigma) * h
Gamma = Gamma / rowSums(Gamma) # normalizing out approximation errors
delta = h * dnorm(bstar, 0, sigma/sqrt(1-phi^2)) # stationary distribution
# approximating state-dependent density based on midpoints
allprobs = t(sapply(y, dnorm, mean = 0, sd = beta * exp(bstar/2)))

probs = stateprobs(delta, Gamma, allprobs)
states = viterbi(delta, Gamma, allprobs)

oldpar = par(mar = c(5,4,3,4.5)+0.1)
plot(y, type = "l", bty = "n", ylim = c(-50,20), yaxt = "n")
# when there are so many states it is not too sensable to only plot the most probable state,
# as its probability might still be very small. Generally, we are approximating continuous 
# distributions, thus it makes sense to plot the entire conditional distribution.
maxprobs = apply(probs, 1, max)
for(t in 1:1000){
  colend = round((probs[t,]/(maxprobs[t]*5))*100)
  colend[which(colend<10)] = paste0("0", colend[which(colend<10)])
  points(rep(t, m), bstar*4-35, col = paste0("#FFA200",colend), pch = 20)
}
# we can add the viterbi decoded volatility levels as a "mean"
lines(bstar[states]*4-35)

axis(side=2, at = seq(-20,20,by=5), labels = seq(-20,20,by=5))
axis(side=4, at = seq(-5,5, by = 5)*4-35, labels = seq(-5,5, by = 5))
mtext("g", side=4, line=3, at = -30)
par(oldpar)

