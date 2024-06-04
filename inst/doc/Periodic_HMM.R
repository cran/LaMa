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
# parameters
mu = c(4, 14)
sigma = c(3, 5)

L = 48 # half-hourly data: 48 observations per day
beta = matrix(c(-1, 1, -1, -1, 1,
                -2, -1, 2, 2, -2), nrow = 2, byrow = TRUE)
Gamma = tpm_p(seq(1, 48, by = 1), L, beta, degree = 2)
Delta = stationary_p(Gamma)

# having a look at the periodically stationary distribution
color = c("orange", "deepskyblue")
plot(Delta[,1], type = "l", lwd = 3, col = color[1], bty = "n", 
     xlab = "time of day", ylab = "Pr(state 1)")
points(Delta[,1], pch = 19, col = color[1])
# only plotting one state, as the other probability is just 1-delta

## ----data---------------------------------------------------------------------
# simulation
z = rep(1:48, 50) # time of day variable, 50 days
n = length(z)
set.seed(123)
s = x = rep(NA, n)
s[1] = sample(1:2, 1, prob = Delta[z[1],])
x[1] = stats::rnorm(1, mu[s[1]], sigma[s[1]])
for(t in 2:n){
  s[t] = sample(1:2, 1, prob = Gamma[s[t-1],,z[t]])
  x[t] = rnorm(1, mu[s[t]], sigma[s[t]])
}

oldpar = par(mfrow = c(1,2))
plot(x[1:400], bty = "n", pch = 20, ylab = "x", 
     col = c(color[1], color[2])[s[1:400]])
boxplot(x ~ z, xlab = "time of day")
# we see a periodic pattern in the data
par(oldpar)

## ----mllk---------------------------------------------------------------------
mllk = function(theta.star, x, z){
  beta = matrix(theta.star[1:10], nrow = 2) # matrix of coefficients
  Gamma = tpm_p(tod = 1:48, L = 48, beta = beta, degree = 2) # calculating all L tpms
  delta = stationary_p(Gamma, t = z[1]) # periodically stationary start
  mu = theta.star[11:12]
  sigma = exp(theta.star[13:14])
  # calculate all state-dependent probabilities
  allprobs = matrix(1, length(x), 2)
  for(j in 1:2){ allprobs[,j] = stats::dnorm(x, mu[j], sigma[j]) }
  # return negative for minimization
  -forward_p(delta, Gamma, allprobs, z)
}

## ----model, warning=FALSE-----------------------------------------------------
theta.star = c(-1,-2, rep(0, 8), # starting values state process
               4, 14 ,log(3),log(5)) # starting values state-dependent process
s = Sys.time()
mod = nlm(mllk, theta.star, x = x, z = z)
Sys.time()-s

## ----visualization------------------------------------------------------------
# transform parameters to working
beta_hat = matrix(mod$estimate[1:10], nrow = 2)
Gamma_hat = tpm_p(tod = 1:48, L = 48, beta = beta_hat, degree = 2)
Delta_hat = stationary_p(Gamma_hat)
mu_hat = mod$estimate[11:12]
sigma_hat = exp(mod$estimate[13:14])

delta_hat = apply(Delta_hat, 2, mean)

oldpar = par(mfrow = c(1,2))
hist(x, prob = TRUE, bor = "white", breaks = 40, main = "")
curve(delta_hat[1]*dnorm(x, mu_hat[1], sigma_hat[1]), add = TRUE, lwd = 2, 
      col = color[1], n=500)
curve(delta_hat[2]*dnorm(x, mu_hat[2], sigma_hat[2]), add = TRUE, lwd = 2, 
      col = color[2], n=500)
curve(delta_hat[1]*dnorm(x, mu_hat[1], sigma_hat[1])+
        delta_hat[2]*dnorm(x, mu[2], sigma_hat[2]),
      add = TRUE, lwd = 2, lty = "dashed", n = 500)
legend("topright", col = c(color[1], color[2], "black"), lwd = 2, bty = "n",
       lty = c(1,1,2), legend = c("state 1", "state 2", "marginal"))

plot(Delta_hat[,1], type = "l", lwd = 3, col = color[1], bty = "n", 
     xlab = "time of day", ylab = "Pr(state 1)")
points(Delta_hat[,1], pch = 19, col = color[1])
par(oldpar)

## ----cSpline------------------------------------------------------------------
nk = 8 # number of basis functions
tod = 1:48
L = 48
k = L * 0:nk / nk # equidistant knots
Z = mgcv::cSplineDes(tod, k) ## cyclic spline design matrix

# plotting the B-Spline basis functions
plot(Z[,1], type = "l", lwd = 2, col = 1, bty = "n",
     xlab = "time of day", ylab = "basis functions", ylim = c(0,0.8))
for(i in 2:nk){
  lines(Z[,i], lwd = 2, col = i)
} 

## ----mllk2--------------------------------------------------------------------
mllk_np = function(theta.star, x, z, Z, lambda){
  beta = matrix(theta.star[1:(2+2*nk)], nrow = 2) # nk params per off-diagonal element
  Gamma = tpm_p(tod = 1:48, L = 48, beta = beta, Z = Z) # calculating all L tpms
  delta = stationary_p(Gamma, t = z[1]) # periodically stationary HMM
  mu = theta.star[2+2*nk + 1:2]
  sigma = exp(theta.star[2+2*nk + 2 + 1:2])
  # calculate all state-dependent probabilities
  allprobs = matrix(1, length(x), 2)
  for(j in 1:2){ allprobs[,j] = stats::dnorm(x, mu[j], sigma[j]) }
  # return negative for minimization
  l = forward_p(delta, Gamma, allprobs, z)
  # penalize curvature
  penalty = sum(diff(beta[1,-1], differences = 2)^2)+
    sum(diff(beta[2,-1], differences = 2)^2)
  return(-l + lambda*penalty)
}

## ----model2, warning=FALSE----------------------------------------------------
theta.star = c(-1,-2, rep(0, 2*nk), # starting values state process
               4, 14 ,log(3),log(5)) # starting values state-dependent process
s = Sys.time()
mod_np = nlm(mllk_np, theta.star, x = x, z = z, Z = Z, lambda = 0)
# in this case we don't seem to need a lot of penalization
Sys.time()-s

## ----visualization2-----------------------------------------------------------
# transform parameters to working
beta_hat_np = matrix(mod_np$estimate[1:(2+2*nk)], nrow = 2)
Gamma_hat_np = tpm_p(tod = 1:48, L = 48, beta = beta_hat_np, Z = Z)
Delta_hat_np = stationary_p(Gamma_hat_np)

# comparing the two fits
plot(Delta_hat_np[,1], type = "l", lwd = 3, col = "purple", bty = "n", 
     xlab = "time of day", ylab = "Pr(state 1)")
# parametric fit
lines(Delta_hat[,1], lwd = 3, col = color[1])

