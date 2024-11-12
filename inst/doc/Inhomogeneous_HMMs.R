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
library(LaMa)

## ----parameters---------------------------------------------------------------
# parameters
mu = c(5, 20)   # state-dependent means
sigma = c(4, 5) # state-dependent standard deviations

# state process regression parameters
beta = matrix(c(-2, -2,       # intercepts
                -1, 0.5,      # linear effects
                0.25, -0.25), # quadratic effects
              nrow = 2)

n = 1000 # number of observations
set.seed(123)
z = rnorm(n) # in practice there will be n covariate values.
# However, we only have n-1 transitions, thererfore we only need n-1 values:
Z = cbind(z, z^2) # quadratic effect of z
Gamma = tpm_g(Z = Z[-1,], beta) # of dimension c(2, 2, n-1)
delta = c(0.5, 0.5) # non-stationary initial distribution

color = c("orange", "deepskyblue")

oldpar = par(mfrow = c(1,2))
zseq = seq(-2,2,by = 0.01)
Gamma_seq = tpm_g(Z = cbind(zseq, zseq^2), beta)
plot(zseq, Gamma_seq[1,2,], type = "l", lwd = 3, bty = "n", ylim = c(0,1), 
     xlab = "z", ylab = "gamma_12", col = color[1])
plot(zseq, Gamma_seq[2,1,], type = "l", lwd = 3, bty = "n", ylim = c(0,1),
     xlab = "z", ylab = "gamma_21", col = color[2])
par(oldpar)

## ----data---------------------------------------------------------------------
s = rep(NA, n)
s[1] = sample(1:2, 1, prob = delta) # sampling first state from initial distr.
for(t in 2:n){
  # sampling next state conditional on previous one with tpm at that time point
  s[t] = sample(1:2, 1, prob = Gamma[s[t-1],,t-1])
}
# sampling observations conditional on the states
x = rnorm(n, mu[s], sigma[s])

plot(x[1:200], bty = "n", pch = 20, ylab = "x", 
     col = c(color[1], color[2])[s[1:200]])

## ----mllk---------------------------------------------------------------------
nll = function(par, x, Z){
  beta = matrix(par[1:6], nrow = 2) # matrix of coefficients
  Gamma = tpm_g(Z[-1,], beta) # excluding the first covariate value -> n-1 tpms
  delta = c(1, exp(par[7]))
  delta = delta / sum(delta)
  mu = par[8:9]
  sigma = exp(par[10:11])
  # calculate all state-dependent probabilities
  allprobs = matrix(1, length(x), 2)
  for(j in 1:2) allprobs[,j] = dnorm(x, mu[j], sigma[j])
  # forward algorithm
  -forward_g(delta, Gamma, allprobs)
}

## ----model, warning=FALSE-----------------------------------------------------
par = c(beta = c(-2, -2, rep(0,4)), # initialising with homogeneous tpm
        logitdelta = 0, # starting value for initial distribution
        mu = c(4, 14), # initial state-dependent means
        sigma = c(log(3),log(5))) # initial state-dependents sds
system.time(
  mod <- nlm(nll, par, x = x, Z = Z)

)

## ----visualization------------------------------------------------------------
# transform parameters to working
beta_hat = matrix(mod$estimate[1:6], nrow = 2)
Gamma_hat = tpm_g(Z = Z[-1,], beta_hat)
delta_hat = c(1, exp(mod$estimate[7]))
delta_hat = delta_hat / sum(delta_hat)
mu_hat = mod$estimate[8:9]
sigma_hat = exp(mod$estimate[10:11])

# we calculate the average state distribution overall all covariate values
zseq = seq(-2, 2, by = 0.01)
Gamma_seq = tpm_g(Z = cbind(zseq, zseq^2), beta_hat)
Prob = matrix(nrow = length(zseq), ncol = 2)
for(i in 1:length(zseq)){ Prob[i,] = stationary(Gamma_seq[,,i]) }
prob = apply(Prob, 2, mean)

hist(x, prob = TRUE, bor = "white", breaks = 20, main = "")
curve(prob[1]*dnorm(x, mu_hat[1], sigma_hat[1]), add = TRUE, lwd = 3, 
      col = color[1], n=500)
curve(prob[2]*dnorm(x, mu_hat[2], sigma_hat[2]), add = TRUE, lwd = 3, 
      col = color[2], n=500)
curve(prob[1]*dnorm(x, mu_hat[1], sigma_hat[1])+
        prob[2]*dnorm(x, mu[2], sigma_hat[2]),
      add = TRUE, lwd = 3, lty = "dashed", n = 500)
legend("topright", col = c(color[1], color[2], "black"), lwd = 3, bty = "n",
       lty = c(1,1,2), legend = c("state 1", "state 2", "marginal"))

oldpar = par(mfrow = c(1,2))
plot(zseq, Gamma_seq[1,2,], type = "l", lwd = 3, bty = "n", ylim = c(0,1), 
     xlab = "z", ylab = "gamma_12_hat", col = color[1])
plot(zseq, Gamma_seq[2,1,], type = "l", lwd = 3, bty = "n", ylim = c(0,1),
     xlab = "z", ylab = "gamma_21_hat", col = color[2])
par(mfrow = c(1,1))
plot(zseq, Prob[,1], type = "l", lwd = 3, bty = "n", ylim = c(0,1), xlab = "z", 
     ylab = "Pr(state 1)", col = color[1])
par(oldpar)

## ----parameters2--------------------------------------------------------------
sigma = c(1, 1) # state-dependent standard deviations (homoscedasticity)

# parameter matrix
# each row contains parameter vector for the corresponding state
beta = matrix(c(8, 10,             # intercepts
                -2, 1, 0.5, -0.5), # slopes
              nrow = 2)
n = 1000 # number of observations
set.seed(123)
z = rnorm(n)
Z = cbind(z, z^2) # quadratic effect of z

Gamma = matrix(c(0.9, 0.1, 0.05, 0.95), 
               nrow = 2, byrow = TRUE) # homogeneous t.p.m.
delta = stationary(Gamma) # stationary Markov chain

## ----data2--------------------------------------------------------------------
s = x = rep(NA, n)
s[1] = sample(1:2, 1, prob = delta)
x[1] = rnorm(1, beta[s[1],]%*%c(1, Z[1,]), # state-dependent regression
                    sigma[s[1]])
for(t in 2:n){
  s[t] = sample(1:2, 1, prob = Gamma[s[t-1],])
  x[t] = rnorm(1, beta[s[t],]%*%c(1, Z[t,]), # state-dependent regression
                      sigma[s[t]])
}

oldpar = par(mfrow = c(1,2))
plot(x[1:400], bty = "n", pch = 20, ylab = "x", 
     col = c(color[1], color[2])[s[1:400]])

plot(z[which(s==1)], x[which(s==1)], pch = 16, col = color[1], bty = "n", 
     ylim = c(0,15), xlab = "z", ylab = "x")
points(z[which(s==2)], x[which(s==2)], pch = 16, col = color[2])
par(oldpar)

## ----mllk3--------------------------------------------------------------------
nllMSR = function(par, x, Z){
  Gamma = tpm(par[1:2]) # homogeneous tpm
  delta = stationary(Gamma) # stationary Markov chain
  beta = matrix(par[2 + 1:(2 + 2*2)], nrow = 2) # parameter matrix
  sigma = exp(par[2 + 2 + 2*2 + 1:2])
  # calculate all state-dependent probabilities
  allprobs = matrix(1, length(x), 2)
  # state-dependent regression
  for(j in 1:2) allprobs[,j] = dnorm(x, cbind(1,Z) %*% beta[j,], sigma[j])
  # forward algorithm
  -forward(delta, Gamma, allprobs)
}

## ----model3, warning=FALSE----------------------------------------------------
par = c(logitgamma = c(-2, -3),      # starting values state process
        beta = c(8, 10, rep(0,4)),   # starting values for regression
        logsigma = c(log(1),log(1))) # starting values for sigma
system.time(
  mod_reg <- nlm(nllMSR, par, x = x, Z = Z)
)

## ----visualization3-----------------------------------------------------------
Gamma_hat_reg = tpm(mod_reg$estimate[1:2]) # calculating all tpms
delta_hat_reg = stationary(Gamma_hat_reg)
beta_hat_reg = matrix(mod_reg$estimate[2+1:(2*2+2)], nrow = 2)
sigma_hat_reg = exp(mod_reg$estimate[2+2*2+2 +1:2])

# we have some label switching
plot(z, x, pch = 16, bty = "n", xlab = "z", ylab = "x", col = color[s])
points(z, x, pch = 20)
curve(beta_hat_reg[1,1] + beta_hat_reg[1,2]*x + beta_hat_reg[1,3]*x^2, 
      add = TRUE, lwd = 4, col = color[2])
curve(beta_hat_reg[2,1] + beta_hat_reg[2,2]*x + beta_hat_reg[2,3]*x^2, 
      add = TRUE, lwd = 4, col = color[1])

