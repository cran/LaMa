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
mu = c(4, 14)   # state-dependent means
sigma = c(3, 5) # state-dependent standard deviations

L = 48 # half-hourly data: 48 observations per day
beta = matrix(c(-1, 1, -1, -1, 1,
                -2, -1, 2, 2, -2), nrow = 2, byrow = TRUE)
Gamma = tpm_p(seq(1, 48, by = 1), L, beta, degree = 2)
Delta = stationary_p(Gamma)

# having a look at the periodically stationary distribution
color = c("orange", "deepskyblue")
plot(Delta[,1], type = "b", lwd = 2, pch = 16, col = color[1], bty = "n", 
     xlab = "time of day", ylab = "Pr(state 1)")
# only plotting one state, as the other probability is just 1-delta

## ----data---------------------------------------------------------------------
# simulation
tod = rep(1:48, 50) # time of day variable, 50 days
n = length(tod)
set.seed(123)
s = rep(NA, n)
s[1] = sample(1:2, 1, prob = Delta[tod[1],]) # initial state from stationary dist
for(t in 2:n){
  # sampling next state conditional on previous one and the periodic t.p.m.
  s[t] = sample(1:2, 1, prob = Gamma[s[t-1],,tod[t]])
}
# sampling observations conditional on the states
x = rnorm(n, mu[s], sigma[s])

oldpar = par(mfrow = c(1,2))
plot(x[1:400], bty = "n", pch = 20, ylab = "x", 
     col = color[s[1:400]])
boxplot(x ~ tod, xlab = "time of day")
# we see a periodic pattern in the data
par(oldpar)

## ----mllk---------------------------------------------------------------------
nll = function(par, x, tod){
  beta = matrix(par[1:10], nrow = 2) # matrix of coefficients
  Gamma = tpm_p(tod = 1:48, L = 48, beta = beta, degree = 2) # calculating all L tpms
  delta = stationary_p(Gamma, t = tod[1]) # periodically stationary start
  mu = par[11:12]
  sigma = exp(par[13:14])
  # calculate all state-dependent probabilities
  allprobs = matrix(1, length(x), 2)
  for(j in 1:2) allprobs[,j] = dnorm(x, mu[j], sigma[j])
  # return negative for minimization
  -forward_p(delta, Gamma, allprobs, tod)
}

## ----model, warning=FALSE-----------------------------------------------------
par = c(beta = c(-1,-2, rep(0, 8)), # starting values state process
        mu = c(4, 14), # initial state-dependent means
        logsigma = c(log(3),log(5))) # initial state-dependent sds
system.time(
  mod <- nlm(nll, par, x = x, tod = tod)
)

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

plot(Delta_hat[,1], type = "b", lwd = 2, pch = 16, col = color[1], bty = "n", 
     xlab = "time of day", ylab = "Pr(state 1)")
par(oldpar)

## ----cosinor------------------------------------------------------------------
tod = 1:24 # cyclic time of day variable
Z = cosinor(tod, period = c(24, 12)) # design matrix
Z = cbind(intercept = 1, Z)
head(Z, 2)

## ----cosinor_form-------------------------------------------------------------
data = data.frame(tod = rep(1:24, 2), 
                  temp = rnorm(48, 20, 5))
modmat = make_matrices(~ temp * cosinor(tod, 24), data)
Z = modmat$Z
head(Z, 2)

## ----tpm----------------------------------------------------------------------
# coefficient matrix
(beta = matrix(c(-2,-2, runif(2*(ncol(Z)-1))), nrow = 2))
# constructing t.p.m.s
Gamma = tpm_p(Z = Z, beta = beta) # not first arguments in tpm_p
Gamma = tpm_g(Z, beta) # but first arguments in tpm_g

