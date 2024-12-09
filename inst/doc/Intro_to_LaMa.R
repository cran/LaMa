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

library(LaMa)

## ----tpm----------------------------------------------------------------------
(Gamma = tpm(c(-2, -3))) # 2 states -> 2*(1-2) = 2 off-diagonal entries

## ----stationary---------------------------------------------------------------
(delta = stationary(Gamma))

## ----data---------------------------------------------------------------------
# parameters
mu = c(0, 6)    # state-dependent means
sigma = c(2, 4) # state-dependent standard deviations
Gamma = matrix(c(0.95, 0.05, 0.15, 0.85), # transition probability matrix
               nrow = 2, byrow = TRUE)
delta = stationary(Gamma) # stationary distribution

# simulation
n = 1000
set.seed(123)
s = rep(NA, n)
s[1] = sample(1:2, 1, prob = delta) # sampling first state from delta
for(t in 2:n){
  # drawing the next state conditional on the last one
  s[t] = sample(1:2, 1, prob = Gamma[s[t-1],]) 
}
# drawing the observation conditional on the states
x = rnorm(n, mu[s], sigma[s])

color = c("orange", "deepskyblue")
plot(x[1:200], bty = "n", pch = 20, ylab = "x", 
     col = color[s[1:200]])

## ----mllk---------------------------------------------------------------------
nll = function(par, x){
  # parameter transformations for unconstrained optimisation
  Gamma = tpm(par[1:2]) # multinomial logistic link
  delta = stationary(Gamma) # stationary initial distribution
  mu = par[3:4] # no transformation needed
  sigma = exp(par[5:6]) # strictly positive
  # calculating all state-dependent probabilities outside the forward algorithm
  allprobs = matrix(1, length(x), 2)
  for(j in 1:2) allprobs[,j] = dnorm(x, mu[j], sigma[j])
  # return negative for minimisation
  -forward(delta, Gamma, allprobs)
}

## ----model, warning=FALSE-----------------------------------------------------
par = c(logitGamma = qlogis(c(0.05, 0.05)),
        mu = c(1,4),
        logsigma = c(log(1),log(3)))
# initial transformed parameters: not chosen too well
system.time(
  mod <- nlm(nll, par, x = x)
)

## ----visualization------------------------------------------------------------
# transform parameters to working
Gamma = tpm(mod$estimate[1:2])
delta = stationary(Gamma) # stationary HMM
mu = mod$estimate[3:4]
sigma = exp(mod$estimate[5:6])

hist(x, prob = TRUE, bor = "white", breaks = 40, main = "")
curve(delta[1]*dnorm(x, mu[1], sigma[1]), add = TRUE, lwd = 2, col = "orange", n=500)
curve(delta[2]*dnorm(x, mu[2], sigma[2]), add = TRUE, lwd = 2, col = "deepskyblue", n=500)
curve(delta[1]*dnorm(x, mu[1], sigma[1])+delta[2]*dnorm(x, mu[2], sigma[2]),
      add = TRUE, lwd = 2, lty = "dashed", n=500)
legend("topright", col = c(color, "black"), lwd = 2, bty = "n",
       lty = c(1,1,2), legend = c("state 1", "state 2", "marginal"))

## ----states-------------------------------------------------------------------
allprobs = matrix(1, length(x), 2)
for(j in 1:2) allprobs[,j] = dnorm(x, mu[j], sigma[j])

states = viterbi(delta, Gamma, allprobs)

plot(x, pch = 20, bty = "n", col = color[states])
legend("topright", pch = 20, legend = c("state 1", "state 2"), 
       col = color, box.lwd = 0)

## ----stateprobs---------------------------------------------------------------
probs = stateprobs(delta, Gamma, allprobs)

## ----pseudores----------------------------------------------------------------
pres = pseudo_res(x, # observations
                  "norm", # parametric distribution to use
                  list(mean = mu, sd = sigma), # parameters for that distribution
                  probs) # local state probabilities

oldpar = par(mfrow = c(1,2))
hist(pres, prob = TRUE, bor = "white")
curve(dnorm(x), lty = 2, add = TRUE)
qqnorm(pres, pch = 16, col = "#00000020", bty = "n")
qqline(pres, col = "orange")
par(oldpar)

