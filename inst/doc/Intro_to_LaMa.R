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

## ----setup--------------------------------------------------------------------
# loading the package
library("LaMa")

## ----data---------------------------------------------------------------------
# parameters
mu = c(0, 6)
sigma = c(2, 4)
Gamma = matrix(c(0.95, 0.05, 0.15, 0.85), nrow = 2, byrow = TRUE)
delta = stationary(Gamma) # stationary HMM

# simulation
n = 1000
set.seed(123)
s = x = rep(NA, n)
s[1] = sample(1:2, 1, prob = delta)
x[1] = rnorm(1, mu[s[1]], sigma[s[1]])
for(t in 2:n){
  # we draw the next state conditional on the last one
  s[t] = sample(1:2, 1, prob = Gamma[s[t-1],]) 
  # we draw the observation conditional on the current state
  x[t] = rnorm(1, mu[s[t]], sigma[s[t]])
}

color = c("orange", "deepskyblue")
plot(x[1:200], bty = "n", pch = 20, ylab = "x", 
     col = color[s[1:200]])

## ----mllk---------------------------------------------------------------------
mllk = function(theta.star, x){
  # parameter transformations for unconstraint optimization
  Gamma = LaMa::tpm(theta.star[1:2])
  delta = LaMa::stationary(Gamma) # stationary HMM
  mu = theta.star[3:4]
  sigma = exp(theta.star[5:6])
  # calculate all state-dependent probabilities outside the forward algorithm
  allprobs = matrix(1, length(x), 2)
  for(j in 1:2){ allprobs[,j] = stats::dnorm(x, mu[j], sigma[j]) }
  # return negative for minimization
  -LaMa::forward(delta, Gamma, allprobs)
}

## ----model, warning=FALSE-----------------------------------------------------
theta.star = c(-1,-1,1,4,log(1),log(3)) 
# initial transformed parameters: not chosen too well
s = Sys.time()
mod = nlm(mllk, theta.star, x = x)
Sys.time()-s

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
for(j in 1:2){ allprobs[,j] = dnorm(x, mu[j], sigma[j]) }

states = viterbi(delta, Gamma, allprobs)

plot(x, pch = 20, bty = "n", col = color[states])
legend("topright", pch = 20, legend = c("state 1", "state 2"), 
       col = color, box.lwd = 0)

