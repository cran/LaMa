## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "85%",
  fig.align = "center",
  error = TRUE
)

## ----setup--------------------------------------------------------------------
library(LaMa)

## ----data---------------------------------------------------------------------
head(trex, 5)

## ----parameters---------------------------------------------------------------
par = list(
  logmu = log(c(0.3, 1)),      # initial means for step length (log-transformed)
  logsigma = log(c(0.2, 0.7)), # initial sds for step length (log-transformed)
  logkappa = log(c(0.2, 0.7)), # initial concentration for turning angle (log-transformed)
  eta = rep(-2, 2)             # initial t.p.m. parameters (on logit scale)
  )    

dat = list(
  step = trex$step,   # hourly step lengths
  angle = trex$angle, # hourly turning angles
  N = 2
  )

## ----mllk---------------------------------------------------------------------
nll = function(par) {
  getAll(par, dat) # makes everything contained available without $
  Gamma = tpm(eta) # computes transition probability matrix from unconstrained eta
  delta = stationary(Gamma) # computes stationary distribution
  # exponentiating because all parameters strictly positive
  mu = exp(logmu)
  sigma = exp(logsigma)
  kappa = exp(logkappa)
  # reporting statements for later use
  REPORT(mu); ADREPORT(mu)
  REPORT(sigma); ADREPORT(sigma)
  REPORT(kappa); ADREPORT(kappa)
  # calculating all state-dependent densities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle)) # only for non-NA obs.
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j])*dvm(angle[ind],0,kappa[j])
  }
  -forward(delta, Gamma, allprobs) # simple forward algorithm
}

## ----ADfunction---------------------------------------------------------------
obj = MakeADFun(nll, par, silent = TRUE) # creating the objective function

## ----tmbobject----------------------------------------------------------------
names(obj)

## ----tmbobject2---------------------------------------------------------------
obj$par
obj$fn()
obj$gr()

## ----modelfit-----------------------------------------------------------------
opt = nlminb(obj$par, obj$fn, obj$gr) # optimization

## ----optpar-------------------------------------------------------------------
opt$par
opt$objective

## ----MLe----------------------------------------------------------------------
mod = obj$report() # runs the reporting from the negative log-likelihood once
(delta = mod$delta)
(Gamma = mod$Gamma)
(mu = mod$mu)
(sigma = mod$sigma)
(kappa = mod$kappa)

## ----decoding, fig.width = 7, fig.height = 4----------------------------------
# manually
mod$states = viterbi(mod$delta, mod$Gamma, mod$allprobs)

# or simpler
mod$states = viterbi(mod = mod)

# defining color vector
color = c("orange", "deepskyblue")

plot(trex$step[1:200], type = "h", xlab = "time", ylab = "step length", 
     col = color[mod$states[1:200]], bty = "n")
legend("topright", col = color, lwd = 1, legend = c("state 1", "state 2"), bty = "n")

## ----statedepdist, fig.width = 8, fig.height = 4------------------------------
oldpar = par(mfrow = c(1,2))
hist(trex$step, prob = TRUE, breaks = 40, 
     bor = "white", main = "", xlab = "step length")
for(j in 1:2) curve(delta[j] * dgamma2(x, mu[j], sigma[j]), 
                    lwd = 2, add = T, col = color[j])
curve(delta[1]*dgamma2(x, mu[1], sigma[1]) + delta[2]*dgamma2(x, mu[2], sigma[2]), 
      lwd = 2, lty = 2, add = T)
legend("top", lwd = 2, col = color, legend = c("state 1", "state 2"), bty = "n")

hist(trex$angle, prob = TRUE, breaks = 40, 
     bor = "white", main = "", xlab = "turning angle")
for(j in 1:2) curve(delta[j] * dvm(x, 0, kappa[j]), 
                    lwd = 2, add = T, col = color[j])
curve(delta[1]*dvm(x, 0, kappa[1]) + delta[2]*dvm(x, 0, kappa[2]), 
      lwd = 2, lty = 2, add = T)
par(oldpar) # resetting to default

## ----sdreport-----------------------------------------------------------------
sdr = sdreport(obj)

## ----sdreport2----------------------------------------------------------------
summary(sdr)

## ----sdreport3, eval = F------------------------------------------------------
#  # estimated parameter in list format
#  as.list(sdr, "Estimate")
#  # parameter standard errors in list format
#  as.list(sdr, "Std")

## ----sdreport4, eval = F------------------------------------------------------
#  # adreported parameters as list
#  as.list(sdr, "Estimate", report = TRUE)
#  # their standard errors
#  as.list(sdr, "Std", report = TRUE)

## ----pres, fig.width = 8, fig.height = 4--------------------------------------
pres_step = pseudo_res(trex$step, "gamma2", list(mean = mu, sd = sigma), mod = mod)
pres_angle = pseudo_res(trex$angle, "vm", list(mu = 0, kappa = kappa), mod = mod)

oldpar = par(mfrow = c(1,2))
hist(pres_step, prob = TRUE, breaks = 40, 
     bor = "white", main = "pseudo-residuals", xlab = "step length")
curve(dnorm(x, 0, 1), lwd = 2, add = T, lty = 2)
hist(pres_angle, prob = TRUE, breaks = 40, 
     bor = "white", main = "pseudo-residuals", xlab = "turning angle")
curve(dnorm(x, 0, 1), lwd = 2, add = T, lty = 2)
par(oldpar)

## ----tod----------------------------------------------------------------------
Z = cosinor(1:24, period = c(24, 12))

modmat = make_matrices(~ cosinor(tod, period = c(24, 12)), 
                       data = data.frame(tod = 1:24))
Z = modmat$Z

# only compute the 24 unique values and index later for entire time series
dat$Z = Z # adding design matrix to dat
dat$tod = trex$tod # adding time of day to dat for indexing

## ----todpar-------------------------------------------------------------------
par = list(logmu = log(c(0.3, 1)), 
           logsigma = log(c(0.2, 0.7)),
           logkappa = log(c(0.2, 0.7)),
           beta = matrix(c(rep(-2, 2), 
                           rep(0, 2*4)), nrow = 2)) # 2 times 4+1 matrix
# replacing eta with regression parameter matrix, initializing slopes at zero

## ----mllk2--------------------------------------------------------------------
nll2 = function(par) {
  getAll(par, dat) # makes everything contained available without $
  Gamma = tpm_g(Z, beta) # covariate-dependent tpms (in this case only 24 unique)
  # tpm_g() automatically checks if intercept column is included
  ADREPORT(Gamma) # adreporting
  Delta = stationary_p(Gamma) # all periodically stationary distributions
  ADREPORT(Delta)
  delta = Delta[tod[1],] # initial periodically stationary distribution
  # exponentiating because all parameters strictly positive
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)
  # calculating all state-dependent densities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle)) # only for non-NA obs.
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j])*dvm(angle[ind],0,kappa[j])
  }
  -forward_g(delta, Gamma[,,tod], allprobs) # indexing 24 unique tpms by tod in data
}

## ----modelfit2----------------------------------------------------------------
obj2 = MakeADFun(nll2, par, silent = TRUE) # creating the objective function
opt2 = nlminb(obj2$par, obj2$fn, obj2$gr) # optimisation

## ----MLE2, fig.width = 8, fig.height = 5--------------------------------------
mod2 = obj2$report()

sdr = sdreport(obj2)
Gamma = as.list(sdr, "Estimate", report = TRUE)$Gamma
Gammasd = as.list(sdr, "Std", report = TRUE)$Gamma

Delta = as.list(sdr, "Estimate", report = TRUE)$Delta
Deltasd = as.list(sdr, "Std", report = TRUE)$Delta

tod_seq = seq(0, 24, length = 200) # sequence for plotting
Z_pred = trigBasisExp(tod_seq, degree = 2) # design matrix for prediction

Gamma_plot = tpm_g(Z_pred, mod2$beta) # interpolating transition probs

plot(tod_seq, Gamma_plot[1,2,], type = "l", lwd = 2, ylim = c(0,1),
     xlab = "time of day", ylab = "transition probability", bty = "n")
segments(x0 = 1:24, y0 = Gamma[1,2,]-1.96*Gammasd[1,2,], 
         y1 = Gamma[1,2,]+1.96*Gammasd[1,2,])
segments(x0 = 1:24, y0 = Gamma[2,1,]-1.96*Gammasd[2,1,], 
         y1 = Gamma[2,1,]+1.96*Gammasd[2,1,])
lines(tod_seq, Gamma_plot[2,1,], lwd = 2, lty = 3)
legend("topleft", lwd = 2, lty = c(1,3), bty = "n",
       legend = c(expression(gamma[12]^(t)), expression(gamma[21]^(t))))
plot(Delta[,2], type = "b", lwd = 2, xlab = "time of day", ylab = "Pr(active)", 
     col = "deepskyblue", bty = "n", xaxt = "n")
segments(x0 = 1:24, y0 = Delta[,2]-1.96*Deltasd[,2], lwd = 2,
         y1 = Delta[,2]+1.96*Deltasd[,2], col = "deepskyblue")
axis(1, at = seq(0,24,by=4), labels = seq(0,24,by=4))

## ----error, eval = FALSE------------------------------------------------------
#  stop("Invalid argument to 'advector' (lost class attribute?)")

## ----overloading--------------------------------------------------------------
"[<-" <- ADoverload("[<-")

## ----overloading2-------------------------------------------------------------
"c" <- ADoverload("c")
"diag<-" <- ADoverload("diag<-")

## ----NA, eval = FALSE---------------------------------------------------------
#  X = array(dim = c(1,2,3))
#  # which is the same as
#  X = array(NA, dim = c(1,2,3))

## ----NaN, eval = FALSE--------------------------------------------------------
#  X = array(NaN, dim = c(1,2,3))
#  # or
#  X = array(0, dim = c(1,2,3))

## ----max2, eval = FALSE-------------------------------------------------------
#  max2 = function(x,y){
#    (x + y + abs(x - y)) / 2
#  }

## ----bytecompiler, message = FALSE--------------------------------------------
compiler::enableJIT(0)

