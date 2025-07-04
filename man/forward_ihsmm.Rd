% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forward_algorithms.R
\name{forward_ihsmm}
\alias{forward_ihsmm}
\title{\href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} for hidden semi-Markov models with inhomogeneous state durations and/ or conditional transition probabilities}
\usage{
forward_ihsmm(
  dm,
  omega,
  allprobs,
  trackID = NULL,
  delta = NULL,
  startInd = NULL,
  eps = 1e-10,
  report = TRUE
)
}
\arguments{
\item{dm}{list of length N containing matrices (or vectors) of dwell-time probability mass functions (PMFs) for each state.

If the dwell-time PMFs are constant, the vectors are the PMF of the dwell-time distribution fixed in time. The vector lengths correspond to the approximating state aggregate sizes, hence there should be little probablity mass not covered by these.

If the dwell-time PMFs are inhomogeneous, the matrices need to have n rows, where n is the number of observations. The number of columns again correponds to the size of the approximating state aggregates.

In the latter case, the first \code{max(sapply(dm, ncol)) - 1} observations will not be used because the first approximating transition probability matrix needs to be computed based on the first \code{max(sapply(dm, ncol))} covariate values (represented by \code{dm}).}

\item{omega}{matrix of dimension c(N,N) or array of dimension c(N,N,n) of conditional transition probabilites, also called embedded transition probability matrix.

It contains the transition probabilities given the current state is left. Hence, the diagonal elements need to be zero and the rows need to sum to one. Such a matrix can be constructed using \code{\link{tpm_emb}} and an array using \code{\link{tpm_emb_g}}.}

\item{allprobs}{matrix of state-dependent probabilities/ density values of dimension c(n, N)}

\item{trackID}{trackID optional vector of length n containing IDs

If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
Instead of a single vector \code{delta} corresponding to the initial distribution, a \code{delta} matrix of initial distributions, of dimension c(k,N), can be provided, such that each track starts with it's own initial distribution.}

\item{delta}{optional vector of initial state probabilities of length N

By default, instead of this, the stationary distribution is computed corresponding to the first approximating transition probability matrix of each track is computed. Contrary to the homogeneous case, this is not theoretically motivated but just for convenience.}

\item{startInd}{optional integer index at which the forward algorithm starts. 

When approximating inhomogeneous HSMMs by inhomogeneous HMMs, the first transition probability matrix that can be constructed is at time \code{max(sapply(dm, ncol))} (as it depends on the previous covariate values).
Hence, when not provided, \code{startInd} is chosen to be \code{max(sapply(dm, ncol))}. Fixing \code{startInd} at a value \strong{larger} than max(aggregate sizes) is useful when models with different aggregate sizes are fitted to the same data and are supposed to be compared. In that case it is important that all models use the same number of observations.}

\item{eps}{small value to avoid numerical issues in the approximating transition matrix construction. Usually, this should not be changed.}

\item{report}{logical, indicating whether initial distribution, approximating transition probability matrix and \code{allprobs} matrix should be reported from the fitted model. Defaults to \code{TRUE}.}
}
\value{
log-likelihood for given data and parameters
}
\description{
Calculates the (approximate) log-likelihood of a sequence of observations under an inhomogeneous hidden semi-Markov model using a modified \strong{forward algorithm}.
}
\details{
Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs, where the state duration distribution is explicitly modelled by a distribution on the positive integers. This function can be used to fit HSMMs where the state-duration distribution and/ or the conditional transition probabilities vary with covariates.
For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.

This function is designed to be used with automatic differentiation based on the \code{R} package \code{RTMB}. It will be very slow without it!
}
\examples{
# currently no examples
}
\references{
Koslik, J. O. (2025). Hidden semi-Markov models with inhomogeneous state dwell-time distributions. Computational Statistics & Data Analysis, 209, 108171.
}
\seealso{
Other forward algorithms: 
\code{\link{forward}()},
\code{\link{forward_g}()},
\code{\link{forward_hsmm}()},
\code{\link{forward_p}()},
\code{\link{forward_phsmm}()}
}
\concept{forward algorithms}
