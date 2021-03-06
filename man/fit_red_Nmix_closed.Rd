% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_reduced_count_models.R
\name{fit_red_Nmix_closed}
\alias{fit_red_Nmix_closed}
\title{fit_red_Nmix_closed}
\usage{
fit_red_Nmix_closed(nit, lambda_site_covariates = NULL,
  pdet_site_covariates = NULL, red = 1, K, starts = c(1, 0),
  VERBOSE = FALSE, PARALLELIZE = FALSE, APA = FALSE,
  precBits = 128, tolerance = 10^-6, method = "DFP",
  outFile = NULL, ...)
}
\arguments{
\item{nit}{R by T matrix of full counts with R sites/rows and T sampling occassions/columns.}

\item{lambda_site_covariates}{Either NULL (no lambda site covariates) or a list of vectors of length R, where each vector represents one site covariate, and where the vector entries correspond to covariate values for each site. Note that the covariate structure is assumed to be log(lambda_it) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...}

\item{pdet_site_covariates}{Either NULL (no pdet site covariates) or a list of vectors of length R, where each vector represents one site covariate, and where the vector entries correspond to covariate values for each site. Note that the covariate structure is assumed to be logit(lambda_it) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...}

\item{red}{reduction factor, either a number, or a vector of reduction factors (R sites reductions).}

\item{K}{Upper bound on summations, either a single number, or a matrix of K values, one for each site and time (input the full counts value, eg if K=300 for full counts, K=reduction(300,red) for reduced counts).}

\item{starts}{Vector of starting values for optimize. Has two elements, log(lambda) and logit(pdet).}

\item{VERBOSE}{If true, prints the log likelihood to console at each optim iteration.}

\item{PARALLELIZE}{If true, calculation will be split over threads by sites. Will use as many threads as have been made available (initialize with START_PARALLEL(num_cores)).}

\item{APA}{If true, will use arbitrary precision arithmetic in the likelihood calculations. Use precBits to specify the number of bits of precision.}

\item{precBits}{If APA=TRUE, then this will specify the number of bits of precision.}

\item{tolerance}{specifies tolerance for convergence (defulat is 10^-6), all components of estimated gradient must be less than tolerance for convergence. If APA=TRUE, then tolerance can be made very small (eg 10^-20) using: tolerance=Rmpfr::mpfr(10^-20, precBits=128). NOTE: currently tolerance is only used if method="DFP".#' @param method      Optimization method to use. Default is "DFP", and is the only method implemented for APA. When APA=FALSE, can use method="BFGS" to use optim().}

\item{outFile}{If not NULL, name of file for saving algorithm progress (overwritten at each iteration).}

\item{...}{Additional input for optim.}
}
\description{
Find maximum likelihood estimates for closed population models, with model parameters log(lambda) and logit(pdet).
}
\examples{
START_PARALLEL(num_cores=4)
Y    <- gen_Nmix_closed(8,8,250,0.5)
out  <- fit_red_Nmix_closed(Y$nit, red=10, K=300, starts = c(log(250),boot::logit(0.5)), PARALLELIZE=TRUE)
out2 <- fit_red_Nmix_closed(Y$nit, red=matrix(c(10,10,10,10,20,20,40,40),nrow=8, ncol=8), K=300, starts = c(log(250),boot::logit(0.5)), PARALLELIZE=TRUE)
END_PARALLEL()
}
