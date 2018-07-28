#TODO: sim study to check if threshold is accurate representation of difference in variance of estimates,

#' Recommend a reduction value \eqn{r} based on maximum threshold on increased variance.
#'
#' @param nit Full counts.
#' @param threshold Desired maximum ratio of \eqn{v2/v1}, where \eqn{v2} is the variance of the full counts,
#'                  and \eqn{v1} is the variance of the reduced counts, scaled by \eqn{r^2}.
#' @return Recommended maximum reduction factor r to use when fitting reduced count model with given sample nit.
#' @examples
#' Y <- redNMix::gen_Nmix_closed(num_sites = 5,num_times = 5,lambda = 250,pdet = 0.5)
#' max_r(nit = Y$nit)
#'
#' Y <- redNMix::gen_Nmix_closed(num_sites = 5,num_times = 5,lambda = 2500,pdet = 0.15)
#' z <- lapply(X = 1:200, FUN = redNMix::reduction, x=Y$nit)
#' z2 <- lapply(X = z, FUN = function(x){var(as.numeric(x))})
#' plot(sqrt(unlist(z2))*1:200,
#'      xlab="reduction factor r",
#'      ylab="standard deviation of scaled reduced counts")
#' abline(h = sqrt(1.1*var(as.numeric(unlist(z[[1]])))))
#' title("sd(nit_r) for r 1:200")
#' @export
max_r <- function(nit, threshold=1.10) {
  r <- 1
  v1 <- var(as.numeric(nit))
  v2 <- 0
  while(v2 < threshold*v1) {
    r  <- r+1
    v2 <- r^2*var(as.numeric(reduction(nit,red=r)))
  }
  return(r-1)
}

#' Generate a population/observation pair with the structure of a closed N-mixture model
#'
#' @param num_sites The number of observation sites.
#' @param num_times The number of sampling occasions.
#' @param lambda    The population rate parameter (\eqn{N_i \sim Poisson(\lambda)})
#' @param pdet      The probability of detection \eqn{p} (\eqn{n_{it} \sim Binomial(N_i,p)})
#' @return A list object with two named matrices. Ni contains the total population per site
#'         (each row represents a site, each column a sampling occasion). nit contains the observed
#'         counts (rows=sites, columns=sampling occasions).
#' @examples
#' pop <- gen_Nmix_closed(num_sites = 5, num_times = 10, lambda = 50, pdetection = 0.4)
#' pop
#' @export
gen_Nmix_closed <- function(num_sites,num_times,lambda,pdet) {
  U    = num_sites
  T    = num_times
  lamb = lambda

  Ntemp <- c(rep(rpois(n=U,lambda = lamb),times=T))
  Ni    <- matrix(data=Ntemp, nrow = U, ncol = T)

  nit   <- Ni
  nit[] <- vapply(Ni, function(x) { rbinom(size = x, n = 1, prob = pdet) }, numeric(1))

  return(list(Ni=Ni, nit=nit))
}

# rounding in R rounds to nearest even number (eg 0.5 rounds to 0), this function does
# "normal" rounding (eg 0.5 rounds to 1)
round2 = function(x, n=0) {
  sgn = sign(x)
  z = abs(x)*10^n
  z = trunc(z + 0.5)
  z = z/10^n
  z*sgn
}

#' Generate a population/observation pair with the structure of an open N-mixture model
#'
#' @param num_sites The number of observation sites.
#' @param num_times The number of sampling occasions.
#' @param lambda    The population rate parameter (\eqn{N_i \sim Poisson(\lambda)}).
#' @param pdet      The probability of detection \eqn{p} (\eqn{n_{it} \sim Binomial(N_i,p)}).
#' @param omega     The probability of survival (\eqn{S_{it} \sim Binomial(N_{it}, \omega)}).
#' @param gamma     The recruitment rate parameter (\eqn{G_{it} \sim Poisson(\gamma)}).
#' @return A list object with two named matrices. Nit contains the total population per site
#'         (each row represents a site, each column a sampling occasion). nit contains the observed
#'         counts (rows=sites, columns=sampling occasions).
#' @examples
#' pop <- gen_Nmix_open(num_sites = 5,
#'                      num_times = 10,
#'                      lambda = 50,
#'                      pdet = 0.4,
#'                      omega = 0.8,
#'                      gamma = 4)
#' pop
#' @export
gen_Nmix_open <- function(num_sites,num_times,lambda,pdet,omega,gamma) {
  U    = num_sites
  T    = num_times
  lamb = lambda
  gamm = gamma
  omeg = omega

  Ntemp <- c(rep(rpois(n=U,lambda = lamb),times=T))
  Ni <- matrix(data=Ntemp, nrow = U, ncol = T)
  for(i in 2:T) {
    Ni[,i] <- rbinom(n = U, size = Ni[,i-1], prob = omeg) + rpois(n = U, lambda = gamm)
  }

  nit <- Ni
  nit[] <- vapply(Ni, function(x) { rbinom(size = x, n = 1, prob = pdet) }, numeric(1))

  return(list(Nit=Ni, nit=nit))
}


#' Reduce an integer value using a reduction function \eqn{R(x;r)}. Currently only implements round2(), which is standard rounding (\{x\}.5 rounds to \{x+1\}.0).
#'
#' @param x The integer which is to be reduced.
#' @param red The factor r by which to reduce the input x.
#' @param FUN The reduction function (default is round2(), some alternatives are ceiling() and floor())
#' @return An integer value which is the reduction of integer x by reduction factor red using function FUN.
#' @examples
#' x <- 104
#' xr1 <- reduction(x, 10)
#' xr2 <- reduction(x, 10, ceiling)
#' @export
reduction <- function(x, red, FUN=round2) {
  FUN(x/red)
}

#' Reduced binomial probability distribution function \eqn{rBinomial(x;N,p,R(x;r))}
#'
#' @param x Full count quantile (alternatively input r*x if x is a reduced count quantile).
#' @param size Number of trials.
#' @param prob Probability of success for each trial.
#' @param red The factor r by which to reduce the input x.
#' @return The probability of observing quantile \eqn{R(x;r)}
#' @examples
#' Y <- drbinom(0:200, 200, 0.3, 10)[seq(1,201,10)]
#' plot(Y, xlab="Y=R(X;10), X~Binomial(N=200,p=0.3)", ylab="P[Y=y]")
#' @export
drbinom <- function(x, size, prob, red, log=FALSE) {
  start <- reduction(x,red, FUN=round2)*red-red/2
  end   <- start + red# reduction(x,red, FUN=round2)*red+red/2

  p <- pbinom(end, size, prob) - pbinom(start, size, prob)

  if(log) {
    return(log(p))
  }
  return(p)
}

#' Reduced poisson probability distribution function \eqn{rPoisson(x;\lambda,r)}
#'
#' @param x Full count quantile (alternatively input r*x if x is a reduced count quantile).
#' @param lambda Mean of the full count poisson distribution.
#' @param red The factor r by which to reduce the input x.
#' @return The probability of observing quantile \eqn{R(x;r)}
#' @examples
#' Y <- drpois(seq(1,200,10), 55, 10)
#' plot(Y, xlab="Y=R(X;10)", main="X~Poisson(lambda=55)", ylab="P[Y=y]")
#' @export
drpois <- function(x, lambda, red, log=FALSE) {
  start <- reduction(x,red, FUN=round2)*red-red/2
  end   <- start + red # reduction(x,red, FUN=round2)*red+red/2

  p <- ppois(end, lambda) - ppois(start, lambda)
  if(log) {
    return(log(p))
  }
  return(p)
}

#' Internal function. Used to calculate the negative of the log likelihood.
#' @param par Vector with two elements, logis(pdet) and log(lambda).
#' @param nit R by T matrix of full counts with R sites/rows and T sampling occassions/columns.
#' @param K   Upper bound on summations.
#' @param red reduction factor
#' @param VERBOSE If true, prints the calculated log likelihood to console.
#' @export
red_Like_closed <- function(par, nit, K, red, FUN=round, VERBOSE=FALSE) {
  T <- ncol(nit)
  R <- nrow(nit)

  pdet <- plogis(par[1])
  lamb <- exp(par[2])
  Y <- FUN(nit/red)
  K <- FUN(K/red)

  l <- 0
  for(i in 1:R) {
    li <- 0
    ni <- max(Y[i,]) # ni

    for(Ni in ni:K) {
      lit <- 1
      for(t in 1:T) {
        lit <- lit*drbinom(x = Y[i,t]*red, size = Ni*red, prob = pdet, red=red)
      }
      li <- li + lit*drpois(x = Ni*red, lambda = lamb, red=red) #lit*dpois(x = Ni, lambda = lamb/red) #
    }
    l <- l+log(li)
  }
  if(VERBOSE) {print(paste0("log likelihood: ",l))}
  return(-1*l)
}


#' Internal function. Used to calculate the negative of the log likelihood.
#' @param par Vector with four elements, log(lambda), log(gamma), logis(omega), and logis(pdet).
#' @param nit R by T matrix of full counts with R sites/rows and T sampling occassions/columns.
#' @param K   Upper bound on summations.
#' @param red reduction factor
#' @param VERBOSE if true, prints the log likelihood to console.
#' @details Note that this function is adapted from the negative log likelihood function from the unmarked package, and uses the recursive method of computation described in Web Appendix A of Dail and Madsen 2011: Models for Estimating Abundance from Repeated Counts of an Open Metapopulation, published in Biometrics volume 67, issue 2.
#' @export
red_Like_open <- function(par, nit, K, red, FUN=round, VERBOSE=FALSE) {
  T <- ncol(nit)
  R <- nrow(nit)

  pdet <- plogis(par[4])
  lamb <- exp(par[1])
  omeg <- plogis(par[3])
  gamm <- exp(par[2])

  Y <- FUN(nit/red)
  K <- FUN(K/red)


  # holds g1[k] * g_star[k]
  g1_t_star <- rep(0, times=K+1)

  # g3 is the transition probability matrix
  g3_ <- matrix(0, nrow = K+1, ncol=K+1)
  g3  <- tp_MAT(M = g3_, omeg = omeg, gamm = gamm, red = red)

  ll <- 0
  # loop over sites
  for(i in 1:R) {

    g_star <- rep(1, times=K+1)
    # loop backwards over times t, stopping at t==2
    for(t in T:2) {

      g1_t <- rep(0, times=K+1)
      # loop over possible value of N at time t
      for(k in 1:K+1) {
        g1_t[k] <- drbinom(x = Y[i,t]*red, size = (k-1)*red, prob = pdet, red = red) #Rf_dbinom(y(i,j,t), k, p(i,j,t), true);
        g1_t_star[k] <- g1_t[k] * g_star[k];
      }

      # update g_star
      g_star = g3 %*% g1_t_star
    }

    ll_i <- 0.0
    g1 <- rep(0, times = K+1)
    g2 <- rep(0, times = K+1)

    # loop over possible values of N at time t==1
    for(k in 1:K+1) {
      g1[k] <- drbinom(x = Y[i,1]*red, size = (k-1)*red, prob = pdet, red = red)
      g2[k] <- drpois(x = (k-1)*red, lambda = lamb, red = red)

      # apply recursive definition of likelihood
      ll_i <- ll_i + g1[k] * g2[k] * g_star[k];
    }

    ll   <- ll + log(ll_i + 1e-320)
  }

  if(VERBOSE) {print(paste0("log likelihood: ",ll))}
  return(-1*ll)
}


#' Internal function, calculates transition probabilities from pop size j to pop size k in the open population likelihood.
tp_jk <- function(j,k,omeg,gamm,red) {
  a <- data.frame(matrix(seq(0,min(j,k),1),ncol = 1))
  prob <- sum(
    apply(X = a, MARGIN = 1, FUN = function(c,j,k,omeg,gamm,red) {
      drbinom(x = c, size = j, prob = omeg, red = red) * drpois(x = k-c, lambda = gamm, red = red)
    }, j=j,k=k,omeg=omeg, gamm=gamm,red=red)
  )
  return(prob)
}

#' Internal function, calculates transition probability matrix (transition from row pop to column pop)
tp_MAT <- function(M, omeg, gamm, red) {
  K1 <- nrow(M)

  for(row in 1:K1) {
    for(col in 1:K1) {
      M[row,col] <- tp_jk(row-1, col-1, omeg=omeg, gamm=gamm, red=red)
    }
  }

  return(M)
}


#' Find maximum likelihood estimates for model parameters logit(pdet) and log(lambda). Uses optimr.
#' @param starts Vector of starting values for optimize. Has two elements, logit(pdet) and log(lambda).
#' @param nit    R by T matrix of full counts with R sites/rows and T sampling occassions/columns.
#' @param K      Upper bound on summations.
#' @param red    reduction factor.
#' @param VERBOSE If true, prints the log likelihood to console at each optim iteration.
#' @param ...    Additional input for optimr.
#' @examples
#' Y <- gen_Nmix_closed(1,5,250,0.5)
#' out <- fit_red_Nmix_closed(Y$nit, red=10, K=300, starts = c(boot::logit(0.5),log(250)))
#' @export
fit_red_Nmix_closed <- function(nit, red, K, starts=c(0,1), VERBOSE=FALSE, ...) {
  require(optimr)
  opt <- optimr(par     = starts,
                fn      = red_Like_closed,
                nit     = nit,
                K       = K,
                red     = red,
                VERBOSE = VERBOSE,
                ...)
  return(opt)
}


#' Find maximum likelihood estimates for model parameters log(lambda), log(gamma), logit(omega), and logit(pdet). Uses optimr.
#' @param starts Vector with four elements, log(lambda), log(gamma), logit(omega), and logit(pdet).
#' @param nit    R by T matrix of full counts with R sites/rows and T sampling occassions/columns.
#' @param K      Upper bound on summations.
#' @param red    reduction factor.
#' @param VERBOSE If true, prints the log likelihood to console at each optim iteration.
#' @param ...    Additional input for optimr.
#' @examples
#' Y <- gen_Nmix_open(num_sites = 3, num_times = 4, lambda = 10, pdet = 0.7, omega = 0.7, gamma = 2)
#' out <- fit_red_Nmix_open(nit = Y$nit, red = 1, K = 30, starts = c(0.5, 0.5, 0.5, 0.5))
#' @export
fit_red_Nmix_open <- function(nit, red, K, starts=c(1,1,0,0), VERBOSE=FALSE, ...) {
  #require(optimr)
  opt <- optim(par      = starts,
                fn      = red_Like_open,
                nit     = nit,
                K       = K,
                red     = red,
                VERBOSE = VERBOSE,
                method = "BFGS",
                ...)
  return(opt)
}

#' Plot likelihood given a pdet and range for lambda.
#' @param nit         R by T matrix of full counts with R sites/rows and T sampling occassions/columns.
#' @param startLambda Starting value for lambda.
#' @param endLambda   Ending value for lambda.
#' @param stepsize    Spacing between values of lambda
#' @param pdet        Probability of detection (pdet = 0.5 means 50\% chance of detection)
#' @param red         Reduction factor.
#' @param K           Upper bound on summations.
#' @examples
#' Y <- gen_Nmix_closed(5,5,250,0.5)
#' plot_red_like_closed_lambda(Y=Y$nit, startLambda = 150, endLambda = 350, stepsize=10, pdet = 0.5, red = 10, K = 400)
#' @export
plot_red_like_closed_lambda <- function(nit, startLambda, endLambda, stepsize, pdet, red, K) {

  parB <- startLambda
  parT <- endLambda
  j <- 1
  L <- NULL
  for(parM in seq(parB,parT,stepsize)) {
    L[j] <- -1*red_Like_closed(nit = nit, par = c(boot::logit(pdet),log(parM)), K = K, red = red)
    j <- j+1
  }

  plot(y=L, x=seq(parB,parT,stepsize), xlab="lambda")
}

#' Plot likelihood given a lambda and range for pdet.
#' @param nit         R by T matrix of full counts with R sites/rows and T sampling occassions/columns.
#' @param startPdet   Starting value for pdet.
#' @param endPdet     Ending value for pdet.
#' @param stepsize    Spacing between values of pdet
#' @param lambda      Initial abundance parameter.
#' @param K   Upper bound on summations.
#' @param red reduction factor.
#' @examples
#' Y <- gen_Nmix_closed(5,5,250,0.5)
#' plot_red_like_closed_pdet(Y=Y$nit, startPdet = 0.1, endPdet = 1.0, stepsize=0.1, lambda = 250, red = 10, K = 400)
#' @export
plot_red_like_closed_pdet <- function(nit, startPdet, endPdet, stepsize, lambda, red, K) {

  parB <- startPdet
  parT <- endPdet
  j <- 1
  L <- NULL
  for(parM in seq(parB,parT,stepsize)) {
    L[j] <- -1*red_Like_closed(nit = nit, par = c(boot::logit(parM),log(lambda)), K = K, red = red)
    j <- j+1
  }

  plot(y=L, x=seq(parB,parT,stepsize), xlab="pdet")
}

#' Plot 2D likelihood given a range for pdet and for lambda.
#' @param nit             R by T matrix of full counts with R sites/rows and T sampling occassions/columns.
#' @param startPdet       Starting value for pdet.
#' @param endPdet         Ending value for pdet.
#' @param stepsizePdet    Spacing between values of pdet
#' @param startLambda     Starting value for lambda.
#' @param endLambda       Ending value for lambda.
#' @param stepsizeLambda  Spacing between values of lambda
#' @param K               Upper bound on summations.
#' @param red             Reduction factor.
#' @return                Returns a matrix of likelihoods where rows represent pdet, and columns represent lambda.
#' @examples
#' Y <- gen_Nmix_closed(5,5,250,0.5)
#' plot_2d_red_like_closed(Y$nit, 0.1, 1.0, 0.25, 100, 400, 100, 10, 400)
#' @export
plot_2d_red_like_closed <- function(nit, startPdet, endPdet, stepsizePdet, startLambda, endLambda, stepsizeLambda, red, K) {
  par1B <- startPdet
  par1T <- endPdet
  par2B <- startLambda
  par2T <- endLambda
  prange <- seq(par1B,par1T,stepsizePdet)
  lrange <- seq(par2B,par2T,stepsizeLambda)
  j <- 1
  L <- matrix(nrow = length(prange), ncol = length(lrange))
  for(par1M in prange) {
    k <- 1
    for(par2M in lrange) {
      L[j,k] <- -1*red_Like_closed(nit = nit, par = c(boot::logit(par1M),log(par2M)), K = K, red = red)
      k <- k+1
    }
    j <- j+1
  }

  L2 <- L
  L[which(L==-Inf)] <- 2*min(L[which(L!=-Inf)])

  #require(graphics)
  layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE),
         widths=c(1,1), heights=c(1,1,1))

  par(mar = c(5,5,5,5))
  image(x = prange, y=lrange, L,
        ylab = "lambda", xlab = "pdet",
        main = "Likelihood Contour")
  box()
  image(x = prange, y=lrange, exp(L),
        ylab = "lambda", xlab = "pdet",
        main = "Exponential Likelihood Contour")
  box()

  par(mar = c(1,1,1,1))
  persp(prange, lrange, L, theta = 30, phi = 30, col = "lightblue", shade = 0.25, ticktype = "detailed", xlab = "pdet", ylab="lambda", zlab = "Likelihood")
  persp(prange, lrange, L, theta = 120, phi = 30, col = "lightblue", shade = 0.25, ticktype = "detailed", xlab = "pdet", ylab="lambda", zlab = "Likelihood")
  persp(prange, lrange, L, theta = 210, phi = 30, col = "lightblue", shade = 0.25, ticktype = "detailed", xlab = "pdet", ylab="lambda", zlab = "Likelihood")
  persp(prange, lrange, L, theta = 300, phi = 30, col = "lightblue", shade = 0.25, ticktype = "detailed", xlab = "pdet", ylab="lambda", zlab = "Likelihood")

  return(L2)
}

