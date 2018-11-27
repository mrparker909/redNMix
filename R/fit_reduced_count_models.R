
#' Find the Most Common Divisor of a data set
#' @param nit Full counts
#'
#' @examples
#' mcd(c(10,20,30,35,12,18), howmany = 5)
#'
#' @export
mcd <- function(nit, howmany=10) {
  x <- list()
  for(n in nit) {
    x <- c(x,numbers::divisors(n))
  }
  y <- unlist(x)
  howmany <- min(howmany, length(table(y)))
  sort(table(y),decreasing=TRUE)[1:howmany]
}

#' @title  lurf: Find the Least Unfavourable Reduction Factors
#' @description Tally the remainders upon division for \eqn{r \in [rmin,rmax]}, return a data.frame df sorted by best r. df has two columns, red and rem. red is the reduction factor and rem is the total of the remainders upon reduction by red.
#' @param rmin  minimum reduction factor to calculate (default is 1)
#' @param rmax  maximum reduction factor to calculate (default is 25)
#' @param rstep step size between r values (default is 1)
#' @param nit   count data
#' @examples
#' data <- c(215,309,116,157,213,248,112)
#' # Suppose you want the reduction factor which loses the least information, but which is at least 10:
#' lurf(data, rmin=10, plothist=TRUE) # the histogram shows patterns in changing r
#' # examining the data.frame for smallest SS shows that 11 and 14 will give better results than 10.
#'
#' # we can compare the reduced counts from these three options (r=10,11,14) to the original data:
#' data
#' 10*reduction(data, red = 10)
#' 11*reduction(data, red = 11)
#' 14*reduction(data, red = 14)
#' # this makes it clear why 11 and 14 are better choices of r than 10 for this data (even though they are larger reductions!)
#' @export
lurf <- function(nit, rmin=1, rmax=25, rstep=1, plothist=FALSE) {
  redseq <- seq(rmin,rmax,1)
  remseq <- sapply(X = redseq, FUN = function(X, nit) {
    sum(nit%%X)
  }, nit=nit)
  SSseq <- sapply(X = redseq, FUN = function(X, nit) {
    sum((nit%%X)^2)
  }, nit=nit)
  df <- data.frame(red=redseq[order(SSseq)], rem=remseq[order(SSseq)], SS=sqrt(SSseq)[order(SSseq)])
  if(plothist) {
    plot(y=df$rem, x=df$red)
  }
  return(df)
}

#' Recommend a reduction value \eqn{r} based on maximum threshold on increased variance.
#'
#' @param nit Full counts.
#' @param threshold Desired maximum ratio of \eqn{v1/v2}, where \eqn{v2} is the variance of the full counts,
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
round2_ = function(x, n=0) {
  sgn = sign(x)
  z = abs(x)*10^n
  z = trunc(z + 0.5)
  z = z/10^n
  z*sgn
}
round2 <- compiler::cmpfun(round2_)

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


#' Reduce an integer value using a reduction function \eqn{R(x;r)}.
#' Currently only implements round2(), which is standard rounding (\{x\}.5 rounds to \{x+1\}.0).
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


drbinom_ <- function(x, size, prob, red, log=FALSE) {
  start <- x*red - red/2
  end   <- start + red

  p <- pbinom(end, size*red, prob) - pbinom(start, size*red, prob)
  #if(p < 0) { print(paste("drbinom < 0!", "x=",x,"size=",size, "prob=",prob, "red=",red))}
  if(log) {
    return(log(p))
  }
  return(p)
}
#' Reduced binomial probability distribution function \eqn{rBinomial(x;N,p,R(x;r))},
#' takes reduced quantiles rather than full quantiles (use drbinom2 for full quantiles).
#'
#' @param x Reduced count quantile (alternatively input reduction(x,r) if x is a full count quantile).
#' @param size Number of trials.
#' @param prob Probability of success for each trial.
#' @param red The factor r by which x has been reduced.
#' @return The probability of observing quantile \eqn{x}.
#' @examples
#' Y <- drbinom(0:20, 20, 0.3, 10)
#' plot(Y, xlab="Y~rBinom(N=20,p=0.3,r=10)", ylab="P[Y=y]")
#' @export
drbinom <- compiler::cmpfun(drbinom_)



#' internal function
drpois_1_  <- function(x, lambda, red, log=FALSE) {
  start <- x*red-red/2
  end   <- start + red # reduction(x,red, FUN=round2)*red+red/2
  p <- ppois(end, lambda) - ppois(start, lambda) #sum(dpois(start:end, lambda)) #
  #p[which(p < 0)] <- 0
  if(log) {
    return(log(p))
  }
  return(p)
}
drpois_1 <- compiler::cmpfun(drpois_1_)


#' Reduced poisson probability distribution function \eqn{rPoisson(x;\lambda,r)}, takes reduced quantiles (use drpois2 for full quantiles).
#'
#' @param x Reduced count quantile (alternatively input reduction(x,r) if x is a full count quantile).
#' @param lambda Mean of the full count poisson distribution.
#' @param red The factor r by which x was reduced.
#' @return The probability of observing quantile \eqn{x}
#' @examples
#' # probability of observing 105 from pois(lam=90)
#' x <- dpois(x = 105, lambda = 90)
#' # probability of observing reduction(105, 10) from rpois(lam=90, r=10) is larger since multiple values of x map to the same value of y
#' y <- drpois(x = reduction(105, 10), lambda = 90, red = 10)
#'
#' Y <- drpois(seq(1,20,1), 55, 10)
#' plot(Y, xlab="Y=rPois(lambda=55, r=10)", main="X~Poisson(lambda=55)", ylab="P[Y=y]")
#' @export
drpois <- drpois_1



red_Like_closed_ <- function(par, nit, K, red, FUN=round, VERBOSE=FALSE, PARALLELIZE=FALSE) {
  T <- length(unique(nit$time))
  R <- length(unique(nit$site))

  pdet <- plogis(par[2])
  lamb <- exp(par[1])
  Y    <- nit
  Y_df <- nit

  # TODO: can this be made more efficient using matrix multiplications?
  l <- 0
  if(PARALLELIZE) {
   li <- foreach(i=1:R) %dopar% {
     Yi_df <- Y_df[which(Y_df$site==i),]
     li <- 0
     ni <- max(Yi_df$count)

     for(Ni in ni:K[i,1]) {
       lit <- with(data = Yi_df, expr = drbinom(x = count, size = Ni, prob = pdet, red=reduc))

       lit_ <- prod(lit)
       li <- li + lit_*drpois(x = Ni, lambda = lamb, red=red[i,1])
     }
     return(log(li))
    }
    l <- sum(unlist(li))
    ###
  } else {

    for(i in 1:R) {
      Yi_df <- Y_df[which(Y_df$site==i),]
      li <- 0
      ni <- max(Yi_df$count)

      for(Ni in ni:K[i,1]) {
        lit <- with(data = Yi_df, expr = drbinom(x = count, size = Ni, prob = pdet, red=reduc))

        lit_ <- prod(lit)
        li <- li + lit_*drpois(x = Ni, lambda = lamb, red=red[i,1])
      }
      l <- l+log(li)
    }
  }
  if(VERBOSE) {print(paste0("log likelihood: ",l))}
  return(-1*l)
}
#' Internal function. Used to calculate the negative of the log likelihood.
#' @param par Vector with two elements, log(lambda) and logis(pdet).
#' @param nit data.frame with R*T rows, and 3 columns: "site", "time", and "count". Deprecated: R by T matrix of reduced counts with R sites/rows and T sampling occassions/columns.
#' @param K   Upper bound on summations (input reduced count upper bound).
#' @param red reduction factor matrix (R by T, with R sites/rows and T sampling occassions/columns)
#' @param VERBOSE If true, prints the calculated log likelihood to console.
#' @param PARALLELIZE If true, will use as many threads as have been made available (initialize with START_PARALLEL(num_cores)).
#' @export
red_Like_closed <- compiler::cmpfun(red_Like_closed_)



red_Like_open_ <- function(par, nit, K, red, FUN=round, VERBOSE=FALSE, PARALLELIZE=FALSE) {
  T <- ncol(nit)
  R <- nrow(nit)

  pdet <- plogis(par[4])
  lamb <- exp(par[1])
  omeg <- plogis(par[3])
  gamm <- exp(par[2])

  Y <- nit

  # g1_t_star[k] holds g1[k] * g_star[k]
  # allocate memory for vectors

  # g1_t_star <- rep(0, times=K+1)
  # g1_t      <- numeric(K+1)
  # g1        <- numeric(K+1)
  # g2        <- numeric(K+1)
  # g_star    <- rep(1, times=K+1)

  # g3 is the transition probability matrix
  # allocate memory for matrix
  #g3 <- matrix(0, nrow = K+1, ncol=K+1)
  #g3 <- tp_MAT(M = g3, omeg = omeg, gamm = gamm, red = red, PARALLELIZE=PARALLELIZE)

  g1_t_star <- list()
  g1_t      <- list()
  g1        <- list()
  g2        <- list()
  g_star    <- list()

  g3 <- list()

  if(var(as.vector(red))==0) {
    tempMat <- matrix(0, nrow = K[1]+1, ncol=K[1]+1)
    g3[[1]] <- tp_MAT(M = tempMat, omeg = omeg, gamm = gamm, red = red[1,1], PARALLELIZE=PARALLELIZE)
    for(i in 1:R) {
      g3[[i]] <- g3[[1]]
    }
  } else {
    if(PARALLELIZE) {
      g3 <- foreach(i=1:R, .packages = c("redNMix","foreach")) %dopar% {
        tempMat        <- matrix(0, nrow = K[i]+1, ncol=K[i]+1)
        tp_MAT(M = tempMat, omeg = omeg, gamm = gamm, red = red[i,1], PARALLELIZE=TRUE)
      }
      # for(i in 1:R) {
      #   tempMat        <- matrix(0, nrow = K[i]+1, ncol=K[i]+1)
      #   g3[[i]]        <- tp_MAT(M = tempMat, omeg = omeg, gamm = gamm, red = red[i,1], PARALLELIZE=PARALLELIZE)
      # }
    } else {
      for(i in 1:R) {
        tempMat        <- matrix(0, nrow = K[i]+1, ncol=K[i]+1)
        g3[[i]]        <- tp_MAT(M = tempMat, omeg = omeg, gamm = gamm, red = red[i,1], PARALLELIZE=PARALLELIZE)
      }
    }

  }

  for(i in 1:R) {
    g1_t_star[[i]] <- rep(0, times=K[i]+1)
    g1_t[[i]]      <- numeric(K[i]+1)
    g1[[i]]        <- numeric(K[i]+1)
    g2[[i]]        <- numeric(K[i]+1)
    g_star[[i]]    <- rep(1, times=K[i]+1)
  }

  # apply over sites (1 to R), this is a prime candidate for parallel processing! since each site i is independent
  ll_i  <- vapply(X = 1:R, FUN = function(i, K, T, Y, lamb, pdet, red, g3, g1_t_star, g1_t,g1,g2, g_star){
    g3        <- g3[[i]]
    g1_t_star <- g1_t_star[[i]]
    g1_t      <- g1_t[[i]]
    g1        <- g1[[i]]
    g2        <- g2[[i]]
    g_star    <- g_star[[i]]
    K         <- K[i]

    # loop backwards over times t, stopping at t==2
    for(t in T:2) {
      # size takes possible value of N (0 to K) at time t (for t > 1)
      g1_t <- drbinom(x = Y[i,t], size = (0:K), prob = pdet, red = red[i,t])
      g1_t_star <- g1_t * g_star

      # update g_star
      g_star = g3 %*% g1_t_star
    }

    # size takes possible values of N (0 to K) at time t==1
    g1 <- drbinom(x = Y[i,1], size = 0:K, prob = pdet, red = red[i,1])
    g2 <- drpois(x = 0:K, lambda = lamb, red = red[i,1])

    # apply recursive definition of likelihood
    return( log(sum(g1 * g2 * g_star)) ) # + 1e-320
  }, FUN.VALUE = numeric(1), K=K, T=T, Y=Y, lamb=lamb, pdet=pdet, red=red, g3=g3, g1_t_star=g1_t_star, g1_t=g1_t,g1=g1,g2=g2,g_star=g_star)

  ll <- sum(unlist(ll_i))
  ###

  if(VERBOSE) { print(paste0("log likelihood: ",ll)) }
  return(-1*ll)
}
#' Internal function. Used to calculate the negative of the log likelihood.
#' @param par     Vector with four elements, log(lambda), log(gamma), logis(omega), and logis(pdet).
#' @param nit     R by T matrix of reduced counts with R sites/rows and T sampling occassions/columns.
#' @param K       Upper bound on summations (reduced counts upper bound).
#' @param red     Reduction factor
#' @param VERBOSE If true, prints the log likelihood to console.
#' @param PARALLELIZE If true, will use as many threads as have been made available (initialize with START_PARALLEL(num_cores)).
#' @details Note that this function is adapted from the negative log likelihood function from the unmarked package, and uses the recursive method of computation described in Web Appendix A of Dail and Madsen 2011: Models for Estimating Abundance from Repeated Counts of an Open Metapopulation, published in Biometrics volume 67, issue 2.
#' @export
red_Like_open <- compiler::cmpfun(red_Like_open_)

#' Not Vectorized.
tp_jk <- function(j,k,omeg,gamm,red) {
  p  <- 0
  rb <- drbinom(x = 0:min(j,k), size = j, prob = omeg, red = red)
  rp <- drpois(x = k-0:min(j,k), lambda = gamm, red = red)

  p <- sum(rb * rp)

  return(p)
}

tp_jk_VV <- Vectorize(tp_jk, vectorize.args = c("j", "k"))

#' Internal function, calculates transition probabilities from pop size j to pop size k in the open population likelihood.
tp_jk_V_ <- function(j_vec,k_vec,omeg,gamm,red) {

  p <- mapply(FUN = function(j,k,omeg,gamm,red) {
    rb <- drbinom(x = 0:min(j,k), size = j, prob = omeg, red = red)
    rp <- drpois_1(x = k-0:min(j,k), lambda = gamm, red = red)

    return(sum(rb * rp))
  }, j=j_vec, k=k_vec, omeg=omeg, gamm=gamm,red=red)


  return(p)
}
tp_jk_V <- compiler::cmpfun(tp_jk_V_)


#' Internal function, calculates transition probability matrix (transition from row pop to column pop)
tp_MAT_ <- function(M, omeg, gamm, red, PARALLELIZE=FALSE) {
  K1 <- 1:(nrow(M))
  if(PARALLELIZE) {
    M <- foreach(a = K1-1, .combine = rbind) %dopar% {
        tp_jk_V(j_vec = a, k_vec = 1:nrow(M)-1, omeg = omeg, gamm = gamm, red = red)
    }
  } else {
    M <- outer(X = K1-1,Y = K1-1, FUN = tp_jk_V, omeg, gamm, red)
  }
  return(M)
}
tp_MAT <- compiler::cmpfun(tp_MAT_)


#' Find maximum likelihood estimates for model parameters log(lambda) and logit(pdet). Uses optim.
#' @param starts Vector of starting values for optimize. Has two elements, log(lambda) and logit(pdet).
#' @param nit    R by T matrix of full counts with R sites/rows and T sampling occassions/columns.
#' @param K      Upper bound on summations, either a single number, or a vector of K values, one for each site (full count value, eg if K=300 for full counts, K=reduction(300,red) for reduced counts).
#' @param red    reduction factor, either a number, or a vector of reduction factors (R sites reductions).
#' @param VERBOSE If true, prints the log likelihood to console at each optim iteration.
#' @param ...    Additional input for optim.
#' @examples
#' START_PARALLEL(num_cores=4)
#' Y    <- gen_Nmix_closed(8,8,250,0.5)
#' out  <- fit_red_Nmix_closed(Y$nit, red=10, K=300, starts = c(log(250),boot::logit(0.5)), PARALLELIZE=TRUE)
#' out2 <- fit_red_Nmix_closed(Y$nit, red=c(10,10,10,10,20,20,40,40), K=300, starts = c(log(250),boot::logit(0.5)), PARALLELIZE=TRUE)
#' END_PARALLEL()
#' @export
fit_red_Nmix_closed <- function(nit, red, K, starts=c(1,0), VERBOSE=FALSE, PARALLELIZE=FALSE, method="BFGS", ...) {

  Y_m <- nit
  row.names(Y_m) <- 1:nrow(nit)
  colnames(Y_m) <- 1:ncol(nit)

  Y_df <- reshape2::melt(Y_m)
  colnames(Y_df) <- c("site", "time", "count")


  if(length(red)==1) {
    red <- rep(x = red, times = nrow(Y_m))
  }

  if(length(red)!=nrow(Y_m)) { stop("reduction must be either one number, or a vector with length equal to nrow(nit)") }

  red <- matrix(red, nrow=nrow(Y_m), ncol=ncol(Y_m))

  redu <- numeric(nrow(Y_df))
  temp <- numeric(nrow(Y_df))
  for(i in 1:nrow(Y_df)) {
    redu[i] <- red[Y_df$site[i],Y_df$time[i]]
    temp[i] <- reduction(x = Y_df$count[i], red = redu[i])
  }
  Y_df$count <- temp
  Y_df$reduc <- redu

  red_K  <- reduction(x = matrix(K, nrow=nrow(red), ncol=ncol(red)), red = red)

  opt <- optim(par      = starts,
                fn      = red_Like_closed,
                nit     = Y_df,
                K       = red_K,
                red     = red,
                VERBOSE = VERBOSE,
                method  = method,
                PARALLELIZE = PARALLELIZE,
                ...)
  return(opt)
}

#' Used to initialize parallel computing (useful for likelihood calculations which can be computationally intensive).
#' @param num_cores Number of cores to use for parallel processing.
#' @export
START_PARALLEL <- function(num_cores) {
  require(doParallel)
  require(foreach)
  registerDoParallel(cores=num_cores)
}

#' Used to end parallel computing.
#' @export
END_PARALLEL <- function() {
  require(doParallel)
  stopImplicitCluster()
}

#' Find maximum likelihood estimates for model parameters log(lambda), log(gamma), logit(omega), and logit(pdet). Uses optim.
#' @param starts Vector with four elements, log(lambda), log(gamma), logit(omega), and logit(pdet).
#' @param nit    R by T matrix of full counts with R sites/rows and T sampling occassions/columns.
#' @param K      Upper bound on summations.
#' @param red    reduction factor, either a number or a vector of length R.
#' @param VERBOSE If true, prints the log likelihood to console at each optim iteration.
#' @param ...    Additional input for optim.
#' @examples
#' START_PARALLEL(num_cores=4)
#' Y   <- gen_Nmix_open(num_sites = 5, num_times = 5, lambda = 20, pdet = 0.7, omega = 0.7, gamma = 2)
#' out <- fit_red_Nmix_open(nit = Y$nit, red = c(1), K = 40, starts = c(0.5, 0.5, 0.5, 0.5), PARALLELIZE=TRUE)
#' END_PARALLEL()
#' @export
fit_red_Nmix_open <- function(nit, red, K, starts=c(1,1,0,0), VERBOSE=FALSE, PARALLELIZE=FALSE, method="BFGS", ...) {
   if(length(red)==1) {
     red <- rep(red, times=nrow(nit))
   }
   red    <- matrix(red, nrow=nrow(nit), ncol=ncol(nit))
   red_K  <- reduction(x = matrix(K, nrow=nrow(red), ncol=ncol(red)), red = red)

   opt <- optim(par      = starts,
                fn      = red_Like_open,
                nit     = reduction(x = nit, red = red),
                K       = red_K,
                red     = red,
                VERBOSE = VERBOSE,
                method  = method,
                PARALLELIZE = PARALLELIZE,
                ...)
   return(opt)
}

#' Plot likelihood given a pdet and range for lambda.
#' @param nit         R by T matrix of reduced counts with R sites/rows and T sampling occassions/columns.
#' @param startLambda Starting value for lambda.
#' @param endLambda   Ending value for lambda.
#' @param stepsize    Spacing between values of lambda
#' @param pdet        Probability of detection (pdet = 0.5 means 50\% chance of detection)
#' @param red         Reduction factor.
#' @param K           Upper bound on summations (reduced counts upper bound).
#' @examples
#' Y <- gen_Nmix_closed(5,5,250,0.5)
#' plot_red_like_closed_lambda(Y=reduction(Y$nit,10), startLambda = 150, endLambda = 350, stepsize=10, pdet = 0.5, red = 10, K = reduction(400,10))
#' @export
plot_red_like_closed_lambda <- function(nit, startLambda, endLambda, stepsize, pdet, red, K) {

  parB <- startLambda
  parT <- endLambda
  j <- 1
  L <- NULL
  for(parM in seq(parB,parT,stepsize)) {
    L[j] <- -1*red_Like_closed(nit = nit, par = c(log(parM), boot::logit(pdet)), K = K, red = red)
    j <- j+1
  }

  plot(y=L, x=seq(parB,parT,stepsize), xlab="lambda")
}

#' Plot likelihood given a lambda and range for pdet.
#' @param nit         R by T matrix of reduced counts with R sites/rows and T sampling occassions/columns.
#' @param startPdet   Starting value for pdet.
#' @param endPdet     Ending value for pdet.
#' @param stepsize    Spacing between values of pdet
#' @param lambda      Initial abundance parameter.
#' @param K   Upper bound on summations (reduced counts upper bound).
#' @param red reduction factor.
#' @examples
#' Y <- gen_Nmix_closed(5,5,250,0.5)
#' plot_red_like_closed_pdet(Y=reduction(Y$nit,10), startPdet = 0.1, endPdet = 1.0, stepsize=0.1, lambda = 250, red = 10, K = reduction(400,10))
#' @export
plot_red_like_closed_pdet <- function(nit, startPdet, endPdet, stepsize, lambda, red, K) {

  parB <- startPdet
  parT <- endPdet
  j <- 1
  L <- NULL
  for(parM in seq(parB,parT,stepsize)) {
    L[j] <- -1*red_Like_closed(nit = nit, par = c(log(lambda), boot::logit(parM)), K = K, red = red)
    j <- j+1
  }

  plot(y=L, x=seq(parB,parT,stepsize), xlab="pdet")
}

#' Plot 2D likelihood given a range for pdet and for lambda.
#' @param nit             R by T matrix of reduced counts with R sites/rows and T sampling occassions/columns.
#' @param startPdet       Starting value for pdet.
#' @param endPdet         Ending value for pdet.
#' @param stepsizePdet    Spacing between values of pdet
#' @param startLambda     Starting value for lambda.
#' @param endLambda       Ending value for lambda.
#' @param stepsizeLambda  Spacing between values of lambda
#' @param K               Upper bound on summations (reduced counts upper bound).
#' @param red             Reduction factor.
#' @return                Returns a matrix of likelihoods where rows represent pdet, and columns represent lambda.
#' @examples
#' Y <- gen_Nmix_closed(5,5,250,0.5)
#' plot_2d_red_like_closed(reduction(Y$nit,10), 0.1, 1.0, 0.25, 100, 400, 100, 10, reduction(400,10))
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
      L[j,k] <- -1*red_Like_closed(nit = nit, par = c(log(par2M),boot::logit(par1M)), K = K, red = red)
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
        main = "Log Likelihood Contour")
  box()
  image(x = prange, y=lrange, exp(L),
        ylab = "lambda", xlab = "pdet",
        main = "Likelihood Contour")
  box()

  par(mar = c(1,1,1,1))
  persp(prange, lrange, L, theta = 30, phi = 30, col = "lightblue", shade = 0.25, ticktype = "detailed", xlab = "pdet", ylab="lambda", zlab = "Log Likelihood")
  persp(prange, lrange, L, theta = 120, phi = 30, col = "lightblue", shade = 0.25, ticktype = "detailed", xlab = "pdet", ylab="lambda", zlab = "Log Likelihood")
  persp(prange, lrange, L, theta = 210, phi = 30, col = "lightblue", shade = 0.25, ticktype = "detailed", xlab = "pdet", ylab="lambda", zlab = "Log Likelihood")
  persp(prange, lrange, L, theta = 300, phi = 30, col = "lightblue", shade = 0.25, ticktype = "detailed", xlab = "pdet", ylab="lambda", zlab = "Log Likelihood")

  return(L2)
}

