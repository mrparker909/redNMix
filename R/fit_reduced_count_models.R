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
#' @description Tally the remainders upon division for \eqn{r \in [rmin,rmax]}, return a data.frame df
#' sorted by best r (smallest sum of remainders after division by r). df has three
#' columns, red, rem, and ss. red is the reduction factor and rem is the sum of the remainders upon
#' reduction by red. ss is the absolute change in standard deviation of the reduced counts compared
#' with the full counts (abs(r*ss_red-ss_full)).
#' abs(full counts std dev - r* reduced counts std dev).
#' @param rmin  minimum reduction factor to calculate (default is 1)
#' @param rmax  maximum reduction factor to calculate (default is 25)
#' @param rstep step size between r values (default is 1)
#' @param nit   count data
#' @examples
#' data <- c(215,309,116,157,213,248,112)
#' # Suppose you want the reduction factor which loses the least information, but which is at least 10:
#' lurf(data, rmin=10, plot=TRUE) # the plot shows patterns in changing r
#' # examining the data.frame for smallest SS shows that 11 and 14 will give better results than 10.
#'
#' # we can compare the reduced counts from these three options (r=10,11,14) to the original data:
#' data
#' 10*reduction(data, red = 10)
#' 11*reduction(data, red = 11)
#' 14*reduction(data, red = 14)
#' # this makes it clear why 11 and 14 are better choices of r than 10 for this data (even though they are larger reductions!)
#' @export
lurf <- function(nit, rmin=1, rmax=25, rstep=1, plot=FALSE) {
  redseq <- seq(rmin,rmax,1)
  remseq <- sapply(X = redseq, FUN = function(X, nit) {
    sum(nit%%X)
  }, nit=nit)
  SSseq <- sapply(X = redseq, FUN = function(X, nit) {
    abs(sqrt(var(X*(reduction(nit,red = X)))) - sqrt(var(nit)))
  }, nit=nit)
  df <- data.frame(red=redseq[order(remseq)], rem=remseq[order(remseq)], SS=SSseq[order(remseq)])
  if(plot) {
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
  Ntemp <- c(rep(rpois(n=num_sites,lambda = lambda),times=num_times))
  Ni    <- matrix(data=Ntemp, nrow = num_sites, ncol = num_times)

  nit   <- Ni
  nit[] <- vapply(Ni, function(x) { rbinom(size = x, n = 1, prob = pdet) }, numeric(1))

  return(list(Ni=Ni, nit=nit))
}

# rounding in R rounds to nearest even number (eg 0.5 rounds to 0), this function does
# "normal" rounding (eg 0.5 rounds to 1)
round2 <- function(x, n=0) {
  z <- sign(x)*trunc(abs(x)*10^n + 0.5)/10^n
  return(z)
}

#' Generate a population/observation pair with the structure of an open N-mixture model
#'
#' @param num_sites The number of observation sites.
#' @param num_times The number of sampling occasions.
#' @param lambda    The population rate parameter (\eqn{N_i \sim Poisson(\lambda)}).
#' @param pdet      The probability of detection \eqn{p} (\eqn{n_{it} \sim Binomial(N_i,p)}).
#' @param omega     The probability of survival (\eqn{S_{it} \sim Binomial(N_{it}, \omega)}). Either a single number, or a vector of length t-1 (omega for each sampling occassion after the first).
#' @param gamma     The recruitment rate parameter (\eqn{G_{it} \sim Poisson(\gamma)}). Either a single number, or a vector of length t-1 (gamma for each sampling occassion after the first).
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
  gamm = if(length(gamma)==1) {gamm<-rep(gamma, times=T-1) } else { gamm<-gamma }
  omeg = if(length(omega)==1) {omeg<-rep(omega, times=T-1) } else { gamm<-gamma }

  Ntemp <- c(rep(rpois(n=U,lambda = lamb),times=T))
  Ni <- matrix(data=Ntemp, nrow = U, ncol = T)
  for(i in 2:T) {
    Ni[,i] <- rbinom(n = U, size = Ni[,i-1], prob = omeg[i-1]) + rpois(n = U, lambda = gamm[i-1])
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

#' Reduced binomial probability distribution function \eqn{rBinomial(x;N,p,R(x;r))},
#' takes reduced quantiles rather than full quantiles (use drbinom2 for full quantiles).
#'
#' @param x Reduced count quantile (alternatively input reduction(x,r) if x is a full count quantile).
#' @param size Number of trials (full count size).
#' @param prob Probability of success for each trial.
#' @param red The factor r by which x has been reduced.
#' @return The probability of observing quantile \eqn{x}.
#' @examples
#' Y <- drbinom(0:20, 20, 0.3, 10)
#' plot(Y, xlab="Y~rBinom(N=20,p=0.3,r=10)", ylab="P[Y=y]")
#' @export
drbinom <- function(x, size, prob, red, log=FALSE) {
  # start <- (x-0.5)*red #- red/2
  # end   <- floor(start + red)
  start <- ceiling(x*red - red/2)-1
  end   <- round2(x*red + red/2)-1

  p <- pbinom(end, size, prob) - pbinom(start, size, prob)
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
#' @param precBits Number of bits of precision for arbitrary precision arithmetic.
#' @return The probability of observing quantile \eqn{x}.
#' @examples
#' Y <- drbinomAPA(0:3, 20, 0.3, 10, precBits=64)
#' @export
drbinomAPA <- function(x, size, prob, red, precBits=128, log=FALSE) {
  if(class(prob)!="mpfr") {
    prob <- Rmpfr::mpfr(prob, precBits)
  }
  if(length(x)!=length(size)) {
    if(length(x) < length(size)) {
      x <- rep(x,length(size))
    } else {
      size <- rep(size, length(x))
    }
  }

  size0 <- which(size==0)

  pt <- Rmpfr::mpfrArray(x = 0, precBits = precBits, dim = c(length(x),1))

  i <- 0
  for(X in x) {
    i <- i+1
    # if(class(X)!="mpfr") {
    #   X <- Rmpfr::mpfr(X, precBits)
    # }
    start <- Rmpfr::pmax(0,ceiling(X*red - red/2)-1)
    end   <- Rmpfr::pmin(size[i],round2(X*red + red/2)-1)

    if(start==0 & end < red) {
      pt[i] <- pbinom_APA(x = end, n = size[i], p = prob, precBits=precBits)
    } else {
      pt[i] <- pbinom_APA(x = end, n = size[i], p = prob, precBits=precBits) - pbinom_APA(x = start, n = size[i], p = prob, precBits = precBits)
    }

  }
  pt[size0] <- 0
  if(log) {
    return(log(pt))
  }
  return(pt)
}

pbinom_ <- function(x, n, p, log.p=FALSE) {
  result <- array(0, dim=c(1, length(n)))
  i <- 0
  for(N in n) {
    i <- i+1
    if(x > N) result[i] <- 1
    if(x < 0) result[i] <- 0
    result[i] <- pbeta(q = p, shape1 = x + 1, shape2 = N - x, log.p = log.p, lower.tail = FALSE)
  }

  return(result)
}

pbinom_APA <- function(x, n, p, log.p=FALSE, precBits=128) {
  result <- Rmpfr::mpfrArray(0, precBits = precBits, dim=c(1, length(n)))
  i <- 0
  for(N in n) {
    i <- i+1
    if(x > N)      result[i] <- Rmpfr::mpfr(1, precBits)
    else if(x < 0) result[i] <- Rmpfr::mpfr(0, precBits) else {
      result[i] <- Rmpfr::pbetaI(q = p, shape1 = x + 1, shape2 = N - x, log.p = log.p, lower.tail = FALSE, precBits = precBits)
    }
  }

  return(result)
}


#' Reduced poisson probability distribution function \eqn{rPoisson(x;\lambda,r)}, takes reduced quantiles (use drpois2 for full quantiles).
#'
#' @param x Reduced count quantile (alternatively input reduction(x,r) if x is a full count quantile).
#' @param lambda Mean of the full count poisson distribution.
#' @param red The factor r by which x was reduced.
#' @param precBits Number of bits of precision for arbitrary precision arithmetic.
#' @return The probability of observing quantile \eqn{x}
#' @examples
#' # probability of observing 105 from pois(lam=90)
#' x <- drpoisAPA(x = 105, lambda = 90, precBits=64)
#' @export
drpoisAPA  <- function(x, lambda, red, precBits=128, log=FALSE) {
  # pt <- NULL
  # i <- 0
  # for(X in x) {
  #   i <- i + 1
  #   start <- as.integer(X*red - red/2) +1
  #   end   <- (start + red)-1
  #   temp  <- optimizeAPA::dpois_APA((start):(end), lambda, precBits=precBits)
  #   pt[i] <- sum(new("mpfr", unlist(temp)))
  # }
  pt <- NULL
  i <- 0
  for(X in x) {
    i <- i + 1
    start <- ceiling(X*red - red/2)
    end   <- round2(X*red + red/2)-1

    ptj <- NULL
    j <- 0
    for(Y in (start):(end)) {
      j <- j+1
      ptj[[j]] <- optimizeAPA::dpois_APA(Y, lambda, precBits = precBits)
    }
    pt[i] <- sum(methods::new("mpfr", unlist(ptj)))
  }
  p <- methods::new("mpfr", unlist(pt))
  if(log) {
    return(log(p))
  }
  return(p)
}


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
#' Y <- drpois(seq(0,20,1), 55, 10)
#' plot(Y, xlab="Y=rPois(lambda=55, r=10)", main="X~Poisson(lambda=55)", ylab="P[Y=y]")
#' @export
drpois  <- function(x, lambda, red, log=FALSE) {
  start <- ceiling(x*red - red/2)-1
  end   <- round2(x*red + red/2)-1

  p <- ppois(start, lambda, lower.tail = FALSE) - ppois(end, lambda, lower.tail = FALSE) #ppois(end, lambda) - ppois(start, lambda) #sum(dpois(start:end, lambda)) #

  if(log) { return(log(p)) }
  return(p)
}



#' Used to calculate the negative of the log likelihood for closed population models.
#' @param par Vector with two elements, log(lambda) and logis(pdet). If l_s_c is not NULL, need length(par) = length(l_s_c) + 2, for the B0...BK coefficients of lambda (B0 is the constant term coefficient).
#' @param nit R by T matrix of reduced counts with R sites/rows and T sampling occassions/columns.
#' @param l_s_c list of lambda site covariates (list of vectors of length R (number of sites))
#' @param p_s_c list of pdet site covariates (list of vectors of length R (number of sites))
#' @param K   Upper bound on summations (input reduced count upper bound).
#' @param red reduction factor matrix (R by T, with R sites/rows and T sampling occassions/columns)
#' @param VERBOSE If true, prints the calculated log likelihood to console.
#' @param PARALLELIZE If true, calculation will be split over threads by sites. Will use as many threads as have been made available (initialize with START_PARALLEL(num_cores)).
#' @param APA Default is FALSE. If true, will use arbitrary precision arithmetic in calculating the likelihood. Note that APA will be slower, however it is required for site population sizes larger than about 200. If APA = TRUE, then precBits specifies the number of bits of precision to use in the calculations.
#' @param precBits If APA=TRUE, then precBits specifies the number of bits of precision for arbitrary precision arithmetic.
#' @export
red_Like_closed <- function(par, nit, l_s_c, p_s_c, K, red, FUN=round, VERBOSE=FALSE, PARALLELIZE=FALSE, APA=FALSE, precBits=64) {
  T <- ncol(nit)
  R <- nrow(nit)

  if(!APA) { # NOT APA
    # extract lambda estimates from par, setup lambda covariate matrix lamb, and covariate vector B_l
    lamb <- matrix(rep(1,times=R),ncol=1) # covariate coefficients for lambda
    B_l <- par[1] # covariates for lambda
    if(!is.null(l_s_c)) {
      lamb <- cbind(lamb, do.call(cbind, l_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix
      B_l  <- sapply(X = 1:(length(l_s_c)+1), FUN = function(X,par) { # one coeff per covariate +1 for baseline B0
        par[X]
      }, par=par)
    }

    # extract pdet estimates from par, setup pdet covariate matrix pdet, and covariate vector B_p
    pdet <- matrix(rep(1,times=R),ncol=1) # covariate coefficients for lambda
    B_p <- par[length(B_l)+1] # covariates for pdet
    if(!is.null(p_s_c)) {
      pdet <- cbind(pdet, do.call(cbind, p_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix
      B_p  <- sapply(X = 1:(length(p_s_c)+1)+length(B_l), FUN = function(X,par) { # one coeff per covariate +1 for baseline B0
        par[X]
      }, par=par)
    }

    Y <- nit

    l <- 0
    if(PARALLELIZE) { # PARALLEL and NOT APA
      li <- foreach(i=1:R) %dopar% {
        li <- 0
        ni <- max(Y[i,])
        for(Ni in ni:K[i,1]) {
          lit <- 1
          pt <- plogis(sum(B_p*pdet[i,]))
          for(t in 1:T) {
            lit <- lit*drbinom(x = Y[i,t], size = Ni*red[i,1], prob = pt, red=red[i,1])
          }
          li <- li + lit*drpois(x = Ni, lambda = exp(sum(B_l*lamb[i,])), red=red[i,1])
        }
        return(log(li))
      }
      l <- sum(unlist(li))
      ###
    } else { # NOT PARALLEL and NOT APA:

      for(i in 1:R) {
        li <- 0
        ni <- max(Y[i,])
        for(Ni in ni:K[i,1]) {
          lit <- 1
          pt <- plogis(sum(B_p*pdet[i,]))
          for(t in 1:T) {
            lit <- lit*drbinom(x = Y[i,t], size = Ni*red[i,1], prob = pt, red=red[i,1])
            #print(paste("lit:", as.numeric(lit)))
          }
          li <- li + lit*drpois(x = Ni, lambda = exp(sum(B_l*lamb[i,])), red=red[i,1])
          #print(paste("li: ",as.numeric(li)))
        }
        l <- l+log(li)
      }
    }
  } else { # DO APA CALCULATIONS
    if(!require(Rmpfr)) { stop("Error, cannot load package Rmpfr") }

    # extract lambda estimates from par, setup lambda covariate matrix lamb, and covariate vector B_l
    lamb <- Rmpfr::mpfrArray(matrix(rep(1,times=R), ncol=1), dim = c(R,1), precBits=precBits) # covariate coefficients for lambda
    B_l <- Rmpfr::mpfr(par[1], precBits=precBits) # covariates for lambda
    if(!is.null(l_s_c)) {
      lamb <- cbind(lamb, do.call(cbind, l_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix
      B_l  <- sapply(X = 1:(length(l_s_c)+1), FUN = function(X,par) { # one coeff per covariate +1 for baseline B0
        Rmpfr::mpfr(par[X], precBits=precBits)
      }, par=par)
    }

    # extract pdet estimates from par, setup pdet covariate matrix pdet, and covariate vector B_p
    pdet <- Rmpfr::mpfrArray(matrix(rep(1,times=R), ncol=1), dim = c(R,1), precBits=precBits) # covariate coefficients for lambda
    B_p <- Rmpfr::mpfr(par[length(B_l)+1], precBits=precBits) # covariates for pdet
    if(!is.null(p_s_c)) {
      pdet <- cbind(pdet, do.call(cbind, p_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix
      B_p  <- sapply(X = 1:(length(p_s_c)+1)+length(B_l), FUN = function(X,par) { # one coeff per covariate +1 for baseline B0
        Rmpfr::mpfr(par[X], precBits=precBits)
      }, par=par)
    }

    Y <- nit

    l <- Rmpfr::mpfr(0, precBits=precBits)
    if(PARALLELIZE) { # PARALLEL and APA
      li <- foreach(i=1:R) %dopar% {
        li <- Rmpfr::mpfr(0, precBits=precBits)
        ni <- max(Y[i,])
        for(Ni in ni:K[i,1]) {
          lit <- Rmpfr::mpfr(1,precBits=precBits)
          pt <- optimizeAPA::plogis_APA(sum(B_p*pdet[i,]), precBits = precBits)
          for(t in 1:T) {
            lit <- lit*(drbinomAPA(x    = Y[i,t],
                                   size = Ni,
                                   prob = pt,
                                   red  = red[i,1],
                                   precBits = precBits))
          }
          li <- li + lit*drpoisAPA(x = Ni, lambda = exp(sum(B_l*lamb[i,])), red=red[i,1], precBits = precBits)

        }
        return(log(li))
      }

      l <- sum(methods::new("mpfr", unlist(li)))
      ###
    } else { # NOT PARALLEL and APA
      for(i in 1:R) {
        li <- Rmpfr::mpfr(0, precBits=precBits)
        ni <- max(Y[i,])
        for(Ni in ni:K[i,1]) {
          lit <- Rmpfr::mpfr(1,precBits=precBits)
          pit <- optimizeAPA::plogis_APA(sum(B_p*pdet[i,]), precBits = precBits)
          for(t in 1:T) {
            lit <- lit*(drbinomAPA(x    = Y[i,t],
                                   size = Ni,
                                   prob = pit,
                                   red  = red[i,1],
                                   precBits = precBits))
          }
          li <- li + lit*drpoisAPA(x = Ni, lambda = exp(sum(B_l*lamb[i,])), red=red[i,1], precBits = precBits)
        }
        l <- l+log(li)
      }
    }
  }

  if(VERBOSE) {
    if(!APA) {
      print(paste0("log likelihood: ", l))
      print(paste0("parameters: ", par))
    }
    if(APA) {
      print(paste0("log likelihood: ", Rmpfr::formatMpfr(l)))
      print(paste0("parameters: ", Rmpfr::formatMpfr(par)))
    }
  }
  return(-1*l)
}


#' Used to calculate the negative of the log likelihood for open population models.
#' @param par     Vector with four elements, log(lambda), log(gamma), logis(omega), and logis(pdet).
#' @param nit     R by T matrix of reduced counts with R sites/rows and T sampling occassions/columns.
#' @param l_s_c   list of lambda site covariates (list of vectors of length R (number of sites)).
#' @param g_s_c   list of gamma site covariates (list of vectors of length R (number of sites))
#' @param g_t_c   list of gamma time covariates (list of vectors of length T (number of sampling occasions))
#' @param o_s_c   list of omega site covariates (list of vectors of length R (number of sites))
#' @param o_t_c   list of omega time covariates (list of vectors of length T (number of sampling occasions))
#' @param p_s_c   list of pdet site covariates (list of vectors of length R (number of sites))
#' @param p_t_c   list of pdet time covariates (list of vectors of length T (number of sampling occasions))
#' @param K       Upper bound on summations (reduced counts upper bound).
#' @param red     Reduction factor
#' @param VERBOSE If true, prints the log likelihood to console.
#' @param PARALLELIZE If true, calculation will be split over threads by sites. Will use as many threads as have been made available (initialize with START_PARALLEL(num_cores)).
#' @param APA Default is FALSE. If TRUE, will use arbitrary precision arithmetic in calculating the likelihood. Note that APA will be slower, however it is required for site population sizes larger than about 200. If APA = TRUE, then precBits specifies the number of bits of precision to use in the calculations.
#' @param precBits If APA=TRUE, then precBits specifies the number of bits of precision for arbitrary precision arithmetic.
#' @details Note that this function is adapted from the negative log likelihood function from the unmarked package, and uses the recursive method of computation described in Web Appendix A of Dail and Madsen 2011: Models for Estimating Abundance from Repeated Counts of an Open Metapopulation, published in Biometrics volume 67, issue 2.
#' @export
red_Like_open <- function(par, nit, l_s_c, g_s_c, g_t_c, o_s_c, o_t_c, p_s_c, p_t_c, K, red, VERBOSE=FALSE, PARALLELIZE=FALSE, APA=FALSE, precBits=128) {
  T <- ncol(nit)
  R <- nrow(nit)

  if(!APA) { # NOT APA
    # extract lambda estimates from par, setup lambda covariate matrix lamb, and covariate vector B_l
    lamb <- matrix(rep(1,times=R),ncol=1) # covariate coefficients for lambda
    B_l <- par[1] # covariates for lambda
    if(!is.null(l_s_c)) {
      lamb <- cbind(lamb, do.call(cbind, l_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix

      B_l <- sapply(X = 1:(length(l_s_c)+1), FUN = function(X,par) { # one coeff per covariate +1 for baseline B0
        par[X]
      }, par=par)
    }

    # extract gamma estimates from par, setup gamma covariate matrix gamm, and covariate vector B_g
    gamm <- matrix(rep(1,times=R),ncol=1) # site covariate coefficients for gamma
    B_g <- par[length(B_l)+1] # covariates for gamma
    if(!is.null(g_s_c)) {
      gamm <- cbind(gamm, do.call(cbind, g_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix

      B_g <- sapply(X = 1:(length(g_s_c)+1), FUN = function(X,par) { # one coeff per covariate +1 for baseline B0
        par[length(B_l)+X]
      }, par=par)
    }

    # extract gamma estimates from par, setup gamma covariate matrix gamt, and covariate vector B_gt
    gamt <- matrix(rep(0,times=T),ncol=1) # time covariate coefficients for gamma
    B_gt <- NULL # covariates for gamma
    if(!is.null(g_t_c)) {
      gamt <- do.call(cbind, g_t_c) # rows are times, cols are covariate values, here we are creating the design matrix
      #gamt <- rbind(gamt,rep(0, times=length(g_t_c)))

      B_gt <- sapply(X = 1:length(g_t_c), FUN = function(X,par) { # one coeff per covariate
        par[length(B_l)+length(B_g)+X]
      }, par=par)
    }

    # extract omega estimates from par, setup omega covariate matrix omeg, and covariate vector B_o
    omeg <- matrix(rep(1,times=R),ncol=1) # site covariate coefficients for omega
    B_o <- par[length(B_l)+length(B_g)+length(B_gt)+1] # covariates for omega
    if(!is.null(o_s_c)) {
      omeg <- cbind(omeg, do.call(cbind, o_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix

      B_o <- sapply(X = 1:(length(o_s_c)+1), FUN = function(X,par) { # one coeff per covariate +1 for baseline B0
        par[length(B_l)+length(B_g)+length(B_gt)+X]
      }, par=par)
    }

    # extract omega estimates from par, setup omega covariate matrix omet, and covariate vector B_ot
    omet <- matrix(rep(0,times=T),ncol=1) # time covariate coefficients for omega
    B_ot <- NULL # covariates for omega
    if(!is.null(o_t_c)) {
      omet <- do.call(cbind, o_t_c) # rows are times, cols are covariate values, here we are creating the design matrix

      B_ot <- sapply(X = 1:length(o_t_c), FUN = function(X,par) { # one coeff per covariate
        par[length(B_l)+length(B_g)+length(B_gt)+length(B_o)+X]
      }, par=par)
    }

    # extract pdet estimates from par, setup pdet covariate matrix pdet, and covariate vector B_p
    pdet <- matrix(rep(1,times=R),ncol=1) # site covariate coefficients for pdet
    B_p <- par[length(B_l)+length(B_g)+length(B_gt)+length(B_o)+length(B_ot)+1] # covariates for pdet
    if(!is.null(p_s_c)) {
      pdet <- cbind(pdet, do.call(cbind, p_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix

      B_p <- sapply(X = 1:(length(p_s_c)+1), FUN = function(X,par) { # one coeff per covariate +1 for baseline B0
        par[length(B_l)+length(B_g)+length(B_gt)+length(B_o)+length(B_ot)+X]
      }, par=par)
    }


    # extract pdet estimates from par, setup pdet covariate matrix pt, and covariate vector B_pt
    pt   <- matrix(rep(0,times=T),ncol=1) # time covariate coefficients for pdet
    B_pt <- NULL # covariates for pdet
    if(!is.null(p_t_c)) {
      pt   <- do.call(cbind, p_t_c) # rows are times, cols are covariate values, here we are creating the design matrix

      B_pt <- sapply(X = 1:length(p_t_c), FUN = function(X,par) { # one coeff per covariate
        par[length(B_l)+length(B_g)+length(B_gt)+length(B_o)+length(B_ot)+length(B_p)+X]
      }, par=par)
    }


    Y <- nit

    g1_t_star <- list()
    g1_t      <- list()
    g1        <- list()
    g2        <- list()
    g_star    <- list()
    g3        <- vector(length = R, mode = "list")#list()
    g3_temp   <- list()

    if(PARALLELIZE) {
      if(!is.null(g_t_c) | !is.null(o_t_c)) { # TIME DEPENDENT
        if(!is.null(g_s_c) | !is.null(o_s_c) | var(as.vector(red))!=0) { # SITE DEPENDENT
          g3_temp <- foreach(i=1:R, .packages = c("redNMix","foreach")) %:%
            foreach(t=1:T, .packages = c("redNMix","foreach")) %dopar% {
              tempMat        <- matrix(0, nrow = K[i]+1, ncol=K[i]+1)
              if(is.null(o_t_c)) tempMat <- tp_MAT(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=0, B_ot=0, gamm = gamm[i,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[i,1], PARALLELIZE=TRUE)
              else if(is.null(g_t_c)) tempMat <- tp_MAT(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm = gamm[i,], B_g=B_g, gamt=0, B_gt=0, red = red[i,1], PARALLELIZE=TRUE)
              else tempMat <- tp_MAT(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm = gamm[i,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[i,1], PARALLELIZE=TRUE)
              return(tempMat)
            }
          for(i in 1:R) {
            g3[[i]]<- g3_temp[[i]]
          }
        } else {# NOT SITE DEPENDENT
          g3_temp <- foreach(t=1:T, .packages = c("redNMix","foreach")) %dopar% {
            tempMat        <- matrix(0, nrow = K[1]+1, ncol=K[1]+1)
            if(is.null(o_t_c)) tempMat <- tp_MAT(M = tempMat, omeg = omeg[1,], B_o=B_o, omet=0, B_ot=0, gamm = gamm[1,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[1,1], PARALLELIZE=TRUE)
            else if(is.null(g_t_c)) tempMat <- tp_MAT(M = tempMat, omeg = omeg[1,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm = gamm[1,], B_g=B_g, gamt=0, B_gt=0, red = red[1,1], PARALLELIZE=TRUE)
            else tempMat <- tp_MAT(M = tempMat, omeg = omeg[1,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm = gamm[1,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[1,1], PARALLELIZE=TRUE)
            return(tempMat)
          }
          for(i in 1:R) {
            for(t in 1:T) {
              g3[[i]][[t]] <- g3_temp[[t]]
            }
          }
        }
      } else {  # NOT TIME DEPENDENT
        g3_temp <- foreach(i=1:R, .packages=c("redNMix","foreach")) %dopar% {
          tempMat        <- matrix(0, nrow = K[i]+1, ncol=K[i]+1)
          tp_MAT(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=0, B_ot=0, gamm = gamm[i,], B_g=B_g, gamt=0, B_gt=0, red = red[i,1], PARALLELIZE=TRUE)
        }
        for(i in 1:R) {
          for(t in 1:T) {
            g3[[i]][[t]] <- g3_temp[[i]]
          }
        }
      }
    } else { # NOT PARALLELIZED
      g3_temp <- list()
      if(var(as.vector(red))==0) { # all reduction factors are the same
        if(!is.null(g_t_c) | !is.null(o_t_c)) { # ANY covariates are time dependent
          if(is.null(g_s_c) & is.null(o_s_c)) { # ALL covariates are constant across sites
            # !r & t & !s
            for(t in 1:T) {
              tempMat <- matrix(0, nrow = K[1]+1, ncol=K[1]+1)

              if(is.null(g_t_c)) g3_temp <- tp_MAT(M = tempMat, omeg = omeg[1,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm=gamm[1,], B_g=B_g, gamt=0, B_gt=0, red = red[1,1], PARALLELIZE=FALSE)
              else if(is.null(o_t_c)) g3_temp <- tp_MAT(M = tempMat, omeg = omeg[1,], B_o=B_o, omet=0, B_ot=0, gamm=gamm[1,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[1,1], PARALLELIZE=FALSE)
              else g3_temp <- tp_MAT(M = tempMat, omeg = omeg[1,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm=gamm[1,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[1,1], PARALLELIZE=FALSE)

              for(i in 1:R) {
                g3[[i]][[t]] <- g3_temp
              }
            }
          } else {                   # at least ONE covariate is site dependent
            # !r & t & s
            for(i in 1:R) {
              for(t in 1:T) {
                tempMat      <- matrix(0, nrow = K[i]+1, ncol=K[i]+1)
                if(is.null(g_t_c)) g3_temp <- tp_MAT(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm=gamm[i,], B_g=B_g, gamt=0, B_gt=0, red = red[1,1], PARALLELIZE=FALSE)
                else if(is.null(o_t_c)) g3_temp <- tp_MAT(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=0, B_ot=0, gamm=gamm[i,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[1,1], PARALLELIZE=FALSE)
                else g3_temp <- tp_MAT(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm=gamm[i,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[1,1], PARALLELIZE=FALSE)

                g3[[i]][[t]] <- g3_temp
              }
            }
          }
        } else {                     # NO covariates are time dependent
          if(is.null(g_s_c) & is.null(o_s_c)) { # ALL covariates are constant across sites
            # !r & !t & !s
            tempMat <- matrix(0, nrow = K[1]+1, ncol=K[1]+1)
            g3_temp <- tp_MAT(M = tempMat, omeg = omeg[1,1], B_o=B_o, omet=0, B_ot=0, gamm = gamm[1,1], B_g=B_g, gamt=0, B_gt=0, red = red[1,1], PARALLELIZE=FALSE)
            for(i in 1:R) {
              for(t in 1:T) {
                g3[[i]][[t]] <- g3_temp
              }
            }
          } else {                     # at least ONE covariate is site dependent
            # !r & !t & s
            for(i in 1:R) {
              tempMat <- matrix(0, nrow = K[i]+1, ncol=K[i]+1)
              g3_temp <- tp_MAT(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=0, B_ot=0, gamm = gamm[i,], B_g=B_g, gamt=0, B_gt=0, red = red[1,1], PARALLELIZE=FALSE)
              for(t in 1:T) {
                g3[[i]][[t]] <- g3_temp
              }
            }
          }
        }
      } else {                     # some reduction factors different
        if(!is.null(g_t_c) | !is.null(o_t_c)) {        # ANY covariates are time dependent
          if(is.null(g_s_c) & is.null(o_s_c)) { # ALL covariates are constant across sites
            # r & t & !s
            # note this is the same as r & t & s since r requires site dependent calc of tp_MAT
            for(i in 1:R) {
              for(t in 1:T) {
                tempMat      <- matrix(0, nrow = K[i]+1, ncol=K[i]+1)
                if(is.null(g_t_c)) g3_temp <- tp_MAT(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm=gamm[i,], B_g=B_g, gamt=0, B_gt=0, red = red[i,1], PARALLELIZE=FALSE)
                else if(is.null(o_t_c)) g3_temp <- tp_MAT(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=0, B_ot=0, gamm=gamm[i,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[i,1], PARALLELIZE=FALSE)
                else g3_temp <- tp_MAT(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm=gamm[i,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[i,1], PARALLELIZE=FALSE)

                g3[[i]][[t]] <- g3_temp
              }
            }
          } else {                   # at least ONE covariate is site dependent
            # r & t & s
            for(i in 1:R) {
              for(t in 1:T) {
                tempMat      <- matrix(0, nrow = K[i]+1, ncol=K[i]+1)
                if(is.null(g_t_c)) g3_temp <- tp_MAT(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm=gamm[i,], B_g=B_g, gamt=0, B_gt=0, red = red[1,1], PARALLELIZE=FALSE)
                else if(is.null(o_t_c)) g3_temp <- tp_MAT(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=0, B_ot=0, gamm=gamm[i,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[1,1], PARALLELIZE=FALSE)
                else g3_temp <- tp_MAT(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm=gamm[i,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[1,1], PARALLELIZE=FALSE)

                g3[[i]][[t]] <- g3_temp
              }
            }
          }
        } else {                     # NO covariates are time dependent
          if(is.null(g_s_c) & is.null(o_s_c)) { # ALL covariates are constant across sites
            # r & !t & !s
            for(i in 1:R) {
              tempMat <- matrix(0, nrow = K[i]+1, ncol=K[i]+1)
              g3_temp <- tp_MAT(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=0, B_ot=0, gamm = gamm[i,], B_g=B_g, gamt=0, B_gt=0, red = red[i,], PARALLELIZE=FALSE)
              for(t in 1:T) {
                g3[[i]][[t]] <- g3_temp
              }
            }
          } else {                     # at least ONE covariate is site dependent
            # r & !t & s
            # note: same as # r & !t & !s since r requires site dependent calc of tp_MAT
            for(i in 1:R) {
              tempMat <- matrix(0, nrow = K[i]+1, ncol=K[i]+1)
              g3_temp <- tp_MAT(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=0, B_ot=0, gamm = gamm[i,], B_g=B_g, gamt=0, B_gt=0, red = red[i,], PARALLELIZE=FALSE)
              for(t in 1:T) {
                g3[[i]][[t]] <- g3_temp
              }
            }
          }
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

    # apply over sites (1 to R), TODO: this is a prime candidate for parallel since each site i is independent
    ll_i  <- vapply(X = 1:R, FUN = function(i, K, T, Y, lamb, B_l, pdet, B_p, pt, B_pt, red, g3, g1_t_star, g1_t,g1,g2, g_star) {
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
        g1_t <- drbinom(x = Y[i,t], size = (0:K)*red[i,1], prob = plogis(sum(pdet[i,] * B_p)+sum(pt[t,] * B_pt)), red = red[i,t])
        g1_t_star <- g1_t * g_star

        # update g_star
        g_star = g3[[t]] %*% g1_t_star
      }

      # size takes possible values of N (0 to K) at time t==1
      g1 <- drbinom(x = Y[i,1], size = (0:K)*red[i,1], prob = plogis(sum(pdet[i,] * B_p)+sum(pt[t,] * B_pt)), red = red[i,1])
      g2 <- drpois(x = 0:K, lambda = exp(sum(lamb[i,] * B_l)), red = red[i,1])


      # apply recursive definition of likelihood
      return( log(sum(g1 * g2 * g_star)) ) # + 1e-320
    }, FUN.VALUE = numeric(1), K=K, T=T, Y=Y, lamb=lamb, B_l=B_l, pdet=pdet, B_p=B_p, pt=pt, B_pt=B_pt, red=red, g3=g3, g1_t_star=g1_t_star, g1_t=g1_t,g1=g1,g2=g2,g_star=g_star)

    ll <- sum(unlist(ll_i))
    ###########################
    ###########################
    ###########################
  } else { # APA = TRUE
    # extract lambda estimates from par, setup lambda covariate matrix lamb, and covariate vector B_l
    lamb <- Rmpfr::mpfrArray(x = 1, precBits = precBits, dim = c(R,1)) # covariate coefficients for lambda
    B_l <- Rmpfr::mpfr(par[1], precBits=precBits) # covariates for lambda
    if(!is.null(l_s_c)) {
      lamb <- cbind(lamb, do.call(cbind, l_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix
      B_l <- sapply(X = 1:(length(l_s_c)+1), FUN = function(X,par) { # one coeff per covariate +1 for baseline B0
        Rmpfr::mpfr(par[X],precBits=precBits)
      }, par=par)
      B_l <- methods::new("mpfr", unlist(B_l))
    }

    # extract gamma estimates from par, setup gamma covariate matrix gamm, and covariate vector B_g
    gamm <- Rmpfr::mpfrArray(x = 1, precBits = precBits, dim = c(R,1)) # site covariate coefficients for gamma
    B_g <- Rmpfr::mpfr(par[length(B_l)+1], precBits=precBits) # covariates for gamma
    if(!is.null(g_s_c)) {
      gamm <- cbind(gamm, do.call(cbind, g_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix

      B_g <- sapply(X = 1:(length(g_s_c)+1), FUN = function(X,par) { # one coeff per covariate +1 for baseline B0
        Rmpfr::mpfr(par[length(B_l)+X],precBits=precBits)
      }, par=par)
      B_g <- methods::new("mpfr", unlist(B_g))
    }

    # extract gamma estimates from par, setup gamma covariate matrix gamt, and covariate vector B_gt
    gamt <- Rmpfr::mpfrArray(x = 0, precBits = precBits, dim = c(T,1)) # time covariate coefficients for gamma
    B_gt <- NULL # covariates for gamma
    if(!is.null(g_t_c)) {
      gamt <- do.call(cbind, g_t_c) # rows are times, cols are covariate values, here we are creating the design matrix
      #gamt <- rbind(gamt,rep(0, times=length(g_t_c)))

      B_gt <- sapply(X = 1:length(g_t_c), FUN = function(X,par) { # one coeff per covariate
        Rmpfr::mpfr(par[length(B_l)+length(B_g)+X],precBits=precBits)
      }, par=par)
      B_gt <- methods::new("mpfr", unlist(B_gt))
    }

    # extract omega estimates from par, setup omega covariate matrix omeg, and covariate vector B_o
    omeg <- Rmpfr::mpfrArray(x = 1, precBits = precBits, dim = c(R,1)) # site covariate coefficients for omega
    B_o  <- Rmpfr::mpfr(par[length(B_l)+length(B_g)+length(B_gt)+1], precBits=precBits) # covariates for omega
    if(!is.null(o_s_c)) {
      omeg <- cbind(omeg, do.call(cbind, o_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix

      B_o <- sapply(X = 1:(length(o_s_c)+1), FUN = function(X,par) { # one coeff per covariate +1 for baseline B0
        Rmpfr::mpfr(par[length(B_l)+length(B_g)+length(B_gt)+X], precBits=precBits)
      }, par=par)
      B_o <- methods::new("mpfr", unlist(B_o))
    }

    # extract omega estimates from par, setup omega covariate matrix omet, and covariate vector B_ot
    omet <- Rmpfr::mpfrArray(x = 0, precBits = precBits, dim = c(T,1)) # time covariate coefficients for omega
    B_ot <- NULL # covariates for omega
    if(!is.null(o_t_c)) {
      omet <- do.call(cbind, o_t_c) # rows are times, cols are covariate values, here we are creating the design matrix

      B_ot <- sapply(X = 1:length(o_t_c), FUN = function(X,par) { # one coeff per covariate
        Rmpfr::mpfr(par[length(B_l)+length(B_g)+length(B_gt)+length(B_o)+X], precBits=precBits)
      }, par=par)
      B_ot <- methods::new("mpfr", unlist(B_ot))
    }

    # extract pdet estimates from par, setup pdet covariate matrix pdet, and covariate vector B_p
    pdet <- Rmpfr::mpfrArray(x = 1, precBits = precBits, dim = c(R,1)) # site covariate coefficients for pdet
    B_p  <- Rmpfr::mpfr(par[length(B_l)+length(B_g)+length(B_gt)+length(B_o)+length(B_ot)+1], precBits=precBits) # covariates for pdet
    if(!is.null(p_s_c)) {
      pdet <- cbind(pdet, do.call(cbind, p_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix

      B_p <- sapply(X = 1:(length(p_s_c)+1), FUN = function(X,par) { # one coeff per covariate +1 for baseline B0
        Rmpfr::mpfr(par[length(B_l)+length(B_g)+length(B_gt)+length(B_o)+length(B_ot)+X], precBits=precBits)
      }, par=par)
      B_p <- methods::new("mpfr", unlist(B_p))
    }


    # extract pdet estimates from par, setup pdet covariate matrix pt, and covariate vector B_pt
    pt   <- Rmpfr::mpfrArray(x = 0, precBits = precBits, dim = c(T,1)) # time covariate coefficients for pdet
    B_pt <- NULL # covariates for pdet
    if(!is.null(p_t_c)) {
      pt   <- do.call(cbind, p_t_c) # rows are times, cols are covariate values, here we are creating the design matrix

      B_pt <- sapply(X = 1:length(p_t_c), FUN = function(X,par) { # one coeff per covariate
        Rmpfr::mpfr(par[length(B_l)+length(B_g)+length(B_gt)+length(B_o)+length(B_ot)+length(B_p)+X], precBits=precBits)
      }, par=par)
      B_pt <- methods::new("mpfr", unlist(B_pt))
    }

    Y <- nit

    g1_t_star <- list()
    g1_t      <- list()
    g1        <- list()
    g2        <- list()
    g_star    <- list()
    g3        <- vector(length = R, mode = "list")#list()
    g3_temp   <- list()

    if(PARALLELIZE) {
      if(!is.null(g_t_c) | !is.null(o_t_c)) { # TIME DEPENDENT
        if(!is.null(g_s_c) | !is.null(o_s_c) | var(as.vector(red))!=0) { # SITE DEPENDENT
          g3_temp <- foreach(i=1:R, .packages = c("redNMix","foreach","optimizeAPA","Rmpfr")) %:%
            foreach(t=1:T, .packages = c("redNMix","foreach","optimizeAPA","Rmpfr")) %dopar% {
              tempMat <- Rmpfr::mpfrArray(0, precBits = precBits, dim = c(K[i]+1,K[i]+1))

              if(is.null(o_t_c)) tempMat <- tp_MAT_APA(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=0, B_ot=0, gamm = gamm[i,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[i,1], PARALLELIZE=TRUE, precBits = precBits)
              else if(is.null(g_t_c)) tempMat <- tp_MAT_APA(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm = gamm[i,], B_g=B_g, gamt=0, B_gt=0, red = red[i,1], PARALLELIZE=TRUE, precBits = precBits)
              else tempMat <- tp_MAT_APA(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm = gamm[i,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[i,1], PARALLELIZE=TRUE, precBits = precBits)
              return(tempMat)
            }
          for(i in 1:R) {
            g3[[i]]<- g3_temp[[i]]
          }
        } else {# NOT SITE DEPENDENT
          g3_temp <- foreach(t=1:T, .packages = c("redNMix","foreach","optimizeAPA","Rmpfr")) %dopar% {
            tempMat <- Rmpfr::mpfrArray(0, precBits = precBits, dim = c(K[i]+1,K[i]+1))
            if(is.null(o_t_c)) tempMat <- tp_MAT_APA(M = tempMat, omeg = omeg[1,], B_o=B_o, omet=0, B_ot=0, gamm = gamm[1,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[1,1], PARALLELIZE=TRUE, precBits = precBits)
            else if(is.null(g_t_c)) tempMat <- tp_MAT_APA(M = tempMat, omeg = omeg[1,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm = gamm[1,], B_g=B_g, gamt=0, B_gt=0, red = red[1,1], PARALLELIZE=TRUE, precBits = precBits)
            else tempMat <- tp_MAT_APA(M = tempMat, omeg = omeg[1,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm = gamm[1,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[1,1], PARALLELIZE=TRUE, precBits = precBits)
            return(tempMat)
          }
          for(i in 1:R) {
            for(t in 1:T) {
              g3[[i]][[t]] <- g3_temp[[t]]
            }
          }
        }
      } else {  # NOT TIME DEPENDENT
        g3_temp <- foreach(i=1:R, .packages=c("redNMix","foreach","optimizeAPA","Rmpfr")) %dopar% {
          tempMat <- Rmpfr::mpfrArray(0, precBits = precBits, dim = c(K[i]+1,K[i]+1))
          tp_MAT_APA(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=0, B_ot=0, gamm = gamm[i,], B_g=B_g, gamt=0, B_gt=0, red = red[i,1], PARALLELIZE=TRUE, precBits = precBits)
        }
        for(i in 1:R) {
          for(t in 1:T) {
            g3[[i]][[t]] <- g3_temp[[i]]
          }
        }
      }
      ######################
      ######################
      ######################
    } else { # NOT PARALLELIZED
      g3_temp <- list()
      if(var(as.vector(red))==0) { # all reduction factors are the same
        if(!is.null(g_t_c) | !is.null(o_t_c)) { # ANY covariates are time dependent
          if(is.null(g_s_c) & is.null(o_s_c)) { # ALL covariates are constant across sites
            # !r & t & !s
            for(t in 1:T) {


              tempMat <- Rmpfr::mpfrArray(0, precBits = precBits, dim = c(K[1]+1,K[1]+1))

              if(is.null(g_t_c)) g3_temp <- tp_MAT_APA(M = tempMat, omeg = omeg[1,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm=gamm[1,], B_g=B_g, gamt=0, B_gt=0, red = red[1,1], PARALLELIZE=FALSE, precBits = precBits)
              else if(is.null(o_t_c)) g3_temp <- tp_MAT_APA(M = tempMat, omeg = omeg[1,], B_o=B_o, omet=0, B_ot=0, gamm=gamm[1,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[1,1], PARALLELIZE=FALSE, precBits = precBits)
              else g3_temp <- tp_MAT_APA(M = tempMat, omeg = omeg[1,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm=gamm[1,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[1,1], PARALLELIZE=FALSE, precBits = precBits)

              for(i in 1:R) {
                g3[[i]][[t]] <- g3_temp
              }
            }
          } else {                   # at least ONE covariate is site dependent
            # !r & t & s
            for(i in 1:R) {
              for(t in 1:T) {
                tempMat <- Rmpfr::mpfrArray(0, precBits = precBits, dim = c(K[i]+1,K[i]+1))
                if(is.null(g_t_c)) g3_temp <- tp_MAT_APA(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm=gamm[i,], B_g=B_g, gamt=0, B_gt=0, red = red[1,1], PARALLELIZE=FALSE, precBits = precBits)
                else if(is.null(o_t_c)) g3_temp <- tp_MAT_APA(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=0, B_ot=0, gamm=gamm[i,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[1,1], PARALLELIZE=FALSE, precBits = precBits)
                else g3_temp <- tp_MAT_APA(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm=gamm[i,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[1,1], PARALLELIZE=FALSE, precBits = precBits)

                g3[[i]][[t]] <- g3_temp
              }
            }
          }
        } else {                     # NO covariates are time dependent
          if(is.null(g_s_c) & is.null(o_s_c)) { # ALL covariates are constant across sites
            # !r & !t & !s
            tempMat <- Rmpfr::mpfrArray(0, precBits = precBits, dim = c(K[1]+1,K[1]+1))
            g3_temp <- tp_MAT_APA(M = tempMat, omeg = omeg[1,1], B_o=B_o, omet=0, B_ot=0, gamm = gamm[1,1], B_g=B_g, gamt=0, B_gt=0, red = red[1,1], PARALLELIZE=FALSE, precBits = precBits)
            for(i in 1:R) {
              for(t in 1:T) {
                g3[[i]][[t]] <- g3_temp
              }
            }
          } else {                     # at least ONE covariate is site dependent
            # !r & !t & s
            for(i in 1:R) {
              tempMat <- Rmpfr::mpfrArray(0, precBits = precBits, dim = c(K[i]+1,K[i]+1))
              g3_temp <- tp_MAT_APA(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=0, B_ot=0, gamm = gamm[i,], B_g=B_g, gamt=0, B_gt=0, red = red[1,1], PARALLELIZE=FALSE, precBits = precBits)
              for(t in 1:T) {
                g3[[i]][[t]] <- g3_temp
              }
            }
          }
        }
      } else {                     # some reduction factors different
        if(!is.null(g_t_c) | !is.null(o_t_c)) {        # ANY covariates are time dependent
          if(is.null(g_s_c) & is.null(o_s_c)) { # ALL covariates are constant across sites
            # r & t & !s
            # note this is the same as r & t & s since r requires site dependent calc of tp_MAT
            for(i in 1:R) {
              for(t in 1:T) {
                tempMat <- Rmpfr::mpfrArray(0, precBits = precBits, dim = c(K[i]+1,K[i]+1))
                if(is.null(g_t_c)) g3_temp <- tp_MAT_APA(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm=gamm[i,], B_g=B_g, gamt=0, B_gt=0, red = red[i,1], PARALLELIZE=FALSE, precBits = precBits)
                else if(is.null(o_t_c)) g3_temp <- tp_MAT_APA(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=0, B_ot=0, gamm=gamm[i,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[i,1], PARALLELIZE=FALSE, precBits = precBits)
                else g3_temp <- tp_MAT_APA(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm=gamm[i,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[i,1], PARALLELIZE=FALSE, precBits = precBits)

                g3[[i]][[t]] <- g3_temp
              }
            }
          } else {                   # at least ONE covariate is site dependent
            # r & t & s
            for(i in 1:R) {
              for(t in 1:T) {
                tempMat <- Rmpfr::mpfrArray(0, precBits = precBits, dim = c(K[i]+1,K[i]+1))
                if(is.null(g_t_c)) g3_temp <- tp_MAT_APA(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm=gamm[i,], B_g=B_g, gamt=0, B_gt=0, red = red[1,1], PARALLELIZE=FALSE, precBits = precBits)
                else if(is.null(o_t_c)) g3_temp <- tp_MAT_APA(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=0, B_ot=0, gamm=gamm[i,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[1,1], PARALLELIZE=FALSE, precBits = precBits)
                else g3_temp <- tp_MAT_APA(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm=gamm[i,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[1,1], PARALLELIZE=FALSE, precBits = precBits)

                g3[[i]][[t]] <- g3_temp
              }
            }
          }
        } else {                     # NO covariates are time dependent
          if(is.null(g_s_c) & is.null(o_s_c)) { # ALL covariates are constant across sites
            # r & !t & !s
            for(i in 1:R) {
              tempMat <- Rmpfr::mpfrArray(0, precBits = precBits, dim = c(K[i]+1,K[i]+1))
              g3_temp <- tp_MAT_APA(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=0, B_ot=0, gamm = gamm[i,], B_g=B_g, gamt=0, B_gt=0, red = red[i,], PARALLELIZE=FALSE, precBits = precBits)
              for(t in 1:T) {
                g3[[i]][[t]] <- g3_temp
              }
            }
          } else {                     # at least ONE covariate is site dependent
            # r & !t & s
            # note: same as # r & !t & !s since r requires site dependent calc of tp_MAT
            for(i in 1:R) {
              tempMat <- Rmpfr::mpfrArray(0, precBits = precBits, dim = c(K[i]+1,K[i]+1))
              g3_temp <- tp_MAT_APA(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=0, B_ot=0, gamm = gamm[i,], B_g=B_g, gamt=0, B_gt=0, red = red[i,], PARALLELIZE=FALSE, precBits = precBits)
              for(t in 1:T) {
                g3[[i]][[t]] <- g3_temp
              }
            }
          }
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

    # apply over sites (1 to R), TODO: this is a prime candidate for parallel since each site i is independent
    ll_i  <- Rmpfr::sapplyMpfr(X = 1:R, FUN = function(i, K, T, Y, lamb, B_l, pdet, B_p, pt, B_pt, red, g3, g1_t_star, g1_t,g1,g2, g_star, precBits) {
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
        p_temp <- ifelse(is.null(B_pt), 0, pt[t,] * B_pt)
        if(class(p_temp[1])=="numeric") {
          p_temp <- unlist(p_temp)
        } else {
          p_temp <- methods::new("mpfr", unlist(p_temp))
        }
        p <- optimizeAPA::plogis_APA(sum(pdet[i,] * B_p)+sum(p_temp), precBits = precBits)
        g1_t <- drbinomAPA(x = Y[i,t], size = (0:K)*red[i,1], prob = p, red = red[i,t], precBits = precBits)
        g1_t_star <- g1_t * g_star

        # update g_star
        g_star = g3[[t]] %*% g1_t_star
      }

      Bp_temp  <- Rmpfr::mpfr2array(ifelse(is.null(B_p), Rmpfr::mpfr(0,precBits), sum(pdet[i,] * B_p)),dim = c(length(pdet[i,]),1))
      Bpt_temp <- Rmpfr::mpfr2array(ifelse(is.null(B_pt), Rmpfr::mpfr(0,precBits), sum(pt[t,] * B_pt)),dim = c(length(pt[t,]),1))
      # size takes possible values of N (0 to K) at time t==1
      p <- optimizeAPA::plogis_APA(Bp_temp+Bpt_temp, precBits = precBits)
      g1 <- drbinomAPA(x = Y[i,1], size = (0:K)*red[i,1], prob = p, red = red[i,1], precBits = precBits)
      g2 <- drpoisAPA(x = 0:K, lambda = exp(sum(lamb[i,] * B_l)), red = red[i,1], precBits = precBits)

      # apply recursive definition of likelihood
      return( Rmpfr::mpfr(log(sum(g1 * g2 * g_star)), precBits) ) # + 1e-320
    }, K=K, T=T, Y=Y, lamb=lamb, B_l=B_l, pdet=pdet, B_p=B_p, pt=pt, B_pt=B_pt, red=red, g3=g3, g1_t_star=g1_t_star, g1_t=g1_t,g1=g1,g2=g2,g_star=g_star,precBits=precBits)

    ll <- sum(ll_i)
  }

  if(VERBOSE) { print(paste0("log likelihood: ",as.numeric(ll))) }

  if(VERBOSE) {
    if(!APA) {
      print(paste0("log likelihood: ", ll))
      print(paste0("parameters: ", par))
    }
    if(APA) {
      print(paste0("log likelihood: ", Rmpfr::formatMpfr(ll)))
      print(paste0("parameters: ", Rmpfr::formatMpfr(par)))
    }
  }

  return(-1*ll)
}


#' Not Vectorized.
tp_jk <- function(j,k,omeg,gamm,red) {
  p  <- 0
  rb <- drbinom(x = 0:min(j,k), size = j*red, prob = omeg, red = red)
  rp <- drpois(x = k-0:min(j,k), lambda = gamm, red = red)

  p <- sum(rb * rp)

  return(p)
}

#' Internal function, calculates transition probabilities from pop size j to pop size k in the open population likelihood.
tp_jk_V_ <- function(j_vec,k_vec,omeg,gamm,red) {

  p <- mapply(FUN = function(j,k,omeg,gamm,red) {
    rb <- drbinom(x = 0:min(j,k), size = j*red, prob = omeg, red = red)
    rp <- drpois(x = k-0:min(j,k), lambda = gamm, red = red)

    return(sum(rb * rp))
  }, j=j_vec, k=k_vec, omeg=omeg, gamm=gamm,red=red)


  return(p)
}
tp_jk_V <- compiler::cmpfun(tp_jk_V_)


#' Internal function, calculates transition probability matrix (transition from row pop to column pop)
tp_MAT_ <- function(M, omeg, B_o, omet, B_ot, gamm, B_g, gamt, B_gt, red, PARALLELIZE=FALSE) {
  K1 <- 1:(nrow(M))
  if(PARALLELIZE) {
    M <- foreach(a = K1-1, .combine = rbind) %dopar% {
        tp_jk_V(j_vec = a, k_vec = 1:nrow(M)-1, omeg = plogis(sum(omeg * B_o)+sum(B_ot*omet)), gamm = exp(sum(B_g*gamm)+sum(B_gt*gamt)), red = red)
    }
  } else {
    M <- outer(X = K1-1,Y = K1-1, FUN = tp_jk_V, omeg = plogis(sum(omeg * B_o)+sum(B_ot*omet)), gamm = exp(sum(B_g*gamm)+sum(B_gt*gamt)), red)
  }
  return(M)
}
tp_MAT <- compiler::cmpfun(tp_MAT_)


#' Internal function, calculates transition probabilities from pop size j to pop size k in the open population likelihood.
tp_jk_V_APA <- function(j_vec,k_vec,omeg,gamm,red,precBits=53) {

  p <- mapply(FUN = function(j,k,omeg,gamm,red,precBits) {
    rb <- drbinomAPA(x = 0:min(j,k), size = j*red, prob = omeg, red = red, precBits = precBits)
    rp <- drpoisAPA(x = k-0:min(j,k), lambda = gamm, red = red, precBits = precBits)

    return(sum(rb * rp))
  }, j=j_vec, k=k_vec, omeg=omeg, gamm=gamm,red=red,precBits=precBits)

  p <- Rmpfr::mpfr2array(p, dim=c(1,min(length(j_vec),length(k_vec))))
  return(p)
}

#' Internal function, calculates transition probability matrix (transition from row pop to column pop)
tp_MAT_APA <- function(M, omeg, B_o, omet, B_ot, gamm, B_g, gamt, B_gt, red, PARALLELIZE=FALSE, precBits=128) {
  K1 <- 1:(nrow(M))
  o_temp <- ifelse(is.null(B_ot), 0, B_ot*omet)
  g_temp <- ifelse(is.null(B_gt), 0, B_gt*gamt)
  if(class(o_temp[1])=="numeric") {
    o_temp <- unlist(o_temp)
  } else {
    o_temp <- methods::new("mpfr", unlist(o_temp))
  }
  if(class(g_temp[1])=="numeric") {
    g_temp <- unlist(g_temp)
  } else {
    g_temp <- methods::new("mpfr", unlist(g_temp))
  }
  if(class(B_g[1])=="numeric") {
    B_g <- unlist(B_g)
  } else {
    B_g <- methods::new("mpfr", unlist(B_g))
  }
  if(class(B_o[1])=="numeric") {
    B_o <- unlist(B_o)
  } else {
    B_o <- methods::new("mpfr", unlist(B_o))
  }
  if(PARALLELIZE) {
    M <- foreach(a = K1-1, .combine = rbind, .packages = "Rmpfr", "optimizeAPA") %dopar% {
      tp_jk_V_APA(j_vec = a, k_vec = 1:nrow(M)-1, omeg = optimizeAPA::plogis_APA(sum(omeg * B_o)+sum(o_temp), precBits = precBits), gamm = exp(sum(B_g*gamm)+sum(g_temp)), red = red, precBits = precBits)
    }
    #M <- Rmpfr::mpfr2array(M, dim=c(length(K1),length(K1)))
  } else {
    M <- outer(X = K1-1,Y = K1-1, FUN = tp_jk_V_APA, omeg = optimizeAPA::plogis_APA(sum(omeg * B_o)+sum(o_temp), precBits = precBits), gamm = exp(sum(B_g*gamm)+sum(g_temp)), red, precBits = precBits)
    #M <- Rmpfr::mpfr2array(M, dim=c(length(K1),length(K1)))
  }
  return(M)
}

#' Find maximum likelihood estimates for model parameters log(lambda) and logit(pdet). Uses optim.
#' @param starts Vector of starting values for optimize. Has two elements, log(lambda) and logit(pdet).
#' @param nit    R by T matrix of full counts with R sites/rows and T sampling occassions/columns.
#' @param lambda_site_covariates Either NULL (no lambda site covariates) or a list of vectors of length R, where each vector represents one site covariate, and where the vector entries correspond to covariate values for each site. Note that the covariate structure is assumed to be log(lambda_it) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param pdet_site_covariates   Either NULL (no pdet site covariates) or a list of vectors of length R, where each vector represents one site covariate, and where the vector entries correspond to covariate values for each site. Note that the covariate structure is assumed to be logit(lambda_it) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param K      Upper bound on summations, either a single number, or a vector of K values, one for each site (full count value, eg if K=300 for full counts, K=reduction(300,red) for reduced counts).
#' @param red    reduction factor, either a number, or a vector of reduction factors (R sites reductions).
#' @param VERBOSE If true, prints the log likelihood to console at each optim iteration.
#' @param PARALLELIZE If true, calculation will be split over threads by sites. Will use as many threads as have been made available (initialize with START_PARALLEL(num_cores)).
#' @param APA    If true, will use arbitrary precision arithmetic in the likelihood calculations. Use precBits to specify the number of bits of precision.
#' @param precBits If APA=TRUE, then this will specify the number of bits of precision.
#' @param tolerance specifies tolerance for convergence (defulat is 10^-6), all components of estimated gradient must be less than tolerance for convergence. If APA=TRUE, then tolerance can be made very small (eg 10^-20) using: tolerance=Rmpfr::mpfr(10^-20, precBits=128). NOTE: currently tolerance is only used if method="DFP".
#' @param ...    Additional input for optim.
#' @examples
#' START_PARALLEL(num_cores=4)
#' Y    <- gen_Nmix_closed(8,8,250,0.5)
#' out  <- fit_red_Nmix_closed(Y$nit, red=10, K=300, starts = c(log(250),boot::logit(0.5)), PARALLELIZE=TRUE)
#' out2 <- fit_red_Nmix_closed(Y$nit, red=c(10,10,10,10,20,20,40,40), K=300, starts = c(log(250),boot::logit(0.5)), PARALLELIZE=TRUE)
#' END_PARALLEL()
#' @export
fit_red_Nmix_closed <- function(nit, lambda_site_covariates=NULL, pdet_site_covariates=NULL, red, K, starts=c(1,0), VERBOSE=FALSE, PARALLELIZE=FALSE, method="BFGS", APA=FALSE, precBits=128, tolerance=10^-6, ...) {

  Y_m <- nit
  row.names(Y_m) <- 1:nrow(nit)
  colnames(Y_m) <- 1:ncol(nit)


  lamb_names  <- c("B_l_0")
  pdet_names  <- c("B_p_0")

  if(!is.null(pdet_site_covariates)) {
    if(!is.list(pdet_site_covariates)) {stop("invalid pdet_site_covariates - must be either NULL or a list of vectors of length R (number of sites).")}
    if(any( !unlist(lapply(X = pdet_site_covariates, FUN = function(X) {is.vector(X) && length(X)==nrow(nit)})) )) {
      stop("invalid pdet_site_covariates - must be either NULL or a list of vectors of length R (number of sites).")
    }
    # update default starting values
    if(identical(starts,c(1,0))) starts <- c(starts, rep(0, times=length(pdet_site_covariates)))
    numCov      <- length(pdet_site_covariates)
    pdet_names  <- c(pdet_names, paste0("B_p_s_",1:numCov))
  }

  if(!is.null(lambda_site_covariates)) {
    if(!is.list(lambda_site_covariates)) {stop("invalid lambda_site_covariates - must be either NULL or a list of vectors of length R (number of sites).")}
    if(any( !unlist(lapply(X = lambda_site_covariates, FUN = function(X) {is.vector(X) && length(X)==nrow(nit)})) )) {
      stop("invalid lambda_site_covariates - must be either NULL or a list of vectors of length R (number of sites).")
    }
    # update default starting values
    if(identical(starts,c(1,0))) starts <- c(rep(1, times=length(lambda_site_covariates)), starts)
    numCov      <- length(lambda_site_covariates)
    lamb_names  <- c(lamb_names, paste0("B_l_s_",1:numCov))
  }

  # Y_df <- reshape2::melt(Y_m)
  # colnames(Y_df) <- c("site", "time", "count")


  if(length(red)==1) {
    red <- rep(x = red, times = nrow(Y_m))
  }

  if(length(red)!=nrow(Y_m)) { stop("reduction must be either one number, or a vector with length equal to nrow(nit)") }

  red <- matrix(red, nrow=nrow(Y_m), ncol=ncol(Y_m))

  # redu <- numeric(nrow(Y_df))
  # temp <- numeric(nrow(Y_df))
  # for(i in 1:nrow(Y_df)) {
  #   redu[i] <- red[Y_df$site[i],Y_df$time[i]]
  #   temp[i] <- reduction(x = Y_df$count[i], red = redu[i])
  # }
  # Y_df$count <- temp
  # Y_df$reduc <- redu

  red_K  <- reduction(x = matrix(K, nrow=nrow(red), ncol=ncol(red)), red = red)

  for(i in 1:nrow(nit)) {
    Y_m[i,] <- reduction(nit[i,], red[i,1])
  }

  opt <- NULL
  if(!APA) {
    if(method=="DFP") {
      opt <- optimizeAPA::optim_DFP_NAPA(starts     = starts,
                           func        = red_Like_closed,
                           nit         = Y_m,
                           l_s_c       = lambda_site_covariates,
                           p_s_c       = pdet_site_covariates,
                           K           = red_K,
                           red         = red,
                           VERBOSE     = VERBOSE,
                           PARALLELIZE = PARALLELIZE,
                           APA         = FALSE,
                           tolerance   = tolerance,
                           ...)

      if(length(opt$x)==length(c(lamb_names, pdet_names))) {
        rownames(opt$x) <- c(lamb_names, pdet_names)
      } else {
        for(i in 1:length(opt$x)) {
          rownames(opt$x[[i]]) <- c(lamb_names, pdet_names)
        }
      }
    } else {
      opt <- optim(par      = starts,
                   fn       = red_Like_closed,
                   nit      = Y_m,
                   l_s_c    = lambda_site_covariates,
                   p_s_c    = pdet_site_covariates,
                   K        = red_K,
                   red      = red,
                   VERBOSE  = VERBOSE,
                   method   = method,
                   PARALLELIZE = PARALLELIZE,
                   control = list(reltol=tolerance),
                   ...)
      names(opt$par) <- c(lamb_names, pdet_names)
    }
  } else {
    opt <- optimizeAPA::optim_DFP_APA(starts      = starts,
                         func        = red_Like_closed,
                         nit         = Y_m,
                         l_s_c       = lambda_site_covariates,
                         p_s_c       = pdet_site_covariates,
                         K           = red_K,
                         red         = red,
                         VERBOSE     = VERBOSE,
                         PARALLELIZE = PARALLELIZE,
                         APA         = TRUE,
                         tolerance   = tolerance,
                         precBits    = precBits,
                         ...)

    if(length(opt$x)==length(c(lamb_names, pdet_names))) {
      rownames(opt$x) <- c(lamb_names, pdet_names)
    } else {
      for(i in 1:length(opt$x)) {
        rownames(opt$x[[i]]) <- c(lamb_names, pdet_names)
      }
    }
  }

  return(opt)
}

#' Used to initialize parallel computing (useful for likelihood calculations which can be computationally intensive).
#' @param num_cores Number of cores to use for parallel processing.
#' @export
START_PARALLEL <- function(num_cores) {
  library(doParallel)
  library(foreach)
  doParallel::registerDoParallel(cores=num_cores)
}

#' Used to end parallel computing.
#' @export
END_PARALLEL <- function() {
  doParallel::stopImplicitCluster()
}

#' Find maximum likelihood estimates for model parameters log(lambda), log(gamma), logit(omega), and logit(pdet). Uses optim.
#' @param starts Vector with four elements (if no covariates), log(lambda), log(gamma), logit(omega), and logit(pdet). When there are X covariates B_x for a parameter, will need X+1 starting values for that parameter (+1 for the constant term B0).
#' @param nit    R by T matrix of full counts with R sites/rows and T sampling occassions/columns.
#' @param lambda_site_covariates Either NULL (no lambda site covariates) or a list of vectors of length R, where each vector represents one site covariate, and where the vector entries correspond to covariate values for each site. Note that the covariate structure is assumed to be log(lambda_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param gamma_site_covariates  Either NULL (no gamma site covariates) or a list of vectors of length R, where each vector represents one site covariate, and where the vector entries correspond to covariate values for each site. Note that the covariate structure is assumed to be log(gamma_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param omega_site_covariates  Either NULL (no omega site covariates) or a list of vectors of length R, where each vector represents one site covariate, and where the vector entries correspond to covariate values for each site. Note that the covariate structure is assumed to be logit(omega_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param pdet_site_covariates   Either NULL (no pdet site covariates) or a list of vectors of length R, where each vector represents one site covariate, and where the vector entries correspond to covariate values for each site. Note that the covariate structure is assumed to be logit(pdet_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param gamma_time_covariates  Either NULL (no gamma time covariates) or a list of vectors of length T, where each vector represents one time covariate, and where the vector entries correspond to covariate values for each time. Note that the covariate structure is assumed to be log(gamma_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param omega_time_covariates  Either NULL (no omega time covariates) or a list of vectors of length T, where each vector represents one time covariate, and where the vector entries correspond to covariate values for each time. Note that the covariate structure is assumed to be logit(omega_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param pdet_time_covariates   Either NULL (no pdet time covariates) or a list of vectors of length T, where each vector represents one time covariate, and where the vector entries correspond to covariate values for each time. Note that the covariate structure is assumed to be logit(omega_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param K           Upper bound on summations, will be reduced by reduction factor red.
#' @param red         reduction factor, either a number or a vector of length R.
#' @param VERBOSE     If TRUE, prints the log likelihood to console at each optim iteration.
#' @param PARALLELIZE If TRUE, calculation will be split over threads by sites and times. This will not improve computation time if there are no site or time covariates. Will use as many threads as have been made available (initialize with START_PARALLEL(num_cores)).
#' @param APA         If true, will use arbitrary precision arithmetic in the likelihood calculations. Use precBits to specify the number of bits of precision.
#' @param precBits    If APA=TRUE, then this will specify the number of bits of precision.
#' @param tolerance   specifies tolerance for convergence (defulat is 10^-6), all components of estimated gradient must be less than tolerance for convergence. If APA=TRUE, then tolerance can be made very small (eg 10^-20) using: tolerance=Rmpfr::mpfr(10^-20, precBits=128). NOTE: currently tolerance is only used if method="DFP".
#' @param ...         Additional input for optim.
#' @examples
#'
#' Y   <- gen_Nmix_open(num_sites = 5, num_times = 5, lambda = 20, pdet = 0.7, omega = 0.7, gamma = 2)
#' out <- fit_red_Nmix_open(nit = Y$nit, red = c(1), K = 40, starts = c(0.5, 0.5, 0.5, 0.5))
#'
#'
#' # example with site covariates:
#' Y1 <- gen_Nmix_open(num_sites = 2, num_times = 5, lambda = 10, gamma = 5, omega = 0.50, pdet = 0.75)
#' Y2 <- gen_Nmix_open(num_sites = 3, num_times = 5, lambda = 5, gamma = 10, omega = 0.75, pdet = 0.50)
#' Y  <- rbind(Y1$nit, Y2$nit)
#' START_PARALLEL(num_cores=5)
#' mod1 <- fit_red_Nmix_open(nit = Y,
#'                           lambda_site_covariates = list(l1=c(0,0,1,1,1)),
#'                           gamma_site_covariates  = list(gs=c(0,0,1,1,1)),
#'                           gamma_time_covariates  = NULL,
#'                           omega_site_covariates  = list(os=c(0,0,1,1,1)),
#'                           omega_time_covariates  = NULL,
#'                           pdet_site_covariates   = list(ps=c(0,0,1,1,1)),
#'                           pdet_time_covariates   = NULL,
#'                           red = 4,
#'                           K   = 50,
#'                           starts  = NULL,
#'                           method  = "BFGS",
#'                           VERBOSE = FALSE,
#'                           PARALLELIZE = TRUE)
#' END_PARALLEL()
#' # lambda sites 1 and 2 estimate:
#' exp(mod1$par[1])
#' # lambda sites 3, 4, and 5 estimate:
#' exp(sum(mod1$par[1:2]))
#' # gamma sites 1 and 2 estimate:
#' exp(mod1$par[3])
#' # gamma sites 3, 4, and 5 estimate:
#' exp(sum(mod1$par[3:4]))
#' # omega sites 1 and 2 estimate:
#' plogis(mod1$par[5])
#' # omega sites 3, 4, and 5 estimate:
#' plogis(sum(mod1$par[5:6]))
#' # pdet sites 1 and 2 estimate:
#' plogis(mod1$par[7])
#' # pdet sites 3, 4, and 5 estimate:
#' plogis(sum(mod1$par[7:8]))
#' @export
fit_red_Nmix_open <- function(nit, lambda_site_covariates=NULL, gamma_site_covariates=NULL, omega_site_covariates=NULL, pdet_site_covariates=NULL, gamma_time_covariates=NULL, omega_time_covariates=NULL, pdet_time_covariates=NULL, red, K, starts=NULL, VERBOSE=FALSE, PARALLELIZE=FALSE, APA=FALSE, precBits=128, tolerance=10^-6, method="BFGS", ...) {
   if(length(red)==1) {
     red <- rep(red, times=nrow(nit))
   }
   red    <- matrix(red, nrow=nrow(nit), ncol=ncol(nit))
   red_K  <- reduction(x = matrix(K, nrow=nrow(red), ncol=ncol(red)), red = red)

   lamb_starts <- c(1)
   gamm_starts <- c(1)
   omeg_starts <- c(0)
   pdet_starts <- c(0)

   lamb_names  <- c("B_l_0")
   gamm_names  <- c("B_g_0")
   omeg_names  <- c("B_o_0")
   pdet_names  <- c("B_p_0")

   if(!is.null(lambda_site_covariates)) {
     if(!is.list(lambda_site_covariates)) {stop("invalid lambda_site_covariates - must be either NULL or a list of vectors of length R (number of sites).")}
     if(any( !unlist(lapply(X = lambda_site_covariates, FUN = function(X) {is.vector(X) && length(X)==nrow(nit)})) )) {
       stop("invalid lambda_site_covariates - must be either NULL or a list of vectors of length R (number of sites).")
     }
     # update default starting values
     lamb_starts <- rep(1, times=length(lambda_site_covariates)+1)
     numCov      <- length(lamb_starts)-1
     lamb_names  <- c(lamb_names, paste0("B_l_s_",1:numCov))
   }

   if(!is.null(gamma_site_covariates)) {
     if(!is.list(gamma_site_covariates)) {stop("invalid gamma_site_covariates - must be either NULL or a list of vectors of length R (number of sites).")}
     if(any( !unlist(lapply(X = gamma_site_covariates, FUN = function(X) {is.vector(X) && length(X)==nrow(nit)})) )) {
       stop("invalid gamma_site_covariates - must be either NULL or a list of vectors of length R (number of sites).")
     }
     # update default starting values
     gamm_starts <- rep(1, times=length(gamma_site_covariates)+1)
     numCov      <- length(gamm_starts)-1
     gamm_names  <- c(gamm_names, paste0("B_g_s_",1:numCov))
   }

   if(!is.null(gamma_time_covariates)) {
     if(!is.list(gamma_time_covariates)) {stop("invalid gamma_time_covariates - must be either NULL or a list of vectors of length T (number of sampling occasions, last entry will be ignored).")}
     if(any( !unlist(lapply(X = gamma_time_covariates, FUN = function(X) {is.vector(X) && length(X)==(ncol(nit))})) )) {
       stop("invalid gamma_time_covariates - must be either NULL or a list of vectors of length T (number of sampling occasions, last entry will be ignored).")
     }
     # update default starting values
     gamm_starts <- c(gamm_starts, rep(1, times=length(gamma_time_covariates)))
     numCov      <- length(gamm_starts)-length(gamm_names)
     gamm_names  <- c(gamm_names, paste0("B_g_t_",1:numCov))
   }

   if(!is.null(omega_site_covariates)) {
     if(!is.list(omega_site_covariates)) {stop("invalid omega_site_covariates - must be either NULL or a list of vectors of length R (number of sites).")}
     if(any( !unlist(lapply(X = omega_site_covariates, FUN = function(X) {is.vector(X) && length(X)==nrow(nit)})) )) {
       stop("invalid omega_site_covariates - must be either NULL or a list of vectors of length R (number of sites).")
     }
     # update default starting values
     omeg_starts <- rep(0, times=length(omega_site_covariates)+1)
     numCov      <- length(omeg_starts)-1
     omeg_names  <- c(omeg_names, paste0("B_o_s_",1:numCov))
   }

   if(!is.null(omega_time_covariates)) {
     if(!is.list(omega_time_covariates)) {stop("invalid omega_time_covariates - must be either NULL or a list of vectors of length T (number of sampling occasions, last entry will be ignored).")}
     if(any( !unlist(lapply(X = omega_time_covariates, FUN = function(X) {is.vector(X) && length(X)==(ncol(nit))})) )) {
       stop("invalid omega_time_covariates - must be either NULL or a list of vectors of length T (number of sampling occasions, last entry will be ignored).")
     }
     # update default starting values
     omeg_starts <- c(omeg_starts, rep(0, times=length(omega_time_covariates)))
     numCov      <- length(omeg_starts)-length(omeg_names)
     omeg_names  <- c(omeg_names, paste0("B_o_t_",1:numCov))
   }

   if(!is.null(pdet_site_covariates)) {
     if(!is.list(pdet_site_covariates)) {stop("invalid pdet_site_covariates - must be either NULL or a list of vectors of length R (number of sites).")}
     if(any( !unlist(lapply(X = pdet_site_covariates, FUN = function(X) {is.vector(X) && length(X)==nrow(nit)})) )) {
       stop("invalid pdet_site_covariates - must be either NULL or a list of vectors of length R (number of sites).")
     }
     # update default starting values
     pdet_starts <- rep(0, times=length(pdet_site_covariates)+1)
     numCov      <- length(pdet_starts)-1
     pdet_names  <- c(pdet_names, paste0("B_p_s_",1:numCov))
   }


   if(!is.null(pdet_time_covariates)) {
     if(!is.list(pdet_time_covariates)) {stop("invalid pdet_time_covariates - must be either NULL or a list of vectors of length T (number of sampling occasions, last entry will be ignored).")}
     if(any( !unlist(lapply(X = pdet_time_covariates, FUN = function(X) {is.vector(X) && length(X)==(ncol(nit))})) )) {
       stop("invalid pdet_time_covariates - must be either NULL or a list of vectors of length T (number of sampling occasions, last entry will be ignored).")
     }
     # update default starting values
     pdet_starts <- c(pdet_starts, rep(0, times=length(pdet_time_covariates)))
     numCov      <- length(pdet_starts)-length(pdet_names)
     pdet_names  <- c(pdet_names, paste0("B_p_t_",1:numCov))
   }


   if(is.null(starts)) {
     starts <- c(lamb_starts, gamm_starts, omeg_starts, pdet_starts)
   }

   Y_m <- matrix(0, nrow=nrow(nit), ncol=ncol(nit))
   for(i in 1:nrow(nit)) {
     Y_m[i,] <- reduction(nit[i,], red[i,1])
   }

   NAMES <- c(lamb_names, gamm_names, omeg_names, pdet_names)

   opt <- NULL
   if(!APA) {
     if(method=="DFP"){
       opt <- optimizeAPA::optim_DFP_NAPA(starts  = starts,
                                          func    = red_Like_open,
                                          nit     = Y_m,
                                          l_s_c   = lambda_site_covariates,
                                          g_s_c   = gamma_site_covariates,
                                          g_t_c   = gamma_time_covariates,
                                          o_s_c   = omega_site_covariates,
                                          o_t_c   = omega_time_covariates,
                                          p_s_c   = pdet_site_covariates,
                                          p_t_c   = pdet_time_covariates,
                                          K           = red_K,
                                          red         = red,
                                          VERBOSE     = VERBOSE,
                                          PARALLELIZE = PARALLELIZE,
                                          APA         = FALSE,
                                          precBits    = precBits,
                                          tolerance   = tolerance,
                                          ...)
       if(length(opt$x)==length(NAMES)) {
         rownames(opt$x) <- NAMES
       } else {
         for(i in 1:length(opt$x)) {
           rownames(opt$x[[i]]) <- NAMES
         }
       }
     } else {
       opt <- optim(par     = starts,
                    fn      = red_Like_open,
                    nit     = Y_m,
                    l_s_c   = lambda_site_covariates,
                    g_s_c   = gamma_site_covariates,
                    g_t_c   = gamma_time_covariates,
                    o_s_c   = omega_site_covariates,
                    o_t_c   = omega_time_covariates,
                    p_s_c   = pdet_site_covariates,
                    p_t_c   = pdet_time_covariates,
                    K       = red_K,
                    red     = red,
                    VERBOSE = VERBOSE,
                    method  = method,
                    PARALLELIZE = PARALLELIZE,
                    ...)
       names(opt$par) <- NAMES
     }
   } else { # APA
     opt <- optimizeAPA::optim_DFP_APA(starts  = starts,
                                       func    = red_Like_open,
                                       nit     = Y_m,
                                       l_s_c   = lambda_site_covariates,
                                       g_s_c   = gamma_site_covariates,
                                       g_t_c   = gamma_time_covariates,
                                       o_s_c   = omega_site_covariates,
                                       o_t_c   = omega_time_covariates,
                                       p_s_c   = pdet_site_covariates,
                                       p_t_c   = pdet_time_covariates,
                                       K           = red_K,
                                       red         = red,
                                       VERBOSE     = VERBOSE,
                                       PARALLELIZE = PARALLELIZE,
                                       APA         = TRUE,
                                       precBits    = precBits,
                                       tolerance   = tolerance,
                                       ...)
     if(length(opt$x)==length(NAMES)) {
       rownames(opt$x) <- NAMES
     } else {
       for(i in 1:length(opt$x)) {
         rownames(opt$x[[i]]) <- NAMES
       }
     }
   }

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
  red   <- matrix(red, nrow=nrow(nit), ncol=ncol(nit))
  K     <- reduction(x = matrix(K, nrow=nrow(red), ncol=ncol(red)), red = red)

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
      L[j,k] <- -1*red_Like_closed(nit = nit, par = c(log(par2M),boot::logit(par1M)), K = K, red = red, l_s_c = NULL, p_s_c = NULL)
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

