
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

drbinomAPA_ <- function(x, size, prob, red, precBits=128, log=FALSE) {
  require(Rmpfr)

  pt <- NULL
  i <- 0
  for(X in x) {
    i <- i + 1
    start <- (X*red - red/2)
    end   <- (start + red)

    temp  <- optimizeAPA::dbinom_APA((start+1):end, size*red, prob=prob, precBits=precBits)
    pt[i] <- sum(new("mpfr", unlist(temp)))
  }
  p <- new("mpfr", unlist(pt))
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
#' Y <- drbinomAPA(0:20, 20, 0.3, 10, precBits=64)
#' @export
drbinomAPA <- compiler::cmpfun(drbinomAPA_)


#' internal function
drpoisAPA_1  <- function(x, lambda, red, precBits=128, log=FALSE) {
  require(Rmpfr)

  pt <- NULL
  i <- 0
  for(X in x) {
    i <- i + 1
    start <- as.integer(X*red - red/2) +1
    end   <- (start + red)-1
    temp  <- optimizeAPA::dpois_APA((start):(end), lambda, precBits=precBits)
    pt[i] <- sum(new("mpfr", unlist(temp)))
  }
  p <- new("mpfr", unlist(pt))
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
#' @param precBits Number of bits of precision for arbitrary precision arithmetic.
#' @return The probability of observing quantile \eqn{x}
#' @examples
#' # probability of observing 105 from pois(lam=90)
#' x <- drpoisAPA(x = 105, lambda = 90, precBits=64)
#' @export
drpoisAPA <- compiler::cmpfun(drpoisAPA_1)


#' internal function
drpois_1  <- function(x, lambda, red, log=FALSE) {
  start <- x*red-red/2
  end   <- start + red # reduction(x,red, FUN=round2)*red+red/2

  p <- ppois(start, lambda, lower.tail = FALSE) - ppois(end, lambda, lower.tail = FALSE) #ppois(end, lambda) - ppois(start, lambda) #sum(dpois(start:end, lambda)) #

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
#' Y <- drpois(seq(1,20,1), 55, 10)
#' plot(Y, xlab="Y=rPois(lambda=55, r=10)", main="X~Poisson(lambda=55)", ylab="P[Y=y]")
#' @export
drpois <- compiler::cmpfun(drpois_1)



red_Like_closed_ <- function(par, nit, l_s_c, p_s_c, K, red, FUN=round, VERBOSE=FALSE, PARALLELIZE=FALSE, APA=FALSE, precBits=64) {
  T <- length(unique(nit$time))
  R <- length(unique(nit$site))

  if(!APA) {
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

    Y    <- nit
    Y_df <- nit

    l <- 0
    if(PARALLELIZE) {
      li <- foreach(i=1:R) %dopar% {
        Yi_df <- Y_df[which(Y_df$site==i),]
        li <- 0
        ni <- max(Yi_df$count)

        for(Ni in ni:K[i,1]) {
          lit <- with(data = Yi_df, expr = drbinom(x = count, size = Ni, prob = plogis(sum(B_p*pdet[i,])), red=reduc))

          lit_ <- prod(lit)

          li <- li + lit_*drpois(x = Ni, lambda = exp(sum(B_l*lamb[i,])), red=red[i,1])
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
          lit <- with(data = Yi_df, expr = drbinom(x = count, size = Ni, prob = plogis(sum(B_p*pdet[i,])), red=reduc))

          lit_ <- prod(lit)
          li <- li + lit_*drpois(x = Ni, lambda = exp(sum(B_l*lamb[i,])), red=red[i,1])
        }
        l <- l+log(li)
      }
    }
  } else { # DO APA CALCULATIONS

    # extract lambda estimates from par, setup lambda covariate matrix lamb, and covariate vector B_l
    lamb <- matrix(rep(1,times=R),ncol=1) # covariate coefficients for lambda
    B_l <- Rmpfr::mpfr(par[1], precBits=precBits) # covariates for lambda
    if(!is.null(l_s_c)) {
      lamb <- cbind(lamb, do.call(cbind, l_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix
      B_l  <- sapply(X = 1:(length(l_s_c)+1), FUN = function(X,par) { # one coeff per covariate +1 for baseline B0
        Rmpfr::mpfr(par[X], precBits=precBits)
      }, par=par)
    }

    # extract pdet estimates from par, setup pdet covariate matrix pdet, and covariate vector B_p
    pdet <- matrix(rep(1,times=R),ncol=1) # covariate coefficients for lambda
    B_p <- par[length(B_l)+1] # covariates for pdet
    if(!is.null(p_s_c)) {
      pdet <- cbind(pdet, do.call(cbind, p_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix
      B_p  <- sapply(X = 1:(length(p_s_c)+1)+length(B_l), FUN = function(X,par) { # one coeff per covariate +1 for baseline B0
        Rmpfr::mpfr(par[X], precBits=precBits)
      }, par=par)
    }

    Y    <- nit
    Y_df <- nit

    l <- Rmpfr::mpfr(0, precBits=precBits)
    if(PARALLELIZE) {
      li <- foreach(i=1:R) %dopar% {
        Yi_df <- Y_df[which(Y_df$site==i),]
        li <- Rmpfr::mpfr(0, precBits=precBits)
        ni <- max(Yi_df$count)

        for(Ni in ni:K[i,1]) {
          lit <- with(data = Yi_df, expr = drbinomAPA(x = count, size = Ni, prob = plogisAPA(sum(B_p*pdet[i,]), precBits=precBits), red=reduc, precBits=precBits))

          lit_ <- prod(lit)

          li <- li + lit_*drpoisAPA(x = Ni, lambda = exp(sum(B_l*lamb[i,])), red=red[i,1], precBits=precBits)
        }
        return(log(li))
      }
      l <- sum(unlist(li))
      ###
    } else {

      for(i in 1:R) {
        Yi_df <- Y_df[which(Y_df$site==i),]
        li <- Rmpfr::mpfr(0, precBits=precBits)
        ni <- max(Yi_df$count)

        for(Ni in ni:K[i,1]) {
          lit <- with(data = Yi_df, expr = drbinomAPA(x = count, size = Ni, prob = plogisAPA(sum(B_p*pdet[i,]), precBits=precBits), red=reduc, precBits=precBits))

          lit_ <- prod(lit)
          li <- li + lit_*drpoisAPA(x = Ni, lambda = exp(sum(B_l*lamb[i,])), red=red[i,1], precBits=precBits)
        }
        l <- l+log(li)
      }
    }
  }
  if(VERBOSE) {print(paste0("log likelihood: ",format(l)))}
  return(-1*l)
}

#' Internal function. Used to calculate the negative of the log likelihood.
#' @param par Vector with two elements, log(lambda) and logis(pdet). If l_s_c is not NULL, need length(par) = length(l_s_c) + 2, for the B0...BK coefficients of lambda (B0 is the constant term coefficient).
#' @param nit data.frame with R*T rows, and 3 columns: "site", "time", and "count". Deprecated: R by T matrix of reduced counts with R sites/rows and T sampling occassions/columns.
#' @param l_s_c list of lambda site covariates (list of vectors of length R (number of sites))
#' @param p_s_c list of pdet site covariates (list of vectors of length R (number of sites))
#' @param K   Upper bound on summations (input reduced count upper bound).
#' @param red reduction factor matrix (R by T, with R sites/rows and T sampling occassions/columns)
#' @param VERBOSE If true, prints the calculated log likelihood to console.
#' @param PARALLELIZE If true, calculation will be split over threads by sites. Will use as many threads as have been made available (initialize with START_PARALLEL(num_cores)).
#' @param APA Default is FALSE. If true, will use arbitrary precision arithmetic in calculating the likelihood. Note that APA will be slower, however it is required for site population sizes larger than about 200. If APA = TRUE, then precBits specifies the number of bits of precision to use in the calculations.
#' @param precBits If APA=TRUE, then precBits specifies the number of bits of precision for arbitrary precision arithmetic.
#' @export
red_Like_closed <- compiler::cmpfun(red_Like_closed_)



red_Like_open_ <- function(par, nit, l_s_c, g_s_c, g_t_c, o_s_c, o_t_c, p_s_c, p_t_c, K, red, VERBOSE=FALSE, PARALLELIZE=FALSE) {
  T <- ncol(nit)
  R <- nrow(nit)

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
            else tp_MAT(M = tempMat <- tempMat, omeg = omeg[i,], B_o=B_o, omet=omet[t,], B_ot=B_ot, gamm = gamm[i,], B_g=B_g, gamt=gamt[t,], B_gt=B_gt, red = red[i,1], PARALLELIZE=TRUE)
            return(tempMat)
          }
        for(i in 1:R) {
          g3[[i]]<- g3_temp[[i]]
        }
      } else {# NOT SITE DEPENDENT
        g3_temp <- foreach(t=1:T, .packages = c("redNMix","foreach")) %dopar% {
          tempMat        <- matrix(0, nrow = K[i]+1, ncol=K[i]+1)
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
              tempMat      <- matrix(0, nrow = K[1]+1, ncol=K[1]+1)
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
            tempMat <- matrix(0, nrow = K[1]+1, ncol=K[1]+1)
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
              tempMat      <- matrix(0, nrow = K[1]+1, ncol=K[1]+1)
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
              tempMat      <- matrix(0, nrow = K[1]+1, ncol=K[1]+1)
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
      g1_t <- drbinom(x = Y[i,t], size = (0:K), prob = plogis(sum(pdet[i,] * B_p)+sum(pt[t,] * B_pt)), red = red[i,t])
      g1_t_star <- g1_t * g_star

      # update g_star
      g_star = g3[[t]] %*% g1_t_star
    }

    # size takes possible values of N (0 to K) at time t==1
    g1 <- drbinom(x = Y[i,1], size = 0:K, prob = plogis(sum(pdet[i,] * B_p)+sum(pt[t,] * B_pt)), red = red[i,1])
    g2 <- drpois(x = 0:K, lambda = exp(sum(lamb[i,] * B_l)), red = red[i,1])


    # apply recursive definition of likelihood
    return( log(sum(g1 * g2 * g_star)) ) # + 1e-320
  }, FUN.VALUE = numeric(1), K=K, T=T, Y=Y, lamb=lamb, B_l=B_l, pdet=pdet, B_p=B_p, pt=pt, B_pt=B_pt, red=red, g3=g3, g1_t_star=g1_t_star, g1_t=g1_t,g1=g1,g2=g2,g_star=g_star)

  ll <- sum(unlist(ll_i))
  ###

  if(VERBOSE) { print(paste0("log likelihood: ",ll)) }
  return(-1*ll)
}

#' Internal function. Used to calculate the negative of the log likelihood.
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


#' Find maximum likelihood estimates for model parameters log(lambda) and logit(pdet). Uses optim.
#' @param starts Vector of starting values for optimize. Has two elements, log(lambda) and logit(pdet).
#' @param nit    R by T matrix of full counts with R sites/rows and T sampling occassions/columns.
#' @param lambda_site_covariates Either NULL (no lambda site covariates) or a list of vectors of length R, where each vector represents one site covariate, and where the vector entries correspond to covariate values for each site. Note that the covariate structure is assumed to be log(lambda_it) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param pdet_site_covariates   Either NULL (no pdet site covariates) or a list of vectors of length R, where each vector represents one site covariate, and where the vector entries correspond to covariate values for each site. Note that the covariate structure is assumed to be logit(lambda_it) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param K      Upper bound on summations, either a single number, or a vector of K values, one for each site (full count value, eg if K=300 for full counts, K=reduction(300,red) for reduced counts).
#' @param red    reduction factor, either a number, or a vector of reduction factors (R sites reductions).
#' @param VERBOSE If true, prints the log likelihood to console at each optim iteration.
#' @param PARALLELIZE If true, calculation will be split over threads by sites. Will use as many threads as have been made available (initialize with START_PARALLEL(num_cores)).
#' @param ...    Additional input for optim.
#' @examples
#' START_PARALLEL(num_cores=4)
#' Y    <- gen_Nmix_closed(8,8,250,0.5)
#' out  <- fit_red_Nmix_closed(Y$nit, red=10, K=300, starts = c(log(250),boot::logit(0.5)), PARALLELIZE=TRUE)
#' out2 <- fit_red_Nmix_closed(Y$nit, red=c(10,10,10,10,20,20,40,40), K=300, starts = c(log(250),boot::logit(0.5)), PARALLELIZE=TRUE)
#' END_PARALLEL()
#' @export
fit_red_Nmix_closed <- function(nit, lambda_site_covariates=NULL, pdet_site_covariates=NULL, red, K, starts=c(1,0), VERBOSE=FALSE, PARALLELIZE=FALSE, method="BFGS", ...) {

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
    numCov      <- length(pdet_site_covariates)-1
    pdet_names  <- c(pdet_names, paste0("B_p_s_",1:numCov))
  }

  if(!is.null(lambda_site_covariates)) {
    if(!is.list(lambda_site_covariates)) {stop("invalid lambda_site_covariates - must be either NULL or a list of vectors of length R (number of sites).")}
    if(any( !unlist(lapply(X = lambda_site_covariates, FUN = function(X) {is.vector(X) && length(X)==nrow(nit)})) )) {
      stop("invalid lambda_site_covariates - must be either NULL or a list of vectors of length R (number of sites).")
    }
    # update default starting values
    if(identical(starts,c(1,0))) starts <- c(rep(1, times=length(lambda_site_covariates)), starts)
    numCov      <- length(lambda_site_covariates)-1
    lamb_names  <- c(lamb_names, paste0("B_l_s_",1:numCov))
  }

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
                l_s_c   = lambda_site_covariates,
                p_s_c   = pdet_site_covariates,
                K       = red_K,
                red     = red,
                VERBOSE = VERBOSE,
                method  = method,
                PARALLELIZE = PARALLELIZE,
                ...)
  names(opt$par) <- c(lamb_names, pdet_names)
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
#' @param starts Vector with four elements (if no covariates), log(lambda), log(gamma), logit(omega), and logit(pdet). When there are X covariates B_x for a parameter, will need X+1 starting values for that parameter (+1 for the constant term B0).
#' @param nit    R by T matrix of full counts with R sites/rows and T sampling occassions/columns.
#' @param lambda_site_covariates Either NULL (no lambda site covariates) or a list of vectors of length R, where each vector represents one site covariate, and where the vector entries correspond to covariate values for each site. Note that the covariate structure is assumed to be log(lambda_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param gamma_site_covariates  Either NULL (no gamma site covariates) or a list of vectors of length R, where each vector represents one site covariate, and where the vector entries correspond to covariate values for each site. Note that the covariate structure is assumed to be log(gamma_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param omega_site_covariates  Either NULL (no omega site covariates) or a list of vectors of length R, where each vector represents one site covariate, and where the vector entries correspond to covariate values for each site. Note that the covariate structure is assumed to be logit(omega_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param pdet_site_covariates   Either NULL (no pdet site covariates) or a list of vectors of length R, where each vector represents one site covariate, and where the vector entries correspond to covariate values for each site. Note that the covariate structure is assumed to be logit(pdet_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param gamma_time_covariates  Either NULL (no gamma time covariates) or a list of vectors of length T, where each vector represents one time covariate, and where the vector entries correspond to covariate values for each time. Note that the covariate structure is assumed to be log(gamma_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param omega_time_covariates  Either NULL (no omega time covariates) or a list of vectors of length T, where each vector represents one time covariate, and where the vector entries correspond to covariate values for each time. Note that the covariate structure is assumed to be logit(omega_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param pdet_time_covariates   Either NULL (no pdet time covariates) or a list of vectors of length T, where each vector represents one time covariate, and where the vector entries correspond to covariate values for each time. Note that the covariate structure is assumed to be logit(omega_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param K      Upper bound on summations, will be reduced by reduction factor red.
#' @param red    reduction factor, either a number or a vector of length R.
#' @param VERBOSE If TRUE, prints the log likelihood to console at each optim iteration.
#' @param PARALLELIZE If TRUE, calculation will be split over threads by sites and times. This will not improve computation time if there are no site or time covariates. Will use as many threads as have been made available (initialize with START_PARALLEL(num_cores)).
#' @param ...    Additional input for optim.
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
fit_red_Nmix_open <- function(nit, lambda_site_covariates=NULL, gamma_site_covariates=NULL, omega_site_covariates=NULL, pdet_site_covariates=NULL, gamma_time_covariates=NULL, omega_time_covariates=NULL, pdet_time_covariates=NULL, red, K, starts=NULL, VERBOSE=FALSE, PARALLELIZE=FALSE, method="BFGS", ...) {
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


   opt <- optim(par     = starts,
                fn      = red_Like_open,
                nit     = reduction(x = nit, red = red),
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
   names(opt$par) <- c(lamb_names, gamm_names, omeg_names, pdet_names)
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

