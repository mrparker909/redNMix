
#' @title drbinom
#' @description Reduced binomial probability distribution function \eqn{rBinomial(x;N,p,R(x;r))}, takes reduced quantiles rather than full quantiles.
#' @param x    Reduced count quantile (alternatively input reduction(x,red) if x is a full count quantile).
#' @param size Number of trials (full count size).
#' @param prob Probability of success for each trial.
#' @param red  The factor r by which x has been reduced.
#' @return The probability of observing quantile \eqn{x}.
#' @examples
#' Y <- drbinom(0:20, 20, 0.3, 2)
#' plot(Y, xlab="Y~rBinom(N=20,p=0.3,r=2)", ylab="P[Y=y]")
#' @export
drbinom <- function(x, size, prob, red, log=FALSE) {
  start <- ceiling(x*red - red/2)-1
  end   <- round2(x*red + red/2)-1
  p <- pbinom(end, size, prob) - pbinom(start, size, prob)
  if(log) { return(log(p)) }
  return(p)
}

#' @title drbinomAPA
#' @description Arbitrary precision reduced binomial probability distribution function \eqn{rBinomial(x;N,p,R(x;r))}, takes reduced quantiles rather than full quantiles.
#' @param x Reduced count quantile (alternatively input reduction(x,r) if x is a full count quantile).
#' @param size Number of trials.
#' @param prob Probability of success for each trial.
#' @param red The factor r by which x has been reduced.
#' @param precBits Number of bits of precision for arbitrary precision arithmetic.
#' @return The probability of observing quantile \eqn{x}.
#' @examples
#' Y <- drbinomAPA(0:3, 20, 0.3, 10, precBits=64)
#' Y
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
    start <- Rmpfr::pmax(0,round2(X*red - red/2))
    end   <- Rmpfr::pmin(size[i],round2(X*red + red/2)-1)
    if(length(end)==0) end = round2(X*red + red/2)-1
    if(start > end) {
      pt[i] = 0
    } else {
      pt[i] = sum(Rmpfr::dbinom(x = start:end, size = size[i], prob = prob))
    }
  }
  pt[size0] = ifelse(x[size0]==0, 1, 0)
  if(log) { return(log(pt)) }
  return(pt)
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
