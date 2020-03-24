
#' @title red_Like_closed
#' @description Used to calculate the negative of the log likelihood for closed population models.
#' @param par Vector with two elements, log(lambda) and logis(pdet). If l_s_c is not NULL, need length(par) = length(l_s_c) + 2, for the B0...BK coefficients of lambda (B0 is the constant term coefficient).
#' @param nit R by T matrix of reduced counts with R sites/rows and T sampling occassions/columns.
#' @param l_s_c list of lambda site covariates (list of vectors of length R (number of sites), one vector per covariate)
#' @param p_s_c list of pdet site covariates (list of vectors of length R (number of sites), one vector per covariate)
#' @param K   Upper bound on summations (input the reduced count upper bound). An R by T matrix, eg: K=matrix(50, nrow=R, ncol=T).
#' @param red reduction factor matrix (R by T, with R sites/rows and T sampling occassions/columns)
#' @param VERBOSE If TRUE, prints the calculated log likelihood to console.
#' @param PARALLELIZE If TRUE, calculations will be done in parallel. Will use as many cores as have been made available (initialize with START_PARALLEL(num_cores)).
#' @param APA Default is FALSE. If TRUE, will use arbitrary precision arithmetic in calculating the likelihood. Note that APA will be slower, however it is required for site population sizes larger than about 300. If APA = TRUE, then precBits specifies the number of bits of precision to use in the calculations.
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
    if(PARALLELIZE) { # PARALLEL
      li <- foreach(i=1:R) %dopar% {
        li <- 0
        ni <- max(Y[i,])
        for(Ni in ni:max(K[i,])) {
          lit <- 1
          pt <- plogis(sum(B_p*pdet[i,]))
          for(t in 1:T) {
            if(Ni > K[i,t]) { next() }
            lit <- lit*drbinom(x    = Y[i,t],
                               size = Ni*red[i,t],
                               prob = pt,
                               red  = red[i,t])
          }
          li <- li + lit*drpois(x = Ni, lambda = exp(sum(B_l*lamb[i,])), red=red[i,1])
        }
        return(log(li))
      }
      l <- sum(unlist(li))
    } else { # NOT PARALLEL
      for(i in 1:R) {
        li <- 0
        ni <- max(Y[i,])
        for(Ni in ni:max(K[i,])) {
          lit <- 1
          pt <- plogis(sum(B_p*pdet[i,]))
          for(t in 1:T) {
            if(Ni > K[i,t]) { next() }
            lit <- lit*drbinom(x    = Y[i,t],
                               size = Ni*red[i,t],
                               prob = pt,
                               red  = red[i,t])
          }
          li <- li + lit*drpois(x = Ni, lambda = exp(sum(B_l*lamb[i,])), red=red[i,1])
        }
        l <- l+log(li)
      }
    }
  } else { # DO APA CALCULATIONS
    if(!require(Rmpfr)) { stop("Error, cannot load package Rmpfr") }

    # extract lambda estimates from par, setup lambda covariate matrix lamb, and covariate vector B_l
    lamb <- matrix(rep(1,times=R),ncol=1) #Rmpfr::mpfrArray(matrix(rep(1,times=R), ncol=1), dim = c(R,1), precBits=precBits) # covariate coefficients for lambda
    B_l <- Rmpfr::mpfr(par[1], precBits=precBits) # covariates for lambda
    if(!is.null(l_s_c)) {
      lamb <- cbind(lamb, do.call(cbind, l_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix
      B_l  <- sapply(X = 1:(length(l_s_c)+1), FUN = function(X,par) { # one coeff per covariate +1 for baseline B0
        Rmpfr::mpfr(par[X], precBits=precBits)
      }, par=par)
      B_l <- methods::new("mpfr", unlist(B_l))
    }
    l_temp <- exp(B_l %*% t(lamb))

    # extract pdet estimates from par, setup pdet covariate matrix pdet, and covariate vector B_p
    pdet <- matrix(rep(1,times=R),ncol=1) #Rmpfr::mpfrArray(matrix(rep(1,times=R), ncol=1), dim = c(R,1), precBits=precBits) # covariate coefficients for lambda
    B_p <- Rmpfr::mpfr(par[length(B_l)+1], precBits=precBits) # covariates for pdet
    if(!is.null(p_s_c)) {
      pdet <- cbind(pdet, do.call(cbind, p_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix
      B_p  <- sapply(X = 1:(length(p_s_c)+1)+length(B_l), FUN = function(X,par) { # one coeff per covariate +1 for baseline B0
        Rmpfr::mpfr(par[X], precBits=precBits)
      }, par=par)
      B_p <- methods::new("mpfr", unlist(B_p))
    }
    p_temp <- B_p %*% t(pdet)

    Y <- nit

    l <- Rmpfr::mpfr(0, precBits=precBits)
    if(PARALLELIZE) { # PARALLEL and APA
      li <- foreach(i=1:R) %dopar% {
        li <- Rmpfr::mpfr(0, precBits=precBits)
        ni <- max(Y[i,])
        for(Ni in ni:max(K[i,])) {
          lit <- Rmpfr::mpfr(1,precBits=precBits)
          pt <- optimizeAPA::plogis_APA(p_temp[i], precBits = precBits)
          for(t in 1:T) {
            if(Ni > K[i,t]) { next() }
            lit <- lit*(drbinomAPA(x    = Y[i,t],
                                   size = Ni*red[i,t],
                                   prob = pt,
                                   red  = red[i,t],
                                   precBits = precBits))
          }
          li <- li + lit*drpoisAPA(x = Ni, lambda = l_temp[i], red=red[i,1], precBits = precBits)
        }
        return(log(li))
      }

      l <- sum(methods::new("mpfr", unlist(li)))
      ###
    } else { # APA and NOT PARALLEL
      for(i in 1:R) {
        li <- Rmpfr::mpfr(0, precBits=precBits)
        ni <- max(Y[i,])
        pit <- optimizeAPA::plogis_APA(sum(B_p*pdet[i,]), precBits = precBits)
        for(Ni in ni:max(K[i,])) {
          lit <- Rmpfr::mpfr(1,precBits=precBits)
          for(t in 1:T) {
            if(Ni > K[i,t]) { next() }
            lit <- lit*(drbinomAPA(x    = Y[i,t],
                                   size = Ni*red[i,t],
                                   prob = pit,
                                   red  = red[i,t],
                                   precBits = precBits))
          }
          li <- li + lit*drpoisAPA(x = Ni, lambda = l_temp[i], red=red[i,1], precBits = precBits)
        }
        l <- l+log(li)
      }
    }
  }

  if(is.nan(l)) { l <- -Inf }

  if(VERBOSE) {
    if(typeof(l)=="double") {
      try(print(paste0("log likelihood: ", l)))
      try(print(paste0("parameters: ", par)))
    } else if(APA) {
      try(print(paste0("log likelihood: ", Rmpfr::formatMpfr(l))))
      try(print(paste0("parameters: ", Rmpfr::formatMpfr(par))))
    }
  }

  return(-1*l)
}
