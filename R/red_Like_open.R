red_Like_open_unfixed <- function(par, omega=1,gamma=1, omeg_index, gamm_index, nit, l_s_c, g_s_c, g_t_c, o_s_c, o_t_c, p_s_c, p_t_c, K, red, VERBOSE=FALSE, PARALLELIZE=FALSE, APA=FALSE, precBits=128) {
   ll <- red_Like_open(par=par, nit=nit,
                      l_s_c=l_s_c, g_s_c=g_s_c,
                      g_t_c=g_t_c, o_s_c=o_s_c,
                      o_t_c=o_t_c, p_s_c=p_s_c,
                      p_t_c=p_t_c, K=K,
                      red=red, VERBOSE=VERBOSE,
                      PARALLELIZE=PARALLELIZE, APA=APA,
                      precBits=precBits)
  return(ll)
}

red_Like_open_fixed_omega <- function(par, omega=1,gamma, omeg_index, gamm_index, nit, l_s_c, g_s_c, g_t_c, o_s_c, o_t_c, p_s_c, p_t_c, K, red, VERBOSE=FALSE, PARALLELIZE=FALSE, APA=FALSE, precBits=128) {
  par2 = numeric(1+length(par))
  if(APA) {
    par2 = Rmpfr::mpfrArray(x = 0, precBits = precBits, dim = c(1,1+length(par)))
  }

  for(i in 1:length(par2)) {
    if(i == omeg_index) {
      par2[i] = omega
    } else {
      j=i
      if(j >= omeg_index) j = j-1
      par2[i] = par[j]
    }
  }

  ll <- red_Like_open(par=par2, nit=nit,
                      l_s_c=l_s_c, g_s_c=g_s_c,
                      g_t_c=g_t_c, o_s_c=o_s_c,
                      o_t_c=o_t_c, p_s_c=p_s_c,
                      p_t_c=p_t_c, K=K,
                      red=red, VERBOSE=VERBOSE,
                      PARALLELIZE=PARALLELIZE, APA=APA,
                      precBits=precBits)
  return(ll)
}



red_Like_open_fixed_gamma <- function(par, omega=1,gamma=0, omeg_index, gamm_index, nit, l_s_c, g_s_c, g_t_c, o_s_c, o_t_c, p_s_c, p_t_c, K, red, VERBOSE=FALSE, PARALLELIZE=FALSE, APA=FALSE, precBits=128) {
  par2 = numeric(1+length(par))
  if(APA) {
    par2 = Rmpfr::mpfrArray(x = 0, precBits = precBits, dim = c(1,1+length(par)))
  }

  for(i in 1:length(par2)) {
    if(i == gamm_index) {
      par2[i] = gamma
    } else {
      j=i
      if(j >= gamm_index) j = j-1
      par2[i] = par[j]
    }
  }

  ll <- red_Like_open(par=par2, nit=nit,
                      l_s_c=l_s_c, g_s_c=g_s_c,
                      g_t_c=g_t_c, o_s_c=o_s_c,
                      o_t_c=o_t_c, p_s_c=p_s_c,
                      p_t_c=p_t_c, K=K,
                      red=red, VERBOSE=VERBOSE,
                      PARALLELIZE=PARALLELIZE, APA=APA,
                      precBits=precBits)
  return(ll)
}


#' @title red_Like_open
#' @description   Used to calculate the negative of the log likelihood for open population models.
#' @param par     Vector with four elements (if no covariates), log(lambda), log(gamma), logis(omega), and logis(pdet). If there are covariates, include a starting value for each covariate. Order is: lambda, lambda site covariates, gamma, gamma site covariates, gamma time covariates, , omega site covariates, omega time covariates, pdet site covariates, pdet time covariates.
#' @param nit     R by T matrix of reduced counts with R sites/rows and T sampling occassions/columns.
#' @param l_s_c   list of lambda site covariates (list of vectors of length R (number of sites), each vector is a covariate).
#' @param g_s_c   list of gamma site covariates (list of vectors of length R (number of sites), each vector is a covariate)
#' @param g_t_c   list of gamma time covariates (list of vectors of length T (number of sampling occasions), each vector is a covariate)
#' @param o_s_c   list of omega site covariates (list of vectors of length R (number of sites), each vector is a covariate)
#' @param o_t_c   list of omega time covariates (list of vectors of length T (number of sampling occasions), each vector is a covariate)
#' @param p_s_c   list of pdet site covariates (list of vectors of length R (number of sites), each vector is a covariate)
#' @param p_t_c   list of pdet time covariates (list of vectors of length T (number of sampling occasions), each vector is a covariate)
#' @param K       Upper bound on summations (reduced counts upper bound). Currently only a vector of K values (one entry for each site) is possible, eg: K=rep(100,times=R).
#' @param red     Reduction factor
#' @param VERBOSE If TRUE, prints the log likelihood to console.
#' @param PARALLELIZE If TRUE, calculations will be done in parallel. Will use as many cores as have been made available (initialize with START_PARALLEL(num_cores)).
#' @param APA     Default is FALSE. If TRUE, will use arbitrary precision arithmetic in calculating the likelihood. Note that APA will be slower, however it is required for site population sizes larger than about 300. If APA = TRUE, then precBits specifies the number of bits of precision to use in the calculations.
#' @param precBits If APA=TRUE, then precBits specifies the number of bits of precision for arbitrary precision arithmetic.
#' @details Note that this function is adapted from the negative log likelihood function from the R package unmarked (Fiske and Chandler 2019), and uses the recursive method of computation described in Web Appendix A of Dail and Madsen 2011: Models for Estimating Abundance from Repeated Counts of an Open Metapopulation, published in Biometrics volume 67, issue 2.
#' @references Fiske, I., Chandler, R., Miller, D., Royle, A., Kery, M., Hostetler, J., Hutchinson, R., Smith, A., & Kellner, K. (2019). unmarked: Models for Data from Unmarked Animals (Version 0.13-1) [Computer software]. https://CRAN.R-project.org/package=unmarked
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
    g3        <- vector(length = R, mode = "list")
    g3_temp   <- list()

    if(is.null(g_t_c)) {
      t_gamt=matrix(0,nrow=T,ncol=1)
      t_B_gt=0
    } else {
      t_gamt=gamt
      t_B_gt=B_gt
    }
    if(is.null(o_t_c)) {
      t_omet=matrix(0,nrow=T,ncol=1)
      t_B_ot=0
    } else {
      t_omet=omet
      t_B_ot=B_ot
    }

    if(!PARALLELIZE) {
      for(i in 1:R) {
        # check if tpMAT is time dependent
        time_dependent = TRUE
        if(is.null(g_t_c) && is.null(o_t_c) && var(red[i,])==0) {
          time_dependent = FALSE
        }

        tempMat <- matrix(0, nrow = K[i]+1, ncol = K[i]+1)
        for(t in 1:T) {
          if(time_dependent | t == 1) {
            g3_temp <- tp_MAT(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=t_omet[t,], B_ot=t_B_ot, gamm=gamm[i,], B_g=B_g, gamt=t_gamt[t,], B_gt=t_B_gt, red = red[i,t], PARALLELIZE=FALSE)
            g3[[i]][[t]] <- g3_temp
          } else {
            g3[[i]][[t]] <- g3[[i]][[1]]
          }
        }
      }
    } else {
      g3_temp <- foreach(i=1:R, .packages = c("redNMix","foreach")) %dopar% {
        # check if tpMAT is time dependent
        time_dependent = TRUE
        if(is.null(g_t_c) && is.null(o_t_c) && var(red[i,])==0) {
          time_dependent = FALSE
        }

        g3_list = list()
        tempMat <- matrix(0, nrow = K[i]+1, ncol = K[i]+1)
        for(t in 1:T) {
          if(time_dependent | t == 1) {
            g3_temp <- tp_MAT(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=t_omet[t,], B_ot=t_B_ot, gamm=gamm[i,], B_g=B_g, gamt=t_gamt[t,], B_gt=t_B_gt, red = red[i,t], PARALLELIZE=FALSE)
            g3_list[[t]] <- g3_temp
          } else {
            g3_list[[t]] <- g3_list[[1]]
          }
        }
        return(g3_list)
      }
      for(i in 1:R) {
        g3[[i]]<- g3_temp[[i]]
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
        g1_t <- drbinom(x = Y[i,t], size = (0:K)*red[i,t], prob = plogis(sum(pdet[i,] * B_p)+sum(pt[t,] * B_pt)), red = red[i,t])
        g1_t_star <- g1_t * g_star

        # update g_star
        g_star = g3[[t]] %*% g1_t_star
      }

      # size takes possible values of N (0 to K) at time t==1
      g1 <- drbinom(x = Y[i,1], size = (0:K)*red[i,1], prob = plogis(sum(pdet[i,] * B_p)+sum(pt[t,] * B_pt)), red = red[i,1])
      g2 <- drpois(x = 0:K, lambda = exp(sum(lamb[i,] * B_l)), red = red[i,1])

      # apply recursive definition of likelihood
      return( log(sum(g1 * g2 * g_star)) )
    }, FUN.VALUE = numeric(1), K=K, T=T, Y=Y, lamb=lamb, B_l=B_l, pdet=pdet, B_p=B_p, pt=pt, B_pt=B_pt, red=red, g3=g3, g1_t_star=g1_t_star, g1_t=g1_t,g1=g1,g2=g2,g_star=g_star)

    ll <- sum(unlist(ll_i))
    ###########################
  } else { # APA = TRUE
    # extract lambda estimates from par, setup lambda covariate matrix lamb, and covariate vector B_l
    lamb <- matrix(rep(1,times=R),ncol=1)
    B_l <- Rmpfr::mpfr(par[1], precBits=precBits) # covariates for lambda
    if(!is.null(l_s_c)) {
      lamb <- cbind(lamb, do.call(cbind, l_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix
      B_l <- sapply(X = 1:(length(l_s_c)+1), FUN = function(X,par) { # one coeff per covariate +1 for baseline B0
        Rmpfr::mpfr(par[X],precBits=precBits)
      }, par=par)
      B_l <- methods::new("mpfr", unlist(B_l))
    }

    # extract gamma estimates from par, setup gamma covariate matrix gamm, and covariate vector B_g
    gamm <- matrix(rep(1,times=R),ncol=1)
    B_g <- Rmpfr::mpfr(par[length(B_l)+1], precBits=precBits) # covariates for gamma
    if(!is.null(g_s_c)) {
      gamm <- cbind(gamm, do.call(cbind, g_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix

      B_g <- sapply(X = 1:(length(g_s_c)+1), FUN = function(X,par) { # one coeff per covariate +1 for baseline B0
        Rmpfr::mpfr(par[length(B_l)+X],precBits=precBits)
      }, par=par)
      B_g <- methods::new("mpfr", unlist(B_g))
    }

    # extract gamma estimates from par, setup gamma covariate matrix gamt, and covariate vector B_gt
    gamt <- matrix(rep(0,times=T),ncol=1)
    B_gt <- NULL # covariates for gamma
    if(!is.null(g_t_c)) {
      gamt <- do.call(cbind, g_t_c) # rows are times, cols are covariate values, here we are creating the design matrix

      B_gt <- sapply(X = 1:length(g_t_c), FUN = function(X,par) { # one coeff per covariate
        Rmpfr::mpfr(par[length(B_l)+length(B_g)+X],precBits=precBits)
      }, par=par)
      B_gt <- methods::new("mpfr", unlist(B_gt))
    }

    # extract omega estimates from par, setup omega covariate matrix omeg, and covariate vector B_o
    omeg <- matrix(rep(1,times=R),ncol=1)
    B_o  <- Rmpfr::mpfr(par[length(B_l)+length(B_g)+length(B_gt)+1], precBits=precBits) # covariates for omega
    if(!is.null(o_s_c)) {
      omeg <- cbind(omeg, do.call(cbind, o_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix

      B_o <- sapply(X = 1:(length(o_s_c)+1), FUN = function(X,par) { # one coeff per covariate +1 for baseline B0
        Rmpfr::mpfr(par[length(B_l)+length(B_g)+length(B_gt)+X], precBits=precBits)
      }, par=par)
      B_o <- methods::new("mpfr", unlist(B_o))
    }

    # extract omega estimates from par, setup omega covariate matrix omet, and covariate vector B_ot
    omet <- matrix(rep(0,times=T),ncol=1)
    B_ot <- NULL # covariates for omega
    if(!is.null(o_t_c)) {
      omet <- do.call(cbind, o_t_c) # rows are times, cols are covariate values, here we are creating the design matrix

      B_ot <- sapply(X = 1:length(o_t_c), FUN = function(X,par) { # one coeff per covariate
        Rmpfr::mpfr(par[length(B_l)+length(B_g)+length(B_gt)+length(B_o)+X], precBits=precBits)
      }, par=par)
      B_ot <- methods::new("mpfr", unlist(B_ot))
    }

    # extract pdet estimates from par, setup pdet covariate matrix pdet, and covariate vector B_p
    pdet <- matrix(rep(1,times=R),ncol=1)
    B_p  <- Rmpfr::mpfr(par[length(B_l)+length(B_g)+length(B_gt)+length(B_o)+length(B_ot)+1], precBits=precBits) # covariates for pdet
    if(!is.null(p_s_c)) {
      pdet <- cbind(pdet, do.call(cbind, p_s_c)) # rows are sites, cols are covariate values, here we are creating the design matrix

      B_p <- sapply(X = 1:(length(p_s_c)+1), FUN = function(X,par) { # one coeff per covariate +1 for baseline B0
        Rmpfr::mpfr(par[length(B_l)+length(B_g)+length(B_gt)+length(B_o)+length(B_ot)+X], precBits=precBits)
      }, par=par)
      B_p <- methods::new("mpfr", unlist(B_p))
    }

    # extract pdet estimates from par, setup pdet covariate matrix pt, and covariate vector B_pt
    pt   <- matrix(rep(0,times=T),ncol=1)
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
    g3        <- vector(length = R, mode = "list")
    g3_temp   <- list()

    if(is.null(g_t_c)) {
      t_gamt=matrix(0,nrow=T,ncol=1)
      t_B_gt=0
    } else {
      t_gamt=gamt
      t_B_gt=B_gt
    }
    if(is.null(o_t_c)) {
      t_omet=matrix(0,nrow=T,ncol=1)
      t_B_ot=0
    } else {
      t_omet=omet
      t_B_ot=B_ot
    }

    if(!PARALLELIZE) {
      for(i in 1:R) {
        # check if tpMAT is time dependent
        time_dependent = TRUE
        if(is.null(g_t_c) && is.null(o_t_c) && var(red[i,])==0) {
          time_dependent = FALSE
        }

        tempMat <- Rmpfr::mpfrArray(0, precBits = precBits, dim = c(K[i]+1,K[i]+1))
        for(t in 1:T) {
          if(time_dependent | t == 1) {
            g3_temp <- tp_MAT_APA(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=t_omet[t,], B_ot=t_B_ot, gamm=gamm[i,], B_g=B_g, gamt=t_gamt[t,], B_gt=t_B_gt, red = red[i,t], PARALLELIZE=FALSE, precBits = precBits)
            g3[[i]][[t]] <- g3_temp
          } else {
            g3[[i]][[t]] <- g3[[i]][[1]]
          }
        }
      }
    } else { # Run in Parallel
      g3_temp <- foreach(i=1:R, .packages = c("redNMix","foreach")) %dopar% {
        # check if tpMAT is time dependent
        time_dependent = TRUE
        if(is.null(g_t_c) && is.null(o_t_c) && var(red[i,])==0) {
          time_dependent = FALSE
        }

        g3_list = list()
        tempMat <- Rmpfr::mpfrArray(0, precBits = precBits, dim = c(K[i]+1,K[i]+1))
        for(t in 1:T) {
          if(time_dependent | t == 1) {
            g3_temp <- tp_MAT_APA(M = tempMat, omeg = omeg[i,], B_o=B_o, omet=t_omet[t,], B_ot=t_B_ot, gamm=gamm[i,], B_g=B_g, gamt=t_gamt[t,], B_gt=t_B_gt, red = red[i,t], PARALLELIZE=FALSE, precBits = precBits)
            g3_list[[t]] <- g3_temp
          } else {
            g3_list[[t]] <- g3_list[[1]]
          }
        }
        return(g3_list)
      }
      for(i in 1:R) {
        g3[[i]]<- g3_temp[[i]]
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
        g1_t <- drbinomAPA(x = Y[i,t], size = (0:K)*red[i,t], prob = p, red = red[i,t], precBits = precBits)
        g1_t_star <- g1_t * g_star

        # update g_star
        g_star = g3[[t]] %*% g1_t_star
      }

      Bp_temp  <- Rmpfr::mpfr2array(ifelse(is.null(B_p), Rmpfr::mpfr(0,precBits), sum(pdet[i,] * B_p)),dim = c(1,1))[1]
      Bpt_temp <- Rmpfr::mpfr2array(ifelse(is.null(B_pt), Rmpfr::mpfr(0,precBits), sum(pt[t,] * B_pt)),dim = c(1,1))[1]

      # size takes possible values of N (0 to K) at time t==1
      p <- optimizeAPA::plogis_APA(Bp_temp+Bpt_temp, precBits = precBits)
      g1 <- drbinomAPA(x = Y[i,1], size = (0:K)*red[i,1], prob = p, red = red[i,1], precBits = precBits)
      g2 <- drpoisAPA(x = 0:K, lambda = exp(sum(lamb[i,] * B_l)), red = red[i,1], precBits = precBits)

      # apply recursive definition of likelihood
      return( Rmpfr::mpfr(log(sum(g1 * g2 * g_star)), precBits) ) # + 1e-320
    }, K=K, T=T, Y=Y, lamb=lamb, B_l=B_l, pdet=pdet, B_p=B_p, pt=pt, B_pt=B_pt, red=red, g3=g3, g1_t_star=g1_t_star, g1_t=g1_t,g1=g1,g2=g2,g_star=g_star,precBits=precBits)

    ll <- sum(ll_i)
  }

  if(is.nan(ll)) { ll <- -Inf }

  if(VERBOSE) {
    if(typeof(ll)=="double") {
      try(print(paste0("log likelihood: ", ll)))
      try(print(paste0("parameters: ", par)))
    } else if(APA) {
      try(print(paste0("log likelihood: ", Rmpfr::formatMpfr(ll))))
      try(print(paste0("parameters: ", Rmpfr::formatMpfr(par))))
    }
  }

  return(-1*ll)
}
