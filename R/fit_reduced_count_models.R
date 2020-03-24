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
  } else {
    M <- outer(X = K1-1,Y = K1-1, FUN = tp_jk_V_APA, omeg = optimizeAPA::plogis_APA(sum(omeg * B_o)+sum(o_temp), precBits = precBits), gamm = exp(sum(B_g*gamm)+sum(g_temp)), red, precBits = precBits)
  }
  return(M)
}

#' @title        fit_red_Nmix_closed
#' @description  Find maximum likelihood estimates for closed population models, with model parameters log(lambda) and logit(pdet).
#' @param starts Vector of starting values for optimize. Has two elements, log(lambda) and logit(pdet).
#' @param nit    R by T matrix of full counts with R sites/rows and T sampling occassions/columns.
#' @param lambda_site_covariates Either NULL (no lambda site covariates) or a list of vectors of length R, where each vector represents one site covariate, and where the vector entries correspond to covariate values for each site. Note that the covariate structure is assumed to be log(lambda_it) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param pdet_site_covariates   Either NULL (no pdet site covariates) or a list of vectors of length R, where each vector represents one site covariate, and where the vector entries correspond to covariate values for each site. Note that the covariate structure is assumed to be logit(lambda_it) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param K      Upper bound on summations, either a single number, or a matrix of K values, one for each site and time (input the full counts value, eg if K=300 for full counts, K=reduction(300,red) for reduced counts).
#' @param red    reduction factor, either a number, or a vector of reduction factors (R sites reductions).
#' @param VERBOSE If true, prints the log likelihood to console at each optim iteration.
#' @param PARALLELIZE If true, calculation will be split over threads by sites. Will use as many threads as have been made available (initialize with START_PARALLEL(num_cores)).
#' @param APA    If true, will use arbitrary precision arithmetic in the likelihood calculations. Use precBits to specify the number of bits of precision.
#' @param precBits If APA=TRUE, then this will specify the number of bits of precision.
#' @param tolerance specifies tolerance for convergence (defulat is 10^-6), all components of estimated gradient must be less than tolerance for convergence. If APA=TRUE, then tolerance can be made very small (eg 10^-20) using: tolerance=Rmpfr::mpfr(10^-20, precBits=128). NOTE: currently tolerance is only used if method="DFP".#' @param method      Optimization method to use. Default is "DFP", and is the only method implemented for APA. When APA=FALSE, can use method="BFGS" to use optim().
#' @param outFile If not NULL, name of file for saving algorithm progress (overwritten at each iteration).
#' @param ...    Additional input for optim.
#' @examples
#' START_PARALLEL(num_cores=4)
#' Y    <- gen_Nmix_closed(8,8,250,0.5)
#' out  <- fit_red_Nmix_closed(Y$nit, red=10, K=300, starts = c(log(250),boot::logit(0.5)), PARALLELIZE=TRUE)
#' out2 <- fit_red_Nmix_closed(Y$nit, red=matrix(c(10,10,10,10,20,20,40,40),nrow=8, ncol=8), K=300, starts = c(log(250),boot::logit(0.5)), PARALLELIZE=TRUE)
#' END_PARALLEL()
#' @export
fit_red_Nmix_closed <- function(nit, lambda_site_covariates=NULL, pdet_site_covariates=NULL, red=1, K, starts=c(1,0), VERBOSE=FALSE, PARALLELIZE=FALSE, APA=FALSE, precBits=128, tolerance=10^-6, method="DFP", outFile=NULL, ...) {

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

  if(length(red)==1) {
    red <- matrix(red, nrow = nrow(Y_m), ncol=ncol(Y_m))
  }
  if(any(dim(red) != dim(Y_m))) {
    stop("reduction 'red' must be either one number, or a matrix with dimension nrow(nit) x ncol(nit)")
  }

  if(length(K)==1) {
    K <- matrix(K, nrow = nrow(Y_m), ncol=ncol(Y_m))
  }
  if(any(dim(K) != dim(Y_m))) {
    stop("K must be either one number, or a matrix with dimension nrow(nit) x ncol(nit)")
  }

  red_K  <- reduction(x = matrix(K, nrow=nrow(red), ncol=ncol(red)), red = red)

  Y_m <- reduction(nit, red)

  opt <- NULL
  if(!APA) {
    if(method=="DFP") {
      opt <- optimizeAPA::optim_DFP_NAPA(
                           starts      = starts,
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
    if(method != "DFP") { warning("USING DFP METHOD, only the DFP optimization method is currently implemented for APA")}
    opt <- optimizeAPA::optim_DFP_APA(
                         starts      = starts,
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
                         outFile     = outFile,
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

#' @title START_PARALLEL
#' @description Used to initialize parallel computing (useful for likelihood calculations which can be computationally intensive).
#' @param num_cores Number of cores to use for parallel processing.
#' @examples
#' START_PARALLEL(4)
#' Y    <- gen_Nmix_closed(16,8,50,0.5)
#' out  <- fit_red_Nmix_closed(Y$nit, red=10, K=100, starts = c(log(40),boot::logit(0.5)), PARALLELIZE=TRUE)
#' END_PARALLEL()
#'
#' out
#' @export
START_PARALLEL <- function(num_cores) {
  require(doParallel)
  require(foreach)
  doParallel::registerDoParallel(cores=num_cores)
}

#' @title END_PARALLEL
#' @description Used to end parallel computing.
#' @examples
#' START_PARALLEL(4)
#' Y    <- gen_Nmix_closed(16,8,50,0.5)
#' out  <- fit_red_Nmix_closed(Y$nit, red=10, K=100, starts = c(log(40),boot::logit(0.5)), PARALLELIZE=TRUE)
#' END_PARALLEL()
#'
#' out
#' @export
END_PARALLEL <- function() {
  doParallel::stopImplicitCluster()
}

#' @title        fit_red_Nmix_open
#' @description  Find maximum likelihood estimates for model parameters log(lambda), log(gamma), logit(omega), and logit(pdet).
#' @param starts Vector with four elements (if no covariates), log(lambda), log(gamma), logit(omega), and logit(pdet). When there are X covariates B_x for a parameter, will need X+1 starting values for that parameter (+1 for the constant term B0).
#' @param nit    R by T matrix of full counts with R sites/rows and T sampling occassions/columns.
#' @param lambda_site_covariates Either NULL (no lambda site covariates) or a list of vectors of length R, where each vector represents one site covariate, and where the vector entries correspond to covariate values for each site. Note that the covariate structure is assumed to be log(lambda_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param gamma_site_covariates  Either NULL (no gamma site covariates) or a list of vectors of length R, where each vector represents one site covariate, and where the vector entries correspond to covariate values for each site. Note that the covariate structure is assumed to be log(gamma_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param omega_site_covariates  Either NULL (no omega site covariates) or a list of vectors of length R, where each vector represents one site covariate, and where the vector entries correspond to covariate values for each site. Note that the covariate structure is assumed to be logit(omega_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param pdet_site_covariates   Either NULL (no pdet site covariates) or a list of vectors of length R, where each vector represents one site covariate, and where the vector entries correspond to covariate values for each site. Note that the covariate structure is assumed to be logit(pdet_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param gamma_time_covariates  Either NULL (no gamma time covariates) or a list of vectors of length T, where each vector represents one time covariate, and where the vector entries correspond to covariate values for each time. Note that the covariate structure is assumed to be log(gamma_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param omega_time_covariates  Either NULL (no omega time covariates) or a list of vectors of length T, where each vector represents one time covariate, and where the vector entries correspond to covariate values for each time. Note that the covariate structure is assumed to be logit(omega_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param pdet_time_covariates   Either NULL (no pdet time covariates) or a list of vectors of length T, where each vector represents one time covariate, and where the vector entries correspond to covariate values for each time. Note that the covariate structure is assumed to be logit(omega_i) = B0 + B1 &ast; V1_i + B2 &ast; V2_i + ...
#' @param K           Upper bound on summations, will be reduced by reduction factor red. Either a single number, or a vector of length R (specifying K for each site).
#' @param red         Reduction factor r, either a number or a vector of length R (specifying reduction factors for each site), or an R by T matrix of reduction factors for each observation.
#' @param VERBOSE     If TRUE, prints the log likelihood to console at each iteration of the optimization (unless optim is being used).
#' @param PARALLELIZE If TRUE, calculation will be split over threads by sites and times. This will not improve computation time if there are no site or time covariates. Will use as many threads as have been made available (initialize with START_PARALLEL(num_cores)).
#' @param APA         If TRUE, will use arbitrary precision arithmetic in the likelihood calculations. Use precBits to specify the number of bits of precision.
#' @param precBits    If APA=TRUE, then this will specify the number of bits of precision for arbitrary precision arithmetic.
#' @param tolerance   Specifies tolerance for convergence (defulat is 10^-6), all components of estimated gradient must be less than tolerance for convergence. If APA=TRUE, then tolerance can be made very small (eg 10^-20) using: tolerance=Rmpfr::mpfr(10^-20, precBits=128). NOTE: currently tolerance is only used if method="DFP".
#' @param method      Optimization method to use. Default is "DFP", and is the only method implemented for APA. When APA=FALSE, can use method="BFGS" to use optim().
#' @param outFile     If not NULL, name of file for saving algorithm progress (overwritten at each iteration).
#' @param ...         Additional input for optimization algorithm.
#' @examples
#' Y   <- gen_Nmix_open(num_sites = 5, num_times = 5, lambda = 20, pdet = 0.7, omega = 0.7, gamma = 2)
#' out <- fit_red_Nmix_open(nit = Y$nit, red = 3, K = 40, starts = c(0.5, 0.5, 0.5, 0.5))
#'
#' # lambda estimate:
#' exp(out$par[1])
#' # gamma estimate:
#' exp(out$par[2])
#' # omega estimate:
#' plogis(out$par[3])
#' # pdet estimate:
#' plogis(out$par[4])
#'
#'
#' # example with site covariates:
#' Y1 <- gen_Nmix_open(num_sites = 4, num_times = 5, lambda = 10, gamma = 5, omega = 0.50, pdet = 0.75)
#' Y2 <- gen_Nmix_open(num_sites = 4, num_times = 5, lambda = 5, gamma = 10, omega = 0.75, pdet = 0.50)
#' Y  <- rbind_pops(Y1, Y2)
#' START_PARALLEL(num_cores=4)
#' mod1 <- fit_red_Nmix_open(nit = Y$nit,
#'                           lambda_site_covariates = list(l1=c(0,0,0,0,1,1,1,1)),
#'                           gamma_site_covariates  = list(gs=c(0,0,0,0,1,1,1,1)),
#'                           gamma_time_covariates  = NULL,
#'                           omega_site_covariates  = list(os=c(0,0,0,0,1,1,1,1)),
#'                           omega_time_covariates  = NULL,
#'                           pdet_site_covariates   = list(ps=c(0,0,0,0,1,1,1,1)),
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
fit_red_Nmix_open <- function(nit, lambda_site_covariates=NULL, gamma_site_covariates=NULL, omega_site_covariates=NULL, pdet_site_covariates=NULL, gamma_time_covariates=NULL, omega_time_covariates=NULL, pdet_time_covariates=NULL, red=1, K, starts=NULL, VERBOSE=FALSE, PARALLELIZE=FALSE, APA=FALSE, precBits=128, tolerance=10^-6, method="DFP", outFile=NULL, ...) {
  if(length(red)==1 | length(red) == nrow(nit)) {
    red <- matrix(red, nrow = nrow(nit), ncol=ncol(nit))
  }
  if(any(dim(red) != dim(nit))) {
    stop("reduction 'red' must be either one number, or a vector of length R. Matrix K currently not supported.")
  }

  if(length(K)!=1 & length(K) != nrow(nit)) {
    stop("K must be one number, or a vector of length R. Matrix K currently not supported.")
  }

  K <- matrix(K, nrow = nrow(nit), ncol=ncol(nit))
  red_K  <- reduction(K, red)

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

   Y_m <- reduction(nit, red)

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
                                          outFile     = outFile,
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
     if(method != "DFP") { warning("USING DFP METHOD, only the DFP optimization method is currently implemented for APA")}
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
