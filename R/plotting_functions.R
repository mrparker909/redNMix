
#' @title             plot_red_like_closed_lambda
#' @description       Plot likelihood given a pdet and range for lambda.
#' @param nit         R by T matrix of reduced counts with R sites/rows and T sampling occassions/columns.
#' @param startLambda Starting value for lambda.
#' @param endLambda   Ending value for lambda.
#' @param stepsize    Spacing between values of lambda
#' @param pdet        Probability of detection (pdet = 0.5 means 50\% chance of detection)
#' @param red         Reduction factor.
#' @param K           Upper bound on summations (reduced counts upper bound).
#' @examples
#' Y <- gen_Nmix_closed(5,5,250,0.5)
#' plot_red_like_closed_lambda(nit=reduction(Y$nit,10), startLambda = 150, endLambda = 350, stepsize=10, pdet = 0.5, red = 10, K = reduction(400,10))
#' @export
plot_red_like_closed_lambda <- function(nit, startLambda, endLambda, stepsize, pdet, red, K) {

  parB <- startLambda
  parT <- endLambda
  j <- 1
  L <- NULL
  for(parM in seq(parB,parT,stepsize)) {
    L[j] <- -1*red_Like_closed(nit = nit, par = c(log(parM), boot::logit(pdet)), K = matrix(K, nrow=nrow(nit),ncol=ncol(nit)), red = matrix(red, nrow=nrow(nit),ncol=ncol(nit)), l_s_c = NULL, p_s_c = NULL)
    j <- j+1
  }

  plot(y=L, x=seq(parB,parT,stepsize), xlab="lambda")
}

#' @title             plot_red_like_closed_pdet
#' @description       Plot likelihood given a lambda and range for pdet.
#' @param nit         R by T matrix of reduced counts with R sites/rows and T sampling occassions/columns.
#' @param startPdet   Starting value for pdet.
#' @param endPdet     Ending value for pdet.
#' @param stepsize    Spacing between values of pdet
#' @param lambda      Initial abundance parameter.
#' @param K           Upper bound on summations (reduced counts upper bound).
#' @param red         Reduction factor r.
#' @examples
#' Y <- gen_Nmix_closed(5,5,250,0.5)
#' plot_red_like_closed_pdet(nit=reduction(Y$nit,10), startPdet = 0.1, endPdet = 1.0, stepsize=0.1, lambda = 250, red = 10, K = reduction(400,10))
#' @export
plot_red_like_closed_pdet <- function(nit, startPdet, endPdet, stepsize, lambda, red, K) {

  parB <- startPdet
  parT <- endPdet
  j <- 1
  L <- NULL
  for(parM in seq(parB,parT,stepsize)) {
    L[j] <- -1*red_Like_closed(nit = nit, par = c(log(lambda), boot::logit(parM)), K = matrix(K, nrow=nrow(nit),ncol=ncol(nit)), red = matrix(red, nrow=nrow(nit),ncol=ncol(nit)), l_s_c = NULL, p_s_c = NULL)
    j <- j+1
  }

  plot(y=L, x=seq(parB,parT,stepsize), xlab="pdet")
}

#' @title                 plot_2d_red_like_closed
#' @description           Plot 2D closed likelihood given a range for pdet and for lambda.
#' @param nit             R by T matrix of reduced counts with R sites/rows and T sampling occassions/columns.
#' @param startPdet       Starting value for pdet.
#' @param endPdet         Ending value for pdet.
#' @param stepsizePdet    Spacing between values of pdet
#' @param startLambda     Starting value for lambda.
#' @param endLambda       Ending value for lambda.
#' @param stepsizeLambda  Spacing between values of lambda
#' @param K               Upper bound on summations (reduced counts upper bound).
#' @param red             Reduction factor r.
#' @return                Returns a matrix of likelihoods where rows represent pdet, and columns represent lambda.
#' @examples
#' set.seed(12345)
#' Y <- gen_Nmix_closed(15,15,150,0.5)
#' L <- plot_2d_red_like_closed(reduction(Y$nit,10), 0.05, 1.0, 0.05, 10, 250, 10, 10, reduction(400,10))
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
