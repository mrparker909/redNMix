
#' @title gen_Nmix_open
#' @description Generate a population/observation pair with the structure of an open N-mixture model. Use cbind_pop() and rbind_pop() functions to combine sub-pops to easily add covariate structures.
#' @param num_sites The number of observation sites.
#' @param num_times The number of sampling occasions.
#' @param lambda    The population rate parameter (\eqn{N_i \sim Poisson(\lambda)}).
#' @param pdet      The probability of detection \eqn{p} (\eqn{n_{it} \sim Binomial(N_i,p)}).
#' @param omega     The probability of survival (\eqn{S_{it} \sim Binomial(N_{it}, \omega)}). Either a single number, or a vector of length t-1 (omega for each sampling occassion after the first).
#' @param gamma     The recruitment rate parameter (\eqn{G_{it} \sim Poisson(\gamma)}). Either a single number, or a vector of length t-1 (gamma for each sampling occassion after the first).
#' @param starts    If not NULL, a vector of starting population values (used instead of lambda).
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
#'
#' # example with time covariate for omega:
#' pop_a <- gen_Nmix_open(2,4,50,0.75,0.9,1) # survival=90%
#' pop_b <- gen_Nmix_open(2,4,50,0.75,0.1,1, starts=pop_a$Ni[,4]) # survival=10%, use last time observation from pop_a to start pop_b, then drop first time observation from pop_b (this makes a seamless transition between pop_a and pop_b)
#' pop_b <- list(Ni=pop_b$Ni[,-1], nit=pop_b$nit[,-1])
#' pop   <- cbind_pops(pop_a,pop_b) # combine to produce one population list
#' pop
#' @export
gen_Nmix_open <- function(num_sites,num_times,lambda,pdet,omega,gamma, starts=NULL) {
  U    = num_sites
  T    = num_times
  lamb = lambda
  gamm = if(length(gamma)==1) {gamm<-rep(gamma, times=T-1) } else { gamm<-gamma }
  omeg = if(length(omega)==1) {omeg<-rep(omega, times=T-1) } else { gamm<-gamma }

  Ntemp <- c(rep(rpois(n=U,lambda = lamb),times=T))
  if(!is.null(starts)) {
    Ntemp = starts
    if(length(starts) != U) { stop("length of starts must equal num_sites") }
  }

  Ni <- matrix(data=Ntemp, nrow = U, ncol = T)
  for(i in 2:T) {
    Ni[,i] <- rbinom(n = U, size = Ni[,i-1], prob = omeg[i-1]) + rpois(n = U, lambda = gamm[i-1])
  }

  nit <- Ni
  nit[] <- vapply(Ni, function(x) { rbinom(size = x, n = 1, prob = pdet) }, numeric(1))

  return(list(Nit=Ni, nit=nit))
}



#' @title gen_Nmix_closed
#' @description Generate a population/observation pair with the structure of a closed N-mixture model.
#' @param num_sites The number of observation sites.
#' @param num_times The number of sampling occasions.
#' @param lambda    The population rate parameter (\eqn{N_i \sim Poisson(\lambda)})
#' @param pdet      The probability of detection \eqn{p} (\eqn{n_{it} \sim Binomial(N_i,p)})
#' @return A list object with two named matrices. Ni contains the total population per site
#'         (each row represents a site, each column a sampling occasion). nit contains the observed
#'         counts (rows=sites, columns=sampling occasions).
#' @examples
#' pop <- gen_Nmix_closed(num_sites = 5, num_times = 10, lambda = 50, pdet = 0.4)
#' pop
#'
#' # population with a site covariate for lambda, B1=c(0,0,0,1,1):
#' pop_a <- gen_Nmix_closed(num_sites = 3, num_times = 4, lambda = 10, pdet = 0.4)
#' pop_b <- gen_Nmix_closed(num_sites = 2, num_times = 4, lambda = 40, pdet = 0.4)
#' pop   <- rbind_pops(pop_a,pop_b)
#' pop
#'
#' # population with a time covariate for lambda, Bt1=c(0,0,0,0,1,1):
#' pop_a <- gen_Nmix_closed(num_sites = 2, num_times = 4, lambda = 10, pdet = 0.20)
#' pop_b <- gen_Nmix_closed(num_sites = 2, num_times = 2, lambda = 10, pdet = 0.85)
#' pop   <- cbind_pops(pop_a,pop_b)
#' pop
#'
#' @export
gen_Nmix_closed <- function(num_sites,num_times,lambda,pdet) {
  Ntemp <- c(rep(rpois(n=num_sites,lambda = lambda),times=num_times))
  Ni    <- matrix(data=Ntemp, nrow = num_sites, ncol = num_times)

  nit   <- Ni
  nit[] <- vapply(Ni, function(x) { rbinom(size = x, n = 1, prob = pdet) }, numeric(1))

  return(list(Ni=Ni, nit=nit))
}

#' @title rbind_pops
#' @description rbind two population lists (as generated from either gen_Nmix_closed or gen_Nmix_open)
#' @param pop1 first population list, must have named matrices Ni and nit as list elements.
#' @param pop2 second population list, must have named matrices Ni and nit with same number of columns (sampling occasions) as pop1.
#' @return Single population list with sites from pop1 and pop2.
rbind_pops <- function(pop1, pop2) {
  return(list(Ni = rbind(pop1$Ni, pop2$Ni), nit = rbind(pop1$nit, pop2$nit)))
}

#' @title cbind_pops
#' @description cbind two population lists (as generated from either gen_Nmix_closed or gen_Nmix_open)
#' @param pop1 first population list, must have named matrices Ni and nit as list elements.
#' @param pop2 second population list, must have named matrices Ni and nit with same number of rows (sampling sites) as pop1.
#' @return Single population list with sites from pop1 and pop2.
cbind_pops <- function(pop1, pop2) {
  return(list(Ni = cbind(pop1$Ni, pop2$Ni), nit = cbind(pop1$nit, pop2$nit)))
}
