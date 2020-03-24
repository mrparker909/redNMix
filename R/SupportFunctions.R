#' @title LR
#' @description Calculate maximum loss of precision from considering reduced counts over full counts.
#' @param fullCounts An R by T matrix of full counts
#' @param r The reduction factor to reduce the full counts matrix by
#' @examples
#' nit = matrix(c(114,117,122,131,112,100,114,138,113,114,
#'                116, 96,105,101,109,107,101,111,102, 98,
#'                 99, 90, 92, 93,103, 82,100, 78,100, 81,
#'                 92,130,112,110,114,114,113, 99,109,108,
#'                 99,102,110,105,106,103,101,104, 92,115),
#'              nrow=5, ncol=10, byrow=T)
#' LR(nit, r=10)
LR <- function(fullCounts, r=1) {
  rCounts = reduction(x = fullCounts, red = r)

  LR = max(abs(fullCounts-r*rCounts))
  return(LR)
}

#' @title mcd
#' @description Find the Most Common Divisor of a data set
#' @param nit Full counts
#' @param howmany Maximum number of common divisors to return.
#' @examples
#' mcd(c(10,20,30,35,12,18), howmany = 5)
#'
#' @export
mcd <- function(nit, howmany=10) {
  if(!require(numbers)) { stop("requires package \"numbers\".") }
  x <- list()
  for(n in nit) {
    x <- c(x,numbers::divisors(n))
  }
  y <- unlist(x)
  howmany <- min(howmany, length(table(y)))
  sort(table(y),decreasing=TRUE)[1:howmany]
}

#' @title round2
#' @description rounding in R rounds to nearest even number (eg 0.5 rounds to 0), this function does "normal" rounding (eg 0.5 rounds to 1)
#' @param x number to be rounded
#' @param n decimal place to round to (n=0 rounds to nearest whole number, n=1 rounds to the nearest tenth, n=-1 rounds to the nearest tens digit, etc).
round2 <- function(x, n=0) {
  z <- sign(x)*trunc(abs(x)*10^n + 0.5)/10^n
  return(z)
}

#' @title reduction
#' @description Reduce an integer value using a reduction function \eqn{R(x;r)}.
#' @param x   The integer which is to be reduced.
#' @param red The factor r by which to reduce the input x.
#' @param FUN The reduction function (default is round2(), some alternatives are ceiling() and floor())
#' @return An integer value which is the reduction of integer x by reduction factor red using function FUN.
#' @examples
#' x <- 104
#' xr1 <- reduction(x, 10)
#' xr2 <- reduction(x, 10, ceiling)
#'
#' x <- matrix(rbinom(9, 51, 0.25), nrow=3, ncol=3)
#' reduction(x, 5, FUN=ceiling)
#' reduction(x, 5)
#' reduction(x, 5, FUN=floor)
#' @export
reduction <- function(x, red, FUN=round2) {
  FUN(x/red)
}

