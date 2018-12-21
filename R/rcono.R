#' Random sample from contaminated normal distribution
#' 
#' @description Generates a random sample from a normal mixture 
#' \eqn{p*N(0,c^2)+(1-p)*N(0,1)}, where the densities are indicated and
#' \eqn{c} is the scale contamination \eqn{c>1} and \eqn{p} is the probability
#' of contamination which is chosen so both distributions contribute equally
#' to the variance. The variable is then standardized to have unit variance.
#' 
#' @param n Sample size.
#' @param c Scale inflation parameter.
#'
#' @return Vector of length n.
#' 
#' @details With probability \code{p} a normal random variable with mean 0
#' and standard deviation \code{c} is generated, where 
#' \eqn{p=1/(1+c^2))}, and with probability \eqn{1-p} a N(0,1) random
#' variable is generated. The resulting variable is standardized by dividing
#' by \eqn{\sqrt{1 - p + p*c^2}}, so the variance is 1.0. For a specified
#' parameter \code{c} the choice of \code{p} given by \eqn{p=1/(1+c^2))}
#' is determined by the requirement that both distributions in the mixture
#' contribute equally the variance and hence this is may be regarded as the
#' worst case.
#' 
#' @author A.I. McLeod
#' 
#' @references
#' Tukey, J. W. (1960).  A survey of sampling from contaminated distributions.  
#' In Contributions to Probability and Statistics: Essays in Honor of 
#' Harold Hotelling, Edited by I. Olkin, S. G. Ghurye, W. Hoeffding, 
#' W. G. Madow and H. B. Mann. Stanford University Press.
#' 
#' 
#' @examples
#' #there are 4 values from the contaminated distribution with this seed!
#'  set.seed(775511)
#'  qqnorm(harper:::rcono(100,5))
#' #
#' #Compute probability at least one value from the contaminated distribution
#' probAtLeastOne <- function(n, c) {
#'   p <- 1/(1+c^2)
#'   1-(1-p)^n 
#' }
#' probAtLeastOne(100, 5)
#' probAtLeastOne(10, 5)
#' 
#'                                                        
#' @keywords internal

rcono <- function(n, c) {
  stopifnot(c>1)
  p <- 1/(1+c^2)
  Z <- c*rbinom(n, size=1, prob=p)
  Z[Z==0] <- 1
  Z*rnorm(n)/sqrt(1-p+p*c^2)
}
