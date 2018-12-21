#' Convert harmonic regression parameters to R*cos(2*pi*f*t+phi) form
#' 
#' @description 
#' The harmonic regression is usually fitted as 
#'    \eqn{z(t) = mu + A*cos(2*pi*f*t) + B*sin(2*pi*f*t)}
#' but is more meaningfully expressed as
#' \eqn{z(t) = mu + R*cos(2*pi*f*t + phi)},
#' where R and phi are the amplitude and phase.
#' This function converts the parameters in the R-phase form to the AB form.
#' 
#' @param RPhi Vector containing the R and phi coefficients respectively.
#'
#' @return vector c(A, B)
#'
#' @author Yuanhao Lai and A. I. McLeod
#' 
#' @examples 
#' #check by back-transform
#' harper:::RPhaseToAB(harper:::ABToRPhase(c(-4.205,1.075)))
#'                                                         
#' @keywords internal

RPhaseToAB <- function(RPhi) { 
 c(RPhi[1]*cos(2*pi*RPhi[2]), -RPhi[1]*sin(2*pi*RPhi[2]))
}