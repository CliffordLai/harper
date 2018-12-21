#' Convert harmonic regression parameters to A*cos+B*sin form
#' 
#' @description 
#' The harmonic regression is usually fitted as 
#'    \eqn{z(t) = mu + A*cos(2*pi*f*t) + B*sin(2*pi*f*t)}
#' but is more meaningfully expressed as
#' \eqn{z(t) = mu + R*cos(2*pi*f*t + 2*pi*phi)},
#' where R and \eqn{-0.5<phi=0.5} are the amplitude and phase.
#' This function converts the parameters to the amplitude-phase 
#' parameterization.
#' 
#' @param AB Vector containing the A and B coefficients respectively.
#'
#' @return vector c(R, phi)
#'
#' @author Yuanhao Lai and A. I. McLeod
#' 
#' @examples 
#' #Bloomfield ozone example, p.17
#' harper:::ABToRPhase(c(-4.205,1.075))
#'                                                         
#' @keywords internal

ABToRPhase <- function(AB) { #Bloomfield, p.13
 R <- sqrt(sum(AB^2))
 if (AB[1]>0)          return(c(R, atan(-AB[2]/AB[1])/(2*pi)))
 if (AB[1]<0&AB[2]>0)  return(c(R, atan(-AB[2]/AB[1])/(2*pi) - 0.5))
 if (AB[1]<0&AB[2]<=0) return(c(R, atan(-AB[2]/AB[1])/(2*pi) + 0.5))
 if (AB[1]==0&AB[2]>0) return(c(R, -0.25))
 if (AB[1]==0&AB[2]<0) return(c(R, 0.25))
 if (all(AB==0))       return(c(R, 0))
}