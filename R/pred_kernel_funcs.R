#' Lognormal predation kernel
#' 
#' This is the most commonly-used predation kernel. The log of the predator/prey
#' mass ratio is normally distributed.
#' 
#' Writing the predator mass as \eqn{w} and the prey mass as \eqn{w_p},
#' the feeding kernel is given as
#' \deqn{\phi_i(w, w_p) = 
#' \exp \left[ \frac{-(\ln(w / w_p / \beta_i))^2}{2\sigma_i^2} \right]
#' }{\phi_i(w/w_p) = exp(-(ln(w/w_p/\beta_i))^2/(2\sigma_i^2))}
#' if \eqn{w/w_p} is larger than 1 and zero otherwise. Here \eqn{\beta_i} is the
#' preferred predator-prey mass ratio and \eqn{\sigma_i} determines the width of
#' the kernel. These two parameters need to be given in the species parameter
#' dataframe in the columns \code{beta} and \code{sigma}.
#' 
#' This function is called from [setPredKernel()] to set up the
#' predation kernel slots in a MizerParams object. 
#' 
#' @param ppmr A vector of predator/prey size ratios
#' @param beta The preferred predator/prey size ratio
#' @param sigma The width parameter of the log-normal kernel
#' 
#' @return A vector giving the value of the predation kernel at each of the
#'   predator/prey mass ratios in the \code{ppmr} argument.
#' @export
#' @family predation kernel
lognormal_pred_kernel <- function(ppmr, beta, sigma) {
    Beta <- log(beta)
    phi <- exp(-(log(ppmr) - Beta)^2 / (2 * sigma^2))
    return(phi)
}

#' Truncated lognormal predation kernel
#' 
#' This is like the [lognormal_pred_kernel()] but with an imposed maximum
#' predator/prey mass ratio
#' 
#' Writing the predator mass as \eqn{w} and the prey mass as \eqn{w_p},
#' the feeding kernel is given as
#' \deqn{\phi_i(w, w_p) = 
#' \exp \left[ \frac{-(\ln(w / w_p / \beta_i))^2}{2\sigma_i^2} \right]
#' }{\phi_i(w/w_p) = exp(-(ln(w/w_p/\beta_i))^2/(2\sigma_i^2))}
#' if \eqn{w/w_p} is between 1 and 
#' \eqn{\beta_i\exp(3\sigma_i)}{\beta_i exp(3\sigma_i)}
#' and zero otherwise. Here \eqn{\beta_i} is the preferred predator-prey mass
#' ratio and \eqn{\sigma_i} determines the width of the kernel. These two
#' parameters need to be given in the species parameter dataframe in the columns
#' `beta` and `sigma`.
#' 
#' This function is called from [setPredKernel()] to set up the
#' predation kernel slots in a MizerParams object. 
#' 
#' @param ppmr A vector of predator/prey size ratios
#' @param beta The preferred predator/prey size ratio
#' @param sigma The width parameter of the log-normal kernel
#' 
#' @return A vector giving the value of the predation kernel at each of the
#'   predator/prey mass ratios in the `ppmr` argument.
#' @export
#' @family predation kernel
truncated_lognormal_pred_kernel <- function(ppmr, beta, sigma) {
    Beta <- log(beta)
    phi <- exp(-(log(ppmr) - Beta)^2 / (2 * sigma^2))
    # rr is the maximal log predator/prey mass ratio
    rr <- exp(Beta + 3 * sigma)
    phi[ppmr > rr] <- 0
    return(phi)
}

#' Box predation kernel
#' 
#' A predation kernel where the predator/prey mass ratio is uniformly
#' distributed on an interval.
#' 
#' Writing the predator mass as \eqn{w} and the prey mass as \eqn{w_p}, the
#' feeding kernel is 1 if \eqn{w/w_p} is between `ppmr_min` and
#' `ppmr_max` and zero otherwise. The parameters need to be given in the
#' species parameter dataframe in the columns `ppmr_min` and
#' `ppmr_max`.
#' 
#' @param ppmr A vector of predator/prey size ratios
#' @param ppmr_min Minimum predator/prey mass ratio
#' @param ppmr_max Maximum predator/prey mass ratio
#' 
#' @return A vector giving the value of the predation kernel at each of the
#'   predator/prey mass ratios in the `ppmr` argument.
#' @export
#' @family predation kernel
box_pred_kernel <- function(ppmr, ppmr_min, ppmr_max) {
    assert_that(ppmr_min < ppmr_max)
    phi <- rep(1, length(ppmr))
    phi[ppmr > ppmr_max] <- 0
    phi[ppmr < ppmr_min] <- 0
    return(phi)
}

#' Power-law predation kernel
#' 
#' This predation kernel is a power-law, with sigmoidal cut-offs at large and
#' small predator/prey mass ratios.
#' 
#' The return value is calculated as
#' 
#' \code{
#' ppmr^kernel_exp /
#'   (1 + (exp(kernel_l_l) / ppmr)^kernel_u_l) /
#'   (1 + (ppmr / exp(kernel_l_r))^kernel_u_r) 
#' }
#'
#' The parameters need to be given as columns in the species parameter
#' dataframe.
#' 
#' @param ppmr A vector of predator/prey size ratios at which to evaluate the
#'   predation kernel.
#' @param kernel_exp The exponent of the power law
#' @param kernel_l_l The location of the left, rising sigmoid
#' @param kernel_u_l The shape of the left, rising sigmoid
#' @param kernel_l_r The location of the right, falling sigmoid
#' @param kernel_u_r The shape of the right, falling sigmoid
#' 
#' @return A vector giving the value of the predation kernel at each of the
#'   predator/prey mass ratios in the `ppmr` argument.
#' @export
#' @family predation kernel
power_law_pred_kernel <- function(ppmr, kernel_exp,
                                  kernel_l_l, kernel_u_l,
                                  kernel_l_r, kernel_u_r) {
    ppmr^kernel_exp /
        (1 + (exp(kernel_l_l) / ppmr)^kernel_u_l) /
        (1 + (ppmr / exp(kernel_l_r))^kernel_u_r) 
}
