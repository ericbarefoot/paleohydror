#' Global empirical relation method for paleoslope inversion
#'
#' \code{trampush_slp} calculates the slope of an ancient river, provided you know the flow depth and the grain size of the bed material. Implemented from Trampush et al. (2014)
#' @param Dbed The median grain size of bedload sediments in meters.
#' @param H The paleoflow depth in meters.
#' @return Returns a vector \code{S}, the dimensionless slope.
#' @export

## Inversion for Paleoslope based on Trampush 2014
## Eric Barefoot
## April 2018

# This function provides functionality for computing paleoslope based on a method described by Sheila Trampush in a 2014 paper.

# describe inputs and other stuff needed

# Trampush et al. use a Bayesian approach for paleohydraulic inversion, meaning that the coefficients that they describe for the empirical fit each have a distribution. They provide percentiles for the distribution in their Table 1.
# In this way, we can ascribe ranges for the coefficients, and potentially provide a range of estimates of paleoslope accordingly. I think this means we can have confidence intervals, but I'm not sure about that.
# Initially, the user is allowed to select a percentile, and that percentile will be chosen for each parameter.

# In the future, users should be allowed to vary each parameter individually.

trampush_slp = function(Dbed, H, perc = 50) {

  a0sp = smooth.spline(x = c(-2.14, -2.10, -2.08, -2.05, -2.01), y = c(2.5,25,50,75,97.5))
  a1sp = smooth.spline(x = c(0.222, 0.244, 0.254, 0.266, 0.287), y = c(2.5,25,50,75,97.5))
  a1sp = smooth.spline(x = c(-1.18, -1.12, -1.09, -1.06, -1.00), y = c(2.5,25,50,75,97.5))

  logS = predict(a0sp, perc) + predict(a1sp, perc) * log(Dbed) + predict(a2sp, perc) * log(H)

  return(S = 10^logS)

}
