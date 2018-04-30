#' Functions for estimating bedform migration rates in modern systems indirectly.
#'
#' \code{mahon_modern} calculates the slope of an ancient river, provided you know the flow depth and the grain size of the bed material. Implemented from Mahon et al. (2018)
#' @param Dbed The median grain size of bedload sediments in meters.
#' @param H The paleoflow depth in meters.
#' @return Returns a vector \code{S}, the dimensionless slope.
#' @export

## Models for bedform migration rates in modern systems
## Eric Barefoot
## May 2018

# TODO #
# document all the sub-functions

# This function provides functionality for calculating bedform migration rates using methods put out by Rob Mahon and Brandon McElroy in 2018.

# describe inputs and other stuff needed

# Mahon et al. use a Bayesian approach for paleohydraulic inversion, meaning that the coefficients that they describe for the empirical fit each have a distribution. They provide uncertainties for the distribution in their Table 1.
# In this way, we can ascribe ranges for the coefficients, and potentially provide a range of estimates of paleoslope accordingly. I think this means we can have confidence intervals, but I'm not sure about that.
# Initially, the user is allowed to select a percentile, and that percentile will be chosen for each parameter.

# In the future, users should be allowed to vary each parameter individually.

# get libraries

library(tibble)

## Set values and lookup table for empirical constants.

param_names = c('alpha0', 'alpha1', 'alpha2', 'beta0', 'beta1', 'gamma', 'delta0', 'delta1')

params_vals = c(-2.08, 0.254, -1.09, 0.6113, 1.305, 2.9, 0.0513, 0.7744)

params_uncert = c(0.036, 0.016, 0.044, 0.144, 0.0515, 0.7, 0.0017, 0.0123)

pars = tibble(name = param_names, vals = params_vals, uncert = params_uncert)

## Define functions from Mahon et al. 2018

# a function for table lookup. just feed it the name of the constant you want.
fun_paramVal = function(pars_list, param_name, zscore = 1) {

  therow = pars_list[which(pars_list$name == param_name), ]
  val = as.numeric(therow[2])
  dev = as.numeric(therow[3] * zscore)
  return(list(val = val, dev = dev))

}

# equation 1 from the paper.
bedform_bedload = function(hd, V, phi = 0.36) {

  qs = (1 - phi) * (hd * V) / 2

  return(qs)

}

# equation 2 from the paper.
bedform_velocity = function(beta0, beta1, S) {

  logV = beta0 + beta1 * log10(S)

  V = 10^(logV)

  return(V)

}

# equation 3 from the paper
height_length = function(delta0, delta1, L) {

  hd = delta0 * L ^ delta1

  return(hd)

}

# equation 4 from the paper
paleo_slope = function(alpha0, alpha1, alpha2, D50, Hbf) {

  logS = alpha0 + alpha1 * log10(D50) + alpha2 * log10(Hbf)

  S = 10^logS

  return(S)

}

# equation 5 from the paper
dune_height = function(gamma, Tsm) {

  hd = gamma * Tsm

  return(hd)

}

# self-contained function for calculating dune migration rates in the modern.
# uses measured dune heights
modern_1 = function(hd, S, uncert = NULL) { # if bedforme height is known

  V_est = bedform_velocity(beta0 = fun_paramVal(pars, 'beta0')$val, beta1 = fun_paramVal(pars, 'beta1')$val, S)
  qs_est = bedform_bedload(hd = hd, V = V_est)

  if (is.null(uncert)) {

    return(list(V = V_est, qs = qs_est))

  } else {

    zscore = abs(qnorm(uncert/2))
    env = c(-zscore, zscore)

    V_env = c(
      bedform_velocity(fun_paramVal(pars, 'beta0', zscore = env[1]), fun_paramVal(pars, 'beta1', zscore = env[1]), S),
      bedform_velocity(fun_paramVal(pars, 'beta0', zscore = env[2]), fun_paramVal(pars, 'beta1', zscore = env[2]), S)
    )

    qs_env = c(
      bedform_bedload(hd, V_env[1]),
      bedform_bedload(hd, V_env[2])
    )

    return(list(V = V_est, qs = qs_est, V_env = V_env, qs_env = qs_env))

  }

}

# self-contained function for calculating dune migration rates in the modern.
# uses measured dune lengths
modern_2 = function(L, S, uncert = NULL) { # if bedforme height is unknown, but length is known

  V_est = bedform_velocity(beta0 = fun_paramVal(pars, 'beta0')$val, beta1 = fun_paramVal(pars, 'beta1')$val, S = S)
  hd_est = height_length(delta0 = fun_paramVal(pars, 'delta0')$val, delta1 = fun_paramVal(pars, 'delta1')$val, L = L)
  qs_est = bedform_bedload(hd = hd_est, V = V_est)

  if (is.null(uncert)) {

    return(list(V = V_est, qs = qs_est))

  } else {

    zscore = abs(qnorm(uncert/2))
    env = c(-zscore, zscore)

    V_env = c(
      bedform_velocity(fun_paramVal(pars, 'beta0', zscore = env[1]), fun_paramVal(pars, 'beta1', zscore = env[1]), S),
      bedform_velocity(fun_paramVal(pars, 'beta0', zscore = env[2]), fun_paramVal(pars, 'beta1', zscore = env[2]), S)
    )

    qs_env = c(
      bedform_bedload(hd, V_env[1]),
      bedform_bedload(hd, V_env[2])
    )

    return(list(V = V_est, qs = qs_est, V_env = V_env, qs_env = qs_env))

  }

}
