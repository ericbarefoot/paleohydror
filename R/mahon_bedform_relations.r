## Inversion for bedform migration rates and sediment flux in the ancient
## Eric Barefoot
## May 2018

## TODO
# clean up documentation more.
# complete the uncertainty routines.
# add the ability to also account for variance in data.

# This collection of functions provides functionality for calculating bedform migration rates using methods put out by Rob Mahon and Brandon McElroy in 2018

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
fun_paramVal = function(pars_list, param_name, zscore = 0) {

  therow = pars_list[which(pars_list$name == param_name), ] # lookup the constant

  val = as.numeric(therow[2]) # extract the value
  dev = as.numeric(therow[3] * zscore) # extract the magnitude of deviation from mean at zscore.

  return(list(val = val, dev = dev))

}

# equation 1 from the paper.

#' Bedform-bedload equation
#'
#' \code{bedform_bedload} calculates the sediment flux through bedload in a river based on the height and migration rate of bedforms. From Simons et al. (1965), Implemented from Mahon et al. (2018)
#' @param hd The mean height of bedforms in meters.
#' @param V The dune migration rate in meters per second.
#' @param phi The porosity of the riverbed material. Default is 0.36.
#' @return Returns a vector \code{qs}, the sediment flux in m^3/s.
#' @export

bedform_bedload = function(hd, V, phi = 0.36) {

  qs = (1 - phi) * (hd * V) / 2

  return(qs)

}

# equation 2 from the paper.

#' Bedform velocity relation
#'
#' \code{bedform_velocity} estimates the velocity of a migrating bedform using a power law fit for slope. Implemented from Mahon et al. (2018)
#' @param S The dimensionless river slope.
#' @param beta0 Empirical constant. Value and associated uncertainty from Table 1 in Mahon + McElroy 2018.
#' @param beta1 Empirical constant. Value and associated uncertainty from Table 1 in Mahon + McElroy 2018.
#' @return Returns a vector \code{V}, the bedform migration rate in m/s.
#' @export

bedform_velocity = function(beta0, beta1, S) {

  logV = beta0 + beta1 * log10(S)

  V = 10^(logV)

  return(V)

}

# equation 3 from the paper

#' Bedform length to height relation
#'
#' \code{height_length} estimates the mean height of bedforms based on the mean bedform length. Implemented from Mahon et al. (2018)
#' @param L The mean bedform length.
#' @param delta0 Empirical constant. Value and associated uncertainty from Table 1 in Mahon + McElroy 2018.
#' @param delta1 Empirical constant. Value and associated uncertainty from Table 1 in Mahon + McElroy 2018.
#' @return Returns a vector \code{hd}, the bedform migration rate in m/s.
#' @export

height_length = function(delta0, delta1, L) {

  hd = delta0 * L ^ delta1

  return(hd)

}

# equation 4 from the paper

#' Paleoslope relation
#'
#' \code{paleo_slope} estimates the slope of an ancient river. Originally developed in Trmapush et al. 2014. Implemented from Mahon et al. (2018)
#' @param D50 The median bedload grain size. Determined from samples.
#' @param Hbf The mean bankfull flow depth. Determined from preserved barforms.
#' @param alpha0 Empirical constant. Value and associated uncertainty from Table 1 in Mahon + McElroy 2018.
#' @param alpha1 Empirical constant. Value and associated uncertainty from Table 1 in Mahon + McElroy 2018.
#' @param alpha2 Empirical constant. Value and associated uncertainty from Table 1 in Mahon + McElroy 2018.
#' @return Returns a vector \code{S}, the dimensionless paleoslope of the ancient river.
#' @export

paleo_slope = function(alpha0, alpha1, alpha2, D50, Hbf) {

  logS = alpha0 + alpha1 * log10(D50) + alpha2 * log10(Hbf)

  S = 10^logS

  return(S)

}

# equation 5 from the paper

#' Dune height from cosets relation
#'
#' \code{dune_height} estimates the height of dunes based on cross-stratification. Originally developed by Leclair and Bridge (2001). Implemented from Mahon et al. (2018)
#' @param Tsm The mean thickness of cross-bedded cosets. (see publication for more information)
#' @param gamma Empirical constant. Value and associated uncertainty from Table 1 in Mahon + McElroy 2018.
#' @return Returns a vector \code{hd}, the dune height in the ancient river.
#' @export

dune_height = function(gamma, Tsm) {

  hd = gamma * Tsm

  return(hd)

}

# self-contained function for calculating sediment flux rates from outcrops.
# have to have grain size measurements and bankfull depth estimates

#' Estimating sediment flux in ancient rivers.
#'
#' \code{ancient_1} estimates the sediment flux (\code{qs}) of ancient river deposits based on dune cross-stratification, bankfull depth reconstruction, and grain size. Implemented from Mahon et al. (2018). Original relations referenced therein.
#' @param Tsm The mean thickness of cross-bedded cosets. (see publication for more information)
#' @param D50 The median bedload grain size. Determined from samples.
#' @param Hbf The mean bankfull flow depth. Determined from preserved barforms.
#' @param uncert The uncertainty window for empirical parameters. The user should put in a confidence interval as a decimal. A lower confidence interval will have a narrower error envelope, but will require more tenuous conclusion.
#' @return Returns a list with elements:
#' \code{S}, the river slope (dimensionless).
#' \code{hd}, the dune height in the ancient river (m).
#' \code{V}, the dune migration rate in the ancient river (m/s).
#' \code{qs}, the bedload flux in the ancient river (m^3/s).
#' If uncertainty is also requested, the output will include upper and lower bounds for confidence, \code{V_env} and \code{qs_env}
#' @export

ancient_1 = function(Tsm, D50, Hbf, uncert = NULL) {

  S_est = paleo_slope(
    alpha0 = fun_paramVal(pars, 'alpha0')$val,
    alpha1 = fun_paramVal(pars, 'alpha1')$val,
    alpha2 = fun_paramVal(pars, 'alpha2')$val,
    D50 = D50,
    Hbf = Hbf
  )
  hd_est = dune_height(gamma = fun_paramVal(pars, 'gamma')$val, Tsm)
  V_est = bedform_velocity(beta0 = fun_paramVal(pars, 'beta0')$val, beta1 = fun_paramVal(pars, 'beta1')$val, S_est)
  qs_est = bedform_bedload(hd = hd_est, V = V_est)

  if (is.null(uncert)) {

    return(list(S = S_est, hd = hd_est, V = V_est, qs = qs_est))

  } else {

    zscore = abs(qnorm(uncert/2))
    env = c(-zscore, zscore)

    V_env = c(
      bedform_velocity(beta0 = fun_paramVal(pars, 'beta0', zscore = env[1])$val, beta1 = fun_paramVal(pars, 'beta1', zscore = env[1])$val, S = S_est),
      bedform_velocity(beta0 = fun_paramVal(pars, 'beta0', zscore = env[2])$val, beta1 = fun_paramVal(pars, 'beta1', zscore = env[2])$val, S = S_est)
      # bedform_velocity(fun_paramVal(pars, 'beta0', zscore = env[2]), fun_paramVal(pars, 'beta1', zscore = env[2]), S_est)
    )

    qs_env = c(
      bedform_bedload(hd_est, V_env[1]),
      bedform_bedload(hd_est, V_env[2])
    )

    return(list(S = S_est, V = V_est, qs = qs_est, V_env = V_env, qs_env = qs_env))

  }

}

# self-contained function for calculating dune migration rates in the modern.
# uses measured dune heights

#' Estimating sediment flux in modern rivers if bedform height is known.
#'
#' \code{modern_1} estimates the sediment flux (\code{qs}) of modern rivers based on bedform height and river slope. Implemented from Mahon et al. (2018).
#' @param hd The mean height of migrating bedforms.
#' @param S The reach-average slope of the river.
#' @param uncert The uncertainty window for empirical parameters. The user should put in a confidence interval. A lower confidence interval will have a narrower error envelope, but will require more tenuous conclusions.
#' @return Returns a list with elements:
#' \code{V}, the dune migration rate in the river (m/s).
#' \code{qs}, the bedload flux in the river (m^3/s).
#' If uncertainty is also requested, the output will include upper and lower bounds for confidence, \code{V_env} and \code{qs_env}
#' @export

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
#' Estimating sediment flux in modern rivers if bedform height is unknown, but bedform length is known.
#'
#' \code{modern_2} estimates the sediment flux (\code{qs}) of modern rivers based on bedform length and river slope. Implemented from Mahon et al. (2018).
#' @param L The mean length of migrating bedforms.
#' @param S The reach-average slope of the river.
#' @param uncert The uncertainty window for empirical parameters. The user should put in a confidence interval. A lower confidence interval will have a narrower error envelope, but will require more tenuous conclusions.
#' @return Returns a list with elements:
#' \code{V}, the dune migration rate in the river (m/s).
#' \code{qs}, the bedload flux in the river (m^3/s).
#' If uncertainty is also requested, the output will include upper and lower bounds for confidence, \code{V_env} and \code{qs_env}
#' @export

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
