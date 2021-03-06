% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mahon_bedform_relations.r
\name{modern_2}
\alias{modern_2}
\title{Estimating sediment flux in modern rivers if bedform height is unknown, but bedform length is known.}
\usage{
modern_2(L, S, uncert = NULL)
}
\arguments{
\item{L}{The mean length of migrating bedforms.}

\item{S}{The reach-average slope of the river.}

\item{uncert}{The uncertainty window for empirical parameters. The user should put in a confidence interval. A lower confidence interval will have a narrower error envelope, but will require more tenuous conclusions.}
}
\value{
Returns a list with elements:
\code{V}, the dune migration rate in the river (m/s).
\code{qs}, the bedload flux in the river (m^3/s).
If uncertainty is also requested, the output will include upper and lower bounds for confidence, \code{V_env} and \code{qs_env}
}
\description{
\code{modern_2} estimates the sediment flux (\code{qs}) of modern rivers based on bedform length and river slope. Implemented from Mahon et al. (2018).
}
