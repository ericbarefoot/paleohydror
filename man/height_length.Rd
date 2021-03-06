% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mahon_bedform_relations.r
\name{height_length}
\alias{height_length}
\title{Bedform length to height relation}
\usage{
height_length(delta0, delta1, L)
}
\arguments{
\item{delta0}{Empirical constant. Value and associated uncertainty from Table 1 in Mahon + McElroy 2018.}

\item{delta1}{Empirical constant. Value and associated uncertainty from Table 1 in Mahon + McElroy 2018.}

\item{L}{The mean bedform length.}
}
\value{
Returns a vector \code{hd}, the bedform migration rate in m/s.
}
\description{
\code{height_length} estimates the mean height of bedforms based on the mean bedform length. Implemented from Mahon et al. (2018)
}
