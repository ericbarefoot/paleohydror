#' paleohydror: A package for implementing published routines in paleohydrologic reconstruction.
#'
#' This package offers several functions and routines for deriving estimates of paleohydrology from field observations of rock outcrops. Data such as grain size are important too.
#'
#' @section paleohydror functions:
#'
#' The functions here make up a mix of tools to calculate and estimate modern suficial processes and back-calculate those processes in ancient rivers.
#'
#' \code{settle} calculates the settling velocity of particles following the formulation by Dietrich (1982).
#'
#' \code{lynds_one} estimates ancient river slope using the first method described by Lynds et al. (2014)
#'
#' \code{lynds_two} estimates ancient river slope using the second method described by Lynds et al. (2014)
#'
#' \code{lynds_three} estimates ancient river slope using the third method described by Lynds et al. (2014)
#'
#' \code{xset2H} estimates ancient river dune heights based on cross-set thicknesses. Methods from Leclair and Paola.
#'
#' \code{H2h} estimates ancient river depths based on dune heights. Methods from Bradley and Venditti.
#'
#' \code{xset2depthSimple} estimates ancient river depths directly from cross-set thicknesses using the simplest assumptions. Methods from Leclair, Paola, and Bradley and Venditti.
#'
#' \code{azStringConvert} Azimuth conversion for 'quadrant' based strike/trend measurements.
#'
#' \code{ancient_1} Estimating sediment flux in ancient rivers following Mahon + McElroy (2018).
#'
#' \code{modern_1} Estimating sediment flux in modern rivers if bedform height is known following Mahon + McElroy (2018).
#'
#' \code{modern_2} Estimating sediment flux in modern rivers if bedform length is known following Mahon + McElroy (2018).
#'
#' \code{bedform_velocity} Estimating bedform migration rate as a function of reach-averaged slope.
#'
#' \code{bedform_bedload} Estimating bedload as a function of bedform height and migration rate. After Simons et al. (1965)
#'
#' \code{height_length} Estimating bedform height if bedform length is known.
#'
#' \code{trampush_slp} estimates paleoslope using an empirical relationship with median bedload grain size and bankfull flow depth from Trampush (2014).
#'
#'
#'
#'
#'
#' @docType package
#' @name paleohydror
NULL
