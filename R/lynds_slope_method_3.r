#' Slackwater/suspension method to calculate paleoslope from preserved sediments
#'
#' \code{lynds_three} calculates the slope of an ancient river, provided you know the flow depth, the grain size of the bed material, and the maximum grain size of suspended sediment. Implemented from Lynds et al. (2014) 
#' @param ws_mode Allows you to choose one of two different methods for calculating the settling velocity of the bedload particles. 'simple' ignores the Corey Shape Factor and the Powers Index. See documentation for \code{settle} for more information. For now, this feature is disabled.
#' @param Dbed The median grain size of bedload sediments in meters.
#' @param Dsusp The maximum grain size of suspended sediments in meters.
#' @param H The paleoflow depth in meters.
#' @param rhos The density of the sediment in kg/m^3
#' @param rhof The density of the fluid in kg/m^3
#' @return Returns a list with two elements: \code{S}, the dimensionless slope, and \code{F}, the ratio of total boundary shear stress to skin friction shear stress.  
#' @export

#	Third method to reconstruct paleoslope from Lynds et al 2014
#	Eric Barefoot
#	Nov 2017

#	define needed variables and data types

#source(paste0(fd,'settle.r'))

lynds_three = function(Dbed, Dsusp, H, rhos = 2650, rhof = 1000, ws_mode = 'simple') {
	
	#	D is the D50b or median bedload grainsize
	#	H is Hbf or the bankfull flow depth
	
	#	constants
	#	NOTE: a number of these relationships depend on older dune scaling relationships. may be worthwhile upgrading to something as recommended by venditti 2017
	
	
	g = 9.8
	Cd = 0.21			#	drag coefficient for separated flow
	k = 0.407			#	von Karman constant
	hd = 0.3 * H		#	height of dunes. *** may also be able to calculate a priori from leclair?? ***
	lam = hd / 0.063 	#	for ratio of hd to lambda
	zo = 0.056 * 2 * Dbed	#	boundary layer thickness
	
	wsmax = settle('simple', Dsusp, constants = list(g = g, rhos = rhos, rhof = rhof, nuu = 1.004e-6))$ws
	
	F = 1 + (Cd / (2 * k^2)) * (hd / lam) * (log(hd / zo) - 1)^2
	
	#	final slope equation
	
	S = (F * wsmax^2) / (g * H)
	
	return(list(S = S, F = F))
	
}