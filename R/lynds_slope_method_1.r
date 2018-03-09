#' Simple method to calculate paleoslope from preserved sediments
#'
#' \code{lynds_one} calculates the slope of an ancient river, provided you know the flow depth, the grain size of the bed material, and make some assumptions about shear stress in the river. Implemented from Lynds et al. (2014) 
#' @param tau_mode Denotes whether to use a fixed value for the shields number of the flow, or to use a logarithmically-spaced range going from one to ten. 
#' @param D The median grain size of bedload sediments.
#' @param rhos The density of the sediment in kg/m^3
#' @param rhof The density of the fluid in kg/m^3
#' @return Returns a single dimensionless number giving the slope.  
#' @export

#	First method to reconstruct paleoslope from Lynds et al 2014
#	Eric Barefoot
#	Nov 2017

#	define needed variables and data types

lynds_one = function(D, H, tau_mode = 'fixed', rhos = 2650, rhof = 1000) {
	
	#	D is the D50b or median bedload grainsize
	#	H is Hbf or the bankfull flow depth
	
	rhos = rhos
	rhof = rhof
	
	#	key assumption is that shields number is equal to one. But paper varies it by up to ten. Provide here a choice for selecting either the recommended value of 1 or a range of estimates from one to ten spaced evenly in log space. 
	
	tau_bf_50 = switch(tau_mode, fixed = 1, range = round(log10(seq(1.258935,10,length = 10)) * 10, 4))
	
	#	R is assumed in the paper, but we will calculate it here, and provide defaults above.
	
	R = rhos/rhof - 1
	
	#	now we calculate the slope based on equation (3) in Lynds et al 2014
	
	S = (tau_bf_50 * R * D)	/ H
	
	return(S)
	
}