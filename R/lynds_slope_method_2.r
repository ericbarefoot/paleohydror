#' Variable Shields number method to calculate paleoslope from preserved sediments
#'
#' \code{lynds_two} calculates the slope of an ancient river, provided you know the flow depth, the grain size of the bed material, and make some assumptions about shear stress in the river. The important part here is that the function implements a different shields number depending on whether the bed material is fine or coarse. Implemented from Lynds et al. (2014) 
#' @param ws_mode Allows you to choose one of two different methods for calculating the settling velocity of the bedload particles. 'simple' ignores the Corey Shape Factor and the Powers Index. See documentation for \code{settle} for more information. For now, this feature is disabled.
#' @param D The median grain size of bedload sediments in meters.
#' @param H The paleoflow depth in meters.
#' @param rhos The density of the sediment in kg/m^3
#' @param rhof The density of the fluid in kg/m^3
#' @return Returns a single dimensionless number giving the slope.  
#' @export

#	Second method to reconstruct paleoslope from Lynds et al. (2014)
#	Eric Barefoot
#	Nov 2017

#	function takes two main arguments, the particle size, and the depth, both in meters. Then it uses them to calculate a slope, using a raft of assumptions.
#	to use a more complicated way to calculate settling velocity, use ws_mode = 'complex'

#	requires the 'settle' function, which lives in the same package

#source(file.path(fd,'settle.r'))

lynds_two = function(D, H, ws_mode = 'simple', rhos = 2650, rhof = 1000) {
	
	#	NOTES:
	#	in the future, need to allow options for csf and PI, but implement more smoothly. 
	
	#	D is the D50b or median bedload grainsize IN METERS
	#	H is Hbf or the bankfull flow depth IN METERS
	
	#	constants
	
	g = 9.8				#	gravity m/s^2
	nuu  = 1.004e-6		#	kinematic viscocity in Pa s at 20C
	rhos = rhos			# 	in kg/m^3
	rhof = rhof			#	in kg/m^3
	R = rhos/rhof - 1	#	should be close to 1650 for qtz sand
	csf = 0.7			# 	some csf for grains 
	PI = 3.5			# 	average for natural grains. Varies from 0-6.

	#	key is to find the settling velocity of particles. there exist a wide variety of schemas for this. lynds uses dietrich 1982, realized in function 'settle'.
	
	#	function giving the particle reynolds number for a size class. 
	
	Rep = function(D, R, g, nuu) {
		Rep = (D * sqrt(D * R * g)) / nuu
		return(Rep)
	}
	
	#	now determine Wstar
	
	Wstar = settle(ws_mode = 'simple', D, constants = list(g = g, nuu = nuu, rhos = rhos, rhof = rhof))$wstar
	
	#	now determine Rep
	
	Rp = Rep(D, R, g, nuu)
	
	#	determine the size class that we are in. divide for fine sand or medium sand
	
	if(D >= 0.25e-3 & D < 0.5e-3){szcls = 'medium'}
	else if (D < 0.25e-3 & D >= 0.0625e-3){szcls = 'fine'}
	else if (D >= 0.5e-3) {szcls = 'medium'; warning('coarse sand, choosing medium value')}
	else if (D < 0.0625e-3) {szcls = 'fine'; warning('very fine silt, choosing fine value')}
	else {stop('something is wrong with the ustar/ws values')}
	
	#	switch defines value for ustar/ws ratio
	
	u_w = switch(szcls,
				 medium = 1.6,
				 fine = 3.1
				)
	
	#	now solve equation 9 in Lynds 2014 for slope
	
	S = ((R * D) / H) * (u_w * (Wstar / Rp)^(1/3))^2
	
	return(S)
	
}