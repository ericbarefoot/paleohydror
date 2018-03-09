#' Calculate settling velocity for a particle.
#'
#' \code{settle} calculates the settling velocity for a particle, given some parameters.
#' @param ws_mode Denotes whether to use a 'simple' equation, neglecting the shape and roundness of particles, or to use the more 'complex' formulation
#' @param D The grain size you wish to calculate settling velocity for.
#' @param constants A list of constants for calculation, need in order (g) gravitational constant, (rhos) sediment density, (rhof) fluid density, (nuu) kinematic viscosity, (csf) corey shape factor, and (PI) Powers roundness index. For the simple formulation, neglect the last two. Defaults to constants = list(g = 9.8, rhos = 2650, rhof = 1000, nuu  = 1.004e-6, csf = 0.7, PI = 3.5)
#' @return Returns a list with two values \code{$ws} and \code{$wstar} which are the settling velocity and dimensionless settling velocity, respectively. 
#' @export
#' @examples
#' settle('simple', 120e-6)
#' settle('complex', 120e-6, list(g = 9.8, rhos = 2650, rhof = 1000, nuu  = 1.004e-6, csf = 0.7, PI = 3.5)

#	function to replicate two implementations of Dietrich 1982 settling velocity measurement for a given grain size and density.
#	Eric Barefoot
#	Nov 2017

#	takes constants as a list with names g, rhos, rhof, nuu, csf, PI with an option to leave off csf and PI for a simpler case

settle = function(ws_mode, D, constants) {
	
	if(missing(constants)) {
		g = 9.8			# m/s^2
		rhos = 2650		# assuming ~quartz sand.
		rhof = 1000		# assuming water
		nuu  = 1.004e-6	# kinematic viscocity of water at 20C in Pa s
		csf = 0.7		# given csf for grains 
		PI = 3.5		# average for natural grains. Varies from 0-6.
	}			
	else if (ws_mode == 'simple' & length(constants) != 4) {
		warning('incorrect number of constants; using defaults for the simple method')
		g = 9.8			# *m/s^2
		rhos = 2650		# *assuming ~quartz sand.
		rhof = 1000		# *assuming water
		nuu  = 1.004e-6}# *kinematic viscocity of water at 20C in Pa s}
	else if (ws_mode == 'complex' & length(constants) != 6) {
		warning('incorrect number of constants; using defaults for the complex method')
		g = 9.8			# m/s^2
		rhos = 2650		# assuming ~quartz sand.
		rhof = 1000		# assuming water
		nuu  = 1.004e-6	# kinematic viscocity of water at 20C in Pa s
		csf = 0.7		# given csf for grains 
		PI = 3.5}		# average for natural grains. Varies from 0-6.
	else if (ws_mode == 'simple') {
		g = constants$g			# m/s^2
		rhos = constants$rhos	# assuming ~quartz sand.
		rhof = constants$rhof	# assuming water
		nuu  = constants$nuu	# kinematic viscocity of water at 20C in Pa s
	}	
	else if (ws_mode == 'complex') {
		g = constants$g			# m/s^2
		rhos = constants$rhos	# assuming ~quartz sand.
		rhof = constants$rhof	# assuming water
		nuu  = constants$nuu	# kinematic viscocity of water at 20C in Pa s
		csf = constants$csf		# given csf for grains 
		PI = constants$PI		# average for natural grains. Varies from 0-6.
	}		

	#	calculate dimensionless grain size

	Dstar = (rhos - rhof) * ((g * D^3) / (rhof * nuu^2)) 

	#	simple case
	
	if (ws_mode == 'simple') {

		Dsl = log10(Dstar)
		
		R1 = -3.76715 + 1.92944*Dsl - 0.09815*Dsl^2 - 0.00575*Dsl^3 + 0.00056*Dsl^4
		
		wstar = 10^(R1)
	
	}
	
	if (ws_mode == 'complex') {
	#	some simplifying common terms
		
		Dsl = log10(Dstar)
		tantan = tanh(Dsl - 4.6)
		csf_t = 1 - csf
		
		#	now the other parameters for settling velocity 
		#	NOTE: lynds et al neglect the other terms, so theres a switch, one with the shape factor and powers index, and one without. default is simple

		R1 = -3.76715 + 1.92944*Dsl - 0.09815*Dsl^2 - 0.00575*Dsl^3 + 0.00056*Dsl^4

		R2 = (log10(1 - (csf_t/0.85)) - csf_t^2.3*tantan + 0.3*(0.5 - csf)*csf_t^2*(Dsl - 4.6))

		R3 = (0.65 - ((csf/2.83)*tantan))^(1+(3.5-PI)/2.5)
		
		wstar = R3 * 10^(R1 + R2)
	}

	#	conversion back to dimensional settling velocity

	ws = (wstar * (rhos - rhof) * g * nuu / rhof)^(1/3)

	return(list(ws = ws, wstar = wstar))

}