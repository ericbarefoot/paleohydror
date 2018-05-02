#	A function to convert 'quadrant' based directional measurements into more easily manipulated 'azimuth' based measurements.
#	Eric Barefoot
#	Nov 2017

#' Azimuth conversion for 'quadrant' based strike/trend measurements
#'
#' Takes measurements of strike or trend like N45E and converts it to an azimuth measurement with N = 0 and E = 90.
#'
#' @param string	A string measurement with format [N/S][##][E/W]. Can be either X00X or X0X. Make sure you're using the right hand rule to determine how quadrant measurements are made, otherwise you're screwed.
#' @return A numeric azimuth trend.
#' @export

azStringConvert = function(string) {
#	must be a string with [N/S][##][E/W]
	if (nchar(string) == 4) {
		ns = substr(string, 1, 1)
		ew = substr(string, 4, 4)
		rotate = as.numeric(substr(string, 2, 3))
	} else if (nchar(string) == 3) {
		ns = substr(string, 1, 1)
		ew = substr(string, 3,3)
		rotate = as.numeric(substr(string, 2, 2))
	}
	if (ns == 'N' & ew == 'E') {
		az = rotate
	} else if (ns == 'N' & ew == 'W') {
		az = 360 - rotate
	} else if (ns == 'S' & ew == 'E') {
		az = 180 - rotate
	} else if (ns == 'S' & ew == 'W') {
		az = 180 + rotate
	} else {
		stop('Measurement format is incorrect.')
	}
	return(az)
}
