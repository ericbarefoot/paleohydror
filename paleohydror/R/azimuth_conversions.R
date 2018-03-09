#	A set of functions to convert 'quadrant' based directional measurements into more easily manipulated 'azimuth' based measurements. Beware, some are more easy to use than others. I recommend 'azStringConvert' since it is the only one with documentation.
#	Eric Barefoot
#	Nov 2017

azconvert = function(datalist) {
	# converts all measurements to hand-corrected azimuths with dip
	# type 1 = strike to the NW and dip to the right
	# type 2 = strike to the NW and dip to the left
	# type 3 = strike to the NE and dip to the right
	# type 4 = strike to the NE and dip to the left
	# datalist is of the form [[1:4]]	
	strike = datalist[[1]]
	dip = datalist[[2]]
	type = datalist[[3]]
	for (i in 1:length(strike)) {
		if (type[i] == 1) {
			strike[i] = 360 - strike[i]
		} else if(type[i] == 2) {
			strike[i] = 180 - strike[i]
		} else if (type[i] == 3) {
			strike[i] = strike[i]
		} else if (type[i] == 4) {
			strike[i] = 180 + strike[i]
		} else {
			stop(paste('You forgot to denote type for datapoint', i, '.', sep = ' '))
		}
	}
	out = cbind(strike, dip)
	return(out)
}
    
lineconvert = function(data) {
    trend = data[[1]]
    plunge = data[[2]]
    type = data[[3]]
    for (i in 1:length(trend)) {
        if(type[i] == 1) {
            trend[i] = 360 - trend[i]
        } else if (type[i] == 2) {
            trend[i] = trend[i]
        } else if (type[i] == 3) {
            trend[i] = 180 + trend[i]
        } else if (type[i] == 4) {
            trend[i] = 180 - trend[i]
        } else {
            stop('you fukked up')
        }
    }
	out = cbind(trend, plunge)
	return(out)
}

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
