#	functions to calculate paleoflow depth from preserved cross-sets.
#	Eric Barefoot
#	Nov 2017
#
#	largely taken from work by Suzanne Leclair, Chris Paola and methods in a review paper by Bradley and Venditti

#	first step is to relate the cross-set (xset) thickness to dune height (H)

#	require packages

require(stats4)

#' Dune height calculation from cross-set thicknesses
#'
#' \code{xset2H} takes measured cross-set thicknesses (\code{y}) and converts it to dune heights (H).
#'
#' @param y 	A vector of measured cross-set thicknesses in meters. 
#' @param mode	Takes two values: \code{simple} and \code{pdfFit}. \code{simple} uses a simple approximation to determine the beta value in Leclair's equation relating x-sets to dune heights. \code{pdfFit} uses a more complicated procedure, fitting a PDF from Paola  + Borgman to get the beta value. \code{pdfFit} is better suited when there are adequate data to estimate the distribution.
#' @param output 	Takes two values: \code{vanilla} and \code{rich}. \code{vanilla} only returns the estimated value of \code{H}, or dune height. \code{rich} returns \code{H}, as well as \code{beta} and either \code{startVal} which is the inital guess of a for the MLE algorithm in \code{pdfFit}, or a span of \code{H} values giving a reasonable range of one particular assumption. 
#' @return Either a single numeric if \code{output = vanilla}, or a list if \code{output = rich}.
#' @export

xset2H = function(y, mode = 'simple', output = 'vanilla') {
	
	#	function takes two modes: 'simple' and 'pdfFit'
	#	y is a vector of cross-set thicknesses -- MUST BE IN METERS
	
	y = y * 1000 # 	convert to millimeters
	
	#	gives output in meters
	
	n = length(y)
	
	#	first a function to estimate the average dune height via fitting a PDF of x-set thicknesses to find the right parameters
	pdfFit = function(y) {
		
		#	define the probability density function from Paola and Borgman
		pfun = function(a,s) {
			p = a * exp(-a*s) * (exp(-a*s) + a * s - 1) / (1 - exp(-a*s))^2
			return(p)
		}
		
		#	construct a likelihood function which takes data as an external argument and the parameter as the internal argument.
		loglikely = function(data) {
			function(a) {
				R = R = pfun(a, data)
				-sum(log(R))
			}
		}
		
		#	construct likelihood function for these data. 
		LL = loglikely(y)
		
		#	default starting a value for fitting
		#	and error handling sequence for finding an inital guess value for that works and doesnt crash the MLE.
		worked = FALSE
		startVal = 1
		try(stop(''), silent = T)
		j = 1
		while (!worked) {

			try_MLE = try(coef(mle(minuslogl = LL, start = list(a = startVal))), silent = T)

			if (class(try_MLE) == 'try-error' & j < 1000) {

				options(warn = -1)
				msg = as.character(try_MLE)

				if (grepl('initial value', msg)) {
					startVal = startVal / 2
					j = j + 1
#					message(paste(startVal))
				} else {
					stop('Something else happened, its real bad')
				}
				
			} else {
				
				worked = TRUE
			
			}
		}
		
		#	estimate parameter using MLE function
		a = coef(mle(minuslogl = LL, start = list(a = startVal)))
		
		#	get beta a la Leclair et al. 0.9 is for dune height/scour correction
		beta = 0.9/a
		
		#	get dune height from beta from Leclair et al. two ways. 
#		H = 2.22 * beta ^ 1.32
		H = 5.3 * beta + 0.001 * beta ^ 2
		
		list(H = as.numeric(H) / 1000, beta = as.numeric(beta), startVal = startVal)
	}
	
	#	Here, a simpler function using only the empirical estimations made by Leclair et al.
	simple = function(y) {
		sm = mean(y, na.rm = T) #	find average set thickness
		hsd_tssd = seq(from = 0.7, to = 1.1, by = 0.1) #	find ratio of scour deviation to height deviation. supposedly an average of 0.9, but ranges randomly from 0.7 to 1.1 this function spits out a range.
		beta = sm / (1.64493 / hsd_tssd) #	use that to calculate beta
#		beta = sm / 1.8
		
#		H = 2.22 * beta ^ 1.32
		H = 5.3 * beta + 0.001 * beta ^ 2 #	use same equation as in pdfFit to get the average height
		
		list(H = mean(H) / 1000, beta = as.numeric(beta), vecH = H / 1000)
	}
	
	#	offer some recommendations to the user. 
	if(n >= 25 & mode == 'simple') {
		
		message('this dataset may benefit from running a more complicated PDF based function.')
		
		H = simple(y)
		
	} else if (n < 25 & mode == 'simple' ) {
		
		H = simple(y)
		
	} else if (n < 25 & mode == 'pdfFit') {
	
		message('this dataset is small, and may benefit from running a less complicated empirical.')
		
		H = pdfFit(y)
		
	} else if (n >= 25 & mode == 'pdfFit') {
				
		H = pdfFit(y)
		
	}
	
	switch(output , 
		   vanilla = return(H$H),
		   rich = return(H)
		  )
}

#' Paleo-depth calculation from dune height
#'
#' \code{H2h} takes estimated dune heights (H) and calculates a flow depth.
#'
#' @param H 	A vector of measured cross-set thicknesses in meters. 
#' @param CL	A pre-specified confidence level. If you pick a low confidence (50\%), the range in estimated flow depth will be smaller. 
#' @param output 	Takes two values: \code{vanilla} and \code{rich}. \code{vanilla} only returns the estimated value of \code{h}. \code{rich} returns \code{h}, as well as \code{CI}, the confidence interval around the \code{h} estimate. 
#' @return Either a single numeric if \code{output = vanilla}, or a list if \code{output = rich}.
#' @export

#	second step is to relate the dune height to the flow depth, following Bradley and Venditti

H2h = function(H, CL = 50, output = 'vanilla') {
	
	#	H is the height of a dune in meters. 
	#	or if H is a vector, it is a set of dune heights
	#	CL is the confidence level in percent
	
	#	confidence ranges
	
	CR = data.frame(CL = c(50, 80, 90, 95), upper = c(10.1, 14.6, 23.9, 29.5), lower = c(4.4, 3.1, 2.8, 2.7))
	
	Cind = which.min(abs(CR$CL - CL)) # match closest one from the table
	
	h = 6.7 * H	#	using the median parameter from Bradley and Venditti
	
	CI = list(upper = H * CR$upper[Cind], lower = H * CR$lower[Cind]) # and based on Confidence Level, assign bounds of confidence
	
	returnval = switch(output ,
					   vanilla = return(h),
					   rich = return(list(h = h, CI = CI))
					  )	
}

#	and now a function that combines both into one function, so you can go straight from xsets to flow depth with no extra frills, or depth

#' Paleo-depth calculation from cross-set thicknesses
#'
#' \code{xset2depthSimple} uses the simplest implementation of both \code{xset2H} and \code{H2h}. It uses the vanilla options of each one and spits out an estimate of paleo-depth based on the defaults. Look in the code for the other functions to see what those defaults are. 
#' @param y 	A vector of measured cross-set thicknesses in meters. 
#' @return The paleo-depth estimate in meters. 
#' @export

xset2depthSimple = function(y) {
	
	#	y is a vector of thicknesses of xsets in meters
	
	H = xset2H(y)
	
	h = H2h(H)
	
	return(h)
	
}

