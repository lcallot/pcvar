#'
#' @name trace.test
#' @aliases trace.test
#' @title Computes the trace test for a pcvar model.  
#' @author Laurent Callot \email{l.callot@@vu.nl}
#' 
#'
#'
#' @param Y A molten data set containing the 4 columns. indiv: the individual (country, sector,...) name as factor levels. Year: the year of observation. Period: the period (day, week, month, quarter) index as consecutive integers. Variable: the name of the variable. value: the value of the variable at the given time for a given individual. 
#' @param W NULL, A matrix, 3 dimensional array with first and second dimension equal to the number of individuals. Columns must sum to one, and values on the diagonal must be zero. In case the input is an array, its third dimension should be equal to the number of observations. In case of NULL, equal weights are used.   
#' @param lags An integer indicating the number of lags of the VAR in level.
#' @param dett An integer indicating the type of deterministics to use, following the typology by Johansen 1988.
#' @param Y0 Optional: a molten set of global exogenous data to use. Its columns must be identical to those of Y, and have the same number of observations. 
#'
#'
#' @return A list with 'est' the individual CVAR and 'trace' the trace test statistics for each model and each hypothesis.  
#'
#'
#'
#'
#'
#'
#' @export
trace.test <- function(Y,W=NULL,lags=2,ddet=1,Y0=NULL,seasonal=FALSE){
	pcvar.est <- list()
	pcvar.trace <- NULL
	ev <- NULL
	evc <- NULL
	season.freq <-ifelse(seasonal, length(unique(Y$Period)) , 1)

	# Construct the weighted averages.
	Ystar <- .mkYstar(Y,W)

	#Push in a // function
	for(i in unique(Y$indiv)){
		# Get data for indiv i
		Yi <- subset(Y,indiv==i)
		YSi <- subset(Ystar,indiv==i)

		# Cast the individual data to matrices:
		aYi <- acast(Yi,  Year + Period ~ Variable)
		aXi <- acast(YSi,  Year + Period ~ Variable)
		if(!is.null(Y0))aXi <- cbind(aXi,acast(Y0, Year + Period ~ Variable))

		# Get the data-elements of the CVAR:
		cvar.dat <- .mkCVAR(aYi,aXi,lags)

		# RRR (eigenvecs in columns)
		eig <- .RRR(cvar.dat,ddet,seasonal,season.freq)
		# SR for each possible cointegration rank 
		cvar.sr <- .cpSR(cvar.dat,eig,ddet,lags,seasonal,season.freq)
		pcvar.est[[i]] <- cvar.sr

		# Rk test stats
		evr <- sort(sort(abs(eig$values),decreasing=TRUE)[1:ncol(aYi)])
		trt <- sort(cumsum(log(1-evr)),decreasing=TRUE)
		trace.stat <- -nrow(aYi) * trt
		pcvar.trace <- rbind(pcvar.trace,c(0,trace.stat))

		# Storing the eigenvalues
		ev <- rbind(ev,sort(evr,decreasing=TRUE))
		# Storing the eigenvectors
		evc <- cbind(evc,eig$vectors[,order(abs(eig$values),decreasing=TRUE)[1:ncol(aYi)]])

	}
	return(list('est'=pcvar.est,'trace'=pcvar.trace,'eigenvalues'=ev,'eigenvectors'=evc))
}
