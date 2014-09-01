#' @description bla!
#' 
#' @details bla, bla?
#'
#'
#' @name meltdat
#' @aliases meltdat
#' @title Transforms a matrix of data to a molten data frame for use in rank.test..  
#' @author Laurent Callot \email{l.callot@@vu.nl}
#' 
#' 
#' @param Ymat A matrix of data with variables in columns. The variables are assumed grouped by countries.
#' @param indiv.names The names of the individuals
#' @param var.names The names of the variables
#' @param time.index A time index, default to 1:T.
#' @param freq The frequency of the seasons, default 1 if not seasonal. 
#'
#' @return A molten dataframe ready for use in rank test. 
#'
#'
#'
#'
#'
#'
#' @export
meltdat <- function(Ymat,indiv.names,var.names,time.index=1:nrow(Ymat),freq=1){

	N <- length(indiv.names)
	p <- length(var.names)
	obs <- nrow(Ymat)

	mY <- data.frame('value'=matrix(Ymat))
	mY$Year <- rep(rep(time.index,each=freq),N*p)
	mY$Period <- rep(1:freq,N*p*obs)
	mY$Variable <- rep(rep(var.names,each=obs),N)
	mY$indiv <- rep(indiv.names,each=obs*p)

	return(mY)
}
