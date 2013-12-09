#'
#' @name rank.test
#' @aliases rank.test
#' @title Applies the sequential panel rank test procedure of Callot (2013) to a pcvar model.  
#' @author Laurent Callot \email{l.callot@@vu.nl}
#' 
#'
#'
#' @param Y A molten data set containing the 5 columns. *indiv*: the individual (country, sector,...) name as factor levels. *Year*: the year of observation. *Period*: the period (day, week, month, quarter) index as consecutive integers. *Variable*: the name of the variable. *value*: the value of the variable at the given time for a given individual. 
#' @param W either NULL, A matrix, 3 dimensional array with first and second dimension equal to the number of individuals. Columns must sum to one, and values on the diagonal must be zero. In case the input is an array, its third dimension should be equal to the number of observations. In case of NULL, equal weights are used.   
#' @param lags An integer indicating the number of lags of the VAR in level.
#' @param dett An integer indicating the type of deterministics to use, following the typology by Johansen 1988.
#' @param BS The number of bootstrap iterations for the test. 
#' @param Y0 A molten set of global exogenous data to use. Its columns must be identical to those of Y, and have the same number of observations. 
#' @param bs.method 'resample' for iid resampling, 'wild' for gaussian wild bootstrap.
#' @param seasonal A boolean indicating whether seasonal dummies should be used. Default FALSE
#' @param cores An integer greater or equal to 1. if it is greater than 1 the bootstrap is parallelized with mclapply. default to 1. 
#'
#'
#' @return A list with *stat* the panel trace test statistic for each hypothesis. *pval* the corresponding p-value. A number or other useful statistics are returned. 
#'
#'
#'
#'
#'
#'
#' @export
rank.test <- function(Y,W=NULL,lags=2,ddet=0,BS=99,Y0=NULL,bs.method='resample',seasonal=FALSE,cores=1){

	
	# 0/ Input checks
	if(!(bs.method%in% c('resample','wild')))stop('invalid bootstrap method selected')
	# Check the format of Y


	# Estimating the CVAR and computing the rank test for each individuals.
	model <- trace.test(Y,W,lags,ddet,Y0,seasonal)

	# Build the W matrices.
	bigW <- .mkW(W,N=length(unique(Y$indiv)),nvar=length(unique(Y$Variable)))

	# Build the PCVAR for each rk assumption.
	pcvar.coef <- .mkpcvar(model$est,ddet,lags)

	# Get the original residuals.
	res <- .mkresid(model$est)
	# Transform the pcvar to levels, and transform the residuals as well.
	pcvar.lev <- .mk.pcvar.level(pcvar.coef,res,lags,bigW)
	# Check the roots of the pcvar.
	rts <- .roots(pcvar.lev,lags)
	boum <- 0
	for(ro in rts)boum <- boum+sum(abs(ro)>1.1)


	if(boum==0){
		# Bootstrap the co-integration rank test.
		bs.rk <- .bs.rank(pcvar.lev,BS,model,bs.method,lags,Y,W,ddet,rts,cores)

		# Compute the panel statistics.
		rktmp <- rowMeans(apply(bs.rk$trace,3,function(x,bmk){return(x>=bmk)},bs.rk$trace[,,1]))
		panel.rk <- rev(-2*colSums(log(matrix(rktmp,nrow=length(unique(Y$indiv))))))
		panel.pv <- 1-pchisq(panel.rk,df=2*length(unique(Y$indiv))) 

		# Some extra info on the BS test
		avg.lr <- rev(apply(bs.rk$trace,c(2),mean))
		orig.lr <- rev(colMeans(model$trace))
		avg.ipv <- rev(colMeans(matrix(rktmp,nrow=length(unique(Y$indiv)))))

		# Store and return. 
		bsl <- list('stat'=panel.rk,'pval'=panel.pv,'avg.lr'=avg.lr,'orig.lr'=orig.lr,'avg.ipv'=avg.ipv,'bs.rk'=bs.rk,'level'=pcvar.lev,'residuals'=res,'coefficients'=pcvar.coef,'roots'=rts,xplo=FALSE)
		}
	else{warning('The estimated PCVAR is explosive, try again.');bsl <- list('roots'=rts,xplo=TRUE)}


	return(bsl)
}
