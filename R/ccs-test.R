#' @description bla!
#' 
#' @details bla, bla?
#'
#' @import  Matrix SparseM geigen reshape2 mvtnorm, parallel
#' @name ccs.test
#' @aliases ccs.test
#' @title Estimates a common cointegration space as in Callot (2014).  
#' @author Laurent Callot \email{l.callot@@vu.nl}
#' 
#' @param Y A molten data set containing the 5 columns. *indiv*: the individual (country, sector,...) name as factor levels. *Year*: the year of observation. *Period*: the period (day, week, month, quarter) index as consecutive integers. *Variable*: the name of the variable. *value*: the value of the variable at the given time for a given individual. 
#' @param W either NULL, A matrix, 3 dimensional array with first and second dimension equal to the number of individuals. Columns must sum to one, and values on the diagonal must be zero. In case the input is an array, its third dimension should be equal to the number of observations. In case of NULL, equal weights are used.   
#' @param rank An integer greater or equal to 1 indicating the co-integration rank of the individual systems. 
#' @param lags An integer indicating the number of lags of the VAR in level.
#' @param BS The number of bootstrap iteration. Default: 99.
#' @param ddet An integer indicating the type of deterministics to use, following the typology by Johansen 1988.
#' @param Y0 A molten set of global exogenous data to use. Its columns must be identical to those of Y, and have the same number of observations. 
#' @param seasonal A boolean indicating whether seasonal dummies should be used. Default FALSE
#' @param ccsinit Optional, a matrix with *rank* columns containing the inital value. If NULL, default to the average of the unconstrained estimator. 
#' @param cores An integer greater or equal to 1. if it is greater than 1 the bootstrap is parallelized with mclapply. default to 1. 
#' @param bs.method The type of bootstrap to use: _resample_ or _wild_.
#'
#' @return A list with *stat* the panel trace test statistic for each hypothesis. *pval* the corresponding p-value. A number or other useful statistics are returned. 
#'
#'
#'
#'
#'
#'
#' @export
ccs.test <- function(Y,W=NULL,rank=1,lags=2,BS=99,ddet=0,Y0=NULL,seasonal=FALSE,ccsinit=NULL,cores=1,bs.method='resample'){

  # Initialization
	pcvar.est <- list()
	pcvar.trace <- NULL
	ev <- NULL
	evc <- NULL
	season.freq <-ifelse(seasonal, length(unique(Y$Period)) , 1)
	N <- length(unique(Y$indiv))
  
	# Build the W matrices.
	bigW <- .mkW(W,N=N,nvar=length(unique(Y$Variable)))
	
  # Unrestricted model
  # Estimating the CVAR and computing the rank test for each individuals.
	model <- trace.test(Y,W,lags,ddet,Y0,seasonal)
	# Build the PCVAR for each rk assumption.
	pcvar.coef <- .mkpcvar(model$est,ddet,lags)
	# Get the unrestricted residuals.
	resu <- .mkresid(model$est,rank)
  
  # initial value if not provided. Average of the normalized unrestricted betas.
  if(W=='zero') exo <- FALSE
  else exo <- TRUE
  if(is.null(ccsinit))ccsinit<-.ccsinit(Y,pcvar.coef,rank,exo)
  
  # Calling the core function that calls the optmizer
  ccscore <- .ccs.core(Y,W,rank,lags,ddet,Y0,seasonal,season.freq,ccsinit,exo)
  ccsopt  <- ccscore$ccsopt
  cvarlst <- ccscore$cvarlst
  rm(ccscore)
  gc()
  
  # Constructing and checking the CCS PCVAR
	# Computing the indiviual SR params 
  cvarccs <- list()
	for(i in 1:length(cvarlst))
	{
	  CVAR    <- cvarlst[[i]]
	  cvsr    <- .ccsSR(CVAR,rbind(ccsopt$ccs,0*ccsopt$ccs),ddet,lags,seasonal,season.freq,exo)
    cvarccs[[i]] <- list()
    cvarccs[[i]][[rank+1]] <- cvsr 
	}  
  # Make the PCVAR under CCS
	pcvar.ccs <- .mkpcvar(cvarccs,ddet,lags,rank)
  
	# Get the CCS residuals.
	resccs <- .mkresid(cvarccs,rank)
	
  # Transform the pcvar to levels, and transform the residuals as well.
	pcvar.ccs.lev <- .mk.pcvar.level(pcvar.ccs,resccs,lags,bigW,rank,exo)
  
  # Check the roots of the pcvar.
	rts <- .roots(pcvar.ccs.lev,lags,rank)
	boum <- 0
	for(ro in rts[[rank+1]])boum <- boum+sum(abs(ro)>1.1)  
  
  # storage
  ccs <- list('coefficients'=ccsopt$ccs,'residuals'=resccs,'init'=ccsinit,'rank'=rank,'roots'=rts[[rank+1]],'xplo'=(boum>0))
    
  # Some convergence info
  ccs$iter <- ccsopt$count[1] 
  ccs$cv <- ccsopt$convergence  
  
  # Likelihoods:
  ccs$llc <- .mkll(resccs,N,rank)
  ccs$llu <- .mkll(resu,N,rank)
  # LR test:
  ccs$lr  <- -2*(ccs$llc-ccs$llu)
  ccs$df  <- ncol(ccs$init)*(nrow(ccs$init)-ncol(ccs$init))
  ccs$pv  <- 1-pchisq(ccs$lr,ccs$df)
  
  # Panel stats (std normal)
  ccs$ppv <- sum(-2*log(ccs$pv)-2)/sqrt(4*N) # pooled p-values, reject for large 
  ccs$plr <- sum(ccs$lr-ccs$df)/sqrt(2*ccs$df*N)# pooled LR, reject for large
  
  # Bootstrap
  bsccs <- .bs.ccs(pcvar.ccs.lev,ccs,BS,bs.method='resample',rank,lags,Y,W,ddet,Y0,rts,cores,seasonal,season.freq,exo)  
  
  # Bootstrap statistics
  ccs$bsipv  <- rowMeans(1*(bsccs$lr>=ccs$lr))
  ccs$bspv   <- sum(-2*log(ccs$bsipv)-2)/sqrt(4*N)
  # bs pooled LR  
  ccs$bslr   <- colSums(bsccs$lr-ccs$df)/sqrt(2*N*ccs$df)
  ccs$bslrpv <- mean(1*(ccs$bslr>= ccs$plr))

	return(ccs)
}

