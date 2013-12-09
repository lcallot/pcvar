#'
#' @name gen.pcvar 
#' @aliases gen.pcvar
#' @title Generates data from a PCVAR model.  
#' @author Laurent Callot \email{l.callot@@vu.nl}
#' 
#' 
#'
#' @param Y A molten data set containing the 4 columns. indiv: the individual (country, sector,...) name as factor levels. Year: the year of observation. Period: the period (day, week, month, quarter) index as consecutive integers. Variable: the name of the variable. value: the value of the variable at the given time for a given individual. 
#' @param bs.method 'resample' for iid resampling, 'wild' for gaussian wild bootstrap.
#'
#' @return A list. 
#'
#'
#'
#'
#'
#'
#' @export
gen.pcvar <- function(obs,N,nvar,W=NULL,Alpha=NULL,Beta=NULL,Lambda0,Gammal=NULL,Omega,err.dist='gaussian',burn.smpl=10,t.df=1,det.type=1,cdet.load=c(1,1),res.dep,garchspec){

	# Construc the long run parameter matrices.
	if(!is.null(Alpha)){
		if(!is.list(Alpha))ALPHA <- bdiag(rep(list(Alpha),N))
		if(is.list(Alpha))ALPHA <- bdiag(Alpha)
		}
	if(!is.null(Beta)){
		if(!is.list(Beta))BETA <- bdiag(rep(list(Beta),N))
		if(is.list(Beta))BETA <- bdiag(Beta)
		}	

	# Contemporanous dependence matrix
	if(!is.null(Lambda0)) Lambda0 <- matrix(0,nvar,nvar)
	if(!is.list(Lambda0))LAMBDA0 <- bdiag(rep(list(Lambda0),N))
	if(is.list(Lambda0))LAMBDA0 <- bdiag(Lambda0)

	# Storing
	pcvar.coef <- list(list('ALPHA'=ALPHA,'BETA'=BETA,'LAMBDA0'=LAMBDA0))

	# Constructing the lagged difference matrices if needed
	lags <- 1
	if(!is.null(Gammal)){
		lags <- length(Gammal)+1
		for(l in 1:length(Gammal)){
			if(!is.list(Gammal[[l]]))assign(paste('LAG',l,sep='_'),bdiag(rep(list(Gammal[[l]]),N)))
			if(is.list(Gammal[[l]]))assign(paste('LAG',l,sep='_'), bdiag(Gammal[[l]]))
			pcvar.coef[[1]][[paste('LAG',l,sep='_')]] <- get(paste('LAG',l,sep='_'))
			}
		}

	# Construct the covariance matrix	
	if(is.null(Omega)) { Omega <-diag(nvar)}

	if(is.list(Omega)) { OMEGA <- bdiag(Omega) }
	else { if(is.matrix(Omega)){	
			if(nrow(Omega)==nvar)OMEGA <- bdiag(rep(list(Omega),N))
			if(nrow(Omega)==nvar*N)OMEGA <- Omega
			}
		}
	
	
	# iid residuals
	if(res.dep=='iid'){
		# Construct the residual matrix
		if(err.dist=='gaussian')pcvar.resid <- rmvnorm(obs+burn.smpl,mean=rep(0,nrow(OMEGA)),sigma=as.matrix(OMEGA)) 
		if(err.dist=='t')pcvar.resid <- rmvt(obs+burn.smpl,sigma=as.matrix(OMEGA),df=t.df) 
	}


	# Garch residuals, garchspec must be a valid garchSpec object from fGarch. 
	if(res.dep=='garch'){
		pcvar.resid <- do.call(cbind,lapply(1:(N*nvar),garchSim,spec=garchspec,n=obs+burn.smpl,extended=FALSE))		
		pcvar.resid <- data.frame(pcvar.resid)
	}


	# Compute the big weight matrices
	bigW <- .mkW(W=W,N,nvar)

	# Makes the pcvar in level
	pcvar.lev <- .mk.pcvar.level(pcvar.coef,list(pcvar.resid),lags,bigW)
	
	# Check the roots of the pcvar.
	roots <- .roots(pcvar.lev,lags)


	# 3/ generate pseudo data
	Ygen <- matrix(0,ncol=obs+burn+lags,nrow=N*nvar)
	for(o in (lags+1):(obs+lags+burn)){
		# Summing over the lags:
		for(l in 1:lags)Ygen[,o] <- Ygen[,o]+t(as.matrix(pcvar.lev$coefficients[[1]][[paste('L',l,sep='')]]))%*%Ygen[,o-l] 
		# Adding noise:
		Ygen[,o] <-Ygen[,o]+ t(pcvar.resid)[,o-lags]
		# Add deterministics
		}

	return(list('roots'=roots,'Y'=t(Ygen[,-c(1:(lags+burn))])))
}
