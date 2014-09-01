#' @description bla!
#' 
#' @details bla, bla?
#'
#'
#' @name mc.ccs.test
#' @aliases mc.ccs.test
#' @title Performs a Monte Carlo experiment on the performance of the panel Common Co-integration space estimator and test of Callot (2014).  
#' @author Laurent Callot \email{l.callot@@vu.nl}
#' 
#'
#'
#' @param obs An integer, the number of observations to generate.
#' @param N an integer, the number of cross section units.
#' @param nvar The number of variables for each cross section unit.
#' @param rank The cointegration rank of the system. Default 1.
#' @param BS The number of bootstrap iterations for the test. Default: 99.
#' @param MC The number of monte carlo iterations. Default: 1.
#' @param Alpha Adjustment matrix (dimension nvar * rank) or list (length N) of parameter matrices.
#' @param Beta Co-integration matrix (dimension 2nvar * rank) or list (length N) of parameter matrices.
#' @param Lambda0 Contemporaneous dependency matrix (dimension nvar  *  nvar) or list (length N) of parameter matrices.
#' @param Gammal List (length is number of lagged first differences) of matrices (N  *  N) or lists (length N) of matrices.
#' @param Omega The covariance matrix.
#' @param err.dist The distribution of the innovations, _gaussian_ or _t_  
#' @param t.df If the innovations are _t_ distributed, the number of degrees of freedom. default 3.
#' @param det.type An integer indicating the type of deterministics to use, following the typology by Johansen 1988.
#' @param burn.smpl The Number of burned observations used to generate the data. 
#' @param res.dep Dependency of the residuals, _iid_ (default) or _garch_.
#' @param garchspec See fGarch package. 
#' @param cdet.load The loadings on the deterministics.
#' @param W The weighting scheme. Default: equal.
#' @param bs.method 'resample' for iid resampling, 'wild' for gaussian wild bootstrap.
#' @param ncores The number of cores, default 1.
#' @param llest Should the Larsson and Lyhagen 2007 JBES CCS estimator be computed? Default FALSE.
#'
#'
#' @return A list. 
#'
#'
#'
#'
#'
#'
#' @export
mc.ccs.test <- function(obs,N,nvar,rank=1,BS=99,MC=1,W='equal',Alpha=NULL,Beta=NULL,Lambda0=NULL,Gammal=NULL,Omega,err.dist='gaussian',det.type=1,cdet.load=c(1,1),bs.method='resample',t.df=3,burn.smpl=10,res.dep='iid',garchspec=NULL,ncores=1,llest=TRUE){

	# 0/ Chk the arguments

	# 1/ Initialisation
	lags <- 1+length(Gammal)
	bsmc <- list()
	if(is.null(W))W <- 'equal'

	# 2/ Chk that the DGP is not explosive and the inputs valid before starting MC.

	# Generate a set of data and hope for no crash.
	Ysim <- gen.pcvar(obs=obs,N=N,nvar=nvar,W=W,Alpha=Alpha,Beta=Beta,Lambda0=Lambda0,Gammal=Gammal,Omega=Omega,err.dist=err.dist,t.df=t.df,burn.smpl=burn,res.dep=res.dep,garchspec=garchspec)

	# check the roots
	if(sum(abs(Ysim$roots[[1]])>1.01)>0)stop(paste('The input data generating process is explosive. Max root: ',max(abs(Ysim$roots[[1]])),'.',sep=''))
	bsmc$dgp.roots <- Ysim$roots

	# 3/ MC loop
	# Parallel 
	mcccs <- mclapply(1:MC,.mcccs.i,MC,obs,N,nvar,BS,rank,
		    W,Alpha,Beta,Lambda0,Gammal,Omega,err.dist,
		    det.type,cdet.load,bs.method,t.df,res.dep,garchspec,burn.smpl,llest,mc.cores=ncores,mc.preschedule=FALSE)

	#mcccs <- lapply(1:MC,.mcccs.i,MC,obs,N,nvar,BS,rank,
	#	    W,Alpha,Beta,Lambda0,Gammal,Omega,err.dist,
	#	    det.type,cdet.load,bs.method,t.df,res.dep,garchspec,burn.smpl,llest)
	#cat('\n')

  mcagg <- .agg.mc.ccs(mcccs,llest)
  
	return(mcagg)
}



