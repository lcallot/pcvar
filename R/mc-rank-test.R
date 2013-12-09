#'
#' @name mc.rank.test
#' @aliases mc.rank.test
#' @title Performs a Monte Carlo experiment on the performance of the panel rank test procedure of Callot (2013) to a pcvar model.  
#' @author Laurent Callot \email{l.callot@@vu.nl}
#' 
#'
#'
#' @param Y A molten data set containing the 4 columns. indiv: the individual (country, sector,...) name as factor levels. Year: the year of observation. Period: the period (day, week, month, quarter) index as consecutive integers. Variable: the name of the variable. value: the value of the variable at the given time for a given individual. 
#' @param W NULL, A matrix, 3 dimensional array with first and second dimension equal to the number of individuals. Columns must sum to one, and values on the diagonal must be zero. In case the input is an array, its third dimension should be equal to the number of observations. In case of NULL, equal weights are used.   
#' @param lags An integer indicating the number of lags of the VAR in level.
#' @param dett An integer indicating the type of deterministics to use, following the typology by Johansen 1988.
#' @param BS The number of bootstrap iterations for the test. 
#' @param Y0 A molten set of global exogenous data to use. Its columns must be identical to those of Y, and have the same number of observations. 
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
mc.rank.test <- function(obs,N,nvar,BS=99,MC=1,W='equal',Alpha=NULL,Beta=NULL,Lambda0=NULL,Gammal=NULL,Omega,err.dist='gaussian',det.type=1,cdet.load=c(1,1),bs.method='resample',t.df=1,burn.smpl=10,res.dep='iid',garchspec=NULL,ncores=1){

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

	# Status indicator (thanks stack overflow!)
	# fifo method for counter.
	# txtProgress bar for output.
	prog.indic <- local({ #evaluates in local environment only
		# Using fifo connection
		f <- fifo(tempfile(), open="w+b", blocking=T) # open fifo connection
		assign(x='f',value=f,envir=.GlobalEnv)
		pb <- msgProgressBar(min=0, max=MC,style=3,file=stderr())

		if (inherits(fork(), "masterProcess")) { #progress tracker
			# Child
			progress <- 0.0
			while (progress < MC && !isIncomplete(f)) {
				msg <- readBin(f, "double")
			    	progress <- progress + as.numeric(msg)

				# Updating the progress bar.
				setMsgProgressBar(pb,progress)
				} 
			exit()
	    		}


		# Parallel 
		mcrk <- mclapply(1:MC,.mcrk.i,MC,obs,N,nvar,BS,
			    W,Alpha,Beta,Lambda0,Gammal,Omega,err.dist,
			    det.type,cdet.load,bs.method,t.df,res.dep,garchspec,burn.smpl,mc.cores=ncores)

		#mcrk <- lapply(1:MC,.mcrk.i,MC,obs,N,nvar,BS,
		#	    W,Alpha,Beta,Lambda0,Gammal,Omega,err.dist,
		#	    det.type,cdet.load,bs.method,t.df,res.dep,garchspec,burn.smpl)
		cat('\n')
		assign(x='mcrk',value=mcrk,envir=.GlobalEnv)
		close.msgProgressBar(pb)
		close(f)
		})


	# Not parallel
	#mcrk <- lapply(1:MC,.mcrk.i,obs,N,nvar,BS,
	#	    W,Alpha,Beta,Lambda0,Gammal,Omega,err.dist,
	#	    det.type,cdet.load,bs.method,t.df,burn.smpl)

	# 4/ Aggregate rank test outcome.
	aggmc <- .agg.mc(mcrk)
	aggmc$dgproots <- Ysim$roots[[1]]

	return(aggmc)
}


# Aggregate the MC results
.agg.mc <- function(mcrk){

	MC <- length(mcrk)

	# The storage list
	aggmc <- list('stat'=NULL,'pval'=NULL,'avg.lr'=NULL,'avg.ipv'=NULL,'orig.lr'=NULL,
		      'trace'=array(NA,dim=c(dim(mcrk[[1]]$trace),MC)),
		      'roots'=array(NA,dim=c(dim(mcrk[[1]]$roots),MC))
		      ,'bad'=NULL,'residuals'=list(),'coefficients'=list())

	for(mc in 1:MC){
		# getting iteration mc
		m <- mcrk[[mc]]
		
		# aggregating the aggregates
		for(s in c('stat','pval','avg.lr','avg.ipv','orig.lr'))aggmc[[s]] <- rbind(aggmc[[s]],m[[s]]) 
		
		# agg the bs trace stats and bs dgp roots
		aggmc$trace[,,,mc] <- m$trace
		aggmc$roots[,,mc] <- m$roots

		# number of bad tries
		aggmc$bad <- c(aggmc$bad,m$bad)

		# List of coefficients and residuals
		aggmc$coefficients[[mc]] <- m$coefficients
		aggmc$residuals[[mc]] <- m$residuals
		}

	return(aggmc)
}



# Monte Carlo internals:
.mcrk.i <- function(mc,MC,obs,N,nvar,BS,
		    W,Alpha=NULL,Beta=NULL,Lambda0=NULL,Gammal=NULL,
		    Omega,err.dist='gaussian',
		    det.type=1,cdet.load=c(1,1)
		    ,bs.method='resample',t.df=1,res.dep='iid',garchspec=NULL,burn.smpl=10){

	# Initialisation
	set.seed(mc*123)
	good.rk <- FALSE
	bad.rk <- 0
	lags <- length(Gammal)+1

	# Working until we get a non-explosive bootstrap DGP
	while(!good.rk){
		# a/ Generate data
		Ysim <- gen.pcvar(obs=obs,N=N,nvar=nvar,W=W,Alpha=Alpha,Beta=Beta,Lambda0=Lambda0,Gammal=Gammal,Omega=Omega,err.dist=err.dist,t.df=t.df,burn.smpl=burn,res.dep=res.dep,garchspec=garchspec)
		mYs <- melt.dat(Ysim$Y,indiv.names=indiv.names,var.names=var.names,time.index=1:obs,freq=1)
		
		# b/ Rank test
		bsrk <- rank.test(mYs,lags=lags,BS=BS,ddet=det.type,bs.method=bs.method)

		# c/ Check if rk test worked out, if not back to b.
		if(bsrk$xplo!=TRUE){good.rk <- TRUE}
		else{bad.rk <- bad.rk+1}
	}

	# For the FIFO progress tracker
	writeBin(1,f) 



	# Sort out the output, add nbr bad ranks.
	return(list('stat'=bsrk$stat,'pval'=bsrk$pval,
		    'avg.lr'=bsrk$avg.lr,'avg.ipv'=bsrk$avg.ipv,'orig.lr'=bsrk$orig.lr,
		    'trace'=bsrk$bs.rk$trace,'roots'=do.call(cbind,bsrk$roots),
		    'bad'=bad.rk,
		    'residuals'=bsrk$residuals,'coefficients'=bsrk$coefficients))
}
