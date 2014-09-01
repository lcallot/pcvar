
# Monte Carlo internals:
.mcccs.i <- function(mc,MC,obs,N,nvar,BS,rank,
  	    W,Alpha=NULL,Beta=NULL,Lambda0=NULL,Gammal=NULL,
		    Omega,err.dist='gaussian',
		    det.type=1,cdet.load=c(1,1)
		    ,bs.method='resample',t.df=1,res.dep='iid',garchspec=NULL,burn.smpl=10,llest=FALSE){

	# Initialisation
	set.seed(mc*123)
	good.ccs <- FALSE
	bad.ccs <- 0
	lags <- length(Gammal)+1
  
	# Working until we get a non-explosive bootstrap DGP
	while(!good.ccs){
		# a/ Generate data
		Ysim <- gen.pcvar(obs=obs,N=N,nvar=nvar,W=W,Alpha=Alpha,Beta=Beta,Lambda0=Lambda0,Gammal=Gammal,Omega=Omega,err.dist=err.dist,t.df=t.df,burn.smpl=burn,res.dep=res.dep,garchspec=garchspec)
		mYs <- meltdat(Ysim$Y,indiv.names=indiv.names,var.names=var.names,time.index=1:obs,freq=1)
		
		# b/ CCS test
		bsccs <- ccs.test(mYs,W=W,rank=rank,lags=lags,BS=BS,ddet=det.type,bs.method=bs.method)
    
    if(obs< N*lags*nvar+2){llest <- FALSE;warning('Not enough observations, Larson-Lyhagen estimator not computed.')}
    
    if(llest){
      bsll <- ll.test(mYs,rank=rank,lags=lags,ddet=det.type)
      bsccs$ll <- bsll
    }
    else bsccs$ll <- NULL
    
		# c/ Check if ccs test worked out, if not back to b.
		if(bsccs$xplo!=TRUE){good.ccs <- TRUE}
		else{bad.ccs <- bad.ccs+1}
    good.ccs<-TRUE
	}

  # Store number of bad models
  bsccs$bad <- bad.ccs

  return(bsccs)  
}






# Aggregate the MC results
.agg.mc.ccs <- function(mcccs,llest){

  MC <- length(mcccs)

  # The storage list
	aggmc <- list('ppv'=NULL,'plr'=NULL,
      		      'bspv'=NULL,'bsipv'=NULL,'pv'=NULL,
      		      'bad'=NULL,'residuals'=list(),'coefficients'=list(),'llu'=NULL,'llc'=NULL,'ll'=NULL,'cv'=NULL,'iter'=NULL)
  if(llest) ll <- list('lr'=NULL,'pv'=NULL,
                       'homogenous'=list('cv'=NULL,'loglik'=NULL,'coefficients'=list()),
                       'heterogenous'=list('cv'=NULL,'loglik'=NULL,'coefficients'=list()),
                       'betaU'=list())
  
	for(mc in 1:MC){
		# getting iteration mc
		m <- mcccs[[mc]]
		
		# aggregating the aggregates
		for(s in c('bslr','bspv','ppv','plr','bad','pv','bsipv','llu','llc','cv','iter','bslrpv')) aggmc[[s]] <- rbind(aggmc[[s]],m[[s]]) 
		
		# List of coefficients and residuals (which ones? only CCS?)
		aggmc$coefficients[[mc]] <- m$coefficients
		aggmc$residuals[[mc]] <- m$residuals
    
    # Storing the larsson lyhagen results
    if(llest){
      ll[['betaU']][[mc]] <-m[['ll']][['betaU']]
      for(s in c('lr','pv')){ll[[s]] <- c(ll[[s]],m[['ll']][[s]])}      
      
      for(h in c('homogenous','heterogenous')){
        for(s in c('cv','loglik')){
          ll[[h]][[s]] <- c(ll[[h]][[s]],m[['ll']][[h]][[s]])
        }
        ll[[h]][['coefficients']][[mc]]<- m[['ll']][[h]][['coefficients']]
      }
    }      
	}

  if(llest)aggmc$ll <- ll
  
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
		mYs <- meltdat(Ysim$Y,indiv.names=indiv.names,var.names=var.names,time.index=1:obs,freq=1)
		
		# b/ Rank test
		bsrk <- rank.test(mYs,lags=lags,BS=BS,ddet=det.type,bs.method=bs.method)

		# c/ Check if rk test worked out, if not back to b.
		if(bsrk$xplo!=TRUE){good.rk <- TRUE}
		else{bad.rk <- bad.rk+1}
	}


	# Sort out the output, add nbr bad ranks.
	return(list('stat'=bsrk$stat,'pval'=bsrk$pval,
		    'avg.lr'=bsrk$avg.lr,'avg.ipv'=bsrk$avg.ipv,'orig.lr'=bsrk$orig.lr,
		    'trace'=bsrk$bs.rk$trace,'roots'=do.call(cbind,bsrk$roots),
		    'bad'=bad.rk,
		    'residuals'=bsrk$residuals,'coefficients'=bsrk$coefficients))
}
