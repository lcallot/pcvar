
.ccs.core <- function(Y,W,rank,lags=2,ddet=0,Y0=NULL,seasonal=FALSE,season.freq,ccsinit,exo=TRUE){
  
  # Construct the weighted averages.
	Ystar <- .mkYstar(Y,W)
  
  #List of individual data matrices
  cvarlst <- list()
  for(i in unique(Y$indiv)){
    # Get data for indiv i
    Yi <- subset(Y,indiv==i)
    YSi <- subset(Ystar,indiv==i)
    
    # Cast the individual data to matrices:
    aYi <- acast(Yi,  Year + Period ~ Variable)
    aXi <- acast(YSi,  Year + Period ~ Variable)
    if(!is.null(Y0))aXi <- cbind(aXi,acast(Y0, Year + Period ~ Variable))
    
    # Get the data-elements of the CVAR:
    cvarlst[[i]] <- .mkCVAR(aYi,aXi,lags)
  }
  
  # Calling the optimizer
  ccsopt <- .ccsest(cvarlst,rank,lags,ddet,seasonal,season.freq,ccsinit,exo)
  
  return(list('ccsopt'=ccsopt,'cvarlst'=cvarlst))
}



# Get initial values, average of normalized unretricted bet 
.ccsinit <- function(Y,pcvar.coef,rank,exo=TRUE){
  bhat <- pcvar.coef[[rank+1]][['BETA']] 
  
  npar <- (1+exo)*length(unique(Y$Variable))
  
  
  ccsinit <- matrix(0,ncol=rank,nrow=npar)
  for(i in 1:length(unique(Y$indiv)) )
  {
    # extract the blocks
    bi <- bhat[ (1+npar*(i-1)):(npar*i) , (1+rank*(i-1)):(i*rank) ]  
    # normalize
    
    if(rank==1) bnorm <- solve(bi[1])
    if(rank>1) bnorm <- solve(bi[1:rank,1:rank])
    binorm <- bi%*%bnorm    
    # sum
    ccsinit <- ccsinit + binorm  
  }
  
  ccsinit <- ccsinit / length(unique(Y$indiv))
  
  return(ccsinit)
}

.ccsest <- function(cvarlst,rank,lags=2,ddet=1,seasonal=FALSE,season.freq=NULL,ccsinit=NULL,exo=TRUE){
   
    # Vectorizing the initial values
    ccsinit <- ccsinit[-c(1:rank),]
    ccsinit <- matrix(ccsinit,ncol=1)
    
    # Optmizer controls
    optctrl <- list('maxit' = 200)
    
    # Optmizer call
    ccsopt <- optim(par=ccsinit,fn=.optfun,cvarlst=cvarlst,rank=rank,ddet=ddet,lags=lags,seasonal=seasonal,season.freq=season.freq,exo=exo,
                    method = 'BFGS',control=optctrl) 
    
    #cat('\n Iterations: ', ccsopt$counts[1],'. Convergence: ',ccsopt$convergence,'.',sep='')
    
    optccs <- matrix(ccsopt$par,ncol=rank)
    #optccs <- rbind(matrix(1,nrow=1,ncol=rank),optccs)
    optccs <- rbind(diag(rank),optccs)
    
    ccsopt$ccs   <- optccs
    ccsopt$roots <- list()
    ccsopt$xplo  <- FALSE
    
  return(ccsopt) 
}


# Function to be passed to the optimizer.
# Takes a parameter value and the list of individual data.
# Returns the panel's log-likelihood
.optfun <- function(optpar,cvarlst,rank,ddet,lags,seasonal,season.freq,exo){
  
  optpar <- matrix(optpar,ncol=rank)
  optpar <- rbind(diag(rank),optpar)
  #optpar <- rbind(matrix(1,nrow=1,ncol=rank),optpar)
  
  #padding with exo zeros
  optpar <- rbind(optpar,matrix(0,nrow=nrow(optpar),ncol=rank))
  
  loglik <- 0
  # Computing the indiviual SR params + residuals, the the log likelihood
  for(i in 1:length(cvarlst))
  {
    CVAR    <- cvarlst[[i]]
    cvsr    <- .ccsSR(CVAR,optpar,ddet,lags,seasonal,season.freq,exo)
    loglik  <- loglik - (nrow(cvsr$residuals)/2)*log(det((var(cvsr$residuals))))
  }

  return(-loglik)
}


# Computes the short run parameters of the CVAR under each rank assumption.
.ccsSR <- function(CVAR,optpar,ddet,lags,seasonal,season.freq,exo=TRUE){

  #initialize storage
	sr <- list()

	# Constructing the long run deterministics
	cdet <- NULL
	if(ddet==1)cdet <- matrix(1,nrow=nrow(CVAR$DY),ncol=1)
	if(ddet==3)cdet <- matrix(1:nrow(CVAR$DY),ncol=1)

	# Constructing the short run deterministics
	srdet <- NULL
	if(ddet==2)srdet <- matrix(1,nrow=nrow(CVAR$DY),ncol=1)

	# Constructing the seasonal dummies
	if(seasonal&season.freq>1){
		obs		 <- nrow(CVAR$DY) 
		ncycle		 <- ceiling(obs/season.freq)
		seasD		 <- diag(season.freq)
		seasD[1,1]	 <- 0
		# constructing the (potentially too long) seasonal dummy matrix 
		long.dum	 <- matrix(1,nrow=ncycle,ncol=1)%x%seasD
		#trimming
		seas.dum	 <- long.dum[1:obs,-1]
		# adding to the short run det
		srdet <- cbind(srdet,seas.dum)
	}
  
  r <- ncol(optpar)

	# Compute the Error correction vector(s).
	ECM <- NULL
  if(!exo) optpar <- matrix(optpar[1:ncol(cbind(cdet,CVAR$LZ)),],ncol=r)
  
	if(r>0)ECM <- cbind(cdet,CVAR$LZ)%*%optpar

	# Plain old ols.
	x <- cbind(ECM,srdet,CVAR$DX,CVAR$LDZ)
	if(!is.null(x))sr.cvar <- lm.fit(y=CVAR$DY,x=x)
	else sr.cvar <- list()
  
  pvec <- sr.cvar$coefficients

	# sorting out the parameters:
	coef.tmp <- list('alpha'=NULL,'beta'=NULL,'rho'=NULL,'lambda0'=NULL,'psi'=NULL)

	# Storing the cointegration parameters
	if(r>0){
		coef.tmp$alpha <- (-1)*matrix(t(pvec[1:r,]),ncol=r)
			pvec <- pvec[-c(1:r),]
			
			if(ddet%in%c(1,3)){
			   coef.tmp$beta <- optpar[-1,1:r]
			   coef.tmp$rho  <- optpar[1,1:r]
				}
			else {coef.tmp$beta <- optpar[,1:r]} 
			}
		# Storing the short run deterministics
		if(!is.null(srdet)){
		  coef.tmp$psi <- pvec[1:ncol(srdet),]
			pvec <- pvec[-c(1:ncol(srdet)),]
			}
		# Storing the contemporanous first diff of the exo variables
    if(!is.null(CVAR$DX)){
  		coef.tmp$lambda0 <- t(pvec[1:ncol(CVAR$DX),])
  		pvec <- pvec[-c(1:ncol(CVAR$DX)),]
    }

		# storing the lagged first differences
		if(lags>1){for(l in 2:lags){
			coef.tmp[[paste('lag',l-1,sep='_')]] <- t(pvec[1:(ncol(CVAR$LDZ)/(lags-1)),])
			pvec <- pvec[-c(1:(ncol(CVAR$LDZ)/(lags-1))),]
			}}

		if(length(pvec)>0)stop('there are parameters left!')

		# Storing
		sr <- list('rank'=r,'ECM'=ECM,'coefficients'=coef.tmp,'y'=CVAR$DY,'x'=x,'residuals'=sr.cvar$residuals)
	
	return(sr)
}

