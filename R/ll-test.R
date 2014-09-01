#' @description bla!
#' 
#' @details bla, bla?
#'
#'
#' @name ll.test
#' @aliases ll.test
#' @title Computes the common cointegration space estimator of Larsson and Lyhagen 2007 JBES  
#' @author Laurent Callot \email{l.callot@@vu.nl}
#' 
#'
#'
#' @param Y A molten data set containing the 4 columns. indiv: the individual (country, sector,...) name as factor levels. Year: the year of observation. Period: the period (day, week, month, quarter) index as consecutive integers. Variable: the name of the variable. value: the value of the variable at the given time for a given individual. 
#' @param rank The co-integration rank of the system, default 1.
#' @param lags An integer indicating the number of lags of the VAR in level.
#' @param ddet An integer indicating the type of deterministics to use, following the typology by Johansen 1988.
#' @param Y0 Optional: a molten set of global exogenous data to use. Its columns must be identical to those of Y, and have the same number of observations. 
#' @param seasonal Optional, should seasonal dummies be included. Default FALSE. 
#'
#' @return A list with 'est' the individual CVAR and 'trace' the trace test statistics for each model and each hypothesis.  
#'
#'
#'
#'
#'
#'
#' @export
ll.test <- function(Y,rank=1,lags=2,ddet=1,Y0=NULL,seasonal=FALSE){
  
  # init
  N <- length(unique(Y$indiv))
  p <- length(unique(Y$Variable))
  # seasonal frequency 
	season.freq <-ifelse(seasonal, length(unique(Y$Period)) , 1)
 
  
  
  # Initialization of the panel estimators
	# Cast the individual data to matrices:
	aY <- acast(Y, Year + Period ~ indiv + Variable)
  aX <- NULL
	if(!is.null(Y0))aX <- cbind(aX,acast(Y0, Year + Period ~ indiv + Variable))
  nobs <- nrow(aY)
	# Get the data-elements of the CVAR:
	ll.cvar <- .mkCVAR(aY,aX,lags)
	# RRR (eigenvecs in columns)
	ll.eig <- .RRR(ll.cvar,ddet,seasonal,season.freq)
  #Concentrated 1st diff and levels  
  DYll <- ll.eig$R0
  LYll <- ll.eig$R1
  
  
  # Switching algorithm stoping criterions
  cvcrit <- 0.001
  maxit <- 200
  
  # individual model (non-panel)
  llid <- .ll.betaID(DYll,LYll,Y,nobs,lags,rank,ddet,p,seasonal,season.freq)
    
  # Hetero co-integrationt test
  llhe <- .ll.betaHE(DYll,LYll,nobs,N,p,rank,llid$coefficients,cvcrit,maxit)
  
  # Homogenous Co-integration vectors
  llho <- .ll.betaHO(DYll,LYll,nobs,N,p,rank,llhe$coefficients,cvcrit,maxit)
  
  
  # Testing 
  LRccs <- -2*(llho$loglik - llhe$loglik)
  pvccs <- 1- pchisq(LRccs,df= (N-1)*rank*(p-rank))
  
  # returning
  ll <- list('individual'=llid,'homogenous'=llho,'heterogenous'=llhe,'lr'=LRccs,'pv'=pvccs,'df'= (N-1)*rank*(p-rank))
	return(ll)
}



# Homogenous co-integration vectors
.ll.betaHE <- function(DYll,LYll,nobs,N,p,rank,betaU,cvcrit,maxit){
  
  #init
  betaHE<- betaU
  count <- 0
  llold <- dll <- 999
  
  #1 Loop unit convergence
  while( ((count < 2)|(abs(dll) > cvcrit)) & (count < maxit))
  {
    #I Loop around individuals
    for(i in 1:N){
      
      #A Construct the matrix of concentration parameters
      bconc <- as.matrix(bdiag(betaHE))
      bconc <- bconc[,-c((1+(i-1)*rank):(i*rank))] 
      
      #B Concentrate co-integration vectors
      dc <- lm.fit(y=DYll,x=LYll%*%bconc)$residuals
      lc <- lm.fit(y=LYll[, c((1+(i-1)*p):(i*p))],x=LYll%*%bconc)$residuals
      
      nobs <- length(lc)
      #C estimate beta i
      S01  <- t(dc)%*%lc/nobs
      S00  <- t(dc)%*%dc/nobs
      S11  <- t(lc)%*%lc/nobs
      S00i <- solve(S00)
      
      # Solving
      eig  <- geigen(A=t(S01) %*% S00i %*% (S01),B=S11)  
      evo <- order(abs(eig$values),decreasing=TRUE)
      betai<- eig$vectors[ ,evo[1:rank]]  
      #D normalize betai and update betaHE
      # normalize beta
      if(rank>1) bnorm <- matrix(rep(betai[1,],each=p),ncol=rank,nrow=p)
      if(rank==1) bnorm <- rep(betai[1],each=p)
      betai <- betai/bnorm     
      
      betaHE[[i]] <- betai
    }
    #II Compute likelihood  
    reg <- lm.fit(y=DYll , x=LYll %*% as.matrix(bdiag(betaHE)))
    llnew <- -(nobs/2)*log(det((var(reg$residuals))))
    
    dll <- llnew - llold
    llold <- llnew
    count <- count+1
  }
  
  #2  Convergence status
  cvhe <- 'nocv'
  if(abs(dll)>cvcrit){if(count>=maxit) cvhe <- 'maxit'}
  if(abs(dll)<cvcrit) cvhe <- 'cv'
  
  return(list('count'=count,'cv'=cvhe,'loglik'=llnew,'coefficients'=betaHE,'alpha'=reg$coefficients))
}
  

# Homogenous co-integration vectors
.ll.betaHO <- function(DYll,LYll,nobs,N,p,rank,betaHE,cvcrit,maxit){
  
  # structure matrix
  H <- NULL
  d <- matrix(0,nrow=N,ncol=1)
  for(i in 1:N){
    di <- d
    di[i] <- 1
    h <- kronecker(kronecker(diag(p),di),diag(rank))
    H <- rbind(H,h)
  }  
  
  # initial value 
  betainit <- betaHO <- Reduce('+',betaHE) / length(betaHE)
    
  # algo initialization
  Dll <- 1
  llp <- 0
  count <- 0
  
  # Switching algo
  while( ((abs(Dll) > cvcrit) | (count < 2)) & (count < maxit) )
  {
    # Regression conditional on beta
    reg1 <- lm.fit(y=DYll,x=LYll%*%kronecker(diag(N),betaHO))
    ll1  <- -(nobs/2)*log(det((var(reg1$residuals))))
      
    # preping the second regression
    omegah <- var(reg1$residuals)
    sqOinv <- (chol(solve(omegah)))
    alphat <- reg1$coefficients%*%sqOinv
    
    # Second regression 
    vecDY <- matrix(DYll%*%sqOinv,ncol=1,byrow = TRUE)
    reg2 <- lm.fit(y=vecDY,x= kronecker(LYll,t(alphat))%*% H)
    res2 <- matrix(reg2$residuals,ncol=N*p,byrow = TRUE) %*% chol(omegah)
    ll2  <- -(nobs/2)*log(det((var(res2))))    
    betaHO  <- matrix(reg2$coefficients,ncol=rank)
    
    # normalize beta
    if(rank>1) bnorm <- solve(betaHO[1:rank,1:rank])
    if(rank==1) bnorm <- solve(betaHO[1])
    betaHO <- betaHO%*%bnorm   
    
    
    Dll <- llp - ll2
    count <- count + 1
    llp <- ll2
  }
  
  # Regression conditional on beta
  reg1 <- lm.fit(y=DYll,x=LYll%*%kronecker(diag(N),betaHO))
  llf  <- -(nobs/2)*log(det((var(reg1$residuals))))
  
  # Strorage
  llccs <- list('coefficients'=betaHO,'count'=count,'loglik'=llf,'alpha' = reg1$coefficients,'init'=betainit)
  
  if(abs(Dll)<cvcrit) llccs$cv <- 'cv'
  else if(count>=maxit) llccs$cv <- 'maxit'
  
  #2  Convergence status
  llccs$cv <- 'nocv'
  if(abs(Dll)>cvcrit){if(count>=maxit) llccs$cv <- 'maxit'}
  if(abs(Dll)<cvcrit) llccs$cv <- 'cv'
  
  return(llccs)
}  

.ll.betaID <- function(DYll,LYll,Y,nobs,lags,rank,ddet,p,seasonal,season.freq,Y0=NULL){
  
  #List of individual data matrices
  betaU <- list()
  llid <- 0
  alphaid <- list()
  icount <- 0
  for(i in unique(Y$indiv)){
    icount <- icount + 1
    # Get data for indiv i
    Yi <- subset(Y,indiv==i)
    
    # Cast the individual data to matrices:
    aYi <- acast(Yi,  Year + Period ~ Variable)
    aXi <- NULL
    if(!is.null(Y0))aXi <- cbind(aXi,acast(Y0, Year + Period ~ Variable))
    
    # Get the data-elements of the CVAR:
    cvari <- .mkCVAR(aYi,aXi,lags)
    # RRR (eigenvecs in columns)
    ll.eig <- .RRR(cvari,ddet,seasonal,season.freq)
    evo <- order(abs(ll.eig$values),decreasing=TRUE)
    betai<- ll.eig$vectors[ ,evo[1:rank]]  
    
      #D normalize betai and update betaHE
    # normalize beta
    if(rank>1) bnorm <- solve(betai[1:rank,1:rank])
    if(rank==1) bnorm <- solve(betai[1])
    betai <- betai%*%bnorm   
    
    betaU[[i]] <- betai
    
    # Regression conditional on beta
    reg1 <- lm.fit(y=DYll[,(1 + p*(icount-1)) :(icount*p)],x=LYll[,(1 + p*(icount-1)) :(icount*p)]%*%betai)
    llid  <- llid + (-nobs/2)*log(det((var(reg1$residuals))))
    alphaid[[i]] <- reg1$coefficients
  }
  
   
  
  unrestricted <- list('coefficients' = betaU,'loglik'= llid,'alpha'=bdiag(alphaid))
  return(unrestricted)
}


