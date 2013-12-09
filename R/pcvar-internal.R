


# Constructs the weighted average of Y using W
.mkYstar <- function(Y,W){

	# The number of countries.
	N  <- length(unique(Y$indiv))

	# Creating equal weights if missing W.
	if(is.null(W)){W <- matrix(1/(N-1),N,N)
		diag(W) <- 0}

	# Initializing Ystar.
	Ystar <- Y
	Ystar$value <- NA

	# In case W is a matrix (or was missing.)
	if(length(dim(W))==2){
		for(vn in unique(Y$Variable)){
			# Selecting the data for variable vn and casting as a matrix. 
			Yvn <- subset(Y,Variable==vn)
			aYvn<- acast(Yvn, Period + Year ~ indiv)

			# Constructing the weighted averages.
			aYvnstar  <- aYvn %*% W
			Yvnstar  <- melt(aYvnstar)

			# Storing 
			Ystar$value[Ystar$Variable==vn] <- Yvnstar$value
		}
	}


	return(Ystar)
}



# Returns the data elements of the individual CVAR (contemporaneous first diff, diff of exo vars, lagged levels, lagged diff of all vars

.mkCVAR <- function(aYi,aXi,lags){

	nr <- nrow(aYi)

	# 1st diff of the dep var (LHS)
	DY <- aYi[-1,]-aYi[-nr,]

	# 1st diff of the exo vars
	DX <- aXi[-1,]-aXi[-nr,]
	
	# Lagged levels
	LZ <- cbind(aYi[-1,],aXi[-1,])

	# Lagged 1st diffs
	LDZ <- NULL
	if(lags>1)for(l in 2:lags){
		LDZl <- cbind(DY,DX)
		LDZl <- LDZl[-(nrow(LDZl)+c((-l+2):0)),] #trimming up here
		if(l<lags)LDZl <- LDZl[-c(1:(lags-l)),] #trimming down there
		LDZ <- cbind(LDZ,LDZl) # merging
	}

	# Trimming the rest:
	if(lags>1){
		DY <- DY[-c(1:(lags-1)),]
		DX <- DX[-c(1:(lags-1)),]
		LZ <- LZ[-c(1:(lags-1)),]		
	}

	return(list('DY'=DY,'DX'=DX,'LZ'=LZ,'LDZ'=LDZ))
}




# Computes a reduced rank regression, resturns the eigenvalue and eigenvectors solving the EV problem. 
.RRR <- function(CVAR,ddet,seasonal,season.freq){
	obs <- nrow(CVAR$DY) 
	
	# Short run dynamics
	Z2 <- cbind(CVAR$DX,CVAR$LDZ)
	# Unrestricted deterministics
	if(ddet%in%c(2,3))Z2 <- cbind(1,Z2)

	# Constructing the seasonal dummies
	if(seasonal&season.freq>1){
		ncycle		 <- ceiling(obs/season.freq)
		seasD		 <- diag(season.freq)
		seasD[1,1]	 <- 0
		# constructing the (potentially too long) seasonal dummy matrix 
		long.dum	 <- matrix(1,nrow=ncycle,ncol=1)%x%seasD
		#trimming
		seas.dum	 <- long.dum[1:obs,-1]
		# adding to Z2
		Z2 <- cbind(Z2,seas.dum)
	}

	R0 <- lm.fit(y=CVAR$DY,x=Z2)$residuals
	R1 <- lm.fit(y=CVAR$LZ,x=Z2)$residuals
	colnames(R1) <- c(colnames(R0),paste(colnames(R0),'s',sep=''))
	if(ddet==1)R1 <- cbind(1,R1)
	if(ddet==3)R1 <- cbind(1:obs,R1)


	S01 <- t(R0)%*%R1/obs
	S00 <- t(R0)%*%R0/obs
	S11 <- t(R1)%*%R1/obs

	S00i <- solve(S00)

	eig <- geigen(A=t(S01) %*% S00i %*% (S01),B=S11)
	eig.orig <- eig
	

	# Take the modulus and sort the vectors
	eig$values <- Mod(eig$values)
	eig$values <- sort(abs(eig$values),decreasing=TRUE)

	# Eigen vectors: in columns.
	eig$vectors <- eig$vectors[,order(abs(eig.orig$values),decreasing=TRUE)]

	return(eig)
}


# Computes the short run parameters of the CVAR under each rank assumption.
.cpSR <- function(CVAR,eig,ddet,lags,seasonal,season.freq){

	#initialize storage
	sr.allrk <- list()

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
	#Loop over the possible ranks. 
	for(r in 0:ncol(CVAR$DY)){

		# Compute the Error correction vector(s).
		ECM <- NULL
		if(r>0)ECM <- cbind(cdet,CVAR$LZ)%*%eig$vector[,1:r]

		# Plain old ols.
		x <- cbind(ECM,srdet,CVAR$DX,CVAR$LDZ)
		sr.cvar <- lm.fit(y=CVAR$DY,x=x)
		pvec <- sr.cvar$coefficients

		# sorting out the parameters:
		coef.tmp <- list('alpha'=NULL,'beta'=NULL,'rho'=NULL,'lambda0'=NULL,'psi'=NULL)

		# Storing the cointegration parameters
		if(r>0){
			coef.tmp$alpha <- (-1)*matrix(t(pvec[1:r,]),ncol=r)
			pvec <- pvec[-c(1:r),]
			
			if(ddet%in%c(1,3)){
			   coef.tmp$beta <- eig$vectors[-1,1:r]
			   coef.tmp$rho<- eig$vectors[1,1:r]
				}
			else {coef.tmp$beta <- eig$vectors[,1:r]} 
			}
		# Storing the short run deterministics
		if(!is.null(srdet)){
		   	coef.tmp$psi <- pvec[1:ncol(srdet),]
			pvec <- pvec[-c(1:ncol(srdet)),]
			}
		# Storing the contemporanous first diff of the exo variables
		coef.tmp$lambda0 <- t(pvec[1:ncol(CVAR$DX),])
		pvec <- pvec[-c(1:ncol(CVAR$DX)),]

		# storing the lagged first differences
		if(lags>1){for(l in 2:lags){
			coef.tmp[[paste('lag',l-1,sep='_')]] <- t(pvec[1:(ncol(CVAR$LDZ)/(lags-1)),])
			pvec <- pvec[-c(1:(ncol(CVAR$LDZ)/(lags-1))),]
			}}

		if(length(pvec)>0)stop('there are parameters left!')

		# Storing
		sr.allrk[[r+1]] <- list('rank'=r,'ECM'=ECM,'coefficients'=coef.tmp,'y'=CVAR$DY,'x'=x,'residuals'=sr.cvar$residuals)
		if(r>0)sr.allrk[[r+1]] <- c(sr.allrk[[r+1]],list('eval'=eig$values[1:r],'evec'=eig$vectors[,1:r]))
	}

	return(sr.allrk)
}


# Takes the raw estimated models, formats and returns the residuals.
.mkresid <- function(pcvar.est){

	res <- list()
	for(r in 1:length(pcvar.est[[1]])){
		# Storing the residuals for rank r
		res[[r]] <- pcvar.est[[1]][[r]]$residuals
		for(i in 2:length(pcvar.est)){res[[r]] <- cbind(res[[r]],pcvar.est[[i]][[r]]$residuals)}
		}
return(res)
}


# Constructs the large weighting matrices to link the panel from the original W matrix. 
.mkW <- function(W,N,nvar){

	# Creating equal weights if missing W.
	if(is.null(W))W='equal'
	if(W=='equal'){
		W <- matrix(1/(N-1),N,N)
		diag(W) <- 0}

	# Constructing the W matrix
	w1 <- matrix(0,nrow=2*N,ncol=N)
	w1[cbind(2*(1:N)-1,1:N)] <- 1
	w1[2*(1:N),] <- W 
	W1 <- kronecker(w1,diag(nvar))

	# Constructing W0
	W0 <- kronecker(W,diag(nvar)) 

	return(list('W0'=W0,'W1'=W1))
}



# Creates the full model parameter matrices from the estimated parameters and the big weight matrices
.mkpcvar <- function(pcvar.est,ddet,lags)
{
	pcvar.coef <- list()
	p <- length(pcvar.est[[1]])-1
	N <- length(pcvar.est)

	# A vec of parameter matrices names
	par.names <- c('alpha','beta','rho','lambda0','psi')
	if(lags>1)par.names <- c(par.names,paste('lag',1:(lags-1),sep='_'))	

	for(r in 1:length(pcvar.est[[1]])){
		pcvar.coef[[r]] <- list()
		for(pn in par.names){
			# create a temporary list with the 'pn' parameter for all countries.
			plst <- list()
			for(i in 1:length(pcvar.est)){plst[[i]] <- pcvar.est[[i]][[r]][['coefficients']][[pn]]}
			if(length(plst)>0){
				assign(toupper(pn),bdiag(plst))
				# assigning it to a block diagonal sparse matrix called PN.
				#if(!(pn%in%c('rho','psi')))assign(toupper(pn),bdiag(plst))
				# Deterministics stacked, not block diagonal.
				#if((pn%in%c('rho','psi')))assign(toupper(pn),as.vector(do.call('cbind',plst)))
				}
			else {assign(toupper(pn),NULL)}
			pcvar.coef[[r]][[toupper(pn)]] <- get(toupper(pn))
			}
		}	
return(pcvar.coef)
}



# Computes the roots of the full pcvar
.roots <- function(pcvar.lev,lags){

	Np <- nrow(pcvar.lev$coefficients[[1]][[1]])

	roots <- list()
	r <- 1

	for(cl in pcvar.lev$coefficients){	
		companion <- matrix(0, nrow = Np*lags, ncol = Np*lags)
		companion[1:Np, 1:Np] <- as.matrix(cl$L1)
		if(lags > 1){
			for(l in 2:lags){companion[1:Np,(1+Np*(l-1)):(Np*l)] <- as.matrix(cl[[paste('L',l,sep='')]])}
			j <- 0
			for( i in (Np+1):(lags*Np)){
				j <- j + 1
				companion[i, j] <- 1
				}
			}
		roots[[r]] <-  Mod(eigen(companion,only.values=TRUE)$values)[1:Np]
		r <- r+1
		}
	return(roots)
}

# Tranform the estimated parameters of the PCVAR oto the model in level. used to compute roots and generate pseudo data. 
.mk.pcvar.level <- function(pcvar.coef,res,lags,bigW){

	W0 <- bigW$W0
	W1 <- bigW$W1
	Np <- nrow(W0)

	coef.lev <- list()
	res.lev <- list()

	r <- 1
	for(cf in pcvar.coef){
		coefL <- list()
		# constructing the A matrix (eq 10 rank paper)
		A <- diag(Np)-(as.matrix(cf$LAMBDA0))%*%W0
		Ai <- solve(A)
		# The parameter matrix for the first lag
		L1 <- diag(Np) # if the rank is 0. 
		if(!is.null(cf$BETA)) L1 <- Ai%*%(A+as.matrix(cf$ALPHA) %*% t(as.matrix(cf$BETA)) %*% W1)

		coefL$L1 <- L1

		# Constructing the following lag matrices if needed.
		if(lags>1){
			for(l in 2:lags){
				coefL[[paste('L',l-1,sep='')]] <- coefL[[paste('L',l-1,sep='')]] + Ai %*% (as.matrix(cf[[paste('LAG_',l-1,sep='')]])) %*% W1
				coefL[[paste('L',l,sep='')]] <-  Ai %*% -(as.matrix(cf[[paste('LAG_',l-1,sep='')]])) %*% W1
				}
			}

		if(!is.null(cf$RHO))coefL$RHOdet <- Ai %*% cf$ALPHA %*% (cf$RHO) 
		if(!is.null(cf$PSI))coefL$PHIdet <- Ai %*% t(as.matrix(cf$PSI)) 

		coef.lev[[r]] <- coefL
		res.lev[[r]] <- Ai%*%t(res[[r]])
		r <- r+1
		}
	return(list('coefficients'=coef.lev,'residuals'=res.lev))
}


# Bootstrap the cointegration rank test:
.bs.rank <- function(pcvar.lev,BS=99,model,bs.method='resample',lags,Y,W,ddet,rts,cores){

	# Storage
	bs <- list()
	bs$trace <- array(NA,c(dim(model$trace),BS+1))
	bs$trace[,,1] <- model$trace
	bs$eval <- array(NA,c(dim(model$eigenvalues),length(pcvar.lev$residuals),BS+1))

	# 0/Loop over the ranks:
	nrk <- length(pcvar.lev$residuals)
	for(r in 1:nrk){
		bs$eval[,,r,1] <- model$eigenvalues

		# 1/ Recenter the residuals
		res <- pcvar.lev$residuals[[r]]-colMeans(pcvar.lev$residuals[[r]])
		obs <- ncol(res)
		Np <- nrow(res)

		# Creating a bootstrap function, will be parallelized or not below.
		for(b in 1:BS){
			# 2/ make pseudo residuals
			if(bs.method=='resample')bsres <- res[,ceiling(obs*runif(obs+lags))]
			if(bs.method=='wild')	bsres <- cbind(matrix(0,ncol=lags,nrow=nrow(res)),res*matrix(rnorm(length(res)),nrow=nrow(res),ncol=ncol(res)))

			# 3/ generate pseudo data
			Y.bs <- matrix(0,ncol=obs+2*lags,nrow=Np)
			for(o in (lags+1):(obs+2*lags)){

				# Summing over the lags:
				for(l in 1:lags)Y.bs[,o] <- Y.bs[,o]+as.matrix(pcvar.lev$coefficients[[r]][[paste('L',l,sep='')]]%*%Y.bs[,o-l]) 
				# Adding noise:
				Y.bs[,o] <-Y.bs[,o]+ bsres[,o-lags]
				# No deterministics, no intial values as per the paper
			
			}
			# 4/ format the data matrix to molten format
			mYbs <- Y
			mYbs$value <- melt(t(Y.bs[,-c(1:lags)]))$value 

			# 5/ compute the rank test on the pseudo data
			bs.trtest	 <- trace.test(mYbs,W,lags,ddet)
			bs$trace[,nrk-r+1,b+1] <- bs.trtest$trace[,nrk-r+1]
			bs$eval[,,1+nrk-r,b+1] <- bs.trtest$eigenvalues

		}
	}
	# 6/ return trace test for each iteration

	return(bs)
}








msgProgressBar <- function (min = 0, max = 1, initial = 0, char = "=", width = NA, title, label, style = 1, file = "") 
{
	if (!identical(file, "") && !(inherits(file, "connection") && isOpen(file))) stop("'file' must be \"\" or an open connection object")
    
	if (!style %in% 1L:3L) style <- 1

    	.val	 <- initial
    	.killed	 <- FALSE
    	.nb	 <- 0L
    	.pc	 <- -1L
    	nw	 <- nchar(char, "w")
    
    if (is.na(width)) {
        width <- getOption("width")
        if (style == 3L) width <- width - 10L
        width <- trunc(width/nw)
    }

    if (max <= min)    stop("must have 'max' > 'min'")

    # 1 up function for each of the three styles.

    up1 <- function(value) {
        if (!is.finite(value) || value < min || value > max) 
            return()
        .val <<- value
        nb <- round(width * (value - min)/(max - min))
        if (.nb < nb) {
            (cat(paste(rep.int(char, nb - .nb), collapse = ""), file = file))
            flush.console()
        }
        else if (.nb > nb) {
            (cat("\r", paste(rep.int(" ", .nb * nw), collapse = ""), 
                "\r", paste(rep.int(char, nb), collapse = ""), 
                sep = "", file = file))
            flush.console()
        }
        .nb <<- nb
    }

    up2 <- function(value) {
        if (!is.finite(value) || value < min || value > max) 
            return()
        .val <<- value
        nb <- round(width * (value - min)/(max - min))
        if (.nb <= nb) {
            (cat("\r", paste(rep.int(char, nb), collapse = ""), 
                sep = "", file = file))
            flush.console()
        }
        else {
            (cat("\r", paste(rep.int(" ", .nb * nw), collapse = ""), 
                "\r", paste(rep.int(char, nb), collapse = ""), 
                sep = "", file = file))
            flush.console()
        }
        .nb <<- nb
    }


    up3 <- function(value) {
        if (!is.finite(value) || value < min || value > max) return()
        .val <<- value
        nb <- round(width * (value - min)/(max - min))
        pc <- round(100 * (value - min)/(max - min))
        if (nb == .nb && pc == .pc) return()
        (cat(paste(c("\r  -", rep.int(" ", nw * width + 6)), collapse = ""), 
            file = file))
        (cat(paste(c("\r  |", rep.int(char, nb),
			    rep.int(" ", nw * (width - nb)),
			    sprintf("| %3d%%", pc)), collapse = ""), file = file))
        flush.console()
        .nb <<- nb
        .pc <<- pc
    	}

    getVal <- function() .val

    kill <- function() if (!.killed) {
        message(cat("\n", file = file))
        flush.console()
        .killed <<- TRUE
    }
    up <- switch(style, up1, up2, up3)
    up(initial)
    structure(list(getVal = getVal, up = up, kill = kill), class = "msgProgressBar")
}

setMsgProgressBar <- function (pb, value, title = NULL, label = NULL) 
{
    if (!inherits(pb, "msgProgressBar")) 
        stop(gettextf("'pb' is not from class %s", dQuote('msgProgressBar')), 
            domain = NA)
    oldval <- pb$getVal()
    pb$up(value)
    invisible(oldval)
}

close.msgProgressBar <- function (con, ...) 
{
    con$kill()
    invisible(NULL)
}

