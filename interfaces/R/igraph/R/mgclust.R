gclust.rsvt <- function(glist,r=1,maxsvt=10,nmfout=FALSE,maxit=10000,nmfmethod='lee') 
{
    # rsvt -- repeated singular value thresholding 
    # input could be a list of graphs (igraph model) 
    # or a matrix, where each column corresponds to 
    # a vectorized version of a graph 

	if (is.list(glist)) {
	    Time <- length(glist)
	    if (class(glist[[1]])=="igraph") {
	        glist <- lapply(glist,get.adjacency)
	    }
	    Xorigin <- sapply(glist,as.vector)
	} else {
	    Xorigin = glist
	}
	
	Xraw = Xorigin 
	Xorigin = Xorigin %*% solve(diag(colSums(Xorigin)))

	mysvd = tryCatch(irlba(Xorigin,r,r,maxit=maxit), error=function(e) svd(Xorigin,r,r))
	U = mysvd$u[,1:r,drop=FALSE]
	V = mysvd$v[,1:r,drop=FALSE]
	if(r == 1) {
	    S = matrix(mysvd$d[1:r],1,1)
	} else {
	    S = diag(mysvd$d[1:r])
	}
	
	maxsvt = ifelse(maxsvt > 0,maxsvt, 1)
	for(itr in 1:maxsvt) {
	    X = U %*% S %*% t(V)
	    X[X<0]=0
	    mysvd = tryCatch(irlba(X,r,r,maxit=maxit), error=function(e) svd(X,r,r))
	    UU = mysvd$u[,1:r,drop=FALSE]
	    VV = mysvd$v[,1:r,drop=FALSE]
	    if(r == 1) {
	        SS = matrix(mysvd$d[1:r],1,1)
	    } else {
	        SS = diag(mysvd$d[1:r])
	    }
	    if(norm(UU-U) + norm(VV-V) + norm(SS-S) < 1e-12) {
	        break
	    }
	    U = UU;
	    V = VV;
	    S = SS;
	}
	
	XX = X %*% solve(diag(colSums(X)))
	stash = which(rowSums(XX)==0)
	if(length(stash)>0){
	    XX =  XX[-stash,]
	}

	mynmf = nmf(XX,rank=r,method=nmfmethod)
	WW = matrix(0,nrow=nrow(X),ncol=r)
	if(length(stash)>0){
	    WW[-stash,] = basis(mynmf)
	} else {
	    WW = basis(mynmf)
	}
	
	HH = coef(mynmf)
	if(r > 1) {
	    HH = diag(colSums(WW)) %*% HH
	    WW = round(WW %*% solve(diag(colSums(WW))),12)
	    HH = round(HH %*% solve(diag(colSums(HH))),12)   
	} else {
	    HH = colSums(WW) %*% HH
	    WW = round(WW/colSums(WW),12)
	    HH = round(HH %*% solve(diag(colSums(HH))),12)   
	}
	
	if(nmfout)
	    return(list(nmf=mynmf, W=WW, H=HH, Xorigin=Xraw))
	else 
	    return(list(nmf=NULL,  W=WW, H=HH, Xorigin=Xraw))
}

gclust.app <- function(glist, r=1, nmfout=FALSE, maxit=10000, nmfmethod='lee') 
{
    # app -- apparent (clustering), i.e. no rsvt
    # input could be a list of graphs (igraph model) 
    # or a matrix, where each column corresponds to 
    # a vectorized version of a graph 

	if (is.list(glist)) {
	    Time <- length(glist)
	    if (class(glist[[1]])=="igraph") {
	        glist <- lapply(glist,get.adjacency)
	    }
	    Xorigin <- sapply(glist,as.vector)
	} else {
	    Xorigin = glist
	}
	
	Xraw = Xorigin 
	XX = Xorigin %*% solve(diag(colSums(Xorigin)))
	
	stash = which(rowSums(XX)==0)
	if(length(stash)>0){
	    XX =  XX[-stash,]
	}

	mynmf = nmf(XX,rank=r,method=nmfmethod)
	WW = matrix(0,nrow=nrow(Xraw),ncol=r)
	if(length(stash)>0){
	    WW[-stash,] = cbind(basis(mynmf))
	} else {
	    WW = basis(mynmf)
	}
	HH = coef(mynmf)
	if(r > 1) {
	    HH = diag(colSums(WW)) %*% HH
	    WW = round(WW %*% solve(diag(colSums(WW))),12)
	    HH = round(HH %*% solve(diag(colSums(HH))),12)   
	} else {
	    HH = colSums(WW) %*% HH
	    WW = round(WW /colSums(WW),12)
	    HH = round(HH %*% solve(diag(colSums(HH))),12)   
	}
	
	if(nmfout)
	    return(list(nmf=mynmf,W=WW,H=HH, Xorigin=Xraw))
	else 
	    return(list(nmf=NULL,W=WW,H=HH, Xorigin=Xraw))
}

getAICc <- function(gfit) {

    # input is the output from a call to gclust.rsvt or gclust.app
    
	Xorigin = gfit$Xorigin
	Xmean = gfit$W %*% gfit$H %*% diag(colSums(Xorigin))
	Xmean[Xmean<1e-12] = 0
	
	Phat = gfit$W %*% gfit$H
	Phat[Phat < 1e-12] = 0
	if(ncol(gfit$W) > 1) { 
	    nparams = colSums((1.0*(gfit$W > 0) ))
	    nparams = nparams %*% diag(1/colSums(Xorigin %*% t(gfit$H)))
	    nparams = sum(nparams)
	} else {
	    nparams = colSums((1.0*(gfit$W > 0) ))
	    nparams = nparams/colSums(Xorigin %*% t(gfit$H))
	    nparams = sum(nparams)
	}
	
	zeroprob = which(Phat<1e-12)
	if(any(Xorigin[zeroprob]> 1e-12)) {
	    retval = Inf
	} else {
	    if(length(zeroprob) > 0) {
	        Phat = Phat[-zeroprob]
	    } 
	    retval = 2*(nparams) - 2*sum(Phat*log(Phat))
	}
	
	data.frame(nclust=ncol(gfit$W),  negloglikpart=-2*sum(Phat*log(Phat)),  parampart=2*nparams, AIC=retval) 
}
