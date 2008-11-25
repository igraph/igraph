
"scg" <- function(X, ev, nt, matrix=c("symmetric","laplacian","stochastic"),
						algo=c("optimum","interv_km","interv","exact_scg"),
						norm="row", sparse="TRUE", 						direction=c("default","left","right"), evec=NULL, p=NULL,
						gr=NULL, use.arpack=FALSE,maxiter=100,
						output=c("default","matrix","graph"),
						semproj=FALSE, epairs=FALSE, stat.prob=FALSE,...)
{
	if( is.matrix(X) || is(X,"Matrix") ){
		if(nrow(X)!=ncol(X)) stop("'X' must be square")
		if(!is.numeric(X) && (!is(X,"Matrix")) ) stop("'X' must be of numberic mode")
		n <- ncol(X)
	}
	else{
		if(is.igraph(X)) n <- vcount(X)
		else stop("'X' must be of type matrix, Matrix or igraph")
	}
	ev <- as.integer(unique(ev))
	if(any(ev<=0))
		stop("'ev' must be a vector of positive integers")
	if(max(ev)>n)
		stop("'max(ev)' must be smaller than vcount(X) or nrow(X)")
	nt <- as.integer(nt)
	nnt <- length(nt)
	nev <- length(ev)
   if (nnt != nev && nnt != 1) 
        stop("expected length(nt)=1 or length(nt)=length(ev)")
    if (any(nt >= n) || any(nt <= 1)) 
        stop("'nt' must be a vector of integers smaller than n and greater than 1")
    if (nnt == 1) 
        nt <- rep(nt, nev)
	matrix <- match.arg(matrix)
	nmatrix <- switch(matrix, "symmetric" = 1, "laplacian" = 2, "stochastic" = 3)
	if(!is.null(evec)){
		evec <- as.matrix(evec)
		if(ncol(evec)!=nev && nrow(evec)!=n)
			stop("'evec' must have length(ev) columns and a number of rows equal to vcount(X) or nrow(X)")
	}
	if( nmatrix==3 && length(p)!=n && !is.null(p) ) 
		stop("'p' must be null or length(p) = vcount(X) or length(p) = nrow(X)")
		
	if(!is.null(gr)) gr <- renumber(gr)
	if( length(gr)!=n && !is.null(gr) ) 
		stop("'gr' must be null or length(gr) = vcount(X) or length(gr) = nrow(X)")
		
#------retrieve the matrix to be coarse-grained-----
	isGraphX <- is.igraph(X)
	if(isGraphX)
		switch(nmatrix,{
			 if(isGraphX)
			 X <- adjacencyMatrix(X,...)},
			{if((isGraphX && output=="default") || output=="graph"){
			 	X <- adjacencyMatrix(X,...)
			 	#D needed to make the coarse-grained graph
			   	if(norm=="row") D <- Diagonal(x=rowSums(X))
			   	else D <- Diagonal(x=colSums(X))
			   	X <- laplacianMatrix(X,norm=norm)
			 	}
			 else X <- laplacianMatrix(X,norm=norm,...)},
			{X <- stochasticMatrix(X,norm=norm,...)})
	sym <- isSymmetric(X)
#------compute eigenpairs if not supplied----------
	isTransposedX <- direction=="left"
	isTransposedX <- isTransposedX||(direction=="default" && (nmatrix==2||nmatrix==3) && norm=="col")
	if(isTransposedX) X <- t(X)
	if(is.null(gr)){
		eval <- NULL
		if(is.null(evec)){
			if(use.arpack){
				f <- function(v, extra=NULL) as.vector(X%*%v)
				ev.pos <- cut(ev,breaks = c(1,n/2,n),include.lowest=TRUE,labels=FALSE)
				ev1 <- ev[ev.pos==1]
				ev2 <- ev[ev.pos==2]
				if(length(ev1)>0){
					if(sym)
						opt <- list(n=n, which="LA", nev=max(ev1), ncv=2*max(ev1)+1)
					else
						opt <- list(n=n, which="LM", nev=max(ev1), ncv=2*max(ev1)+2)
					a <- arpack(f, options = opt, sym = sym)
					eval <- a$values
					if(sym)
						oo <- order(eval,decreasing=TRUE)
					else
						oo <- order(Mod(eval),decreasing=TRUE)
					eval <- eval[oo][ev1]
					evec <- a$vectors[,oo][,ev1]
					evec <- as.matrix(evec)
					names(eval) <- ev1
					colnames(evec) <- ev1
				}
				if(length(ev2)>0){
					if(sym)
						opt <- list(n=n, which="SA", nev=max(n-ev2+1), ncv=2*max(n-ev2+1)+1)
					else
						opt <- list(n=n, which="SM", nev=max(n-ev2+1), ncv=2*max(n-ev2+1)+2)
					a <- arpack(f, options = opt , sym = sym)
					if(sym)
						oo <- order(eval)
					else
						oo <- order(Mod(eval))
					eval <- cbind(eval,a$values[oo][n-ev2+1])
					evec <- cbind(evec,a$vectors[,oo][,n-ev2+1])
					evec <- as.matrix(evec)
					names(eval) <- c(names(eval),ev2)
					colnames(evec) <- c(colnames(evec),ev2)
				}
			}
			else{
				e <- eigen(X,symmetric=sym)
				eval <- e$values[ev]
				evec <- as.matrix(e$vectors[,ev])
				names(eval) <- ev
				colnames(evec) <- ev
			}
		}
		else{
			evec <- as.matrix(evec)
			names(evec) <- ev
		}
	}
#------handle complex eigenvectors if any------------
#------for now real and immaginary parts-------------
#------are partitioned the same way------------------
	isCmplxEvec <- apply(evec,2,function(z)any(Im(z)!=0))
	if(any(isCmplxEvec)){
		ncc <- ncol(evec) + sum(isCmplxEvec)
		evec.cmplx <- matrix(nrow=nrow(evec), ncol=ncc)
		nt.cmplx <- numeric(ncc)
		k <- 1
		for(i in 1:ncol(evec)){
			if(isCmplxEvec[i]){
				evec.cmplx[,k] <- Re(evec[,i])
				evec.cmplx[,k+1] <- Im(evec[,i])
				nt.cmplx[k] <- nt[i]
				nt.cmplx[k+1] <- nt[i]
				k <- k+2
			}
			else{
				evec.cmplx[,k] <- Re(evec[,i])
				nt.cmplx[k] <- nt[i]
				k <- k+1
			}
		}
	evec.cmplx <- Re(evec.cmplx)
	}
#------compute the stationary probability p---------
#------if not supplied (stochastic case)------------
  if(is.null(p) && nmatrix==3){
		X <- t(X)
		f <- function(v, extra=NULL) as.vector(X%*%v)
		which <- ifelse(sym,"LA","LM")
		p <- arpack(f,options=list(n=n,which=which,nev=1,ncv=5),sym=sym)$vectors
		p <- as.vector(abs(Re(p)))
		p <- p/sum(p)
		X <- t(X)
	}		
#------work out the groups if not supplied---------
	if(is.null(gr)){
		if(any(isCmplxEvec))
			gr <- grouping(V=evec.cmplx, nt=nt.cmplx, matrix=matrix, algo=algo, p=p, maxiter=maxiter)
		else
			gr <- grouping(V=evec, nt=nt, matrix=matrix, algo=algo, p=p, maxiter=maxiter)
	}
	gr <- as.factor(gr)
#------perform the coarse graining-----------------
	if(isTransposedX) X <- t(X)
	lr <- semiProjectors(gr, matrix=matrix, norm=norm, sparse=sparse, p=p)
	output <- match.arg(output)
	#computes a coarse-grained graph
	if( (isGraphX && output=="default") || output=="graph"){
		switch(nmatrix,
			{Xt <- lr$L %*% tcrossprod(X,lr$R)},
			{Xt <- lr$L %*% tcrossprod(D-X,lr$R)},
			{p.sum <- as.vector(tapply(p,gr,sum))
				as.vector
			 if(norm=="row")
				Xt <- p.sum * lr$L%*%tcrossprod(X,lr$R)
			 else
				Xt <- t(p.sum * tcrossprod(R,lr$L%*%X))})
				
		mode <- ifelse(isSymmetric(Xt),"undirected","directed")
		colnames(Xt) <- NULL
		rownames(Xt) <- NULL
		Xt <- graph.adjacency(Xt,weighted=TRUE,mode=mode)
		#self-loops with zero weight may appear
		Xt <- delete.edges(Xt,which(E(Xt)$weight==0)-1)
		V(Xt)$comp <- groupComposition(gr,0)  
	}
#-------computes a coarse-grained matrix
	else{
		Xt <- lr$L %*% tcrossprod(X,lr$R)
		colnames(Xt) <- NULL
		rownames(Xt) <- NULL
	}
	colnames(lr$L) <- colnames(lr$R) <- NULL
	rownames(lr$L) <- rownames(lr$R) <- NULL
	return(list(
			Xt = Xt,
			gr = as.integer(gr),
			L = if(semproj) lr$L,
			R = if(semproj) lr$R,
			values = if(epairs) eval, 
			vectors = if(epairs) evec,
			p = if(stat.prob) p))
}
