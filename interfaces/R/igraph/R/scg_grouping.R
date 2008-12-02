
"grouping" <- function(V, nt, matrix = c("symmetric","laplacian","stochastic"),
                       algo = c("optimum","interv_km","interv","exact_scg"), p=NULL, maxiter=100)
{
	V <- as.matrix(V,mode = "numeric")
	n <- nrow(V)
	nev <- ncol(V)
	if(n<1 || nev<1 || n<=nev)
		stop("invalid 'V'")
	V <- as.vector(V)
	nt <- as.integer(nt)
	nnt <- length(nt)
	if(nnt!=nev && nnt!=1)
		stop("length(nt) and ncol(V) expected to be the same")
	if(any(nt>=n) || any(nt<=1))
		stop("'nt' must be a vector of integers smaller than n and greater than 1")
	if(nnt==1) nt <- rep(nt,nev)
	matrix <- match.arg(matrix)
	nmatrix <- switch(matrix,"symmetric" = 1,"laplacian" = 2,"stochastic" = 3)
	algo <- match.arg(algo)
	nalgo <- switch(algo,"optimum" = 1,"interv_km" = 2,"interv" = 3,"exact_scg" = 4)
	#for now 'p' only required with algo="optimum", but this might change in the future 
	if(nmatrix==3 && length(p)!=n) 
		stop("invalid 'p'")
	#if(!is.loaded("scg_r_wrapper"))
	#	dyn.load("~/Desktop/SCG_project/SCG/src/scg_r_wrapper.so")
	out <- .C("scg_r_wrapper", as.numeric(V), gr=integer(n),as.integer(n),
				as.integer(nt), as.integer(nev),as.integer(nmatrix),
				as.integer(nalgo), as.numeric(p),as.integer(abs(maxiter)),PACKAGE = "SCG")
	return(gr=out$gr)
}
