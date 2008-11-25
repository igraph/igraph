
"adjacencyMatrix" <- function(X,sparse = TRUE,...)
{
	if(!sparse) return(get.adjacency(X,sparse = sparse,...))
	#to correct a bug in get.adjacency when sparse=TRUE
	A <- get.adjacency(X,sparse = sparse,...)
	diag(A) <- diag(A)/2
	return(A)
}

"laplacianMatrix" <- function(X,norm = c("row","col"),sparse = TRUE,
										normalized = FALSE,...)
{
	if(!is.matrix(X) && !is(X,"Matrix") && !is.igraph(X))
		stop("'X' must be of type matrix, Matrix or igraph")
	if(is.igraph(X))
		X <- get.adjacency(X,type = "both",sparse=sparse,...)
	if(nrow(X)!=ncol(X))
		stop("'X' must be square")
	if(any(X<0))
		stop("'X' must have non-negative entries")
	norm <- match.arg(norm)
	if(normalized){
		if(norm == "row") y <- rowSums(X)
		if(norm == "col") y <- colSums(X)
		#by convention
		y[is.infinite(1/y)] <- 0
		X <- X / sqrt(y)
		X <- t(t(X) / sqrt(y))
		return( Diagonal(nrow(X),1)-X )
	}	
	if(norm == "row")
		return(Diagonal(x=rowSums(X))-X)
	if(norm == "col")
		return(Diagonal(x=colSums(X))-X)
}

"stochasticMatrix" <- function(X,norm = c("row","col"),sparse = TRUE,...)
{
	if( !is.matrix(X) && !is(X,"Matrix") && !is.igraph(X) )
		stop("'X' must be of type matrix, Matrix or igraph")
	if(is.igraph(X))
		X <- get.adjacency(X,type = "both",sparse=sparse,...)
	if(nrow(X)!=ncol(X))
		stop("'X' must be square")
	if(any(X<0))
		stop("'X' must have non-negative entries")
	norm <- match.arg(norm)
	if(norm == "row")
	{
		y <- rowSums(X)
		if(any(y==0))
			stop("zero rows in 'X' not allowed")
		return( X/y )
	}
	if(norm == "col")
	{
		y <- colSums(X)
		if(any(y==0))
			stop("zero columns in 'X not allowed'")
		return( t(t(X) / y) )
	}
}
