
"semiProjectors" <- function(gr,matrix=c("symmetric","laplacian","stochastic"),
								norm="row",sparse=TRUE,p=NULL)
{
	gr <- renumber(gr)
	#same skeleton for every semi-projector
	L <- as(as.factor(gr),"sparseMatrix")
	if(norm !="row" && norm !="col")
		stop("'norm' expects either 'row' or 'col'")

	matrix <- match.arg(matrix)
	nmatrix <- switch(matrix, "symmetric" = 1, "laplacian" = 2, "stochastic" = 3)
	if(nmatrix==3 && length(p)!=length(gr))
		stop("when 'stochastic' is required provide the stationary probability 'p'")
	switch(nmatrix,
	{	
		L <- L/sqrt(as.vector(table(gr)))
		R <- L	
	},
	{
		R <- L
		if(norm == "row")
			L <- L/as.vector(table(gr))
		else
			R <- R/as.vector(table(gr))
	},
	{
		gr.ind <- split(1:length(gr),gr)
		p.sum <- tapply(p,gr,sum)
		k = 1
		for(gr in gr.ind){
			p[gr] <- p[gr]/p.sum[k]
			k <- k+1
		}
		R <- L
		if(norm == "row")
			L <- t( p*t(L) ) 
		else
			R <- t( p*t(R) ) 
	})
	if(sparse)
			return( list(L = L, R = R) )
	 return( list(L = as.matrix(L), R = as.matrix(R)) )
}
