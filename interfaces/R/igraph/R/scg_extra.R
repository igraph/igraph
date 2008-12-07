
## "normEps" <- function(V, gr, matrix = c("symmetric","laplacian","stochastic"), 
## 							norm = "row", p = NULL)
## {
## 	V <- as.matrix(V,mode="numeric")
## 	gr <- renumber(gr)
## 	if(length(gr)!=nrow(V))
## 		stop("length(gr) and nrow(V) expected to be equal")
## 	lr <- semiProjectors(gr,matrix=matrix,norm=norm,sparse=TRUE,p=p)
## 	P <- crossprod(lr$R,lr$L)
## 	normEpsCalc <- numeric(ncol(V))
## 	for(k in 1:ncol(V))
## 			normEpsCalc[k] <- sqrt(sum( (V[,k] - P%*%V[,k])^2 ))
## 	return(normEpsCalc)
## }

## "groupComposition" <- function(gr, vertex.smallest = 1)
## {
## 	gr <- renumber(gr)
## 	n <- length(gr)
## 	if(n<1) stop("'gr' empty")
## 	return(split(seq(vertex.smallest,n+vertex.smallest-1),as.factor(gr)))
## }	

## "renumber" <- function(x, smallest = 1, order = FALSE)
## {
## 	x <- as.vector(x)
## 	n <- length(x)
## 	if(n==1) return(smallest)
		
## 	x <- as.integer(as.factor(x))
## 	if(!order) return(x)
	
## 	y <- integer(n)
## 	k <- smallest-1
## 	for(i in 1:n){
## 		if(y[i]>0)next
## 		k = k+1
## 		y[x==x[i]] <- k	
## 	}
## 	return(y)
## }
