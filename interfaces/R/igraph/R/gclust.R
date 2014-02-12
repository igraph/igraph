require(igraph)
require(Matrix)
require(NMF)
require(irlba)




########### source code for gclust 
gclust <- function(glist,r=1,maxsvt=10) 
{
    # Model Based Graph Clutering Approach    
    if (is.list(glist)) {
        Time <- length(glist)
        if (class(glist[[1]])=="igraph") {
            glist <- lapply(glist,get.adjacency)
        }
        Xorigin <- sapply(glist,as.vector)
    }

    mysvd = irlba(Xorigin,r,r)
    U = mysvd$u[,1:r,drop=FALSE]
    V = mysvd$v[,1:r,drop=FALSE]
    if(r == 1) {
        S = matrix(mysvd$d[1:r],1,1)
    } else {
        S = diag(mysvd$d[1:r])
    }

    maxsvt = ifelse(maxsvt > 0,maxsvt,1)

    for(itr in 1:maxsvt) {
        X = U %*% S %*% t(V)
        X[X<0]=0
        mysvd = irlba(X,r,r)
        UU = mysvd$u[,1:r,drop=FALSE]
        VV = mysvd$v[,1:r,drop=FALSE]
        if(r == 1) {
            SS = matrix(mysvd$d[1:r],1,1)
        } else {
            SS = diag(mysvd$d[1:r])
        }
        if(norm(UU-U) + norm(VV-V) + norm(SS-S) <
           1e-12) {
            break
        }
        U = UU;
        V = VV;
        S = SS;
    }        

    npairs = nrow(X)
    ntimes = ncol(X)

    # Poisson Approximation Based Specification (`a la Le Cam's
    # Theorem)
    # N.B. this is subject to modification
    nparams = (npairs-1) * r + (r-1) * ntimes  
    Xmean = X
    Xmean[X>1] = 1
    cXmean = 1 - Xmean
    Z = (Xmean^Xorigin) * (cXmean^(1-Xorigin))
    myAIC = 2*nparams - 2*log(prod(Z))

    XX = X %*% solve(diag(colSums(X)))
    stash = which(rowSums(XX)==0)
    if(length(stash)>0){
        XX =  XX[-stash,]
    }
    mynmf = nmf(XX,rank=r)
    WW = matrix(0,nrow=nrow(X),ncol=r)
    if(length(stash)>0){
        WW[-stash,] = basis(mynmf)
    } else {
        WW = basis(mynmf)
    }
    HH = coef(mynmf)
    HH = diag(colSums(WW)) %*% HH
    WW = round(WW %*%
               solve(diag(colSums(WW))),12)
    HH = round(HH %*%
               solve(diag(colSums(HH))),12)   

    membership <- apply(HH,2,which.max)

    return(list(nmf=mynmf,W=WW,H=HH,membership=membership,
                AIC=myAIC))
}

