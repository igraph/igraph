
#Algorithm 1:
#D is NxN non-negative weighted matrix
#level is a number between minimum and maximum element in D
#puts a Threshold on D at level and returns all maximal cliques of the binary network

threshold.net=function(D,level)
{
N=dim(D)[1]

D.t=data.matrix(D>=level)

El=which(D.t>0,arr.ind=TRUE)
Dth=graph.edgelist(El, directed=FALSE)

clqt=maximal.cliques(Dth)
pobt=length(clqt)

Bc=matrix(0,N,pobt)
for(j in 1:pobt){Bc[clqt[[j]],j]=1}

W=matrix(0,1,pobt)
for(j in 1:pobt){W[j]=sum(Bc[,j])}
WS=sort(W,decreasing = TRUE,index.return = TRUE)
Bc=Bc[,WS$ix[1:length(WS$ix)]]
if(pobt==1){Bc=t(Bc)}
;Bc
}


#Algorithm 2:
#D is NxN non negative weighted matrix
# projects the network on the basis elements(Bc) and finds Mu

project.net=function(D,Bc,BB,Mu,iter)
{
K=dim(Bc)[2]
Qt=list
a=0
for (j in 1:K)
        {
        a[j]=sum(Bc[,j])^2
        }

for (i in 1:iter)
{
Dh=Bc%*%diag(Mu)%*%t(Bc)
        for (j in 1:K)
        {
        Qt=BB[which(BB[,2]==j),1]
        Mu[j]=Mu[j]*sum(D[Qt,Qt]/(Dh[Qt,Qt]+.0001))/a[j]
        }
}
;Mu
}


#Main:
#Finds all basis(Bc) for the network D (using algorithm 1), projects on the Bc and finds Muc (using algorithm 2). Input D must be symettric with all nonnegative elements

bsd.net=function(D,iter)
{  
print("Summary of the network:")
Summary.net(D)
print("Bit String Decomposition:")
N=dim(D)[1]
Bc=matrix(0,N,1)
ls=sort(D,decreasing=TRUE)
ldif=diff(ls)<0
L=ls[which(ldif)]
L=c(L,0.0001)
for(level in L)
{
Bt=threshold.net(D,level)
pobt=dim(Bt)[2]
W=matrix(0,1,pobt)
for(j in 1:pobt){W[j]=sum(Bt[,j])}
WS=sort(W,decreasing = TRUE,index.return = TRUE)
Bt=Bt[,WS$ix[1:length(WS$ix)]]
#Bt=Bt[,which(WS$x>1)]

Wc=matrix(0,1,dim(Bc)[2])
for(j in 1:dim(Bc)[2]){Wc[j]=sum(Bc[,j]*(1:N))}

indic=length(Bt)/N

if(indic>1)
{
for(j in 1:dim(Bt)[2])
{
flag=0
for(iij in which(Wc==sum(Bt[,j]*(1:N))))
{
if (sum(Bc[,iij]==Bt[,j])==N){flag=1}
}
if((flag==0) & (sum(Bt[,j])>1) ){Bc=matrix(c(Bc,Bt[,j]),nrow=N)}
}
}

}
Bc=Bc[,2:dim(Bc)[2]]
BB=which(Bc==1,arr.ind=TRUE)
Mu=array(1,dim(Bc)[2])
#iter=10000
Mu=project.net(D,Bc,BB,Mu,iter)
Smb=sort(Mu,decreasing=TRUE,index=TRUE)
re=list(Bc[,Smb$ix],Mu[Smb$ix])
names(re)<-c("Bc","Muc")
;re

}
