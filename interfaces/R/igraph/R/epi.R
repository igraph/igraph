
#   IGraph R package
#   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
#   334 Harvard street, Cambridge, MA 02139 USA
#   
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
#   02110-1301 USA
#
###################################################################

# (c) Eric Kolaczyk 2014

# R function to simulate (N_S,N_I,N_R) curves on a network.
# Inputs: gg   := network graph
#            beta     := infection rate
#            gamma := recovery rate 
#               ns := number of epidemic processes to simulate                   
# Outputs: sim.res := list of simulation results, each element of which is composed of 
#         times := vector of times at which transitions occur
#         NS   := vector of values N_S(t), at time points in tvec
#         NI    := vector of values N_I(t), at time points in tvec
#         NR   := vector of values N_R(t), at time points in tvec
#
# Notes: Requires igraph.
simul.net.epi <- function(gg,beta,gamma,ns){
  nv <- vcount(gg)
  sim.res <- vector("list",ns)

  for(j in (1:ns)){
    # Randomly place an infected and initialize status vector.
    initialinf <- sample(seq(1,nv),1)
    status <- rep("S",nv)
    status[initialinf] <- "I"
    NSvec <- nv-1
    NIvec <- 1
    NRvec <- 0
    tvec <- 0
 
    # Set up initial rate vector.
    lamvec <- numeric(nv)
    lamvec[initialinf] <- gamma
    lamvec[V(gg)[nei(initialinf)]] <- beta
    lam <- sum(lamvec)
  
    # Simulate SIR process on graph.
    NS <- nv-1
    NI <- 1
    NR <- 0
    cumtt <- 0
    while(lam > 0){
      # Simulate an event.
      tt <- rexp(1,lam)
      pvec <- lamvec/lam
      vchange <- sample(nv,1,replace=T,prob=pvec)
      if(status[vchange]=="R"){print("WHAT THE HECK? SHOULD NOT HAPPEN!")}
      if(status[vchange]=="I"){status[vchange] <- "R"
                                         NI <- NI - 1
                                         NR <- NR + 1}
      if(status[vchange]=="S"){status[vchange] <- "I"
                                          NS <- NS-1
                                          NI <- NI+1}
      cumtt <- cumtt + tt
    
      # Update time series.
      NSvec <- c(NSvec,NS)
      NIvec <- c(NIvec,NI)
      NRvec <- c(NRvec,NR)
      tvec <- c(tvec,cumtt)
    
      # Update vector of rates and total rate.
      lamvec[status=="I"] <- gamma
      lamvec[status=="R"] <- 0
      susceps <- seq(1,nv)[status=="S"]
      infec.ind <- as.numeric(status=="I")
      if(length(susceps)>0){
        for(i in seq(1,length(susceps))){
          snbrs <- V(gg)[nei(susceps[i])]
	  if(length(snbrs)>0){lamtmp <- beta*sum(infec.ind[snbrs])}
          else{lamtmp <- 0}
          lamvec[susceps[i]] <- lamtmp
        }
      }
      lam <- sum(lamvec)  
    }
    SIRtmp <- list(tvec,NSvec,NIvec,NRvec)
    names(SIRtmp) <- c("times","NS","NI","NR")
    sim.res[[j]] <- SIRtmp
  }
  return(sim.res)
}

# R function to plot compartment total curves from simul.net.epi .
# Inputs:  sim.res :=  list of simulated network SIR processes
#             comp := compartment (i.e., "NS", "NI", or "NR")
#                q := vector of lower and upper quantiles, resp
#             cols := char vector of colors for lines, median, and quantiles, resp.
# Outputs:  None.  Just produces the plot of all compartment curves, 
#           with median and quantiles.

plot.net.epi <- function(sim.res,comp,q,cols){
  
  if(missing(comp))
     comp <- "NI"
  if(missing(q))
     q <- c(0.1,0.9)
  if((missing(cols)) & (comp=="NI"))
     cols <- c("skyblue","blue","blue")
  if((missing(cols)) & (comp=="NS"))
     cols <- c("pink","red","red")
  if((missing(cols)) & (comp=="NR"))
     cols <- c("palegoldenrod","gold","gold")

  ns <- length(sim.res)
  max.x <- max(unlist(lapply(sim.res,function(x){max(x$times)})))
  max.y <- max(unlist(lapply(sim.res,function(x){max(x[[comp]])})))

  ## Work out median and quantile vectors.
  # Create single huge vector of times and compartment counts.
  big.time <- unlist(sapply(sim.res,function(x){x$times}))
  big.N <- unlist(sapply(sim.res,function(x){x[[comp]]}))
  # Adhoc use of Freedman-Diaconis binwidth; rescale time accordingly.
  w <- 1/(2*(quantile(big.time,0.75)-quantile(big.time,0.25))/(100^(1/3)))
  time.bin <- floor(big.time*w)
  
  # Compute median and lower/upper quantiles.
  my.m <- tapply(big.N,time.bin,median)
  my.ql <- tapply(big.N,time.bin,function(x){quantile(x,prob=q[1])})
  my.qu <- tapply(big.N,time.bin,function(x){quantile(x,prob=q[2])})

  
  ## Generate the plot, first with individual curves, and then 
  ## adding median and quantile curves.
  #Clean the frame, set up axes, and plots individual curves.
  par(new=F)
  par(oma=c(1,1.2,1,1))

  bin.vals <- as.numeric(levels(as.factor(time.bin)))/w

  if(comp=="NI"){my.ylab <- expression(N[I](t))}
  if(comp=="NR"){my.ylab <- expression(N[R](t))}
  if(comp=="NS"){my.ylab <- expression(N[S](t))}

  # Plot the stochastic curves individually.
  plot(0,0,type="n",xlim=c(0,max.x),ylim=c(0,max.y),
       xlab="Time",ylab=my.ylab)
  lapply(seq_along(sim.res), 
       function(i) lines(sim.res[[i]]$time,sim.res[[i]][[comp]],col=cols[1]))

  # Plot the median and quantiles.
  par(new=T)
  plot(bin.vals,my.m,type="l",lwd=2,col=cols[2],
       xlim=c(0,max.x),ylim=c(0,max.y),xlab="",ylab="")
  par(new=T)
  plot(bin.vals,my.ql,type="l",lty=3,lwd=2,col=cols[3],
       xlim=c(0,max.x),ylim=c(0,max.y),xlab="",ylab="")
  par(new=T)
  plot(bin.vals,my.qu,type="l",lty=3,lwd=2,col=cols[3],
       xlim=c(0,max.x),ylim=c(0,max.y),xlab="",ylab="")

}

function() {
  ## Code to illustrate epidemic processes, w/out and w/ networks.
  
  library(igraph)
  source("simul.net.epi.R")
  source("plot.net.epi.R")
  
  my.g <- barabasi.game(100,directed=FALSE)
  my.g <- erdos.renyi.game(100,100,type=c("gnm"),directed=FALSE)
  
  ## Simulate network-based SIR process.
  beta <- 5
  gamma <- 1
  ntrials <- 100
  sim.res <- simul.net.epi(my.g,beta,gamma,ntrials)
  plot.net.epi(sim.res)
}
