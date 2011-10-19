
#   IGraph R package
#   Copyright (C) 2005  Gabor Csardi <csardi@rmki.kfki.hu>
#   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
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

######################################################
####  Cohesive Blocking Algorithm and Utilities   ####
######################################################
####           Written by Peter McMahan           ####
######################################################

######################################################
## Main cohesive blocking function. Feed it a graph ##
## of class igraph and it returns an object of      ##
## class bgraph. Definitely requires digest, and    ##
## will use RSQLite for large graphs if it is avai- ##
## lable.                                           ##
######################################################

cohesive.blocks <- function(graph, db=NULL,
                            useDB=(vcount(graph)>400 && require(RSQLite)),
                            cutsetHeuristic=TRUE,
                            verbose=igraph.par("verbose")) {

    if(useDB && !require(RSQLite)) stop("package `RSQLite` required")
    if(!require(digest)) stop("package `digest` required")
    if(!is.igraph(graph)) stop("`graph' must be an igraph object")
    verbose <- as.logical(verbose)
    
    Gin <- graph
    graph <- simplify(as.undirected2(graph))
    V(graph)$cbid <- as.numeric(V(graph))
    infCohesion <- vcount(graph)+1 ## for fully connected subcomponents
    g <- graph
    v <- sort(as.numeric(V(g)))
    
    ## set up database connection and schema
    if(is.null(db) && useDB){
        ## initialize the database and define schema
        dbfile <- tempfile("cohesive.blocks")
        con <- dbConnect(dbDriver("SQLite"), dbname=dbfile)
        ## branchMembership
        dbSendQuery(con, "create table branchMembership (vertexId integer, branchId integer)")
        dbSendQuery(con, "create index index_membership on branchMembership (vertexId, branchId)")
        dbSendQuery(con, "create index index_vertex on branchMembership (vertexId)")
        dbSendQuery(con, "create index index_branch on branchMembership (branchId)")
        ## branches (unexamined branches given cohesion of (-1))
        dbSendQuery(con, "create table branches (branchId integer primary key, membershipHash text, branchCohesion integer default -1, maxAncestorCohesion integer, isCohesiveBlock integer default 0)")
        dbSendQuery(con, "create index index_hash on branches (membershipHash)")
        ## subBranches
        dbSendQuery(con, "create table subBranches (parentId integer, childId integer)")
        dbSendQuery(con, "create index index_parentId on subBranches (parentId)")
        
        ##
        ## Add the root branch Ñ the whole graph
        ##
        ## create the hash
        v <- sort(as.numeric(V(g)))
        thisHash <- digest(paste(paste(v, collapse=" "), 0))
        
        dbBeginTransaction(con)
        dbSendQuery(con, paste("insert into branches (branchId,membershipHash,maxAncestorCohesion,isCohesiveBlock) values(", 0, ",'", thisHash, "',", 0, ",", 1, "); ", sep="")) 
        for(thisV in v){
            dbSendQuery(con, paste("insert into branchMembership values(", thisV, ",", 0, "); ", sep=""))
        }
        dbCommit(con)
        maxId <- 0
        
    }else if(useDB){
        ## connect to existing db
        dbfile <- db
        con <- dbConnect(dbDriver("SQLite"), dbname=dbfile)
        maxId <- dbGetQuery(con, "select max(branchId) from branches")
    } else {
        branchMembership <- data.frame(vertexId=sort(as.numeric(V(graph))), branchId=0)
        thisHash <- digest(paste(paste(v, collapse=" "), 0))
        branches <- data.frame(branchId=0, membershipHash=thisHash, branchCohesion=-1, maxAncestorCohesion=0, isCohesiveBlock=1)
        subBranches <- data.frame(parentId=NA, childId=NA)
    }
    
    ## now that setup's done, start looping
    notDone <- TRUE
    while(notDone){
        ## which branch are we working on this time
        ## DATA ACCESS:
        if(useDB){
            branchId <- as.numeric(dbGetQuery(con, "select branchId from branches where branchCohesion=-1 limit 1"))
        } else {
            branchId <- branches$branchId[branches$branchCohesion==-1][1]
        }
        if (verbose) cat("Branch ", branchId, ":  ", sep="")
        
        ## get data and make g
        ## DATA ACCESS:
        if(useDB){
            v <- sort(as.numeric(unlist(dbGetQuery(con, paste("select vertexId from branchMembership where branchId=", branchId, sep=""))[[1]])))
            maxId <- as.numeric(dbGetQuery(con, "select max(branchId) from branches"))
        } else {
            v <- sort(branchMembership$vertexId[branchMembership$branchId==branchId])
            maxId <- max(branches$branchId)
        }

        cbid <- NULL # To fix a check NOTE
        g <- subgraph(graph, V(graph)[cbid %in% v])

        ## check if trivial or fully connected, else treat normally
        if(vcount(g) < 3){ #trivial
            if (verbose) cat("Branch single or binary -- no cohesive subgroups.\n")
            k <- 0
            if(useDB){
                dbSendQuery(con, paste("update branches set branchCohesion=", k, " where branchId=", branchId, ";", sep=""))
                notDone <- dbGetQuery(con, "select count(*) from branches where branchCohesion=-1")>0
                #dbSendQuery(con, paste("update branches set branchCohesion=0 where branchId=", branchId, ";", sep=""))
            } else {
                branches$branchCohesion[branches$branchId==branchId] <- k
                notDone <- any(branches$branchCohesion==-1)
            }
        }else if(all(degree(g)>=vcount(g))){ #fully  connected
            if (verbose) cat("Branch fully connected -- setting cohesion to ", infCohesion, ".\n", sep="")
            k <- infCohesion
            if(useDB){
                dbSendQuery(con, paste("update branches set branchCohesion=", k, " where branchId=", branchId, ";", sep=""))
                notDone <- dbGetQuery(con, "select count(*) from branches where branchCohesion=-1")>0
                #dbSendQuery(con, paste("update branches set branchCohesion=", infCohesion, ",  isCohesiveBlock=1 where branchId=", branchId, ";", sep=""))
            } else {
                branches$branchCohesion[branches$branchId==branchId] <- k
                notDone <- any(branches$branchCohesion==-1)
            }
        } else { #keep looking
            ## find cohesion of g
            if(!is.connected(g)){
                k <- 0
            }else if(min(degree(g))==1){
                k <- 1
            } else {
                k <- graph.cohesion(g)
            }
            if (verbose) cat(k, "-cohesive; ", sep="")
            
            ## trim it down:
            if(useDB){
                mac <- max(k, unlist(dbGetQuery(con, paste("select maxAncestorCohesion from branches where branchId =", branchId, ";"))))
            } else {
                mac <- max(k, unlist(branches$maxAncestorCohesion[branches$branchId==branchId]))
            }
            if (verbose) cat(length(v))
            while(any(degree(g)<=mac)){ ## trimming goes on here
                g <- subgraph(g, V(g)[degree(g)>mac])
            }
            if (verbose) cat("(", vcount(g), ") vertices; ", sep="")
            
            ## find k-components
            ## For nastier graphs:
            ## check if cohesion after trimming is greater than k. If so,  kcomp <- list(as.numeric(V(g)))
            if(vcount(g)>100 && k>1){ ## only try speedup for relatively complicated kComponents() calls
                if(!is.connected(g)){
                    newk <- 0
                }else if(min(degree(g))==1){
                    newk <- 1
                } else {
                    newk <- graph.cohesion(g)
                }
                if(newk>k){
                    kcomp <- list(as.numeric(V(g)))
                } else {
                    kcomp <- kComponents(g, k, verbose=verbose)
                }
            } else { ## otherwise just use kComponents()
                kcomp <- kComponents(g, k, cutsetHeuristic, verbose=verbose)
            }
            
            kcomp <- lapply(kcomp, function(thisV){V(g)[thisV]$cbid})
            if (verbose) cat(length(kcomp), " sub-branches,  ", sep="")
            ## which have been searched?
            theseHashes <- lapply(kcomp, function(x){digest(paste(paste(sort(unlist(x)), collapse=" "), mac))})
            
            ## DATA ACCESS:
            if(useDB){
                searched <- dbGetQuery(con, paste("select branchId, membershipHash from branches where membershipHash in ('", paste(theseHashes, collapse="','"), "');", sep=""))
            } else {
                searched <- branches[c("branchId", "membershipHash")][branches$membershipHash %in% theseHashes, ]
            }
            searchedCount <- nrow(searched)
            
            ## assign appropriate branchIds
            branchIds <- as.numeric(rep(NA, length(kcomp)))
            isNew <- rep(TRUE, length(kcomp))
            if(nrow(searched)>0){
                branchIds[ind <- match(searched$membershipHash, theseHashes)] <- searched$branchId
                branchIds[-ind] <- (maxId+1):(maxId+length(branchIds)-length(ind))
                isNew[ind] <- FALSE
            }else if(length(branchIds)>0){
                branchIds <- (maxId+1):(maxId+length(branchIds))
            }
            
            ## DATA ACCESS:
            ## add all branches to `subBranches` and only new branches to `branches` and `branchMembership`
            if(useDB){
                dbBeginTransaction(con)
                if(length(kcomp)>0){for(i in 1:length(kcomp)){
                    dbSendQuery(con, paste("insert into subBranches values(", branchId, ",", branchIds[i], ");", sep=""))
                    if(isNew[i]){ ## if new insert into branches branchMembership
                        dbSendQuery(con, paste("insert into branches (branchId,membershipHash,maxAncestorCohesion) values(", branchIds[i], ",'", theseHashes[i], "',", mac, ");", sep=""))
                        for(thisV in sort(unlist(kcomp[[i]]))){
                            dbSendQuery(con, paste("insert into branchMembership values(", thisV, ",", branchIds[i], "); ", sep=""))
                        }
                    }
                }}
                dbSendQuery(con, paste("update branches set branchCohesion=", k, " where branchId=", branchId, ";", sep=""))
                dbCommit(con)
                notDone <- dbGetQuery(con, "select count(*) from branches where branchCohesion=-1")>0
            } else {
                if(length(kcomp)>0){
                    subBranches <- rbind(subBranches, data.frame(parentId=rep(branchId, length(branchIds)), childId=branchIds))
                    branches <- rbind(branches, data.frame(branchId=branchIds, membershipHash=unlist(theseHashes), branchCohesion=-1, maxAncestorCohesion=mac, isCohesiveBlock=0)[isNew, ])
                    branchMembership <- rbind(branchMembership, data.frame(vertexId=unlist(kcomp[isNew]), branchId=rep(branchIds[isNew], lapply(kcomp[isNew], length))))
                }
                
                branches$branchCohesion[branches$branchId==branchId] <- k
                notDone <- any(branches$branchCohesion==-1)
            }
            
            if (verbose) cat("(", searchedCount, " complete);\n", sep="")
        }
    }
    
    ## mark the cohesive blocks
    if(useDB){
        dbSendQuery(con, "update branches set isCohesiveBlock=1 where maxAncestorCohesion < branchCohesion")
    } else {
        branches$isCohesiveBlock[branches$maxAncestorCohesion < branches$branchCohesion] <- 1
    }
    
    ## and extract the data we need
    if(useDB){
        cBlocks <- dbGetQuery(con, "select * from branches where isCohesiveBlock>0")
        cbMembership <- dbGetQuery(con, paste("select * from branchMembership where branchId in (", paste(cBlocks$branchId, collapse=", "), ");", sep=""))
        subBranches <- dbReadTable(con, "subBranches")
        lapply(dbListResults(con), dbClearResult) ## just-in-case cleanup
        dbDisconnect(con)
    } else {
        cBlocks <- branches[branches$isCohesiveBlock>0, ]
        cbMembership <- branchMembership[branchMembership$branchId %in% cBlocks$branchId, ]
        subBranches <- subBranches[!(is.na(subBranches$parentId) | is.na(subBranches$childId)), ]
    }
    
    ## get block members and cohesion
    blocks <- list()
    block.cohesion <- numeric()
    for(i in cBlocks$branchId){
        blocks <- c(blocks, list(cbMembership$vertexId[cbMembership$branchId==i]))
        block.cohesion <- c(block.cohesion, cBlocks$branchCohesion[cBlocks$branchId==i])
    }
    
    ## build the block hierarchy graph. first the vertices and attributes:
    tree <- graph.empty(nrow(cBlocks))
    V(tree)$branchId <- cBlocks$branchId
    V(tree)$nodes <- blocks
    V(tree)$chnodes <- lapply(blocks, paste, collapse=",")
    V(tree)$cohesion <- block.cohesion
    ## and now the edges:
    for(v in V(tree)){
        p <- find.cohesive.parents(V(tree)[v]$branchId, cBlocks$branchId, subBranches)
        e <- rep(v, times=2*length(p))
        e[(1:length(p))*2-1] <- V(tree)[match(p, V(tree)$branchId)-1]
        tree <- add.edges(tree, e)
    }
    
    ## put it all together and return
    Gin <- set.graph.attribute(Gin, "tree", tree)
    Gin <- set.graph.attribute(Gin, "blocks", blocks)
    Gin <- set.graph.attribute(Gin, "block.cohesion", block.cohesion)
    if (useDB) {
      Gin <- set.graph.attribute(Gin, "data", dbfile)
    } else {
      data <- list(branches=branches, branchMembership=branchMembership,
                   subBranches=subBranches)
      Gin <- set.graph.attribute(Gin, "data", data)
    }

    class(Gin) <- c("bgraph", "igraph")
    return(Gin)
}

structurally.cohesive.blocks <- cohesive.blocks

######################################################
## A small helper function to identify the cohesive ##
## parents of a given block inside cohesive.blocks. ##
######################################################

find.cohesive.parents <- function(id, blockIds, subBranches){
    res <- numeric()
    Q <- subBranches$parentId[subBranches$childId==id]
    while(length(Q)>0){
        thisId <- Q[1]
        Q <- Q[-1]
        if(is.na(thisId)) stop("NA in subBranches. something's wrong")
        if(thisId %in% blockIds){
            res <- c(res, thisId)
        } else {
            Q <- c(Q, subBranches$parentId[subBranches$childId==thisId])
        }
    }
    return(unique(res))
}

find.all.min.cutsets <- function(g, k=NULL){
    ## implements (mostly) algorithm from:
    ## Kanevsky,  Arkady 1990 "On the Number of Minimum Size
    ## Separating Vertex Sets in a Graph and How to Find
    ## All of Them"
    
    #- preliminary setup and checks
    if (!is.igraph(g)) {
        stop("Not an igraph object")
    }

    #- find connectivity if necessary:
    if(is.null(k)){
        if(!is.connected(g)){
            k <- 0
        }else if(min(degree(g))<=1){
            k <- 1
        } else {
            k <- vertex.connectivity(g)
        }
    }
    
    if (k==0){
        #warning("Graph not connected. Minimal cutset is empty.")
        return(numeric())
    }
    res <- list()
    v <- as.numeric(V(g))
    
    
    if(k==1){ # for 1-connected graphs
        return(articulation.points(g))
    } else { # else generic cutset finder
        #- take a k-set X
        x <- v[order(degree(g), decreasing=TRUE)][1:k]
    
        #- check x for cutsetness
        if(!is.connected(subgraph(g, setdiff(v, x)))){
            res <- c(res, list(x))
        }
    
        #- for each nonadjacent pair (x, v) in (X, V) with connectivity k
        for(i in x){for(j in v){if(length(get.shortest.paths(g, i, j)[[1]])>2){
            #cat(i, j, "\n")
        
        #   - grab the maxflow f between these
            g.r <- simplify(g)
            phi <- numeric()
            done <- FALSE
            while(!done){
                p <- get.shortest.paths(g.r, i, j, mode="out")[[1]]
                if(length(p)<1) done <- TRUE
                phi <- c(phi, p)
                e <- E(g.r, path=p)
                g.r <- delete.edges(g.r, e)
                if(vertex.disjoint.paths(g.r, i, j)==0) done <- TRUE
            }
            phi <- phi[phi!=i]
            phi <- phi[phi!=j]
        #   - try every k-set of vertices in f
        #   - add the cutsets to result set
            if(length(phi)>0){
                ksets <- combn(phi, k)
                for(thisSet in 1:ncol(ksets)){
                    if(!is.connected(subgraph(g, setdiff(v, ksets[, thisSet])))){
                        if(!(list(ksets[, thisSet]) %in% res)){
                            res <- c(res, list(ksets[, thisSet]))
                        }
                    }
                }
            }
            
        #   - add edge (x, v) to E
            add.edges(g, c(i, j))
        }}}
        if(length(res)>0){
            for(i in 1:length(res)){
                res[[i]] <- sort(res[[i]])
            }
        }
        return(unique(res))
    }
}

kComponents <- function(g, k=NULL, useHeuristic=TRUE, verbose=igraph.par("verbose")){
    if(vcount(g)<1) return(list())
    V(g)$csid <- as.numeric(V(g))
    cs <- if(useHeuristic){find.all.min.cutsets(g, k)} else {kCutsets(g, k)}
    theseBlocks <- list()
    if(length(cs)==0){## not connected
        cls <- clusters(g)
        for(cl in unique(cls$membership)){
            theseBlocks <- c(theseBlocks, list(V(subgraph(g, which(cls$membership==cl)-1))$csid))
        }
        return(theseBlocks)
    }
    for(thisCS in cs){
        if(is.list(thisCS)){thisCS <- unlist(thisCS)}
        gprime <- subgraph(g, V(g)[!(V(g) %in% thisCS)])
        gclusts <- clusters(gprime)
        for(i in ((1:length(gclusts$csize))-1)){
            theseIDs <- V(gprime)[V(gprime) %in% (which(gclusts$membership==i)-1)]$csid
            thisBlockBig <- list(sort(c(theseIDs, thisCS)))
            if(!(thisBlockBig %in% theseBlocks) && length(theseIDs)>0){
                theseBlocks <- c(theseBlocks, thisBlockBig)
            }
        }
    }
    return(theseBlocks)
}

is.cutset <- function(v, g){ ## does removal of `v` disconnect `g`?
    return(!is.connected(subgraph(g, setdiff(V(g), v))))
}

etReduction <- function(g){
    G <- graph.empty(vcount(g)*2, directed=TRUE)
    
    el1 <- el2 <- get.edgelist(g)
    el1[, 1] <- el1[, 1] + vcount(g)
    el1 <- as.numeric(t(el1)) ## u' -> v'' edges (external)
    el2[, 2] <- el2[, 2] + vcount(g)
    el2 <- el2[, 2:1]
    el2 <- as.numeric(t(el2)) ## v' -> u'' edges (external)
    el3 <- as.numeric(rbind(0:(vcount(g)-1), vcount(g):(2*vcount(g)-1))) ## v' -> v'' edges (internal)
    
    G <- add.edges(G, c(el1, el2), capacity=Inf) ## (external)
    G <- add.edges(G, el3, capacity=1) ## (internal)
    
    return(simplify(G))
}

# Not actually needed

Phi <- function(g, i, j){
    g.r <- simplify(g)
    phi <- numeric()
    done <- FALSE
    while(!done){
        p <- get.shortest.paths(g.r, i, j, mode="out")[[1]]
        if(length(p)<1) done <- TRUE
        phi <- c(phi, p)
        e <- E(g.r, path=p)
        g.r <- delete.edges(g.r, e)
        if(vertex.disjoint.paths(g.r, i, j)==0) done <- TRUE
    }
    phi <- phi[phi!=i]
    phi <- phi[phi!=j]
    return(phi)
}    

## a small fix for as.undirected so it won't strip vertex attributes:
as.undirected2 <- function(g){
    ## store the attributes:
    tempAtt <- list()
    for(a in list.vertex.attributes(g)){
        tempAtt <- c(tempAtt, list(get.vertex.attribute(g, a)))
        names(tempAtt)[length(tempAtt)] <- a
    }
    
    ## convert the graph:
    res <- as.undirected(g)
    
    ## restore the attributes:
    for(a in names(tempAtt)){
        set.vertex.attribute(res, a, value=tempAtt[[a]])
    }
    return(res)
}

######################################################
## The following are specific utility functions for ##
## dealing with the bgraphs output by the cohesive  ##
## blocking algorithm. They include functions to    ##
## plot, recalculate and save into different        ##
## formats. (pretty messy, sorry)                   ##
## The most important are:                          ##
## - plot.bgraph()                                  ##
## - maxcohesion()                                 ##
## - write.pajek.bgraph()                           ##
######################################################

is.bgraph <- function(graph) {
  "bgraph" %in% class(graph)
}

recalc.tree.bgraph <- function(graph){

    if (!is.bgraph(graph)) {
      stop("Not a bgraph object")
    }
  
    ## make the block-hierarchy tree:
    blocks <- get.graph.attribute(graph,  "blocks")
    block.cohesion <- get.graph.attribute(graph,  "block.cohesion")
    tree <- graph.empty(n=length(blocks), directed=TRUE)
    block.parents <- rep(0, length(block.cohesion))
    for(i in 1:length(blocks)){
        do.contain <- rep(FALSE, length(blocks))
        for(j in 1:length(blocks)){
            do.contain[j] <- all(blocks[[i]] %in% blocks[[j]])
        }
        this.parent <- 0
        ml <- Inf
        for(j in which(do.contain)){
            if(length(blocks[[j]])<ml && !(length(blocks[[i]]) == length(blocks[[j]]))){
                ml <- length(blocks[[j]])
                this.parent <- j
            }
        }
        if(this.parent>0){
            tree <- add.edges(tree, c(this.parent-1, i-1))
            block.parents[i] <- this.parent-1
        }
        V(tree)$nodes[i] <- list(blocks[[i]])
        V(tree)$chnodes[i] <- paste(blocks[[i]], collapse=", ")
        V(tree)$cohesion[i] <- block.cohesion[i]
    }
    graph <- set.graph.attribute(graph, "tree", tree)
    return(graph)
}

print.bgraph <- function(x, ...){

    if (!is.bgraph(x)) {
      stop("Not a bgraph object")
    }
    cat("Graph:\n\n")
    print.igraph(x, ...)
    cat("\nBlock Hierarchy Tree:\n\n")
    print.igraph(get.graph.attribute(x, "tree"), ...)
    invisible(x)
}

maxcohesion <- function(graph){

    if (!is.bgraph(graph)) {
      stop("Not a bgraph object")
    }

    bc <- get.graph.attribute(graph, "block.cohesion")
    bco <- order(bc)
    bc <- bc[bco]
    b <- get.graph.attribute(graph, "blocks")[bco]
    mc <- numeric(vcount(graph))
    for(i in 1:length(bc)){
        mc[b[[i]]+1] <- bc[[i]]
    }
    return(mc)
}

layout.svd.bgraph <- function(graph, d=shortest.paths(graph), s=.6, ...) {

  if (!is.bgraph(graph)) {
    stop("Not a bgraph object")
  }

  for(b in graph$blocks){
    d[b+1, b+1] <- d[b+1, b+1]*s
  }

  layout.svd.igraph(graph, d=d)
}

layout.mds.bgraph <- function(graph, d=shortest.paths(graph), s=.6, ...) {

  if (!is.bgraph(graph)) {
    stop("Not a bgraph object")
  }

  for(b in graph$blocks){
    d[b+1, b+1] <- d[b+1, b+1]*s
  }

  layout.mds.igraph(graph, d=d)
}

plot.bgraph <- function(x, mc=NULL, vertex.size=3, colpal=NULL, emph=NULL, ...){

    if (!is.bgraph(x)) {
      stop("Not a bgraph object")
    }
  
    layout(matrix(c(1, 1, 2), nrow=1))

    g <- x
    if(is.null(colpal)){
        block.cohesion <- get.graph.attribute(g, "block.cohesion")
        palette(rev(rainbow(max(block.cohesion)+1, end=5/6)))
    } else {
        palette(colpal)   
    }
    
    if(is.null(mc)){
      if("mc" %in% list.vertex.attributes(g)){
        mc <- V(g)$mc
      } else {
        mc <- maxcohesion(g)
      }
    }
    
    bord <- rep("black", times=vcount(g))
    vertex.size <- rep(vertex.size, times=vcount(g))
##     vlabels <- rep("", times=vcount(g))
    if(!is.null(emph)){
        for(i in emph){
            emphme <- get.graph.attribute(g, "blocks")[[i]]+1
            bord[emphme] <- "white"
            vertex.size[emphme] <- vertex.size[1]*2.5
##             vlabels[emphme] <- "*"
        }
    }
    
    plot.igraph(g, vertex.size=vertex.size, vertex.color=mc+1,
                vertex.frame.color=bord, vertex.label.color="white", ...)

    tree <- get.graph.attribute(g, "tree")
    r <- which.min(V(tree)$cohesion)-1
    l <- layout.reingold.tilford(tree, root=r)
    l[, 2] <- -l[, 2]
    V(tree)$label <- V(tree)$cohesion
    plot.igraph(tree, layout=l, vertex.color=V(tree)$cohesion+1)
    
    layout(1)
    invisible(NULL)
}

plotkCore.bgraph <- function(graph, vertex.size=3, colpal=NULL, emph=NULL,
                              remove.multiple=TRUE, remove.loops=TRUE, ...){

    if (!is.bgraph(graph)) {
      stop("Not a bgraph object")
    }

    g <- graph
    kc <- graph.coreness(simplify(g, remove.multiple, remove.loops))
    
    if(is.null(colpal)){
        palette(rev(rainbow(max(kc)+1, end=5/6)))
    } else {
        palette(colpal)   
    }
    
    bord <- rep("black", times=vcount(g))
    vertex.size <- rep(vertex.size, times=vcount(g))
##     vlabels <- rep("", times=vcount(g))
    if(!is.null(emph)){
        for(i in emph){
            emphme <- get.graph.attribute(g, "blocks")[[i]]+1
            bord[emphme] <- "white"
            vertex.size[emphme] <- vertex.size[1]*2.5
##             vlabels[emphme] <- "*"
        }
    }
    
    plot.igraph(g, vertex.size=vertex.size, vertex.color=kc+1,
                vertex.frame.color=bord, vertex.label.color="white", ...)
    
    kcores <- sort(unique(kc))
    legend("topright", legend=kcores, pch=".", pt.cex=10, col=kcores+1)
}

## writes out a series of Pajek files for supplied bgraph object
## filename should not include an extension
write.pajek.bgraph <- function(graph, filename, hierarchy=FALSE){
  if (!is.bgraph(graph)) {
    stop("Not a bgraph object")
  }
  
  g <- graph
  ## .net file:
  write.graph(g, file=paste(filename, ".net", sep=""), format="pajek")
  
  ## .clu with max cohesions
  mc <- if("mc" %in% list.vertex.attributes(g)){
    V(g)$mc
  } else {
    maxcohesion(g)
  }
  cat(paste("*Vertices", vcount(g)), mc, "\r", sep="\r\n",
      file=paste(filename, ".clu", sep=""))
  
  ## .clu for each block giving binary membership (toggled with blockfiles argument)
  if(hierarchy){
    for(i in 1:length(g$blocks)){
      thisblock <- rep(0, vcount(g))
      thisblock[g$blocks[[i]]+1] <- 1
      cat(paste("*Vertices", vcount(g)), thisblock, "\r", sep="\r\n",
          file=paste(filename, "_block", i, "(", g
            $block.cohesion[i], ").clu", sep=""))
    }
    V(g$tree)$id <- as.numeric(V(g$tree))+1
    write.graph(g$tree, file=paste(filename, "_blocktree.net",
                          sep=""), format="pajek")
  }
}

maxcohesion <- function(graph){

  if (!is.bgraph(graph)) {
    stop("Not a bgraph object")
  }
  
  bc <- get.graph.attribute(graph, "block.cohesion")
  bco <- order(bc)
  bc <- bc[bco]
  b <- get.graph.attribute(graph, "blocks")[bco]
  mc <- numeric(vcount(graph))
  for(i in 1:length(bc)){
    mc[b[[i]]+1] <- bc[[i]]
  }
  return(mc)
}

## same as write.pajek.bgraph() but uses k-coreness rather than cohesion for clusters
write.pajek.kCore.bgraph <- function(graph, filename, remove.multiple=TRUE, remove.loops=TRUE){ # filename should not include an extension
    if (!is.bgraph(graph)) {
      stop("Not a bgraph object")
    }

    g <- graph
    ## .net file:
    write.graph(g, file=paste(filename, ".net", sep=""), format="pajek")
    
    ## .clu with kcoreness
    kc <- graph.coreness(simplify(g, remove.multiple, remove.loops))

    cat(paste("*Vertices", vcount(g)), kc, "\r", sep="\r\n", file=paste(filename, ".clu", sep=""))
}

######################################################
## kCutsets uses Kanevsy's algorithm to find all    ##
## cutsets of size k in the given graph.            ##
######################################################

kCutsets <- function(g, k=NULL){
    if(!is.igraph(g)){
        stop("g must be an igraph object")
    }
    if(is.null(k)){
        if(!is.connected(g)){
            k <- 0
        }else if(min(degree(g))<=1){
            k <- 1
        } else {
            k <- vertex.connectivity(g)
        }
    }
    
    if(k==0 || vcount(g) <=2) return(list(numeric()))
    
    if(k==1){ # for 1-connected graphs
        return(articulation.points(g))
    }
    
    #############
    ## begin Kanevsky's (parallel) algorithm:
    #############
    v <- as.numeric(V(g))
    
    K <- v[order(degree(g), decreasing=TRUE)][1:k] #2a
    if(is.cutset(K, g)) res <- c(res, list(K))
    
    P <- expand.grid(K, v) #2b
    
    G <- list(g) #3
    for(i in 2:nrow(P)){
        G[[i]] <- add.edges(g, as.numeric(t(P[1:(i-1), ])))
    }
    
    Gbar <- lapply(G, etReduction) #4 (sequential)
    mflow <- unlist(lapply(1:nrow(P), function(i){graph.maxflow(Gbar[[i]], P[i, 1]+vcount(g), P[i, 2])}))
    
    Gbar <- Gbar[mflow==k] #5 (sequential)
    P <- P[mflow==k, ]
        
    if(nrow(P)<1) return(list())
    GbarRes <- lapply(1:nrow(P), function(i){graph.residual(Gbar[[i]], P[i, 1]+vcount(g), P[i, 2])}) #6 (sequential)
    comp <- lapply(GbarRes, clusters, mode="strong")
    L <- lapply(1:length(comp), function(i){graph.shrink(GbarRes[[i]], comp=comp[[i]])})
    
    L <- lapply(L, graph.antichains) #7 (sequential)
    
    ## convert antichains in L to vertex sets in g
    res <- lapply(1:length(L), 
                  function(i){
                    gc()
                    l <- L[[i]]
                    cmp <- comp[[i]]
                    thisres <- list()
                    for(cn in l){
                      thisres <- unique(c(thisres, list(unique(sort((which(cmp$membership %in% cn)-1) %% vcount(g))))))
                    }
                    return(thisres)
                  }
                  )
    res <- unlist(res, recursive=FALSE)
        
    ## and see which are cutsets
    res <- unique(res[unlist(lapply(res, function(x){length(x)==k && is.cutset(x, g)}))])
    
    return(res)
}

graph.residual <- function(g, i, j){ ## NOT GENERIC: deletes edges once they get to nonpositive capacity
    gc()
    if(!("capacity" %in% list.edge.attributes(g))) E(g)$capacity <- 1
    g.r <- simplify(g)
    done <- FALSE
    while(!done){
        p <- get.shortest.paths(g.r, i, j, mode="out")[[1]]
        if(length(p)<1) done <- TRUE
        E(g.r, path=p)$capacity <- E(g.r, path=p)$capacity - 1
        capacity <- NULL # eliminate a check NOTE
        g.r <- delete.edges(g.r, E(g.r)[capacity<=0])
        if(vertex.disjoint.paths(g.r, i, j)==0) done <- TRUE
    }
    return(g.r)
}

graph.shrink <- function(g, comp=clusters(g, mode="strong"), remove.multiple=TRUE, remove.loops=TRUE){
    gc()
    res <- graph.empty(length(comp$csize), directed=TRUE)
    el <- get.edgelist(g)
    el.comp <- rbind(comp$membership[el[, 1]+1], comp$membership[el[, 2]+1])
    res <- add.edges(res, el.comp)
    res <- simplify(res, remove.multiple=remove.multiple, remove.loops=remove.loops)
    return(res)
}

graph.antichains <- function(g){ ## g must be a directed acyclic graph Ñ not checked!
    ## add edges from every vertex to each of its descendants
    gc()
    for(i in as.numeric(V(g))){
        el <- unlist(neighborhood(g, vcount(g), i, "out"))
        el <- el[el != i]
        if(length(el)>0){
            el <- as.numeric(rbind(i, el))
            g <- add.edges(g, el)
        }
    }
    
    ## make into an adjacency matrix
    adj <- get.adjacency(simplify(as.undirected(g)))
    
    ## iterate through combinations of antichains to look for new antichains:
    res <- as.list(as.numeric(V(g)))
    acount <- 0
    while(length(res) != acount){
        acount <- length(res)
        newchains <- combn(1:acount, 2, 
                           function(x){
                             thisChain <- unique(unlist(res[x]))
                             if(all(adj[thisChain+1, thisChain+1]==0)) return(sort(thisChain))
                             return(NA)
                           }, simplify=FALSE)
        newchains <- newchains[!is.na(newchains)]
        res <- unique(c(res, newchains))
    }
    
    return(res)
}
