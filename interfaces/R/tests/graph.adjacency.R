
# sink("graph.adjacency.Rout.save")

library(igraph)

M1 <- rbind(c(0,0,1,1),
            c(1,0,0,0),
            c(0,1,0,1),
            c(1,0,0,1))
g1 <- graph.adjacency(M1)
el1 <- get.edgelist(g1)
el1[order(el1[,1], el1[,2]),]

M2 <- rbind(c(0,1,1,1),
            c(1,0,0,0),
            c(1,0,0,1),
            c(1,0,1,0))
g2 <- graph.adjacency(M2, mode="undirected")
el2 <- get.edgelist(g2)
el2[order(el2[,1], el2[,2]),]

M3 <- rbind(c(0,1,1,2),
            c(1,0,0,0),
            c(1,0,0,0),
            c(1,0,1,0))
g3 <- graph.adjacency(M3, mode="min")
el3 <- get.edgelist(g3)
el3[order(el3[,1], el3[,2]),]

M4 <- rbind(c(0,1,1,2),
            c(1,0,0,0),
            c(1,0,0,0),
            c(1,0,1,0))
g4 <- graph.adjacency(M4, mode="max")
el4 <- get.edgelist(g4)
el4[order(el4[,1], el4[,2]),]

M5 <- rbind(c(0,1,1,2),
            c(1,0,0,0),
            c(1,0,0,0),
            c(1,0,1,0))
g5 <- graph.adjacency(M5, mode="upper")
el5 <- get.edgelist(g5)
el5[order(el5[,1], el5[,2]),]

M6 <- rbind(c(0,1,1,2),
            c(1,0,0,0),
            c(1,0,0,0),
            c(1,0,1,0))
g6 <- graph.adjacency(M6, mode="lower")
el6 <- get.edgelist(g6)
el6[order(el6[,1], el6[,2]),]

M7 <- rbind(c(0,1,1,2),
            c(1,0,0,0),
            c(1,0,0,0),
            c(1,0,1,0))
g7 <- graph.adjacency(M7, mode="plus")
el7 <- get.edgelist(g7)
el7[order(el7[,1], el7[,2]),]

M8 <- rbind(c(0,1,1,0.5),
            c(1,0,0,0),
            c(1,0,0,0),
            c(1,0,2,0))
g8 <- graph.adjacency(M8, mode="directed", weighted=TRUE)
el8 <- cbind(get.edgelist(g8), E(g8)$weight)
el8[order(el8[,1], el8[,2]),]

M9 <- rbind(c(0,1,1,3),
            c(1,0,0,0),
            c(1,0,0,2),
            c(3,0,2,0))
g9 <- graph.adjacency(M9, mode="undirected", weighted=TRUE)
el9 <- cbind(get.edgelist(g9), E(g9)$weight)
el9[order(el9[,1], el9[,2]),]

M10 <- rbind(c(0,1,1,0.5),
            c(1,0,0,0),
            c(1,0,0,0),
            c(1,0,2,0))
g10 <- graph.adjacency(M10, mode="max", weighted=TRUE)
el10 <- cbind(get.edgelist(g10), E(g10)$weight)
el10[order(el10[,1], el10[,2]),]

M11 <- rbind(c(0,1,1,0.5),
            c(1,0,0,0),
            c(1,0,0,0),
            c(1,0,2,0))
g11 <- graph.adjacency(M11, mode="min", weighted=TRUE)
el11 <- cbind(get.edgelist(g11), E(g11)$weight)
el11[order(el11[,1], el11[,2]),]

M12 <- rbind(c(0,1,1,0.5),
            c(1,0,0,0),
            c(1,0,0,0),
            c(1,0,2,0))
g12 <- graph.adjacency(M12, mode="lower", weighted=TRUE)
el12 <- cbind(get.edgelist(g12), E(g12)$weight)
el12[order(el12[,1], el12[,2]),]

M13 <- rbind(c(0,1,1,0.5),
            c(1,0,0,0),
            c(1,0,0,0),
            c(1,0,2,0))
g13 <- graph.adjacency(M13, mode="upper", weighted=TRUE)
el13 <- cbind(get.edgelist(g13), E(g13)$weight)
el13[order(el13[,1], el13[,2]),]

M14 <- rbind(c(0,1,1,0.5),
            c(1,0,0,0),
            c(1,0,0,0),
            c(1,0,2,0))
g14 <- graph.adjacency(M14, mode="plus", weighted=TRUE)
el14 <- cbind(get.edgelist(g14), E(g14)$weight)
el14[order(el14[,1], el14[,2]),]

# sink(NULL)
