
library(igraph)

pattern <- graph.formula(1:2:3:4:5,
                         1 - 2:5, 2 - 1:5:3, 3 - 2:4, 4 - 3:5, 5 - 4:2:1)
target <- graph.formula(1:2:3:4:5:6:7:8:9,
                        1 - 2:5:7, 2 - 1:5:3, 3 - 2:4, 4 - 3:5:6:8:9,
                        5 - 1:2:4:6:7, 6 - 7:5:4:9, 7 - 1:5:6,
                        8 - 4:9, 9 - 6:4:8)
domains <- list(`1` = c(1,3,9), `2` = c(5,6,7,8), `3` = c(2,4,6,7,8,9),
                `4` = c(1,3,9), `5` = c(2,4,8,9))
graph.subisomorphic.lad(pattern, target, all.maps=TRUE)
graph.subisomorphic.lad(pattern, target, induced=TRUE, all.maps=TRUE)
graph.subisomorphic.lad(pattern, target, domains=domains, all.maps=TRUE)

