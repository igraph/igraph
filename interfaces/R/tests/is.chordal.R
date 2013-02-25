
library(igraph)

## The examples from the Tarjan-Yannakakis paper
g1 <- graph.formula(A-B:C:I, B-A:C:D, C-A:B:E:H, D-B:E:F,
                    E-C:D:F:H, F-D:E:G, G-F:H, H-C:E:G:I,
                    I-A:H)
g1 <- simplify(g1)
maximum.cardinality.search(g1)
is.chordal(g1, fillin=TRUE)

g2 <- graph.formula(A-B:E, B-A:E:F:D, C-E:D:G, D-B:F:E:C:G,
                    E-A:B:C:D:F, F-B:D:E, G-C:D:H:I, H-G:I:J,
                    I-G:H:J, J-H:I)
g2 <- simplify(g2)
maximum.cardinality.search(g2)
is.chordal(g2, fillin=TRUE)
