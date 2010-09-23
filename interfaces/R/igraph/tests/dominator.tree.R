
library(igraph)

g <- graph.formula(R-+A:B:C, A-+D, B-+A:D:E, C-+F:G, D-+L,
                   E-+H, F-+I, G-+I:J, H-+E:K, I-+K, J-+I,
                   K-+I:R, L-+H)
dtree <- dominator.tree(g, root="R")

dtree$dom <- V(g)$name[ dtree$dom ]
dtree$leftout <- V(g)$name[ dtree$leftout ]
dtree
