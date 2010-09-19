
library(igraph)

g <- graph.formula( a -+ b -+ c -+ d -+ e )
all.st.cuts(g, source="a", target="e")

g2 <- graph.formula( s -+ a:b -+ t, a -+ 1:2:3 -+ b )
all.st.cuts(g2, source="s", target="t")

g3 <- graph.formula( s -+ a:b -+ t, a -+ 1:2:3:4:5 -+ b )
all.st.mincuts(g2, source="s", target="t")

