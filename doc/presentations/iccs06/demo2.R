library(igraph)

exps <- seq(0.5, 1.5, length=16)

par(mar=c(2,2,2,2))
layout( matrix(1:16, nr=4, byrow=TRUE))
layout.show(16)

maxdeg <- numeric()
for (ex in exps) \{
  g <- \emph{barabasi.game}(100000, power=ex)
  maxdeg <- c(maxdeg, max(\emph{degree}(g, mode="in")))  
  plot(\emph{degree.distribution}(g, mode="in"), log="xy", xlab=NA, ylab=NA)
\}

x11()
plot(exps, maxdeg, type="b")

