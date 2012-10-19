library(igraph)

g &lt;- barabasi.game(100000)
d &lt;- degree(g, mode="in")
dd &lt;- degree.distribution(g, mode="in", cumulative=TRUE)
alpha &lt;- power.law.fit(d, xmin=20)
plot(dd, log="xy", xlab="degree", ylab="cumulative frequency",
     col=1, main="Nonlinear preferential attachment")
lines(10:500, 10*(10:500)^(-coef(alpha)+1))

powers &lt;- c(0.9, 0.8, 0.7, 0.6)
for (p in seq(powers)) {
  g &lt;- barabasi.game(100000, power=powers[p])
  dd &lt;- degree.distribution(g, mode="in", cumulative=TRUE)
  points(dd, col=p+1, pch=p+1)
}

legend(1, 1e-5, c(1,powers), col=1:5, pch=1:5, ncol=1, yjust=0, lty=0)
