
url <- "http://people.maths.ox.ac.uk/~porterm/data/facebook5.zip"
urlfull <- "http://people.maths.ox.ac.uk/~porterm/data/facebook100.zip"

tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "facebook5.zip")
destfull <- paste(sep="", tmp, "/", "facebook100.zip")
download.file(url, dest)
download.file(urlfull, destfull)

system(paste("cd", tmp, ";", "unzip", destfull))

library(R.matlab)

matfiles <- grep("\\.mat$", list.files(paste(sep="", tmp, "/facebook100")),
                 value=TRUE)

nullna <- function(x) {
  ifelse(x==0, NA, x)
}

dofile <- function(mfile) {
  u <- sub("[0-9]+\\.mat", "", mfile)
  M <- readMat(paste(sep="", tmp, "/facebook100/", mfile))
  g <- graph.adjacency(M[[1]], mode="undirected")
  V(g)$SFStatus   <- nullna(M[[2]][,1])
  V(g)$Gender     <- nullna(M[[2]][,2])
  V(g)$Major      <- nullna(M[[2]][,3])
  V(g)$Minor      <- nullna(M[[2]][,4])
  V(g)$Dorm       <- nullna(M[[2]][,5])
  V(g)$Year       <- nullna(M[[2]][,6])
  V(g)$HighSchool <- nullna(M[[2]][,7])
  g$name <- paste("Facebook network,", u)
  g$Author <- "Mason A. Porter"
  g$Citation <- "Amanda L. Traud, Eric D. Kelsic, Peter J. Mucha, and Mason A. Porter,Comparing Community Structure to Characteristics in Online Collegiate Social Networks, SIAM Review, in press (arXiv:0809.0690)."
  g$URL <- "http://www.insna.org/member/profiles/356.html"

  g
}

