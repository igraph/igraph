

source("sgb-coappearance.R")
huckleberry <- 
  getnet("http://mirror.switch.ch/ftp/mirror/tex/support/graphbase/huck.dat")

huckleberry$name <- "Huckleberry Finn coappearance network"
save(huckleberry, file="/tmp/huckleberry.Rdata")

