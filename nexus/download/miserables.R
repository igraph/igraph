
source("sgb-coappearance.R")
tmp <- tempdir()
miserables <-
  getnet("http://mirror.switch.ch/ftp/mirror/tex/support/graphbase/jean.dat")

miserables$name <- "Les Miserables coappearance network"
save(miserables, file="miserables.Rdata")

