
source("sgb-coappearance.R")
copperfield <-
  getnet("http://mirror.switch.ch/ftp/mirror/tex/support/graphbase/david.dat")

copperfield$name <- "David Copperfield coappearance network"
save(copperfield, file="/tmp/copperfield.Rdata")
