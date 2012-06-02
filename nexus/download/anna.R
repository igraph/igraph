
source("sgb-coappearance.R")

tmp <- tempdir()
annakarenina <-
  getnet("http://mirror.switch.ch/ftp/mirror/tex/support/graphbase/anna.dat")

annakarenina$name <- "Anna Karenina coappearance network"
save(annakarenina, file="/tmp/annakarenina.Rdata")
