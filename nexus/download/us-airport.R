
url <- "http://sites.google.com/site/cxnets/US_largest500_airportnetwork.txt"

tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "usairport.txt")

download.file(url, dest)

library(igraph)
usairport <- read.graph(dest, format="ncol")
usairport$name <- "US air transportation network"
usairport$Author <- "V. Colizza, R. Pastor-Satorras and A. Vespignani"
usairport$Citation <- "V. Colizza, R. Pastor-Satorras and A. Vespignani. Reaction-diffusion processes and metapopulation models in heterogeneous networks. Nature Physics 3, 276-282 (2007)."
usairport$URL <- "http://sites.google.com/site/cxnets/usairtransportationnetwork"

save(air, file="/tmp/usairport.Rdata")
