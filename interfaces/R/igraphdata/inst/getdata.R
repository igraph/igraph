
#   igraphdata R package
#   Copyright (C) 2010  Gabor Csardi <csardi.gabor@gmail.com>
#   Rue de l'Industrie 5, Lausanne 1005, Switzerland
#   
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
#   02110-1301 USA
#
###################################################################

library(igraph)

tmp <- tempdir()

#####################################################################
## Foodwebs

u1 <- "http://vlado.fmf.uni-lj.si/pub/networks/data/bio/FoodWeb/Webs_paj.zip"
u2 <- "http://vlado.fmf.uni-lj.si/pub/networks/data/bio/FoodWeb/ATLSS_paj.zip"
foodzip <- paste(tmp, sep="/", c("f1.zip", "f2.zip"))
download.file(url=u1, destfile=foodzip[1])
download.file(url=u2, destfile=foodzip[2])
unlink(paste(tmp, sep="/", "paj"), recursive=TRUE)
system(paste("cd", tmp, ";", "unzip", foodzip[1]))
system(paste("cd", tmp, ";", "unzip", foodzip[2]))
system(paste("cd", tmp, ";", "mv *.paj paj/"))

pajfiles <- list.files(paste(tmp, sep="/", "paj"), full.names=TRUE)

readpaj <- function(filename) {
  lines <- readLines(filename)
  lines <- grep("^%", lines, invert=TRUE, value=TRUE)       # comments
  lines <- grep("^[ \t]*$", lines, invert=TRUE, value=TRUE) # empty lines

  eco <- lines[grep("^\\*partitio", lines)[1]:
               (grep("^\\*network", lines)[1]-1)]
  net <- lines[grep("^\\*network", lines)[1]:(grep("^\\*vector", lines)[1]-1)]
  bim <- lines[grep("^\\*vector", lines)[1]:length(lines)]

  tf <- tempfile()
  cat(net, file=tf, sep="\n")
  G <- read.graph(tf, format="pajek")
  V(G)$name <- V(G)$id
  G <- remove.vertex.attribute(G, "id")
  V(G)$ECO <- as.numeric(eco[-(1:2)])
  V(G)$Biomass <- as.numeric(bim[-(1:2)])
  G
}

foodwebs <- lapply(pajfiles, readpaj)
names(foodwebs) <- sub("\\.paj$", "", basename(pajfiles))

foodwebs <- foodwebs[setdiff(names(foodwebs), c("Everglades", "Florida"))]

authors <- c("ChesLower", "Hagy, J.D.",
             "ChesMiddle", "Hagy, J.D.",
             "ChesUpper", "Hagy, J.D.",
             "Chesapeake", "Baird, D. and R.E. Ulanowicz",
             "CrystalC", "Homer, M. and W.M. Kemp",
             "CryslalD", "Homer, M. and W.M. Kemp",
             "Maspalomas",
             "Almunia, J., G. Basterretxea, J. Aristegui, and R.E. Ulanowicz",
             "Michigan", "Krause, A. and D. Mason",
             "Mondego", "Patricio, J.",
             "Narragan", "Monaco, M.E. and R.E. Ulanowicz",
             "Rhode", "Correll, D",
             "StMarks", "Baird, D., J. Luczkovich and R. R. Christian",
             "baydry",
             "Ulanowicz, R. E., C. Bondavalli, and M. S. Egnotovich",
             "baywet",
             "Ulanowicz, R. E., C. Bondavalli, and M. S. Egnotovich",
             "cypdry",
             "Ulanowicz, R. E., C. Bondavalli, and M. S. Egnotovich",
             "cypwet",
             "Ulanowicz, R. E., C. Bondavalli, and M. S. Egnotovich",
             "gramdry",
             "Ulanowicz, R. E., C. Bondavalli, and M. S. Egnotovich",
             "gramwet",
             "Ulanowicz, R. E., C. Bondavalli, and M. S. Egnotovich",
             "mangdry",
             "Ulanowicz, R. E., C. Bondavalli, J. J. Heymans, and M. S. Egnotovich",
             "mangwet",
             "Ulanowicz, R. E., C. Bondavalli, J. J. Heymans, and M. S. Egnotovich"
             )                           

Authors <- matrix(authors, nc=2, byrow=TRUE)

citations <- 
"ChesLower,ChesMiddle,ChesUpper| Hagy, J.D. (2002) Eutrophication, hypoxia
                   and trophic transfer efficiency in Chesapeake Bay PhD
                   Dissertation, University of Maryland at College
                   Park (USA), 446 pp.
Chesapeake|        Baird D. & Ulanowicz R.E. (1989) The seasonal dynamics
                   of the Chesapeake Bay ecosystem. Ecological Monographs
                   59:329-364.
CrystalC,CrystalD| Homer, M. and W.M. Kemp. Unpublished Ms. See also
                   Ulanowicz, R.E. 1986. Growth and Development:
                   Ecosystems Phenomenology. Springer, New York. pp 69-79. 
Maspalomas|        Almunia, J., G. Basterretxea, J. Aristegui, and R.E.
                   Ulanowicz. (1999) Benthic- Pelagic switching in a coastal
                   subtropical lagoon. Estuarine, Coastal and Shelf
                   Science 49:363-384.
Michigan|          Krause, A. and D. Mason. (In preparation.) A. Krause,
                   PhD. Dissertation, Michigan State University.
                   Ann Arbor, MI. USA 
Mondego|           Patricio, J. (In Preparation) Master's Thesis.
                   University of Coimbra, Coimbra, Portugal. 
Narragan|          Monaco, M.E. and R.E. Ulanowicz. (1997) Comparative
                   ecosystem trophic structure of three U.S. Mid-Atlantic
                   estuaries. Mar. Ecol. Prog. Ser. 161:239-254.
Rhode|             Correll, D. (Unpublished manuscript) Smithsonian
                   Institute, Chesapeake Bay Center for Environmental
                   Research, Edgewater, Maryland 21037-0028 USA. 
StMarks|           Baird, D., J. Luczkovich and R. R. Christian. (1998)
                   Assessment of spatial and temporal variability in
                   ecosystem attributes of the St Marks National Wildlife
                   Refuge, Apalachee Bay, Florida. Estuarine, Coastal, and
                   Shelf Science 47: 329-349.
baydry,baywet|     Ulanowicz, R. E., C. Bondavalli, and M. S. Egnotovich.
                   1998. Network analysis of trophic dynamics in South
                   Florida ecosystems, FY 97: the Florida Bay ecosystem.
                   Annual Report to the United States Geological Service
                   Biological Resources Division, University of Miami Coral
                   Gables, [UMCES] CBL 98-123, Maryland System Center for
                   Environmental Science, Chesapeake Biological Laboratory,
                   Maryland, USA.
cypdry,cypwet|     Ulanowicz, R. E., C. Bondavalli, and M. S. Egnotovich.
                   1997. Network analysis of trophic dynamics in South
                   Florida ecosystems, FY 96: the cypress wetland ecosystem.
                   Annual Report to the United States Geological Service
                   Biological Resources Division, University of Miami Coral
                   Gables, [UM-CES] CBL 97-075, Maryland System Center for
                   Environmental Science, Chesapeake Biological Laboratory.
gramdry,gramwet|   Ulanowicz, R. E., J. J. Heymans, and M. S. Egnotovich.
                   2000. Network analysis of trophic dynamics in South
                   Florida ecosystems, FY 99: the graminoid ecosystem.
                   Technical Report TS-191-99, Maryland System Center for
                   Environmental Science, Chesapeake Biological Laboratory,
                   Maryland, USA.
mangdry,mangwet|   Ulanowicz, R. E., C. Bondavalli, J. J. Heymans, and
                   M. S. Egnotovich. 1999. Network analysis of trophic
                   dynamics in South Florida ecosystems, FY 98: the mangrove
                   ecosystem. Technical Report TS-191-99, Maryland System
                   Center for Environmental Science, Chesapeake Biological
                   Laboratory, Maryland, USA."


Citations <- readLines(textConnection(citations))
Citations2 <- Citations[1]
for (i in 2:length(Citations)) {
  if (grepl("^[ ]", Citations[i])) {
    Citations2[length(Citations2)] <- paste(Citations2[length(Citations2)],
                                            sub("^[ ]*", "", Citations[i]))
  } else {
    Citations2 <- c(Citations2, Citations[i])
  }
}
Citations2 <- strsplit(Citations2, split="|", fixed=TRUE)
ids <- lapply(Citations2, function(x) strsplit(x[1], ",")[[1]])
cits <- sub("^[ ]*", "", sapply(Citations2, "[[", 2))
Citations2 <- cbind(unlist(ids), rep(cits, sapply(ids, length)))

url <- "http://vlado.fmf.uni-lj.si/pub/networks/data/bio/foodweb/foodweb.htm"

name <-
"ChesLower|         Lower Chesapeake Bay in Summer
ChesMiddle|         Middle Chesapeake Bay in Summer
ChesUpper|          Upper Chesapeake Bay in Summer
Chesapeake|         Chesapeake Bay Mesohaline Network
CrystalC|           Crystal River Creek (Control)
CrystalD|           Crystal River Creek (Delta Temp)
Maspalomas|         Charca de Maspalomas
Michigan|           Lake Michigan Control network
Mondego|            Mondego Estuary - Zostrea site
Narragan|           Narragansett Bay Model
Rhode|              Rhode River Watershed - Water Budget
StMarks|            St. Marks River (Florida) Flow network
baydry|             Florida Bay Trophic Exchange Matrix, dry season
baywet|             Florida Bay Trophic Exchange Matrix, wet season
cypdry|             Cypress Dry Season
cypwet|             Cypress Wet Season
gramdry|            Everglades Graminoids - Dry Season
gramwet|            Everglades Graminoids - Wet Season
mangdry|            Mangrove Estuary, Dry Season
mangwet|            Mangrove Estuary, Wet Season"

Name <- read.delim(textConnection(name), sep="|", header=FALSE)
Name[,2] <- sub("^[ ]*", "", Name[,2])

for (n in names(foodwebs)) {
  foodwebs[[n]]$Citation <- Citations2[,2][match(n, Citations2[,1])]
  foodwebs[[n]]$Author   <- Authors[,2][match(n, Authors[,1])]
  foodwebs[[n]]$URL      <- url
  foodwebs[[n]]$name     <- Name[,2][match(n, Name[,1])]
}
  
save(foodwebs, file="/tmp/foodwebs.rda")

#####################################################################

## Konigsberg

library(igraph)

edges <- '
from,to,Euler_letter,name
Altstadt-Loebenicht,Kneiphof,a,Kraemer Bruecke
Altstadt-Loebenicht,Kneiphof,b,Schmiedebruecke
Altstadt-Loebenicht,Lomse,f,Holzbruecke
Kneiphof,Lomse,e,Honigbruecke
Vorstadt-Haberberg,Lomse,g,Hohe Bruecke
Vorstadt-Haberberg,Kneiphof,c,Gruene Bruecke
Vorstadt-Haberberg,Kneiphof,d,Koettelbruecke'

vertices <- "
name,Euler_letter
Altstadt-Loebenicht,B
Kneiphof,A
Vorstadt-Haberberg,C
Lomse,D"

Koenigsberg <- graph.data.frame(read.csv(textConnection(edges)),
                                vertices=read.csv(textConnection(vertices)),
                                directed=FALSE)

Koenigsberg$name <- "The seven bidges of Koenigsberg"

save(Koenigsberg, file="/tmp/Koenigsberg.rda")

