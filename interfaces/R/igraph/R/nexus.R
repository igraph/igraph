
#   IGraph R package
#   Copyright (C) 2011  Gabor Csardi <csardi.gabor@gmail.com>
#   Rue de l'Industrie 5, 1005 Lausanne, Switzerland
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

makeNexusDatasetInfo <- function(entries) {
  dsi <- lapply(entries, "[", 2)
  nam <- sapply(entries, "[", 1)

  attr <- nam=="attribute"
  myattr <- unlist(dsi[attr])
  dsi <- dsi[!attr]
  nam <- nam[!attr]
  names(dsi) <- nam
  class(dsi) <- "nexusDatasetInfo"

  if (length(myattr) != 0) {
    myattr <- strsplit(myattr, "\n", fixed=TRUE)
    attrdat <- lapply(myattr, function(x) strsplit(x[1], " ")[[1]])
    myattr <- sapply(myattr, "[", 2)
    dsi$attributes <- mapply(attrdat, myattr, SIMPLIFY=FALSE,
                             FUN=function(dat, desc) {
                               list(type=dat[1], datatype=dat[2], name=dat[3],
                                    description=desc)
                             })
  }
  
  dsi$id <- as.numeric(dsi$id)
  dsi$tags <- strsplit(dsi$tags, ";", fixed=TRUE)[[1]]
  
  dsi
}

print.nexusDatasetInfo <- function(x, ...) {
  ve <- strsplit(parseVE(x$`vertices/edges`), "/")[[1]]
  nc <- c("U", "-", "-", "-")
  if ("directed" %in% x$tags && "undirected" %in% x$tags) {
    nc[1] <- "B"
  } else if ("directed" %in% x$tags) {
    nc[1] <- "D"
  }
  if (is.null(x$attributes)) {
    nc[2] <- "?"
  } else if (any(sapply(x$attributes,
                        function(X) X$name=="name" && X$type=="vertex"))) {
    nc[2] <- "N"
  } 
  if ("weighted" %in% x$tags) {
    nc[3] <- "W"
  }
  if ("bipartite" %in% x$tags) {
    nc[4] <- "B"
  }
  nc <- paste(nc, collapse="")
  head <- paste(sep="", "NEXUS ", nc, " ", ve[1], " ", ve[2], " -- #",
                x$id, " ", x$sid, " ", x$name)
  if (nchar(head) > getOption("width")) {
    head <- paste(sep="", substr(head, 1, getOption("width")-1), "+")
  }
  cat(head, sep="", "\n")
  if (length(x$tags) != 0) {
    tt <- strwrap(paste(sep="", "+ tags:", paste(x$tags, collapse="; ")),
                  initial="", prefix="  ")
    cat(tt, sep="\n")
  }
  if ("networks" %in% names(x)) {
    nets <- strsplit(x$networks, " ")[[1]]
    nn <- strwrap(paste(sep="", "+ nets:", paste(nets, collapse="; ")),
                  initial="", prefix="  ")
    cat(nn, sep="\n")
  }
  attr <- x[["attributes"]]
  printed <- c("id", "sid", "vertices/edges", "name", "tags", "networks",
               "attributes")
  x <- x[ setdiff(names(x), printed) ]
  if (length(attr)>0) {
    dcode <- function(d) {
      if (d=="numeric") return("n")
      if (d=="string") return("c")
      "x"
    }
    cat("+ attr: ")
    astr <- sapply(attr, function(a) {
      paste(sep="", a$name, " (", substr(a$type, 1, 1), "/",
            dcode(a$datatype), ")")
    })
    cat(strwrap(paste(astr, collapse=", "), exdent=2), "\n")
  }
  for (i in names(x)) {
    xx <- strsplit(x[[i]], "\n")[[1]]
    ff <- strwrap(paste(sep="", "+ ", i, ": ", xx[1]), initial="",
                  prefix="  ")
    xx <- unlist(sapply(xx[-1], strwrap, prefix="  "))
    cat(ff, sep="\n")
    if (length(xx)>0) {
      cat(xx, sep="\n")
    }
  }
  invisible(x)
}

summary.nexusDatasetInfoList <- function(object, ...) {
  o <- as.numeric(attr(object, "offset"))
  s <- as.numeric(attr(object, "size"))
  t <- as.numeric(attr(object, "totalsize"))
  n <- attr(object, "name")
  cat(sep="", "NEXUS ", o+1, "-", o+s, "/", t, " -- ", n, "\n")
  invisible(object)
}

parseVE <- function(ve) {
  if (length(ve)==0) { return(character(0)) }
  ve <- strsplit(unname(ve), " ")
  ve <- lapply(ve, strsplit, "/")
  v <- lapply(ve, function(x) sapply(x, "[", 1))
  e <- lapply(ve, function(x) sapply(x, "[", 2))
  int <- function(x) {
    if (length(unique(x))==1) {
      as.character(x[1])
    } else {
      paste(sep="", min(x), "-", max(x))
    }
  }
  v <- sapply(v, int)
  e <- sapply(e, int)
  paste(v, sep="/", e)
}

print.nexusDatasetInfoList <- function(x, ...) {
  summary(x)
  
  if (length(x)==0) { return(invisible(x)) }
  
  ve <- parseVE(unname(sapply(x, "[[", "vertices/edges")))
  nets <- sapply(x, function(y) length(strsplit(y$networks, " ")[[1]]))
  sid <- sapply(x, "[[", "sid")
  if (any(nets>1)) {
    sid[nets > 1] <- paste(sep="", sid[nets>1], ".", nets[nets>1])
  }
  df <- data.frame(no=paste(sep="", "[", format(seq_along(x)), "] "),
                   sid=format(sid),
                   size=paste(sep="", " ", format(ve)),
                   id=paste(sep="", " #", format(sapply(x, "[[", "id")), " "),
                   name=sapply(x, "[[", "name"))
  out <- do.call(paste, c(as.list(df), sep=""))
  long <- nchar(out) > getOption("width")  
  out <- paste(sep="", substr(out, 1, getOption("width")-1),
               ifelse(long, "+", ""))
  cat(out, sep="\n")
  invisible(x)
}

nexus.format.result <- function(l, name="") {
  
  if (length(l)==0) {
    res <- list()
    class(res) <- "nexusDatasetInfoList"
    return(res)
  }
  
  l <- lapply(l, function(x) c(sub("[ ]*:[^:]*$", "", x),
                               sub("^[^:]*:[ ]*", "", x)))
  spos <- which(sapply(l, function(x) x[1]=="id"))
  epos <- c((spos-1), length(l))
  ehead <- epos[1]
  epos <- epos[-1]
  
  res <- mapply(spos, epos, SIMPLIFY=FALSE, FUN=function(s, e)
                makeNexusDatasetInfo(l[s:e]))
  class(res) <- "nexusDatasetInfoList"

  for (h in 1:ehead) {
    attr(res, l[[h]][1]) <- l[[h]][2]
    attr(res, "name") <- name
  }
  
  res
}

nexus.list <- function(tags=NULL, offset=0, limit=10,
                       operator=c("or", "and"),
                       order=c("date", "name", "popularity"),
                       nexus.url=getIgraphOpt("nexus.url")) {

  operator=igraph.match.arg(operator)
  order=igraph.match.arg(order)
  
  if (is.null(tags)) {
    u <- paste(sep="", nexus.url, "/api/dataset_info?format=text",
               "&offset=", offset, "&limit=", limit, "&order=", order)
    name <- "data set list"
  } else {
    tags <- paste(tags, collapse="|")
    u <- paste(sep="", nexus.url, "/api/dataset_info?tag=", tags,
               "&operator=", operator, "&format=text",
               "&offset=", offset, "&limit=", limit, "&order=", order)
    name <- paste("tags:", gsub("|", "; ", tags, fixed=TRUE))
  }
  f <- url(URLencode(u))
  l <- readLines(f)
  close(f)

  nexus.format.result(l, name)
}

nexus.info <- function(id, nexus.url=getIgraphOpt("nexus.url")) {

  if (inherits(id, "nexusDatasetInfo")) {
    id <- id$id
  } else if (inherits(id, "nexusDatasetInfoList")) {
    rid <- sapply(id, "[[", "id")
    res <- lapply(rid, nexus.info, nexus.url=nexus.url)
    class(res) <- class(id)
    attributes(res) <- attributes(id)
    return(res)
  }  
  
  u <- paste(sep="", nexus.url, "/api/dataset_info?format=text&id=", id)
  f <- url(URLencode(u))
  l <- readLines(f)
  close(f)
  l2 <- character()
  for (i in seq_along(l)) {
    if (!grepl("^  ", l[i])) {
      l2 <- c(l2, l[i])
    } else {
      l2[length(l2)] <- paste(sep="\n", l2[length(l2)],
                              sub("  ", "", l[i], fixed=TRUE))
    }
  }
  l2 <- lapply(l2, function(x)
               c(sub("[ ]*:.*$", "", x), sub("^[^:]*:[ ]*", "", x)))
  res <- makeNexusDatasetInfo(l2)
  if (! "attributes" %in% names(res)) { res$attributes <- list() }
  return(res)
}  

nexus.get <- function(id, offset=0,
                      order=c("date", "name", "popularity"),
                      nexus.url=getIgraphOpt("nexus.url")) {

  order=igraph.match.arg(order)

  if (inherits(id, "nexusDatasetInfo")) {
    id <- id$id
  } else if (inherits(id, "nexusDatasetInfoList")) {
    id <- sapply(id, "[[", "id")
    return(lapply(id, nexus.get, nexus.url=nexus.url))
  }
  
  u <- paste(sep="", nexus.url, "/api/dataset?id=", id, "&format=R-igraph")
  env <- new.env()
  rdata <- url(URLencode(u))
  load(rdata, envir=env)
  close(rdata)
  return(get(ls(env)[1], env))
}

nexus.search <- function(q, offset=0, limit=10,
                         order=c("date", "name", "popularity"),
                         nexus.url=getIgraphOpt("nexus.url")) {

  order=igraph.match.arg(order)

  u <- paste(sep="", nexus.url, "/api/search?q=", q,
             "&format=text","&offset=", offset, "&limit=", limit,
             "&order=", order)
  f <- url(URLencode(u))
  l <- readLines(f)
  close(f)

  if (length(l)==0) {
    res <- list()
    class(res) <- "nexusDatasetInfoList"
    return(res)
  }

  nexus.format.result(l, name=paste("q:", q))
}

`[.nexusDatasetInfoList` <- function(x, i) {
  res <- unclass(x)[i]
  class(res) <- class(x)
  attributes(res) <- attributes(x)
  res
}

'
DATA SET LIST:
--------------

NEXUS 1-10/18 -- data set list
[ 1] kaptail.4         #18 39/109-223   Kapferer tailor shop
[ 2] condmatcollab2003 #17 31163/120029 Condensed matter collaborations, 2003
[ 3] condmatcollab     #16 16726/47594  Condensed matter collaborations, 1999
[ 4] powergrid         #15 4941/6594    Western US power grid
[ 5] celegansneural    #14 297/2359     C. Elegans neural network
[ 6] polblogs          #13 1490/19090   US political blog network
[ 7] dolphins          #12 62/159       Dolphin social network
[ 8] football          #11 115/616      Network of American college ...
[ 9] adjnoun           #10 112/425      Word adjacencies from David ...
[10] huckleberry       # 9 74/301       Coappearance network from ...


TAG SEARCH:
-----------

NEXUS 1-4/4 -- tags: directed
[1] kaptail.4 #18 39/109-223 Kapferer tailor shop
[2] polblogs  #13 1490/19090 US political blog network
[3] macaque   # 4 45/463     Macaque visuotactile brain areas
[4] UKfaculty # 2 81/817     UK faculty social network


FULL TEXT SEARCH:
-----------------
  
NEXUS 1-2/2 -- q: US
[1] powergrid #15 4941/6594  Western US power grid
[2] polblogs  #13 1490/19090 US political blog network


DATA SET SUMMARY:
-----------------
  
NEXUS B--- 39 109-223 -- #18 Kapferer tailor shop
+ tags: directed; social network; undirected
+ networks: 1/KAPFTI2; 2/KAPFTS2; 3/KAPFTI1; 4/KAPFTS1

NEXUS U--- 4941 6594 -- #15 Western US power grid
+ tags: technology

DATA SET INFO:
--------------

NEXUS B--- 39 109-223 -- #18 Kapferer tailor shop
+ tags: directed; social network; undirected
+ attr: name (v/c) [Actor names]
+ networks: 1/KAPFTI2; 2/KAPFTS2; 3/KAPFTI1; 4/KAPFTS1
+ nets: #1 KAPFTI2; #2 KAPFTS2; #3 KAPFTI1; #4 KAPFTS1
+ date: 2011-01-23
+ licence: Creative Commons by-sa 3.0
+ licence url: http://creativecommons.org/licenses/by-sa/3.0/
+ summary: Interactions in a tailor shop in Zambia (then
  Northern Rhodesia) over a period of ten months.
+ details: Bruce Kapferer (1972) observed interactions in a tailor
  shop in Zambia (then Northern Rhodesia) over a period of ten months.
  His focus was the changing patterns of alliance among workers during
  extended negotiations for higher wages. . The matrices represent two
  different types of interaction, recorded at two different times
  (seven months apart) over a period of one month. TI1 and TI2 record
  the "instrumental" (work- and assistance-related) interactions at the
  two times; TS1 and TS2 the "sociational" (friendship, socioemotional)
  interactions. . The data are particularly interesting since an
  abortive strike occurred after the first set of observations, and a
  successful strike took place after the second.
+ formats: Pajek; R-igraph
+ citation: Kapferer B. (1972). Strategy and transaction in an African
  factory. Manchester: Manchester University Press. 
'
