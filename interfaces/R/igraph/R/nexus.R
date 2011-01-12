
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

  attr <- nam=="Attribute"
  myattr <- unlist(dsi[attr])
  dsi <- dsi[!attr]
  nam <- nam[!attr]
  names(dsi) <- nam
  class(dsi) <- "nexusDatasetInfo"

  if (length(myattr) != 0) {
    myattr <- strsplit(myattr, "\n", fixed=TRUE)
    attrdat <- lapply(myattr, function(x) strsplit(x[1], " ")[[1]])
    myattr <- sapply(myattr, "[", 2)
    dsi$Attributes <- mapply(attrdat, myattr, SIMPLIFY=FALSE,
                             FUN=function(dat, desc) {
                               list(Type=dat[1], Datatype=dat[2], Name=dat[3],
                                    Description=desc)
                             })
  }
  
  dsi$Id <- as.numeric(dsi$Id)
  dsi$Vertices <- as.numeric(dsi$Vertices)
  dsi$Edges <- as.numeric(dsi$Edges)
  dsi$Tags <- strsplit(dsi$Tags, ";", fixed=TRUE)[[1]]
  
  dsi
}

print.nexusDatasetInfo <- function(x, ...) {
  cat(sep="", "[Nexus dataset # ", x$Id, ", ", x$Vertices, "/",
      x$Edges, ", ", x$Name, "]\n")
  if (length(x$Tags) != 0) {
    cat("Tags: ")
    cat(x$Tags, sep="; ")
    cat("\n")
  }
  attr <- x[["Attributes"]]
  printed <- c("Id", "Vertices", "Edges", "Name", "Tags", "Attributes")
  x <- x[ setdiff(names(x), printed) ]
  for (a in attr) {
    cat(paste(sep="", "Attribute '", a$Name, "', ", a$Type,
              ", ", a$Datatype, ":"), "\n")
    cat(strwrap(a$Description, initial="  ", prefix="  "), sep="\n")
  }
  for (i in names(x)) {
    cat(strwrap(paste(sep="", i, ": ", x[[i]]), initial="", prefix="  "),
        sep="\n")
  }
  invisible(x)
}

summary.nexusDatasetInfoList <- function(object, ...) {
  cat(paste(sep="", "Nexus data set information list, results ",
            as.numeric(attr(object, "Offset")) + 1, "-",
            as.numeric(attr(object, "Offset")) +
            as.numeric(attr(object, "Size")),
            " out of ", attr(object, "Totalsize"), ".\n\n"))
  invisible(object)
}

print.nexusDatasetInfoList <- function(x, ...) {
  summary(x)
  attributes(x) <- NULL
  print(x)
}

nexus.list <- function(tags=NULL, offset=0, limit=10,
                       operator=c("or", "and"),
                       order=c("date", "name", "popularity"),
                       nexus.url=getIgraphOpt("nexus.url")) {

  operator=match.arg(operator)
  order=match.arg(order)
  
  if (is.null(tags)) {
    u <- paste(sep="", nexus.url, "/api/dataset_info?format=text",
               "&offset=", offset, "&limit=", limit, "&order=", order)
  } else {
    tags <- paste(tags, collapse="|")
    u <- paste(sep="", nexus.url, "/api/dataset_info?tag=", tags,
               "&operator=", operator, "&format=text",
               "&offset=", offset, "&limit=", limit, "&order=", order)
  }
  f <- url(URLencode(u))
  l <- readLines(f)
  close(f)

  if (length(l)==0) {
    return(list())
  }

  l <- lapply(l, function(x) c(sub("[ ]*:[^:]*$", "", x),
                               sub("^[^:]*:[ ]*", "", x)))
  spos <- which(sapply(l, function(x) x[1]=="Id"))
  epos <- c((spos-1), length(l))
  ehead <- epos[1]
  epos <- epos[-1]
  
  res <- mapply(spos, epos, SIMPLIFY=FALSE, FUN=function(s, e)
                makeNexusDatasetInfo(l[s:e]))
  class(res) <- "nexusDatasetInfoList"

  for (h in 1:ehead) {
    attr(res, l[[h]][1]) <- l[[h]][2]
  }
  
  res
}

nexus.info <- function(id, nexus.url=getIgraphOpt("nexus.url")) {
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
  makeNexusDatasetInfo(l2)
}  

nexus.get <- function(id, nexus.url=getIgraphOpt("nexus.url")) {
  u <- paste(sep="", nexus.url, "/api/dataset?id=", id, "&format=R-igraph")
  env <- new.env()
  rdata <- url(URLencode(u))
  load(rdata, envir=env)
  close(rdata)
  return(get(ls(env)[1], env))
}
