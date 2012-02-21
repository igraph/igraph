
urls <- paste(sep="",
              "http://w3.usf.edu/FreeAssociation/AppendixA/Cue_Target_Pairs.",
              c("A-B", "C", "D-F", "G-K", "L-O", "P-R", "S", "T-Z"))

tmp <- tempdir()
for (u in seq_along(urls)) {
  dest <- paste(sep="", tmp, "/", u, ".dat")
  download.file(urls[u], dest)
}

tables <- lapply(seq_along(urls), function(u) {
  dest <- paste(sep="", tmp, "/", u, ".dat")
  lines <- gsub("\xa5", "", readLines(dest))
  read.csv(textConnection(lines), comment.char="<")
})

## Columns:
## CUE      edge source
## TARGET   edge target
## NORMED.  YES/NO
## X.G      number of participant serving in the group norming the word
## X.P      number of participants producing the response
## FSG      forward strength, X.P / X.G
## BSG      backward strength
## MSG      mediated strength (two-step strength)
## OSG      overlapping strength
## X.M      number of mediated connections linking cue and target
## MMIAS    number of potential mediators that were not normed
## X.O      number of overlapping associates shared
## OMIAS    overlaps not normed
## QSS      number of neighbors of cue
## QFR      printed frequency of the cue
## QCON     concreteness rating, 1-7, many missing
## QH       homograph, according to various sources
## QPS      part of the speech classification
## QMC      mean connectivity among the associates of the normed word
## QPR      probability that each associate in the set produces
##          the normed cue as an associate
## QRSG     resonance strength of the cue
## QUC      Use Code value for the cue, whether all important associated were
##          normed
## TSS      The rest is the same, but for the target word
## TFR
## TCON
## TH
## TPS
## TMC
## TPR
## TRSG

## First put the tables together

for (i in seq_along(tables)) {
  colnames(tables[[i]]) <- colnames(tables[[1]])
}
tables <- do.call(rbind, tables)

edges <- data.frame(from=tables$CUE, to=sub("^[ ]+", "", tables$TARGET),
                    Normed=as.numeric(tables$NORMED.=="YES"),
                    GroupSize=tables$X.P,
                    NormingSize=tables$X.G)

vertices1 <- data.frame(id=tables$CUE, PrintedFreq=tables$QFR,
                        Concreteness=tables$QCON,
                        Homograph=sub("^[ ]+", "", tables$QH),
                        PartOfSpeech=sub("^[ ]+", "", tables$QPS))
vertices1 <- vertices1[ !duplicated(vertices1), ]

vertices2 <- data.frame(id=sub("^[ ]+", "", tables$TARGET),
                        PrintedFreq=tables$TFR,
                        Concreteness=tables$TCON,
                        Homograph=sub("^[ ]+", "", tables$TH),
                        PartOfSpeech=sub("^[ ]+", "", tables$TPS))
vertices2 <- vertices2[ !duplicated(vertices2), ]                        

vertices <- rbind(vertices1, vertices2)
vertices <- vertices[ !duplicated(vertices), ]
vertices <- vertices[ order(vertices$id), ]

## Fix errors

eqorna <- function(x, y) {
  if (is.na(x) && is.na(y)) {
    TRUE
  } else if ( (is.na(x) && !is.na(y)) || (!is.na(x) && is.na(y)) ) {
    FALSE
  } else {
    x==y
  }
}

## Word is twice, because can be different parts of speech

dupids <- vertices$id[duplicated(vertices$id)]
for (id in dupids) {
  rr <- vertices[vertices$id==id,]
  nps <- setdiff(colnames(vertices), "PartOfSpeech")
  eq <- sapply(nps, function(x) eqorna(rr[[x]][1], rr[[x]][2]))
  if (all(eq)) {
    vertices[vertices$id==id,]$PartOfSpeech <- rr$PartOfSpeech[1]
  }
}
vertices <- vertices[ !duplicated(vertices), ]

## Homograph entry is different

hord <- c("N", "P", "W", "T", "G", "C", "")
dupids <- vertices$id[duplicated(vertices$id)]
for (id in dupids) {
  rr <- vertices[vertices$id==id,]
  nh <- setdiff(colnames(vertices), "Homograph")
  eq <- sapply(nh, function(x) eqorna(rr[[x]][1], rr[[x]][2]))
  if (all(eq)) {
    h <- rr$Homograph    
    vertices[vertices$id==id,]$Homograph <- hord[which(hord %in% h)[1]]
  }
}
vertices <- vertices[ !duplicated(vertices), ]

## BANDANNA

vertices[ vertices$id=="BANDANNA", ]$PrintedFreq <- 0
vertices[ vertices$id=="BANDANNA", ]$PartOfSpeech <- "N"
vertices <- vertices[ !duplicated(vertices), ]

## Final check

dupids <- vertices$id[duplicated(vertices$id)]
if (length(dupids) != 0) { stop("GEBASZ") }

library(igraph)
freeassoc <- graph.data.frame(edges, vertices=vertices, directed=TRUE)

freeassoc$name <- "University of South Florida Free Association Norms"
freeassoc$Author <- "Douglas L. Nelson and Cathy L. McEvoy"
freeassoc$Citation <- "Nelson, D. L., McEvoy, C. L., & Schreiber, T. A. (1998). The University of South Florida word association, rhyme, and word fragment norms. http://www.usf.edu/FreeAssociation/."
freeassoc$URL <- "http://w3.usf.edu/FreeAssociation"

save(freeassoc, file="/tmp/freeassoc.Rdata")

