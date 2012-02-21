
base <- "~/Downloads/enron/maildir"

## Check that every maildir has a _sent-mail folder, with at least
## one message

sentFiles <- c("sent", "sent_items", "_sent_mail", "chris_stokley/sent")

checkSentMail <- function(x) {
  d <- paste(sep="", base, "/", x)
  f <- list.files(d)
  sf <- ifelse(grepl("/", sentFiles), dirname(sentFiles), sentFiles)
  any(sf %in% f)
}

dirs <- list.files(base)
sentmail <- sapply(dirs, checkSentMail)

## One person did not have a sent-mail folder, (s)he did not send
## any mail, apparently

getEmailAddress <- function(x) {
  d <- paste(sep="", base, "/", x)

  for (s in sentFiles) {
    f <- paste(sep="", d, "/", s)
    if (file.exists(f)) {
      mailfile <- list.files(f, full.names=TRUE)[1]
      l <- readLines(mailfile, warn=FALSE)
      from <- l[grep("^From:", l)[1]]
      return(sub(".*[ :]([^: ]+@[^ >]+).*", "\\1", from))
    }
  }
  as.character(NA)
}

email <- sapply(dirs, getEmailAddress)
email[ is.na(email) ] <- "steven.harris@enron.com"

getAllEdges <- function(x) {
  print(x)
  folders <- list.files(paste(sep="", base, "/", x), full.names=TRUE)
  e <- list()
  for (f in folders) {
    mailfiles <- list.files(f, full.names=TRUE)
    e <- c(e, lapply(mailfiles, getEdges))
  }
  e
}

getEdges <- function(x) {
  cat('.')

  addl <- function(i) {
    r <- l[i]
    i <- i+1
    while (grepl("^[ \t]", l[i])) {
      r <- c(r, l[i])
      i <- i+1
    }
    r <- paste(r, collapse=" ")
    r <- strsplit(r, ",?[ \t]+")[[1]][-1]
  }
  
  l <- readLines(x, warn=FALSE)
  from <- addl(grep("^From:", l)[1])
  to <- addl(grep("^To:", l)[1])
  cc <- addl(grep("^Cc:", l)[1])
  bcc <- addl(grep("^Bcc:", l)[1])
  date <- l[grep("^Date:", l)[1]]
  
  bcc <- setdiff(bcc, c(from, to, cc))
  cc <- setdiff(cc, c(from, to))
  to <- setdiff(to, from)
  
  list(from=from, date=date, to=to, cc=cc, bcc=bcc)
}

edges <- lapply(dirs, getAllEdges)
edges <- unlist(edges, recursive=FALSE)

## Try to filter out messages appearing multiple times...

encode <- function(e) {
  f <- e$from
  d <- e$date
  t <- paste(sort(e$to), collapse=";")
  paste(f, d, t, sep=";")
}

edgecodes <- sapply(edges, encode)
dup <- duplicated(edgecodes)

uniEdges <- edges[!dup]

people <- unique(c(unique(unlist(sapply(uniEdges, "[[", "from"))),
                   unique(unlist(sapply(uniEdges, "[[", "to"))),
                   unique(unlist(sapply(uniEdges, "[[", "cc")))))

people <- sort(people)

## This is a big mess, it should be cleaned up, but
## of course I don't have the time for this....
