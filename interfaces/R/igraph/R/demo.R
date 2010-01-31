
igraphdemo <- function(which) {
  require(igraph)
  require(tcltk)
  
  if (missing(which)) {
    demodir <- system.file("demo", package="igraph")
    if (demodir=="") {
      stop("Could not find igraph demos, broken igraph installation?")
    }
    return( sub("\\.R$", "", list.files(demodir)) )
  }

  if (!grepl("\\.R$", which)) {
    which <- paste(which, sep=".", "R")
  }

  if (!file.exists(which) && ! grepl("^/", which)) {
    which <- system.file( paste("demo", sep="/", which), package="igraph" )
  }

  if (which=="" || !file.exists(which)) {
    stop("Could not find demo file")
  }
  
  .igraphdemo.next <- function(top, txt) {
    act <- as.character(tktag.nextrange(txt, "active", "0.0"))
    if (length(act)==0) {
      return()
    }

    options(keep.source=TRUE)
    
    text <- tclvalue(tkget(txt, act[1], act[2]))
    cat("=======================================================\n");

    expr <- parse(text=text)
    for (i in seq_along(expr)) {
      co <- as.character(attributes(expr)$srcref[[i]])
      co[1] <- paste("> ", sep="", co[1])
      if (length(co)>1) {
        co[-1] <- paste(" +", sep="", co[-1])
      }
      cat(co, sep="\n")
      res <- withVisible(eval(expr[[i]], envir=.GlobalEnv))
      if (res$visible) {
        print(res$value)
      }
    }
    cat("> -------------------------------------------------------\n");
    cat(options()$prompt)
    
    tktag.remove(txt, "activechunk", act[1], act[2])
    tktag.remove(txt, "active", act[1], act[2])

    nex <- as.character(tktag.nextrange(txt, "activechunk", act[1]))
    if (length(nex)!=0) {
      tktag.add(txt, "active", nex[1], nex[2])
      tksee(txt, paste(sep="", as.numeric(nex[2]), ".0"))
      tksee(txt, paste(sep="", as.numeric(nex[1]), ".0"))
    }
  }
  
  .igraphdemo.close <- function(top) {
    tkdestroy(top)
  }
  
  .igraphdemo.reset <- function(top, txt, which) {
    demolines <- readLines(which)
    demolines <- demolines[!grepl("^pause\\(\\)$", demolines)]
    demolines <- paste(" ", sep="", demolines)

    ch <- grep("^[ ]*###", demolines)
    ch <- c(ch, length(demolines)+1)
    if (length(ch)==1) {
      warning("Demo source file does not contain chunks")
    } else {
      demolines <- demolines[ch[1]:length(demolines)]
      ch <- grep("^[ ]*###", demolines)
      ch <- c(ch, length(demolines)+1)
    }
    
    tkconfigure(txt, state="normal")
    tkdelete(txt, "0.0", "end")
    tkinsert(txt, "insert", paste(demolines, collapse="\n"))
    tkconfigure(txt, state="disabled")

    for (i in seq_along(ch[-1])) {
      from <- paste(sep="", ch[i], ".0")
      to <- paste(sep="", ch[i+1]-1, ".0")
      tktag.add(txt, "chunk", from, to)
      tktag.add(txt, "activechunk", from, to)
    }
    tktag.configure(txt, "chunk", "-borderwidth", "1")
    tktag.configure(txt, "chunk", "-relief", "sunken")
    if (length(ch) >= 2) {
      tktag.add(txt, "active", paste(sep="", ch[1], ".0"),
                paste(sep="", ch[2]-1, ".0"))
      tktag.configure(txt, "active", "-foreground", "red")
      tktag.configure(txt, "active", "-background", "lightgrey")
    }

    comm <- grep("^#", demolines)
    for (i in comm) {
      tktag.add(txt, "comment", paste(sep="", i, ".0"),
                paste(sep="", i, ".end"))
    }
    tktag.configure(txt, "comment", "-font", "bold")
    tktag.configure(txt, "comment", "-foreground", "darkolivegreen")
  }

  top <- tktoplevel(background="lightgrey")
  tktitle(top) <- paste("igraph demo:", which)
  
  main.menu <- tkmenu(top)
  tkadd(main.menu, "command", label="Close", command=function()
        .igraphdemo.close(top))
  tkadd(main.menu, "command", label="Reset", command=function()
        .igraphdemo.reset(top, txt, which))
  tkconfigure(top, "-menu", main.menu)

  scr <- tkscrollbar(top, repeatinterval=5,
                     command=function(...) tkyview(txt,...))
  txt <- tktext(top, yscrollcommand=function(...) tkset(scr, ...),
                width=80, height=40)
  but <- tkbutton(top, text="Next", command=function()
                  .igraphdemo.next(top, txt))
  
  tkpack(but, side="bottom", fill="x", expand=0)
  tkpack(scr, side="right", fill="y", expand=0)
  tkpack(txt, side="left", fill="both", expand=1)

  .igraphdemo.reset(top, txt, which)
  
  invisible()
}
