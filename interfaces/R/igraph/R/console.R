
#   IGraph R package
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

.igraph.pb <- NULL

.igraph.progress <- function(percent, message, clean=FALSE) {
  if (clean) {
    if (!is.null(.igraph.pb)) { close(.igraph.pb) }
    return(invisible())
  }
  type <- getIgraphOpt("verbose")
  if (is.logical(type) && type) {
    .igraph.progress.txt(percent, message)
  } else {
    switch (type,
            "tk"=.igraph.progress.tk(percent, message),
            "tkconsole"=.igraph.progress.tkconsole(percent, message),
            stop("Cannot interpret 'verbose' option, this should not happen"))
  }
}

.igraph.status <- function(message) {
  type <- getIgraphOpt("verbose")
  if (is.logical(type) && type) {
    message(message, appendLF=FALSE)
  } else {
    switch(type,
           "tk"=message(message, appendLF=FALSE),
           "tkconsole"=.igraph.progress.tkconsole.message(message, start=TRUE),
           stop("Cannot interpret 'verbose' option, this should not happen"))
  }
  0L
}

.igraph.progress.txt <- function(percent, message) {
  pb <- get(".igraph.pb", asNamespace("igraph"))
  if (percent==0) {
    if (!is.null(pb)) { close(pb) }
    cat(sep="", "  ", message, "\n")
    pb <- txtProgressBar(min=0, max=100, style=3)
  }
  setTxtProgressBar(pb, percent)
  if (percent==100) {
    close(pb);
    pb <- NULL
  }
  assign(".igraph.pb", pb, env=asNamespace("igraph"))
  0L
}

.igraph.progress.tk <- function(percent, message) {
  pb <- get(".igraph.pb", asNamespace("igraph"))
  if (percent==0) {
    if (!is.null(pb)) { close(pb) }
    pb <- tkProgressBar(min=0, max=100, title=message, label="0 %")
  }
  setTkProgressBar(pb, percent, label=paste(percent, "%"))
  if (percent==100) {
    close(pb);
    pb <- NULL
  }
  assign(".igraph.pb", pb, env=asNamespace("igraph"))
  0L
}

.igraph.progress.tkconsole <- function(percent, message) {
  pb <- get(".igraph.pb", asNamespace("igraph"))
  startmess <- FALSE

  ## Open the console, if it is not open
  if (is.null(pb)) {
    startmess <- TRUE
    pb <- .igraph.progress.tkconsole.create()
  }

  ## Update progress bar
  pb$pb$set(pb$pb$widget, percent)
  tkconfigure(pb$pb$label, text=substr(message, 1, 20))
  tcl("update", "idletasks")
  
  ## Done
  assign(".igraph.pb", pb, env=asNamespace("igraph"))  
  if (startmess) .igraph.progress.tkconsole.message("Console started.\n")
  0L
}

.igraph.progress.tkconsole.create <- function() {
  console <- tktoplevel()
  tktitle(console) <- "igraph console"

  fn <- tkfont.create(family="courier", size=8)

  lfr <- tkframe(console)
  image <- tkimage.create("photo", "img", format="gif",
                          file=system.file("igraph2.gif", package="igraph"))
  logo <- tklabel(lfr, relief="flat", padx=10, pady=10, image=image)
  
  scr <- tkscrollbar(console, repeatinterval=5,
                     command=function(...) tkyview(txt, ...)) 
  txt <- tktext(console, yscrollcommand=function(...) tkset(scr, ...),
                width=60, height=7, font=fn)
  tkconfigure(txt, state="disabled")
  pbar <- .igraph.progress.tkconsole.pbar(console)

  bclear <- tkbutton(lfr, text="Clear", command=function() {
    tkconfigure(txt, state="normal")
    tkdelete(txt, "0.0", "end")
    tkconfigure(txt, state="disabled")
  })
  bstop  <- tkbutton(lfr, text="Stop",  command=function() {})
  bclose <- tkbutton(lfr, text="Close", command=function() {
    tkdestroy(console) })

  tkpack(logo, side="top", fill="none", expand=0, anchor="n",
         ipadx=10, ipady=10)
  tkpack(bclear, side="top", fill="x", expand=0, padx=10)
  ## tkpack(bstop, side="top", fill="x", expand=0, padx=10)
  tkpack(bclose, side="top", fill="x", expand=0, padx=10)  
  
  tkpack(lfr, side="left", fill="none", expand=0, anchor="n")
  tkpack(pbar$frame, side="bottom", fill="x", expand=0)
  tkpack(scr, side="right", fill="y", expand=0)
  tkpack(txt, side="left", fill="both", expand=1)

  tkbind(console, "<Destroy>", function() {
    assign(".igraph.pb", NULL, env=asNamespace("igraph"))
  })
  
  res <- list(top=console, txt=txt, pb=pbar$pb)
  class(res) <- "igraphconsole"
  res
}

.igraph.progress.tkconsole.message <- function(message, start=FALSE) {
  txt <- get(".igraph.pb", asNamespace("igraph"))$txt
  if (is.null(txt)) {
    if (start) {
      pb <- .igraph.progress.tkconsole.create()
      assign(".igraph.pb", pb, env=asNamespace("igraph"))
      txt <- pb$txt
    } else { 
      return()
    }
  }
  tkconfigure(txt, state="normal")
  now <- paste(sep="", substr(date(), 5, 19), ": ")
  s1 <- grepl("^ ", message)
  if (!s1) { tkinsert(txt, "insert", now) }
  tkinsert(txt, "insert", message)
  tksee(txt, "end")
  tkconfigure(txt, state="disabled")
  tcl("update", "idletasks")
}

close.igraphconsole <- function(con, ...) {
  invisible()
}

## Much of this is from tkProgressbar

.igraph.progress.tkconsole.pbar <- function(top) {
  useText <- FALSE
  have_ttk <- as.character(tcl("info", "tclversion")) >= "8.5"
  if (!have_ttk && as.character(tclRequire("PBar")) == "FALSE") 
    useText <- TRUE
  fn <- tkfont.create(family = "helvetica", size = 10)
  frame <- tkframe(top)
  if (useText) {
    .lab <- tklabel(frame, text = " ", font = fn, anchor="w",
                    padx = 20)
    tkpack(.lab, side = "left", anchor="w", padx=5)
    fn2 <- tkfont.create(family = "helvetica", size = 12)
    .vlab <- tklabel(frame, text = "0%", font = fn2, padx = 20)
    tkpack(.vlab, side = "right")
  } else {
    .lab <- tklabel(frame, text = " ", font = fn, anchor="w",
                    pady = 5)
    tkpack(.lab, side = "top", anchor="w", padx=5)
    tkpack(tklabel(frame, text = "", font = fn), side = "bottom")
    .val <- tclVar()
    pBar <- if (have_ttk) {
      ttkprogressbar(frame, length = 300, variable=.val)
    } else {
      tkwidget(frame, "ProgressBar", width = 300, variable=.val)
    }
    tkpack(pBar, side = "bottom", anchor="w", padx=5)
  }
  get <- function(w) { return(tclvalue(.val)); }
  set <- function(w, val) { tclvalue(.val) <<- val }
  pb <- list(widget=pBar, get=get, set=set, label=.lab)
  list(frame=frame, pb=pb)
}
