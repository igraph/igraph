
#   SimpleGraph R package
#   Copyright (C) 2003, 2004  Gabor Csardi <csardi@rmki.kfki.hu>
#   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
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
#   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
###################################################################

###################################################################
# Internal variables
###################################################################

# the environment containing all the plots
if (!exists(".tkplot.env")) {
  .tkplot.env <- new.env()
  assign(".next", 1, .tkplot.env)
}

###################################################################
# Main function
###################################################################

tkplot <- function(graph, layout=layout.random, layout.par=list(),
                   labels=NULL, label.color="darkblue",
                   label.font=NULL, label.degree=-pi/4, label.dist=0,
                   vertex.color="SkyBlue2", vertex.size=15,
                   edge.color="darkgrey", edge.width=1) {

  # Libraries
  require(tcltk) || stop("tcl/tk library not available")

  # Layout
  if (is.function(layout)) {
    layout <- layout(graph, layout.par)
  } else if (is.character(layout) && length(layout)==1 &&
             substr(layout, 1, 2)=="a:") {
    layout <- matrix(unlist(get.vertex.attribute(graph, substring(layout,3))),
                     nr=vcount(graph), byrow=TRUE)[,1:2]
  }

  # Vertex color
  if (length(vertex.color)==1 && substr(vertex.color, 1, 2)=="a:") {
    vertex.color <- as.character(get.vertex.attribute
                                 (graph, substring(vertex.color,3)))
  }

  # Vertex size
  if (is.character(vertex.size) &&
      length(vertex.size)==1 && substr(vertex.size, 1, 2)=="a:") {
    vertex.size <- as.numeric(get.vertex.attribute
                                 (graph, substring(vertex.size,3)))
  }

  # Edge color
  if (length(edge.color)==1 && substr(edge.color, 1, 2)=="a:") {
    edge.color <- as.character(get.edge.attribute
                               (graph, substring(edge.color,3)))
  }

  # Edge width
  if (is.character(edge.width) &&
      length(edge.width)==1 && substr(edge.width, 1, 2)=="a:") {
    edge.width <- as.character(get.edge.attribute
                               (graph, substring(edge.width,3)))
  }

  # Label degree
  if (is.character(label.degree) &&
      length(label.degree)==1 && substr(label.degree, 1, 2)=="a:") {
    label.degree <- as.numeric(get.vertex.attribute
                               (graph, substring(label.degree,3)))
  }
  
  # Label font
  if (is.null(label.font)) {
    label.font=tkfont.create(family="helvetica", size=16, weight="bold")
  }
  
  # Create window & canvas
  top <- tktoplevel(background="lightgrey")
  canvas <- tkcanvas(top, relief="raised", width=450, height=450,
                     borderwidth=2)
  tkpack(canvas, fill="both", expand=1)

  # Create parameters
  params <- list(vertex.color=vertex.color, vertex.size=vertex.size,
                 edge.color=edge.color, label.color=label.color,
                 labels.state=1, edge.width=edge.width, padding=30,
                 grid=0, label.font=label.font, label.degree=label.degree,
                 label.dist=label.dist)

  # The popup menu
  popup.menu <- tkmenu(canvas)
  tkadd(popup.menu, "command", label="Fit to screen", command=function() {
    tkplot.fit.to.screen(tkp.id)})  
  
  # Create plot object
  tkp <- list(top=top, canvas=canvas, graph=graph, coords=layout,
              labels=labels, params=params, popup.menu=popup.menu)
  tkp.id <- .tkplot.new(tkp)
  tktitle(top) <- paste("Graph plot", as.character(tkp.id))

  # The main pull-down menu
  main.menu <- tkmenu(top)
  tkadd(main.menu, "command", label="Quit", command=function() {
    tkplot.close(tkp.id, TRUE)})
  layout.menu <- .tkplot.layout.menu(tkp.id, main.menu)
  tkadd(main.menu, "cascade", label="Layout", menu=layout.menu)
  view.menu <- tkmenu(main.menu)
  tkadd(main.menu, "cascade", label="View", menu=view.menu)
  tkadd(view.menu, "command", label="Fit to screen", command=function() {
    tkplot.fit.to.screen(tkp.id)})
  tkadd(view.menu, "command", label="Center on screen", command=function() {
    tkplot.center(tkp.id)})
  tkadd(view.menu, "separator")
  view.menu.labels <- tclVar(1)
  view.menu.grid <- tclVar(0)
  tkadd(view.menu, "checkbutton", label="Labels",
        variable=view.menu.labels, command=function() {
          .tkplot.toggle.labels(tkp.id)})
# grid canvas object not implemented in tcltk (?) :(
#   tkadd(view.menu, "checkbutton", label="Grid",
#         variable=view.menu.grid, command=function() {
#           .tkplot.toggle.grid(tkp.id)})
  tkadd(view.menu, "separator")
  rotate.menu <- tkmenu(view.menu)
  tkadd(view.menu, "cascade", label="Rotate", menu=rotate.menu)
  sapply(c(-90,-45,-15,-5,-1,1,5,15,45,90),
         function(deg) {
           tkadd(rotate.menu, "command",
                 label=paste(deg, "degree"), command=function() {
                   tkplot.rotate(tkp.id, deg=deg)
                 })
         })
  export.menu <- tkmenu(main.menu)
  tkadd(main.menu, "cascade", label="Export", menu=export.menu)
  tkadd(export.menu, "command", label="Postscript", command=function() {
    tkplot.export.postscript(tkp.id)})
  tkconfigure(top, "-menu", main.menu)
  
  # plot it
  .tkplot.create.edges(tkp.id)
  .tkplot.create.vertices(tkp.id)
  # we would need an update here
  tkplot.fit.to.screen(tkp.id, 450, 450)

  # Kill myself if window was closed
  tkbind(top, "<Destroy>", function() tkplot.close(tkp.id, FALSE))

###################################################################
# The callbacks for interactive editing
###################################################################  
  
  tkitembind(canvas, "vertex||label", "<1>", function(x, y) {
    tkp <- .tkplot.get(tkp.id)
    canvas <- .tkplot.get(tkp.id, "canvas")
    tkdtag(canvas, "selected")
    tkaddtag(canvas, "selected", "withtag", "current")
                                        # get the id
    tags <- as.character(tkgettags(tkp$canvas, "selected"))
    id <- as.numeric(strsplit(tags[pmatch("v-", tags)],
                              "-", fixed=TRUE)[[1]][2])
    if (is.na(id)) { return() }
    tkitemraise(canvas, paste(sep="", "v-", id))
  })
  tkitembind(canvas, "vertex||label", "<ButtonRelease-1>", function(x, y) {
    canvas <- .tkplot.get(tkp.id, "canvas")
    tkdtag(canvas, "selected")
  })
  tkitembind(canvas, "vertex||edge||label", "<Shift-Alt-1>", function(x, y) {
    canvas <- .tkplot.get(tkp.id, "canvas")
    tkitemlower(canvas, "current")
  })
  tkitembind(canvas, "vertex||edge||label", "<Alt-1>", function(x, y) {
    canvas <- .tkplot.get(tkp.id, "canvas")
    tkitemraise(canvas, "current")
  })
  tkbind(top, "<3>", function(x, y) {
    menu <- .tkplot.get(tkp.id, "popup.menu")
    canvas <- .tkplot.get(tkp.id, "canvas")
    x <- as.integer(x) + as.integer(tkwinfo("rootx", canvas))
    y <- as.integer(y) + as.integer(tkwinfo("rooty", canvas))
    .Tcl(paste("tk_popup", .Tcl.args(menu, x, y)))
  })
  tkitembind(canvas, "vertex", "<B1-Motion>", function(x, y) {
    tkp <- .tkplot.get(tkp.id)
    x <- as.numeric(x)
    y <- as.numeric(y)
    width <- as.numeric(tkwinfo("width", tkp$canvas))
    height <- as.numeric(tkwinfo("height", tkp$canvas))
    if (x < 10) { x <- 10 }
    if (x > width-10) { x <- width-10 }
    if (y < 10) { y <- 10 }
    if (y > height-10) { y <- height-10 }
    
                                        # get the id
    tags <- as.character(tkgettags(tkp$canvas, "selected"))
    id <- as.numeric(strsplit(tags[pmatch("v-", tags)],
                              "-", fixed=TRUE)[[1]][2])
    if (is.na(id)) { return() }
                                        # move the vertex
    .tkplot.set.vertex.coords(tkp.id, id, x, y)
    .tkplot.update.vertex(tkp.id, id, x, y)
  })  
  tkitembind(canvas, "label", "<B1-Motion>", function(x,y) {
    tkp <- .tkplot.get(tkp.id)
    x <- as.numeric(x)
    y <- as.numeric(y)
                                        # get the id
    tags <- as.character(tkgettags(tkp$canvas, "selected"))
    id <- as.numeric(strsplit(tags[pmatch("v-", tags)],
                              "-", fixed=TRUE)[[1]][2])
    if (is.na(id)) { return() }
    phi <- pi+atan2(tkp$coords[id,2]-y, tkp$coords[id,1]-x)
    .tkplot.set.label.degree(tkp.id, id, phi)
    .tkplot.update.label(tkp.id, id, tkp$coords[id,1], tkp$coords[id,2])
  })
  
  # We don't need these any more, they are stored in the environment
  rm(tkp, params, layout, vertex.color, edge.color, top, canvas,
     main.menu, layout.menu, view.menu, export.menu, label.font, label.degree)
  
  tkp.id
}

###################################################################
# Internal functions handling data about layouts for the GUI
###################################################################

.tkplot.addlayout <- function(name, layout.data) {
  if (!exists(".layouts", envir=.tkplot.env)) {
    assign(".layouts", list(), .tkplot.env)
  }

  assign("tmp", layout.data, .tkplot.env)
  cmd <- paste(sep="", ".layouts[[\"", name, "\"]]", " <- tmp")
  eval(parse(text=cmd), .tkplot.env)
  rm("tmp", envir=.tkplot.env)
}

.tkplot.getlayout <- function(name) {
  cmd <- paste(sep="", ".layouts[[\"", name, "\"]]")
  eval(parse(text=cmd), .tkplot.env)
}

.tkplot.layouts.newdefaults <- function(name, defaults) {
  assign("tmp", defaults, .tkplot.env)
  for (i in seq(along=defaults)) {
    cmd <- paste(sep="", '.layouts[["', name, '"]]$params[[', i,
                 ']]$default <- tmp[[', i, ']]')
    eval(parse(text=cmd), .tkplot.env)
  }
}

.tkplot.getlayoutlist <- function() {
  eval(parse(text="names(.layouts)"), .tkplot.env)
}

.tkplot.getlayoutname <- function(name) {
  cmd <- paste(sep="", '.layouts[["', name, '"]]$name')
  eval(parse(text=cmd), .tkplot.env)
}

.tkplot.addlayout("random",
                  list(name="Random", f=layout.random, params=list()))
.tkplot.addlayout("circle",
                  list(name="Circle", f=layout.circle, params=list()))
.tkplot.addlayout("fruchterman.reingold",
                  list(name="Fruchterman-Reingold",
                       f=layout.fruchterman.reingold,
                       params=list(
                         niter=list(name="Number of iterations",
                           type="numeric",
                           default=100),
                         coolexp=list(name="Cooling exponent",
                           type="numeric",
                           default=1.5),
                         frame=list(name="Frame shape",
                           type="choice",
                           default="rectangle",
                           values=c("rectangle", "circle")),
                         initial=list(name="Use current as initial layout",
                           type="initial",
                           default=TRUE),
                         temp=list(name="Initial temperature",
                           type="numeric",
                           default=1/10)
                         )
                       )
                  )
.tkplot.addlayout("kamada.kawai",
                  list(name="Kamada-Kawai",
                       f=layout.kamada.kawai,
                       params=list(
                         niter=list(name="Number of iterations",
                           type="numeric",
                           default=1000),
                         initemp=list(name="Initial temperature",
                           type="numeric",
                           default=10),
                         coolexp=list(name="Cooling exponent",
                           type="numeric",
                           default=0.99)
                         )
                       )
                  )
.tkplot.addlayout("spring",
                  list(name="Spring Embedder",
                       f=layout.spring,
                       params=list(
                         mass=list(names="The vertex mass",
                           type="numeric",
                           default=0.1),
                         equil=list(names="The equilibrium spring extension",
                           type="numeric",
                           default=1),
                         k=list(names="The spring coefficient",
                           type="numeric",
                           default=0.001),
                         repeqdis=list(names="Repulsion balance point",
                           type="numeric",
                           default=0.1),
                         kfr=list(names="Friction base coefficient",
                           type="numeric",
                           default=0.01),
                         repulse=list(names="Use repulsion",
                           type="logical",
                           default=FALSE)
                         )
                       )
                  )
                         
                       
###################################################################
# Other public functions, misc.
###################################################################

tkplot.close <- function(tkp.id, window.close=TRUE) {
  if (window.close) {
    cmd <- paste(sep="", "tkp.", tkp.id, "$top")
    top <- eval(parse(text=cmd), .tkplot.env)
    tkbind(top, "<Destroy>", "")
    tkdestroy(top)
  }
  cmd <- paste(sep="", "tkp.", tkp.id)
  rm(list=cmd, envir=.tkplot.env)
}

tkplot.off <- function() {
  eapply(.tkplot.env, function(tkp) { tkdestroy(tkp$top) })
  rm(list=ls(.tkplot.env), envir=.tkplot.env)
  TRUE
}

tkplot.fit.to.screen <- function(tkp.id, width=NULL, height=NULL) {
  tkp <- .tkplot.get(tkp.id)
  if (is.null(width)) {
    width  <- as.numeric(tkwinfo("width", tkp$canvas))
  }
  if (is.null(height)) {
    height <- as.numeric(tkwinfo("height", tkp$canvas))
  }
  coords <- .tkplot.get(tkp.id, "coords")
  # Shift to zero
  coords[,1] <- coords[,1]-min(coords[,1])
  coords[,2] <- coords[,2]-min(coords[,2])
  # Scale
  coords[,1] <- coords[,1] / max(coords[,1]) * (width-(tkp$params$padding*2))
  coords[,2] <- coords[,2] / max(coords[,2]) * (height-(tkp$params$padding*2))
  # Padding
  coords <- coords+tkp$params$padding
  # Store
  .tkplot.set(tkp.id, "coords", coords)
  # Update
  .tkplot.update.vertices(tkp.id)
  invisible(TRUE)
}

tkplot.center <- function(tkp.id) {
  tkp <- .tkplot.get(tkp.id)
  width  <- as.numeric(tkwinfo("width", tkp$canvas))
  height <- as.numeric(tkwinfo("height", tkp$canvas))
  coords <- .tkplot.get(tkp.id, "coords")
  canvas.center.x <- width/2
  canvas.center.y <- height/2
  coords <- .tkplot.get(tkp.id, "coords")
  r1 <- range(coords[,1])
  r2 <- range(coords[,2])
  coords.center.x <- (r1[1]+r1[2])/2
  coords.center.y <- (r2[1]+r2[2])/2
  # Shift to center
  coords[,1] <- coords[,1]+canvas.center.x-coords.center.x
  coords[,2] <- coords[,2]+canvas.center.y-coords.center.y
  # Store
  .tkplot.set(tkp.id, "coords", coords)
  # Update
  .tkplot.update.vertices(tkp.id)
  invisible(TRUE)
}
  
tkplot.reshape <- function(tkp.id, newlayout, ...) {
  tkp <- .tkplot.get(tkp.id)
  .tkplot.set(tkp.id, "coords", newlayout(tkp$graph, ...))
  tkplot.fit.to.screen(tkp.id)
  .tkplot.update.vertices(tkp.id)
}

tkplot.export.postscript <- function(tkp.id) {

  tkp <- .tkplot.get(tkp.id)

  filename <- tkgetSaveFile(initialfile="Rplots.eps",
                            defaultextension="eps",
                            title="Export graph to PostScript file")
  tkpostscript(tkp$canvas, file=filename)
}

tkplot.getcoords <- function(tkp.id, norm=FALSE) {
  coords <- .tkplot.get(tkp.id, "coords")
  if (norm) {
    # Shift
    coords[,1] <- coords[,1]-min(coords[,1])
    coords[,2] <- coords[,2]-min(coords[,2])
    # Scale
    coords[,1] <- coords[,1] / max(coords[,1])-0.5
    coords[,2] <- coords[,2] / max(coords[,2])-0.5
  }
  coords
}

tkplot.rotate <- function(tkp.id, degree=NULL, rad=NULL) {
  coords <- .tkplot.get(tkp.id, "coords")

  if (is.null(degree) && is.null(rad)) {
    rad <- pi/2
  } else if (is.null(rad) && !is.null(degree)) {
    rad <- degree/180*pi
  }
  
  center <- c(mean(range(coords[,1])), mean(range(coords[,2])))
  phi <- atan2(coords[,2]-center[2], coords[,1]-center[1])
  r   <- sqrt((coords[,1]-center[1])**2 + (coords[,2]-center[2])**2)

  phi <- phi + rad

 coords[,1] <- r * cos(phi)
  coords[,2] <- r * sin(phi)
  
  .tkplot.set(tkp.id, "coords", coords)
  tkplot.center(tkp.id)
  invisible(TRUE)
}

###################################################################
# Internal functions, handling the internal environment
###################################################################

.tkplot.new <- function(tkp) {
  id <- get(".next", .tkplot.env)
  assign(".next", id+1, .tkplot.env)
  assign("tmp", tkp, .tkplot.env)
  cmd <- paste("tkp.", id, "<- tmp", sep="")
  eval(parse(text=cmd), .tkplot.env)
  rm("tmp", envir=.tkplot.env)
  id
}

.tkplot.get <- function(tkp.id, what=NULL) {
  if (is.null(what)) {
    get(paste("tkp.", tkp.id, sep=""), .tkplot.env)
  } else {
    cmd <- paste("tkp.", tkp.id, "$", what, sep="")
    eval(parse(text=cmd), .tkplot.env)
  }
}

.tkplot.set <- function(tkp.id, what, value) {
  assign("tmp", value, .tkplot.env)
  cmd <- paste(sep="", "tkp.", tkp.id, "$", what, "<-tmp")
  eval(parse(text=cmd), .tkplot.env)
  rm("tmp", envir=.tkplot.env)
  TRUE
}

.tkplot.set.params <- function(tkp.id, what, value) {
  assign("tmp", value, .tkplot.env)
  cmd <- paste(sep="", "tkp.", tkp.id, "$params$", what, "<-tmp")
  eval(parse(text=cmd), .tkplot.env)
  rm("tmp", envir=.tkplot.env)
  TRUE
}

.tkplot.set.vertex.coords <- function(tkp.id, id, x, y) {
  cmd <- paste(sep="", "tkp.", tkp.id, "$coords[",id,",]<-c(",x,",",y,")")
  eval(parse(text=cmd), .tkplot.env)
  TRUE
}

.tkplot.set.label.degree <- function(tkp.id, id, phi) {
  tkp <- .tkplot.get(tkp.id)
  
  if (length(tkp$params$label.degree)==1) {
    label.degree <- rep(tkp$params$label.degree, times=vcount(tkp$graph))
    label.degree[id] <- phi
    assign("tmp", label.degree, .tkplot.env)
    cmd <- paste(sep="", "tkp.", tkp.id, "$params$label.degree <- tmp")
    eval(parse(text=cmd), .tkplot.env)
    rm("tmp", envir=.tkplot.env)
  } else {
    cmd <- paste(sep="", "tkp.", tkp.id, "$params$label.degree[", id,
                 "] <- ", phi)
    eval(parse(text=cmd), .tkplot.env)
  }
  TRUE
}   

###################################################################
# Internal functions, creating and updating canvas objects
###################################################################

# Creates a new vertex tk object
.tkplot.create.vertex <- function(tkp.id, id, label, x=0, y=0) {
  tkp <- .tkplot.get(tkp.id)
  vertex.size <- ifelse(length(tkp$params$vertex.size)>1,
                        tkp$params$vertex.size[id],
                        tkp$params$vertex.size)
  vertex.color <- ifelse(length(tkp$params$vertex.color)>1,
                         tkp$params$vertex.color[id],
                         tkp$params$vertex.color)
  label.degree <- ifelse(length(tkp$params$label.degree)>1,
                         tkp$params$label.degree[id],
                         tkp$params$label.degree)
  label.dist <- tkp$params$label.dist
  label.x <- x+label.dist*cos(label.degree)*
    (vertex.size+6+4*(ceiling(log10(id+1))))
  label.y <- y+label.dist*sin(label.degree)*
    (vertex.size+6+4*(ceiling(log10(id+1))))
  item <- tkcreate(tkp$canvas, "oval", x-vertex.size, y-vertex.size,
                   x+vertex.size, y+vertex.size, width=1,
                   outline="black",  fill=vertex.color,
                   activefill="red", activeoutline="yellow", activewidth=2)
  tkaddtag(tkp$canvas, "vertex", "withtag", item)
  tkaddtag(tkp$canvas, paste("v-", id, sep=""), "withtag", item)
  litem <- tkcreate(tkp$canvas, "text", label.x, label.y,
                    text=as.character(label), state="normal",
                    fill=tkp$params$label.color, activefill="red",
                    font=tkp$params$label.font)
  tkaddtag(tkp$canvas, "label", "withtag", litem)
  tkaddtag(tkp$canvas, paste("v-", id, sep=""), "withtag", litem)
  item
}

# Create all vertex objects and move them into correct position
.tkplot.create.vertices <- function(tkp.id) {
  tkp <- .tkplot.get(tkp.id)
  n <- vcount(tkp$graph)

  # Labels
  if (is.null(tkp$labels)) {
    labels <- 1:n
  } else if (is.character(tkp$labels) && length(tkp$labels)==1 &&
             substr(tkp$labels, 1, 2)=="a:") {
    labels <- get.vertex.attribute(tkp$graph, substring(tkp$labels, 3))
  } else {
    labels <- tkp$labels
  }
  
  mapply(function(v, l, x, y) .tkplot.create.vertex(tkp.id, v, l, x, y),
         1:n, labels, tkp$coords[,1], tkp$coords[,2])
}

.tkplot.update.label <- function(tkp.id, id, x, y) {
  tkp <- .tkplot.get(tkp.id)
  vertex.size <- ifelse(length(tkp$params$vertex.size)>1,
                        tkp$params$vertex.size[id],
                        tkp$params$vertex.size)
  label.degree <- ifelse(length(tkp$params$label.degree)>1,
                         tkp$params$label.degree[id],
                         tkp$params$label.degree)
  label.dist <- tkp$params$label.dist
  label.x <- x+label.dist*cos(label.degree)*
    (vertex.size+6+4*(ceiling(log10(id+1))))
  label.y <- y+label.dist*sin(label.degree)*
    (vertex.size+6+4*(ceiling(log10(id+1))))
  tkcoords(tkp$canvas, paste("label&&v-", id, sep=""),
           label.x, label.y)
}

.tkplot.update.vertex <- function(tkp.id, id, x, y) {
  tkp <- .tkplot.get(tkp.id)
  vertex.size <- ifelse(length(tkp$params$vertex.size)>1,
                        tkp$params$vertex.size[id],
                        tkp$params$vertex.size)
  # Vertex
  tkcoords(tkp$canvas, paste("vertex&&v-", id, sep=""),
           x-vertex.size, y-vertex.size,
           x+vertex.size, y+vertex.size)
  # Label
  .tkplot.update.label(tkp.id, id, x, y)
  
  # Edges
  edge.from.ids <- as.numeric(tkfind(tkp$canvas, "withtag",
                                     paste("from-", id, sep="")))
  edge.to.ids <- as.numeric(tkfind(tkp$canvas, "withtag",
                                   paste("to-", id, sep="")))
  for (i in seq(along=edge.from.ids)) {
    .tkplot.update.edge(tkp.id, edge.from.ids[i])
  }
  for (i in seq(along=edge.to.ids)) {
    .tkplot.update.edge(tkp.id, edge.to.ids[i])
  }
}

.tkplot.update.vertices <- function(tkp.id) {
  tkp <- .tkplot.get(tkp.id)
  n <- vcount(tkp$graph)
  mapply(function(v, x, y) .tkplot.update.vertex(tkp.id, v, x, y), 1:n,
         tkp$coords[,1], tkp$coords[,2])
}

# Creates tk object for edge 'id'
.tkplot.create.edge <- function(tkp.id, from, to, id) {
  tkp <- .tkplot.get(tkp.id)  
  arrow <- ifelse(is.directed(tkp$graph), "last", "none")
  from.c <- tkp$coords[from,]
  to.c   <- tkp$coords[to,]
  edge.color <- ifelse(length(tkp$params$edge.color)>1,
                       tkp$params$edge.color[id],
                       tkp$params$edge.color)
  edge.width <- ifelse(length(tkp$params$edge.width)>1,
                       tkp$params$edge.width[id],
                       tkp$params$edge.width)
  tkcreate(tkp$canvas, "line", from.c[1], from.c[2], to.c[1], to.c[2],
           width=edge.width, activewidth=2*edge.width,
           arrow=arrow, arrowshape=c(10, 10, 5),
           fill=edge.color, activefill="red",           
           tags=c("edge", paste(sep="", "edge-", id),
             paste(sep="", "from-", from),
             paste(sep="", "to-", to)))
}

# Creates all edges
.tkplot.create.edges <- function(tkp.id) {
  tkp <- .tkplot.get(tkp.id)
  n <- ecount(tkp$graph)
  edgematrix <- get.edgelist(tkp$graph)
  mapply(function(from, to, id) .tkplot.create.edge(tkp.id, from, to, id),
         edgematrix[,1],
         edgematrix[,2], 1:nrow(edgematrix))
}

# Update an edge with given itemid (not edge id!)
.tkplot.update.edge <- function(tkp.id, itemid) {
  tkp <- .tkplot.get(tkp.id)
  tags <- as.character(tkgettags(tkp$canvas, itemid))
  from <- as.numeric(substring(grep("from-", tags, value=TRUE, fixed=TRUE),6))
  to <- as.numeric(substring(grep("to-", tags, value=TRUE, fixed=TRUE),4))
  from.c <- tkp$coords[from,]
  to.c <- tkp$coords[to,]
  if (is.directed(tkp$graph)) {
    phi <- atan2(to.c[2]-from.c[2], to.c[1]-from.c[1])
    r <- sqrt( (to.c[1]-from.c[1])^2 + (to.c[2]-from.c[2])^2 )
    vertex.size <- ifelse(length(tkp$params$vertex.size)>1,
                          tkp$params$vertex.size[to],
                          tkp$params$vertex.size)
    to.c[1] <- from.c[1] + (r-vertex.size)*cos(phi)
    to.c[2] <- from.c[2] + (r-vertex.size)*sin(phi)
  }
  tkcoords(tkp$canvas, itemid, from.c[1], from.c[2], to.c[1], to.c[2])
}

.tkplot.toggle.labels <- function(tkp.id) {
  .tkplot.set.params(tkp.id, "labels.state",
                    1 - .tkplot.get(tkp.id, "params")$labels.state)
  tkp <- .tkplot.get(tkp.id)
  state <- ifelse(tkp$params$labels.state==1, "normal", "hidden")
  tkitemconfigure(tkp$canvas, "label", "-state", state)  
}

.tkplot.toggle.grid <- function(tkp.id) {
  .tkplot.set.params(tkp.id, "grid",
                    1 - .tkplot.get(tkp.id, "params")$grid)
  tkp <- .tkplot.get(tkp.id)
  state <- ifelse(tkp$params$grid==1, "normal", "hidden")
  if (state=="hidden") {
    tkdelete(tkp$canvas, "grid")
  } else {
    tkcreate(tkp$canvas, "grid", 0, 0, 10, 10, tags=c("grid"))
  }
}

###################################################################
# Internal functions: manipulating the UI
###################################################################

.tkplot.layout.menu <- function(tkp.id, main.menu) {
  layout.menu <- tkmenu(main.menu)
  
  sapply(.tkplot.getlayoutlist(), function(n) {
    tkadd(layout.menu, "command", label=.tkplot.getlayoutname(n),
          command=function() {
            .tkplot.layout.dialog(tkp.id, n)
          })
  })
  
  layout.menu
}

.tkplot.layout.dialog <- function(tkp.id, layout.name) {
  layout <- .tkplot.getlayout(layout.name)

  # No parameters
  if (length(layout$params)==0) {
    return(tkplot.reshape(tkp.id, layout$f, params=list()))
  }
  
  submit <- function() {
    realparams <- params <- vector(mode="list", length(layout$params))
    names(realparams) <- names(params) <- names(layout$params)
    for (i in seq(along=layout$params)) {
      realparams[[i]] <-
        params[[i]] <- switch(layout$params[[i]]$type,
                              "numeric"=as.numeric(tkget(values[[i]])),
                              "character"=as.character(tkget(values[[i]])),
                              "logical"=as.logical(tclvalue(values[[i]])),
                              "choice"=as.character(tclvalue(values[[i]])),
                              "initial"=as.logical(tclvalue(values[[i]]))
                              )
      if (layout$params[[i]]$type=="initial" &&
          params[[i]]) {
        realparams[[i]] <- tkplot.getcoords(tkp.id, norm=TRUE)
      }
    }
    if (as.logical(tclvalue(save.default))) {
      .tkplot.layouts.newdefaults(layout.name, params)
    }
    tkdestroy(dialog)
    tkplot.reshape(tkp.id, layout$f, params=realparams)
  }
  
  dialog <- tktoplevel(.tkplot.get(tkp.id, "top"), padx=10, pady=10)
  
  tkwm.title(dialog, paste("Layout parameters for graph plot", tkp.id))
  tkwm.transient(dialog, .tkplot.get(tkp.id, "top"))

  tkgrid(tklabel(dialog, text=paste(layout$name, "layout"),
                 font=tkfont.create(family="helvetica",size=20,weight="bold")),
         row=0, column=0, columnspan=2, padx=10, pady=10)
                 
  row <- 1
  values <- list()
  for (i in seq(along=layout$params)) {
    
    tkgrid(tklabel(dialog, text=paste(sep="", layout$params[[i]]$name, ":")),
                   row=row, column=0, sticky="ne", padx=5, pady=5)
    
    if (layout$params[[i]]$type %in% c("numeric", "character")) {
      values[[i]] <- tkentry(dialog)
      tkinsert(values[[i]], 0, as.character(layout$params[[i]]$default))
      tkgrid(values[[i]], row=row, column=1, sticky="nw", padx=5, pady=5)
    } else if (layout$params[[i]]$type=="logical") {
      values[[i]] <- tclVar(as.character(layout$params[[i]]$default))
      tmp <- tkcheckbutton(dialog, onvalue="TRUE", offvalue="FALSE",
                           variable=values[[i]])
      tkgrid(tmp, row=row, column=1, sticky="nw", padx=5, pady=5)
    } else if (layout$params[[i]]$type=="choice") {
      tmp.frame <- tkframe(dialog)
      tkgrid(tmp.frame, row=row, column=1, sticky="nw", padx=5, pady=5)
      values[[i]] <- tclVar(layout$params[[i]]$default)
      for (j in 1:length(layout$params[[i]]$values)) {
        tmp <- tkradiobutton(tmp.frame, variable=values[[i]],
                             value=layout$params[[i]]$values[j],
                             text=layout$params[[i]]$values[j])
        tkpack(tmp, anchor="nw")
      }
    } else if (layout$params[[i]]$type=="initial") {
      values[[i]] <- tclVar(as.character(layout$params[[i]]$default))
      tkgrid(tkcheckbutton(dialog, onvalue="TRUE", offvalue="FALSE",
                           variable=values[[i]]),
             row=row, column=1, sticky="nw", padx=5, pady=5)
    }

    row <- row + 1
    
  } # for along layout$params

  tkgrid(tklabel(dialog, text="Set these as defaults"), sticky="ne",
         row=row, column=0, padx=5, pady=5)
  save.default <- tclVar("FALSE")
  tkgrid(tkcheckbutton(dialog, onvalue="TRUE", offvalue="FALSE",
                       variable=save.default, text=""), row=row,
         column=1, sticky="nw", padx=5, pady=5)
  row <- row + 1
  
  tkgrid(tkbutton(dialog, text="OK", command=submit), row=row, column=0)
  tkgrid(tkbutton(dialog, text="Cancel",
                  command=function() { tkdestroy(dialog); invisible(TRUE) }),
         row=row, column=1)
}

