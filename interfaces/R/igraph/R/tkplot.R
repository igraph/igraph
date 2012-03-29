
#   IGraph R package
#   Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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
#   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
#   02110-1301 USA
#
###################################################################

###################################################################
# Internal variables
###################################################################

# the environment containing all the plots
.tkplot.env <- new.env()
assign(".next", 1, .tkplot.env)

###################################################################
# Main function
###################################################################

tkplot <- function(graph, canvas.width=450, canvas.height=450, ...) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  
  # Libraries
  require(tcltk) || stop("tcl/tk library not available")

  # Visual parameters
  params <- i.parse.plot.params(graph, list(...))
  labels <- params("vertex", "label")
  label.color <- .tkplot.convert.color(params("vertex", "label.color"))
  label.font <- .tkplot.convert.font(params("vertex", "label.font"),
                                     params("vertex", "label.family"),
                                     params("vertex", "label.cex"))
  label.degree <- params("vertex", "label.degree")
  label.dist <- params("vertex", "label.dist")
  vertex.color <- .tkplot.convert.color(params("vertex", "color"))
  vertex.size <- params("vertex", "size")
  vertex.frame.color <- .tkplot.convert.color(params("vertex", "frame.color"))

  edge.color <- .tkplot.convert.color(params("edge", "color"))
  edge.width <- params("edge", "width")
  edge.labels <- params("edge", "label")
  edge.lty <- params("edge", "lty")
  loop.angle <- params("edge", "loop.angle")
  arrow.mode <- params("edge", "arrow.mode")
  edge.label.font <- .tkplot.convert.font(params("edge", "label.font"),
                                          params("edge", "label.family"),
                                          params("edge", "label.cex"))
  edge.label.color <- params("edge", "label.color")
  arrow.size  <- params("edge", "arrow.size")[1]
  curved <- params("edge", "curved")
  curved <- rep(curved, length=ecount(graph))
  
  layout <- params("plot", "layout")
  layout[,2] <- -layout[,2]
  margin <- params("plot", "margin")
  margin <- rep(margin, length=4)

  # the new style parameters can't do this yet
  arrow.mode         <- i.get.arrow.mode(graph, arrow.mode)
  
  # Edge line type
  edge.lty <- i.tkplot.get.edge.lty(edge.lty)

  # Create window & canvas
  top <- tktoplevel(background="lightgrey")
  canvas <- tkcanvas(top, relief="raised",
                     width=canvas.width, height=canvas.height,
                     borderwidth=2)
  tkpack(canvas, fill="both", expand=1)

  # Create parameters
  params <- list(vertex.color=vertex.color, vertex.size=vertex.size,
                 edge.color=edge.color, label.color=label.color,
                 labels.state=1, edge.width=edge.width,
                 padding=margin*300+max(vertex.size)+5,
                 grid=0, label.font=label.font, label.degree=label.degree,
                 label.dist=label.dist, edge.labels=edge.labels,
                 vertex.frame.color=vertex.frame.color,
                 loop.angle=loop.angle, edge.lty=edge.lty, arrow.mode=arrow.mode,
                 edge.label.font=edge.label.font,
                 edge.label.color=edge.label.color, arrow.size=arrow.size,
                 curved=curved)

  # The popup menu
  popup.menu <- tkmenu(canvas)
  tkadd(popup.menu, "command", label="Fit to screen", command=function() {
    tkplot.fit.to.screen(tkp.id)})  

  # Different popup menu for vertices
  vertex.popup.menu <- tkmenu(canvas)
  tkadd(vertex.popup.menu, "command", label="Vertex color",
        command=function() {
          tkp <- .tkplot.get(tkp.id)
          vids <- .tkplot.get.selected.vertices(tkp.id)
          if (length(vids)==0) return(FALSE)
          
          initialcolor <- ifelse(length(tkp$params$vertex.color)>1,
                                 tkp$params$vertex.color[vids[1]+1],
                                 tkp$params$vertex.color)
          color <- .tkplot.select.color(initialcolor)
          if (color=="") return(FALSE) # Cancel
          
          .tkplot.update.vertex.color(tkp.id, vids, color)
        })
  tkadd(vertex.popup.menu, "command", label="Vertex size",
        command=function() {
          tkp <- .tkplot.get(tkp.id)
          vids <- .tkplot.get.selected.vertices(tkp.id)
          if (length(vids)==0) return(FALSE)

          initialsize <- ifelse(length(tkp$params$vertex.size)>1,
                                tkp$params$vertex.size[vids[1]+1],
                                tkp$params$vertex.size)
          size <- .tkplot.select.number("Vertex size", initialsize, 1, 20)
          if (is.na(size)) return(FALSE)

          .tkplot.update.vertex.size(tkp.id, vids, size)
        })
  
  # Different popup menu for edges
  edge.popup.menu <- tkmenu(canvas)
  tkadd(edge.popup.menu, "command", label="Edge color",
        command=function() {
          tkp <- .tkplot.get(tkp.id)
          eids <- .tkplot.get.selected.edges(tkp.id)
          if (length(eids)==0) return(FALSE)
          
          initialcolor <- ifelse(length(tkp$params$edge.color)>1,
                                 tkp$params$edge.color[eids[1]],
                                 tkp$params$edge.color)
          color <- .tkplot.select.color(initialcolor)
          if (color=="") return(FALSE) # Cancel

          .tkplot.update.edge.color(tkp.id, eids, color)          
        })
  tkadd(edge.popup.menu, "command", label="Edge width",
        command=function() {
          tkp <- .tkplot.get(tkp.id)
          eids <- .tkplot.get.selected.edges(tkp.id)
          if (length(eids)==0) return(FALSE)
          
          initialwidth <- ifelse(length(tkp$params$edge.width)>1,
                                 tkp$params$edge.width[eids[1]],
                                 tkp$params$edge.width)
          width <- .tkplot.select.number("Edge width", initialwidth, 1, 10)
          if (is.na(width)) return(FALSE) # Cancel

          .tkplot.update.edge.width(tkp.id, eids, width)
        })
          
  
  # Create plot object
  tkp <- list(top=top, canvas=canvas, graph=graph, coords=layout,
              labels=labels, params=params, popup.menu=popup.menu,
              vertex.popup.menu=vertex.popup.menu,
              edge.popup.menu=edge.popup.menu)
  tkp.id <- .tkplot.new(tkp)
  tktitle(top) <- paste("Graph plot", as.character(tkp.id))

  # The main pull-down menu
  main.menu <- tkmenu(top)
  tkadd(main.menu, "command", label="Close", command=function() {
    tkplot.close(tkp.id, TRUE)})
  select.menu <- .tkplot.select.menu(tkp.id, main.menu)
  tkadd(main.menu, "cascade", label="Select", menu=select.menu)  
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
                   tkplot.rotate(tkp.id, degree=deg)
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
  tkplot.fit.to.screen(tkp.id, canvas.width, canvas.height)

  # Kill myself if window was closed
  tkbind(top, "<Destroy>", function() tkplot.close(tkp.id, FALSE))

###################################################################
# The callbacks for interactive editing
###################################################################  

  tkitembind(canvas, "vertex||label||edge", "<1>", function(x, y) {
    tkp <- .tkplot.get(tkp.id)
    canvas <- .tkplot.get(tkp.id, "canvas")
    .tkplot.deselect.all(tkp.id)
    .tkplot.select.current(tkp.id)
#     tkitemraise(canvas, "current")
  })
  tkitembind(canvas, "vertex||label||edge", "<Control-1>", function(x,y) {
    canvas <- .tkplot.get(tkp.id, "canvas")
    curtags <- as.character(tkgettags(canvas, "current"))
    seltags <- as.character(tkgettags(canvas, "selected"))
    if ("vertex" %in% curtags && "vertex" %in% seltags) {
      if ("selected" %in% curtags) {
        .tkplot.deselect.current(tkp.id)
      } else {
        .tkplot.select.current(tkp.id)
      }
    } else if ("edge" %in% curtags && "edge" %in% seltags) {
      if ("selected" %in% curtags) {
        .tkplot.deselect.current(tkp.id)
      } else {
        .tkplot.select.current(tkp.id)
      }
    } else if ("label" %in% curtags && "vertex" %in% seltags) {
      vtag <- curtags[pmatch("v-", curtags)]
      tkid <- as.numeric(tkfind(canvas, "withtag",
                                paste(sep="", "vertex&&", vtag)))
      vtags <- as.character(tkgettags(canvas, tkid))
      if ("selected" %in% vtags) {
        .tkplot.deselect.vertex(tkp.id, tkid)
      } else {
        .tkplot.select.vertex(tkp.id, tkid)
      }
    } else {
      .tkplot.deselect.all(tkp.id)
      .tkplot.select.current(tkp.id)
    }
  })
  tkitembind(canvas, "vertex||edge||label", "<Shift-Alt-1>", function(x, y) {
    canvas <- .tkplot.get(tkp.id, "canvas")
    tkitemlower(canvas, "current")
  })
  tkitembind(canvas, "vertex||edge||label", "<Alt-1>", function(x, y) {
    canvas <- .tkplot.get(tkp.id, "canvas")
    tkitemraise(canvas, "current")
  })  
  tkbind(canvas, "<3>", function(x, y) {
    canvas <- .tkplot.get(tkp.id, "canvas")
    tags <- as.character(tkgettags(canvas, "current"))
    if ("label" %in% tags) {
      vtag <- tags[ pmatch("v-", tags) ]
      vid <- as.character(tkfind(canvas, "withtag",
                                 paste(sep="", "vertex&&", vtag)))
      tags <- as.character(tkgettags(canvas, vid))
    }
    if ("selected" %in% tags) {
      # The selection is active
    } else {
      # Delete selection, single object
      .tkplot.deselect.all(tkp.id)
      .tkplot.select.current(tkp.id)
    }
    tags <- as.character(tkgettags(canvas, "selected"))
    ## TODO: what if different types of objects are selected
    if ("vertex" %in% tags || "label" %in% tags) {
      menu <- .tkplot.get(tkp.id, "vertex.popup.menu")
    } else if ("edge" %in% tags) {
      menu <- .tkplot.get(tkp.id, "edge.popup.menu")
    } else {
      menu <- .tkplot.get(tkp.id, "popup.menu")
    }
    x <- as.integer(x) + as.integer(tkwinfo("rootx", canvas))
    y <- as.integer(y) + as.integer(tkwinfo("rooty", canvas))
    .Tcl(paste("tk_popup", .Tcl.args(menu, x, y)))
  })
  if (tkp$params$label.dist==0) tobind <- "vertex||label"
  else tobind <- "vertex"
  tkitembind(canvas, tobind, "<B1-Motion>", function(x, y) {
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
  if (tkp$params$label.dist!=0) {
    tkitembind(canvas, "label", "<B1-Motion>", function(x,y) {
      tkp <- .tkplot.get(tkp.id)
      x <- as.numeric(x)
      y <- as.numeric(y)
                                        # get the id
      tags <- as.character(tkgettags(tkp$canvas, "selected"))
      id <- as.numeric(strsplit(tags[pmatch("v-", tags)],
                                "-", fixed=TRUE)[[1]][2])
      if (is.na(id)) { return() }
      phi <- pi+atan2(tkp$coords[id+1,2]-y, tkp$coords[id+1,1]-x)
      .tkplot.set.label.degree(tkp.id, id, phi)
      .tkplot.update.label(tkp.id, id, tkp$coords[id+1,1], tkp$coords[id+1,2])
    })
  }
  
  # We don't need these any more, they are stored in the environment
  rm(tkp, params, layout, vertex.color, edge.color, top, canvas,
     main.menu, layout.menu, view.menu, export.menu, label.font, label.degree,
     vertex.frame.color)
  
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
                           default=500),
                         maxdelta=list(name="Maximum change (n)",
                           type="expression",
                           default=expression(vcount(.tkplot.g))),
                         area=list(name="Area parameter (n^2)",
                           type="expression",
                           default=expression(vcount(.tkplot.g)^2)),
                         coolexp=list(name="Cooling exponent",
                           type="numeric",
                           default=3),                         
                         repulserad=list(name="Cancellation radius (n^3)",
                           type="expression",
                           # FIXME: this should be area * n, but parameters
                           # can't depend on each other....
                           default=expression(vcount(.tkplot.g)^3))
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
.tkplot.addlayout("reingold.tilford",
                  list(names="Reingold-Tilford",
                       f=layout.reingold.tilford,
                       params=list(
                         root=list(name="Root vertex",
                           type="numeric",
                           default=0)
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
  invisible(NULL)
}

tkplot.off <- function() {
  eapply(.tkplot.env, function(tkp) { tkdestroy(tkp$top) })
  rm(list=ls(.tkplot.env), envir=.tkplot.env)
  invisible(NULL)
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
  coords[,1] <- coords[,1] / max(coords[,1]) *
    (width-(tkp$params$padding[2]+tkp$params$padding[4]))
  coords[,2] <- coords[,2] / max(coords[,2]) *
    (height-(tkp$params$padding[1]+tkp$params$padding[3]))
  # Padding
  coords[,1] <- coords[,1]+tkp$params$padding[2]
  coords[,2] <- coords[,2]+tkp$params$padding[3]
  # Store
  .tkplot.set(tkp.id, "coords", coords)
  # Update
  .tkplot.update.vertices(tkp.id)
  invisible(NULL)
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
  invisible(NULL)
}
  
tkplot.reshape <- function(tkp.id, newlayout, ...) {
  tkp <- .tkplot.get(tkp.id)
  .tkplot.set(tkp.id, "coords", newlayout(tkp$graph, ...))
  tkplot.fit.to.screen(tkp.id)
  .tkplot.update.vertices(tkp.id)
  invisible(NULL)
}

tkplot.export.postscript <- function(tkp.id) {

  tkp <- .tkplot.get(tkp.id)

  filename <- tkgetSaveFile(initialfile="Rplots.eps",
                            defaultextension="eps",
                            title="Export graph to PostScript file")
  tkpostscript(tkp$canvas, file=filename)
  invisible(NULL)
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
  invisible(NULL)
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
  cmd <- paste(sep="", "tkp.", tkp.id, "$coords[",id+1,",]<-c(",x,",",y,")")
  eval(parse(text=cmd), .tkplot.env)
  TRUE
}

.tkplot.set.label.degree <- function(tkp.id, id, phi) {
  tkp <- .tkplot.get(tkp.id)
  
  if (length(tkp$params$label.degree)==1) {
    label.degree <- rep(tkp$params$label.degree, times=vcount(tkp$graph))
    label.degree[id+1] <- phi
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
                        tkp$params$vertex.size[id+1],
                        tkp$params$vertex.size)
  vertex.color <- ifelse(length(tkp$params$vertex.color)>1,
                         tkp$params$vertex.color[id+1],
                         tkp$params$vertex.color)
  vertex.frame.color <- ifelse(length(tkp$params$vertex.frame.color)>1,
                               tkp$params$vertex.frame.color[id+1],
                               tkp$params$vertex.frame.color)
  item <- tkcreate(tkp$canvas, "oval", x-vertex.size, y-vertex.size,
                   x+vertex.size, y+vertex.size, width=1,
                   outline=vertex.frame.color,  fill=vertex.color)
  tkaddtag(tkp$canvas, "vertex", "withtag", item)
  tkaddtag(tkp$canvas, paste("v-", id, sep=""), "withtag", item)
  if (!is.na(label)) {
    label.degree <- ifelse(length(tkp$params$label.degree)>1,
                           tkp$params$label.degree[id+1],
                           tkp$params$label.degree)
    label.color <- if (length(tkp$params$label.color)>1) {
      tkp$params$label.color[id+1]
    } else {
      tkp$params$label.color
    }
    label.dist <- tkp$params$label.dist
    label.x <- x+label.dist*cos(label.degree)*
      (vertex.size+6+4*(ceiling(log10(id+1))))
    label.y <- y+label.dist*sin(label.degree)*
      (vertex.size+6+4*(ceiling(log10(id+1))))
    if (label.dist==0)
      { afill <- label.color }
    else
      { afill <- "red" }
    litem <- tkcreate(tkp$canvas, "text", label.x, label.y,
                      text=as.character(label), state="normal",
                      fill=label.color, activefill=afill,
                      font=tkp$params$label.font)
    tkaddtag(tkp$canvas, "label", "withtag", litem)
    tkaddtag(tkp$canvas, paste("v-", id, sep=""), "withtag", litem)
  }
  item
}

# Create all vertex objects and move them into correct position
.tkplot.create.vertices <- function(tkp.id) {
  tkp <- .tkplot.get(tkp.id)
  n <- vcount(tkp$graph)

  # Labels
  labels <- i.get.labels(tkp$graph, tkp$labels)

  mapply(function(v, l, x, y) .tkplot.create.vertex(tkp.id, v, l, x, y),
         0:(n-1), labels, tkp$coords[,1], tkp$coords[,2])
}

.tkplot.update.label <- function(tkp.id, id, x, y) {
  tkp <- .tkplot.get(tkp.id)
  vertex.size <- ifelse(length(tkp$params$vertex.size)>1,
                        tkp$params$vertex.size[id+1],
                        tkp$params$vertex.size)
  label.degree <- ifelse(length(tkp$params$label.degree)>1,
                         tkp$params$label.degree[id+1],
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
                        tkp$params$vertex.size[id+1],
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
  mapply(function(v, x, y) .tkplot.update.vertex(tkp.id, v, x, y), 0:(n-1),
         tkp$coords[,1], tkp$coords[,2])
}

# Creates tk object for edge 'id'
.tkplot.create.edge <- function(tkp.id, from, to, id) {
  tkp <- .tkplot.get(tkp.id)  
  from.c <- tkp$coords[from+1,]
  to.c   <- tkp$coords[to+1,]
  edge.color <- ifelse(length(tkp$params$edge.color)>1,
                       tkp$params$edge.color[id],
                       tkp$params$edge.color)
  edge.width <- ifelse(length(tkp$params$edge.width)>1,
                       tkp$params$edge.width[id],
                       tkp$params$edge.width)
  edge.lty <- ifelse(length(tkp$params$edge.lty)>1,
                       tkp$params$edge.lty[[id]],
                       tkp$params$edge.lty)
  arrow.mode <- ifelse(length(tkp$params$arrow.mode)>1,
                       tkp$params$arrow.mode[[id]],
                       tkp$params$arrow.mode)
  arrow.size <- tkp$params$arrow.size
  curved <- tkp$params$curved[[id]]
  arrow <- c("none", "first", "last", "both")[arrow.mode+1]
  
  if (from != to) {
    ## non-loop edge
    if (is.logical(curved)) curved <- curved * 0.5
    if (curved != 0) {
      smooth <- TRUE
      midx <- (from.c[1]+to.c[1])/2
      midy <- (from.c[2]+to.c[2])/2        
      spx <- midx - curved * 1/2 * (from.c[2]-to.c[2])
      spy <- midy + curved * 1/2 * (from.c[1]-to.c[1])
      coords <- c(from.c[1], from.c[2], spx, spy, to.c[1], to.c[2])
    } else {
      smooth <- FALSE
      coords <- c(from.c[1], from.c[2], to.c[1], to.c[2])
    }
    args <- c(list(tkp$canvas, "line"),
              coords, 
              list(width=edge.width, activewidth=2*edge.width,
                   arrow=arrow, arrowshape=arrow.size * c(10, 10, 5),
                   fill=edge.color, activefill="red", dash=edge.lty,
                   tags=c("edge", paste(sep="", "edge-", id),
                     paste(sep="", "from-", from),
                     paste(sep="", "to-", to))), smooth=smooth)
    do.call(tkcreate, args)
  } else {
    ## loop edge
    ## the coordinates are not correct but we will call update anyway...
    tkcreate(tkp$canvas, "line", from.c[1], from.c[2],
             from.c[1]+20, from.c[1]-10, from.c[2]+30, from.c[2],
             from.c[1]+20, from.c[1]+10, from.c[1], from.c[2],
             width=edge.width, activewidth=2*edge.width,
             arrow=arrow, arrowshape=arrow.size * c(10,10,5), dash=edge.lty,
             fill=edge.color, activefill="red", smooth=TRUE,
             tags=c("edge", "loop", paste(sep="", "edge-", id),
               paste(sep="", "from-", from),
               paste(sep="", "to-", to)))
    
  }

  edge.label <- ifelse(length(tkp$params$edge.labels)>1,
                       tkp$params$edge.labels[id],
                       tkp$params$edge.labels)
  if (!is.na(edge.label)) {
    label.color <- ifelse(length(tkp$params$edge.label.color)>1,
                          tkp$params$edge.label.color[id],
                          tkp$params$edge.label.color)
    ## not correct for loop edges but we will update anyway...
    label.x <- (to.c[1]+from.c[1])/2
    label.y <- (to.c[2]+from.c[2])/2
    litem <- tkcreate(tkp$canvas, "text", label.x, label.y,
                      text=as.character(edge.label), state="normal",
                      fill=label.color,
                      font=tkp$params$edge.label.font)
    tkaddtag(tkp$canvas, "label", "withtag", litem)
    tkaddtag(tkp$canvas, paste(sep="", "edge-", id), "withtag", litem)
  }
}

# Creates all edges
.tkplot.create.edges <- function(tkp.id) {
  tkp <- .tkplot.get(tkp.id)
  n <- ecount(tkp$graph)
  edgematrix <- get.edgelist(tkp$graph, names=FALSE)
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
  from.c <- tkp$coords[from+1,]
  to.c <- tkp$coords[to+1,]

  edgeid <- as.numeric(substring(tags[ pmatch("edge-", tags) ], 6))

  if (from != to) {
    phi <- atan2(to.c[2]-from.c[2], to.c[1]-from.c[1])
    r <- sqrt( (to.c[1]-from.c[1])^2 + (to.c[2]-from.c[2])^2 )
    vertex.size <- ifelse(length(tkp$params$vertex.size)>1,
                          tkp$params$vertex.size[to+1],
                          tkp$params$vertex.size)
    vertex.size2 <- ifelse(length(tkp$params$vertex.size)>1,
                           tkp$params$vertex.size[from+1],
                           tkp$params$vertex.size)
    curved <- tkp$params$curved[[edgeid]]
    to.c[1] <- from.c[1] + (r-vertex.size)*cos(phi)
    to.c[2] <- from.c[2] + (r-vertex.size)*sin(phi)
    from.c[1] <- from.c[1] + vertex.size2*cos(phi)
    from.c[2] <- from.c[2] + vertex.size2*sin(phi)
    if (is.logical(curved)) curved <- curved * 0.5
    if (curved == 0) {
      tkcoords(tkp$canvas, itemid, from.c[1], from.c[2], to.c[1], to.c[2])
    } else {
      midx <- (from.c[1]+to.c[1])/2
      midy <- (from.c[2]+to.c[2])/2        
      spx <- midx - curved * 1/2 * (from.c[2]-to.c[2])
      spy <- midy + curved * 1/2 * (from.c[1]-to.c[1])
      tkcoords(tkp$canvas, itemid, from.c[1], from.c[2], spx, spy,
               to.c[1], to.c[2])
    }
  } else {
    vertex.size <- ifelse(length(tkp$params$vertex.size)>1,
                          tkp$params$vertex.size[to+1],
                          tkp$params$vertex.size)
    loop.angle <- ifelse(length(tkp$param$loop.angle)>1,
                         tkp$params$loop.angle[edgeid+1],
                         tkp$params$loop.angle)
    xx <- from.c[1] + cos(loop.angle/180*pi)*vertex.size
    yy <- from.c[2] + sin(loop.angle/180*pi)*vertex.size
    cc <- matrix(c(xx,yy, xx+20,yy-10, xx+30,yy, xx+20,yy+10, xx,yy),
                 ncol=2, byrow=TRUE)

    phi <- atan2(cc[,2]-yy, cc[,1]-xx)
    r <- sqrt((cc[,1]-xx)**2 + (cc[,2]-yy)**2)
    phi <- phi+loop.angle/180*pi
    cc[,1] <- xx+r*cos(phi)
    cc[,2] <- yy+r*sin(phi)
    tkcoords(tkp$canvas, itemid, cc[1,1], cc[1,2], cc[2,1], cc[2,2],
             cc[3,1], cc[3,2], cc[4,1], cc[4,2], cc[5,1]+0.001, cc[5,2]+0.001)
  }

  edge.label <- ifelse(length(tkp$params$edge.labels)>1,
                       tkp$params$edge.labels[edgeid],
                       tkp$params$edge.labels)
  if (!is.na(edge.label)) {
    if (from != to) {
      label.x <- (to.c[1]+from.c[1])/2
      label.y <- (to.c[2]+from.c[2])/2
    } else {
      ## loops
      label.x <- xx+cos(loop.angle/180*pi)*30
      label.y <- yy+sin(loop.angle/180*pi)*30
    }
    litem <- as.numeric(tkfind(tkp$canvas, "withtag",
                               paste(sep="", "label&&edge-", edgeid)))
    tkcoords(tkp$canvas, litem, label.x, label.y)
  }
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

.tkplot.update.vertex.color <- function(tkp.id, vids, newcolor) {
  tkp <- .tkplot.get(tkp.id)
  colors <- tkp$params$vertex.color
  if (length(colors)==1 && length(vids)==vcount(tkp$graph)) {
    ## Uniform color -> uniform color
    .tkplot.set(tkp.id, "params$vertex.color", newcolor)
  } else if (length(colors)==1) {
    ## Uniform color -> nonuniform color
    colors <- rep(colors, vcount(tkp$graph))
    colors[vids+1] <- newcolor
    .tkplot.set(tkp.id, "params$vertex.color", colors)
  } else if (length(vids)==vcount(tkp$graph)) {
    ## Non-uniform -> uniform
    .tkplot.set(tkp.id, "params$vertex.color", newcolor)
  } else {
    ## Non-uniform -> non-uniform
    colors[vids+1] <- newcolor
    .tkplot.set(tkp.id, "params$vertex.color", colors)
  }

  tkitemconfigure(tkp$canvas, "selected&&vertex", "-fill", newcolor)
}

.tkplot.update.edge.color <- function(tkp.id, eids, newcolor) {
  tkp <- .tkplot.get(tkp.id)
  colors <- tkp$params$edge.color
  if (length(colors)==1 && length(eids)==ecount(tkp$graph)) {
    ## Uniform color -> uniform color
    .tkplot.set(tkp.id, "params$edge.color", newcolor)
  } else if (length(colors)==1) {
    ## Uniform color -> nonuniform color
    colors <- rep(colors, ecount(tkp$graph))
    colors[eids] <- newcolor
    .tkplot.set(tkp.id, "params$edge.color", colors)
  } else if (length(eids)==ecount(tkp$graph)) {
    ## Non-uniform -> uniform
    .tkplot.set(tkp.id, "params$edge.color", newcolor)
  } else {
    ## Non-uniform -> non-uniform
    colors[eids] <- newcolor
    .tkplot.set(tkp.id, "params$edge.color", colors)
  }

  tkitemconfigure(tkp$canvas, "selected&&edge", "-fill", newcolor)
}

.tkplot.update.edge.width <- function(tkp.id, eids, newwidth) {
  tkp <- .tkplot.get(tkp.id)
  widths <- tkp$params$edge.width
  if (length(widths)==1 && length(eids)==ecount(tkp$graph)) {
    ## Uniform width -> uniform width
    .tkplot.set(tkp.id, "params$edge.width", newwidth)
  } else if (length(widths)==1) {
    ## Uniform width -> nonuniform width
    widths <- rep(widths, ecount(tkp$graph))
    widths[eids] <- newwidth
    .tkplot.set(tkp.id, "params$edge.width", widths)
  } else if (length(eids)==ecount(tkp$graph)) {
    ## Non-uniform -> uniform
    .tkplot.set(tkp.id, "params$edge.width", newwidth)
  } else {
    ## Non-uniform -> non-uniform
    widths[eids] <- newwidth
    .tkplot.set(tkp.id, "params$edge.width", widths)
  }

  tkitemconfigure(tkp$canvas, "selected&&edge", "-width", newwidth)
}
  

.tkplot.update.vertex.size <- function(tkp.id, vids, newsize) {
  tkp <- .tkplot.get(tkp.id)
  sizes <- tkp$params$vertex.size
  if (length(sizes)==1 && length(vids)==vcount(tkp$graph)) {
    ## Uniform -> uniform
    .tkplot.set(tkp.id, "params$vertex.size", newsize)
  } else if (length(sizes)==1) {
    ## Uniform size -> nonuniform size
    sizes <- rep(sizes, vcount(tkp$graph))
    sizes[vids+1] <- newsize
    .tkplot.set(tkp.id, "params$vertex.size", sizes)
  } else if (length(vids)==vcount(tkp$graph)) {
    ## Non-uniform -> uniform
    .tkplot.set(tkp.id, "params$vertex.size", newsize)
  } else {
    ## Non-uniform -> non-uniform
    sizes[vids+1] <- newsize
    .tkplot.set(tkp.id, "params$vertex.size", sizes)
  }

  sapply(vids, function(id) {
    .tkplot.update.vertex(tkp.id, id, tkp$coords[id+1,1], tkp$coords[id+1,2])
  })
}

.tkplot.get.numeric.vector <- function(...) {
  labels <- list(...)
  if (length(labels)==0) return(FALSE)
  
  answers <- as.list(rep("", length(labels)))
  dialog <- tktoplevel()
  vars <- lapply(answers, tclVar)

  retval <- list()

  OnOK <- function() {
    retval <<- lapply(vars, tclvalue)
    tkdestroy(dialog)
  }
  
  OK.but <- tkbutton(dialog, text="   OK   ", command=OnOK)
  for (i in seq(along=labels)) {
    tkgrid(tklabel(dialog, text=labels[[i]]))
    tmp <- tkentry(dialog, width="40",textvariable=vars[[i]])
    tkgrid(tmp)
    tkbind(tmp, "<Return>", OnOK)
  }
  tkgrid(OK.but)
  tkwait.window(dialog)

  retval <- lapply(retval, function(v)
                   { eval(parse(text=paste("c(", v, ")"))) })
  return (retval)
}

.tkplot.select.number <- function(label, initial, low=1, high=100) {
  dialog <- tktoplevel()
  SliderValue <- tclVar(as.character(initial))
  SliderValueLabel <- tklabel(dialog,text=as.character(tclvalue(SliderValue)))
  tkgrid(tklabel(dialog,text=label), SliderValueLabel)
  tkconfigure(SliderValueLabel, textvariable=SliderValue)
  slider <- tkscale(dialog, from=high, to=low,
                    showvalue=F, variable=SliderValue,
                    resolution=1, orient="horizontal")  
  OnOK <- function() {
    SliderValue <<- as.numeric(tclvalue(SliderValue))
    tkdestroy(dialog)
  }
  OnCancel <- function() {
    SliderValue <<- NA
    tkdestroy(dialog)
  }
  OK.but <- tkbutton(dialog, text="   OK   ", command=OnOK)
  cancel.but <- tkbutton(dialog, text=" Cancel ", command=OnCancel)
  tkgrid(slider)
  tkgrid(OK.but, cancel.but)
  
  tkwait.window(dialog)
  return(SliderValue)
}
  
###################################################################
# Internal functions, vertex and edge selection
###################################################################

.tkplot.deselect.all <- function(tkp.id) {
  canvas <- .tkplot.get(tkp.id, "canvas")
  ids <- as.numeric(tkfind(canvas, "withtag", "selected"))
  for (i in ids) {
    .tkplot.deselect.this(tkp.id, i)
  }
}

.tkplot.select.all.vertices <- function(tkp.id) {
  canvas <- .tkplot.get(tkp.id, "canvas")
  vertices <- as.numeric(tkfind(canvas, "withtag", "vertex"))
  for (i in vertices) {
    .tkplot.select.vertex(tkp.id, i)
  }
}

.tkplot.select.some.vertices <- function(tkp.id, vids) {
  canvas <- .tkplot.get(tkp.id, "canvas")
  vids <- unique(vids)
  for (i in vids) {
    tkid <- as.numeric(tkfind(canvas, "withtag",
                              paste(sep="", "vertex&&v-", i)))
    .tkplot.select.vertex(tkp.id, tkid)
  }
}

.tkplot.select.all.edges <- function(tkp.id, vids) {
  canvas <- .tkplot.get(tkp.id, "canvas")
  edges <- as.numeric(tkfind(canvas, "withtag", "edge"))
  for (i in edges) {
    .tkplot.select.edge(tkp.id, i)
  }
}

.tkplot.select.some.edges <- function(tkp.id, from, to) {
  canvas <- .tkplot.get(tkp.id, "canvas")
  fromtags <- sapply(from, function(i) { paste(sep="", "from-", i) })
  totags <- sapply(from, function(i) { paste(sep="", "to-", i) })
  edges <- as.numeric(tkfind(canvas, "withtag", "edge"))
  for (i in edges) {
    tags <- as.character(tkgettags(canvas, i))
    ftag <- tags[ pmatch("from-", tags) ]
    ttag <- tags[ pmatch("to-", tags) ]
    if (ftag %in% fromtags && ttag %in% totags) {
      .tkplot.select.edge(tkp.id, i)
    }
  }
}

.tkplot.select.vertex <- function(tkp.id, tkid) {
  canvas <- .tkplot.get(tkp.id, "canvas")
  tkaddtag(canvas, "selected", "withtag", tkid)
  tkitemconfigure(canvas, tkid, "-outline", "red",
                  "-width", 2)
}

.tkplot.select.edge <- function(tkp.id, tkid) {
  canvas <- .tkplot.get(tkp.id, "canvas")
  tkaddtag(canvas, "selected", "withtag", tkid)
  tkitemconfigure(canvas, tkid, "-dash", "-")
}

.tkplot.select.label <- function(tkp.id, tkid) {
  canvas <- .tkplot.get(tkp.id, "canvas")
  tkaddtag(canvas, "selected", "withtag", tkid)
}  

.tkplot.deselect.vertex <- function(tkp.id, tkid) {
  canvas <- .tkplot.get(tkp.id, "canvas")
  tkdtag(canvas, tkid, "selected")
  tkp <- .tkplot.get(tkp.id)
  tags <- as.character(tkgettags(canvas, tkid))
  id <- as.numeric(substring(tags[pmatch("v-", tags)], 3))
  vertex.frame.color <- ifelse(length(tkp$params$vertex.frame.color)>1,
                               tkp$params$vertex.frame.color[id+1],
                               tkp$params$vertex.frame.color)  
  tkitemconfigure(canvas, tkid, "-outline", vertex.frame.color,
                  "-width", 1)
}

.tkplot.deselect.edge <- function(tkp.id, tkid) {
  canvas <- .tkplot.get(tkp.id, "canvas")
  tkdtag(canvas, tkid, "selected")
  tkp <- .tkplot.get(tkp.id)
  tags <- as.character(tkgettags(canvas, tkid))
  id <- as.numeric(substring(tags[pmatch("edge-", tags)], 6))
  edge.lty <- ifelse(length(tkp$params$edge.lty)>1,
                     tkp$params$edge.lty[[id]],
                     tkp$params$edge.lty)
  tkitemconfigure(canvas, tkid, "-dash", edge.lty)
}

.tkplot.deselect.label <- function(tkp.id, tkid) {
  canvas <- .tkplot.get(tkp.id, "canvas")
  tkdtag(canvas, tkid, "selected")
}

.tkplot.select.current <- function(tkp.id) {
  canvas <- .tkplot.get(tkp.id, "canvas")
  tkid <- as.numeric(tkfind(canvas, "withtag", "current"))
  .tkplot.select.this(tkp.id, tkid)
}

.tkplot.deselect.current <- function(tkp.id) {
  canvas <- .tkplot.get(tkp.id, "canvas")
  tkid <- as.numeric(tkfind(canvas, "withtag", "current"))
  .tkplot.deselect.this(tkp.id, tkid)
}

.tkplot.select.this <- function(tkp.id, tkid) {
  canvas <- .tkplot.get(tkp.id, "canvas")
  tags <- as.character(tkgettags(canvas, tkid))
  if ("vertex" %in% tags) {
    .tkplot.select.vertex(tkp.id, tkid)
  } else if ("edge" %in% tags) {
    .tkplot.select.edge(tkp.id, tkid)
  } else if ("label" %in% tags) {
    tkp <- .tkplot.get(tkp.id)
    if (tkp$params$label.dist == 0) {
      id <- tags[pmatch("v-", tags)]
      tkid <- as.character(tkfind(canvas, "withtag",
                                  paste(sep="", "vertex&&", id)))
      .tkplot.select.vertex(tkp.id, tkid)
    } else {
      .tkplot.select.label(tkp.id, tkid)
    }
  } 
}

.tkplot.deselect.this <- function(tkp.id, tkid) {
  canvas <- .tkplot.get(tkp.id, "canvas")
  tags <- as.character(tkgettags(canvas, tkid))
  if ("vertex" %in% tags) {
    .tkplot.deselect.vertex(tkp.id, tkid)
  } else if ("edge" %in% tags) {
    .tkplot.deselect.edge(tkp.id, tkid)
  } else if ("label" %in% tags) {
    tkp <- .tkplot.get(tkp.id)
    if (tkp$params$label.dist == 0) {
      id <- tags[pmatch("v-", tags)]
      tkid <- as.character(tkfind(canvas, "withtag",
                                  paste(sep="", "vertex&&", id)))
      .tkplot.deselect.vertex(tkp.id, tkid)
    } else {
      .tkplot.deselect.label(tkp.id, tkid)
    }
  }
}

.tkplot.get.selected.vertices <- function(tkp.id) {
  canvas <- .tkplot.get(tkp.id, "canvas")
  tkids <- as.numeric(tkfind(canvas, "withtag", "vertex&&selected"))

  ids <- sapply(tkids, function(tkid) {
    tags <- as.character(tkgettags(canvas, tkid))
    id <- as.numeric(substring(tags [pmatch("v-", tags)], 3))
    id})

  ids
}

.tkplot.get.selected.edges <- function(tkp.id) {
  canvas <- .tkplot.get(tkp.id, "canvas")
  tkids <- as.numeric(tkfind(canvas, "withtag", "edge&&selected"))

  ids <- sapply(tkids, function(tkid) {
    tags <- as.character(tkgettags(canvas, tkid))
    id <- as.numeric(substring(tags [pmatch("edge-", tags)], 6))
    id})

  ids
}

###################################################################
# Internal functions: manipulating the UI
###################################################################

.tkplot.select.menu <- function(tkp.id, main.menu) {
  select.menu <- tkmenu(main.menu)

  tkadd(select.menu, "command", label="Select all vertices",
        command=function() {
          .tkplot.deselect.all(tkp.id)
          .tkplot.select.all.vertices(tkp.id)
        })
  tkadd(select.menu, "command", label="Select all edges",
        command=function() {
          .tkplot.deselect.all(tkp.id)
          .tkplot.select.all.edges(tkp.id)
        })
  tkadd(select.menu, "command", label="Select some vertices...",
        command=function() {
          vids <- .tkplot.get.numeric.vector("Select vertices")
          .tkplot.select.some.vertices(tkp.id, vids[[1]])
        })
  tkadd(select.menu, "command", label="Select some edges...",
        command=function() {
          fromto <- .tkplot.get.numeric.vector("Select edges from vertices",
                                               "to vertices")
          .tkplot.select.some.edges(tkp.id, fromto[[1]], fromto[[2]])
        })
  tkadd(select.menu, "separator")
  tkadd(select.menu, "command", label="Deselect everything",
        command=function() { .tkplot.deselect.all(tkp.id) })

  select.menu
}

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
  
  dialog <- tktoplevel(.tkplot.get(tkp.id, "top"))
  
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
    } else if (layout$param[[i]]$type=="expression") {
      values[[i]] <- tkentry(dialog)
      .tkplot.g <- .tkplot.get(tkp.id, "graph")
      tkinsert(values[[i]], 0, as.character(eval(layout$params[[i]]$default)))
      tkgrid(values[[i]], row=row, column=1, sticky="nw", padx=5, pady=5)
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

.tkplot.select.color <- function(initialcolor) {
  
  color <- tclvalue(tcl("tk_chooseColor", initialcolor=initialcolor,
                        title="Choose a color"))
  return(color);
}

###################################################################
# Internal functions: other
###################################################################

.tkplot.convert.color <- function(col) {
  if (is.numeric(col)) {
    ## convert numeric color based on current palette
    p <- palette()
    col <- col %% length(p)
    col[col==0] <- length(p)
    col <- palette()[col]
  } else if (is.character(col) && any(substr(col,1,1)=="#" & nchar(col)==9)) {
    ## drop alpha channel, tcltk doesn't support it
    idx <- substr(col,1,1)=="#" & nchar(col)==9
    col[idx] <- substr(col[idx],1,7)
  }
  
  col
}

.tkplot.convert.font <- function(font, family, cex) {
  tk.fonts <- as.character(tkfont.names())
  if (as.character(font) %in% tk.fonts) {
    ## already defined Tk font
    as.character(font)
  } else {
    ## we create a font from familiy, font & cex
    font <- as.numeric(font)
    family <- as.character(family)
    cex <- as.numeric(cex)    

    ## set slant & weight
    if (font==2) {
      slant <- "roman"
      weight <- "bold"
    } else if (font==3) {
      slant <- "italic"
      weight <- "normal"
    } else if (font==4) {
      slant <- "italic"
      weight <- "bold"
    } else {
      slant <- "roman"
      weight <- "normal"
    }

    ## set tkfamily
    if (family=="symbol" || font==5) {
      tkfamily <- "symbol"
    } else if (family=="serif") {
      tkfamily <- "Times"
    } else if (family=="sans") {
      tkfamily <- "Helvetica"
    } else if (family=="mono") {
      tkfamily <- "Courier"
    } else {
      ## pass the family and see what happens
      tkfamily <- family
    }
    
    newfont <- tkfont.create(family=tkfamily, slant=slant, weight=weight,
                             size=12*cex)
    as.character(newfont)
  }
}

i.tkplot.get.edge.lty <- function(edge.lty) {

  if (is.numeric(edge.lty)) {
    lty <- c( " ", "", "-", ".", "-.", "--", "--.")
    edge.lty <- lty[edge.lty %% 7 + 1]
  } else if (is.character(edge.lty)) {
    wh <- edge.lty %in% c("blank", "solid", "dashed", "dotted", "dotdash",
                          "longdash", "twodash")
    lty <- c( " ", "", "-", ".", "-.", "--", "--.")
    names(lty) <- c("blank", "solid", "dashed", "dotted", "dotdash",
                    "longdash", "twodash")
    edge.lty[wh] <- lty[ edge.lty[wh] ]
  }
  edge.lty
}
