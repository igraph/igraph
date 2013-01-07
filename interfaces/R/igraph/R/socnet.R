
# TODO LIST:
#   * adding edges to a graph
#   * exporting graphics
#   * scroll bar for the graph list area   == IMPOSSIBLE right now, should be a list
#   * window title in the error dialog
#   * keyboard shortcuts
#   * implement min & max in .tkigraph.dialog

.tkigraph.env <- new.env()

tkigraph <- function() {

  require(tcltk) || stop("tcl/tk library not available")

  options(scipen=10000)
  
  if (!exists("window", envir=.tkigraph.env, inherits=FALSE)) {
    assign("window", TRUE, envir=.tkigraph.env)
    assign("graphs", list(), envir=.tkigraph.env)
    assign("selected", list(), envir=.tkigraph.env)
    assign("tklines", list(), envir=.tkigraph.env)
  } else {
    stop("tkigraph window is already open!")
  }
  
  # Create top window
  top <- tktoplevel(background="lightgrey", width=700, height=400)
  tktitle(top) <- "iGraph GUI (Social Network Basics)"
  topframe <- tkframe(top, relief="sunken", borderwidth=1)
  scr <- tkscrollbar(top, repeatinterval=5,
                     command=function(...) tkyview(topframe))
  tkplace(topframe, x=0, y=0, relwidth=1.0)
  
  # Store myself in the environment if needed
  if (!exists("top", envir=.tkigraph.env, inherits=FALSE)) {
    assign("top", top, envir=.tkigraph.env)
    assign("topframe", topframe, envir=.tkigraph.env)
  }
  
  # kill myself if window was closed
  tkbind(top, "<Destroy>", function() .tkigraph.close())

  # pull-down menu
  main.menu <- tkmenu(top)

  graph.menu <- tkmenu(main.menu)
  
  create.menu <- tkmenu(main.menu)
  tkadd(create.menu, "command", label="By hand", command=function() {
    .tkigraph.by.hand()
  })
  tkadd(create.menu, "separator")
  tkadd(create.menu, "command", label="Ring", command=function() {
    .tkigraph.ring()
  })
  tkadd(create.menu, "command", label="Tree", command=function() {
    .tkigraph.tree()
  })
  tkadd(create.menu, "command", label="Lattice", command=function() {
    .tkigraph.lattice()
  })
  tkadd(create.menu, "command", label="Star", command=function() {
    .tkigraph.star()
  })
  tkadd(create.menu, "command", label="Full", command=function() {
    .tkigraph.full()
  })
  tkadd(create.menu, "separator")
  tkadd(create.menu, "command", label="Graph atlas...", command=function() {
    .tkigraph.atlas()
  })
  tkadd(create.menu, "separator")
  tkadd(create.menu, "command", label="Moody-White network", command=function() {
    g <- graph.adjacency(.tkigraph.net.moody.white, mode="undirected")
    g <- set.graph.attribute(g, "name", "Moody-White network")
    .tkigraph.add.graph(g)
  })
  
  tkadd(create.menu, "separator") 
  tkadd(create.menu, "command", label="Random (Erdos-Renyi G(n,p))",
        command=function() {
          .tkigraph.erdos.renyi.game()
        })
  tkadd(create.menu, "command", label="Random (Erdos-Renyi G(n,m))",
        command=function() { .tkigraph.erdos.renyi.gnm.game() })
  tkadd(create.menu, "command", label="Random (Barabasi-Albert)",
        command=function() {
          .tkigraph.barabasi.game()
        })
  tkadd(create.menu, "command", label="Random (Configuration model)",
        command=function() {
          .tkigraph.degree.sequence.game()
        })
  tkadd(create.menu, "command", label="Watts-Strogatz random graph",
        command=function() {
          .tkigraph.watts.strogatz()
        })
  tkadd(create.menu, "separator")
  tkadd(create.menu, "command", label="Simplify", command=function() {
    .tkigraph.simplify()
  })
  tkadd(graph.menu, "cascade", label="Create", menu=create.menu)
  tkadd(graph.menu, "command", label="Delete", command=function() {
    .tkigraph.delete() })  
  tkadd(graph.menu, "separator")
  tkadd(graph.menu, "command", label="Show graph",
        command=function() { .tkigraph.show() })
  tkadd(graph.menu, "command", label="Basic statistics",
        command=function() { .tkigraph.stat() })
  tkadd(graph.menu, "separator")
  tkadd(graph.menu, "command", label="Import session", command=function() {
    .tkigraph.load()
  })
#  tkadd(graph.menu, "command", label="Load from the Web", command=function() {
#    .tkigraph.load.online()
#  })
  tkadd(graph.menu, "command", label="Export session", command=function() {
    .tkigraph.save()
  })
  tkadd(graph.menu, "separator")
  tkadd(graph.menu, "command", label="Import adjacency matrix",
        command=function() .tkigraph.import.adjacency())
  tkadd(graph.menu, "command", label="Import edge list",
        command=function() .tkigraph.import.edgelist())
  tkadd(graph.menu, "command", label="Import Pajek file",
        command=function() .tkigraph.import.pajek())
  tkadd(graph.menu, "command", label="Export adjacency matrix",
        command=function() .tkigraph.export.adjacency())
  tkadd(graph.menu, "command", label="Export edge list",
        command=function() .tkigraph.export.edgelist())
  tkadd(graph.menu, "command", label="Export Pajek file",
        command=function() .tkigraph.export.pajek())
  tkadd(main.menu, "cascade", label="Graph", menu=graph.menu)

  plot.menu <- tkmenu(main.menu)
  tkadd(plot.menu, "command", label="Simple", command=function() {
    .tkigraph.plot(simple=TRUE)
  })
  tkadd(plot.menu, "command", label="Advanced", command=function() {
    .tkigraph.plot(simple=FALSE)
  })
  tkadd(main.menu, "cascade", label="Draw", menu=plot.menu)
  
  centrality.menu <- tkmenu(main.menu)
  tkadd(centrality.menu, "command", label="Degree (out)", command=function() {
    .tkigraph.degree("out")
  })
  tkadd(centrality.menu, "command", label="Degree (in)", command=function() {
    .tkigraph.degree("in")
  })
  tkadd(centrality.menu, "command", label="Degree (total)",
        command=function() {
          .tkigraph.degree("total")
        })
  tkadd(centrality.menu, "command", label="Plot log-log degree distribution",
        command=function() {
          .tkigraph.degree.dist(power=FALSE)
        })
  tkadd(centrality.menu, "command", label="Fit a power-law to degree distribution",
        command=function() {
          .tkigraph.degree.dist(power=TRUE)
        })
  tkadd(centrality.menu, "separator")
  tkadd(centrality.menu, "command", label="Closeness", command=function() {
    .tkigraph.closeness()
  })
  tkadd(centrality.menu, "command", label="Betweenness", command=function() {
    .tkigraph.betweenness()
  })
  tkadd(centrality.menu, "command", label="Burt's constraint", command=function() {
    .tkigraph.constraints()
  })
  tkadd(centrality.menu, "command", label="Page rank", command=function() {
    .tkigraph.page.rank()
  })
  tkadd(centrality.menu, "separator")  
  tkadd(centrality.menu, "command", label="Edge betweenness",
        command=function() {
          .tkigraph.edge.betweenness()
        })
  tkadd(main.menu, "cascade", label="Centrality", menu=centrality.menu)

  distances.menu <- tkmenu(main.menu)
  tkadd(distances.menu, "command", label="Distance matrix",
        command=function() { .tkigraph.dist.matrix() })
  tkadd(distances.menu, "command", label="Distances from/to vertex",
        command=function() { .tkigraph.distance.tofrom() })
  tkadd(distances.menu, "command", label="Diameter",
        command=function() { .tkigraph.diameter() })
  tkadd(distances.menu, "command", label="Draw diameter",
        command=function() { .tkigraph.plot.diameter(simple=FALSE) })
  tkadd(distances.menu, "command", label="Average path length",
        command=function() { .tkigraph.diameter(mode="path") })
  tkadd(main.menu, "cascade", label="Distances", menu=distances.menu)

  component.menu <- tkmenu(main.menu)
  tkadd(component.menu, "command", label="Show components",
        command=function() { .tkigraph.clusters() })
  tkadd(component.menu, "command", label="Show membership",
        command=function() { .tkigraph.clusters.membership() })
  tkadd(component.menu, "command", label="Calculate component sizes",
        command=function() { .tkigraph.calculate.clusters() })
  tkadd(component.menu, "command", label="Draw components",
        command=function() { .tkigraph.plot.comp(simple=FALSE) })
  tkadd(component.menu, "command", label="Create graph from giant component",
        command=function() { .tkigraph.create.giantcomp() })
  tkadd(component.menu, "command", label="Create graph from component of a vertex",
        command=function() { .tkigraph.create.mycomp() })
  tkadd(component.menu, "command", label="Create graph from a component",
        command=function() { .tkigraph.create.comp() })

  community.menu <- tkmenu(main.menu)
  tkadd(community.menu, "command", label="Spinglass algorithm",
        command=function() { .tkigraph.spinglass() })
  tkadd(community.menu, "command", label="Spinglass algorithm, single vertex",
        command=function() { .tkigraph.my.spinglass() })

  cohesion.menu <- tkmenu(main.menu)
  tkadd(cohesion.menu, "command", label="Cohesion of all components",
        command=function() { .tkigraph.cohesion() })
  
  subgraph.menu <- tkmenu(main.menu)
  tkadd(subgraph.menu, "cascade", label="Components", menu=component.menu)
  tkadd(subgraph.menu, "cascade", label="Communities", menu=community.menu)
  tkadd(subgraph.menu, "cascade", label="Cohesion", menu=cohesion.menu)
  
  tkadd(main.menu, "cascade", label="Subgraphs", menu=subgraph.menu)
  
  motif.menu <- tkmenu(main.menu)
  tkadd(motif.menu, "command", label="Draw motifs", command=function() {
    .tkigraph.motifs.draw()
  })
  tkadd(motif.menu, "command", label="Find motifs", command=function() {
    .tkigraph.motifs.find()
  })
  tkadd(main.menu, "cascade", label="Motifs", menu=motif.menu)

  help.menu <- tkmenu(main.menu)
  tkadd(help.menu, "command", label="Contents", command=function() { .tkigraph.help() })
  tkadd(help.menu, "command", label="In external browser",
        command=function() { .tkigraph.help.external() })
  tkadd(help.menu, "separator")
  tkadd(help.menu, "command", label="About", command=function() { .tkigraph.about() })
  tkadd(main.menu, "cascade", label="Help", menu=help.menu)
  
  tkadd(main.menu, "command", label="Quit", command=.tkigraph.close)
  
  tkconfigure(top, "-menu", main.menu)

  # Set up the main area
  tkgrid(tklabel(top, text=""),
         tklabel(top, text="#", justify="center", relief="raised"),
         tklabel(top, text="Name", width=50, relief="raised",
                 justify="left"),
         tklabel(top, text="|V|", width=6, relief="raised",
                 justify="left"),
         tklabel(top, text="|E|", width=6, relief="raised",
                 justify="left"),
         tklabel(top, text="Dir.", width=6, relief="raised",
                 justify="left"),
         sticky="nsew", "in"=topframe)
  tkgrid.columnconfigure(topframe, 2, weight=1)

  invisible(NULL)
}

.tkigraph.close <- function() {
  message <- "Are you sure?"
  yesno <- tkmessageBox(message=message, icon="question", type="yesno",
                        default="yes")
  if (as.character(yesno) == "no") { return() }
  top <- get("top", .tkigraph.env)
  tkbind(top, "<Destroy>", "")
  tkdestroy(top)
  rm(list=ls(envir=.tkigraph.env), envir=.tkigraph.env)
}

.tkigraph.get.selected <- function() {
  gnos <- get("selected", .tkigraph.env)
  which(as.logical(sapply(gnos, tclvalue)))
}

.tkigraph.error <- function(message) {
  tkmessageBox(message=message, icon="error", type="ok")
}

.tkigraph.warning <- function(message) {
  tkmessageBox(message=message, icon="warning", type="ok")
}

.tkigraph.dialogbox <- function(TITLE="Setup parameters", ...) {
  
  params <- list(...)
  answers <- lapply(params, "[[", "default")
  dialog <- tktoplevel()
  frame <- tkframe(dialog)
  tkgrid(frame)
  tktitle(dialog) <- TITLE
  vars <- lapply(answers, tclVar)
  retval <- list()
  widgets <- list()

  OnOK <- function() {
    retval <<- lapply(vars, tclvalue)
    for (i in seq(along=params)) {
      if (params[[i]]$type == "listbox") {
        retval[[i]] <<- as.numeric(tclvalue(tkcurselection(widgets[[i]])))
      }
    }
    tkdestroy(dialog)    
  }

  tkgrid(tklabel(dialog, text=TITLE,
                 font=tkfont.create(family="times", size="16", weight="bold")),
         columnspan=2, sticky="nsew", "in"=frame, padx=10, pady=10)
  
  OK.but <- tkbutton(dialog, text="   OK    ", command=OnOK)
  for (i in seq(along=params)) {
    tkgrid(tklabel(dialog, text=params[[i]]$name),
           column=0, row=i, sticky="nw", padx=10, "in"=frame)
    if (params[[i]]$type == "numeric" || params[[i]]$type == "text") {
      tmp <- tkentry(dialog, width="10", textvariable=vars[[i]])
      tkgrid(tmp, column=1, row=i, sticky="nsew", padx=10, "in"=frame)
      tkbind(tmp, "<Return>", OnOK)
    } else if (params[[i]]$type == "boolean") {
      b <- tkcheckbutton(dialog, onvalue="TRUE", offvalue="FALSE",
                           variable=vars[[i]])
      if (params[[i]]$default == "TRUE") { tkselect(b) }
      tkgrid(b, column=1, row=i, sticky="w", padx=10, "in"=frame)
    } else if (params[[i]]$type == "listbox") {
      f <- tkframe(dialog)
      tkgrid(f, "in"=frame, padx=10, sticky="nsew", column=1, row=i)
      scr <- tkscrollbar(f, repeatinterval=5)
      fun <- eval(eval(substitute(expression(function(...) tkset(scr,...)),
                                  list(scr=scr))))
      lb <- tklistbox(f, selectmode="single", exportselection=FALSE,
                      height=3, yscrollcommand=fun)
      fun <- eval(eval(substitute(expression(function(...) tkyview(lb, ...)),
          list(lb=lb))))
      tkconfigure(scr, "-command", fun)
      tkselection.set(lb, as.numeric(params[[i]]$default)+1)
      lapply(params[[i]]$values, function(l) tkinsert(lb, "end", l))
      tkselection.set(lb, as.numeric(params[[i]]$default))
      tkgrid(lb, scr, sticky="nsew", "in"=f)
      tkgrid.configure(scr, sticky="nsw")
      tkgrid.columnconfigure(f, 0, weight=1)
      widgets[[i]] <- lb
    }
  }
  tkgrid(OK.but, column=0, columnspan=2, sticky="nsew", "in"=frame, pady=10,
         padx=10)
  tkgrid.columnconfigure(frame, 1, weight=1)
  tkwait.window(dialog)

  for (i in seq(retval)) {
    if (params[[i]]$type == "numeric") {
      retval[[i]] <- eval(parse(text=retval[[i]]))
    } else if (params[[i]]$type == "text") {
      retval[[i]] <- eval(retval[[i]])
    } else if (params[[i]]$type == "boolean") {
      if (retval[[i]] == "FALSE") {
        retval[[i]] <- FALSE
      } else {
        retval[[i]] <- TRUE
      }
    } else if (params[[i]]$type == "listbox") {
      ## nothing to do
    }
  }
  names(retval) <- names(params)
  return (retval)  
}

.tkigraph.add.graph <- function(g) {
  top <- get("top", .tkigraph.env)
  topframe <- get("topframe", .tkigraph.env)

  ## add 'name' attribute if not present
  if (!"name" %in% list.vertex.attributes(g)) {
    V(g)$name <- as.integer(seq(vcount(g)))
  }
  if (!"name" %in% list.edge.attributes(g)) {
    E(g)$name <- as.integer(seq(ecount(g)))
  }
  
  graphs <- get("graphs", .tkigraph.env)
  selected <- get("selected", .tkigraph.env)
  assign("graphs", append(graphs, list(g)), .tkigraph.env)
  no <- length(graphs)+1

  selected[[no]] <- tclVar("FALSE")
  assign("selected", selected, .tkigraph.env)

  name <- get.graph.attribute(g, "name")
  tmpvar <- tclVar(as.character(name))
  but <- tkcheckbutton(top, onvalue="TRUE", offvalue="FALSE",
                       variable=selected[[no]])
  lab <- tklabel(top, text=as.character(no), width=2)
  ent <- tkentry(top, width=30, textvariable=tmpvar)
  lab2 <- tklabel(top, text=as.character(vcount(g)),
                  justify="right", padx=2)
  lab3 <- tklabel(top, text=as.character(ecount(g)), justify="right",
                  padx=2)
  lab4 <- tklabel(top, text=if (is.directed(g)) "YES" else "NO")
  tkgrid(but, lab, ent, lab2, lab3, lab4, "in"=topframe, sticky="nsew")

  tklines <- get("tklines", .tkigraph.env)
  tklines[[no]] <- list(but, lab, ent, lab2, lab3, lab4)
  assign("tklines", tklines, .tkigraph.env)
}

.tkigraph.delete <- function() {
  gnos <- .tkigraph.get.selected()  

  if (length(gnos) == 0) { return() }
  if (length(gnos) > 1) {
    message <- paste("Are you sure to delete", length(gnos), "graphs?")
  } else {
    message <- paste("Are you sure to delete graph #", gnos, "?")
  }
  yesno <- tkmessageBox(message=message, icon="question", type="yesno",
                        default="yes")
  if (as.character(yesno) == "no") { return() }

  ## remove from the screen
  graphs <- get("graphs", .tkigraph.env)
  topframe <- get("topframe", .tkigraph.env)
  todel <- get("tklines", .tkigraph.env)[gnos]
  todel <- unlist(recursive=FALSE, todel)
  for (i in todel) {
    tkgrid.remove(topframe, i)
  }
  ## delete the graphs
  graphs[gnos] <- NA
  assign("graphs", graphs, .tkigraph.env)
  selected <- get("selected", .tkigraph.env)
  for (i in gnos) { 
    selected[[i]] <- tclVar("FALSE")
  }
  assign("selected", selected, .tkigraph.env)
}

.tkigraph.load <- function() {
  filename <- tkgetOpenFile(defaultextension="Rdata",
                            title="Load graphs")
  env <- new.env()
  load(paste(as.character(filename), collapse=" "), envir=env)
  .tkigraph.graphs <- get("graphs", envir=env) 
  for (i in seq(.tkigraph.graphs)) {
    .tkigraph.add.graph(.tkigraph.graphs[[i]])
  }
  if (".tkigraph.graphs" %in% ls(all.names=TRUE)) {
    rm(.tkigraph.graphs)
  }
}

.tkigraph.load.online <- function() {
  ## TODO
}

.tkigraph.save <- function() {
  graphs <- get("graphs", .tkigraph.env)
  topframe <- get("topframe", .tkigraph.env)
  for (i in seq(graphs)) {
    if (is.na(graphs)[i]) { next }
    entry <- tkgrid.slaves(topframe, row=i, col=2)
    graphs[[i]] <- set.graph.attribute(graphs[[i]], "name",
                                       as.character(tcl(entry, "get")))
  }
  graphs <- graphs[ !is.na(graphs) ]
  filename <- tkgetSaveFile(initialfile="graphs.Rdata",
                            defaultextension="Rdata",
                            title="Save graphs")
  save(graphs, file=paste(as.character(filename), collapse=" "))
}

.tkigraph.import.adjacency <- function() {
  filename <- tkgetOpenFile(defaultextension="adj",
                            title="Import adjacency matrix")
  filename <- paste(as.character(filename), collapse=" ")
  if (filename=="") { return() }
  tab <- read.table(filename)
  tab <- as.matrix(tab)
  if (ncol(tab) != nrow(tab)) {
    .tkigraph.error("Cannot interpret as adjacency matrix")
    return()
  }
  dir <- if (all(t(tab)==tab)) "undirected" else "directed"
  if (all(unique(tab) %in% c(0,1))) {
    weighted <- NULL
  } else {
    weighted <- "weight"
  }
  g <- .tkigraph.graph.adjacency(tab, mode=dir, weighted=weighted)
  g <- set.graph.attribute(g, "name", "Imported adjacency matrix")
  .tkigraph.add.graph(g)
}

.tkigraph.graph.adjacency <- function(adjmatrix, mode, weighted) {
  if (is.null(weighted)) {
    g <- graph.adjacency(adjmatrix, mode=mode)
  } else {
    ## there is bug in the currect igraph version, this is a workaround
    if (mode=="undirected") {
      adjmatrix[ lower.tri(adjmatrix) ] <- 0
    }
    g <- graph.adjacency(adjmatrix, mode=mode, weighted=weighted)
  }
  g
}

.tkigraph.import.edgelist <- function() {
  filename <- tkgetOpenFile(defaultextension="el",
                            title="Import edge list")
  filename <- paste(as.character(filename), collapse=" ")
  if (filename=="") { return() }
  tab <- read.table(filename,
                    colClasses="character")
  cn <- rep("", ncol(tab))
  if (ncol(tab)>=3) { cn[3] <- "weight" }
  colnames(tab) <- cn
  read <- .tkigraph.dialogbox(TITLE="Importing an edge list",
                              directed=list(name="Directed", type="boolean",
                                default="FALSE"))
  g <- graph.data.frame(tab, directed=read$directed)
  g <- set.graph.attribute(g, "name", "Imported edge list")
  .tkigraph.add.graph(g)
}

.tkigraph.import.pajek <- function() {
  filename <- tkgetOpenFile(defaultextension="net",
                            title="Import Pajek file")
  filename <- paste(as.character(filename), collapse=" ")
  if (filename=="") { return() }
  g <- read.graph(file=filename, format="pajek")
  color <- NULL # To eliminate a check NOTE
  if ("color" %in% list.vertex.attributes(g)) { V(g)[ color=="" ]$color <- "black" }
  if ("color" %in% list.edge.attributes(g)) { E(g)[ color=="" ]$color <- "black" }
  g <- set.graph.attribute(g, "name", "Imported Pajek fie")
  .tkigraph.add.graph(g)
}

.tkigraph.export.adjacency <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos)!=1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graph <- get("graphs", .tkigraph.env)[[gnos]]
  if ("weight" %in% list.graph.attributes(graph)) {
    tab <- get.adjacency(graph, attr="weight", names=FALSE)
  } else {
    tab <- get.adjacency(graph, names=FALSE)
  }
  filename <- tkgetSaveFile(initialfile="graph.adj",
                            defaultextension="adj",
                            title="Export adjacency matrix")
  filename <- paste(as.character(filename), collapse=" ")
  if (filename=="") { return() }
  write.table(tab, file=filename, row.names=FALSE, col.names=FALSE)
}

.tkigraph.export.edgelist <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos)!=1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graph <- get("graphs", .tkigraph.env)[[gnos]]
  el <- get.edgelist(graph)
  if ("weight" %in% list.edge.attributes(graph)) {
    el <- cbind(el, E(graph)$weight)
  }
  filename <- tkgetSaveFile(initialfile="graph.el",
                            defaultextension="el",
                            title="Export edge list")
  filename <- paste(as.character(filename), collapse=" ")
  if (filename=="") { return() }
  write.table(el, file=filename,
              row.names=FALSE, col.names=FALSE)
}

.tkigraph.export.pajek <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos)!=1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graph <- get("graphs", .tkigraph.env)[[gnos]]
  filename <- tkgetSaveFile(initialfile="pajek.net",
                            defaultextension="net",
                            title="Export Pajek file")
  filename <- paste(as.character(filename), collapse=" ")
  if (filename=="") { return() }
  write.graph(graph, file=filename, format="pajek")
}

.tkigraph.show <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos)!=1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graphs <- get("graphs", .tkigraph.env)
  el <- get.edgelist(graphs[[gnos]])
  el <- data.frame(from=el[,1], to=el[,2])
#  if (any(V(graphs[[gnos]])$name != seq(length=vcount(graphs[[gnos]])))) {
#    el2 <- get.edgelist(graphs[[gnos]], names=FALSE)
#    el <- cbind(el, el2)
#  }
  if ("weight" %in% list.edge.attributes(graphs[[gnos]])) {
    el <- cbind(el, value=E(graphs[[gnos]])$weight)
  }

  .tkigraph.showData(el, title=paste(sep="", "Graph #", gnos), right=FALSE)
}

.tkigraph.stat <- function() {
  gnos <- .tkigraph.get.selected()

  if (length(gnos) == 0) {
    .tkigraph.error("Please select some graphs")
    return()
  }

  read <- .tkigraph.dialogbox(TITLE="Choose statistics",
                              vertices=list(name="Vertices", type="boolean",
                                default="FALSE"),
                              edges=list(name="Edges", type="boolean",
                                default="FALSE"),
                              recip=list(name="Reciprocity", type="boolean",
                                default="FALSE"),
                              dens=list(name="Density", type="boolean",
                                default="FALSE"),
                              trans=list(name="Transitivity (global)",
                                type="boolean", default="FALSE"),
                              ltrans=list(name="Mean local transitivity",
                                type="boolean", default="FALSE"),
                              deg=list(name="Average degree",
                                type="boolean", default="FALSE"),
                              maxdeg=list(name="Maximum degree (total)",
                                type="boolean", default="FALSE"),
                              maxindeg=list(name="Maximum degree (in)",
                                type="boolean", default="FALSE"),
                              maxoutdeg=list(name="Maximum degree (out)",
                                type="boolean", default="FALSE"),
                              mindeg=list(name="Minimum degree (total)",
                                type="boolean", default="FALSE"),
                              minindeg=list(name="Minimum degree (in)",
                                type="boolean", default="FALSE"),
                              minoutdeg=list(name="Minimum degree (out)",
                                type="boolean", default="FALSE")
                              )
  
  graphs <- get("graphs", .tkigraph.env)[gnos]

  v <- e <- recip <- dens <- trans <- ltrans <- 
    deg <- maxdeg <- maxindeg <- maxoutdeg <-
      mindeg <- minindeg <- minoutdeg <- numeric()
  for (i in seq(along=gnos)) {
    if (read$vertices) {
      v[i] <- vcount( graphs[[ i ]] )
    }
    if (read$edges) {
      e[i] <- ecount( graphs[[ i ]] )
    }
    if (read$recip) {
      recip[i] <- reciprocity( graphs[[ i ]] )
    }
    if (read$dens) {
      dens[i] <- graph.density( graphs[[ i ]] )
    }
    if (read$trans) {
      trans[i] <- transitivity( graphs[[ i ]], type="global")
    }
    if (read$ltrans) {
      ltrans[i] <- transitivity( graphs[[ i ]], type="localaverage")
    }
    if (read$deg) {
      deg[i] <- mean(degree( graphs[[ i ]], mode="total"))
    }
    if (read$maxdeg) {
      maxdeg[i] <- max(degree( graphs[[ i ]], mode="total"))
    }
    if (read$maxindeg) {
      maxindeg[i] <- max(degree( graphs[[ i ]], mode="in"))
    }
    if (read$maxoutdeg) {
      maxoutdeg[i] <- max(degree( graphs[[ i ]], mode="out"))
    }
    if (read$mindeg) {
      mindeg[i] <- min(degree( graphs[[ i ]], mode="total"))
    }
    if (read$minindeg) {
      minindeg[i] <- min(degree( graphs[[ i ]], mode="in"))
    }
    if (read$minoutdeg) {
      minoutdeg[i] <- min(degree( graphs[[ i ]], mode="out"))
    }
  }

  value <- numeric()
  cn <- character()
  if (read$vertices) {
    value <- cbind(value, v)
    cn <- c(cn, "Vertices")
  }
  if (read$edges) {
    value <- cbind(value, e)
    cn <- c(cn, "Edges")
  }
  if (read$recip) {
    value <- cbind(value, recip)
    cn <- c(cn, "Reciprocity")
  }
  if (read$dens) {
    value <- cbind(value, dens)
    cn <- c(cn, "Density")
  }
  if (read$trans) {
    value <- cbind(value, trans)
    cn <- c(cn, "Transitivity")
  }
  if (read$ltrans) {
    value <- cbind(value, ltrans)
    cn <- c(cn, "Local trans.")
  }
  if (read$deg) {
    value <- cbind(value, deg)
    cn <- c(cn, "Mean degree")
  }
  if (read$maxdeg) {
    value <- cbind(value, maxdeg)
    cn <- c(cn, "Max. degree")
  }
  if (read$maxindeg) {
    value <- cbind(value, maxindeg)
    cn <- c(cn, "Max. in-deg.")
  }
  if (read$maxoutdeg) {
    value <- cbind(value, maxoutdeg)
    cn <- c(cn, "Max. out-deg.")
  }
  if (read$mindeg) {
    value <- cbind(value, mindeg)
    cn <- c(cn, "Min. deg.")
  }
  if (read$minindeg) {
    value <- cbind(value, minindeg)
    cn <- c(cn, "Min. in-deg.")
  }
  if (read$minoutdeg) {
    value <- cbind(value, minoutdeg)
    cn <- c(cn, "Min. out-deg.")
  }

  value <- t(value)
  rownames(value) <- cn
  colnames(value) <- gnos
  .tkigraph.showData(value, title="Graphs properties", sort.button=FALSE)  
}

.tkigraph.plot <- function(simple=TRUE, gnos=NULL, ...) {

  if (is.null(gnos)) {
    gnos <- .tkigraph.get.selected()
  }
  graphs <- get("graphs", .tkigraph.env)

  if (length(gnos)==0) {
    return (.tkigraph.error("Please select one or more graphs to draw."))
  }

  max.vcount <- max(sapply(graphs[gnos], vcount))
  if (max.vcount > 5000) {
    vertex.size <- 1
  } else if (max.vcount > 30) {
    vertex.size <- 3
  } else {
    vertex.size <- 15
  }

  if (!simple) {
    read <- .tkigraph.dialogbox(TITLE="Drawing graphs",
                                interactive=list(name="Interactive",
                                  type="boolean", default="FALSE"),
                                vertex.size=list(name="Vertex size",
                                  type="numeric", default=vertex.size),
                                labels=list(name="Vertex labels", type="listbox",
                                  default="3", values=c("None", "IDs", "Names",
                                                 "Labels")),
                                elabels=list(name="Edge labels", type="listbox",
                                  default="0", values=c("None", "IDs", "Names",
                                                 "Values")),
                                layout=list(name="Layout",
                                  type="listbox", default="0",
                                  values=c("Default", "Force-based (KK)",
                                    "Force-based (FR)", "Tree (RT)",
                                    "Circle", "Random")))
  } else {
    read <- list(interactive=FALSE,
                 vertex.size=vertex.size,
                 labels=3,              # labels
                 elabels=0,             # none
                 layout=0)
  }
    
  if (!read$interactive) {
    fun <- function(...) { x11() ; plot.igraph(...) }
  } else {
    fun <- tkplot
  }
  
  layout.default <- function(graph, layout.par) {
    if ("x" %in% list.vertex.attributes(graph) &&
        "y" %in% list.vertex.attributes(graph)) {
      cbind( V(graph)$x , V(graph)$y )
    } else if ("layout" %in% list.graph.attributes(graph)) {
      l <- get.graph.attribute(graph, "layout")
      if (is.function(l)) {
        l(graph)
      } else {
        l
      }
    } else if (vcount(graph) < 300 && is.connected(graph)) {
      layout.kamada.kawai(graph)
    } else if (vcount(graph) < 1000) {
      layout.fruchterman.reingold(graph)
    } else {
      layout.circle(graph)
    }
  }
  
  layouts <- list(layout.default, layout.kamada.kawai,
                  layout.fruchterman.reingold,
                  layout.reingold.tilford, layout.circle, layout.random)

  if (read$vertex.size < 10) {
    label.dist <- 0.4
  } else {
    label.dist <- 0
  }
  
  for (i in gnos) {

    if (read$labels == "0") {
      labels <- NA
    } else if (read$labels == "1") {
      labels <- seq(vcount(graphs[[i]]))
    } else if (read$labels == "2") {
      labels <- V(graphs[[i]])$name
    } else if (read$labels == "3") {
      if ("label" %in% list.vertex.attributes(graphs[[i]])) {
        labels <- V(graphs[[i]])$label
      } else {
        labels <- V(graphs[[i]])$name
      }
    }

    if (read$elabels == "0") {
      elabels <- NA
    } else if (read$labels == "1") {
      elabels <- seq(ecount(graphs[[i]]))
    } else if (read$labels == "2") {
      elabels <- E(graphs[[i]])$name
    } else if (read$labels == "3") {
      if ("weight" %in% list.edge.attributes(graphs[[i]])) {
        elabels <- E(graphs[[i]])$weight
      } else {
        .tkigraph.warning("No edge weights, not a valued graph");
        elabels <- NA
      }
    }    

    if (vcount(graphs[[i]]) > 10) {
      eas <- 0.5
    } else {
      eas <- 1
    }
    
    g <- graphs[[i]]
    g <- remove.vertex.attribute(g, "name")
    fun(g, layout=layouts[[ read$layout+1 ]],
        vertex.size=read$vertex.size, ## vertex.color=read$vertex.color,
        vertex.label=labels, vertex.label.dist=label.dist,
        edge.label=elabels, edge.arrow.size=eas,
        ...)
  }
}

.tkigraph.by.hand <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos) > 1) {
    .tkigraph.error("Please select zero or one graph")
    return()
  }

  if (length(gnos)==0) {
    newdf <- edit(data.frame(list(from=character(), to=character())))
    if (ncol(newdf) > 2) {
      colnames(newdf) <- c("from", "to", "weight")
    }
    read <- .tkigraph.dialogbox(TITLE="Creating a graph by hand",
                                directed=list(name="Directed", type="boolean",
                                  default="FALSE"))
    g <- graph.data.frame(newdf, directed=read$directed)
    g <- set.graph.attribute(g, "name", "New graph")
    .tkigraph.add.graph(g)
  } else {
    graphs <- get("graphs", .tkigraph.env)
    df <- get.edgelist(graphs[[gnos]])
    colnames <- c("from", "to")    
    if ("weight" %in% list.edge.attributes(graphs[[gnos]])) {
      df <- cbind(df, E(g)$weight)
      colnames <- c("from", "to", "weight")
    }
    df <- as.data.frame(df)
        
    colnames(df) <- colnames
    df <- edit(df)
    if (ncol(df) > 2) {
      colnames(df) <- c("from", "to", "weight")
    }
    graphs[[gnos]] <- graph.data.frame(df,
                                       directed=is.directed(graphs[[gnos]]))
    assign("graphs", graphs, .tkigraph.env)
  }
  invisible(NULL)
}

.tkigraph.tree <- function() {
  read <- .tkigraph.dialogbox(TITLE="Regular tree",
                              n=list(name="Vertices", type="numeric",
                                default=63, min=0),
                              b=list(name="Branches", type="numeric",
                                default=2, min=1),
                              mode=list(name="Mode", type="listbox",
                                values=c("Directed (out)", "Directed (in)",
                                  "Undirected"), default="2"))
  read$mode <- c("out", "in", "undirected")[read$mode]
  g <- graph.tree(n=read$n, children=read$b, mode=read$mode)
  lay <- layout.reingold.tilford(g, root=0, mode="all")
  g <- set.graph.attribute(g, "layout", lay)
  g <- set.graph.attribute(g, "name", "Regular tree")
  .tkigraph.add.graph(g)
}

.tkigraph.ring <- function() {
  read <- .tkigraph.dialogbox(TITLE="Regular ring",
                              n=list(name="Vertices", type="numeric",
                                default=100, min=0))
  g <- graph.ring(n=read$n)
  g <- set.graph.attribute(g, "layout", layout.circle)
  g <- set.graph.attribute(g, "name", "Regular ring")
  .tkigraph.add.graph(g)
}

.tkigraph.lattice <- function() {
  read <- .tkigraph.dialogbox(TITLE="Regular lattice",
                              dim=list(name="Dimensions", type="numeric",
                                default=2, min=1, max=5),
                              s1=list(name="Size 1", type="numeric",
                                default=10, min=1),
                              s2=list(name="Size 2", type="numeric",
                                default=10, min=1),
                              s3=list(name="Size 3", type="numeric",
                                default=10, min=1),
                              s4=list(name="Size 4", type="numeric",
                                default=10, min=1),
                              s5=list(name="Size 5", type="numeric",
                                default=10, min=1))
  if (read$dim > 5) { read$dim <- 5 }
  dimv <- c(read$s1, read$s2, read$s3, read$s4, read$s5)[1:read$dim]
  g <- graph.lattice(dimvector=dimv)
  g <- set.graph.attribute(g, "name", "Regular Lattice")
  .tkigraph.add.graph(g)
}

.tkigraph.star <- function() {
  read <- .tkigraph.dialogbox(TITLE="Star graph",
                              n=list(name="Vertices", type="numeric",
                                default=100, min=0),
                              mode=list(name="Mode", type="listbox",
                                values=c("Directed (out)", "Directed (in)",
                                  "Undirected"), default="2"))
  read$mode <- c("out", "in", "undirected")[read$mode]
  g <- graph.star(read$n, mode=read$mode)
  g <- set.graph.attribute(g, "name", "Star graph")
  .tkigraph.add.graph(g)
}

.tkigraph.full <- function() {
  read <- .tkigraph.dialogbox(TITLE="Full graph",
                              n=list(name="Vertices", type="numeric",
                                default=30, min=0),
                              directed=list(name="Directed", type="boolean",
                                default="FALSE"),
                              loops=list(name="Loops", type="boolean",
                                default="FALSE"))
  g <- graph.full(read$n, read$directed, read$loops)
  g <- set.graph.attribute(g, "name", "Full graph")
  .tkigraph.add.graph(g)
}                             

.tkigraph.atlas <- function() {
  read <- .tkigraph.dialogbox(TITLE="Graph Atlas",
                              n=list(name="Number", type="numeric",
                                default=sample(0:1252, 1), min=0, max=1252))
  g <- graph.atlas(read$n)
  g <- set.graph.attribute(g, "name",
                           paste("Graph Atlas #", read$n))
  .tkigraph.add.graph(g)
}
  
.tkigraph.erdos.renyi.game <- function() {
  read <- .tkigraph.dialogbox(TITLE="Erdos-Renyi random graph, G(n,p)",
                              n=list(name="Vertices", type="numeric",
                                default=100, min=0),
                              p=list(name="Connection probability",
                                type="numeric", default=0.02, min=0, max=1),
                              directed=list(name="Directed",
                                type="boolean", default="FALSE"))
  
  g <- erdos.renyi.game(read$n,read$p,directed=read$directed)
  g <- set.graph.attribute(g, "name", "Random graph (Erdos-Renyi G(n,p))")
  .tkigraph.add.graph(g)
}

.tkigraph.erdos.renyi.gnm.game <- function() {
  read <- .tkigraph.dialogbox(TITLE="Erdos-Renyi random graph, G(n,m)",
                              n=list(name="Vertices", type="numeric",
                                default=100, min=0),
                              m=list(name="Edges", type="numeric",
                                default=200, min=0),
                              directed=list(name="Directed",
                                type="boolean", default="FALSE"))

  g <- erdos.renyi.game(read$n, read$m, type="gnm", directed=read$directed)
  g <- set.graph.attribute(g, "name", "Random graph (Erdos-Renyi G(n,m))")
  .tkigraph.add.graph(g)
}

.tkigraph.barabasi.game <- function() {
  read <- .tkigraph.dialogbox(TITLE="Scale Free graph",
                              n=list(name="Vertices", type="numeric",
                                default=100, min=0),
                              m=list(name="Edges per time step",
                                type="numeric", default=1, min=0),
                              directed=list(name="Directed",
                                type="boolean", default="TRUE"))
  g <- barabasi.game(n=read$n, m=read$m, directed=read$directed)
  g <- set.graph.attribute(g, "name", "Scale-free random graph")
  .tkigraph.add.graph(g)
}

.tkigraph.degree.sequence.game <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos) == 0) {
    .tkigraph.error("Please select at least one graph")
    return()
  }
  graphs <- get("graphs", .tkigraph.env)

  for (i in gnos) {
    if (is.directed(graphs[[i]])) {
      indeg <- degree(graphs[[i]], mode="in")
      outdeg <- degree(graphs[[i]], mode="out")
      g <- degree.sequence.game(out.deg=outdeg, in.deg=indeg)
    } else {
      deg <- degree(graphs[[i]])
      g <- degree.sequence.game(deg)
    }
    g <- set.graph.attribute(g, "name",
                             paste(sep="", "Configuration model (#", i,")"))
    .tkigraph.add.graph(g)
  }
}

.tkigraph.watts.strogatz <- function() {
  read <- .tkigraph.dialogbox(TITLE="Watts-Strogatz graph",
                              dim=list(name="Dimensions", type="numeric",
                                default=1, min=1),
                              size=list(name="Lattice size", type="numeric",
                                default=1000, min=1),
                              nei=list(name="Neighborhood", type="numeric",
                                default=5, min=1),
                              p=list(name="Rewiring probability",
                                type="numeric", default=0.01, min=0, max=1))
  g <- watts.strogatz.game(dim=read$dim, size=read$size, nei=read$nei,
                           p=read$p)
  g <- set.graph.attribute(g, "name", "Watts-Strogatz small-world graph")
  if (read$dim == 1) { 
    g <- set.graph.attribute(g, "layout", layout.circle)
  }
  .tkigraph.add.graph(g)
}

.tkigraph.simplify <- function() {
  gnos <- .tkigraph.get.selected()
  graphs <- get("graphs", .tkigraph.env)

  for (i in gnos) {
    g <- simplify(graphs[[i]])
    g <- set.graph.attribute(g, "name",
                             paste(sep="", "Simplification of #", i))
    .tkigraph.add.graph(g)
  }
}

#####################################################

.tkigraph.degree <- function(mode) {
  gnos <- .tkigraph.get.selected()
  if (length(gnos)!=1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graphs <- get("graphs", .tkigraph.env)
  deg <- degree(graphs[[gnos]], mode=mode)
  value <- data.frame(Vertex=V(graphs[[gnos]])$name, deg)
  colnames(value) <- c("Vertex", paste(sep="","Degree (", mode, ")"))
  plot.command <- function() {
    read <- .tkigraph.dialogbox(TITLE="Plot degree distribution",
                                logx=list(name="Logarithmic `X' axis",
                                  type="boolean", default="FALSE"),
                                logy=list(name="Logarithmic `Y' axis",
                                  type="boolean", default="FALSE"),
                                hist=list(name="Histogram",
                                  type="boolean", default="FALSE"))

    if (!read$hist) {
      h <- hist(value[,2], -1:max(value[,2]), plot=FALSE)$density  
      log <- ""
      if (read$logx) { log <- paste(sep="", log, "x") }
      if (read$logy) { log <- paste(sep="", log, "y") }
      x11()
      plot(0:max(value[,2]), h, xlab="Degree", ylab="Relative frequency",
           type="b", main="Degree distribution", log=log)
    } else {
      x11()
      hist(value[,2], main="Degree distribution", xlab="Degree")
    }
  }
  value <- value[ order(value[,2], decreasing=TRUE), ]
  mv <- paste("Mean degree:", round(mean(deg), 2))
  .tkigraph.showData(value,
                     title=paste(sep="", "Degree for graph #", gnos),
                     plot.text="Plot distribution",
                     plot.command=plot.command,
                     showmean=mv)
}

.tkigraph.degree.dist <- function(power=FALSE) {
  gnos <- .tkigraph.get.selected()
  if (length(gnos)!=1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graphs <- get("graphs", .tkigraph.env)
  read <- .tkigraph.dialogbox(TITLE="Choose degree type",
                              type=list(name="Degree type",
                                type="listbox", default="0",
                                values=c("Out", "In", "Total")))
  mode <- c("out", "in", "all")[read$type+1]
  deg <- degree(graphs[[gnos]], mode=mode)
  x11()
  h <- hist(deg, -1:max(deg), plot=FALSE)$density
  plot(0:max(deg), h, xlab="Degree", ylab="Relative frequency",
       type="b", main="Degree distribution", log="xy")

  if (power) {
    if (max(deg)<10) {
      .tkigraph.error("Degrees are too small for a power-law fit")
      return()
    }
    fit <- power.law.fit(deg, xmin=10)
    lines(0:max(deg), (0:max(deg))^(-coef(fit)), col="red")
    legend("topright", c(paste("exponent:", round(coef(fit), 2)),
                         paste("standard error:", round(sqrt(vcov(fit)), 2))),
           bty="n", cex=1.5)
  }
}

.tkigraph.closeness <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos)!=1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graphs <- get("graphs", .tkigraph.env)
  cl <- closeness(graphs[[gnos]], mode="out")
  value <- data.frame(Vertex=V(graphs[[gnos]])$name, Closeness=cl)
  value <- value[ order(value[,2], decreasing=TRUE), ]
  mv <- paste("Mean value:", round(mean(cl),2))
  .tkigraph.showData(value,
                     title=paste(sep="", "Closeness for graph #", gnos),
                     showmean=mv)
}

.tkigraph.betweenness <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos)!=1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graphs <- get("graphs", .tkigraph.env)
  btw <- betweenness(graphs[[gnos]])
  vc <- vcount(graphs[[gnos]])
  m <- (vc-1)*(vc-2)
  nbtw <- btw/m
  value <- data.frame(V(graphs[[gnos]])$name, btw, nbtw)
  colnames(value) <- c("Vertex", "Betweenness", "Normalized Betweenness")
  value <- value[ order(value[,2], decreasing=TRUE), ]
  mv <- paste("Mean value:", round(mean(btw),2), "&", round(mean(nbtw),5))
  .tkigraph.showData(value,
                     title=paste(sep="", "Betweenness for graph #", gnos),
                     showmean=mv)
}

.tkigraph.constraints <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos)!=1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graphs <- get("graphs", .tkigraph.env)
  const <- constraint(graphs[[gnos]])
  value <- data.frame(V(graphs[[gnos]])$name, const)
  colnames(value) <- c("Vertex", "Constraint")
  value <- value[ order(value[,2], decreasing=TRUE), ]
  mv <- paste("Mean value:", round(mean(const),2))
  .tkigraph.showData(value,
                     title=paste(sep="", "Constraint for graph #", gnos),
                     showmean=mv)
}  

.tkigraph.power.centrality <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos)!=1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graphs <- get("graphs", .tkigraph.env)
  bp <- bonpow(graphs[[gnos]])
  value <- data.frame(V(graphs[[gnos]])$name, bp)
  colnames(value) <- c("Vertex", "Power centrality")
  value <- value[ order(value[,2], decreasing=TRUE), ]
  mv <- paste("Mean value:", round(mean(bp),2))
  .tkigraph.showData(value,
                     title=paste(sep="", "Power centrality for graph #", gnos),
                     showmean=mv)
}    

.tkigraph.page.rank <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos)!=1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graphs <- get("graphs", .tkigraph.env)
  bp <- page.rank(graphs[[gnos]])$vector
  value <- data.frame(V(graphs[[gnos]])$name, bp)
  colnames(value) <- c("Vertex", "Page rank")
  value <- value[ order(value[,2], decreasing=TRUE), ]
  mv <- paste("Mean value:", round(mean(bp),2))
  .tkigraph.showData(value,
                     title=paste(sep="", "Page rank centrality for graph #", gnos),
                     showmean=mv)
}    

.tkigraph.edge.betweenness <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos)!=1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graphs <- get("graphs", .tkigraph.env)
  ebtw <- edge.betweenness(graphs[[gnos]])
  el <- get.edgelist(graphs[[gnos]])
  value <- data.frame(E(graphs[[gnos]])$name, el[,1], el[,2], ebtw)
  colnames(value) <- c("Edge", "From", "To", "Betweenness")
  value <- value[ order(value[,4], decreasing=TRUE), ]
  mv <- paste("Mean value:", round(mean(ebtw),2))
  .tkigraph.showData(value,
                     title=paste(sep="", "Edge betweenness for graph #",gnos),
                     showmean=mv)
}

#####################################################

.tkigraph.dist.matrix <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos)!=1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }

  graph <- get("graphs", .tkigraph.env)[[gnos]]
  if (vcount(graph) > 100) {
    .tkigraph.error("Graphs is too large to do this")
    return()
  }

  value <- shortest.paths(graph, mode="out")
  rownames(value) <- colnames(value) <- V(graph)$name
  .tkigraph.showData(value, sort.button=FALSE, 
                     title=paste(sep="", "Distance matrix for graph #", gnos))
}

.tkigraph.distance.tofrom <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos)!=1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }

  graph <- get("graphs", .tkigraph.env)[[gnos]]
  read <- .tkigraph.dialogbox(TITLE="Distance from a vertex",
                              v=list(name="Vertex ID", type="numeric",
                                default=1, min=1, max=vcount(graph)))
  if (read$v < 1 || read$v > vcount(graph)) {
    .tkigraph.error("Invalid vertex ID")
    return()
  }

  value <- shortest.paths(graph, read$v-1, mode="out")
  dim(value) <- NULL
  value <- data.frame( V(graph)$name, value)
  colnames(value) <- c("Vertex", "Distance")
  mv <- paste("Mean distance:", round(mean(value),2))
  .tkigraph.showData(value,
                     title=paste("Distance from vertex", read$v, "in graph #",
                       gnos), showmean=mv)
}

.tkigraph.diameter <- function(mode="dia") {
  gnos <- .tkigraph.get.selected()
  if (length(gnos)==0) {
    .tkigraph.error("Please select one or more graphs")
    return()
  }

  isconn <- logical()
  dia <- numeric()
  graphs <- get("graphs", .tkigraph.env)
  for (i in seq(along=gnos)) {
    if (mode=="dia") {
      dia[i] <- diameter(graphs[[ gnos[i] ]])
    } else if (mode=="path") {
      dia[i] <- average.path.length(graphs[[ gnos[i] ]], directed=FALSE)
    }
    isconn[i] <- is.connected(graphs[[ gnos[i] ]])
  }

  value <- data.frame( gnos, isconn, dia)
  if (mode=="dia") {
    title <- "Diameter"
    colnames(value) <- c("Graph #", "Connected", "Diameter")
  } else if (mode=="path") {
    title <- "Average path length"
    colnames(value) <- c("Graph #", "Connected", "Mean path length")  }
  title <- paste(title, "of graph")
  if (length(gnos) > 1) {
    title <- paste(sep="", title, "s")
  }
  .tkigraph.showData(value, title=title)
  
}

.tkigraph.plot.diameter <- function(simple=FALSE) {
  gnos <- .tkigraph.get.selected()
  if (length(gnos)!=1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }

  graph <- get("graphs", .tkigraph.env)[[gnos]]
  edges <- E(graph, path=get.diameter(graph, directed=FALSE), directed=FALSE)
  color <- rep("black", ecount(graph))
  color[edges+1] <- "red"
  width <- rep(1, ecount(graph))
  width[edges+1] <- 2
  .tkigraph.plot(gnos=gnos, simple=simple, edge.color=color, edge.width=width)
}

.tkigraph.clusters <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos) != 1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graph <- get("graphs", .tkigraph.env)[[gnos]]
  comm <- clusters(graph)
  members <- sapply(sapply(seq(along=comm$csize)-1,
                           function(i) which(comm$membership==i)),
                    paste, collapse=", ")
  value <- data.frame("Component"=seq(along=comm$csize), "Members"=members)
  .tkigraph.showData(value, title=paste("Components of graph #",
                       gnos), right=FALSE)
}

.tkigraph.clusters.membership <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos) != 1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graph <- get("graphs", .tkigraph.env)[[gnos]]
  comm <- clusters(graph)
  value <- data.frame("Vertex"=seq(along=comm$membership),
                 "Component"=comm$membership+1)
  .tkigraph.showData(value, title=paste("Components of graph #", gnos))
}

.tkigraph.calculate.clusters <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos) != 1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graph <- get("graphs", .tkigraph.env)[[gnos]]
  cs <- clusters(graph)$csize
  value <- data.frame(seq(along=cs), cs)
  colnames(value) <- c("Cluster #", "Size")

  plot.command <- function() {
    read <- .tkigraph.dialogbox(TITLE="Plot degree distribution",
                                logx=list(name="Logarithmic `X' axis",
                                  type="boolean", default="FALSE"),
                                logy=list(name="Logarithmic `Y' axis",
                                  type="boolean", default="FALSE"),
                                hist=list(name="Histogram",
                                  type="boolean", default="FALSE"))

    if (!read$hist) {
      h <- hist(value[,2], 0:max(value[,2]), plot=FALSE)$density  
      log <- ""
      if (read$logx) { log <- paste(sep="", log, "x") }
      if (read$logy) { log <- paste(sep="", log, "y") }
      x11()
      plot(1:max(value[,2]), h, xlab="Component size",
           ylab="Relative frequency",
           type="b", main="Component size distribution", log=log)
    } else {
      x11()
      hist(value[,2], main="Component size distribution", xlab="Degree")
    }
  }
  value <- value[ order(value[,2], decreasing=TRUE), ]
  mv <- paste("Mean component size:", round(mean(cs),2))
  .tkigraph.showData(value,
                     title=paste(sep="", "Component sizes, graph #", gnos),
                     plot.text="Plot distribution",
                     plot.command=plot.command, showmean=mv)
}

.tkigraph.plot.comp <- function(simple=FALSE) {
  gnos <- .tkigraph.get.selected()
  if (length(gnos) != 1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graph <- get("graphs", .tkigraph.env)[[gnos]]
  clu <- clusters(graph)
  colbar <- rainbow(length(clu$csize)*2)
  vertex.color <- colbar[ clu$membership+1 ]
  .tkigraph.plot(gnos=gnos, simple=simple, vertex.color=vertex.color)
}

.tkigraph.create.giantcomp <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos) != 1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graph <- get("graphs", .tkigraph.env)[[gnos]]
  clu <- clusters(graph)
  v <- which(clu$membership+1 == which.max(clu$csize)) - 1
  g <- subgraph(graph, v)
  .tkigraph.add.graph(g)
}

.tkigraph.create.mycomp <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos) != 1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graph <- get("graphs", .tkigraph.env)[[gnos]]
  read <- .tkigraph.dialogbox(TITLE="Component of a vertex",
                              vertex=list(name="Vertex", type="numeric",
                                default=1, min=1, max=vcount(graph)))
  if (read$vertex<1 || read$vertex >vcount(graph)) {
    .tkigraph.error("Invalid vertex id")
    return()
  }

  g <- subgraph(graph, subcomponent(graph, read$vertex-1))
  .tkigraph.add.graph(g)
}

.tkigraph.create.comp <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos) != 1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graph <- get("graphs", .tkigraph.env)[[gnos]]
  read <- .tkigraph.dialogbox(TITLE="Graph from component",
                              comp=list(name="Component id", type="numeric",
                                default=1, min=1))
  clu <- clusters(graph)
  if (read$comp<1 || read$comp > length(clu$csize)) {
    .tkigraph.error("Invalid component id")
    return()
  }
  
  v <- which(clu$membership==read$comp-1)-1
  g <- subgraph(graph, v)
  .tkigraph.add.graph(g)  
}

.tkigraph.motifs.draw <- function() {
  read <- .tkigraph.dialogbox(TITLE="Draw all motifs",
                              size=list(name="Size", type="numeric",
                                default=3, min=3, max=4),
                              dir=list(name="Directed", type="boolean",
                                default="FALSE"))

  if (read$size < 3 || read$size > 4) {
    .tkigraph.error("Invalid motif size, should be 3 or 4")
    return()
  }

  if (read$size == 3) {
    co <- matrix( c(1,1, 0,0, 2,0), ncol=2, byrow=TRUE)
  } else {
    co <- matrix( c(0,1, 1,1, 0,0, 1,0), ncol=2, byrow=TRUE)
  }

  if (read$size == 3 && read$dir) {
    no <- 16
    rows <- cols <- 4
  } else if (read$size == 3 && !read$dir) {
    no <- 4
    rows <- cols <- 2
  } else if (read$size == 4 && read$dir) {
    no <- 216
    rows <- cols <- 15
  } else if (read$size == 4 && !read$dir) {
    no <- 11
    rows <- 4
    cols <- 3
  }
  names <- as.character(seq(no)-1)
  x11()
  layout( matrix(1:(rows*cols), nrow=rows, byrow=TRUE) )
  layout.show(rows*cols)
  for (i in seq(no)-1) {
    g <- graph.isocreate(read$size, i, directed=read$dir)
    par(mai=c(0,0,0,0), mar=c(0,0,0,0))
    par(cex=2)
    plot(g, layout=co, vertex.color="red", vertex.label=NA, frame=TRUE,
         edge.color="black", margin=0.1)
    text(0,0, names[i+1], col="blue")
  }
  
}

.tkigraph.motifs.find <- function() {

  gnos <- .tkigraph.get.selected()
  if (length(gnos)!=1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  
  read <- .tkigraph.dialogbox(TITLE="Find motifs",
                              size=list(name="Size", type="numeric",
                                default=3, min=3, max=4))
  if (read$size < 3 || read$size > 4) {
    .tkigraph.error("Invalid motif size, should be 3 or 4")
    return()
  }
  
  graphs <- get("graphs", .tkigraph.env)
  motifs <- graph.motifs(graphs[[gnos]], size=read$size)

  x11()
  barplot(motifs)

  if (read$size == 3) {
    co <- matrix( c(1,1, 0,0, 2,0), ncol=2, byrow=TRUE)
  } else {
    co <- matrix( c(0,1, 1,1, 0,0, 1,0), ncol=2, byrow=TRUE)
  }

  if (read$size == 3 && is.directed(graphs[[gnos]])) {
    no <- 16
    rows <- cols <- 4
  } else if (read$size == 3 && !is.directed(graphs[[gnos]])) {
    no <- 4
    rows <- cols <- 2
  } else if (read$size == 4 && is.directed(graphs[[gnos]])) {
    no <- 216
    rows <- cols <- 15
  } else if (read$size == 4 && !is.directed(graphs[[gnos]])) {
    no <- 11
    rows <- 4
    cols <- 3
  }
  names <- as.character(seq(no)-1)
  x11()
  layout( matrix(1:(rows*cols), nrow=rows, byrow=TRUE) )
  layout.show(rows*cols)
  for (i in seq(no)-1) {
    g <- graph.isocreate(read$size, i, directed=is.directed(graphs[[gnos]]))
    par(mai=c(0,0,0,0), mar=c(0,0,0,0))
    par(cex=2)
    plot(g, layout=co, vertex.color="red", vertex.label=NA, frame=TRUE,
         edge.color="black", margin=0.1)
    text(0,0, motifs[i+1], col="green")
  }  
  
}

.tkigraph.spinglass <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos)!=1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graph <- get("graphs", .tkigraph.env)[[gnos]]

  if (!is.connected(graph)) {
    .tkigraph.error("Graph is not connected")
    return()
  }

  weights <- if ("weight" %in% list.edge.attributes(graph)) "TRUE" else "FALSE"
  read <- .tkigraph.dialogbox(TITLE="Spinglass community structure",
                              gamma=list(name="Gamma parameter",
                                type="numeric", default=1),
                              weights=list(name="Use edge weights",
                                type="boolean", default=weights),
                              spins=list(name="Number of spins",
                                type="numeric", default=25),
                              parupdate=list(name="Parallel update",
                                type="boolean", default="FALSE"),
                              update.rule=list(name="Update rule",
                                type="listbox", default="1",
                                values=c("Simple", "Configuration model")),
                              start.temp=list(name="Start temperature",
                                type="numeric", default=1),
                              stop.temp=list(name="Stop temperature",
                                type="numeric", default=0.1),
                              cool.fact=list(name="Cooling factor",
                                type="numeric", default=0.99))

  read$update.rule <- c("simple", "config")[read$update.rule+1]
  if (read$weights) {
    if (!"weight" %in% list.edge.attributes(graph)) {
      .tkigraph.warning("This graphs is not weighted")
      read$weights <- NULL
    } else {
      read$weights <- E(graph)$weight
    }
  } else {
    read$weights <- NULL
  }
  comm <- spinglass.community(graph, weights=read$weights, spins=read$spins,
                              parupdate=read$parupdate, start.temp=read$start.temp,
                              stop.temp=read$stop.temp, cool.fact=read$cool.fact,
                              update.rule=read$update.rule, gamma=read$gamma)

  .tkigraph.spinglass.community.dialog(comm, read, gnos)
}

.tkigraph.spinglass.community.dialog <- function(comm, read, gnos) {
  dialog <- tktoplevel()
  frame <- tkframe(dialog)
  tkgrid(frame)
  tktitle(dialog) <- "Spinglass community structure algorithm results"

  read$update.rule <- if (read$update.rule=="simple") "Simple" else "Configuration model"
  tkgrid(tklabel(dialog, text="Spinglass community structure algorithm results",
                 font=tkfont.create(family="times", size=16, weight="bold")),
         columnspan=3, sticky="nsew", "in"=frame, padx=10, pady=10)
  tkgrid(txt <- tktext(dialog), columnspan=1, rowspan=5, sticky="nsew",
         "in"=frame, padx=10, pady=10)
  tkconfigure(txt, height=15)
  tkinsert(txt, "end", "Parameters were:\n")
  tkinsert(txt, "end", paste("  Gamma=", read$gamma, "\n"))
  tkinsert(txt, "end", if (is.null(read$weights)) "  Weights were not used.\n" else
           "  Weights were used.\n")
  tkinsert(txt, "end", paste("  Number of spins=", read$spins, "\n"))
  tkinsert(txt, "end", if (read$parupdate) "  Parallel updating.\n" else
           "  Sequential updating.\n")
  tkinsert(txt, "end", paste("  Update rule:", read$update.rule, "\n"))
  tkinsert(txt, "end", paste("  Start temperature was", read$start.temp, "\n"))
  tkinsert(txt, "end", paste("  Stop temperaure was", read$stop.temp, "\n"))
  tkinsert(txt, "end", paste("  Cooling factor was", read$cool.fact, "\n"))

  tkinsert(txt, "end", "\nResults:\n")
  tkinsert(txt, "end", paste("  Number of communities found:", length(comm$csize),
                             "\n"))
  tkinsert(txt, "end", paste("  Modularity of the result:", comm$modularity, "\n"))
  tkinsert(txt, "end", paste("  Stopped at temperature:", comm$temperature, "\n"))
  tkconfigure(txt, state="disabled")

  show.communities <- function() {
    members <- sapply(sapply(seq(along=comm$csize)-1,
                             function(i) which(comm$membership==i)),
                      paste, collapse=", ")
    value <- data.frame("Community"=seq(along=comm$csize), "Members"=members)
    .tkigraph.showData(value,
                       title=paste("Communities, spinglass algorithm on graph #",
                         gnos), right=FALSE)
  }
  show.membership <- function() {
    value <- data.frame("Vertex"=seq(along=comm$membership),
                   "Community"=comm$membership+1)
    .tkigraph.showData(value,
                       title=paste("Communities, spinglass algorithm on graph #",
                         gnos))
  }
  show.csize <- function() {
    value <- data.frame("Comm. #"=seq(along=comm$csize), "Size"=comm$csize)
    value <- value[ order(value[,2], decreasing=TRUE), ]
    .tkigraph.showData(value,
                       title=paste("Communities, spinglass algorithm on graph #",
                         gnos))                       
  }
  plot.communities <- function(simple=FALSE) {
    colbar <- rainbow(length(comm$csize)*2)
    vertex.color=colbar[ comm$membership+1 ]
    .tkigraph.plot(gnos=gnos, simple=simple, vertex.color=vertex.color)
  }
  create.subgraph <- function() {
    ## TODO
  }
  
  tkgrid(tkbutton(dialog, text="Show communities", command=show.communities),
         "in"=frame, sticky="ew", column=1, row=1, padx=10, pady=10)
  tkgrid(tkbutton(dialog, text="Show membership", command=show.membership),
         "in"=frame, sticky="ew", column=1, row=2, padx=10, pady=10)
  tkgrid(tkbutton(dialog, text="Show community sizes", command=show.csize),
         "in"=frame, sticky="ew", column=1, row=3, padx=10, pady=10)
  tkgrid(tkbutton(dialog, text="Draw communities",
                  command=function() plot.communities(simple=FALSE)),
         "in"=frame, sticky="ew", column=1, row=4, padx=10, pady=10)
##   tkgrid(tkbutton(dialog, text="Create subgraph", command=create.subgraph),
##          "in"=frame, sticky="nsew", column=1, row=6, padx=10, pady=10)
  
  tkgrid(tkbutton(dialog, text="Close", command=function() tkdestroy(dialog)),
         "in"=frame, sticky="nsew", columnspan=2, padx=10, pady=10)
}

.tkigraph.my.spinglass <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos)!=1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graph <- get("graphs", .tkigraph.env)[[gnos]]

  if (!is.connected(graph)) {
    .tkigraph.error("Graph is not connected")
    return()
  }

  weights <- if ("weight" %in% list.edge.attributes(graph)) "TRUE" else "FALSE"
  read <- .tkigraph.dialogbox(TITLE="Spinglass community of a vertex",
                              vertex=list(name="Vertex", type="numeric",
                                default=1, min=1, max=vcount(graph)),
                              gamma=list(name="Gamma parameter",
                                type="numeric", default=1),
                              weights=list(name="Use edge weights",
                                type="boolean", default=weights),
                              spins=list(name="Number of spins",
                                type="numeric", default=25),
                              update.rule=list(name="Update rule",
                                type="listbox", default="1",
                                values=c("Simple", "Configuration model")))
                              
  if (read$vertex<1 || read$vertex > vcount(graph)) {
    .tkigraph.error("Invalid vertex id")
    return()
  }
  
  read$update.rule <- c("simple", "config")[read$update.rule+1]
  if (read$weights) {
    if (!"weight" %in% list.edge.attributes(graph)) {
      .tkigraph.warning("This graphs is not weighted")
      read$weights <- NULL
    } else {
      read$weights <- E(graph)$weight
    }
  } else {
    read$weights <- NULL
  }
  comm <- spinglass.community(graph, vertex=read$vertex,
                              weights=read$weights, spins=read$spins,
                              update.rule=read$update.rule, gamma=read$gamma)
  .tkigraph.spinglass.mycommunity.dialog(comm, read, gnos)
}

.tkigraph.spinglass.mycommunity.dialog <- function(comm, read, gnos) {
  dialog <- tktoplevel()
  frame <- tkframe(dialog)
  tkgrid(frame)
  tktitle(dialog) <- "Spinglass community of a single vertex"

  scr <- tkscrollbar(dialog, repeatinterval=5,
                     command=function(...) tkyview(txt,...))
  
  read$update.rule <- if (read$update.rule=="simple") "Simple" else "Configuration model"
  tkgrid(tklabel(dialog, text="Spinglass community of a single vertex",
                 font=tkfont.create(family="times", size=16, weight="bold")),
         columnspan=3, sticky="nsew", "in"=frame, padx=10, pady=10)
  tkgrid(txt <- tktext(dialog, yscrollcommand=function(...) tkset(scr,...)),
         columnspan=1, rowspan=3, sticky="nsew", "in"=frame, padx=10, pady=10)
  tkconfigure(txt, height=17)
  tkgrid(scr, row=1, column=1, rowspan=3, sticky="ns", "in"=frame, pady=10)
  tkinsert(txt, "end", "Parameters were:\n")
  tkinsert(txt, "end", paste("  Vertex:", read$vertex, "\n"));
  tkinsert(txt, "end", paste("  Gamma=", read$gamma, "\n"))
  tkinsert(txt, "end", if (is.null(read$weights)) "  Weights were not used.\n" else
           "  Weights were used.\n")
  tkinsert(txt, "end", paste("  Number of spins=", read$spins, "\n"))
  tkinsert(txt, "end", paste("  Update rule:", read$update.rule, "\n"))

  tkinsert(txt, "end", "\nResults:\n")
  tkinsert(txt, "end", paste("  Size of the community:", length(comm$community),
                             "\n"))
  tkinsert(txt, "end", paste("  Cohesion:", comm$cohesion, "\n"))
  tkinsert(txt, "end", paste("  Adhesion:", comm$adhesion, "\n"))
  tkinsert(txt, "end", paste("  Inner links:", comm$inner.links, "\n"))
  tkinsert(txt, "end", paste("  Outer links:", comm$outer.links, "\n"))

  tkinsert(txt, "end", "\nThe community:\n")
  .tkigraph.ttt <- NULL
  con <- textConnection(".tkigraph.ttt", open="w")
  cat(sort(comm$community), file=con, fill=TRUE, sep=", ")
  close(con)
  tkinsert(txt, "end", .tkigraph.ttt)
  tkconfigure(txt, state="disabled")

  plot.communities <- function(simple=FALSE) {
    graph <- get("graphs", .tkigraph.env)[[gnos]]
    color <- rep("skyblue2", vcount(graph))
    color[ comm$community+1 ] <- "red"
    .tkigraph.plot(gnos=gnos, simple=simple, vertex.color=color)
  }
  
  create.graph <- function() {
    graph <- get("graphs", .tkigraph.env)[[gnos]]
    g <- subgraph(graph, comm$community)
    .tkigraph.add.graph(g)
  }
  
  tkgrid(tkbutton(dialog, text="Draw community",
                  command=function() plot.communities(simple=FALSE)),
         "in"=frame, sticky="ew", column=2, row=1, padx=10, pady=10)  
  tkgrid(tkbutton(dialog, text="Create graph from community",
                  command=create.graph),
         "in"=frame, sticky="ew", column=2, row=2, padx=10, pady=10)         
  
  tkgrid(tkbutton(dialog, text="Close", command=function() tkdestroy(dialog)),
         "in"=frame, sticky="nsew", columnspan=3, padx=10, pady=10)  
}

.tkigraph.cohesion <- function() {
  gnos <- .tkigraph.get.selected()
  if (length(gnos) != 1) {
    .tkigraph.error("Please select exactly one graph")
    return()
  }
  graphs <- decompose.graph(get("graphs", .tkigraph.env)[[gnos]])
  coh <- sapply(graphs, graph.cohesion)
  value <- data.frame("Component"=seq(length=length(graphs)), "Cohesion"=coh)
  .tkigraph.showData(value, title=paste("Cohesion of components in graph #",
                              gnos), right=FALSE)
}  

.tkigraph.help <- function(page="index.html") {
  dialog <- tktoplevel()
  tktitle(dialog) <- "Help (main page)"

  close <- function() {
    tkdestroy(dialog)
  }

  scr <- tkscrollbar(dialog, repeatinterval=5,
                     command=function(...) tkyview(txt,...))
  txt <- tktext(dialog, yscrollcommand=function(...) tkset(scr, ...),
                width=80, height=40)

  main.menu <- tkmenu(dialog)
  tkadd(main.menu, "command", label="Back", command=function() {
    tcl("render_back", txt)
  })
  tkadd(main.menu, "command", label="Forw", command=function() {
    tcl("render_forw", txt)
  })
  tkadd(main.menu, "command", label="Home", command=function() {
    tcl("render", txt, "index.html"); return()
  })
  tkadd(main.menu, "command", label="Close", command=function() {
    tkdestroy(dialog); return()
  })
  tkconfigure(dialog, "-menu", main.menu)
  
  tkpack(scr, side="right", fill="y", expand=0)
  tkpack(txt, side="left", fill="both", expand=1)

  browser.button <- tkbutton(dialog, command=function() {
    browseURL(tclvalue("browser_url"))
  })
  
  tcl("global", "tkigraph_help_root", "tkigraph_help_history",
      "tkigraph_help_history_pos", "browser_button", "browser_url")  
  tcl("set", "tkigraph_help_root",
      system.file("tkigraph_help", package="igraph0"))
  tcl("set", "browser_button", browser.button)
  
  tcl("source", system.file("html_library.tcl", package="igraph0"))
  tcl("source", system.file("my_html_library.tcl", package="igraph0"))
  tcl("HMinit_win", txt)
  tcl("start_history", txt)
  tcl("render", txt, "index.html")
  
  tkconfigure(txt, state="disabled")
}

.tkigraph.help.external <- function(page="index.html") {
  f <- system.file("tkigraph_help/index.html", package="igraph0")
  browseURL(f)
}

.tkigraph.about <- function() {
  dialog <- tktoplevel()
  tktitle(dialog) <- "About tkigraph"
  image <-tkimage.create("photo", "img", format="gif",
                         file=system.file("igraph.gif", package="igraph0"))
  logo <- tklabel(dialog, relief="flat", padx=10, pady=10, image=image)
  label <- tklabel(dialog, padx=30, pady=10,
                   text=paste(sep="", "tkigraph (c) 2009 Gabor Csardi\n",
                     "igraph (c) 2003-2009 Gabor Csardi and Tamas Nepusz\n\n",
                     "This is igraph version ",
                     packageDescription("igraph0")$Version, " and\n",
                     R.version$version.string))
  close <- tkbutton(dialog, text="Close", command=function() {
    tkdestroy(dialog); return()
  })

  tkpack(logo, side="top", anchor="c", expand=0)
  tkpack(label, side="top", anchor="c", expand=0)
  tkpack(close, side="bottom", anchor="c", expand=0)
}

#####################################################
# This is from the 'relimp' package by David Firth, thanks

.tkigraph.showData <-
    function (dataframe,
              colname.bgcolor = "grey50",
              rowname.bgcolor = "grey50",
              body.bgcolor = "white",
              colname.textcolor = "white",
              rowname.textcolor = "white",
              body.textcolor = "black",
              font = "Courier 12",
              maxheight = 30,
              maxwidth = 80,
              title = NULL,
              rowname.bar = "left",
              colname.bar = "top",
              rownumbers = FALSE,
              placement = "-20-40",
              plot.text="Plot",
              plot.command=NULL,
              suppress.X11.warnings = FALSE,
              right=TRUE,
              showmean=NULL,
              sort.button=TRUE,
              inthis=NULL)
{
    if (suppress.X11.warnings) { ## as in John Fox's Rcmdr package
        messages.connection <- textConnection(".messages", open = "w",
                                              local = TRUE)
        sink(messages.connection, type = "message")
        on.exit({
            sink(type="message")
            close(messages.connection)
        })
    }
    object.name <- deparse(substitute(dataframe))
    if (!is.data.frame(dataframe)){
        temp <- try(dataframe <- as.data.frame(dataframe), silent = FALSE)
        if (inherits(temp, "try-error")) {
            stop(paste(object.name, "cannot be coerced to a data frame"))
        }
        object.name <- paste("as.data.frame(", object.name, ")", sep = "")
    }

    if (is.numeric(rownumbers) &&
        length(rownumbers) != nrow(dataframe))
        stop("rownumbers argument must be TRUE, FALSE or have length nrow(dataframe)")
    oldwidth <- unlist(options("width"))
    options(width = 10000)
    conn <- textConnection(".tkigraph.ttt", open="w")
    sink(conn)
    options(max.print=10000000)
    print(dataframe, right=right)
    sink()
    close(conn)
    zz <- strsplit(.tkigraph.ttt, "\n", fixed=TRUE)
    if (length(zz) > 1 + nrow(dataframe)) stop(
       "data frame too wide")
    options(width = oldwidth)
    if (is.null(inthis)) {
      base <- tktoplevel()
      tkwm.geometry(base, placement)
      tkwm.title(base, {
        if (is.null(title))
          object.name
        else title
      })
    } else { base <- inthis }      
    nrows <- length(zz) - 1
    if (is.numeric(rownumbers))
        rowname.text <- paste(rownumbers, row.names(dataframe))
    else if (rownumbers)
        rowname.text <- paste(1:nrows, row.names(dataframe))
    else rowname.text <- row.names(dataframe)
    namewidth = max(nchar(rowname.text))
    yy <- substring(zz, 2 + max(nchar(row.names(dataframe))))
    datawidth <- max(nchar(yy))
    winwidth <- min(1 + datawidth, maxwidth)
    hdr <- tktext(base,
                  bg = colname.bgcolor,
                  fg = colname.textcolor,
                  font = font,
                  height = 1,
                  width = winwidth,
                  takefocus = TRUE)
    ftr <- tktext(base,
                  bg = colname.bgcolor,
                  fg = colname.textcolor,
                  font = font,
                  height = 1,
                  width = winwidth,
                  takefocus = TRUE)
    textheight <- min(maxheight, nrows)
    txt <- tktext(base,
                  bg = body.bgcolor,
                  fg = body.textcolor,
                  font = font,
                  height = textheight,
                  width = winwidth,
                  setgrid = 1,
                  takefocus = TRUE)
     lnames <- tktext(base,
                     bg = rowname.bgcolor,
                     fg = rowname.textcolor,
                     font = font,
                     height = textheight,
                     width = namewidth,
                     takefocus = TRUE)
    rnames <- tktext(base,
                     bg = rowname.bgcolor,
                     fg = rowname.textcolor,
                     font = font,
                     height = textheight,
                     width = namewidth,
                     takefocus = TRUE)
    xscroll <- tkscrollbar(base,
                           orient = "horizontal",
                           repeatinterval = 1,
                           command = function(...) {
                               tkxview(txt, ...)
                               tkxview(hdr, ...)
                               tkxview(ftr, ...)
                           })
    string.to.vector <- function(string.of.indices) {
        string.of.indices <- tclvalue(string.of.indices)
        as.numeric(strsplit(string.of.indices, split = " ")[[1]])
    }
    tkconfigure(txt, xscrollcommand = function(...) {
        tkset(xscroll, ...)
        xy <- string.to.vector(tkget(xscroll))
        tkxview.moveto(hdr, xy[1])
        tkxview.moveto(ftr, xy[1])
    })
    tkconfigure(hdr, xscrollcommand = function(...) {
        tkset(xscroll, ...)
        xy <- string.to.vector(tkget(xscroll))
        tkxview.moveto(txt, xy[1])
        tkxview.moveto(ftr, xy[1])
    })
    tkconfigure(ftr, xscrollcommand = function(...) {
        tkset(xscroll, ...)
        xy <- string.to.vector(tkget(xscroll))
        tkxview.moveto(hdr, xy[1])
        tkxview.moveto(txt, xy[1])
    })
    yscroll <- tkscrollbar(base,
                           orient = "vertical",
                           repeatinterval = 1,
                           command = function(...) {
                               tkyview(txt, ...)
                               tkyview(lnames, ...)
                               tkyview(rnames, ...)
                           })
    tkconfigure(txt, yscrollcommand = function(...) {
        tkset(yscroll, ...)
        xy <- string.to.vector(tkget(yscroll))
        tkyview.moveto(lnames, xy[1])
        tkyview.moveto(rnames, xy[1])
    })
    tkconfigure(lnames, yscrollcommand = function(...) {
        tkset(yscroll, ...)
        xy <- string.to.vector(tkget(yscroll))
        tkyview.moveto(txt, xy[1])
        tkyview.moveto(rnames, xy[1])
    })
    tkconfigure(rnames, yscrollcommand = function(...) {
        tkset(yscroll, ...)
        xy <- string.to.vector(tkget(yscroll))
        tkyview.moveto(txt, xy[1])
        tkyview.moveto(lnames, xy[1])
    })
    tkbind(txt, "<B2-Motion>", function(x, y) {
        tkscan.dragto(txt, x, y)
    })
## The next block just enables copying from the text boxes
{
    copyText.hdr <- function(){
        tcl("event", "generate",
              .Tk.ID(hdr),
              "<<Copy>>")}
    tkbind(hdr, "<Button-1>", function() tkfocus(hdr))
    editPopupMenu.hdr <- tkmenu(hdr, tearoff = FALSE)
    tkadd(editPopupMenu.hdr, "command", label = "Copy <Ctrl-C>",
              command = copyText.hdr)
    RightClick.hdr <- function(x,y) # x and y are the mouse coordinates
    {
        rootx <- as.integer(tkwinfo("rootx", hdr))
        rooty <- as.integer(tkwinfo("rooty", hdr))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", editPopupMenu.hdr, xTxt, yTxt)
    }
    tkbind(hdr, "<Button-3>", RightClick.hdr)
    tkbind(hdr, "<Control-KeyPress-c>", copyText.hdr)
    ##
    copyText.ftr <- function(){
        tcl("event", "generate",
              .Tk.ID(ftr),
              "<<Copy>>")}
    tkbind(ftr, "<Button-1>", function() tkfocus(ftr))
    editPopupMenu.ftr <- tkmenu(ftr, tearoff = FALSE)
    tkadd(editPopupMenu.ftr, "command", label = "Copy <Ctrl-C>",
              command = copyText.ftr)
    RightClick.ftr <- function(x,y) # x and y are the mouse coordinates
    {
        rootx <- as.integer(tkwinfo("rootx", ftr))
        rooty <- as.integer(tkwinfo("rooty", ftr))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", editPopupMenu.ftr, xTxt, yTxt)
    }
    tkbind(ftr, "<Button-3>", RightClick.ftr)
    tkbind(ftr, "<Control-KeyPress-c>", copyText.ftr)
    ##
    copyText.txt <- function(){
        tcl("event", "generate",
              .Tk.ID(txt),
              "<<Copy>>")}
    tkbind(txt, "<Button-1>", function() tkfocus(txt))
    editPopupMenu.txt <- tkmenu(txt, tearoff = FALSE)
    tkadd(editPopupMenu.txt, "command", label = "Copy <Ctrl-C>",
              command = copyText.txt)
    RightClick.txt <- function(x,y) # x and y are the mouse coordinates
    {
        rootx <- as.integer(tkwinfo("rootx", txt))
        rooty <- as.integer(tkwinfo("rooty", txt))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", editPopupMenu.txt, xTxt, yTxt)
    }
    tkbind(txt, "<Button-3>", RightClick.txt)
    tkbind(txt, "<Control-KeyPress-c>", copyText.txt)
    ##
    copyText.lnames <- function(){
        tcl("event", "generate",
              .Tk.ID(lnames),
              "<<Copy>>")}
    tkbind(lnames, "<Button-1>", function() tkfocus(lnames))
    editPopupMenu.lnames <- tkmenu(lnames, tearoff = FALSE)
    tkadd(editPopupMenu.lnames, "command", label = "Copy <Ctrl-C>",
              command = copyText.lnames)
    RightClick.lnames <- function(x,y) # x and y are the mouse coordinates
    {
        rootx <- as.integer(tkwinfo("rootx", lnames))
        rooty <- as.integer(tkwinfo("rooty", lnames))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", editPopupMenu.lnames, xTxt, yTxt)
    }
    tkbind(lnames, "<Button-3>", RightClick.lnames)
    tkbind(lnames, "<Control-KeyPress-c>", copyText.lnames)
    ##
        copyText.rnames <- function(){
        tcl("event", "generate",
              .Tk.ID(rnames),
              "<<Copy>>")}
    tkbind(rnames, "<Button-1>", function() tkfocus(rnames))
    editPopupMenu.rnames <- tkmenu(rnames, tearoff = FALSE)
    tkadd(editPopupMenu.rnames, "command", label = "Copy <Ctrl-C>",
              command = copyText.rnames)
    RightClick.rnames <- function(x,y) # x and y are the mouse coordinates
    {
        rootx <- as.integer(tkwinfo("rootx", rnames))
        rooty <- as.integer(tkwinfo("rooty", rnames))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", editPopupMenu.rnames, xTxt, yTxt)
    }
    tkbind(rnames, "<Button-3>", RightClick.rnames)
    tkbind(rnames, "<Control-KeyPress-c>", copyText.rnames)
}

    tktag.configure(hdr, "notwrapped", wrap = "none")
    tktag.configure(ftr, "notwrapped", wrap = "none")
    tktag.configure(txt, "notwrapped", wrap = "none")
    tktag.configure(lnames, "notwrapped", wrap = "none")
    tktag.configure(rnames, "notwrapped", wrap = "none")
    tkinsert(txt, "end", paste(paste(yy[-1], collapse = "\n"),
                               sep = ""), "notwrapped")

    tkgrid(txt, row = 1, column = 1, sticky = "nsew")
    if ("top" %in% colname.bar) {
        tkinsert(hdr, "end", paste(yy[1], sep = ""), "notwrapped")
        tkgrid(hdr, row = 0, column = 1, sticky = "ew")
    }
    if ("bottom" %in% colname.bar) {
        tkinsert(ftr, "end", paste(yy[1], sep = ""), "notwrapped")
        tkgrid(ftr, row = 2, column = 1, sticky = "ew")
    }
    if ("left" %in% rowname.bar) {
        tkinsert(lnames, "end",
                 paste(rowname.text, collapse = "\n"),
                 "notwrapped")
        tkgrid(lnames, row = 1, column = 0, sticky = "ns")
    }
    if ("right" %in% rowname.bar) {
        tkinsert(rnames, "end",
                 paste(rowname.text, collapse = "\n"),
                 "notwrapped")
        tkgrid(rnames, row = 1, column = 2, sticky = "ns")
    }
    tkconfigure(hdr, state = "disabled")
    tkconfigure(ftr, state = "disabled")
    tkconfigure(txt, state = "disabled")
    tkconfigure(lnames, state = "disabled")
    tkconfigure(rnames, state = "disabled")
    if (maxheight < nrows) {
        tkgrid(yscroll, row = 1, column = 3, sticky = "ns")
    }
    if (maxwidth < datawidth) {
        tkgrid(xscroll, row = 3, column = 1, sticky = "ew")
    }

    sortColumn <- function(n, decreasing=FALSE) {
      dataframe <<- dataframe[ order(dataframe[[n]], decreasing=decreasing), ]
      rownames(dataframe) <- seq(length=nrow(dataframe))
      .tkigraph.showData(dataframe,
                         colname.bgcolor = colname.bgcolor,
                         rowname.bgcolor = rowname.bgcolor,
                         body.bgcolor = body.bgcolor,
                         colname.textcolor = colname.textcolor,
                         rowname.textcolor = rowname.textcolor,
                         body.textcolor = body.textcolor,
                         font = font,
                         maxheight = maxheight,
                         maxwidth = maxwidth,
                         title = title,
                         rowname.bar = rowname.bar,
                         colname.bar = colname.bar,
                         rownumbers = rownumbers,
                         placement = placement,
                         plot.text=plot.text,
                         plot.command=plot.command,
                         suppress.X11.warnings = suppress.X11.warnings,
                         right=right,
                         showmean=showmean,
                         sort.button=sort.button,
                         inthis=base)
    }

    pf <- tkframe(base)
    if (is.null(inthis)) { tkgrid(pf, column=5, row=0, rowspan=10, sticky="new") }

    if (!is.null(showmean) && is.null(inthis)) {
      for (i in seq(along=showmean)) {
        tkgrid(tklabel(base, text=showmean[1]), sticky="nsew",
               column=0, padx=1, pady=1, columnspan=4)
      }
    }

    sortBut <- tkbutton(base, text="Sort otherwise", command=function() {})

    sortPopup <- function() {
      sortMenu <- tkmenu(base, tearoff=FALSE)
      sapply(seq(along=colnames(dataframe)),
             function(n) {
               tkadd(sortMenu, "command", label=colnames(dataframe)[n],
                     command=function() sortColumn(colnames(dataframe)[n]))
               label <- paste(colnames(dataframe)[n], "decreasing", sep=", ")
               tkadd(sortMenu, "command", label=label,
                     command=function() sortColumn(colnames(dataframe)[n],
                       decreasing=TRUE))
             })
      rootx <- as.integer(tkwinfo("rootx", sortBut))
      rooty <- as.integer(tkwinfo("rooty", sortBut))
      tkpopup(sortMenu, rootx, rooty)
    }

    if (!is.null(plot.command)) {
      but <- tkbutton(base, text=plot.text, command=plot.command)
      tkgrid(but, "in"=pf, sticky="ew", column=10, row=1, padx=1, pady=1)
    }

    if (sort.button) { tkgrid(sortBut, "in"=pf, sticky="ew", column=10, row=2,
                              padx=1, pady=1) }
    tkconfigure(sortBut, command=sortPopup)

    savebut <- tkbutton(base, text="Export table to file", command=function() {
      filename <- tkgetSaveFile(initialfile="data.txt",
                                defaultextension="txt",
                                title="Export as table")
      filename <- paste(as.character(filename), collapse=" ")
      write.table(dataframe, file=filename, row.names=FALSE, col.names=FALSE)
    })
    tkgrid(savebut, "in"=pf, sticky="ew", column=10, row=3, padx=1, pady=1)

    but <- tkbutton(base, text="Close", command=function() tkdestroy(base))
    tkgrid(but, "in"=pf, sticky="ew", column=10, row=4, padx=1, pady=1) 
    tkgrid.columnconfigure(pf, 0, weight=1)

    tkgrid.rowconfigure(base, 1, weight = 1)
    tkgrid.columnconfigure(base, 1, weight = 1)
    tkwm.maxsize(base, 2 + datawidth, nrows)
    tkwm.minsize(base, 2 + nchar(names(dataframe)[1]), 1)
invisible(NULL)
}

.tkigraph.net.moody.white <-
matrix( c(0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          1,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          1,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          1,1,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          1,1,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
          1,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,1,1,1,1,1,0,1,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,
          0,0,0,0,0,0,1,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,1,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,
          0,0,0,0,0,0,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,
          0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,1,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,0,
          0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,1,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0), nrow=23, ncol=23)


