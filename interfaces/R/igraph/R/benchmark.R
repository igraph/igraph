
benchmark_context <- function(id, version) { }

benchmark_this <- function(label, ..., init={}, init_each={}, replications=100,
                           environment=new.env(), child=FALSE) {

  init <- substitute(init)
  init_each <- substitute(init_each)
  
  arguments <- match.call()[-1]
  parameters <- names(arguments)
  if (is.null(parameters)) {
    parameters = as.character(arguments)
  } else {
    keep <- ! parameters %in% names(formals(benchmark_this))
    arguments <- arguments[keep]
    parameters <- parameters[keep]
  }

  if (length(arguments) == 0) {
    warning("benchmark_this called with no expressions to be evaluated.")
    return(invisible(data.frame()))
  }

  
  n <- list(tests=length(arguments), replications=length(replications))
  replications <- rep(replications, each=n$tests)
  labels <- rep(ifelse(parameters=='', as.character(arguments), parameters),
    n$replications)

  eval(init, envir=environment)

  replicator <- function(test, repl) {
    res <- replicate(repl, {
      eval(init_each, environment)
      st <- try(system.time(eval(test, environment)), silent=TRUE)
      if (inherits(st, "try-error")) {
        st <- system.time(eval({}, environment))
        st[] <- NA
      }
      if (!child) { st <- st[!grepl("\\.child$", names(st))] }
      st
    })
    rm <- rowMeans(res)
    rs <- apply(res, 1, sd)
    structure(c(rm, rs), names=c(paste("mean", sep=".", names(rm)),
                           paste("sd", sep=".", names(rs))))
  }
  
  timings <- mapply(replicator, arguments, replications)
  
  result <- data.frame(row.names=NULL, test=labels,
                       replications=as.integer(replications), t(timings))

  result[order(result[["mean.user.self"]]), , drop=FALSE]
}

