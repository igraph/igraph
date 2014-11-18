
generate_uuid <- function(use_time = NA) {
  .Call("UUID_gen", as.logical(use_time), PACKAGE="igraph")
}
