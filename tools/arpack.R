#! /usr/bin/env r
# You must set the ARPACK environment variable, and it must
# point to the root directory of your (Fortran) ARPACK source

files <- list.files("src/arpack", pattern=".*\\.[fh]$", ignore.case=TRUE,
                    full.names=TRUE)

arpack <- Sys.getenv("ARPACK")
if (arpack=="") { stop("ARPACK environment variable is not set") }

cat("Copying files from ARPACK.\n")

for (f in basename(files)) {
  cat(f, "...")
  ff <- paste(arpack, "BLAS/", f, sep="")
  ff1 <- paste(arpack, "LAPACK/", f, sep="")
  ff2 <- paste(arpack, "UTIL/", f, sep="")
  ff3 <- paste(arpack, "SRC/", f, sep="")
  if (!file.exists(ff)) { ff <- ff1 }
  if (!file.exists(ff)) { ff <- ff2 }
  if (!file.exists(ff)) { ff <- ff3 }
  if (!file.exists(ff)) {
    stop( paste("File does not exists:", ff) )
  }
  tt <- paste("src/arpack", sep="/", f)
  file.copy(ff, tt, overwrite=TRUE)
  cat("OK\n")
}

## Correct etime_

cat("Correcting 'etime_' in second.f ...")
ll <- readLines("src/arpack/second.f")
etime <- grep("EXTERNAL[ \t]*ETIME", ll)
ll <- ll[-etime]
writeLines(ll, "src/arpack/second.f")
cat("OK\n")

cat("Rewriting with 'igrapharpack' prefix\n")

files <- grep("\\.f$", files, value=TRUE, ignore.case=TRUE)
oldnames <- toupper(sub("\\.f$", "", basename(files),
                        fixed=FALSE, ignore.case=TRUE))
newnames <- paste("IGRAPHARPACK", sep="", oldnames)

print(oldnames)
print(newnames)

for (f in files) {
  cat("Converting", f, "...")
  ll <- readLines(f)
  ig <- grep("igraph", ll)
  if (length(ig) != 0) {
    cat("skipping...\n")
    next
  }
  for (i in seq_along(oldnames)) {
    ll <- gsub(paste(sep="", oldnames[i], "\\b"),
               paste(sep="", newnames[i]), ll, perl=TRUE,
               ignore.case=TRUE)
  }
  writeLines(ll, con=f)
  cat("OK.\n")
}
