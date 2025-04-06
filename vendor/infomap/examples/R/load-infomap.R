loadInfomap <- function() {
	message("Loading Infomap...")
	dyn.load(paste("infomap", .Platform$dynlib.ext, sep=""))
	source("infomap.R")
	# The cacheMetaData(1) will cause R to refresh its object tables. Without it, inheritance of wrapped objects may fail.
	cacheMetaData(1)
}


# Only load if not already loaded
if (!exists("InfomapWrapper")) {
	loadInfomap()
} else {
	message("Infomap already loaded!")
}