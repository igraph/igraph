require("reticulate")
require("jsonlite")
require("readtext")

source_python("LaplaceDynamic.py",envir=globalenv())

resultsfile = 'Results_Laplace.JSON'

#read data function
#uses Python Networkx Package
main(inputfile = 'Laplace/test/dummy_unweighted_t/as19971108.txt', sequencedirectory = 'Laplace/test/dummy_unweighted_t/sequences', outputfile = resultsfile, weighted_graph = FALSE , norm = TRUE)

results_dataframe <- as.data.frame(unlist(fromJSON(resultsfile)), optional = TRUE, col.names = c("Laplacian Centrality"))

results_dataframe <- results_dataframe[order(results_dataframe[,1],decreasing = TRUE),]

results_dataframe

