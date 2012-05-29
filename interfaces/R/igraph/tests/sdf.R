
library(igraph)

sdf <- igraph:::sdf(id=1:10, color="black")
as.data.frame(sdf)

## access

sdf[1,"id"]
sdf[1:4, "id"]
sdf[, "id"]

sdf[1, "color"]
sdf[1:4, "color"]
sdf[, "color"]

## set

sdf2 <- sdf
sdf2[5, "id"] <- 100
as.data.frame(sdf2)

sdf2 <- sdf
sdf2[, "id"] <- 0
as.data.frame(sdf2)

sdf2 <- sdf
sdf2[2:10, "id"] <- 1
as.data.frame(sdf2)

sdf2 <- sdf
sdf2[, "color"] <- "white"
as.data.frame(sdf2)

sdf2 <- sdf
sdf2[5:6, "color"] <- "white"
as.data.frame(sdf2)

