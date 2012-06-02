
#   IGraph R package
#   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
#   334 Harvard street, Cambridge, MA 02139 USA
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

# This file contains code, that was written by Magnus Torfason

# This function checks some very basic pitfalls that can cause
# incorrect merge results, and then merges using the merge function.
#
# Additionally, if indicate.origin=TRUE is specified, it will add 
# variables called was.in.x and was.in.y indicating if an entry was
# in the x/y part of the merge.
safer.merge = function(
        x, y, by = intersect(names(x), names(y)),
        by.x = by, by.y = by, all = FALSE, all.x = all, all.y = all,
        suffixes = c(".x",".y"), incomparables = NULL, 
        indicate.origin=FALSE, allow.duplicates=FALSE, 
        allow.duplicates.x=allow.duplicates, allow.duplicates.y=allow.duplicates, 
        ...)
{
    if ( !allow.duplicates.x & (sum(duplicated(x[,by.x]))!=0) )
    {
        stop("safer.merge: Duplicates found in x, but allow.duplicates.x is FALSE" )
    }
    if ( !allow.duplicates.y & (sum(duplicated(y[,by.y]))!=0) )
    {
        stop("safer.merge: Duplicates found in y, but allow.duplicates.y is FALSE" )
    }
    stopifnot( sum(is.na(x[,by.x])) == 0 )
    stopifnot( sum(is.na(y[,by.y])) == 0 )
    
    if ( indicate.origin )
    {
        stopifnot( !("was.in" %in% names(x) ) )
        stopifnot( !("was.in" %in% names(y) ) )
        x$was.in=TRUE
        y$was.in=TRUE
    }
    
    result = merge(x, y, by=by, 
        by.x=by.x, by.y=by.y, all=all, all.x=all.x, all.y=all.y,
        sort=TRUE, suffixes=suffixes, incomparables=incomparables, ...)
    
    if ( indicate.origin )
    {
        was.in.x = paste("was.in", suffixes[1], sep="")
        was.in.y = paste("was.in", suffixes[2], sep="")
        result[ is.na(result[,was.in.x]) , was.in.x] = FALSE
        result[ is.na(result[,was.in.y]) , was.in.y] = FALSE
    }
    
    return(result)
}
    
# Return graph vertices as a data.frame.
# If keep.attributes == TRUE, the data.frame contains all vertex attributes.
get.vertices.as.data.frame = function(graph, keep.attributes=FALSE)
{
    stopifnot( length(V(graph)$name) == length(unique(V(graph)$name))) # Vertex names must be unique
    d = data.frame(V=V(graph)$name)
    if ( keep.attributes )
    {
        # Bind attributes if requested
        attr.names = list.vertex.attributes(graph) 
        d.a = as.data.frame(lapply( attr.names, 
                      function(attr, g){get.vertex.attribute(g,attr)},
                      graph))
        names(d.a) = attr.names
        d = cbind(d,d.a)
    }
    return(d)
}

# Return graph edges as a mergable data.frame
# If the graph is undirected, the result will have all V1 <= V2
# If keep.attributes == TRUE, the data.frame contains all edge attributes.
get.edges.as.data.frame = function(graph, keep.attributes=FALSE)
{
    stopifnot( length(V(graph)$name) == length(unique(V(graph)$name))) # Vertex names must be unique
    m = get.edgelist(graph)
    if ( !is.directed(graph) )
    {
        ix.rev        = m[,1] > m[,2]       # Choose which pairs to switch
        m[ix.rev,1:2] = m[ix.rev,2:1]       # Switch
    }
    d = as.data.frame(m)
    if ( keep.attributes )
    {
        # Bind attributes if requested
        attr.names = list.edge.attributes(graph) 
        d.a = as.data.frame(lapply( attr.names, 
                      function(attr, g){get.edge.attribute(g,attr)},
                      graph))
        names(d.a) = attr.names
        d = cbind(d,d.a)
    }
    return(d)
}

# Create an intersection of two graphs.
# This function correctly intersects the graphs based
# on the name attributes, such that:
#   Any vertex that is in g1 and g2 is in the result.
#   Any edge that is in g1 and g2 is in the result.
# This function preserves some attributes, but not all.
graph.intersection.by.name = function(g1, g2, 
        keep.x.vertices          = FALSE,
        keep.x.vertex.attributes = FALSE,
        keep.x.edge.attributes   = FALSE)
{
    # Only undirected graphs are supported
    stopifnot(!is.directed(g1) & !is.directed(g2))

    # Construct data.frames of nodes and edges
    dv1 = get.vertices.as.data.frame(g1, keep.attributes=keep.x.vertex.attributes)
    dv2 = get.vertices.as.data.frame(g2)
    de1 = get.edges.as.data.frame(g1, keep.attributes=keep.x.edge.attributes)
    de2 = get.edges.as.data.frame(g2)

    # Merge data.frames and construct result
    if ( keep.x.vertices )
        dv = dv1
    else
        dv = safer.merge(dv1, dv2)
    de = safer.merge(de1, de2)
    g  = graph.data.frame(de, directed=FALSE, vertices=dv)
    return(g)
}

# Create union of two graphs.
# This function correctly unions the graphs based
# on the name attributes, such that:
#   Any vertex that is in g1 or g2 is in the result.
#   Any edge that is in g1 or g2 is in the result.
# However, this function does NOT preserve other attributes.
graph.union.by.name = function(g1, g2)
{
    # Only undirected graphs are supported
    stopifnot(!is.directed(g1) & !is.directed(g2))

    # Construct data.frames of nodes and edges
    dv1 = get.vertices.as.data.frame(g1)
    dv2 = get.vertices.as.data.frame(g2)
    de1 = get.edges.as.data.frame(g1)
    de2 = get.edges.as.data.frame(g2)

    # Merge data.frames and construct result
    dv = safer.merge(dv1, dv2, all=TRUE)
    de = safer.merge(de1, de2, all=TRUE)
    g  = graph.data.frame(de, directed=FALSE, vertices=dv)
    return(g)
}

# Create difference of two graphs.
# This function correctly diffs the graphs based
# on the name attributes. Note that this difference applies
# only to the edges, not the vertices, so that:
#   Any vertex that is in g1 is in the result.
#   Any edge that is in g1 but not in g2 is in the result.
# However, this function does NOT preserve other attributes.
graph.difference.by.name = function(g1, g2,
            keep.x.vertex.attributes = FALSE,
            keep.x.edge.attributes   = FALSE)
{
    # Only undirected graphs are supported
    stopifnot(!is.directed(g1) & !is.directed(g2))

    # Construct data.frames of nodes and edges
    dv = get.vertices.as.data.frame(g1, keep.attributes=keep.x.vertex.attributes)
    de1 = get.edges.as.data.frame(g1, keep.attributes=keep.x.edge.attributes)
    de2 = get.edges.as.data.frame(g2)

    # Merge data.frames and construct result
    de = safer.merge(de1, de2, all=TRUE, indicate.origin=TRUE) # All edges
    de = de[ (de$was.in.x & !de$was.in.y) ,                    # Remove g2 edges
             setdiff(names(de),c("was.in.x","was.in.y")) ]     # Keep attributes from g1     
    g  = graph.data.frame(de, directed=FALSE, vertices=dv)     # The new graph
    return(g)
}


# Test/demonstration code
if (FALSE)
{

# Load library and set parameters for better printing
library(igraph)
igraph.options(print.vertex.attributes=TRUE)
igraph.options(print.edge.attributes=TRUE)

# Define input graphs
g1 <- graph.formula(a-b-c)
V(g1)$v.attr=c(1,2,3)
E(g1)$e.attr=c(5,7)
g2 <- graph.formula(b-c-d)

# Test the functions
graph.intersection.by.name(g1,g2) # Vertices are intersected as well
graph.union.by.name(g1,g2)        # Vertices are unioned as well
graph.difference.by.name(g1,g2)   # Vertices from x (g1) are used

# graph.intersection.by.name() has some extra parameters
graph.intersection.by.name(g1,g2,
    keep.x.vertices          = TRUE) # Keep all x vertices (only intersect edges)

graph.intersection.by.name(g1,g2,
    keep.x.vertices          = FALSE, # Keep all x vertices (only intersect edges)
    keep.x.vertex.attributes = TRUE, # Don't throw away V(g1) attributes
    keep.x.edge.attributes   = TRUE) # Don't throw away E(g1) attributes

# graph.difference.by.name() has some extra parameters
graph.difference.by.name(g1,g2,
    keep.x.vertex.attributes = TRUE, # Don't throw away V(g1) attributes
    keep.x.edge.attributes   = TRUE) # Don't throw away E(g1) attributes

}
