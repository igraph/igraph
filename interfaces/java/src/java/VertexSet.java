/* 
   IGraph library Java interface.
   Copyright (C) 2007  Tamas Nepusz <ntamas@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA 

*/

/*

ATTENTION: This is a highly experimental, proof-of-concept Java interface.
Its main purpose was to convince me that it can be done in finite time :)
The interface is highly incomplete, at the time of writing even some
essential functions (e.g. addEdges) are missing. Since I don't use Java
intensively, chances are that this interface gets finished only if there
is substantial demand for it and/or someone takes the time to send patches
or finish it completely.

*/

package net.sf.igraph;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;
import net.sf.igraph.util.LongRangeIterator;

/// Class representing a vertex (sub)set of a graph
public class VertexSet implements Iterable<Long> {
    /// List of vertex IDs in this vertex (sub)set
    protected List<Long> ids = null;

    /// Graph object we are attached to (if any)
    protected Graph graph = null;

    /// Static instance variable that refers to a vertex set with all the vertices in a graph
    static final VertexSet ALL = new VertexSet();

    /// Constructor that creates a VertexSet referring to all the vertices with no attached graph
    public VertexSet() {}

    /// Constructor that creates a VertexSet referring to all the vertices attached to a given graph
    public VertexSet(Graph graph) {
        this.graph = graph;
    }

    /// Constructor that creates a VertexSet referring to a single vertex attached to a given graph
    public VertexSet(Long vertexID, Graph graph) {
        this(graph);
        this.ids = new Vector<Long>(1);
        this.ids.add(vertexID);
    }

    /// Constructor that creates a VertexSet referring to a single vertex with no attached graph
    public VertexSet(Long vertexID) {
        this(vertexID, null);
    }

    /// Constructor that creates a VertexSet referring to a single vertex attached to a given graph
    public VertexSet(Integer vertexID, Graph graph) {
		this(Long.valueOf(vertexID), graph);
    }

    /// Constructor that creates a VertexSet referring to a single vertex with no attached graph
    public VertexSet(Integer vertexID) {
        this(vertexID, null);
    }

    /// Constructor that creates a VertexSet referring to a single vertex attached to a given graph
    public VertexSet(long vertexID, Graph graph) {
		this(Long.valueOf(vertexID), graph);
    }

    /// Constructor that creates a VertexSet referring to a single vertex with no attached graph
    public VertexSet(long vertexID) {
        this(vertexID, null);
    }

    /// Constructor that creates a VertexSet referring to a predefined set of vertices attached to a graph
    public VertexSet(Collection<? extends Long> vertexIDs, Graph graph) {
        this(graph);
        this.ids = new Vector<Long>(vertexIDs);
    }

    /// Constructor that creates a VertexSet referring to a predefined set of vertices with no attached graph
    public VertexSet(Collection<? extends Long> vertexIDs) {
        this(vertexIDs, null);
    }

    /// Constructor that creates a VertexSet referring to a predefined set of vertices attached to a graph
    public VertexSet(long[] vertexIDs, Graph graph) {
        this(graph);
        this.ids = new Vector<Long>(vertexIDs.length);
        for (int i = 0; i < vertexIDs.length; i++)
            this.ids.add(Long.valueOf(vertexIDs[i]));
    }

    /// Constructor that creates a VertexSet referring to a predefined set of vertices with no attached graph
    public VertexSet(long[] vertexIDs) {
        this(vertexIDs, null);
    }

    /// Returns an iterator for the vertex IDs (may change in the future!)
    public Iterator<Long> iterator() {
        if (ids == null) {
            /* No IDs specified. If we are assigned to a graph, get the vertex
             * count and iterate over all the vertices */
            if (graph == null)
                throw new UnsupportedOperationException("VertexSet is not assigned to a graph");

            return new LongRangeIterator(0L, graph.vcount());
        }

        return this.ids.iterator();
    }

	/// Returns an array for the vertex IDs
	public long[] getIdArray() {
		if (this.ids == null) {
            /* No IDs specified. If we are assigned to a graph, get the vertex
             * count and return all the IDs */
            if (graph == null)
                throw new UnsupportedOperationException("VertexSet is not assigned to a graph");

			long[] result = new long[(int)graph.vcount()];
			for (int i = 0; i < result.length; i++)
				result[i] = i;

            return result;
		}

		long[] result = new long[this.ids.size()];
		for (int i = 0; i < result.length; i++)
			result[i] = this.ids.get(i);
		return result;
	}

    /**
     * Returns the type hint of this vertex set.
     *
     * This variable tells which igraph_vs_t constructor should be used in the C core
     * of igraph. It is a non-negative integer with the following meanings:
     *
     * - 0 = igraph_vs_all
     * - 1 = igraph_vs_1
     * - 2 = igraph_vs_vector
     *
     * Other igraph_vs_t constructors are not supported yet.
     */
    public int getTypeHint() {
        if (ids == null)
            return 0;

        if (ids.size() == 1)
            return 1;

        return 2;
    }

    /// Get the graph we are attached to
    public Graph getGraph() { return this.graph; }

    /// Set the graph we are attached to
    public void setGraph(Graph graph) { this.graph = graph; }

    /// Set the graph we are attached to
    public void attach(Graph graph) { setGraph(graph); }

    /// Detach from the graph we are attached to
    public void detach() { this.graph = null; }
}
