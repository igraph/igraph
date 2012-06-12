/* 
   IGraph library Java interface.
   Copyright (C) 2006-2012  Tamas Nepusz <ntamas@gmail.com>
   Pázmány Péter sétány 1/a, 1117 Budapest, Hungary
   
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

import org.junit.*;
import static org.junit.Assert.*;

public class BasicGraphTests {
    @Test
    public void testEmptyGraph() {
        Graph graph = Graph.Empty(10, false);
        assertEquals("empty graph should have ten vertices", 10, graph.vcount());
        assertEquals("empty graph should have no edges", 0, graph.ecount());
        assertFalse(graph.isDirected());

        graph = Graph.Empty(5, true);
        assertEquals("empty graph should have five vertices", 5, graph.vcount());
        assertEquals("empty graph should have no edges", 0, graph.ecount());
        assertTrue(graph.isDirected());
    }

    @Test(expected=CoreException.class)
    public void testEmptyGraphException() {
        Graph graph = Graph.Empty(-1, true);
    }

    @Test
    public void testFullGraph() {
        Graph graph = Graph.Full(10, false, false);
        assertEquals("full graph should have 10 vertices", 10, graph.vcount());
        assertEquals("full graph should have 45 edges", 45, graph.ecount());

        graph = Graph.Full(10, true, false);
        assertEquals("full directed graph should have 10 vertices", 10, graph.vcount());
        assertEquals("full directed graph should have 90 edges", 90, graph.ecount());

        graph = Graph.Full(10, false, true);
        assertEquals("full graph with loops should have 10 vertices", 10, graph.vcount());
        assertEquals("full graph with loops should have 45 edges", 55, graph.ecount());
    }

    @Test(expected=CoreException.class)
    public void testFullGraphException() {
        Graph graph = Graph.Full(-1, true, true);
    }

    @Test
    public void testNeighbors() {
        Graph graph = Graph.Full(10, false, false);
        double[] neighbors = graph.neighbors(3, NeighborMode.OUT);
        double[] expectedNeighbors = { 0, 1, 2, 4, 5, 6, 7, 8, 9 };

        for (int i = 0; i < expectedNeighbors.length; i++)
            assertEquals("neighbor list element "+i+" invalid for undirected full graph",
                    expectedNeighbors[i], neighbors[i], 0.0);
    }

    @Test(expected=CoreException.class)
    public void testNeighborsInvalidVertexID() {
        Graph graph = Graph.Full(10, false, false);
        double[] neighbors = graph.neighbors(-1, NeighborMode.OUT);
    }

    @Test(expected=NullPointerException.class)
    public void testNeighborsInvalidNeighborMode() {
        Graph graph = Graph.Empty(3, false);
        double[] neighbors = graph.neighbors(2, null);
    }

    @Test
    public void testGetEids() {
        Graph graph = Graph.Full(10, false, false);
        double[] pairs = {2, 7, 5, 1, 3, 4};
        double[] eids = graph.getEids(pairs, false);
        double[] expectedEids = { 21, 12, 24 };

        for (int i = 0; i < expectedEids.length; i++)
            assertEquals("edge ID list element "+i+" invalid for undirected full graph",
                    expectedEids[i], eids[i], 0.0);
    }

    @Test
    public void testDegreeAll1() {
        Graph graph = Graph.Full(10, false, true);
        double[] expectedDegrees = { 9, 9, 9, 9, 9, 9, 9, 9, 9, 9 };
        double[] degrees = graph.degree(null, NeighborMode.OUT, false);

        for (int i = 0; i < expectedDegrees.length; i++) {
            assertEquals("degree list element "+i+" invalid for undirected full graph (no loops)",
                    expectedDegrees[i], degrees[i], 0.0);
		}
	}

    @Test
    public void testDegreeAll2() {
        Graph graph = Graph.Full(10, false, true);
        double[] expectedDegrees = { 9, 9, 9, 9, 9, 9, 9, 9, 9, 9 };
        double[] degrees = graph.degree(VertexSet.ALL, NeighborMode.OUT, false);

        for (int i = 0; i < expectedDegrees.length; i++) {
            assertEquals("degree list element "+i+" invalid for undirected full graph (no loops)",
                    expectedDegrees[i], degrees[i], 0.0);
		}
	}

	@Test
	public void testDegreeSingle() {
        Graph graph = Graph.Full(10, false, true);
        double[] expectedDegrees = { 9 };
        double[] degrees = graph.degree(new VertexSet(2), NeighborMode.OUT, false);
        for (int i = 0; i < expectedDegrees.length; i++)
            assertEquals("degree list element "+i+" invalid for undirected full graph (no loops)",
                    expectedDegrees[i], degrees[i], 0.0);
    }

	@Test
	public void testDegreeMultiple() {
        Graph graph = Graph.Full(10, false, true);
		long[] ids = { 2, 3, 7 };
        double[] expectedDegrees = { 9, 9, 9 };
        double[] degrees = graph.degree(new VertexSet(ids), NeighborMode.OUT, false);
        for (int i = 0; i < expectedDegrees.length; i++)
            assertEquals("degree list element "+i+" invalid for undirected full graph (no loops)",
                    expectedDegrees[i], degrees[i], 0.0);
    }
};
