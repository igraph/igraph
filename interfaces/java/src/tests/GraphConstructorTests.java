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

import org.junit.*;
import static org.junit.Assert.*;

public class GraphConstructorTests {
    @Test
    public void testStarGraph() {
		Graph graph = Graph.Star(10, StarMode.OUT, 2);
        assertEquals("star graph must have ten vertices", 10, graph.vcount());
        assertEquals("star graph must have nine edges", 9, graph.ecount());
        assertTrue("star graph must be directed", graph.isDirected());

        double[] expectedEdgelist = { 2, 0, 2, 1, 2, 3, 2, 4, 2, 5, 2, 6, 2, 7, 2, 8, 2, 9 };
        double[] edgelist = graph.getEdgelist(false);

        for (int i = 0; i < expectedEdgelist.length; i++) {
            assertEquals("edge list element "+i+" mismatch", expectedEdgelist[i], edgelist[i], 0);
        }
	}

	@Test(expected=CoreException.class)
	public void testInvalidStarGraphException1() {
		Graph graph = Graph.Star(-1, StarMode.OUT, 0);
	}

	@Test(expected=NullPointerException.class)
	public void testInvalidStarGraphException2() {
		Graph graph = Graph.Star(10, null, 0);
	}

	@Test(expected=CoreException.class)
	public void testInvalidStarGraphException3() {
		Graph graph = Graph.Star(10, StarMode.IN, 15);
	}
};
