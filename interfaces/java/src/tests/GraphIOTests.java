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

import java.io.File;
import java.io.IOException;

import org.junit.*;
import static org.junit.Assert.*;

public class GraphIOTests {
    @Test
    public void testReadWritePajek() throws IOException {
		Graph graph = Graph.Star(10, StarMode.OUT, 2);
        Graph graph2;

        File f = File.createTempFile("igraph", null);

        graph.writePajek(f);
        graph2 = Graph.ReadPajek(f);
        assertEquals("star graph from Pajek file must have ten vertices", 10, graph2.vcount());
        assertEquals("star graph from Pajek file must have nine edges", 9, graph2.ecount());
        assertTrue("star graph from Pajek file must be directed", graph2.isDirected());
        assertTrue("graph loaded from Pajek file should be isomorphic to original",
                graph2.isIsomorphic(graph));
	}
    
    @Test(expected=IOException.class)
    public void testWriteIOError() throws IOException {
		Graph graph = Graph.Star(10, StarMode.OUT, 2);

        File f = new File("////this/is/almost surely an invalid filename");
        graph.writePajek(f);
	}
    
    @Test(expected=IOException.class)
    public void testReadIOError() throws IOException {
        File f = new File("////this/is/almost surely a nonexistent filename");
		Graph graph = Graph.ReadPajek(f);
	}
};
