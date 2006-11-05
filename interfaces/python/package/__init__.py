#   IGraph library.
#   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
#   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
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
from igraph._igraph import *
from igraph._igraph import __version__, __build_date__

import os
import gzip
from tempfile import mkstemp

class Graph(_igraph.Graph):
    def degree_distribution(self):
	pass

    def write_graphmlz(self, f, compresslevel=9):
	"""write_graphmlz(f, compresslevel=9)
    
	Writes the graph to a zipped GraphML file. The library uses
	the gzip compression algorithm, so the resulting file can be
	unzipped with regular gzip uncompression (like C{gunzip} or C{zcat}
	from Unix command line) or the Python C{gzip} module.

	Uses a temporary file to store intermediate GraphML data, so
	make sure you have enough free space to store the unzipped
	GraphML file as well.

	@param f: the name of the file to be written.
	@param compresslevel: the level of compression. 1 is fastest and
 	  produces the least compression, and 9 is slowest and produces
	  the most compression."""
	tmpfilename=None
	try:
	    tmpfileno, tmpfilename = mkstemp(text=True)
	    os.close(tmpfileno)
	    self.write_graphml(tmpfilename)
	    outf = gzip.GzipFile(f, "wb", compresslevel)
	    inf = open(tmpfilename)
	    for line in inf: outf.write(line)
	    inf.close()
	    outf.close()
	finally:
	    if tmpfilename is not None: os.unlink(tmpfilename)

    def Read_GraphMLz(cls, f, *params, **kwds):
	"""Read_GraphMLz(f, directed=True, index=0) -> Graph

	Reads a graph from a zipped GraphML file.

	@param f: the name of the file
	@param index: if the GraphML file contains multiple graphs,
	  specified the one that should be loaded. Graph indices
	  start from zero, so if you want to load the first graph,
	  specify 0 here.
        @return: the loaded graph object"""
	tmpfilename=None
	try:
	    tmpfileno, tmpfilename = mkstemp(text=True)
	    os.close(tmpfileno)
	    inf = gzip.GzipFile(f, "rb")
	    outf = open(tmpfilename, "wt")
	    for line in inf: outf.write(line)
	    inf.close()
	    outf.close()
	    result=cls.Read_GraphML(tmpfilename)
	finally:
	    if tmpfilename is not None: os.unlink(tmpfilename)
	return result
    Read_GraphMLz = classmethod(Read_GraphMLz)
