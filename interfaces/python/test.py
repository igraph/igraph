import sys;
import gc;
import tempfile;
import os;
import types;
import weakref;
from time import clock;

sys.path.insert(0, 'interfaces/python/.libs');
import igraph;

longtests=False;

# Recursively expand slist's objects
# into olist, using seen to track
# already processed objects.
def _getr(slist, olist, seen):
    for e in slist:
	if id(e) in seen:
	    continue
	seen[id(e)] = None
	olist.append(e)
	tl = gc.get_referents(e)
	if tl:
	    _getr(tl, olist, seen)
					      
# The public function.
def get_all_objects(gcl = None):
    """Return a list of all live Python
    objects, not including the list itself."""
    if gcl==None:
	gcl = gc.get_objects()
    olist = []
    seen = {}
    # Just in case
    seen[id(gcl)] = None
    seen[id(olist)] = None
    seen[id(seen)] = None
    # _getr does the real work.
    _getr(gcl, olist, seen)
    return olist
									  
######################## Testing framework routines ########################
def section(msg):
    global lastsection;
    lastsection=msg;
    sys.stdout.write("\n"+msg+"\n");
    sys.stdout.write("="*len(msg)+"\n\n");
    
def start(msg):
    global lasttest;
    lasttest=msg;
    sys.stdout.write(msg+"... ");
    sys.stdout.flush();
    
def ok():
    global passed;
    passed=passed+1;
    print "ok.";

def fail():
    global failed, lastsection, lasttest;
    failed.append((lastsection, lasttest));
    print "FAILED.";

def skip():
    global skipped;
    skipped=skipped+1;
    print "result validity not checked.";
    
def test(f):
    if f:
	ok()
    else:
	fail()

def results():
    global passed, failed
    print "========================================================"
    print "Passed: %d, failed: %d, skipped: %d, percentage: %f" % (passed, len(failed), skipped, float(passed)/(passed+len(failed))*100)
    if len(failed)>0:
	print "SOME TESTS FAILED!"
	print "Check the implementation of the Python module."
	print "Failed test cases are:"
	for (s, f) in failed:
	    print " *", f
	    print "   in section:", s
    else:
	print "Everything went OK."

def test_leaks(objects):
    for o in objects:
	if (type(o) == types.IntType): continue
	rc=sys.getrefcount(o)
	refs=gc.get_referrers(o)
	if rc>len(refs)+1:
	    print o, hex(id(o)), type(o), sys.getrefcount(o), len(refs)
	    print "Known referrers:"
	    print "================"
	    for i in refs:
		print hex(id(i)), type(i)
		if type(i) == types.ListType: print i
	    print "======================================================="
	    print "The current stack frame was:", sys._getframe()
	del refs

def unique(l):
    l2=l
    l2.sort()
    l3=[]
    prev=None
    for item in l:
	if prev==item:
	    pass
	else:
	    l3.append(item)
	prev=item
    return l3
	
passed=0; failed=[]; skipped=0;

section("General tests")

start("Testing igraph.InternalError exception")
try:
    raise igraph.InternalError
except igraph.InternalError, e:
    ok();
else:
    fail();
    
start("Creating default igraph.Graph object")
g=igraph.Graph(); ok();

start("Destroying igraph.Graph object")
g=0; ok();

start("Creating weak reference to an igraph.Graph object")
g=igraph.Graph(); h=weakref.ref(g);
test(h() == g);

start("Destroying graph and checking the weak reference")
del g
test(h() == None);
del h

start("Creating igraph.Graph object with 5 vertices")
g=igraph.Graph(5); test(g.vcount() == 5);

start("Creating directed igraph.Graph object from edge list")
g=igraph.Graph(edges=[(0,1), (1,2), (2,3), (3,0)], directed=True);
test(g.vcount() == 4 and g.ecount() == 4);

start("Creating igraph.Graph object with -2 vertices (should throw exception)")
try:
    g=igraph.Graph(-2);
except AssertionError, e:
    ok();
    print "Expected exception arrived: ", e
else:
    fail();
    
start("Representing a directed igraph.Graph object with 5 vertices as a string")
g=igraph.Graph(n=5, directed=True); print g; ok();

start("Checking whether it was really directed")
test(g.is_directed());

start("Extending this graph with two and three more vertices (in two steps)")
g.add_vertices(2);
g.add_vertices(3);
test(g.vcount() == 10);

start("Extending this graph with one and four more vertices (in one step)")
g.add_vertices(1).add_vertices(4);
test(g.vcount() == 15);

start("Deleting the 3rd vertex from the graph");
g.delete_vertices(2);
test(g.vcount() == 14);

start("Deleting the 1st and 7th vertices from the graph in two steps");
g.delete_vertices([0]).delete_vertices([6]);
test(g.vcount() == 12);

start("Trying to delete with erroneous vertex list");
try:
    g.delete_vertices([-3]);
except TypeError, e:
    ok();
    print "Expected exception arrived: ", e
else:
    fail();

start("Creating a new graph, extending it to 45 vertices and adding an edge");
g=igraph.Graph();
g.add_vertices(44);
g.add_edges([(0, 1)]);
test(g.vcount() == 45 and g.ecount() == 1);

start("Recreating a graph with 5 nodes and a couple of edges");
g=igraph.Graph(n=5);
g.add_edges([(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (4, 1), (0, 3), (2, 2)]);
test(g.vcount() == 5 and g.ecount() == 8);

start("Testing node degrees without self-loops");
degs=g.degree([0, 1, 2, 3, 4]);
test(degs == [3, 3, 2, 3, 3]);

start("Testing node degrees with self-loops");
degs=g.degree([0,1,2], igraph.ALL, True);
test(degs == [3, 3, 4]);

start("Testing special case: retrieving node degrees for no nodes");
degs=g.degree([]);
test(degs == []);

start("Testing the case of invalid node ID");
try:
    degs=g.degree(666);
except igraph.InternalError, e:
    ok();
    print "Expected exception arrived: ", e
else:
    fail();

start("Testing neighbor retrieval #1");
neis=map(lambda x: g.neighbors(x), range(5))
for n in neis: n.sort()
test(neis == [[1,3,4], [0,2,4], [1,2,2,3], [0,2,4], [0,1,3]])
start("Testing neighbor retrieval #2");
neis=g.successors(0)
neis.sort()
test(neis == [1,3,4])

start("Testing edge deletion")
g.delete_edges([(0, 4), (2, 2)]);
test(g.vcount() == 5 and g.ecount() == 6 and g.degree() == [2, 3, 2, 3, 2])

start("Calculating graph diameter")
test(g.diameter(True, True) == 2)

section("Graph generators")

start("Generating a graph from the Graph Atlas")
g=igraph.Graph.Atlas(4)
ok()

start("Generating a graph using the Barabasi-Albert model, n=1000, m=10")
g=igraph.Graph.Barabasi(1000, 10)
ok()

if longtests:
    start("Calculating graph diameter (should take several seconds)")
    t1=clock()
    sys.stdout.write(str(g.diameter())+", took ")
    t2=clock()
    sys.stdout.write(str(t2-t1)+"s, ")
    ok()

start("Generating a graph using the Erdos-Renyi model, n=1000, m=500")
g=igraph.Graph.Erdos_Renyi(n=1000, m=500)
test(g.vcount() == 1000 and g.ecount() == 500)

start("Generating a graph using the Erdos-Renyi model, n=1000, p=0.1")
g=igraph.Graph.Erdos_Renyi(n=1000, p=0.1)
skip() # test(g.vcount() == 1000)

start("Trying ill-formed arguments for the Erdos-Renyi model (missing n and p)")
try:
    g=igraph.Graph.Erdos_Renyi(n=1000)
except TypeError, e:
    ok();
    print "Expected exception arrived: ", e
else:
    fail();
    
start("Trying ill-formed arguments for the Erdos-Renyi model (p out of domain)")
try:
    g=igraph.Graph.Erdos_Renyi(n=1000, p=2.0)
except ValueError, e:
    ok();
    print "Expected exception arrived: ", e
else:
    fail();

start("Generating full directed graph with 10 nodes and self-loops")
g=igraph.Graph.Full(10, True, True)
test(g.vcount() == 10 and g.ecount() == 100)

start("Generating growing random graph with 10 nodes and 2 new edges per step")
g=igraph.Graph.Growing_Random(10, 2)
test(g.vcount() == 10 and g.ecount() == 18)

start("Generating directed outgoing star graph with 10 nodes and a center of node 5")
g=igraph.Graph.Star(10, igraph.STAR_OUT, 5)
neis=g.neighbors(5)
neis.sort()
test(g.vcount() == 10 and g.ecount() == 9 and neis == [0,1,2,3,4,6,7,8,9])

start("Generating directed closed ring graph with 10 nodes and mutual edges")
g=igraph.Graph.Ring(10, directed=True, mutual=True)
test(g.vcount() == 10 and g.ecount() == 20)

start("Generating an undirected graph with a given degree sequence")
g=igraph.Graph.Degree_Sequence([3, 2, 2, 1])
test(g.vcount() == 4 and g.degree(loops=True) == [3, 2, 2, 1])

start("Generating a ternary directed tree with 40 nodes")
g=igraph.Graph.Tree(40, 3, type=igraph.TREE_OUT)
test(g.vcount() == 40 and g.ecount() == 39)

section("Structural properties")

start("Checking whether this tree is weakly connected")
test(g.is_connected(igraph.WEAK))

start("Checking whether this tree is strongly connected")
test(g.is_connected(igraph.STRONG) == False)

start("Checking whether there's an edge from vertex 0 to 1")
test(g.are_connected(0,1))

start("Calculating average path length, considering paths as undirected")
sys.stdout.write(str(g.average_path_length(directed=False))+", ")
skip()

start("Calculating average path length, considering paths as directed (1)")
sys.stdout.write(str(g.average_path_length(directed=True))+", ")
skip()

start("Calculating average path length, considering paths as directed (1)")
sys.stdout.write(str(g.average_path_length(directed=True, unconn=False))+", ")
skip()

start("Calculating node betweennesses, considering paths as undirected")
sys.stdout.write(str(g.betweenness(directed=False))+", ")
skip()

start("Calculating node betweennesses")
sys.stdout.write(str(g.betweenness())+", ")
skip()

start("Generating a directed graph using the Barabasi-Albert model, n=50, m=10")
g=igraph.Graph.Barabasi(50, 10, directed=True)
ok()

start("Calculating closeness centralities, mode=IN")
sys.stdout.write(str(g.closeness(mode=igraph.IN))+", ")
skip()

start("Calculating closeness centralities, mode=OUT")
sys.stdout.write(str(g.closeness(mode=igraph.OUT))+", ")
skip()

start("Calculating closeness centralities, mode=ALL, several nodes")
sys.stdout.write(str(g.closeness([4,2,7]))+", ")
skip()

start("Calculating clusters in an undirected Erdos-Renyi graph, n=1000, m=1000")
comps=igraph.Graph.Erdos_Renyi(n=1000, m=2000).components()
compsizes=[0]*(max(comps)+1)
for c in comps:
    compsizes[c]=compsizes[c]+1
sys.stdout.write(str(compsizes)+", ")
skip()

start("Calculating diameter for an unconnected graph")
g2=igraph.Graph(n=3, edges=[(1,2)])
test(g2.diameter(unconn=False)==3)

start("Calculating edge betweennesses in a graph")
g=igraph.Graph(edges=[(0,1), (2,1), (1,3), (3,4), (3,5)], directed=True)
test(g.edge_betweenness(directed=False) == [5,5,9,5,5])

start("Calculating shortest paths to node 4 in the graph")
test(g.get_shortest_paths(4, mode=igraph.IN) == [[4,3,1,0], [4,3,1], [4,3,1,2], [4,3], [4], []])

start("Adding some more edges (to make sense for a spanning tree calculation)")
g.add_edges([(0,2), (4, 5), (0, 4), (2, 5)])
test(g.vcount() == 6 and g.ecount() == 9)

start("Calculating all of the shortest paths from node 0 in the graph")
test(g.get_all_shortest_paths(0, mode=igraph.OUT) == [[0], [0,1], [0,2], [0,1,3], [0,4], [0,2,5], [0,4,5]])

start("Calculating unweighted spanning tree")
g2=g.spanning_tree()
test(g2.vcount() == 6 and g2.ecount() == 5)

start("Calculating weighted spanning tree")
g2=g.spanning_tree([1,1,1,1,1,2,2,2,2])
test(g2.vcount() == 6 and g2.ecount() == 5)

start("Adding some duplicate edges and loops to the spanning tree and simplifying")
g2.add_edges([(2,2), (3,3), (0,2), (0,2), (1,3)])
g2.simplify()
test(g2.vcount() == 6 and g2.ecount() == 6)

start("Testing igraph_subcomponent on a graph, mode=IN")
g=igraph.Graph(edges=[(0,1), (1,2), (2,0), (3,4), (4,5), (5,3), (0,3)], directed=True)
comps=g.subcomponent(2, mode=igraph.IN)
test(unique(comps)==[0,1,2])

start("Testing igraph_subcomponent on a graph, mode=OUT")
comps=g.subcomponent(3, mode=igraph.OUT)
test(unique(comps)==[3,4,5])

start("Testing igraph_subcomponent on a graph, mode=ALL")
comps=g.subcomponent(3, mode=igraph.ALL)
test(unique(comps)==[0,1,2,3,4,5])

start("Testing igraph_subgraph");
g2=g.subgraph(g.subcomponent(2, mode=igraph.IN))
test(g2.vcount() == 3 and g2.ecount() == 3)

start("Testing igraph_get_edgelist");
test(g2.get_edgelist()==[(0,1), (1,2), (2,0)]);

start("Retrieving adjacency matrix");
test(g2.get_adjacency()==[[0,1,0], [0,0,1], [1,0,0]]);

start("Retrieving adjacency matrix, only upper half, undirected graph");
g=igraph.Graph(edges=[(0,1), (1,2), (2,0), (3,0)])
test(g.get_adjacency(igraph.GET_ADJACENCY_UPPER)==[[0,1,1,1], [0,0,1,0], [0,0,0,0], [0,0,0,0]]);

start("Calculating bibliographic coupling")
g=igraph.Graph(edges=[(0,1), (2,1), (2,0), (3,0)], directed=True)
test(g.bibcoupling()==[[0,0,1,0], [0,0,0,0], [1,0,0,1], [0,0,1,0]])

start("Calculating cocitation scores")
test(g.cocitation()==[[0,1,0,0], [1,0,0,0], [0,0,0,0], [0,0,0,0]])

start("Calculating outgoing shortest paths")
test(g.shortest_paths(mode=igraph.OUT)==[[0,1,4,4], [4,0,4,4], [1,1,0,4], [1,2,4,0]])

g=igraph.Graph(edges=[(0,1), (1,2), (2,0), (0,2), (3,2)], directed=True)
start("Calculating Google PageRank values in a graph")
l=g.pagerank()
expected=[1.48, 0.78, 1.57, 0.15]
diff=map(lambda x: abs(round(l[x], 2)-expected[x]), range(len(l)))
print l
test(sum(diff)<0.001)

g=igraph.Graph.Barabasi(100, 3, directed=True)
start("Randomly rewiring a Barabasi-Albert graph while preserving degrees")
ddist=[(g.degree(i, type=igraph.OUT), g.degree(i, type=igraph.IN)) for i in range(g.vcount())]
g.rewire()
ddist2=[(g.degree(i, type=igraph.OUT), g.degree(i, type=igraph.IN)) for i in range(g.vcount())]
test(ddist == ddist2)
del ddist
del ddist2

"""
g=igraph.Graph(edges=[(0,1), (1,2), (2,0), (3,4), (4,5), (5,3), (6,4)])
start("Decomposing a graph into components")
l=g.decompose()
test(len(l) == 2 and l[0].vcount() == 3 and l[0].ecount() == 3 and l[1].vcount() == 4 and l[1].ecount() == 4)
"""

section("Layout algorithms")

start("Testing layout of vertices in a circle or sphere")
g.layout_circle()
g.layout_sphere()
skip()

start("Testing layout of vertices randomly either in 2D or 3D")
g.layout_random()
g.layout_random_3d()
skip()

start("Testing layout of vertices according to the Kamada-Kawai layout")
g.layout_kamada_kawai()
g.layout_kamada_kawai_3d()
skip()

start("Testing layout of vertices according to the Fruchterman-Reingold layout")
g.layout_fruchterman_reingold()
g.layout_fruchterman_reingold_3d()
skip()

start("Testing layout of vertices according to the Fruchterman-Reingold grid layout")
g.layout_grid_fruchterman_reingold()
skip()

start("Testing layout of vertices according to the Large Graph Layout")
g.layout_lgl()
skip()

start("Testing layout of vertices according to the Reingold-Tilford tree layout")
g.layout_reingold_tilford(0)
skip()

section("Import and export functions")

start("Trying to load edge list from nonexistent file")
try:
    g=igraph.Graph.Read_Edgelist("nonexistent")
except IOError, e:
    ok()
    print "Expected exception arrived: ", e
else:
    fail()

start("Trying to load edge list from existing file")
(fd,fname)=tempfile.mkstemp()
f=os.fdopen(fd, "w")
f.write("1 2\n")
f.write("2 3\n")
f.write("3 4\n")
f.write("1 3\n")
f.close()
g=igraph.Graph.Read_Edgelist(fname, directed=False)
test(g.vcount() == 5 and g.ecount() == 4)
os.unlink(fname)

start("Trying to save graph edge list")
(fd,fname)=tempfile.mkstemp()
os.close(fd)
g.write_edgelist(fname)
os.unlink(fname)
skip()

start("Trying to load nonexistent NCOL file")
try:
    g=igraph.Graph.Read_Ncol("nonexistent")
except IOError, e:
    ok()
    print "Expected exception arrived: ", e
else:
    fail()

start("Trying to load graph from existing NCOL file")
(fd,fname)=tempfile.mkstemp()
f=os.fdopen(fd, "w")
f.write("a b 2.0\n")
f.write("b c\n")
f.write("c d 3.0\n")
f.write("a c\n")
f.close()
g=igraph.Graph.Read_Ncol(fname)
test(g.vcount() == 4 and g.ecount() == 4)
os.unlink(fname)

start("Trying to save graph edge list as NCOL file")
(fd,fname)=tempfile.mkstemp()
os.close(fd)
g.write_ncol(fname)
os.unlink(fname)
skip()

start("Trying to save graph edge list as NCOL file without weights and names")
(fd,fname)=tempfile.mkstemp()
os.close(fd)
g.write_ncol(fname, None, None)
os.unlink(fname)
skip()

start("Trying to load Pajek file")
g=igraph.Graph.Read_Pajek("examples/simple/LINKS.NET")
test(g.vcount() == 4 and g.ecount() == 7)

start("Trying to load GraphML file")
g=igraph.Graph.Read_GraphML("examples/simple/test.gxl")
ok()
#test(g.vcount() == 11 and g.ecount() == 22)

start("Trying to load GraphML file as undirected")
g=igraph.Graph.Read_GraphML("examples/simple/test.gxl", directed=False)
ok()
#test(g.vcount() == 11 and g.ecount() == 12)

start("Trying to load GraphML file with invalid index")
try:
    g=igraph.Graph.Read_GraphML("examples/simple/test.gxl", index=2)
except igraph.InternalError, e:
    ok()
    print "Expected exception arrived:", e
else:
    fail()

"""
g=igraph.Graph.Full(5)
start("Trying to save graph not having name/weight attributes as NCOL file")
(fd,fname)=tempfile.mkstemp()
os.close(fd)
g.write_ncol(fname)
os.unlink(fname)
skip()
"""

start("Trying to load nonexistent LGL file")
try:
    g=igraph.Graph.Read_Lgl("nonexistent")
except IOError, e:
    ok()
    print "Expected exception arrived: ", e
else:
    fail()

start("Trying to load graph from existing LGL file")
(fd,fname)=tempfile.mkstemp()
f=os.fdopen(fd, "w")
f.write("# a\n")
f.write("b 2.0\n")
f.write("c\n")
f.write("# b\n")
f.write("c\n")
f.write("# c\n")
f.write("d 3.0\n")
f.write("# e\n")
f.close()
g=igraph.Graph.Read_Lgl(fname)
test(g.vcount() == 5 and g.ecount() == 4)
os.unlink(fname)

start("Trying to save graph edge list as LGL file")
(fd,fname)=tempfile.mkstemp()
os.close(fd)
g.write_lgl(fname)
os.unlink(fname)
skip()

start("Trying to save graph edge list as LGL file without weights and names and isolated vertices")
(fd,fname)=tempfile.mkstemp()
os.close(fd)
g.write_lgl(fname, None, None, False)
os.unlink(fname)
skip()

section("Vertex sequence of a graph")

g=igraph.Graph.Full(5)

start("Testing whether it's possible to reference the vertex sequence")
test(g.vs)

start("Testing whether it's possible to reference a vertex in the vertex sequence")
test(g.vs[0] != None)

start("Testing whether it's possible to get the sequence length")
test(len(g.vs) == 5)

start("Iterating through the sequence")
for v in g.vs:
    str(v)
ok()

start("Weak reference testing: exporting the vertex sequence to a variable and destroying the graph")
vs=g.vs
del g
try:
    vs[0]
except TypeError, e:
    print "Expected exception arrived:", e
    ok()
else:
    fail()
del vs

start("Weak reference testing: referencing a vertex of the destroyed graph")
try:
    print v
except TypeError, e:
    print "Expected exception arrived:", e
    ok()
else:
    fail()
del v

section("Edge sequence of a graph")

g=igraph.Graph.Full(5)

start("Testing whether it's possible to reference the edge sequence")
test(g.es)

start("Testing whether it's possible to reference an edge in the edge sequence")
test(g.es[0] != None)

start("Testing whether it's possible to get the sequence length")
test(len(g.es) == 10)

start("Iterating through the sequence")
for ed in g.es:
    str(ed)
ok()

start("Getting the source/target node of an edge")
test(g.es[0].source == 0 and g.es[0].target == 1)

start("Weak reference testing: exporting the edge sequence to a variable and destroying the graph")
es=g.es
del g
try:
    es[0]
except TypeError, e:
    print "Expected exception arrived:", e
    ok()
else:
    fail()
del es

start("Weak reference testing: referencing an edge of the destroyed graph")
try:
    print ed
except TypeError, e:
    print "Expected exception arrived:", e
    ok()
else:
    fail()
del ed

section("Graph, vertex and edge attributes")

g=igraph.Graph.Full(5)

l=[g, g.vs[0], g.es[0]]
for i in l:
    start("Trying to get a nonexistent attribute")
    try:
	print i['aaa']
    except KeyError, e:
	ok()
	print "Expected exception arrived: ", e
    else:
	fail()

    start("Setting some string attributes and reading them back")
    i["name"]="Test"
    i["date"]="9:37 11.10.2005"
    test(i["name"] == "Test" and i["date"] == "9:37 11.10.2005")

    start("Overwriting a string attribute")
    s="11:48 11.11.2005"
    i["date"]=s
    s=""
    test(i["date"] == "11:48 11.11.2005")

    start("Setting a numeric attribute from integer")
    i["size"]=2
    test(i["size"] == 2)

    start("Setting a numeric attribute from float")
    i["size"]=2.5
    test(i["size"] == 2.5)

    start("Getting list of attributes")
    l=i.attributes()
    l.sort()
    test(l == ["date", "name", "size"])

    start("Overwriting a numeric attribute with a string")
    try:
	i["size"]="6000 bytes"
    except TypeError, e:
	del i["size"]
	i["size"]="6000 bytes"
    test(i["size"] == "6000 bytes")
del i
del l

section("Graph operators")

start("Disjoint union of graphs")
g=igraph.Graph.Tree(10, 3) + igraph.Graph(edges=[(0, 1), (2, 1)])
l=g.get_edgelist()
l.sort()
test(l == [(1, 0), (2, 0), (3, 0), (4, 1), (5, 1), (6, 1), (7, 2), (8, 2), (9, 2), (10, 11), (12, 11)])

start("Disjoint union of graph and number (should fail)")
try:
    g=g + 2
except TypeError:
    ok()
else:
    fail()

start("In-place disjoint union")
g = igraph.Graph.Tree(10, 3)
g += igraph.Graph(edges=[(0, 1), (2, 1)])
l=g.get_edgelist()
l.sort()
test(l == [(1, 0), (2, 0), (3, 0), (4, 1), (5, 1), (6, 1), (7, 2), (8, 2), (9, 2), (10, 11), (12, 11)])

start("Multiple disjoint union")
g = igraph.Graph.Tree(10, 3).disjoint_union([igraph.Graph.Tree(8, 2), igraph.Graph.Full(5)])
skip()

start("Union of graphs")
g=igraph.Graph.Tree(10, 3) | igraph.Graph(edges=[(0, 1), (2, 1)])
l=g.get_edgelist()
l.sort()
test(l == [(0, 1), (0, 2), (0, 3), (1, 2), (1, 4), (1, 5), (1, 6), (2, 7), (2, 8), (2, 9)])

start("Union of graph and number (should fail)")
try:
    g=g | 2
except TypeError:
    ok()
else:
    fail()

start("Multiple union")
g = igraph.Graph.Tree(10, 3).union([igraph.Graph.Tree(8, 2), igraph.Graph.Full(5)])
skip()

start("Intersection of graphs")
g=igraph.Graph.Tree(10, 3) & igraph.Graph(edges=[(0, 1), (2, 1), (5, 6), (4, 1)])
l=g.get_edgelist()
l.sort()
test(l == [(0, 1), (1, 4)])

start("Difference of graphs")
g=igraph.Graph.Tree(10, 3) - igraph.Graph(edges=[(0, 1), (2, 1), (5, 6), (4, 1)])
l=g.get_edgelist()
l.sort()
test(l == [(0, 2), (0, 3), (1, 5), (1, 6), (2, 7), (2, 8), (2, 9)])

start("Complementer of graph")
l=(~igraph.Graph.Tree(3, 2, igraph.TREE_OUT)).get_edgelist()
l.sort()
test(l == [(1, 0), (1, 2), (2, 0), (2, 1)])

start("Complementer of graph with loop edges")
l=igraph.Graph.Tree(3, 2, igraph.TREE_OUT).complementer(True).get_edgelist()
l.sort()
test(l == [(0, 0), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)])

start("Composition of two graphs")
l=igraph.Graph(4, [(0, 1), (1, 2)]).compose(igraph.Graph(4, [(2, 3), (2, 1)])).get_edgelist()
l.sort()
test(l == [(0, 2), (1, 1), (1, 3), (2, 2)])

section("Miscellaneous functions")

vs = [[3, 2], [5, 1], [4, 4], [6, 4], [4, 3], [2, 5], [1, 3], [2, 4], [6, 3], [9, 2]];
start("Convex hull test #1")
test(igraph.convex_hull(vs) == [1, 6, 5, 3, 9])
start("Convex hull test #2")
print(igraph.convex_hull(vs, True))
test(igraph.convex_hull(vs, True) == [[5, 1], [1, 3], [2, 5], [6, 4], [9, 2]])

results()

print "\nTesting for possible leaks..."

# This is a list containing all possible leaky objects.
# Don't EVER use a list here. Lists don't get deallocated immediately even
# when using their __del__ method. Python keeps track of the last 80 lists
# and doesn't free them in order to speed up new list allocation. So if
# we use a real list here, our graphs won't get deallocated even when we
# __del__ them because of their references in the list. The real
# deallocation will take place only when Python destroys the list at
# last, which mostly happens at the end of the testing framework, after
# the last print message.
rootnames=["g", "g2", "degs", "neis"]
roots=[degs, neis]
objs=get_all_objects(roots)
test_leaks(roots)

print "Invoking garbage collector..."
gc.set_debug(gc.DEBUG_UNCOLLECTABLE | gc.DEBUG_INSTANCES | gc.DEBUG_SAVEALL)
x=gc.collect()
if x>0:
    print x, "uncollectable objects found."
    for o in gc.garbage:
	rs=gc.get_referrers(o)
	print type(o), hex(id(o))
	# print type(o), hex(id(o)), "referred by", map(lambda x: (type(x), hex(id(x))), rs)
	# print type(o), hex(id(o)), "referred by", map(lambda x: repr(x), rs)
	del rs
print "Garbage collected."

for o in rootnames:
    exec "del "+o
del roots
del rootnames

print "Testing finished. No deallocation should happen in the igraph module from now on."
