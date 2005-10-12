import sys;
import gc;
sys.path.append('interfaces/python/.libs');
import igraph;
from time import clock;

longtests=False;

def start(msg):
    sys.stdout.write(msg+"... ");
    sys.stdout.flush();
    
def ok():
    global passed;
    passed=passed+1;
    print "ok.";

def fail():
    global failed;
    failed=failed+1;
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
    print "Passed: %d, failed: %d, skipped: %d, percentage: %f" % (passed, failed, skipped, float(passed)/(passed+failed)*100)
    if failed>0:
	print "SOME TESTS FAILED!"
	print "Check the implementation of the Python module."
    else:
	print "Everything went OK."

def test_leaks(objects):
    for o in objects:
	rc=sys.getrefcount(o)
	refs=len(gc.get_referrers(o))
	if rc>refs+1:
	    print hex(id(o)), type(o), sys.getrefcount(o), len(gc.get_referrers(o))

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
	
passed=0; failed=0; skipped=0;

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

start("Testing neighbor retrieval");
neis=map(lambda x: g.neighbors(x), range(5))
for n in neis: n.sort()
test(neis == [[1,3,4], [0,2,4], [1,2,2,3], [0,2,4], [0,1,3]])

start("Testing edge deletion")
g.delete_edges([(0, 4), (2, 2)]);
test(g.vcount() == 5 and g.ecount() == 6 and g.degree() == [2, 3, 2, 3, 2])

start("Calculating graph diameter")
test(g.diameter(True, True) == 2)

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

start("Generating a ternary directed tree with 40 nodes")
g=igraph.Graph.Tree(40, 3, type=igraph.TREE_OUT)
test(g.vcount() == 40 and g.ecount() == 39)

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

start("Calculating shortest paths from/to node 1 in the graph")
test(g.get_shortest_paths(4, mode=igraph.IN) == [[4,3,1,0], [4,3,1], [4,3,1,2], [4,3], [4], []])

start("Adding some more edges (to make sense for a spanning tree calculation)")
g.add_edges([(0,2), (4, 5), (0, 4), (2, 5)])
test(g.vcount() == 6 and g.ecount() == 9)

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

start("Testing layout of vertices in a circle")
print g.layout_circle()
skip()

start("Testing layout of vertices randomly")
print g.layout_random()
skip()

results()
test_leaks([g, g2, degs, neis, comps])

gc.set_debug(gc.DEBUG_LEAK)
gc.collect()
