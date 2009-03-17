.. include:: include/global.rst

.. Tutorial

========
Tutorial
========

This chapter contains a short overview of |igraph|'s capabilities. It is highly recommended
to read it at least once if you are new to |igraph|. I assume that you have already installed
|igraph|; if you did not, see :ref:`installing-igraph` first. Familiarity with the Python
language is also assumed; if this is the first time you are trying to use Python, there are
many good Python tutorials on the Internet to get you started. Mark Pilgrim's
`Dive Into Python <http://www.diveintopython.org>`_ is one that I personally suggest.
If this is the first time you ever try to use a programming language,
`A Byte of Python <http://swaroopch.com/notes/Python>`_ is even better. If
you already have a stable programming background in other languages and you just want a
quick overview of Python, `Learn Python in 10 minutes
<http://www.poromenos.org/tutorials/python>`_ is probably your best bet.


Starting |igraph|
=================

|igraph| is a Python module, hence it can be imported exactly the same way as any other
ordinary Python module at the Python prompt::

  $ python
  Python 2.5.1 (r251:54863, Apr 15 2008, 22:57:26)
  [GCC 4.0.1 (Apple Inc. build 5465)] on darwin
  Type "help", "copyright", "credits" or "license" for more information.
  >>> import igraph

This imports |igraph|'s objects and methods inside an own namespace called :mod:`igraph`. Whenever
you would like to call any of |igraph|'s methods, you will have to provide the appropriate
namespace-qualification. E.g., to check which |igraph| version you are using, you could do the
following:

>>> import igraph
>>> print igraph.__version__
0.6

Another way to make use of |igraph| is to import all its objects and methods into the main
Python namespace (so you do not have to type the namespace-qualification every time).
This is fine as long as none of your own objects and methods do not conflict with the ones
provided by |igraph|:

>>> from igraph import *

The third way to start |igraph| is to simply call the startup script that was supplied with
the |igraph| package you installed. Not too surprisingly, the script is called :command:`igraph`,
and provided that the script is on your path in the command line of your operating system
(which is almost surely the case on Linux and OS X), you can simply type :command:`igraph` at the
command line. Windows users will find the script inside the :file:`scripts` subdirectory of Python
and you may have to add it manually to your path in order to be able to use the script from
the command line without typing the whole path.

When you start the script, you will see something like this::

  $ igraph
  No configuration file, using defaults
  igraph 0.6 running inside Python 2.5.1 (r251:54863, Apr 15 2008, 22:57:26)
  Type "copyright", "credits" or "license" for more information.
  >>>

The command-line startup script imports all of |igraph|'s methods and objects into the main
namespace, so it is practically equivalent to ``from igraph import *``. The difference between
the two approaches (apart from saving some typing) is that the command-line script checks
whether you have any of Python's more advanced shells installed and uses that instead of the
standard Python shell. Currently the module looks for `IPython <http://ipython.scipy.org>`_ and
IDLE (the Tcl/Tk-based graphical shell supplied with Python). If neither IPython nor IDLE is
installed, the startup script launches the default Python shell. You can also modify the
order in which these shells are searched by tweaking |igraph|'s configuration file
(see :ref:`configuring-igraph`).

In general, it is advised to use the command line startup script when using |igraph|
interactively (i.e., when you just want to quickly load or generate some graphs, calculate
some basic properties and save the results somewhere). For non-disposable graph analysis
routines that you intend to re-run from time to time, you should write a script separately
in a ``.py`` source file and import |igraph| using one of the above methods at the start of
the script, then launch the script using the Python interpreter.

From now on, every example in the documentation will assume that |igraph|'s objects and
methods are imported into the main namespace (i.e., we used ``from igraph import *``
instead of ``import igraph``). If you let |igraph| take its own namespace, please adjust
all the examples accordingly.


Creating a graph from scratch
=============================

Assuming that you have started |igraph| successfully, it is time to create your first
|igraph| graph. This is pretty simple:

>>> g = Graph(1)

The above statement created an undirected graph with a single vertex and assigned it to
the variable `g`. To confirm that it's really an |igraph| graph, we can
print it:

>>> g
<igraph.Graph object at 0x4c87a0>

This tells us that `g` is an instance of |igraph|'s :class:`Graph` class and that it is currently
living at the memory address ``0x4c87a0`` (the exact output will almost surely be different
for your platform). To obtain a more user-friendly output, we can try to print the graph
using Python's :keyword:`print` statement:

>>> print(g)
Undirected graph (|V| = 1, |E| = 0)

This is not too exciting so far; a graph with a single vertex and no edges is not really useful
for us. Let's add some vertices first!

>>> g.add_vertices(2)
<igraph.Graph object at 0x4c87a0>

:meth:`Graph.add_vertices` (i.e., the :meth:`~Graph.add_vertices` method of the :class:`Graph`
class) adds the given number of vertices to the graph and returns the graph itself (hence the
output; you can see that the memory address of the :class:`Graph` object that was returned
is exactly the same).

Now our graph has three vertices but no edges, so let's add some edges as well! You can
add edges by calling :meth:`Graph.add_edges` - but in order to add edges, you have to refer to
existing vertices somehow. |igraph| uses integer vertex IDs starting from zero, thus the
first vertex of your graph has index zero, the second vertex has index 1 and so on.
Edges are specified by pairs of integers, so ``[(0,1), (1,2)]`` denotes a list of two
edges: one between the first and the second, and the other one between the second and the
third vertices of the graph. Passing this list to :meth:`Graph.add_edges` adds these two edges
to your graph:

>>> g.add_edges([(0,1), (1,2)])
<igraph.Graph object at 0x4c87a0>

:meth:`~Graph.add_edges` is clever enough to figure out what you want to do in most of the
cases: if you supply a single pair of integers, it will automatically assume that you want
to add a single edge. However, if you try to add edges to vertices with invalid IDs (i.e.,
you try to add an edge to vertex 5 when you only have three edges), you will get an
exception:

>>> g.add_edges((5, 0))
Traceback (most recent call last):
  File "<stdin>", line 6, in <module>
igraph.core.InternalError: Error at ../../src/type_indexededgelist.c:245: cannot add edges, invalid vertex id

Most |igraph| functions will raise an :exc:`igraph.core.InternalError` if
something goes wrong. The message corresponding to the exception gives you a
short textual explanation of what went wrong (``cannot add edges, invalid
vertex id``) along with the corresponding line in the C source where the error
occurred. The exact filename and line number may not be too informative to you,
but it is invaluable for |igraph| developers if you think you found an error in
|igraph| and you want to report it.

You may be wondering why it is useful to return the graph itself when adding
vertices or edges. The reason is that you can conveniently chain your calls to
:meth:`~Graph.add_vertices` and :meth:`~Graph.add_edges`. Let us go on with our
graph ``g`` and add some more vertices and edges to it:

>>> g.add_edges((2,0)).add_vertices(3).add_edges([(2,3),(3,4),(4,5),(5,3)])
<igraph.Graph object at 0x4c87a0>
>>> print g
Undirected graph (|V| = 6, |E| = 7)

Now, this is better. We have an undirected graph with six vertices and seven edges.
Edges also have IDs, similarly to vertices; they also start from zero and edges that
were added later have higher IDs than edges that were added earlier. Vertex and
edge IDs are always *continuous*, and a direct consequence of this fact is that
if you happen to delete an edge, chances are that some (or all) of the edges will
be renumbered. Moreover, if you delete a vertex, even the vertex IDs will change.
Edges can be deleted by :meth:`~Graph.delete_edges` and it requires a list of edge IDs
to be deleted (or a single edge ID). Vertices can be deleted by :meth:`~Graph.delete_vertices`
and you may have already guessed that it requires a list of vertex IDs to be deleted
(or a single vertex ID). If you do not know the ID of an edge you wish to delete,
but you know the IDs of the vertices at its two endpoints, you can use :meth:`~Graph.get_eid`
to get the edge ID. Remember, all these are *methods* of the :class:`Graph` class and
you must call them on the appropriate :class:`Graph` instance!

>>> g.get_eid(2,3)
3
>>> g.delete_edges(3)
<igraph.Graph object at 0x4c87a0>
>>> summary(g)
6 nodes, 6 edges, undirected
Number of components: 2
Diameter: 1
Density: 0.4000
Average path length: 1.0000

:meth:`summary` is a new command that you haven't seen before; it is a member of |igraph|'s
own namespace and it can be used to get an overview of a given graph object. It lists the
number of nodes and edges, checks whether the graph is directed, counts the connected components,
calculates the graph diameter, the edge density and the average path lengths. All of these
informations can be calculated separately by the appropriate methods of :class:`Graph` of course;
we will talk about these later in the reference manual. In general, :meth:`summary` is primarily
meant for smaller graph objects as calculating some of these properties on a really large graph
would take a lot of time.


Generating graphs
=================

|igraph| includes a large set of graph generators which can be divided into two groups:
deterministic and stochastic graph generators. Deterministic generators produce the same
graph if you call them with exactly the same parameters, while stochastic generators
produce a different graph every time. Deterministic generators include methods for
creating trees, regular lattices, rings, extended chordal rings, several famous graphs
and so on, while stochastic generators are used to create Erdős-Rényi random networks,
Barabási-Albert networks, geometric random graphs and such. |igraph| has too many
generators to cover them all in this tutorial, so we will only show an example for both:

>>> g = Graph.Tree(127, 2)
>>> summary(g)
127 nodes, 126 edges, undirected
Number of components: 1
Diameter: 12
Density: 0.0157
Average path length: 8.3510

:meth:`Graph.Tree` generates a regular tree graph. The one that we generated has 127
vertices and each vertex (apart from the root) has two children (and of course one
parent). No matter how many times you call :meth:`Graph.Tree`, the generated graph will
always be the same if you use the same parameters:

>>> g2 = Graph.Tree(127, 2)
>>> g2.get_edgelist() == g.get_edgelist()
True

The above code snippet also shows you that the :meth:`~Graph.get_edgelist()` method
of :class:`Graph` graph objects return a list that contains pairs of integers, one for
each edge. The first member of the pair is the source vertex ID and the second member
is the target vertex ID of the corresponding edge. This list is too long, so let's
just print the first 10 elements!

>>> g2.get_edgelist()[0:10]
[(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6), (3, 7), (3, 8), (4, 9), (4, 10)]

Let's do the same with a stochastic generator!

>>> g = Graph.GRG(100, 0.2)
>>> summary(g)
100 nodes, 524 edges, undirected
Number of components: 1
Diameter: 9
Density: 0.1059
Average path length: 3.7701

:meth:`Graph.GRG` generates a geometric random graph: *n* points are chosen randomly and
uniformly inside the unit square and pairs of points closer to each other than a predefined
distance *d* are connected by an edge. In our case, *n* is 100 and *d* is 0.2. Due to
the random nature of the algorithm, chances are that the exact graph you got is different
from the one that was generated when I wrote this tutorial, hence the values above in the
summary will not match the ones you got. This is normal and expected. Even if you generate
two geometric random graphs on the same machine, they will be different for the same parameter
set:

>>> g2 = Graph.GRG(100, 0.2)
>>> g.get_edgelist() == g2.get_edgelist()
False
>>> g.isomorphic(g2)
False

:meth:`~Graph.isomorphic()` tells you whether two graphs are isomorphic or not. In general,
it might take quite a lot of time, especially for large graphs, but in our case, the
answer can quickly be given by checking the degree sequences of the two graphs.


Setting and retrieving attributes
=================================

|igraph| uses vertex and edge IDs in its core. These IDs are integers, starting from zero,
and they are always continuous at any given time instance during the lifetime of the graph.
This means that whenever vertices and edges are deleted, a large set of edge and possibly
vertex IDs will be renumbered to ensure the continuiuty. Now, let us assume that our graph
is a social network where vertices represent people and edges represent social connections
between them. One way to maintain the association between vertex IDs and say, the corresponding
names is to have an additional Python list that maps from vertex IDs to names. The drawback
of this approach is that this additional list must be maintained in parallel to the
modifications of the original graph. Luckily, |igraph| knows the concept of *attributes*,
i.e., auxiliary objects associated to a given vertex or edge of a graph, or even to the
graph as a whole. Every |igraph| :class:`Graph`, vertex and edge behaves as a standard
Python dictionary in some sense: you can add key-value pairs to any of them, with the key
representing the name of your attribute (the only restriction is that it must be a string)
and the value representing the attribute itself.

.. warning::
   Attributes can be arbitrary Python objects, but if you are saving graphs to a
   file, only string and numeric attributes will be kept. See the :mod:`pickle` module in
   the standard Python library if you are looking for a way to save other attribute types.
   You can either pickle your attributes individually, store them as strings and save them,
   or you can pickle the whole :class:`Graph` if you know that you want to load the graph
   back into Python only.


Let us create a simple imaginary social network the usual way by hand.

>>> g = Graph([(0,1), (0,2), (2,3), (3,4), (4,2), (2,5), (5,0), (6,3), (5,6)])

Now, let us assume that we want to store the names, ages and genders of people in this network as
vertex attributes, and for every connection, we want to store whether this is an informal
friendship tie or a formal tie. Every :class:`Graph` object contains two special members
called :attr:`~Graph.vs` and :attr:`~Graph.es`, standing for the sequence of all vertices
and all edges, respectively. If you try to use :attr:`~Graph.vs` or :attr:`~Graph.es` as
a Python dictionary, you will manipulate the attribute storage area of the graph:

>>> g.vs
<igraph.VertexSeq object at 0x1b23b90>
>>> g.vs["name"] = ["Alice", "Bob", "Claire", "Dennis", "Esther", "Frank", "George"]
>>> g.vs["age"] = [25, 31, 18, 47, 22, 23, 50]
>>> g.vs["gender"] = ["f", "m", "f", "m", "f", "m", "m"]
>>> g.es["formal"] = [False, False, True, True, True, False, True, False, False]

Whenever you use :attr:`~Graph.vs` or :attr:`~Graph.es` as a dictionary, you are assigning
attributes to *all* vertices/edges of the graph. However, you can simply alter the attributes
of vertices and edges individually by *indexing* :attr:`~Graph.vs` or :attr:`~Graph.es`
with integers as if they were lists (remember, they are sequences, they contain all the
vertices or all the edges). When you index them, you obtain a :class:`Vertex` or
:class:`Edge` object, which refers to (I am sure you already guessed that) a single vertex
or a single edge of the graph. :class:`Vertex` and :class:`Edge` objects can also be used
as dictionaries to alter the attributes of that single vertex or edge:

>>> g.es[0]
igraph.Edge(<igraph.Graph object at 0x4c87a0>,0,{'formal': False})
>>> g.es[0].attributes()
{'formal': False}
>>> g.es[0]["formal"] = True
igraph.Edge(<igraph.Graph object at 0x4c87a0>,0,{'formal': True})

The above snippet illustrates that indexing an :class:`EdgeSeq` object returns
:class:`Edge` objects; the representation above shows the graph the object belongs to,
the edge ID (zero in our case) and the dictionary of attributes assigned to that edge.
:class:`Edge` objects have some useful attributes, too: the :attr:`~Edge.source` property
gives you the source vertex of that edge, :attr:`~Edge.target` gives you the target vertex,
:attr:`~Edge.index` gives you the corresponding edge ID, :attr:`~Edge.tuple` gives you a
tuple containing the source and target vertices and :meth:`~Edge.attributes` gives you
a dictionary containing the attributes of this edge. :class:`Vertex` instances only have
:attr:`~Vertex.index` and :meth:`~Vertex.attributes`.

Since :attr:`Graph.es` always represents all the edges in a graph, indexing it by
*i* will always return the edge with ID *i*, and of course the same applies
to :attr:`Graph.vs`. However, keep in mind that an :class:`EdgeSeq` object *in general*
does not necessarily represent the whole edge sequence of a graph;
:ref:`later in this tutorial <querying_vertices_and_edges>`
we will see methods that can filter :class:`EdgeSeq` objects and return other
:class:`EdgeSeq` objects that are restricted to a subset of edges, and of course the same
applies to :class:`VertexSeq` objects. But before we dive into that, let's see how we
can assign attributes to the whole graph. Not too surprisingly, :class:`Graph` objects
themselves can also behave as dictionaries:

>>> g["date"] = "2009-01-10"
>>> print g["date"]
2009-01-10

Finally, it should be mentioned that attributes can be deleted by the Python keyword
:keyword:`del` just as you would do with any member of an ordinary dictionary:

>>> g.vs[3]["foo"] = "bar"
>>> g.vs["foo"]
[None, None, None, 'bar', None, None, None]
>>> del g.vs["foo"]
>>> g.vs["foo"]
Traceback (most recent call last):
  File "<stdin>", line 25, in <module>
KeyError: 'Attribute does not exist'


Structural properties of graphs
===============================

Besides the simple graph and attribute manipulation routines described above,
|igraph| provides a large set of methods to calculate various structural properties
of graphs. It is beyond the scope of this tutorial to document all of them, hence
this section will only introduce a few of them for illustrative purposes.
We will work on the small social network we built in the previous section.

Probably the simplest property one can think of is the :dfn:`vertex degree`. The
degree of a vertex equals the number of edges adjacent to it. In case of directed
networks, we can also define :dfn:`in-degree` (the number of edges pointing towards
the vertex) and :dfn:`out-degree` (the number of edges originating from the vertex).
|igraph| is able to calculate all of them using a simple syntax:

>>> g.degree()
[3, 1, 4, 3, 2, 3, 2]

If the graph was directed, we would have been able to calculate the in- and out-degrees
separately using ``g.degree(type="in")`` and ``g.degree(type="out")``. You can
also pass a single vertex ID or a list of vertex IDs to :meth:`~Graph.degree` if you
want to calculate the degrees for only a subset of vertices:

>>> g.degree(6)
2
>>> g.degree([2,3,4])
[4, 3, 2]

This calling convention applies to most of the structural properties |igraph| can
calculate. For vertex properties, the methods accept a vertex ID or a list of vertex IDs
(and if they are omitted, the default is the set of all vertices). For edge properties,
the methods accept a single edge ID or a list of edge IDs. Instead of a list of IDs,
you can also supply a :class:`VertexSeq` or an :class:`EdgeSeq` instance appropriately.
Later in the :ref:`next chapter <querying_vertices_and_edges>`, you will learn how to
restrict them to exactly the vertices or edges you want.

.. note::

  For some measures, it does not make sense to calculate them only for a few vertices
  or edges instead of the whole graph, as it would take the same time anyway. In this
  case, the methods won't accept vertex or edge IDs, but you can still restrict the
  resulting list later using standard list indexing and slicing operators. One such
  example is eigenvector centrality (:meth:`Graph.evcent()`).

Besides degree, |igraph| includes built-in routines to calculate many other centrality
properties, including node and edge betweenness (:meth:`Graph.betweenness`,
:meth:`Graph.edge_betweenness`) or Google's PageRank (:meth:`Graph.pagerank`)
just to name a few. Here we just illustrate edge betweenness:

>>> g.edge_betweenness()
[6.0, 6.0, 4.0, 2.0, 4.0, 3.0, 4.0, 3.0. 4.0]

Now we can also figure out which connections have the highest betweenness centrality
with some Python magic:

>>> ebs = g.edge_betweenness()
>>> max_eb = max(ebs)
>>> [g.es[idx].tuple for idx, eb in enumerate(ebs) if eb == max_eb]
[(0, 1), (0, 2)]

.. _querying_vertices_and_edges:

Querying vertices and edges based on attributes
===============================================

Imagine that in a given social network, you would like to find out who has the largest
degree or betweenness centrality. You can do that with the tools presented so far and
some basic Python knowledge, but since it is a common task to select vertices and edges
based on attributes or structural properties, |igraph| gives you an easier way to do that:

>>> g.vs.select(_degree = g.maxdegree())["name"]
["Alice", "Bob"]

The syntax may seem a little bit awkward for the first sight, so let's try to interpret
it step by step. :meth:`~VertexSeq.select` is a method of :class:`VertexSeq` and its
sole purpose is to filter a :class:`VertexSeq` based on the properties of individual
vertices. The way it filters the vertices depends on its positional and keyword
arguments and positional arguments (the ones without an explicit name like
``_degree`` above) are always processed before keyword arguments as follows:

- If the first positional argument is :keyword:`None`, an empty sequence (containing no
  vertices) is returned:

  >>> seq = g.vs.select(None)
  >>> len(seq)
  0

- If the first positional argument is a callable object (i.e., a function, a bound
  method or anything that behaves like a function), the object will be called for
  every vertex that's currently in the sequence. If the function returns :keyword:`True`,
  the vertex will be included, otherwise it will be excluded:

  >>> graph = Graph.Full(10)
  >>> only_odd_vertices = graph.select(lambda x: x % 2 == 1)
  >>> len(only_odd_vertices)
  5

- If the first positional argument is an iterable (i.e., a list, a generator or
  anything that can be iterated over), it *must* return integers and these integers
  will be considered as indices into the current vertex set (which is *not* necessarily
  the whole graph). Only those vertices that match the given indices will be included
  in the filtered vertex set. Floats, strings, invalid vertex IDs will silently be
  ignored:

  >>> seq = graph.vs.select([2, 3, 7])
  >>> len(seq)
  3
  >>> [v.index for v in seq]
  [2, 3, 7]
  >>> seq = seq.select([0, 2])         # filtering an existing vertex set
  >>> [v.index for v in seq]
  [2, 7]
  >>> seq = graph.vs.select([2, 3, 7, "foo", 3.5])
  >>> len(seq)
  3

- If the first positional argument is an integer, all remaining arguments are also
  expected to be integers and they are interpreted as indices into the current vertex
  set. This is just syntactical sugar, you could achieve an equivalent effect by
  passing a list as the first positional argument, but this way you can omit the
  square brackets:

  >>> seq = graph.vs.select(2, 3, 7)
  >>> len(seq)
  3

Keyword arguments can be used to filter the vertices based on their attributes
or their structural properties. The name of each keyword argument should consist
of at most two parts: the name of the attribute or structural property and the
filtering operator. The operator can be omitted; in that case, we automatically
assume the equality operator. The possibilities are as follows (where
*name* denotes the name of the attribute or property):

# TODO: table here

- ``name_eq``: the attribute/property value must be equal to the value of the
  keyword argument

- ``name_ne``: the attribute/property value must not be equal to the value of the
  keyword argument

- ``name_lt``: the attribute/property value must be less than the value of the
  keyword argument

- ``name_le``: the attribute/property value must be less than or equal to the
  value of the keyword argument

- ``name_gt``: the attribute/property value must be greater than the value of the
  keyword argument

- ``name_ge``: the attribute/property value must be greater than or equal to the
  value of the keyword argument

- ``name_in``: the attribute/property value must be included in the
  value of the keyword argument, which must be a sequence in this case
  
- ``name_notin``: the attribute/property value must not be included in the value
  of the the keyword argument, which must be a sequence in this case

For instance, the following command gives you people younger than 30 years in
our imaginary social network:

>>> g.vs.select(age_lt=30)

.. note:
   Due to the syntactical constraints of Python, you cannot use the admittedly
   simpler syntax of ``g.vs.select(age < 30)`` as only the equality operator is
   allowed to appear in an argument list in Python.

To save you some typing, you can even omit the :meth:`~VertexSeq.select` method if
you wish:

>>> g.vs(age_lt=30)

Theoretically, it can happen that there exists an attribute and a structural property
with the same name (e.g., you could have a vertex attribute named ``degree``). In that
case, we would not be able to decide whether the user meant ``degree`` as a structural
property or as a vertex attribute. To resolve this ambiguity, structural property names
*must* always be preceded by an underscore (``_``) when used for filtering. For example, to
find vertices with degree larger than 2:

>>> g.vs(_degree_gt=2)



Layouts and plotting
====================

TODO

>>> layout = g.layout("kamada_kawai")
>>> plot(g)
>>> plot(g, "graph.png")
>>> plot(g, "graph.png", vertex_color = g.vs["colors"])


|igraph| and the outside world
==============================

TODO

>>> g = load("karate.net")
>>> g.write_graphml("karate.graphml")
>>> g.save("karate.graphml")


Where to go next
================


