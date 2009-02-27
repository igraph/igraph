.. include:: include/global.rst

.. Tutorial

========
Tutorial
========

This chapter contains a short overview of |igraph|'s capabilities. It is highly recommended
to read it at least once if you are new to |igraph|.

Starting |igraph|
=================

|igraph| is a Python module, hence it can be imported exactly the same way as any other
ordinary Python module:

>>> import igraph

This imports |igraph|'s objects and methods inside an own namespace called `igraph`. Whenever
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
the |igraph| package you installed. Not too surprisingly, the script is called `igraph`,
and providing that the script is on your path in the command line of your operating system
(which is almost surely the case on Linux and OS X), you can simply type `igraph` at the
command line. Windows users will find the script inside the ``scripts`` subdirectory of Python
and you may have to add it manually to your path in order to be able to use the script from
the command line without typing the whole path.

When you start the script, you will see something like this::

  $ igraph
  Using configuration from /home/ntamas/.igraphrc
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
(see `Configuring igraph`_).

In general, it is advised to use the command line startup script when using |igraph|
interactively (i.e., when you just want to quickly load or generate some graphs, calculate
some basic properties and saving the results somewhere). For non-disposable graph analysis
routines that you intend to re-run from time to time, you should write a script separately
in a `.py` source file and import |igraph| using one of the above methods at the start of
the script, then launch the script using the Python interpreter.

From now on, every example in the documentation will assume that |igraph|'s objects and
methods are imported into the main namespace (i.e., we used ``from igraph import *``
instead of ``import igraph``). If you let |igraph| take its own namespace, please adjust
all the examples accordingly.


Creating a graph from scratch
=============================


