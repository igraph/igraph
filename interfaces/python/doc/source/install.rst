.. include:: include/global.rst

.. Installing igraph

.. _installing-igraph:

===================
Installing |igraph|
===================

This chapter describes how to install the C core of |igraph| and its Python bindings
on various operating systems.

Which |igraph| is right for you?
================================

|igraph| is primarily a library written in C. It is *not* a standalone program, nor it is
a Python package that you can just drop on your Python path to start using it. Therefore,
if you would like to exploit |igraph|'s functionality in Python, you must compile and install
a few packages. Do not worry, though, there are precompiled packages for the major operating
systems, so you will not have to compile |igraph| from source unless you use an esoteric
operating system or you have specific requirements (i.e., adding a custom patch to |igraph|'s
C core). Precompiled packages are often called *binary packages*, while the raw source code
is usually referred to as the *source package*.

In general, you should almost always opt for the binary package unless a binary package is not
available for your platform or you have some local modifications that you want to incorporate
into |igraph|'s source. `Installation from a binary package`_ tells you how to install |igraph|
from a precompiled binary package on various platforms. `Compiling igraph from source`_ tells
you how to compile |igraph| from the source package.

Installation from a binary package
==================================

|igraph| on Windows
-------------------

There is a Windows installer for |igraph|'s Python interface on the
`Python Package Index <http://pypi.python.org/pypi/python-igraph>`_.
Download the one that is suitable for your Python version (currently
there are binary packages for Python 2.4 and Python 2.5, though it
might change in the future). To test the installed package, launch
your favourite Python IDE and type the following:

  >>> import igraph.test
  >>> igraph.test.test()

The above commands run the bundled test cases to ensure that everything
is fine with your |igraph| installation.

Graph plotting in |igraph| on Windows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Graph plotting in |igraph| is implemented using a third-party package
called `Cairo <http://www.cairographics.org>`_. If you want to create
publication-quality plots in |igraph| on Windows, you must also install
Cairo and its Python bindings. The Cairo project does not provide
pre-compiled binaries for Windows, but other projects depending on
Cairo do, so the preferred way to install Cairo on Windows along with
its Python bindings is as follows:

1. Get the Cairo DLLs from http://www.gtk.org/download-windows.
   You will need ``cairo_1.8.6-1_win32.zip`` and the binary version
   of ``libpng`` and ``zlib``. From ``cairo_1.8.6-1_win32.zip``,
   you need ``libcairo-2.dll``; from ``libpng``, you need
   ``libpng13.dll``; from ``zlib``, you need ``zlib1.dll``.
   Version numbers may vary, so be adaptive! Put these DLLs to
   somewhere where Windows can find them, e.g., into
   ``C:\Windows\System32``.

2. Get the latest PyCairo for Windows installer from
   http://ftp.gnome.org/pub/GNOME/binaries/win32/pycairo/1.4. Make
   sure you grab the one that matches your Python version.

3. Run the PyCairo installer.

Having done that, you can launch Python again and check if it worked:

  >>> from igraph import * 
  >>> g = Graph.Famous("petersen")
  >>> plot(g)

|igraph| on Linux
-----------------

|igraph| on Debian GNU/Linux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

|igraph| on RedHat Linux
^^^^^^^^^^^^^^^^^^^^^^^^

|igraph| on other Linux distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

|igraph| on Mac OS X
--------------------

|igraph| on other operating systems
-----------------------------------

Compiling |igraph| from source
==============================

Summary
=======
