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
if you would like to exploit |igraph|'s functionality in Python, you must install
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

1. Get the latest PyCairo for Windows installer from
   http://ftp.gnome.org/pub/GNOME/binaries/win32/pycairo/1.8. Make sure you
   grab the one that matches your Python version. At the time of writing,
   the above folder contained installers for Python 2.6 only. You may
   also try and go one level up, then down then 1.4 subfolder -- these are
   older versions, but they work with Python 2.5 and Python 2.6 as well.

2. Install PyCairo using the installer. The installer extracts the necessary
   files into ``Lib\site-packages\cairo`` within the folder where Python is
   installed. Unfortunately there are some extra DLLs which are required to
   make Cairo work, so we have to get these as well.

3. Head to http://www.gtk.org/download-windows and get the binary versions of
   Cairo (``cairo_1.8.10-3_win32.zip`` at the time of writing), Fontconfig
   (``fontconfig_2.8.0-2_win32.zip``), Expat (``expat_2.0.1-1_win32.zip``),
   ``libpng`` (``libpng_1.4.0-1_win32.zip``) and ``zlib`` (``zlib_1.2.4-2_win32.zip``).
   Version numbers may vary, so be adaptive! Each ZIP file will contain a
   ``bin`` subfolder with a DLL file in it. Put the following DLLs in
   ``Lib\site-packages\cairo`` within your Python installation:

   - ``libcairo-2.dll`` (from ``cairo_1.8.10-3_win32.zip``)
   - ``libexpat-1.dll`` (from ``expat_2.0.1-1_win32.zip``)
   - ``libfontconfig-1.dll`` (from ``fontconfig_2.8.0-2_win32.zip``)
   - ``libpng14-14.dll`` (from ``libpng_1.4.0-1_win32.zip``)
   - ``zlib1.dll`` (from ``zlib_1.2.4-2_win32.zip``).

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

There is a Mac OS X installer for |igraph|'s Python interface on the
`Python Package Index <http://pypi.python.org/pypi/python-igraph>`_
which works for Intel-based Macs running OS X Leopard. The default
Python version in Leopard is Python 2.5, so the package is compiled
for this specific version. PowerPC users should compile the package
themselves (see `Compiling igraph from source`_). To test the
installed package, launch your favourite Python IDE or the default
command line interpreter and type the following:

  >>> import igraph.test
  >>> igraph.test.test()

The above commands run the bundled test cases to ensure that everything
is fine with your |igraph| installation.

Graph plotting in |igraph| on Mac OS X
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Graph plotting in |igraph| is implemented using a third-party package
called `Cairo <http://www.cairographics.org>`_. If you want to create
publication-quality plots in |igraph| on Mac OS X, you must also install
Cairo and its Python bindings. The Cairo project does not provide
pre-compiled binaries for Mac OS X, but `MacPorts <http://www.macports.org>`_
and `Fink <http://www.finkproject.org>`_ does, so you can use them to
install Cairo. The `Cairo homepage <http://www.cairographics.org>` gives
you some installation instructions. However, this is only one half of the
job, you will also need the Python bindings of Cairo from the
`PyCairo homepage <http://www.cairographics.org/pycairo>`. At the moment
there are no precompiled PyCairo packages for Mac OS X either.

TODO: detailed compilation instructions for PyCairo

|igraph| on other operating systems
-----------------------------------

Compiling |igraph| from source
==============================

Summary
=======
