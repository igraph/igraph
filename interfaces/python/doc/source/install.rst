.. include:: include/global.rst

.. Installing igraph

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
