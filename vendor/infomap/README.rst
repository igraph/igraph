.. image:: https://github.com/mapequation/infomap/actions/workflows/build.yml/badge.svg

Infomap
=======

Infomap is a network clustering algorithm based on the `Map equation`_.

For detailed documentation, see `mapequation.org/infomap`_.

For a list of recent changes, see `CHANGELOG.md`_ in the source directory.

.. _Map equation: https://www.mapequation.org/publications.html#Rosvall-Axelsson-Bergstrom-2009-Map-equation
.. _`mapequation.org/infomap`: https://www.mapequation.org/infomap
.. _`CHANGELOG.md`: https://github.com/mapequation/infomap/blob/master/CHANGELOG.md

Getting started
---------------

Infomap can be installed either from `PyPI`_ using ``pip`` or by
compiling from source.

An experimental Javascript version for browsers is available on `NPM`_.

.. _PyPI: https://pypi.org/project/infomap/

Using pip
---------

A pre-compiled version is available for macOS users.

Installing on other operating systems requires a
working ``gcc`` or ``clang`` compiler.

To install, run::

    pip install infomap


To upgrade, run::

    pip install --upgrade infomap


When the Python package is installed, an executable called
``infomap`` (with lowercase i) is available from any directory.

To get started, read `Infomap Python API`_.

.. _`Infomap Python API`: https://mapequation.github.io/infomap/python/

Using Docker
------------

There are currently two Docker images available on `Docker Hub`_.

- ``mapequation/infomap``
- ``mapequation/infomap:notebook`` based on ``jupyter/scipy-notebook``

The image ``mapequation/infomap`` can be started with

.. code-block:: bash

    docker run -it --rm \
        -v `pwd`:/data \
        mapequation/infomap
        [infomap arguments]

You can also use the supplied `docker-compose.yml`_:

.. code-block:: bash

    docker-compose run --rm infomap

The image ``mapequation/infomap:notebook`` can be started with

.. code-block:: bash

    docker run \
        -v `pwd`:/home/jovyan/work \
        -p 8888:8888 \
        mapequation/infomap:notebook \
        start.sh jupyter lab

Or similarly, using docker-compose:

.. code-block:: bash

    docker-compose up notebook

.. _`Docker Hub`: https://hub.docker.com/r/mapequation/infomap
.. _`docker-compose.yml`: https://github.com/mapequation/infomap/blob/master/docker-compose.yml

Compiling from source
---------------------

Installing Infomap from source requires a working ``gcc`` or ``clang`` compiler.

To download and compile the newest version from `Github`_, clone the repository
by running

.. code-block:: shell

    git clone git@github.com:mapequation/infomap.git
    cd infomap
    make

This creates the binary ``Infomap``, run it using::

    ./Infomap [options] network_data destination

For a list of options, run::

    ./Infomap --help

Read `the documentation`_ to learn more about the different options.

.. _Github: https://www.github.com/mapequation/infomap
.. _the documentation: https://www.mapequation.org/infomap

Npm package
-----------

An experimental Javascript web worker is available on `NPM`_.

To install it, run

.. code-block:: shell

    npm install @mapequation/infomap

.. _NPM: https://www.npmjs.com/package/@mapequation/infomap

Feedback
--------

If you have any questions, suggestions or issues regarding the software,
please add them to `GitHub issues`_.

.. _Github issues: http://www.github.com/mapequation/infomap/issues

Authors
-------

Daniel Edler, Anton Holmgren, Martin Rosvall

For contact information, see `mapequation.org/about.html`_.

.. _`mapequation.org/about.html`: https://www.mapequation.org/about.html

Terms of use
------------

Infomap is released under a dual licence.

To give everyone maximum freedom to make use of Infomap
and derivative works, we make the code open source under
the GNU General Public License version 3 or any
later version (see `LICENSE_GPLv3.txt`_).

For a non-copyleft license, please contact us.

.. _LICENSE_GPLv3.txt: https://github.com/mapequation/infomap/blob/master/LICENSE_GPLv3.txt
