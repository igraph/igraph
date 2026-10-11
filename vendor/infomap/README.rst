|ci|

Infomap
=======

Infomap is a network clustering algorithm based on the `Map equation`_.
This repository contains the native CLI, the Python package, the R package,
the JavaScript browser worker and Node.js package, the Docker images, the
tutorial notebooks, and the source for the published Python documentation.

Start with `mapequation.org/infomap/`_ for the user guide, the
`Infomap Python API`_ for Python examples and tutorial notebooks, and
`CHANGELOG.md`_ for release notes.

For contributing, security reporting, and maintainer workflows, see
`CONTRIBUTING.md`_, `SECURITY.md`_, `BUILD.md`_, `ARCHITECTURE.md`_, and
`AGENTS.md`_.

.. |ci| image:: https://github.com/mapequation/infomap/actions/workflows/ci.yml/badge.svg
   :target: https://github.com/mapequation/infomap/actions/workflows/ci.yml
   :alt: CI

.. _Map equation: https://www.mapequation.org/publications.html?utm_source=infomap&utm_medium=readme&utm_campaign=infomap#Rosvall-Axelsson-Bergstrom-2009-Map-equation
.. _`mapequation.org/infomap/`: https://www.mapequation.org/infomap/?utm_source=infomap&utm_medium=readme&utm_campaign=infomap
.. _`CHANGELOG.md`: https://github.com/mapequation/infomap/blob/master/CHANGELOG.md
.. _`CONTRIBUTING.md`: https://github.com/mapequation/infomap/blob/master/CONTRIBUTING.md
.. _`SECURITY.md`: https://github.com/mapequation/infomap/blob/master/SECURITY.md
.. _`BUILD.md`: https://github.com/mapequation/infomap/blob/master/BUILD.md
.. _`ARCHITECTURE.md`: https://github.com/mapequation/infomap/blob/master/ARCHITECTURE.md
.. _`AGENTS.md`: https://github.com/mapequation/infomap/blob/master/AGENTS.md

Install
-------

Python package
^^^^^^^^^^^^^^

Install from `PyPI`_:

.. code-block:: bash

    pip install infomap

Install optional integrations for common Python graph and analysis workflows:

.. code-block:: bash

    pip install "infomap[networkx]"
    pip install "infomap[igraph]"
    pip install "infomap[pandas]"

Upgrades use the usual ``pip`` flow:

.. code-block:: bash

    pip install --upgrade infomap

The package also installs the ``infomap`` CLI entry point.
The Python API reference lives at `Infomap Python API`_.

Quick start with Python:

.. code-block:: python

    import networkx as nx
    import infomap

    graph = nx.karate_club_graph()
    result = infomap.run(graph, seed=123, num_trials=20)

    print(result.num_top_modules, result.codelength)
    print(result.modules())  # {node_id: module_id}

``infomap.run`` accepts a NetworkX or igraph graph, a SciPy sparse matrix, a
``(2, E)`` edge index, a network file path, or an iterable of links, and returns
an immutable ``Result``. If you only need the communities in the graph's own
node labels, use ``infomap.find_communities(graph, seed=123, num_trials=20)``,
which returns a NetworkX-style ``list`` of ``set``\ s of node labels (its igraph
counterpart ``infomap.find_igraph_communities`` returns an
``igraph.VertexClustering``).

For incremental construction -- adding nodes and links one at a time -- build a
``Network`` and run it, reading results off the returned ``Result``:

.. code-block:: python

    from infomap import Network, run

    net = Network()
    net.add_link(0, 1)
    net.add_link(1, 2)
    result = run(net, two_level=True, num_trials=20, seed=123)

    print(result.num_top_modules, result.codelength)
    print(result.to_dataframe(columns=["node_id", "module_id", "flow"], index="node_id"))

To keep one configured engine and run it repeatedly -- or to maintain code
written against the original API -- the stateful ``Infomap`` class works the same
way (``im = Infomap(...); im.add_link(...); result = im.run()``).

For Jupyter, start with the
`quickstart notebook <https://github.com/mapequation/infomap/blob/master/examples/notebooks/quickstart.ipynb>`_.
It shows the Infomap result summary, dataframe inspection, a copyable static
network partition helper, and export paths for further analysis.

.. _PyPI: https://pypi.org/project/infomap/
.. _`Infomap Python API`: https://mapequation.org/infomap-python-docs/

R package
^^^^^^^^^

Pre-built binaries are published on `r-universe`_; this is the recommended path:

.. code-block:: r

    install.packages(
      "infomap",
      repos = c("https://mapequation.r-universe.dev", "https://cloud.r-project.org")
    )

Quick start with R:

.. code-block:: r

    library(infomap)

    im <- Infomap(silent = TRUE, two_level = TRUE, num_trials = 20)
    im$add_link(0, 1)
    im$add_link(1, 2)
    im$run()

    print(im$num_top_modules)
    print(im$codelength)

See ``?Infomap`` for the user-facing constructor plus the ``InfomapClass``
method and active-binding reference. The R-specific source README lives at
`interfaces/R/infomap/README.md`_.

.. _r-universe: https://mapequation.r-universe.dev
.. _`interfaces/R/infomap/README.md`: https://github.com/mapequation/infomap/blob/master/interfaces/R/infomap/README.md

Homebrew CLI
^^^^^^^^^^^^

If you want the native CLI without the Python package, install the tap and
formula with:

.. code-block:: bash

    brew tap mapequation/infomap
    brew install infomap

Or install directly in one command:

.. code-block:: bash

    brew install mapequation/infomap/infomap

Upgrade the CLI with the normal Homebrew flow:

.. code-block:: bash

    brew upgrade infomap

The Homebrew formula installs Bash and Zsh completion files into Homebrew's
standard completion directories.

JavaScript package
^^^^^^^^^^^^^^^^^^

The package is published on `NPM`_ and provides a browser web worker, a Node.js
module, and an ``infomap`` command line tool:

.. code-block:: bash

    npm install @mapequation/infomap

Quick start in Node.js with the ``@mapequation/infomap/node`` entry point:

.. code-block:: javascript

    import { run } from "@mapequation/infomap/node";

    const network = "0 1\n0 2\n0 3\n1 2\n3 4\n3 5\n4 5";
    const result = await run(network, { args: ["-o", "tree,json", "-2"] });

    console.log(result.json.codelength);
    console.log(result.tree);

Installing the package also provides an ``infomap`` command that behaves like
the native binary:

.. code-block:: bash

    npx @mapequation/infomap network.net . --tree

Browser worker and React usage are documented on the `NPM`_ package page.

.. _NPM: https://www.npmjs.com/package/@mapequation/infomap

Docker
^^^^^^

Multi-arch images are published to `GHCR`_ for ``linux/amd64`` and
``linux/arm64``:

- ``ghcr.io/mapequation/infomap:latest``
- ``ghcr.io/mapequation/infomap:X.Y.Z``
- ``ghcr.io/mapequation/infomap:notebook``
- ``ghcr.io/mapequation/infomap:notebook-X.Y.Z``

Run the CLI image with:

.. code-block:: bash

    docker run -it --rm \
        -v "$(pwd)":/data \
        ghcr.io/mapequation/infomap:latest \
        [infomap arguments]

Start the notebook image with:

.. code-block:: bash

    docker run --rm \
        -p 8888:8888 \
        ghcr.io/mapequation/infomap:notebook \
        start.sh jupyter lab

The notebook image includes the survey companion notebooks from
``examples/notebooks`` and opens in that workspace by default. To keep local
copies or outputs, mount a host directory as a separate workspace path:

.. code-block:: bash

    docker run --rm \
        -v "$(pwd)":/home/jovyan/work/local \
        -p 8888:8888 \
        ghcr.io/mapequation/infomap:notebook \
        start.sh jupyter lab

The Dockerfiles in this repository are also smoke-tested in CI and can be
built locally:

.. code-block:: bash

    docker build -f docker/infomap.Dockerfile -t infomap:local .
    docker build -f docker/notebook.Dockerfile -t infomap:notebook-local .

Or use the local Compose file:

.. code-block:: bash

    docker compose run --rm infomap

.. _GHCR: https://github.com/mapequation/infomap/pkgs/container/infomap

Build from source
-----------------

Building locally requires a working ``gcc`` or ``clang`` toolchain.

.. code-block:: bash

    git clone git@github.com:mapequation/infomap.git
    cd infomap
    make build-native

On macOS, the default OpenMP-enabled build may require Homebrew ``libomp``.
If OpenMP is unavailable, use:

.. code-block:: bash

    make build-native OPENMP=0

This creates the ``Infomap`` binary in the repository root.
Show the available CLI options with:

.. code-block:: bash

    ./Infomap --help

Install shell completion scripts manually with:

.. code-block:: bash

    mkdir -p ~/.zfunc
    ./Infomap --completion zsh > ~/.zfunc/_Infomap

    mkdir -p ~/.local/share/bash-completion/completions
    ./Infomap --completion bash > ~/.local/share/bash-completion/completions/infomap

For Zsh, make sure ``~/.zfunc`` is in ``fpath`` and ``compinit`` is loaded from
``~/.zshrc``. For Bash, make sure ``bash-completion`` is sourced from
``~/.bashrc``.

See ``BUILD.md`` for platform-specific maintainer build details.

Maintainers should use:

- ``BUILD.md`` for local build and verification commands
- ``RELEASING.md`` for the release flow
- ``ARCHITECTURE.md`` for ownership and source-of-truth rules
- ``AGENTS.md`` for repo-local maintenance guidance
- ``CONTRIBUTING.md`` for pull request and contributor guidance
- ``SECURITY.md`` for vulnerability reporting

Agent skill
-----------

This repository includes an Infomap agent skill in ``skills/infomap/`` for
reproducible CLI, Python, R, and notebook research workflows.

Feedback
--------

Usage questions and setup help belong in `GitHub Discussions`_.
Bug reports and feature requests belong in `GitHub issues`_.

.. _`GitHub Discussions`: https://github.com/mapequation/infomap/discussions
.. _`GitHub issues`: https://github.com/mapequation/infomap/issues

Authors
-------

Daniel Edler, Anton Holmgren, Martin Rosvall

For contact information, see `mapequation.org/about.html`_.

.. _`mapequation.org/about.html`: https://www.mapequation.org/about.html

Terms of use
------------

Infomap is released under a dual license.

The code is available under the GNU General Public License version 3 or any
later version; see `LICENSE_GPLv3.txt`_.
For a non-copyleft license, please contact us.

.. _`LICENSE_GPLv3.txt`: https://github.com/mapequation/infomap/blob/master/LICENSE_GPLv3.txt
