# Versioning and stability

This document is provided for informational purposes only, to help igraph users understand how igraph's programming interface is evolving, and what stability guaratees you can rely on.  It concerns the igraph C library only.  igraph's high level interfaces (R, Python, Mathematica) have their own separate versioning schemes and compatibility policies.

## Versioning

Starting with version 1.0, igraph follows the spirit of semantic versioning, with some differences, as described below.  The version number consists of three parts in the `MAJOR.MINOR.PATCH` format.

- `MAJOR` is incremented after making _incompatible changes_ to the stable programming interface.  Major releases are intended to be infrequent, and are accompanied by release notes where we make the effort to provide detailed guidance on adapting to incompatible changes.
- `MINOR` is incremented after making _additions_ to the stable programming interface.  Minor releases are issued regularly.
- `PATCH` is incremented when making changes that do not affect compatibility (usually bug fixes, documentation improvements, or build systems changes).

The three version parts are available at compile time as the macros `IGRAPH_VERSION_MAJOR`, `IGRAPH_VERSION_MINOR` and `IGRAPH_VERSION_PATCH`, or at runtime through the `igraph_version()` function.

The majority of public, documented functions are considered to be part of the _stable programming interface_, but there are some notable exceptions:

- **Experimental functions** may change at any time without notice.  These are clearly marked in their documentation ([example](https://igraph.org/c/html/0.10.13/igraph-Generators.html#igraph_chung_lu_game)).  They are also marked in igraph's header files with the `IGRAPH_EXPERIMENTAL` macro ([example](https://github.com/igraph/igraph/blob/3629c46b2784cb10fc27fc6e9fab4404a13d031c/include/igraph_cycles.h#L49-L53)).  Most newly added functions start out as _experimental_, and stay in this state until we are confident in their design, typically for one or two minor releases.  User feedback about experimental functions is particularly welcome.  We make the effort to avoid changes to experimental functions in patch releases, but do not guarantee this.
- **Internal and undocumented functions** are not part of the stable programming interface, not even if they are present in public headers.  They may change at any time. The names of internal functions usually start with the prefix `igraph_i_` (capitalized for macros), while public functions start with `igraph_`.

## Symbol lifecycle

Most new symbols start out as _experimental_ in minor releases.  Eventually, their API is declared stable, and the experimental marker is removed from their declaration and documentation in an upcoming minor release.

Symbols go through a deprecation phase before they are removed.  Deprecated symbols are marked in their documentation ([example](https://igraph.org/c/html/0.10.13/igraph-Structural.html#igraph_clusters)), and functions are prefixed with `IGRAPH_DEPRECATED` in the public headers ([example](https://github.com/igraph/igraph/blob/997f59ad742892fff199824a248fab382b40f526/include/igraph_components.h#L45-L47)). With GCC-compatible compilers, use the `-Wdeprecated` flag to get a warning for the use of deprecated functions, but keep in mind that deprecation warnings are not supported for all symbol types (e.g. macros) with all compilers.  The ultimate reference for deprecations is the [changelog][1].

## Stability of behaviour

Whether changes in function behaviour are considered breaking is somewhat subjective, and is decided on a case-by-case basis.  Expect some changes within minor releases.  For example, a function that ignored edge multiplicities may gain support for multigraphs in a new minor release.  Do not rely on details of behaviour that are not explicitly documented.

A notable case is stochastic functions: we do not guarantee that the same output is returned across different releases (even patch releases) for the same random seed.  We only guarantee the same statistical properties.

## Advice to users and package maintainers

**Users:**

For as long as you don't use _experimental_ functions, you can be confident that your code will continue to work with future releases having the same major version.  If you do use _experimental_ functions, it is your responsibility to check the [release notes][1] of each igraph release and adapt accordingly.  The use of _internal_ functions is completely unsupported: if you feel you need them, please talk to us.

If you do use _experimental_ functions, make this clear in your README file for the benefit of package maintainers.

While igraph comes with multiple header files, only `#include <igraph.h>` is supported. The rest of the headers exist for internal organizational purposes only, and may change without notice.

**Package maintainers:**

Software that does not use experimental functions from igraph can safely link to future igraph versions with the same major version.  Ask the developer of any software you are packaging if they are using experimental igraph functions.

The high-level interfaces of igraph do use both experimental and internal functions. Each high-level interface release is only guaranteed to be compatible with one specific release of C/igraph. As of this writing, this is a concern only for the Python interface, as the other interfaces (R and Mathematica) cannot link dynamically to C/igraph.

We provide the `IGRAPH_WARN_EXPERIMENTAL` compile-time macro to help maintainers in determining whether a piece of software uses experimental igraph functions or not. Compilers that support `__attribute__((__warning(...)))` clauses will issue a warning when `IGRAPH_WARN_EXPERIMENTAL` is defined to a non-zero value at compile-time and an experimental function is used somewhere in the code.

## Notes

For the purposes of this document, "API compatibility" means that the same sources can be compiled using headers from different igraph versions.  "ABI compatibility" means that a program that only uses stable API can be linked to a different version of the igraph shared library than what it was compiled with.

We strive to maintain both API and ABI compatibility.

However, it must be pointed out that we do not support manipulating the same in-memory igraph data structures with different igraph versions (for example, if two libraries that exchange data are each statically linked to different igraph versions).

 [1]: https://github.com/igraph/igraph/blob/main/CHANGELOG.md
