This file lists changes that were made to the original Bliss package (version 0.75) to integrate it into igraph.

Exclude `CMakeLists.txt`, `Doxyfile`, `Makefile-manual`, `readme.txt`. Make sure not to accidentally overwrite igraph's own `bliss/CMakeLists.txt`.

Removed `bliss.cc`, `bliss_C.cc`, `bliss_C.h`.

Remove `timer.hh`. Remove references to `timer.hh` and `Timer` class in `graph.cc`.

Replace `#pragma once` by traditional header guards in all headers.

### In `bignum.hh`:

Replace `#include <gmp.h>` by `#include "internal/gmp_internal.h"`.

At the beginning, add `#define BLISS_USE_GMP`. Verify that this macro is only used in this file.

### In `defs.cc` and `defs.hh`:

Remove the `...` argument from `fatal_error` for simplicity, and make the function simply invoke `IGRAPH_FATAL`.

### In `graph.cc`:

Define `_INTERNAL_ERROR` in terms of `IGRAPH_FATAL`.

### MSVC compatibility

Bliss uses `and`, `or`, etc. instead of `&&`, `||`, etc. These are not supported by MSVC by default. Bliss 0.74 uses the `/permissive` option to enable support in MSVC, but this option is only supported wit VS2019. Instead, in igraph we add the following where relevant:

```
/* Allow using 'and' instead of '&&' with MSVC */
#if _MSC_VER
#include <ciso646>
#endif
```
