Modifying from https://github.com/qhull/qhull on commit d1c2fc0caa5f644f3a0f220290d4a868c68ed4f6

Changed `userprintf_r.c` to not print when passed a null pointer as file handle.

Changed `usermem_r.c` to not reference `exit()`. Even though this `exit()` is never called, it is necessary to keep igraph free of any references to it to satisfy CRAN requirements and comply with Debian guidelines.

Changed the following settings in `user_r.h`:

 - `#define qh_KEEPstatistics 0`, do not incude statistics gathering code
 - `#define qh_NOtrace`, do not include tracing code

Patched with https://github.com/qhull/qhull/pull/165.
