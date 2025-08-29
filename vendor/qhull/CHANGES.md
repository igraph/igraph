Modifying from https://github.com/qhull/qhull on commit c814cb78883e6bb63a3b09a8a981df61d3ecaabb

Changed `userprintf_r.c` to not print when passed a null pointer as file handle.

Changes the following settings in `user_r.h`:

 - `#define qh_KEEPstatistics 0`, do not incude statistics gathering code
 - `#define qh_NOtrace`, do not include tracing code

Patched with https://github.com/qhull/qhull/pull/165 and https://github.com/qhull/qhull/pull/167.
