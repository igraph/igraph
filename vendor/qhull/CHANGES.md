Modifying from https://github.com/qhull/qhull on commit d1c2fc0caa5f644f3a0f220290d4a868c68ed4f6

Changed `userprintf_r.c` to not print when passed a null pointer as file handle.

Changes the following settings in `user_r.h`:

 - `#define qh_KEEPstatistics 0`, do not incude statistics gathering code
 - `#define qh_NOtrace`, do not include tracing code

Patched with https://github.com/qhull/qhull/pull/165.
