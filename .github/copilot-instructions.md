
# Code review guidelines

When performing a code review, keep in mind the following.

## Coding standards

 - Check that all values of integer type are represented as `igraph_int_t`, except in justified special cases. Such special cases are:
    * The interruption counter used with `IGRAPH_ALLOW_INTERRUPTION_LIMITED`, typically called `iter`, which should be an `int`.
    * When interfacing with a library that uses a different type. Before passing igraph data to this library, it must be checked that the values will not overflow.

 - The `long` and `long long` types must never be used in native igraph code, as their sizes are inconsistent between platforms.

 - Check that the parameter ordering of newly added public functions conforms to the guidelines in our wiki at https://github.com/igraph/igraph/wiki/Guidelines-for-function-argument-ordering

 - Check that input parameters that are passed by reference (i.e. as a pointer) are marked as `const`. When a parameter that is passed by reference to a public function is modified, and therefore cannot be `const`, the documentation of that public function should make this clear. Parameters that are passed by value should not be marked as `const`.

 - Check that header includes are ordered according to the following guidelines:
    * The first included header must be the one that declares the functions *defined* in this file (optionally followed by an empty line).
    * Next come public igraph headers (those in `/include` and having names starting with `igraph_`), ordered alphabetically.
    * Next come private igraph headers (those in subdirectories of `/src`), ordered alphabetically.
    * Next come headers of igraph's dependencies, if any are used.
    * Standard C library headers must come last.

 - Check that all newly added public and private functions follow these conventions:
    * Public functions must have a name starting with `igraph_`.
    * Private functions that can be used from multiple translation units must have a name starting with `igraph_i_`, and must be included in a private header.
    * Private functions that are local to their translation unit must be marked as `static` and should not use the above two prefixes in their name.

 - Public igraph functions that are written in C++ (i.e. defined in `.cpp` files) must catch all exceptions before returning. This can be done using the `IGRAPH_HANDLE_EXCEPTIONS_BEGIN;` and `IGRAPH_HANDLE_EXCEPTIONS_END;` macros.

 - Check that `for` loops that iterate through vectors do not repeatedly call the `igraph_vector_size()` (or equivalent) function, unless the vector size changes during the iteration. Instead, the vector size should be saved into a variable before entering the `for` loop.

 - Check that the formatting of the code is consistent with the rest of the igraph codebase. If the formatting is considerably different, suggest reformatting using `astyle` and the `.astylerc` file at the root of the repository.

 - Check that C code follows the C99 standard while staying compatible with recent MSVC versions. VLAs and array designators must not be used.

## Memory management and error handling

 - Check that in the library implementation, all calls to functions that return `igraph_error_t` are protected with `IGRAPH_CHECK`, or that there is a comment such as `/* reserved */` that makes it clear why such protection is not needed. This ensures that allocated resources will be freed if an error occurs. This is not necessary for tests, benchmarks and examples.

 - Check that whenever a new data structure is initialized, it is registered on igraph's "finally stack" using `IGRAPH_FINALLY()`, or by using a convenience macro that combines initialization with registration on the "finally stack" (e.g. `IGRAPH_VECTOR_INT_INIT_FINALLY()`). The `IGRAPH_FINALLY()` call must follow the initialization call directly.

 - Check that when data structures are destroyed, they are also removed from the "finally stack" using `IGRAPH_FINALLY_CLEAN()`. The parameter passed to `IGRAPH_FINALLY_CLEAN()` indicates the number of items to remove. Thus n calls to destructors (usually named as `igraph_..._destroy()`) must be followed by `IGRAPH_FINALLY_CLEAN(n)`.

 - Check that data structures are destroyed in a reverse order compared to their initialization. This ensures that the "finally stack" is kept in a consistent state.

 - Suggest replacing pairs of a data structure initialization and `IGRAPH_FINALLY` calls with a single convenience macro, when such a macro is available. For example, `IGRAPH_CHECK(igraph_vector_init(&v, 0)); IGRAPH_FINALLY(igraph_vector_destroy, &v);` can be replaced by `IGRAPH_VECTOR_INIT_FINALLY(&v, 0);`.

 - Check that `calloc()` and `realloc()` are not used directly. Instead, the igraph-specific `IGRAPH_CALLOC()` and `IGRAPH_REALLOC()` macros should be used.

## Input validation

 - Check that newly added public functions validate their input.

 - Check that integer overflow checks are performed when there is a risk of overflow. `math/safe_intop.h` contains helper macros and functions for overflow-safe arithmetic.

 - Check that edge and vertex weights, as well as spatial coordinates, are protected against invalid values such as NaN, or that there is a comment explaining why such a check is not present.

## Testing

 - Check that all newly added public functions have tests. If they do not, include a link in your response to the testing guide in our wiki at https://github.com/igraph/igraph/wiki/How-to-write-unit-tests%3F. Check that edge cases, such as the null graph and singleton graph (i.e. graphs with 0 and 1 vertices) are included in tests when relevant.

 - Tests for public igraph functions must only include the main igraph header `<igraph.h>` (with angle brackets), not any of igraph's sub-headers such as `"igraph_interface.h"`.

 - Tests must contain `#include "test_utilities.h"` at the top. Whenever possible, tests should use the helper functions in this header (such as `print_vector()`) for printing output. Tests must use the `VERIFY_FINALLY_STACK();` macro to check the consistency of igraph's "finally stack". This macro is typically called before the `main()` function returns.

 - Ensure that important edge cases are tested. What these are differs from situation to situation, but typical cases are making sure that graphs with zero vertices (null graph) or zero edges (empty graph) are handled well.

 - It is not necessary to use `IGRAPH_CHECK`, `IGRAPH_FINALLY` and related functions in tests. Do not suggest using these.

## Documentation and error messages

 - Check that newly added public functions are documented in detail. Check that the documentation is clear, concise, complete and accurate. When the documentation of public functions is missing or incomplete, point this out.

 - Check the spelling and grammar of documentation and error/warning messages.

 - Check that American spelling is used in all documentation, error messages and symbol names. This instruction does not apply to code comments.

 - Check that the documentation of newly added or updated public functions describes all function parameters, as well as the return value. The parameters must be documented in the same order, and must have the same name as in the function signature.

 - Check that the documentation of newly added public functions has a note on their time complexity.

 - Check that within the main description of the function each paragraph, except the first, is preceded by `</para><para>` on a separate line.

 - Check that any scientific articles referenced in the documentation are accompanied by a weblink, ideally a DOI link of the form `https://www.doi.org/...`. References to books do not need a weblink.

 - Check that newly added public functions are marked as experimental by including `\experimental` in their documentation. `\experimental` must appear on a separate line, directly following the `\brief` line.

 - Check that the documentation of newly added public functions is referenced in the appropriate `.xxml` file in the `/doc` subdirectory of the repo.

 - Check that the documentation of public functions conforms to the guidelines in our wiki at https://github.com/igraph/igraph/wiki/Documentation-guidelines. If they do not, include a link in your response to these guidelines.

 - Check that error and warning messages use sentence case and end in a fullstop. Refer to https://github.com/igraph/igraph/wiki/Error-reporting-guidelines for more details, and link this page when commenting on error or warning messages.

 - Check that error messages do not include indices (such as vertex or edge IDs), as igraph has high-level interfaces to languages using both 0-based and 1-based indexing.

## Review hints

 - Keep in mind that igraph vectors and matrices are always initialized to null values. However, in contrast to initialization, *resizing* a vector to a larger size does not ensure that newly added values are null.

## Miscellaneous points to check during reviews

 - Check that newly added public functions that are intended for use from igraph's high-level interfaces have an interface description in `interfaces/functions.yaml`, and that the interface description uses the appropriate types. Note that the Stimulus CI check verifies that the interface description matched the C function prototype.

 - Check that no irrelevant changes have been made that could cause merge conflicts in the future. Examples of inappropriate changes include reformatting code that is not relevant to the pull request, or changing whitespace in unrelated code.

 - Check that newly added public functions follow the guidelines within this checklist in our wiki at https://github.com/igraph/igraph/wiki/Checklist-for-new-(and-old)-functions. If they do not, include a link in your response to this checklist.

 - Newly added files must contain an appropriate copyright header.

 - If `CHANGELOG.md` was updated, check that the description of changes is complete and accurate. If `CHANGELOG.md` was not updated, do not suggest doing so.

 - When performing a code review on a pull request that adds new features, you may comment on adherence to the general design principles of igraph explained at https://github.com/igraph/igraph/wiki/Design-principles-for-functions
