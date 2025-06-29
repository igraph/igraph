
# Code review guidelines

When performing a code review, keep in mind the following.

## Coding standards

 - Check that all calls to functions that return `igraph_error_t` are protected with `IGRAPH_CHECK`, or that there is a comment such as `/* reserved */` that makes it clear why such protection is not needed.

 - Check that all values of integer type are represented as `igraph_integer_t`, except in justified special cases such as the interruption counter used with `IGRAPH_ALLOW_INTERRUPTION_LIMITED`, which may be an `int`. The `long` and `long long` types must never be used.

 - Check that the parameter ordering of newly added public functions conforms to the guidelines in our wiki at https://github.com/igraph/igraph/wiki/Guidelines-for-function-argument-ordering

 - Check that input parameters that are passed by reference (i.e. as a pointer) are marked as `const`. When a parameter that is passed by reference to a public function is modified, and therefore cannot be const, the documentation of that public function should make this clear.

 - Suggest replacing pairs of a data structure initialization and `IGRAPH_FINALLY` calls with a single convenience macro, when such a macro is available. For example, `IGRAPH_CHECK(igraph_vector_init(&v, 0)); IGRAPH_FINALLY(igraph_vector_destroy, &v);` can be replaced by `IGRAPH_VECTOR_INIT_FINALLY(&v, 0);`.

 - Check that header includes are ordered according to the following guidelines:
    * The first included header must be the one that declares the functions *defined* in this file (optionally followed by an empty line).
    * Next come public igraph headers (those in `/include` and having names starting with `igraph_`), ordered alphabetically.
    * Next come private igraph headers (those in subdirectories of `/src`), ordered alphabetically.
    * Next come headers of igraph's dependencies, if any are used.
    * Standard C library headers must come last.

 - Public igraph functions that are written in C++ (i.e. defined in `.cpp` files) must catch all exceptions before returning. This can be done using the `IGRAPH_HANDLE_EXCEPTIONS_BEGIN;` and `IGRAPH_HANDLE_EXCEPTIONS_END;` macros.

 - Check that the formatting of the code is consistent with the rest of the igraph codebase. If the formatting is considerably different, suggest reformatting using `astyle` and the `.astylerc` file at the root of the repository.

## Input validation

 - Check that newly added public functions validate their input.

 - Check that integer overflow checks are performed when there is a risk of overflow. `math/safe_intop.h` contains helper macros and functions to perform overflow checks.

## Testing

 - Check that all newly added public functions have tests. If they do not, include a link in your response to the testing guide in our wiki at https://github.com/igraph/igraph/wiki/How-to-write-unit-tests%3F. Check that edge cases, such as the null graph and singleton graph (i.e. graph with 0 and 1 vertices) are included in tests when relevant.

 - Check that all newly added public and private functions follow these conventions:
    * Public functions must have a name starting with `igraph_`.
    * Private functions that can be used from multiple translation units must have a name starting with `igraph_i_`, and must be included in a private header.
    * Private functions that are local to their translation unit must be marked as `static` and should not use the above two prefixes in their name.

 - Test for public igraph functions must only include the main igraph header `<igraph.h>`, not any of igraph's sub-headers such as `igraph_interface.h`.

 - Tests must contain `#include "test_utilities.h"` at the top. Whenever possible, tests should use the helper functions in this header (such as `print_vector()`) for printing output. Tests must use the `VERIFY_FINALLY_STACK();` macro to check the consistency of igraph's "finally stack". This macro is typically called before the `main()` function returns.

## Documentation and error messages

 - Check that newly added public functions are documented in detail, that their documentation is included in the appropriate `.xxml` file in the `doc` subdirectory of the repo.

 - Check that the documentation of public functions conforms to the guidelines in our wiki at https://github.com/igraph/igraph/wiki/Documentation-guidelines. If they do not, include a link in your response to these guidelines.

 - When performing a code review, check that error and warning messages use sentence case and end in a fullstop. Messages should not include indices (such as vertex or edge IDs), as igraph has high-level interfaces to languages using both 0-based and 1-based indexing. Refer to https://github.com/igraph/igraph/wiki/Error-reporting-guidelines for more details, and link this page when commenting on error or warning messages.

## Review hints

 - Keep in mind that igraph vectors and matrices are always initialized to null values. However, in contrast to initialization, *resizing* a vector to a larger size does not ensure that newly added values are null.

## Miscellaneous points to check during reviews

 - Check that newly added public functions that are intended for use from igraph's high-level interfaces have an interface description in `interfaces/functions.yaml`, and that the interface description uses the appropriate types. Note that the Stimulus CI check verifies that the interface description matched the C function prototype.

 - Check that no irrelevant changes have been made that could cause merge conflicts in the future. Examples of inappropriate changes include reformatting code that is not relevant to the pull request, or changing whitespace in unrelated code.

 - Check that newly added public functions follow the guidelines within this checklist in our wiki at https://github.com/igraph/igraph/wiki/Checklist-for-new-(and-old)-functions. If they do not, include a link in your response to this checklist.

 - Newly added files must contain an appropriate copyright header.

 - When performing a code review on a pull request that adds new features, you may comment on adherence to the general design principles of igraph explained at https://github.com/igraph/igraph/wiki/Design-principles-for-functions

