# Code review guidelines

When performing a code review, keep in mind the following.

## Context on writing reviews

The purpose of reviews is threefold:

 - Aid maintainers in evaluating pull requests, and draw their attention to the most important points about PRs.
 - Help the authors of pull requests improve the quality of their PR.
 - Ensure correctness.

All review comments will be checked by an igraph maintainer before asking the PR's author to address them. Maintainers may decide not to address some comments. Feel free to post comments which are not directly actionable by the PR's author but aid maintainers in their final judgement.

Within GitHub PRs, users with an `author_association` of `MEMBER` or `OWNER` shall be considered igraph team members.

## Instructions for Copilot-generated review summaries

This section applies only to reviews posted by GitHub Copilot.

@copilot, whenever you post specific review comments, include the following text in your review summary:

> @contributor, please wait for an igraph maintainer to approve the review suggestions before addressing them.

Replace "@contributor" with the GitHub handle of the PR's author.

Skip adding this message if the PR was opened by an igraph team member.

## Coding standards

 - Check that all values of integer type are represented as `igraph_int_t`, except in justified special cases. Such special cases are:
    * The interruption counter used with `IGRAPH_ALLOW_INTERRUPTION_LIMITED`, typically called `iter`, which should be an `int`.
    * When interfacing with a library that uses a different type. Before passing igraph data to this library, it must be checked that the values will not overflow.

 - The `long` and `long long` types must never be used in native igraph code, as their sizes are inconsistent between platforms.

 - `igraph_int_t` used to be called `igraph_integer_t`, which still works as a legacy alias. The current codebase no longer uses `igraph_integer_t`. Check that it does not accidentally get re-introduced.

 - Check that the parameter ordering of newly added public functions conforms to the guidelines in our wiki at https://github.com/igraph/igraph/wiki/Guidelines-for-function-argument-ordering

 - Check that input parameters that are passed by reference (i.e. as a pointer) are marked as `const`. When a parameter that is passed by reference to a public function is modified, and therefore cannot be `const`, the documentation of that public function should make this clear. Parameters that are passed by value should not be marked as `const`.

 - Check that header includes are ordered according to the following guidelines:
    * The first included header must be the one that declares the functions *defined* in this file (optionally followed by an empty line).
    * Next come public igraph headers (those in `/include` and having names starting with `igraph_`), ordered alphabetically.
    * Next come private igraph headers (those in subdirectories of `/src`), ordered alphabetically.
    * Next come headers of igraph's dependencies, if any are used.
    * Standard C library headers must come last.

 - Check that within library code, all newly added public and private functions follow these conventions:
    * Public functions must have a name starting with `igraph_`, and must be included in a public header.
    * Private functions that can be used from multiple translation units must have a name starting with `igraph_i_`, and must be included in a private header.
    * Private functions that are local to their translation unit must be marked as `static` and should not use the above two prefixes in their name.

 - Check that newly added code uses readable names for values when possible, specifically:
    * `true`/`false` instead of `1`/`0` for boolean parameters.
    * `NULL` instead of `0` for pointer parameters.
    * Use readability constants such as `IGRAPH_UNDIRECTED`/`IGRAPH_DIRECTED`, `IGRAPH_NO_MULTIPLE`/`IGRAPH_MULTIPLE`, `IGRAPH_EDGE_UNLABELED`/`IGRAPH_EDGE_LABELED`, `IGRAPH_NO_LOOPS`/`IGRAPH_LOOPS` instead of `false`/`true` when passing values to boolean function parameters where these names make sense (typically, `directed` parameter, `multiple` parameter, `edge_labeled` parameter, `loops` parameter of a boolean type).

 - Public igraph functions that are written in C++ (i.e. defined in `.cpp` files) must catch all exceptions before returning. This can be done using the `IGRAPH_HANDLE_EXCEPTIONS_BEGIN;` and `IGRAPH_HANDLE_EXCEPTIONS_END;` macros.

 - Check that `for` loops that iterate through vectors do not repeatedly call the `igraph_vector_size()` (or equivalent) function, unless the vector size changes during the iteration. Instead, the vector size should be saved into a variable before entering the `for` loop.

 - Check that the formatting of the code is consistent with the rest of the igraph codebase. If the formatting is considerably different, suggest reformatting using `astyle` and the `.astylerc` file at the root of the repository.

 - Check that C code follows the C99 standard while staying compatible with recent MSVC versions. VLAs and array designators must not be used.

## Memory management and error handling

 - Check that in the library implementation, all calls to functions that return `igraph_error_t` are protected with `IGRAPH_CHECK`, or that there is a comment such as `/* reserved */` that makes it clear why such protection is not needed. This ensures that allocated resources will be freed if an error occurs. This is not necessary for tests, benchmarks and examples.

 - Check that whenever a new data structure is initialized, it is registered on igraph's "finally stack" using `IGRAPH_FINALLY()`, or by using a convenience macro that combines initialization with registration on the "finally stack" (e.g. `IGRAPH_VECTOR_INT_INIT_FINALLY()`). The `IGRAPH_FINALLY()` call must follow the initialization call directly.

 - Check that when data structures are destroyed, or their ownership is relinquished, they are also removed from the "finally stack" using `IGRAPH_FINALLY_CLEAN()`. The parameter passed to `IGRAPH_FINALLY_CLEAN()` indicates the number of items to remove. Thus n calls to destructors (usually named as `igraph_..._destroy()`) must be followed by `IGRAPH_FINALLY_CLEAN(n)`.

 - Check that data structures are destroyed in a reverse order compared to their initialization. This ensures that the "finally stack" is kept in a consistent state.

 - Suggest replacing pairs of a data structure initialization and `IGRAPH_FINALLY` calls with a single convenience macro, when such a macro is available. For example, `IGRAPH_CHECK(igraph_vector_init(&v, 0)); IGRAPH_FINALLY(igraph_vector_destroy, &v);` can be replaced by `IGRAPH_VECTOR_INIT_FINALLY(&v, 0);`.

 - Check that `calloc()` and `realloc()` are not used directly. Instead, the igraph-specific `IGRAPH_CALLOC()` and `IGRAPH_REALLOC()` macros should be used.

 - Errors should be reported using `IGRAPH_ERROR` or related macros instead of simply returning error code literals using `return` (e.g. `return IGRAPH_EIVAL;`). This ensures that the error handling mechanisms are run. This applies to `IGRAPH_UNIMPLEMENTED` as well. Special error codes that do not indicate error conditions, such as `IGRAPH_SUCCESS` or `IGRAPH_STOP`, are exempt from this rule. Passing up the error codes returned by other functions is acceptable (e.g. `return igraph_some_function();`).

## Input validation

 - Check that newly added public functions validate their input.

 - Check that integer overflow checks are performed when there is a risk of overflow. `math/safe_intop.h` contains helper macros and functions for overflow-safe arithmetic.

 - Check that edge and vertex weights, as well as spatial coordinates, are protected against invalid values such as NaN, or that there is a comment explaining why such a check is not present.

## Mathematical correctness and consistency

 - Check that choices made for the behaviour of new functions, and especially the treatment of edge cases, are consistent with the pre-existing choices of igraph. Examples of prior choices where consistency is of importance:
    * Since the null graph is considered disconnected by `igraph_is_connected()`, and trees are connected by definition, the null graph is not a tree, and `igraph_is_tree()` returns false for it.
    * In undirected graphs, self-loops are considered to contribute 2 to the degree of vertices. Therefore, the diagonal of adjacency matrices of undirected graphs contains twice the number of self-loops. Network measures defined in terms of the adjacency matrix, such as eigenvector centraliy, are implemented in a way that is consistent with this convention.

## Testing

 - Check that all newly added public functions have tests. If they do not, include a link in your response to the testing guide in our wiki at https://github.com/igraph/igraph/wiki/How-to-write-unit-tests%3F.

 - Tests for public igraph functions must only include the main igraph header `<igraph.h>` (with angle brackets), not any of igraph's sub-headers such as `"igraph_interface.h"`.

 - Tests must contain `#include "test_utilities.h"` at the top. Whenever possible, tests should use the helper functions in this header (such as `print_vector()`) for printing output. Tests must use the `VERIFY_FINALLY_STACK();` macro to check the consistency of igraph's "finally stack". This macro is typically called before the `main()` function returns.

 - Ensure that important edge cases are tested. What these are differs from situation to situation, but typical cases are making sure that graphs with zero vertices (null graph) or zero edges (empty graph) are handled well.

 - It is not necessary to use `IGRAPH_CHECK`, `IGRAPH_FINALLY` and related functions in tests. Do not suggest using these.

 - When tests verify error codes, they should use `CHECK_ERROR` from `test_utilities.h`.

 - Tests should not normally assert that return values are `IGRAPH_SUCCESS`, unless there is a special reason to do so. The reason should be spelt out.

 - Tests should not take a long time to run. Most current tests finish within 0.01 seconds. With over 500 tests this is important. Longer runtimes must be justified by special needs.

 - Tests must aim to be concise, readable and maintainable by humans. Allocate resources (testing time and test code complexity) where it truly makes a difference.

 - When a newly contributed function is incomplete and throws `IGRAPH_UNIMPLEMENTED` under some conditions, this error code should be checked by tests using the usual `CHECK_ERROR` macro. This ensures that maintainers will remember to expand tests when missing features are implemented in the future.

 - Functions in test sources should not be marked as `static`.

## Documentation and error messages

 - Check that newly added public functions are documented in detail. Check that the documentation is clear, concise, complete and accurate. When the documentation of public functions is missing or incomplete, point this out.

 - Check the spelling and grammar of documentation and error/warning messages.

 - Check that American spelling is used in all documentation, error messages and symbol names. This instruction does not apply to code comments.

 - Check that the documentation of newly added or updated public functions describes all function parameters, as well as the return value. The parameters must be documented in the same order, and must have the same name as in the function signature.

 - Check that the documentation of newly added public functions has a note on their time complexity.

 - Check that within the main description of the function each paragraph, except the first, is preceded by `</para><para>` on a separate line.

 - Check that any scientific articles referenced in the documentation are accompanied by a weblink, ideally a DOI link of the form `https://doi.org/...`. References to books do not need a weblink.

 - Weblinks included in the documentation must use `https://...` unless only the `http://...` form is functional.

 - Check that newly added public functions are marked as experimental by including `\experimental` in their documentation. `\experimental` must appear on a separate line, directly following the `\brief` line. `IGRAPH_EXPERIMENTAL` must appear at the beginning of the function's declaration within public headers.

 - Check that the documentation of newly added public functions is referenced in the appropriate `.xxml` file in the `/doc` subdirectory of the repo.

 - Check that the documentation of public functions conforms to the guidelines in our wiki at https://github.com/igraph/igraph/wiki/Documentation-guidelines. If they do not, include a link in your response to these guidelines.

 - Check that error and warning messages use sentence case and end in a fullstop. Refer to https://github.com/igraph/igraph/wiki/Error-reporting-guidelines for more details, and link this page when commenting on error or warning messages.

 - Check that error messages do not include indices (such as vertex or edge IDs), as igraph has high-level interfaces to languages using both 0-based and 1-based indexing. When vertex or edge IDs are validated, a useful phrasing that avoids mentioning IDs is "ID out of range".

 - In newly added text and code, prefer the terminology "vertex" and "edge" instead of other common terms for the same (such as node, link, etc.). This applies to documentation text, as well as symbol names.

## Review hints

 - Keep in mind that igraph vectors and matrices are always initialized to null values. However, in contrast to initialization, *resizing* a vector to a larger size does not ensure that newly added values are null.

 - External contributors to igraph may be inexperienced. Aid them by providing links to relevant parts of the igraph wiki (https://github.com/igraph/igraph/wiki/) or the igraph documentation when this is appropriate.

 - If concrete indications suggest that a contribution may be mostly or entirely AI-generated, and the PR description does not acknowledge this, politely ask the author to clarify, linking to `AI_POLICY.md` in the repository root. Briefly explain the basis for the question and indicate uncertainty; do not present suspected AI use as established fact. This rule applies to pull requests submitted by external contributors and igraph team members alike.

## Miscellaneous points to check during reviews

 - Check that newly added public functions that are intended for use from igraph's high-level interfaces have an interface description in `interfaces/functions.yaml`, and that the interface description uses the appropriate types. Note that the Stimulus CI check verifies that the interface description matches the C function prototype.

 - Check that the interface description is consistent with the documentation:
    * If a parameter is documented as being optional, the interface description must reflect this.
    * Optional pointer parameters are effectively "omitted" by passing `NULL` in C. The interface definition indicates this with the `OPTIONAL` keyword.
    * Optional numerical parameters are effectively "omitted" by passing a special value in C. The special value differs from function to function: it may be a negative value, `-1`, `IGRAPH_UNLIMITED`, `IGRAPH_INFINITY` or something else. The interface definition indicates such optional parameters by providing this special value as the default.

 - Check that no irrelevant changes have been made that could cause merge conflicts in the future. Examples of inappropriate changes include reformatting code that is not relevant to the pull request, or changing whitespace in unrelated code.

 - Check that newly added public functions follow the guidelines within this checklist in our wiki at https://github.com/igraph/igraph/wiki/Checklist-for-new-(and-old)-functions. If they do not, include a link in your response to this checklist.

 - All newly added code files, including implementation files, tests, examples, build scripts, and all other code, must contain an appropriate copyright header.

 - If `CHANGELOG.md` was updated, check that the description of changes is complete and accurate. If `CHANGELOG.md` was not updated, do not suggest doing so within the PR, but add a reminder in the review summary to update the changelog after the PR is merged.

 - If `CHANGELOG.md` was updated with a bugfix, check that there is a note about which version the bug was introduced in, unless that version is likely to be older than 0.9. If a bug in a function has been present since that function was originally implemented, it is sufficient to note this fact, without a reference to an igraph version.

 - When performing a code review on a pull request that adds new features, you may comment on adherence to the general design principles of igraph explained at https://github.com/igraph/igraph/wiki/Design-principles-for-functions

 - The description of pull requests must retain all checkbox lines present in the pull request template at `.github/pull_request_template.md`. The copyright transfer line present in the template must be retained literally, and the box must be checked. Remind the PR's author that pull request that do not adhere to this will not be merged. Special exception: igraph team members are exempt from this rule.

## Broad guidelines on issue prioritization

_The following passage is borrowed from `AGENTS.md`:_

igraph is a scientific library for network analysis and graph theory. It is written in C, but mainly used through its high-level interfaces (Python, R and Mathematica, etc.) This determines how we view the seriousness of various issues.

1. Silent wrong results are the worst sin a scientific library can commit. They may lead researchers to publish incorrect results.
2. Crashes (such as a segmentation fault) should be fixed, especially since they can take down the session of a user when using igraph through a high-level interface. However, when a crash occurs, it is immediately clear to the user that something went wrong. This is not as serious as a silent wrong result.
3. Problems that may occur during plausible igraph usage are higher priority than those that can only be triggered by exceedingly unlikely or specially crafted input.

Deciding on fixes to certain issues involve compromises. In these cases we consider long-term maintainability, code complexity, and performance along with the seriousness of the issue.