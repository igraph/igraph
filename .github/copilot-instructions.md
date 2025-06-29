When performing a code review, check that all newly added public functions have tests. If they do not, include a link in your response to the testing guide in our wiki at https://github.com/igraph/igraph/wiki/How-to-write-unit-tests%3F.

When performing a code review, check that newly added public functions are documented in detail, that their documentation is included in the appropriate `.xxml` file in the `doc` subdirectory of the repo.

When performing a code review, check that the documentation of public functions conforms to the guidelines in our wiki at https://github.com/igraph/igraph/wiki/Documentation-guidelines. If they do not, include a link in your response to these guidelines.

When performing a code review, check that the parameter ordering of newly added public functions conforms to the guidelines in our wiki at https://github.com/igraph/igraph/wiki/Guidelines-for-function-argument-ordering

When performing a code review, check that newly added public functions that are intended for use from igraph's high-level interfaces have an interface description in `interfaces/functions.yaml`, and that the interface description uses the appropriate types.

When performing a code review, keep in mind that igraph vectors and matrices are always initialized to null values. However, in contrast to initialization, *resizing* a vector to a larger size does not ensure that newly added values are null.

When performing a code review, check that all values of integer type are represented as `igraph_integer_t`, except in justified special cases such as the interruption counter used with `IGRAPH_ALLOW_INTERRUPTION_LIMITED`, which may be an `int`. The `long` and `long long` types must never be used.

When performing a code review, make sure that all calls to functions that return `igraph_error_t` are protected with `IGRAPH_CHECK`, or that there is a comment such as `/* reserved */` that makes it clear why such protection is not needed.

When performing a code review, check that newly added public functions follow the guidelines within this checklist in our wiki at https://github.com/igraph/igraph/wiki/Checklist-for-new-(and-old)-functions. If they do not, include a link in your response to this checklist.

When performing a code review, check that the formatting of the code is consistent with the igraph codebase.

When performing a code review, check that no irrelevant changes have been made that could cause merge conflicts in the future (such as reformatting unrelated code or adding whitespace to unrelated parts of the code).
