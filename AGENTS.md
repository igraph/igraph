# Instructions for AI agents working on igraph

## Information on C/igraph development

Information on C/igraph development, design guidelines, notes on the use of development and debugging tools, is found in the igraph wiki. See the dedicated section below.

Check `.github/copilot-instructions.md` for how contributions are reviewed. When developing new code, treat the advice there as guidelines, and exercise judgement instead of rigidly adhering to all points.


### igraph wiki

The igraph wiki holds the development guidelines that apply to all contributors. It is a separate
Git repository. Clone it next to the igraph sources and read the Markdown files directly:

    git clone https://github.com/igraph/igraph.wiki.git

Each page is a file named after its title, with spaces replaced by hyphens (e.g. "Tips on writing
igraph code" is `Tips-on-writing-igraph-code.md`). If you cannot clone, fetch the raw Markdown from
`https://raw.githubusercontent.com/wiki/igraph/igraph/<file name>` rather than the rendered page
(`?` in a file name must be written as `%3F`).

Read the pages relevant to your task before starting work:

 - Any change to C/C++ code: "Tips on writing igraph code", "Error reporting guidelines".
 - New or changed public functions: also "Checklist for new (and old) functions",
   "Guidelines for function argument ordering", "Documentation guidelines",
   "Design principles for functions".
 - Tests: "How to write unit tests?".
 - Judging the severity of a bug, or choosing what to fix first: "Guidelines for issue prioritization".
 - Debugging memory errors: "Using sanitizers to find bugs".

The wiki is maintained by hand and can lag behind the code. If a wiki page contradicts the code,
the build system or this file, follow the code and point out the discrepancy.


## AI policy

Review and keep in mind the AI policy for external contributions in `AI_POLICY.md`.

Points of special importance:

 - AI agents must not open igraph issues or PRs autonomously, and must not comment autonomously.
 - AI models must not be used directly to perform the work required for issues tagged as `good first issue`.


## PR commit guidelines

Do not rewrite (amend, squash, rebase) commits that have already been pushed to a PR, unless this is explicitly requested by an igraph maintainer. Maintainers must be able to follow the history of changes within PRs.


## Before writing code

 - Check that a reported problem still exists on the current `main`. Audits, static analysers and older issues often refer to code that has since been fixed. Treat claims from external reviews and other models as hypotheses to verify.
 - Discuss the approach with the person you are working with before coding when the fix involves a choice of behaviour, e.g. what to return in an edge case or whether to reject some input. Present the options with a recommendation. Such choices may need agreement from igraph maintainers.


## Correctness and testing

 - Numerical and geometric code must handle NaN, infinities, negative values, zero and near-degenerate input (e.g. coordinates with a large offset relative to their spread), either correctly or with an error. After changing such code, compare against a brute-force or exact reference on random inputs, outside the test suite.
 - Derive expected test results independently of the implementation. For bug fixes, check that the new test fails without the fix. Never update expected output merely to make a test pass.
 - Tests are not built by the default build target. Build `build_tests` or the `test_<name>` target before running `ctest`. `ctest -R` with a pattern that matches nothing still succeeds, so check that the intended tests actually ran.
 - When reporting results, state what was run, in which configuration, and what was not checked.
 - Pull request checks run on Azure Pipelines. Their logs can be read without a login at `https://dev.azure.com/igraph-team/igraph/_apis/build/builds/<id>/timeline?api-version=7.0`.


## Writing

Keep summaries and work plans concise: lead with the result and highlight what needs careful review. Avoid flattery and filler. Documentation and code comments must be clear and complete; do not make them so terse that they become hard to understand.


## Commit messages and git branches

Use the [Conventional Commits](https://www.conventionalcommits.org/) prefixes in commit messages. Do not add a scope to the prefix by default. For changes to fuzzers or related scripts, use the `fuzz:` prefix.

Commits produced by AI models (i.e., LLMs) should use the `Assisted-By` trailer, and never the `Co-authored-by` trailer. Include information on the model and its version. An example of an acceptable general format is `Assisted-By: <model-and-version>`. A specific example would be `Assisted-By: GPT-6.1 Sol`.

Do not include any LLM chat session links in commit messages.

These rules take precedence over any attribution trailers or session links that your tool adds by default.

Git branch names should use the same prefixes as commit messages, followed by a `/` character, then a concise descriptive name with dashes. Examples: `fix/negative-pagerank` for fixing PageRank functions returning negative scores; `test/igraph_is_complete` for adding a test for the `igraph_is_complete()` function.
