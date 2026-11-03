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


## More is not always better

When producing code, documentation text, or work plans and summaries, keep readability by humans in mind. More is not always better.

Keep code your write maintainable and readable. Assume that human developers with limited time and capacity will be working on the code.

Keep summaries complete but concise. Highlight changes that require careful review.

Avoid flattery, avoid stating the obvious, avoid filler phrases, prioritize technical clarity.


## Commit messages and git branches

Use the [Conventional Commits](https://www.conventionalcommits.org/) prefixes in commit messages. Do not add a scope to the prefix by default. For changes to fuzzers or relates scripts, use the `fuzz:` prefix.

Commits produced by AI models (i.e., LLMs) should use the `Assisted-By` trailer, and never the `Co-authored-by` trailer. Include information on the model and its version. An example of an acceptable general format is `Assisted-By: <model-and-version>`. A specific example would be `Assisted-By: GPT-6.1 Sol`.

Do not include any LLM chat session links in commit messages.

Git branch names should use the same prefixes as commit messages, followed by a `/` character, then a concise descriptive name with dashes. Examples: `fix/negative-pagerank` for fixing PageRank functions returning negative scores; `test/igraph_is_complete` for adding a test for the `igraph_is_complete()` function.
