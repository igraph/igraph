# Instructions for AI agents working on igraph

## Information on C/igraph development

Information on C/igraph development, design guidelines, notes on the use of development and debugging tools, is found in the igraph wiki at https://github.com/igraph/igraph/wiki/

In particular, make sure you read in full: Checklist for functions, Documentation guidelines, Error reporting guidelines, Guidelines for function argument ordering, How to write unit tests?

Check `.github/copilot-instructions.md` for how contributions are reviewed. When developing new code, treat the advice there as guidelines, and exercise judgement instead of rigidly adhering to all points.


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
