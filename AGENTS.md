# Instructions for AI agents working on igraph

### Information on C/igraph development

Information on C/igraph development, design guidelines, notes on the use of development and debugging tools, is found in the igraph wiki at https://github.com/igraph/igraph/wiki/

In particular, make sure you read in full: Checklist for functions, Documentation guidelines, Error reporting guidelines, Guidelines for function argument ordering, How to write unit tests?

Check `.github/copilot-instructions.md` for how contributions are reviewed. When developing new code, treat the advice there as guidelines, and exercise judgement instead of rigidly adhering to all points.


### AI policy

Review and keep in mind the AI policy for external contributions in `AI_POLICY.md`.

Points of special importance:

 - AI agents must not open igraph issues or PRs autonomously, and must not comment autonomously.
 - AI models must not be used directly to perform the work required for issues tagged as `good first issue`.


### Guidelines on issue prioritization

igraph is a scientific library for network analysis and graph theory. It is written in C, but mainly used through its high-level interfaces (Python, R and Mathematica, etc.) This determines how we view the seriousness of various issues.

1. Silent wrong results are the worst sin a scientific library can commit. They may lead researchers to publish incorrect results.
2. Crashes (such as a segmentation fault) should be fixed, especially since they can take down the session of a user when using igraph through a high-level interface. However, when a crash occurs, it is immediately clear to the user that something went wrong. This is not as serious as a silent wrong result.
3. Problems that may occur during plausible igraph usage are higher priority than those that can only be triggered by exceedingly unlikely or specially crafted input.

Deciding on fixes to certain issues involve compromises. In these cases we consider long-term maintainability, code complexity, and performance along with the seriousness of the issue.


### PR commit guidelines

Do not rewrite (amend, squash, rebase) commits that have already been pushed to a PR, unless this is explicitly requested by a maintainer with commmit access to the igraph repository. Maintainers must be able to follow the history of changes within PRs.


### More is not always better

When producing code, documentation text, or work plans and summaries, keep readability by humans in mind. More is not always better.

Keep code your write maintainable and readable. Assume that human developers with limited time and capacity will be working on the code.

Keep summaries complete but concise. Highlight changes that require careful review.

Avoid flattery, avoid stating the obvious, avoid filler phrases, prioritize technical clarity.


### Commit message standards for AI agents

Use standard prefixes in commit messages, such as `fix:`, `feat:`, `docs:`, `test:`, etc., following the convention set by prior commits.

Commits produced by AI models (i.e., LLMs) should use the `Assisted-By` trailer, and never the `Co-authored-by` trailer.

Include information on the model, its version, any any relevant analysis or refactoring tools that were used (such as clang-tidy, astyle, libfuzzer, etc.). Do not include basic development tools (such as editors, make, gcc, git, etc.)

An example of an acceptable general format is the following:

```
Assisted-By: <model-and-version>, <tools-used>
```

A specific example would be:

```
Assisted-By: GPT-6.1 Sol, clang-tidy
```

Do not include any LLM chat session links in commit messages.

These standards are based on those used by the Linux kernel, see https://docs.kernel.org/process/coding-assistants.html
