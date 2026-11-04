# Instructions for AI agents working on igraph

igraph is a C library for network analysis and graph theory. All work goes on a branch based on `main`; never commit to or otherwise modify `main` itself. `vendor/` holds third-party code and must not be modified.

## Design context

igraph is used mostly through its high-level interfaces (R, Python, Mathematica and others); the C library is their shared core, and its design follows from that:

 - Errors and warnings are translated into the host language's own mechanisms. Report them only through `IGRAPH_ERROR`, `IGRAPH_WARNING` and related macros, never by printing or aborting. Write messages that make sense from any host: describe the problem in graph terms, and do not quote vertex or edge IDs (some hosts index from 1) or C parameter names.
 - igraph runs in interactive sessions and notebooks. Long computations must be interruptible (`IGRAPH_ALLOW_INTERRUPTION()`) and must release all resources on interruption or error.
 - Interface maintainers rely on the function documentation to choose defaults and to know which parameters are optional. State these explicitly for new or changed public functions, and describe the function in `interfaces/functions.yaml`.
 - For a scientific library, a silently wrong result is the worst failure. When it is unclear how a function should behave on some input, reject the input with an error and raise the question, rather than guess.
 - Only documented behaviour is a promise. Undocumented details, such as the order of results or the output of a stochastic function for a given seed, may change when there is a reason to. Say so when a change affects them, and do not lock them into tests.

## Read first

 - `.github/copilot-instructions.md` is igraph's coding conventions, written as a review checklist. Read it before changing or reviewing C code, and apply it with judgement rather than rigidly.
 - The igraph wiki holds the development guidelines. Clone it with `git clone https://github.com/igraph/igraph.wiki.git` and read the Markdown files, named after the page titles with hyphens for spaces; if you cannot clone, fetch `https://raw.githubusercontent.com/wiki/igraph/igraph/<file name>` (`?` written as `%3F`). Read what applies to your task:

    - Always, for any C/C++ change (both short): "Tips on writing igraph code", "Error reporting guidelines".
    - New or changed public function: also "Contributing your first function to igraph", "Checklist for new (and old) functions", "Guidelines for function argument ordering", "Design principles for functions", "Documentation guidelines".
    - Tests: "How to write unit tests?". Changelog entries: "Changelog guidelines".
    - Deciding how serious a bug is: "Guidelines for issue prioritization".
    - Memory errors or crashes: "Using sanitizers to find bugs". Other build targets and CMake flags: "Developer documentation". Work in `fuzzing/`: "Fuzzing". Adding a benchmark to `tests/benchmarks/`: "Benchmarking igraph" and the existing files there.

   The wiki is maintained by hand and can lag behind the code. If it contradicts the code, the build system or this file, follow the code and point out the discrepancy.
 - The documentation lives in the sources: function documentation in comments above each definition in `src/`, chapters such as installation and prerequisites in `doc/*.xml`. Read it there; there is no need to build it.

## Build and test

Requirements: CMake 3.18 or later, C and C++ compilers, `bison`, `flex` and `git`. `libxml2` is optional (GraphML support); every other dependency is bundled, but the first build is much faster with external libraries, which CMake uses automatically when they are installed. On Debian or Ubuntu:

    apt-get install bison flex ninja-build ccache libxml2-dev libgmp-dev libglpk-dev libarpack2-dev libopenblas-dev

The version number is derived from `git describe`, so a shallow clone without tags fails to configure. Check with `git rev-parse --is-shallow-repository` and fix with `git fetch --unshallow --tags` before building. Full history is also needed before using `git log`, `blame` or `merge-base` to reason about what changed.

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DUSE_SANITIZER="Address;Undefined"
    cmake --build build --target build_tests
    ctest --test-dir build -R <name> --output-on-failure

 - The default target builds only the library; build `build_tests`, or `test_<name>` for one test, before `ctest`. A `ctest -R` pattern that matches nothing still succeeds, so check that the intended tests ran.
 - Warnings are errors by default. Fix the warning rather than turning this off.
 - Use `build` or `build-*` as the build directory name; these are ignored by Git. Re-run the configure step after adding a `.out` file.
 - Timings from a Debug or sanitizer build are meaningless. Measure in a separate build with `-DCMAKE_BUILD_TYPE=Release -DIGRAPH_ENABLE_LTO=ON -DBUILD_SHARED_LIBS=ON` and no `USE_SANITIZER`. (LTO matters for some functions; the shared library keeps the benchmark's own code out of the optimisation and keeps link times short.) For a before/after comparison, write the benchmark first, build its `benchmark_<name>` target in a second build directory such as `build-base` before changing any code, and leave that directory untouched; after the change, build the same target in your working build directory and run both executables on the same machine. Temporary measurement programs are not committed; long-lived benchmarks go in `tests/benchmarks/`.
 - Pull request checks run on Azure Pipelines and GitHub Actions. Look at the `linux_static_vendored` Azure job first; it is the most stringent (sanitizers, locale tests). `linux_static_vendored_32` uses 32-bit `igraph_int_t`; the Stimulus action checks `interfaces/functions.yaml` against the C prototypes. Azure results are readable without a login: `https://dev.azure.com/igraph-team/igraph/_apis/build/builds?api-version=7.0&branchName=refs/pull/<PR>/merge` lists a PR's build, and `.../builds/<id>/timeline?api-version=7.0` lists its jobs and logs.

## Working with people

`AI_POLICY.md` is addressed to human contributors; following it is their responsibility, not yours. Two points affect how you work:

 - Leave communication to the person you are working with. Do not open issues or pull requests, post comments, or push commits unless they ask you to for that specific case; prepare the material and hand it to them, and they decide what is submitted, in their own words. A review that a maintainer has requested from a review tool is such a case.
 - If asked to work on an issue labelled `good first issue`, say before starting that the policy asks contributors to do these themselves as a learning exercise, then continue unless told otherwise.

Do not rewrite (amend, squash, rebase) commits that have already been pushed to a PR unless an igraph maintainer explicitly asks for it. Maintainers must be able to follow the history of changes within PRs.

## Working on a change

 - Work happens in forks. Before starting, fetch `main` from `https://github.com/igraph/igraph` and tell the person if the branch you are on is behind it; do not merge or rebase unless asked.
 - Reproduce a reported problem on current `main` before fixing it. Old issues and analyser reports often describe code that has since changed. Report a suspected bug only with a minimal reproduction and an independent argument for what the correct result is.
 - When a fix involves a choice of behaviour, e.g. what to return in an edge case or whether to reject some input, present the options with a recommendation before coding; such choices may need agreement within the igraph team. If you cannot ask, take the option most consistent with existing igraph behaviour and list every such choice in your report.
 - Match the formatting of the surrounding code; do not run formatters or reformat code you did not change.
 - Derive expected test results independently of the implementation. For bug fixes, check that the new test fails without the fix. Never update expected output merely to make a test pass.
 - Propose a `CHANGELOG.md` entry under `[main]` for every user-visible change, following "Changelog guidelines" in the wiki.
 - When asked to review a pull request, fetch it with `git fetch https://github.com/igraph/igraph pull/<N>/head`, apply the review checklist, and report findings to the person in order of importance, naming the checklist item or wiki page each one rests on.
 - Run the full test suite before you report done. When reporting results, state what was run, in which configuration, and what was not checked.

## Writing

Keep summaries and work plans concise: lead with the result and highlight what needs careful review. Avoid flattery and filler. Documentation and code comments must be clear and complete; do not make them so terse that they become hard to understand.

## Commit messages and git branches

Use [Conventional Commits](https://www.conventionalcommits.org/) prefixes without a scope; igraph additionally uses `fuzz:` for fuzzers and related scripts, `bench:` for benchmarks, and `chore:` for changelog updates. Make one logical change per commit; keep bug fixes separate from performance work and refactoring.

The commit author is the person you are working with, never an AI model, even if your tool configures a model identity by default. If no Git identity is configured, say so in your report rather than inventing one.

Commits produced with the help of AI models must carry an `Assisted-by` trailer naming the model and its version as reported by your tool, e.g. `Assisted-by: GPT-6.1 Sol`. Do not use the `Co-authored-by` trailer for AI models, and do not include LLM chat session links in commit messages. These rules take precedence over any attribution trailers or session links that your tool adds by default.

Branch names use the same prefixes followed by `/` and a concise dash-separated name: `fix/negative-pagerank`, `test/igraph_is_complete`.
