# Contributing to this project

Thank you for being interested in contributing to `igraph`! We need the help of
volunteers to keep the package going, so every little bit is welcome. You can help out
the project in several different ways.

This repository only hosts the C code of the `igraph` project. Even if you are not so
experienced with C, you can contribute in a number of ways:

1. Respond to user questions on our [support forum](https://igraph.discourse.group/).
2. Correct or improve our [documentation](https://igraph.org/c/html/latest/).
3. Go over [open issues](https://github.com/igraph/igraph/issues):
   - Are some older issues still relevant in the most recent version? If not, write a
     comment to the issue stating that you feel that the issue is not relevant any more.
   - Can you reproduce some of the bugs that are reported? If so, write a comment to
     the issue stating that this is still a problem in version X.
   - Some [issues point out problems with the documentation](https://github.com/igraph/igraph/labels/documentation);
     perhaps you could help correct these?
   - Some [issues require clarifying a mathematical problem, or some literature research](https://github.com/igraph/igraph/labels/theory),
     before any programming can begin. Can you contribute through your theoretical expertise?
   - Looking to contribute code? Take a look at some [good first issues](https://github.com/igraph/igraph/labels/good%20first%20issue).

## Using the issue tracker

- The issue tracker is the preferred channel for [bug reports](#bugs),
  [feature requests](#features) and [submitting pull requests](#pull-requests).

- Do you have a question? Please use our [igraph support forum](https://igraph.discourse.group)
  for support requests.

- Please keep the discussion on topic and respect the opinions of others, and
  adhere to our [Code of Conduct](https://igraph.org/code-of-conduct.html).

<a name="bugs"></a>
## Bug reports

A bug is a _demonstrable problem_ that is caused by the code in the repository.
Good bug reports are extremely helpful &mdash; thank you for reporting!

Guidelines for bug reports:

1. **Make sure that the bug is in the C code of igraph and not in one of the
   higher level interfaces** &mdash; if you are using igraph from R, Python
   or Mathematica, consider submitting your issue in
   [igraph/rigraph](https://github.com/igraph/rigraph/issues/new),
   [igraph/python-igraph](https://github.com/igraph/python-igraph/issues/new)
   or [szhorvat/IGraphM](https://github.com/szhorvat/IGraphM/issues/new)
   instead. If you are unsure whether your issue is in the C layer, submit
   a bug report in the repository of the higher level interface &mdash;
   we will transfer the issue here if it indeed affects the C layer.

2. **Use the GitHub issue search** &mdash; check if the issue has already been
   reported.

3. **Check if the issue has been fixed** &mdash; try to reproduce it using the
   latest `master` or development branch in the repository.

4. **Isolate the problem** &mdash; create a [short, self-contained, correct
   example](http://sscce.org/).

Please try to be as detailed as possible in your report and provide all
necessary information. What is your environment? What steps will reproduce the
issue? What would you expect to be the outcome? All these details will help us
to fix any potential bugs.

Example:

> Short and descriptive example bug report title
>
> A summary of the issue and the compiler/OS environment in which it occurs. If
> suitable, include the steps required to reproduce the bug.
>
> 1. This is the first step
> 2. This is the second step
> 3. Further steps, etc.
>
> `<url>` - a link to the reduced test case
>
> Any other information you want to share that is relevant to the issue being
> reported. This might include the lines of code that you have identified as
> causing the bug, and potential solutions (and your opinions on their
> merits).


<a name="features"></a>
## Feature requests

Feature requests are always welcome. First, take a moment to find out whether your
idea fits with the scope and aims of the project. Please provide as much detail
and context as possible, and where possible, references to relevant literature.
Having said that, implementing new features can be quite time consuming, and as
such they might not be implemented quickly. In addition, the development team
might decide not to implement a certain feature. It is up to you to make a case
to convince the project's developers of the merits of this feature.

<a name="pull-requests"></a>
## Pull requests

_**Note:** The wiki has a lot of useful information for newcomers, as well as a
[quick start guide](https://github.com/igraph/igraph/wiki/Quickstart-for-new-contributors)!_

Good pull requests - patches, improvements, new features - are a fantastic help.
They should remain focused in scope and avoid containing unrelated commits.
Please also take a look at our [tips on writing igraph code](#tips) before
getting your hands dirty.

**Please ask first** before embarking on any significant pull request (e.g.
implementing features, refactoring code, porting to a different language),
otherwise you risk spending a lot of time working on something that the
project's developers might not want to merge into the project.

Please adhere to the coding conventions used throughout a project (indentation,
accurate comments, etc.) and any other requirements (such as test coverage).

Follow the following steps if you would like to make a new pull request:

1. [Fork](http://help.github.com/fork-a-repo/) the project, clone your fork,
   and configure the remotes:

   ```bash
   # Clone your fork of the repo into the current directory
   git clone https://github.com/<your-username>/<repo-name>
   # Navigate to the newly cloned directory
   cd <repo-name>
   # Assign the original repo to a remote called "upstream"
   git remote add upstream https://github.com/<upstream-owner>/<repo-name>
   ```

2. Please checkout the section on [branching](#branching) to see whether you
   need to branch off from the `master` branch or the `develop` branch.

   If you cloned a while ago, get the latest changes from upstream:

   ```bash
   git checkout <dev-branch>
   git pull --rebase upstream <dev-branch>
   ```

3. Create a new topic branch (off the targeted branch, see
   [branching](#branching) section) to contain your feature, change, or fix:

   ```bash
   git checkout -b <topic-branch-name>
   ```

4. Please commit your changes in logical chunks, and try to provide clear commit
   messages. It helps us during the review process if we can follow your thought
   process during the implementation. If you hit a dead end, use `git revert`
   to revert your commits or just go back to an earlier commit with `git checkout`
   and continue your work from there.

5. We have a [checklist for new igraph functions](https://github.com/igraph/igraph/wiki/Checklist-for-new-(and-old)-functions).
   If you have added any new functions to igraph, please go through the
   checklist to ensure that your functions play nicely with the rest of the
   library.

6. Make sure that your PR is based off the latest code and locally merge (or
   rebase) the upstream development branch into your topic branch:

   ```bash
   git pull [--rebase] upstream <dev-branch>
   ```

   Rebasing is preferable over merging as you do not need to deal with merge
   conflicts; however, if you already have many commits, merging the upstream
   development branch may be faster.

7. WHen your topic branch is up-to-date with the upstream development branch, you can
   push your topic branch up to your fork:

   ```bash
   git push origin <topic-branch-name>
   ```

8. [Open a pull request](https://help.github.com/articles/using-pull-requests/)
    with a clear title and description.

**IMPORTANT**: By submitting a pull request, you agree to allow the project
owner to license your work under the same license as that used by the project,
see also [Legal Stuff](#legal).

<a name="branching"></a>
### Branching

`igraph` is committed to [semantic versioning](https://semver.org/). We are
currently still in the development release (0.x), which in principle is a mark
that the public API is not yet stable. Regardless, we try to maintain semantic
versioning also for the development releases. We do so as follows. Any released
minor version (0.x.z) will be API backwards-compatible with any previous release
of the *same* minor version (0.x.y, with y < z). This means that *if* there is
an API incompatible change, we will increase the minor version. For example,
release 0.8.1 is API backwards-compatible with release 0.8.0, while release
0.9.0 might be API incompatible with version 0.8.1. Note that this only concerns
the *public* API, internal functions may change also within a minor version.

There will always be two versions of `igraph`: the most recent released version,
and the next upcoming minor release, which is by definition not yet released.
The most recent release version is in the `master` branch, while the next
upcoming minor release is in the `develop` branch. If you make a change that is
API incompatible with the most recent release, it **must** be merged to
the `develop` branch. If the change is API backwards-compatible, it **can** be
merged to the `master` branch. It is possible that you build on recent
improvements in the `develop` branch, in which case your change should of course
target the `develop` branch. If you only add new functionality, but do not
change anything of the existing API, this should be backwards-compatible, and
can be merged in the `master` branch.

When you make a new pull request, please specify the correct target branch. The
maintainers of `igraph` may decide to retarget your pull request to the correct
branch. Retargeting you pull request may result in merge conflicts, so it is
always good to decide **before** starting to work on something whether you
should start from the `master` branch or from the `develop` branch. In most
cases, changes in the `master` branch will also be merged to the `develop`
branch by the maintainers.

If you are unsure about the branch to target, open an issue about your proposed
feature and we can discuss the appropriate target branch in the issue before
you send a PR.

<a name="tips"></a>
## Writing igraph Code

[Some tips on writing igraph code](https://github.com/igraph/igraph/wiki/Tips-on-writing-igraph-code).

## Ask Us!

In general, if you are not sure about something, please ask! You can
open an issue on GitHub, open a thread in our
[igraph support forum](https://igraph.discourse.group), or write to
[@ntamas](https://github.com/ntamas), [@vtraag](https://github.com/vtraag),
[@szhorvat](https://github.com/szhorvat), [@iosonofabio](https://github.com/iosonofabio) or
[@gaborcsardi](https://github.com/gaborcsardi).
We prefer open communication channels, because others can then learn from it
too.

<a name="legal"></a>
## Legal Stuff

This is a pain to deal with, but we can't avoid it, unfortunately.

`igraph` is licensed under the "General Public License (GPL) version 2, or
later". The igraph manual is licensed under the "GNU Free Documentation
License". By submitting a patch or pull request, you agree to allow the project
owner to license your work under the same license as that used by the project.
