
# How to Contribute

## TL;DR

Send your PR! Thanks!

## Slightly more Details

You want to contribute? That's very generous of you and we are truly
grateful! Small changes, like fixing typos in documentation are
completely fine and also most welcome. For bigger changes, we suggest
that you open an issue before you start coding, so that we can maximize
the probability that we can successfully merge in your code.

## Making Trivial Changes

* Fork the igraph repository at github: https://github.com/igraph/igraph
* Please always use the `develop` branch. Choose this branch in your
  fork.
* Then look for the file you want to modify. C sources and docs are in `src`.
  Python sources and docs are in `interfaces/python`. R sources and docs are in
  `interfaces/R`. 
* Click on the edit symbol (pen) on the upper right corner of the file view.
* Make your edits.
* Write a short commit message, less than 65 characters. E.g.
  "Fix C docs typo" or "Fix degree bug for loops". If needed, elaborate
  your changes below in the "extended description" field.
* Commit your changes.
* Go back to the start page of *your* forked repository. It is at
  `https://github.com/<username>/igraph`.
* Click on the green button before the branch name to create a pull request.
* Click on "Create pull request".
* Provide a more detailed description if you like. Please also indicate that
  you are fine with licensing your contribution under igraph's license (see
  Legal Stuff below).
* Click on "Create pull request".
* That's it! It is probably a good idea to keep your forked repository
  until the change is accepted into igraph, in case you need to modify it.
* Now you need to wait for us, unfortunately. Please ping us, if it takes
  long to respond. E.g. a week is considered to be long.
* Once your pull request is accepted, you can delete your forked repository.

## Making More Involved Changes

This is mostly the same as for trivial changes, but you probably want to
edit the sources on your computer, instead of online on Github.

* Open an issue in the issue tracker about the proposed changes.
  This is not required for smaller things, but I suggest you do it
  for others. Just in case somebody is already working on the same thing,
  or it is something we don't want in igraph.
* Fork the repository, and clone it to the machine you'll work on.
* We usually build igraph on OSX, so the `develop` branch is usally fine
  on that platform. It might have problems on other systems. If this
  happens, please open an issue and tell us.
* Make sure you work on the `develop` branch.
* Once ready with your changes, build igraph, and run the tests. For the C
  library this means

  ```sh
  ./bootstrap.sh
  ./configure
  make
  make check
  ```
  
  For the R package you need to do:
  
  ```sh
  cd interfaces/R
  make
  R CMD INSTALL igraph
  Rscript -e 'library(igraph); igraph_test()'
  ```

  For Python, you need to do:

  ```sh
  cd interfaces/python
  python setup.py test
  ```

  Note that if you modify C code that is used in R/Python (most C code is!), then
  you need to test these, too.
* Submit your pull request.
* Now you need to wait for us, unfortunately. Please ping us, if it takes
  long to respond. E.g. a week is considered to be long.

## Writing igraph Code 

Some tips on writing igraph code. In general, look at how things are done,
and try to do them similarly. (Unless you think they are not done well, in which
case please tell us.)

### R vs Python vs C

An enhancement that is specific to R or Python or C (or a subset) is perfectly
fine. You are not obliged to write wrappers in both R or Python. E.g. updating
plotting code in R or Python is also fine.

### Code Formatting

Look at the style (indentation, braces, etc.) of some recently committed
bigger change, and try to mimic that. The code style within igraph is not
stricly the same, but we want to keep it reasonably similar.

### C vs. C++

Try to use C, unless you are updating already existing C++ code, or
you have other good reason for C++ (but then maybe ask us first).

### Data types

Please try to use igraph's data types for vectors, matrices, stacks, etc.
If they lack some functionality you need, please tell us.

### Memory Allocation, Error Handling

Please use igraph's memory allocation functions. Please also use the
`FINALLY` stack: `IGRAPH_FINALLY`, `IGRAPH_FINALLY_CLEAN`, etc. See examples
in the C code.

### Random Numbers

Please look at how random numbers are generated in any function in `src/games.c`.
Do the same. I.e. use `RNG_BEGIN`, `RNG_END`, and igraph's RNG calls. Do
not use the libc RNGs or other RNGs.

### Documentation

Please document your new functions. The C documentation is included in the C
source code. The R documentation is now generated with `roxygen2`, and
it is included in the R source code. Python docs are in the Python source
files.

### Test Cases

Unless you change something trivial, please consider adding test cases,
either in C or R or Python (or all of them, of you contribute to all three).
This is important!

### Your name in the list of AUTHORS

There is a list of authors in the `AUTHORS` file in the R package. If you
add C or R code, please feel free to add your name here. Note that there might
be two copies of the same file, at least a the time of writing there are.

### Ask Us!

In general, if you are not sure about something, please ask! You can
open an issue on Github, write to the igraph-help mailing list (see the
homepage at http://igraph.org), or write to Tamás and Gábor. We preper
the public forums, though, because then others can learn from it, too.

## Legal Stuff

This is a pain to deal with, but we can't avoid it, unfortunately.
So, igraph is licensed under the "General Public License (GPL) version 2,
or later". The igraph manual is licensed under the "GNU Free Documentation License".
If your contribution is bigger than a typo fix, then please
indicate that you are fine with releasing your code/text under these licenses.
E.g. adding a sentence that reads as "I'm fine with GPL 2 or later and FDL."
is perfectly enough.
