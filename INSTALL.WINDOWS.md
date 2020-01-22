There are two ways to compile igraph on Windows. You can either build from the released source code, or you can build from the `git` repository. Compiling from the released source code is sometimes easier, as it already provides a project file that can be used to compile with Microsoft Visual C++ (MSVC). Of course, compilation from the released source code is limited to released versions, while compilation from the `git` repository allows you to also compile development versions. Compilation can always be done using MinGW or `cygwin`, but if you want to use `igraph` with Python (e.g. as an external library for `python-igraph`) you have to use MSVC. Each Python version requires a specific version of MSVC, please see the Python [website](https://wiki.python.org/moin/WindowsCompilers) for more details. Preparing a project file to compile with MSVC is again done using MinGW (`msys2`) or `cygwin`.

In summary, there are two routes to compiling `igraph`:

1. Compile using `msys2`/`cygwin`
2. Compile using Microsoft Visual C++
   - If no project file is provided, this should be generated using `msys2`/`cygwin`.

More detailed instructions are provided below.

# Compilation using `msys2`/`cygwin`

## Configuration of `msys2`

We recommend to use [`msys2`](https://www.msys2.org/) to compile `igraph`, because it includes [MinGW-w64](http://mingw-w64.org/), which supports compiling both for 32- and 64-bit targets. The old [MinGW](http://mingw.org/) project only provides a 32-bit version.

The instructions below assume that you want to compile for a 64-bit target.

After installing `msys2`, first make sure that it is up-to-date. You can do so by opening an `msys2` `bash` terminal and run:

```
pacman -Syuu
```

You may be requested to close and restart the `msys2` `bash` terminal at several points. Please repeat executing `pacman -Syuu` until it says there are no longer any packages to update. More details regarding the installation of `msys2` are available from their [Wiki](https://github.com/msys2/msys2/wiki/MSYS2-installation).

Compilation of `igraph` requires the installation of a number of packages, namely:

- `autoconf`
- `automake`
- `bison`
- `flex`
- `libtool`
- `make`
- `mingw-w64-x86_64-toolchain`
- `mingw-w64-x86_64-libxml2`
- `git`
- `zip`

The `mingw-w64-x86_64-toolchain` also includes the optional [Gnu Multiple Precision](https://gmplib.org/) package `mingw-w64-x86_64-gmp`. This enables extended functionality in igraph.

You can install these packages by running the following command in an `msys2` terminal:

```
pacman -S autoconf automake bison flex libtool make mingw-w64-x86_64-toolchain mingw-w64-x86_64-libxml2 mingw-w64-x86_64-gmp git zip
```

Once everything is installed you need to open the `msys2` "MinGW 64-bit" terminal (`mingw64.exe`).

## Configuration of `cygwin`

You can use [`cygwin`](https://www.cygwin.com/) to compile `igraph`. The 32-bits version should not be used (as recommended by `cygwin` also), so make sure to use the 64-bits version. The compilation of `igraph` requires the installation of a number of packages, namely:

- `autoconf2.5`
- `automake`
- `bison`
- `flex`
- `gcc-core`
- `gcc-g++`
- `git`
- `libgmp-devel`
- `libtool`
- `libxml2-devel`
- `make`
- `util-linux`
- `zip`

You can install `cygwin` including these packages by running the following command from a Windows Command Prompt (make sure that you are in the directory that contains `setup-x86_64.exe`.):

```
setup-x86_64.exe -q -P autoconf2.5,automake,bison,flex,gcc-core,gcc-g++,git,libgmp-devel,libtool,libxml2-devel,make,util-linux,zip
```

You will be asked to choose a mirror, simply choose one, preferably close by. You will then be given an overview of packages to install, you can simply press Next (two times).

## Compilation

You can either download the [source code archive of the latest release](https://igraph.org/c/#downloads) or clone the [`git` repository](https://github.com/igraph/igraph) using `git clone https://github.com/igraph/igraph.git`. When cloning from GitHub, please ensure that you run `git config --global core.autocrlf false` before cloning the repository. Otherwise, you might run into problems with line endings.

After extracting the source code, open a `bash` terminal (from either `msys2` or `cygwin`) and change to the directory in which the source code is located. When using `msys2`, use the "MSYS2 MinGW 64-bit" terminal, and *not* the "MSYS2 MSYS" one.

If you have cloned the source code from the `git` repository, you must first create the build scripts by running

```
./bootstrap.sh
```

If you are compiling using MinGW (using `msys2`), do not forget to set the `CPPFLAGS` before running `./configure`:
```
export CPPFLAGS="$CPPFLAGS -DMSDOS"
```

Compilation is done with the usual steps of

```
./configure
make
```

If you have multiple cores, for example 4, you can run `make -j4` to speed up the compilation substantially.

You can install the library by running `make install`. If you prefer to install it in a specific location, you can specify `--prefix=[DIRECTORY]` to `./configure`, in which case the library will be installed into the `bin`, `lib` and `include` subdirectories of `[DIRECTORY]`. If you just need the built `dll` file, this is located in `src/.libs/libigraph-0.dll`.

You may optionally run the test suite using `make check`. Before doing this, ensure that the directory containing `libigraph-0.dll` is in the `PATH`. This can be achieved for example using the command

```
export PATH=$PATH:`realpath src/.libs`
```

## Preparing the MSVC project file

If you intend to compile with Microsoft Visual C++, you will need a project file for Microsoft Visual C++. This can be done by following the above steps, except for the last step, you create the MSVC project file using

```
./configure
make msvc
```

The necessary source code for compilation with MSVC is then contained in the newly created `igraph-[VERSION]-msvc` directory, where `[VERSION]` refers to the version that you are compiling.

# Compilation using MSVC

Compilation using MSVC is primarily necessary for use with Python, in particular for `python-igraph`. Note that different versions of Python require different versions of MSVC. See [the Python website](https://wiki.python.org/moin/WindowsCompilers) for more details.

The provided project file is still provided in an older MSVC format. To use it with more recent versions of MSVC, you will have to upgrade the project file. Note that this can only be done if you have the the full Microsoft Visual Studio install. You can install the Community edition, which is [freely available](https://visualstudio.microsoft.com/downloads/). The Build Tools that are available do not allow you to upgrade the project file, which is necessary to be able to compile `igraph` correctly.

Compilation can be done in two ways. The easiest way for most users is probably simply to open the project solution file (`igraph.sln`) in the Visual Studio IDE. It will ask you to convert it to the current format, and you should be able to compile it, once that is done.

Alternatively, you can do this via the command line. In that case, open the "Developer Command Prompt", and depending on your version of MSVC, you should do the following

MSVC  9.0
```
vcbuild.exe /upgrade
vcbuild.exe igraph.vcproj "Release|x64"
```
MSVC 10.0
```
VCUpgrade.exe /overwrite igraph.vcproj
msbuild.exe igraph.vcxproj
```

MSVC 14.x
```
devenv /upgrade igraph.vcproj
msbuild.exe igraph.vcxproj
```

By default, compilation on Windows targets a static library (`.lib`), which is necessary for the compilation of `python-igraph`. To import the static library, you have to make sure to specify
```
#define IGRAPH_NO_IMPORTS 1
```
before you include
```
#include "igraph.h"
```

Alternatively, you can also reconfigure the MSVC project to compile a dynamic library (`.dll`), in which case the `IGRAPH_NO_IMPORTS` definition should *not* be included.

A small test project is included in the `igraphtest` directory. It does not (yet) include the `#define IGRAPH_NO_IMPORTS 1` so you have to include that yourself.

# Compilation problems

If you have any problems with compilation on Windows, please post a message on [the support forum](https://igraph.discourse.group/) or file an issue at https://github.com/igraph/igraph.
