There are two ways to compile igraph on Windows. You can either build from the released source code, or you can build from the `git` repository. Compiling from the released source code is sometimes easier, as it already provides a project file that can be used to compile with Microsoft Visual C++ (MSVC). Of course, compilation from the release source code is limited to released versions, while compilation from the `git` repository allows you to also compile development versions. Compilation can always be done using MinGW, but if you want to use `igraph` with Python (e.g. as an external library for `python-igraph`) you have to use MSVC. Each Python version requires a specific version of MSVC, please see the Python [website](https://wiki.python.org/moin/WindowsCompilers) for more details. Preparing a project file to compile with MSVC is again done using MinGW.

In summary, there are two ways of compiling `igraph`:

1. Compile using `msys2`/`cygwin`
2. Compile using Microsoft Visual C++
   - If no project file is provided, this should be generated using `msys2`/`cygwin`.

More detailed instructions are provided below.

# Compilation using `msys2`/`cygwin`

## Configuration of `msys2`

We recommend to use [`msys2`](https://www.msys2.org/) to compile `igraph`, because it support both 32-bits and a [64-bits MinGW](http://mingw-w64.org/), whereas the old [MinGW](http://mingw.org/) project only provides a 32-bits version. After installation, the compilation of `igraph` requires the installation of a number of packages, namely:

- `autoconf`
- `automake`
- `bison`
- `flex`
- `libtool`
- `libxml2-devel`
- `make`
- `zip`
- `mingw-w64-x86_64-toolchain`

You can install these packages by running the following commands in an open `msys2` `bash` terminal:

```
pacman --needed --noconfirm -Sy pacman-mirrors
pacman --noconfirm -Sy
pacman --noconfirm -S autoconf automake bison flex libtool libxml2-devel make zip mingw-w64-x86_64-toolchain
```

Running these commands may require administrator rights. In that case, be sure to run the `msys2` `bash` terminal as administrator.

Finally, the `gcc` compiler from `minGW` may possibly not be found correctly, in which case the path needs to be corrected, which can be done by running `export PATH="/mingw64/bin:$PATH"`.
`

## Configuration of `cygwin`

You can use [`cygwin`](https://www.cygwin.com/) to compile `igraph`. The compilation of `igraph` requires the installation of a number of packages, namely:

- `automake`
- `autoconf2.5`
- `libtool`
- `flex`
- `bison`
- `libxml2-devel`
- `gcc-g++`
- `gcc-core`
- `make`
- `libgmp-devel`
- `util-linux`
- `git`
- `zip`

You can install `cygwin` including these packages by running the following command from a Windows Command Prompt:

```
setup-x86_64.exe -q -P automake,autoconf2.5,libtool,flex,bison,libxml2-devel,gcc-g++,gcc-core,make,libgmp-devel,util-linux,git,zip
```

Make sure that you are in the directory that contains `setup-x86_64.exe`.

## Compilation

You can either download the [released source code](https://igraph.org/c/#downloads) or clone  the [`git` repository](https://github.com/igraph/igraph). When cloning from GitHub, please ensure that you run `git config --global core.autocrlf false` before cloning the repository. Otherwise, you might run into problems with line endings. After extracting the source code you should open the `bash` terminal (from either `msys2` or `cygwin`) and change to the directory in which the source code is located.

If you have cloned the source code from the `git` repository, you first have create the build scripts by running

```
./bootstrap.sh
```

Compilation is done with the typical steps of

```
./configure
make
```

This built `igraph` library is then located at XXX.

## Preparing the MSVC project file

If you intend to compile with Microsoft Visual C++, you will need a project file for Microsoft Visual C++. This can be done by following the above steps, except for the last step, you create the MSVC project file using

```
./configure
make msvc
```

The necessary source code for compilation with MSVC is then contained in the newly created `igraph-[VERSION]-msvc` directory, where `[VERSION]` refers to the version that you are compiling.

# Compilation using MSVC

Compilation using MSVC is primarily necessary for use with Python. Note that different versions of Python require different versions of MSVC, see the Python [website](https://wiki.python.org/moin/WindowsCompilers) for more details.

The provided project file is still provided in an older MSVC format. To use it with more recent versions of MSVC, you will have to upgrade the project file. This can be done in various ways. If you dispose of the full Microsoft Visual Studio, you can simply open the project file in the Visual Studio IDE. It will ask you to convert it to a modern format, and you should be able to compile it.

If you do not have the full Microsoft Visual Studio, but only the standalone tools (as indicated on the Python [website](https://wiki.python.org/moin/WindowsCompilers)), you should follow a slightly different route. In that case, open the "Developer Command Prompt", and depending on your version of MSVC, you should do the following

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

# Compilation problems

If you have any problems with compilation on Windows, please post a message on XXX or file an issue at https://github.com/igraph/igraph.