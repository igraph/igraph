#!/usr/bin/env python
try:
    from setuptools import setup
    from setuptools.command.build_ext import build_ext
    build_py = None
except ImportError:
    from distutils.core import setup
    try:
        from distutils.command.build_py import build_py_2to3 as build_py
    except ImportError:
        from distutils.command.build_py import build_py

import atexit
import distutils.ccompiler
import glob
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib2

from distutils.core import Extension
from distutils.util import get_platform
from select import select
from subprocess import Popen, PIPE
from textwrap import dedent
from urllib import urlretrieve

# Global version number
VERSION = "0.7"

# Check Python's version info and exit early if it is too old
if sys.version_info < (2, 5):
    print("This module requires Python >= 2.5")
    sys.exit(0)

###########################################################################

LIBIGRAPH_FALLBACK_INCLUDE_DIRS = ['/usr/include/igraph', '/usr/local/include/igraph']
LIBIGRAPH_FALLBACK_LIBRARIES = ['igraph']
LIBIGRAPH_FALLBACK_LIBRARY_DIRS = []

###########################################################################

def cleanup_tmpdir(dirname):
    """Removes the given temporary directory if it exists."""
    if dirname is not None and os.path.exists(dirname):
        shutil.rmtree(dirname)

def find_static_library(library_name, library_path):
    """Given the raw name of a library in `library_name`, tries to find a
    static library with this name in the given `library_path`. `library_path`
    is automatically extended with common library directories on Linux and Mac
    OS X."""

    variants = ["lib{0}.a", "{0}.a", "{0}.lib", "lib{0}.lib"]
    if is_unix_like():
        extra_libdirs = ["/usr/local/lib64", "/usr/local/lib",
                "/usr/lib64", "/usr/lib", "/lib64", "/lib"]
    else:
        extra_libdirs = []

    for path in extra_libdirs:
        if path not in library_path and os.path.isdir(path):
            library_path.append(path)

    for path in library_path:
        for variant in variants:
            full_path = os.path.join(path, variant.format(library_name))
            if os.path.isfile(full_path):
                return full_path

def first(iterable):
    """Returns the first element from the given iterable."""
    for item in iterable:
        return item
    raise ValueError("iterable is empty")

def get_output(command):
    """Returns the output of a command returning a single line of output."""
    p = Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
    p.stdin.close()
    p.stderr.close()
    line=p.stdout.readline().strip()
    p.wait()
    if type(line).__name__ == "bytes":
        line = str(line, encoding="utf-8")
    return line, p.returncode

def http_url_exists(url):
    """Returns whether the given HTTP URL 'exists' in the sense that it is returning
    an HTTP error code or not. A URL is considered to exist if it does not return
    an HTTP error code."""
    class HEADRequest(urllib2.Request):
        def get_method(self):
            return "HEAD"
    try:
        response = urllib2.urlopen(HEADRequest(url))
        return True
    except urllib2.URLError:
        return False

def is_unix_like(platform=None):
    """Returns whether the given platform is a Unix-like platform with the usual
    Unix filesystem. When the parameter is omitted, it defaults to ``sys.platform``
    """
    platform = platform or sys.platform
    platform = platform.lower()
    return platform.startswith("linux") or platform.startswith("darwin") or \
            platform.startswith("cygwin")

def preprocess_args():
    """Preprocesses the command line options before they are passed to
    setup.py"""
    # Yes, this is ugly, but we don't want to interfere with setup.py's own
    # option handling
    opts_to_remove = []
    for idx, option in enumerate(sys.argv):
        if not option.startswith("--"):
            continue
        if option == "--static":
            opts_to_remove.append(idx)
            buildcfg.static_extension = True
        elif option == "--no-download":
            opts_to_remove.append(idx)
            buildcfg.download_igraph_if_needed = False
        elif option == "--no-pkg-config":
            opts_to_remove.append(idx)
            buildcfg.use_pkgconfig = False
        elif option.startswith("--c-core-version"):
            opts_to_remove.append(idx)
            if option == "--c-core-version":
                value = sys.argv[idx+1]
                opts_to_remove.append(idx+1)
            else:
                value = option.split("=", 1)[1]
            buildcfg.c_core_versions = [value]

    for idx in reversed(opts_to_remove):
        sys.argv[idx:(idx+1)] = []

def preprocess_fallback_config():
    """Preprocesses the fallback include and library paths depending on the
    platform."""
    global LIBIGRAPH_FALLBACK_INCLUDE_DIRS
    global LIBIGRAPH_FALLBACK_LIBRARY_DIRS
    global LIBIGRAPH_FALLBACK_LIBRARIES

    if os.name == 'nt' and distutils.ccompiler.get_default_compiler() == 'msvc':
        # if this setup is run in the source checkout *and* the igraph msvc was build,
        # this code adds the right library and include dir
        all_msvc_dirs = glob.glob(os.path.join('..', '..', 'igraph-*-msvc'))
        if len(all_msvc_dirs) > 0:
            if len(all_msvc_dirs) > 1:
                print("More than one MSVC build directory (..\\..\\igraph-*-msvc) found!")
                print("It could happen that setup.py uses the wrong one! Please remove all but the right one!\n\n")

            msvc_builddir = all_msvc_dirs[-1]
            if not os.path.exists(os.path.join(msvc_builddir, "Release")):
                print("There is no 'Release' dir in the MSVC build directory\n(%s)" % msvc_builddir)
                print("Please build the MSVC build first!\n")
            else:
                print("Using MSVC build dir as a fallback: %s\n\n" % msvc_builddir)
                LIBIGRAPH_FALLBACK_INCLUDE_DIRS = [os.path.join(msvc_builddir, "include")]
                LIBIGRAPH_FALLBACK_LIBRARY_DIRS = [os.path.join(msvc_builddir, "Release")]

def version_variants(version):
    """Given an igraph version number, returns a list of possible version
    number variants to try when looking for a suitable nightly build of the
    C core to download from igraph.org."""

    result = [version]

    # Add trailing ".0" as needed to ensure that we have at least
    # major.minor.patch
    parts = version.split(".")
    while len(parts) < 3:
        parts.append("0")
        result.append(".".join(parts))

    return result

###########################################################################

class IgraphCCoreBuilder(object):
    """Class responsible for downloading and building the C core of igraph
    if it is not installed yet."""

    def __init__(self, versions_to_try):
        self.versions_to_try = versions_to_try
        self._builddir = None
        self._tmpdir = None

    @property
    def builddir(self):
        """The directory in which igraph is built."""
        if self._builddir is None:
            self._builddir = os.path.join(self.tmpdir, "igraph-%s" % self.version)
        return self._builddir

    @property
    def tmpdir(self):
        """The temporary directory in which igraph is downloaded and extracted."""
        if self._tmpdir is None:
            self._tmpdir = tempfile.mkdtemp(prefix="igraph.")
            atexit.register(cleanup_tmpdir, self._tmpdir)
        return self._tmpdir

    def download_and_compile(self):
        """Downloads and compiles the C core of igraph."""

        def _progress_hook(count, block_size, total_size):
            if total_size < 0:
                sys.stdout.write("\rDownloading %s... please wait." % local_file)
            else:
                percentage = count * block_size * 100.0 / total_size
                percentage = min(percentage, 100.0)
                sys.stdout.write("\rDownloading %s... %.2f%%" % (local_file, percentage))
            sys.stdout.flush()

        # Determine the remote URL
        version, remote_url = self.find_first_version()
        if not version:
            print("Version %s of the C core of igraph is not found among the "
                    "nightly builds." % self.versions_to_try[0])
            print("Use the --c-core-version switch to try a different version.")
            print("")
            return False

        # Now determine the local filename of the C core
        local_file = "igraph-%s.tar.gz" % version
        local_file_full_path = os.path.join(self.tmpdir, local_file)

        # Download the C core
        urlretrieve(remote_url, local_file_full_path, reporthook=_progress_hook)
        print("")

        # Extract it in the temporary directory
        print("Extracting %s..." % local_file)
        archive = tarfile.open(local_file_full_path, "r:gz")
        archive.extractall(self.tmpdir)

        # Try to compile
        cwd = os.getcwd()
        try:
            print("Configuring igraph...")
            os.chdir(self.builddir)
            retcode = subprocess.call("./configure", shell=True)
            if retcode:
                return False

            retcode = subprocess.call("make", shell=True)
            if retcode:
                return False

            libraries = []
            for line in open(os.path.join(self.builddir, "igraph.pc")):
                if line.startswith("Libs: ") or line.startswith("Libs.private: "):
                    words = line.strip().split()
                    libraries.extend(word[2:] for word in words if word.startswith("-l"))

            if not libraries:
                # Educated guess
                libraries = ["igraph"]

            self.include_dirs = [os.path.join(self.builddir, "include")]
            self.library_dirs = [os.path.join(self.builddir, "src", ".libs")]
            self.libraries = libraries
            return True

        finally:
            os.chdir(cwd)

        return False

    def find_first_version(self):
        """Finds the first version of igraph that exists in the nightly build
        repo from the version numbers provided in ``self.versions_to_try``."""
        for version in self.versions_to_try:
            remote_url = self.get_download_url(version=version)
            if http_url_exists(remote_url):
                return version, remote_url
        return None, None

    def get_download_url(self, version):
        return "http://igraph.org/nightly/get/c/igraph-%s.tar.gz" % version

    def run(self):
        return self.download_and_compile()


class BuildConfiguration(object):
    def __init__(self):
        global VERSION
        self.c_core_versions = version_variants(VERSION)
        self.include_dirs = []
        self.library_dirs = []
        self.libraries = []
        self.extra_compile_args = []
        self.extra_link_args = []
        self.extra_objects = []
        self.static_extension = False
        self.download_igraph_if_needed = True
        self.use_pkgconfig = True
        self._has_pkgconfig = None

    @property
    def has_pkgconfig(self):
        """Returns whether ``pkg-config`` is available on the current system
        and it knows about igraph or not."""
        if self._has_pkgconfig is None:
            if self.use_pkgconfig:
                line, exit_code = get_output("pkg-config igraph")
                self._has_pkgconfig = (exit_code == 0)
            else:
                self._has_pkgconfig = False
        return self._has_pkgconfig

    @property
    def build_ext(self):
        """Returns a class that can be used as a replacement for the
        ``build_ext`` command in ``distutils`` and that will download and
        compile the C core of igraph if needed."""
        try:
            from setuptools.command.build_ext import build_ext
        except ImportError:
            from distutils.command.build_ext import build_ext

        buildcfg = self
        class custom_build_ext(build_ext):
            def run(self):
                # Print a warning if pkg-config is not available or does not know about igraph
                if buildcfg.use_pkgconfig:
                    detected = buildcfg.detect_from_pkgconfig()
                else:
                    detected = False

                # Download and compile igraph if the user did not disable it and
                # we do not know the libraries from pkg-config yet
                if not detected:
                    if buildcfg.download_igraph_if_needed and is_unix_like():
                        detected = buildcfg.download_and_compile_igraph()

                # Fall back to an educated guess if everything else failed
                if not detected:
                    buildcfg.use_educated_guess()

                # Replaces library names with full paths to static libraries
                # where possible
                if buildcfg.static_extension:
                    buildcfg.replace_static_libraries()

                # Prints basic build information
                buildcfg.print_build_info()

                ext = first(extension for extension in self.extensions
                        if extension.name == "igraph._igraph")
                buildcfg.configure(ext)

                # Run the original build_ext command
                build_ext.run(self)

        return custom_build_ext

    def configure(self, ext):
        """Configures the given Extension object using this build configuration."""
        ext.include_dirs = self.include_dirs
        ext.library_dirs = self.library_dirs
        ext.libraries = self.libraries
        ext.extra_compile_args = self.extra_compile_args
        ext.extra_link_args = self.extra_link_args
        ext.extra_objects = self.extra_objects

    def detect_from_pkgconfig(self):
        """Detects the igraph include directory, library directory and the
        list of libraries to link to using ``pkg-config``."""
        if not buildcfg.has_pkgconfig:
            print("Cannot find the C core of igraph on this system using pkg-config.")
            return False

        cmd = "pkg-config igraph --cflags --libs"
        if self.static_extension:
            cmd += " --static"
        line, exit_code = get_output(cmd)
        if exit_code > 0 or len(line) == 0:
            return False

        opts = line.strip().split()
        self.libraries = [opt[2:] for opt in opts if opt.startswith("-l")]
        self.library_dirs = [opt[2:] for opt in opts if opt.startswith("-L")]
        self.include_dirs = [opt[2:] for opt in opts if opt.startswith("-I")]
        return True

    def download_and_compile_igraph(self):
        """Downloads and compiles the C core of igraph."""
        print("We will now try to download and compile the C core from scratch.")
        print("Version number of the C core: %s" % self.c_core_versions[0])
        if len(self.c_core_versions) > 1:
            print("We will also try: %s" % ", ".join(self.c_core_versions[1:]))
        print("")

        igraph_builder = IgraphCCoreBuilder(self.c_core_versions)
        if igraph_builder.run():
            self.include_dirs = igraph_builder.include_dirs
            self.library_dirs = igraph_builder.library_dirs
            self.libraries = igraph_builder.libraries
            self.static_extension = True
            return True
        else:
            print("Could not download and compile the C core of igraph.")
            print("")
            return False

    def print_build_info(self):
        """Prints the include and library path being used for debugging purposes."""
        if self.static_extension:
            build_type = "static extension"
        else:
            build_type = "dynamic extension"
        print("Build type: %s" % build_type)
        print("Include path: %s" % " ".join(self.include_dirs))
        print("Library path: %s" % " ".join(self.library_dirs))
        print("Linked dynamic libraries: %s" % " ".join(self.libraries))
        print("Linked static libraries: %s" % " ".join(self.extra_objects))
        print("Extra compiler options: %s" % " ".join(self.extra_compile_args))
        print("Extra linker options: %s" % " ".join(self.extra_link_args))

    def replace_static_libraries(self):
        """Replaces references to libraries with full paths to their static
        versions if the static version is to be found on the library path."""
        if "-static" in self.extra_link_args:
            return

        self.extra_link_args.append("-static")
        for library_name in self.libraries[:]:
            static_lib = find_static_library(library_name, self.library_dirs)
            if static_lib:
                self.libraries.remove(library_name)
                self.extra_objects.append(static_lib)

    def use_educated_guess(self):
        """Tries to guess the proper library names, include and library paths
        if everything else failed."""

        preprocess_fallback_config()

        global LIBIGRAPH_FALLBACK_LIBRARIES
        global LIBIGRAPH_FALLBACK_INCLUDE_DIRS
        global LIBIGRAPH_FALLBACK_LIBRARY_DIRS

        print("WARNING: we were not able to detect where igraph is installed on")
        print("your machine (if it is installed at all). We will use the fallback")
        print("library and include pathss hardcoded in setup.py and hope that the")
        print("C core of igraph is installed there.")
        print("")
        print("If the compilation fails and you are sure that igraph is installed")
        print("on your machine, adjust the following two variables in setup.py")
        print("accordingly and try again:")
        print("")
        print("- LIBIGRAPH_FALLBACK_INCLUDE_DIRS")
        print("- LIBIGRAPH_FALLBACK_LIBRARY_DIRS")
        print("")

        seconds_remaining = 10
        while seconds_remaining > 0:
            if seconds_remaining > 1:
                plural = "s"
            else:
                plural = ""

            sys.stdout.write("\rContinuing in %2d second%s; press Enter to continue "
                    "immediately. " % (seconds_remaining, plural))
            sys.stdout.flush()

            rlist, _, _ = select([sys.stdin], [], [], 1)
            if rlist:
                sys.stdin.readline()
                break

            seconds_remaining -= 1
        sys.stdout.write("\r" + " "*65 + "\r")

        self.libraries = LIBIGRAPH_FALLBACK_LIBRARIES[:]
        if self.static_extension:
            self.libraries.extend(["xml2", "z", "m"])
        self.include_dirs = LIBIGRAPH_FALLBACK_INCLUDE_DIRS[:]
        self.library_dirs = LIBIGRAPH_FALLBACK_LIBRARY_DIRS[:]

buildcfg = BuildConfiguration()

###########################################################################

# Process command line options
preprocess_args()

# Define the extension
sources=glob.glob(os.path.join('src', '*.c'))
igraph_extension = Extension('igraph._igraph', sources)
        # library_dirs=library_dirs,
        # libraries=libraries,
        # include_dirs=include_dirs,
        # extra_objects=extra_objects,
        # extra_link_args=extra_link_args
       
description = """Python interface to the igraph high performance graph
library, primarily aimed at complex network research and analysis.

Graph plotting functionality is provided by the Cairo library, so make
sure you install the Python bindings of Cairo if you want to generate
publication-quality graph plots.

See the `Cairo homepage <http://cairographics.org/pycairo>`_ for details.

From release 0.5, the C core of the igraph library is **not** included
in the Python distribution - you must compile and install the C core
separately. Windows installers already contain a compiled igraph DLL,
so they should work out of the box. Linux users should refer to the
`igraph homepage <http://igraph.org>`_ for
compilation instructions (but check your distribution first, maybe
there are pre-compiled packages available). OS X users may
benefit from the disk images in the Python Package Index.

Unofficial installers for 64-bit Windows machines and/or different Python
versions can also be found `here <http://www.lfd.uci.edu/~gohlke/pythonlibs>`_.
Many thanks to the maintainers of this page!
"""

options = dict(
    name = 'python-igraph',
    version = VERSION,
    description = 'High performance graph data structures and algorithms',
    long_description = description,
    license = 'GNU General Public License (GPL)',

    author = 'Tamas Nepusz',
    author_email = 'tamas@cs.rhul.ac.uk',

    ext_modules = [igraph_extension],
    package_dir = {'igraph': 'igraph'},
    packages = ['igraph', 'igraph.test', 'igraph.app', 'igraph.drawing',
        'igraph.remote', 'igraph.vendor'],
    scripts = ['scripts/igraph'],
    test_suite = "igraph.test.suite",

    headers = ['src/igraphmodule_api.h'],

    platforms = 'ALL',
    keywords = ['graph', 'network', 'mathematics', 'math', 'graph theory', 'discrete mathematics'],
    classifiers = [
      'Development Status :: 4 - Beta',
      'Intended Audience :: Developers',
      'Intended Audience :: Science/Research',
      'Operating System :: OS Independent',
      'Programming Language :: C',
      'Programming Language :: Python',
      'Topic :: Scientific/Engineering',
      'Topic :: Scientific/Engineering :: Information Analysis',
      'Topic :: Scientific/Engineering :: Mathematics',
      'Topic :: Scientific/Engineering :: Physics',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
      'Topic :: Software Development :: Libraries :: Python Modules'
    ],

    cmdclass = {
        "build_ext": buildcfg.build_ext
    }
)

if "macosx" in get_platform() and "bdist_mpkg" in sys.argv:
    # OS X specific stuff to build the .mpkg installer
    options["data_files"] = [ \
            ('/usr/local/lib', [os.path.join('..', '..', 'fatbuild', 'libigraph.0.dylib')])
    ]

if sys.version_info > (3, 0):
    if build_py is None:
        options["use_2to3"] = True
    else:
        options["cmdclass"]["build_py"] = build_py

setup(**options)
