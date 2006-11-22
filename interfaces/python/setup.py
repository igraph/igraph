#!/usr/bin/env python

from distutils.core import setup, Extension
from distutils.file_util import copy_file
from sys import version, exit
import glob
import os.path
from os import popen, mkdir
from shutil import copy2

LIBXML2_FALLBACK_INCLUDE_DIRS = ['/usr/include/libxml2']
LIBXML2_FALLBACK_LIBRARIES = ['xml2', 'iconv']
LIBXML2_FALLBACK_LIBRARY_DIRS = []

if version < '2.3':
    print "This module requires Python >= 2.3"
    exit(0)
    
def get_output(command):
    """Returns the output of a command returning a single line of output"""
    f=popen(command)
    line=f.readline().strip()
    exit_code=f.close()
    return line, exit_code
    
def detect_libxml2_include_dirs(default = LIBXML2_FALLBACK_INCLUDE_DIRS):
    """Tries to detect the libxml2 include directory"""
    line, exit_code = get_output("xml2-config --cflags")
    if exit_code>0: return default
    opts=line.split()
    return [opt[2:] for opt in opts if opt[0:2]=="-I"]

def detect_libxml2_libraries(default = LIBXML2_FALLBACK_LIBRARIES):
    """Tries to detect the libraries that libxml2 uses"""
    line, exit_code = get_output("xml2-config --libs")
    if exit_code>0: return default
    opts=line.split()
    return [opt[2:] for opt in opts if opt[0:2]=="-l"]
    
def detect_libxml2_library_dirs(default = LIBXML2_FALLBACK_LIBRARY_DIRS):
    """Tries to detect the libxml2 library directory"""
    line, exit_code = get_output("xml2-config --libs")
    if exit_code>0: return default
    opts=line.split()
    return [opt[2:] for opt in opts if opt[0:2]=="-L"]

def detect_igraph_source():
    """Tries to detect the igraph sources and copy it to igraph/ if necessary"""
    if not os.path.isdir('igraph'): os.mkdir('igraph')
    if os.path.isfile(os.path.join('..', '..', 'src', 'igraph.h')):
	files_to_copy = ['*.c', '*.h', '*.y']
	src_files = [os.path.join('..', '..', 'config.h')]
	for wildcard in files_to_copy:
	    src_files.extend(glob.glob(os.path.join('..', '..', 'src', wildcard)))
	for src_file in src_files:
	    copy_file(src_file, 'igraph', update=1)
	
	return
	
def igraph_version():
    """Returns igraph version number"""
    try:
	f=open(os.path.join('igraph', 'config.h'))
	while True:
	    line=f.readline()
	    if not line: break
	    line=line.strip()
	    if line[0:16] == "#define VERSION ":
		return line[16:].replace('"', '')
	    
	f.close()
	raise Exception, "Can't detect igraph version number from source code!"
    except:
	raise Exception, "igraph source code not found"
    
try:
    detect_igraph_source()
except:
    print "An error happened while trying to find igraph source!"
    print "Consider downloading the C source code of the igraph library and "
    print "put the contents of the src subdirectory and config.h in "
    print "subdirectory called igraph."
    exit(1)

module_sources=glob.glob(os.path.join('src', '*.c'))
sources=glob.glob(os.path.join('igraph', '*.c'))

sources.extend(module_sources)

include_dirs=['igraph', '.']
library_dirs=[]
libraries=[]

line, exit_code = get_output("xml2-config --version")
if exit_code>0:
    print "An error happened while trying to get libxml2 compilation options!"
    print "Maybe you don't have libxml2 installed?"
    print "Falling back to reasonable defaults."
    print "If the compilation fails, please edit the LIBXML2_FALLBACK_* variables"
    print "in setup.py to point to the correct directories and libraries."
    print

include_dirs.extend(detect_libxml2_include_dirs())
library_dirs.extend(detect_libxml2_library_dirs())
libraries.extend(detect_libxml2_libraries())

igraph_extension = Extension('igraph._igraph', sources, \
  library_dirs=library_dirs, libraries=libraries, \
  include_dirs=include_dirs)
       
description = """This module extends Python with a Graph class which is capable
of handling arbitrary directed and undirected graphs with thousands of nodes and
millions of edges. Since the module makes use of the open source igraph library
written in almost 100% pure C, it is blazing fast and outperforms most other
pure Python-based graph packages around."""

setup(name = 'igraph',
      version = igraph_version(),
      description = 'High performance graph data structures and algorithms',
      long_description = description,
      license = 'GNU General Public License (GPL)',
      author = 'Tamas Nepusz',
      author_email = 'ntamas@rmki.kfki.hu',
      url = 'http://cneurocvs.rmki.kfki.hu/igraph/',
      download_url = 'http://cneurocvs.rmki.kfki.hu/igraph/',
      ext_modules = [igraph_extension],
      package_dir = {'igraph': 'package'},
      packages = ['igraph', 'igraph.test'],
      platforms = 'ALL',
      keywords = ['graph', 'network', 'mathematics', 'math', 'graph theory', 'discrete mathematics'],
      classifiers = [
        'Development Status :: 4 - Beta',
	'Intended Audience :: Developers',
	'Intended Audience :: Science/Research',
	'License :: OSI Approved :: GNU General Public License',
	'Operating System :: OS Independent',
	'Programming Language :: C',
	'Programming Language :: Python',
	'Topic :: Scientific/Engineering',
	'Topic :: Scientific/Engineering :: Information Analysis',
	'Topic :: Scientific/Engineering :: Mathematics',
	'Topic :: Scientific/Engineering :: Physics',
	'Topic :: Scientific/Engineering :: Bio-Informatics',
	'Topic :: Software Development :: Libraries :: Python Modules'
      ])
