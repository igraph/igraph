#!/usr/bin/env python

import distutils.util
import os.path
import sys

build_dir = os.path.join('build', 'lib.'+distutils.util.get_platform()+'-'+sys.version[0:3])

sys.path.insert(0, build_dir)

import igraph.test
igraph.test.test()
