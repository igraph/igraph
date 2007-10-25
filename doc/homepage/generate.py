#!/usr/bin/env python
"""
Homepage generator for igraph.

Usage: ./generate.py version
"""

import sys
import string
import glob

if len(sys.argv)<2:
    print __doc__
    sys.exit(1)

version = sys.argv[1]

template = [string.Template(line) for line in open("template.html")]
tokens = { "VERSION": version }

files = glob.glob('*.html.in')

for file in files:
    f = open(file[:-3], "w")
    content = open(file).read()
    tokens["CONTENT"] = string.Template(content).safe_substitute(tokens)
    for tmpl in template:
	print >>f, tmpl.safe_substitute(tokens),
	
    f.close()
