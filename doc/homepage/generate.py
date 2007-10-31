#!/usr/bin/env python
"""
Homepage generator for igraph.

Usage: ./generate.py version [newsfile]
"""

import sys
import string
import glob
import os.path
import re

def detect_newsfile():
    """Tries to detect the full path of the NEWS file"""
    dirs = [".", "..", os.path.join("..", "..")]
    for d in dirs:
	fname = os.path.join(d, "NEWS")
	if os.path.exists(fname): return fname
    raise ValueError, "NEWS file not found"

def process_newsfile(fname):
    """Processes the NEWS file and formats its content in HTML"""
    news = open(fname).read()
    
    # Titles are wrapped in <h3>
    reg = re.compile(r"^=====*\n(?P<title>[^\n]*)\n====*\n", re.MULTILINE | re.DOTALL)
    news = reg.sub("<h3>\g<title></h3>", news)
    
    # Lists are wrapped in <ul>
    reg = re.compile(r"^[\s]*\n(?P<lines>(- [^\n]\n)+)[\s]*\n", re.MULTILINE)
    news = reg.sub("<ul>\g<lines></ul>", news)
    
    # List items are wrapped in <li>
    reg = re.compile(r"^- (?P<entry>.*?)(?=(^$)|(^-[ ])|(\Z))", re.MULTILINE | re.DOTALL)
    news = reg.sub("<li>\g<entry></li>", news)
    
    # Paragraphs are separated by empty lines
    reg = re.compile(r"^[\s]*$", re.MULTILINE)
    news = reg.sub("<p></p>", news)
    
    # Replace URLs
    reg = re.compile(r"(?P<url>http://[^\n ]+)")
    news = reg.sub("<a href=\"\g<url>\">\g<url></a>", news)
    
    # Add subheadings in <h4> tags
    reg = re.compile(r"^(?P<text>.*)\n-------*$", re.MULTILINE)
    news = reg.sub("<h4>\g<text></h4>", news)
    
    return news

if len(sys.argv)<2:
    print __doc__
    sys.exit(1)

version = sys.argv[1]
if len(sys.argv)>= 3:
    newsfile = sys.argv[2]
else:
    newsfile = detect_newsfile()
    
template = [string.Template(line) for line in open("template.html")]
tokens = { "VERSION": version, "NEWS": process_newsfile(newsfile) }

files = glob.glob('*.html.in')

for file in files:
    f = open(file[:-3], "w")
    content = open(file).read()
    tokens["CONTENT"] = string.Template(content).safe_substitute(tokens)
    for tmpl in template:
	print >>f, tmpl.safe_substitute(tokens),
	
    f.close()
