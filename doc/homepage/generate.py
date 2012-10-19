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
import HTMLParser

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
    
    # Add some <code> tags
    reg = re.compile(r"(?P<func>\w*_\w*)")
    news = reg.sub("<code>\g<func></code>", news)

    # Titles are wrapped in <h3>
    reg = re.compile(r"^=====*\n(?P<title>[^\n]*)\n====*\n", re.MULTILINE | re.DOTALL)
    news = reg.sub("</p><h3>\g<title></h3><p class=\"news\">", news)
    
    # Lists are wrapped in <ol>
    reg = re.compile(r"\n(?P<lines>(^\s*\n|^(- |  )[^\n]*\n)+)[\s]*\n(?P<n>[^-])", re.MULTILINE)
    news = reg.sub("\n</p><ul class=\"newslist\">\n\g<lines></ul>\n<p class=\"news\">\n\g<n>", news)
    
    # List items are wrapped in <li>
    reg = re.compile(r"^- (?P<entry>.*?)(?=(</ul>)|(^$)|(^-[ ])|(\Z))", re.MULTILINE | re.DOTALL)
    news = reg.sub("<li>\g<entry></li>", news)
    
    # Paragraphs are separated by empty lines
    reg = re.compile(r"^[\s]*$", re.MULTILINE)
    news = reg.sub("</p><p class=\"news\">", news)
    news = "<p>" + news + "</p>"
    
    # Replace URLs
    reg = re.compile(r"(?P<url>https?://[^\n ]+)")
    news = reg.sub("<a href=\"\g<url>\">\g<url></a>", news)
    
    # Add subheadings in <h4> tags
    reg = re.compile(r"^(?P<text>.*)\n-------*$", re.MULTILINE)
    news = reg.sub("<h4 class=\"news\">\g<text></h4>", news)
    
    return news

def addExamples(text):
    ex=re.search("<example[^>]*>", text)
    while ex:
        exfile=re.search(' file="([^"]*)"', ex.group()).groups()[0]
        extext=open("examples/" + exfile).read()
        text = text[0:ex.start()] + extext + text[ex.end():]
        ex=re.search("<example[^>]*>", text)
    return text

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
    content = addExamples(content)
    tokens["CONTENT"] = string.Template(content).safe_substitute(tokens)
    tokens["PAGENAME"] = file[:-8]
    for tmpl in template:
	print >>f, tmpl.safe_substitute(tokens),
	
    f.close()
