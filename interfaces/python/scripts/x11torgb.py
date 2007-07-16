#!/usr/bin/env python
"""Converts X11 color files to RGB formatted Python dicts"""
import sys
import pprint

if len(sys.argv)<2:
    print "Usage: %s filename" % sys.argv[0]
    sys.exit(1)

colors = {
  "black":    (0.  , 0.  , 0.  ),
  "silver":   (0.75, 0.75, 0.75),
  "gray":     (0.5 , 0.5 , 0.5 ),
  "white":    (1.  , 1.  , 1.  ),
  "maroon":   (0.5 , 0.  , 0.  ),
  "red":      (1.  , 0.  , 0.  ),
  "purple":   (0.5 , 0.  , 0.5 ),
  "fuchsia":  (1.  , 0.  , 1.  ),
  "green":    (0.  , 0.5 , 0.  ),
  "lime":     (0.  , 1.  , 0.  ),
  "olive":    (0.5 , 0.5 , 0.  ),
  "yellow":   (1.  , 1.  , 0.  ),
  "navy":     (0.  , 0.  , 0.5 ),
  "blue":     (0.  , 0.  , 1.  ),
  "teal":     (0.  , 0.5 , 0.5 ),
  "aqua":     (0.  , 1.  , 1.  ),
}

f = open(sys.argv[1])
for line in f:
    if line[0] == '!': continue
    parts = line.strip().split(None, 3)
    for x in xrange(3):
	parts[x] = float(parts[x])/255.
    colors[parts[3].lower()] = tuple(parts[0:3])

pp = pprint.PrettyPrinter(indent=4)
pp.pprint(colors)

