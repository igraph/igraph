#! /usr/bin/env python

import sys
import os.path
import re

from subprocess import check_output

# Some notes:
# - we have some sources with .cc extensions, these are marked as type 0
# - we have some non-standard header files, e.g. .hh, .pmt, these are type 2

srctext = """                       <File
                               RelativePath=".\src\%s"
                               FileType="0">
                       </File>"""

headtext = """                       <File
                               RelativePath=".\include\%s"
                               FileType="2">
                       </File>"""
headptext = """                       <File
                               RelativePath=".\src\%s"
                               FileType="2">
                       </File>"""

def runmake(makefile, target):
    out = check_output("make -q -C " + os.path.dirname(makefile) + 
                       " -f " + os.path.basename(makefile) + 
                       " " + target, shell=True)
    out = out.decode(sys.stdout.encoding or "utf8")
    out = out.replace("/", "\\")
    out = re.sub("make.*'.*'", "", out) # msys2 make adds "make[x]: entering '<dir>'""
    return out

def rreplace(s, old, new, occurrence):
    li = s.rsplit(old, occurrence)
    return new.join(li)

def main():
    if len(sys.argv) != 4:
        print( "Error: need three arguments")
        sys.exit(1)
    package = sys.argv[1]
    projectfile = sys.argv[2]
    makefile = sys.argv[3]

    proj = open(projectfile).read()
    sources = runmake(makefile, "echosources").split()
    headers = runmake(makefile, "echoheaders").split()
    headersprivate = runmake(makefile, "echoheadersprivate").split()

    # lex and bison stuff
    headers2 = [ rreplace(s, ".y", ".h",1) for s in sources if s[-2:]==".y" ]
    sources = [ rreplace(s, ".l", ".c", 1) for s in sources ]
    sources = [ rreplace(s, ".y", ".c", 1) for s in sources ]

    stext = "\n".join([ srctext % s for s in sources ])
    htext = "\n".join([ headtext % s for s in headers ])
    hptext = "\n".join([ headptext % s for s in headersprivate + headers2 ])

    proj = proj.replace("<!-- SOURCE-FILES -->", stext)
    proj = proj.replace("<!-- HEADER-FILES -->", htext + "\n" + hptext)

    out_file = open(package + "/igraph.vcproj", "w")
    out_file.write(proj)
    out_file.close()
    
if __name__ == "__main__":
    main()

