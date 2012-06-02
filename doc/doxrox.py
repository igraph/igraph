#! /usr/bin/env python

#   IGraph R package
#   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
#   334 Harvard street, Cambridge, MA 02139 USA
#   
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
#   02110-1301 USA
#
###################################################################

import sys
import getopt
import re
import string

#################
# constants, these might turn to parameters some time
#################
doxhead='\/\*\*'

#################
# global variables
#################
verbose=False
cutit=False

#########################################################################
# The main function
#########################################################################

def main():

    global verbose, cutit
    
    # get command line arguments
    try:
        optlist, args = getopt.getopt(sys.argv[1:], 't:e:o:hvc', ['help'])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)

    # handle command line arguments
    templatefile=regexfile=outputfile=""
    verbose=False
    
    for o, a in optlist:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o == "-t":
            templatefile = a
        if o == "-e":
            regexfile = a
        if o == "-o":
            outputfile = a
        if o == "-v":
            verbose = True
        if o == "-c":
            cutit = True
        
    if templatefile == "" or regexfile == "" or outputfile == "":
        print("Error, some special file is not given")
        usage()
        sys.exit(2)

    if templatefile in args or regexfile in args or outputfile in args:
        print("Error, special file is also used as an input file")
        usage()
        sys.exit(2)

    if templatefile == regexfile or templatefile == outputfile or \
        regexfile == outputfile:
        print('Error, some special files are the same')
        usage()
        sys.exit(2)

    # get all regular expressions
    if verbose:
        print 'Reading regular expressions...',
    regexlist=readregex(regexfile)
    if verbose:
        print("done, "+str(len(regexlist))+" rules read.")

    # parse all input files and extract chunks, apply rules
    docchunks=dict()
    for ifile in args:
        if verbose:
            print 'Parsing input file '+ifile+'...',
        try:
            f=open(ifile, 'r')
            strinput=f.read()
            f.close()
        except IOError:
            print("Error reading input file: "+ifile)
            sys.exit(3)
        parsestring(strinput, regexlist, docchunks)
        if verbose:
            print('done, '+str(len(docchunks))+" chunks read.")

    # substitute the template file
    try:
        if verbose:
            print "Reading template file...",
        tfile=open(templatefile, 'r')
        tstring=tfile.read()
        tfile.close()
        if verbose:
            print('done.')
    except IOError:
        print("Error reading the template file: "+templatefile)
        sys.exit(7)
    if verbose:
        print "Substituting template file...",        
    chunkit=re.finditer(r'<!--\s*doxrox-include\s+(\w+)\s+-->', tstring)
    outstring=""
    last=0
    for chunk in chunkit:
        outstring=outstring+tstring[last:chunk.start()]+\
                   docchunks[chunk.group(1)]
        last=chunk.end()
    outstring=outstring+tstring[last:]
    if verbose:
        print "done."

    # write output file
    try:
        if verbose:
            print "Writing output file...",
        ofile=open(outputfile, 'w')
        ofile.write(outstring)
        ofile.close()
    except IOError:
        print("Error writing output file:"+outputfile)
        sysexit(8)
    if verbose:
        print "done."
    
#########################################################################
# End of the main function
#########################################################################

#################
# read the regular expressions
#################
def readregex(regexfile):
    lines=[]
    mode="empty"
    actreplace=""
    actwith=""
    acttype=""
    lineno=1
    try:
        f=open(regexfile, "r")
        for line in f:
            # a new pattern block starts
            if line[0:7]=="REPLACE":
                if mode not in ("empty","with"):
                    print("Parse error in regex file ("+regexfile+"), line "+
                          str(lineno))
                    sys.exit(4)
                else:
                    if (actreplace != ""):
                        readregexappend(lines, actreplace, actwith, acttype)
                    actreplace=actwith=""
                    mode="replace"
            # the second half of the pattern block starts
            elif line[0:4]=="WITH" or line[0:3]=="RUN":
                if mode != "replace":
                    print("Parse error in regex file ("+regexfile+"), line "+
                          str(lineno))
                    sys.exit(4)
                else:
                    mode="with"
                if line[0:4]=="WITH":
                    acttype="with"
                else:
                    acttype="run"
            # empty line, do nothing
            elif re.match("^\s*$", line):
                1==1
            # normal line, append
            else:
                if mode=="replace":
                    actreplace=actreplace+line
                elif mode=="with":
                    actwith=actwith+line
                else:
                    print("Parse error in regex file ("+regexfile+"), line "+
                          str(lineno))
                    sys.exit(4)
            lineno=lineno+1

        if actreplace != "":
            readregexappend(lines, actreplace, actwith, acttype)
        f.close()
    except IOError:
        print("Error reading regex file: "+regexfile)
        sys.exit(4)
    return (lines)
        
def readregexappend(lines, actreplace, actwith, acttype):
    compactreplace=re.compile(actreplace,re.VERBOSE|re.MULTILINE|re.DOTALL)
    actwith=actwith[:(len(actwith)-1)]
    lines.append( (compactreplace, actwith, acttype) )
        
#################
# parse an input file string
#################
def parsestring(strinput, regexlist, docchunks):
    global cutit
    # split the file
    chunks=re.split(doxhead, strinput)
    chunks=chunks[1:]
    # apply all rules to the chunks
    for ch in chunks:
        if cutit:
            ch=ch.split("/*")[0]
        actch=ch
        name=""
        for reg in regexlist:
            matched=reg[0].match(actch)
            if name=="" and matched != None:
                try:
                    name=matched.group('name')
                except IndexError:
                    name=""
            if reg[2]=="with":
                try:
                    actch=reg[0].sub(reg[1], actch)
                except IndexError:
                    print("Index error:"+ch[0:60]+"...")
                    print("Pattern:\n"+reg[0].pattern)
                    print("Current state:"+actch[0:60]+"...")
                    sys.exit(6)
            elif reg[2]=="run":
                exec(reg[1])
        if name=="":
            print("Chunk without a name ignored:"+ch[0:60]+"...")
            continue
        if docchunks.has_key(name):
            print("Multiple defined name: "+name)
            sys.exit(6)
        if verbose:
            print name,
        docchunks[name]=string.strip(actch)
    return(docchunks)

#################
# print out some help
#################
def usage():
    print("Usage: " + sys.argv[0] + " [-vh] -t template-file -e regex-file -o output-file\n"+
          "                   [--help] [inputfile]...")


if __name__ == "__main__":
    main()
