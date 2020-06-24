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

"""DocBook XML generator for igraph.

"""

import sys
import re

from argparse import ArgumentParser
from contextlib import contextmanager

#################
# constants, these might turn to parameters some time
#################
doxhead = r"/\*\*"

#################
# global variables
#################
verbose = False
cutit = False

#########################################################################
# The main function
#########################################################################


def main():

    global verbose, cutit

    # get command line arguments
    parser = create_argument_parser()
    arguments = parser.parse_args()

    templatefile = arguments.template_file
    regexfile = arguments.rules_file
    outputfile = arguments.output_file
    verbose = arguments.verbose
    cutit = arguments.cut_it
    args = arguments.inputs

    if templatefile in args or regexfile in args or outputfile in args:
        print("Error, special file is also used as an input file")
        parser.print_help()
        sys.exit(2)

    if (
        templatefile == regexfile
        or templatefile == outputfile
        or regexfile == outputfile
    ):
        print("Error, some special files are the same")
        parser.print_help()
        sys.exit(2)

    # get all regular expressions
    with operation("Reading regular expressions...") as op:
        regexlist = readregex(regexfile)
        op("{0} rules read".format(len(regexlist)))

    # parse all input files and extract chunks, apply rules
    docchunks = dict()
    for ifile in args:
        with operation("Parsing input file {0}...".format(ifile)) as op:
            try:
                with open(ifile, "r") as f:
                    strinput = f.read()
            except IOError:
                print("Error reading input file: " + ifile)
                sys.exit(3)
            parsestring(strinput, regexlist, docchunks)
            op("{0} chunks read".format(len(docchunks)))

    # substitute the template file
    with operation("Reading template file..."):
        try:
            with open(templatefile, "r") as tfile:
                tstring = tfile.read()
        except IOError:
            print("Error reading the template file: " + templatefile)
            sys.exit(7)

    with operation("Substituting template file..."):
        chunkit = re.finditer(r"<!--\s*doxrox-include\s+(\w+)\s+-->", tstring)
        outstring = ""
        last = 0
        for chunk in chunkit:
            outstring = (
                outstring + tstring[last : chunk.start()] + docchunks[chunk.group(1)]
            )
            last = chunk.end()
        outstring = outstring + tstring[last:]

    # write output file
    with operation("Writing output file..."):
        try:
            with open(outputfile, "w") as ofile:
                ofile.write(outstring)
        except IOError:
            print("Error writing output file:" + outputfile)
            sys.exit(8)


#########################################################################
# Argument parser
#########################################################################


def create_argument_parser():
    parser = ArgumentParser(description=sys.modules[__name__].__doc__.strip())

    parser.add_argument(
        "-t",
        "--template",
        metavar="FILE",
        dest="template_file",
        required=True,
        help="template file to process",
    )
    parser.add_argument(
        "-e",
        "--rules",
        metavar="FILE",
        dest="rules_file",
        required=True,
        help="file containing matching and replacement rules",
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="FILE",
        dest="output_file",
        required=True,
        help="name of the output file",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        dest="verbose",
        help="enable verbose output",
    )
    parser.add_argument(
        "-c", "--cut", action="store_true", default=False, dest="cut_it"
    )

    parser.add_argument(
        "inputs", metavar="INPUT", nargs="*", help="input files to process"
    )

    return parser


#################
# read the regular expressions
#################
def readregex(regexfile):
    lines = []
    mode = "empty"
    actreplace = ""
    actwith = ""
    acttype = ""
    lineno = 1
    try:
        with open(regexfile, "r") as f:
            for line in f:
                # a new pattern block starts
                if line[0:7] == "REPLACE":
                    if mode not in ("empty", "with"):
                        print(
                            "Parse error in regex file ("
                            + regexfile
                            + "), line "
                            + str(lineno)
                        )
                        sys.exit(4)
                    else:
                        if actreplace != "":
                            readregexappend(lines, actreplace, actwith, acttype)
                        actreplace = actwith = ""
                        mode = "replace"
                # the second half of the pattern block starts
                elif line[0:4] == "WITH" or line[0:3] == "RUN":
                    if mode != "replace":
                        print(
                            "Parse error in regex file ("
                            + regexfile
                            + "), line "
                            + str(lineno)
                        )
                        sys.exit(4)
                    else:
                        mode = "with"
                    if line[0:4] == "WITH":
                        acttype = "with"
                    else:
                        acttype = "run"
                # empty line, do nothing
                elif re.match(r"^\s*$", line):
                    pass
                # normal line, append
                else:
                    if mode == "replace":
                        actreplace = actreplace + line
                    elif mode == "with":
                        actwith = actwith + line
                    else:
                        print(
                            "Parse error in regex file ("
                            + regexfile
                            + "), line "
                            + str(lineno)
                        )
                        sys.exit(4)
                lineno = lineno + 1

            if actreplace != "":
                readregexappend(lines, actreplace, actwith, acttype)
    except IOError:
        print("Error reading regex file: " + regexfile)
        sys.exit(4)
    return lines


def readregexappend(lines, actreplace, actwith, acttype):
    compactreplace = re.compile(actreplace, re.VERBOSE | re.MULTILINE | re.DOTALL)
    actwith = actwith[: (len(actwith) - 1)]
    lines.append((compactreplace, actwith, acttype))


#################
# parse an input file string
#################
def parsestring(strinput, regexlist, docchunks):
    global cutit
    # split the file
    chunks = re.split(doxhead, strinput)
    chunks = chunks[1:]
    # apply all rules to the chunks
    for ch in chunks:
        if cutit:
            ch = ch.split("/*")[0]
        actch = ch
        name = ""
        for reg in regexlist:
            matched = reg[0].match(actch)
            if name == "" and matched is not None:
                try:
                    name = matched.group("name")
                except IndexError:
                    name = ""
            if reg[2] == "with":
                try:
                    actch = reg[0].sub(reg[1], actch)
                except IndexError:
                    print("Index error:" + ch[0:60] + "...")
                    print("Pattern:\n" + reg[0].pattern)
                    print("Current state:" + actch[0:60] + "...")
                    sys.exit(6)
            elif reg[2] == "run":
                exec(reg[1])
        if name == "":
            print("Chunk without a name ignored:" + ch[0:60] + "...")
            continue
        if name in docchunks:
            print("Multiple defined name: " + name)
            sys.exit(6)
        docchunks[name] = actch.strip()
    return docchunks


@contextmanager
def operation(message):
    """Helper function to show progress messages for a potentially long-running
    operation in verbose mode.

    Parameters:
        message (str): the message to show
    """
    global verbose

    if verbose:
        print(message, end="")

    result = [None]

    def set_result(obj):
        result[0] = obj

    success = False
    try:
        yield set_result
        success = True
    finally:
        if verbose and success:
            if result[0] is None:
                print(" done.")
            else:
                print(" done, {0}.".format(result[0]))


if __name__ == "__main__":
    main()
