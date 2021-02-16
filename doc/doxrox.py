#! /usr/bin/env python3

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

The generator parses one or more input files for documentation chunks
(embedded in the source code as Doxygen-style comments), and processes
them with a set of regex-based rules. The processed chunks are then
substituted into a template file containing <!-- doxrox-include -->
directives.

When a template file is not provided, the generator will read the input
files, process them with the ruleset and save a dictionary mapping chunk
names to the corresponding processed chunks into a Python pickle. This
can be used to speed up the processing of multiple input files as you can
generate the chunks once and then re-use them for multiple input files.
"""

import sys
import re

from argparse import ArgumentParser
from contextlib import contextmanager
from pickle import dump, load

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

    outputfile = arguments.output_file
    verbose = arguments.verbose
    cutit = arguments.cut_it
    args = arguments.inputs

    if (
        arguments.template_file in args
        or arguments.rules_file in args
        or outputfile in args
    ):
        print("Error, special file is also used as an input file")
        parser.print_help()
        sys.exit(2)

    # get all regular expressions
    if arguments.rules_file:
        with operation("Reading regular expressions...") as op:
            regexlist = readregex(arguments.rules_file)
            op("{0} rules read".format(len(regexlist)))
    else:
        regexlist = []

    # parse all input files and extract chunks, apply rules
    if arguments.chunk_file:
        with operation("Reading pickled chunks...") as op:
            try:
                with open(arguments.chunk_file, "rb") as f:
                    all_chunks = load(f)
            except IOError:
                print("Error reading chunk file: " + arguments.chunk_file)
                sys.exit(9)
            op("{0} chunks read".format(len(all_chunks)))
    else:
        all_chunks = dict()

    for ifile in args:
        with operation("Parsing input file {0}...".format(ifile)) as op:
            try:
                with open(ifile, "r") as f:
                    strinput = f.read()
            except IOError:
                print("Error reading input file: " + ifile)
                sys.exit(3)

            num_new_chunks = collect_chunks_from_input_file(
                strinput, regexlist, all_chunks
            )
            op("{0} chunks read".format(num_new_chunks))

    if arguments.template_file:
        # substitute the template file
        with operation("Reading template file..."):
            try:
                with open(arguments.template_file, "r") as tfile:
                    tstring = tfile.read()
            except IOError:
                print("Error reading the template file: " + arguments.template_file)
                sys.exit(7)

        with operation("Substituting template file..."):
            chunk_iterator = re.finditer(
                r"<!--\s*doxrox-include\s+(\w+)\s+-->", tstring
            )
            outstring = []
            last = 0
            for match in chunk_iterator:
                try:
                    chunk = all_chunks[match.group(1)]
                except KeyError:
                    print("Chunk not found: {0}".format(match.group(1)))
                    sys.exit(4)
                outstring.append(tstring[last : match.start()])
                outstring.append(chunk)
                last = match.end()
            outstring.append(tstring[last:])
            outstring = "".join(outstring)

        # write output file
        with operation("Writing output file..."):
            try:
                with open(outputfile, "w") as ofile:
                    ofile.write(outstring)
            except IOError:
                print("Error writing output file:" + outputfile)
                sys.exit(8)
    else:
        # no template file given so just save the chunks as a pickle into the
        # output file
        with operation("Writing output file..."):
            try:
                with open(outputfile, "wb") as ofile:
                    dump(all_chunks, ofile)
            except IOError:
                print("Error writing output file:" + outputfile)
                sys.exit(5)


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
        help="template file to process",
    )
    parser.add_argument(
        "-e",
        "--rules",
        metavar="FILE",
        dest="rules_file",
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
        "--chunks",
        dest="chunk_file",
        metavar="FILE",
        help="name of a previously saved chunk file",
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
def collect_chunks_from_input_file(strinput, regexlist, all_chunks):
    global cutit

    num_new_chunks = 0

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
            # print("Chunk without a name ignored:" + ch[0:60] + "...")
            continue

        if name in all_chunks:
            print("Multiple defined name: " + name)
            sys.exit(6)

        all_chunks[name] = actch.strip()
        num_new_chunks += 1
    return num_new_chunks


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
