#! /usr/bin/env python3

#   IGraph library
#   Copyright (C) 2005-2021  The igraph development team
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
from collections import defaultdict, namedtuple
from contextlib import contextmanager
from hashlib import sha1
from operator import itemgetter
from pathlib import Path
from pickle import dump, load
from time import time

#: Constant indicating the start of a comment that doxrox.py will process
DOXHEAD = r"/\*\*"

#: Stores whether we want verbose output
verbose = False


def fatal(message, code = 1):
    """Prints an error message and exits the program with the given error code."""
    print(message, file=sys.stderr)
    sys.exit(code)


#########################################################################
# The main function
#########################################################################


def main():
    """Main entry point of the script."""

    global verbose

    # get command line arguments
    parser = create_argument_parser()
    arguments = parser.parse_args()

    outputfile = arguments.output_file
    verbose = arguments.verbose
    inputs = arguments.inputs

    if (
        arguments.template_file in inputs
        or arguments.rules_file in inputs
        or outputfile in inputs
    ):
        fatal("Special file is also used as an input file", 2)

    # open the cache file if needed
    cache = ChunkCache(arguments.cache_file) if arguments.cache_file else None

    # get all regular expressions
    if arguments.rules_file:
        with operation("Reading regular expressions...") as op:
            rules = read_regex_rules_file(arguments.rules_file)
            op("{0} rules read".format(len(rules)))
    else:
        rules = []

    # parse all input files and extract chunks, apply rules
    if arguments.chunk_file:
        with operation("Reading pickled chunks...") as op:
            try:
                with open(arguments.chunk_file, "rb") as f:
                    all_chunks = load(f)
            except IOError:
                fatal("Error reading chunk file: " + arguments.chunk_file, 9)
            op("{0} chunks read".format(len(all_chunks)))
    else:
        all_chunks = {}

    rule_timings = defaultdict(list)

    for ifile in inputs:
        with operation("Parsing input file {0}...".format(ifile)) as op:
            try:
                with open(ifile, "r") as f:
                    contents = f.read()
            except IOError:
                fatal("Error reading input file: " + ifile, 3)

            if cache:
                key = cache.key_of(contents)
                chunks = cache.get(key)
            else:
                key, chunks = None, None

            if chunks is not None:
                op("{0} chunks read from cache".format(len(chunks)))
            else:
                chunks = collect_chunks_from_input_file(contents, rules, rule_timings)
                op("{0} chunks parsed".format(len(chunks)))
                if key:
                    cache.put(key, chunks)

            for name, chunk in chunks.items():
                if name in all_chunks:
                    fatal("Multiple files provide chunks for {0!r}".format(name), code=4)
                all_chunks[name] = chunk

    if arguments.timing_stats and rule_timings:
        rule_timings = {name: sum(dts) / len(dts) for name, dts in rule_timings.items()}
        for name, dt in sorted(rule_timings.items(), key=itemgetter(1), reverse=True):
            print("{0}: {1:.3f}us".format(name, dt))
        print("======")

    if cache:
        cache.close()

    if arguments.template_file:
        # substitute the template file
        with operation("Reading template file..."):
            try:
                with open(arguments.template_file, "r") as tfile:
                    tstring = tfile.read()
            except IOError:
                fatal("Error reading the template file: " + arguments.template_file, 7)

        with operation("Substituting template file..."):
            chunk_iterator = re.finditer(
                r"<!--\s*doxrox-include\s+(\w+)\s*-->", tstring
            )
            outstring = []
            last = 0
            for match in chunk_iterator:
                try:
                    chunk = all_chunks[match.group(1)]
                except KeyError:
                    fatal("Chunk not found: {0}".format(match.group(1)), code=4)
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
                fatal("Error writing output file:" + outputfile, 8)
    else:
        # no template file given so just save the chunks as a pickle into the
        # output file
        with operation("Writing output file..."):
            try:
                with open(outputfile, "wb") as ofile:
                    dump(all_chunks, ofile)
            except IOError:
                fatal("Error writing output file:" + outputfile, 5)


#########################################################################
# Argument parser
#########################################################################


def create_argument_parser():
    """Creates the command line argument parser that the script uses."""
    parser = ArgumentParser(description=sys.modules[__name__].__doc__.strip())

    parser.add_argument(
        "--cache", metavar="FILE",
        dest="cache_file",
        help="optional cache file to store chunks from already processed files"
    )
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
        "--chunks",
        dest="chunk_file",
        metavar="FILE",
        help="name of a previously saved chunk file",
    )
    parser.add_argument(
        "--timing-stats",
        dest="timing_stats",
        action="store_true",
        default=False,
        help="print the average time it takes to process regex rules from the rules file"
    )
    parser.add_argument(
        "inputs", metavar="INPUT", nargs="*", help="input files to process"
    )

    return parser


#################
# classes and functions to read the regular expression rules
#################


Rule = namedtuple("Rule", "regex replacement name type")


def read_regex_rules_file(filename):
    """Parses the file containing the regex-based rules that we use to chop
    up the input source files into chunks that can later be fed into a
    DocBook document.

    Parameters:
        filename: name of the input file

    Returns:
        the rules that were parsed from the input file
    """

    rules = []

    def parse_error(lineno):
        """Helper function to indicate a parse error at the given line."""
        fatal(
            "Parse error in regex file ({0}), line {1}".format(regexfile, lineno),
            code = 4
        )

    def store(rule, replacement, rule_name, rule_type):
        """Helper function to append the current rule to the result."""
        regex = re.compile("".join(rule), re.VERBOSE | re.MULTILINE | re.DOTALL)
        replacement = "".join(replacement)[:-1]
        rules.append(Rule(regex, replacement, rule_name, rule_type))

    mode = "empty"
    regex, replacement = [], []
    rule_name, rule_type = None, ""

    try:
        with open(filename, "r") as f:
            for lineno, line in enumerate(f, 1):
                if line.startswith("REPLACE"):
                    # a new pattern block starts
                    if mode not in ("empty", "with"):
                        parse_error(lineno)
                    else:
                        if regex:
                            store(regex, replacement, rule_name, rule_type)
                        regex.clear()
                        replacement.clear()
                        mode = "replace"

                        match = re.match(r"^REPLACE\s+-+\s+(?P<name>.*)\s+-", line)
                        rule_name = match.group("name") if match else None

                elif line.startswith("WITH") or line.startswith("RUN"):
                    # the second half of the pattern block starts
                    if mode != "replace":
                        parse_error(lineno)
                    else:
                        mode = "with"
                    rule_type = "with" if line.startswith("WITH") else "run"

                elif re.match(r"^\s*$", line):
                    # empty line, do nothing
                    pass

                else:
                    # normal line, append
                    if mode == "replace":
                        regex.append(line)
                    elif mode == "with":
                        replacement.append(line)
                    else:
                        parse_error(lineno)

            if regex != "":
                store(regex, replacement, rule_name, rule_type)

    except IOError:
        fatal("Error reading regex file: " + regexfile, code = 4)

    return rules


#################
# parse an input file string
#################
def collect_chunks_from_input_file(strinput, rules, rule_timings):
    result = {}

    # split the file
    chunks = re.split(DOXHEAD, strinput)
    chunks = chunks[1:]

    # apply all rules to the chunks
    for chunk in chunks:
        name = None

        for index, rule in enumerate(rules):
            start = time()

            if not name and "name" in rule.regex.groupindex:
                # The regex might provide us with a chunk name so try figuring
                # out what the "name" group might match to
                matched = rule.regex.search(chunk)
                if matched:
                    try:
                        name = matched.group("name")
                    except IndexError:
                        name = ""

            if rule.type == "with":
                # This is a simple regex replacement rule
                try:
                    chunk = rule.regex.sub(rule.replacement, chunk)
                except IndexError:
                    print("Index error:" + ch[0:60] + "...")
                    print("Pattern:\n" + rule.regex.pattern)
                    print("Current state:" + chunk[0:60] + "...")
                    fatal(code=6)
            elif rule.type == "run":
                # This is a piece of Python code that has to be executed on
                # the part that matched
                matched = rule.regex.search(chunk)
                if matched:
                    exec(rule.replacement)
            else:
                fatal("Invalid rule type: {0!r}".format(rule.type), code=6)

            rule_timings[rule.name].append((time() - start) * 1000000)

        if not name:
            # print("Chunk without a name ignored:" + ch[0:60] + "...")
            continue

        result[name] = chunk.strip()

    return result


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


class ChunkCache:
    """Simple on-disk cache that stores SHA256 hashes of files along with the
    DocBook documentation chunks that were parsed from them.
    """

    def __init__(self, filename, hash=sha1):
        """Constructor.

        Parameters:
            filename: name of the file on the disk where the cache resides
            hash: the hash function to use
        """
        self._data = None
        self._dirty = False
        self._hash = hash
        self._path = Path(filename)

    def _load(self):
        """Populates the in-memory copy of the cache from the disk."""
        if self._path.exists():
            try:
                with self._path.open("rb") as fp:
                    self._data = load(fp)
            except (IOError, EOFError):
                # cache corrupted
                self._data = {}
        else:
            self._data = {}
        self._dirty = False

    def close(self):
        """Closes the cache and flushes its contents back to the disk if it
        changed recently.
        """
        if self._dirty:
            self.flush()

    def flush(self):
        """Flushes the contents of the cache back to the disk."""
        with self._path.open("wb") as fp:
            dump(self._data, fp)
        self._dirty = False

    def get(self, key):
        """Returns the chunks associated to the file with the given key, or
        `None` if the key is not in the cache.
        """
        if self._data is None:
            self._load()
        return self._data.get(key)

    def key_of(self, contents, encoding = "utf-8"):
        """Returns the hash key corresponding to the file with the given
        contents.
        """
        if not isinstance(contents, bytes):
            contents = contents.encode("utf-8")

        key = self._hash()
        key.update(contents)
        return key.hexdigest()

    def put(self, key, chunks):
        """Stores some chunks associated to the file with the given key."""
        self._data[key] = chunks
        self._dirty = True


if __name__ == "__main__":
    main()
