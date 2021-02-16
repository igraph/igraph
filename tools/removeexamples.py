#!/usr/bin/env python3

import sys
from xml.etree.ElementTree import ElementTree


def usage():
    print(sys.argv[0], "<infile> <outfile>")


def main():
    if len(sys.argv) != 3:
        usage()
        sys.exit(2)

    # Read in
    tree = ElementTree()
    tree.parse(sys.argv[1])

    # Remove examples
    examples = tree.findall(".//example")
    for ex in examples:
        prog = ex.find("programlisting")
        ex.remove(prog)

    # Write result
    tree.write(sys.argv[2])


if __name__ == "__main__":
    main()
