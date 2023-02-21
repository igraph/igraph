#!/usr/bin/env python3
"""Strips the license notices from the bundled igraph examples so they are
not visible in the documentation.
"""

import argparse
import sys

from pathlib import Path
from typing import IO


def strip_license_notice(infp: IO[str], outfp: IO[str]) -> None:
    seen_code = False

    for line in infp:
        if not seen_code:
            # This is an approximation; we consider the first line that starts
            # with a non-whitespace character, * or / as "real code", and we strip
            # everything before that
            if line and line[0] not in ' \n\t\r/*':
                seen_code = True
                outfp.write(line)
        else:
            outfp.write(line)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o", "--out-dir", help="output directory to put the stripped files in",
    )
    parser.add_argument('files', help="input files to process", nargs="*")
    options = parser.parse_args()

    if not options.out_dir:
        # We must have one or two args; one arg means that we are printing to
        # stdout
        if len(options.files) == 0:
            return 0
        elif len(options.files) > 2:
            parser.error(
                "At most two files (one input, one output) must be given "
                "when -o is not used"
            )
            return 1

        with Path(options.files[0]).open("r") as infp:
            if len(options.files) > 1:
                with Path(options.files[1]).open("w") as outfp:
                    strip_license_notice(infp, outfp)
            else:
                strip_license_notice(infp, sys.stdout)

    else:
        # We have an output dir so we can handle an arbitrary number of input
        # files
        for filename in options.files:
            path = Path(filename)
            with path.open("r") as infp:
                with (Path(options.out_dir) / path.name).open("w") as outfp:
                    strip_license_notice(infp, outfp)

    return 0

if __name__ == "__main__":
    sys.exit(main())
