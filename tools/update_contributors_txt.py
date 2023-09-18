#!/usr/bin/env python3
"""Update the contents of CONTRIBUTORS.txt based on .all-contributorsrc."""

from json import load
from os import chdir
from pathlib import Path


HEADER = """\
Thanks goes to these wonderful people:

"""

FOOTER = """\

This project follows the [all-contributors][1] specification. Contributions of any kind welcome!

This file is an automatically generated, plain-text version of CONTRIBUTORS.md.

[1]: https://github.com/all-contributors/all-contributors
"""

def main():
    root_dir = Path(__file__).parent.parent
    chdir(root_dir)

    with (root_dir / ".all-contributorsrc").open("r") as fp:
        contributors = load(fp)["contributors"]

    with (root_dir / "CONTRIBUTORS.txt").open("w") as fp:
        fp.write(HEADER)
        for c in contributors:
            if c["name"]:
                fp.write(f"{c['name']} (@{c['login']})\n")
            else:
                fp.write(f"@{c['login']}\n")
        fp.write(FOOTER)


if __name__ == "__main__":
    main()
