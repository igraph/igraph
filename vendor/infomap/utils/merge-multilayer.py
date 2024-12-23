#!/usr/bin/env python3
import sys
from collections import defaultdict


def merge_multilayer(input_files, out_file):
    nodes = {}
    links = defaultdict(list)

    for layer, file in enumerate(input_files):
        print(f"Reading {file}...")
        with open(file) as f:
            lines = f.readlines()

        context = None
        for line in map(str.strip, lines):
            if line.startswith("#"):
                continue
            if line.startswith("*"):
                context = line.split(maxsplit=1)[0].lower()
                continue

            if context in ("*vertices", "*nodes"):
                node_id, name = line.split(maxsplit=1)
                if node_id in nodes and nodes[node_id] != name:
                    print(f"Warning: node {node_id} has a different names ({nodes[node_id]} and {name})")
                nodes[node_id] = name
            elif context in ("*links", "*edges") or context is None:
                try:
                    source, target, weight = line.split()
                except ValueError:
                    source, target = line.split()
                    weight = 1

                links[layer + 1].append((source, target, weight))
            else:
                print(f"Unknown heading {context}")
                sys.exit(1)

    print(f"Writing {out_file}")
    with open(outfile, "w") as f:
        f.write(f"# {' '.join(sys.argv)}\n")
        if len(nodes) != 0:
            f.write(f"*Vertices {len(nodes)}\n")
            f.writelines(f"{node} {name}\n" for node, name in nodes.items())
        f.write("*Intra\n")
        for layer, layer_links in links.items():
            f.writelines(f"{layer} {source} {target} {weight}\n" for source, target, weight in layer_links)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Merge pajek network files into a multilayer network with intra-layer links.")
        print("Usage:")
        print("\tmerge-multilayer.py pajek_file1 [pajek_file2 pajek_fileN] outfile")
        sys.exit(0)

    input_files = sys.argv[1:-1]
    outfile = sys.argv[-1]
    merge_multilayer(input_files, outfile)
