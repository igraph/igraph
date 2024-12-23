#!/usr/bin/env python3
import sys
from collections import namedtuple, defaultdict


def read_network(lines, directed=False):
    nodes = {}
    links = []
    # for undirected
    link_set = defaultdict(float)

    context = None
    for line in map(str.strip, lines):
        if line.startswith("#"):
            continue
        if line.startswith("*"):
            context = line.split(maxsplit=1)[0].lower()
            continue

        if context in ("*vertices", "*nodes"):
            node_id, name = line.split(maxsplit=1)
            nodes[int(node_id)] = name
        elif context in ("*links", "*edges") or context is None:
            try:
                source, target, weight = line.split()
            except ValueError:
                source, target = line.split()
                weight = 1

            source, target, weight = int(source), int(target), float(weight)
            links.append((source, target, weight))
            if not directed:
                link_set[source, target] += weight
        else:
            print(f"Unknown heading {context}")
            sys.exit(1)

    if not directed:
        for (source, target), weight in link_set.items():
            if (target, source) not in link_set:
                links.append((target, source, weight))

    return nodes, links


Node = namedtuple("Node", "id, name, out_links")


def create_state_network(network):
    nodes, links = network
    states = {}
    state_links = []
    max_id = 1 << 16

    original_num_nodes = len(nodes)
    if original_num_nodes == 0:
        for source, target, _ in links:
            nodes[source] = None
            nodes[target] = None

    for node_id in nodes:
        if node_id > max_id:
            raise RuntimeError(f"Max number of nodes ({1 << 16}) exceeded")
        name = nodes[node_id]
        nodes[node_id] = Node(node_id, name, out_links=[])

    for source, target, weight in links:
        nodes[source].out_links.append((target, weight))

    def create_state_id(prev, current):
        return prev << 16 | current

    num_links = len(links)

    for i, (prev, current, curr_weight) in enumerate(links):
        prev_state_id = create_state_id(prev, prev)
        states[prev_state_id] = prev

        current_state_id = create_state_id(prev, current)
        states[current_state_id] = current

        state_links.append((prev_state_id, current_state_id, curr_weight))

        for next_, next_weight in nodes[current].out_links:
            next_state_id = create_state_id(current, next_)
            states[next_state_id] = next_

            state_links.append((current_state_id, next_state_id, next_weight))

        print(f"\r{int(100 * (i + 1) / num_links)}%", end="")

    print()

    # restore nodes
    if original_num_nodes == 0:
        nodes = {}
    else:
        for node_id in nodes:
            node = nodes[node_id]
            nodes[node_id] = node.name

    return nodes, states, state_links


def write_states(state_network, f):
    nodes, states, links = state_network
    f.write(f"# {' '.join(sys.argv)}\n")

    if len(nodes) != 0:
        f.write(f"*Vertices {len(nodes)}\n")
        f.write("# node_id name\n")
        f.writelines(f"{node} {name}\n" for node, name in nodes.items())

    f.write(f"*States {len(states)}\n")
    f.write("# state_id node_id\n")
    f.writelines(f"{state} {node}\n" for state, node in states.items())

    f.write(f"*Links {len(links)}\n")
    f.write("# source target weight\n")
    f.writelines(f"{source} {target} {weight}\n" for source, target, weight in links)


def network_to_states(network_file, states_file, directed=False):
    print(f"Reading network {network_file}...")
    with open(network_file) as f:
        lines = f.readlines()

    print("Creating state network...")
    network = read_network(lines, directed=directed)
    state_network = create_state_network(network)
    print(f"Created {len(state_network[1])} state nodes and {len(state_network[2])} links.")

    print(f"Writing state network {states_file}...")
    with open(states_file, "w") as f:
        write_states(state_network, f)


if __name__ == "__main__":
    directed = len(sys.argv) != 1 and sys.argv[1] == "-d"

    if len(sys.argv) < 3 or (len(sys.argv) == 3 and directed):
        print("Create second-order state network from first-order network.")
        print("Usage:")
        print("\tnetwork-to-states.py [-d] network_file states_file\n")
        print("\t-d\tInput network is directed")
        sys.exit(0)

    arg_index = 2 if directed else 1
    network_file, states_file = sys.argv[arg_index:]

    network_to_states(network_file, states_file, directed)
