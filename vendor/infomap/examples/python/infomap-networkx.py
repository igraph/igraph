import pathlib

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import networkx as nx

from infomap import Infomap

"""
Generate and draw a network with NetworkX, colored
according to the community structure found by Infomap.
"""


def draw_network(G):
    # position map
    pos = nx.spring_layout(G)
    # community index
    communities = [c - 1 for c in nx.get_node_attributes(G, "community").values()]
    num_communities = max(communities) + 1

    # color map from http://colorbrewer2.org/
    cmap_light = colors.ListedColormap(
        ["#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6"],
        "indexed",
        num_communities,
    )
    cmap_dark = colors.ListedColormap(
        ["#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a"],
        "indexed",
        num_communities,
    )

    # edges
    nx.draw_networkx_edges(G, pos)

    # nodes
    node_collection = nx.draw_networkx_nodes(
        G, pos=pos, node_color=communities, cmap=cmap_light
    )

    # set node border color to the darker shade
    dark_colors = [cmap_dark(v) for v in communities]
    node_collection.set_edgecolor(dark_colors)

    # Print node labels separately instead
    for n in G.nodes:
        plt.annotate(
            n,
            xy=pos[n],
            textcoords="offset points",
            horizontalalignment="center",
            verticalalignment="center",
            xytext=[0, 2],
            color=cmap_dark(communities[n]),
        )

    plt.axis("off")
    pathlib.Path("output").mkdir(exist_ok=True)
    print("Writing network figure to output/karate.png")
    plt.savefig("output/karate.png")
    # plt.show()


G = nx.karate_club_graph()

print("Building Infomap network from a NetworkX graph...")
im = Infomap(two_level=True, silent=True, num_trials=20)
im.add_networkx_graph(G)

print("Find communities with Infomap...")
im.run()

print(f"Found {im.num_top_modules} modules with codelength {im.codelength:.8f} bits")

nx.set_node_attributes(G, im.get_modules(), "community")

draw_network(G)
