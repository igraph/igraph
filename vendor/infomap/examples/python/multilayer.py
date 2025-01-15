from infomap import Infomap


def print_result(im):
    print(
        f"Found {im.num_top_modules} modules with codelength {im.codelength:.8f} bits"
    )

    print("#layer_id node_id module_id:")
    for node in im.nodes:
        print(f"{node.layer_id} {node.node_id} {node.module_id}")


im = Infomap(silent=True)
print("\nAdding multilayer network...")

# from (layer1, node1) to (layer2, node2) with optional weight
im.add_multilayer_link((2, 1), (1, 2), 1.0)
im.add_multilayer_link((1, 2), (2, 1), 1.0)
im.add_multilayer_link((3, 2), (2, 3), 1.0)

im.run()

print_result(im)

# Add only intra-layer links and let Infomap provide
# inter-layer links by relaxing the random walker's
# constraint to its current layer
im = Infomap(silent=True)
print("\nAdding intra-layer network...")

# Add layer_id, source_node_id, target_node_id and optional weight
im.add_multilayer_intra_link(1, 1, 2)
im.add_multilayer_intra_link(1, 2, 3)
im.add_multilayer_intra_link(2, 1, 3)
im.add_multilayer_intra_link(2, 3, 4)
# Optionally add inter-layer links but that will prevent Infomap
# from simulating those
# im.add_multilayer_inter_link(1, 1, 2)

im.run()

print_result(im)
