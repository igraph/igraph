from infomap import Infomap

im = Infomap(two_level=True, silent=True)

im.set_name(1, "PRE")
im.set_name(2, "SCIENCE")
im.set_name(3, "PRL")
im.set_name(4, "BIO")

im.add_state_node(0, 1)
im.add_state_node(1, 2)
im.add_state_node(2, 3)
im.add_state_node(3, 2)
im.add_state_node(4, 2)
im.add_state_node(5, 4)

im.add_link(0, 1)
im.add_link(1, 2)
im.add_link(3, 2)
im.add_link(4, 5)

im.run()

print(f"Found {im.num_top_modules} modules with codelength {im.codelength:.8f} bits")

print("\n#node_id module")
for node, module in im.get_modules(states=True).items():
    print(node, module)

print("\nState nodes:")
print("#state_id node_id module_id")
for node in im.nodes:
    print(node.state_id, node.node_id, node.module_id)

print(
    "\nPhysical nodes (merging state nodes with same physical node id within modules):"
)
print("#node_id module_id")
for node in im.physical_tree:
    if node.is_leaf:
        print(node.node_id, node.module_id)

# for node in im.physical_nodes:
#     print(node.node_id, node.module_id)
