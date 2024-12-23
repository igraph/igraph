from infomap import Infomap

im = Infomap(two_level=True, silent=True)

# Set the start id for bipartite nodes
im.bipartite_start_id = 5

# Add weight as an optional third argument
im.add_link(5, 0)
im.add_link(5, 1)
im.add_link(5, 2)
im.add_link(6, 2, 0.5)
im.add_link(6, 3)
im.add_link(6, 4)

im.run()

print(f"Found {im.num_top_modules} modules with codelength {im.codelength:.8f} bits")

print("\n#node module:")
for node, module in im.modules:
    print(node, module)
