import pathlib

from infomap import Infomap

im = Infomap(silent=True)

name = "ninetriangles"
filename = f"../networks/{name}.net"

# You can read a network with the method read_file,
# which by default will accumulate to existing network data
im.read_file(filename, accumulate=False)

im.run(num_trials=5)

print(
    f"Found {im.max_depth} levels with {im.num_leaf_modules} leaf modules in {im.num_top_modules} top modules and codelength: {im.codelength:.8f} bits"
)
print(f"All codelengths: {im.codelengths}")

print("Tree:\n# path node_id module_id flow")
for node in im.nodes:
    print(f"{node.path} {node.node_id} {node.module_id} {node.flow:.8f}")

for module_level in range(1, im.max_depth):
    print(
        f"Modules at level {module_level}: {tuple(im.get_modules(module_level).values())}"
    )

print("\nModules at all levels:")
for node_id, modules in im.get_multilevel_modules().items():
    print(f"{node_id}: {modules}")

pathlib.Path("output").mkdir(exist_ok=True)
print(f"Writing top level modules to output/{name}.clu...")
im.write(f"output/{name}.clu")

print(f"Writing second level modules to output/{name}_level2.clu...")
im.write(f"output/{name}_level2.clu", depth_level=2)

print(f"Writing bottom level modules to output/{name}_level-1.clu...")
im.write(f"output/{name}_level-1.clu", depth_level=-1)

print(f"Writing tree to output/{name}.tree...")
im.write(f"output/{name}.tree")

print("Read back .clu file and only calculate codelength...")
im2 = Infomap(
    silent=True, two_level=True, no_infomap=True, cluster_data=f"output/{name}.clu"
)
im2.read_file(filename)
im2.run()
print(
    f"Found {im2.max_depth} levels with {im2.num_top_modules} top modules and codelength: {im2.codelength:.8f} bits"
)

print("Done!")
