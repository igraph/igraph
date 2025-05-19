from infomap import Infomap

im = Infomap(silent=True)
im.read_file("../networks/states.net")
im.run()

print("source target weight")
for source, target, weight in im.get_links():  # or im.links:
    print(source, target, weight)

print("source target flow")
for source, target, flow in im.get_links(data="flow"):  # or im.flow_links
    print(source, target, f"{flow:.4f}")
