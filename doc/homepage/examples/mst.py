from igraph import *

def distance(p1, p2):
    return ((p1[0]-p2[0]) ** 2 + (p1[1]-p2[1]) ** 2) ** 0.5
    
g = Graph.GRG(100, 0.2)
layout = Layout(zip(g.vs["x"], g.vs["y"]))

weights = [distance(layout[edge.source], layout[edge.target]) for edge in g.es]
max_weight = max(weights)
g.es["width"] = [6 - 5*weight/max_weight for weight in weights]
mst = g.spanning_tree(weights)

fig = Plot()
fig.add(g, layout=layout, opacity=0.25, vertex_label=None)
fig.add(mst, layout=layout, edge_color="red", vertex_label=None)
fig.show()
