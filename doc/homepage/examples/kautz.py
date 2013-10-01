from igraph import *

g = Graph.Kautz(m=3, n=2)
adj = g.get_adjacency()
fig = Plot(bbox=(480, 480))
fig.add(g, layout="fr", vertex_label=None)
fig.add(adj, bbox=(360, 0, 480, 120), grid_width=0, opacity=0.7)
fig.show()
