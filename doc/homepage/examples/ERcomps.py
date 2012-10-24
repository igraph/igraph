from igraph import *

g = Graph.Erdos_Renyi(n=300, m=250)
colors = ["lightgray", "cyan", "magenta", "yellow", "blue", "green", "red"]
components = g.components()
for component in components:
  color = colors[min(6, len(component)-1)]
  g.vs[component]["color"] = color
  
plot(g, layout="fr", vertex_label=None)
