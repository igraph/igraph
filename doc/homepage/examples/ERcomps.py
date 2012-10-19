from igraph import *

g = Graph.Erdos_Renyi(n=300, m=250)
colors = ["lightgray", "cyan", "magenta", "yellow", "blue", "green", "red"]
for component in g.components():
  color = colors[min(6, len(component)-1)]
  for vidx in component: g.vs[vidx]["color"] = color
  
plot(g, layout="fr", vertex_label=None)
