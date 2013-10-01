import unittest

from igraph import *

class HomepageExampleTests(unittest.TestCase):
    """Smoke tests for the Python examples found on the homepage to ensure
    that they do not break."""

    def testErdosRenyiComponents(self):
        g = Graph.Erdos_Renyi(n=300, m=250)
        colors = ["lightgray", "cyan", "magenta", "yellow", "blue", "green", "red"]
        components = g.components()
        for component in components:
            color = colors[min(6, len(components)-1)]
            g.vs[component]["color"] = color

        # No plotting here, but we calculate the FR layout
        fr = g.layout("fr")

    def testKautz(self):
        g = Graph.Kautz(m=3, n=2)
        adj = g.get_adjacency()
        # Plotting omitted

    def testMSTofGRG(self):
        def distance(p1, p2):
            return ((p1[0]-p2[0]) ** 2 + (p1[1]-p2[1]) ** 2) ** 0.5

        g = Graph.GRG(100, 0.2)
        layout = Layout(zip(g.vs["x"], g.vs["y"]))

        weights = [distance(layout[edge.source], layout[edge.target]) \
                for edge in g.es]
        max_weight = max(weights)
        g.es["width"] = [6 - 5*weight / max_weight for weight in weights]
        mst = g.spanning_tree(weights)
        # Plotting omitted

def suite():
    homepage_example_suite = unittest.makeSuite(HomepageExampleTests)
    return unittest.TestSuite([homepage_example_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

        
