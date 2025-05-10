from infomap import Infomap
from operator import itemgetter


def test_iter_physical_on_physical():
    im = Infomap(num_trials=10, silent=True)
    im.read_file("examples/networks/twotriangles.net")
    im.run()

    modules = sorted([(node.node_id, node.module_id) for node in im.physical_nodes], key=itemgetter(0))
    assert modules == [(1, 1), (2, 1), (3, 1), (4, 2), (5, 2), (6, 2)]

def test_iter_physical_on_states():
    im = Infomap(num_trials=10, silent=True)
    im.read_file("examples/networks/states.net")
    im.run()

    modules = sorted([(node.state_id, node.node_id, node.module_id) for node in im.physical_nodes], key=itemgetter(0))
    assert modules == [(1, 1, 1), (2, 2, 1), (3, 3, 1), (4, 1, 2), (5, 4, 2), (6, 5, 2)]


def test_iter_physical_reliability():
    
    for _ in range(100):
        im = Infomap(num_trials=10, silent=True)
        im.read_file("examples/networks/states.net")
        im.run()
        
        modules = [(node.node_id, node.module_id) for node in im.physical_nodes]
        assert modules == [(1, 1), (2, 1), (3, 1), (1, 2), (4, 2), (5, 2)]

def test_multilevel_modules_on_states():
    im = Infomap(silent=True)
    im.read_file("examples/networks/states.net")
    im.run()
    modules = [(node, modules) for node, modules in im.get_multilevel_modules(states=True).items()]
    assert modules == [(1, (1,)), (2, (1,)), (3, (1,)), (4, (2,)), (5, (2,)), (6, (2,))]

def test_multilevel_modules_on_physical():
    im = Infomap(silent=True)
    im.read_file("examples/networks/states.net")
    im.run()
    # todo: Use unit test function to assert exception is raised
    # RuntimeError: Cannot get multilevel modules on higher-order network without states.
    modules = im.get_multilevel_modules(states=False)



if __name__ == "__main__":
    test_iter_physical_on_physical()
    test_iter_physical_on_states()
    test_iter_physical_reliability()
    test_multilevel_modules_on_states()
