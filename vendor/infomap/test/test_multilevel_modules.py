from collections import defaultdict
from infomap import Infomap
import pytest


def multilevel_modules(im, states):
  modules = defaultdict(list)

  for level in range(0, im.num_levels - 1):
    module_id = 1
    prev_path = None
    for node in im.get_nodes(states=states):
      if prev_path is None:
        prev_path = node.path
      index = min(level, len(node.path) - 2) + 1
      if node.path[:index] != prev_path[:index]:
        module_id += 1
      node_id = node.state_id if states else node.node_id
      modules[node_id].append(module_id)
      prev_path = node.path
      print(f"{node_id=}, {module_id=}, {node.path=}")
  
  return {node: tuple(m) for node, m in modules.items()}


def test_multilevel_modules_states():
  for _ in range(1):
    im = Infomap(num_trials=5, silent=True)
    im.read_file("test/data/multilayer.net")
    im.run()

    assert im.num_top_modules == 7
    assert im.num_levels == 4

    assert im.get_multilevel_modules(states=True) == multilevel_modules(im, states=True)

    with pytest.raises(RuntimeError):
      im.get_multilevel_modules(states=False)


def test_multilevel_modules():
  im = Infomap(num_trials=10, silent=True)
  im.read_file("test/data/ninetriangles.net")
  im.run()

  assert im.num_top_modules == 5
  assert im.num_levels == 3

  assert sorted(im.get_multilevel_modules(states=False).values()) == sorted(multilevel_modules(im, states=False).values())
