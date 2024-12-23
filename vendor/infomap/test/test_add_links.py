from infomap import Infomap


def test_add_link():
    im = Infomap(silent=True)
    links = [
        [1, 2],
        [2, 3],
        [3, 1],
        [3, 4],
        [4, 5],
        [5, 6],
        [6, 4],
    ]

    for link in links:
        im.add_link(*link)

    im.run()

    assert im.num_top_modules == 2


def test_add_links_unweighted():
    im = Infomap(silent=True)
    links = [
        [1, 2],
        [2, 3],
        [3, 1],
        [3, 4],
        [4, 5],
        [5, 6],
        [6, 4],
    ]

    im.add_links(links)

    im.run()

    assert im.num_top_modules == 2


def test_add_links_weighted():
    im = Infomap(silent=True)
    links = [
        [1, 2, 0.5],
        [2, 3, 0.5],
        [3, 1, 0.5],
        [3, 4, 0.5],
        [4, 5, 0.5],
        [5, 6, 0.5],
        [6, 4, 0.5],
    ]

    im.add_links(links)

    im.run()

    assert im.num_top_modules == 2
