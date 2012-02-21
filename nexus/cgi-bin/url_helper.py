"""
Helper functions for generating URLs and links.
"""

# List of symbols to be exported from this module
__all__ = [
    "link_to", "link_to_dataset", "link_to_dataset_id", "link_to_tag",
    "link_to_licence", "url_for_dataset", "url_for_tag", "url_for_licence"
]

from web import websafe

def link_to(url, caption):
    """Generates a link to the given URL. `caption` is optional and may
    specify the caption of the URL. If `caption` is ``None``, the URL itself
    is used."""
    if caption is None:
        caption = url
    return '<a href="%s">%s</a>' % (websafe(url), websafe(caption))

def link_to_dataset(base, dataset, caption=None, bullet=""):
    """Generates a link to the given dataset. `caption` is optional and may
    specify the caption of the URL. If `caption` is ``None``, the name of
    the dataset is used. `bullet` specifies the bullet to be shown before
    the caption."""
    if caption is None:
        caption = dataset.name
    if bullet:
        caption = u"%s%s" % (bullet, caption)
    return link_to(url_for_dataset(base, dataset.id), caption)

def link_to_licence(base, licence, caption=None, bullet=""):
    if caption is None:
        caption = ""
    if bullet:
        caption = u"%s%s" % (bullet, caption)
    return link_to(url_for_licence(base, licence), caption)

def link_to_dataset_id(base, dataset_id, caption):
    """Generates a link to the dataset with the given ID. `caption` specifies
    the caption of the URL."""
    return link_to(url_for_dataset(base, dataset_id), caption)

def link_to_tag(base, tag, caption=None):
    """Generates a link to the page of the given tag. `caption` is optional
    and may specify the caption of the URL. If `caption` is ``None``, the
    tag itself is used."""
    if caption is None:
        caption = tag.tag
    caption = websafe(caption)
    return link_to(url_for_tag(base, tag), caption)

def url_for_dataset(base, dataset, format="html"):
    """Returns the URL where a given dataset is to be found.
    `dataset` is the ID of the dataset or the dataset itself, `format`
    is the format of the dataset."""
    if not isinstance(dataset, int):
        dataset = dataset.id
    return "%s/api/dataset_info?id=%s&format=%s" % (base, dataset, format)

def url_for_licence(base, licence, format="html"):
    return getbase() + "/api/licence?id=%s&format=%s" % (base, 
                                                         licence, format)

def url_for_tag(base, tag, format="html"):
    """Returns the URL where the list of datasets with a given tag are
    to be found. `tag` is the tag itself or any object with a ``tag``
    attribute). `format` is the desired format of the result."""
    if hasattr(tag, "tag"):
        tag = tag.tag
    return "%s/api/dataset_info?tag=%s&format=%s" % (base, tag, format)

