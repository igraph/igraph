"""Utility functions that cannot be categorised anywhere else"""

from contextlib import contextmanager

import os
import tempfile

__license__ = "GPL"

__all__ = ["named_temporary_file"]
__docformat__ = "restructuredtext en"

@contextmanager
def named_temporary_file(*args, **kwds):
    """Context manager that creates a named temporary file and
    returns its name.

    All parameters are passed on to `tempfile.mkstemp`, see
    its documentation for more info.
    """
    handle, tmpfile = tempfile.mkstemp(*args, **kwds)
    os.close(handle)
    try:
        yield tmpfile
    finally:
        os.unlink(tmpfile)

def consecutive_pairs(iterable, circular=False):
    """Returns consecutive pairs of items from the given iterable.

    When `circular` is ``True``, the pair consisting of the last
    and first elements is also returned.

    Example:
        
        >>> list(consecutive_pairs(range(5)))
        [(0, 1), (1, 2), (2, 3), (3, 4)]
        >>> list(consecutive_pairs(range(5), circular=True))
        [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
        >>> list(consecutive_pairs([]))
        []
        >>> list(consecutive_pairs([], circular=True))
        []
        >>> list(consecutive_pairs([0]))
        []
        >>> list(consecutive_pairs([0], circular=True))
        [(0, 0)]
    """
    it = iter(iterable)

    try:
        prev = it.next()
    except StopIteration:
        return
    first = prev

    for item in it:
        yield prev, item
        prev = item

    if circular:
        try:
            yield item, first
        except UnboundLocalError:
            yield first, first

