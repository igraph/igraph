"""Utility functions that cannot be categorised anywhere else"""

from contextlib import contextmanager
from itertools import chain

import os
import tempfile

__license__ = "GPL"

__all__ = ["rescale", "safemin", "safemax"]
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

def rescale(values, out_range = (0., 1.), in_range = None, clamp = False):
    """Rescales a list of numbers into a given range.

    `out_range` gives the range of the output values; by default, the minimum
    of the original numbers in the list will be mapped to the first element
    in the output range and the maximum will be mapped to the second element.
    Elements between the minimum and maximum values in the input list will be
    interpolated linearly between the first and second values of the output
    range.

    `in_range` may be used to override which numbers are mapped to the first
    and second values of the output range. This must also be a tuple, where
    the first element will be mapped to the first element of the output range
    and the second element to the second.

    If `clamp` is ``True``, elements which are outside the given `out_range`
    after rescaling are clamped to the output range to ensure that no number
    will be outside `out_range` in the result.
    
    Examples:
        
        >>> rescale(range(5), (0, 8))
        [0.0, 2.0, 4.0, 6.0, 8.0]
        >>> rescale(range(5), (2, 10))
        [2.0, 4.0, 6.0, 8.0, 10.0]
        >>> rescale(range(5), (0, 4), (1, 3))
        [-2.0, 0.0, 2.0, 4.0, 6.0]
        >>> rescale(range(5), (0, 4), (1, 3), clamp=True)
        [0.0, 0.0, 2.0, 4.0, 4.0]
        >>> rescale([0]*5, (1, 3))
        [2.0, 2.0, 2.0, 2.0, 2.0]
    """
    if in_range is None:
        mi, ma = min(values), max(values)
    else:
        mi, ma = in_range

    ratio = float(ma - mi)
    if not ratio:
        return [(out_range[0] + out_range[1]) / 2.] * len(values)

    min_out, max_out = map(float, out_range)
    ratio = (max_out - min_out) / ratio
    result = [(x - mi) * ratio + min_out for x in values]

    if clamp:
        return [max(min(x, max_out), min_out) for x in result]
    else:
        return result

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

def safemax(iterable, default=0):
    """Safer variant of C{max()} that returns a default value if the iterable
    is empty.
    
    Example:
        
        >>> safemax([-5, 6, 4])
        4
        >>> safemax([])
        0
        >>> safemax((), 2)
        2
    """
    it = iter(iterable)
    try:
        first = it.next()
    except StopIteration:
        return default
    else:
        return max(chain([first], it))

def safemin(iterable, default=0):
    """Safer variant of C{min()} that returns a default value if the iterable
    is empty.
    
    Example:
        
        >>> safemin([-5, 6, 4])
        -5
        >>> safemin([])
        0
        >>> safemin((), 2)
        2
    """
    it = iter(iterable)
    try:
        first = it.next()
    except StopIteration:
        return default
    else:
        return min(chain([first], it))
