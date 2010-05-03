"""Utility functions that cannot be categorised anywhere else"""

from contextlib import contextmanager
from tempfile import mkstemp

import os

__license__ = "GPL"

__all__ = ["named_temporary_file"]


@contextmanager
def named_temporary_file(*args, **kwds):
    """Context manager that creates a named temporary file and
    returns its name.

    All parameters are passed on to L{tempfile.mkstemp}, see
    its documentation for more info.

    @see: tempfile.mkstemp
    """
    handle, tmpfile = mkstemp(*args, **kwds)
    os.close(handle)
    try:
        yield tmpfile
    finally:
        os.unlink(tmpfile)

