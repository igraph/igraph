"""Utility functions for unit testing."""

import functools
import sys
import types

__all__ = ["skip", "skipIf"]

def _id(obj):
    return obj

try:
    from unittest import skip
except ImportError:
    # Provide basic replacement for unittest.skip
    def skip(reason):
        """Unconditionally skip a test."""
        def decorator(test_item):
            if not isinstance(test_item, (type, types.ClassType)):
                @functools.wraps(test_item)
                def skip_wrapper(*args, **kwds):
                    if reason:
                        sys.stderr.write("skipped, %s ... " % reason)
                    else:
                        sys.stderr.write("skipped, ")
                    return
                test_item = skip_wrapper
            return test_item
        return decorator

try:
    from unittest import skipIf
except ImportError:
    # Provide basic replacement for unittest.skipIf
    def skipIf(condition, reason):
        """Skip a test if the condition is true."""
        if condition:
            return skip(reason)
        return _id
