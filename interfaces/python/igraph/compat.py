# vim:ts=4:sw=4:sts=4:et
"""
Compatibility methods and backported versions of newer Python features
to enable igraph to run on Python 2.5.
"""

import sys

#############################################################################
# Providing @property.setter syntax for Python 2.5

if sys.version_info < (2, 6):
    _property = property
    class property(property):
        def __init__(self, fget, *args, **kwds):
            self.__doc__ = fget.__doc__
            super(property, self).__init__(fget, *args, **kwds)

        def setter(self, fset):
            cls_ns = sys._getframe(1).f_locals
            for k, v in cls_ns.iteritems():
                if v == self:
                    propname = k
                    break
            cls_ns[propname] = property(self.fget, fset, self.fdel,
                    self.__doc__)
            return cls_ns[propname]
else:
    property = __builtins__["property"]


