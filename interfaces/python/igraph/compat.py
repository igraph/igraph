# vim:ts=4:sw=4:sts=4:et
# -*- coding: utf-8 -*-
"""
Compatibility methods and backported versions of newer Python features
to enable igraph to run on Python 2.5.
"""

import sys

__license__ = u"""\
Copyright (C) 2006-2012  Tamás Nepusz <ntamas@gmail.com>
Pázmány Péter sétány 1/a, 1117 Budapest, Hungary

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
02110-1301 USA
"""

#############################################################################
# Simulating math.isnan

try:
    from math import isnan
except ImportError:
    def isnan(num):
        return num != num

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


