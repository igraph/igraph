"""
Layout-related code in the IGraph library.

This package contains the implementation of the L{Layout} object.
"""

from copy import copy
from igraph.statistics import RunningMean
import math

__license__ = """
Copyright (C) 2006-2007  Gabor Csardi <csardi@rmki.kfki.hu>,
Tamas Nepusz <ntamas@rmki.kfki.hu>

MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary

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

class Layout(object):
    """Represents the layout of a graph.

    A layout is practically a list of coordinates in an n-dimensional
    space. This class is generic in the sense that it can store coordinates
    in any n-dimensional space.

    Layout objects are not associated directly with a graph. This is deliberate:
    there were times when I worked with almost identical copies of the same
    graph, the only difference was that they had different colors assigned to
    the vertices. It was particularly convenient for me to use the same layout
    for all of them, especially when I made figures for a paper. However,
    C{igraph} will of course refuse to draw a graph with a layout that has
    less coordinates than the node count of the graph.

    Layouts behave exactly like lists when they are accessed using the item
    index operator (C{[...]}). They can even be iterated through. Items
    returned by the index operator are only copies of the coordinates,
    but the stored coordinates can be modified by directly assigning to
    an index.

        >>> layout = Layout([(0, 1), (0, 2)])
        >>> coords = layout[1]
        >>> print coords
        [0, 2]
        >>> coords = (0, 3)
        >>> print layout[1]
        [0, 2]
        >>> layout[1] = coords
        >>> print layout[1]
        [0, 3]
    """
    
    def __init__(self, coords=[], dim=None):
        """Constructor.

        @param coords: the coordinates to be stored in the layout.
        @param dim: the number of dimensions. If C{None}, the number of
          dimensions is determined automatically from the length of the
          first item of the coordinate list. An exception is thrown if
          the coordinate list is empty and this parameter is C{None}.
          Generally, this should be given if the length of the coordinate
          list is zero, otherwise it should be left as is.
        @raise ValueError: if the coordinate list is empty and the
          number of dimensions is not given.
        """
        self._coords = map(list, coords)
        if dim is None:
            if len(self._coords) == 0:
                raise ValueError, "the number of dimensions must be given if the coordinate list is empty"
            else:
                self._dim = len(self._coords[0])
        else:
            self._dim = int(dim)
            for row in self._coords:
                if len(row) != self._dim:
                    raise ValueError, "all items in the coordinate list must have a length of %d" % self._dim

    def __len__(self): return len(self._coords)
    def __getitem__(self, idx): return copy(self._coords[idx])
    def __setitem__(self, idx, value):
        if len(value) != self._dim:
            raise ValueError, "assigned item must have %d elements" % self._dim
        self._coords[idx] = list(value)
    def __delitem__(self, idx): del self._coords[idx]
    def __copy__(self):
        return self.__class__(copy(self._coords), self._dim)

    def _get_dim(self): return self._dim
    dim = property(_get_dim, "the number of dimensions")
    def _get_coords(self): return [copy(row) for row in self._coords]
    coords = property(_get_coords, "the coordinates as a list of lists")

    def append(self, value):
        """Appends a new point to the layout"""
        if len(value) < self._dim:
            raise ValueError, "appended item must have %d elements" % self._dim
        self._coords.append(map(float, value[0:self._dim]))

    def mirror(self, dim):
        """Mirrors the layout along the given dimension(s)

        @param dim: the list of dimensions or a single dimension
        """
        if isinstance(dim, int): dim = [dim]
        else: dim = map(int, dim)

        _n = self._dim
        vec = [1]*_n
        for d in dim: vec[d] *= -1

        self._coords = [[row[i]*vec[i] for i in xrange(_n)] for row in self._coords]


    def rotate(self, degree, dim1=0, dim2=1):
        """Rotates the layout by the given degrees along the plane defined by
        the given two dimensions"""
        ca, sa = math.cos(degree * math.pi / 180.), math.sin(degree * math.pi / 180.)
        
        for idx, row in enumerate(self._coords): 
            new_row = list(row)
            new_row[dim1] = ca*row[dim1] - sa*row[dim2]
            new_row[dim2] = sa*row[dim1] + ca*row[dim2]
            self._coords[idx] = new_row


    def scale(self, *args, **kwds):
        """Scales the layout.

        Scaling parameters can be provided either through the C{scale} keyword
        argument or through plain unnamed arguments. If a single integer or
        float is given, it is interpreted as a uniform multiplier to be applied
        on all dimensions. If it is a list or tuple, its length must be equal to
        the number of dimensions in the layout, and each element must be an
        integer or float describing the scaling coefficient in one of the
        dimensions.

        @keyword scale: scaling coefficients (integer, float, list or tuple)
        @keyword origin: the origin of scaling (this point will stay in place).
          Optional, defaults to the origin of the coordinate system being used.
        """
        origin = list(kwds.get("origin", [0.]*self._dim))
        if len(origin) != self._dim:
            raise ValueError, "origin must have %d dimensions" % self._dim

        scaling = kwds.get("scale") or args
        if type(scaling) == int or type(scaling) == float: scaling = [scaling]
        if len(scaling) == 0:
            raise ValueError, "scaling factor must be given"
        elif len(scaling) == 1:
            if type(scaling[0]) == int or type(scaling[0]) == float:
                scaling = scaling*self._dim
            else:
                scaling = scaling[0]
        if len(scaling) != self._dim:
            raise ValueError, "scaling factor list must have %d elements" % self._dim

        for idx, row in enumerate(self._coords):
            self._coords[idx] = [(row[d]-origin[d])*scaling[d]+origin[d] \
                                 for d in xrange(self._dim)]

    def translate(self, *args, **kwds):
        """Translates the layout.

        The translation vector can be provided either through the C{v} keyword
        argument or through plain unnamed arguments. If unnamed arguments are
        used, the vector can be supplied as a single list (or tuple) or just as
        a series of arguments. In all cases, the translation vector must have
        the same number of dimensions as the layout.

        @keyword v: the translation vector
        """
        v = kwds.get("v") or args
        if len(v) == 0:
            raise ValueError, "translation vector must be given"
        elif len(v) == 1 and type(v[0]) != int and type(v[0]) != float:
            v = v[0]
        if len(v) != self._dim:
            raise ValueError, "translation vector must have %d dimensions" % self._dim

        for idx, row in enumerate(self._coords):
            self._coords[idx] = [row[d]+v[d] for d in xrange(self._dim)]


    def to_radial(self, min_angle = 100, max_angle = 80, \
        min_radius=0.0, max_radius=1.0):
        """Converts a planar layout to a radial one

        This method applies only to 2D layouts. The X coordinate of the
        layout is transformed to an angle, with min(x) corresponding to
        the parameter called I{min_angle} and max(y) corresponding to
        I{max_angle}. Angles are given in degrees, zero degree corresponds
        to the direction pointing upwards. The Y coordinate is
        interpreted as a radius, with min(y) belonging to the minimum and
        max(y) to the maximum radius given in the arguments.

        This is not a fully generic polar coordinate transformation, but
        it is fairly useful in creating radial tree layouts from ordinary
        top-down ones (that's why the Y coordinate belongs to the radius).
        It can also be used in conjunction with the Fruchterman-Reingold
        layout algorithm via its I{miny} and I{maxy} parameters (see
        L{Graph.layout_fruchterman_reingold}) to produce radial layouts
        where the radius belongs to some property of the vertices.

        @param min_angle: the angle corresponding to the minimum X value
        @param max_angle: the angle corresponding to the maximum X value
        @param min_radius: the radius corresponding to the minimum Y value
        @param max_radius: the radius corresponding to the maximum Y value
        """
        bbox = self.bounding_box()
        if len(bbox) != 4: raise TypeError, "implemented only for 2D layouts"

        while min_angle > max_angle: max_angle += 360
        while min_angle > 360:
            min_angle -= 360
            max_angle -= 360
        while min_angle < 0:
            min_angle += 360
            max_angle += 360

        minx, miny, maxx, maxy = bbox
        rx = (max_angle - min_angle) / (maxx-minx)
        rx *= math.pi / 180.
        min_angle *= math.pi / 180.
        ry = (max_radius - min_radius) / (maxy-miny)
        for idx, (x, y) in enumerate(self._coords):
            alpha, r = (x-minx) * rx + min_angle, (y-miny) * ry + min_radius
            self._coords[idx] = math.cos(alpha)*r, -math.sin(alpha)*r


    def transform(self, function, *args, **kwds):
        """Performs an arbitrary transformation on the layout

        Additional positional and keyword arguments are passed intact to
        the given function.

        @param function: a function which receives the coordinates as a
          tuple and returns the transformed tuple.
        """
        self._coords = [list(function(tuple(row), *args, **kwds)) \
            for row in self._coords]


    def centroid(self):
        """Returns the centroid of the layout.

        The centroid of the layout is the arithmetic mean of the points in
        the layout.
        
        @return: the centroid as a list of floats"""
        centroid = [RunningMean() for d in xrange(self._dim)]
        for row in self._coords:
            for d in xrange(self._dim):
                centroid[d] << row[d]
        return [rm.mean for rm in centroid]

    def bounding_box(self, border=0):
        """Returns the bounding box of the layout.

        The bounding box of the layout is the smallest box enclosing all the
        points in the layout.

        @param border: this value gets subtracted from the minimum bounds
          and gets added to the maximum bounds before returning the coordinates
          of the box. Defaults to zero.
        @return: the coordinates of the lower left and the upper right corner
          of the box. "Lower left" means the minimum coordinates and "upper right"
          means the maximum."""
        mins, maxs = [], []
        for d in xrange(self._dim):
            col = [row[d] for row in self._coords]
            mins.append(min(col)-border)
            maxs.append(max(col)+border)
        mins.extend(maxs)
        return tuple(mins)

    def center(self, *args, **kwds):
        """Centers the layout around the given point.

        The point itself can be supplied as multiple unnamed arguments, as a
        simple unnamed list or as a keyword argument. This operation moves
        the centroid of the layout to the given point. If no point is supplied,
        defaults to the origin of the coordinate system.

        @keyword p: the point where the centroid of the layout will be after
          the operation."""
        center = kwds.get("p") or args
        if len(center) == 0:
            center = [0.] * self._dim
        elif len(center) == 1 and type(center[0]) != int \
            and type(center[0]) != float:
            center = center[0]
        if len(center) != self._dim:
            raise ValueError, "the given point must have %d dimensions"%self._dim
        centroid = self.centroid()
        vec = [center[d]-centroid[d] for d in xrange(self._dim)]
        self.translate(vec)


    def copy(self):
        """Creates an exact copy of the layout."""
        return copy(self)

