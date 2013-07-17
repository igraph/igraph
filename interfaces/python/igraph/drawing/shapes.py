# vim:ts=4:sw=4:sts=4:et
# -*- coding: utf-8 -*-
"""
Shape drawing classes for igraph

Vertex shapes in igraph are usually referred to by short names like
C{"rect"} or C{"circle"}. This module contains the classes that
implement the actual drawing routines for these shapes, and a
resolver class that determines the appropriate shape drawer class
given the short name.

Classes that are derived from L{ShapeDrawer} in this module are
automatically registered by L{ShapeDrawerDirectory}. If you
implement a custom shape drawer, you must register it in
L{ShapeDrawerDirectory} manually if you wish to refer to it by a
name in the C{shape} attribute of vertices.
"""

from __future__ import division

__all__ = ["ShapeDrawerDirectory"]

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

from math import atan2, cos, pi, sin
import sys

from igraph.drawing.baseclasses import AbstractCairoDrawer
from igraph.drawing.utils import Point
from igraph.utils import consecutive_pairs

class ShapeDrawer(object):
    """Static class, the ancestor of all vertex shape drawer classes.
    
    Custom shapes must implement at least the C{draw_path} method of the class.
    The method I{must not} stroke or fill, it should just set up the current
    Cairo path appropriately."""

    @staticmethod
    def draw_path(ctx, center_x, center_y, width, height=None):
        """Draws the path of the shape on the given Cairo context, without
        stroking or filling it.

        This method must be overridden in derived classes implementing custom shapes
        and declared as a static method using C{staticmethod(...)}.

        @param ctx: the context to draw on
        @param center_x: the X coordinate of the center of the object
        @param center_y: the Y coordinate of the center of the object
        @param width: the width of the object
        @param height: the height of the object. If C{None}, equals to the width.
        """
        raise NotImplementedError("abstract class")

    # pylint: disable-msg=W0613
    @staticmethod
    def intersection_point(center_x, center_y, source_x, source_y, \
            width, height=None):
        """Determines where the shape centered at (center_x, center_y)
        intersects with a line drawn from (source_x, source_y) to
        (center_x, center_y).

        Can be overridden in derived classes. Must always be defined as a static
        method using C{staticmethod(...)}

        @param width: the width of the shape
        @param height: the height of the shape. If C{None}, defaults to the width
        @return: the intersection point (the closest to (source_x, source_y) if
            there are more than one) or (center_x, center_y) if there is no
            intersection
        """
        return center_x, center_y


class NullDrawer(ShapeDrawer):
    """Static drawer class which draws nothing.

    This class is used for graph vertices with unknown shapes"""
    names = ["null", "none", "empty", "hidden", ""]

    @staticmethod
    def draw_path(ctx, center_x, center_y, width, height=None):
        """Draws nothing."""
        pass


class RectangleDrawer(ShapeDrawer):
    """Static class which draws rectangular vertices"""
    names = "rectangle rect rectangular square box"

    @staticmethod
    def draw_path(ctx, center_x, center_y, width, height=None):
        """Draws a rectangle-shaped path on the Cairo context without stroking
        or filling it.
        @see: ShapeDrawer.draw_path"""
        height = height or width
        ctx.rectangle(center_x - width/2., center_y - height/2.,
                width, height)

    # pylint: disable-msg=C0103, R0911
    # R0911: too many return statements
    @staticmethod
    def intersection_point(center_x, center_y, source_x, source_y, \
            width, height=None):
        """Determines where the rectangle centered at (center_x, center_y)
        having the given width and height intersects with a line drawn from
        (source_x, source_y) to (center_x, center_y).

        @see: ShapeDrawer.intersection_point"""
        height = height or width
        delta_x, delta_y = center_x-source_x, center_y-source_y

        if delta_x == 0 and delta_y == 0:
            return center_x, center_y

        if delta_y > 0 and delta_x <= delta_y and delta_x >= -delta_y:
            # this is the top edge
            ry = center_y - height/2.
            ratio = (height/2.) / delta_y
            return center_x-ratio*delta_x, ry

        if delta_y < 0 and delta_x <= -delta_y and delta_x >= delta_y:
            # this is the bottom edge
            ry = center_y + height/2.
            ratio = (height/2.) / -delta_y
            return center_x-ratio*delta_x, ry

        if delta_x > 0 and delta_y <= delta_x and delta_y >= -delta_x:
            # this is the left edge
            rx = center_x - width/2.
            ratio = (width/2.) / delta_x
            return rx, center_y-ratio*delta_y

        if delta_x < 0 and delta_y <= -delta_x and delta_y >= delta_x:
            # this is the right edge
            rx = center_x + width/2.
            ratio = (width/2.) / -delta_x
            return rx, center_y-ratio*delta_y

        if delta_x == 0:
            if delta_y > 0:
                return center_x, center_y - height/2.
            return center_x, center_y + height/2.

        if delta_y == 0:
            if delta_x > 0:
                return center_x - width/2., center_y
            return center_x + width/2., center_y


class CircleDrawer(ShapeDrawer):
    """Static class which draws circular vertices"""
    names = "circle circular"

    @staticmethod
    def draw_path(ctx, center_x, center_y, width, height=None):
        """Draws a circular path on the Cairo context without stroking or
        filling it.

        Height is ignored, it is the width that determines the diameter of the circle.

        @see: ShapeDrawer.draw_path"""
        ctx.arc(center_x, center_y, width/2., 0, 2*pi)

    @staticmethod
    def intersection_point(center_x, center_y, source_x, source_y, \
            width, height=None):
        """Determines where the circle centered at (center_x, center_y)
        intersects with a line drawn from (source_x, source_y) to
        (center_x, center_y).

        @see: ShapeDrawer.intersection_point"""
        height = height or width
        angle = atan2(center_y-source_y, center_x-source_x)
        return center_x-width/2. * cos(angle), \
               center_y-height/2.* sin(angle)


class UpTriangleDrawer(ShapeDrawer):
    """Static class which draws upright triangles"""
    names = "triangle triangle-up up-triangle arrow arrow-up up-arrow"

    @staticmethod
    def draw_path(ctx, center_x, center_y, width, height=None):
        """Draws an upright triangle on the Cairo context without stroking or
        filling it.
        
        @see: ShapeDrawer.draw_path"""
        height = height or width
        ctx.move_to(center_x-width/2., center_y+height/2.)
        ctx.line_to(center_x, center_y-height/2.)
        ctx.line_to(center_x+width/2., center_y+height/2.)
        ctx.close_path()

    @staticmethod
    def intersection_point(center_x, center_y, source_x, source_y, \
            width, height=None):
        """Determines where the triangle centered at (center_x, center_y)
        intersects with a line drawn from (source_x, source_y) to
        (center_x, center_y).

        @see: ShapeDrawer.intersection_point"""
        # TODO: finish it properly
        height = height or width
        return center_x, center_y

class DownTriangleDrawer(ShapeDrawer):
    """Static class which draws triangles pointing down"""
    names = "down-triangle triangle-down arrow-down down-arrow"

    @staticmethod
    def draw_path(ctx, center_x, center_y, width, height=None):
        """Draws a triangle on the Cairo context without stroking or
        filling it.
        
        @see: ShapeDrawer.draw_path"""
        height = height or width
        ctx.move_to(center_x-width/2., center_y-height/2.)
        ctx.line_to(center_x, center_y+height/2.)
        ctx.line_to(center_x+width/2., center_y-height/2.)
        ctx.close_path()

    @staticmethod
    def intersection_point(center_x, center_y, source_x, source_y, \
            width, height=None):
        """Determines where the triangle centered at (center_x, center_y)
        intersects with a line drawn from (source_x, source_y) to
        (center_x, center_y).

        @see: ShapeDrawer.intersection_point"""
        # TODO: finish it properly
        height = height or width
        return center_x, center_y

class DiamondDrawer(ShapeDrawer):
    """Static class which draws diamonds (i.e. rhombuses)"""
    names = "diamond rhombus"

    @staticmethod
    def draw_path(ctx, center_x, center_y, width, height=None):
        """Draws a rhombus on the Cairo context without stroking or
        filling it.
        
        @see: ShapeDrawer.draw_path"""
        height = height or width
        ctx.move_to(center_x-width/2., center_y)
        ctx.line_to(center_x, center_y+height/2.)
        ctx.line_to(center_x+width/2., center_y)
        ctx.line_to(center_x, center_y-height/2.)
        ctx.close_path()

    @staticmethod
    def intersection_point(center_x, center_y, source_x, source_y, \
            width, height=None):
        """Determines where the rhombus centered at (center_x, center_y)
        intersects with a line drawn from (source_x, source_y) to
        (center_x, center_y).

        @see: ShapeDrawer.intersection_point"""
        # TODO: finish it properly
        height = height or width
        return center_x, center_y

#####################################################################

class PolygonDrawer(AbstractCairoDrawer):
    """Class that is used to draw polygons.
    
    The corner points of the polygon can be set by the C{points}
    property of the drawer, or passed at construction time. Most
    drawing methods in this class also have an extra C{points}
    argument that can be used to override the set of points in the
    C{points} property."""

    def __init__(self, context, bbox=(1, 1), points = []):
        """Constructs a new polygon drawer that draws on the given
        Cairo context.

        @param  context: the Cairo context to draw on
        @param  bbox:    ignored, leave it at its default value
        @param  points:  the list of corner points
        """
        super(PolygonDrawer, self).__init__(context, bbox)
        self.points = points

    def draw_path(self, points=None, corner_radius=0):
        """Sets up a Cairo path for the outline of a polygon on the given
        Cairo context.

        @param points: the coordinates of the corners of the polygon,
          in clockwise or counter-clockwise order, or C{None} if we are
          about to use the C{points} property of the class.
        @param corner_radius: if zero, an ordinary polygon will be drawn.
          If positive, the corners of the polygon will be rounded with
          the given radius.
        """
        if points is None:
            points = self.points

        self.context.new_path()

        if len(points) < 2:
            # Well, a polygon must have at least two corner points
            return

        ctx = self.context
        if corner_radius <= 0:
            # No rounded corners, this is simple
            ctx.move_to(*points[-1])
            for point in points:
                ctx.line_to(*point)
            return

        # Rounded corners. First, we will take each side of the
        # polygon and find what the corner radius should be on
        # each corner. If the side is longer than 2r (where r is
        # equal to corner_radius), the radius allowed by that side
        # is r; if the side is shorter, the radius is the length
        # of the side / 2. For each corner, the final corner radius
        # is the smaller of the radii on the two sides adjacent to
        # the corner.
        points = [Point(*point) for point in points]
        side_vecs = [v-u for u, v in consecutive_pairs(points, circular=True)]
        half_side_lengths = [side.length() / 2 for side in side_vecs]
        corner_radii = [corner_radius] * len(points)
        for idx in xrange(len(corner_radii)):
            prev_idx = -1 if idx == 0 else idx - 1
            radii = [corner_radius, half_side_lengths[prev_idx],
                     half_side_lengths[idx]]
            corner_radii[idx] = min(radii)

        # Okay, move to the last corner, adjusted by corner_radii[-1]
        # towards the first corner
        ctx.move_to(*(points[-1].towards(points[0], corner_radii[-1])))
        # Now, for each point in points, draw a line towards the
        # corner, stopping before it in a distance of corner_radii[idx],
        # then draw the corner
        u = points[-1]
        for idx, (v, w) in enumerate(consecutive_pairs(points, True)):
            radius = corner_radii[idx]
            ctx.line_to(*v.towards(u, radius))
            aux1 = v.towards(u, radius / 2.)
            aux2 = v.towards(w, radius / 2.)
            ctx.curve_to(aux1.x, aux1.y, aux2.x, aux2.y,
                         *v.towards(w, corner_radii[idx]))
            u = v

    def draw(self, points=None):
        """Draws the polygon using the current stroke of the Cairo context.

        @param points: the coordinates of the corners of the polygon,
          in clockwise or counter-clockwise order, or C{None} if we are
          about to use the C{points} property of the class.
        """
        self.draw_path(points)
        self.context.stroke()

#####################################################################

class ShapeDrawerDirectory(object):
    """Static class that resolves shape names to their corresponding
    shape drawer classes.
        
    Classes that are derived from L{ShapeDrawer} in this module are
    automatically registered by L{ShapeDrawerDirectory} when the module
    is loaded for the first time.
    """

    known_shapes = {}

    @classmethod
    def register(cls, drawer_class):
        """Registers the given shape drawer class under the given names.

        @param drawer_class: the shape drawer class to be registered
        """
        names = drawer_class.names
        if isinstance(names, (str, unicode)):
            names = names.split()

        for name in names:
            cls.known_shapes[name] = drawer_class

    @classmethod
    def register_namespace(cls, namespace):
        """Registers all L{ShapeDrawer} classes in the given namespace

        @param namespace: a Python dict mapping names to Python objects."""
        for name, value in namespace.iteritems():
            if name.startswith("__"):
                continue
            if isinstance(value, type):
                if issubclass(value, ShapeDrawer) and value != ShapeDrawer:
                    cls.register(value)

    @classmethod
    def resolve(cls, shape):
        """Given a shape name, returns the corresponding shape drawer class
        
        @param shape: the name of the shape
        @return: the corresponding shape drawer class

        @raise ValueError: if the shape is unknown
        """
        try:
            return cls.known_shapes[shape]
        except KeyError:
            raise ValueError("unknown shape: %s" % shape)

    @classmethod
    def resolve_default(cls, shape, default=NullDrawer):
        """Given a shape name, returns the corresponding shape drawer class
        or the given default shape drawer if the shape name is unknown.
        
        @param shape: the name of the shape
        @param default: the default shape drawer to return when the shape
          is unknown
        @return: the shape drawer class corresponding to the given name or
          the default shape drawer class if the name is unknown
        """
        return cls.known_shapes.get(shape, default)

ShapeDrawerDirectory.register_namespace(sys.modules[__name__].__dict__)

