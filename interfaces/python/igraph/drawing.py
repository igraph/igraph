"""
Drawing and plotting routines for IGraph.

Plotting is dependent on the C{pycairo} library which provides Python bindings
to the popular U{Cairo library<http://www.cairographics.org>}. This means that
if you don't have U{pycairo<http://www.cairographics.org/pycairo>} installed,
you won't be able to use the plotting capabilities. However, you can still use
L{Graph.write_svg} to save the graph to an SVG file and view it from
U{Mozilla Firefox<http://www.mozilla.org/firefox>} (free) or edit it in
U{Inkscape<http://www.inkscape.org>} (free), U{Skencil<http://www.skencil.org>}
(formerly known as Sketch, also free) or Adobe Illustrator (not free, therefore
I'm not linking to it :)).
"""

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

from warnings import warn
from operator import itemgetter

import math
import os
import platform
import time

import igraph.colors as colors

from igraph.configuration import Configuration

__all__ = ["BoundingBox", "DefaultGraphDrawer", "Plot", "Point", "plot"]

try:
    import cairo
except ImportError:
    # No cairo support is installed. Create a fake module
    class FakeModule(object):
        def __getattr__(self, a):
            raise TypeError("plotting not available")
        def __call__(self, a):
            raise TypeError("plotting not available")
        def __setattr__(self, k, v):
            raise TypeError("plotting not available")
    cairo=FakeModule()


class Point(tuple):
    """Class representing a point on the 2D plane."""
    __slots__ = ()
    _fields = ('x', 'y')

    def __new__(_cls, x, y):
        """Creates a new point with the given coordinates"""
        return tuple.__new__(_cls, (x, y))

    @classmethod
    def _make(cls, iterable, new = tuple.__new__, len = len):
        """Creates a new point from a sequence or iterable"""
        result = new(cls, iterable)
        if len(result) != 2:
            raise TypeError('Expected 2 arguments, got %d' % len(result))
        return result

    def __repr__(self):
        """Returns a nicely formatted representation of the point"""
        return 'Point(x=%r, y=%r)' % self

    def _asdict(self):
        """Returns a new dict which maps field names to their values"""
        return dict(zip(self._fields, self))

    def _replace(self, **kwds):
        """Returns a new point object replacing specified fields with new
        values"""
        result = self._make(map(kwds.pop, ('x', 'y'), self))
        if kwds:
            raise ValueError('Got unexpected field names: %r' % kwds.keys())
        return result

    def __getnewargs__(self):
        """Return self as a plain tuple. Used by copy and pickle."""
        return tuple(self)

    x = property(itemgetter(0), doc="Alias for field number 0")
    y = property(itemgetter(1), doc="Alias for field number 1")

    def __add__(self, other):
        """Adds the coordinates of a point to another one"""
        return self.__class__(x = self.x + other.x, y = self.y + other.y)

    def __sub__(self, other):
        """Subtracts the coordinates of a point to another one"""
        return self.__class__(x = self.x - other.x, y = self.y - other.y)

    def __mul__(self, scalar):
        """Multiplies the coordinates by a scalar"""
        return self.__class__(x = self.x * scalar, y = self.y * scalar)

    def __div__(self, scalar):
        """Divides the coordinates by a scalar"""
        return self.__class__(x = self.x / scalar, y = self.y / scalar)

    def interpolate(self, other, ratio = 0.5):
        """Linearly interpolates between the coordinates of this point and
        another one.

        @param  other:  the other point
        @param  ratio:  the interpolation ratio between 0 and 1. Zero will
          return this point, 1 will return the other point.
        """
        r = float(ratio)
        return Point(x = self.x * (1.0 - r) + other.x * r, \
                     y = self.y * (1.0 - r) + other.y * r)

    def length(self):
        """Returns the length of the vector pointing from the origin to this point."""
        return (self.x ** 2 + self.y ** 2) ** 0.5

    def sq_length(self):
        """Returns the squared length of the vector pointing from the origin to this point."""
        return (self.x ** 2 + self.y ** 2)


class BoundingBox(object):
    """Class representing a bounding box (a rectangular area)."""

    def __init__(self, *args):
        """Creates a bounding box.

        The corners of the bounding box can be specified by either a tuple
        (four items, two for each corner, respectively), four separate numbers
        (X and Y coordinates for each corner) or two separate numbers (width
        and height, the upper left corner is assumed to be at (0,0))"""
        coords = None
        if len(args) == 1:
            if isinstance(args[0], BoundingBox):
                coords = args[0].coords
            elif len(args[0]) >= 4:
                coords = tuple(args[0])[0:4]
            elif len(args[0]) == 2:
                coords = (0, 0, args[0][0], args[0][1])
        elif len(args) == 4:
            coords = tuple(args)
        elif len(args) == 2:
            coords = (0, 0, args[0], args[1])
        if coords is None:
            raise ValueError("invalid coordinate format")
        try:
            coords = tuple(map(float, coords))
        except ValueError:
            raise ValueError("invalid coordinate format, numbers expected")
        self.coords = coords

    def _set_coords(self, coords):
        """Sets the coordinates of the corners.

        @param coords: a 4-tuple with the coordinates of the corners
        """
        if coords[0] > coords[2]:
            coords = (coords[2], coords[1], coords[0], coords[3])
        if coords[1] > coords[3]:
            coords = (coords[0], coords[3], coords[2], coords[1])
        self._coords = coords

    def _get_coords(self):
        """Returns the coordinates of the corners."""
        return self._coords

    coords = property(_get_coords, _set_coords,
        doc="Sets or returns the coordinates of the corners")

    def _get_width(self): return self._coords[2]-self._coords[0]
    def _get_height(self): return self._coords[3]-self._coords[1]
    def _get_left(self): return self._coords[0]
    def _get_right(self): return self._coords[2]
    def _get_top(self): return self._coords[1]
    def _get_bottom(self): return self._coords[3]
    def _get_shape(self):
        return self._coords[2]-self._coords[0], self._coords[3]-self._coords[1]
    width = property(_get_width, doc="Gets the width of the bounding box")
    height = property(_get_height, doc="Gets the height of the bounding box")
    left = property(_get_left, doc="X coordinate of the left side of the box")
    right = property(_get_right, doc="X coordinate of the right side of the box")
    top = property(_get_top, doc="Y coordinate of the top of the box")
    bottom = property(_get_bottom, doc="Y coordinate of the bottom of the box")
    shape = property(_get_shape, doc="Gets the shape of the bounding box (width, height)")

    def contract(self, margins):
        """Contracts the bounding box by the given margins.

        @return: a new L{BoundingBox} object.
        """
        if isinstance(margins, int) or isinstance(margins, float):
            margins = [float(margins)] * 4
        if len(margins) != 4:
            raise ValueError("margins must be a 4-tuple or a single number")
        nx1, ny1 = self._coords[0]+margins[0], self._coords[1]+margins[1]
        nx2, ny2 = self._coords[2]-margins[2], self._coords[3]-margins[3]
        if nx1 > nx2:
            nx1 = (nx1+nx2)/2.
            nx2 = nx1
        if ny1 > ny2:
            ny1 = (ny1+ny2)/2.
            ny2 = ny1
        return BoundingBox(nx1, ny1, nx2, ny2)

    def __repr__(self):
        return "%s(%s, %s, %s, %s)" % (self.__class__.__name__, \
            self._coords[0], self._coords[1], self._coords[2], \
            self._coords[3])

    def __eq__(self, other): return self.coords == other.coords
    def __ne__(self, other): return self.coords != other.coords
    def __hash__(self): return hash(self.coords)



class Plot(object):
    """Class representing an arbitrary plot

    Every plot has an associated surface object where the plotting is done. The
    surface is an instance of C{cairo.Surface}, a member of the C{pycairo}
    library. The surface itself provides a unified API to various plotting
    targets like SVG files, X11 windows, PostScript files, PNG files and so on.
    C{igraph} usually does not know on which surface it is plotting right now,
    since C{pycairo} takes care of the actual drawing. Everything that's supported
    by C{pycairo} should be supported by this class as well.

    Current Cairo surfaces that I'm aware of are:

      - C{cairo.GlitzSurface} -- OpenGL accelerated surface for the X11
        Window System.

      - C{cairo.ImageSurface} -- memory buffer surface. Can be written to a
        C{PNG} image file.

      - C{cairo.PdfSurface} -- PDF document surface.

      - C{cairo.PsSurface} -- PostScript document surface.

      - C{cairo.Win32Surface} -- Microsoft Windows screen rendering.

      - C{cairo.XlibSurface} -- X11 Window System screen rendering.

    If you create a C{Plot} object with a string given as the target surface,
    the string will be treated as a filename, and its extension will decide
    which surface class will be used. Please note that not all surfaces might
    be available, depending on your C{pycairo} installation.

    A C{Plot} has an assigned default palette (see L{colors.Palette}) which
    is used for plotting objects.

    A C{Plot} object also has a list of objects to be plotted with their
    respective bounding boxes, palettes and opacities. Palettes assigned
    to an object override the default palette of the plot. Objects can be
    added by the L{Plot.add} method and removed by the L{Plot.remove} method.
    """
    
    def __init__(self, target=None, bbox=None, palette=None):
        """Creates a new plot.

        @param target: the target surface to write to. It can be one of the
          following types:

            - C{None} -- an appropriate surface will be created and the object
              will be plotted there.

            - C{cairo.Surface} -- the given Cairo surface will be used.

            - C{string} -- a file with the given name will be created and an
              appropriate Cairo surface will be attached to it.

        @param bbox: the bounding box of the surface. It is interpreted
          differently with different surfaces: PDF and PS surfaces will
          treat it as points (1 point = 1/72 inch). Image surfaces will
          treat it as pixels. SVG surfaces will treat it as an abstract
          unit, but it will mostly be interpreted as pixels when viewing
          the SVG file in Firefox.

        @param palette: the palette primarily used on the plot if the
          added objects do not specify a private palette. Must be either
          a L{colors.Palette} object or a string referring to a valid
          key of C{colors.palettes} (see module L{colors}) or C{None}.
          In the latter case, the default palette given by the configuration
          key C{plotting.palette} is used.
        """
        self._filename=None
        self._surface_was_created=not isinstance(target, cairo.Surface)
        self._tmpfile=False
        self._tmpfile_name=None
        self._filename=None
        self._bgcolor=None

        # Several Windows-specific hacks will be used from now on, thanks
        # to Dale Hunscher for debugging and fixing all that stuff
        self._windows_hacks=("Windows" in platform.platform())

        if bbox is None:
            bbox = BoundingBox(600, 600)
        elif isinstance(bbox, tuple) or isinstance(bbox, list):
            bbox = BoundingBox(bbox)

        if palette is None:
            config = Configuration.instance()
            palette = config["plotting.palette"]
        if not isinstance(palette, colors.Palette):
            palette = colors.palettes[palette]
        self._palette = palette

        if target is None:
            self._tmpfile=True
            self._surface=cairo.ImageSurface(cairo.FORMAT_ARGB32, \
                int(bbox.width), int(bbox.height))
            self._bgcolor=(1., 1., 1.)
        elif isinstance(target, cairo.Surface):
            self._surface = target
        else:
            import os.path
            self._filename = target
            fname, ext = os.path.splitext(target)
            ext=ext.lower()
            if ext == ".pdf":
                self._surface = cairo.PDFSurface(target, bbox.width, bbox.height)
            elif ext == ".ps":
                self._surface = cairo.PSSurface(target, bbox.width, bbox.height)
            elif ext == ".png":
                self._surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, \
                    int(bbox.width), int(bbox.height))
                self._bgcolor = (1., 1., 1.)
            elif ext == ".svg":
                self._surface = cairo.SVGSurface(target, bbox.width, bbox.height)

        self._width = bbox.width
        self._height = bbox.height
        self._ctx = cairo.Context(self._surface)
        self._objects = []
        self._is_dirty = False

    def add(self, object, bbox=None, palette=None, opacity=1.0, *args, **kwds):
        """Adds an object to the plot.

        Arguments not specified here are stored and passed to the object's
        plotting function when necessary. Since you are most likely interested
        in the arguments acceptable by graphs, see L{Graph.__plot__} for more
        details.

        @param object: the object to be added
        @param bbox: the bounding box of the object. If C{None}, the object
          will fill the entire area of the plot.
        @param palette: the color palette used for drawing the object. If the
          object tries to get a color assigned to a positive integer, it
          will use this palette. If C{None}, defaults to the global palette
          of the plot.
        @param opacity: the opacity of the object being plotted, in the range
          0.0-1.0

        @see: Graph.__plot__
        """
        if opacity<0.0 or opacity>1.0:
            raise ValueError("opacity must be between 0.0 and 1.0")
        bbox = bbox or (0,0,self._width,self._height)
        if not isinstance(bbox, BoundingBox): bbox = BoundingBox(bbox)
        self._objects.append((object, bbox, palette, opacity, args, kwds))
        self.mark_dirty()

    def remove(self, object, bbox=None, idx=1):
        """Removes an object from the plot.

        If the object has been added multiple times and no bounding box
        was specified, it removes the instance which occurs M{idx}th
        in the list of identical instances of the object.

        @param object: the object to be removed
        @param bbox: optional bounding box specification for the object.
          If given, only objects with exactly this bounding box will be
          considered.
        @param idx: if multiple objects match the specification given by
          M{object} and M{bbox}, only the M{idx}th occurrence will be removed.
        @return: C{True} if the object has been removed successfully,
          C{False} if the object was not on the plot at all or M{idx}
          was larger than the count of occurrences
        """
        for i, (o, b, _, _, _, _) in enumerate(self._objects):
            if o is object and (bbox is None or b == bbox):
                idx -= 1
                if idx == 0:
                    self._objects[i:(i+1)] = []
                    self.mark_dirty()
                    return True
        return False

    def mark_dirty(self):
        """Marks the plot as dirty (should be redrawn)"""
        self._is_dirty = True

    def redraw(self, context=None):
        """Redraws the plot"""
        ctx = context or self._ctx
        if self._bgcolor is not None:
            ctx.set_source_rgb(*self._bgcolor)
            ctx.rectangle(0, 0, self._width, self._height)
            ctx.fill()

        for obj, bbox, palette, opacity, args, kwds in self._objects:
            palette = palette or self._palette
            plotter = getattr(obj, "__plot__", None)
            if plotter is None:
                warn("%s does not support plotting" % obj)
            else:
                if opacity < 1.0: ctx.push_group()
                else: ctx.save()
                plotter(ctx, bbox, palette, *args, **kwds)
                if opacity < 1.0:
                    ctx.pop_group_to_source()
                    ctx.paint_with_alpha(opacity)
                else:
                    ctx.restore()

        self._is_dirty = False

    def _create_tmpfile(self):
        from tempfile import mkstemp
        handle, self._tmpfile_name = mkstemp(prefix="igraph", suffix=".png")
        os.close(handle)
        return self._tmpfile_name

    def _close_tmpfile(self):
        if self._tmpfile_name:
            os.unlink(self._tmpfile_name)
            self._tmpfile_name = None

    def save(self, fname=None):
        """Saves the plot.

        @param fname: the filename to save to. It is ignored if the surface
          of the plot is not an C{ImageSurface}.
        """
        if self._is_dirty: self.redraw()
        if isinstance(self._surface, cairo.ImageSurface):
            if self._tmpfile: self._create_tmpfile()
            fname = fname or self._filename or self._tmpfile_name
            if fname is None:
                raise ValueError("no file name is known for the surface and none given")
            result = self._surface.write_to_png(fname)
            if self._tmpfile: self._close_tmpfile()
            if not self._tmpfile: return result
        else:
            if fname is not None:
                warn("filename is ignored for surfaces other than ImageSurface")

            self._ctx.show_page()
            self._surface.finish()

    def show(self):
        """Saves the plot to a temporary file and shows it."""
        if not isinstance(self._surface, cairo.ImageSurface):
            sur=cairo.ImageSurface(cairo.FORMAT_ARGB32,int(self._width),
                int(self._height))
            ctx=cairo.Context(sur)
            self.redraw(ctx)
        else:
            sur=self._surface
            ctx=self._ctx
            if self._is_dirty: self.redraw(ctx)

        self._create_tmpfile()
        sur.write_to_png(self._tmpfile_name)

        config = Configuration.instance()
        imgviewer = config["apps.image_viewer"]
        if not imgviewer:
            # No image viewer was given and none was detected. This
            # should only happen on unknown platforms.
            raise NotImplementedError("showing plots is not implemented on this platform: %s" % platform.system())
        else:
            os.system("%s %s" % (imgviewer, self._tmpfile_name))
            if platform.system() == "Darwin" or self._windows_hacks:
                # On Mac OS X and Windows, launched applications are likely to
                # fork and give control back to Python immediately.
                # Chances are that the temporary image file gets removed
                # before the image viewer has a chance to open it, so
                # we wait here a little bit. Yes, this is quite hackish :(
                time.sleep(5)
        self._close_tmpfile()

#####################################################################

class ShapeDrawer(object):
    """Static class, the ancestor of all vertex shape drawer classes.
    
    Custom shapes must implement at least the C{draw_path} method of the class.
    The method I{must not} stroke or fill, it should just set up the current
    Cairo path appropriately."""
    
    def draw_path(ctx, cx, cy, w, h=None):
        """Draws the path of the shape on the given Cairo context, without stroking
        or filling it.

        This method must be overridden in derived classes implementing custom shapes
        and declared as a static method using C{staticmethod(...)}.

        @param ctx: the context to draw on
        @param cx: the X coordinate of the center of the object
        @param cy: the Y coordinate of the center of the object
        @param w: the width of the object
        @param h: the height of the object. If C{None}, equals to the width.
        """
        raise TypeError("abstract class")
    draw_path=staticmethod(draw_path)


    def intersection_point(cx, cy, sx, sy, w, h=None):
        """Determines where the shape centered at (cx, cy) intersects with a
        line drawn from (sx, sy) to (cx, cy).

        Can be overridden in derived classes. Must always be defined as a static
        method using C{staticmethod(...)}

        @param w: the width of the shape
        @param h: the height of the shape. If C{None}, defaults to the width
        @return: the intersection point (the closest to (sx, sy) if there are
            more than one) or (cx, cy) if there is no intersection
        """
        return cx, cy
    intersection_point=staticmethod(intersection_point)

class NullDrawer(ShapeDrawer):
    """Static drawer class which draws nothing.

    This class is used for graph vertices with unknown shapes"""
    def draw_path(ctx, cx, cy, w, h=None): pass
    draw_path=staticmethod(draw_path)


class RectangleDrawer(ShapeDrawer):
    """Static class which draws rectangular vertices"""

    def draw_path(ctx, cx, cy, w, h=None):
        """Draws a rectangle-shaped path on the Cairo context without stroking or
        filling it.
        @see: ShapeDrawer.draw_path"""
        h = h or w
        ctx.rectangle(cx-w/2., cy-h/2., w, h)
    draw_path=staticmethod(draw_path)

    def intersection_point(cx, cy, sx, sy, w, h=None):
        h = h or w
        dx, dy = cx-sx, cy-sy
        if dx == 0 and dy == 0: return cx, cy
        if dy>0 and dx<=dy and dx>=-dy: # top edge
            ry = cy - h/2.
            ratio = (h/2.) / dy
            return cx-ratio*dx, ry
        if dy<0 and dx<=-dy and dx>=dy: # bottom edge
            ry = cy + h/2.
            ratio = (h/2.) / -dy
            return cx-ratio*dx, ry
        if dx>0 and dy<=dx and dy>=-dx: # left edge
            rx = cx - w/2.
            ratio = (w/2.) / dx
            return rx, cy-ratio*dy
        if dx<0 and dy<=-dx and dy>=dx: # right edge
            rx = cx + w/2.
            ratio = (w/2.) / -dx
            return rx, cy-ratio*dy
        if dx == 0:
            if dy>0: return cx, cy - h/2.
            return cx, cy + h/2.
        if dy == 0:
            if dx>0: return cx - w/2., cy
            return cx + w/2., cy
    intersection_point=staticmethod(intersection_point)


class CircleDrawer(ShapeDrawer):
    """Static class which draws circular vertices"""
    
    def draw_path(ctx, cx, cy, w, h=None):
        """Draws a circular path on the Cairo context without stroking or filling.

        Height is ignored, it is the width that determines the diameter of the circle.

        @see: ShapeDrawer.draw_path"""
        ctx.arc(cx, cy, w/2., 0, 2*math.pi)
    draw_path=staticmethod(draw_path)

    def intersection_point(cx, cy, sx, sy, w, h=None):
        h = h or w
        angle = math.atan2(cy-sy, cx-sx)
        return cx-w/2.*math.cos(angle), cy-h/2.*math.sin(angle)
    intersection_point=staticmethod(intersection_point)


class UpTriangleDrawer(ShapeDrawer):
    """Static class which draws upright triangles"""
    def draw_path(ctx, cx, cy, w, h=None):
        """Draws an upright triangle on the Cairo context without stroking or filling."""
        h = h or w
        ctx.move_to(cx-w/2., cy+h/2.)
        ctx.line_to(cx, cy-h/2.)
        ctx.line_to(cx+w/2., cy+h/2.)
        ctx.line_to(cx-w/2., cy+h/2.)
    draw_path=staticmethod(draw_path)

    def intersection_point(cx, cy, sx, sy, w, h=None):
        # TODO: finish it properly
        h = h or w
        return cx, cy
    intersection_point=staticmethod(intersection_point)

class DownTriangleDrawer(ShapeDrawer):
    """Static class which draws triangles pointing down"""
    def draw_path(ctx, cx, cy, w, h=None):
        """Draws a triangle on the Cairo context without stroking or filling."""
        h = h or w
        ctx.move_to(cx-w/2., cy-h/2.)
        ctx.line_to(cx, cy+h/2.)
        ctx.line_to(cx+w/2., cy-h/2.)
        ctx.line_to(cx-w/2., cy-h/2.)
    draw_path=staticmethod(draw_path)

    def intersection_point(cx, cy, sx, sy, w, h=None):
        # TODO: finish it properly
        h = h or w
        return cx, cy
    intersection_point=staticmethod(intersection_point)

def draw_shape_path(shape, ctx, cx, cy, w, h=None):
    """Draws a path of a shape on the given Cairo context.

    @param shape: the shape to be drawn
    @param ctx: the context to draw on
    @param cx: X coordinate of the center of the shape
    @param cy: Y coordinate of the center of the shape
    @param w: desired width of the shape
    @param h: desired height of the shape. If omitted, defaults to the width.
    """
    global known_shapes
    try:
        drawer = known_shapes[shape]
    except:
        raise ValueError("unknown shape: %s" % shape)
    drawer.draw_path(ctx, cx, cy, w, h)

known_shapes = {
    "rectangle": RectangleDrawer,
    "rect": RectangleDrawer,
    "rectangular": RectangleDrawer,
    "square": RectangleDrawer,
    "box": RectangleDrawer,

    "circle": CircleDrawer,
    "circular": CircleDrawer,

    "null": NullDrawer,
    "": NullDrawer,
    "empty": NullDrawer,
    "hidden": NullDrawer,

    "triangle": UpTriangleDrawer,
    "triangle-up": UpTriangleDrawer,
    "arrow-up": UpTriangleDrawer,
    "up-arrow": UpTriangleDrawer,
    "triangle-down": DownTriangleDrawer,
    "arrow-down": DownTriangleDrawer,
    "down-arrow": DownTriangleDrawer,
}

#####################################################################

class AbstractDrawer(object):
    """Abstract class that serves as a base class for anything that
    draws on a Cairo context within a given bounding box."""

    def __init__(self, context, bbox):
        """Constructs the drawer and associates it to the given
        Cairo context and the given L{BoundingBox}.

        @param context  the context on which we will draw
        @param bbox     the bounding box within which we will draw.
                        Can be anything accepted by the constructor
                        of L{BoundingBox} (i.e., a 2-tuple, a 4-tuple
                        or a L{BoundingBox} object).
        """
        self.context = context
        if not isinstance(bbox, BoundingBox):
            bbox = BoundingBox(bbox)
        self.bbox = bbox

    def draw(self):
        """Abstract method, must be implemented in derived classes."""
        raise NotImplementedError("abstract class")


#####################################################################

class AbstractGraphDrawer(AbstractDrawer):
    """Abstract class that serves as a base class for anything that
    draws a graph on a Cairo context within a given bounding box.
    
    This class is primarily used to collect routines that can be
    potentially useful in different kinds of graph drawers."""

    def collect_attributes(self, name, alt_name, kwds, vs, default, transform=None, config=None):
        """Collects graph visualization attributes from various sources.

        This method can be used to collect the attributes required for graph
        visualization from various sources. Attribute value sources are:

          - A specific value of a Python dict belonging to a given key. This dict
            is given by the argument M{kwds}, and the name of the key is determined
            by the argument M{name}.

          - A vertex or edge sequence of a graph, given in M{vs}

          - The global configuration, given in M{config}

          - A default value when all other sources fail to provide the value.
            given in M{default}

        Attribute sources are considered exactly in the order mentioned above.
        Optionally, the retrieved value is passed through an arbitrary
        transformation.

        @param  name:      the name of the attribute when it is coming from a
                           list of Python keyword arguments
        @param  alt_name:  the name of the attribute when it is coming from the
                           graph attributes directly
        @param  kwds:      a Python dict of keyword arguments that will be
                           indexed by C{name}
        @param  vs:        a L{VertexSeq} or L{EdgeSeq} that will be indexed
                           by C{alt_name}
        @param  default:   the default value of the attribute
        @param  transform: optional callable transformation to call on the values.
                           This can be used to ensure that the attributes are of
                           a given type.
        @param  config:    a L{Configuration} object to be used for determining the
                           defaults if all else fails. If C{None}, the global
                           igraph configuration will be used
        @return: the collected attributes
        """
        n = len(vs)

        if config is None:
            config = Configuration.instance()

        try:
            attrs = vs[alt_name]
        except:
            attrs = None

        result = kwds.get(name, None)
        if attrs:
            if not result:
                result = attrs
            else:
                if isinstance(result, str): result = [result] * n
                try:
                    len(result)
                except TypeError:
                    result = [result] * n
                result = [result[idx] or attrs[idx] for idx in xrange(len(result))]

        if isinstance(result, str): result = [result] * n
        try:
            m = len(result)
        except TypeError:
            result = [result] * n

        if not hasattr(result, "extend"): result = list(result)
        m = len(result)
        while len(result) < n:
            if len(result) <= n/2:
                result.extend(result)
            else:
                result.extend(result[0:(n-len(result))])

        # By now, the length of the result vector should be n as requested
        try:
            conf_def = config["plotting.%s" % name]
        except:
            conf_def = None

        if conf_def and None in result:
            result = [result[idx] or conf_def for idx in xrange(len(result))]

        if None in result:
            result = [result[idx] or default for idx in xrange(len(result))]

        if transform is not None:
            result = [transform(x) for x in result]

        return result


#####################################################################

class DefaultGraphDrawer(AbstractGraphDrawer):
    """Class implementing the default visualisation of a graph.

    The default visualisation of a graph draws the nodes on a 2D plane
    according to a given L{Layout}, then draws a straight or curved
    edge between nodes connected by edges. This is the visualisation
    used when one invokes the L{plot()} function on a L{Graph} object.

    See L{Graph.__plot__()} for the keyword arguments understood by
    this drawer."""

    def __init__(self, context, bbox):
        """Constructs the graph drawer and associates it to the given
        Cairo context and the given L{BoundingBox}.

        @param context  the context on which we will draw
        @param bbox     the bounding box within which we will draw.
                        Can be anything accepted by the constructor
                        of L{BoundingBox} (i.e., a 2-tuple, a 4-tuple
                        or a L{BoundingBox} object).
        """
        AbstractGraphDrawer.__init__(self, context, bbox)

    def draw(self, graph, palette, *args, **kwds):
        from igraph.layout import Layout

        vcount, ecount, directed = graph.vcount(), graph.ecount(), graph.is_directed()
        context = self.context

        margin = kwds.get("margin", [0., 0., 0., 0.])
        try:
            margin = list(margin)
        except TypeError:
            margin = [margin]
        while len(margin)<4: margin.extend(margin)
        margin = tuple(map(float, margin[:4]))

        vertex_colors = self.collect_attributes("vertex_color", \
            "color", kwds, graph.vs, "red", palette.get)
        vertex_sizes = self.collect_attributes("vertex_size", \
            "size", kwds, graph.vs, 10, float)
        vertex_shapes = [known_shapes.get(x, NullDrawer) \
            for x in self.collect_attributes("vertex_shape", \
            "shape", kwds, graph.vs, "circle")]

        max_vertex_size = max(vertex_sizes)

        layout = kwds.get("layout", None)
        if isinstance(layout, Layout):
            layout = Layout(layout.coords)
        elif isinstance(layout, str) or layout is None:
            layout = graph.layout(layout)
        else:
            layout = Layout(layout)

        margin = [x + max_vertex_size/2. for x in margin]
        bbox = self.bbox.contract(margin)
        layout.fit_into(bbox, keep_aspect_ratio=False)

        context.set_line_width(1)

        edge_colors = self.collect_attributes("edge_color", \
            "color", kwds, graph.es, "black", palette.get)
        edge_widths = self.collect_attributes("edge_width", \
            "width", kwds, graph.es, 1, float)
        edge_arrow_sizes = self.collect_attributes( \
            "edge_arrow_size", "arrow_size", kwds, graph.es, 1, float)
        edge_arrow_widths = self.collect_attributes( \
            "edge_arrow_width", "arrow_width", kwds, graph.es, 1, float)

        # Draw the edges
        for idx, e in enumerate(graph.es):
            context.set_source_rgb(*edge_colors[idx])
            context.set_line_width(edge_widths[idx])

            src, tgt = e.tuple
            if src == tgt:
                # Loop edge
                r = vertex_sizes[src]*2
                cx, cy = layout[src][0]+math.cos(math.pi/4)*r/2, \
                  layout[src][1]-math.sin(math.pi/4)*r/2
                context.arc(cx, cy, r/2., 0, math.pi*2)
            else:
                # Determine where the edge intersects the circumference of the
                # vertex shape. TODO: theoretically this need not to be done
                # if there are no arrowheads on the edge, but maybe it's not worth
                # testing for
                p1 = vertex_shapes[src].intersection_point( \
                    layout[src][0], layout[src][1], layout[tgt][0], layout[tgt][1],
                    vertex_sizes[src])
                p2 = vertex_shapes[tgt].intersection_point( \
                    layout[tgt][0], layout[tgt][1], layout[src][0], layout[src][1],
                    vertex_sizes[tgt])
                context.move_to(*p1)
                context.line_to(*p2)
            context.stroke()

            if directed and src != tgt:
                # Draw an arrowhead
                angle = math.atan2(p2[1]-p1[1], p2[0]-p1[0])
                arrow_size = 15.*edge_arrow_sizes[idx]
                arrow_width = 10./edge_arrow_widths[idx]
                a1 = (p2[0]-arrow_size*math.cos(angle-math.pi/arrow_width),
                  p2[1]-arrow_size*math.sin(angle-math.pi/arrow_width))
                a2 = (p2[0]-arrow_size*math.cos(angle+math.pi/arrow_width),
                  p2[1]-arrow_size*math.sin(angle+math.pi/arrow_width))
                context.move_to(*p2)
                context.line_to(*a1)
                context.line_to(*a2)
                context.line_to(*p2)
                context.fill()

        del edge_colors
        del edge_widths

        # Draw the vertices
        context.set_line_width(1)
        for idx in xrange(vcount):
            vertex_shapes[idx].draw_path(context, layout[idx][0], layout[idx][1], vertex_sizes[idx])
            context.set_source_rgb(*vertex_colors[idx])
            context.fill_preserve()
            context.set_source_rgb(0., 0., 0.)
            context.stroke()
        del vertex_colors
        del vertex_shapes

        # Draw the vertex labels
        if "vertex_label" not in kwds and "label" not in graph.vs.attribute_names():
            vertex_labels = [str(i) for i in xrange(vcount)]
        elif "vertex_label" in kwds and kwds["vertex_label"] is None:
            vertex_labels = [""] * vcount
        else:
            vertex_labels = self.collect_attributes("vertex_label", \
                "label", kwds, graph.vs, None)
        vertex_dists = self.collect_attributes("vertex_label_dist", \
            "label_dist", kwds, graph.vs, 1.6, float)
        vertex_degrees = self.collect_attributes(\
            "vertex_label_angle", "label_angle", kwds, graph.vs, \
            -math.pi/2, float)
        vertex_label_colors = self.collect_attributes(\
            "vertex_label_color", "label_color", kwds, graph.vs, \
            "black", palette.get)
        vertex_label_sizes = self.collect_attributes(\
            "vertex_label_size", "label_size", kwds, graph.vs, \
            14, float)

        context.select_font_face("sans-serif", cairo.FONT_SLANT_NORMAL, \
            cairo.FONT_WEIGHT_BOLD)
        
        for idx in xrange(vcount):
            xb, yb, w, h = context.text_extents(vertex_labels[idx])[:4]
            cx, cy = layout[idx]
            si, co = math.sin(vertex_degrees[idx]), math.cos(vertex_degrees[idx])
            cx += co * vertex_dists[idx] * vertex_sizes[idx] / 2.
            cy += si * vertex_dists[idx] * vertex_sizes[idx] / 2.
            cx += (co - 1) * w/2. + xb
            cy += (si + 1) * h/2.
            context.move_to(cx, cy)
            context.set_font_size(vertex_label_sizes[idx])
            context.set_source_rgb(*vertex_label_colors[idx])
            context.text_path(vertex_labels[idx])
            context.fill()


#####################################################################

class CoordinateSystem(AbstractDrawer):
    """Class implementing a coordinate system object.

    Coordinate system objects are used when drawing plots which
    2D or 3D coordinate system axes. This is an abstract class
    which must be extended in order to use it. In general, you'll
    only need the documentation of this class if you intend to
    implement an own coordinate system not present in igraph yet.
    """

    def __init__(self, context, bbox):
        """Initializes the coordinate system.

        @param context: the context on which the coordinate system will
          be drawn.
        @param bbox: the bounding box that will contain the coordinate
          system.
        """
        AbstractDrawer.__init__(self, context, bbox)

    def draw(self):
        """Draws the coordinate system.

        This method must be overridden in derived classes.
        """
        raise NotImplementedError("abstract class")

    def local_to_context(self, *args):
        """Converts local coordinates to the context coordinate system (given
        by the bounding box).
        
        This method must be overridden in derived classes."""
        raise NotImplementedError("abstract class")


class DescartesCoordinateSystem(object):
    """Class implementing a 2D Descartes coordinate system object."""

    def __init__(self, context, bbox, bounds):
        """Initializes the coordinate system.

        @param context: the context on which the coordinate system will
          be drawn.
        @param bbox: the bounding box that will contain the coordinate
          system.
        @param bounds: minimum and maximum X and Y values in a 4-tuple.
        """
        self.context = context
        self._bounds = None
        self.bbox = bbox
        self.bounds = bounds
        self._recalc_scale_factors()

    def _get_bbox(self): return BoundingBox(self._bbox.coords)
    def _set_bbox(self, bbox):
        self._bbox = bbox
        self._recalc_scale_factors()
    bbox = property(_get_bbox, _set_bbox, doc="The bounding box of the coordinate system")

    def _get_bounds(self): return self._bounds.coords
    def _set_bounds(self, bounds):
        self._bounds = BoundingBox(bounds)
        self._recalc_scale_factors()
    bounds = property(_get_bounds, _set_bounds, doc="The lower and upper bounds of the X and Y values")

    def _recalc_scale_factors(self):
        if self._bounds is None: return
        self._sx = self._bbox.width / self._bounds.width
        self._sy = self._bbox.height / self._bounds.height
        self._ox = self._bounds.left
        self._oy = self._bounds.top
        self._ox2 = self._bbox.left
        self._oy2 = self._bbox.bottom

    def draw(self):
        """Draws the coordinate system."""
        # Draw the frame
        coords = self.bbox.coords
        self.context.set_source_rgb(0., 0., 0.)
        self.context.set_line_width(1)
        self.context.rectangle(coords[0], coords[1], \
            coords[2]-coords[0], coords[3]-coords[1])
        self.context.stroke()

    def local_to_context(self, x, y):
        """Converts local coordinates to the context coordinate system (given
        by the bounding box).
        """
        return (x-self._ox)*self._sx+self._ox2, self._oy2-(y-self._oy)*self._sy

#####################################################################

def plot(obj, target=None, bbox=(0, 0, 600, 600), *args, **kwds):
    """Plots the given object to the given target.

    Positional and keyword arguments not explicitly mentioned here will be
    passed down to the C{__plot__} method of the object being plotted.
    Since you are most likely interested in the keyword arguments available
    for graph plots, see L{Graph.__plot__} as well.

    @param obj: the object to be plotted
    @param target: the target where the object should be plotted. It can be one
      of the following types:
      
        - C{None} -- an appropriate surface will be created and the object will
          be plotted there.

        - C{cairo.Surface} -- the given Cairo surface will be used. This can
          refer to a PNG image, an arbitrary window, an SVG file, anything that
          Cairo can handle.

        - C{string} -- a file with the given name will be created and an
          appropriate Cairo surface will be attached to it.
          
    @param bbox: the bounding box of the plot. It must be a tuple with four
      integers, the first two denoting the X and Y coordinates of a corner
      and the latter two denoting the X and Y coordinates of the opposite
      corner. It can also be a L{BoundingBox} object.

    @keyword opacity: the opacity of the object being plotted. It can be
      used to overlap several plots of the same graph if you use the same
      layout for them -- for instance, you might plot a graph with opacity
      0.5 and then plot its spanning tree over it with opacity 0.1. To
      achieve this, you'll need to modify the L{Plot} object returned with
      L{Plot.add}.
    @return: an appropriate L{Plot} object.

    @see: Graph.__plot__
    """
    if not isinstance(bbox, BoundingBox): bbox=BoundingBox(bbox)
    result = Plot(target, bbox)
    cw, ch = bbox.width/60., bbox.height/60.
    result.add(obj, bbox.contract((cw, ch, cw, ch)), *args, **kwds)
    if target is None: result.show()
    if isinstance(target, basestring): result.save()
    return result

#####################################################################

