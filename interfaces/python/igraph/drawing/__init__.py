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

from __future__ import with_statement

from itertools import izip
from operator import itemgetter
from warnings import warn

import math
import os
import platform
import time

from ConfigParser import NoOptionError

from igraph.configuration import Configuration
from igraph.drawing.colors import Palette, palettes
from igraph.drawing.edge import ArrowEdgeDrawer
from igraph.drawing.shapes import ShapeDrawerDirectory
from igraph.drawing.metamagic import AttributeCollectorBase, \
                                     AttributeSpecification
from igraph.utils import named_temporary_file

__all__ = ["BoundingBox", "DefaultGraphDrawer", "Plot", "Point", "plot"]

__license__ = "GPL"

try:
    import cairo
except ImportError:
    # No cairo support is installed. Create a fake module
    # pylint: disable-msg=R0903
    # R0903: too few public methods
    class FakeModule(object):
        """Fake module that raises an exception for everything"""

        def __getattr__(self, _):
            raise TypeError("plotting not available")
        def __call__(self, _):
            raise TypeError("plotting not available")
        def __setattr__(self, key, value):
            raise TypeError("plotting not available")

    # pylint: disable-msg=C0103
    # C0103: invalid name
    cairo = FakeModule()


class Point(tuple):
    """Class representing a point on the 2D plane."""
    __slots__ = ()
    _fields = ('x', 'y')

    def __new__(cls, x, y):
        """Creates a new point with the given coordinates"""
        return tuple.__new__(cls, (x, y))

    # pylint: disable-msg=W0622
    # W0622: redefining built-in 'len'
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

    # pylint: disable-msg=W0141
    # W0141: used builtin function 'map'
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
        ratio = float(ratio)
        return Point(x = self.x * (1.0 - ratio) + other.x * ratio, \
                     y = self.y * (1.0 - ratio) + other.y * ratio)

    def length(self):
        """Returns the length of the vector pointing from the origin to this
        point."""
        return (self.x ** 2 + self.y ** 2) ** 0.5

    def sq_length(self):
        """Returns the squared length of the vector pointing from the origin
        to this point."""
        return (self.x ** 2 + self.y ** 2)


class BoundingBox(object):
    """Class representing a bounding box (a rectangular area)."""

    __slots__ = ("_coords", )

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
            coords = tuple(float(coord) for coord in coords)
        except ValueError:
            raise ValueError("invalid coordinate format, numbers expected")

        self._coords = None    # to make pylint happy
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

    @property
    def width(self):
        """Returns the width of the bounding box"""
        return self._coords[2]-self._coords[0]

    @property
    def height(self):
        """Returns the height of the bounding box"""
        return self._coords[3]-self._coords[1]

    @property
    def left(self):
        """Returns the X coordinate of the left side of the box"""
        return self._coords[0]

    @property
    def right(self):
        """Returns the X coordinate of the right side of the box"""
        return self._coords[2]

    @property
    def top(self):
        """Returns the Y coordinate of the top edge of the box"""
        return self._coords[1]

    @property
    def bottom(self):
        """Returns the Y coordinate of the bottom edge of the box"""
        return self._coords[3]

    @property
    def shape(self):
        """Returns the shape of the bounding box (width, height)"""
        return self._coords[2]-self._coords[0], self._coords[3]-self._coords[1]

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

    def __eq__(self, other):
        return self.coords == other.coords
    def __ne__(self, other):
        return self.coords != other.coords
    def __hash__(self):
        return hash(self.coords)



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

    A C{Plot} has an assigned default palette (see L{igraph.colors.Palette}) which
    is used for plotting objects.

    A C{Plot} object also has a list of objects to be plotted with their
    respective bounding boxes, palettes and opacities. Palettes assigned
    to an object override the default palette of the plot. Objects can be
    added by the L{Plot.add} method and removed by the L{Plot.remove} method.
    """

    # pylint: disable-msg=E1101
    # E1101: Module 'cairo' has no 'foo' member - of course it has! :)
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
          an L{igraph.colors.Palette} object or a string referring to a valid
          key of C{igraph.colors.palettes} (see module L{igraph.colors}) or C{None}.
          In the latter case, the default palette given by the configuration
          key C{plotting.palette} is used.
        """
        self._filename = None
        self._surface_was_created = not isinstance(target, cairo.Surface)
        self._need_tmpfile = False
        self._bgcolor = None

        # Several Windows-specific hacks will be used from now on, thanks
        # to Dale Hunscher for debugging and fixing all that stuff
        self._windows_hacks = "Windows" in platform.platform()

        if bbox is None:
            self.bbox = BoundingBox(600, 600)
        elif isinstance(bbox, tuple) or isinstance(bbox, list):
            self.bbox = BoundingBox(bbox)
        else:
            self.bbox = bbox

        if palette is None:
            config = Configuration.instance()
            palette = config["plotting.palette"]
        if not isinstance(palette, Palette):
            palette = palettes[palette]
        self._palette = palette

        if target is None:
            self._need_tmpfile = True
            self._surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, \
                int(self.bbox.width), int(self.bbox.height))
            self._bgcolor = (1., 1., 1.)
        elif isinstance(target, cairo.Surface):
            self._surface = target
        else:
            self._filename = target
            _, ext = os.path.splitext(target)
            ext = ext.lower()
            if ext == ".pdf":
                self._surface = cairo.PDFSurface(target, self.bbox.width, \
                                                 self.bbox.height)
            elif ext == ".ps":
                self._surface = cairo.PSSurface(target, self.bbox.width, \
                                                self.bbox.height)
            elif ext == ".png":
                self._surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, \
                    int(self.bbox.width), int(self.bbox.height))
                self._bgcolor = (1., 1., 1.)
            elif ext == ".svg":
                self._surface = cairo.SVGSurface(target, self.bbox.width, \
                                                 self.bbox.height)
            else:
                raise ValueError("image format not handled by Cairo: %s" % ext)

        self._ctx = cairo.Context(self._surface)
        self._objects = []
        self._is_dirty = False

    def add(self, obj, bbox=None, palette=None, opacity=1.0, *args, **kwds):
        """Adds an object to the plot.

        Arguments not specified here are stored and passed to the object's
        plotting function when necessary. Since you are most likely interested
        in the arguments acceptable by graphs, see L{Graph.__plot__} for more
        details.

        @param obj: the object to be added
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
        if opacity < 0.0 or opacity > 1.0:
            raise ValueError("opacity must be between 0.0 and 1.0")
        bbox = bbox or self.bbox
        if not isinstance(bbox, BoundingBox):
            bbox = BoundingBox(bbox)
        self._objects.append((obj, bbox, palette, opacity, args, kwds))
        self.mark_dirty()

    def remove(self, obj, bbox=None, idx=1):
        """Removes an object from the plot.

        If the object has been added multiple times and no bounding box
        was specified, it removes the instance which occurs M{idx}th
        in the list of identical instances of the object.

        @param obj: the object to be removed
        @param bbox: optional bounding box specification for the object.
          If given, only objects with exactly this bounding box will be
          considered.
        @param idx: if multiple objects match the specification given by
          M{obj} and M{bbox}, only the M{idx}th occurrence will be removed.
        @return: C{True} if the object has been removed successfully,
          C{False} if the object was not on the plot at all or M{idx}
          was larger than the count of occurrences
        """
        for i in xrange(len(self._objects)):
            current_obj, current_bbox = self._objects[i][0:2]
            if current_obj is obj and (bbox is None or current_bbox == bbox):
                idx -= 1
                if idx == 0:
                    self._objects[i:(i+1)] = []
                    self.mark_dirty()
                    return True
        return False

    def mark_dirty(self):
        """Marks the plot as dirty (should be redrawn)"""
        self._is_dirty = True

    # pylint: disable-msg=W0142
    # W0142: used * or ** magic
    def redraw(self, context=None):
        """Redraws the plot"""
        ctx = context or self._ctx
        if self._bgcolor is not None:
            ctx.set_source_rgba(*self._bgcolor)
            ctx.rectangle(0, 0, self.bbox.width, self.bbox.height)
            ctx.fill()

        for obj, bbox, palette, opacity, args, kwds in self._objects:
            palette = palette or self._palette
            plotter = getattr(obj, "__plot__", None)
            if plotter is None:
                warn("%s does not support plotting" % obj)
            else:
                if opacity < 1.0:
                    ctx.push_group()
                else:
                    ctx.save()
                plotter(ctx, bbox, palette, *args, **kwds)
                if opacity < 1.0:
                    ctx.pop_group_to_source()
                    ctx.paint_with_alpha(opacity)
                else:
                    ctx.restore()

        self._is_dirty = False

    def save(self, fname=None):
        """Saves the plot.

        @param fname: the filename to save to. It is ignored if the surface
          of the plot is not an C{ImageSurface}.
        """
        if self._is_dirty:
            self.redraw()
        if isinstance(self._surface, cairo.ImageSurface):
            if self._need_tmpfile:
                with named_temporary_file(prefix="igraph", suffix=".png") as fname:
                    self._surface.write_to_png(fname)
                    return None

            fname  = fname or self._filename
            if fname is None:
                raise ValueError("no file name is known for the surface " + \
                                 "and none given")
            return self._surface.write_to_png(fname)

        if fname is not None:
            warn("filename is ignored for surfaces other than ImageSurface")

        self._ctx.show_page()
        self._surface.finish()


    def show(self):
        """Saves the plot to a temporary file and shows it."""
        if not isinstance(self._surface, cairo.ImageSurface):
            sur = cairo.ImageSurface(cairo.FORMAT_ARGB32,
                    int(self.bbox.width), int(self.bbox.height))
            ctx = cairo.Context(sur)
            self.redraw(ctx)
        else:
            sur = self._surface
            ctx = self._ctx
            if self._is_dirty:
                self.redraw(ctx)

        with named_temporary_file(prefix="igraph", suffix=".png") as tmpfile:
            sur.write_to_png(tmpfile)
            config = Configuration.instance()
            imgviewer = config["apps.image_viewer"]
            if not imgviewer:
                # No image viewer was given and none was detected. This
                # should only happen on unknown platforms.
                plat = platform.system()
                raise NotImplementedError("showing plots is not implemented " + \
                                          "on this platform: %s" % plat)
            else:
                os.system("%s %s" % (imgviewer, tmpfile))
                if platform.system() == "Darwin" or self._windows_hacks:
                    # On Mac OS X and Windows, launched applications are likely to
                    # fork and give control back to Python immediately.
                    # Chances are that the temporary image file gets removed
                    # before the image viewer has a chance to open it, so
                    # we wait here a little bit. Yes, this is quite hackish :(
                    time.sleep(5)


    @property
    def bounding_box(self):
        """Returns the bounding box of the Cairo surface as a
        L{BoundingBox} object"""
        return BoundingBox(self.bbox)

    @property
    def height(self):
        """Returns the height of the Cairo surface on which the plot
        is drawn"""
        return self.bbox.height

    @property
    def surface(self):
        """Returns the Cairo surface on which the plot is drawn"""
        return self._surface

    @property
    def width(self):
        """Returns the width of the Cairo surface on which the plot
        is drawn"""
        return self.bbox.width

#####################################################################

# pylint: disable-msg=R0903
# R0903: too few public methods
class AbstractDrawer(object):
    """Abstract class that serves as a base class for anything that
    draws on a Cairo context within a given bounding box."""

    def __init__(self, context, bbox):
        """Constructs the drawer and associates it to the given
        Cairo context and the given L{BoundingBox}.

        @param context: the context on which we will draw
        @param bbox:    the bounding box within which we will draw.
                        Can be anything accepted by the constructor
                        of L{BoundingBox} (i.e., a 2-tuple, a 4-tuple
                        or a L{BoundingBox} object).
        """
        self.context = context
        if not isinstance(bbox, BoundingBox):
            bbox = BoundingBox(bbox)
        self.bbox = bbox

    def draw(self, *args, **kwds):
        """Abstract method, must be implemented in derived classes."""
        raise NotImplementedError("abstract class")


#####################################################################

# pylint: disable-msg=R0922
# R0922: Abstract class is only referenced 1 times
class AbstractGraphDrawer(AbstractDrawer):
    """Abstract class that serves as a base class for anything that
    draws a graph on a Cairo context within a given bounding box.
    
    This class is primarily used to collect routines that can be
    potentially useful in different kinds of graph drawers."""

    def draw(self, *args, **kwds):
        """Abstract method, must be implemented in derived classes."""
        raise NotImplementedError("abstract class")


#####################################################################

class DefaultGraphDrawer(AbstractGraphDrawer):
    """Class implementing the default visualisation of a graph.

    The default visualisation of a graph draws the nodes on a 2D plane
    according to a given L{Layout}, then draws a straight or curved
    edge between nodes connected by edges. This is the visualisation
    used when one invokes the L{plot()} function on a L{Graph} object.

    See L{Graph.__plot__()} for the keyword arguments understood by
    this drawer."""

    def __init__(self, context, bbox, \
                 edge_drawer_factory = ArrowEdgeDrawer):
        """Constructs the graph drawer and associates it to the given
        Cairo context and the given L{BoundingBox}.

        @param context: the context on which we will draw
        @param bbox:    the bounding box within which we will draw.
                        Can be anything accepted by the constructor
                        of L{BoundingBox} (i.e., a 2-tuple, a 4-tuple
                        or a L{BoundingBox} object).
        @param edge_drawer_factory: a factory method that returns an
                        L{AbstractEdgeDrawer} instance bound to a
                        given Cairo context. You can use any of the
                        actual L{AbstractEdgeDrawer} implementations
                        here to control the style of edges drawn by
                        igraph. The default edge drawer is
                        L{ArrowEdgeDrawer}.
        """
        AbstractGraphDrawer.__init__(self, context, bbox)
        self.edge_drawer_factory = edge_drawer_factory

    # pylint: disable-msg=W0142,W0221,E1101
    # W0142: Used * or ** magic
    # W0221: argument number differs from overridden method
    # E1101: Module 'cairo' has no 'foo' member - of course it does :)
    def draw(self, graph, palette, *args, **kwds):
        from igraph.layout import Layout

        # Some abbreviations for sake of simplicity
        directed = graph.is_directed()
        context = self.context

        # Calculate/get the layout of the graph
        layout = kwds.get("layout", None)
        if isinstance(layout, Layout):
            layout = Layout(layout.coords)
        elif isinstance(layout, str) or layout is None:
            layout = graph.layout(layout)
        else:
            layout = Layout(layout)

        # Determine the size of the margin on each side
        margin = kwds.get("margin", [0., 0., 0., 0.])
        try:
            margin = list(margin)
        except TypeError:
            margin = [margin]
        while len(margin)<4:
            margin.extend(margin)
        margin = [x + 5. for x in margin[:4]]

        # Contract the drawing area by the margin and fit the layout
        bbox = self.bbox.contract(margin)
        layout.fit_into(bbox, keep_aspect_ratio=False)

        # Construct the visual vertex/edge builders
        class VisualVertexBuilder(AttributeCollectorBase):
            """Collects some visual properties of a vertex for drawing"""
            _kwds_prefix = "vertex_"
            color = ("red", palette.get)
            label = None
            label_angle = -math.pi/2
            label_dist  = 1.6
            label_color = ("black", palette.get)
            label_size  = 14.0
            position = dict(func=layout.__getitem__)
            shape = ("circle", ShapeDrawerDirectory.resolve_default)
            size  = 10.0

        class VisualEdgeBuilder(AttributeCollectorBase):
            """Collects some visual properties of an edge for drawing"""
            _kwds_prefix = "edge_"
            arrow_size = 1.0
            arrow_width = 1.0
            color = ("#444", palette.get)
            width = 1.0

        vertex_builder = VisualVertexBuilder(graph.vs, kwds)
        edge_builder = VisualEdgeBuilder(graph.es, kwds)

        # Draw the edges
        edge_drawer = self.edge_drawer_factory(context)
        if directed:
            drawer_method = edge_drawer.draw_directed_edge
        else:
            drawer_method = edge_drawer.draw_undirected_edge
        for edge, visual_edge in izip(graph.es, edge_builder):
            src, dest = edge.tuple
            src_vertex, dest_vertex = vertex_builder[src], vertex_builder[dest]
            drawer_method(visual_edge, src_vertex, dest_vertex)

        # Draw the vertices
        context.set_line_width(1)
        for vertex, coords in izip(vertex_builder, layout):
            vertex.shape.draw_path(context, \
                    coords[0], coords[1], vertex.size)
            context.set_source_rgba(*vertex.color)
            context.fill_preserve()
            context.set_source_rgb(0., 0., 0.)
            context.stroke()

        # Draw the vertex labels
        context.select_font_face("sans-serif", cairo.FONT_SLANT_NORMAL, \
            cairo.FONT_WEIGHT_BOLD)
        
        for vertex, coords in izip(vertex_builder, layout):
            if vertex.label is None:
                continue
            xb, _, w, h = context.text_extents(vertex.label)[:4]
            cx, cy = coords
            si = math.sin(vertex.label_angle)
            co = math.cos(vertex.label_angle)
            cx += co * vertex.label_dist * vertex.size / 2.
            cy += si * vertex.label_dist * vertex.size / 2.
            cx += (co - 1) * w/2. + xb
            cy += (si + 1) * h/2.
            context.move_to(cx, cy)
            context.set_font_size(vertex.label_size)
            context.set_source_rgba(*vertex.label_color)
            context.text_path(vertex.label)
            context.fill()


#####################################################################

# pylint: disable-msg=R0922
# R0922: Abstract class is only referenced 1 times
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

    def local_to_context(self, x, y):
        """Converts local coordinates to the context coordinate system (given
        by the bounding box).
        
        This method must be overridden in derived classes."""
        raise NotImplementedError("abstract class")


class DescartesCoordinateSystem(CoordinateSystem):
    """Class implementing a 2D Descartes coordinate system object."""

    def __init__(self, context, bbox, bounds):
        """Initializes the coordinate system.

        @param context: the context on which the coordinate system will
          be drawn.
        @param bbox: the bounding box that will contain the coordinate
          system.
        @param bounds: minimum and maximum X and Y values in a 4-tuple.
        """
        self._bounds, self._bbox = None, None
        self._sx, self._sy = None, None
        self._ox, self._oy, self._ox2, self._oy2 = None, None, None, None

        CoordinateSystem.__init__(self, context, bbox)
        self.bbox = bbox
        self.bounds = bounds

    def _get_bbox(self):
        """Returns the bounding box of the coordinate system"""
        return BoundingBox(self._bbox.coords)
    def _set_bbox(self, bbox):
        """Sets the bounding box of the coordinate system"""
        self._bbox = bbox
        self._recalc_scale_factors()
    bbox = property(_get_bbox, _set_bbox)

    def _get_bounds(self):
        """Returns the lower and upper bounds of the X and Y values"""
        return self._bounds.coords
    def _set_bounds(self, bounds):
        """Sets the lower and upper bounds of the X and Y values"""
        self._bounds = BoundingBox(bounds)
        self._recalc_scale_factors()
    bounds = property(_get_bounds, _set_bounds)

    def _recalc_scale_factors(self):
        """Recalculates some cached scale factors used within the class"""
        if self._bounds is None:
            return
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
    if not isinstance(bbox, BoundingBox):
        bbox = BoundingBox(bbox)

    result = Plot(target, bbox)
    contract_w, contract_h = bbox.width/60., bbox.height/60.
    bbox = bbox.contract((contract_w, contract_h, contract_w, contract_h))

    result.add(obj, bbox, *args, **kwds)
    if target is None:
        result.show()

    if isinstance(target, basestring):
        result.save()

    return result

#####################################################################

