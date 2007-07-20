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
import os
import platform
import colors
import time

__all__ = ["BoundingBox", "Plot", "plot"]

try:
    import cairo
except ImportError:
    # No cairo support is installed. Drop a warning and create a fake module
    class FakeModule(object):
        def __getattr__(self, a):
            raise TypeError, "plotting not available"
        def __call__(self, a):
            raise TypeError, "plotting not available"
        def __setattr__(self, k, v):
            raise TypeError, "plotting not available"
    cairo=FakeModule()
    warn("Cairo support not installed, no plotting functionality will be available")

class BoundingBox(object):
    """Class representing a bounding box (a rectangular area)."""

    def __init__(self, *args):
        """Creates a bounding box.

        The corners of the bounding box can be specified by either a tuple
        (four items, two for each corner, respectively), four separate numbers
        (X and Y coordinates for each corner) or two separate numbers (width
        and height, the upper left corner is assumed to be at (0,0))"""
        if len(args) == 1:
            coords = tuple(args[0])[0:4]
        elif len(args) == 4:
            coords = tuple(args)
        elif len(args) == 2:
            coords = (0, 0, args[0], args[1])
        else:
            raise ValueError, "invalid coordinate format"
        try:
            coords = tuple(map(float, coords))
        except ValueError:
            raise ValueError, "invalid coordinate format, numbers expected"
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

    def _get_width(self):
        return self._coords[2]-self._coords[0]
    def _get_height(self):
        return self._coords[3]-self._coords[1]
    def _get_shape(self):
        return self._coords[2]-self._coords[0], self._coords[3]-self._coords[1]
    width = property(_get_width, "Gets the width of the bounding box")
    height = property(_get_height, "Gets the height of the bounding box")
    shape = property(_get_shape, "Gets the shape of the bounding box (width, height)")


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

      * C{cairo.GlitzSurface} -- OpenGL accelerated surface for the X11
        Window System.
      * C{cairo.ImageSurface} -- memory buffer surface. Can be written to a
        C{PNG} image file.
      * C{cairo.PdfSurface} -- PDF document surface.
      * C{cairo.PsSurface} -- PostScript document surface.
      * C{cairo.Win32Surface} -- Microsoft Windows screen rendering.
      * C{cairo.XlibSurface} -- X11 Window System screen rendering.

    If you create a C{Plot} object with a string given as the target surface,
    the string will be treated as a filename, and its extension will decide
    which surface class will be used. Please note that not all surfaces might
    be available, depending on your C{pycairo} installation.

    A C{Plot} has an assigned default palette (see L{colors.Palette}) which
    is used for plottingo bjects.

    A C{Plot} object also has a list of objects to be plotted with their
    respective bounding boxes, palettes and opacities. Palettes assigned
    to an object override the default palette of the plot. Objects can be
    added by the L{Plot.add} method and removed by the L{Plot.remove} method.
    """
    
    def __init__(self, target=None, bbox=None, palette=None):
        """Creates a new plot.

        @param target: the target surface to write to. It can be one of the
          following types:

          * C{None} -- an appropriate surface will be created and the object
            will be plotted there.
          * C{cairo.Surface} -- the given Cairo surface will be used.
          * C{string} -- a file with the given name will be created and an
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
        self._tmpfile_obj=None
        self._filename=None
        self._bgcolor=None

        if bbox is None:
            bbox = BoundingBox(600, 600)

        if palette is None:
            from igraph import config
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
            raise ValueError, "opacity must be between 0.0 and 1.0"
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
            if o == object and (bbox is None or b == bbox):
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
                ctx.push_group()
                plotter(ctx, bbox, palette, *args, **kwds)
                ctx.pop_group_to_source()
                ctx.paint_with_alpha(opacity)

        self._is_dirty = False

    def _create_tmpfile(self):
        from tempfile import NamedTemporaryFile
        self._tmpfile_obj=NamedTemporaryFile(prefix="igraph", suffix=".png")
        return self._tmpfile_obj

    def _close_tmpfile(self):
        if self._tmpfile_obj:
            self._tmpfile_obj.close()
            self._tmpfile_obj=None

    def save(self, fname=None):
        """Saves the plot.

        @param fname: the filename to save to. It is ignored if the surface
          of the plot is not an C{ImageSurface}.
        """
        if self._is_dirty: self.redraw()
        if isinstance(self._surface, cairo.ImageSurface):
            if self._tmpfile: self._create_tmpfile()
            fname = fname or self._filename or self._tmpfile_obj.name
            if fname is None:
                raise ValueError, "no file name is known for the surface and none given"
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
        sur.write_to_png(self._tmpfile_obj.name)
        from igraph import config # Can only be imported here
        imgviewer = config["apps.image_viewer"]
        if not imgviewer:
            # No image viewer was given and none was detected. This
            # should only happen on unknown platforms.
            raise NotImplementedError, "showing plots is not implemented on this platform: %s" % platform.system()
        else:
            os.system("%s %s" % (imgviewer, self._tmpfile_obj.name))
            if platform.system() == "Darwin":
                # On Mac OS X, launched applications are likely to
                # fork and give control back to Python immediately.
                # Chances are that the temporary image file gets removed
                # before the image viewer has a chance to open it, so
                # we wait here a little bit. Yes, this is quite hackish :(
                time.sleep(5)
        self._close_tmpfile()


def plot(obj, target=None, bbox=(0, 0, 600, 600), *args, **kwds):
    """Plots the given object to the given target.

    Positional and keyword arguments not explicitly mentioned here will be
    passed down to the C{__plot__} method of the object being plotted.
    Since you are most likely interested in the keyword arguments available
    for graph plots, see L{Graph.__plot__} as well.

    @param obj: the object to be plotted
    @param target: the target where the object should be plotted. It can be one
      of the following types:
      
        * C{None} -- an appropriate surface will be created and the object will
          be plotted there.
        * C{cairo.Surface} -- the given Cairo surface will be used. This can
          refer to a PNG image, an arbitrary window, an SVG file, anything that
          Cairo can handle.
        * C{string} -- a file with the given name will be created and an
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
    result.add(obj, bbox, *args, **kwds)
    if target is None: result.show()
    if isinstance(target, basestring): result.save()
    return result

def collect_attributes(n, name, alt_name, kwds, vs, config, default, transform=None):
    """Collects graph visualization attributes from various sources.

    This method is used by L{Graph.__plot__} to collect the attributes required
    for graph visualization from various sources. Attribute value sources are:

      * A specific value of a Python dict belonging to a given key. This dict
        is given by the argument M{kwds}, and the name of the key is determined
        by the argument M{name}.

      * A vertex or edge sequence of a graph, given in M{vs}

      * The global configuration, given in M{config}

      * A default value when all other sources fail to provide the value.
        given in M{default}

    Attribute sources are considered exactly in the order mentioned above.
    Optionally, the retrieved value is passed through an arbitrary
    transformation.

    @return: the collected attributes
    """
    try:
        attrs = vs[alt_name]
    except:
        attrs = None

    result = kwds.get(name, None)
    if attrs:
        if not result:
            result = attrs
        else:
            result = [result[idx] or attrs[idx] for idx in xrange(len(result))]

    try:
        m = len(result)
    except TypeError:
        result = [result] * n
    
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

