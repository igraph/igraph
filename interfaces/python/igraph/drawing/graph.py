"""
Drawing routines to draw graphs.

This module contains routines to draw graphs on:

- Cairo surfaces (L{DefaultGraphDrawer})
- UbiGraph displays (L{UbiGraphDrawer}, see L{http://ubietylab.net/ubigraph})

It also contains routines to send an igraph graph directly to Cytoscape
(L{http://www.cytoscape.org}) using the CytoTalk plugin.
"""

from itertools import izip
from math import cos, pi, sin

from igraph.drawing.baseclasses import AbstractDrawer, AbstractCairoDrawer
from igraph.drawing.edge import ArrowEdgeDrawer
from igraph.drawing.metamagic import AttributeCollectorBase, \
                                     AttributeSpecification
from igraph.drawing.shapes import ShapeDrawerDirectory
from igraph.layout import Layout

__all__ = ["DefaultGraphDrawer", "UbiGraphDrawer"]
__license__ = "GPL"

try:
    import cairo
except ImportError:
    # No cairo support is installed. Create a fake module
    # pylint: disable-msg=C0103
    from igraph.drawing.utils import FakeModule
    cairo = FakeModule()

#####################################################################

# pylint: disable-msg=R0903
# R0903: too few public methods
class AbstractGraphDrawer(AbstractDrawer):
    """Abstract class that serves as a base class for anything that
    draws an igraph.Graph."""

    # pylint: disable-msg=W0221
    # W0221: argument number differs from overridden method
    # E1101: Module 'cairo' has no 'foo' member - of course it does :)
    def draw(self, graph, *args, **kwds):
        """Abstract method, must be implemented in derived classes."""
        raise NotImplementedError("abstract class")

#####################################################################

class DefaultGraphDrawer(AbstractGraphDrawer, AbstractCairoDrawer):
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
        AbstractCairoDrawer.__init__(self, context, bbox)
        AbstractGraphDrawer.__init__(self)
        self.edge_drawer_factory = edge_drawer_factory

    # pylint: disable-msg=W0142,W0221,E1101
    # W0142: Used * or ** magic
    # W0221: argument number differs from overridden method
    # E1101: Module 'cairo' has no 'foo' member - of course it does :)
    def draw(self, graph, palette, *args, **kwds):
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
            label_angle = -pi/2
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
            si = sin(vertex.label_angle)
            co = cos(vertex.label_angle)
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

class UbiGraphDrawer(AbstractGraphDrawer):
    """Graph drawer that draws a given graph on an UbiGraph display
    using the XML-RPC API of UbiGraph.
    """

    def __init__(self, url="http://localhost:20738/RPC2"):
        """Constructs an UbiGraph drawer using the display at the given
        URL."""
        from urlparse import urlparse, urlunparse
        import xmlrpclib
        import re

        url_parts = urlparse(url)
        hostname = url_parts.netloc
        if not re.match("[0-9.:]+", hostname):
            # hostname is not an IP address, look it up first
            from socket import gethostbyname
            if ":" in hostname:
                hostname = hostname[0:hostname.index(":")]
            hostname = gethostbyname(hostname)
            if url_parts.port is not None:
                hostname = "%s:%d" % (hostname, url_parts.port)
            url_parts = list(url_parts)
            url_parts[1] = hostname
            url = urlunparse(url_parts)

        self.server = xmlrpclib.ServerProxy(url)
        
    def draw(self, graph, *args, **kwds):
        """Draws the given graph on an UbiGraph display.
        
        @keyword clear: whether to clear the current UbiGraph display before
                        plotting. Default: C{True}."""
        display = self.server.ubigraph

        # Clear the display
        if kwds.get("clear", True):
            display.clear()

        # Add the vertices
        n = graph.vcount()
        vertex_ids = [None] * n
        for i in xrange(n):
            vertex_ids[i] = display.new_vertex()

        # Add the edges
        for edge in graph.es:
            v1, v2 = vertex_ids[edge.source], vertex_ids[edge.target]
            display.new_edge(v1, v2)

#####################################################################


