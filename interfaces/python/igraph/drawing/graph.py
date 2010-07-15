"""
Drawing routines to draw graphs.

This module contains routines to draw graphs on:

- Cairo surfaces (L{DefaultGraphDrawer})
- UbiGraph displays (L{UbiGraphDrawer}, see L{http://ubietylab.net/ubigraph})

It also contains routines to send an igraph graph directly to Cytoscape
(L{http://www.cytoscape.org}) using the CytoscapeRPC plugin
(L{http://gforge.nbic.nl/projects/cytoscaperpc/})
"""

from itertools import izip
from math import cos, pi, sin

from igraph.core import convex_hull
from igraph.drawing.baseclasses import AbstractDrawer, AbstractCairoDrawer, \
                                       AbstractXMLRPCDrawer
from igraph.drawing.edge import ArrowEdgeDrawer
from igraph.drawing.text import TextDrawer
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

        # Draw the highlighted groups (if any)
        if "mark_groups" in kwds:
            mark_groups = kwds["mark_groups"]

            # Figure out what to do with mark_groups in order to be able to
            # iterate over it and get color-memberlist pairs
            if isinstance(mark_groups, dict):
                group_iter = mark_groups.iteritems()
            elif hasattr(mark_groups, "__iter__"):
                if hasattr(mark_groups, "next"):
                    # Already an iterator, let's hope it works
                    group_iter = mark_groups
                else:
                    # Lists, tuples etc
                    group_iter = enumerate(mark_groups)
            else:
                # False
                group_iter = {}.iteritems()

            # Iterate over color-memberlist pairs
            for color_id, group in group_iter:
                color = palette.get(color_id)

                if not hasattr(group, "__iter__"):
                    raise TypeError("group membership list must be iterable")

                polygon = convex_hull([layout[idx] for idx in group], True)

                context.set_source_rgba(color[0], color[1], color[2], color[3]*0.25)
                context.move_to(*polygon[-1])
                for point in polygon:
                    context.line_to(*point)
                context.fill_preserve()
                context.set_source_rgba(*color)
                context.stroke()

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

        label_drawer = TextDrawer(context, halign="center")
        for vertex, coords in izip(vertex_builder, layout):
            if vertex.label is None:
                continue

            cx, cy = coords
            radius = vertex.label_dist * vertex.size / 2.
            cx = coords[0] + cos(vertex.label_angle) * radius
            cy = coords[1] + sin(vertex.label_angle) * radius
            # cx += (co - 1) * w/2. + xb
            # cy += (si + 1) * h/2.

            context.set_font_size(vertex.label_size)
            context.set_source_rgba(*vertex.label_color)
            label_drawer.set_text(vertex.label)
            label_drawer.draw(cx, cy)

#####################################################################

class UbiGraphDrawer(AbstractXMLRPCDrawer, AbstractGraphDrawer):
    """Graph drawer that draws a given graph on an UbiGraph display
    using the XML-RPC API of UbiGraph.
    """

    def __init__(self, url="http://localhost:20738/RPC2"):
        """Constructs an UbiGraph drawer using the display at the given
        URL."""
        super(UbiGraphDrawer, self).__init__(url, "ubigraph")
        
    def draw(self, graph, *args, **kwds):
        """Draws the given graph on an UbiGraph display.
        
        @keyword clear: whether to clear the current UbiGraph display before
                        plotting. Default: C{True}."""
        display = self.service

        # Clear the display
        if kwds.get("clear", True):
            display.clear()
            display.set_vertex_style_attribute(0, "color", "#ff0000")

        # Add the vertices
        n = graph.vcount()
        new_vertex = display.new_vertex
        vertex_ids = [new_vertex() for _ in xrange(n)]

        # Add the edges
        new_edge = display.new_edge
        eids = [new_edge(vertex_ids[edge.source], vertex_ids[edge.target]) \
                for edge in graph.es]

        # Add arrowheads if needed
        display.set_edge_style_attribute(0, "arrow", "true")

#####################################################################


