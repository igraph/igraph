"""
Drawing routines to draw the vertices of graphs.

This module provides implementations of vertex drawers, i.e. drawers that the
default graph drawer will use to draw vertices.
"""

from igraph.drawing.baseclasses import AbstractDrawer, AbstractCairoDrawer
from igraph.drawing.metamagic import AttributeCollectorBase
from igraph.drawing.shapes import ShapeDrawerDirectory
from math import pi

__all__ = ["AbstractVertexDrawer", "AbstractCairoVertexDrawer", \
        "DefaultVertexDrawer"]
__license__ = "GPL"

class AbstractVertexDrawer(AbstractDrawer):
    """Abstract vertex drawer object from which all concrete vertex drawer
    implementations are derived."""

    def __init__(self, palette, layout):
        """Constructs the vertex drawer and associates it to the given
        palette.

        @param palette: the palette that can be used to map integer
                        color indices to colors when drawing vertices
        @param layout:  the layout of the vertices in the graph being drawn
        """
        self.layout = layout
        self.palette = palette

    def draw(self, visual_vertex, vertex, coords):
        """Draws the given vertex.
        
        @param visual_vertex: object specifying the visual properties of the
            vertex. Its structure is defined by the VisualVertexBuilder of the
            L{DefaultGraphDrawer}; see its source code.
        @param vertex: the raw igraph vertex being drawn
        @param coords: the X and Y coordinates of the vertex as specified by the
            layout algorithm, scaled into the bounding box.
        """
        raise NotImplementedError("abstract class")

class AbstractCairoVertexDrawer(AbstractVertexDrawer, AbstractCairoDrawer):
    """Abstract base class for vertex drawers that draw on a Cairo canvas."""

    def __init__(self, context, bbox, palette, layout):
        """Constructs the vertex drawer and associates it to the given
        Cairo context and the given L{BoundingBox}.

        @param context: the context on which we will draw
        @param bbox:    the bounding box within which we will draw.
                        Can be anything accepted by the constructor
                        of L{BoundingBox} (i.e., a 2-tuple, a 4-tuple
                        or a L{BoundingBox} object).
        @param palette: the palette that can be used to map integer
                        color indices to colors when drawing vertices
        @param layout:  the layout of the vertices in the graph being drawn
        """
        AbstractCairoDrawer.__init__(self, context, bbox)
        AbstractVertexDrawer.__init__(self, palette, layout)

class DefaultVertexDrawer(AbstractCairoVertexDrawer):
    """The default vertex drawer implementation of igraph."""

    def __init__(self, context, bbox, palette, layout):
        AbstractCairoVertexDrawer.__init__(self, context, bbox, palette, layout)
        self.VisualVertexBuilder = self._construct_visual_vertex_builder()

    def _construct_visual_vertex_builder(self):
        class VisualVertexBuilder(AttributeCollectorBase):
            """Collects some visual properties of a vertex for drawing"""
            _kwds_prefix = "vertex_"
            color = ("red", self.palette.get)
            frame_color = ("black", self.palette.get)
            frame_width = 1.0
            label = None
            label_angle = -pi/2
            label_dist  = 0.0
            label_color = ("black", self.palette.get)
            label_size  = 14.0
            position = dict(func=self.layout.__getitem__)
            shape = ("circle", ShapeDrawerDirectory.resolve_default)
            size  = 20.0
        return VisualVertexBuilder

    def draw(self, visual_vertex, vertex, coords):
        context = self.context

        visual_vertex.shape.draw_path(context, \
                coords[0], coords[1], visual_vertex.size)
        context.set_source_rgba(*visual_vertex.color)
        context.fill_preserve()
        context.set_source_rgba(*visual_vertex.frame_color)
        context.set_line_width(visual_vertex.frame_width)
        context.stroke()


