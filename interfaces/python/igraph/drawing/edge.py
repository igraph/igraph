"""
Drawers for various edge styles in graph plots.
"""

__all__ = ["AbstractEdgeDrawer", "AlphaVaryingEdgeDrawer",
           "ArrowEdgeDrawer", "DarkToLightEdgeDrawer",
           "LightToDarkEdgeDrawer", "TaperedEdgeDrawer"]

__license__ = "GPL"

try:
    from cairo import LinearGradient
except ImportError:
    # No cairo support is installed. Don't worry, there will
    # be a fake Cairo module in igraph.drawing
    pass

from igraph.drawing.colors import clamp
from igraph.drawing.metamagic import AttributeCollectorBase
from igraph.drawing.text import TextAlignment
from math import atan2, cos, pi, sin

class AbstractEdgeDrawer(object):
    """Abstract edge drawer object from which all concrete edge drawer
    implementations are derived."""

    def __init__(self, context, palette):
        """Constructs the edge drawer.

        @param context: a Cairo context on which the edges will be drawn.
        @param palette: the palette that can be used to map integer
                        color indices to colors when drawing edges
        """
        self.context = context
        self.palette = palette
        self.VisualEdgeBuilder = self._construct_visual_edge_builder()

    @staticmethod
    def _curvature_to_float(value):
        """Converts values given to the 'curved' edge style argument
        in plotting calls to floating point values."""
        if value is None or value is False:
            return 0.0
        if value is True:
            return 0.5
        return float(value)

    def _construct_visual_edge_builder(self):
        """Construct the visual edge builder that will collect the visual
        attributes of an edge when it is being drawn."""
        class VisualEdgeBuilder(AttributeCollectorBase):
            """Builder that collects some visual properties of an edge for
            drawing"""
            _kwds_prefix = "edge_"
            arrow_size  = 1.0
            arrow_width = 1.0
            color       = ("#444", self.palette.get)
            curved      = (0.0, self._curvature_to_float)
            label       = None
            label_color = ("black", self.palette.get)
            label_size  = 12.0
            width       = 1.0
        return VisualEdgeBuilder

    def draw_directed_edge(self, edge, src_vertex, dest_vertex):
        """Draws a directed edge.

        @param edge: the edge to be drawn. Visual properties of the edge
          are defined by the attributes of this object.
        @param src_vertex: the source vertex. Visual properties are given
          again as attributes.
        @param dest_vertex: the target vertex. Visual properties are given
          again as attributes.
        """
        raise NotImplementedError()

    def draw_loop_edge(self, edge, vertex):
        """Draws a loop edge.

        The default implementation draws a small circle.

        @param edge: the edge to be drawn. Visual properties of the edge
          are defined by the attributes of this object.
        @param vertex: the vertex to which the edge is attached. Visual
          properties are given again as attributes.
        """
        ctx = self.context
        ctx.set_source_rgba(*edge.color)
        ctx.set_line_width(edge.width)
        radius = vertex.size * 1.5
        center_x = vertex.position[0] + cos(pi/4) * radius / 2.
        center_y = vertex.position[1] - sin(pi/4) * radius / 2.
        ctx.arc(center_x, center_y, radius/2., 0, pi * 2)
        ctx.stroke()

    def draw_undirected_edge(self, edge, src_vertex, dest_vertex):
        """Draws an undirected edge.

        The default implementation of this method draws undirected edges
        as straight lines. Loop edges are drawn as small circles.

        @param edge: the edge to be drawn. Visual properties of the edge
          are defined by the attributes of this object.
        @param src_vertex: the source vertex. Visual properties are given
          again as attributes.
        @param dest_vertex: the target vertex. Visual properties are given
          again as attributes.
        """
        if src_vertex == dest_vertex:    # TODO
            return self.draw_loop_edge(edge, src_vertex)

        ctx = self.context
        ctx.set_source_rgba(*edge.color)
        ctx.set_line_width(edge.width)
        ctx.move_to(*src_vertex.position)

        if edge.curved:
            (x1, y1), (x2, y2) = src_vertex.position, dest_vertex.position
            aux1 = (2*x1+x2) / 3.0 - edge.curved * 0.5 * (y2-y1), \
                   (2*y1+y2) / 3.0 + edge.curved * 0.5 * (x2-x1)
            aux2 = (x1+2*x2) / 3.0 - edge.curved * 0.5 * (y2-y1), \
                   (y1+2*y2) / 3.0 + edge.curved * 0.5 * (x2-x1)
            ctx.curve_to(aux1[0], aux1[1], aux2[0], aux2[1], *dest_vertex.position)
        else:
            ctx.line_to(*dest_vertex.position)

        ctx.stroke()

    def get_label_position(self, edge, src_vertex, dest_vertex):
        """Returns the position where the label of an edge should be drawn. The
        default implementation returns the midpoint of the edge and an alignment
        that tries to avoid overlapping the label with the edge.
        
        @param edge: the edge to be drawn. Visual properties of the edge
          are defined by the attributes of this object.
        @param src_vertex: the source vertex. Visual properties are given
          again as attributes.
        @param dest_vertex: the target vertex. Visual properties are given
          again as attributes.
        @return: a tuple containing two more tuples: the desired position of the
          label and the desired alignment of the label, where the position is
          given as C{(x, y)} and the alignment is given as C{(horizontal, vertical)}.
          Members of the alignment tuple are taken from constants in the
          L{TextAlignment} class.
        """
        # Determine the angle of the line
        dx = dest_vertex.position[0] - src_vertex.position[0]
        dy = dest_vertex.position[1] - src_vertex.position[1]
        if dx != 0 or dy != 0:
            # Note that we use -dy because the Y axis points downwards
            angle = atan2(-dy, dx) % (2*pi)
        else:
            angle = None

        # Determine the midpoint
        pos = ((src_vertex.position[0] + dest_vertex.position[0]) / 2., \
                (src_vertex.position[1] + dest_vertex.position[1]) / 2)

        # Determine the alignment based on the angle
        pi4 = pi / 4
        if angle is None:
            halign, valign = TextAlignment.CENTER, TextAlignment.CENTER
        else:
            index = int((angle / pi4) % 8)
            halign = [TextAlignment.RIGHT, TextAlignment.RIGHT,
                    TextAlignment.RIGHT, TextAlignment.RIGHT,
                    TextAlignment.LEFT, TextAlignment.LEFT,
                    TextAlignment.LEFT, TextAlignment.LEFT][index]
            valign = [TextAlignment.BOTTOM, TextAlignment.CENTER,
                    TextAlignment.CENTER, TextAlignment.TOP,
                    TextAlignment.TOP, TextAlignment.CENTER,
                    TextAlignment.CENTER, TextAlignment.BOTTOM][index]

        return pos, (halign, valign)


class ArrowEdgeDrawer(AbstractEdgeDrawer):
    """Edge drawer implementation that draws undirected edges as
    straight lines and directed edges as arrows.
    """

    def draw_directed_edge(self, edge, src_vertex, dest_vertex):
        if src_vertex == dest_vertex:    # TODO
            return self.draw_loop_edge(edge, src_vertex)

        ctx = self.context
        (x1, y1), (x2, y2) = src_vertex.position, dest_vertex.position


        # Draw the edge
        ctx.set_source_rgba(*edge.color)
        ctx.set_line_width(edge.width)
        ctx.move_to(x1, y1)

        if edge.curved:
            # Calculate the curve
            aux1 = (2*x1+x2) / 3.0 - edge.curved * 0.5 * (y2-y1), \
                   (2*y1+y2) / 3.0 + edge.curved * 0.5 * (x2-x1)
            aux2 = (x1+2*x2) / 3.0 - edge.curved * 0.5 * (y2-y1), \
                   (y1+2*y2) / 3.0 + edge.curved * 0.5 * (x2-x1)
            ctx.curve_to(aux1[0], aux1[1], aux2[0], aux2[1], x2, y2)
            x1, y1 = aux2
        else:
            # Draw the line
            ctx.line_to(x2, y2)

        # Determine where the edge intersects the circumference of the
        # vertex shape.
        x2, y2 = dest_vertex.shape.intersection_point(
                x2, y2, x1, y1, dest_vertex.size)

        ctx.stroke()

        # Draw the arrowhead
        angle = atan2(y2-y1, x2-x1)
        arrow_size  = 15. * edge.arrow_size
        arrow_width = 10. / edge.arrow_width
        aux_points = [
            (x2 - arrow_size * cos(angle - pi/arrow_width),
             y2 - arrow_size * sin(angle - pi/arrow_width)),
            (x2 - arrow_size * cos(angle + pi/arrow_width),
             y2 - arrow_size * sin(angle + pi/arrow_width)),
        ]
        ctx.move_to(x2, y2)
        ctx.line_to(*aux_points[0])
        ctx.line_to(*aux_points[1])
        ctx.line_to(x2, y2)
        ctx.fill()


class TaperedEdgeDrawer(AbstractEdgeDrawer):
    """Edge drawer implementation that draws undirected edges as
    straight lines and directed edges as tapered lines that are
    wider at the source and narrow at the destination.
    """

    def draw_directed_edge(self, edge, src_vertex, dest_vertex):
        if src_vertex == dest_vertex:    # TODO
            return self.draw_loop_edge(edge, src_vertex)

        # Determine where the edge intersects the circumference of the
        # destination vertex.
        src_pos, dest_pos = src_vertex.position, dest_vertex.position
        dest_pos = dest_vertex.shape.intersection_point(
                dest_pos[0], dest_pos[1], src_pos[0], src_pos[1],
                dest_vertex.size
        )

        ctx = self.context

        # Draw the edge
        ctx.set_source_rgba(*edge.color)
        ctx.set_line_width(edge.width)
        angle = atan2(dest_pos[1]-src_pos[1], dest_pos[0]-src_pos[0])
        arrow_size  = src_vertex.size / 4.
        aux_points = [
            (src_pos[0] + arrow_size * cos(angle + pi/2),
             src_pos[1] + arrow_size * sin(angle + pi/2)),
            (src_pos[0] + arrow_size * cos(angle - pi/2),
             src_pos[1] + arrow_size * sin(angle - pi/2))
        ]
        ctx.move_to(*dest_pos)
        ctx.line_to(*aux_points[0])
        ctx.line_to(*aux_points[1])
        ctx.line_to(*dest_pos)
        ctx.fill()


class AlphaVaryingEdgeDrawer(AbstractEdgeDrawer):
    """Edge drawer implementation that draws undirected edges as
    straight lines and directed edges by varying the alpha value
    of the specified edge color between the source and the destination.
    """

    def __init__(self, context, alpha_at_src, alpha_at_dest):
        super(AlphaVaryingEdgeDrawer, self).__init__(context)
        self.alpha_at_src = (clamp(float(alpha_at_src), 0., 1.), )
        self.alpha_at_dest = (clamp(float(alpha_at_dest), 0., 1.), )

    def draw_directed_edge(self, edge, src_vertex, dest_vertex):
        if src_vertex == dest_vertex:    # TODO
            return self.draw_loop_edge(edge, src_vertex)

        src_pos, dest_pos = src_vertex.position, dest_vertex.position
        ctx = self.context

        # Set up the gradient
        lg = LinearGradient(src_pos[0], src_pos[1], dest_pos[0], dest_pos[1])
        edge_color = edge.color[:3] + self.alpha_at_src
        edge_color_end = edge_color[:3] + self.alpha_at_dest
        lg.add_color_stop_rgba(0, *edge_color)
        lg.add_color_stop_rgba(1, *edge_color_end)

        # Draw the edge
        ctx.set_source(lg)
        ctx.set_line_width(edge.width)
        ctx.move_to(*src_pos)
        ctx.line_to(*dest_pos)
        ctx.stroke()


class LightToDarkEdgeDrawer(AlphaVaryingEdgeDrawer):
    """Edge drawer implementation that draws undirected edges as
    straight lines and directed edges by using an alpha value of
    zero (total transparency) at the source and an alpha value of
    one (full opacity) at the destination. The alpha value is
    interpolated in-between.
    """

    def __init__(self, context):
        super(LightToDarkEdgeDrawer, self).__init__(context, 0.0, 1.0)


class DarkToLightEdgeDrawer(AlphaVaryingEdgeDrawer):
    """Edge drawer implementation that draws undirected edges as
    straight lines and directed edges by using an alpha value of
    one (full opacity) at the source and an alpha value of zero
    (total transparency) at the destination. The alpha value is
    interpolated in-between.
    """

    def __init__(self, context):
        super(DarkToLightEdgeDrawer, self).__init__(context, 1.0, 0.0)

