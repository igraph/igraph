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
from math import atan2, cos, pi, sin

class AbstractEdgeDrawer(object):
    """Abstract edge drawer object from which all concrete edge drawer
    implementations are derived."""

    def __init__(self, context):
        """Constructs the edge drawer.

        @param context: a Cairo context on which the edges will be drawn.
        """
        self.context = context

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
        radius = src_vertex.size * 2
        center_x = src_vertex.position[0] + cos(pi/4) * radius / 2.
        center_y = src_vertex.position[0] - sin(pi/4) * radius / 2.
        ctx.arc(center_x, center_y, radius/2., 0, pi * 2)

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
        ctx.line_to(*dest_vertex.position)
        ctx.stroke()


class ArrowEdgeDrawer(AbstractEdgeDrawer):
    """Edge drawer implementation that draws undirected edges as
    straight lines and directed edges as arrows.
    """

    def draw_directed_edge(self, edge, src_vertex, dest_vertex):
        if src_vertex == dest_vertex:    # TODO
            return self.draw_loop_edge(edge, src_vertex)

        # Determine where the edge intersects the circumference of the
        # vertex shape.
        src_pos, dest_pos = src_vertex.position, dest_vertex.position
        dest_pos = dest_vertex.shape.intersection_point(
                dest_pos[0], dest_pos[1], src_pos[0], src_pos[1],
                dest_vertex.size
        )

        ctx = self.context

        # Draw the edge
        ctx.set_source_rgba(*edge.color)
        ctx.set_line_width(edge.width)
        ctx.move_to(*src_pos)
        ctx.line_to(*dest_pos)
        ctx.stroke()

        # Draw the arrowhead
        angle = atan2(dest_pos[1]-src_pos[1], dest_pos[0]-src_pos[0])
        arrow_size  = 15. * edge.arrow_size
        arrow_width = 10. / edge.arrow_width
        aux_points = [
            (dest_pos[0] - arrow_size * cos(angle - pi/arrow_width),
             dest_pos[1] - arrow_size * sin(angle - pi/arrow_width)),
            (dest_pos[0] - arrow_size * cos(angle + pi/arrow_width),
             dest_pos[1] - arrow_size * sin(angle + pi/arrow_width)),
        ]
        ctx.move_to(*dest_pos)
        ctx.line_to(*aux_points[0])
        ctx.line_to(*aux_points[1])
        ctx.line_to(*dest_pos)
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

