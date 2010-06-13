"""
Drawers for various edge styles in graph plots.
"""

__all__ = ["AbstractEdgeDrawer", "ArrowEdgeDrawer"]

__license__ = "GPL"


import math

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
        ctx.set_source_rgb(*edge.color)
        ctx.set_line_width(edge.width)
        radius = src_vertex.size * 2
        center_x = src_vertex.position[0] + math.cos(math.pi/4) * radius / 2.
        center_y = src_vertex.position[0] - math.sin(math.pi/4) * radius / 2.
        ctx.arc(center_x, center_y, radius/2., 0, math.pi * 2)

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
        ctx.set_source_rgb(*edge.color)
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
        src_point = src_vertex.shape.intersection_point(
                src_pos[0], src_pos[1], dest_pos[0], dest_pos[1],
                src_vertex.size
        )
        dest_point = dest_vertex.shape.intersection_point(
                dest_pos[0], dest_pos[1], src_pos[0], src_pos[1],
                dest_vertex.size
        )

        ctx = self.context

        # Draw the edge
        ctx.set_source_rgb(*edge.color)
        ctx.set_line_width(edge.width)
        ctx.move_to(*src_point)
        ctx.line_to(*dest_point)
        ctx.stroke()

        # Draw the arrowhead
        angle = math.atan2(dest_point[1]-src_point[1], dest_point[0]-src_point[0])
        arrow_size  = 15. * edge.arrow_size
        arrow_width = 10. / edge.arrow_width
        aux_points = [
            (dest_point[0] - arrow_size * math.cos(angle - math.pi/arrow_width),
             dest_point[1] - arrow_size * math.sin(angle - math.pi/arrow_width)),
            (dest_point[0] - arrow_size * math.cos(angle + math.pi/arrow_width),
             dest_point[1] - arrow_size * math.sin(angle + math.pi/arrow_width)),
        ]
        ctx.move_to(*dest_point)
        ctx.line_to(*aux_points[0])
        ctx.line_to(*aux_points[1])
        ctx.line_to(*dest_point)
        ctx.fill()

