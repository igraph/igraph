"""
Drawers for various edge styles in graph plots.
"""

__all__ = ["AbstractEdgeDrawer", "AlphaVaryingEdgeDrawer",
           "ArrowEdgeDrawer", "DarkToLightEdgeDrawer",
           "LightToDarkEdgeDrawer", "TaperedEdgeDrawer"]

__license__ = "GPL"

from igraph.drawing.colors import clamp
from igraph.drawing.metamagic import AttributeCollectorBase
from igraph.drawing.text import TextAlignment
from igraph.drawing.utils import find_cairo
from math import atan2, cos, pi, sin, sqrt

cairo = find_cairo()

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
        (x_src, y_src), (x_dest, y_dest) = src_vertex.position, dest_vertex.position


        def bezier_cubic(x0,y0, x1,y1, x2,y2, x3,y3, t):
            """ Computes the Bezier curve from point (x0,y0) to (x3,y3) 
            via control points (x1,y1) and (x2,y2) with parameter t.
            """
            xt = (1.0 - t) ** 3 * x0 + 3. *t * (1.0 - t) ** 2 * x1 + 3. * t**2 * (1. - t) * x2 + t**3 * x3 
            yt = (1.0 - t) ** 3 * y0 + 3. *t * (1.0 - t) ** 2 * y1 + 3. * t**2 * (1. - t) * y2 + t**3 * y3 
            return xt,yt
        
        def euclidean_distance(x1,y1,x2,y2):
            """ Computes the Euclidean distance between points (x1,y1) and (x2,y2).
            """
            return sqrt( (1.0*x1-x2) **2 + (1.0*y1-y2) **2 )

        def intersect_bezier_circle(x0,y0, x1,y1, x2,y2, x3,y3, radius):
            """ Binary search solver for finding the intersection of a Bezier curve
            and a circle centered at the curve's end point.
            Returns the x,y of the intersection point.
            TODO: implement safeguard to ensure convergence in ALL possible cases.
            """
            precision = radius / 20.0
            source_target_distance = euclidean_distance(x0,y0,x3,y3)
            radius = float(radius)
            t0 = 1.0
            t1 = 1.0 - radius / source_target_distance
            
            xt0, yt0 = x3, y3
            xt1, yt1 = bezier_cubic(x0,y0, x1,y1, x2,y2, x3,y3, t1)
            
            distance_t0 = 0
            distance_t1 = euclidean_distance(x3,y3, xt1,yt1)
            counter = 0
            while abs(distance_t1 - radius) > precision:
                if ((distance_t1-radius) > 0) !=  ((distance_t0-radius) > 0):
                    t_new = (t0 + t1)/2.0
                else:
                    if (abs(distance_t1 - radius) < abs(distance_t0 - radius)):
                        # If t1 gets us closer to the circumference step in the same direction
                        t_new = t1 + (t1 - t0)/ 2.0
                    else:
                        t_new = t1 - (t1 - t0)
                t_new = 1 if t_new > 1 else (0 if t_new < 0 else t_new)        
                t0,t1 = t1,t_new
                distance_t0 = distance_t1
                xt1, yt1 = bezier_cubic(x0,y0, x1,y1, x2,y2, x3,y3, t1)
                distance_t1 = euclidean_distance(x3,y3, xt1,yt1)
                counter += 1
            return bezier_cubic(x0,y0, x1,y1, x2,y2, x3,y3, t1)



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

            # Coordinates of the control points of the Bezier curve
            xc1, yc1 = aux1 
            xc2, yc2 = aux2

            # Determine where the edge intersects the circumference of the
            # vertex shape: Tip of the arrow
            x2, y2 = intersect_bezier_circle(x_src,y_src, xc1,yc1, xc2,yc2, x_dest,y_dest, dest_vertex.size/2.0)

            # Calculate the arrow head coordinates
            angle = atan2(y_dest - y2, x_dest - x2) # navid
            arrow_size  = 15. * edge.arrow_size
            arrow_width = 10. / edge.arrow_width
            aux_points = [
                        (x2 - arrow_size * cos(angle - pi/arrow_width),
                         y2 - arrow_size * sin(angle - pi/arrow_width)),
                        (x2 - arrow_size * cos(angle + pi/arrow_width),
                         y2 - arrow_size * sin(angle + pi/arrow_width)),
                    ]

            # Midpoint of the base of the arrow triangle
            x_arrow_mid , y_arrow_mid = (aux_points [0][0] + aux_points [1][0]) / 2.0, (aux_points [0][1] + aux_points [1][1]) / 2.0   

            # Vector representing the base of the arrow triangle
            x_arrow_base_vec, y_arrow_base_vec = (aux_points [0][0] - aux_points [1][0]) , (aux_points [0][1] - aux_points [1][1])   

            # Recalculate the curve such that it lands on the base of the arrow triangle
            aux1 = (2*x_src+x_arrow_mid) / 3.0 - edge.curved * 0.5 * (y_arrow_mid-y_src), \
                   (2*y_src+y_arrow_mid) / 3.0 + edge.curved * 0.5 * (x_arrow_mid-x_src)
            aux2 = (x_src+2*x_arrow_mid) / 3.0 - edge.curved * 0.5 * (y_arrow_mid-y_src), \
                   (y_src+2*y_arrow_mid) / 3.0 + edge.curved * 0.5 * (x_arrow_mid-x_src)

            # Offset the second control point (aux2) such that it falls precisely on the normal to the arrow base vector
            # Strictly speaking, offset_length is the offset length divided by the length of the arrow base vector.
            offset_length = (x_arrow_mid - aux2[0]) * x_arrow_base_vec + (y_arrow_mid - aux2[1]) * y_arrow_base_vec 
            offset_length /= euclidean_distance(0,0, x_arrow_base_vec, y_arrow_base_vec) ** 2
            
            aux2 = aux2[0] + x_arrow_base_vec * offset_length, \
                   aux2[1] + y_arrow_base_vec * offset_length

            # Draw tthe curve from the first vertex to the midpoint of the base of the arrow head
            ctx.curve_to(aux1[0], aux1[1], aux2[0], aux2[1], x_arrow_mid, y_arrow_mid)
        else:
            # Determine where the edge intersects the circumference of the
            # vertex shape.
            x2, y2 = dest_vertex.shape.intersection_point(
                    x2, y2, x1, y1, dest_vertex.size)

            # Draw the arrowhead
            angle = atan2(y_dest - y2, x_dest - x2) 
            arrow_size  = 15. * edge.arrow_size
            arrow_width = 10. / edge.arrow_width
            aux_points = [
                (x2 - arrow_size * cos(angle - pi/arrow_width),
                 y2 - arrow_size * sin(angle - pi/arrow_width)),
                (x2 - arrow_size * cos(angle + pi/arrow_width),
                 y2 - arrow_size * sin(angle + pi/arrow_width)),
            ]

            # Midpoint of the base of the arrow triangle
            x_arrow_mid , y_arrow_mid = (aux_points [0][0] + aux_points [1][0]) / 2.0, (aux_points [0][1] + aux_points [1][1]) / 2.0   
            # Draw the line
            ctx.line_to(x_arrow_mid, y_arrow_mid)

        # Draw the edge
        ctx.stroke()


        # Draw the arrow head
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
        lg = cairo.LinearGradient(src_pos[0], src_pos[1], dest_pos[0], dest_pos[1])
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

