"""
Drawing routines to draw graphs.

This module contains routines to draw graphs on:

  - Cairo surfaces (L{DefaultGraphDrawer})
  - UbiGraph displays (L{UbiGraphDrawer}, see U{http://ubietylab.net/ubigraph})

It also contains routines to send an igraph graph directly to
(U{Cytoscape<http://www.cytoscape.org>}) using the
(U{CytoscapeRPC plugin<http://gforge.nbic.nl/projects/cytoscaperpc/>}), see
L{CytoscapeGraphDrawer}. L{CytoscapeGraphDrawer} can also fetch the current
network from Cytoscape and convert it to igraph format.
"""

from collections import defaultdict
from itertools import izip
from math import cos, pi, sin

from igraph._igraph import convex_hull, VertexSeq
from igraph.configuration import Configuration
from igraph.drawing.baseclasses import AbstractDrawer, AbstractCairoDrawer, \
                                       AbstractXMLRPCDrawer
from igraph.drawing.colors import color_to_html_format, color_name_to_rgb
from igraph.drawing.edge import ArrowEdgeDrawer
from igraph.drawing.text import TextDrawer
from igraph.drawing.metamagic import AttributeCollectorBase, \
                                     AttributeSpecification
from igraph.drawing.shapes import ShapeDrawerDirectory, PolygonDrawer
from igraph.drawing.utils import Point
from igraph.layout import Layout

__all__ = ["DefaultGraphDrawer", "UbiGraphDrawer", "CytoscapeGraphDrawer"]
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

    def ensure_layout(self, layout, graph = None):
        """Helper method that ensures that I{layout} is an instance
        of L{Layout}. If it is not, the method will try to convert
        it to a L{Layout} according to the following rules:

          - If I{layout} is a string, it is assumed to be a name
            of an igraph layout, and it will be passed on to the
            C{layout} method of the given I{graph} if I{graph} is
            not C{None}.

          - If I{layout} is C{None}, the C{layout} method of
            I{graph} will be invoked with no parameters, which
            will call the default layout algorithm.

          - Otherwise, I{layout} will be passed on to the constructor
            of L{Layout}. This handles lists of lists, lists of tuples
            and such.

        If I{layout} is already a L{Layout} instance, it will still
        be copied and a copy will be returned. This is because graph
        drawers are allowed to transform the layout for their purposes,
        and we don't want the transformation to propagate back to the
        caller.
        """
        if isinstance(layout, Layout):
            layout = Layout(layout.coords)
        elif isinstance(layout, str) or layout is None:
            layout = graph.layout(layout)
        else:
            layout = Layout(layout)
        return layout

#####################################################################

class AbstractCairoGraphDrawer(AbstractGraphDrawer, AbstractCairoDrawer):
    """Abstract base class for graph drawers that draw on a Cairo canvas.
    """

    def __init__(self, context, bbox):
        """Constructs the graph drawer and associates it to the given
        Cairo context and the given L{BoundingBox}.

        @param context: the context on which we will draw
        @param bbox:    the bounding box within which we will draw.
                        Can be anything accepted by the constructor
                        of L{BoundingBox} (i.e., a 2-tuple, a 4-tuple
                        or a L{BoundingBox} object).
        """
        AbstractCairoDrawer.__init__(self, context, bbox)
        AbstractGraphDrawer.__init__(self)

#####################################################################

class DefaultGraphDrawer(AbstractCairoGraphDrawer):
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
        AbstractCairoGraphDrawer.__init__(self, context, bbox)
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
        layout = self.ensure_layout(kwds.get("layout", None), graph)

        # Determine the size of the margin on each side
        margin = kwds.get("margin", [0., 0., 0., 0.])
        try:
            margin = list(margin)
        except TypeError:
            margin = [margin]
        while len(margin)<4:
            margin.extend(margin)
        margin = [x + 20. for x in margin[:4]]

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
            label_dist  = 0.0
            label_color = ("black", palette.get)
            label_size  = 14.0
            position = dict(func=layout.__getitem__)
            shape = ("circle", ShapeDrawerDirectory.resolve_default)
            size  = 20.0

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
            # iterate over it and get memberlist-color pairs
            if isinstance(mark_groups, dict):
                group_iter = mark_groups.iteritems()
            elif hasattr(mark_groups, "__iter__"):
                # Lists, tuples, iterators etc
                group_iter = iter(mark_groups)
            else:
                # False
                group_iter = {}.iteritems()

            # We will need a polygon drawer to draw the convex hulls
            polygon_drawer = PolygonDrawer(context, bbox)

            # Iterate over color-memberlist pairs
            for group, color_id in group_iter:
                if not group or color_id is None:
                    continue

                color = palette.get(color_id)

                if isinstance(group, VertexSeq):
                    group = [vertex.index for vertex in group]
                if not hasattr(group, "__iter__"):
                    raise TypeError("group membership list must be iterable")

                # Get the vertex indices that constitute the convex hull
                hull = [group[i] for i in convex_hull([layout[idx] for idx in group])]

                # Calculate the preferred rounding radius for the corners
                corner_radius = 1.25 * max(vertex_builder[idx].size for idx in hull)

                # Construct the polygon
                polygon = [layout[idx] for idx in hull]

                if len(polygon) == 2:
                    # Expand the polygon (which is a flat line otherwise)
                    a, b = Point(*polygon[0]), Point(*polygon[1])
                    c = corner_radius * (a-b).normalized()
                    n = Point(-c[1], c[0])
                    polygon = [a + n, b + n, b - c, b - n, a - n, a + c]
                else:
                    # Expand the polygon around its center of mass
                    center = Point(*[sum(coords) / float(len(coords))
                                      for coords in zip(*polygon)])
                    polygon = [Point(*point).towards(center, -corner_radius)
                               for point in polygon]

                # Draw the hull
                context.set_source_rgba(color[0], color[1], color[2],
                                        color[3]*0.25)
                polygon_drawer.draw_path(polygon, corner_radius=corner_radius)
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
            cairo.FONT_WEIGHT_NORMAL)

        wrap = kwds.get("wrap_labels", None)
        if wrap is None:
            wrap = Configuration.instance()["plotting.wrap_labels"]
        else:
            wrap = bool(wrap)

        label_drawer = TextDrawer(context, halign="center", valign="center")
        for vertex, coords in izip(vertex_builder, layout):
            if vertex.label is None:
                continue

            context.set_font_size(vertex.label_size)
            context.set_source_rgba(*vertex.label_color)
            label_drawer.text = vertex.label

            if vertex.label_dist:
                # Label is displaced from the center of the vertex
                cx, cy = coords
                radius = vertex.label_dist * vertex.size / 2.
                cx = coords[0] + cos(vertex.label_angle) * radius
                cy = coords[1] + sin(vertex.label_angle) * radius
                label_drawer.draw_at(cx, cy, wrap=wrap)
            else:
                # Label is exactly in the center of the vertex
                cx, cy = coords
                half_size = vertex.size / 2.
                label_drawer.bbox = (cx - half_size, cy - half_size,
                                     cx + half_size, cy + half_size)
                label_drawer.draw(wrap=wrap)

#####################################################################

class UbiGraphDrawer(AbstractXMLRPCDrawer, AbstractGraphDrawer):
    """Graph drawer that draws a given graph on an UbiGraph display
    using the XML-RPC API of UbiGraph.

    The following vertex attributes are supported: C{color}, C{label},
    C{shape}, C{size}. See the Ubigraph documentation for supported shape
    names. Sizes are relative to the default Ubigraph size.

    The following edge attributes are supported: C{color}, C{label},
    C{width}. Edge widths are relative to the default Ubigraph width.

    All color specifications supported by igraph (e.g., color names,
    palette indices, RGB triplets, RGBA quadruplets, HTML format)
    are understood by the Ubigraph graph drawer.

    The drawer also has two attributes, C{vertex_defaults} and
    C{edge_defaults}. These are dictionaries that can be used to
    set default values for the vertex/edge attributes in Ubigraph.
    """

    def __init__(self, url="http://localhost:20738/RPC2"):
        """Constructs an UbiGraph drawer using the display at the given
        URL."""
        super(UbiGraphDrawer, self).__init__(url, "ubigraph")
        self.vertex_defaults = dict(
            color="#ff0000",
            shape="cube",
            size=1.0
        )
        self.edge_defaults = dict(
            color="#ffffff",
            width=1.0
        )

    def draw(self, graph, *args, **kwds):
        """Draws the given graph on an UbiGraph display.
        
        @keyword clear: whether to clear the current UbiGraph display before
                        plotting. Default: C{True}."""
        display = self.service

        # Clear the display and set the default visual attributes
        if kwds.get("clear", True):
            display.clear()

            for k, v in self.vertex_defaults.iteritems():
                display.set_vertex_style_attribute(0, k, str(v))
            for k, v in self.edge_defaults.iteritems():
                display.set_edge_style_attribute(0, k, str(v))

        # Custom color converter function
        def color_conv(color):
            return color_to_html_format(color_name_to_rgb(color))

        # Construct the visual vertex/edge builders
        class VisualVertexBuilder(AttributeCollectorBase):
            """Collects some visual properties of a vertex for drawing"""
            _kwds_prefix = "vertex_"
            color = (str(self.vertex_defaults["color"]), color_conv)
            label = None
            shape = str(self.vertex_defaults["shape"])
            size  = float(self.vertex_defaults["size"])

        class VisualEdgeBuilder(AttributeCollectorBase):
            """Collects some visual properties of an edge for drawing"""
            _kwds_prefix = "edge_"
            color = (str(self.edge_defaults["color"]), color_conv)
            label = None
            width = float(self.edge_defaults["width"])

        vertex_builder = VisualVertexBuilder(graph.vs, kwds)
        edge_builder = VisualEdgeBuilder(graph.es, kwds)

        # Add the vertices
        n = graph.vcount()
        new_vertex = display.new_vertex
        vertex_ids = [new_vertex() for _ in xrange(n)]

        # Add the edges
        new_edge = display.new_edge
        eids = [new_edge(vertex_ids[edge.source], vertex_ids[edge.target]) \
                for edge in graph.es]

        # Add arrowheads if needed
        if graph.is_directed():
            display.set_edge_style_attribute(0, "arrow", "true")

        # Set the vertex attributes
        set_attr = display.set_vertex_attribute
        vertex_defaults = self.vertex_defaults
        for vertex_id, vertex in izip(vertex_ids, vertex_builder):
            if vertex.color != vertex_defaults["color"]:
                set_attr(vertex_id, "color", vertex.color)
            if vertex.label:
                set_attr(vertex_id, "label", str(vertex.label))
            if vertex.shape != vertex_defaults["shape"]:
                set_attr(vertex_id, "shape", vertex.shape)
            if vertex.size != vertex_defaults["size"]:
                set_attr(vertex_id, "size", str(vertex.size))

        # Set the edge attributes
        set_attr = display.set_edge_attribute
        edge_defaults = self.edge_defaults
        for edge_id, edge in izip(eids, edge_builder):
            if edge.color != edge_defaults["color"]:
                set_attr(edge_id, "color", edge.color)
            if edge.label:
                set_attr(edge_id, "label", edge.label)
            if edge.width != edge_defaults["width"]:
                set_attr(edge_id, "width", str(edge.width))

#####################################################################

class CytoscapeGraphDrawer(AbstractXMLRPCDrawer, AbstractGraphDrawer):
    """Graph drawer that sends/receives graphs to/from Cytoscape using
    CytoscapeRPC.
    
    This graph drawer cooperates with U{Cytoscape<http://www.cytoscape.org>}
    using U{CytoscapeRPC<http://wiki.nbic.nl/index.php/CytoscapeRPC>}.
    You need to install the CytoscapeRPC plugin first and start the
    XML-RPC server on a given port (port 9000 by default) from the
    appropriate Plugins submenu in Cytoscape.

    Graph, vertex and edge attributes are transferred to Cytoscape whenever
    possible (i.e. when a suitable mapping exists between a Python type
    and a Cytoscape type). If there is no suitable Cytoscape type for a
    Python type, the drawer will use a string attribute on the Cytoscape
    side and invoke C{str()} on the Python attributes.

    If an attribute to be created on the Cytoscape side already exists with
    a different type, an underscore will be appended to the attribute name
    to resolve the type conflict.

    You can use the C{network_id} attribute of this class to figure out the
    network ID of the last graph drawn with this drawer.
    """

    def __init__(self, url="http://localhost:9000/Cytoscape"):
        """Constructs a Cytoscape graph drawer using the XML-RPC interface
        of Cytoscape at the given URL."""
        super(CytoscapeGraphDrawer, self).__init__(url, "Cytoscape")
        self.network_id = None

    def draw(self, graph, name = "Network from igraph", *args, **kwds):
        """Sends the given graph to Cytoscape as a new network.
        
        @param name: the name of the network in Cytoscape."""
        from xmlrpclib import Fault

        cy = self.service

        # Create the network
        network_id = cy.createNetwork(name)
        self.network_id = network_id

        # Create the nodes
        node_ids = [str(idx) for idx in xrange(graph.vcount())]
        cy.createNodes(network_id, node_ids)

        # Create the edges
        edgelists = [[], []]
        for v1, v2 in graph.get_edgelist():
            edgelists[0].append(node_ids[v1])
            edgelists[1].append(node_ids[v2])
        edge_ids = cy.createEdges(network_id,
                edgelists[0], edgelists[1],
                ["unknown"] * graph.ecount(),
                [graph.is_directed()] * graph.ecount(),
                False
        )

        if "layout" in kwds:
            # Calculate/get the layout of the graph
            layout = self.ensure_layout(kwds["layout"], graph)
            size = 100 * graph.vcount() ** 0.5
            layout.fit_into((size, size), keep_aspect_ratio=True)
            layout.translate(-size/2., -size/2.)
            cy.setNodesPositions(network_id,
                    node_ids, *zip(*list(layout)))
        else:
            # Ask Cytoscape to perform the default layout so the user can
            # at least see something in Cytoscape while the attributes are
            # being transferred
            cy.performDefaultLayout(network_id)

        # Send the network attributes
        attr_names = set(cy.getNetworkAttributeNames())
        for attr in graph.attributes():
            cy_type, value = self.infer_cytoscape_type([graph[attr]])
            value = value[0]
            if value is None:
                continue

            # Resolve type conflicts (if any)
            try:
                while attr in attr_names and \
                      cy.getNetworkAttributeType(attr) != cy_type:
                    attr += "_"
            except Fault:
                # getNetworkAttributeType is not available in some older versions
                # so we simply pass here
                pass
            cy.addNetworkAttributes(attr, cy_type, {network_id: value})

        # Send the node attributes
        attr_names = set(cy.getNodeAttributeNames())
        for attr in graph.vertex_attributes():
            cy_type, values = self.infer_cytoscape_type(graph.vs[attr])
            values = dict(pair for pair in izip(node_ids, values)
                    if pair[1] is not None)
            # Resolve type conflicts (if any)
            while attr in attr_names and \
                  cy.getNodeAttributeType(attr) != cy_type:
                attr += "_"
            # Send the attribute values
            cy.addNodeAttributes(attr, cy_type, values, True)

        # Send the edge attributes
        attr_names = set(cy.getEdgeAttributeNames())
        for attr in graph.edge_attributes():
            cy_type, values = self.infer_cytoscape_type(graph.es[attr])
            values = dict(pair for pair in izip(edge_ids, values)
                    if pair[1] is not None)
            # Resolve type conflicts (if any)
            while attr in attr_names and \
                  cy.getEdgeAttributeType(attr) != cy_type:
                attr += "_"
            # Send the attribute values
            cy.addEdgeAttributes(attr, cy_type, values)

    def fetch(self, name = None, directed = False, keep_canonical_names = False):
        """Fetches the network with the given name from Cytoscape.
        
        When fetching networks from Cytoscape, the C{canonicalName} attributes
        of vertices and edges are not converted by default. Use the
        C{keep_canonical_names} parameter to retrieve these attributes as well.

        @param name: the name of the network in Cytoscape.
        @param directed: whether the network is directed.
        @param keep_canonical_names: whether to keep the C{canonicalName}
            vertex/edge attributes that are added automatically by Cytoscape
        @return: an appropriately constructed igraph L{Graph}."""
        from igraph import Graph

        cy = self.service

        # Check the version number. Anything older than 1.3 is bad.
        if map(int, cy.version().split(".")[:2]) < (1, 3):
            raise NotImplementedError("CytoscapeGraphDrawer requires "
                                      "Cytoscape-RPC 1.3 or newer")

        # Find out the ID of the network we are interested in
        if name is None:
            network_id = cy.getNetworkID()
        else:
            network_id = [k for k, v in cy.getNetworkList().iteritems()
                          if v == name]
            if not network_id:
                raise ValueError("no such network: %r" % name)
            elif len(network_id) > 1:
                raise ValueError("more than one network exists with name: %r" % name)
            network_id = network_id[0]

        # Fetch the list of all the nodes and edges
        vertices = cy.getNodes(network_id)
        edges = cy.getEdges(network_id)
        n, m = len(vertices), len(edges)

        # Fetch the graph attributes
        graph_attrs = cy.getNetworkAttributes(network_id)

        # Fetch the vertex attributes
        vertex_attr_names = cy.getNodeAttributeNames()
        vertex_attrs = {}
        for attr_name in vertex_attr_names:
            if attr_name == "canonicalName" and not keep_canonical_names:
                continue
            has_attr = cy.nodesHaveAttribute(attr_name, vertices)
            filtered = [idx for idx, ok in enumerate(has_attr) if ok]
            values = cy.getNodesAttributes(attr_name,
                    [name for name, ok in izip(vertices, has_attr) if ok]
            )
            attrs = [None] * n
            for idx, value in izip(filtered, values):
                attrs[idx] = value
            vertex_attrs[attr_name] = attrs

        # Fetch the edge attributes
        edge_attr_names = cy.getEdgeAttributeNames()
        edge_attrs = {}
        for attr_name in edge_attr_names:
            if attr_name == "canonicalName" and not keep_canonical_names:
                continue
            has_attr = cy.edgesHaveAttribute(attr_name, edges)
            filtered = [idx for idx, ok in enumerate(has_attr) if ok]
            values = cy.getEdgesAttributes(attr_name,
                    [name for name, ok in izip(edges, has_attr) if ok]
            )
            attrs = [None] * m
            for idx, value in izip(filtered, values):
                attrs[idx] = value
            edge_attrs[attr_name] = attrs

        # Create a vertex name index
        vertex_name_index = dict((v, k) for k, v in enumerate(vertices))
        del vertices

        # Remap the edges list to numeric IDs
        edge_list = []
        for edge in edges:
            parts = edge.split()
            edge_list.append((vertex_name_index[parts[0]], vertex_name_index[parts[2]]))
        del edges

        return Graph(n, edge_list, directed=directed,
                graph_attrs=graph_attrs, vertex_attrs=vertex_attrs,
                edge_attrs=edge_attrs)

    @staticmethod
    def infer_cytoscape_type(values):
        """Returns a Cytoscape type that can be used to represent all the
        values in `values` and an appropriately converted copy of `values` that
        is suitable for an XML-RPC call.  Note that the string type in
        Cytoscape is used as a catch-all type; if no other type fits, attribute
        values will be converted to string and then posted to Cytoscape.
        
        ``None`` entries are allowed in `values`, they will be ignored on the
        Cytoscape side.
        """
        types = [type(value) for value in values if value is not None]
        if all(t == bool for t in types):
            return "BOOLEAN", values
        if all(issubclass(t, (int, long)) for t in types):
            return "INTEGER", values
        if all(issubclass(t, float) for t in types):
            return "FLOATING", values
        return "STRING", [
                str(value) if not isinstance(value, basestring) else value
                for value in values
        ]

#####################################################################


