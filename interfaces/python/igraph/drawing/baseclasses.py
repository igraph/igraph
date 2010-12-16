"""
Abstract base classes for the drawing routines.
"""

from igraph.compat import property
from igraph.drawing.utils import BoundingBox

#####################################################################

# pylint: disable-msg=R0903
# R0903: too few public methods
class AbstractDrawer(object):
    """Abstract class that serves as a base class for anything that
    draws an igraph object."""

    def draw(self, *args, **kwds):
        """Abstract method, must be implemented in derived classes."""
        raise NotImplementedError("abstract class")

#####################################################################

# pylint: disable-msg=R0903
# R0903: too few public methods
class AbstractCairoDrawer(AbstractDrawer):
    """Abstract class that serves as a base class for anything that
    draws on a Cairo context within a given bounding box.
    
    A subclass of L{AbstractCairoDrawer} is guaranteed to have an
    attribute named C{context} that represents the Cairo context
    to draw on, and an attribute named C{bbox} for the L{BoundingBox}
    of the drawing area.
    """

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
        self._bbox = None
        self.bbox = bbox

    @property
    def bbox(self):
        """The bounding box of the drawing area where this drawer will
        draw."""
        return self._bbox

    @bbox.setter
    def bbox(self, bbox):
        """Sets the bounding box of the drawing area where this drawer
        will draw."""
        if not isinstance(bbox, BoundingBox):
            self._bbox = BoundingBox(bbox)
        else:
            self._bbox = bbox

    def draw(self, *args, **kwds):
        """Abstract method, must be implemented in derived classes."""
        raise NotImplementedError("abstract class")

#####################################################################

class AbstractXMLRPCDrawer(AbstractDrawer):
    """Abstract drawer that uses a remote service via XML-RPC
    to draw something on a remote display.
    """

    def __init__(self, url, service=None):
        """Constructs an abstract drawer using the XML-RPC service
        at the given URL.
        
        @param url: the URL where the XML-RPC calls for the service should
          be addressed to.
        @param service: the name of the service at the XML-RPC address. If
          C{None}, requests will be directed to the server proxy object
          constructed by L{xmlrpclib.ServerProxy}; if not C{None}, the
          given attribute will be looked up in the server proxy object.
        """
        import xmlrpclib
        url = self._resolve_hostname(url)
        self.server = xmlrpclib.ServerProxy(url)
        if service is None:
            self.service = self.server
        else:
            self.service = getattr(self.server, service)

    @staticmethod
    def _resolve_hostname(url):
        """Parses the given URL, resolves the hostname to an IP address
        and returns a new URL with the resolved IP address. This speeds
        up things big time on Mac OS X where an IP lookup would be
        performed for every XML-RPC call otherwise."""
        from urlparse import urlparse, urlunparse
        import re

        url_parts = urlparse(url)
        hostname = url_parts.netloc
        if re.match("[0-9.:]+$", hostname):
            # the hostname is already an IP address, possibly with a port
            return url

        from socket import gethostbyname
        if ":" in hostname:
            hostname = hostname[0:hostname.index(":")]
        hostname = gethostbyname(hostname)
        if url_parts.port is not None:
            hostname = "%s:%d" % (hostname, url_parts.port)
        url_parts = list(url_parts)
        url_parts[1] = hostname
        return urlunparse(url_parts)

