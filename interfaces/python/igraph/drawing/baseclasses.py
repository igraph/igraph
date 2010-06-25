"""
Abstract base classes for the drawing routines.
"""

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
        if not isinstance(bbox, BoundingBox):
            bbox = BoundingBox(bbox)
        self.bbox = bbox

    def draw(self, *args, **kwds):
        """Abstract method, must be implemented in derived classes."""
        raise NotImplementedError("abstract class")


