"""
Utility classes for drawing routines.
"""

from operator import itemgetter

__all__ = ["BoundingBox", "FakeModule", "Point"]
__license__ = "GPL"

#####################################################################

class BoundingBox(object):
    """Class representing a bounding box (a rectangular area)."""

    __slots__ = ("_coords", )

    def __init__(self, *args):
        """Creates a bounding box.

        The corners of the bounding box can be specified by either a tuple
        (four items, two for each corner, respectively), four separate numbers
        (X and Y coordinates for each corner) or two separate numbers (width
        and height, the upper left corner is assumed to be at (0,0))"""
        coords = None
        if len(args) == 1:
            if isinstance(args[0], BoundingBox):
                coords = args[0].coords
            elif len(args[0]) >= 4:
                coords = tuple(args[0])[0:4]
            elif len(args[0]) == 2:
                coords = (0, 0, args[0][0], args[0][1])
        elif len(args) == 4:
            coords = tuple(args)
        elif len(args) == 2:
            coords = (0, 0, args[0], args[1])
        if coords is None:
            raise ValueError("invalid coordinate format")

        try:
            coords = tuple(float(coord) for coord in coords)
        except ValueError:
            raise ValueError("invalid coordinate format, numbers expected")

        self._coords = None    # to make pylint happy
        self.coords = coords

    def _set_coords(self, coords):
        """Sets the coordinates of the corners.

        @param coords: a 4-tuple with the coordinates of the corners
        """
        if coords[0] > coords[2]:
            coords = (coords[2], coords[1], coords[0], coords[3])
        if coords[1] > coords[3]:
            coords = (coords[0], coords[3], coords[2], coords[1])
        self._coords = coords

    def _get_coords(self):
        """Returns the coordinates of the corners."""
        return self._coords

    coords = property(_get_coords, _set_coords,
        doc="Sets or returns the coordinates of the corners")

    @property
    def width(self):
        """Returns the width of the bounding box"""
        return self._coords[2]-self._coords[0]

    @property
    def height(self):
        """Returns the height of the bounding box"""
        return self._coords[3]-self._coords[1]

    @property
    def left(self):
        """Returns the X coordinate of the left side of the box"""
        return self._coords[0]

    @property
    def right(self):
        """Returns the X coordinate of the right side of the box"""
        return self._coords[2]

    @property
    def top(self):
        """Returns the Y coordinate of the top edge of the box"""
        return self._coords[1]

    @property
    def bottom(self):
        """Returns the Y coordinate of the bottom edge of the box"""
        return self._coords[3]

    @property
    def shape(self):
        """Returns the shape of the bounding box (width, height)"""
        return self._coords[2]-self._coords[0], self._coords[3]-self._coords[1]

    def contract(self, margins):
        """Contracts the bounding box by the given margins.

        @return: a new L{BoundingBox} object.
        """
        if isinstance(margins, int) or isinstance(margins, float):
            margins = [float(margins)] * 4
        if len(margins) != 4:
            raise ValueError("margins must be a 4-tuple or a single number")
        nx1, ny1 = self._coords[0]+margins[0], self._coords[1]+margins[1]
        nx2, ny2 = self._coords[2]-margins[2], self._coords[3]-margins[3]
        if nx1 > nx2:
            nx1 = (nx1+nx2)/2.
            nx2 = nx1
        if ny1 > ny2:
            ny1 = (ny1+ny2)/2.
            ny2 = ny1
        return BoundingBox(nx1, ny1, nx2, ny2)

    def __repr__(self):
        return "%s(%s, %s, %s, %s)" % (self.__class__.__name__, \
            self._coords[0], self._coords[1], self._coords[2], \
            self._coords[3])

    def __eq__(self, other):
        return self.coords == other.coords
    def __ne__(self, other):
        return self.coords != other.coords
    def __hash__(self):
        return hash(self.coords)

#####################################################################

# pylint: disable-msg=R0903
# R0903: too few public methods
class FakeModule(object):
    """Fake module that raises an exception for everything"""

    def __getattr__(self, _):
        raise TypeError("plotting not available")
    def __call__(self, _):
        raise TypeError("plotting not available")
    def __setattr__(self, key, value):
        raise TypeError("plotting not available")

#####################################################################

class Point(tuple):
    """Class representing a point on the 2D plane."""
    __slots__ = ()
    _fields = ('x', 'y')

    def __new__(cls, x, y):
        """Creates a new point with the given coordinates"""
        return tuple.__new__(cls, (x, y))

    # pylint: disable-msg=W0622
    # W0622: redefining built-in 'len'
    @classmethod
    def _make(cls, iterable, new = tuple.__new__, len = len):
        """Creates a new point from a sequence or iterable"""
        result = new(cls, iterable)
        if len(result) != 2:
            raise TypeError('Expected 2 arguments, got %d' % len(result))
        return result

    def __repr__(self):
        """Returns a nicely formatted representation of the point"""
        return 'Point(x=%r, y=%r)' % self

    def _asdict(self):
        """Returns a new dict which maps field names to their values"""
        return dict(zip(self._fields, self))

    # pylint: disable-msg=W0141
    # W0141: used builtin function 'map'
    def _replace(self, **kwds):
        """Returns a new point object replacing specified fields with new
        values"""
        result = self._make(map(kwds.pop, ('x', 'y'), self))
        if kwds:
            raise ValueError('Got unexpected field names: %r' % kwds.keys())
        return result

    def __getnewargs__(self):
        """Return self as a plain tuple. Used by copy and pickle."""
        return tuple(self)

    x = property(itemgetter(0), doc="Alias for field number 0")
    y = property(itemgetter(1), doc="Alias for field number 1")

    def __add__(self, other):
        """Adds the coordinates of a point to another one"""
        return self.__class__(x = self.x + other.x, y = self.y + other.y)

    def __sub__(self, other):
        """Subtracts the coordinates of a point to another one"""
        return self.__class__(x = self.x - other.x, y = self.y - other.y)

    def __mul__(self, scalar):
        """Multiplies the coordinates by a scalar"""
        return self.__class__(x = self.x * scalar, y = self.y * scalar)

    def __div__(self, scalar):
        """Divides the coordinates by a scalar"""
        return self.__class__(x = self.x / scalar, y = self.y / scalar)

    def interpolate(self, other, ratio = 0.5):
        """Linearly interpolates between the coordinates of this point and
        another one.

        @param  other:  the other point
        @param  ratio:  the interpolation ratio between 0 and 1. Zero will
          return this point, 1 will return the other point.
        """
        ratio = float(ratio)
        return Point(x = self.x * (1.0 - ratio) + other.x * ratio, \
                     y = self.y * (1.0 - ratio) + other.y * ratio)

    def length(self):
        """Returns the length of the vector pointing from the origin to this
        point."""
        return (self.x ** 2 + self.y ** 2) ** 0.5

    def sq_length(self):
        """Returns the squared length of the vector pointing from the origin
        to this point."""
        return (self.x ** 2 + self.y ** 2)

