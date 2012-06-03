"""
Utility classes for drawing routines.
"""

from igraph.compat import property
from itertools import izip
from math import atan2, cos, sin
from operator import itemgetter

__all__ = ["BoundingBox", "FakeModule", "Point", "Rectangle"]
__license__ = "GPL"

#####################################################################

class Rectangle(object):
    """Class representing a rectangle."""

    __slots__ = ("_left", "_top", "_right", "_bottom")

    def __init__(self, *args):
        """Creates a rectangle.

        The corners of the rectangle can be specified by either a tuple
        (four items, two for each corner, respectively), four separate numbers
        (X and Y coordinates for each corner) or two separate numbers (width
        and height, the upper left corner is assumed to be at (0,0))"""
        coords = None
        if len(args) == 1:
            if isinstance(args[0], Rectangle):
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

        self.coords = coords

    @property
    def coords(self):
        """The coordinates of the corners.
        
        The coordinates are returned as a 4-tuple in the following order:
        left edge, top edge, right edge, bottom edge.
        """
        return self._left, self._top, self._right, self._bottom

    @coords.setter
    def coords(self, coords):
        """Sets the coordinates of the corners.

        @param coords: a 4-tuple with the coordinates of the corners
        """
        self._left, self._top, self._right, self._bottom = coords
        if self._left > self._right:
            self._left, self._right = self._right, self._left
        if self._top > self._bottom:
            self._bottom, self._top = self._top, self._bottom

    @property
    def width(self):
        """The width of the rectangle"""
        return self._right - self._left

    @width.setter
    def width(self, value):
        """Sets the width of the rectangle by adjusting the right edge."""
        self._right = self._left + value

    @property
    def height(self):
        """The height of the rectangle"""
        return self._bottom - self._top

    @height.setter
    def height(self, value):
        """Sets the height of the rectangle by adjusting the bottom edge."""
        self._bottom = self._top + value

    @property
    def left(self):
        """The X coordinate of the left side of the box"""
        return self._left

    @left.setter
    def left(self, value):
        """Sets the X coordinate of the left side of the box"""
        self._left = float(value)
        self._right = max(self._left, self._right)

    @property
    def right(self):
        """The X coordinate of the right side of the box"""
        return self._right

    @right.setter
    def right(self, value):
        """Sets the X coordinate of the right side of the box"""
        self._right = float(value)
        self._left = min(self._left, self._right)

    @property
    def top(self):
        """The Y coordinate of the top edge of the box"""
        return self._top

    @top.setter
    def top(self, value):
        """Sets the Y coordinate of the top edge of the box"""
        self._top = value
        self._bottom = max(self._bottom, self._top)

    @property
    def bottom(self):
        """The Y coordinate of the bottom edge of the box"""
        return self._bottom

    @bottom.setter
    def bottom(self, value):
        """Sets the Y coordinate of the bottom edge of the box"""
        self._bottom = value
        self._top = min(self._bottom, self._top)

    @property
    def midx(self):
        """The X coordinate of the center of the box"""
        return (self._left + self._right) / 2.0

    @midx.setter
    def midx(self, value):
        """Moves the center of the box to the given X coordinate"""
        dx = value - (self._left + self._right) / 2.0
        self._left += dx
        self._right += dx

    @property
    def midy(self):
        """The Y coordinate of the center of the box"""
        return (self._top + self._bottom) / 2.0

    @midy.setter
    def midy(self, value):
        """Moves the center of the box to the given Y coordinate"""
        dy = value - (self._top + self._bottom) / 2.0
        self._top += dy
        self._bottom += dy

    @property
    def shape(self):
        """The shape of the rectangle (width, height)"""
        return self._right - self._left, self._bottom - self._top

    @shape.setter
    def shape(self, shape):
        """Sets the shape of the rectangle (width, height)."""
        self.width, self.height = shape

    def contract(self, margins):
        """Contracts the rectangle by the given margins.

        @return: a new L{Rectangle} object.
        """
        if isinstance(margins, int) or isinstance(margins, float):
            margins = [float(margins)] * 4
        if len(margins) != 4:
            raise ValueError("margins must be a 4-tuple or a single number")
        nx1, ny1 = self._left+margins[0], self._top+margins[1]
        nx2, ny2 = self._right-margins[2], self._bottom-margins[3]
        if nx1 > nx2:
            nx1 = (nx1+nx2)/2.
            nx2 = nx1
        if ny1 > ny2:
            ny1 = (ny1+ny2)/2.
            ny2 = ny1
        return self.__class__(nx1, ny1, nx2, ny2)

    def expand(self, margins):
        """Expands the rectangle by the given margins.

        @return: a new L{Rectangle} object.
        """
        if isinstance(margins, int) or isinstance(margins, float):
            return self.contract(-float(margins))
        return self.contract([-float(margin) for margin in margins])

    def isdisjoint(self, other):
        """Returns ``True`` if the two rectangles have no intersection.
        
        Example::
            
            >>> r1 = Rectangle(10, 10, 30, 30)
            >>> r2 = Rectangle(20, 20, 50, 50)
            >>> r3 = Rectangle(70, 70, 90, 90)
            >>> r1.isdisjoint(r2)
            False
            >>> r2.isdisjoint(r1)
            False
            >>> r1.isdisjoint(r3)
            True
            >>> r3.isdisjoint(r1)
            True
        """
        return self._left > other._right or self._right < other._left \
                or self._top > other._bottom or self._bottom < other._top

    def isempty(self):
        """Returns ``True`` if the rectangle is empty (i.e. it has zero
        width and height).

        Example::

            >>> r1 = Rectangle(10, 10, 30, 30)
            >>> r2 = Rectangle(70, 70, 90, 90)
            >>> r1.isempty()
            False
            >>> r2.isempty()
            False
            >>> r1.intersection(r2).isempty()
            True
        """
        return self._left == self._right and self._top == self._bottom

    def intersection(self, other):
        """Returns the intersection of this rectangle with another.
        
        Example::
            
            >>> r1 = Rectangle(10, 10, 30, 30)
            >>> r2 = Rectangle(20, 20, 50, 50)
            >>> r3 = Rectangle(70, 70, 90, 90)
            >>> r1.intersection(r2)
            Rectangle(20.0, 20.0, 30.0, 30.0)
            >>> r2 & r1
            Rectangle(20.0, 20.0, 30.0, 30.0)
            >>> r2.intersection(r1) == r1.intersection(r2)
            True
            >>> r1.intersection(r3)
            Rectangle(0.0, 0.0, 0.0, 0.0)
        """
        if self.isdisjoint(other):
            return Rectangle(0, 0, 0, 0)
        return Rectangle(max(self._left, other._left),
                max(self._top, other._top),
                min(self._right, other._right),
                min(self._bottom, other._bottom))
    __and__ = intersection

    def translate(self, dx, dy):
        """Translates the rectangle in-place.

        Example:

            >>> r = Rectangle(10, 20, 50, 70)
            >>> r.translate(30, -10)
            >>> r
            Rectangle(40.0, 10.0, 80.0, 60.0)

        @param dx: the X coordinate of the translation vector
        @param dy: the Y coordinate of the translation vector
        """
        self._left += dx
        self._right += dx
        self._top += dy
        self._bottom += dy

    def union(self, other):
        """Returns the union of this rectangle with another.
        
        The resulting rectangle is the smallest rectangle that contains both
        rectangles.

        Example::

            >>> r1 = Rectangle(10, 10, 30, 30)
            >>> r2 = Rectangle(20, 20, 50, 50)
            >>> r3 = Rectangle(70, 70, 90, 90)
            >>> r1.union(r2)
            Rectangle(10.0, 10.0, 50.0, 50.0)
            >>> r2 | r1
            Rectangle(10.0, 10.0, 50.0, 50.0)
            >>> r2.union(r1) == r1.union(r2)
            True
            >>> r1.union(r3)
            Rectangle(10.0, 10.0, 90.0, 90.0)
        """
        return Rectangle(min(self._left, other._left),
                min(self._top, other._top),
                max(self._right, other._right),
                max(self._bottom, other._bottom))
    __or__ = union

    def __ior__(self, other):
        """Expands this rectangle to include itself and another completely while
        still being as small as possible.

        Example::

            >>> r1 = Rectangle(10, 10, 30, 30)
            >>> r2 = Rectangle(20, 20, 50, 50)
            >>> r3 = Rectangle(70, 70, 90, 90)
            >>> r1 |= r2
            >>> r1
            Rectangle(10.0, 10.0, 50.0, 50.0)
            >>> r1 |= r3
            >>> r1
            Rectangle(10.0, 10.0, 90.0, 90.0)
        """
        self._left   = min(self._left,   other._left)
        self._top    = min(self._top,    other._top)
        self._right  = max(self._right,  other._right)
        self._bottom = max(self._bottom, other._bottom)
        return self

    def __repr__(self):
        return "%s(%s, %s, %s, %s)" % (self.__class__.__name__, \
            self._left, self._top, self._right, self._bottom)

    def __eq__(self, other):
        return self.coords == other.coords

    def __ne__(self, other):
        return self.coords != other.coords

    def __bool__(self):
        return self._left != self._right or self._top != self._bottom

    def __nonzero__(self):
        return self._left != self._right or self._top != self._bottom

    def __hash__(self):
        return hash(self.coords)

#####################################################################

class BoundingBox(Rectangle):
    """Class representing a bounding box (a rectangular area) that
    encloses some objects."""

    def __ior__(self, other):
        """Replaces this bounding box with the union of itself and
        another.

        Example::
            
            >>> box1 = BoundingBox(10, 20, 50, 60)
            >>> box2 = BoundingBox(70, 40, 100, 90)
            >>> box1 |= box2
            >>> print(box1)
            BoundingBox(10.0, 20.0, 100.0, 90.0)
        """
        self._left   = min(self._left, other._left)
        self._top    = min(self._top, other._top)
        self._right  = max(self._right, other._right)
        self._bottom = max(self._bottom, other._bottom)
        return self

    def __or__(self, other):
        """Takes the union of this bounding box with another.

        The result is a bounding box which encloses both bounding
        boxes.
        
        Example::
            
            >>> box1 = BoundingBox(10, 20, 50, 60)
            >>> box2 = BoundingBox(70, 40, 100, 90)
            >>> box1 | box2
            BoundingBox(10.0, 20.0, 100.0, 90.0)
        """
        return self.__class__(
                       min(self._left, other._left),
                       min(self._top, other._top),
                       max(self._right, other._right),
                       max(self._bottom, other._bottom)
        )


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
    __rmul__ = __mul__

    def __div__(self, scalar):
        """Divides the coordinates by a scalar"""
        return self.__class__(x = self.x / scalar, y = self.y / scalar)

    def as_polar(self):
        """Returns the polar coordinate representation of the point.

        @return: the radius and the angle in a tuple.
        """
        return len(self), atan2(self.y, self.x)

    def distance(self, other):
        """Returns the distance of the point from another one.
        
        Example:
            
            >>> p1 = Point(5, 7)
            >>> p2 = Point(8, 3)
            >>> p1.distance(p2)
            5.0
        """
        dx, dy = self.x - other.x, self.y - other.y
        return (dx * dx + dy * dy) ** 0.5
    
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

    def normalized(self):
        """Normalizes the coordinates of the point s.t. its length will be 1
        after normalization. Returns the normalized point."""
        len = self.length()
        if len == 0:
            return Point(x = self.x, y = self.y)
        return Point(x = self.x / len, y = self.y / len)

    def sq_length(self):
        """Returns the squared length of the vector pointing from the origin
        to this point."""
        return (self.x ** 2 + self.y ** 2)

    def towards(self, other, distance = 0):
        """Returns the point that is at a given distance from this point
        towards another one."""
        if not distance:
            return self

        angle = atan2(other.y - self.y, other.x - self.x)
        return Point(self.x + distance * cos(angle),
                     self.y + distance * sin(angle))

    @classmethod
    def FromPolar(cls, radius, angle):
        """Constructs a point from polar coordinates.

        `radius` is the distance of the point from the origin; `angle` is the
        angle between the X axis and the vector pointing to the point from
        the origin.
        """
        return cls(radius * cos(angle), radius * sin(angle))

