# vim:ts=4:sw=4:sts=4:et
# -*- coding: utf-8 -*-
"""Additional auxiliary data types"""

from itertools import islice

__license__ = """\
Copyright (C) 2006-2012  Tamás Nepusz <ntamas@gmail.com>
Pázmány Péter sétány 1/a, 1117 Budapest, Hungary

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
02110-1301 USA
"""

class Matrix(object):
    """Simple matrix data type.

    Of course there are much more advanced matrix data types for Python (for
    instance, the C{ndarray} data type of Numeric Python) and this implementation
    does not want to compete with them. The only role of this data type is to
    provide a convenient interface for the matrices returned by the C{Graph}
    object (for instance, allow indexing with tuples in the case of adjacency
    matrices and so on).
    """

    def __init__(self, data=None):
        """Initializes a matrix.

        @param data: the elements of the matrix as a list of lists, or C{None} to
          create a 0x0 matrix.
        """
        self._nrow, self._ncol, self._data = 0, 0, []
        self.data = data

    # pylint: disable-msg=C0103
    @classmethod
    def Fill(cls, value, *args):
        """Creates a matrix filled with the given value

        @param value: the value to be used
        @keyword shape: the shape of the matrix. Can be a single integer,
          two integers or a tuple. If a single integer is
          given here, the matrix is assumed to be square-shaped.
        """
        if len(args) < 1:
            raise TypeError("expected an integer or a tuple")
        if len(args) == 1:
            if hasattr(args[0], "__len__"):
                height, width = int(args[0][0]), int(args[0][1])
            else:
                height, width = int(args[0]), int(args[0])
        else:
            height, width = int(args[0]), int(args[1])
        mtrx = [[value]*width for _ in xrange(height)]
        return cls(mtrx)

    # pylint: disable-msg=C0103
    @classmethod
    def Zero(cls, *args):
        """Creates a matrix filled with zeros.

        @keyword shape: the shape of the matrix. Can be a single integer,
          two integers or a tuple. If a single integer is
          given here, the matrix is assumed to be square-shaped.
        """
        result = cls.Fill(0, *args)
        return result

    # pylint: disable-msg=C0103
    @classmethod
    def Identity(cls, *args):
        """Creates an identity matrix.

        @keyword shape: the shape of the matrix. Can be a single integer,
          two integers or a tuple. If a single integer is
          given here, the matrix is assumed to be square-shaped.
        """
        # pylint: disable-msg=W0212
        result = cls.Fill(0, *args)
        for i in xrange(min(result.shape)):
            result._data[i][i] = 1
        return result

    def _set_data(self, data=None):
        """Sets the data stored in the matrix"""
        if data is not None:
            self._data = [list(row) for row in data]
            self._nrow = len(self._data)
            if self._nrow > 0:
                self._ncol = max(len(row) for row in self._data)
            else:
                self._ncol = 0
            for row in self._data:
                if len(row) < self._ncol:
                    row.extend([0]*(self._ncol-len(row)))

    def _get_data(self):
        """Returns the data stored in the matrix as a list of lists"""
        return [list(row) for row in self._data]
    data = property(_get_data, _set_data)

    @property
    def shape(self):
        """Returns the shape of the matrix as a tuple"""
        return self._nrow, self._ncol

    def __add__(self, other):
        """Adds the given value to the matrix.

        @param other: either a scalar or a matrix. Scalars will
          be added to each element of the matrix. Matrices will
          be added together elementwise.
        @return: the result matrix
        """
        if isinstance(other, Matrix):
            if self.shape != other.shape:
                raise ValueError("matrix shapes do not match")
            return self.__class__([
                [a+b for a, b in izip(row_a, row_b)]
                for row_a, row_b in izip(self, other)
            ])
        else:
            return self.__class__([
                [item+other for item in row] for row in self])

    def __eq__(self, other):
        """Checks whether a given matrix is equal to another one"""
        return isinstance(other, Matrix) and \
                self._nrow == other._nrow and \
                self._ncol == other._ncol and \
                self._data == other._data

    def __getitem__(self, i):
        """Returns a single item, a row or a column of the matrix

        @param i: if a single integer, returns the M{i}th row as a list. If a
          slice, returns the corresponding rows as another L{Matrix} object. If
          a 2-tuple, the first element of the tuple is used to select a row and
          the second is used to select a column.
        """
        if isinstance(i, int):
            return list(self._data[i])
        elif isinstance(i, slice):
            return self.__class__(self._data[i])
        elif isinstance(i, tuple):
            try:
                first = i[0]
            except IndexError:
                first = slice(None)
            try:
                second = i[1]
            except IndexError:
                second = slice(None)
            if type(first) == slice and type(second) == slice:
                return self.__class__(row[second] for row in self._data[first])
            elif type(first) == slice:
                return [row[second] for row in self._data[first]]
            else:
                return self._data[first][second]
        else:
            raise IndexError("invalid matrix index")

    def __hash__(self):
        """Returns a hash value for a matrix."""
        return hash(self._nrow, self._ncol, self._data)

    def __iadd__(self, other):
        """In-place addition of a matrix or scalar."""
        if isinstance(other, Matrix):
            if self.shape != other.shape:
                raise ValueError("matrix shapes do not match")
            for row_a, row_b in izip(self._data, other):
                for i in xrange(len(row_a)):
                    row_a[i] += row_b[i]
        else:
            for row in self._data:
                for i in xrange(len(row)):
                    row[i] += other
        return self

    def __isub__(self, other):
        """In-place subtraction of a matrix or scalar."""
        if isinstance(other, Matrix):
            if self.shape != other.shape:
                raise ValueError("matrix shapes do not match")
            for row_a, row_b in izip(self._data, other):
                for i in xrange(len(row_a)):
                    row_a[i] -= row_b[i]
        else:
            for row in self._data:
                for i in xrange(len(row)):
                    row[i] -= other
        return self

    def __ne__(self, other):
        """Checks whether a given matrix is not equal to another one"""
        return not self == other

    def __setitem__(self, i, value):
        """Sets a single item, a row or a column of the matrix

        @param i: if a single integer, sets the M{i}th row as a list. If a
          slice, sets the corresponding rows from another L{Matrix} object.
          If a 2-tuple, the first element of the tuple is used to select a row
          and the second is used to select a column.
        @param value: the new value
        """
        if isinstance(i, int):
            # Setting a row
            if len(value) != len(self._data[i]):
                raise ValueError("new value must have %d items" % self._ncol)
            self._data[i] = list(value)
        elif isinstance(i, slice):
            # Setting multiple rows
            if len(value) != len(self._data[i]):
                raise ValueError("new value must have %d items" % self._ncol)
            if any(len(row) != self._ncol for row in value):
                raise ValueError("rows of new value must have %d items" % \
                        self._ncol)
            self._data[i] = [list(row) for row in value]
        elif isinstance(i, tuple):
            try:
                first = i[0]
            except IndexError:
                first = slice(None)
            try:
                second = i[1]
            except IndexError:
                second = slice(None)
            if type(first) == slice and type(second) == slice:
                # Setting a submatrix
                # TODO
                raise NotImplementedError
            elif type(first) == slice:
                # Setting a submatrix
                raise NotImplementedError
            else:
                # Setting a single element
                self._data[first][second] = value
        else:
            raise IndexError("invalid matrix index")

    def __sub__(self, other):
        """Subtracts the given value from the matrix.

        @param other: either a scalar or a matrix. Scalars will
          be subtracted from each element of the matrix. Matrices will
          be subtracted together elementwise.
        @return: the result matrix
        """
        if isinstance(other, Matrix):
            if self.shape != other.shape:
                raise ValueError("matrix shapes do not match")
            return self.__class__([
                [a-b for a, b in izip(row_a, row_b)]
                for row_a, row_b in izip(self, other)
            ])
        else:
            return self.__class__([
                [item-other for item in row] for row in self])

    def __repr__(self):
        class_name = self.__class__.__name__
        rows = ("[%s]" % ", ".join(repr(item) for item in row) for row in self)
        return "%s([%s])" % (class_name, ", ".join(rows))

    def __str__(self):
        rows = ("[%s]" % ", ".join(repr(item) for item in row) for row in self)
        return "[%s]" % "\n ".join(rows)

    def __iter__(self):
        """Support for iteration.

        This is actually implemented as a generator, so there is no need for a
        separate iterator class. The generator returns I{copies} of the rows in
        the matrix as lists to avoid messing around with the internals. Feel
        free to do anything with the copies, the changes won't be reflected in
        the original matrix."""
        return (list(row) for row in self._data)

    def __plot__(self, context, bbox, palette, **kwds):
        """Plots the matrix to the given Cairo context in the given box

        Besides the usual self-explanatory plotting parameters (C{context},
        C{bbox}, C{palette}), it accepts the following keyword arguments:

          - C{style}: the style of the plot. C{boolean} is useful for plotting
            matrices with boolean (C{True}/C{False} or 0/1) values: C{False}
            will be shown with a white box and C{True} with a black box.
            C{palette} uses the given palette to represent numbers by colors,
            the minimum will be assigned to palette color index 0 and the maximum
            will be assigned to the length of the palette. C{None} draws transparent
            cell backgrounds only. The default style is C{boolean} (but it may
            change in the future). C{None} values in the matrix are treated
            specially in both cases: nothing is drawn in the cell corresponding
            to C{None}.

          - C{square}: whether the cells of the matrix should be square or not.
            Default is C{True}.

          - C{grid_width}: line width of the grid shown on the matrix. If zero or
            negative, the grid is turned off. The grid is also turned off if the size
            of a cell is less than three times the given line width. Default is C{1}.
            Fractional widths are also allowed.

          - C{border_width}: line width of the border drawn around the matrix.
            If zero or negative, the border is turned off. Default is C{1}.

          - C{row_names}: the names of the rows

          - C{col_names}: the names of the columns.

          - C{values}: values to be displayed in the cells. If C{None} or
            C{False}, no values are displayed. If C{True}, the values come
            from the matrix being plotted. If it is another matrix, the
            values of that matrix are shown in the cells. In this case,
            the shape of the value matrix must match the shape of the
            matrix being plotted.

          - C{value_format}: a format string or a callable that specifies how
            the values should be plotted. If it is a callable, it must be a
            function that expects a single value and returns a string.
            Example: C{"%#.2f"} for floating-point numbers with always exactly
            two digits after the decimal point. See the Python documentation of
            the C{%} operator for details on the format string. If the format
            string is not given, it defaults to the C{str} function.

        If only the row names or the column names are given and the matrix
        is square-shaped, the same names are used for both column and row
        names.
        """
        # pylint: disable-msg=W0142
        # pylint: disable-msg=C0103
        grid_width = float(kwds.get("grid_width", 1.))
        border_width = float(kwds.get("border_width", 1.))
        style = kwds.get("style", "boolean")
        row_names = kwds.get("row_names")
        col_names = kwds.get("col_names", row_names)
        values = kwds.get("values")
        value_format = kwds.get("value_format", str)

        # Validations
        if style not in ("boolean", "palette", "none", None):
            raise ValueError("invalid style")
        if style == "none":
            style = None
        if row_names is None and col_names is not None:
            row_names = col_names
        if row_names is not None:
            row_names = [str(name) for name in islice(row_names, self._nrow)]
            if len(row_names) < self._nrow:
                row_names.extend([""]*(self._nrow-len(row_names)))
        if col_names is not None:
            col_names = [str(name) for name in islice(col_names, self._ncol)]
            if len(col_names) < self._ncol:
                col_names.extend([""]*(self._ncol-len(col_names)))
        if values == False:
            values = None
        if values == True:
            values = self
        if isinstance(values, list):
            values = Matrix(list)
        if values is not None and not isinstance(values, Matrix):
            raise TypeError("values must be None, False, True or a matrix")
        if values is not None and values.shape != self.shape:
            raise ValueError("values must be a matrix of size %s" % self.shape)

        # Calculate text extents if needed
        if row_names is not None or col_names is not None:
            te = context.text_extents
            space_width = te(" ")[4]
            max_row_name_width = max([te(s)[4] for s in row_names])+space_width
            max_col_name_width = max([te(s)[4] for s in col_names])+space_width
        else:
            max_row_name_width, max_col_name_width = 0, 0

        # Calculate sizes
        total_width = float(bbox.width)-max_row_name_width
        total_height = float(bbox.height)-max_col_name_width
        dx = total_width / self.shape[1]
        dy = total_height / self.shape[0]
        if kwds.get("square", True):
            dx, dy = min(dx, dy), min(dx, dy)
        total_width, total_height = dx*self.shape[1], dy*self.shape[0]
        ox = bbox.left + (bbox.width - total_width - max_row_name_width) / 2.0
        oy = bbox.top + (bbox.height - total_height - max_col_name_width) / 2.0
        ox += max_row_name_width
        oy += max_col_name_width

        # Determine rescaling factors for the palette if needed
        if style == "palette":
            mi, ma = self.min(), self.max()
            color_offset = mi
            color_ratio = (len(palette)-1) / float(ma-mi)

        # Validate grid width
        if dx < 3*grid_width or dy < 3*grid_width:
            grid_width = 0.
        if grid_width > 0:
            context.set_line_width(grid_width)
        else:
            # When the grid width is zero, we will still stroke the
            # rectangles, but with the same color as the fill color
            # of the cell - otherwise we would get thin white lines
            # between the cells as a drawing artifact
            context.set_line_width(1)

        # Draw row names (if any)
        context.set_source_rgb(0., 0., 0.)
        if row_names is not None:
            x, y = ox, oy 
            for heading in row_names:
                _, _, _, h, xa, _ = context.text_extents(heading)
                context.move_to(x-xa-space_width, y + (dy+h)/2.)
                context.show_text(heading)
                y += dy

        # Draw column names (if any)
        if col_names is not None:
            context.save()
            context.translate(ox, oy)
            context.rotate(-1.5707963285)   # pi/2
            x, y = 0., 0.
            for heading in col_names:
                _, _, _, h, _, _ = context.text_extents(heading)
                context.move_to(x+space_width, y + (dx+h)/2.)
                context.show_text(heading)
                y += dx
            context.restore()

        # Draw matrix
        x, y = ox, oy
        if style is None:
            fill = lambda: None
        else:
            fill = context.fill_preserve
        for row in self:
            for item in row:
                if item is None:
                    x += dx
                    continue
                if style == "boolean":
                    if item:
                        context.set_source_rgb(0., 0., 0.)
                    else:
                        context.set_source_rgb(1., 1., 1.)
                elif style == "palette":
                    cidx = int((item-color_offset)*color_ratio)
                    if cidx < 0:
                        cidx = 0
                    context.set_source_rgba(*palette.get(cidx))
                context.rectangle(x, y, dx, dy)
                if grid_width > 0:
                    fill()
                    context.set_source_rgb(0.5, 0.5, 0.5)
                    context.stroke()
                else:
                    fill()
                    context.stroke()
                x += dx
            x, y = ox, y+dy

        # Draw cell values
        if values is not None:
            x, y = ox, oy
            context.set_source_rgb(0., 0., 0.)
            for row in values.data:
                if hasattr(value_format, "__call__"):
                    values = [value_format(item) for item in row]
                else:
                    values = [value_format % item for item in row]
                for item in values:
                    th, tw = context.text_extents(item)[3:5]
                    context.move_to(x+(dx-tw)/2., y+(dy+th)/2.)
                    context.show_text(item)
                    x += dx
                x, y = ox, y+dy

        # Draw borders
        if border_width > 0:
            context.set_line_width(border_width)
            context.set_source_rgb(0., 0., 0.)
            context.rectangle(ox, oy, dx*self.shape[1], dy*self.shape[0])
            context.stroke()


    def min(self, dim=None):
        """Returns the minimum of the matrix along the given dimension

        @param dim: the dimension. 0 means determining the column minimums, 1 means
          determining the row minimums. If C{None}, the global minimum is
          returned.
        """
        if dim == 1:
            return [min(row) for row in self._data]
        if dim == 0:
            return [min(row[idx] for row in self._data) \
                        for idx in xrange(self._ncol)]
        return min(min(row) for row in self._data)

    def max(self, dim=None):
        """Returns the maximum of the matrix along the given dimension

        @param dim: the dimension. 0 means determining the column maximums, 1 means
          determining the row maximums. If C{None}, the global maximum is
          returned.
        """
        if dim == 1:
            return [max(row) for row in self._data]
        if dim == 0:
            return [max(row[idx] for row in self._data) \
                        for idx in xrange(self._ncol)]
        return max(max(row) for row in self._data)


class DyadCensus(tuple):
    """Dyad census of a graph.

    This is a pretty simple class - basically it is a tuple, but it allows
    the user to refer to its individual items by the names C{mutual} (or
    C{mut}), C{asymmetric} (or C{asy} or C{asym} or C{asymm}) and C{null}.

    Examples:

      >>> from igraph import Graph
      >>> g=Graph.Erdos_Renyi(100, 0.2, directed=True)
      >>> dc=g.dyad_census()
      >>> print dc.mutual             #doctest:+SKIP
      179
      >>> print dc["asym"]            #doctest:+SKIP
      1609
      >>> print tuple(dc), list(dc)   #doctest:+SKIP
      (179, 1609, 3162) [179, 1609, 3162]
      >>> print sorted(dc.as_dict().items())   #doctest:+ELLIPSIS
      [('asymmetric', ...), ('mutual', ...), ('null', ...)]

    @undocumented: _remap
    """
    _remap = {"mutual": 0, "mut": 0, "sym": 0, "symm": 0,
        "asy": 1, "asym": 1, "asymm": 1, "asymmetric": 1, "null": 2}

    def __getitem__(self, idx):
        return tuple.__getitem__(self, self._remap.get(idx, idx))

    def __getattr__(self, attr):
        if attr in self._remap:
            return tuple.__getitem__(self, self._remap[attr])
        raise AttributeError("no such attribute: %s" % attr)

    def __repr__(self):
        return "DyadCensus((%d, %d, %d))" % self

    def __str__(self):
        return "%d mutual, %d asymmetric, %d null dyads" % self

    def as_dict(self):
        """Converts the dyad census to a dict using the known dyad names."""
        return {"mutual": self[0], "asymmetric": self[1], "null": self[2]}


class TriadCensus(tuple):
    """Triad census of a graph.

    This is a pretty simple class - basically it is a tuple, but it allows
    the user to refer to its individual items by the following triad names:

      - C{003} -- the empty graph
      - C{012} -- a graph with a single directed edge (C{A --> B, C})
      - C{102} -- a graph with a single mutual edge (C{A <-> B, C})
      - C{021D} -- the binary out-tree (C{A <-- B --> C})
      - C{021U} -- the binary in-tree (C{A --> B <-- C})
      - C{021C} -- the directed line (C{A --> B --> C})
      - C{111D} -- C{A <-> B <-- C}
      - C{111U} -- C{A <-> B --> C}
      - C{030T} -- C{A --> B <-- C, A --> C}
      - C{030C} -- C{A <-- B <-- C, A --> C}
      - C{201} -- C{A <-> B <-> C}
      - C{120D} -- C{A <-- B --> C, A <-> C}
      - C{120U} -- C{A --> B <-- C, A <-> C}
      - C{120C} -- C{A --> B --> C, A <-> C}
      - C{210C} -- C{A --> B <-> C, A <-> C}
      - C{300} -- the complete graph (C{A <-> B <-> C, A <-> C})

    Attribute and item accessors are provided. Due to the syntax of Python,
    attribute names are not allowed to start with a number, therefore the
    triad names must be prepended with a lowercase C{t} when accessing
    them as attributes. This is not necessary with the item accessor syntax.

    Examples:

      >>> from igraph import Graph
      >>> g=Graph.Erdos_Renyi(100, 0.2, directed=True)
      >>> tc=g.triad_census()
      >>> print tc.t003                     #doctest:+SKIP
      39864
      >>> print tc["030C"]                  #doctest:+SKIP
      1206
    """
    _remap = {"003": 0, "012": 1, "102": 2, "021D": 3, "021U": 4, "021C": 5, \
        "111D": 6, "111U": 7, "030T": 8, "030C": 9, "201": 10, "120D": 11, \
        "120U": 12, "120C": 13, "210": 14, "300": 15}

    def __getitem__(self, idx):
        if isinstance(idx, basestring):
            idx = idx.upper()
        return tuple.__getitem__(self, self._remap.get(idx, idx))

    def __getattr__(self, attr):
        if isinstance(attr, basestring) and attr[0] == 't' \
                and attr[1:].upper() in self._remap:
            return tuple.__getitem__(self, self._remap[attr[1:].upper()])
        raise AttributeError("no such attribute: %s" % attr)

    def __repr__(self):
        return "TriadCensus((%s))" % ", ".join(str(item) for item in self)

    def __str__(self):
        maxidx = len(self)
        maxcount = max(self)
        numwidth = len(str(maxcount))
        captionwidth = max(len(key) for key in self._remap)
        colcount = 4

        rowcount = maxidx / colcount
        if rowcount * colcount < maxidx:
            rowcount += 1

        invmap = dict((v, k) for k, v in self._remap.iteritems())
        result, row, idx = [], [], 0
        for _ in xrange(rowcount):
            for _ in xrange(colcount):
                if idx >= maxidx:
                    break 
                row.append("%-*s: %*d" % (captionwidth, invmap.get(idx, ""),
                  numwidth, self[idx]))
                idx += 1
            result.append(" | ".join(row))
            row = []

        return "\n".join(result)


class UniqueIdGenerator(object):
    """A dictionary-like class that can be used to assign unique IDs to
    names (say, vertex names).

    Usage:
    
    >>> gen = UniqueIdGenerator()
    >>> gen["A"]
    0
    >>> gen["B"]
    1
    >>> gen["C"]
    2
    >>> gen["A"]      # Retrieving already existing ID
    0
    >>> gen.add("D")  # Synonym of gen["D"]
    3
    >>> len(gen)      # Number of already used IDs
    4
    >>> "C" in gen
    True
    >>> "E" in gen
    False
    """

    def __init__(self, id_generator=None, initial=None):
        """Creates a new unique ID generator. `id_generator` specifies how do we
        assign new IDs to elements that do not have an ID yet. If it is `None`,
        elements will be assigned integer identifiers starting from 0. If it is
        an integer, elements will be assigned identifiers starting from the given
        integer. If it is an iterator or generator, its `next` method will be
        called every time a new ID is needed."""
        if id_generator is None:
            id_generator = 0
        if isinstance(id_generator, int):
            import itertools
            self._generator = itertools.count(id_generator)
        else:
            self._generator = id_generator
        self._ids = {}
        if initial:
            for value in initial:
                self.add(value)

    def __contains__(self, item):
        """Checks whether `item` already has an ID or not."""
        return item in self._ids

    def __getitem__(self, item):
        """Retrieves the ID corresponding to `item`. Generates a new ID for
        `item` if it is the first time we request an ID for it."""
        try:
            return self._ids[item]
        except KeyError:
            self._ids[item] = self._generator.next()
            return self._ids[item]

    def __setitem__(self, item, value):
        """Overrides the ID for `item`."""
        self._ids[item] = value

    def __len__(self):
        """"Returns the number of items"""
        return len(self._ids)

    def reverse_dict(self):
        """Returns the reverse mapping, i.e., the one that maps from generated
        IDs to their corresponding objects"""
        return dict((v, k) for k, v in self._ids.iteritems())

    def values(self):
        """Returns the values stored so far. If the generator generates items
        according to the standard sorting order, the values returned will be
        exactly in the order they were added. This holds for integer IDs for
        instance (but for many other ID generators as well)."""
        return sorted(self._ids.keys(), key = self._ids.__getitem__)

    add = __getitem__


