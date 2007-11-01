"""Additional auxiliary data types"""

class Matrix(object):
    """Simple matrix data type.

    Of course there are much more advanced matrix data types for Python (for instance,
    the C{ndarray} data type of Numeric Python) and this implementation does not want
    to compete with them. The only role of this data type is to provide a convenient
    interface for the matrices returned by the C{Graph} object (for instance, allow
    indexing with tuples in the case of adjacency matrices and so on).
    """

    def __init__(self, data=None):
        """Initializes a matrix.

        @param data: the elements of the matrix as a list of lists, or C{None} to
          create a 0x0 matrix.
        """
        self._data = []
        self.data = data

    def _set_data(self, data=None):
        """Sets the data stored in the matrix"""
        if data is not None:
            self._data = map(list, data)
            self._nrow = len(data)
            self._ncol = max(map(len, data))
            for row in self._data:
                if len(row) < self._ncol: row.extend([0]*(self._ncol-len(row)))

    def _get_data(self):
        """Returns the data stored in the matrix as a list of lists"""
        return [list(row) for row in self._data]
    data = property(_get_data, _set_data, doc="Elements of the matrix as a list of lists")

    def _get_shape(self):
        """Returns the shape of the matrix"""
        return self._nrow, self._ncol
    shape = property(_get_shape, doc="Shape of the matrix as a tuple")

    def __getitem__(self, i):
        """Returns a single item, a row or a column of the matrix

        @param i: if a single integer, returns the M{i}th row as a list. If a slice,
          returns the corresponding rows as another L{Matrix} object. If a 2-tuple,
          the first element of the tuple is used to select a row and the second is
          used to select a column.
        """
        if isinstance(i, int):
            return list(self._data[i])
        elif isinstance(i, slice):
            return self.__class__(self._data[i])
        elif isinstance(i, tuple):
            try:
                first = i[0]
            except:
                first = slice(None)
            try:
                second = i[1]
            except:
                second = slice(None)
            if type(first) == slice and type(second) == slice:
                return self.__class__([row[second] for row in self._data[first]])
            elif type(first) == slice:
                return [row[second] for row in self._data[first]]
            else:
                return self._data[first][second]
        else:
            raise IndexError, "invalid matrix index"


    def __setitem__(self, i, value):
        """Sets a single item, a row or a column of the matrix

        @param i: if a single integer, returns the M{i}th row as a list. If a slice,
          returns the corresponding rows as another L{Matrix} object. If a 2-tuple,
          the first element of the tuple is used to select a row and the second is
          used to select a column.
        @param value: the new value
        """
        if isinstance(i, int):
            # Setting a row
            if len(value) != len(self._data[i]):
                raise ValueError, "new value must have %d items" % len(self._data[i])
            self._data[i] = list(value)
        elif isinstance(i, slice):
            # Setting multiple rows
            if len(value) != len(self._data[i]):
                raise ValueError, "new value must have %d items" % len(self._data[i])
            for j in xrange(len(value)):
                if len(value[j]) != len(self._data[0]):
                    raise ValueError, "rows of new value must have %d items" % len(self._data[0])
            self._data[i] = list(map(list, value))
        elif isinstance(i, tuple):
            try:
                first = i[0]
            except:
                first = slice(None)
            try:
                second = i[1]
            except:
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
            raise IndexError, "invalid matrix index"


    def __repr__(self):
        return "Matrix([%s])" % "\n ".join(["[%s]" % (", ".join(map(repr, row))) for row in self])

    def __str__(self):
        return "[%s]" % "\n ".join(["[%s]" % (", ".join(map(str, row))) for row in self])

    def __iter__(self):
        """Support for iteration.

        This is actually implemented as a generator, so there is no need for a
        separate iterator class. The generator returns I{copies} of the rows in
        the matrix as lists to avoid messing around with the internals. Feel
        free to do anything with the copies, the changes won't be reflected in
        the original matrix."""
        for row in self._data: yield list(row)

    def __plot__(self, context, bbox, palette, *args, **kwds):
        """Plots the matrix to the given Cairo context in the given box

        Besides the usual self-explanatory plotting parameters (C{context},
        C{bbox}, C{palette}), it accepts the following keyword arguments:

          - C{style}: the style of the plot. C{boolean} is useful for plotting
            matrices with boolean (C{True}/C{False} or 0/1) values: C{False}
            will be shown with a white box and C{True} with a black box.
            C{palette} uses the given palette to represent numbers by colors,
            the minimum will be assigned to palette color index 0 and the maximum
            will be assigned to the length of the palette. The default style is
            C{boolean} (but it may change in the future). C{None} values in
            the matrix are treated specially in both cases: nothing is drawn in
            the cell corresponding to C{None}.

          - C{square}: whether the cells of the matrix should be square or not.
            Default is C{True}.

          - C{grid_width}: line width of the grid shown on the matrix. If zero or
            negative, the grid is turned off. The grid is also turned off if the size
            of a cell is less than three times the given line width. Default is C{1}.
            Fractional widths are also allowed.

          - C{border_width}: line width of the border shown around the matrix.
            If zero or negative, the border is turned off. Default is C{1}.
        """
        import colors
        grid_width = float(kwds.get("grid_width", 1.))
        border_width = float(kwds.get("border_width", 1.))
        style = kwds.get("style", "boolean")
        if style not in ("boolean", "palette"): raise ValueError, "invalid style"

        # Calculate sizes
        dx = float(bbox.width) / self.shape[1]
        dy = float(bbox.height) / self.shape[0]
        if kwds.get("square", True): dx, dy = min(dx, dy), min(dx, dy)
        ox = bbox.left + (bbox.width - dx*self.shape[1]) / 2.0
        oy = bbox.top + (bbox.height - dy*self.shape[0]) / 2.0

        # Determine rescaling factors for the palette if needed
        if style == "palette":
            mi = self.min()
            ma = self.max()
            color_offset = mi
            color_ratio = (len(palette)-1) / float(ma-mi)

        # Validate grid width
        if dx < 3*grid_width or dy < 3*grid_width: grid_width = 0.
        if grid_width>0: context.set_line_width(grid_width)

        x, y = ox, oy
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
                    if cidx < 0: cidx = 0
                    context.set_source_rgb(*palette.get(cidx))
                context.rectangle(x, y, dx, dy)
                if grid_width>0:
                    context.fill_preserve()
                    context.set_source_rgb(0.5, 0.5, 0.5)
                    context.stroke()
                else:
                    context.fill()
                x += dx
            x, y = ox, y+dy

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
        if dim == 1: return map(min, self._data)
        if dim == 0: return map(min, [[row[idx] for row in self._data] for idx in xrange(self._ncol)])
        return min(map(min, self._data))

    def max(self, dim=None):
        """Returns the maximum of the matrix along the given dimension

        @param dim: the dimension. 0 means determining the column maximums, 1 means
          determining the row maximums. If C{None}, the global maximum is
          returned.
        """
        if dim == 1: return map(max, self._data)
        if dim == 0: return map(max, [[row[idx] for row in self._data] for idx in xrange(self._ncol)])
        return max(map(max, self._data))


class DyadCensus(tuple):
    """Dyad census of a graph.

    This is a pretty simple class - basically it is a tuple, but it allows
    the user to refer to its individual items by the names C{mutual} (or
    C{mut}), C{asymmetric} (or C{asy} or C{asym} or C{asymm}) and C{null}.

    Examples:

      >>> g=Graph.Erdos_Renyi(100, 0.2, directed=True)
      >>> dc=g.dyad_census()
      >>> print dc.mutual
      179
      >>> print dc["asym"]
      1609
      >>> print tuple(dc), list(dc)
      (179, 1609, 3162) [179, 1609, 3162]
      >>> print dc.as_dict()
      {"mutual": 179, "asymmetric": 1609, "null": 3162}
    """
    _remap = {"mutual": 0, "mut": 0, "sym": 0, "symm": 0,
        "asy": 1, "asym": 1, "asymm": 1, "asymmetric": 1, "null": 2}

    def __getitem__(self, idx):
        return tuple.__getitem__(self, self._remap.get(idx, idx))

    def __getattr__(self, attr):
        if attr in self._remap.keys():
            return tuple.__getitem__(self, self._remap[attr])
        raise AttributeError, "no such attribute: %s" % attr

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

      >>> g=Graph.Erdos_Renyi(100, 0.2, directed=True)
      >>> tc=g.triad_census()
      >>> print tc.t003
      39864
      >>> print tc["030C"]
      1206
    """
    _remap = {"003": 0, "012": 1, "102": 2, "021D": 3, "021U": 4, "021C": 5, \
        "111D": 6, "111U": 7, "030T": 8, "030C": 9, "201": 10, "120D": 11, \
        "120U": 12, "120C": 13, "210": 14, "300": 15}

    def __getitem__(self, idx):
        if isinstance(idx, basestring): idx=idx.upper()
        return tuple.__getitem__(self, self._remap.get(idx, idx))

    def __getattr__(self, attr):
        if isinstance(attr, basestring) and attr[0] == 't' \
            and attr[1:].upper() in self._remap.keys():
                return tuple.__getitem__(self, self._remap[attr[1:].upper()])
        raise AttributeError, "no such attribute: %s" % attr

    def __repr__(self):
        return "TriadCensus((%s))" % ", ".join(map(str, self))

    def __str__(self):
        # TODO
        return "%d mutual, %d asymmetric, %d null dyads" % self


