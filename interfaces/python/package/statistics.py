"""
Statistics related stuff in igraph
"""

__license__ = """
Copyright (C) 2006-2007  Gabor Csardi <csardi@rmki.kfki.hu>,
Tamas Nepusz <ntamas@rmki.kfki.hu>

MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary

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

class Interval:
    """A class representing an interval over the real numbers"""

    def __init__(self, left, right):
        """Constructor.

        @param left: the left edge of the interval
        @param right: the right edge of the interval"""
        if left == right:
            raise ArgumentError, "Interval can't have zero length"""
        self._left = min(left, right)
        self._right = max(left, right)

    def __contains__(self, x):
        """Returns C{True} if x is in the interval, C{False} otherwise"""
        return (x >= self._left and x < self._right)

    def __eq__(self, x):
        """Tests for the equality of two intervals"""
        return (isinstance(x, Interval) and \
                x.left == self._left and x.right == self._right)
 
    def __cmp__(self, x): return cmp(self._left, x.left)
    def __getattr__(self, attr):
        if attr == "left": return self._left
        if attr == "right": return self._right
        return object.__getattr__(self, attr)

    def __hash__(self): return hash(self._left) | hash(self._right)

    def __str__(self): return "[%f, %f)" % (self._left, self._right)

class Histogram(object):
    """Generic histogram class for real numbers
    
    Example:
        
        >>> h = Histogram(5)     # Initializing, bin width = 5
        >>> h << [2,3,2,7,8,5,5,0,7,9]     # Adding more items
        >>> print h
        N = 10, mean +- sd: 4.8000 +- 2.9740
        [ 0.000,  5.000): ****
        [ 5.000, 10.000): ******
    """

    def __init__(self, bin_width = 1, data = []):
        """Initializes the histogram with the given data set.

        @param bin_width: the bin width of the histogram.
        @param data: the data set to be used. Must contain real numbers.
        """
        self._bin_width = float(bin_width)
        self.clear()
        self.add_many(data)

    def _get_bin(self, num, create = False):
        """Returns the bin corresponding to the given number.

        @param num: the number for which the bin is being sought
        @param create: whether to create a new bin if no bin exists yet.
        @return: the range of the bin or C{None} if no bin exists yet and
          {create} is C{False}."""
        for bin in self._bins.keys():
            if num in bin: return bin

        if create:
            left = (num // self._bin_width) * self._bin_width
            right = left + self._bin_width
            bin = Interval(left, right)
            self._bins[bin] = 0

            if self._min is None or left < self._min: self._min = left
            if self._max is None or right > self._max: self._max = right

            return bin

        return None

    def _get_n(self): return self._n
    def _get_mean(self): return self._mean
    def _get_sd(self): return self._sd
    def _get_var(self): return self._sd ** 2
    n = property(_get_n, doc="Number of elements in the histogram")
    mean = property(_get_mean, doc="Mean of the elements")
    sd = property(_get_sd, doc="Standard deviation of the elements")
    var = property(_get_var, doc="Variance of the elements")

    def add(self, num):
        """Adds a single number to the histogram.
        
        @param num: the number to be added"""
        bin = self._get_bin(num, True)
        self._bins[bin] += 1
        self._n += 1
        delta = num - self._mean
        self._mean += (delta / self._n)
        self._s += delta * (num - self._mean)
        if self._n > 1:
            self._sd = (self._s / (self._n-1)) ** 0.5

    def add_many(self, data):
        """Adds a single number or elements of an iterable to the histogram.

        @param data: the data to be added"""
        try:
            it = iter(data)
        except:
            it = iter([data])
        for x in it: self.add(x)
    __lshift__ = add_many

    def clear(self):
        """Clears the collected data"""
        self._bins = {}
        self._min = None
        self._max = None
        self._n = 0
        self._mean = 0.0
        self._s = 0.0
        self._sd = 0.0

    def bins(self):
        """Generator returning the bins of the histogram in increasing order
        
        @return: a tuple with the following elements: left bound, right bound,
          number of elements in the bin"""
        x = self._min
        while x < self._max:
            bin = self._get_bin(x)
            if bin is None:
                yield (x, x+self._bin_width, 0)
            else:
                yield (x, x+self._bin_width, self._bins[bin])
            x += self._bin_width

    def __str__(self):
        """Returns the string representation of the histogram"""
        if self._min is None or self._max is None: return str()
        num_length = max(len("%.3f" % self._min), \
                         len("%.3f" % self._max))
        format_string = "[%%%d.3f, %%%d.3f): %%s" % (num_length, num_length)
        #bins = self._bins.keys()
        #bins.sort()
        maxval = max(self._bins.itervalues())
        scale = maxval // (70-2*num_length)
        if scale<1: scale = 1

        result=["N = %d, mean +- sd: %.4f +- %.4f " % \
            (self.n, self.mean, self.sd)]

        if scale>1: result.append("Each * represents %d items" % scale)

        for left, right, cnt in self.bins():
            cnt //= scale
            result.append(format_string % (left, right, '*'*cnt))

        return "\n".join(result)


