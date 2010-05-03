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
import math

__all__ = ["Histogram", "RunningMean", "power_law_fit"]

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

    def __init__(self, bin_width = 1, data = None):
        """Initializes the histogram with the given data set.

        @param bin_width: the bin width of the histogram.
        @param data: the data set to be used. Must contain real numbers.
        """
        self._bin_width = float(bin_width)
        self._bins = None
        self._min, self._max = None, None
        self._running_mean = RunningMean()
        self.clear()

        if data:
            self.add_many(data)

    def _get_bin(self, num, create = False):
        """Returns the bin index corresponding to the given number.

        @param num: the number for which the bin is being sought
        @param create: whether to create a new bin if no bin exists yet.
        @return: the index of the bin or C{None} if no bin exists yet and
          {create} is C{False}."""
        if len(self._bins) == 0:
            if not create:
                result = None
            else: 
                self._min = int(num/self._bin_width)*self._bin_width
                self._max = self._min+self._bin_width
                self._bins = [0]
                result = 0
            return result

        if num >= self._min:
            binidx = int((num-self._min)/self._bin_width)
            if binidx < len(self._bins):
                return binidx
            if not create:
                return None
            extra_bins = binidx-len(self._bins)+1
            self._bins.extend([0]*extra_bins)
            self._max = self._min + len(self._bins)*self._bin_width
            return binidx

        if not create:
            return None

        extra_bins = int(math.ceil((self._min-num)/self._bin_width))
        self._bins[0:0] = [0]*extra_bins
        self._min -= extra_bins*self._bin_width
        self._max = self._min + len(self._bins)*self._bin_width
        return 0

    @property
    def n(self):
        """Returns the number of elements in the histogram"""
        return len(self._running_mean)

    @property
    def mean(self):
        """Returns the mean of the elements in the histogram"""
        return self._running_mean.mean

    # pylint: disable-msg=C0103
    @property
    def sd(self):
        """Returns the standard deviation of the elements in
        the histogram"""
        return self._running_mean.sd

    @property
    def var(self):
        """Returns the variance of the elements in the histogram"""
        return self._running_mean.var

    def add(self, num, repeat=1):
        """Adds a single number to the histogram.
        
        @param num: the number to be added
        @param repeat: number of repeated additions
        """
        num = float(num)
        binidx = self._get_bin(num, True)
        self._bins[binidx] += repeat 
        self._running_mean.add(num, repeat)

    def add_many(self, data):
        """Adds a single number or the elements of an iterable to the histogram.

        @param data: the data to be added"""
        try:
            iterator = iter(data)
        except TypeError:
            iterator = iter([data])
        for x in iterator:
            self.add(x)
    __lshift__ = add_many

    def clear(self):
        """Clears the collected data"""
        self._bins = []
        self._min, self._max = None, None
        self._running_mean = RunningMean()

    def bins(self):
        """Generator returning the bins of the histogram in increasing order
        
        @return: a tuple with the following elements: left bound, right bound,
          number of elements in the bin"""
        x = self._min
        for elem in self._bins:
            yield (x, x+self._bin_width, elem)
            x += self._bin_width

    def __plot__(self, context, bbox, _, **kwds):
        """Plotting support"""
        import igraph.drawing
        coord_system = igraph.drawing.DescartesCoordinateSystem(context, bbox, \
            (kwds.get("min", self._min), 0, \
             kwds.get("max", self._max), kwds.get("max_value", max(self._bins))
            ))

        # Draw the boxes
        context.set_line_width(1)
        context.set_source_rgb(1., 0., 0.)
        x = self._min
        for value in self._bins:
            top_left_x, top_left_y = coord_system.local_to_context(x, value)
            x += self._bin_width
            bottom_right_x, bottom_right_y = coord_system.local_to_context(x, 0)
            context.rectangle(top_left_x, top_left_y, \
                              bottom_right_x - top_left_x, \
                              bottom_right_y - top_left_y)
            context.fill()

        # Draw the axes
        coord_system.draw()

    def __str__(self):
        """Returns the string representation of the histogram"""
        if self._min is None or self._max is None:
            return "N = 0"
        num_length = max(len("%.3f" % self._min), \
                         len("%.3f" % self._max))
        format_string = "[%%%d.3f, %%%d.3f): %%s" % (num_length, num_length)
        maxval = max(self._bins)
        scale = maxval // (70-2*num_length)
        if scale < 1:
            scale = 1

        result = ["N = %d, mean +- sd: %.4f +- %.4f " % \
            (self.n, self.mean, self.sd)]

        if scale > 1:
            result.append("Each * represents %d items" % scale)

        for left, right, cnt in self.bins():
            cnt //= scale
            result.append(format_string % (left, right, '*'*cnt))

        return "\n".join(result)



class RunningMean(object):
    """Running mean calculator.
    
    This class can be used to calculate the mean of elements from a
    list, tuple, iterable or any other data source. The mean is
    calculated on the fly without explicitly summing the values,
    so it can be used for data sets with arbitrary item count. Also
    capable of returning the standard deviation (also calculated on
    the fly)
    """

    # pylint: disable-msg=C0103
    def __init__(self, n = 0.0, mean = 0.0, sd = 0.0):
        """RunningMean(n=0.0, mean=0.0, sd=0.0)
        
        Initializes the running mean calculator. Optionally the
        number of already processed elements and an initial mean
        can be supplied if we want to continue an interrupted
        calculation.

        @param n: the initial number of elements already processed
        @param mean: the initial mean
        @param sd: the initial standard deviation"""
        self._nitems = float(n)
        self._mean = float(mean)
        if n > 1:
            self._sqdiff = float(sd) ** 2 * float(n-1)
            self._sd = float(sd)
        else:
            self._sqdiff = 0.0
            self._sd = 0.0
        
    def add(self, value, repeat=1):
        """RunningMean.add(value, repeat=1)
        
        Adds the given value to the elements from which we calculate
        the mean and the standard deviation.

        @param value: the element to be added
        @param repeat: number of repeated additions
        @return: the new mean and standard deviation as a tuple"""
        repeat = int(repeat)
        self._nitems += repeat
        delta = value - self._mean
        self._mean += (repeat*delta / self._nitems)
        self._sqdiff += (repeat*delta) * (value - self._mean)
        if self._nitems > 1:
            self._sd = (self._sqdiff / (self._nitems-1)) ** 0.5
        return self._mean, self._sd

    def add_many(self, values):
        """RunningMean.add(values)
        
        Adds the values in the given iterable to the elements from
        which we calculate the mean. Can also accept a single number.
        The left shift (C{<<}) operator is aliased to this function,
        so you can use it to add elements as well:
            
          >>> rm=RunningMean()
          >>> rm << [1,2,3,4]
          (2.5, 1.6666666666667)
        
        @param values: the element(s) to be added
        @type values: iterable
        @return: the new mean"""
        try:
            iterator = iter(values)
        except TypeError:
            iterator = iter([values])
        for value in iterator:
            self.add(value)
        return self._mean, self._sd

    @property
    def result(self):
        """Returns the current mean and standard deviation as a tuple"""
        return self._mean, self._sd

    @property
    def mean(self):
        """Returns the current mean"""
        return self._mean

    @property
    def sd(self):
        """Returns the current standard deviation"""
        return self._sd

    @property
    def var(self):
        """Returns the current variation"""
        return self._sd ** 2

    def __str__(self):
        return "Running mean (N=%d, %f +- %f)" % \
            (self._nitems, self._mean, self._sd)
    
    __lshift__ = add_many
    
    def __float__(self):
        return float(self._mean)

    def __int__(self):
        return int(self._mean)

    def __long__(self):
        return long(self._mean)

    def __complex__(self):
        return complex(self._mean)

    def __len__(self):
        return self._nitems

def power_law_fit(x, xmin=None, method="discrete_approx"):
    """Fitting a power-law distribution to empirical data

    @param x: the data to fit, a list containing integer values
    @param xmin: the lower bound for fitting the power-law. If C{None}, the
      smallest x value is used. This argument makes it possible to fit
      only the tail of the distribution.
    @param method: the fitting method to use. The following methods are
      implemented so far:

        - C{continuous}, C{hill}: exact maximum likelihood estimation
          when the input data comes from a continuous scale. This is
          known as the Hill estimator. The statistical error of
          this estimator is M{(alpha-1) / sqrt(n)}, where alpha is the
          estimated exponent and M{n} is the number of data points above
          M{xmin}. The estimator is known to exhibit a small finite
          sample-size bias of order M{O(n^-1)}, which is small when
          M{n > 100}.

        - C{discrete_approx}: approximation of the maximum likelihood
          estimation in discrete case (see Clauset et al among the
          references). This is said to produce quite results provided
          M{xmin} >= 6 (approx.).

    @return: the estimated power-law exponent
    
    @newfield ref: Reference
    @ref: MEJ Newman: Power laws, Pareto distributions and Zipf's law.
      Contemporary Physics 46, 323-351 (2005)
    @ref: A Clauset, CR Shalizi, MEJ Newman: Power-law distributions
      in empirical data. E-print (2007). arXiv:0706.1062"""
    real_xmin = float(min(x))
    if xmin is not None:
        real_xmin = float(max(xmin, real_xmin))
    filtered_xs = [x1 for x1 in x if x1 >= xmin]

    method = method.lower()

    if method in ("continuous", "hill", "discrete_approx"):
        if method == "discrete_approx":
            real_xmin -= 0.5
        log_sum = sum(math.log(x/real_xmin) for x in filtered_xs)
        if log_sum == 0:
            raise ValueError("lower bound too high")
        return 1.0+len(filtered_xs)/log_sum

    raise ValueError("unknown method: %s" % method)

