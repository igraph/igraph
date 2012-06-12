/* 
   IGraph library Java interface.
   Copyright (C) 2006-2012  Tamas Nepusz <ntamas@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA 

*/

/*

ATTENTION: This is a highly experimental, proof-of-concept Java interface.
Its main purpose was to convince me that it can be done in finite time :)
The interface is highly incomplete, at the time of writing even some
essential functions (e.g. addEdges) are missing. Since I don't use Java
intensively, chances are that this interface gets finished only if there
is substantial demand for it and/or someone takes the time to send patches
or finish it completely.

*/

package net.sf.igraph.util;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Iterator generating numbers in increasing order from min to max with a given step size.
 *
 * The interval is inclusive from the left (i.e. min will be included) but exclusive from
 * the right (i.e. max will not be included).
 */
public class LongRangeIterator implements Iterator<Long> {
    /// The next value to be returned
    private long nextValue;

    /// The maximum value
    private final long max;

    /// The step size
    private final long step;

    /// Constructor that creates a range iterator in [min..max) with a given step size
    public LongRangeIterator(long min, long max, long step) {
        if (min > max) {
            throw new IllegalArgumentException("min must be <= max");
        }
        this.nextValue = min;
        this.max = max;
        this.step = step;
    }

    /// Constructor that creates a range iterator in [min..max) with step size 1
    public LongRangeIterator(long min, long max) {
        this(min, max, 1);
    }

    /// Constructor that creates an unlimited range iterator starting from min with step size 1
    public LongRangeIterator(long min) {
        this(min, Long.MAX_VALUE, 1);
    }

    /// Checks whether there are more elements left
    public boolean hasNext() {
        return nextValue < max;
    }

    /// Returns the next element
    public Long next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        Long result = Long.valueOf(nextValue);
        nextValue = nextValue + step;
        return result;
    }

    /// Removes the current element - not implemented of course
    public void remove() {
        throw new UnsupportedOperationException();
    }
}
