
# -*- coding: utf-8 -*-

u'''Handle 2-d NumPy or other arrays as I{LatLon}s, psuedo-x, -ys.

NumPy arrays are assumed to contain rows of points with a lat-, a
longitude -and possibly other- values in different columns.  Iterating
over the array rows, creates an instance of a given I{LatLon} class
"on-the-fly" for each row with the row's lat- and longitude.  The
array is read-accessed only and never duplicated, except to create
a I{subset}.

For example, to process a NumPy array, wrap the array by instantiating
class L{Numpy2LatLon} and specifying the index for the lat- and
longitude column in each row.  Then, pass the L{Numpy2LatLon} instance
to any L{pygeodesy} function or method as the I{points} argument.

Tested with 64-bit Python 2.6.9 (and numpy 1.6.2), 2.7.13 (and numpy 1.13.1),
3.5.3 and 3.6.2 on macOS 10.12.6 Sierra, with 64-bit Intel-Python 3.5.3 (and
numpy 1.11.3) on macOS 10.12.5 Sierra and with Pythonista 3.1 using 64-bit
Python 2.7.12 and 3.5.1 (both with numpy 1.8.0) on iOS 10.3.3.

@newfield example: Example, Examples
'''
from utils import fdot, fsum, isint, issequence, iStr, \
                  polygon, wrap90, wrap180
try:
    from collections import Sequence as _Sequence  # immutable
except ImportError:
    _Sequence = object
from inspect import isclass
from math import radians

__all__ = ('LatLon2psxy', 'Numpy2LatLon',  # class
           'bounds', 'isclockwise', 'isconvex')
__version__ = '17.08.10'


class _Basequence(_Sequence):  # immutable, on purpose
    '''(INTERNAL) Base class.
    '''
    _array = []
    _itemname = 'point'

    def _contains(self, point):
        '''(INTERNAL) Check for a matching point.
        '''
        for _ in self._findall(point, ()):
            return True
        return False

    def _count(self, point):
        '''(INTERNAL) Count the number of matching points.
        '''
        n = 0
        for _ in self._findall(point, ()):
            n += 1
        return n

    def _find(self, point, start_end):
        '''(INTERNAL) Find the first matching point index.
        '''
        for i in self._findall(point, start_end):
            return i
        return -1

    def _findall(self, unused, start_end):  # PYCHOK unused
        '''Must be overloaded.
        '''
        raise NotImplementedError('method: %s' % ('_findall',))

    def _getitem(self, index):
        '''(INTERNAL) Return point [index] or return a slice.
        '''
        # Luciano Ramalho, "Fluent Python", page 290+, O'Reilly, 2016
        if isinstance(index, slice):
            # XXX an numpy.array slice is a view, not a copy
            return self.__class__(self._array[index], **self._slicekwds())
        else:
            return self.point(self._array[index])

    def _index(self, point, start_end):
        '''(INTERNAL) Find the first matching point index.
        '''
        for i in self._findall(point, start_end):
            return i
        raise IndexError('%s not found: %r' % (self._itemname, point))

    def _iter(self):
        '''(INTERNAL) Yield all points.
        '''
        for row in self._array:
            yield self.point(row)

    def point(self, unused):  # PYCHOK unused
        '''Must be overloaded.
        '''
        raise NotImplementedError('method: %s' % ('point',))

    def _range(self, start=None, end=None, step=1):
        '''(INTERNAL) Return the range.
        '''
        if step > 0:
            if start is None:
                start = 0
            if end is None:
                end = len(self)
        elif step < 0:
            if start is None:
                start = len(self) - 1
            if end is None:
                end = -1
        else:
            raise ValueError('%s invalid: %r' % ('step', step))
        return range(start, end, step)

    def _repr(self):
        '''(INTERNAL) Return a string representation.
        '''
        # XXX use Python 3+ reprlib.repr
        t = repr(self._array[:1])  # only first row
        t = '%s, ...%s[%s]' % (t[:-1], t[-1:], len(self))
        t = ' '.join(t.split())  # coalesce spaces
        return iStr(self, t, **self._slicekwds())

    def _reversed(self):  # PYCHOK false
        '''(INTERNAL) Yield all points in reverse order.
        '''
        for row in reversed(self._array):
            yield self.point(row)

    def _rfind(self, point, start_end):
        '''(INTERNAL) Find the last matching point index.
        '''
        def _rt3(start=None, end=None):
            return (start, end, -1)

        for i in self._findall(point, _rt3(*start_end)):
            return i
        return -1

    def _slicekwds(self):
        '''(INTERNAL) Should be overloaded.
        '''
        return {}


class LatLon2psxy(_Basequence):
    '''Wrapper for I{LatLon} points to pseudo-x- and -y-coordinates.
    '''
    _closed   = False
    _len      = 0
    _deg2m    = None  # default, keep degrees
    _radius   = None
    _wrap     = True

    def __init__(self, latlons, closed=False, radius=None, wrap=True):
        '''Handle I{LatLon} points as pseudo-x, -y coordinates.

           @param latlons: Points list, sequence, set, tuple, etc. (LatLon[]).
           @keyword closed: Points are a closed polygon (bool).
           @keyword radius: Optional, mean earth radius (meter).
           @keyword wrap: Optionally, wrap90(lat) and wrap180(lon) (bool).

           @raise TypeError: Some points are not I{LatLon}.

           @raise ValueError: Too few points.
        '''
        if closed:
            self._closed = True
        self._len, self._array = polygon(latlons, closed=closed)
        if radius:
            self._radius = radius
            self._deg2m = radius * radians(1.0)
        if not wrap:
            self._wrap = False

    def __contains__(self, xy):
        '''Check for a matching point.

           @param xy: Point (I{LatLon}) or 2-tuple (x, y).

           @return: True if present, False otherwise.

           @raise TypeError: Invalid xy.
        '''
        return self._contains(xy)

    def __getitem__(self, index):
        '''Return the pseudo-x, -y or return a L{LatLon2psxy} slice.
        '''
        return self._getitem(index)

    def __iter__(self):
        '''Yield all pseudo-x, -y's.
        '''
        return self._iter()

    def __len__(self):
        '''Return the number of pseudo-x and -y's.
        '''
        return self._len

    def __repr__(self):
        '''Return a string representation.
        '''
        return self._repr()

    def __reversed__(self):  # PYCHOK false
        '''Yield all pseudo-x, -y's in reverse order.
        '''
        return self._reversed()

    __str__ = __repr__

    def count(self, xy):
        '''Count the number of matching points.

           @param xy: Point (I{LatLon}) or 2-tuple (x, y).

           @return: Count (integer).

           @raise TypeError: Invalid xy.
        '''
        return self._count(xy)

    def find(self, xy, *start_end):
        '''Find the first matching point.

           @param xy: Point (I{LatLon}) or 2-tuple (x, y).
           @param start_end: Optional [start [, end]] index (integers).

           @return: Index or -1 if not found (integer).

           @raise TypeError: Invalid xy.
        '''
        return self._find(xy, start_end)

    def _findall(self, xy, start_end):
        '''(INTERNAL) Yield indices of all matching points.
        '''
        try:
            x, y = xy.lon, xy.lat

            def _3xyll(ll):  # match LatLon
                return ll.lon, ll.lat, ll

        except AttributeError:
            try:
                x, y = xy
            except (TypeError, ValueError):
                raise TypeError('%s invalid: %r' % ('xy', xy))

            def _3xyll(ll):  # PYCHOK expected
                return self.point(ll)

        for i in self._range(*start_end):
            xi, yi, _ = _3xyll(self._array[i])
            if xi == x and yi == y:
                yield i

    def findall(self, xy, *start_end):
        '''Yield indices of all matching points.

           @param xy: Point (I{LatLon}) or 2-tuple (x, y).
           @param start_end: Optional [start [, end]] index (integers).

           @return: Indices (iterator).

           @raise TypeError: Invalid xy.
        '''
        return self._findall(xy, start_end)

    def index(self, xy, *start_end):  # PYCHOK Python 2- issue
        '''Find the first matching point.

           @param xy: Point (I{LatLon}) or 2-tuple (x, y).
           @param start_end: Optional [start [, end]] index (integers).

           @return: Index (integer).

           @raise IndexError: Point not found.

           @raise TypeError: Invalid xy.
        '''
        return self._index(xy, start_end)

    @property
    def isPoints2(self):
        '''Is this a Points2 wrapper/converter?
        '''
        return True  # isinstance(self, Points))

#   next = __iter__

    def point(self, ll):
        '''Create a pseudo-x and -y.

           @param ll: Point (I{LatLon}).

           @return: 3-Tuple (x, y, ll) of (float, float, I{ll}).
        '''
        x, y = ll.lon, ll.lat  # note, x, y = lon, lat
        if self._wrap:
            x, y = wrap180(x), wrap90(y)
        if self._deg2m:  # convert degrees to meter
            x, y = x * self._deg2m, y * self._deg2m
        return x, y, ll

    def rfind(self, xy, *start_end):
        '''Find the last matching point.

           @param xy: Point (I{LatLon}) or 2-tuple (x, y).
           @param start_end: Optional [start [, end]] index (integers).

           @return: Index or -1 if not found (integer).

           @raise TypeError: Invalid xy.
        '''
        return self._rfind(xy, start_end)

    def _slicekwds(self):
        '''(INTERNAL) Slice kwds.
        '''
        return dict(closed=False, radius=self._radius, wrap=self._wrap)


class _LatLon(object):
    '''(INTERNAL) L{Numpy2LatLon} helper'
    '''
    __slots__ = ('lat', 'lon')

    def __init__(self, lat, lon):
        self.lat = lat
        self.lon = lon

    def __eq__(self, other):
        return isinstance(other, _LatLon) and \
                          other.lat == self.lat and \
                          other.lon == self.lon

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return iStr(self, lat=self.lat, lon=self.lon)

    __str__ = __repr__

    def to2ab(self):
        '''Return lat-/longitude in radians.

           @return: 2-Tuple (lat, lon) in (radians, radians).
        '''
        return radians(self.lat), radians(self.lon)


class Numpy2LatLon(_Basequence):  # immutable, on purpose
    '''Wrapper for NumPy arrays as "on-the-fly" I{LatLon} points.
    '''
    _ilat = 0  # row column index
    _ilon = 0  # row column index
    _LatLon = _LatLon  # default
    _shape = ()

    def __init__(self, array, ilat=0, ilon=1, LatLon=None):
        '''Handle a NumPy array as I{simplify...}-compatible I{LatLon} points.

           @param array: NumPy array (I{numpy.array}).
           @keyword ilat: index of the latitudes column (integer).
           @keyword ilon: index of the longitudes column (integer).
           @keyword LatLon: I{LatLon} class to use (internal).

           @raise IndexError: If I{array.shape} is not (1+, 2+).

           @raise TypeError: If I{array} is not a NumPy array or
                             I{LatLon} is not a class with I{ilat}
                             and I{ilon} attributes.

           @raise ValueError: If the I{ilat} and/or I{ilon} values are
                              out of range or the same.

           @example:

           >>> type(array)
           <type 'numpy.ndarray'>  # <class ...> in Python 3+
           >>> points = Numpy2LatLon(array, lat=0, lon=1)
           >>> simply = simplifyRDP(points, ...)
           >>> type(simply)
           <type 'numpy.ndarray'>  # <class ...> in Python 3+
           >>> sliced = points[1:-1]
           >>> type(sliced)
           <class '...Numpy2LatLon'>
        '''
        ais = [('ilat', ilat), ('ilon', ilon)]

        try:  # get shape and check some other numpy.array attrs
            s, _, _ = array.shape, array.nbytes, array.ndim  # PYCHOK expected
        except AttributeError:
            raise TypeError('%s not NumPy: %s' % ('array', type(array)))
        if len(s) != 2 or s[0] < 1 or s[1] < len(ais):
            raise IndexError('%s shape invalid: %r' % ('array', s))
        self._array = array
        self._shape = s

        # check the point class
        if LatLon is not None:
            if isclass(LatLon) and all(hasattr(LatLon, a) for a in _LatLon.__slots__):
                self._LatLon = LatLon
            else:
                raise TypeError('%s invalid: %r' % ('LatLon', LatLon))

        # check the attr indices
        for n, (ai, i) in enumerate(ais):
            if not isint(i):
                raise TypeError('%s invalid: %r' % (ai, i))
            i = int(i)
            if not 0 <= i < s[1]:
                raise ValueError('%s invalid: %s' % (ai, i))
            for aj, j in ais[:n]:
                if int(j) == i:
                    raise ValueError('%s == %s == %s' % (ai, aj, i))
            setattr(self, '_' + ai, i)

    def __contains__(self, latlon):
        '''Check for a specific lat-/longitude.

           @param latlon: Point (I{LatLon}) or 2-tuple (lat, lon).

           @return: True if present, False otherwise.

           @raise TypeError: Invalid latlon.
        '''
        return self._contains(latlon)

    def __getitem__(self, index):
        '''Return row[index] as I{LatLon} or return a L{Numpy2LatLon} slice.
        '''
        return self._getitem(index)

    def __iter__(self):
        '''Yield rows as I{LatLon}.
        '''
        return self._iter()

    def __len__(self):
        '''Return the number of rows.
        '''
        return self._shape[0]

    def __repr__(self):
        '''Return a string representation.
        '''
        return self._repr()

    def __reversed__(self):  # PYCHOK false
        '''Yield rows as I{LatLon} in reverse order.
        '''
        return self._reversed()

    __str__ = __repr__

    def count(self, latlon):
        '''Count the number of rows with a specific lat-/longitude.

           @param latlon: Point (I{LatLon}) or 2-tuple (lat, lon)

           @return: Count (integer).

           @raise TypeError: Invalid latlon.
        '''
        return self._count(latlon)

    def find(self, latlon, *start_end):
        '''Find the first row with a specific lat-/longitude.

           @param latlon: Point (I{LatLon}) or 2-tuple (lat, lon).
           @param start_end: Optional [start [, end]] index (integers).

           @return: Index or -1 if not found (integer).

           @raise TypeError: Invalid latlon.
        '''
        return self._find(latlon, start_end)

    def _findall(self, latlon, start_end):
        '''(INTERNAL) Yield indices of all matching rows.
        '''
        try:
            lat, lon = latlon.lat, latlon.lon
        except AttributeError:
            try:
                lat, lon = latlon
            except (TypeError, ValueError):
                raise TypeError('%s invalid: %r' % ('latlon', latlon))

        for i in self._range(*start_end):
            row = self._array[i]
            if row[self._ilat] == lat and \
               row[self._ilon] == lon:
                yield i

    def findall(self, latlon, *start_end):
        '''Yield indices of all rows with a specific lat-/longitude.

           @param latlon: Point (I{LatLon}) or 2-tuple (lat, lon).
           @param start_end: Optional [start [, end]] index (integers).

           @return: Indices (iterator).

           @raise TypeError: Invalid latlon.
        '''
        return self._findall(latlon, start_end)

    def index(self, latlon, *start_end):  # PYCHOK Python 2- issue
        '''Find index of the first row with a specific lat-/longitude.

           @param latlon: Point (I{LatLon}) or 2-tuple (lat, lon).
           @param start_end: Optional [start [, end]] index (integers).

           @return: Index (integer).

           @raise IndexError: Point not found.

           @raise TypeError: Invalid latlon.
        '''
        return self._index(latlon, start_end)

    @property
    def ilat(self):
        '''Get the latitudes column index (integer).
        '''
        return self._ilat

    @property
    def ilon(self):
        '''Get the longitudes column index (integer).
        '''
        return self._ilon

    @property
    def isNumpy2(self):
        '''Is this a Numpy2 wrapper?
        '''
        return True  # isinstance(self, (Numpy2LatLon, ...))

#   next = __iter__

    def point(self, row):
        '''Instantiate a point I{LatLon}.

           @param row: Array row (numpy.array).

           @return: Point (I{LatLon}).
        '''
        return self._LatLon(row[self._ilat], row[self._ilon])

    def rfind(self, latlon, *start_end):
        '''Find the last row with a specific lat-/longitude.

           @param latlon: Point (I{LatLon}) or 2-tuple (lat, lon).
           @param start_end: Optional [start [, end]] index (integers).

           @note: Keyword order, first stop, then start.

           @return: Index or -1 if not found (integer).

           @raise TypeError: Invalid latlon.
        '''
        return self._rfind(latlon, start_end)

    def _slicekwds(self):
        '''(INTERNAL) Slice kwds.
        '''
        return dict(ilat=self._ilat, ilon=self._ilon)

    @property
    def shape(self):
        '''Get the shape of the NumPy array (2-tuple).
        '''
        return self._shape

    def subset(self, indices):
        '''Return a subset of the NumPy array.

           @param indices: Row indices (ints).

           @note: A I{subset} is different from a I{slice} in 2 ways:
                  (a) the I{subset} is typically specified as a list of
                  (un-)ordered indices and (b) the I{subset} allocates
                  a new, separate NumPy array while a I{slice} is just
                  an other I{view} of the original NumPy array.

           @return: Sub-array (numpy.array).

           @raise TypeError: If indices is not a I{range}
                             or I{list} of ints.

           @raise ValueError: Out of range indices value.
        '''
        if not issequence(indices, tuple):  # NO tuple, only list
            # and range work properly to get Numpy array sub-sets
            raise TypeError('%s invalid: %s' % ('indices', type(indices)))

        n = len(self)
        for i, v in enumerate(indices):
            if not isint(v):
                raise TypeError('%s[%s] invalid: %r' % ('indices', i, v))
            elif not 0 <= v < n:
                raise ValueError('%s[%s] invalid: %r' % ('indices', i, v))

        return self._array[indices]


def bounds(points, radius=None, wrap=True):
    '''Determine the lower-left and upper-right corners of a polygon
       defined by a list, sequence, set or tuple of I{LatLon} points.

       @param points: The points defining the polygon (I{LatLon}[]).
       @keyword radius: Optional, mean earth radius (meter).
       @keyword wrap: Optionally, wrap90(lat) and wrap180(lon) (bool).

       @return: 4-Tuple (lolat, lolon, hilat, hilon) corners (degrees).

       @raise TypeError: Some points are not I{LatLon}.

       @raise ValueError: Too few points.

       @example:

       >>> b = LatLon(45,1), LatLon(45,2), LatLon(46,2), LatLon(46,1)
       >>> bounds(b)  # False
       >>> 45.0, 1.0, 46.0, 2.0
    '''
    pts = LatLon2psxy(points, closed=False, radius=radius, wrap=wrap)

    lox, loy, _ = hix, hiy, _ = pts[0]

    for x, y, _ in pts:  # [1:]
        if lox > x:
            lox = x
        elif hix < x:
            hix = x
        if loy > y:
            loy = y
        elif hiy < y:
            hiy = y

    return loy, lox, hiy, hix


def isclockwise(points, radius=None, wrap=True):
    '''Determine the direction of a polygon defined by a list,
       sequence, set or tuple of I{LatLon} points.

       @param points: The points defining the polygon (I{LatLon}[]).
       @keyword radius: Optional, mean earth radius (meter).
       @keyword wrap: Optionally, wrap90(lat) and wrap180(lon) (bool).

       @return: True if clockwise, False otherwise.

       @raise TypeError: Some points are not I{LatLon}.

       @raise ValueError: Too few points or zero area polygon.

       @example:

       >>> f = LatLon(45,1), LatLon(45,2), LatLon(46,2), LatLon(46,1)
       >>> isclockwise(f)  # False

       >>> t = LatLon(45,1), LatLon(46,1), LatLon(46,2)
       >>> isclockwise(t)  # True
    '''
    pts = LatLon2psxy(points, closed=True, radius=radius, wrap=wrap)

    if len(pts):

        def _sareas(pts):  # iterate
            x1, y1, _ = pts[-1]
            for x2, y2, _ in pts:
                yield (x2 - x1) * (y2 + y1)  # segment pseudo-area
                x1, y1 = x2, y2

        pa = fsum(_sareas(pts))  # signed pseudo-area
        if pa > 0:
            return True
        elif pa < 0:
            return False

    raise ValueError('zero area: %r' % (points[:3],))


def isconvex(points, radius=None, wrap=True):
    '''Determine whether a polygon defined by a list, sequence,
       set or tuple of I{LatLon} points is convex.

       @param points: The points defining the polygon (I{LatLon}[]).
       @keyword radius: Optional, mean earth radius (meter).
       @keyword wrap: Optionally, wrap90(lat) and wrap180(lon) (bool).

       @return: True if convex, False otherwise.

       @raise TypeError: Some points are not I{LatLon}.

       @raise ValueError: Too few points or a colinear point.

       @example:

       >>> t = LatLon(45,1), LatLon(46,1), LatLon(46,2)
       >>> isconvex(t)  # True

       >>> f = LatLon(45,1), LatLon(46,2), LatLon(45,2), LatLon(46,1)
       >>> isconvex(f)  # False
    '''
    pts = LatLon2psxy(points, closed=True, radius=radius, wrap=wrap)

    s = None
    x1, y1, _ = pts[-2]
    x2, y2, _ = pts[-1]
    for x3, y3, ll in pts:
        # get the sign of the distance from point
        # x3, y3 to the line from x1, y1 to x2, y2
        # <http://wikipedia.org/wiki/Distance_from_a_point_to_a_line>
        s3 = fdot((x3, y3, x1, y1), y2 - y1, x1 - x2, -y2, x2)
        if s3 > 0:  # x3, y3 on the left
            if s is None:
                s = True
            elif not s:  # different side
                return False
        elif s3 < 0:  # x3, y3 on the right
            if s is None:
                s = False
            elif s:  # different side
                return False
        elif fdot((x3 - x2, y1 - y2), y3 - y2, x1 - x2) < 0:
            # colinear u-turn: x3, y3 not on the
            # opposite side of x2, y2 as x1, y1
            raise ValueError('colinear: %r' % (ll,))
        x1, y1, x2, y2 = x2, y2, x3, y3

    return True  # all points on the same side

# **) MIT License
#
# Copyright (C) 2016-2017 -- mrJean1 at Gmail dot com
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
