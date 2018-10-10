
# -*- coding: utf-8 -*-

u'''Handle 2-d U{NumPy<http://www.NumPy.org>}, I{arrays} or tuples
as I{LatLon}s or as I{pseudo-x/-y} pairs.

C{NumPy} arrays are assumed to contain rows of points with a lat-, a
longitude -and possibly other- values in different columns.  While
iterating over the array rows, create an instance of a given I{LatLon}
class "on-the-fly" for each row with the row's lat- and longitude.

The original C{NumPy} array is read-accessed only and never duplicated,
except to create a I{subset} of the original array.

For example, to process a C{NumPy} array, wrap the array by instantiating
class L{Numpy2LatLon} and specifying the column index for the lat- and
longitude in each row.  Then, pass the L{Numpy2LatLon} instance to any
L{pygeodesy} function or method accepting a I{points} argument.

Similarly, class L{Tuple2LatLon} is used to instantiate a I{LatLon}
for each 2+tuple in a list, tuple or sequence of such 2+tuples from
the index for the lat- and longitude index in each 2+tuple.

@newfield example: Example, Examples
'''

from bases import classname, inStr
from fmath import EPS, favg, fdot, fStr, Fsum, fsum, isint, scalar
from utily import R_M, equirectangular_, issequence, points2, \
                  property_RO, unroll180, unrollPI, wrap90, wrap180
from vector3d import CrossError, crosserrors

try:
    from collections import Sequence as _Sequence  # immutable
except ImportError:
    _Sequence = object  # XXX or tuple
from inspect import isclass
from math import cos, hypot, radians, sin

__all__ = ('LatLon_',  # classes
           'LatLon2psxy', 'Numpy2LatLon', 'Tuple2LatLon',
           'areaOf',  # functions
           'bounds',
           'isclockwise', 'isconvex',
           'isenclosedBy', 'isenclosedby',  # DEPRECATED
           'nearestOn3',
           'perimeterOf')
__version__ = '18.10.10'


class LatLon_(object):
    '''Low-overhead I{LatLon} class for L{Numpy2LatLon} or L{Tuple2LatLon}'
    '''
    # __slots__ efficiency is voided if the __slots__ class attribute
    # is used in a subclass of a class with the traditional __dict__,
    # see <http://docs.Python.org/2/reference/datamodel.html#slots>
    # and __slots__ must be repeated in sub-classes, see "Problems
    # with __slots__" in Luciano Ramalho, "Fluent Python", page
    # 276+, O'Reilly, 2016, also at <http://Books.Google.ie/
    #   books?id=bIZHCgAAQBAJ&lpg=PP1&dq=fluent%20python&pg=
    #   PT364#v=onepage&q=“Problems%20with%20__slots__”&f=false>
    __slots__ = ('lat', 'lon', 'name')

    def __init__(self, lat, lon, name=''):
        '''Creat a new, mininal, low-overhead L{LatLon_} instance,
           without heigth and datum.

           @param lat: Latitude (C{degrees}).
           @param lon: Longitude (C{degrees}).
           @keyword name: Optional name (C{str}).

           @note: The lat- and longitude are taken as-given,
                  un-clipped and un-validated.
        '''
        self.lat = float(lat)
        self.lon = float(lon)
        self.name = str(name)

    def __eq__(self, other):
        return isinstance(other, LatLon_) and \
                          other.lat == self.lat and \
                          other.lon == self.lon

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return self.toStr2()

    def __str__(self):
        return self.toStr()

    def others(self, other, name='other'):
        '''Check this and an other instance for type compatiblility.

           @param other: The other instance (any C{type}).
           @keyword name: Optional, name for other (C{str}).

           @return: C{None}.

           @raise TypeError: Incompatible I{other} C{type}.
        '''
        if not (isinstance(other, self.__class__) or
                (hasattr(other, 'lat') and hasattr(other, 'lon'))):
            raise TypeError('type(%s) mismatch: %s vs %s' % (name,
                             classname(other), classname(self)))

    def points(self, points, closed=False, base=None):
        return points2(points, closed=closed, base=base)

    points.__doc__ = points2.__doc__

    def to2ab(self):
        '''Return the lat- and longitude in radians.

           @return: 2-Tuple (lat, lon) in (C{radians}, C{radians}).
        '''
        return radians(self.lat), radians(self.lon)

    def toStr(self, **kwds):
        '''This L{LatLon_} as a string "<degrees>, <degrees>".

           @keyword kwds: Optional, keyword arguments.

           @return: Instance (C{str}).
        '''
        t = [fStr(getattr(self, _)) for _ in self.__slots__]
        if kwds:
            t += ['%s=%s' % _ for _ in sorted(kwds.items())]
        return ', '.join(t)

    def toStr2(self, **kwds):
        '''This L{LatLon_} as a string "class(<degrees>, ...)".

           @keyword kwds: Optional, keyword arguments.

           @return: Class instance (C{str}).
        '''
        return '%s(%s)' % (classname(self), self.toStr(**kwds))


class _Basequence(_Sequence):  # immutable, on purpose
    '''(INTERNAL) Base class.
    '''
    _array    = []
    _epsilon  = EPS
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
        raise ValueError('%s not found: %r' % (self._itemname, point))

    def _iter(self):
        '''(INTERNAL) Yield all points.
        '''
        for i in range(len(self)):
            yield self.point(self._array[i])

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
        return inStr(self, t, **self._slicekwds())

    def _reversed(self):  # PYCHOK false
        '''(INTERNAL) Yield all points in reverse order.
        '''
        for i in range(len(self) - 1, -1):
            yield self.point(self._array[i])

    def _rfind(self, point, start_end):
        '''(INTERNAL) Find the last matching point index.
        '''
        def _r3(start=None, end=None, step=-1):
            return (start, end, step)

        for i in self._findall(point, _r3(*start_end)):
            return i
        return -1

    def _slicekwds(self):
        '''(INTERNAL) Should be overloaded.
        '''
        return {}

    def _zeros(self, *zeros):
        '''(INTERNAL) Check for near-zero values.
        '''
        return all(abs(z) <= self._epsilon for z in zeros)

    @property
    def epsilon(self):
        '''Get the tolerance for equality tests (C{float}).
        '''
        return self._epsilon

    @epsilon.setter  # PYCHOK setter!
    def epsilon(self, tol):
        '''Set the tolerance for equality tests.

           @param tol: New tolerance (C{scalar}).

           @raise TypeError: Non-scalar I{tol}.

           @raise ValueError: Out-of-bounds I{tol}.
        '''
        self._epsilon = scalar(tol, 0.0, name='tolerance')

    @property_RO
    def isNumpy2(self):
        '''Is this a Numpy2 wrapper?
        '''
        return False  # isinstance(self, (Numpy2LatLon, ...))

    @property_RO
    def isPoints2(self):
        '''Is this a LatLon2 wrapper/converter?
        '''
        return False  # isinstance(self, (LatLon2psxy, ...))

    @property_RO
    def isTuple2(self):
        '''Is this a Tuple2 wrapper?
        '''
        return False  # isinstance(self, (Tuple2LatLon, ...))


class _Array2LatLon(_Basequence):  # immutable, on purpose
    '''Base class for Numpy2LatLon or Tuple2LatLon.
    '''
    _array  = ()
    _ilat   = 0  # row column index
    _ilon   = 0  # row column index
    _LatLon = LatLon_  # default
    _shape  = ()

    def __init__(self, array, ilat=0, ilon=1, LatLon=None, shape=()):
        '''Handle a C{NumPy} or C{Tuple} array as a sequence of I{LatLon} points.
        '''
        ais = ('ilat', ilat), ('ilon', ilon)

        if len(shape) != 2 or shape[0] < 1 or shape[1] < len(ais):
            raise IndexError('%s shape invalid: %r' % ('array', shape))
        self._array = array
        self._shape = shape

        # check the point class
        if LatLon is not None:
            if isclass(LatLon) and all(hasattr(LatLon, a) for a in LatLon_.__slots__):
                self._LatLon = LatLon
            else:
                raise TypeError('%s invalid: %r' % ('LatLon', LatLon))

        # check the attr indices
        for n, (ai, i) in enumerate(ais):
            if not isint(i):
                raise TypeError('%s invalid: %r' % (ai, i))
            i = int(i)
            if not 0 <= i < shape[1]:
                raise ValueError('%s invalid: %s' % (ai, i))
            for aj, j in ais[:n]:
                if int(j) == i:
                    raise ValueError('%s == %s == %s' % (ai, aj, i))
            setattr(self, '_' + ai, i)

    def __contains__(self, latlon):
        '''Check for a specific lat-/longitude.

           @param latlon: Point (I{LatLon}) or 2-tuple (lat, lon).

           @return: C{True} if I{latlon} is present, C{False} otherwise.

           @raise TypeError: Invalid I{latlon}.
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

           @return: Count (C{int}).

           @raise TypeError: Invalid I{latlon}.
        '''
        return self._count(latlon)

    def find(self, latlon, *start_end):
        '''Find the first row with a specific lat-/longitude.

           @param latlon: Point (I{LatLon}) or 2-tuple (lat, lon).
           @param start_end: Optional I{[start[, end]]} index (integers).

           @return: Index or -1 if not found (C{int}).

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
            if self._zeros(row[self._ilat] - lat,
                           row[self._ilon] - lon):
                yield i

    def findall(self, latlon, *start_end):
        '''Yield indices of all rows with a specific lat-/longitude.

           @param latlon: Point (I{LatLon}) or 2-tuple (lat, lon).
           @param start_end: Optional I{[start[, end]]} index (C{int}s).

           @return: Indices (C{iterable}).

           @raise TypeError: Invalid I{latlon}.
        '''
        return self._findall(latlon, start_end)

    def index(self, latlon, *start_end):  # PYCHOK Python 2- issue
        '''Find index of the first row with a specific lat-/longitude.

           @param latlon: Point (I{LatLon}) or 2-tuple (lat, lon).
           @param start_end: Optional I{[start[, end]]} index (C{int}s).

           @return: Index (C{int}).

           @raise TypeError: Invalid I{latlon}.

           @raise ValueError: Point not found.
        '''
        return self._index(latlon, start_end)

    @property_RO
    def ilat(self):
        '''Get the latitudes column index (C{int}).
        '''
        return self._ilat

    @property_RO
    def ilon(self):
        '''Get the longitudes column index (C{int}).
        '''
        return self._ilon

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
           @param start_end: Optional I{[start[, end]]} index (C{int}s).

           @note: Keyword order, first stop, then start.

           @return: Index or -1 if not found (C{int}).

           @raise TypeError: Invalid I{latlon}.
        '''
        return self._rfind(latlon, start_end)

    def _slicekwds(self):
        '''(INTERNAL) Slice kwds.
        '''
        return dict(ilat=self._ilat, ilon=self._ilon)

    @property_RO
    def shape(self):
        '''Get the shape of the C{NumPy} array or the C{Tuples} as
           2-tuple of I{(number of rows, number of colums)}.
        '''
        return self._shape

    def _subset(self, indices):  # PYCHOK unused
        '''Must be overloaded.
        '''
        raise NotImplementedError('method: %s' % ('_subset',))

    def subset(self, indices):
        '''Return a subset of the C{NumPy} array.

           @param indices: Row indices (C{int}s).

           @note: A I{subset} is different from a I{slice} in 2 ways:
                  (a) the I{subset} is typically specified as a list of
                  (un-)ordered indices and (b) the I{subset} allocates
                  a new, separate C{NumPy} array while a I{slice} is
                  just an other I{view} of the original C{NumPy} array.

           @return: Sub-array (C{numpy.array}).

           @raise IndexError: Out-of-range I{indices} value.

           @raise TypeError: If I{indices} is not a C{range} or a
                             C{list} of C{int}s.
        '''
        if not issequence(indices, tuple):  # NO tuple, only list
            # and range work properly to get Numpy array sub-sets
            raise TypeError('%s invalid: %s' % ('indices', type(indices)))

        n = len(self)
        for i, v in enumerate(indices):
            if not isint(v):
                raise TypeError('%s[%s] invalid: %r' % ('indices', i, v))
            elif not 0 <= v < n:
                raise IndexError('%s[%s] invalid: %r' % ('indices', i, v))

        return self._subset(indices)


class LatLon2psxy(_Basequence):
    '''Wrapper for I{LatLon} points as "on-the-fly" pseudo-xy coordinates.
    '''
    _closed = False
    _len    = 0
    _deg2m  = None  # default, keep degrees
    _radius = None
    _wrap   = True

    def __init__(self, latlons, closed=False, radius=None, wrap=True):
        '''Handle I{LatLon} points as pseudo-xy coordinates.

           @note: The I{LatLon}'s latitude is considered the pseudo-y
                  and longitude the pseudo-x coordinate.  Similarly,
                  2-tuples (x, y) are (longitude, latitude).

           @param latlons: Points C{list}, C{sequence}, C{set}, C{tuple},
                           etc. (I{LatLon[]}).
           @keyword closed: Optionally, close the polygon (C{bool}).
           @keyword radius: Optional, mean earth radius (C{meter}).
           @keyword wrap: Wrap lat- and longitudes (C{bool}).

           @raise TypeError: Some I{latlons} are not I{LatLon}.

           @raise ValueError: Insufficient number of I{latlons}.
        '''
        self._closed = closed
        self._len, self._array = points2(latlons, closed=closed)
        if radius:
            self._radius = radius
            self._deg2m = radius * radians(1.0)
        self._wrap = wrap

    def __contains__(self, xy):
        '''Check for a matching point.

           @param xy: Point (I{LatLon}) or 2-tuple (x, y) in (C{degrees}).

           @return: C{True} if I{xy} is present, C{False} otherwise.

           @raise TypeError: Invalid I{xy}.
        '''
        return self._contains(xy)

    def __getitem__(self, index):
        '''Return the pseudo-xy or return a L{LatLon2psxy} slice.
        '''
        return self._getitem(index)

    def __iter__(self):
        '''Yield all pseudo-xy's.
        '''
        return self._iter()

    def __len__(self):
        '''Return the number of pseudo-xy's.
        '''
        return self._len

    def __repr__(self):
        '''Return a string representation.
        '''
        return self._repr()

    def __reversed__(self):  # PYCHOK false
        '''Yield all pseudo-xy's in reverse order.
        '''
        return self._reversed()

    __str__ = __repr__

    def count(self, xy):
        '''Count the number of matching points.

           @param xy: Point (I{LatLon}) or 2-tuple (x, y) in (C{degrees}).

           @return: Count (C{int}).

           @raise TypeError: Invalid I{xy}.
        '''
        return self._count(xy)

    def find(self, xy, *start_end):
        '''Find the first matching point.

           @param xy: Point (I{LatLon}) or 2-tuple (x, y) in (C{degrees}).
           @param start_end: Optional I{[start[, end]]} index (C{int}s).

           @return: Index or -1 if not found (C{int}).

           @raise TypeError: Invalid I{xy}.
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
            if self._zeros(xi - x, yi - y):
                yield i

    def findall(self, xy, *start_end):
        '''Yield indices of all matching points.

           @param xy: Point (I{LatLon}) or 2-tuple (x, y) in (C{degrees}).
           @param start_end: Optional I{[start[, end]]} index (C{int}s).

           @return: Indices (C{iterator}).

           @raise TypeError: Invalid I{xy}.
        '''
        return self._findall(xy, start_end)

    def index(self, xy, *start_end):  # PYCHOK Python 2- issue
        '''Find the first matching point.

           @param xy: Point (I{LatLon}) or 2-tuple (x, y) in (C{degrees}).
           @param start_end: Optional I{[start[, end]]} index (C{int}s).

           @return: Index (C{int}).

           @raise TypeError: Invalid C{xy}.

           @raise ValueError: Point not found.
        '''
        return self._index(xy, start_end)

    @property_RO
    def isPoints2(self):
        '''Is this a LatLon2 wrapper/converter?
        '''
        return True  # isinstance(self, (LatLon2psxy, ...))

#   next = __iter__

    def point(self, ll):
        '''Create a pseudo-xy.

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

           @param xy: Point (I{LatLon}) or 2-tuple (x, y) in (C{degrees}).
           @param start_end: Optional I{[start[, end]]} index (C{int}s).

           @return: Index or -1 if not found (C{int}).

           @raise TypeError: Invalid I{xy}.
        '''
        return self._rfind(xy, start_end)

    def _slicekwds(self):
        '''(INTERNAL) Slice kwds.
        '''
        return dict(closed=self._closed, radius=self._radius, wrap=self._wrap)


class Numpy2LatLon(_Array2LatLon):  # immutable, on purpose
    '''Wrapper for C{NumPy} arrays as "on-the-fly" I{LatLon} points.
    '''
    def __init__(self, array, ilat=0, ilon=1, LatLon=None):
        '''Handle a C{NumPy} array as a sequence of I{LatLon} points.

           @param array: C{NumPy} array (I{numpy.array}).
           @keyword ilat: Optional index of the latitudes column (C{int}).
           @keyword ilon: Optional index of the longitudes column (C{int}).
           @keyword LatLon: Optional I{LatLon} (sub-)class to use (L{LatLon_}).

           @raise IndexError: If I{array.shape} is not (1+, 2+).

           @raise TypeError: If I{array} is not a C{NumPy} array or
                             I{LatLon} is not a class with I{lat}
                             and I{lon} attributes.

           @raise ValueError: If the I{ilat} and/or I{ilon} values are
                              the same or out of range.

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
        try:  # get shape and check some other numpy.array attrs
            s, _, _ = array.shape, array.nbytes, array.ndim  # PYCHOK expected
        except AttributeError:
            raise TypeError('%s not NumPy: %s' % ('array', type(array)))

        _Array2LatLon.__init__(self, array, ilat=ilat, ilon=ilon,
                                     LatLon=LatLon, shape=s)

    @property_RO
    def isNumpy2(self):
        '''Is this a Numpy2 wrapper?
        '''
        return True  # isinstance(self, (Numpy2LatLon, ...))

    def _subset(self, indices):
        return self._array[indices]  # NumPy special


class Tuple2LatLon(_Array2LatLon):
    '''Wrapper for tuple sequences as "on-the-fly" I{LatLon} points.
    '''
    def __init__(self, tuples, ilat=0, ilon=1, LatLon=None):
        '''Handle a list of tuples, each containing a lat- and longitude
           and perhaps other values as a sequence of I{LatLon} points.

           @param tuples: The C{list}, C{tuple} or C{sequence} of tuples (C{tuple}[]).
           @keyword ilat: Optional index of the latitudes value (C{int}).
           @keyword ilon: Optional index of the longitudes value (C{int}).
           @keyword LatLon: Optional I{LatLon} (sub-)class to use (L{LatLon_}).

           @raise IndexError: If I{(len(tuples), min(len(t) for t in tuples))}
                              is not (1+, 2+).

           @raise TypeError: If I{tuples} is not a C{list}, C{tuple} or
                             C{sequence} or if I{LatLon} is not a C{class}
                             with I{lat} and I{lon} attributes.

           @raise ValueError: If the I{ilat} and/or I{ilon} values are
                              the same or out of range.

           @example:

           >>> tuples = [(0, 1), (2, 3), (4, 5)]
           >>> type(tuples)
           <type 'list'>  # <class ...> in Python 3+
           >>> points = Tuple2LatLon(tuples, lat=0, lon=1)
           >>> simply = simplifyRW(points, 0.5, ...)
           >>> type(simply)
           <type 'list'>  # <class ...> in Python 3+
           >>> simply
           [(0, 1), (4, 5)]
           >>> sliced = points[1:-1]
           >>> type(sliced)
           <class '...Tuple2LatLon'>
           >>> sliced
           ...Tuple2LatLon([(2, 3), ...][1], ilat=0, ilon=1)

           >>> closest, _ = nearestOn2(LatLon_(2, 1), points, adjust=False)
           >>> closest
           LatLon_(lat=1.0, lon=2.0)

           >>> closest, _ = nearestOn2(LatLon_(3, 2), points)
           >>> closest
           LatLon_(lat=2.001162, lon=3.001162)
        '''
        if isinstance(tuples, (list, tuple)):
            s = len(tuples), min(len(_) for _ in tuples)
        else:
            TypeError('%s not sequence: %s' % ('tuple', type(tuples)))
        _Array2LatLon.__init__(self, tuples, ilat=ilat, ilon=ilon,
                                     LatLon=LatLon, shape=s)

    @property_RO
    def isTuple2(self):
        '''Is this a Tuple2 wrapper?
        '''
        return True  # isinstance(self, (Tuple2LatLon, ...))

    def _subset(self, indices):
        return type(self._array)(self._array[i] for i in indices)


def _areaOf(points, adjust, wrap):
    # return the signed area in radians squared

    # setting radius=1 converts degrees to radians
    pts = LatLon2psxy(points, closed=True, radius=1, wrap=wrap)

    def _rads2(n, pts):  # trapezoidal areas in rads**2
        x1, y1, _ = pts[n-1]
        for i in range(n):
            x2, y2, _ = pts[i]
            w, x2 = unrollPI(x1, x2, wrap=wrap)
            # approximate trapezoid by a rectangle, adjusting
            # the top width by the cosine of the latitudinal
            # average and bottom width by some fudge factor
            h = (y2 + y1) * 0.5
            if adjust:
                w *= (cos(h) + 1.2876) * 0.5
            yield h * w  # signed trapezoidal area

            x1, y1 = x2, y2

    return fsum(_rads2(len(pts), pts))


def areaOf(points, adjust=True, radius=R_M, wrap=True):
    '''Approximate the area of a polygon.

       @param points: The polygon points (I{LatLon}s).
       @keyword adjust: Adjust the wrapped, unrolled longitudinal delta
                        by the cosine of the mean latitude (C{bool}).
       @keyword radius: Optional, mean earth radius (C{meter}).
       @keyword wrap: Wrap lat-, wrap and unroll longitudes (C{bool}).

       @return: Approximate area (C{meter}, same units as I{radius}, squared).

       @raise TypeError: Some I{points} are not I{LatLon}.

       @raise ValueError: Insufficient number of I{points}.

       @note: This is an area approximation with limited accuracy,
       ill-suited for regions exceeding several hundred Km or Miles
       or with near-polar latitudes.

       @see: L{sphericalNvector.areaOf}, L{sphericalTrigonometry.areaOf}
             and L{ellipsoidalKarney.areaOf}.
    '''
    return abs(_areaOf(points, adjust, wrap)) * float(radius)**2


def bounds(points, wrap=True, LatLon=None):
    '''Determine the lower-left and upper-right corners of a polygon.

       @param points: The polygon points (I{LatLon}s).
       @keyword wrap: Wrap lat- and longitudes (C{bool}).
       @keyword LatLon: Optional (sub-)class to use to return I{bounds}
                        (I{LatLon}) or C{None}.

       @return: 2-tuple (loLatLon, hiLatLon) of I{LatLon}s for the
                lower-left respectively upper-right corners or 4-Tuple
                (loLat, loLon, hiLat, hiLon) of bounds (C{degrees}) if
                I{LatLon} is C{None}.

       @raise TypeError: Some I{points} are not I{LatLon}.

       @raise ValueError: Insufficient number of I{points}.

       @example:

       >>> b = LatLon(45,1), LatLon(45,2), LatLon(46,2), LatLon(46,1)
       >>> bounds(b)  # False
       >>> 45.0, 1.0, 46.0, 2.0
    '''
    pts = LatLon2psxy(points, closed=False, radius=None, wrap=wrap)

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

    if LatLon is None:
        b = loy, lox, hiy, hix
    else:
        b = LatLon(loy, lox), LatLon(hiy, hix)
    return b


def _imdex2(closed, n):  # imported by sphericalNvector, -Trigonometry
    '''(INTERNAL) Return first and second index.
    '''
    return (n-1, 0) if closed else (0, 1)


def isclockwise(points, adjust=False, wrap=True):
    '''Determine the direction of a polygon.

       @param points: The polygon points (I{LatLon}s).
       @keyword adjust: Adjust the wrapped, unrolled longitudinal delta
                        by the cosine of the mean latitude (C{bool}).
       @keyword wrap: Wrap lat-, wrap and unroll longitudes (C{bool}).

       @return: C{True} if I{points} are clockwise, C{False} otherwise.

       @raise TypeError: Some I{points} are not I{LatLon}.

       @raise ValueError: Insufficient number of I{points} or zero
                          area polygon.

       @example:

       >>> f = LatLon(45,1), LatLon(45,2), LatLon(46,2), LatLon(46,1)
       >>> isclockwise(f)  # False

       >>> t = LatLon(45,1), LatLon(46,1), LatLon(46,2)
       >>> isclockwise(t)  # True
    '''
    a = _areaOf(points, adjust, wrap) or 0
    if a > 0:
        return True
    elif a < 0:
        return False

    raise ValueError('zero area: %r' % (points[:3],))


def isconvex(points, adjust=False, wrap=True):
    '''Determine whether a polygon is convex.

       @param points: The polygon points (I{LatLon}s).
       @keyword adjust: Adjust the wrapped, unrolled longitudinal delta
                        by the cosine of the mean latitude (C{bool}).
       @keyword wrap: Wrap lat-, wrap and unroll longitudes (C{bool}).

       @return: C{True} if I{points} are convex, C{False} otherwise.

       @raise CrossError: Some I{points} are colinear.

       @raise TypeError: Some I{points} are not I{LatLon}.

       @raise ValueError: Insufficient number of I{points}.

       @example:

       >>> t = LatLon(45,1), LatLon(46,1), LatLon(46,2)
       >>> isconvex(t)  # True

       >>> f = LatLon(45,1), LatLon(46,2), LatLon(45,2), LatLon(46,1)
       >>> isconvex(f)  # False
    '''
    def _unroll_adjust(x1, y1, x2, y2):
        x21, x2 = unroll180(x1, x2, wrap=wrap)
        if adjust:
            x21 *= cos(radians(y1 + y2) * 0.5)
        return x21, x2

    pts = LatLon2psxy(points, closed=True, radius=None, wrap=wrap)
    c, n, s = crosserrors(), len(pts), None

    x1, y1, _ = pts[n-2]
    x2, y2, _ = pts[n-1]
    x21, x2 = _unroll_adjust(x1, y1, x2, y2)

    for i in range(n):
        x3, y3, ll = pts[i]
        x32, x3 = _unroll_adjust(x2, y2, x3, y3)

        # get the sign of the distance from point
        # x3, y3 to the line from x1, y1 to x2, y2
        # <http://WikiPedia.org/wiki/Distance_from_a_point_to_a_line>
        s3 = fdot((x3, y3, x1, y1), y2 - y1, -x21, -y2, x2)
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

        elif c and fdot((x32, y1 - y2), y3 - y2, -x21) < 0:
            # colinear u-turn: x3, y3 not on the
            # opposite side of x2, y2 as x1, y1
            raise CrossError('%s %s: %r' % ('colinear', 'point', ll))

        x1, y1, x2, y2, x21 = x2, y2, x3, y3, x32

    return True  # all points on the same side


def isenclosedBy(latlon, points, wrap=False):  # MCCABE 14
    '''Determine whether a point is enclosed by a polygon.

       @param latlon: The point (I{LatLon} or 2-tuple (lat, lon)).
       @param points: The polygon points (I{LatLon}s).
       @keyword wrap: Wrap lat-, wrap and unroll longitudes (C{bool}).

       @return: C{True} if I{latlon} is inside the polygon, C{False}
                otherwise.

       @raise TypeError: Some I{points} are not I{LatLon}.

       @raise ValueError: Insufficient number of I{points} or invalid
                          I{latlon}.

       @see: L{sphericalNvector.LatLon.isEnclosedBy},
             L{sphericalTrigonometry.LatLon.isEnclosedBy} and
             U{MultiDop GeogContainPt<http://GitHub.com/NASA/MultiDop>}
             (U{Shapiro et al. 2009, JTECH
             <http://journals.AMetSoc.org/doi/abs/10.1175/2009JTECHA1256.1>}
             and U{Potvin et al. 2012, JTECH
             <http://journals.AMetSoc.org/doi/abs/10.1175/JTECH-D-11-00019.1>}).
    '''
    pts = LatLon2psxy(points, closed=True, radius=None, wrap=wrap)

    def _xy(i):
        x, y, _ = pts[i]
        if not wrap:
            x %= 360.0
            if x < (x0 - 180):
                x += 360
            elif x >= (x0 + 180):
                x -= 360
        return x, y

    try:
        y0, x0 = latlon.lat, latlon.lon
    except AttributeError:
        try:
            y0, x0 = latlon[:2]
        except (IndexError, TypeError, ValueError):
            raise ValueError('%s invalid: %r' % ('latlon', latlon))

    if wrap:
        x0, y0 = wrap180(x0), wrap90(y0)
    else:
        x0 %= 360.0

    n = len(pts)
    e = m = False
    s = Fsum()

    x1, y1 = _xy(n-1)
    for i in range(n):
        x2, y2 = _xy(i)
        dx, x2 = unroll180(x1, x2, wrap=wrap)
        # determine if polygon edge (x1, y1)..(x2, y2) straddles
        # point (lat, lon) or is on boundary, but do not count
        # edges on boundary as more than one crossing
        if abs(dx) < 180 and (x1 < x0 <= x2 or x2 < x0 <= x1):
            m = not m
            dy = (x0 - x1) * (y2 - y1) - (y0 - y1) * dx
            if (dy > 0 and dx >= 0) or (dy < 0 and dx <= 0):
                e = not e

        s.fadd(sin(radians(y2)))
        x1, y1 = x2, y2

    # An odd number of meridian crossings means, the polygon
    # contains a pole.  Assume that is the pole on the hemisphere
    # containing the polygon mean point.  If the polygon does
    # contain the North Pole, flip the result.
    if m and s.fsum() > 0:
        e = not e

    return e


def isenclosedby(latlon, points, wrap=False):
    '''DEPRECATED, use function L{isenclosedBy}.
    '''
    return isenclosedBy(latlon, points, wrap=wrap)


def nearestOn3(point, points, closed=False, **options):
    '''Locate the point on a polygon closest to an other point.

       If the given point is within the extent of a polygon edge,
       the closest point is on that edge, otherwise the closest
       point is the nearest of that edge's end points.

       Distances are approximated by function L{equirectangular_},
       subject to the supplied I{options}.

       @param point: The other, reference point (C{LatLon}).
       @param points: The polygon points (C{LatLon}s).
       @keyword closed: Optionally, close the polygon (C{bool}).
       @keyword options: Optional keyword arguments for function
                         L{equirectangular_}.

       @return: 3-Tuple (lat, lon, I{distance}) all in C{degrees}.
                The I{distance} is the L{equirectangular_} distance
                between the closest and the reference I{point} in
                C{degrees}.

       @raise LimitError: Lat- and/or longitudinal delta exceeds
                          I{limit}, see function L{equirectangular_}.

       @raise TypeError: Some I{points} are not I{LatLon}.

       @raise ValueError: Insufficient number of I{points}.

       @see: Function L{degrees2m} to convert C{degrees} to C{meter}.
    '''
    def _d2yx(p2, p1, u):
        # equirectangular_ returns a 4-tuple (distance in
        # degrees squared, delta lat, delta lon, p2.lon
        # unroll/wrap); the previous p2.lon unroll/wrap
        # is also applied to the next edge's p1.lon
        return equirectangular_(p1.lat, p1.lon + u,
                                p2.lat, p2.lon, **options)

    # point (x, y) on axis rotated ccw by angle a:
    #   x' = y * sin(a) + x * cos(a)
    #   y' = y * cos(a) - x * sin(a)
    #
    # distance (w) along and perpendicular (h) to
    # a line thru point (dx, dy) and the origin:
    #   w = (y * dy + x * dx) / hypot(dx, dy)
    #   h = (y * dx - x * dy) / hypot(dx, dy)
    #
    # closest point on that line thru (dx, dy):
    #   xc = dx * w / hypot(dx, dy)
    #   yc = dy * w / hypot(dx, dy)
    # or
    #   xc = dx * f
    #   yc = dy * f
    # with
    #   f = w / hypot(dx, dy)
    # or
    #   f = (y * dy + x * dx) / (dx**2 + dy**2)

    n, points = points2(points, closed=closed)

    i, m = _imdex2(closed, n)
    p2 = c = points[i]
    u2 = u = 0
    d, dy, dx, _ = _d2yx(point, p2, u2)
    for i in range(m, n):
        p1, u1, p2 = p2, u2, points[i]
        # iff wrapped, unroll lon1 (actually previous
        # lon2) like function unroll180/-PI would've
        d21, y21, x21, u2 = _d2yx(p2, p1, u1)
        if d21 > EPS:
            # distance point to p1
            d2, y01, x01, _ = _d2yx(point, p1, u1)
            if d2 > EPS:
                w2 = y01 * y21 + x01 * x21
                if w2 > 0:
                    if w2 < d21:
                        # closest is between p1 and p2, use
                        # original delta's, not y21 and x21
                        f = w2 / d21
                        p1 = LatLon_(favg(p1.lat, p2.lat, f=f),
                                     favg(p1.lon, p2.lon + u2, f=f))
                        u1 = 0
                    else:  # p2 is closest
                        p1, u1 = p2, u2
                    d2, y01, x01, _ = _d2yx(point, p1, u1)
            if d2 < d:  # p1 is closer
                c, u, d, dy, dx = p1, u1, d2, y01, x01
    return c.lat, c.lon + u, hypot(dy, dx)


def perimeterOf(points, closed=False, adjust=True, radius=R_M, wrap=True):
    '''Approximate the perimeter of a polygon.

       @param points: The polygon points (I{LatLon}s).
       @keyword closed: Optionally, close the polygon (C{bool}).
       @keyword adjust: Adjust the wrapped, unrolled longitudinal delta
                        by the cosine of the mean latitude (C{bool}).
       @keyword radius: Optional, mean earth radius (C{meter}).
       @keyword wrap: Wrap lat-, wrap and unroll longitudes (C{bool}).

       @return: Approximate perimeter (C{meter}, same units as I{radius}).

       @raise TypeError: Some I{points} are not I{LatLon}.

       @raise ValueError: Insufficient number of I{points}.

       @note: This perimeter is based on the L{equirectangular_}
              distance approximation and is ill-suited for regions
              exceeding several hundred Km or Miles or with
              near-polar latitudes.

       @see: L{sphericalTrigonometry.perimeterOf} and
             L{ellipsoidalKarney.perimeterOf}.
    '''
    pts = LatLon2psxy(points, closed=closed, radius=None, wrap=False)

    def _degs(n, pts, closed):  # angular edge lengths in degrees
        i, m = _imdex2(closed, n)
        x1, y1, _ = pts[i]
        u = 0  # previous x2's unroll/wrap
        for i in range(m, n):
            x2, y2, _ = pts[i]
            # apply previous x2's unroll/wrap to new x1
            _, dy, dx, u = equirectangular_(y1, x1 + u, y2, x2,
                                            adjust=adjust,
                                            limit=None,
                                            wrap=wrap)
            yield hypot(dx, dy)
            x1, y1 = x2, y2

    d = fsum(_degs(len(pts), pts, closed))
    return radians(d) * float(radius)

# **) MIT License
#
# Copyright (C) 2016-2018 -- mrJean1 at Gmail dot com
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
