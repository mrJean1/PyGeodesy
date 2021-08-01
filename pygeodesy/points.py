
# -*- coding: utf-8 -*-

u'''Utilities for point lists, tuples, etc.

Functions to handle collections and sequences of C{LatLon} points
specified as 2-d U{NumPy<https://www.NumPy.org>}, C{arrays} or tuples as
C{LatLon} or as C{pseudo-x/-y} pairs.

C{NumPy} arrays are assumed to contain rows of points with a lat-, a
longitude -and possibly other- values in different columns.  While
iterating over the array rows, create an instance of a given C{LatLon}
class "on-the-fly" for each row with the row's lat- and longitude.

The original C{NumPy} array is read-accessed only and never duplicated,
except to create a I{subset} of the original array.

For example, to process a C{NumPy} array, wrap the array by instantiating
class L{Numpy2LatLon} and specifying the column index for the lat- and
longitude in each row.  Then, pass the L{Numpy2LatLon} instance to any
L{pygeodesy} function or method accepting a I{points} argument.

Similarly, class L{Tuple2LatLon} is used to instantiate a C{LatLon}
for each 2+tuple in a list, tuple or sequence of such 2+tuples from
the index for the lat- and longitude index in each 2+tuple.
'''

from pygeodesy.basics import isclass, isint, isnear0, isscalar, \
                             issequence, issubclassof, map1, \
                             _Sequence, _xcopy, _xinstanceof
from pygeodesy.datums import _spherical_datum
from pygeodesy.dms import F_D, latDMS, lonDMS, parseDMS2
from pygeodesy.errors import CrossError, crosserrors, _IndexError, \
                            _IsnotError, _TypeError, _ValueError, \
                            _xkwds, _xkwds_pop
from pygeodesy.fmath import favg, fdot, Fsum, fsum, hypot
from pygeodesy.formy import _bearingTo2, equirectangular_, latlon2n_xyz
from pygeodesy.interns import EPS, EPS1, NN, PI_2, R_M, _angle_, \
                             _colinear_, _COMMASPACE_, _DASH_, _DEQUALSPACED_, \
                             _distance_, _ELLIPSIS_, _height_, _lat_, _lon_, \
                             _near_, _not_, _point_, _SPACE_, _UNDER_, _valid_, \
                             _0_0, _0_5, _1_0, _3_0, _90_0, _180_0, _360_0
from pygeodesy.iters import LatLon2PsxyIter, PointsIter, points2
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import classname, nameof, notImplemented, notOverloaded, \
                           _NamedTuple, _xnamed, _xother3, _xotherError
from pygeodesy.namedTuples import Bounds2Tuple, Bounds4Tuple, \
                                  LatLon2Tuple, NearestOn3Tuple, \
                                  PhiLam2Tuple, Point3Tuple, Vector4Tuple
from pygeodesy.nvectorBase import NvectorBase, _N_vector_
from pygeodesy.props import deprecated_method, Property_RO, property_doc_, \
                            property_RO
from pygeodesy.streprs import Fmt, hstr, instr, pairs
from pygeodesy.units import Degrees, Lat, Lon, Number_, Radius, Scalar_
from pygeodesy.utily import atan2b, degrees90, degrees180, degrees2m, \
                            unroll180, unrollPI, wrap90, wrap180

from math import cos, fmod, radians, sin

__all__ = _ALL_LAZY.points
__version__ = '21.07.31'

_fin_   = 'fin'
_ilat_  = 'ilat'
_ilon_  = 'ilon'
_ncols_ = 'ncols'
_nrows_ = 'nrows'


class LatLon_(object):  # XXX imported by heights._HeightBase.height
    '''Low-overhead C{LatLon} class for L{Numpy2LatLon} and L{Tuple2LatLon}.
    '''
    # __slots__ efficiency is voided if the __slots__ class attribute
    # is used in a subclass of a class with the traditional __dict__,
    # see <https://docs.Python.org/2/reference/datamodel.html#slots>
    # and __slots__ must be repeated in sub-classes, see "Problems
    # with __slots__" in Luciano Ramalho, "Fluent Python", page
    # 276+, O'Reilly, 2016, also at <https://Books.Google.ie/
    #   books?id=bIZHCgAAQBAJ&lpg=PP1&dq=fluent%20python&pg=
    #   PT364#v=onepage&q=“Problems%20with%20__slots__”&f=false>
    #
    # __slots__ = (_lat_, _lon_, _height_, _datum_, _name_)
    # Property_RO = property_RO  # no __dict__ with __slots__
    #
    # However, sys.getsizeof(LatLon_(1, 2)) is 72-88 with __slots__
    # and only 48-64 bytes without in Python 2.7.18+ and Python 3+.

    def __init__(self, lat, lon, name=NN, height=0, datum=None):
        '''Creat a new, mininal, low-overhead L{LatLon_} instance,
           without height and datum.

           @arg lat: Latitude (C{degrees}).
           @arg lon: Longitude (C{degrees}).
           @kwarg name: Optional name (C{str}).
           @kwarg height: Optional height (C{float} or C{int}).
           @kwarg datum: Optional datum (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2}, L{a_f2Tuple} or I{scalar}
                         radius) or C{None}.

           @raise TypeError: Invalid B{C{datum}}.

           @note: The lat- and longitude are taken as-given,
                  un-clipped and un-validated .
        '''
        try:  # most common use case
            self.lat, self.lon = float(lat), float(lon)
        except (TypeError, ValueError):
            self.lat, self.lon = parseDMS2(lat, lon, clipLat=0, clipLon=0)  # PYCHOK LatLon2Tuple
        self.name   = str(name)
        self.height = height
        self.datum  = datum if datum is None else \
                     _spherical_datum(datum, name=self.name)

    def __eq__(self, other):
        return isinstance(other, LatLon_) and \
                          other.lat == self.lat and \
                          other.lon == self.lon

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return self.toRepr()

    def __str__(self):
        return self.toStr()

    def classof(self, *args, **kwds):
        '''Instantiate this very class.

           @arg args: Optional, positional arguments.
           @kwarg kwds: Optional, keyword arguments.

           @return: New instance (C{self.__class__}).
        '''
        return _xnamed(self.__class__(*args, **kwds), self.name)

    def copy(self, deep=False):
        '''Make a shallow or deep copy of this instance.

           @kwarg deep: If C{True} make a deep, otherwise a
                        shallow copy (C{bool}).

           @return: The copy (C{This class} or subclass thereof).
        '''
        return _xcopy(self, deep=deep)

    def heightStr(self, prec=-2):
        '''Return a string for the height B{C{height}}.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).

           @see: Function L{hstr}.
        '''
        return hstr(self.height, prec=prec)

    def isequalTo(self, other, eps=None):
        '''Compare this point with an other point, I{ignoring} height.

           @arg other: The other point (C{LatLon}).
           @kwarg eps: Tolerance for equality (C{degrees}).

           @return: C{True} if both points are identical,
                    I{ignoring} height, C{False} otherwise.

           @raise UnitError: Invalid B{C{eps}}.
        '''
        self.others(other)

        if eps:
            return max(abs(self.lat - other.lat),
                       abs(self.lon - other.lon)) < Scalar_(eps=eps)
        else:
            return self.lat == other.lat and \
                   self.lon == other.lon

    @Property_RO
    def latlon(self):
        '''Get the lat- and longitude in C{degrees} (L{LatLon2Tuple}C{(lat, lon)}).
        '''
        return LatLon2Tuple(self.lat, self.lon, name=self.name)

    @Property_RO
    def latlonheight(self):
        '''Get the lat-, longitude and height (L{LatLon3Tuple}C{(lat, lon, height)}).
        '''
        return self.latlon.to3Tuple(self.height)

    @Property_RO
    def _N_vector(self):
        '''(INTERNAL) Get the minimal, low-overhead (C{nvectorBase._N_vector_})
        '''
        return _N_vector_(*latlon2n_xyz(self.lat, self.lon),
                           h=self.height, name=self.name)

    def others(self, *other, **name_other_up):  # see .named._namedBase.others
        '''Refined class comparison.

           @arg other: The other instance (any C{type}).
           @kwarg name_other_up: Overriding C{name=other} and C{up=1}
                                 keyword arguments.

           @return: The B{C{other}} if compatible.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        other, name, up = _xother3(self, other, **name_other_up)
        if isinstance(other, self.__class__) or (hasattr(other, _lat_)
                                             and hasattr(other, _lon_)):
            return other
        raise _xotherError(self, other, name=name, up=up + 1)

    @Property_RO
    def philam(self):
        '''Get the lat- and longitude in C{radians} (L{PhiLam2Tuple}C{(phi, lam)}).
        '''
        return PhiLam2Tuple(radians(self.lat), radians(self.lon), name=self.name)

    @Property_RO
    def philamheight(self):
        '''Get the lat-, longitude in C{radians} and height (L{PhiLam3Tuple}C{(phi, lam, height)}).
        '''
        return self.philam.to3Tuple(self.height)

    @deprecated_method
    def points(self, points, closed=False, base=None):  # PYCHOK no cover
        '''DEPRECATED, use method C{points2}.'''
        return points2(points, closed=closed, base=base)

    def points2(self, points, closed=False, base=None):
        '''Check a path or polygon represented by points.

           @arg points: The path or polygon points (C{LatLon}[])
           @kwarg closed: Optionally, consider the polygon closed,
                          ignoring any duplicate or closing final
                          B{C{points}} (C{bool}).
           @kwarg base: Optionally, check all B{C{points}} against
                        this base class, if C{None} don't check.

           @return: A L{Points2Tuple}C{(number, points)} with the number
                    of points and the points C{list} or C{tuple}.

           @raise PointsError: Insufficient number of B{C{points}}.

           @raise TypeError: Some B{C{points}} are not B{C{base}}.
        '''
        return points2(points, closed=closed, base=base)

    def PointsIter(self, points, loop=0, dedup=False):
        '''Return a points iterator.

           @arg points: The path or polygon points (C{LatLon}[])
           @kwarg loop: Number of loop-back points (non-negative C{int}).
           @kwarg dedup: Skip duplicate points (C{bool}).

           @return: A new C{PointsIter} iterator.

           @raise PointsError: Insufficient number of B{C{points}}.
        '''
        return PointsIter(points, loop=loop, base=self, dedup=dedup)

    @deprecated_method
    def to2ab(self):  # PYCHOK no cover
        '''DEPRECATED, use property L{philam}.'''
        return self.philam

    def toNvector(self, h=None, Nvector=NvectorBase, **Nvector_kwds):
        '''Convert this point to C{n-vector} (normal to the earth's
           surface) components, I{including height}.

           @kwarg h: Optional height, overriding this point's height
                     (C{meter}).
           @kwarg Nvector: Optional class to return the C{n-vector}
                           components (C{Nvector}) or C{None}.
           @kwarg Nvector_kwds: Optional, additional B{C{Nvector}} keyword
                                arguments, ignored if C{B{Nvector} is None}.

           @return: The C{n-vector} components B{C{Nvector}} or if
                    B{C{Nvector}} is C{None}, a L{Vector4Tuple}C{(x,
                    y, z, h)}.

           @raise TypeError: Invalid B{C{Nvector}} or B{C{Nvector_kwds}}.
        '''
        x, y, z = latlon2n_xyz(self.lat, self.lon)
        r = Vector4Tuple(x, y, z, self.height if h is None else h)
        if Nvector is not None:
            r = Nvector(x, y, z, **_xkwds(Nvector_kwds, h=r.h, ll=self))
        return _xnamed(r, self.name)

    def toRepr(self, **kwds):
        '''This L{LatLon_} as a string "class(<degrees>, ...)".

           @kwarg kwds: Optional, keyword arguments.

           @return: Class instance (C{str}).
        '''
        _ = _xkwds_pop(kwds, std=None)  # PYCHOK std unused
        return Fmt.PAREN(classname(self), self.toStr(**kwds))

    def toStr(self, form=F_D, prec=6, sep=_COMMASPACE_, **kwds):
        '''This L{LatLon_} as a string "<degrees>, <degrees>".

           @kwarg form: Optional format, F_D, F_DM, F_DMS for
                        deg°, deg°min′, deg°min′sec″ (C{str}).
           @kwarg prec: Optional number of decimal digits (0..8 or C{None}).
           @kwarg sep: Optional separator to join (C{str}).
           @kwarg kwds: Optional, keyword arguments.

           @return: Instance (C{str}).
        '''
        t = (latDMS(self.lat, form=form, prec=prec),
             lonDMS(self.lon, form=form, prec=prec))
        if self.height:
            t += (self.heightStr(),)
        if self.name:
            t += (repr(self.name),)
        if kwds:
            t += pairs(kwds, prec=prec)
        return sep.join(t)

    @deprecated_method
    def toStr2(self, **kwds):  # PYCHOK no cover
        '''DEPRECATED, used method L{toRepr}.'''
        return self.toRepr(**kwds)


try:
    _LatLon_attrs = LatLon_.__slots__  # PYCHOK no __slots__
except AttributeError:
    _LatLon_attrs = tuple(LatLon_(0, 0).__dict__.keys())


def _isLatLon_(LL):
    '''(INTERANL) Check attributes of class C{LL}.
    '''
    return issubclassof(LL, LatLon_) or (isclass(LL) and
                                     all(hasattr(LL, a) for a in _LatLon_attrs))


class _Basequence(_Sequence):  # immutable, on purpose
    '''(INTERNAL) Base class.
    '''
    _array    =  []
    _epsilon  =  EPS
    _itemname = _point_

    def _contains(self, point):
        '''(INTERNAL) Check for a matching point.
        '''
        for _ in self._findall(point, ()):
            return True
        return False

    def copy(self, deep=False):
        '''Make a shallow or deep copy of this instance.

           @kwarg deep: If C{True} make a deep, otherwise a
                        shallow copy (C{bool}).

           @return: The copy (C{This class} or subclass thereof).
        '''
        return _xcopy(self, deep=deep)

    def _count(self, point):
        '''(INTERNAL) Count the number of matching points.
        '''
        n = 0
        for _ in self._findall(point, ()):
            n += 1
        return n

    @property_doc_(''' the equality tolerance (C{float}).''')
    def epsilon(self):
        '''Get the tolerance for equality tests (C{float}).
        '''
        return self._epsilon

    @epsilon.setter  # PYCHOK setter!
    def epsilon(self, tol):
        '''Set the tolerance for equality tests.

           @arg tol: New tolerance (C{scalar}).

           @raise TypeError: Non-scalar B{C{tol}}.

           @raise ValueError: Out-of-bounds B{C{tol}}.
        '''
        self._epsilon = Scalar_(tol, name='tolerance')

    def _find(self, point, start_end):
        '''(INTERNAL) Find the first matching point index.
        '''
        for i in self._findall(point, start_end):
            return i
        return -1

    def _findall(self, point, start_end):  # PYCHOK no cover
        '''(INTERNAL) I{Must be implemented/overloaded}.
        '''
        notImplemented(self, point, start_end)

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
        raise _IndexError(self._itemname, point, txt=_not_('found'))

    @property_RO
    def isNumpy2(self):  # PYCHOK no cover
        '''Is this a Numpy2 wrapper?
        '''
        return False  # isinstance(self, (Numpy2LatLon, ...))

    @property_RO
    def isPoints2(self):  # PYCHOK no cover
        '''Is this a LatLon2 wrapper/converter?
        '''
        return False  # isinstance(self, (LatLon2psxy, ...))

    @property_RO
    def isTuple2(self):  # PYCHOK no cover
        '''Is this a Tuple2 wrapper?
        '''
        return False  # isinstance(self, (Tuple2LatLon, ...))

    def _iter(self):
        '''(INTERNAL) Yield all points.
        '''
        for i in range(len(self)):
            yield self.point(self._array[i])

    def point(self, *attrs):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}, see function C{notOverloaded}.

           @arg attrs: Optional arguments.
        '''
        notOverloaded(self, *attrs)

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
            raise _ValueError(step=step)
        return range(start, end, step)

    def _repr(self):
        '''(INTERNAL) Return a string representation.
        '''
        # XXX use Python 3+ reprlib.repr
        t =  repr(self._array[:1])  # only first row
        t = _SPACE_(t[:-1], _ELLIPSIS_, Fmt.SQUARE(t[-1:], len(self)))
        t = _SPACE_.join(t.split())  # coalesce spaces
        return instr(self, t, **self._slicekwds())

    def _reversed(self):  # PYCHOK false
        '''(INTERNAL) Yield all points in reverse order.
        '''
        for i in range(len(self) - 1, -1, -1):
            yield self.point(self._array[i])

    def _rfind(self, point, start_end):
        '''(INTERNAL) Find the last matching point index.
        '''
        def _r3(start=None, end=None, step=-1):
            return (start, end, step)  # PYCHOK returns

        for i in self._findall(point, _r3(*start_end)):
            return i
        return -1

    def _slicekwds(self):  # PYCHOK no cover
        '''(INTERNAL) I{Should be overloaded}.
        '''
        return {}

    def _zeros(self, *zeros):
        '''(INTERNAL) Check for near-zero values.
        '''
        return all(abs(z) <= self._epsilon for z in zeros)


class _Array2LatLon(_Basequence):  # immutable, on purpose
    '''Base class for Numpy2LatLon or Tuple2LatLon.
    '''
    _array  = ()
    _ilat   = 0  # row column index
    _ilon   = 0  # row column index
    _LatLon = LatLon_  # default
    _shape  = ()

    def __init__(self, array, ilat=0, ilon=1, LatLon=None, shape=()):
        '''Handle a C{NumPy} or C{Tuple} array as a sequence of C{LatLon} points.
        '''
        ais = (_ilat_, ilat), (_ilon_, ilon)

        if len(shape) != 2 or shape[0] < 1 or shape[1] < len(ais):
            raise _IndexError('array.shape', shape)

        self._array = array
        self._shape = Shape2Tuple(*shape)

        if LatLon:  # check the point class
            if _isLatLon_(LatLon):
                self._LatLon = LatLon
            else:
                raise _IsnotError(_valid_, LatLon=LatLon)

        # check the attr indices
        for n, (ai, i) in enumerate(ais):
            if not isint(i):
                raise _IsnotError(int.__name__, **{ai: i})
            i = int(i)
            if not 0 <= i < shape[1]:
                raise _ValueError(ai, i)
            for aj, j in ais[:n]:
                if int(j) == i:
                    raise _ValueError(_DEQUALSPACED_(ai, aj, i))
            setattr(self, NN(_UNDER_, ai), i)

    def __contains__(self, latlon):
        '''Check for a specific lat-/longitude.

           @arg latlon: Point (C{LatLon}, L{LatLon2Tuple} or 2-tuple
                        C{(lat, lon)}).

           @return: C{True} if B{C{latlon}} is present, C{False} otherwise.

           @raise TypeError: Invalid B{C{latlon}}.
        '''
        return self._contains(latlon)

    def __getitem__(self, index):
        '''Return row[index] as C{LatLon} or return a L{Numpy2LatLon} slice.
        '''
        return self._getitem(index)

    def __iter__(self):
        '''Yield rows as C{LatLon}.
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
        '''Yield rows as C{LatLon} in reverse order.
        '''
        return self._reversed()

    __str__ = __repr__

    def count(self, latlon):
        '''Count the number of rows with a specific lat-/longitude.

           @arg latlon: Point (C{LatLon}, L{LatLon2Tuple} or 2-tuple
                        C{(lat, lon)}).

           @return: Count (C{int}).

           @raise TypeError: Invalid B{C{latlon}}.
        '''
        return self._count(latlon)

    def find(self, latlon, *start_end):
        '''Find the first row with a specific lat-/longitude.

           @arg latlon: Point (C{LatLon}) or 2-tuple (lat, lon).
           @arg start_end: Optional C{[start[, end]]} index (integers).

           @return: Index or -1 if not found (C{int}).

           @raise TypeError: Invalid B{C{latlon}}.
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
                raise _IsnotError(_valid_, latlon=latlon)

        for i in self._range(*start_end):
            row = self._array[i]
            if self._zeros(row[self._ilat] - lat,
                           row[self._ilon] - lon):
                yield i

    def findall(self, latlon, *start_end):
        '''Yield indices of all rows with a specific lat-/longitude.

           @arg latlon: Point (C{LatLon}, L{LatLon2Tuple} or 2-tuple
                        C{(lat, lon)}).
           @arg start_end: Optional C{[start[, end]]} index (C{int}).

           @return: Indices (C{iterable}).

           @raise TypeError: Invalid B{C{latlon}}.
        '''
        return self._findall(latlon, start_end)

    def index(self, latlon, *start_end):  # PYCHOK Python 2- issue
        '''Find index of the first row with a specific lat-/longitude.

           @arg latlon: Point (C{LatLon}, L{LatLon2Tuple} or 2-tuple
                        C{(lat, lon)}).
           @arg start_end: Optional C{[start[, end]]} index (C{int}).

           @return: Index (C{int}).

           @raise IndexError: Point not found.

           @raise TypeError: Invalid B{C{latlon}}.
        '''
        return self._index(latlon, start_end)

    @Property_RO
    def ilat(self):
        '''Get the latitudes column index (C{int}).
        '''
        return self._ilat

    @Property_RO
    def ilon(self):
        '''Get the longitudes column index (C{int}).
        '''
        return self._ilon

#   next = __iter__

    def point(self, row):  # PYCHOK *attrs
        '''Instantiate a point C{LatLon}.

           @arg row: Array row (numpy.array).

           @return: Point (C{LatLon}).
        '''
        return self._LatLon(row[self._ilat], row[self._ilon])

    def rfind(self, latlon, *start_end):
        '''Find the last row with a specific lat-/longitude.

           @arg latlon: Point (C{LatLon}, L{LatLon2Tuple} or 2-tuple
                        C{(lat, lon)}).
           @arg start_end: Optional C{[start[, end]]} index (C{int}).

           @note: Keyword order, first stop, then start.

           @return: Index or -1 if not found (C{int}).

           @raise TypeError: Invalid B{C{latlon}}.
        '''
        return self._rfind(latlon, start_end)

    def _slicekwds(self):
        '''(INTERNAL) Slice kwds.
        '''
        return dict(ilat=self._ilat, ilon=self._ilon)

    @Property_RO
    def shape(self):
        '''Get the shape of the C{NumPy} array or the C{Tuples} as
           L{Shape2Tuple}C{(nrows, ncols)}.
        '''
        return self._shape

    def _subset(self, indices):  # PYCHOK no cover
        '''(INTERNAL) I{Must be implemented/overloaded}.
        '''
        notImplemented(self, indices)

    def subset(self, indices):
        '''Return a subset of the C{NumPy} array.

           @arg indices: Row indices (C{range} or C{int}[]).

           @note: A C{subset} is different from a C{slice} in 2 ways:
                  (a) the C{subset} is typically specified as a list of
                  (un-)ordered indices and (b) the C{subset} allocates
                  a new, separate C{NumPy} array while a C{slice} is
                  just an other C{view} of the original C{NumPy} array.

           @return: Sub-array (C{numpy.array}).

           @raise IndexError: Out-of-range B{C{indices}} value.

           @raise TypeError: If B{C{indices}} is not a C{range}
                             nor an C{int}[].
        '''
        if not issequence(indices, tuple):  # NO tuple, only list
            # and range work properly to get Numpy array sub-sets
            raise _IsnotError(_valid_, indices=type(indices))

        n = len(self)
        for i, v in enumerate(indices):
            if not isint(v):
                raise _TypeError(Fmt.SQUARE(indices=i), v)
            elif not 0 <= v < n:
                raise _IndexError(Fmt.SQUARE(indices=i), v)

        return self._subset(indices)


class LatLon2psxy(_Basequence):
    '''Wrapper for C{LatLon} points as "on-the-fly" pseudo-xy coordinates.
    '''
    _closed = False
    _len    = 0
    _deg2m  = None  # default, keep degrees
    _radius = None
    _wrap   = True

    def __init__(self, latlons, closed=False, radius=None, wrap=True):
        '''Handle C{LatLon} points as pseudo-xy coordinates.

           @note: The C{LatLon} latitude is considered the I{pseudo-y}
                  and longitude the I{pseudo-x} coordinate, likewise
                  for L{LatLon2Tuple}.  However, 2-tuples C{(x, y)} are
                  considered as I{(longitude, latitude)}.

           @arg latlons: Points C{list}, C{sequence}, C{set}, C{tuple},
                         etc. (C{LatLon[]}).
           @kwarg closed: Optionally, close the polygon (C{bool}).
           @kwarg radius: Mean earth radius (C{meter}).
           @kwarg wrap: Wrap lat- and longitudes (C{bool}).

           @raise PointsError: Insufficient number of B{C{latlons}}.

           @raise TypeError: Some B{C{points}} are not B{C{base}}.
        '''
        self._closed = closed
        self._len, self._array = points2(latlons, closed=closed)
        if radius:
            self._radius = r = Radius(radius)
            self._deg2m  = degrees2m(_1_0, r)
        self._wrap = wrap

    def __contains__(self, xy):
        '''Check for a matching point.

           @arg xy: Point (C{LatLon}, L{LatLon2Tuple} or 2-tuple
                    C{(x, y)}) in (C{degrees}.

           @return: C{True} if B{C{xy}} is present, C{False} otherwise.

           @raise TypeError: Invalid B{C{xy}}.
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

           @arg xy: Point (C{LatLon}, L{LatLon2Tuple} or 2-tuple
                    C{(x, y)}) in (C{degrees}.

           @return: Count (C{int}).

           @raise TypeError: Invalid B{C{xy}}.
        '''
        return self._count(xy)

    def find(self, xy, *start_end):
        '''Find the first matching point.

           @arg xy: Point (C{LatLon}, L{LatLon2Tuple} or 2-tuple
                    C{(x, y)}) in (C{degrees}.
           @arg start_end: Optional C{[start[, end]]} index (C{int}).

           @return: Index or -1 if not found (C{int}).

           @raise TypeError: Invalid B{C{xy}}.
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
                x, y = xy[:2]
            except (IndexError, TypeError, ValueError):
                raise _IsnotError(_valid_, xy=xy)

            def _3xyll(ll):  # PYCHOK expected
                return self.point(ll)

        for i in self._range(*start_end):
            xi, yi, _ = _3xyll(self._array[i])
            if self._zeros(xi - x, yi - y):
                yield i

    def findall(self, xy, *start_end):
        '''Yield indices of all matching points.

           @arg xy: Point (C{LatLon}, L{LatLon2Tuple} or 2-tuple
                    C{(x, y)}) in (C{degrees}.
           @arg start_end: Optional C{[start[, end]]} index (C{int}).

           @return: Indices (C{iterator}).

           @raise TypeError: Invalid B{C{xy}}.
        '''
        return self._findall(xy, start_end)

    def index(self, xy, *start_end):  # PYCHOK Python 2- issue
        '''Find the first matching point.

           @arg xy: Point (C{LatLon}) or 2-tuple (x, y) in (C{degrees}).
           @arg start_end: Optional C{[start[, end]]} index (C{int}).

           @return: Index (C{int}).

           @raise IndexError: Point not found.

           @raise TypeError: Invalid B{C{xy}}.
        '''
        return self._index(xy, start_end)

    @property_RO
    def isPoints2(self):
        '''Is this a LatLon2 wrapper/converter?
        '''
        return True  # isinstance(self, (LatLon2psxy, ...))

    def point(self, ll):  # PYCHOK *attrs
        '''Create a pseudo-xy.

           @arg ll: Point (C{LatLon}).

           @return: An L{Point3Tuple}C{(x, y, ll)}.
        '''
        x, y = ll.lon, ll.lat  # note, x, y = lon, lat
        if self._wrap:
            x, y = wrap180(x), wrap90(y)
        if self._deg2m:  # convert degrees to meter (or radians)
            x *= self._deg2m
            y *= self._deg2m
        return Point3Tuple(x, y, ll)

    def rfind(self, xy, *start_end):
        '''Find the last matching point.

           @arg xy: Point (C{LatLon}, L{LatLon2Tuple} or 2-tuple
                    C{(x, y)}) in (C{degrees}.
           @arg start_end: Optional C{[start[, end]]} index (C{int}).

           @return: Index or -1 if not found (C{int}).

           @raise TypeError: Invalid B{C{xy}}.
        '''
        return self._rfind(xy, start_end)

    def _slicekwds(self):
        '''(INTERNAL) Slice kwds.
        '''
        return dict(closed=self._closed, radius=self._radius, wrap=self._wrap)


class NearestOn5Tuple(_NamedTuple):
    '''5-Tuple C{(lat, lon, distance, angle, height)} all in C{degrees},
       except C{height}.  The C{distance} is the L{equirectangular_}
       distance between the closest and the reference B{C{point}} in
       C{degrees}.  The C{angle} from the reference B{C{point}} to
       the closest point is in compass C{degrees360}, see function
       L{compassAngle}.  The C{height} is the (interpolated) height
       at the closest point in C{meter} or C{0}.
    '''
    _Names_ = (_lat_, _lon_, _distance_, _angle_, _height_)
    _Units_ = ( Lat,   Lon,   Degrees,    Degrees, Degrees)


class Numpy2LatLon(_Array2LatLon):  # immutable, on purpose
    '''Wrapper for C{NumPy} arrays as "on-the-fly" C{LatLon} points.
    '''
    def __init__(self, array, ilat=0, ilon=1, LatLon=None):
        '''Handle a C{NumPy} array as a sequence of C{LatLon} points.

           @arg array: C{NumPy} array (C{numpy.array}).
           @kwarg ilat: Optional index of the latitudes column (C{int}).
           @kwarg ilon: Optional index of the longitudes column (C{int}).
           @kwarg LatLon: Optional C{LatLon} class to use (L{LatLon_}).

           @raise IndexError: If B{C{array.shape}} is not (1+, 2+).

           @raise TypeError: If B{C{array}} is not a C{NumPy} array or
                             C{LatLon} is not a class with C{lat}
                             and C{lon} attributes.

           @raise ValueError: If the B{C{ilat}} and/or B{C{ilon}} values
                              are the same or out of range.

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
            raise _IsnotError('NumPy', array=type(array))

        _Array2LatLon.__init__(self, array, ilat=ilat, ilon=ilon,
                                     LatLon=LatLon, shape=s)

    @property_RO
    def isNumpy2(self):
        '''Is this a Numpy2 wrapper?
        '''
        return True  # isinstance(self, (Numpy2LatLon, ...))

    def _subset(self, indices):
        return self._array[indices]  # NumPy special


class Shape2Tuple(_NamedTuple):
    '''2-Tuple C{(nrows, ncols)}, the number of rows and columns,
       both C{int}.
    '''
    _Names_ = (_nrows_, _ncols_)
    _Units_ = ( Number_, Number_)


class Tuple2LatLon(_Array2LatLon):
    '''Wrapper for tuple sequences as "on-the-fly" C{LatLon} points.
    '''
    def __init__(self, tuples, ilat=0, ilon=1, LatLon=None):
        '''Handle a list of tuples, each containing a lat- and longitude
           and perhaps other values as a sequence of C{LatLon} points.

           @arg tuples: The C{list}, C{tuple} or C{sequence} of tuples (C{tuple}[]).
           @kwarg ilat: Optional index of the latitudes value (C{int}).
           @kwarg ilon: Optional index of the longitudes value (C{int}).
           @kwarg LatLon: Optional C{LatLon} class to use (L{LatLon_}).

           @raise IndexError: If C{(len(B{tuples}), min(len(t) for t
                              in B{tuples}))} is not (1+, 2+).

           @raise TypeError: If B{C{tuples}} is not a C{list}, C{tuple}
                             or C{sequence} or if B{C{LatLon}} is not a
                             C{LatLon} with C{lat}, C{lon} and C{name}
                             attributes.

           @raise ValueError: If the B{C{ilat}} and/or B{C{ilon}} values
                              are the same or out of range.

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
        _xinstanceof(list, tuple, tuples=tuples)
        s = len(tuples), min(len(_) for _ in tuples)
        _Array2LatLon.__init__(self, tuples, ilat=ilat, ilon=ilon,
                                     LatLon=LatLon, shape=s)

    @property_RO
    def isTuple2(self):
        '''Is this a Tuple2 wrapper?
        '''
        return True  # isinstance(self, (Tuple2LatLon, ...))

    def _subset(self, indices):
        return type(self._array)(self._array[i] for i in indices)


def _area2(points, adjust, wrap):
    '''(INTERNAL) Approximate the area in radians squared, I{signed}.
    '''
    if adjust:
        # approximate trapezoid by a rectangle, adjusting
        # the top width by the cosine of the latitudinal
        # average and bottom width by some fudge factor
        def _adjust(w, h):
            c = cos(h) if abs(h) < PI_2 else _0_0
            return w * h * (c + 1.2876) * _0_5
    else:
        def _adjust(w, h):  # PYCHOK expected
            return w * h

    # setting radius=1 converts degrees to radians
    Ps = LatLon2PsxyIter(points, wrap=wrap, radius=_1_0, loop=1)
    x1, y1, ll = Ps[0]
    pts = [ll]

    A2 = Fsum()  # trapezoidal area in radians**2
    for i, p in Ps.enumerate(closed=True):
        x2, y2, ll = p
        if 0 < i < 4:
            pts.append(ll)
        w, x2 = unrollPI(x1, x2, wrap=wrap if i else False)
        A2 += _adjust(w, (y2 + y1) * _0_5)
        x1, y1 = x2, y2

    return A2.fsum(), tuple(pts)


def _areaError(pts, near_=NN):  # imported by .ellipsoidalKarney
    '''(INTERNAL) Area issue.
    '''
    t = _ELLIPSIS_(pts[:3], NN)
    return _ValueError(NN(near_, 'zero or polar area'), txt=t)


def areaOf(points, adjust=True, radius=R_M, wrap=True):
    '''Approximate the area of a polygon.

       @arg points: The polygon points (C{LatLon}[]).
       @kwarg adjust: Adjust the wrapped, unrolled longitudinal delta
                      by the cosine of the mean latitude (C{bool}).
       @kwarg radius: Mean earth radius (C{meter}) or C{None}.
       @kwarg wrap: Wrap lat-, wrap and unroll longitudes (C{bool}).

       @return: Approximate area (I{square} C{meter}, same units as
                B{C{radius}} or C{radians} I{squared} if B{C{radius}}
                is C{None}).

       @raise PointsError: Insufficient number of B{C{points}}

       @raise TypeError: Some B{C{points}} are not C{LatLon}.

       @raise ValueError: Invalid B{C{radius}}.

       @note: This area approximation has limited accuracy and is
              ill-suited for regions exceeding several hundred Km
              or Miles or with near-polar latitudes.

       @see: L{sphericalNvector.areaOf}, L{sphericalTrigonometry.areaOf}
             and L{ellipsoidalKarney.areaOf}.
    '''
    a, _ = _area2(points, adjust, wrap)
    return abs(a) if radius is None else (abs(a) * Radius(radius)**2)


def boundsOf(points, wrap=True, LatLon=None):
    '''Determine the lower-left SW and upper-right NE corners of a
       path or polygon.

       @arg points: The path or polygon points (C{LatLon}[]).
       @kwarg wrap: Wrap lat- and longitudes (C{bool}).
       @kwarg LatLon: Optional class to return the C{bounds}
                      corners (C{LatLon}) or C{None}.

       @return: A L{Bounds2Tuple}C{(latlonSW, latlonNE)} as
                B{C{LatLon}}s if B{C{LatLon}} is C{None} a
                L{Bounds4Tuple}C{(latS, lonW, latN, lonE)}.

       @raise PointsError: Insufficient number of B{C{points}}

       @raise TypeError: Some B{C{points}} are not C{LatLon}.

       @see: Function L{quadOf}.

       @example:

        >>> b = LatLon(45,1), LatLon(45,2), LatLon(46,2), LatLon(46,1)
        >>> boundsOf(b)  # False
        >>> 45.0, 1.0, 46.0, 2.0
    '''
    Ps = LatLon2PsxyIter(points, wrap=wrap, loop=1)
    lox, loy, _ = hix, hiy, _ = Ps[0]

    for x, y, _ in Ps.iterate(closed=False):  # [1:]
        if lox > x:
            lox = x
        elif hix < x:
            hix = x

        if loy > y:
            loy = y
        elif hiy < y:
            hiy = y

    return Bounds4Tuple(loy, lox, hiy, hix) if LatLon is None else \
           Bounds2Tuple(LatLon(loy, lox), LatLon(hiy, hix))  # PYCHOK inconsistent


def centroidOf(points, wrap=True, LatLon=None):
    '''Determine the centroid of a polygon.

       @arg points: The polygon points (C{LatLon}[]).
       @kwarg wrap: Wrap lat-, wrap and unroll longitudes (C{bool}).
       @kwarg LatLon: Optional class to return the centroid (C{LatLon})
                      or C{None}.

       @return: Centroid (B{C{LatLon}}) or a L{LatLon2Tuple}C{(lat, lon)}
                if C{B{LatLon} is None}.

       @raise PointsError: Insufficient number of B{C{points}}

       @raise TypeError: Some B{C{points}} are not C{LatLon}.

       @raise ValueError: The B{C{points}} enclose a pole or
                          near-zero area.

       @see: U{Centroid<https://WikiPedia.org/wiki/Centroid#Of_a_polygon>} and
             Paul Bourke's U{Calculating The Area And Centroid Of A Polygon
             <https://www.SEAS.UPenn.edu/~ese502/lab-content/extra_materials/
             Polygon%20Area%20and%20Centroid.pdf>}, 1988.
    '''
    A, X, Y = Fsum(), Fsum(), Fsum()

    # setting radius=1 converts degrees to radians
    Ps = LatLon2PsxyIter(points, wrap=wrap, radius=_1_0, loop=1)
    x1, y1, ll = Ps[0]
    pts = [ll]  # for _areaError
    for i, p in Ps.enumerate(closed=True):
        x2, y2, ll = p
        if 0 < i < 4:
            pts.append(ll)
        if wrap and i != 0:
            _, x2 = unrollPI(x1, x2, wrap=wrap)
        t  = x1 * y2 - x2 * y1
        A += t
        X += t * (x1 + x2)
        Y += t * (y1 + y2)
        # XXX more elaborately:
        # t1, t2 = x1 * y2, -(x2 * y1)
        # A.fadd_(t1, t2)
        # X.fadd_(t1 * x1, t1 * x2, t2 * x1, t2 * x2)
        # Y.fadd_(t1 * y1, t1 * y2, t2 * y1, t2 * y2)
        x1, y1 = x2, y2

    A = A.fsum() * _3_0  # 6.0 / 2.0
    if isnear0(A):
        raise _areaError(pts, near_=_DASH_(_near_, NN))
    Y, X = degrees90(Y.fsum() / A), degrees180(X.fsum() / A)
    return LatLon2Tuple(Y, X) if LatLon is None else LatLon(Y, X)


def _fractional(points, fi, fin=None):  # wrap [n] to [0]
    '''(INTERNAL) Compute point at L{fractional} index.
    '''
    i = int(fi)
    p = points[i]
    r = fi - float(i)
    if r > EPS:
        i += 1
        q  = points[0 if i == fin else i]
        if r < EPS1:
            p = LatLon2Tuple(favg(p.lat, q.lat, f=r),
                             favg(p.lon, q.lon, f=r))
        else:
            p = q
    return p


def fractional(points, fi, LatLon=None, **LatLon_kwds):
    '''Return the point at a given I{fractional} index.

       @arg points: The points (C{LatLon}[], L{Numpy2LatLon}[],
                    L{Tuple2LatLon}[] or C{other}[]).
       @arg fi: The fractional index (L{FIx}, C{float} or C{int}).
       @kwarg LatLon: Optional class to return the I{intermediate},
                      I{fractional} point (C{LatLon}) or C{None}.
       @kwarg LatLon_kwds: Optional B{C{LatLon}} keyword arguments,
                           ignored of C{B{LatLon} is None}.

       @return: A B{C{LatLon}} or if B{C{LatLon}} is C{None}, a
                L{LatLon2Tuple}C{(lat, lon)} for C{B{points}[B{fi}]} if
                I{fractional} index B{C{fi}} is C{int}, otherwise the
                intermediate point between C{B{points}[int(B{fi})]} and
                C{B{points}[int(B{fi}) + 1]} for C{float} I{fractional}
                index B{C{fi}}.

       @raise IndexError: Fractional index B{C{fi}} invalid or
                          B{C{points}} not subscriptable or not
                          closed.

       @see: Class L{FIx} and method L{FIx.fractional}.
    '''
    try:
        if not isscalar(fi) or fi < 0:
            raise IndexError
        p = _fractional(points, fi, fin=getattr(fi, _fin_, 0))  # see .units.FIx
    except (IndexError, TypeError):
        raise _IndexError(fractional.__name__, fi)

    if LatLon and isinstance(p, LatLon2Tuple):
        p = LatLon(*p, **LatLon_kwds)
    return p


def isclockwise(points, adjust=False, wrap=True):
    '''Determine the direction of a path or polygon.

       @arg points: The path or polygon points (C{LatLon}[]).
       @kwarg adjust: Adjust the wrapped, unrolled longitudinal delta
                      by the cosine of the mean latitude (C{bool}).
       @kwarg wrap: Wrap lat-, wrap and unroll longitudes (C{bool}).

       @return: C{True} if B{C{points}} are clockwise, C{False} otherwise.

       @raise PointsError: Insufficient number of B{C{points}}

       @raise TypeError: Some B{C{points}} are not C{LatLon}.

       @raise ValueError: The B{C{points}} enclose a pole or zero area.

       @example:

        >>> f = LatLon(45,1), LatLon(45,2), LatLon(46,2), LatLon(46,1)
        >>> isclockwise(f)  # False
        >>> isclockwise(reversed(f))  # True
    '''
    a, pts = _area2(points, adjust, wrap)
    if a > 0:
        return True
    elif a < 0:
        return False
    # <https://blog.Element84.com/determining-if-a-spherical-polygon-contains-a-pole.html>
    raise _areaError(pts)


def isconvex(points, adjust=False, wrap=True):
    '''Determine whether a polygon is convex.

       @arg points: The polygon points (C{LatLon}[]).
       @kwarg adjust: Adjust the wrapped, unrolled longitudinal delta
                      by the cosine of the mean latitude (C{bool}).
       @kwarg wrap: Wrap lat-, wrap and unroll longitudes (C{bool}).

       @return: C{True} if B{C{points}} are convex, C{False} otherwise.

       @raise CrossError: Some B{C{points}} are colinear.

       @raise PointsError: Insufficient number of B{C{points}}

       @raise TypeError: Some B{C{points}} are not C{LatLon}.

       @example:

        >>> t = LatLon(45,1), LatLon(46,1), LatLon(46,2)
        >>> isconvex(t)  # True

        >>> f = LatLon(45,1), LatLon(46,2), LatLon(45,2), LatLon(46,1)
        >>> isconvex(f)  # False
    '''
    return bool(isconvex_(points, adjust=adjust, wrap=wrap))


def isconvex_(points, adjust=False, wrap=True):
    '''Determine whether a polygon is convex I{and clockwise}.

       @arg points: The polygon points (C{LatLon}[]).
       @kwarg adjust: Adjust the wrapped, unrolled longitudinal delta
                      by the cosine of the mean latitude (C{bool}).
       @kwarg wrap: Wrap lat-, wrap and unroll longitudes (C{bool}).

       @return: C{+1} if B{C{points}} are convex clockwise, C{-1} for
                convex counter-clockwise B{C{points}}, C{0} otherwise.

       @raise CrossError: Some B{C{points}} are colinear.

       @raise PointsError: Insufficient number of B{C{points}}

       @raise TypeError: Some B{C{points}} are not C{LatLon}.

       @example:

        >>> t = LatLon(45,1), LatLon(46,1), LatLon(46,2)
        >>> isconvex_(t)  # +1

        >>> f = LatLon(45,1), LatLon(46,2), LatLon(45,2), LatLon(46,1)
        >>> isconvex_(f)  # 0
    '''
    def _unroll_adjust(x1, y1, x2, y2, w):
        x21, x2 = unroll180(x1, x2, wrap=w)
        if adjust:
            y = radians(y1 + y2) * _0_5
            x21 *= cos(y) if abs(y) < PI_2 else _0_0
        return x21, x2

    c, s = crosserrors(), 0

    Ps = LatLon2PsxyIter(points, wrap=wrap, loop=2)
    x1, y1, _ = Ps[0]
    x2, y2, _ = Ps[1]
    x21, x2 = _unroll_adjust(x1, y1, x2, y2, False)

    for i, p in Ps.enumerate(closed=True):
        x3, y3, ll = p
        x32, x3 = _unroll_adjust(x2, y2, x3, y3, (wrap if i > 1 else False))

        # get the sign of the distance from point
        # x3, y3 to the line from x1, y1 to x2, y2
        # <https://WikiPedia.org/wiki/Distance_from_a_point_to_a_line>
        s3 = fdot((x3, y3, x1, y1), y2 - y1, -x21, -y2, x2)
        if s3 > 0:  # x3, y3 on the right
            if s < 0:  # non-convex
                return 0
            s = +1

        elif s3 < 0:  # x3, y3 on the left
            if s > 0:  # non-convex
                return 0
            s = -1

        elif c and fdot((x32, y1 - y2), y3 - y2, -x21) < 0:
            # colinear u-turn: x3, y3 not on the
            # opposite side of x2, y2 as x1, y1
            t = Fmt.SQUARE(points=i)
            raise CrossError(t, ll, txt=_colinear_)

        x1, y1, x2, y2, x21 = x2, y2, x3, y3, x32

    return s  # all points on the same side


def isenclosedBy(point, points, wrap=False):  # MCCABE 15
    '''Determine whether a point is enclosed by a polygon.

       @arg point: The point (C{LatLon} or 2-tuple C{(lat, lon)}).
       @arg points: The polygon points (C{LatLon}[]).
       @kwarg wrap: Wrap lat-, wrap and unroll longitudes (C{bool}).

       @return: C{True} if B{C{point}} is inside the polygon, C{False}
                otherwise.

       @raise PointsError: Insufficient number of B{C{points}}

       @raise TypeError: Some B{C{points}} are not C{LatLon}.

       @raise ValueError: Invalid B{C{point}}, lat- or longitude.

       @see: L{sphericalNvector.LatLon.isenclosedBy},
             L{sphericalTrigonometry.LatLon.isenclosedBy} and
             U{MultiDop GeogContainPt<https://GitHub.com/NASA/MultiDop>}
             (U{Shapiro et al. 2009, JTECH
             <https://Journals.AMetSoc.org/doi/abs/10.1175/2009JTECHA1256.1>}
             and U{Potvin et al. 2012, JTECH
             <https://Journals.AMetSoc.org/doi/abs/10.1175/JTECH-D-11-00019.1>}).
    '''
    try:
        y0, x0 = point.lat, point.lon
    except AttributeError:
        try:
            y0, x0 = map1(float, *point[:2])
        except (IndexError, TypeError, ValueError) as x:
            raise _ValueError(point=point, txt=str(x))

    if wrap:
        x0, y0 = wrap180(x0), wrap90(y0)

        def _dxy3(x1, x2, y2, w):
            dx, x2 = unroll180(x1, x2, wrap=w)
            return dx, x2, y2

    else:
        x0 = fmod(x0, _360_0)  # not x0 % 360
        x0_180_ = x0 - _180_0
        x0_180  = x0 + _180_0

        def _dxy3(x1, x, y, unused):  # PYCHOK expected
            x = float(x) % _360_0
            if x < x0_180_:
                x += _360_0
            elif x >= x0_180:
                x -= _360_0
            return (x - x1), x, y

    Ps = LatLon2PsxyIter(points, wrap=wrap, loop=1)
    p  = Ps[0]
    e  = m = False
    S  = Fsum()

    _, x1, y1 = _dxy3(x0, p.x, p.y, False)
    for i, p in Ps.enumerate(closed=True):
        dx, x2, y2 = _dxy3(x1, p.x, p.y, (wrap if i else False))
        # ignore duplicate and near-duplicate pts
        if max(abs(dx), abs(y2 - y1)) > EPS:
            # determine if polygon edge (x1, y1)..(x2, y2) straddles
            # point (lat, lon) or is on boundary, but do not count
            # edges on boundary as more than one crossing
            if abs(dx) < 180 and (x1 < x0 <= x2 or x2 < x0 <= x1):
                m = not m
                dy = (x0 - x1) * (y2 - y1) - (y0 - y1) * dx
                if (dy > 0 and dx >= 0) or (dy < 0 and dx <= 0):
                    e = not e

            S += sin(radians(y2))
            x1, y1 = x2, y2

    # An odd number of meridian crossings means, the polygon
    # contains a pole.  Assume it is the pole on the hemisphere
    # containing the polygon mean point and if the polygon does
    # contain the North Pole, flip the result.
    if m and S.fsum() > 0:
        e = not e
    return e


def ispolar(points, wrap=False):
    '''Check whether a polygon encloses a pole.

       @arg points: The polygon points (C{LatLon}[]).
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).

       @return: C{True} if the polygon encloses a pole, C{False}
                otherwise.

       @raise PointsError: Insufficient number of B{C{points}}

       @raise TypeError: Some B{C{points}} are not C{LatLon} or don't
                         have C{bearingTo2}, C{initialBearingTo}
                         and C{finalBearingTo} methods.
    '''
    def _cds(points, wrap):  # iterate over course deltas
        Ps = PointsIter(points, loop=2)
        p2, p1 = Ps[0:2]
        b1, _ = _bearingTo2(p2, p1, wrap=wrap)

        for p2 in Ps.iterate(closed=True):
            if not p2.isequalTo(p1, EPS):
                b, b2 = _bearingTo2(p1, p2, wrap=wrap)
                yield wrap180(b - b1)  # (b - b1 + 540) % 360 - 180
                yield wrap180(b2 - b)  # (b2 - b + 540) % 360 - 180
                p1, b1 = p2, b2

    # summation of course deltas around pole is 0° rather than normally ±360°
    # <https://blog.Element84.com/determining-if-a-spherical-polygon-contains-a-pole.html>
    s = fsum(_cds(points, wrap))
    # XXX fix (intermittant) edge crossing pole - eg (85,90), (85,0), (85,-90)
    return abs(s) < 90  # "zero-ish"


def luneOf(lon1, lon2, closed=False, LatLon=LatLon_, **LatLon_kwds):
    '''Generate an ellipsoidal or spherical U{lune
       <https://WikiPedia.org/wiki/Spherical_lune>}-shaped path or polygon.

       @arg lon1: Left longitude (C{degrees90}).
       @arg lon2: Right longitude (C{degrees90}).
       @kwarg closed: Optionally, close the path (C{bool}).
       @kwarg LatLon: Class to use (L{LatLon_}).
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}}
                           keyword arguments.

       @return: A tuple of 4 or 5 B{C{LatLon}} instances outlining
                the lune shape.

       @see: U{Latitude-longitude quadrangle
             <https://www.MathWorks.com/help/map/ref/areaquad.html>}.
    '''
    t = (LatLon(  _0_0, lon1, **LatLon_kwds),
         LatLon( _90_0, lon1, **LatLon_kwds),
         LatLon(  _0_0, lon2, **LatLon_kwds),
         LatLon(-_90_0, lon2, **LatLon_kwds))
    if closed:
        t += t[:1]
    return t


def nearestOn5(point, points, closed=False, wrap=False, LatLon=None, **options):
    '''Locate the point on a path or polygon closest to an other point.

       If the given point is within the extent of a polygon edge,
       the closest point is on that edge, otherwise the closest
       point is the nearest of that edge's end points.

       Distances are approximated by function L{equirectangular_},
       subject to the supplied B{C{options}}.

       @arg point: The other, reference point (C{LatLon}).
       @arg points: The path or polygon points (C{LatLon}[]).
       @kwarg closed: Optionally, close the path or polygon (C{bool}).
       @kwarg wrap: Wrap and L{unroll180} longitudes and longitudinal
                    delta (C{bool}) in function L{equirectangular_}.
       @kwarg LatLon: Optional class to return the closest point (C{LatLon})
                      or C{None}.
       @kwarg options: Other keyword arguments for function L{equirectangular_}.

       @return: A L{NearestOn3Tuple}C{(closest, distance, angle)} with the
                {closest} point (B{C{LatLon}}) or if C{B{LatLon} is None},
                a L{NearestOn5Tuple}C{(lat, lon, distance, angle, height)}.
                The C{distance} is the L{equirectangular_} distance between
                the C{closest} and reference B{C{point}} in C{degrees}.  The
                C{angle} from the reference B{C{point}} to the C{closest} is
                in compass C{degrees360}, like function L{compassAngle}.

       @raise LimitError: Lat- and/or longitudinal delta exceeds the
                          B{C{limit}}, see function L{equirectangular_}.

       @raise PointsError: Insufficient number of B{C{points}}

       @raise TypeError: Some B{C{points}} are not C{LatLon}.

       @see: Function L{degrees2m} to convert C{degrees} to C{meter}.
    '''
    def _d2yx(p2, p1, u, w):
        # w = wrap if (not closed or w < (n - 1)) else False
        # equirectangular_ returns a Distance4Tuple(distance
        # in degrees squared, delta lat, delta lon, p2.lon
        # unroll/wrap); the previous p2.lon unroll/wrap
        # is also applied to the next edge's p1.lon
        return equirectangular_(p1.lat, p1.lon + u,
                                p2.lat, p2.lon, wrap=w, **options)

    def _h(p):  # get height or default 0
        return getattr(p, _height_, 0) or 0

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
    #   f = (y * dy + x * dx) / hypot2(dx, dy)

    Ps = PointsIter(points, loop=1)
    p1 = c = Ps[0]
    u1 = u = _0_0
    d, dy, dx, _ = _d2yx(p1, point, u1, False)
    for i, p2 in Ps.enumerate(closed=closed):
        # iff wrapped, unroll lon1 (actually previous
        # lon2) like function unroll180/-PI would've
        w = False if closed and i == 0 else wrap
        d21, y21, x21, u2 = _d2yx(p2, p1, u1, w)
        if d21 > EPS:
            # distance point to p1, y01 and x01 inverted
            d2, y01, x01, _ = _d2yx(point, p1, u1, closed)
            if d2 > EPS:
                w2 = y01 * y21 + x01 * x21
                if w2 > 0:
                    if w2 < d21:
                        # closest is between p1 and p2, use
                        # original delta's, not y21 and x21
                        f = w2 / d21
                        p1 = LatLon_(favg(p1.lat, p2.lat, f=f),
                                     favg(p1.lon, p2.lon + u2, f=f),
                              height=favg(_h(p1), _h(p2), f=f))
                        u1 = _0_0
                    else:  # p2 is closest
                        p1, u1 = p2, u2
                    d2, y01, x01, _ = _d2yx(point, p1, u1, closed)
            if d2 < d:  # p1 is closer, y01 and x01 inverted
                c, u, d, dy, dx = p1, u1, d2, -y01, -x01
        p1, u1 = p2, u2

    d, a, h = hypot(dx, dy), atan2b(dx, dy), _h(c)
    if LatLon is None:
        r = NearestOn5Tuple(c.lat, c.lon + u, d, a, h)
    else:
        r = LatLon(c.lat, c.lon + u, height=h)
        r = NearestOn3Tuple(r, d, a)
    return _xnamed(r, nameof(point))


def perimeterOf(points, closed=False, adjust=True, radius=R_M, wrap=True):
    '''Approximate the perimeter of a path or polygon.

       @arg points: The path or polygon points (C{LatLon}[]).
       @kwarg closed: Optionally, close the path or polygon (C{bool}).
       @kwarg adjust: Adjust the wrapped, unrolled longitudinal delta
                      by the cosine of the mean latitude (C{bool}).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg wrap: Wrap lat-, wrap and unroll longitudes (C{bool}).

       @return: Approximate perimeter (C{meter}, same units as
                B{C{radius}}).

       @raise PointsError: Insufficient number of B{C{points}}

       @raise TypeError: Some B{C{points}} are not C{LatLon}.

       @raise ValueError: Invalid B{C{radius}}.

       @note: This perimeter is based on the L{equirectangular_}
              distance approximation and is ill-suited for regions
              exceeding several hundred Km or Miles or with
              near-polar latitudes.

       @see: L{sphericalTrigonometry.perimeterOf} and
             L{ellipsoidalKarney.perimeterOf}.
    '''
    def _degs(points, closed, wrap):  # angular edge lengths in degrees
        Ps = LatLon2PsxyIter(points, wrap=wrap, loop=1)
        p1, u = Ps[0], _0_0  # previous x2's unroll/wrap
        for i, p2 in Ps.enumerate(closed=closed):
            w = False if closed and i == 0 else wrap
            # apply previous x2's unroll/wrap to new x1
            _, dy, dx, u = equirectangular_(p1.y, p1.x + u, p2.y, p2.x,
                                            adjust=adjust,
                                            limit=None,
                                            wrap=w)  # PYCHOK non-sequence
            yield hypot(dx, dy)
            p1 = p2

    d = fsum(_degs(points, closed, wrap))
    return degrees2m(d, radius=radius)


def quadOf(latS, lonW, latN, lonE, closed=False, LatLon=LatLon_, **LatLon_kwds):
    '''Generate a quadrilateral path or polygon from two points.

       @arg latS: Southernmost latitude (C{degrees90}).
       @arg lonW: Westernmost longitude (C{degrees180}).
       @arg latN: Northernmost latitude (C{degrees90}).
       @arg lonE: Easternmost longitude (C{degrees180}).
       @kwarg closed: Optionally, close the path (C{bool}).
       @kwarg LatLon: Class to use (L{LatLon_}).
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}}
                           keyword arguments.

       @return: Return a tuple of 4 or 5 B{C{LatLon}} instances
                outlining the quadrilateral.

       @see: Function L{boundsOf}.
    '''
    t = (LatLon(latS, lonW, **LatLon_kwds),
         LatLon(latN, lonW, **LatLon_kwds),
         LatLon(latN, lonE, **LatLon_kwds),
         LatLon(latS, lonE, **LatLon_kwds))
    if closed:
        t += t[:1]
    return t

# **) MIT License
#
# Copyright (C) 2016-2021 -- mrJean1 at Gmail -- All Rights Reserved.
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
