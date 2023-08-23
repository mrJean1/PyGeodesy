
# -*- coding: utf-8 -*-

u'''Utilities for point lists, tuples, etc.

Functions to handle collections and sequences of C{LatLon} points
specified as 2-d U{NumPy<https://www.NumPy.org>}, C{arrays} or tuples
as C{LatLon} or as C{pseudo-x/-y} pairs.

C{NumPy} arrays are assumed to contain rows of points with a lat-, a
longitude -and possibly other- values in different columns.  While
iterating over the array rows, create an instance of a given C{LatLon}
class "on-the-fly" for each row with the row's lat- and longitude.

The original C{NumPy} array is read-accessed only and never duplicated,
except to return a I{subset} of the original array.

For example, to process a C{NumPy} array, wrap the array by instantiating
class L{Numpy2LatLon} and specifying the column index for the lat- and
longitude in each row.  Then, pass the L{Numpy2LatLon} instance to any
L{pygeodesy} function or method accepting a I{points} argument.

Similarly, class L{Tuple2LatLon} is used to instantiate a C{LatLon} from
each 2+tuple in a sequence of such 2+tuples using the C{ilat} lat- and
C{ilon} longitude index in each 2+tuple.
'''

from pygeodesy.basics import isclass, isint, isscalar, issequence, \
                             issubclassof, _Sequence, _xcopy, _xdup, \
                             _xinstanceof
from pygeodesy.constants import EPS, EPS1, PI_2, R_M, isnear0, isnear1, \
                               _umod_360, _0_0, _0_5, _1_0, _2_0, _6_0, \
                               _90_0, _N_90_0, _180_0, _360_0
# from pygeodesy.datums import _spherical_datum  # from .formy
from pygeodesy.dms import F_D, latDMS, lonDMS, parseDMS, S_DEG, S_DMS, \
                          S_MIN, S_SEC, S_SEP
from pygeodesy.errors import CrossError, crosserrors, _IndexError, \
                            _IsnotError, _TypeError, _ValueError, \
                            _xattr, _xkwds, _xkwds_pop
from pygeodesy.fmath import favg, fdot, hypot,  Fsum, fsum
# from pygeodesy.fsums import Fsum, fsum  # from .fmath
from pygeodesy.formy import _bearingTo2, equirectangular_, _isequalTo, \
                             isnormal, latlon2n_xyz, normal, _spherical_datum
from pygeodesy.interns import NN, _colinear_, _COMMASPACE_, _composite_, \
                             _DEQUALSPACED_, _ELLIPSIS_, _EW_, _immutable_, \
                             _lat_, _lon_, _LatLon_, _near_, _no_, _not_, \
                             _NS_, _point_, _SPACE_, _UNDER_, _valid_
from pygeodesy.iters import LatLon2PsxyIter, PointsIter, points2
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import classname, nameof, notImplemented, notOverloaded, \
                           _NamedTuple, _xnamed, _xother3, _xotherError
from pygeodesy.namedTuples import Bounds2Tuple, Bounds4Tuple, \
                                  LatLon2Tuple, NearestOn3Tuple, \
                                  NearestOn5Tuple, PhiLam2Tuple, \
                                  Point3Tuple, Vector3Tuple, Vector4Tuple
from pygeodesy.nvectorBase import NvectorBase, _N_vector_
from pygeodesy.props import deprecated_method, Property, Property_RO, \
                            property_doc_, property_RO, _update_all
from pygeodesy.streprs import Fmt, hstr, instr, pairs
from pygeodesy.units import Number_, Radius, Scalar, Scalar_
from pygeodesy.utily import atan2b, degrees90, degrees180, degrees2m, \
                            unroll180, _unrollon, unrollPI, _Wrap, wrap180

from math import cos, fabs, fmod, radians, sin

__all__ = _ALL_LAZY.points
__version__ = '23.08.23'

_ilat_  = 'ilat'
_ilon_  = 'ilon'
_ncols_ = 'ncols'
_nrows_ = 'nrows'


class LatLon_(object):  # XXX in heights._HeightBase.height
    '''Low-overhead C{LatLon} class for L{Numpy2LatLon} and L{Tuple2LatLon}.
    '''
    # __slots__ efficiency is voided if the __slots__ class attribute
    # is used in a subclass of a class with the traditional __dict__,
    # see <https://docs.Python.org/2/reference/datamodel.html#slots>
    # and __slots__ must be repeated in sub-classes, see Luciano
    # Ramalho, "Fluent Python", O'Reilly, 2016 p. 276+ "Problems
    # with __slots__" (also at <https://Books.Google.ie/books?
    #   id=bIZHCgAAQBAJ&lpg=PP1&dq=fluent%20python&pg=PT364#
    #   v=onepage&q=“Problems%20with%20__slots__”&f=false>),
    #   2022 p. 390 "Summarizing the Issues with __slots__".
    #
    # __slots__ = (_lat_, _lon_, _height_, _datum_, _name_)
    # Property_RO = property_RO  # no __dict__ with __slots__!
    #
    # However, sys.getsizeof(LatLon_(1, 2)) is 72-88 with __slots__
    # but only 48-56 bytes without in Python 2.7.18+ and Python 3+.

    def __init__(self, latlonh, lon=None, name=NN, height=0,
                                          datum=None, wrap=False):
        '''Creat a new, mininal, low-overhead L{LatLon_} instance.

           @arg latlonh: Latitude (C{degrees} or DMS C{str} with N or S suffix) or
                         a previous C{LatLon} instance provided C{B{lon}=None}.
           @kwarg lon: Longitude (C{degrees} or DMS C{str} with E or W suffix) or
                       C(None), indicating B{C{latlonh}} is a C{LatLon}.
           @kwarg name: Optional name (C{str}).
           @kwarg height: Optional height (C{meter}, conventionally).
           @kwarg datum: Optional datum (L{Datum}, L{Ellipsoid}, L{Ellipsoid2},
                         L{a_f2Tuple} or I{scalar} radius) or C{None}.
           @kwarg wrap: If C{True}, wrap or I{normalize} B{C{lat}} and B{C{lon}}
                        (C{bool}).

           @raise TypeError: Invalid B{C{datum}} or B{C{latlonh}} not a C{LatLon}.

           @note: The lat- and longitude are taken as-given,
                  un-clipped and un-validated .
        '''
        if lon is None:  # PYCHOK no cover
            try:
                ll = latlonh.lat, latlonh.lon
                height = _xattr(latlonh, height=height)
            except AttributeError:
                raise _IsnotError(_LatLon_, latlonh=latlonh)
            if wrap:
                ll = _Wrap.latlon(*ll)
        elif wrap:  # PYCHOK no cover
            ll = _Wrap.latlonDMS2(latlonh, lon)
        else:  # must be latNS, lonEW
            ll = parseDMS(latlonh, suffix=_NS_), parseDMS(lon, suffix=_EW_)

        self.lat, self.lon = ll
        self.name   = str(name) if name else NN
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

           @return: The copy (C{This} (sub-)class).
        '''
        return _xcopy(self, deep=deep)

    def dup(self, **items):
        '''Duplicate this instance, replacing some items.

           @kwarg items: Attributes to be changed (C{any}).

           @return: The duplicate (C{This} (sub-)class).

           @raise AttributeError: Some B{C{items}} invalid.
        '''
        d = _xdup(self, **items)
        if items:
            _update_all(d)
        return d

    def heightStr(self, prec=-2):
        '''Return a string for the height B{C{height}}.

           @kwarg prec: Number of (decimal) digits, unstripped (C{int}).

           @see: Function L{pygeodesy.hstr}.
        '''
        return hstr(self.height, prec=prec)

    def intermediateTo(self, other, fraction, height=None, wrap=False):
        '''Locate the point at a given fraction between (or along) this
           and an other point.

           @arg other: The other point (C{LatLon}).
           @arg fraction: Fraction between both points (C{float},
                          0.0 for this and 1.0 for the other point).
           @kwarg height: Optional height (C{meter}), overriding the
                          intermediate height.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                        the B{C{other}} point (C{bool}).

           @return: Intermediate point (this C{LatLon}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        f = Scalar(fraction=fraction)
        if isnear0(f):
            r = self
        elif isnear1(f) and not wrap:
            r = self.others(other)
        else:
            p = self.others(other)
            h = favg(self.height, p.height, f=f) if height is None else height
            _, lat, lon = _Wrap.latlon3(self.lon, p.lat, p.lon, wrap=wrap)
            r = self.classof(favg(self.lat, lat, f=f),
                             favg(self.lon, lon, f=f),
                             height=h, datum=self.datum,
                             name=self.intermediateTo.__name__)
        return r

    @Property_RO  # PYCHOK no cover
    def isEllipsoidal(self):
        '''Check whether this point is ellipsoidal (C{bool} or C{None} if unknown).
        '''
        return self.datum.isEllipsoidal if self.datum else None

    @Property_RO  # PYCHOK no cover
    def isEllipsoidalLatLon(self):
        '''Get C{LatLon} base.
        '''
        return False

    def isequalTo(self, other, eps=None):
        '''Compare this point with an other point, I{ignoring} height.

           @arg other: The other point (C{LatLon}).
           @kwarg eps: Tolerance for equality (C{degrees}).

           @return: C{True} if both points are identical,
                    I{ignoring} height, C{False} otherwise.

           @raise UnitError: Invalid B{C{eps}}.
        '''
        return _isequalTo(self, self.others(other), eps=eps)

    @Property_RO
    def isnormal(self):
        '''Return C{True} if this point is normal (C{bool}),
           meaning C{abs(lat) <= 90} and C{abs(lon) <= 180}.

           @see: Methods L{normal}, L{toNormal} and functions
                 L{pygeodesy.isnormal} and L{pygeodesy.normal}.
        '''
        return isnormal(self.lat, self.lon, eps=0)

    @Property_RO
    def isSpherical(self):  # PYCHOK no cover
        '''Check whether this point is spherical (C{bool} or C{None} if unknown).
        '''
        return self.datum.isSpherical if self.datum else None

    @Property_RO
    def lam(self):
        '''Get the longitude (B{C{radians}}).
        '''
        return radians(self.lon)

    @Property
    def latlon(self):
        '''Get the lat- and longitude in C{degrees} (L{LatLon2Tuple}C{(lat, lon)}).
        '''
        return LatLon2Tuple(self.lat, self.lon, name=self.name)

    @latlon.setter  # PYCHOK setter!
    def latlon(self, latlon):
        '''Set the lat- and longitude.

           @arg latlon: New lat- and longitude in C{degrees} (C{2-tuple} or {-list}).
        '''
        if self.latlon != latlon[:2]:
            _update_all(self)
            self.lat, self.lon = latlon[:2]

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

    def normal(self):
        '''Normalize this point I{in-place} to C{abs(lat) <= 90} and
           C{abs(lon) <= 180}.

           @return: C{True} if this point was I{normal}, C{False} if it
                    wasn't (but is now).

           @see: Property L{isnormal} and method L{toNormal}.
        '''
        n = self.isnormal
        if not n:
            self.latlon = normal(*self.latlon)
        return n

    @Property_RO
    def phi(self):
        '''Get the latitude (B{C{radians}}).
        '''
        return radians(self.lat)

    @Property_RO
    def philam(self):
        '''Get the lat- and longitude (L{PhiLam2Tuple}C{(phi, lam)}).
        '''
        return PhiLam2Tuple(self.phi, self.lam, name=self.name)

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
        return PointsIter(points, base=self, loop=loop, dedup=dedup)

    @deprecated_method
    def to2ab(self):  # PYCHOK no cover
        '''DEPRECATED, use property L{philam}.'''
        return self.philam

    def toNormal(self, deep=False, name=NN):
        '''Get a copy of this point normalized to C{abs(lat) <= 90}
           and C{abs(lon) <= 180}.

           @kwarg deep: If C{True} make a deep, otherwise a
                        shallow copy (C{bool}).
           @kwarg name: Optional name of the copy (C{str}).

           @return: A copy of this point, I{normalized} and
                    optionally renamed (C{LatLon}).

           @see: Property L{isnormal}, method L{normal} and
                 function L{pygeodesy.normal}.
        '''
        ll = self.copy(deep=deep)
        _  = ll.normal()
        if name:
            ll.name = str(name)
        return ll

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

           @raise TypeError: Invalid B{C{Nvector}} or B{C{Nvector_kwds}}
                             argument.
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

    def toStr(self, form=F_D, joined=_COMMASPACE_, **prec_sep_s_D_M_S_kwds):
        '''Convert this point to a "lat, lon[, height][, name][, ...]" string,
           formatted in the given C{B{form}at}.

           @kwarg form: The lat-/longitude C{B{form}at} to use (C{str}), see
                        functions L{pygeodesy.latDMS} or L{pygeodesy.lonDMS}.
           @kwarg joined: Separator to join the lat-, longitude, heigth, name
                          and other strings (C{str} or C{None} or C{NN} for
                          non-joined).
           @kwarg prec_sep_s_D_M_S_kwds: Optional C{B{prec}ision}, C{B{sep}arator},
                      B{C{s_D}}, B{C{s_M}}, B{C{s_S}}, B{C{s_DMS}} and possibly
                      other keyword arguments, see functions L{pygeodesy.latDMS}
                      or L{pygeodesy.lonDMS}.

           @return: This point in the specified C{B{form}at}, etc. (C{str} or
                    a 2- or 3+tuple C{(lat_str, lon_str[, height_str][, name_str][,
                    ...])} if C{B{joined}=NN} or C{B{joined}=None} and with the
                    C{height_str} and C{name_str} only included if non-zero
                    respectively non-empty).

           @see: Function L{pygeodesy.latDMS} or L{pygeodesy.lonDMS} for more
                 details about keyword arguments C{B{form}at}, C{B{prec}ision},
                 C{B{sep}arator}, B{C{s_D}}, B{C{s_M}}, B{C{s_S}} and B{C{s_DMS}}.
        '''
        def _t3(prec=None, sep=S_SEP, s_D=S_DEG, s_M=S_MIN, s_S=S_SEC, s_DMS=S_DMS, **kwds):
            return dict(prec=prec, sep=sep, s_D=s_D, s_M=s_M, s_S=s_S, s_DMS=s_DMS), kwds, prec

        prec_sep_s_D_M_S, kwds, prec = _t3(**prec_sep_s_D_M_S_kwds)
        t = (latDMS(self.lat, form=form, **prec_sep_s_D_M_S),
             lonDMS(self.lon, form=form, **prec_sep_s_D_M_S))
        if self.height:
            t += (self.heightStr(),)
        if self.name:
            t += (repr(self.name),)
        if kwds:
            t += pairs(kwds) if prec is None else pairs(kwds, prec=prec)
        return joined.join(t) if joined else t

    @deprecated_method
    def toStr2(self, **kwds):  # PYCHOK no cover
        '''DEPRECATED, used method L{toRepr}.'''
        return self.toRepr(**kwds)


def _isLatLon(inst):
    '''(INTERNAL) Check a C{LatLon} or C{LatLon_} instance.
    '''
    return isinstance(inst, (LatLon_, _MODS.latlonBase.LatLonBase))


def _isLatLon_(LL):
    '''(INTERNAL) Check a (sub-)class of C{LatLon_}.
    '''
    return issubclassof(LL, LatLon_) or (isclass(LL) and
            all(hasattr(LL, a) for a in _ALL_ATTRS_))


# get all pseudo-slots for class C{LatLon_}
_ALL_ATTRS_ = tuple(LatLon_(0, 0).__dict__.keys())


class _Basequence(_Sequence):  # immutable, on purpose
    '''(INTERNAL) Base class.
    '''
    _array    =  []
    _epsilon  =  EPS
    _itemname = _point_

    def _contains(self, point):
        '''(INTERNAL) Check for a matching point.
        '''
        return any(self._findall(point, ()))

    def copy(self, deep=False):  # PYCHOK no cover
        '''Make a shallow or deep copy of this instance.

           @kwarg deep: If C{True} make a deep, otherwise a
                        shallow copy (C{bool}).

           @return: The copy (C{This class} or subclass thereof).
        '''
        return _xcopy(self, deep=deep)

    def _count(self, point):
        '''(INTERNAL) Count the number of matching points.
        '''
        return sum(1 for _ in self._findall(point, ()))  # NOT len()!

    def dup(self, **items):  # PYCHOK no cover
        '''Duplicate this instance, I{without replacing items}.

           @kwarg items: No attributes (I{not allowed}).

           @return: The duplicate (C{This} (sub-)class).

           @raise TypeError: Any B{C{items}} invalid.
        '''
        if items:
            t = _SPACE_(classname(self), _immutable_)
            raise _TypeError(txt=t, this=self, **items)
        return _xdup(self)

    @property_doc_(''' the equality tolerance (C{float}).''')
    def epsilon(self):
        '''Get the tolerance for equality tests (C{float}).
        '''
        return self._epsilon

    @epsilon.setter  # PYCHOK setter!
    def epsilon(self, tol):
        '''Set the tolerance for equality tests (C{scalar}).

           @raise UnitError: Non-scalar or invalid B{C{tol}}.
        '''
        self._epsilon = Scalar_(tolerance=tol)

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
        # Luciano Ramalho, "Fluent Python", O'Reilly, 2016 p. 290+, 2022 p. 405+
        if isinstance(index, slice):
            # XXX an numpy.[nd]array slice is a view, not a copy
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
        _array, _point = self._array, self.point
        for i in range(len(self)):
            yield _point(_array[i])

    def point(self, *attrs):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}, see function C{notOverloaded in}.

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
        _array, point = self._array, self.point
        for i in range(len(self) - 1, -1, -1):
            yield point(_array[i])

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


class _Array2LatLon(_Basequence):  # immutable, on purpose
    '''(INTERNAL) Base class for Numpy2LatLon or Tuple2LatLon.
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
        self._shape = Shape2Tuple(shape)  # *shape

        if LatLon:  # check the point class
            if not _isLatLon_(LatLon):
                raise _IsnotError(_valid_, LatLon=LatLon)
            self._LatLon = LatLon

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

        _ilat,  _ilon = self._ilat,  self._ilon
        _array, _eps  = self._array, self._epsilon
        for i in self._range(*start_end):
            row = _array[i]
            if fabs(row[_ilat] - lat) <= _eps and \
               fabs(row[_ilon] - lon) <= _eps:
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
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                        the B{C{latlons}} points (C{bool}).

           @raise PointsError: Insufficient number of B{C{latlons}}.

           @raise TypeError: Some B{C{points}} are not B{C{base}}.
        '''
        self._closed = closed
        self._len, self._array = points2(latlons, closed=closed)
        if radius:
            self._radius = r = Radius(radius)
            self._deg2m  = degrees2m(_1_0, r)
        if not wrap:
            self._wrap = False

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

            def _x_y_ll3(ll):  # match LatLon
                return ll.lon, ll.lat, ll

        except AttributeError:
            try:
                x, y = xy[:2]
            except (IndexError, TypeError, ValueError):
                raise _IsnotError(_valid_, xy=xy)

            _x_y_ll3 = self.point  # PYCHOK expected

        _array, _eps = self._array, self._epsilon
        for i in self._range(*start_end):
            xi, yi, _ = _x_y_ll3(_array[i])
            if fabs(xi - x) <= _eps and \
               fabs(yi - y) <= _eps:
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
            y, x = _Wrap.latlon(y, x)
        d = self._deg2m
        if d:  # convert degrees to meter (or radians)
            x *= d
            y *= d
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
            c = cos(h) if fabs(h) < PI_2 else _0_0
            return w * h * (c + 1.2876) * _0_5
    else:
        def _adjust(w, h):  # PYCHOK expected
            return w * h

    # setting radius=1 converts degrees to radians
    Ps = LatLon2PsxyIter(points, loop=1, radius=_1_0, wrap=wrap)
    x1, y1, ll = Ps[0]
    pts = [ll]  # for _areaError

    A2 = Fsum()  # trapezoidal area in radians**2
    for p in Ps.iterate(closed=True):
        x2, y2, ll = p
        if len(pts) < 4:
            pts.append(ll)
        w, x2 = unrollPI(x1, x2, wrap=wrap and not Ps.looped)
        A2 += _adjust(w, (y2 + y1) * _0_5)
        x1, y1 = x2, y2

    return A2.fsum(), tuple(pts)


def _areaError(pts, near_=NN):  # imported by .ellipsoidalKarney
    '''(INTERNAL) Area issue.
    '''
    t = _ELLIPSIS_(pts[:3], NN)
    return _ValueError(NN(near_, 'zero or polar area'), txt=t)


def areaOf(points, adjust=True, radius=R_M, wrap=True):
    '''Approximate the area of a polygon or composite.

       @arg points: The polygon points or clips (C{LatLon}[],
                    L{BooleanFHP} or L{BooleanGH}).
       @kwarg adjust: Adjust the wrapped, unrolled longitudinal delta
                      by the cosine of the mean latitude (C{bool}).
       @kwarg radius: Mean earth radius (C{meter}) or C{None}.
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                    the B{C{points}} (C{bool}).

       @return: Approximate area (I{square} C{meter}, same units as
                B{C{radius}} or C{radians} I{squared} if B{C{radius}}
                is C{None}).

       @raise PointsError: Insufficient number of B{C{points}}

       @raise TypeError: Some B{C{points}} are not C{LatLon}.

       @raise ValueError: Invalid B{C{radius}}.

       @note: This area approximation has limited accuracy and is
              ill-suited for regions exceeding several hundred Km
              or Miles or with near-polar latitudes.

       @see: L{sphericalNvector.areaOf}, L{sphericalTrigonometry.areaOf},
             L{ellipsoidalExact.areaOf} and L{ellipsoidalKarney.areaOf}.
    '''
    if _MODS.booleans.isBoolean(points):
        a = points._sum1(areaOf, adjust=adjust, radius=None, wrap=wrap)
    else:
        a, _ = _area2(points, adjust, wrap)
    return fabs(a if radius is None else (Radius(radius)**2 * a))


def boundsOf(points, wrap=False, LatLon=None):  # was=True
    '''Determine the bottom-left SW and top-right NE corners of a
       path or polygon.

       @arg points: The path or polygon points (C{LatLon}[]).
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                    the B{C{points}} (C{bool}).
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
    Ps = LatLon2PsxyIter(points, loop=1, wrap=wrap)
    w, s, _ = e, n, _ = Ps[0]

    v = w
    for x, y, _ in Ps.iterate(closed=False):  # [1:]
        if wrap:
            _, x = unroll180(v, x, wrap=True)
        v = x

        if w > x:
            w = x
        elif e < x:
            e = x

        if s > y:
            s = y
        elif n < y:
            n = y

    return Bounds4Tuple(s, w, n, e) if LatLon is None else \
           Bounds2Tuple(LatLon(s, w), LatLon(n, e))  # PYCHOK inconsistent


def _distanceTo(Error, **name_points):  # .frechet, .hausdorff, .heights
    '''(INTERNAL) Chech callable C{distanceTo} methods.
    '''
    name, ps = name_points.popitem()
    for i, p in enumerate(ps):
        if not callable(_xattr(p, distanceTo=None)):
            n = _distanceTo.__name__[1:]
            t = _SPACE_(_no_, callable.__name__, n)
            raise Error(Fmt.SQUARE(name, i), p, txt=t)
    return ps


def centroidOf(points, wrap=False, LatLon=None):  # was=True
    '''Determine the centroid of a polygon.

       @arg points: The polygon points (C{LatLon}[]).
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                    B{C{points}} (C{bool}).
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
    Ps = LatLon2PsxyIter(points, loop=1, radius=_1_0, wrap=wrap)
    x1, y1, ll = Ps[0]
    pts = [ll]  # for _areaError
    for p in Ps.iterate(closed=True):
        x2, y2, ll = p
        if len(pts) < 4:
            pts.append(ll)
        if wrap and not Ps.looped:
            _, x2 = unrollPI(x1, x2, wrap=True)
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

    a = A.fmul(_6_0).fover(_2_0)
    if isnear0(a):
        raise _areaError(pts, near_=_near_)
    y, x = degrees90(Y.fover(a)), degrees180(X.fover(a))
    return LatLon2Tuple(y, x) if LatLon is None else LatLon(y, x)


def fractional(points, fi, j=None, wrap=None, LatLon=None, Vector=None, **kwds):
    '''Return the point at a given I{fractional} index.

       @arg points: The points (C{LatLon}[], L{Numpy2LatLon}[],
                    L{Tuple2LatLon}[], C{Cartesian}[], C{Vector3d}[],
                    L{Vector3Tuple}[]).
       @arg fi: The fractional index (L{FIx}, C{float} or C{int}).
       @kwarg j: Optionally, index of the other point (C{int}).
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                    B{{points}} (C{bool}) or C{None} for a backward
                    compatible L{LatLon2Tuple} or B{C{LatLon}} with
                    averaged lat- and longitudes.  Use C{True} or
                    C{False} to get the I{fractional} point computed
                    by method C{B{points}[fi].intermediateTo}.
       @kwarg LatLon: Optional class to return the I{intermediate},
                      I{fractional} point (C{LatLon}) or C{None}.
       @kwarg Vector: Optional class to return the I{intermediate},
                      I{fractional} point (C{Cartesian}, C{Vector3d})
                      or C{None}.
       @kwarg kwds: Optional, additional B{C{LatLon}} I{or} B{C{Vector}}
                    keyword arguments, ignored if both C{B{LatLon}} and
                    C{B{Vector}} are C{None}.

       @return: A L{LatLon2Tuple}C{(lat, lon)} if B{C{wrap}}, B{C{LatLon}}
                and B{C{Vector}} all are C{None}, the defaults.

                An instance of B{C{LatLon}} if not C{None} I{or} an instance
                of B{C{Vector}} if not C{None}.

                Otherwise with B{C{wrap}} either C{True} or C{False} and
                B{C{LatLon}} and B{C{Vector}} both C{None}, an instance of
                B{C{points}}' (sub-)class C{intermediateTo} I{fractional}.

                Summarized as follows:

                  >>>  wrap  | LatLon | Vector | returned type/value
                  #   -------+--------+--------+--------------+------
                  #          |        |        | LatLon2Tuple | favg
                  #    None  |  None  |  None  |   or**       |
                  #          |        |        | Vector3Tuple | favg
                  #    None  | LatLon |  None  | LatLon       | favg
                  #    None  |  None  | Vector | Vector       | favg
                  #   -------+--------+--------+--------------+------
                  #    True  |  None  |  None  | points'      | .iTo
                  #    True  | LatLon |  None  | LatLon       | .iTo
                  #    True  |  None  | Vector | Vector       | .iTo
                  #   -------+--------+--------+--------------+------
                  #    False |  None  |  None  | points'      | .iTo
                  #    False | LatLon |  None  | LatLon       | .iTo
                  #    False |  None  | Vector | Vector       | .iTo
                  # _____
                  # favg) averaged lat, lon or x, y, z values
                  # .iTo) value from points[fi].intermediateTo
                  # **) depends on base class of points[fi]

       @raise IndexError: Fractional index B{C{fi}} invalid or B{C{points}}
                          not subscriptable or not closed.

       @raise TypeError: Invalid B{C{LatLon}}, B{C{Vector}} or B{C{kwds}}
                         argument.

       @see: Class L{FIx} and method L{FIx.fractional}.
    '''
    if LatLon and Vector:  # PYCHOK no cover
        kwds = _xkwds(kwds, fi=fi, LatLon=LatLon, Vector=Vector)
        raise _TypeError(txt=fractional.__name__, **kwds)
    w = wrap if LatLon else False  # intermediateTo
    try:
        if not isscalar(fi) or fi < 0:
            raise IndexError
        n = _xattr(fi, fin=0)
        p = _fractional(points, fi, j, fin=n, wrap=w)  # see .units.FIx
        if LatLon:
            p = LatLon(p.lat, p.lon, **kwds)
        elif Vector:
            p = Vector(p.x, p.y, p.z, **kwds)
    except (IndexError, TypeError):
        raise _IndexError(fi=fi, points=points, wrap=w, txt=fractional.__name__)
    return p


def _fractional(points, fi, j, fin=None, wrap=None):  # in .frechet.py
    '''(INTERNAL) Compute point at L{fractional} index C{fi} and C{j}.
    '''
    i = int(fi)
    p = points[i]
    r = fi - float(i)
    if r > EPS:  # EPS0?
        if j is None:  # in .frechet.py
            j = i + 1
            if fin:
                j %= fin
        q = points[j]
        if r >= EPS1:  # PYCHOK no cover
            p = q
        elif wrap is not None:  # in (True, False)
            p = p.intermediateTo(q, r, wrap=wrap)
        elif _isLatLon(p):  # backward compatible default
            p = LatLon2Tuple(favg(p.lat, q.lat, f=r),
                             favg(p.lon, q.lon, f=r),
                             name=fractional.__name__)
        else:  # assume p and q are cartesian or vectorial
            z = p.z if p.z is q.z else favg(p.z, q.z, f=r)
            p = Vector3Tuple(favg(p.x, q.x, f=r),
                             favg(p.y, q.y, f=r), z,
                             name=fractional.__name__)
    return p


def isclockwise(points, adjust=False, wrap=True):
    '''Determine the direction of a path or polygon.

       @arg points: The path or polygon points (C{LatLon}[]).
       @kwarg adjust: Adjust the wrapped, unrolled longitudinal delta
                      by the cosine of the mean latitude (C{bool}).
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                    B{C{points}} (C{bool}).

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
    if a > 0:  # opposite of ellipsoidalExact and -Karney
        return True
    elif a < 0:
        return False
    # <https://blog.Element84.com/determining-if-a-spherical-polygon-contains-a-pole.html>
    raise _areaError(pts)


def isconvex(points, adjust=False, wrap=False):  # was=True
    '''Determine whether a polygon is convex.

       @arg points: The polygon points (C{LatLon}[]).
       @kwarg adjust: Adjust the wrapped, unrolled longitudinal delta
                      by the cosine of the mean latitude (C{bool}).
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                    B{C{points}} (C{bool}).

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


def isconvex_(points, adjust=False, wrap=False):  # was=True
    '''Determine whether a polygon is convex I{and clockwise}.

       @arg points: The polygon points (C{LatLon}[]).
       @kwarg adjust: Adjust the wrapped, unrolled longitudinal delta
                      by the cosine of the mean latitude (C{bool}).
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                    B{C{points}} (C{bool}).

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
    if adjust:
        def _unroll2(x1, x2, w, y1, y2):
            x21, x2 = unroll180(x1, x2, wrap=w)
            y = radians(y1 + y2) * _0_5
            x21 *= cos(y) if fabs(y) < PI_2 else _0_0
            return x21, x2
    else:
        def _unroll2(x1, x2, w, *unused):  # PYCHOK expected
            return unroll180(x1, x2, wrap=w)

    c, s = crosserrors(), 0

    Ps = LatLon2PsxyIter(points, loop=2, wrap=wrap)
    x1, y1, _ = Ps[0]
    x2, y2, _ = Ps[1]

    x21, x2 = _unroll2(x1, x2, False, y1, y2)
    for i, p in Ps.enumerate(closed=True):
        x3, y3, ll = p
        x32, x3 = _unroll2(x2, x3, bool(wrap and not Ps.looped), y2, y3)

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

        elif c and fdot((x32, y1 - y2), y3 - y2, -x21) < 0:  # PYCHOK no cover
            # colinear u-turn: x3, y3 not on the
            # opposite side of x2, y2 as x1, y1
            t = Fmt.SQUARE(points=i)
            raise CrossError(t, ll, txt=_colinear_)

        x1, y1, x2, y2, x21 = x2, y2, x3, y3, x32

    return s  # all points on the same side


def isenclosedBy(point, points, wrap=False):  # MCCABE 15
    '''Determine whether a point is enclosed by a polygon or composite.

       @arg point: The point (C{LatLon} or 2-tuple C{(lat, lon)}).
       @arg points: The polygon points or clips (C{LatLon}[], L{BooleanFHP}
                    or L{BooleanGH}).
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                    B{C{points}} (C{bool}).

       @return: C{True} if the B{C{point}} is inside the polygon or
                composite, C{False} otherwise.

       @raise PointsError: Insufficient number of B{C{points}}

       @raise TypeError: Some B{C{points}} are not C{LatLon}.

       @raise ValueError: Invalid B{C{point}}, lat- or longitude.

       @see: Functions L{pygeodesy.isconvex} and L{pygeodesy.ispolar} especially
             if the B{C{points}} may enclose a pole or wrap around the earth
             I{longitudinally}, methods L{sphericalNvector.LatLon.isenclosedBy},
             L{sphericalTrigonometry.LatLon.isenclosedBy} and U{MultiDop
             GeogContainPt<https://GitHub.com/NASA/MultiDop>} (U{Shapiro et.al. 2009,
             JTECH<https://Journals.AMetSoc.org/doi/abs/10.1175/2009JTECHA1256.1>}
             and U{Potvin et al. 2012, JTECH <https://Journals.AMetSoc.org/doi/abs/
             10.1175/JTECH-D-11-00019.1>}).
    '''
    try:
        y0, x0 = point.lat, point.lon
    except AttributeError:
        try:
            y0, x0 = map(float, point[:2])
        except (IndexError, TypeError, ValueError) as x:
            raise _ValueError(point=point, cause=x)

    if wrap:
        y0, x0 = _Wrap.latlon(y0, x0)

        def _dxy3(x, x2, y2, Ps):
            dx, x2 = unroll180(x, x2, wrap=not Ps.looped)
            return dx, x2, y2

    else:
        x0 = fmod(x0, _360_0)  # not x0 % 360!
        x0_180_ = x0 - _180_0
        x0_180  = x0 + _180_0

        def _dxy3(x1, x, y, unused):  # PYCHOK expected
            x = _umod_360(float(x))
            if x < x0_180_:
                x += _360_0
            elif x >= x0_180:
                x -= _360_0
            return (x - x1), x, y

    if _MODS.booleans.isBoolean(points):
        return points._encloses(y0, x0, wrap=wrap)

    Ps = LatLon2PsxyIter(points, loop=1, wrap=wrap)
    p  = Ps[0]
    e  = m = False
    S  = Fsum()

    _, x1, y1 = _dxy3(x0, p.x, p.y, False)
    for p in Ps.iterate(closed=True):
        dx, x2, y2 = _dxy3(x1, p.x, p.y, Ps)
        # ignore duplicate and near-duplicate pts
        if fabs(dx) > EPS or fabs(y2 - y1) > EPS:
            # determine if polygon edge (x1, y1)..(x2, y2) straddles
            # point (lat, lon) or is on boundary, but do not count
            # edges on boundary as more than one crossing
            if fabs(dx) < 180 and (x1 < x0 <= x2 or x2 < x0 <= x1):
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
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                    the B{C{points}} (C{bool}).

       @return: C{True} if the polygon encloses a pole, C{False}
                otherwise.

       @raise PointsError: Insufficient number of B{C{points}}

       @raise TypeError: Some B{C{points}} are not C{LatLon} or don't
                         have C{bearingTo2}, C{initialBearingTo}
                         and C{finalBearingTo} methods.
    '''
    def _cds(ps, w):  # iterate over course deltas
        Ps     =  PointsIter(ps, loop=2, wrap=w)
        p2, p1 =  Ps[0:2]
        b1, _  = _bearingTo2(p2, p1, wrap=False)
        for p2 in Ps.iterate(closed=True):
            if not p2.isequalTo(p1, EPS):
                if w and not Ps.looped:
                    p2 = _unrollon(p1, p2)
                b,  b2 = _bearingTo2(p1, p2, wrap=False)
                yield wrap180(b - b1)  # (b - b1 + 540) % 360 - 180
                yield wrap180(b2 - b)  # (b2 - b + 540) % 360 - 180
                p1, b1 = p2, b2

    # summation of course deltas around pole is 0° rather than normally ±360°
    # <https://blog.Element84.com/determining-if-a-spherical-polygon-contains-a-pole.html>
    s = fsum(_cds(points, wrap), floats=True)
    # XXX fix (intermittant) edge crossing pole - eg (85,90), (85,0), (85,-90)
    return fabs(s) < 90  # "zero-ish"


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
    t = (LatLon(   _0_0, lon1, **LatLon_kwds),
         LatLon(  _90_0, lon1, **LatLon_kwds),
         LatLon(   _0_0, lon2, **LatLon_kwds),
         LatLon(_N_90_0, lon2, **LatLon_kwds))
    if closed:
        t += t[:1]
    return t


def nearestOn5(point, points, closed=False, wrap=False, adjust=True,
                                            limit=9, **LatLon_and_kwds):
    '''Locate the point on a path or polygon closest to a reference point.

       The closest point on each polygon edge is either the nearest of that
       edge's end points or a point in between.

       @arg point: The reference point (C{LatLon}).
       @arg points: The path or polygon points (C{LatLon}[]).
       @kwarg closed: Optionally, close the path or polygon (C{bool}).
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                    B{C{points}} (C{bool}).
       @kwarg adjust: See function L{pygeodesy.equirectangular_} (C{bool}).
       @kwarg limit: See function L{pygeodesy.equirectangular_} (C{degrees}),
                     default C{9 degrees} is about C{1,000 Kmeter} (for mean
                     spherical earth radius L{R_KM}).
       @kwarg LatLon_and_kwds: Optional, C{B{LatLon}=None} class to use for
                     the closest point and additional B{C{LatLon}} keyword
                     arguments, ignored if C{B{LatLon}=None} or not given.

       @return: A L{NearestOn3Tuple}C{(closest, distance, angle)} with the
                {closest} point (B{C{LatLon}}) or if C{B{LatLon} is None},
                a L{NearestOn5Tuple}C{(lat, lon, distance, angle, height)}.
                The C{distance} is the L{pygeodesy.equirectangular} distance
                between the C{closest} and reference B{C{point}} in C{degrees}.
                The C{angle} from the B{C{point}} to the C{closest} is in
                compass C{degrees}, like function L{pygeodesy.compassAngle}.

       @raise LimitError: Lat- and/or longitudinal delta exceeds the B{C{limit}},
                          see function L{pygeodesy.equirectangular_}.

       @raise PointsError: Insufficient number of B{C{points}}

       @raise TypeError: Some B{C{points}} are not C{LatLon}.

       @note: Distances are I{approximated} by function L{pygeodesy.equirectangular_}.
              For more accuracy use one of the C{LatLon.nearestOn6} methods.

       @see: Function L{pygeodesy.degrees2m}.
    '''
    def _d2yx4(p2, p1, u, alw):
        # w = wrap if (i < (n - 1) or not closed) else False
        # equirectangular_ returns a Distance4Tuple(distance
        # in degrees squared, delta lat, delta lon, p2.lon
        # unroll/wrap'd); the previous p2.lon unroll/wrap'd
        # is also applied to the next edge's p1.lon
        return equirectangular_(p1.lat, p1.lon + u,
                                p2.lat, p2.lon, **alw)

    def _h(p):  # get height or default 0
        return _xattr(p, height=0) or 0

    # 3-D version used in .vector3d._nearestOn2
    #
    # point (x, y) on axis rotated ccw by angle a:
    #   x' = y * sin(a) + x * cos(a)
    #   y' = y * cos(a) - x * sin(a)
    #
    # distance (w) along and (h) perpendicular to
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
    #
    # i.e. no need for sqrt or hypot

    Ps = PointsIter(points, loop=1, wrap=wrap)
    p1 = c = Ps[0]
    u1 = u = _0_0
    kw = dict(adjust=adjust, limit=limit, wrap=False)
    d, dy, dx, _ = _d2yx4(p1, point, u1, kw)
    for p2 in Ps.iterate(closed=closed):
        # iff wrapped, unroll lon1 (actually previous
        # lon2) like function unroll180/-PI would've
        if wrap:
            kw.update(wrap=not (closed and Ps.looped))
        d21, y21, x21, u2 = _d2yx4(p2, p1, u1, kw)
        if d21 > EPS:
            # distance point to p1, y01 and x01 negated
            d2, y01, x01, _ = _d2yx4(point, p1, u1, kw)
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
                    d2, y01, x01, _ = _d2yx4(point, p1, u1, kw)
            if d2 < d:  # p1 is closer, y01 and x01 negated
                c, u, d, dy, dx = p1, u1, d2, -y01, -x01
        p1, u1 = p2, u2

    a =  atan2b(dx, dy)
    d =  hypot(dx, dy)
    h = _h(c)
    n =  nameof(point) or nearestOn5.__name__
    if LatLon_and_kwds:
        kwds = _xkwds(LatLon_and_kwds, height=h, name=n)
        LL   = _xkwds_pop(kwds, LatLon=None)
        if LL is not None:
            r = LL(c.lat, c.lon + u, **kwds)
            return NearestOn3Tuple(r, d, a, name=n)
    return NearestOn5Tuple(c.lat, c.lon + u, d, a, h, name=n)  # PYCHOK expected


def perimeterOf(points, closed=False, adjust=True, radius=R_M, wrap=True):
    '''I{Approximate} the perimeter of a path, polygon. or composite.

       @arg points: The path or polygon points or clips (C{LatLon}[],
                    L{BooleanFHP} or L{BooleanGH}).
       @kwarg closed: Optionally, close the path or polygon (C{bool}).
       @kwarg adjust: Adjust the wrapped, unrolled longitudinal delta
                      by the cosine of the mean latitude (C{bool}).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                    B{C{points}} (C{bool}).

       @return: Approximate perimeter (C{meter}, same units as
                B{C{radius}}).

       @raise PointsError: Insufficient number of B{C{points}}

       @raise TypeError: Some B{C{points}} are not C{LatLon}.

       @raise ValueError: Invalid B{C{radius}} or C{B{closed}=False} with
                          C{B{points}} a composite.

       @note: This perimeter is based on the L{pygeodesy.equirectangular_}
              distance approximation and is ill-suited for regions exceeding
              several hundred Km or Miles or with near-polar latitudes.

       @see: Functions L{sphericalTrigonometry.perimeterOf} and
             L{ellipsoidalKarney.perimeterOf}.
    '''
    def _degs(ps, c, a, w):  # angular edge lengths in degrees
        Ps = LatLon2PsxyIter(ps, loop=1)  # wrap=w
        p1, u = Ps[0], _0_0  # previous x2's unroll/wrap
        for p2 in Ps.iterate(closed=c):
            if w and c:
                w = not Ps.looped
            # apply previous x2's unroll/wrap'd to new x1
            _, dy, dx, u = equirectangular_(p1.y, p1.x + u,
                                            p2.y, p2.x,
                                            adjust=a, limit=None,
                                            wrap=w)  # PYCHOK non-seq
            yield hypot(dx, dy)
            p1 = p2

    if _MODS.booleans.isBoolean(points):
        if not closed:
            notImplemented(None, closed=closed, points=_composite_)
        d = points._sum1(perimeterOf, closed=True, adjust=adjust,
                                      radius=radius, wrap=wrap)
    else:
        d = fsum(_degs(points, closed, adjust, wrap), floats=True)
    return degrees2m(d, radius=radius)


def quadOf(latS, lonW, latN, lonE, closed=False, LatLon=LatLon_, **LatLon_kwds):
    '''Generate a quadrilateral path or polygon from two points.

       @arg latS: Souther-nmost latitude (C{degrees90}).
       @arg lonW: Western-most longitude (C{degrees180}).
       @arg latN: Norther-nmost latitude (C{degrees90}).
       @arg lonE: Eastern-most longitude (C{degrees180}).
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


__all__ += _ALL_DOCS(_Array2LatLon, _Basequence)

# **) MIT License
#
# Copyright (C) 2016-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
