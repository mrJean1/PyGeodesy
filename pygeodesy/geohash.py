
# -*- coding: utf-8 -*-

u'''Classes L{Geohash} and L{GeohashError} and several functions to
encode, decode and inspect I{geohashes}.

Transcribed from JavaScript originals by I{(C) Chris Veness 2011-2015}
and published under the same MIT Licence**, see U{Geohashes
<https://www.Movable-Type.co.UK/scripts/geohash.html>}.

See also U{Geohash<https://WikiPedia.org/wiki/Geohash>},
U{Geohash<https://GitHub.com/vinsci/geohash>},
U{PyGeohash<https://PyPI.org/project/pygeohash>} and
U{Geohash-Javascript<https://GitHub.com/DaveTroy/geohash-js>}.

@newfield example: Example, Examples
@newfield JSname: JS name, JS names
'''

from pygeodesy.basics import isstr, map2
from pygeodesy.dms import parse3llh  # parseDMS2
from pygeodesy.errors import _ValueError, _xkwds
from pygeodesy.fmath import favg
from pygeodesy.formy import equirectangular_ as _equirectangular_, \
                            equirectangular, euclidean, haversine, vincentys
from pygeodesy.interns import EPS, NN, R_M, _COMMA_, _DOT_, _E_, \
                             _floatuple, _N_, _NE_, _NW_, _S_, _SE_, \
                             _SW_, _W_, _0_0, _0_5, _180_0, _360_0
from pygeodesy.interns import _90_0  # PYCHOK used!
from pygeodesy.lazily import _ALL_LAZY, _ALL_OTHER
from pygeodesy.named import _NamedDict, _NamedTuple, nameof, _xnamed
from pygeodesy.namedTuples import Bounds2Tuple, Bounds4Tuple, \
                                  LatLon2Tuple, PhiLam2Tuple
from pygeodesy.props import deprecated_function, deprecated_method, \
                            deprecated_property_RO, Property_RO
from pygeodesy.streprs import fstr
from pygeodesy.units import Degrees_, Int, Lat, Lon, Precision_, Str, \
                           _xStrError

from math import ldexp, log10, radians

__all__ = _ALL_LAZY.geohash
__version__ = '21.04.15'


class _GH(object):
    '''(INTERNAL) Lazily defined constants.
    '''
    def _4d(self, n, e, s, w):  # helper
        return dict(N=(n, e), S=(s, w),
                    E=(e, n), W=(w, s))

    @Property_RO
    def Borders(self):
        return self._4d('prxz', 'bcfguvyz', '028b', '0145hjnp')

    Bounds4 = (-_90_0, -_180_0, _90_0, _180_0)

    @Property_RO
    def DecodedBase32(self):  # inverse GeohashBase32 map
        return dict((c, i) for i, c in enumerate(self.GeohashBase32))

    # Geohash-specific base32 map
    GeohashBase32 = '0123456789bcdefghjkmnpqrstuvwxyz'  # no a, i, j and o

    @Property_RO
    def Neighbors(self):
        return self._4d('p0r21436x8zb9dcf5h7kjnmqesgutwvy',
                        'bc01fg45238967deuvhjyznpkmstqrwx',
                        '14365h7k9dcfesgujnmqp0r2twvyx8zb',
                        '238967debc01fg45kmstqrwxuvhjyznp')

    @Property_RO
    def Sizes(self):  # lat-, lon and radial size (in meter)
        # ... where radial = sqrt(latSize * lonWidth / PI)
        return (_floatuple(20032e3, 20000e3, 11292815.096),  # 0
                _floatuple( 5003e3,  5000e3,  2821794.075),  # 1
                _floatuple(  650e3,  1225e3,   503442.397),  # 2
                _floatuple(  156e3,   156e3,    88013.575),  # 3
                _floatuple(  19500,   39100,    15578.683),  # 4
                _floatuple(   4890,    4890,     2758.887),  # 5
                _floatuple(    610,    1220,      486.710),  # 6
                _floatuple(    153,     153,       86.321),  # 7
                _floatuple(     19.1,    38.2,     15.239),  # 8
                _floatuple(      4.77,    4.77,     2.691),  # 9
                _floatuple(      0.596,   1.19,     0.475),  # 10
                _floatuple(      0.149,   0.149,    0.084),  # 11
                _floatuple(      0.0186,  0.0372,   0.015))  # 12  _MaxPrec

_GH      = _GH()  # PYCHOK singleton
_MaxPrec =  12


def _2bounds(LatLon, LatLon_kwds, s, w, n, e, name=NN):
    '''(INTERNAL) Return SW and NE bounds.
    '''
    if LatLon is None:
        r = Bounds4Tuple(s, w, n, e, name=name)
    else:
        sw = _xnamed(LatLon(s, w, **LatLon_kwds), name)
        ne = _xnamed(LatLon(n, e, **LatLon_kwds), name)
        r  =  Bounds2Tuple(sw, ne, name=name)
    return r  # _xnamed(r, name)


def _2center(bounds):
    '''(INTERNAL) Return the C{bounds} center.
    '''
    return (favg(bounds.latN, bounds.latS),
            favg(bounds.lonE, bounds.lonW))


def _2fll(lat, lon, *unused):
    '''(INTERNAL) Convert lat, lon to 2-tuple of floats.
    '''
    # lat, lon = parseDMS2(lat, lon)
    return (Lat(lat, Error=GeohashError),
            Lon(lon, Error=GeohashError))


def _2Geohash(geohash):
    '''(INTERNAL) Check or create a Geohash instance.
    '''
    return geohash if isinstance(geohash, Geohash) else \
                         Geohash(geohash)


def _2geostr(geohash):
    '''(INTERNAL) Check a geohash string.
    '''
    try:
        if not (0 < len(geohash) <= _MaxPrec):
            raise ValueError
        geostr = geohash.lower()
        for c in geostr:
            if c not in _GH.DecodedBase32:
                raise ValueError
        return geostr
    except (AttributeError, TypeError, ValueError) as x:
        raise GeohashError(Geohash.__name__, geohash, txt=str(x))


class Geohash(Str):
    '''Geohash class, a named C{str}.
    '''
    # no str.__init__ in Python 3
    def __new__(cls, cll, precision=None, name=NN):
        '''New L{Geohash} from an other L{Geohash} instance or C{str}
           or from a C{LatLon} instance or C{str}.

           @arg cll: Cell or location (L{Geohash}, C{LatLon} or C{str}).
           @kwarg precision: Optional, the desired geohash length (C{int}
                             1..12), see function L{geohash.encode} for
                             some examples.
           @kwarg name: Optional name (C{str}).

           @return: New L{Geohash}.

           @raise GeohashError: INValid or non-alphanumeric B{C{cll}}.

           @raise TypeError: Invalid B{C{cll}}.
        '''
        ll = None

        if isinstance(cll, Geohash):
            gh = _2geostr(str(cll))

        elif isstr(cll):
            if _COMMA_ in cll:
                ll = _2fll(*parse3llh(cll))
                gh =  encode(*ll, precision=precision)
            else:
                gh = _2geostr(cll)

        else:  # assume LatLon
            try:
                ll = _2fll(cll.lat, cll.lon)
                gh =  encode(*ll, precision=precision)
            except AttributeError:
                raise _xStrError(Geohash, cll=cll, Error=GeohashError)

        self = Str.__new__(cls, gh, name=name or nameof(cll))
        self._latlon = ll
        return self

    @deprecated_property_RO
    def ab(self):
        '''DEPRECATED, use property C{philam}.'''
        return self.philam

    def adjacent(self, direction, name=NN):
        '''Determine the adjacent cell in given compass direction.

           @arg direction: Compass direction ('N', 'S', 'E' or 'W').
           @kwarg name: Optional name (C{str}), otherwise the name
                        of this cell plus C{.D}irection.

           @return: Geohash of adjacent cell (L{Geohash}).

           @raise GeohashError: Invalid geohash or B{C{direction}}.
        '''
        # based on <https://GitHub.com/DaveTroy/geohash-js>

        D = direction[:1].upper()
        if D not in _GH.Neighbors:
            raise GeohashError(direction=direction)

        e = len(self) & 1  # % 2

        c = self[-1:]  # last hash char
        i = _GH.Neighbors[D][e].find(c)
        if i < 0:
            raise GeohashError(geohash=self)

        p = self[:-1]  # hash without last char
        # check for edge-cases which don't share common prefix
        if p and (c in _GH.Borders[D][e]):
            p = Geohash(p).adjacent(D)

        n = name or self.name
        if n:
            n = _DOT_(n, D)
        # append letter for direction to parent
        return Geohash(p + _GH.GeohashBase32[i], name=n)

    @Property_RO
    def _bounds(self):
        '''(INTERNAL) Cache for L{bounds}.
        '''
        return bounds(self)

    def bounds(self, LatLon=None, **LatLon_kwds):
        '''Return the lower-left SW and upper-right NE bounds of this
           geohash cell.

           @kwarg LatLon: Optional class to return I{bounds} (C{LatLon})
                          or C{None}.
           @kwarg LatLon_kwds: Optional keyword arguments for B{C{LatLon}},
                               ignored if B{C{LatLon}} is C{None}.

           @return: A L{Bounds2Tuple}C{(latlonSW, latlonNE)} of B{C{LatLon}}s
                    or a L{Bounds4Tuple}C{(latS, lonW, latN, lonE)} if
                    C{B{LatLon}=None},
        '''
        r = self._bounds
        return r if LatLon is None else \
              _2bounds(LatLon, LatLon_kwds, *r, name=self.name)

    def _distanceTo(self, func_, other, **kwds):
        '''(INTERNAL) Helper for distances, see C{formy._distanceTo*}.
        '''
        lls = self.latlon + _2Geohash(other).latlon
        return func_(*lls, **kwds)

    def distanceTo(self, other):
        '''Estimate the distance between this and an other geohash
           based the cell sizes.

           @arg other: The other geohash (L{Geohash}, C{LatLon} or C{str}).

           @return: Approximate distance (C{meter}).

           @raise TypeError: The B{C{other}} is not a L{Geohash},
                             C{LatLon} or C{str}.
        '''
        other = _2Geohash(other)

        n = min(len(self), len(other), len(_GH.Sizes))
        if n:
            for n in range(n):
                if self[n] != other[n]:
                    break
        return _GH.Sizes[n][2]

    @deprecated_method
    def distance1To(self, other):  # PYCHOK no cover
        '''DEPRECATED, use method L{distanceTo}.'''
        return self.distanceTo(other)

    distance1 = distance1To

    @deprecated_method
    def distance2To(self, other, radius=R_M, adjust=False, wrap=False):  # PYCHOK no cover
        '''DEPRECATED, use method L{equirectangularTo}.'''
        return self.equirectangularTo(other, radius=radius, adjust=adjust, wrap=wrap)

    distance2 = distance2To

    @deprecated_method
    def distance3To(self, other, radius=R_M, wrap=False):  # PYCHOK no cover
        '''DEPRECATED, use method L{haversineTo}.'''
        return self.haversineTo(other, radius=radius, wrap=wrap)

    distance3 = distance3To

    def equirectangularTo(self, other, radius=R_M, adjust=False, wrap=False):
        '''Approximate the distance between this and an other geohash
           using the L{equirectangular} function.

           @arg other: The other geohash (L{Geohash}, C{LatLon} or C{str}).
           @kwarg radius: Mean earth radius, ellipsoid or datum
                          (C{meter}, L{Ellipsoid}, L{Ellipsoid2},
                          L{Datum} or L{a_f2Tuple}) or C{None}.
           @kwarg adjust: Adjust the wrapped, unrolled longitudinal
                          delta by the cosine of the mean latitude
                          C{bool}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Distance (C{meter}, same units as B{C{radius}} or the
                    ellipsoid or datum axes or C{radians I{squared}} if
                    B{C{radius}} is C{None} or C{0}).

           @raise TypeError: The B{C{other}} is not a L{Geohash}, C{LatLon}
                             or C{str} or invalid B{C{radius}}.

           @see: U{Local, flat earth approximation
                 <https://www.EdWilliams.org/avform.htm#flat>}, functions
        '''
        lls  = self.latlon + _2Geohash(other).latlon
        kwds = dict(adjust=adjust, limit=None, wrap=wrap)
        return equirectangular( *lls, radius=radius, **kwds) if radius else \
              _equirectangular_(*lls, **kwds).distance2

    def euclideanTo(self, other, radius=R_M, adjust=False, wrap=False):
        '''Approximate the distance between this and an other geohash
           using the L{euclidean} function.

           @arg other: The other geohash (L{Geohash}, C{LatLon} or C{str}).
           @kwarg radius: Mean earth radius, ellipsoid or datum
                          (C{meter}, L{Ellipsoid}, L{Ellipsoid2},
                          L{Datum} or L{a_f2Tuple}).
           @kwarg adjust: Adjust the wrapped, unrolled longitudinal
                          delta by the cosine of the mean latitude
                          C{bool}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Distance (C{meter}, same units as B{C{radius}} or the
                    ellipsoid or datum axes).

           @raise TypeError: The B{C{other}} is not a L{Geohash}, C{LatLon}
                             or C{str} or invalid B{C{radius}}.
        '''
        return self._distanceTo(euclidean, other, radius=radius,
                                           adjust=adjust, wrap=wrap)

    def haversineTo(self, other, radius=R_M, wrap=False):
        '''Compute the distance between this and an other geohash using
           the L{haversine} formula.

           @arg other: The other geohash (L{Geohash}, C{LatLon} or C{str}).
           @kwarg radius: Mean earth radius, ellipsoid or datum
                          (C{meter}, L{Ellipsoid}, L{Ellipsoid2},
                          L{Datum} or L{a_f2Tuple}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Distance (C{meter}, same units as B{C{radius}} or the
                    ellipsoid or datum axes).

           @raise TypeError: The B{C{other}} is not a L{Geohash}, C{LatLon}
                             or C{str} or invalid B{C{radius}}.
        '''
        return self._distanceTo(haversine, other, radius=radius, wrap=wrap)

    @Property_RO
    def latlon(self):
        '''Get the lat- and longitude of (the approximate center of)
           this geohash as a L{LatLon2Tuple}C{(lat, lon)} in C{degrees}.

           @example:

            >>> geohash.Geohash('geek').latlon  # 65.478515625, -17.75390625
            >>> geohash.decode('geek')  # '65.48', '-17.75'
        '''
        lat, lon = self._latlon or _2center(self.bounds())
        return LatLon2Tuple(lat, lon, name=self.name)

    @Property_RO
    def neighbors(self):
        '''Get all 8 adjacent cells as a L{Neighbors8Dict}C{(N, NE,
           E, SE, S, SW, W, NW)} of L{Geohash}es.

           @JSname: I{neighbours}.
        '''
        return Neighbors8Dict(N=self.N, NE=self.NE, E=self.E, SE=self.SE,
                              S=self.S, SW=self.SW, W=self.W, NW=self.NW,
                              name=self.name)

    @Property_RO
    def philam(self):
        '''Get the lat- and longitude of (the approximate center of)
           this geohash as a L{PhiLam2Tuple}C{(phi, lam)} in C{radians}.
        '''
        return PhiLam2Tuple(*map2(radians, self.latlon), name=self.name)

    @Property_RO
    def precision(self):
        '''Get this geohash's precision (C{int}).
        '''
        return len(self)

    @Property_RO
    def sizes(self):
        '''Get the lat- and longitudinal size of this cell as
           a L{LatLon2Tuple}C{(lat, lon)} in (C{meter}).
        '''
        z = _GH.Sizes
        n =  min(len(z) - 1, max(self.precision, 1))
        return LatLon2Tuple(*z[n][:2], name=self.name)  # XXX Height, Width?

    def toLatLon(self, LatLon=None, **LatLon_kwds):
        '''Return (the approximate center of) this geohash cell
           as an instance of the supplied C{LatLon} class.

           @arg LatLon: Class to use (C{LatLon}) or C{None}.
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}}
                               keyword arguments, ignored if
                               C{B{LatLon}=None}.

           @return: This geohash location (B{C{LatLon}}) or a
                    L{LatLon2Tuple}C{(lat, lon)} if B{C{LatLon}}
                    is C{None}.

           @raise TypeError: Invalid B{C{LatLon}} or B{C{LatLon_kwds}}.

           @example:

            >>> from sphericalTrigonometry import LatLon
            >>> ll = Geohash('u120fxw').toLatLon(LatLon)
            >>> print(repr(ll))  # LatLon(52°12′17.9″N, 000°07′07.64″E)
            >>> print(ll)  # 52.204971°N, 000.11879°E
        '''
        return self.latlon if LatLon is None else _xnamed(LatLon(
              *self.latlon, **LatLon_kwds), self.name)

    def vincentysTo(self, other, radius=R_M, wrap=False):
        '''Compute the distance between this and an other geohash using
           the L{vincentys} formula.

           @arg other: The other geohash (L{Geohash}, C{LatLon} or C{str}).
           @kwarg radius: Mean earth radius, ellipsoid or datum
                          (C{meter}, L{Ellipsoid}, L{Ellipsoid2},
                          L{Datum} or L{a_f2Tuple}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Distance (C{meter}, same units as B{C{radius}} or the
                    ellipsoid or datum axes).

           @raise TypeError: The B{C{other}} is not a L{Geohash}, C{LatLon}
                             or C{str} or invalid B{C{radius}}.
        '''
        return self._distanceTo(vincentys, other, radius=radius, wrap=wrap)

    @Property_RO
    def N(self):
        '''Get the cell North of this (L{Geohash}).
        '''
        return self.adjacent(_N_)

    @Property_RO
    def S(self):
        '''Get the cell South of this (L{Geohash}).
        '''
        return self.adjacent(_S_)

    @Property_RO
    def E(self):
        '''Get the cell East of this (L{Geohash}).
        '''
        return self.adjacent(_E_)

    @Property_RO
    def W(self):
        '''Get the cell West of this (L{Geohash}).
        '''
        return self.adjacent(_W_)

    @Property_RO
    def NE(self):
        '''Get the cell NorthEast of this (L{Geohash}).
        '''
        return self.N.E

    @Property_RO
    def NW(self):
        '''Get the cell NorthWest of this (L{Geohash}).
        '''
        return self.N.W

    @Property_RO
    def SE(self):
        '''Get the cell SouthEast of this (L{Geohash}).
        '''
        return self.S.E

    @Property_RO
    def SW(self):
        '''Get the cell SouthWest of this (L{Geohash}).
        '''
        return self.S.W


class GeohashError(_ValueError):
    '''Geohash encode, decode or other L{Geohash} issue.
    '''
    pass


class Neighbors8Dict(_NamedDict):
    '''8-Dict C{(N, NE, E, SE, S, SW, W, NW)} of L{Geohash}es,
       providing key I{and} attribute access to the items.
    '''
    _Keys_ = (_N_, _NE_, _E_, _SE_, _S_, _SW_, _W_, _NW_)

    def __init__(self, **kwds):  # PYCHOK no *args
        kwds = _xkwds(kwds, **_Neighbors8Defaults)
        _NamedDict.__init__(self, **kwds)  # name=...


_Neighbors8Defaults = dict(zip(Neighbors8Dict._Keys_, (None,) *
                           len(Neighbors8Dict._Keys_)))  # XXX frozendict


def bounds(geohash, LatLon=None, **LatLon_kwds):
    '''Returns the lower-left SW and upper-right NE corners of a geohash.

       @arg geohash: To be bound (L{Geohash}).
       @kwarg LatLon: Optional class to return the bounds (C{LatLon})
                      or C{None}.
       @kwarg LatLon_kwds: Optional keyword arguments for B{C{LatLon}},
                           ignored if C{B{LatLon}=None}.

       @return: A L{Bounds2Tuple}C{(latlonSW, latlonNE)} of B{C{LatLon}}s
                or if B{C{LatLon}} is C{None}, a L{Bounds4Tuple}C{(latS,
                lonW, latN, lonE)}.

       @raise TypeError: The B{C{geohash}} is not a L{Geohash}, C{LatLon}
                         or C{str} or invalid B{C{LatLon}} or invalid
                         B{C{LatLon_kwds}}.

       @raise GeohashError: Invalid or C{null} B{C{geohash}}.

       @example:

        >>> geohash.bounds('u120fxw')  #  52.20428467, 0.11810303, 52.20565796, 0.11947632
        >>> geohash.decode('u120fxw')  # '52.205',    '0.1188'
    '''
    gh = _2Geohash(geohash)
    if len(gh) < 1:
        raise GeohashError(geohash=geohash)

    s, w, n, e = _GH.Bounds4
    try:
        d = True
        for c in gh.lower():
            i = _GH.DecodedBase32[c]
            for m in (16, 8, 4, 2, 1):
                if d:  # longitude
                    if i & m:
                        w = favg(w, e)
                    else:
                        e = favg(w, e)
                else:  # latitude
                    if i & m:
                        s = favg(s, n)
                    else:
                        n = favg(s, n)
                d = not d
    except KeyError:
        raise GeohashError(geohash=geohash)

    return _2bounds(LatLon, LatLon_kwds, s, w, n, e,
                                   name=nameof(geohash))


def _bounds3(geohash):
    '''(INTERNAL) Return 3-tuple C{(bounds, height, width)}.
    '''
    b = bounds(geohash)
    return b, (b.latN - b.latS), (b.lonE - b.lonW)


def decode(geohash):
    '''Decode a geohash to lat-/longitude of the (approximate
       centre of) geohash cell to reasonable precision.

       @arg geohash: To be decoded (L{Geohash}).

       @return: 2-Tuple C{(latStr, lonStr)}, both C{str}.

       @raise TypeError: The B{C{geohash}} is not a L{Geohash},
                         C{LatLon} or C{str}.

       @raise GeohashError: Invalid or null B{C{geohash}}.

       @example:

        >>> geohash.decode('u120fxw')  # '52.205', '0.1188'
        >>> geohash.decode('sunny')  # '23.708', '42.473'  Saudi Arabia
        >>> geohash.decode('fur')  # '69.6', '-45.7'  Greenland
        >>> geohash.decode('reef')  # '-24.87', '162.95'  Coral Sea
        >>> geohash.decode('geek')  # '65.48', '-17.75'  Iceland
    '''
    b, h, w  = _bounds3(geohash)
    lat, lon = _2center(b)

    # round to near centre without excessive precision to
    # ⌊2-log10(Δ°)⌋ decimal places, strip trailing zeros
    return (fstr(lat, prec=int(2 - log10(h))),
            fstr(lon, prec=int(2 - log10(w))))  # strs!


def decode2(geohash, LatLon=None, **LatLon_kwds):
    '''Decode a geohash to lat-/longitude of the (approximate
       centre of) geohash cell to reasonable precision.

       @arg geohash: To be decoded (L{Geohash}).
       @kwarg LatLon: Optional class to return the location (C{LatLon})
                      or C{None}.
       @kwarg LatLon_kwds: Optional keyword arguments for B{C{LatLon}},
                           ignored if C{B{LatLon}=None}.

       @return: L{LatLon2Tuple}C{(lat, lon)}, both C{degrees} if
                C{B{LatLon}=None}, otherwise a B{C{LatLon}} instance.

       @raise TypeError: The B{C{geohash}} is not a L{Geohash},
                         C{LatLon} or C{str}.

       @raise GeohashError: Invalid or null B{C{geohash}}.
    '''
    t = map2(float, decode(geohash))
    r = LatLon2Tuple(*t) if LatLon is None else LatLon(*t, **LatLon_kwds)
    return _xnamed(r, decode2.__name__)


def decode_error(geohash):
    '''Return the relative lat-/longitude decoding errors for
       this geohash.

       @arg geohash: To be decoded (L{Geohash}).

       @return: A L{LatLon2Tuple}C{(lat, lon)} with the lat- and
                longitudinal errors in (C{degrees}).

       @raise TypeError: The B{C{geohash}} is not a L{Geohash},
                         C{LatLon} or C{str}.

       @raise GeohashError: Invalid or null B{C{geohash}}.

       @example:

        >>> geohash.decode_error('u120fxw')  # 0.00068665, 0.00068665
        >>> geohash.decode_error('fur')  # 0.703125, 0.703125
        >>> geohash.decode_error('fu')  # 2.8125, 5.625
        >>> geohash.decode_error('f')  # 22.5, 22.5
    '''
    _, h, w = _bounds3(geohash)
    return LatLon2Tuple(h * _0_5,  # Height error
                        w * _0_5)  # Width error


def distance_(geohash1, geohash2):
    '''Estimate the distance between two geohash (from the cell sizes).

       @arg geohash1: First geohash (L{Geohash}, C{LatLon} or C{str}).
       @arg geohash2: Second geohash (L{Geohash}, C{LatLon} or C{str}).

       @return: Approximate distance (C{meter}).

       @raise TypeError: If B{C{geohash1}} or B{C{geohash2}} is
                         not a L{Geohash}, C{LatLon} or C{str}.

       @example:

        >>> geohash.distance_('u120fxwsh', 'u120fxws0')  # 15.239
    '''
    return _2Geohash(geohash1).distanceTo(geohash2)


@deprecated_function
def distance1(geohash1, geohash2):
    '''DEPRECATED, used L{geohash.distance_}.'''
    return distance_(geohash1, geohash2)


@deprecated_function
def distance2(geohash1, geohash2):
    '''DEPRECATED, used L{geohash.equirectangular_}.'''
    return equirectangular_(geohash1, geohash2)


@deprecated_function
def distance3(geohash1, geohash2):
    '''DEPRECATED, used L{geohash.haversine_}.'''
    return haversine_(geohash1, geohash2)


def encode(lat, lon, precision=None):
    '''Encode a lat-/longitude as a C{geohash}, either to the specified
       precision or if not provided, to an automatically evaluated
       precision.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).
       @kwarg precision: Optional, the desired geohash length (C{int}
                         1..12).

       @return: The C{geohash} (C{str}).

       @raise GeohashError: Invalid B{C{lat}}, B{C{lon}} or B{C{precision}}.

       @example:

        >>> geohash.encode(52.205, 0.119,   7)  # 'u120fxw'
        >>> geohash.encode(52.205, 0.119,  12)  # 'u120fxwshvkg'
        >>> geohash.encode(52.205, 0.1188, 12)  # 'u120fxws0jre'
        >>> geohash.encode(52.205, 0.1188)      # 'u120fxw'
        >>> geohash.encode(     0, 0)           # 's00000000000'
    '''
    lat, lon = _2fll(lat, lon)

    if precision is None:
        # Infer precision by refining geohash until
        # it matches precision of supplied lat/lon.
        for p in range(1, _MaxPrec + 1):
            gh = encode(lat, lon, p)
            ll = map2(float, decode(gh))
            if abs(lat - ll[0]) < EPS and \
               abs(lon - ll[1]) < EPS:
                return gh
        p = _MaxPrec
    else:
        p = Precision_(precision, Error=GeohashError, low=1, high=_MaxPrec)

    b = i = 0
    d, gh = True, []
    s, w, n, e = _GH.Bounds4

    while p > 0:
        i += i
        if d:  # bisect longitude
            m = favg(e, w)
            if lon < m:
                e = m
            else:
                w = m
                i += 1
        else:  # bisect latitude
            m = favg(n, s)
            if lat < m:
                n = m
            else:
                s = m
                i += 1
        d = not d

        b += 1
        if b == 5:
            # 5 bits gives a character:
            # append it and start over
            gh.append(_GH.GeohashBase32[i])
            b  = i = 0
            p -= 1

    return NN.join(gh)


def equirectangular_(geohash1, geohash2, radius=R_M):
    '''Approximate the distance between two geohashes using the
       L{equirectangular} formula.

       @arg geohash1: First geohash (L{Geohash}, C{LatLon} or C{str}).
       @arg geohash2: Second geohash (L{Geohash}, C{LatLon} or C{str}).
       @kwarg radius: Mean earth radius (C{meter}) or C{None}.

       @return: Approximate distance (C{meter}, same units as
                B{C{radius}}).

       @raise TypeError: If B{C{geohash1}} or B{C{geohash2}} is
                         not a L{Geohash}, C{LatLon} or C{str}.

       @example:

        >>> geohash.equirectangular_('u120fxwsh', 'u120fxws0')  # 19.0879
    '''
    return _2Geohash(geohash1).equirectangularTo(geohash2, radius=radius)


def haversine_(geohash1, geohash2, radius=R_M):
    '''Compute the great-circle distance between two geohashes
       using the L{haversine} formula.

       @arg geohash1: First geohash (L{Geohash}, C{LatLon} or C{str}).
       @arg geohash2: Second geohash (L{Geohash}, C{LatLon} or C{str}).
       @kwarg radius: Mean earth radius (C{meter}).

       @return: Great-circle distance (C{meter}, same units as
                B{C{radius}}).

       @raise TypeError: If B{C{geohash1}} or B{C{geohash2}} is
                         not a L{Geohash}, C{LatLon} or C{str}.

       @example:

        >>> geohash.haversine_('u120fxwsh', 'u120fxws0')  # 11.6978
    '''
    return _2Geohash(geohash1).haversineTo(geohash2, radius=radius)


def neighbors(geohash):
    '''Return the L{Geohash}es for all 8 adjacent cells.

       @arg geohash: Cell for which neighbors are requested
                     (L{Geohash} or C{str}).

       @return: A L{Neighbors8Dict}C{(N, NE, E, SE, S, SW, W, NW)}
                of L{Geohash}es.

       @raise TypeError: The B{C{geohash}} is not a L{Geohash},
                         C{LatLon} or C{str}.

       @JSname: I{neighbours}.
    '''
    return _2Geohash(geohash).neighbors


def precision(res1, res2=None):
    '''Determine the L{Geohash} precisions to meet a or both given
       (geographic) resolutions.

       @arg res1: The required primary I{(longitudinal)} resolution
                  (C{degrees}).
       @kwarg res2: Optional, required secondary I{(latitudinal)}
                    resolution (C{degrees}).

       @return: The L{Geohash} precision or length (C{int}, 1..12).

       @raise GeohashError: Invalid B{C{res1}} or B{C{res2}}.

       @see: C++ class U{Geohash
             <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1Geohash.html>}.
    '''
    r = Degrees_(res1=res1, low=_0_0, Error=GeohashError)
    if res2 is None:
        t = r, r
        for p in range(1, _MaxPrec):
            if resolution2(p, None) <= t:
                return p

    else:
        t = r, Degrees_(res2=res2, low=_0_0, Error=GeohashError)
        for p in range(1, _MaxPrec):
            if resolution2(p, p) <= t:
                return p

    return _MaxPrec


class Resolutions2Tuple(_NamedTuple):
    '''2-Tuple C{(res1, res2)} with the primary I{(longitudinal)} and
       secondary I{(latitudinal)} resolution, both in C{degrees}.
    '''
    _Names_ = ('res1',   'res2')
    _Units_ = ( Degrees_, Degrees_)


def resolution2(prec1, prec2=None):
    '''Determine the (geographic) resolutions of given L{Geohash}
       precisions.

       @arg prec1: The given primary I{(longitudinal)} precision
                   (C{int} 1..12).
       @kwarg prec2: Optional, secondary I{(latitudinal)} precision
                     (C{int} 1..12).

       @return: L{Resolutions2Tuple}C{(res1, res2)} with the
                (geographic) resolutions C{degrees}, where C{res2}
                B{C{is}} C{res1} if no B{C{prec2}} is given.

       @raise GeohashError: Invalid B{C{prec1}} or B{C{prec2}}.

       @see: I{Karney}'s C++ class U{Geohash
             <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1Geohash.html>}.
    '''
    res1, res2 = _360_0, _180_0  # note ... lon, lat!

    if prec1:
        p = 5 * max(0, min(Int(prec1=prec1, Error=GeohashError), _MaxPrec))
        res1 = res2 = ldexp(res1, -(p - p // 2))

    if prec2:
        p = 5 * max(0, min(Int(prec2=prec2, Error=GeohashError), _MaxPrec))
        res2 = ldexp(res2, -(p // 2))

    return Resolutions2Tuple(res1, res2)


def sizes(geohash):
    '''Return the lat- and longitudinal size of this L{Geohash} cell.

       @arg geohash: Cell for which size are required (L{Geohash} or
                     C{str}).

       @return: A L{LatLon2Tuple}C{(lat, lon)} with the latitudinal
                height and longitudinal width in (C{meter}).

       @raise TypeError: The B{C{geohash}} is not a L{Geohash},
                         C{LatLon} or C{str}.
    '''
    return _2Geohash(geohash).sizes


__all__ += _ALL_OTHER(bounds,  # functions
                      decode, decode2, decode_error, distance_,
                      encode, equirectangular_, haversine_,
                      neighbors, precision, resolution2, sizes)

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
