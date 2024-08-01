
# -*- coding: utf-8 -*-

u'''I{Gustavo Niemeyer}’s Geohash en-, decoding and caching.

Classes L{Geohash} and L{GeohashError} and several functions to encode,
decode and inspect I{geohashes} and optional L{Geohashed} caches.

Transcoded from JavaScript originals by I{(C) Chris Veness 2011-2015} and
published under the same MIT Licence**, see
U{Geohashes<https://www.Movable-Type.co.UK/scripts/geohash.html>}.

@see: U{Geohash<https://WikiPedia.org/wiki/Geohash>},
      U{Geohash<https://GitHub.com/vinsci/geohash>},
      U{PyGeohash<https://PyPI.org/project/pygeohash>} and
      U{Geohash-Javascript<https://GitHub.com/DaveTroy/geohash-js>}.
'''

from pygeodesy.basics import isodd, isstr, map2
from pygeodesy.constants import EPS, R_M, _0_0, _0_5, _180_0, _360_0, \
                               _90_0, _N_90_0, _N_180_0  # PYCHOK used!
# from pygeodesy.dms import parse3llh  # _MODS
from pygeodesy.errors import _ValueError, _xkwds, _xStrError
# from pygeodesy import formy as _formy  # _MODS
from pygeodesy.interns import NN, _COMMA_, _DOT_, _E_, _height_, _N_, _NE_, \
                             _NW_, _radius_, _S_, _SE_, _SW_, _W_, _width_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import _name__, _NamedDict, _NamedTuple, nameof, _xnamed
from pygeodesy.namedTuples import Bounds2Tuple, Bounds4Tuple, LatLon2Tuple, \
                                  PhiLam2Tuple
from pygeodesy.props import deprecated_function, deprecated_method, \
                            deprecated_property_RO, Property_RO, \
                            property_RO, property_ROver
# from pygeodesy.streprs import Fmt, fstr  # _MODS
from pygeodesy.units import Degrees_, Int, Lat, Lon, Meter, Precision_, Str

from math import fabs, ldexp, log10, radians

__all__ = _ALL_LAZY.geohash
__version__ = '24.08.01'

_formy   = _MODS.into(formy=__name__)
_MASK5   =  16, 8, 4, 2, 1  # PYCHOK used!
_MaxPrec =  12


def _2avg(a, b):
    '''(INTERNAL) Bisect C{a} to C{b}.
    '''
    return (a + b) * _0_5  # favg


def _2avg_ndigits2(a, b):
    '''(INTERNAL) Return 2-tuple C{(_2avg, ndigits)}.
    '''
    # round to near centre without excessive
    # precision to ⌊2-log10(Δ°)⌋ ndigits
    return _2avg(a, b), int(2 - log10(a - b))


def _2bounds(LatLon, LatLon_kwds, s, w, n, e, **name):
    '''(INTERNAL) Return SW and NE bounds.
    '''
    if LatLon is None:
        r    =  Bounds4Tuple(s, w, n, e, **name)
    else:
        kwds = _xkwds(LatLon_kwds, **name)
        r    =  Bounds2Tuple(LatLon(s, w, **kwds),
                             LatLon(n, e, **kwds), **name)
    return r


def _2center(bounds):
    '''(INTERNAL) Return the C{bounds} center.
    '''
    return (_2avg(bounds.latN, bounds.latS),
            _2avg(bounds.lonE, bounds.lonW))


def _2dab(d, a, b):
    '''(INTERNAL) Get delta from rounded.
    '''
    return fabs(d - round(*_2avg_ndigits2(a, b)))


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


def _2latlon(s, w, n, e, fstr=None):
    '''(INTERNAL) Get the center C{lat, lon}, rounded.
    '''
    lat, a = _2avg_ndigits2(n, s)
    lon, b = _2avg_ndigits2(e, w)
    return (fstr(lat, prec=a), fstr(lon, prec=b)) if fstr else \
           (round(lat, a),     round(lon, b))


def _2Precision(p):
    '''(INTERNAL) Get a valid C{Precision}.
    '''
    return Precision_(p, Error=GeohashError, low=1, high=_MaxPrec)


class _GH(object):
    '''(INTERNAL) Lazily defined constants.
    '''
    def _4d(self, s, w, n, e):  # helper
        return dict(S=(s, w), W=(w, s),
                    N=(n, e), E=(e, n))

    @property_ROver
    def Borders(self):
        return self._4d('028b', '0145hjnp', 'prxz', 'bcfguvyz')

    @property_ROver
    def DecodeB32(self):  # inverse EncodeB32 map
        return dict((c, i) for i, c in enumerate(self.EncodeB32))

    def decode2(self, geohash):
        '''Decode C{geohash} to C{(lat, lon)}.
        '''
        swne = self.swne4(geohash)
        return _2latlon(*swne)

    # Geohash base32 map: no a, i, l and o
    EncodeB32 = '0123456789bcdefghjkmnpqrstuvwxyz'

    def encode(self, *lat_lon_prec_eps):
        '''Encode C{lat, lon} to C{prec}ision or C{eps}.
        '''
        def _encodes(lat, lon, prec, eps=0):
            s, w, n, e = self.SWNE4
            E, d, _avg = self.EncodeB32, True, _2avg
            for _ in range(prec):
                i = 0
                for _ in range(5):  # len(_MASK5)
                    i += i
                    if d:  # bisect longitude
                        a = _avg(e, w)
                        if lon < a:
                            e  = a
                        else:
                            w  = a
                            i += 1
                    else:  # bisect latitude
                        a = _avg(n, s)
                        if lat < a:
                            n  = a
                        else:
                            s  = a
                            i += 1
                    d = not d
                yield E[i]
                if eps > 0:  # infer prec
                    if _2dab(lon, e, w) < eps and \
                       _2dab(lat, n, s) < eps:
                        break

        return NN.join(_encodes(*lat_lon_prec_eps))

    @property_ROver
    def _LatLon2Tuple(self):

        class _LatLon2Tuple(_NamedTuple):
            '''DEPRECATED on 2024.07.28, C{(lat, lon)} in B{C{meter}}, use L{Sizes3Tuple}.'''
            _Names_ = LatLon2Tuple._Names_
            _Units_ = Meter, Meter

        return _LatLon2Tuple

    @property_ROver
    def Neighbors(self):
        return self._4d('14365h7k9dcfesgujnmqp0r2twvyx8zb',
                        '238967debc01fg45kmstqrwxuvhjyznp',
                        'p0r21436x8zb9dcf5h7kjnmqesgutwvy',
                        'bc01fg45238967deuvhjyznpkmstqrwx')

    @property_ROver
    def Sizes(self):  # height, width and radius (in meter)
        # where radius = sqrt(height * width / PI), the
        # radius of a circle with area (height * width)
        T = Sizes3Tuple
        return (T(20000e3, 20032e3, 11292815.096),  # 0
                T( 5000e3,  5003e3,  2821794.075),  # 1
                T(  650e3,  1225e3,   503442.397),  # 2
                T(  156e3,   156e3,    88013.575),  # 3
                T(  19500,   39100,    15578.683),  # 4
                T(   4890,    4890,     2758.887),  # 5
                T(    610,    1220,      486.710),  # 6
                T(    153,     153,       86.321),  # 7
                T(     19.1,    38.2,     15.239),  # 8
                T(      4.77,    4.77,     2.691),  # 9
                T(      0.596,   1.19,     0.475),  # 10
                T(      0.149,   0.149,    0.084),  # 11
                T(      0.0186,  0.0372,   0.015))  # 12  _MaxPrec

    SWNE4 = (_N_90_0, _N_180_0, _90_0, _180_0)

    def swne4(self, geohash, mask5=_MASK5):
        '''Decode C{geohash} into 4-tuple C{(s, w, n, e)}.
        '''
        nc = len(geohash) if isstr(geohash) else 0
        if not (0 < nc <= _MaxPrec):
            raise GeohashError(geohash=geohash, len=nc)
        s, w, n, e = self.SWNE4
        D, d, _avg = self.DecodeB32, True, _2avg
        try:
            for j, c in enumerate(geohash.lower()):
                i = D[c]
                for m in mask5:
                    if d:  # longitude
                        a = _avg(e, w)
                        if (i & m):
                            w = a
                        else:
                            e = a
                    else:  # latitude
                        a = _avg(n, s)
                        if (i & m):
                            s = a
                        else:
                            n = a
                    d = not d
        except KeyError:
            c = _MODS.streprs.Fmt.INDEX(repr(c), j)
            raise GeohashError(geohash=geohash, len=nc, txt=c)
        return s, w, n, e

_GH = _GH()  # PYCHOK singleton


class Geohash(Str):
    '''Geohash class, a named C{str}.
    '''
    # no str.__init__ in Python 3
    def __new__(cls, cll, precision=None, **name):
        '''New L{Geohash} from an other L{Geohash} instance or C{str}
           or from a C{LatLon} instance or C{str}.

           @arg cll: Cell or location (L{Geohash}, C{LatLon} or C{str}).
           @kwarg precision: Optional, the desired geohash length (C{int}
                             1..12).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: New L{Geohash}.

           @raise GeohashError: INValid or non-alphanumeric B{C{cll}}.

           @raise TypeError: Invalid B{C{cll}}.
        '''
        ll = None

        if isinstance(cll, Geohash):
            gh = str(cll)

        elif isstr(cll):  # "lat, lon" or "geohash"
            if _COMMA_ in cll:
                ll = _2fll(*_MODS.dms.parse3llh(cll))  # parseDMS2
                gh =  encode(*ll, precision=precision)
            else:
                gh =  cll.lower()
                _  = _GH.swne4(gh, mask5=())  # validate

        else:  # assume LatLon
            try:
                ll = _2fll(cll.lat, cll.lon)
                gh =  encode(*ll, precision=precision)
            except AttributeError:
                raise _xStrError(Geohash, cll=cll, Error=GeohashError)

        self = Str.__new__(cls, gh, name=_name__(name, _or_nameof=cll))
        self._latlon = ll
        return self

    @deprecated_property_RO
    def ab(self):
        '''DEPRECATED, use property C{philam}.'''
        return self.philam

    def adjacent(self, direction, **name):
        '''Determine the adjacent cell in the given compass direction.

           @arg direction: Compass direction ('N', 'S', 'E' or 'W').
           @kwarg name: Optional C{B{name}=NN} (C{str}) otherwise this
                        cell's name, either extended with C{.D}irection.

           @return: Geohash of adjacent cell (L{Geohash}).

           @raise GeohashError: Invalid geohash or B{C{direction}}.
        '''
        # based on <https://GitHub.com/DaveTroy/geohash-js>

        D = direction[:1].upper()
        if D not in _GH.Neighbors:
            raise GeohashError(direction=direction)

        e = 1 if isodd(len(self)) else 0

        c = self[-1:]  # last hash char
        i = _GH.Neighbors[D][e].find(c)
        if i < 0:
            raise GeohashError(geohash=self)

        p = self[:-1]  # hash without last char
        # check for edge-cases which don't share common prefix
        if p and (c in _GH.Borders[D][e]):
            p = Geohash(p).adjacent(D)

        n = self._name__(name)
        if n:
            n = _DOT_(n, D)
        # append letter for direction to parent
        return Geohash(p + _GH.EncodeB32[i], name=n)

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
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                               arguments, ignored if C{B{LatLon} is None}.

           @return: A L{Bounds2Tuple}C{(latlonSW, latlonNE)} of B{C{LatLon}}s
                    or a L{Bounds4Tuple}C{(latS, lonW, latN, lonE)} if
                    C{B{LatLon} is None},
        '''
        r = self._bounds
        return r if LatLon is None else \
           _2bounds(LatLon, LatLon_kwds, *r, name=self.name)

    def _distanceTo(self, func_, other, **kwds):
        '''(INTERNAL) Helper for distances, see C{.formy._distanceTo*}.
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
        return _GH.Sizes[n].radius

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

    def equirectangularTo(self, other, radius=R_M, **adjust_limit_wrap):
        '''Approximate the distance between this and an other geohash
           using function L{pygeodesy.equirectangular}.

           @arg other: The other geohash (L{Geohash}, C{LatLon} or C{str}).
           @kwarg radius: Mean earth radius, ellipsoid or datum (C{meter},
                          L{Ellipsoid}, L{Ellipsoid2}, L{Datum} or L{a_f2Tuple})
                          or C{None}, see function L{pygeodesy.equirectangular}.
           @kwarg adjust_limit_wrap: Optional keyword arguments for function
                         L{pygeodesy.equirectangular4}, overriding defaults
                         C{B{adjust}=False, B{limit}=None} and C{B{wrap}=False}.

           @return: Distance (C{meter}, same units as B{C{radius}} or the ellipsoid
                    or datum axes or C{radians I{squared}} if B{C{radius} is None}
                    or C{0}).

           @raise TypeError: The B{C{other}} is not a L{Geohash}, C{LatLon} or
                             C{str} or invalid B{C{radius}}.

           @see: U{Local, flat earth approximation
                 <https://www.EdWilliams.org/avform.htm#flat>}, functions
        '''
        lls  =  self.latlon + _2Geohash(other).latlon
        kwds = _xkwds(adjust_limit_wrap, adjust=False, limit=None, wrap=False)
        return _formy.equirectangular( *lls, radius=radius, **kwds) if radius else \
               _formy.equirectangular4(*lls, **kwds).distance2

    def euclideanTo(self, other, **radius_adjust_wrap):
        '''Approximate the distance between this and an other geohash using
           function L{pygeodesy.euclidean}.

           @arg other: The other geohash (L{Geohash}, C{LatLon} or C{str}).
           @kwarg radius_adjust_wrap: Optional keyword arguments for function
                                      L{pygeodesy.euclidean}.

           @return: Distance (C{meter}, same units as B{C{radius}} or the
                    ellipsoid or datum axes).

           @raise TypeError: The B{C{other}} is not a L{Geohash}, C{LatLon}
                             or C{str} or invalid B{C{radius}}.
        '''
        return self._distanceTo(_formy.euclidean, other, **radius_adjust_wrap)

    def haversineTo(self, other, **radius_wrap):
        '''Compute the distance between this and an other geohash using
           the L{pygeodesy.haversine} formula.

           @arg other: The other geohash (L{Geohash}, C{LatLon} or C{str}).
           @kwarg radius_wrap: Optional keyword arguments for function
                               L{pygeodesy.haversine}.

           @return: Distance (C{meter}, same units as B{C{radius}} or the
                    ellipsoid or datum axes).

           @raise TypeError: The B{C{other}} is not a L{Geohash}, C{LatLon}
                             or C{str} or invalid B{C{radius}}.
        '''
        return self._distanceTo(_formy.haversine, other, **radius_wrap)

    @Property_RO
    def latlon(self):
        '''Get the lat- and longitude of (the approximate center of)
           this geohash as a L{LatLon2Tuple}C{(lat, lon)} in C{degrees}.
        '''
        lat, lon = self._latlon or _2center(self.bounds())
        return LatLon2Tuple(lat, lon, name=self.name)

    @Property_RO
    def neighbors(self):
        '''Get all 8 adjacent cells as a L{Neighbors8Dict}C{(N, NE,
           E, SE, S, SW, W, NW)} of L{Geohash}es.
        '''
        return Neighbors8Dict(N=self.N, NE=self.NE, E=self.E, SE=self.SE,
                              S=self.S, SW=self.SW, W=self.W, NW=self.NW,
                              name=self.name)

    @Property_RO
    def philam(self):
        '''Get the lat- and longitude of (the approximate center of)
           this geohash as a L{PhiLam2Tuple}C{(phi, lam)} in C{radians}.
        '''
        return PhiLam2Tuple(map2(radians, self.latlon), name=self.name)  # *map2

    @Property_RO
    def precision(self):
        '''Get this geohash's precision (C{int}).
        '''
        return len(self)

    @deprecated_property_RO
    def sizes(self):
        '''DEPRECATED on 2024.07.28, use property C{Geohash.sizes3}.'''
        t = self.sizes3
        return _GH._LatLon2Tuple(t.height, t.width, name=t.name)

    @Property_RO
    def sizes3(self):
        '''Get the lat-, longitudinal and radial size of this cell as
           a L{Sizes3Tuple}C{(height, width, radius)}, all C{meter}.
        '''
        z = _GH.Sizes
        n =  min(len(z) - 1, max(self.precision, 1))
        return Sizes3Tuple(z[n], name=self.name)

    def toLatLon(self, LatLon=None, **LatLon_kwds):
        '''Return (the approximate center of) this geohash cell
           as an instance of the supplied C{LatLon} class.

           @arg LatLon: Class to use (C{LatLon}) or C{None}.
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                               arguments, ignored if C{B{LatLon} is None}.

           @return: This geohash location (B{C{LatLon}}) or if C{B{LatLon}
                    is None}, a L{LatLon2Tuple}C{(lat, lon)}.

           @raise TypeError: Invalid B{C{LatLon}} or B{C{LatLon_kwds}}.
        '''
        return self.latlon if LatLon is None else _xnamed(LatLon(
              *self.latlon, **LatLon_kwds), self.name)

    def vincentysTo(self, other, **radius_wrap):
        '''Compute the distance between this and an other geohash using
           the L{pygeodesy.vincentys} formula.

           @arg other: The other geohash (L{Geohash}, C{LatLon} or C{str}).
           @kwarg radius_wrap: Optional keyword arguments for function
                               L{pygeodesy.vincentys}.

           @return: Distance (C{meter}, same units as B{C{radius}} or the
                    ellipsoid or datum axes).

           @raise TypeError: The B{C{other}} is not a L{Geohash}, C{LatLon}
                             or C{str} or invalid B{C{radius}}.
        '''
        return self._distanceTo(_formy.vincentys, other, **radius_wrap)

    @Property_RO
    def E(self):
        '''Get the cell East of this (L{Geohash}).
        '''
        return self.adjacent(_E_)

    @Property_RO
    def N(self):
        '''Get the cell North of this (L{Geohash}).
        '''
        return self.adjacent(_N_)

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
    def S(self):
        '''Get the cell South of this (L{Geohash}).
        '''
        return self.adjacent(_S_)

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

    @Property_RO
    def W(self):
        '''Get the cell West of this (L{Geohash}).
        '''
        return self.adjacent(_W_)


class Geohashed(object):
    '''A cache of en- and decoded geohashes of one precision.
    '''
    _nn = None,  # 1-tuple

    def __init__(self, precision, ndigits=None):
        '''New L{Geohashed} cache.

           @arg precision: The geohash encoded length (C{int}, 1..12).
           @kwarg ndigits: Optional number of digits to round C{lat}
                           and C{lon} to cache keys (C{int}, typically
                           C{B{ndigits}=B{precision}}) or C{None} for
                           no rounding.
        '''
        self._p = _2Precision(precision)
        if ndigits is None:
            self._ab2 = self._ab2float
        else:
            self._ab2 = self._ab2round
            n         = Int(ndigits=ndigits)
            self._nn  = n, n
        self.clear()

    def __len__(self):
        '''Return the number of I{unigue} geohashes (C{int}).
        '''
        d = self._d
        d = set(d.keys())
        n = len(d)
        for e in self._e.values():
            e  = set(e.values())
            n += len(e - d)
        return n

    def _ab2(self, *ll):  # overwritten
        '''(INTERNAL) Make encoded keys C{a, b}.
        '''
        return ll

    def _ab2float(self, *ll):
        '''(INTERNAL) Make encoded keys C{a, b}.
        '''
        return map(float, ll)

    def _ab2round(self, *ll):
        '''(INTERNAL) Make encoded keys C{a, b}.
        '''
        return map(round, ll, self._nn)

    def clear(self):
        '''Clear the C{en-} and C{decoded} cache.
        '''
        self._e = {}
        self._d = {}

    def decoded(self, geohash, encoded=False):
        '''Get and cache the C{(lat, lon)} for C{geohash}, see L{decode<pygeodesy.geohash.decode>}.

           @kwarg encoded: If C{True}, cache the result as C{encoded}.

           @return: The C{(lat, lon}) pair for C{geohash}.
        '''
        try:
            ll = self._d[geohash]
        except KeyError:
            self._d[geohash] = ll = _GH.decode2(geohash)
        if encoded:
            a, b = self._ab2(*ll)
            try:
                _ = self._e[b][a]
            except KeyError:
                self._e.setdefault(b, {})[a] = geohash
        return ll

    def encoded(self, lat, lon, decoded=False):
        '''Get and cache the C{geohash} for C{(lat, lon)}, see L{encode<pygeodesy.geohash.encode>}.

           @kwarg decoded: If C{True}, cache the result as C{decoded}.

           @return: The C{geohash} for pair C{(lat, lon}).
        '''
        lat, lon = ll = _2fll(lat, lon)
        a, b = self._ab2(*ll)
        try:
            gh = self._e[b][a]
        except KeyError:
            gh = _GH.encode(lat, lon, self._p, 0)
            self._e.setdefault(b, {})[a] = gh
        if decoded and gh not in self._d:
            self._d[gh] = ll
        return gh

    @property_RO
    def len2(self):
        '''Return 2-tuple C{(lencoded, ldecoded)} with the C{len}gths of the
           C{en-} and C{decoded} cache.
        '''
        return sum(len(e) for e in self._e.values()), len(self._d)

    @Property_RO
    def ndigits(self):
        '''Get the rounding (C{int} or C{None}).
        '''
        return self._nn[0]

    @Property_RO
    def precision(self):
        '''Get the C{precision} (C{int}).
        '''
        return self._p


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


class Resolutions2Tuple(_NamedTuple):
    '''2-Tuple C{(res1, res2)} with the primary I{(longitudinal)} and
       secondary I{(latitudinal)} resolution, both in C{degrees}.
    '''
    _Names_ = ('res1',   'res2')
    _Units_ = ( Degrees_, Degrees_)

    @property_RO
    def lat(self):
        '''Get the secondary, latitudinal resolution (C{degrees}).
        '''
        return self[1]

    @property_RO
    def lon(self):
        '''Get the primary, longitudinal resolution (C{degrees}).
        '''
        return self[0]


class Sizes3Tuple(_NamedTuple):
    '''3-Tuple C{(height, width, radius)} with latitudinal C{height},
       longitudinal C{width} and area C{radius}, all in C{meter}.
    '''
    _Names_ = (_height_, _width_, _radius_)
    _Units_ = ( Meter,    Meter,   Meter)


def bounds(geohash, LatLon=None, **LatLon_kwds):
    '''Returns the lower-left SW and upper-right NE corners of a geohash.

       @arg geohash: To be "bound" (L{Geohash}).
       @kwarg LatLon: Optional class to return the bounds (C{LatLon}) or C{None}.
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword arguments,
                           ignored if C{B{LatLon} is None}.

       @return: A L{Bounds2Tuple}C{(latlonSW, latlonNE)}, each a B{C{LatLon}}
                or if C{B{LatLon} is None}, a L{Bounds4Tuple}C{(latS, lonW,
                latN, lonE)}.

       @raise TypeError: The B{C{geohash}} is not a L{Geohash}, C{LatLon} or
                         C{str} or invalid B{C{LatLon}} or invalid B{C{LatLon_kwds}}.

       @raise GeohashError: Invalid or C{null} B{C{geohash}}.
    '''
    swne = _GH.swne4(geohash)
    return _2bounds(LatLon, LatLon_kwds, *swne,
                            name=nameof(geohash))  # _or_nameof=geohash


def decode(geohash):
    '''Decode a geohash to lat-/longitude of the (approximate
       centre of) geohash cell to reasonable precision.

       @arg geohash: To be decoded (L{Geohash}).

       @return: 2-Tuple C{(latStr, lonStr)}, both C{str}.

       @raise TypeError: The B{C{geohash}} is not a L{Geohash},
                         C{LatLon} or C{str}.

       @raise GeohashError: Invalid or null B{C{geohash}}.
    '''
    # round to near centre without excessive precision to
    # ⌊2-log10(Δ°)⌋ decimal places, strip trailing zeros
    swne = _GH.swne4(geohash)
    return _2latlon(*swne, fstr=_MODS.streprs.fstr)


def decode2(geohash, LatLon=None, **LatLon_kwds):
    '''Decode a geohash to lat-/longitude of the (approximate center
       of) geohash cell to reasonable precision.

       @arg geohash: To be decoded (L{Geohash}).
       @kwarg LatLon: Optional class to return the location (C{LatLon})
                      or C{None}.
       @kwarg LatLon_kwds: Optional, addtional B{C{LatLon}} keyword
                           arguments, ignored if C{B{LatLon} is None}.

       @return: L{LatLon2Tuple}C{(lat, lon)}, both C{degrees} if
                C{B{LatLon} is None}, otherwise a B{C{LatLon}} instance.

       @raise TypeError: The B{C{geohash}} is not a L{Geohash},
                         C{LatLon} or C{str}.

       @raise GeohashError: Invalid or null B{C{geohash}}.
    '''
    ll = _GH.decode2(geohash)
    r  =  LatLon2Tuple(ll) if LatLon is None else \
          LatLon(     *ll,  **LatLon_kwds)
    return _xnamed(r, name__=decode2)


@deprecated_function
def decode_error(geohash):
    '''DEPRECATED on 2024.07.28, use L{geohash.decode_error2}.'''
    return decode_error2(geohash)


def decode_error2(geohash):
    '''Return the lat- and longitude decoding error for a geohash.

       @arg geohash: To be decoded (L{Geohash}).

       @return: A L{LatLon2Tuple}C{(lat, lon)} with the lat- and
                longitudinal errors in (C{degrees}).

       @raise TypeError: The B{C{geohash}} is not a L{Geohash},
                         C{LatLon} or C{str}.

       @raise GeohashError: Invalid or null B{C{geohash}}.
    '''
    s, w, n, e = _GH.swne4(geohash)
    return LatLon2Tuple((n - s) * _0_5,  # lat error
                        (e - w) * _0_5)  # lon error


def distance_(geohash1, geohash2):
    '''Estimate the distance between two geohash (from the cell sizes).

       @arg geohash1: First geohash (L{Geohash}, C{LatLon} or C{str}).
       @arg geohash2: Second geohash (L{Geohash}, C{LatLon} or C{str}).

       @return: Approximate distance (C{meter}).

       @raise TypeError: If B{C{geohash1}} or B{C{geohash2}} is not a
                         L{Geohash}, C{LatLon} or C{str}.
    '''
    return _2Geohash(geohash1).distanceTo(geohash2)


@deprecated_function
def distance1(geohash1, geohash2):
    '''DEPRECATED, use L{geohash.distance_}.'''
    return distance_(geohash1, geohash2)


@deprecated_function
def distance2(geohash1, geohash2):
    '''DEPRECATED, use L{geohash.equirectangular4}.'''
    return equirectangular4(geohash1, geohash2)


@deprecated_function
def distance3(geohash1, geohash2):
    '''DEPRECATED, use L{geohash.haversine_}.'''
    return haversine_(geohash1, geohash2)


def encode(lat, lon, precision=None, eps=EPS):
    '''Encode a lat-/longitude as a C{geohash}, either to the specified
       precision or if not provided, to an inferred precision.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).
       @kwarg precision: Optional, the desired geohash length (C{int}
                         1..12) or C{None} or C{0} for inferred.
       @kwarg eps: Optional inference tolerance (C{degrees}), ignored
                   if B{C{precision}} is not C{None} or C{0}.

       @return: The C{geohash} (C{str}).

       @raise GeohashError: Invalid B{C{lat}}, B{C{lon}} or B{C{precision}}.
    '''
    lat, lon = _2fll(lat, lon)

    if precision:
        p, e = _2Precision(precision), 0
    else:  # infer precision by refining geohash
        p, e = _MaxPrec, max(eps, EPS)
    return _GH.encode(lat, lon, p, e)


def equirectangular4(geohash1, geohash2, radius=R_M):
    '''Approximate the distance between two geohashes using the
       L{pygeodesy.equirectangular} formula.

       @arg geohash1: First geohash (L{Geohash}, C{LatLon} or C{str}).
       @arg geohash2: Second geohash (L{Geohash}, C{LatLon} or C{str}).
       @kwarg radius: Mean earth radius (C{meter}) or C{None}, see method
                      L{Geohash.equirectangularTo}.

       @return: Approximate distance (C{meter}, same units as B{C{radius}}),
                see method L{Geohash.equirectangularTo}.

       @raise TypeError: If B{C{geohash1}} or B{C{geohash2}} is not a
                         L{Geohash}, C{LatLon} or C{str}.
    '''
    return _2Geohash(geohash1).equirectangularTo(geohash2, radius=radius)


def euclidean_(geohash1, geohash2, **radius_adjust_wrap):
    '''Approximate the distance between two geohashes using the
       L{pygeodesy.euclidean} formula.

       @arg geohash1: First geohash (L{Geohash}, C{LatLon} or C{str}).
       @arg geohash2: Second geohash (L{Geohash}, C{LatLon} or C{str}).
       @kwarg radius_adjust_wrap: Optional keyword arguments for function
                                  L{pygeodesy.euclidean}.

       @return: Approximate distance (C{meter}, same units as B{C{radius}}).

       @raise TypeError: If B{C{geohash1}} or B{C{geohash2}} is not a
                         L{Geohash}, C{LatLon} or C{str}.
    '''
    return _2Geohash(geohash1).euclideanTo(geohash2, **radius_adjust_wrap)


def haversine_(geohash1, geohash2, **radius_wrap):
    '''Compute the great-circle distance between two geohashes
       using the L{pygeodesy.haversine} formula.

       @arg geohash1: First geohash (L{Geohash}, C{LatLon} or C{str}).
       @arg geohash2: Second geohash (L{Geohash}, C{LatLon} or C{str}).
       @kwarg radius_wrap: Optional keyword arguments for function
                           L{pygeodesy.haversine}.

       @return: Great-circle distance (C{meter}, same units as
                B{C{radius}}).

       @raise TypeError: If B{C{geohash1}} or B{C{geohash2}} is
                         not a L{Geohash}, C{LatLon} or C{str}.
    '''
    return _2Geohash(geohash1).haversineTo(geohash2, **radius_wrap)


def neighbors(geohash):
    '''Return the L{Geohash}es for all 8 adjacent cells.

       @arg geohash: Cell for which neighbors are requested
                     (L{Geohash} or C{str}).

       @return: A L{Neighbors8Dict}C{(N, NE, E, SE, S, SW, W, NW)}
                of L{Geohash}es.

       @raise TypeError: The B{C{geohash}} is not a L{Geohash},
                         C{LatLon} or C{str}.
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
             <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Geohash.html>}.
    '''
    r = Degrees_(res1=res1, low=_0_0, Error=GeohashError)
    N = res2 is None
    t = r, (r if N else Degrees_(res2=res2, low=_0_0, Error=GeohashError))
    for p in range(1, _MaxPrec):
        if resolution2(p, (None if N else p)) <= t:
            return p
    return _MaxPrec


def resolution2(prec1, prec2=None):
    '''Determine the (geographic) resolutions of given L{Geohash}
       precisions.

       @arg prec1: The given primary I{(longitudinal)} precision
                   (C{int} 1..12).
       @kwarg prec2: Optional, secondary I{(latitudinal)} precision
                     (C{int} 1..12).

       @return: L{Resolutions2Tuple}C{(res1, res2)} with the
                (geographic) resolutions C{degrees}, where C{res2
                B{is} res1} if no B{C{prec2}} is given.

       @raise GeohashError: Invalid B{C{prec1}} or B{C{prec2}}.

       @see: I{Karney}'s C++ class U{Geohash
             <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Geohash.html>}.
    '''
    res1, res2 = _360_0, _180_0  # note ... lon, lat!

    if prec1:
        p = 5 * max(0, min(Int(prec1=prec1, Error=GeohashError), _MaxPrec))
        res1 = res2 = ldexp(res1, -(p - p // 2))

    if prec2:
        p = 5 * max(0, min(Int(prec2=prec2, Error=GeohashError), _MaxPrec))
        res2 = ldexp(res2, -(p // 2))

    return Resolutions2Tuple(res1, res2)


@deprecated_function
def sizes(geohash):
    '''DEPRECATED on 2024.07.28, use function L{pygeodesy.geohash.sizes3}.'''
    t = sizes3(geohash)
    return _GH._LatLon2Tuple(t.height, t.width, name=t.name)


def sizes3(geohash):
    '''Return the lat-, longitudinal and radial size of this L{Geohash} cell.

       @arg geohash: Cell for which size are required (L{Geohash} or C{str}).

       @return: A L{Sizes3Tuple}C{(height, width, radius)}, all C{meter}.

       @raise TypeError: The B{C{geohash}} is not a L{Geohash}, C{LatLon} or C{str}.
    '''
    return _2Geohash(geohash).sizes3


def vincentys_(geohash1, geohash2, **radius_wrap):
    '''Compute the distance between two geohashes using the
       L{pygeodesy.vincentys} formula.

       @arg geohash1: First geohash (L{Geohash}, C{LatLon} or C{str}).
       @arg geohash2: Second geohash (L{Geohash}, C{LatLon} or C{str}).
       @kwarg radius_wrap: Optional keyword arguments for function
                           L{pygeodesy.vincentys}.

       @return: Distance (C{meter}, same units as B{C{radius}}).

       @raise TypeError: If B{C{geohash1}} or B{C{geohash2}} is not a
                         L{Geohash}, C{LatLon} or C{str}.
    '''
    return _2Geohash(geohash1).vincentysTo(geohash2, **radius_wrap)


__all__ += _ALL_DOCS(bounds,  # functions
                     decode, decode2, decode_error2, distance_,
                     encode, equirectangular4, euclidean_, haversine_,
                     neighbors, precision, resolution2, sizes3, vincentys_,
                     decode_error, sizes)  # DEPRECATED

if __name__ == '__main__':

    from pygeodesy.internals import printf, _versions
    from timeit import timeit

    for f, p in (('encode', _MaxPrec), ('infer', None)):

        def _t(prec=p):
            i = 0
            for lat in range(-90, 90, 3):
                for lon in range(-180, 180, 7):
                    _ = encode(lat, lon, prec)
                    i += 1
            return i

        i = _t()  # prime
        n =  10
        t =  timeit(_t, number=n) / (i * n)
        printf('%s %.3f usec, %s', f, t * 1e6, _versions())

# % python3.12 -m pygeodesy.geohash
# encode 9.666 usec, pygeodesy 24.8.4 Python 3.13.0b4 64bit arm64 macOS 14.5
# infer 13.735 usec, pygeodesy 24.8.4 Python 3.13.0b4 64bit arm64 macOS 14.5
# or about 6.88 and 79.75 times faster than pygeodesy 24.7.24 and older:
# encode 66.524 usec, pygeodesy 24.7.24 Python 3.12.4 64bit arm64 macOS 14.5
# infer 1095.386 usec, pygeodesy 24.7.24 Python 3.12.4 64bit arm64 macOS 14.5

# **) MIT License
#
# Copyright (C) 2016-2024 -- mrJean1 at Gmail -- All Rights Reserved.
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
