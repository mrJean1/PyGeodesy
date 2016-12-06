
# -*- coding: utf-8 -*-

# Python implementation of UTM / WGS-84 conversion functions using
# an ellipsoidal earth model.  Transcribed from JavaScript originals
# by (C) Chris Veness 2011-2016 published under the same MIT Licence**,
# see <http://www.movable-type.co.uk/scripts/latlong-utm-mgrs.html>
# and <http://www.movable-type.co.uk/scripts/geodesy/docs/module-utm.html>

from math import asinh, atan, atanh, atan2, cos, cosh, \
                 hypot, sin, sinh, tan, tanh
from operator import mul
from bases import _Base
from datum import Datums
from dms   import S_DEG
from utils import EPS, degrees, degrees90, degrees180, \
                  fdot3, fStr, hypot1, isscalar, len2, map2, \
                  radians, wrap90, wrap180

# The Universal Transverse Mercator (UTM) system is a 2-dimensional
# Cartesian coordinate system providing locations on the surface of
# the Earth.

# UTM is a set of 60 transverse Mercator projections, normally based
# on the WGS-84 ellipsoid.  Within each zone, coordinates are
# represented as eastings and northings, measured in metres.

# This method based on Karney 2011 'Transverse Mercator with an
# accuracy of a few nanometers', building on Krüger 1912 'Konforme
# Abbildung des Erdellipsoids in der Ebene'.

# References <https://arxiv.org/pdf/1002.1417v3.pdf>,
# <http://bib.gfz-potsdam.de/pub/digi/krueger2.pdf>,
# <http://henrik-seidel.gmxhome.de/gausskrueger.pdf> and
# <http://wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system>.

# all public contants, classes and functions
__all__ = ('Utm',  # classes
           'parseUTM', 'toUtm')  # functions
__version__ = '16.12.06'

# Latitude bands C..X of 8° each, covering 80°S to 84°N
_Bands         = 'CDEFGHJKLMNPQRSTUVWXX'  # X repeated for 80-84°N
_FalseEasting  =   500e3  # meter
_FalseNorthing = 10000e3  # meter
_K0            = 0.9996   # UTM scale on the central meridian


class _Ks(object):
    # Alpha or Beta Krüger series summations for x, y,
    # p and q and caching the cos, sin, cosh, and sinh
    # for the given x and y angles (in radians).

    def __init__(self, AB, x, y):
        # AB is a 1-origin, 6th-order
        # Krüger Alpha or Beta series
        n, j2 = len2(range(2, len(AB) * 2, 2))

        self._ab = AB[1:]  # 0-origin
        self._pq = map2(mul, j2, self._ab)

        assert len(self._ab) == len(self._pq) == n

        x2 = map2(mul, j2, (x,) * n)
        self._chx = map2(cosh, x2)
        self._shx = map2(sinh, x2)

        assert len(x2) == len(self._chx) == len(self._shx) == n

        y2 = map2(mul, j2, (y,) * n)
        self._cy = map2(cos, y2)
        self._sy = map2(sin, y2)

        assert len(y2) == len(self._cy) == len(self._sy) == n

    def xs(self, x):
        return fdot3(self._ab, self._cy, self._shx, start=x)

    def ys(self, y):
        return fdot3(self._ab, self._sy, self._chx, start=y)

    def ps(self, p):
        return fdot3(self._pq, self._cy, self._chx, start=p)

    def qs(self, q):
        return fdot3(self._pq, self._sy, self._shx, start=q)


def _toZBL(zone, band, mgrs=False):  # used by mgrs.Mgrs
    # check and return zone, Band and band latitude
    try:
        if isscalar(zone) or zone.isdigit():
            z, B, x = int(zone), str(band), band
        else:
            z, B, x = int(zone[:-1] or 0), zone[-1:], zone

        if 1 > z or z > 60:
            raise ValueError

    except (AttributeError, TypeError, ValueError):
        raise ValueError('%s invalid: %r' % ('zone', zone))

    b = None
    if B:
        b = _Bands.find(B)
        if b < 0:
            raise ValueError('%s invalid: %r' % ('band', x))
        b = (b - 10) * 8
    elif mgrs:
        raise ValueError('%s missing' % ('band',))

    return z, B, b


class Utm(_Base):
    '''UTM coordinate.
    '''
    _band     = ''
    _converge = None
    _datum    = Datums.WGS84
    _easting  = 0
    _hemi     = ''
    _latlon   = None  # also set by ellipsoidalBase._LatLonHeightDatumBase.toUtm.
    _mgrs     = None
    _northing = 0
    _scale    = None
    _zone     = 0

    def __init__(self, zone, hemisphere, easting, northing, band='',
                       datum=Datums.WGS84, convergence=None, scale=None):
        '''Creates a UTM coordinate.

           @param {number|string} zone - UTM 6° longitudinal zone (1..60
                                         covering 180°W..180°E) as number
                                         or 00B zone and band as string.
           @param {string} hemisphere - N for the northern or S for the
                                        southern hemisphere.
           @param {meter} easting - Easting from false easting (-500km
                                    from central meridian).
           @param {meter} northing - Northing from equator (N) or from
                                     false northing -10,000km (S).
           @param {string} [band=''] - Latitudinal band (C..X), optional.
           @param {Datum} [datum=Datums.WGS84] - This coordinate's datum.
           @param {degrees} [convergence=None] - Meridian convergence
                                                 (bearing of grid north,
                                                 clockwise from true
                                                 North) in degrees.
           @param {number} [scale=None] - Grid scale factor.

           @throws {ValueError} Invalid UTM coordinate.

           @example
           import utm
           g = utm.Utm(31, 'N', 448251, 5411932)
        '''
        self._zone, B, _ = _toZBL(zone, band)

        h = str(hemisphere)[:1]
        if h not in 'NnSs':
            raise ValueError('%s invalid: %r' % ('hemisphere', hemisphere))

        e, n = float(easting), float(northing)
        # check easting/northing (with 40km overlap
        # between zones) - is this worthwhile?
        if 120e3 > e or e > 880e3:
            raise ValueError('%s invalid: %r' % ('easting', easting))
        if 0 > n or n > _FalseNorthing:
            raise ValueError('%s invalid: %r' % ('northing', northing))

        self._hemi     = h.upper()
        self._easting  = e
        self._northing = n
        if self._band != B:
            self._band = B
        if self._datum != datum:
            self._datum = datum
        if self._converge != convergence:
            self._converge = convergence
        if self._scale != scale:
            self._scale = scale

    @property
    def band(self):
        '''Return latitudinal band (C..X) if available.'''
        return self._band

    @property
    def convergence(self):
        '''Return convergence in degrees or None.'''
        return self._converge

    @property
    def datum(self):
        '''Return the datum.'''
        return self._datum

    @property
    def easting(self):
        '''Return easting in meter.'''
        return self._easting

    @property
    def hemisphere(self):
        '''Return hemisphere (N|S).'''
        return self._hemi

    @property
    def northing(self):
        '''Return northing in meter.'''
        return self._northing

    def parseUTM(self, strUTM):
        '''Parse a string to a UTM coordinate.

           For details, see function parseUTM
           in this module utm.
        '''
        return parseUTM(strUTM, datum=self.datum)

    @property
    def scale(self):
        '''Return convergence in degrees or None.'''
        return self._scale

    def toLatLon(self, LatLon):
        '''Converts UTM coordinate to an ellipsoidal lat-/longitude.

           @param {LatLon} LatLon - Ellipsoidal LatLon class to use.

           @returns {LatLon} Lat-/longitude of this UTM coordinate.

           @example
           g = Utm(31, 'N', 448251.795, 5411932.678)
           from geodesy import ellipsoidalVincenty as eV
           ll = g.toLatLon(eV.LatLon)  # 48°51′29.52″N, 002°17′40.20″E
        '''
        if self._latlon and self._latlon.__class__ is LatLon \
                        and self._latlon.datum == self._datum:
            return self._latlon  # set below

        E = self._datum.ellipsoid  # XXX vs LatLon.datum.ellipsoid

        x = self._easting - _FalseEasting  # relative to central meridian
        y = self._northing
        if self._hemi == 'S':  # relative to equator
            y -= _FalseNorthing

        # from Karney 2011 Eq 15-22, 36
        A0 = _K0 * E.A
        x /= A0  # η eta
        y /= A0  # ξ ksi

        B6 = _Ks(E.Beta6, x, y)  # 6th-order Krüger series, 1-origin
        y = -B6.ys(-y)  # ξ'
        x = -B6.xs(-x)  # η'

        shx = sinh(x)
        cy, sy = cos(y), sin(y)

        H = hypot(shx, cy)

        T = t0 = sy / H
        q = 1.0 / E.e12
        d = 1
        # note, a relatively large convergence test as d
        # toggles on ±1.12e-16 for eg 31 N 400000 5000000
        while abs(d) > EPS:  # 1e-12
            h = hypot1(T)
            s = sinh(E.e * atanh(E.e * T / h))
            t = T * hypot1(s) - s * h
            d = (t0 - t) / hypot1(t) * (q + T * T) / h
            T += d

        a = atan(T)  # lat
        b = atan2(shx, cy) + radians(self._zone * 6 - 183)  # lon of central meridian
        ll = LatLon(degrees180(a), degrees90(b), datum=self._datum)

        # convergence: Karney 2011 Eq 26, 27
        p = -B6.ps(-1)
        q =  B6.qs(0)
        ll.convergence = degrees(atan(tan(y) * tanh(x)) + atan2(q, p))

        # scale: Karney 2011 Eq 28
        ll.scale = E.e2s2(sin(a)) * hypot1(T) * H * (A0 / E.a / hypot(p, q))

        self._latlon = ll
        return ll

    def toMgrs(self):
        '''Convert this UTM coordinate to an MGRS grid reference.

           See function toMgrs in module mgrs for more details.

           @returns {Mgrs} The MGRS grid reference.
        '''
        if self._mgrs is None:
            from mgrs import toMgrs  # PYCHOK recursive import
            self._mgrs = toMgrs(self)
        return self._mgrs

    def toStr(self, prec=0, sep=' ', B=False, cs=False):  # PYCHOK expected
        '''Returns a string representation of this UTM coordinate.

           To distinguish from MGRS grid zone designators, a
           space is left between the zone and the hemisphere.

           Note that UTM coordinates are rounded, not truncated
           (unlike MGRS grid references).

           @param {number} [prec=0] - Number of decimal, unstripped.
           @param {string} [sep=' '] - Separator to join.
           @param {bool} [B=False] - Include latitudinal band.
           @param {bool} [cs=False] - Include grid convergence and
                                      scale factor.

           @returns {string} UTM as string "00 N|S meter meter"
                             plus "degrees float" if cs is True.

           @example
           u = Utm(3, 'N', 448251, 5411932.0001)
           u.toStr(4)  # 03 N 448251.0 5411932.0001
           u.toStr(sep=', ')  # 03 N, 448251, 5411932
        '''
        b = self._band if B else ''
        t = ['%02d%s %s' % (self._zone, b, self._hemi),
             fStr(self._easting, prec=prec),
             fStr(self._northing, prec=prec)]
        if cs:
            t += ['n/a' if self._converge is None else
                      fStr(self._converge, prec=8, fmt='%+013.*f') + S_DEG,
                  'n/a' if self._scale is None else
                      fStr(self._scale, prec=8)]
        return sep.join(t)

    def toStr2(self, prec=0, fmt='[%s]', sep=', ', B=False, cs=False):  # PYCHOK expected
        '''Returns a string representation of this UTM coordinate.

           Note that UTM coordinates are rounded, not truncated
           (unlike MGRS grid references).

           @param {number} [prec=0] - Number of decimals, unstripped.
           @param {string} [fmt='[%s]'] - Enclosing backets format.
           @param {string} [sep=', '] - Separator between name:values.
           @param {bool} [B=False] - Include latitudinal band.
           @param {bool} [cs=False] - Include grid convergence and
                                      scale factor.

           @returns {string} This Utm as "[Z:00, H:N|S, E:meter, N:meter]"
                             string plus "C:degrees, S:float" if cs is True.
        '''
        t = self.toStr(prec=prec, sep=' ', B=B, cs=cs).split()
        k = 'ZHENCS' if cs else 'ZHEN'
        return fmt % (sep.join('%s:%s' % t for t in zip(k, t)),)

    @property
    def zone(self):
        '''Return longitudal zone (1..60).'''
        return self._zone


def parseUTM(strUTM, datum=Datums.WGS84):
    '''Parse a string representing a UTM coordinate, consisting of
       zone, hemisphere, easting and northing.

       @param {string} strUTM - A UTM coordinate.
       @param {Datum} [datum=Datums.WGS84] - The datum to use.

       @returns {Utm} The UTM coordinate.

       @throws {ValueError} Invalid strUTM.

       @example
       u = parseUTM('31 N 448251 5411932')
       u.toStr2()  # [Z:31, H:N, E:448251, N:5411932]
       u = parseUTM('31 N 448251.8 5411932.7')
       u.toStr()  # 31 N 448252 5411933
    '''
    u = strUTM.strip().replace(',', ' ').split()
    try:
        if len(u) != 4:
            raise ValueError  # caught below
        z, h = u[:2]
        if z.isdigit():
            z = int(z)
        else:
            z = z.upper()
        e, n = map(float, u[2:])
    except ValueError:
        raise ValueError('%s invalid: %r' % ('strUTM', strUTM))

    return Utm(z, h.upper(), e, n, datum=datum)


def toUtm(latlon, lon=None, datum=Datums.WGS84):
    '''Convert lat-/longitude location to a UTM coordinate.

       Implements Karney’s method, using 6-th order Krüger series,
       giving results accurate to 5 nm for distances up to 3900 km
       from the central meridian.

       @param {degrees|LatLon} latlon - Latitude in degrees or a
                                        LatLon instance.
       @param {degrees} [lon=None] - Longitude in degrees or None.
       @param {Datum} [datum=WGS84] - Datum for this UTM coordinate.

       @returns {Utm} The UTM coordinate.

       @throws {TypeError} If latlon is not an ellipsoidal LatLon.
       @throws {ValueError} If latitude value is missing.

       @example
       p = LatLon(48.8582, 2.2945)  # 31 N 448251.8 5411932.7
       u = toUtm(p)  # 31 N 448252 5411933
       p = LatLon(13.4125, 103.8667) # 48 N 377302.4 1483034.8
       u = toUtm(p)  # 48 N 377302 1483035
    '''
    try:
        lat, lon = latlon.lat, latlon.lon
        if not hasattr(latlon, 'toUtm'):
            raise TypeError('%s not ellipsoidal: %r' % ('latlon', latlon))
    except AttributeError:
        if lon is None:
            raise ValueError('%s invalid: %r' % ('latlon', latlon))
        lat = latlon

    lat = wrap90(lat)
    if -80 < lat < 84:  # UTM latitude range
        B = _Bands[int(lat + 80) >> 3]
    else:
        raise ValueError('%s outside UTM: %s' % ('lat', lat))

    lon = wrap180(lon)
    z = int((lon + 180) / 6) + 1  # longitudinal zone
    if lat > 55 and lon > -1:
        if lat < 64:  # southern Norway
            if lon < 3:
                z = 31
            elif lon < 12:
                z = 32
        elif 71 < lat < 84:  # Svalbard
            if lon < 9:
                z = 31
            elif lon < 21:
                z = 33
            elif lon < 33:
                z = 35
            elif lon < 42:
                z = 37

    b = radians(lon - (z * 6) + 183)  # lon off central meridian
    a = radians(lat)  # lat off equator
    h = 'S' if a < 0 else 'N'  # hemisphere

    E = datum.ellipsoid

    # easting, northing: Karney 2011 Eq 7-14, 29, 35
    cb, sb, tb = cos(b), sin(b), tan(b)

    T = tan(a)
    T12 = hypot1(T)
    S = sinh(E.e * atanh(E.e * T / T12))

    T_ = T * hypot1(S) - S * T12
    H = hypot(T_, cb)

    y = atan2(T_, cb)  # ξ' ksi
    x = asinh(sb / H)  # η' eta

    A0 = _K0 * E.A

    A6 = _Ks(E.Alpha6, x, y)  # 6th-order Krüger series, 1-origin
    y = A6.ys(y) * A0  # ξ
    x = A6.xs(x) * A0  # η

    x += _FalseEasting  # make x relative to false easting
    if y < 0:
        y += _FalseNorthing  # y relative to false northing in S

    # convergence: Karney 2011 Eq 23, 24
    p_ = A6.ps(1)
    q_ = A6.qs(0)
    c = degrees(atan(T_ / hypot1(T_) * tb) + atan2(q_, p_))

    # scale: Karney 2011 Eq 25
    k = E.e2s2(sin(a)) * T12 / H * (A0 / E.a * hypot(p_, q_))

    return Utm(z, h, x, y, band=B, datum=datum, convergence=c, scale=k)

# **) MIT License
#
# Copyright (c) 2016-2017 -- mrJean1@Gmail.com
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
