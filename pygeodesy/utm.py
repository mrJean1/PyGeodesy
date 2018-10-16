
# -*- coding: utf-8 -*-

u'''Universal Transverse Mercator (UTM) class L{Utm} and functions
L{parseUTM} and L{toUtm}.

Pure Python implementation of UTM / WGS-84 conversion functions using
an ellipsoidal earth model, transcribed from JavaScript originals by
I{(C) Chris Veness 2011-2016} published under the same MIT Licence**, see
U{UTM<http://www.Movable-Type.co.UK/scripts/latlong-utm-mgrs.html>} and
U{Module utm<http://www.Movable-Type.co.UK/scripts/geodesy/docs/module-utm.html>}.

The UTM system is a 2-dimensional cartesian coordinate system providing
locations on the surface of the earth.

UTM is a set of 60 transverse Mercator projections, normally based on
the WGS-84 ellipsoid.  Within each zone, coordinates are represented
as eastings and northings, measured in metres.

This method based on Karney U{'Transverse Mercator with an accuracy of
a few nanometers'<http://Arxiv.org/pdf/1002.1417v3.pdf>}, 2011 (building
on Krüger U{'Konforme Abbildung des Erdellipsoids in der Ebene'
<http://bib.GFZ-Potsdam.DE/pub/digi/krueger2.pdf>}, 1912).

Other references Seidel U{'Die Mathematik der Gauß-Krueger-Abbildung'
<http://Henrik-Seidel.GMXhome.DE/gausskrueger.pdf>}, 2006,
U{Transverse Mercator Projection<http://GeographicLib.SourceForge.io/tm.html>},
and U{Universal Transverse Mercator coordinate system
<http://WikiPedia.org/wiki/Universal_Transverse_Mercator_coordinate_system>}.

@newfield example: Example, Examples
'''

from bases import Based, _nameof, _xattrs, _xnamed
from datum import Datums
from dms import S_DEG, parseDMS2, RangeError
from ellipsoidalBase import LatLonEllipsoidalBase as _LLEB
from fmath import EPS, fdot3, fStr, Fsum, hypot, hypot1, \
                  isscalar, len2, map2
from utily import degrees90, degrees180, property_RO, wrap90, wrap180

from math import asinh, atan, atanh, atan2, cos, cosh, \
                 degrees, radians, sin, sinh, tan, tanh
from operator import mul

# all public contants, classes and functions
__all__ = ('Utm', 'UTMError',  # classes
           'parseUTM', 'toUtm', 'utmZoneBand2')  # functions
__version__ = '18.10.12'

# Latitude bands C..X of 8° each, covering 80°S to 84°N with X repeated
# for 80-84°N
_Bands         = 'CDEFGHJKLMNPQRSTUVWXX'  #: (INTERNAL) Latitude bands.
_FalseEasting  =   500e3  #: (INTERNAL) False (C{meter}).
_FalseNorthing = 10000e3  #: (INTERNAL) False (C{meter}).
_K0            = 0.9996   #: (INTERNAL) UTM scale central meridian.


class UTMError(ValueError):
    '''UTM parse or other error.
    '''
    pass


class _Kseries(object):
    '''(INTERNAL) Alpha or Beta Krüger series.

       Krüger series summations for I{eta}, I{ksi}, I{p} and I{q}
       while caching the cos, sin, cosh and sinh values for the
       given I{eta} and I{ksi} angles (in C{radians}).
    '''
    def __init__(self, AB, x, y):
        '''(INTERNAL) New Alpha or Beta Krüger series

           @param AB: Krüger Alpha or Beta series coefficients (4-,
                      6- or 8-tuple).
           @param x: Eta angle (C{radians}).
           @param y: Ksi angle (C{radians}).
        '''
        n, j2 = len2(range(2, len(AB) * 2 + 1, 2))

        self._ab = AB
        self._pq = map2(mul, j2, self._ab)
#       assert len(self._ab) == len(self._pq) == n

        x2 = map2(mul, j2, (x,) * n)
        self._chx = map2(cosh, x2)
        self._shx = map2(sinh, x2)
#       assert len(x2) == len(self._chx) == len(self._shx) == n

        y2 = map2(mul, j2, (y,) * n)
        self._cy = map2(cos, y2)
        self._sy = map2(sin, y2)
#       assert len(y2) == len(self._cy) == len(self._sy) == n

    def xs(self, x0):
        '''(INTERNAL) Eta summation (C{float}).
        '''
        return fdot3(self._ab, self._cy, self._shx, start=x0)

    def ys(self, y0):
        '''(INTERNAL) Ksi summation (C{float}).
        '''
        return fdot3(self._ab, self._sy, self._chx, start=y0)

    def ps(self, p0):
        '''(INTERNAL) P summation (C{float}).
        '''
        return fdot3(self._pq, self._cy, self._chx, start=p0)

    def qs(self, q0):
        '''(INTERNAL) Q summation (C{float}).
        '''
        return fdot3(self._pq, self._sy, self._shx, start=q0)


def _cmlon(zone):
    '''(INTERNAL) Central meridian longitude (C{degrees180}).
    '''
    return (zone * 6) - 183


def _toZBab4(lat, lon, cmoff=True):
    '''(INTERNAL) Return zone, Band and central lat- and longitude in radians.

       @param lat: Latitude (C{degrees}).
       @param lon: Longitude (C{degrees}).
       @keyword cmoff: Offset i{lon} from zone's central meridian.

       @return: 4-Tuple (zone, Band, a, b).
    '''
    # return zone, Band and lat- and medidian (in C{radians})
    lat = wrap90(lat)
    if -80 > lat or lat > 84:
        raise RangeError('%s outside UTM: %s' % ('lat', lat))
    B = _Bands[int(lat + 80) >> 3]

    lon = wrap180(lon)
    z = int((lon + 180) / 6) + 1  # longitudinal zone
    if B == 'X':
        x = {32: 9, 34: 21, 36: 33}.get(z, None)
        if x:  # Svalbard
            z += 1 if lon >= x else -1
    elif B == 'V' and z == 31 and lon >= 3:
        z += 1  # southern Norway

    a = radians(lat)  # lat off equator
    if cmoff:  # lon off central meridian
        lon -= _cmlon(z)
    return z, B, a, radians(lon)


def _toZBlat3(zone, band, mgrs=False):  # used by mgrs.Mgrs
    '''(INTERNAL) Check and return zone, Band and band latitude.

       @param zone: Zone number or string.
       @param band: Band letter.
       @param mgrs: Optionally, raise UTMError (C{bool}).

       @return: 3-Tuple (zone, Band, latitude).
    '''
    try:
        if isscalar(zone) or zone.isdigit():
            z, B, x = int(zone), str(band), band
        else:
            z, B, x = int(zone[:-1] or 0), zone[-1:], zone

        if 1 > z or z > 60:
            raise ValueError

    except (AttributeError, TypeError, ValueError):
        raise UTMError('%s invalid: %r' % ('zone', zone))

    b = None
    if B:
        b = _Bands.find(B)
        if b < 0:
            raise UTMError('%s invalid: %r' % ('band', x))
        b = (b << 3) - 80
    elif mgrs:
        raise UTMError('%s missing' % ('band',))

    return z, B, b


class Utm(Based):
    '''Universal Transverse Mercator (UTM) coordinate.
    '''
    _band        = ''    #: (INTERNAL) Latitude band letter ('C..X').
    _convergence = None  #: (INTERNAL) Meridian conversion (C{degrees}).
    _datum       = Datums.WGS84  #: (INTERNAL) L{Datum}.
    _easting     = 0     #: (INTERNAL) Easting from false easting (C{meter}).
    _hemi        = ''    #: (INTERNAL) Hemisphere ('N' or 'S')
    # _latlon also set by ellipsoidalBase.LatLonEllipsoidalBase.toUtm
    _latlon      = None  #: (INTERNAL) toLatLon cache.
    _mgrs        = None  #: (INTERNAL) toMgrs cache.
    _northing    = 0     #: (INTERNAL) Northing from false northing (C{meter}).
    _scale       = None  #: (INTERNAL) Grid scale factor (scalar or C{None}).
    _zone        = 0     #: (INTERNAL) Longitudinal zone (1..60).

    def __init__(self, zone, hemisphere, easting, northing, band='',
                             datum=Datums.WGS84, convergence=None,
                             scale=None, name=''):
        '''New UTM coordinate.

           @param zone: UTM 6° longitudinal zone (C{int}, 1..60 covering
                        180°W..180°E) or '00B' (C{str}) zone and band letter.
           @param hemisphere: N for the northern or S for the southern
                              hemisphere (C{str}).
           @param easting: Easting from false easting (C{meter}), -500km
                           from central meridian.
           @param northing: Northing from equator (C{meter}), N or from
                            false northing -10,000km S.
           @keyword band: Optional, latitudinal band (string, C..X).
           @keyword datum: Optional, this coordinate's datum (L{Datum}).
           @keyword convergence: Optional meridian convergence, bearing
                                 of grid North, clockwise from true
                                 North (C{degrees} or C{None}).
           @keyword scale: Optional grid scale factor (C{scalar} or C{None}).
           @keyword name: Optional name (C{str}).

           @raise RangeError: If I{easting} or I{northing} is outside
                              the valid UTM range.

           @raise UTMError: Invalid I{band}, I{hemishere} or I{zone}.

           @example:

           >>> import pygeodesy
           >>> u = pygeodesy.Utm(31, 'N', 448251, 5411932)
        '''
        if name:
            self.name = name

        self._zone, B, _ = _toZBlat3(zone, band)

        h = str(hemisphere)[:1].upper()
        if not h or h not in ('N', 'S'):
            raise UTMError('%s invalid: %r' % ('hemisphere', hemisphere))

        e, n = float(easting), float(northing)
        # check easting/northing (with 40km overlap
        # between zones) - is this worthwhile?
        if 120e3 > e or e > 880e3:
            raise RangeError('%s invalid: %r' % ('easting', easting))
        if 0 > n or n > _FalseNorthing:
            raise RangeError('%s invalid: %r' % ('northing', northing))

        self._hemi        = h
        self._easting     = e
        self._northing    = n
        self._band        = B
        self._datum       = datum
        self._convergence = convergence
        self._scale       = scale

    def _xcopy(self, *attrs):
        '''(INTERNAL) Make copy with add'l, subclass attributes.
        '''
        return _xattrs(self.classof(self.zone, self.hemisphere,
                                    self.easting, self.northing,
                                    band=self.band, datum=self.datum),
                       self, *attrs)

    @property_RO
    def band(self):
        '''Get the latitudinal band (C..X or '').
        '''
        return self._band

    @property_RO
    def convergence(self):
        '''Get the meridian convergence (C{degrees} or C{None}).
        '''
        return self._convergence

    def copy(self):
        '''Copy this UTM coordinate.

           @return: The copy (L{Utm} or subclass thereof).
        '''
        return self._xcopy()

    @property_RO
    def datum(self):
        '''Get the datum (L{Datum}).
        '''
        return self._datum

    @property_RO
    def easting(self):
        '''Get the easting (C{meter}).'''
        return self._easting

    @property_RO
    def hemisphere(self):
        '''Get the hemisphere (N|S).
        '''
        return self._hemi

    @property_RO
    def northing(self):
        '''Get the northing (C{meter}).
        '''
        return self._northing

    def parseUTM(self, strUTM):
        '''Parse a string to a UTM coordinate.

           For more details, see function L{parseUTM} in
           this module L{utm}.
        '''
        return parseUTM(strUTM, datum=self.datum)

    @property_RO
    def scale(self):
        '''Get the grid scale (C{scalar} or C{None}).
        '''
        return self._scale

    def toLatLon(self, LatLon=None):
        '''Convert this UTM coordinate to an (ellipsoidal) geodetic point.

           @keyword LatLon: Optional, ellipsoidal (sub-)class to use
                            for the point (C{LatLon}) or C{None}.

           @return: This UTM coordinate as (I{LatLon}) or 5-tuple
                    (lat, lon, datum, convergence, scale) if I{LatLon}
                    is C{None}.

           @raise TypeError: If I{LatLon} is not ellipsoidal.

           @raise UTMError: Invalid meridional radius or H-value.

           @example:

           >>> u = Utm(31, 'N', 448251.795, 5411932.678)
           >>> from pygeodesy import ellipsoidalVincenty as eV
           >>> ll = u.toLatLon(eV.LatLon)  # 48°51′29.52″N, 002°17′40.20″E
        '''
        if self._latlon:
            return self._latlon5(LatLon)

        E = self._datum.ellipsoid  # XXX vs LatLon.datum.ellipsoid

        x = self._easting - _FalseEasting  # relative to central meridian
        y = self._northing
        if self._hemi == 'S':  # relative to equator
            y -= _FalseNorthing

        # from Karney 2011 Eq 15-22, 36
        A0 = _K0 * E.A
        if A0 < EPS:
            raise UTMError('%s invalid: %r' % ('meridional', E.A))
        x /= A0  # η eta
        y /= A0  # ξ ksi

        Ks = _Kseries(E.BetaKs, x, y)  # Krüger series
        y = -Ks.ys(-y)  # ξ'
        x = -Ks.xs(-x)  # η'

        shx = sinh(x)
        cy, sy = cos(y), sin(y)

        H = hypot(shx, cy)
        if H < EPS:
            raise UTMError('%s invalid: %r' % ('H', H))

        T = t0 = sy / H  # τʹ
        q = 1.0 / E.e12
        d = 1
        sd = Fsum(T)
        # toggles on +/-1.12e-16 eg. 31 N 400000 5000000
        while abs(d) > EPS:  # 1e-12
            h = hypot1(T)
            s = sinh(E.e * atanh(E.e * T / h))
            t = T * hypot1(s) - s * h
            d = (t0 - t) / hypot1(t) * (q + T**2) / h
            T = sd.fsum_(d)  # τi

        a = atan(T)  # lat
        b = atan2(shx, cy) + radians(_cmlon(self._zone))
        ll = _LLEB(degrees90(a), degrees180(b), datum=self._datum, name=self.name)

        # convergence: Karney 2011 Eq 26, 27
        p = -Ks.ps(-1)
        q =  Ks.qs(0)
        ll._convergence = degrees(atan(tan(y) * tanh(x)) + atan2(q, p))

        # scale: Karney 2011 Eq 28
        ll._scale = E.e2s(sin(a)) * hypot1(T) * H * (A0 / E.a / hypot(p, q))

        self._latlon = ll
        return self._latlon5(LatLon)

    def _latlon5(self, LatLon):
        '''(INTERNAL) Convert cached LatLon
        '''
        ll = self._latlon
        if LatLon is None:
            return ll.lat, ll.lon, ll.datum, ll.convergence, ll.scale
        elif issubclass(LatLon, _LLEB):
            return _xnamed(_xattrs(LatLon(ll.lat, ll.lon, datum=ll.datum),
                                   ll, '_convergence', '_scale'), ll.name)
        raise TypeError('%s not ellipsoidal: %r' % ('LatLon', LatLon))

    def toMgrs(self):
        '''Convert this UTM coordinate to an MGRS grid reference.

           See function L{toMgrs} in module L{mgrs} for more details.

           @return: The MGRS grid reference (L{Mgrs}).
        '''
        if self._mgrs is None:
            from mgrs import toMgrs  # PYCHOK recursive import
            self._mgrs = toMgrs(self, name=self.name)
        return self._mgrs

    def toStr(self, prec=0, sep=' ', B=False, cs=False):  # PYCHOK expected
        '''Return a string representation of this UTM coordinate.

           To distinguish from MGRS grid zone designators, a
           space is left between the zone and the hemisphere.

           Note that UTM coordinates are rounded, not truncated
           (unlike MGRS grid references).

           @keyword prec: Optional number of decimals, unstripped (C{int}).
           @keyword sep: Optional separator to join (C{str}).
           @keyword B: Optionally, include latitudinal band (C{bool}).
           @keyword cs: Optionally, include meridian convergence and
                        grid scale factor (C{bool}).

           @return: This UTM as string "00 N|S meter meter" (C{str})
                    plus " degrees float" if I{cs} is C{True}.

           @example:

           >>> u = Utm(3, 'N', 448251, 5411932.0001)
           >>> u.toStr(4)  # 03 N 448251.0 5411932.0001
           >>> u.toStr(sep=', ')  # 03 N, 448251, 5411932
        '''
        b = self._band if B else ''
        t = ['%02d%s %s' % (self._zone, b, self._hemi),
             fStr(self._easting, prec=prec),
             fStr(self._northing, prec=prec)]
        if cs:
            t += ['n/a' if self._convergence is None else
                      fStr(self._convergence, prec=8, fmt='%+013.*f') + S_DEG,
                  'n/a' if self._scale is None else
                      fStr(self._scale, prec=8)]
        return sep.join(t)

    def toStr2(self, prec=0, fmt='[%s]', sep=', ', B=False, cs=False):  # PYCHOK expected
        '''Return a string representation of this UTM coordinate.

           Note that UTM coordinates are rounded, not truncated
           (unlike MGRS grid references).

           @keyword prec: Optional number of decimals, unstripped (C{int}).
           @keyword fmt: Optional, enclosing backets format (C{str}).
           @keyword sep: Optional separator between name:value pairs (C{str}).
           @keyword B: Optionally, include latitudinal band (C{bool}).
           @keyword cs: Optionally, include meridian convergence and
                        grid scale factor (C{bool}).

           @return: This UTM as "[Z:00, H:N|S, E:meter, N:meter]" (C{str})
                    plus "C:degrees, S:float" if I{cs} is C{True}.
        '''
        t = self.toStr(prec=prec, sep=' ', B=B, cs=cs).split()
        k = 'ZHENCS' if cs else 'ZHEN'
        return fmt % (sep.join('%s:%s' % t for t in zip(k, t)),)

    @property_RO
    def zone(self):
        '''Get the longitudinal zone (1..60).
        '''
        return self._zone


def parseUTM(strUTM, datum=Datums.WGS84, Utm=Utm, name=''):
    '''Parse a string representing a UTM coordinate, consisting
       of zone, hemisphere, easting and northing.

       @param strUTM: A UTM coordinate (C{str}).
       @keyword datum: Optional datum to use (L{Datum}).
       @keyword Utm: Optional (sub-)class to use for the UTM
                     coordinate (L{Utm}) or C{None}.
       @keyword name: Optional I{Utm} name (C{str}).

       @return: The UTM coordinate (L{Utm}) or 4-tuple (zone,
                hemisphere, easting, northing) if I{Utm} is C{None}.

       @raise UTMError: Invalid I{strUTM}.

       @example:

       >>> u = parseUTM('31 N 448251 5411932')
       >>> u.toStr2()  # [Z:31, H:N, E:448251, N:5411932]
       >>> u = parseUTM('31 N 448251.8 5411932.7')
       >>> u.toStr()  # 31 N 448252 5411933
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
        raise UTMError('%s invalid: %r' % ('strUTM', strUTM))

    return (z, h, e, n) if Utm is None else _xnamed(Utm(
            z, h, e, n, datum=datum), name)


def toUtm(latlon, lon=None, datum=None, Utm=Utm, name='', cmoff=True):
    '''Convert a lat-/longitude point to a UTM coordinate.

       @param latlon: Latitude (C{degrees}) or an (ellipsoidal)
                      geodetic C{LatLon} point.
       @keyword lon: Optional longitude (C{degrees} or C{None}).
       @keyword datum: Optional datum for this UTM coordinate,
                       overriding I{latlon}'s datum (C{Datum}).
       @keyword Utm: Optional (sub-)class to use for the UTM
                     coordinate (L{Utm}) or C{None}.
       @keyword name: Optional I{Utm} name (C{str}).
       @keyword cmoff: Offset longitude from zone's central meridian,
                       apply false easting and false northing (C{bool}).

       @return: The UTM coordinate (L{Utm}) or a 6-tuple (zone, easting,
                northing, band, convergence, scale) if I{Utm} is C{None}
                or I{cmoff} is C{False}.

       @raise TypeError: If I{latlon} is not ellipsoidal.

       @raise RangeError: If I{lat} is outside the valid UTM bands or if
                          I{lat} or I{lon} outside the valid range and
                          I{rangerrrors} set to C{True}.

       @raise ValueError: If I{lon} value is missing or if I{latlon}
                          is invalid.

       @note: Implements Karney’s method, using 8-th order Krüger series,
              giving results accurate to 5 nm (or better) for distances
              up to 3900 km from the central meridian.

       @example:

       >>> p = LatLon(48.8582, 2.2945)  # 31 N 448251.8 5411932.7
       >>> u = toUtm(p)  # 31 N 448252 5411933
       >>> p = LatLon(13.4125, 103.8667) # 48 N 377302.4 1483034.8
       >>> u = toUtm(p)  # 48 N 377302 1483035
    '''
    try:
        lat, lon = latlon.lat, latlon.lon
        if not isinstance(latlon, _LLEB):
            raise TypeError('%s not %s: %r' % ('latlon', 'ellipsoidal', latlon))
        if not name:  # use latlon.name
            name = _nameof(latlon) or name  # PYCHOK no effect
        d = datum or latlon.datum
    except AttributeError:
        lat, lon = parseDMS2(latlon, lon)
        d = datum or Datums.WGS84

    E = d.ellipsoid

    z, B, a, b = _toZBab4(lat, lon, cmoff)

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

    Ks = _Kseries(E.AlphaKs, x, y)  # Krüger series
    y = Ks.ys(y) * A0  # ξ
    x = Ks.xs(x) * A0  # η

    if cmoff:
        # C.F.F. Karney, "Test data for the transverse Mercator projection (2009)",
        # <http://GeographicLib.SourceForge.io/html/transversemercator.html> and
        # <http://Zenodo.org/record/32470#.W4LEJS2ZON8>
        x += _FalseEasting  # make x relative to false easting
        if y < 0:
            y += _FalseNorthing  # y relative to false northing in S

    # convergence: Karney 2011 Eq 23, 24
    p_ = Ks.ps(1)
    q_ = Ks.qs(0)
    c = degrees(atan(T_ / hypot1(T_) * tb) + atan2(q_, p_))

    # scale: Karney 2011 Eq 25
    s = E.e2s(sin(a)) * T12 / H * (A0 / E.a * hypot(p_, q_))

    if cmoff and Utm is not None:
        h = 'S' if a < 0 else 'N'  # hemisphere
        return _xnamed(Utm(z, h, x, y, band=B, datum=d,
                                       convergence=c, scale=s), name)
    else:  # zone, easting, northing, band, convergence and scale
        return z, x, y, B, c, s


def utmZoneBand2(lat, lon):
    '''Return the UTM zone number and UTM Band letter for a location.

       @param lat: Latitude (C{degrees}) or string.
       @param lon: Longitude (C{degrees}) or string.

       @return: 2-Tuple (zone, Band) as (int, string).

       @raise RangeError: If I{lat} is outside the valid UTM bands or if
                          I{lat} or I{lon} outside the valid range and
                          I{rangerrrors} set to C{True}.

       @raise ValueError: Invalid I{lat} or I{lon}.
    '''
    return _toZBab4(*parseDMS2(lat, lon))[:2]


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
