
# -*- coding: utf-8 -*-

u'''Universal Transverse Mercator (UTM) class L{Utm} and functions
L{parseUTM} and L{toUtm}.

Pure Python implementation of UTM / WGS-84 conversion functions using
an ellipsoidal earth model, transcribed from JavaScript originals by
I{(C) Chris Veness 2011-2016} published under the same MIT Licence**, see
U{UTM<http://www.movable-type.co.uk/scripts/latlong-utm-mgrs.html>} and
U{Module utm<http://www.movable-type.co.uk/scripts/geodesy/docs/module-utm.html>}.

The UTM system is a 2-dimensional cartesian coordinate system providing
locations on the surface of the earth.

UTM is a set of 60 transverse Mercator projections, normally based on
the WGS-84 ellipsoid.  Within each zone, coordinates are represented
as eastings and northings, measured in metres.

This method based on Karney U{'Transverse Mercator with an accuracy of
a few nanometers'<http://arxiv.org/pdf/1002.1417v3.pdf>}, 2011 (building
on Krüger U{'Konforme Abbildung des Erdellipsoids in der Ebene'
<http://bib.gfz-potsdam.de/pub/digi/krueger2.pdf>}, 1912).

Other references Seidel U{'Die Mathematik der Gauß-Krueger-Abbildung'
<http://henrik-seidel.gmxhome.de/gausskrueger.pdf>}, 2006,
U{Transverse Mercator Projection<http://geographiclib.sourceforge.io/tm.html>},
and U{Universal Transverse Mercator coordinate system
<http://wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system>}.

@newfield example: Example, Examples
'''

from bases import Base
from datum import Datums
from dms import S_DEG, parseDMS2
from ellipsoidalBase import LatLonEllipsoidalBase as _eLLb
from fmath import EPS, fdot3, fStr, Fsum, hypot, hypot1, \
                  isscalar, len2, map2
from utils import degrees, degrees90, degrees180, \
                  radians, wrap90, wrap180

from math import asinh, atan, atanh, atan2, \
                 cos, cosh, sin, sinh, tan, tanh
from operator import mul

# all public contants, classes and functions
__all__ = ('Utm',  # classes
           'parseUTM', 'toUtm')  # functions
__version__ = '18.02.09'

# Latitude bands C..X of 8° each, covering 80°S to 84°N with X repeated
# for 80-84°N
_Bands         = 'CDEFGHJKLMNPQRSTUVWXX'  #: (INTERNAL) Latitude bands.
_FalseEasting  =   500e3  #: (INTERNAL) False (meter).
_FalseNorthing = 10000e3  #: (INTERNAL) False (meter).
_K0            = 0.9996   #: (INTERNAL) UTM scale central meridian.


class _K6(object):
    '''(INTERNAL) Alpha or Beta Krüger series.

       Krüger series summations for I{eta}, I{ksi}, I{p} and I{q}
       while caching the cos, sin, cosh and sinh values for the
       given I{eta} and I{ksi} angles (in radians).
    '''
    def __init__(self, AB, x, y):
        '''(INTERNAL) New Alpha or Beta Krüger series

           @param AB: 6th-order Krüger Alpha or Beta series (1-origin).
           @param x: Eta angle (radians).
           @param y: Ksi angle (radians).
        '''
        n, j2 = len2(range(2, len(AB) * 2, 2))

        self._ab = AB[1:]  # 0-origin
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
        '''(INTERNAL) Eta summations (float).
        '''
        return fdot3(self._ab, self._cy, self._shx, start=x0)

    def ys(self, y0):
        '''(INTERNAL) Ksi summations (float).
        '''
        return fdot3(self._ab, self._sy, self._chx, start=y0)

    def ps(self, p0):
        '''(INTERNAL) P summations (float).
        '''
        return fdot3(self._pq, self._cy, self._chx, start=p0)

    def qs(self, q0):
        '''(INTERNAL) Q summations (float).
        '''
        return fdot3(self._pq, self._sy, self._shx, start=q0)


def _cml(zone):
    '''(INTERNAL) Central meridian longitude (degrees180).
    '''
    return (zone * 6) - 183


def _toZBL(zone, band, mgrs=False):  # used by mgrs.Mgrs
    '''(INTERNAL) Check and return zone, Band and band latitude.

       @param zone: Zone number or string.
       @param band: Band letter.
       @param mgrs: Optionally, raise ValueError (bool).

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
        raise ValueError('%s invalid: %r' % ('zone', zone))

    b = None
    if B:
        b = _Bands.find(B)
        if b < 0:
            raise ValueError('%s invalid: %r' % ('band', x))
        b = (b << 3) - 80
    elif mgrs:
        raise ValueError('%s missing' % ('band',))

    return z, B, b


def _toZBll(lat, lon):
    '''(INTERNAL) Return zone, Band and central lat- and longitude.

       @param lat: Latitude (degrees).
       @param lon: Longitude (degrees).

       @return: 4-Tuple (zone, Band, lat, lon).
    '''
    # return zone, Band and central
    # lat- and longitude (in radians)
    lat = wrap90(lat)
    if -80 > lat or lat > 84:
        raise ValueError('%s outside UTM: %s' % ('lat', lat))
    B = _Bands[int(lat + 80) >> 3]

    lon = wrap180(lon)
    z = int((lon + 180) / 6) + 1  # longitudinal zone
    if B == 'X':
        x = {32: 9, 34: 21, 36: 33}.get(z, None)
        if x:  # Svalbard
            if lon >= x:
                z += 1
            else:
                z -= 1
    elif B == 'V' and z == 31 and lon >=3:
        z += 1  # southern Norway

    b = radians(lon - _cml(z))  # lon off central meridian
    a = radians(lat)  # lat off equator
    return z, B, a, b


class Utm(Base):
    '''Universal Transverse Mercator (UTM) coordinate.
    '''
    _band        = ''    #: (INTERNAL) Latitude band letter ('C..X').
    _convergence = None  #: (INTERNAL) Meridian conversion (degrees).
    _datum       = Datums.WGS84  #: (INTERNAL) L{Datum}.
    _easting     = 0     #: (INTERNAL) Easting from false easting (meter).
    _hemi        = ''    #: (INTERNAL) Hemisphere ('N' or 'S')
    # _latlon also set by ellipsoidalBase.LatLonEllipsoidalBase.toUtm
    _latlon      = None  #: (INTERNAL) toLatLon cache.
    _mgrs        = None  #: (INTERNAL) toMgrs cache.
    _northing    = 0     #: (INTERNAL) Northing from false northing (meter).
    _scale       = None  #: (INTERNAL) Grid scale factor (scalar or None).
    _zone        = 0     #: (INTERNAL) Longitudinal zone (1..60).

    def __init__(self, zone, hemisphere, easting, northing, band='',
                       datum=Datums.WGS84, convergence=None, scale=None):
        '''New UTM coordinate.

           @param zone: UTM 6° longitudinal zone (int 1..60 covering 180°W.
                        180°E) or '00B' zone and band letter (string).
           @param hemisphere: N for the northern or S for the southern
                              hemisphere (string).
           @param easting: Easting from false easting, -500km from
                           central meridian (meter).
           @param northing: Northing from equator N or from false
                            northing -10,000km S (meter).
           @keyword band: Optional, latitudinal band (string, C..X).
           @keyword datum: Optional, this coordinate's datum (L{Datum}).
           @keyword convergence: Optional meridian convergence, bearing
                                 of grid North, clockwise from true
                                 North (degrees or None).
           @keyword scale: Optional grid scale factor (scalar or None).

           @raise ValueError: Invalid I{easting} or I{northing}.

           @example:

           >>> import pygeodesy
           >>> u = pygeodesy.Utm(31, 'N', 448251, 5411932)
        '''
        self._zone, B, _ = _toZBL(zone, band)

        h = str(hemisphere)[:1].upper()
        if not h or h not in ('N', 'S'):
            raise ValueError('%s invalid: %r' % ('hemisphere', hemisphere))

        e, n = float(easting), float(northing)
        # check easting/northing (with 40km overlap
        # between zones) - is this worthwhile?
        if 120e3 > e or e > 880e3:
            raise ValueError('%s invalid: %r' % ('easting', easting))
        if 0 > n or n > _FalseNorthing:
            raise ValueError('%s invalid: %r' % ('northing', northing))

        self._hemi        = h
        self._easting     = e
        self._northing    = n
        self._band        = B
        self._datum       = datum
        self._convergence = convergence
        self._scale       = scale

    @property
    def band(self):
        '''Get the latitudinal band (C..X or '').
        '''
        return self._band

    @property
    def convergence(self):
        '''Get the meridian convergence (degrees or None).
        '''
        return self._convergence

    @property
    def datum(self):
        '''Get the datum (L{Datum}).
        '''
        return self._datum

    @property
    def easting(self):
        '''Get the easting (meter).'''
        return self._easting

    @property
    def hemisphere(self):
        '''Get the hemisphere (N|S).
        '''
        return self._hemi

    @property
    def northing(self):
        '''Get the northing (meter).
        '''
        return self._northing

    def parseUTM(self, strUTM):
        '''Parse a string to a UTM coordinate.

           For more details, see function L{parseUTM} in
           this module L{utm}.
        '''
        return parseUTM(strUTM, datum=self.datum)

    @property
    def scale(self):
        '''Get the grid scale (scalar or None).
        '''
        return self._scale

    def toLatLon(self, LatLon):
        '''Convert this UTM coordinate to an (ellipsoidal) geodetic point.

           @param LatLon: The I{LatLon} class to use for the point
                          (I{LatLon}) or None.

           @return: Point of this UTM coordinate (I{LatLon}) or 5-tuple
                    (lat, lon, datum, convergence, scale) if I{LatLon}
                    is None.

           @raise TypeError: If I{LatLon} is not ellipsoidal.

           @raise ValueError: Invalid meridional radius or H-value.

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
            raise ValueError('%s invalid: %r' % ('meridional', E.A))
        x /= A0  # η eta
        y /= A0  # ξ ksi

        B6 = _K6(E.Beta6, x, y)  # 6th-order Krüger series, 1-origin
        y = -B6.ys(-y)  # ξ'
        x = -B6.xs(-x)  # η'

        shx = sinh(x)
        cy, sy = cos(y), sin(y)

        H = hypot(shx, cy)
        if H < EPS:
            raise ValueError('%s invalid: %r' % ('H', H))

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
            sd.fadd(d)
            T = sd.fsum()  # τi

        a = atan(T)  # lat
        b = atan2(shx, cy) + radians(_cml(self._zone))
        ll = _eLLb(degrees90(a), degrees180(b), datum=self._datum)

        # convergence: Karney 2011 Eq 26, 27
        p = -B6.ps(-1)
        q =  B6.qs(0)
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
        elif issubclass(LatLon, _eLLb):
            ll = LatLon(ll.lat, ll.lon, datum=ll.datum)
            ll._convergence = self._latlon.convergence
            ll._scale = self._latlon.scale
            return ll
        raise TypeError('%s not ellipsoidal: %r' % ('LatLon', LatLon))

    def toMgrs(self):
        '''Convert this UTM coordinate to an MGRS grid reference.

           See function L{toMgrs} in module L{mgrs} for more details.

           @return: The MGRS grid reference (L{Mgrs}).
        '''
        if self._mgrs is None:
            from mgrs import toMgrs  # PYCHOK recursive import
            self._mgrs = toMgrs(self)
        return self._mgrs

    def toStr(self, prec=0, sep=' ', B=False, cs=False):  # PYCHOK expected
        '''Return a string representation of this UTM coordinate.

           To distinguish from MGRS grid zone designators, a
           space is left between the zone and the hemisphere.

           Note that UTM coordinates are rounded, not truncated
           (unlike MGRS grid references).

           @keyword prec: Optional number of decimals, unstripped (int).
           @keyword sep: Optional separator to join (string).
           @keyword B: Optionally, include latitudinal band (bool).
           @keyword cs: Optionally, include meridian convergence and
                        grid scale factor (bool).

           @return: This UTM as string "00 N|S meter meter" plus
                    "degrees float" if I{cs} is True (string).

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

           @keyword prec: Optional number of decimals, unstripped (int).
           @keyword fmt: Optional, enclosing backets format (string).
           @keyword sep: Optional separator between name:value pairs (string).
           @keyword B: Optionally, include latitudinal band (bool).
           @keyword cs: Optionally, include meridian convergence and
                        grid scale factor (bool).

           @return: This UTM as "[Z:00, H:N|S, E:meter, N:meter]"
                    (string) plus "C:degrees, S:float" if I{cs} is True.
        '''
        t = self.toStr(prec=prec, sep=' ', B=B, cs=cs).split()
        k = 'ZHENCS' if cs else 'ZHEN'
        return fmt % (sep.join('%s:%s' % t for t in zip(k, t)),)

    @property
    def zone(self):
        '''Get the longitudinal zone (1..60).
        '''
        return self._zone


def parseUTM(strUTM, datum=Datums.WGS84, Utm=Utm):
    '''Parse a string representing a UTM coordinate, consisting
       of zone, hemisphere, easting and northing.

       @param strUTM: A UTM coordinate (string).
       @keyword datum: Optional datum to use (L{Datum}).
       @keyword Utm: Optional I{Utm} class to use for the UTM
                     coordinate (L{Utm}) or None.

       @return: The UTM coordinate (L{Utm}) or 4-tuple (zone,
                hemisphere, easting, northing) if I{Utm} is None.

       @raise ValueError: Invalid I{strUTM}.

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
        raise ValueError('%s invalid: %r' % ('strUTM', strUTM))

    return (z, h, e, n) if Utm is None else Utm(z, h, e, n, datum=datum)


def toUtm(latlon, lon=None, datum=None, Utm=Utm):
    '''Convert a lat-/longitude point to a UTM coordinate.

       @note: Implements Karney’s method, using 6-th order Krüger
       series, giving results accurate to 5 nm for distances up to
       3900 km from the central meridian.

       @param latlon: Latitude (degrees) or an (ellipsoidal)
                      geodetic I{LatLon} point.
       @keyword lon: Optional longitude (degrees or None).
       @keyword datum: Optional datum for this UTM coordinate,
                       overriding latlon's datum (I{Datum}).
       @keyword Utm: Optional I{Utm} class to usefor the UTM
                     coordinate (L{Utm}).

       @return: The UTM coordinate (L{Utm}).

       @raise TypeError: If I{latlon} is not ellipsoidal.

       @raise ValueError: If I{lon} value is missing, if I{latlon}
                          is not scalar or I{latlon} is outside
                          the valid UTM bands.

       @example:

       >>> p = LatLon(48.8582, 2.2945)  # 31 N 448251.8 5411932.7
       >>> u = toUtm(p)  # 31 N 448252 5411933
       >>> p = LatLon(13.4125, 103.8667) # 48 N 377302.4 1483034.8
       >>> u = toUtm(p)  # 48 N 377302 1483035
    '''
    try:
        lat, lon = latlon.lat, latlon.lon
        if not isinstance(latlon, _eLLb):
            raise TypeError('%s not %s: %r' % ('latlon', 'ellipsoidal', latlon))
        d = datum or latlon.datum
    except AttributeError:
        lat, lon = parseDMS2(latlon, lon)
        d = datum or Datums.WGS84

    E = d.ellipsoid

    z, B, a, b = _toZBll(lat, lon)
    h = 'S' if a < 0 else 'N'  # hemisphere

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

    A6 = _K6(E.Alpha6, x, y)  # 6th-order Krüger series, 1-origin
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
    s = E.e2s(sin(a)) * T12 / H * (A0 / E.a * hypot(p_, q_))

    return Utm(z, h, x, y, band=B, datum=d, convergence=c, scale=s)

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
