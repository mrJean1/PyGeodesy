
# -*- coding: utf-8 -*-

u'''Cassini-Soldner projection class L{CassiniSoldner} requiring
U{GeographicLib<http://PyPI.org/project/geographiclib/>}.

@newfield example: Example, Examples
'''

from bases import _Based, _nameof, _xattrs, _xnamed
from datum import Datums
from ellipsoidalBase import LatLonEllipsoidalBase as _LLEB
from fmath import fStr
from lazily import _ALL_LAZY
from utily import property_RO

# all public contants, classes and functions
__all__ = _ALL_LAZY.css
__version__ = '19.04.09'

_CassiniSoldner0 = None  # default projection


def _CassiniSoldner(cs0):
    '''(INTERNAL) Get/det default projection.
    '''
    if cs0 is None:
        global _CassiniSoldner0
        if _CassiniSoldner0 is None:
            _CassiniSoldner0 = CassiniSoldner(0, 0, name='Default')
        return _CassiniSoldner0

    elif isinstance(cs0, CassiniSoldner):
        return cs0

    raise TypeError('%s not %s: %r' % ('cs0', CassiniSoldner.__name__, cs0))


class CassiniSoldner(_Based):
    '''A Python version of Charles Karney's U{CassiniSoldner
       <http://GeographicLib.SourceForge.io/1.49/classGeographicLib_1_1CassiniSoldner.html>}
       C++ class.
    '''
    _cb0      = 0
    _datum    = Datums.WGS84  #: (INTERNAL) L{Datum}.
    _latlon0  = ()
    _meridian = None
    _sb0      = 0

    def __init__(self, lat0, lon0, datum=Datums.WGS84, name=''):
        '''New L{CassiniSoldner} projection.

           @param lat0: Latitude of center point (C{degrees90}).
           @param lon0: Longitude of center point (C{degrees180}).
           @keyword datum: Optional, the geodesic datum (L{Datum}).
           @keyword name: Optional name (C{str}).

           @raise ImportError: Package U{GeographicLib<http://PyPI.org/
                               project/geographiclib>} missing.

           @example:

           >>> p = CassiniSoldner(48 + 50/60.0, 2 + 20/60.0)  # Paris
           >>> p.forward(50.9, 1.8)  # Calais
           (-37518.854545, 230003.561828)

           >>> p.reverse4(-38e3, 230e3)
           (50.899937, 1.793161, 89.580797, 0.999982)
        '''
        if name:
            self.name = name

        if datum and datum != self._datum:
            self._datum = datum

        self.reset(lat0, lon0)

    def _xcopy(self, *attrs):
        '''(INTERNAL) Make copy with add'l, subclass attributes.
        '''
        return _xattrs(self.classof(self.lat0, self.lon0,
                                    datum=self.datum),
                       self, *attrs)

    def copy(self):
        '''Copy this Cassini-Soldner projection.

           @return: The copy (L{Utm} or subclass thereof).
        '''
        return self._xcopy()

    @property_RO
    def datum(self):
        '''Get the datum (L{Datum}).
        '''
        return self._datum

    @property_RO
    def flattening(self):
        '''Get the geodesic's flattening (C{float}).
        '''
        return self.geodesic.f

    def forward(self, lat, lon):
        '''Convert an (ellipsoidal) geodetic location Cassini-Soldner
           easting and northing.

           @param lat: Latitude of the location (C{degrees90}).
           @param lon: Longitude of the location (C{degrees180}).

           @return: 2-Tuple (C{easting, northing}) with C{easting} and
                    C{northing} in C{meters}.
        '''
        return self.forward4(lat, lon)[:2]

    def forward4(self, lat, lon):
        '''Convert an (ellipsoidal) geodetic location Cassini-Soldner
           easting and northing.

           @param lat: Latitude of the location (C{degrees90}).
           @param lon: Longitude of the location (C{degrees180}).

           @return: 4-Tuple (C{easting, northing, azi, rk}) for the
                    Cassini-Soldner location with C{easting} and
                    C{northing} in C{meters}, the azimuth of easting
                    direction C{azi} and the reciprocal of azimuthal
                    northing scale C{rk} in C{degrees}.
        '''
        g, M = self.datum.ellipsoid._geodesic2

        d = M.AngDiff(self.lon0, lon)[0]  # _2sum
        r = g.Inverse(lat, -abs(d), lat, abs(d))
        z1, a = r['azi1'], (r['a12'] * 0.5)
        z2, s = r['azi2'], (r['s12'] * 0.5)
        if s == 0:
            z = M.AngDiff(z1, z2)[0] * 0.5  # _2sum
            c = -90 if abs(d) > 90 else 90
            z1, z2 = c - z, c + z
        if d < 0:
            a, s, z2 = -a, -s, z1

        # z: azimuth of easting direction
        e, z = s, M.AngNormalize(z2)
        p = g.Line(lat, d, z, g.DISTANCE | g.GEODESICSCALE)
        # rk: reciprocal of azimuthal northing scale
        rk = p.ArcPosition(-a, g.GEODESICSCALE)['M21']
        # rk = p._GenPosition(True, -a, g.DISTANCE)[7]

        # s, c = M.sincosd(p.EquatorialAzimuth())
        s, c = M.sincosd(M.atan2d(p._salp0, p._calp0))
        sb1 = -c if lat < 0 else c
        cb1 = -abs(s) if abs(d) > 90 else abs(s)  # copysign(s, 90 - abs(d))
        d = M.atan2d(sb1 * self._cb0 - cb1 * self._sb0,
                     cb1 * self._cb0 + sb1 * self._sb0)
        n = self._meridian.ArcPosition(d, g.DISTANCE)['s12']
        # n = self._meridian._GenPosition(True, d, g.DISTANCE)[4]

        return e, n, z, rk

    @property_RO
    def geodesic(self):
        '''Get this projection's U{Geodesic
           <http://GeographicLib.SourceForge.io/html/python/code.html>},
           provided package U{geographiclib
           <http://PyPI.org/project/geographiclib>} is installed.
        '''
        return self._datum.ellipsoid.geodesic

    @property_RO
    def lat0(self):
        '''Get the center latitude (C{degrees90}).
        '''
        return self._latlon0[0]

    @property_RO
    def latlon0(self):
        '''Get the center lat- and longitude (C{degrees90}, C{degrees180}).
        '''
        return self._latlon0

    @property_RO
    def lon0(self):
        '''Get the center longitude (C{degrees180}).
        '''
        return self._latlon0[1]

    @property_RO
    def majoradius(self):
        '''Get the geodesic's major (equatorial) radius (C{float}).
        '''
        return self.geodetic.a

    def reset(self, lat0, lon0):
        '''Set the center point of this projection.

           @param lat0: Latitude of center point (C{degrees90}).
           @param lon0: Longitude of center point (C{degrees180}).
        '''
        g, M = self.datum.ellipsoid._geodesic2

        self._meridian = m = g.Line(lat0, lon0, 0.0, g.STANDARD | g.DISTANCE_IN)
        self._latlon0 = m.lat1, m.lon1
        s, c = M.sincosd(m.lat1)  # == self.lat0 == self..LatitudeOrigin()
        self._sb0, self._cb0 = M.norm(s * (1.0 - g.f), c)

    def reverse(self, easting, northing, LatLon=None):
        '''Convert a Cassini-Soldner location to (ellipsoidal) geodetic
           lat- and longitude.

           @param easting: Easting of the location (C{meter}).
           @param northing: Northing of the location (C{meter}).
           @keyword LatLon: Optional, ellipsoidal (sub-)class to return
                            the location as (C{LatLon}) or C{None}.

           @return: Geodetic location as I{LatLon} or 2-tuple (C{lat,
                    lon}) if I{LatLon} is C{None}.

           @raise TypeError: If I{LatLon} is not ellipsoidal.
        '''
        r = self.reverse4(easting, northing)[:2]
        if LatLon is None:
            pass
        elif issubclass(LatLon, _LLEB):
            r = LatLon(r[0], r[1], datum=self.datum, name=self.name)
        else:
            raise TypeError('%s not ellipsoidal: %r' % ('LatLon', LatLon))
        return r

    toLatLon = reverse

    def reverse4(self, easting, northing):
        '''Convert a Cassini-Soldner location to (ellipsoidal) geodetic
           lat- and longitude.

           @param easting: Easting of the location (C{meter}).
           @param northing: Northing of the location (C{meter}).

           @return: 4-Tuple (C{lat, lon, azi, rk}), all in C{degrees}
                    where C{azi} is the azimuth of easting direction
                    and C{rk} the reciprocal of azimuthal northing
                    scale.
        '''
        g = self.geodesic

        r = self._meridian.Position(northing)
        a, b, z = r['lat2'], r['lon2'], r['azi2']
        r = g.Direct(a, b, z + 90, easting, g.STANDARD | g.GEODESICSCALE)
        # include azimuth of easting direction and reciprocal of
        # azimuthal northing scale (see C++ member Direct() 5/6
        # <http://GeographicLib.SourceForge.io/1.49/classGeographicLib_1_1Geodesic.html>)
        return r['lat2'], r['lon2'], r['azi2'], r['M12']

    def toStr(self, prec=6, sep=' '):  # PYCHOK expected
        '''Return a string representation of this projection.

           @keyword prec: Optional number of decimal, unstripped (C{int}).
           @keyword sep: Optional separator to join (C{str}).

           @return: This projection as C{"lat0 lon0"} (C{str}).
        '''
        return fStr(self.latlon0, prec=prec, sep=sep)

    def toStr2(self, prec=6):  # PYCHOK expected
        '''Return a string representation of this projection.

           @keyword prec: Optional number of decimals, unstripped (C{int}).

           @return: This projection as C{"<classname>(lat0, lon0, ...)"}
                    (C{str}).
        '''
        t = self.toStr(prec=prec, sep=', ')
        n = self.name or ''
        if n:
            n = ', name=%r' % (n,)
        return '%s(%s%s)' % (self.classname, t, n)


class Css(_Based):
    '''Cassini-Soldner East-/Northing location.
    '''
    _azi      = None  #: (INTERNAL) azimuth of easting direction (C{degrees})
    _cs0      = None  #: (INTERNAL) projection (L{CassiniSoldner})
    _easting  = 0     #: (INTERNAL) Easting (C{float})
    _height   = 0     #: (INTERNAL) Height (C{meter})
    _latlon   = None  #: (INTERNAL) Geodetic (lat, lon)
    _northing = 0     #: (INTERNAL) Northing (C{float})
    _rk       = None  #: (INTERNAL) reciprocal of azimuthal northing scale (C{float})

    def __init__(self, e, n, h=0, cs0=_CassiniSoldner0, name=''):
        '''New L{Css} position.

           @param e: Easting (C{meter}).
           @param n: Northing (C{meter}).
           @keyword h: Optional height (C{meter}).
           @keyword cs0: Optional, the Cassini-Soldner projection
                         (L{CassiniSoldner}).
           @keyword name: Optional name (C{str}).

           @return: The Cassini-Soldner location (L{Css}).

           @raise ImportError: Package U{GeographicLib<http://PyPI.org/
                               project/geographiclib>} missing.

           @raise TypeError: If I{cs0} is not L{CassiniSoldner}.

           @raise ValueError: If I{e} or I{n} is invalid.

           @example:

           >>> cs = Css(448251, 5411932.0001)
        '''
        self._cs0 = _CassiniSoldner(cs0)
        self._easting  = float(e)
        self._northing = float(n)
        if h:
            self._height = float(h)
        if name:
            self.name = name

    def _reverse4(self):
        '''(INTERNAL) Convert to geodetic location.
        '''
        a, b, z, rk = self.cs0.reverse4(self.easting, self.northing)
        self._latlon = a, b
        self._azi, self._rk = z, rk
        return a, b, z, rk

    def _xcopy(self, *attrs):
        '''(INTERNAL) Make copy with add'l, subclass attributes.
        '''
        return _xattrs(self.classof(self.easting, self.northing,
                                    h=self.height, cs0=self.cs0),
                       self, *attrs)

    @property_RO
    def azi(self):
        '''Get the azimuth of easting direction (C{degrees}).
        '''
        if self._azi is None:
            self._reverse4()
        return self._azi

    @property_RO
    def cs0(self):
        '''Get the projection (L{CassiniSoldner}).
        '''
        return self._cs0

    def copy(self):
        '''Copy this Css location.

           @return: The copy (L{Css} or subclass thereof).
        '''
        return self._xcopy()

    @property_RO
    def easting(self):
        '''Get the easting (C{meter}).
        '''
        return self._easting

    @property_RO
    def height(self):
        '''Get the height (C{meter}).
        '''
        return self._height

    @property_RO
    def latlon(self):
        '''Get the lat- and longitude (C{degrees}, C{degrees}).
        '''
        if self._latlon is None:
            self._reverse4()
        return self._latlon

    @property_RO
    def northing(self):
        '''Get the northing (C{meter}).
        '''
        return self._northing

    @property_RO
    def rk(self):
        '''Get the reciprocal of azimuthal northing scale (L{float}).
        '''
        if self._rk is None:
            self._reverse4()
        return self._rk

    def toLatLon(self, LatLon=None, height=None):
        '''Convert this L{Css} to an (ellipsoidal) geodetic point.

           @keyword LatLon: Optional, ellipsoidal (sub-)class to return
                            the geodetic point (C{LatLon}) or C{None}.
           @keyword height: Optional height for the point, overriding
                            the default height (C{meter}).

           @return: The point (I{LatLon}) or 4-tuple (C{degrees90},
                    C{degrees180}, height, datum) if I{LatLon} is C{None}.

           @raise TypeError: If I{LatLon} or I{datum} is not ellipsoidal.
        '''
        if LatLon and not issubclass(LatLon, _LLEB):
            raise TypeError('%s not %s: %r' % ('LatLon', 'ellipsoidal', LatLon))

        a, b = self.latlon
        d = self.cs0.datum
        h = self.height if height is None else height

        return (a, b, h, d) if LatLon is None else _xnamed(LatLon(
                a, b, height=h, datum=d), self.name)

    def toStr(self, prec=6, sep=' ', m='m'):  # PYCHOK expected
        '''Return a string representation of this L{Css} position.

           @keyword prec: Optional number of decimal, unstripped (C{int}).
           @keyword sep: Optional separator to join (C{str}).
           @keyword m: Optional height units, default C{meter} (C{str}).

           @return: This Css as "easting nothing" C{str} in C{meter} plus
                    " height" and 'm' if heigth is non-zero (C{str}).
        '''
        t = [fStr(self.easting,  prec=prec),
             fStr(self.northing, prec=prec)]
        if self.height:
            t += ['%+.2f%s' % (self.height, m)]
        return sep.join(t)

    def toStr2(self, prec=6, fmt='[%s]', sep=', ', m='m', C=False):  # PYCHOK expected
        '''Return a string representation of this L{Css} position.

           @keyword prec: Optional number of decimals, unstripped (C{int}).
           @keyword fmt: Optional, enclosing backets format (C{str}).
           @keyword sep: Optional separator between name:values (C{str}).
           @keyword m: Optional unit of the height, default meter (C{str}).
           @keyword C: Optionally, include name of projection (C{bool}).

           @return: This Css as "[E:meter, N:meter, H:m, name:'', C:Conic.Datum]"
                    (C{str}).
        '''
        t = self.toStr(prec=prec, sep=' ', m=m).split()
        k = ('E', 'N', 'H')[:len(t)]
        if self.name:
            k += 'name',
            t += [repr(self.name)]
        if C:
            k += 'C',
            t += [repr(self.cs0)]
        return fmt % (sep.join('%s:%s' % t for t in zip(k, t)),)


def toCss(latlon, cs0=_CassiniSoldner0, height=None, Css=Css, name=''):
    '''Convert an (ellipsoidal) geodetic point to a Cassini-Soldner location.

       @param latlon: Ellipsoidal point (C{LatLon}).
       @keyword cs0: Optional, the Cassini-Soldner projection to use
                     (L{CassiniSoldner}).
       @keyword height: Optional height for the point, overriding
                        the default height (C{meter}).
       @keyword Css: Optional (sub-)class to return the location
                     (L{Css}) or C{None}.
       @keyword name: Optional I{Css} name (C{str}).

       @return: The Cassini-Soldner location (L{Css}) or 3-tuple
                (C{easting, northing, height}) if I{Css} is C{None}.

       @raise ImportError: Package U{GeographicLib<http://PyPI.org/
                           project/geographiclib>} missing.

       @raise TypeError: If I{latlon} is not ellipsoidal.
    '''
    if not isinstance(latlon, _LLEB):
        raise TypeError('%s not %s: %r' % ('latlon', 'ellipsoidal', latlon))

    cs = _CassiniSoldner(cs0)

    C, E = cs.datum.ellipsoid, latlon.datum.ellipsoid
    if C.a != E.a or C.f != E.f:
        raise ValueError('%s mistmatch %r vs %r' % ('ellipsoidal', C, E))

    e, n, z, rk = cs.forward4(latlon.lat, latlon.lon)
    h = latlon.height if height is None else height

    if Css is None:
        r = e, n, h
    else:
        r = _xnamed(Css(e, n, h=h, cs0=cs), name or _nameof(latlon))
        r._latlon = latlon.lat, latlon.lon
        r._azi, r._rk = z, rk
    return r

# **) MIT License
#
# Copyright (C) 2016-2019 -- mrJean1 at Gmail dot com
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
