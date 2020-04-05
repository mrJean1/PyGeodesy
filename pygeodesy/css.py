
# -*- coding: utf-8 -*-

u'''Cassini-Soldner projection classes L{CassiniSoldner}, L{Css} and
L{CSSError} requiring I{Charles Karney's} U{geographiclib
<https://PyPI.org/project/geographiclib/>} package to be installed.

@newfield example: Example, Examples
'''

from pygeodesy.basics import _isnotError, issubclassof, property_RO, \
                             _TypeError, _xkwds
from pygeodesy.datum import Datums
from pygeodesy.ellipsoidalBase import LatLonEllipsoidalBase as _LLEB
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import EasNor2Tuple, EasNor3Tuple, EasNorAziRk4Tuple, \
                            LatLon2Tuple, LatLon4Tuple, LatLonAziRk4Tuple, \
                            _NamedBase, nameof, _xnamed
from pygeodesy.streprs import fstr, strs
from pygeodesy.utily import false2f

# all public contants, classes and functions
__all__ = _ALL_LAZY.css
__version__ = '20.04.04'

_CassiniSoldner0 = None  # default projection


def _CassiniSoldner(cs0):
    '''(INTERNAL) Get/set default projection.
    '''
    if cs0 is None:
        global _CassiniSoldner0
        if _CassiniSoldner0 is None:
            _CassiniSoldner0 = CassiniSoldner(0, 0, name='Default')
        cs0 = _CassiniSoldner0
    else:
        _TypeError(CassiniSoldner, cs0=cs0)
    return cs0


class CassiniSoldner(_NamedBase):
    '''Cassini-Soldner projection, a Python version of Karney's C++ class U{CassiniSoldner
       <https://GeographicLib.SourceForge.io/1.49/classGeographicLib_1_1CassiniSoldner.html>}.
    '''
    _cb0      = 0
    _datum    = Datums.WGS84  #: (INTERNAL) L{Datum}.
    _latlon0  = ()
    _meridian = None
    _sb0      = 0

    def __init__(self, lat0, lon0, datum=Datums.WGS84, name=''):
        '''New L{CassiniSoldner} projection.

           @arg lat0: Latitude of center point (C{degrees90}).
           @arg lon0: Longitude of center point (C{degrees180}).
           @kwarg datum: Optional, the geodesic datum (L{Datum}).
           @kwarg name: Optional name (C{str}).

           @raise ImportError: Package U{GeographicLib<https://PyPI.org/
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

    @property_RO
    def datum(self):
        '''Get the datum (L{Datum}).
        '''
        return self._datum

    def _datumatch(self, latlon):
        '''Check for matching datum ellipsoids.

           @raise CSSError: Ellipsoidal mismatch of B{C{latlon}} and this projection.
        '''
        C, E = self.datum.ellipsoid, latlon.datum.ellipsoid
        if E != C:
            raise CSSError('%s mistmatch: %r vs %r' % ('ellipsoids', E, C))

    @property_RO
    def flattening(self):
        '''Get the geodesic's flattening (C{float}).
        '''
        return self.geodesic.f

    def forward(self, lat, lon):
        '''Convert an (ellipsoidal) geodetic location to Cassini-Soldner
           easting and northing.

           @arg lat: Latitude of the location (C{degrees90}).
           @arg lon: Longitude of the location (C{degrees180}).

           @return: An L{EasNor2Tuple}C{(easting, northing)}.

           @see: Methods L{CassiniSoldner.forward4}, L{CassiniSoldner.reverse}
                 and L{CassiniSoldner.reverse4}.
        '''
        r = EasNor2Tuple(*self.forward4(lat, lon)[:2])
        return self._xnamed(r)

    def forward4(self, lat, lon):
        '''Convert an (ellipsoidal) geodetic location to Cassini-Soldner
           easting and northing.

           @arg lat: Latitude of the location (C{degrees90}).
           @arg lon: Longitude of the location (C{degrees180}).

           @return: An L{EasNorAziRk4Tuple}C{(easting, northing,
                    azimuth, reciprocal)}.

           @see: Method L{CassiniSoldner.forward}, L{CassiniSoldner.reverse}
                 and L{CassiniSoldner.reverse4}.
        '''
        g, M = self.datum.ellipsoid._geodesic_Math2

        d = M.AngDiff(self.lon0, lon)[0]  # _2sum
        r = g.Inverse(lat, -abs(d), lat, abs(d))
        z1, a = r.azi1, (r.a12 * 0.5)
        z2, s = r.azi2, (r.s12 * 0.5)
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
        rk = p.ArcPosition(-a, g.GEODESICSCALE).M21
        # rk = p._GenPosition(True, -a, g.DISTANCE)[7]

        # s, c = M.sincosd(p.EquatorialAzimuth())
        s, c = M.sincosd(M.atan2d(p._salp0, p._calp0))
        sb1 = -c if lat < 0 else c
        cb1 = -abs(s) if abs(d) > 90 else abs(s)  # copysign(s, 90 - abs(d))
        d = M.atan2d(sb1 * self._cb0 - cb1 * self._sb0,
                     cb1 * self._cb0 + sb1 * self._sb0)
        n = self._meridian.ArcPosition(d, g.DISTANCE).s12
        # n = self._meridian._GenPosition(True, d, g.DISTANCE)[4]
        r = EasNorAziRk4Tuple(e, n, z, rk)
        return self._xnamed(r)

    @property_RO
    def geodesic(self):
        '''Get this projection's I{wrapped} U{Karney Geodesic
           <https://GeographicLib.SourceForge.io/html/python/code.html>},
           provided package U{geographiclib
           <https://PyPI.org/project/geographiclib>} is installed.
        '''
        return self.datum.ellipsoid.geodesic

    @property_RO
    def lat0(self):
        '''Get the center latitude (C{degrees90}).
        '''
        return self._latlon0.lat

    @property
    def latlon0(self):
        '''Get the center lat- and longitude (L{LatLon2Tuple}C{(lat, lon)}) in (C{degrees90}, (C{degrees180}).
        '''
        return self._latlon0

    @latlon0.setter  # PYCHOK setter!
    def latlon0(self, latlon0):
        '''Set the center lat- and longitude (L{LatLon2Tuple}, ellipsoidal C{LatLon} or L{LatLon4Tuple}).

           @raise CSSError: Ellipsoidal mismatch of B{C{latlon0}} and this projection.

           @raise TypeError: Invalid B{C{latlon0}}.
        '''
        _TypeError(_LLEB, LatLon4Tuple, LatLon2Tuple, latlon0=latlon0)
        if hasattr(latlon0, 'datum'):
            self._datumatch(latlon0)
        self.reset(latlon0.lat, latlon0.lon)

    @property_RO
    def lon0(self):
        '''Get the center longitude (C{degrees180}).
        '''
        return self._latlon0.lon

    @property_RO
    def majoradius(self):
        '''Get the geodesic's major (equatorial) radius (C{float}).
        '''
        return self.geodesic.a

    def reset(self, lat0, lon0):
        '''Set or reset the center point of this Cassini-Soldner projection.

           @arg lat0: Center point latitude (C{degrees90}).
           @arg lon0: Center point longitude (C{degrees180}).
        '''
        g, M = self.datum.ellipsoid._geodesic_Math2

        self._meridian = m = g.Line(lat0, lon0, 0.0, g.STANDARD | g.DISTANCE_IN)
        self._latlon0 = LatLon2Tuple(m.lat1, m.lon1)
        s, c = M.sincosd(m.lat1)  # == self.lat0 == self.LatitudeOrigin()
        self._sb0, self._cb0 = M.norm(s * (1.0 - g.f), c)

    def reverse(self, easting, northing, LatLon=None, **LatLon_kwds):
        '''Convert a Cassini-Soldner location to (ellipsoidal) geodetic
           lat- and longitude.

           @arg easting: Easting of the location (C{meter}).
           @arg northing: Northing of the location (C{meter}).
           @kwarg LatLon: Optional, ellipsoidal class to return the
                          geodetic location as (C{LatLon}) or C{None}.
           @kwarg LatLon_kwds: Optional (C{LatLon}) keyword arguments,
                               ignored if B{C{LatLon=None}}.

           @return: Geodetic location B{C{LatLon}} or if B{C{LatLon}}
                    is C{None}, a L{LatLon2Tuple}C{(lat, lon)}.

           @raise CSSError: Ellipsoidal mismatch of B{C{LatLon}} and this projection.

           @raise TypeError: Invalid B{C{LatLon}} or B{C{LatLon_kwds}}.

           @see: Method L{CassiniSoldner.reverse4}, L{CassiniSoldner.forward}
                 and L{CassiniSoldner.forward4}.
        '''
        r = self.reverse4(easting, northing)
        if LatLon is None:
            r = LatLon2Tuple(r.lat, r.lon)  # PYCHOK expected
        elif issubclassof(LatLon, _LLEB):
            kwds = _xkwds(LatLon_kwds, datum=self.datum)
            r = LatLon(r.lat, r.lon, **kwds)  # PYCHOK expected
            self._datumatch(r)
        else:
            raise _isnotError(_LLEB.__name__, LatLon=LatLon)
        return self._xnamed(r)

    toLatLon = reverse

    def reverse4(self, easting, northing):
        '''Convert a Cassini-Soldner location to (ellipsoidal) geodetic
           lat- and longitude.

           @arg easting: Easting of the location (C{meter}).
           @arg northing: Northing of the location (C{meter}).

           @return: A L{LatLonAziRk4Tuple}C{(lat, lon, azimuth, reciprocal)}.

           @see: Method L{CassiniSoldner.reverse}, L{CassiniSoldner.forward}
                 and L{CassiniSoldner.forward4}.
        '''
        g = self.geodesic

        n = self._meridian.Position(northing)
        r = g.Direct(n.lat2, n.lon2, n.azi2 + 90, easting, g.STANDARD | g.GEODESICSCALE)
        # include azimuth of easting direction and reciprocal of
        # azimuthal northing scale (see C++ member Direct() 5/6
        # <https://GeographicLib.SourceForge.io/1.49/classGeographicLib_1_1Geodesic.html>)
        r = LatLonAziRk4Tuple(r.lat2, r.lon2, r.azi2, r.M12)
        return self._xnamed(r)

    def toStr(self, prec=6, sep=' '):  # PYCHOK expected
        '''Return a string representation of this projection.

           @kwarg prec: Optional number of decimal, unstripped (C{int}).
           @kwarg sep: Optional separator to join (C{str}).

           @return: This projection as C{"lat0 lon0"} (C{str}).
        '''
        return sep.join(strs(self.latlon0, prec=prec))

    def toStr2(self, prec=6):  # PYCHOK expected
        '''Return a string representation of this projection.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).

           @return: This projection as C{"<classname>(lat0, lon0, ...)"}
                    (C{str}).
        '''
        t = self.toStr(prec=prec, sep=', ')
        n = self.name or ''
        if n:
            n = ', name=%r' % (n,)
        return '%s(%s%s)' % (self.classname, t, n)


class CSSError(ValueError):
    '''Cassini-Soldner (CSS) conversion or other L{Css} issue.
    '''
    pass


class Css(_NamedBase):
    '''Cassini-Soldner East-/Northing location.
    '''
    _cs0      = None  #: (INTERNAL) projection (L{CassiniSoldner})
    _easting  = 0     #: (INTERNAL) Easting (C{float})
    _height   = 0     #: (INTERNAL) Height (C{meter})
    _northing = 0     #: (INTERNAL) Northing (C{float})
    _reverse4 = None  #: (INTERNAL) Cached reverse4 (L{LatLonAziRk4Tuple})

    def __init__(self, e, n, h=0, cs0=_CassiniSoldner0, name=''):
        '''New L{Css} Cassini-Soldner position.

           @arg e: Easting (C{meter}).
           @arg n: Northing (C{meter}).
           @kwarg h: Optional height (C{meter}).
           @kwarg cs0: Optional, the Cassini-Soldner projection
                       (L{CassiniSoldner}).
           @kwarg name: Optional name (C{str}).

           @return: The Cassini-Soldner location (L{Css}).

           @raise CSSError: If B{C{e}} or B{C{n}} is invalid.

           @raise ImportError: Package U{GeographicLib<https://PyPI.org/
                               project/geographiclib>} missing.

           @raise TypeError: If B{C{cs0}} is not L{CassiniSoldner}.

           @example:

           >>> cs = Css(448251, 5411932.0001)
        '''
        self._cs0 = _CassiniSoldner(cs0)
        self._easting  = false2f(e, 'easting',  false=False, Error=CSSError)
        self._northing = false2f(n, 'northing', false=False, Error=CSSError)
        if h:
            self._height = float(h)
        if name:
            self.name = name

    @property_RO
    def azi(self):
        '''Get the azimuth of easting direction (C{degrees}).
        '''
        return self.reverse4.azimuth

    azimuth = azi

    @property_RO
    def cs0(self):
        '''Get the projection (L{CassiniSoldner}).
        '''
        return self._cs0

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
        '''Get the lat- and longitude (L{LatLon2Tuple}).
        '''
        r = LatLon2Tuple(self.reverse4.lat, self.reverse4.lon)
        return self._xnamed(r)

    @property_RO
    def northing(self):
        '''Get the northing (C{meter}).
        '''
        return self._northing

    @property_RO
    def reverse4(self):
        '''Get the lat, lon, azimuth and reciprocal (L{LatLonAziRk4Tuple}).
        '''
        if self._reverse4 is None:
            self._reverse4 = self.cs0.reverse4(self.easting, self.northing)
        return self._reverse4

    @property_RO
    def rk(self):
        '''Get the reciprocal of azimuthal northing scale (C{degrees}).
        '''
        return self.reverse4.reciprocal

    reciprocal = rk

    def toLatLon(self, LatLon=None, height=None, **LatLon_kwds):
        '''Convert this L{Css} to an (ellipsoidal) geodetic point.

           @kwarg LatLon: Optional, ellipsoidal class to return the
                          geodetic point (C{LatLon}) or C{None}.
           @kwarg height: Optional height for the point, overriding the
                          default height (C{meter}).
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                               arguments, ignored if B{C{LatLon=None}}.

           @return: The geodetic point (B{C{LatLon}}) or if B{C{LatLon}}
                    is C{None}, a L{LatLon4Tuple}C{(lat, lon, height,
                    datum)}.

           @raise TypeError: If B{C{LatLon}} or B{C{datum}} is not
                             ellipsoidal or invalid B{C{LatLon_kwds}}.
        '''
        if LatLon and not issubclassof(LatLon, _LLEB):
            raise _isnotError(_LLEB.__name__, LatLon=LatLon)

        lat, lon = self.latlon
        d = self.cs0.datum
        h = self.height if height is None else height

        r = LatLon4Tuple(lat, lon, h, d) if LatLon is None else \
                  LatLon(lat, lon, height=h, datum=d, **LatLon_kwds)
        return self._xnamed(r)

    def toStr(self, prec=6, sep=' ', m='m'):  # PYCHOK expected
        '''Return a string representation of this L{Css} position.

           @kwarg prec: Optional number of decimal, unstripped (C{int}).
           @kwarg sep: Optional separator to join (C{str}).
           @kwarg m: Optional height units, default C{meter} (C{str}).

           @return: This position as C{"easting nothing"} C{str} in
                    C{meter} plus C{" height"} and C{'m'} if heigth
                    is non-zero (C{str}).
        '''
        t = [fstr(self.easting,  prec=prec),
             fstr(self.northing, prec=prec)]
        if self.height:  # abs(self.height) > EPS
            t += ['%+.2f%s' % (self.height, m)]
        return sep.join(t)

    def toStr2(self, prec=6, fmt='[%s]', sep=', ', m='m', C=False):  # PYCHOK expected
        '''Return a string representation of this L{Css} position.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg fmt: Optional, enclosing backets format (C{str}).
           @kwarg sep: Optional separator between name:values (C{str}).
           @kwarg m: Optional unit of the height, default meter (C{str}).
           @kwarg C: Optionally, include name of projection (C{bool}).

           @return: This position as C{"[E:meter, N:meter, H:m, name:'',
                    C:Conic.Datum]"} (C{str}).
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
    '''Convert an (ellipsoidal) geodetic point to a Cassini-Soldner
       location.

       @arg latlon: Ellipsoidal point (C{LatLon} or L{LatLon4Tuple}).
       @kwarg cs0: Optional, the Cassini-Soldner projection to use
                   (L{CassiniSoldner}).
       @kwarg height: Optional height for the point, overriding the
                      default height (C{meter}).
       @kwarg Css: Optional class to return the location (L{Css}) or C{None}.
       @kwarg name: Optional B{C{Css}} name (C{str}).

       @return: The Cassini-Soldner location (B{C{Css}}) or an
                L{EasNor3Tuple}C{(easting, northing, height)}
                if B{C{Css}} is C{None}.

       @raise CSSError: Ellipsoidal mismatch of B{C{latlon}} and B{C{cs0}}.

       @raise ImportError: Package U{GeographicLib<https://PyPI.org/
                           project/geographiclib>} missing.

       @raise TypeError: If B{C{latlon}} is not ellipsoidal.
    '''
    _TypeError(_LLEB, LatLon4Tuple, latlon=latlon)

    cs = _CassiniSoldner(cs0)
    cs._datumatch(latlon)

    c = cs.forward4(latlon.lat, latlon.lon)
    h = latlon.height if height is None else height

    if Css is None:
        r = EasNor3Tuple(c.easting, c.northing, h)
    else:
        r = Css(c.easting, c.northing, h=h, cs0=cs)
        r._latlon = LatLon2Tuple(latlon.lat, latlon.lon)
        r._azi, r._rk = c.azimuth, c.reciprocal
    return _xnamed(r, name or nameof(latlon))

# **) MIT License
#
# Copyright (C) 2016-2020 -- mrJean1 at Gmail -- All Rights Reserved.
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
