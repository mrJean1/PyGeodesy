
# -*- coding: utf-8 -*-

u'''Web Mercator (WM) class L{Wm} and functions L{parseWM} and L{toWm}.

Pure Python implementation of a U{Web Mercator<http://WikiPedia.org/wiki/Web_Mercator>}
(aka I{Pseudo-Mercator}) class and conversion functions for spherical and
near-spherical earth models.

References U{Google Maps / Bing Maps Spherical Mercator Projection
<http://alastaira.WordPress.com/2011/01/23/the-google-maps-bing-maps-spherical-mercator-projection>},
U{Geomatics Guidance Note 7, part 2<http://www.EPSG.org/Portals/0/373-07-02.pdf>} and
U{Implementation Practice Web Mercator Map Projection
<http://Earth-Info.NGA.mil/GandG/wgs84/web_mercator/%28U%29%20NGA_SIG_0011_1.0.0_WEBMERC.pdf>}.

@newfield example: Example, Examples
'''

from bases import Based, classname, _nameof, _xnamed
from datum import R_MA
from dms import clipDMS, parseDMS2
from ellipsoidalBase import LatLonEllipsoidalBase as _ELLB
from fmath import EPS, fStr, isscalar, map1
from utily import PI_2, degrees90, degrees180, property_RO

from math import atan, atanh, exp, radians, sin, tanh

# all public contants, classes and functions
__all__ = ('Wm',  # classes
           'parseWM', 'toWm')  # functions
__version__ = '18.10.12'

# _FalseEasting  = 0   #: (INTERNAL) False Easting (C{meter}).
# _FalseNorthing = 0   #: (INTERNAL) False Northing (C{meter}).
_LatLimit = 85.051129  #: (INTERNAL) Latitudinal limit (C{degrees}).
# _LonOrigin     = 0   #: (INTERNAL) Longitude of natural origin (C{degrees}).


class Wm(Based):
    '''Web Mercator (WM) coordinate.
    '''
    _radius = 0  #: (INTERNAL) earth radius (C{meter}).
    _x      = 0  #: (INTERNAL) easting (C{meter}).
    _y      = 0  #: (INTERNAL) northing (C{meter}).

    def __init__(self, x, y, radius=R_MA, name=''):
        '''New Web Mercator (WM) coordinate.

           @param x: Easting from central meridian (C{meter}).
           @param y: Northing from equator (C{meter}).
           @keyword radius: Optional earth radius (C{meter}).
           @keyword name: Optional name (C{str}).

           @raise ValueError: Invalid I{x}, I{y} or I{radius}.

           @example:

           >>> import pygeodesy
           >>> w = pygeodesy.Wm(448251, 5411932)
        '''
        if name:
            self.name = name

        try:
            self._x, self._y, r = map1(float, x, y, radius)
        except (TypeError, ValueError):
            raise ValueError('%s invalid: %r' % (Wm.__name__, (x, y, radius)))

        if r < EPS:  # check radius
            t = '%s.%s' % (classname(self), 'radius')
            raise ValueError('%s invalid: %r' % (t, r))
        self._radius = r

    def copy(self):
        '''Copy this Web Marcator coordinate.

           @return: The copy (L{Wm} or subclass thereof).
        '''
        return self.classof(self.x, self.y, radius=self.radius)

    def parseWM(self, strWM, name=''):
        '''Parse a string to a WM coordinate.

           For more details, see function L{parseWM} in
           this module L{webmercator}.
        '''
        return parseWM(strWM, radius=self.radius, Wm=self.classof, name=name)

    @property_RO
    def radius(self):
        '''Get the earth radius (C{meter}).
        '''
        return self._radius

    def to2ll(self, datum=None):
        '''Convert this WM coordinate to a geodetic lat- and longitude.

           @keyword datum: Optional datum (C{Datum}).

           @return: 2-Tuple (lat, lon) in (C{degrees90}, C{degrees180}).

           @raise TypeError: Non-ellipsoidal I{datum}.

           @raise ValueError: Invalid I{radius}.
        '''
        r = self.radius
        x = self._x / r
        y = 2 * atan(exp(self._y / r)) - PI_2
        if datum:
            E = datum.ellipsoid
            if not E.isEllipsoidal:
                raise TypeError('%s not %s: %r' % ('datum', 'ellipsoidal', datum))
            # <http://Earth-Info.NGA.mil/GandG/wgs84/web_mercator/
            #       %28U%29%20NGA_SIG_0011_1.0.0_WEBMERC.pdf>
            y = y / r
            if E.e:
                y -= E.e * atanh(E.e * tanh(y))
            y *= E.a
            x *= E.a / r
        return degrees90(y), degrees180(x)

    def toLatLon(self, LatLon, datum=None):
        '''Convert this WM coordinate to a geodetic point.

           @param LatLon: Ellipsoidal (sub-)class to use for the
                          point (C{LatLon}).
           @keyword datum: Optional datum for ellipsoidal or C{None} for
                           spherical I{LatLon} (C{Datum}).

           @return: Point of this WM coordinate (I{LatLon}).

           @raise TypeError: If I{LatLon} and I{datum} are not compatible
                             or if I{datum} is not ellipsoidal.

           @raise ValueError: Invalid I{radius}.

           @example:

           >>> w = Wm(448251.795, 5411932.678)
           >>> from pygeodesy import sphericalTrigonometry as sT
           >>> ll = w.toLatLon(sT.LatLon)  # 43°39′11.58″N, 004°01′36.17″E
        '''
        if issubclass(LatLon, _ELLB):
            if datum:
                return _xnamed(LatLon(*self.to2ll(datum=datum), datum=datum),
                               self.name)
        elif datum is None:
            return _xnamed(LatLon(*self.to2ll(datum=datum)), self.name)
        raise TypeError('%r and %s %r' % (LatLon, 'datum', datum))

    def toStr(self, prec=3, sep=' ', radius=False):  # PYCHOK expected
        '''Return a string representation of this WM coordinate.

           @keyword prec: Optional number of decimals, unstripped (C{int}).
           @keyword sep: Optional separator to join (C{str}).
           @keyword radius: Optionally, include radius (C{bool} or C{scalar}).

           @return: This WM as "meter meter" (C{str}) plus " radius"
                    if I{radius} is C{True} or C{scalar}.

           @raise ValueError: Invalid I{radius}.

           @example:

           >>> w = Wm(448251, 5411932.0001)
           >>> w.toStr(4)  # 448251.0 5411932.0001
           >>> w.toStr(sep=', ')  # 448251, 5411932
        '''
        fs = self._x, self._y
        if radius in (False, None):
            pass
        elif radius is True:
            fs += (self._radius,)
        elif isscalar(radius):
            fs += (radius,)
        else:
            raise ValueError('% invalid: %r' % ('radius', radius))
        return fStr(fs, prec=prec, sep=sep)

    def toStr2(self, prec=3, fmt='[%s]', sep=', ', radius=False):  # PYCHOK expected
        '''Return a string representation of this WM coordinate.

           @keyword prec: Optional number of decimals, unstripped (C{int}).
           @keyword fmt: Optional, enclosing backets format (C{str}).
           @keyword sep: Optional separator between name:value pairs (C{str}).
           @keyword radius: Optionally, include radius (C{bool} or C{scalar}).

           @return: This WM as "[x:meter, y:meter]" (C{str}) plus
                    ", radius:meter]" if I{radius} is C{True} or
                    C{scalar}.
        '''
        t = self.toStr(prec=prec, sep=' ', radius=radius).split()
        k = 'x', 'y', 'radius'
        return fmt % (sep.join('%s:%s' % t for t in zip(k, t)),)

    @property_RO
    def x(self):
        '''Get the easting (C{meter}).'''
        return self._x

    @property_RO
    def y(self):
        '''Get the northing (C{meter}).
        '''
        return self._y


def parseWM(strWM, radius=R_MA, Wm=Wm, name=''):
    '''Parse a string representing a WM coordinate, consisting
       of easting, northing and an optional radius.

       @param strWM: A WM coordinate (C{str}).
       @keyword radius: Optional earth radius (C{meter}).
       @keyword Wm: Optional (sub-)class to use (L{Wm}) or C{None}.
       @keyword name: Optional name (C{str}).

       @return: The WM coordinate (L{Wm}) or 3-tuple (easting,
                northing, radius) if I{Wm} is C{None}.

       @raise ValueError: Invalid I{strWM}.

       @example:

       >>> u = parseWM('448251 5411932')
       >>> u.toStr2()  # [E:448251, N:5411932]
    '''
    w = strWM.strip().replace(',', ' ').split()
    try:
        if len(w) == 2:
            w += [radius]
        elif len(w) != 3:
            raise ValueError  # caught below
        x, y, r = map(float, w)

    except (TypeError, ValueError):
        raise ValueError('%s invalid: %r' % ('strWM', strWM))

    return (x, y, r) if Wm is None else _xnamed(Wm(
            x, y, radius=r), name)


def toWm(latlon, lon=None, radius=R_MA, Wm=Wm, name=''):
    '''Convert a lat-/longitude point to a WM coordinate.

       @param latlon: Latitude (C{degrees}) or an (ellipsoidal or
                      spherical) geodetic C{LatLon} point.
       @keyword lon: Optional longitude (C{degrees} or C{None}).
       @keyword radius: Optional earth radius (C{meter}).
       @keyword Wm: Optional (sub-)class for the WM coordinate
                    (L{Wm}) or C{None}.
       @keyword name: Optional name (C{str}).

       @return: The WM coordinate (L{Wm}) or 3-tuple (easting,
                northing, radius) if I{Wm} is C{None}.

       @raise ValueError: If I{lon} value is missing, if I{latlon}
                          is not scalar, if I{latlon} is beyond the
                          valid WM range and L{rangerrors} is set
                          to C{True} or if I{radius} is invalid.

       @example:

       >>> p = LatLon(48.8582, 2.2945)  # 448251.8 5411932.7
       >>> w = toWm(p)  # 448252 5411933
       >>> p = LatLon(13.4125, 103.8667)  # 377302.4 1483034.8
       >>> w = toWm(p)  # 377302 1483035
    '''
    r, e = radius, None
    try:
        lat, lon = latlon.lat, latlon.lon
        if isinstance(latlon, _ELLB):
            r = latlon.datum.ellipsoid.a
            e = latlon.datum.ellipsoid.e
            if not name:  # use latlon.name
                name = _nameof(latlon) or name  # PYCHOK no effect
        lat = clipDMS(lat, _LatLimit)
    except AttributeError:
        lat, lon = parseDMS2(latlon, lon, clipLat=_LatLimit)

    s = sin(radians(lat))
    y = atanh(s)  # == log(tan((90 + lat) / 2)) == log(tanPI_2_2(radians(lat)))
    if e:
        y -= e * atanh(e * s)

    e, n = r * radians(lon), r * y
    return (e, n, r) if Wm is None else _xnamed(Wm(
            e, n, radius=r), name)

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
