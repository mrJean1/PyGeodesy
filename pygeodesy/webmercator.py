
# -*- coding: utf-8 -*-

u'''Web Mercator (WM) class L{Wm} and functions L{parseWM} and L{toWm}.

Pure Python implementation of a U{Web Mercator<http://wikipedia.org/wiki/Web_Mercator>}
(aka Pseudo-Mercator) class and conversion functions for spherical and
near-spherical earth models.

References U{The Google Maps / Bing Maps Spherical Mercator Projection
<http://alastaira.wordpress.com/2011/01/23/the-google-maps-bing-maps-spherical-mercator-projection>},
U{Geomatics Guidance Note 7, part 2<http://www.epsg.org/Portals/0/373-07-2.pdf>} and
U{Implementation Practice Web Mercator Map Projection
<http://earth-info.nga.mil/GandG/wgs84/web_mercator/%28U%29%20NGA_SIG_0011_1.0.0_WEBMERC.pdf>}.

@newfield example: Example, Examples
'''

from bases import Base
from datum import R_MA
from dms import clipDMS, parseDMS2
from ellipsoidalBase import LatLonEllipsoidalBase
from utils import EPS, PI_2, degrees90, degrees180, fStr, \
                  isscalar, map1, radians

from math import atan, atanh, exp, sin, tanh

# all public contants, classes and functions
__all__ = ('Wm',  # classes
           'parseWM', 'toWm')  # functions
__version__ = '18.01.31'

# _FalseEasting  = 0   #: (INTERNAL) False Easting (meter).
# _FalseNorthing = 0   #: (INTERNAL) False Northing (meter).
_LatLimit = 85.051129  #: (INTERNAL) Latitudinal limit (degrees).
# _LonOrigin     = 0   #: (INTERNAL) Longitude of natural origin (degrees).


class Wm(Base):
    '''Web Mercator (WM) coordinate.
    '''
    _radius = 0  #: (INTERNAL) earth radius (meter).
    _x      = 0  #: (INTERNAL) easting (meter).
    _y      = 0  #: (INTERNAL) northing (meter).

    def __init__(self, x, y, radius=R_MA):
        '''New Web Mercator (WM) coordinate.

           @param x: Easting from central meridian (meter).
           @param y: Northing from equator (meter).
           @keyword radius: Optional earth radius (meter).

           @raise ValueError: Invalid I{x}, I{y} or I{radius}.

           @example:

           >>> import pygeodesy
           >>> w = pygeodesy.Wm(448251, 5411932)
        '''
        try:
            self._radius, self._x, self._y = map1(float, radius, x, y)
        except (TypeError, ValueError):
            raise ValueError('%s invalid: %r' % (Wm.__name__, (x, y, radius)))
        if self._radius < EPS:
            raise ValueError('%s invalid: %r' % ('radius', radius))

    def parseWM(self, strWM):
        '''Parse a string to a WM coordinate.

           For more details, see function L{parseWM} in
           this module L{webmercator}.
        '''
        return parseWM(strWM, radius=self._radius)

    @property
    def radius(self):
        '''Get the earth radius (meter).
        '''
        return self._radius

    def toLatLon(self, LatLon, datum=None):
        '''Convert this WM coordinate to a geodetic point.

           @param LatLon: Geodetic class for the point (I{LatLon}).
           @keyword datum: Optional datum for ellipsoidal I{LatLon}
                           or None for spherical I{LatLon} (I{Datum}).

           @return: Point of this WM coordinate (I{LatLon}) or 2-tuple
                    (lat, lon) in (degrees) if I{LatLon} is None.

           @raise TypeError: I{LatLon} and I{datum} not compatible.

           @example:

           >>> w = Wm(448251.795, 5411932.678)
           >>> from pygeodesy import sphericalTrigonometry as sT
           >>> ll = w.toLatLon(sT.LatLon)  # 43°39′11.58″N, 004°01′36.17″E
        '''
        r = self._radius
        x = self._x / r
        y = 2 * atan(exp(self._y / r)) - PI_2
        if datum:
            E = datum.ellipsoid
            # <http://earth-info.nga.mil/GandG/wgs84/web_mercator/
            #         %28U%29%20NGA_SIG_0011_1.0.0_WEBMERC.pdf>
            y = y / r
            if E.e:
                y -= E.e * atanh(E.e * tanh(y))
            y *= E.a
            x *= E.a / r
        x, y = degrees180(x), degrees90(y)
        if not LatLon:
            return y, x

        if issubclass(LatLon, LatLonEllipsoidalBase):
            if datum:
                return LatLon(y, x, datum=datum)
        elif not datum:
            return LatLon(y, x)
        raise TypeError('%r and %s %r' % (LatLon, 'datum', datum))

    def toStr(self, prec=3, sep=' ', radius=False):  # PYCHOK expected
        '''Return a string representation of this WM coordinate.

           @keyword prec: Optional number of decimals, unstripped (int).
           @keyword sep: Optional separator to join (string).
           @keyword radius: Optionally, include radius (bool or scalar).

           @return: This WM as string "meter meter" plus " radius"
                    if I{radius} is True or scalar (string).

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
        return fStr(fs, prec=prec, sep=sep)

    def toStr2(self, prec=3, fmt='[%s]', sep=', ', radius=False):  # PYCHOK expected
        '''Return a string representation of this WM coordinate.

           @keyword prec: Optional number of decimals, unstripped (int).
           @keyword fmt: Optional, enclosing backets format (string).
           @keyword sep: Optional separator between name:value pairs (string).
           @keyword radius: Optionally, include radius (bool or scalar).

           @return: This WM as "[x:meter, y:meter]" string plus
                    ", radius:meter]" if I{radius} is True or scalar
                    (string).
        '''
        t = self.toStr(prec=prec, sep=' ', radius=radius).split()
        k = 'x', 'y', 'radius'
        return fmt % (sep.join('%s:%s' % t for t in zip(k, t)),)

    @property
    def x(self):
        '''Get the easting (meter).'''
        return self._x

    @property
    def y(self):
        '''Get the northing (meter).
        '''
        return self._y


def parseWM(strWM, radius=R_MA):
    '''Parse a string representing a WM coordinate, consisting
       of easting, northing and an optional radius.

       @param strWM: A WM coordinate (string).
       @keyword radius: Optional earth radius (meter).

       @return: The WM coordinate (L{Wm}).

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

    return Wm(x, y, radius=r)


def toWm(latlon, lon=None, radius=R_MA, Wm=Wm):
    '''Convert a lat-/longitude point to a WM coordinate.

       @param latlon: Latitude (degrees) or an (ellipsoidal or
                      spherical) geodetic I{LatLon} point.
       @keyword lon: Optional longitude (degrees or None).
       @keyword radius: Optional earth radius (meter).
       @keyword Wm: Optional Wm class for the WM coordinate (L{Wm}).

       @return: The WM coordinate (L{Wm}).

       @raise ValueError: If I{lon} value is missing, if I{latlon}
                          is not scalar or if I{latlon} is beyond
                          the valid WM range and L{rangerrors} is
                          set to True.

       @example:

       >>> p = LatLon(48.8582, 2.2945)  # 448251.8 5411932.7
       >>> w = toWm(p)  # 448252 5411933
       >>> p = LatLon(13.4125, 103.8667)  # 377302.4 1483034.8
       >>> w = toWm(p)  # 377302 1483035
    '''
    r, e = radius, None
    try:
        lat, lon = latlon.lat, latlon.lon
        if isinstance(latlon, LatLonEllipsoidalBase):
            r = latlon.datum.ellipsoid.a
            e = latlon.datum.ellipsoid.e
        lat = clipDMS(lat, _LatLimit)
    except AttributeError:
        lat, lon = parseDMS2(latlon, lon, clipLat=_LatLimit)

    s = sin(radians(lat))
    y = atanh(s)  # == log(tan(radians((90 + lat) * 0.5)))
    if e:
        y -= e * atanh(e * s)
    return Wm(r * radians(lon), r * y, radius=r)

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
