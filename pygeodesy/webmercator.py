
# -*- coding: utf-8 -*-

u'''Web Mercator (WM) class L{Wm} and functions L{parseWM} and L{toWm}.

Pure Python implementation of a U{Web Mercator<http://wikipedia.org/wiki/Web_Mercator>}
(aka I{Pseudo-Mercator}) class and conversion functions for spherical and
near-spherical earth models.

References U{Google Maps / Bing Maps Spherical Mercator Projection
<http://alastaira.wordpress.com/2011/01/23/the-google-maps-bing-maps-spherical-mercator-projection>},
U{Geomatics Guidance Note 7, part 2<http://www.epsg.org/Portals/0/373-07-2.pdf>} and
U{Implementation Practice Web Mercator Map Projection
<http://earth-info.nga.mil/GandG/wgs84/web_mercator/%28U%29%20NGA_SIG_0011_1.0.0_WEBMERC.pdf>}.

@newfield example: Example, Examples
'''

from bases import Base
from datum import R_MA
from dms import clipDMS, parseDMS2
from ellipsoidalBase import LatLonEllipsoidalBase as _eLLb
from fmath import EPS, fStr, isscalar, map1
from utils import PI_2, classname, degrees90, degrees180, radians

from math import atan, atanh, exp, sin, tanh

# all public contants, classes and functions
__all__ = ('Wm',  # classes
           'parseWM', 'toWm')  # functions
__version__ = '18.02.02'

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
            self._x, self._y, self._radius = map1(float, x, y, radius)
        except (TypeError, ValueError):
            raise ValueError('%s invalid: %r' % (Wm.__name__, (x, y, radius)))
        _ = self.radius  # PYCHOK check radius

    def parseWM(self, strWM):
        '''Parse a string to a WM coordinate.

           For more details, see function L{parseWM} in
           this module L{webmercator}.
        '''
        return parseWM(strWM, radius=self.radius)

    @property
    def radius(self):
        '''Get the earth radius (meter).
        '''
        r = self._radius
        if r < EPS:
            t = '%s.%s' % (classname(self), 'radius')
            raise ValueError('%s invalid: %r' % (t, r))
        return r

    def to2ll(self, datum=None):
        '''Convert this WM coordinate to a geodetic lat- and longitude.

           @keyword datum: Optional datum (I{Datum}).

           @return: 2-Tuple (lat, lon) in (degrees90, degrees180).

           @raise TypeError: If I{datum} is not ellipsoidal.

           @raise ValueError: Invalid I{radius}.
        '''
        r = self.radius
        x = self._x / r
        y = 2 * atan(exp(self._y / r)) - PI_2
        if datum:
            E = datum.ellipsoid
            if not E.isEllipsoidal:
                raise TypeError('%s not %s: %r' % ('datum', 'ellipsoidal', datum))
            # <http://earth-info.nga.mil/GandG/wgs84/web_mercator/
            #         %28U%29%20NGA_SIG_0011_1.0.0_WEBMERC.pdf>
            y = y / r
            if E.e:
                y -= E.e * atanh(E.e * tanh(y))
            y *= E.a
            x *= E.a / r
        return degrees90(y), degrees180(x)

    def toLatLon(self, LatLon, datum=None):
        '''Convert this WM coordinate to a geodetic point.

           @param LatLon: LatLon class to use for the point (I{LatLon}).
           @keyword datum: Optional datum for ellipsoidal or None for
                           spherical I{LatLon} (I{Datum}).

           @return: Point of this WM coordinate (I{LatLon}).

           @raise TypeError: If I{LatLon} and I{datum} are not compatible
                             or if I{datum} is not ellipsoidal.

           @raise ValueError: Invalid I{radius}.

           @example:

           >>> w = Wm(448251.795, 5411932.678)
           >>> from pygeodesy import sphericalTrigonometry as sT
           >>> ll = w.toLatLon(sT.LatLon)  # 43°39′11.58″N, 004°01′36.17″E
        '''
        if issubclass(LatLon, _eLLb):
            if datum:
                return LatLon(*self.to2ll(datum=datum), datum=datum)
        elif datum is None:
            return LatLon(*self.to2ll(datum=datum))
        raise TypeError('%r and %s %r' % (LatLon, 'datum', datum))

    def toStr(self, prec=3, sep=' ', radius=False):  # PYCHOK expected
        '''Return a string representation of this WM coordinate.

           @keyword prec: Optional number of decimals, unstripped (int).
           @keyword sep: Optional separator to join (string).
           @keyword radius: Optionally, include radius (bool or scalar).

           @return: This WM as string "meter meter" plus " radius"
                    if I{radius} is True or scalar (string).

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


def parseWM(strWM, radius=R_MA, Wm=Wm):
    '''Parse a string representing a WM coordinate, consisting
       of easting, northing and an optional radius.

       @param strWM: A WM coordinate (string).
       @keyword radius: Optional earth radius (meter).
       @keyword Wm: Optional Wm class to use (L{Wm}) or None.

       @return: The WM coordinate (L{Wm}) or 3-tuple (easting,
                northing, radius) if I{Wm} is None.

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

    return (x, y, r) if Wm is None else Wm(x, y, radius=r)


def toWm(latlon, lon=None, radius=R_MA, Wm=Wm):
    '''Convert a lat-/longitude point to a WM coordinate.

       @param latlon: Latitude (degrees) or an (ellipsoidal or
                      spherical) geodetic I{LatLon} point.
       @keyword lon: Optional longitude (degrees or None).
       @keyword radius: Optional earth radius (meter).
       @keyword Wm: Optional Wm class for the WM coordinate (L{Wm}).

       @return: The WM coordinate (L{Wm}).

       @raise ValueError: If I{lon} value is missing, if I{latlon}
                          is not scalar, if I{latlon} is beyond
                          the valid WM range and L{rangerrors} is
                          set to True or if I{radius} is invalid.

       @example:

       >>> p = LatLon(48.8582, 2.2945)  # 448251.8 5411932.7
       >>> w = toWm(p)  # 448252 5411933
       >>> p = LatLon(13.4125, 103.8667)  # 377302.4 1483034.8
       >>> w = toWm(p)  # 377302 1483035
    '''
    r, e = radius, None
    try:
        lat, lon = latlon.lat, latlon.lon
        if isinstance(latlon, _eLLb):
            r = latlon.datum.ellipsoid.a
            e = latlon.datum.ellipsoid.e
        lat = clipDMS(lat, _LatLimit)
    except AttributeError:
        lat, lon = parseDMS2(latlon, lon, clipLat=_LatLimit)

    s = sin(radians(lat))
    y = atanh(s)  # == log(tan((90 + lat) / 2)) == log(tanPI_2_2(radians(lat)))
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
