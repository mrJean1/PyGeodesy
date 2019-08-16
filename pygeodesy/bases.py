
# -*- coding: utf-8 -*-

u'''(INTERNAL) Common base classes.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**,
see U{https://www.Movable-Type.co.UK/scripts/latlong.html}
and U{https://www.Movable-Type.co.UK/scripts/latlong-vectors.html}.

@newfield example: Example, Examples
'''
from pygeodesy.dms import F_D, F_DMS, latDMS, lonDMS, parseDMS, parseDMS2
from pygeodesy.fmath import EPS, favg, map1, scalar
from pygeodesy.formy import antipode, compassAngle, equirectangular, \
                            euclidean, haversine, isantipode, vincentys
from pygeodesy.lazily import _ALL_LAZY, _ALL_DOCS
from pygeodesy.named import Bounds2Tuple, LatLon2Tuple, _NamedBase, \
                            PhiLam2Tuple, Vector3Tuple, _xattrs
from pygeodesy.utily import R_M, points2, sincos2

from math import asin, cos, degrees, radians

# XXX the following classes are listed only to get
# Epydoc to include class and method documentation
__all__ = _ALL_LAZY.bases + _ALL_DOCS('_VectorBase')
__version__ = '19.08.14'


class _VectorBase(_NamedBase):
    '''(INTERNAL) Base class for L{Vector3d}.
    '''
    def __init__(self, name='', **unused):
        if name:
            self.name = name


class LatLonHeightBase(_NamedBase):
    '''(INTERNAL) Base class for C{LatLon} points on
       spherical or ellipsiodal earth models.
    '''
    _ab     = ()    #: (INTERNAL) Cache (L{PhiLam2Tuple})
    _datum  = None  #: (INTERNAL) Datum, overriden
    _height = 0     #: (INTERNAL) Height (C{meter})
    _lat    = 0     #: (INTERNAL) Latitude (C{degrees})
    _latlon = None  #: (INTERNAL) Cache (L{LatLon2Tuple})
    _lon    = 0     #: (INTERNAL) Longitude (C{degrees})
    _name   = ''    #: (INTERNAL) name (C{str})

    def __init__(self, lat, lon, height=0, name=''):
        '''New C{LatLon}.

           @param lat: Latitude (C{degrees} or DMS C{str} with N or S suffix).
           @param lon: Longitude (C{degrees} or DMS C{str} with E or W suffix).
           @keyword height: Optional height (C{meter} above or below the earth surface).
           @keyword name: Optional name (C{str}).

           @return: New instance (C{LatLon}).

           @raise RangeError: Value of B{C{lat}} or B{C{lon}} outside the valid
                              range and C{rangerrors} set to C{True}.

           @raise ValueError: Invalid B{C{lat}} or B{C{lon}}.

           @example:

           >>> p = LatLon(50.06632, -5.71475)
           >>> q = LatLon('50°03′59″N', """005°42'53"W""")
        '''
        self._lat, self._lon = parseDMS2(lat, lon)  # PYCHOK LatLon2Tuple
        if height:  # elevation
            self._height = scalar(height, None, name='height')
        if name:
            self.name = name

    def __eq__(self, other):
        return self.isequalTo(other)

    def __ne__(self, other):
        return not self.isequalTo(other)

    def __str__(self):
        return self.toStr(form=F_D, prec=6)

    def _havg(self, other, f=0.5):
        '''(INTERNAL) Weighted, average height.

           @param other: An other point (C{LatLon}).
           @keyword f: Optional fraction (C{float}).

           @return: Average, fractional height (C{float}).
        '''
        return favg(self.height, other.height, f=f)

    def _update(self, updated):
        '''(INTERNAL) Reset caches if updated.
        '''
        if updated:  # reset caches
            self._ab = self._latlon = self._wm = None

    def _xcopy(self, *attrs):
        '''(INTERNAL) Make copy with add'l, subclass attributes.
        '''
        return _xattrs(self.classof(self.lat, self.lon),
                       self, '_height', *attrs)

    def antipode(self, height=None):
        '''Return the antipode, the point diametrically opposite
           to this point.

           @keyword height: Optional height of the antipode, height
                            of this point otherwise (C{meter}).

           @return: The antipodal point (C{LatLon}).
        '''
        a, b = antipode(self.lat, self.lon)  # PYCHOK LatLon2Tuple
        h = self.height if height is None else height
        return self.classof(a, b, height=h)

    def bounds(self, wide, high, radius=R_M):
        '''DEPRECATED, use method C{boundsOf}.
        '''
        return self.boundsOf(wide, high, radius=radius)

    def boundsOf(self, wide, high, radius=R_M):
        '''Return the SE and NW lat-/longitude of a great circle
           bounding box centered at this location.

           @param wide: Longitudinal box width (C{meter}, same units as
                        B{C{radius}} or C{degrees} if B{C{radius}} is C{None}).
           @param high: Latitudinal box height (C{meter}, same units as
                        B{C{radius}} or C{degrees} if B{C{radius}} is C{None}).
           @keyword radius: Optional, mean earth radius (C{meter}).

           @return: A L{Bounds2Tuple}C{(latlonSW, latlonNE)}, the
                    lower-left and upper-right corner (C{LatLon}).

           @see: U{https://www.Movable-Type.co.UK/scripts/latlong-db.html}
        '''
        w = wide * 0.5
        h = high * 0.5
        if radius:
            a, _ = self.to2ab()
            ca = cos(a)
            if ca > EPS:
                w = abs(degrees(asin(w / radius) / ca))
            else:
                w = 0  # XXX
            h = abs(degrees(h / radius))

        r = Bounds2Tuple(self.classof(self.lat - h, self.lon - w, height=self.height),
                         self.classof(self.lat + h, self.lon + w, height=self.height))
        return self._xnamed(r)

    def compassAngle(self, other):
        '''DEPRECATED, use method C{compassAngleTo}.
        '''
        return self.compassAngleTo(other)

    def compassAngleTo(self, other, adjust=True, wrap=False):
        '''Return the angle from North for the direction vector between
           this and an other point.

           Suitable only for short, non-near-polar vectors up to a few
           hundred Km or Miles.  Use method C{initialBearingTo} for
           larger distances.

           @param other: The other point (C{LatLon}).
           @keyword adjust: Adjust the longitudinal delta by the
                            cosine of the mean latitude (C{bool}).
           @keyword wrap: Wrap and L{unroll180} longitudes and longitudinal
                          delta (C{bool}).

           @return: Compass angle from North (C{degrees360}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @note: Courtesy Martin Schultz.

           @see: U{Local, flat earth approximation
                 <https://www.EdWilliams.org/avform.htm#flat>}.
        '''
        self.others(other)
        return compassAngle(self.lat, self.lon, other.lat, other.lon,
                            adjust=adjust, wrap=wrap)

    def copy(self):
        '''Copy this point.

           @return: The copy (C{LatLon} or subclass thereof).
        '''
        return self._xcopy()

    def _distanceTo(self, func, other, radius, **options):
        '''(INTERNAL) Helper for methods C{<func>To}.
        '''
        self.others(other)
        if radius is None:
            radius = self._datum.ellipsoid.R1 if self._datum else R_M
        return func(self.lat, self.lon, other.lat, other.lon,
                                        radius=radius, **options)

    def equals(self, other, eps=None):
        '''DEPRECATED, use method C{isequalTo}.
        '''
        return self.isequalTo(other, eps=eps)

    def equals3(self, other, eps=None):
        '''DEPRECATED, use method C{isequalTo3}.
        '''
        return self.isequalTo3(other, eps=eps)

    def equirectangularTo(self, other, radius=None, **options):
        '''Compute the distance between this and an other point
           using the U{Equirectangular Approximation / Projection
           <https://www.Movable-Type.co.UK/scripts/latlong.html>}.

           Suitable only for short, non-near-polar distances up to a
           few hundred Km or Miles.  Use method C{haversineTo} or
           C{distanceTo*} for more accurate and/or larger distances.

           See function L{equirectangular_} for more details, the
           available B{C{options}} and errors raised.

           @param other: The other point (C{LatLon}).
           @keyword radius: Optional, mean earth radius (C{meter}) or
                            C{None} for the mean radius of this
                            point's datum ellipsoid.
           @keyword options: Optional keyword arguments for function
                             L{equirectangular}.

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{equirectangular}, methods C{euclideanTo},
                 C{distanceTo*}, C{haversineTo} and C{vincentysTo}.
        '''
        return self._distanceTo(equirectangular, other, radius, **options)

    def euclideanTo(self, other, radius=None, **options):
        '''Approximate the C{Euclidian} distance between this and
           an other point.

           See function L{euclidean} for the available B{C{options}}.

           @param other: The other point (C{LatLon}).
           @keyword radius: Optional, mean earth radius (C{meter}) or
                            C{None} for the mean radius of this
                            point's datum ellipsoid.
           @keyword options: Optional keyword arguments for function
                             L{euclidean}.

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{euclidean}, methods C{equirectangularTo},
                 C{distanceTo*}, C{haversineTo} and C{vincentysTo}.
        '''
        return self._distanceTo(euclidean, other, radius, **options)

    def haversineTo(self, other, radius=None, wrap=False):
        '''Compute the distance between this and an other point using the
           U{Haversine<https://www.Movable-Type.co.UK/scripts/latlong.html>}
           formula.

           @param other: The other point (C{LatLon}).
           @keyword radius: Optional, mean earth radius (C{meter}) or
                            C{None} for the mean radius of this
                            point's datum ellipsoid.
           @keyword wrap: Wrap and L{unroll180} longitudes (C{bool}).

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{haversine}, methods C{equirectangularTo},
                 C{euclideanTo}, C{distanceTo*} and C{vincentysTo}.
        '''
        return self._distanceTo(haversine, other, radius, wrap=wrap)

    @property
    def height(self):
        '''Get the height (C{meter}).
        '''
        return self._height

    @height.setter  # PYCHOK setter!
    def height(self, height):
        '''Set the height.

           @param height: New height (C{meter}).

           @raise TypeError: Invalid B{C{height}} C{type}.

           @raise ValueError: Invalid B{C{height}}.
        '''
        h = scalar(height, None, name='height')
        self._update(h != self._height)
        self._height = h

    def isantipodeTo(self, other, eps=EPS):
        '''Check whether this and an other point are antipodal,
           on diametrically opposite sides of the earth.

           @param other: The other point (C{LatLon}).
           @keyword eps: Tolerance for near-equality (C{degrees}).

           @return: C{True} if points are antipodal within the given
                    tolerance, C{False} otherwise.
        '''
        return isantipode(self.lat,  self.lon,
                          other.lat, other.lon, eps=eps)

    def isantipode(self, other, eps=EPS):
        '''DEPRECATED, use method C{isantipodeTo}.
        '''
        return self.isantipodeTo(other, eps=eps)

    def isequalTo(self, other, eps=None):
        '''Compare this point with an other point.

           @param other: The other point (C{LatLon}).
           @keyword eps: Tolerance for equality (C{degrees}).

           @return: C{True} if both points are identical,
                    I{ignoring} height, C{False} otherwise.

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Method L{isequalTo3}.

           @example:

           >>> p = LatLon(52.205, 0.119)
           >>> q = LatLon(52.205, 0.119)
           >>> e = p.isequalTo(q)  # True
        '''
        self.others(other)

        if eps and eps > 0:
            return abs(self.lat - other.lat) < eps and \
                   abs(self.lon - other.lon) < eps
        else:
            return self.lat == other.lat and \
                   self.lon == other.lon

    def isequalTo3(self, other, eps=None):
        '''Compare this point with an other point.

           @param other: The other point (C{LatLon}).
           @keyword eps: Tolerance for equality (C{degrees}).

           @return: C{True} if both points are identical
                    I{including} height, C{False} otherwise.

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Method L{isequalTo}.

           @example:

           >>> p = LatLon(52.205, 0.119, 42)
           >>> q = LatLon(52.205, 0.119)
           >>> e = p.isequalTo3(q)  # False
        '''
        return self.isequalTo(other, eps=eps) and self.height == other.height

    @property
    def lat(self):
        '''Get the latitude (C{degrees90}).
        '''
        return self._lat

    @lat.setter  # PYCHOK setter!
    def lat(self, lat):
        '''Set the latitude.

           @param lat: New latitude (C{str[N|S]} or C{degrees}).

           @raise ValueError: Invalid B{C{lat}}.
        '''
        lat = parseDMS(lat, suffix='NS', clip=90)
        self._update(lat != self._lat)
        self._lat = lat

    @property
    def latlon(self):
        '''Get the lat- and longitude (L{LatLon2Tuple}).
        '''
        if self._latlon is None:
            self._latlon = LatLon2Tuple(self._lat, self._lon)
        return self._xrenamed(self._latlon)

    @latlon.setter  # PYCHOK setter!
    def latlon(self, latlonh):
        '''Set the lat- and longitude and optionally the height.

           @param latlonh: New lat-, longitude and height (2- or
                           3-tuple of C{degrees} and C{meter}).

           @raise TypeError: Height of B{C{latlonh}} not C{scalar} or
                             B{C{latlonh}} not C{list} or C{tuple}.

           @raise ValueError: Invalid B{C{latlonh}} or M{len(latlonh)}.

           @see: Function L{parse3llh} to parse a B{C{latlonh}} string
                 into a 3-tuple (lat, lon, h).
        '''
        if not isinstance(latlonh, (list, tuple)):
            raise TypeError('%s invalid: %r' % ('latlonh', latlonh))

        if len(latlonh) == 3:
            h = scalar(latlonh[2], None, name='latlonh')
        elif len(latlonh) != 2:
            raise ValueError('%s invalid: %r' % ('latlonh', latlonh))
        else:
            h = self._height

        lat, lon = parseDMS2(latlonh[0], latlonh[1])
        self._update(lat != self._lat or
                     lon != self._lon or h != self._height)
        self._lat, self._lon, self._height = lat, lon, h

    def latlon_(self, ndigits=0):
        '''DEPRECATED, use method C{latlon2}.
        '''
        return self.latlon2(ndigits)

    def latlon2(self, ndigits=0):
        '''Return this point's lat- and longitude, rounded.

           @keyword ndigits: Number of decimal digits (C{int}).

           @return: A L{LatLon2Tuple}C{(lat, lon)}, both rounded
                    away from zero.

           @see: Built-in function C{round}.
        '''
        r = LatLon2Tuple(round(self.lat, ndigits),
                         round(self.lon, ndigits))
        return self._xnamed(r)

    def latlon2round(self, ndigits=0):
        '''DEPRECATED, use method C{latlon2}.
        '''
        return self.latlon2(ndigits)

    @property
    def lon(self):
        '''Get the longitude (C{degrees180}).
        '''
        return self._lon

    @lon.setter  # PYCHOK setter!
    def lon(self, lon):
        '''Set the longitude.

           @param lon: New longitude (C{str[E|W]} or C{degrees}).

           @raise ValueError: Invalid B{C{lon}}.
        '''
        lon = parseDMS(lon, suffix='EW', clip=180)
        self._update(lon != self._lon)
        self._lon = lon

    def points(self, points, closed=True):
        '''DEPRECATED, use method C{points2}.
        '''
        return self.points2(points, closed=closed)

    def points2(self, points, closed=True):
        '''Check a polygon represented by points.

           @param points: The polygon points (C{LatLon}[])
           @keyword closed: Optionally, consider the polygon closed,
                            ignoring any duplicate or closing final
                            B{C{points}} (C{bool}).

           @return: 2-Tuple (number, ...) of points (C{int}, C{list} or
                    C{tuple}).

           @raise TypeError: Some B{C{points}} are not C{LatLon}.

           @raise ValueError: Insufficient number of B{C{points}}.
        '''
        return points2(points, closed=closed, base=self)

    def to2ab(self):
        '''Return this point's lat- and longitude in C{radians}.

           @return: A L{PhiLam2Tuple}C{(phi, lambda)}.
        '''
        if not self._ab:
            a, b = map1(radians, self.lat, self.lon)
            self._ab = self._xnamed(PhiLam2Tuple(a, b))
        return self._ab

    def to3llh(self, height=None):
        '''Return this point's lat-, longitude and height.

           @keyword height: Optional height, overriding this
                            point's height (C{meter}).

           @return: A L{LatLon3Tuple}C{(lat, lon, height)}.
        '''
        h = self.height if height is None else height
        return self.latlon._3Tuple(h)

    def to3xyz(self):
        '''Convert this (geodetic) point to (n-)vector (normal
           to the earth's surface) x/y/z components, ignoring
           the height.

           @return: A L{Vector3Tuple}C{(x, y, z)} in C{units},
                    NOT C{meter}.
        '''
        # Kenneth Gade eqn 3, but using right-handed
        # vector x -> 0°E,0°N, y -> 90°E,0°N, z -> 90°N
        a, b = self.to2ab()
        sa, ca, sb, cb = sincos2(a, b)
        r = Vector3Tuple(ca * cb, ca * sb, sa)
        return self._xnamed(r)

    def toStr(self, form=F_DMS, prec=None, m='m', sep=', '):  # PYCHOK expected
        '''Convert this point to a "lat, lon [+/-height]" string,
           formatted in the given form.

           @keyword form: Optional format, F_D, F_DM, F_DMS for
                          deg°, deg°min′, deg°min′sec″ (C{str}).
           @keyword prec: Optional number of decimal digits (0..8 or C{None}).
           @keyword m: Optional unit of the height (C{str}).
           @keyword sep: Optional separator to join (C{str}).

           @return: Point in the specified form (C{str}).

           @example:

           >>> LatLon(51.4778, -0.0016).toStr()  # 51°28′40″N, 000°00′06″W
           >>> LatLon(51.4778, -0.0016).toStr(F_D)  # 51.4778°N, 000.0016°W
           >>> LatLon(51.4778, -0.0016, 42).toStr()  # 51°28′40″N, 000°00′06″W, +42.00m

        '''
        t = [latDMS(self.lat, form=form, prec=prec),
             lonDMS(self.lon, form=form, prec=prec)]
        if self.height:
            t += ['%+.2f%s' % (self.height, m)]
        return sep.join(t)

    def vincentysTo(self, other, radius=None, wrap=False):
        '''Compute the distance between this and an other point using
           U{Vincenty's<https://WikiPedia.org/wiki/Great-circle_distance>}
           spherical formula.

           @param other: The other point (C{LatLon}).
           @keyword radius: Optional, mean earth radius (C{meter}) or
                            C{None} for the mean radius of this
                            point's datum ellipsoid.
           @keyword wrap: Wrap and L{unroll180} longitudes (C{bool}).

           @return: Distance (C{meter}, same units as B{C{radius}}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @see: Function L{vincentys}, methods C{equirectangularTo},
                 C{euclideanTo}, C{distanceTo*} and C{haversineTo}.
        '''
        return self._distanceTo(vincentys, other, radius, wrap=wrap)

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
