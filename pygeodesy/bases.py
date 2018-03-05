
# -*- coding: utf-8 -*-

u'''(INTERNAL) Common base classes.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**,
see U{http://www.movable-type.co.uk/scripts/latlong.html}
and U{http://www.movable-type.co.uk/scripts/latlong-vectors.html}.

@newfield example: Example, Examples
'''
from dms   import F_D, F_DMS, latDMS, lonDMS, parseDMS, parseDMS2
from fmath import EPS, favg, map1, scalar
from utils import R_M, antipode, classname, isantipode, polygon

from math import asin, cos, degrees, radians, sin

# XXX the following classes are listed only to get
# Epydoc to include class and method documentation
__all__ = ('Base', 'LatLonHeightBase', 'Named', 'VectorBase')
__version__ = '18.03.01'


class Base(object):
    '''(INTERNAL) Base class.
    '''
    def __repr__(self):
        return self.toStr2()

    def __str__(self):
        return self.toStr()

    def _update(self, unused):
        '''(INTERNAL) To be overloaded.
        '''
        pass

    def classname(self, *other):
        '''Build module.class name of this object.

           @param other: Optional object other than self (any).

           @return: Name of module and class (string).
        '''
        return classname(self if not other else other[0])

    def classof(self, *args, **kwds):
        '''Instantiate this very class.

           @param args: Optional, positional arguments.
           @keyword kwds: Optional, keyword arguments.

           @return: New instance (I{self.__class__}).
        '''
        return self.__class__(*args, **kwds)

#   def notImplemented(self, attr):
#       '''Raise error for a missing method, function or attribute.
#
#          @param attr: Attribute name (string).
#
#          @raise NotImplementedError: No such attribute.
#       '''
#       c = self.__class__.__name__
#       return NotImplementedError('%s.%s' % (c, attr))

    def others(self, other, name='other'):
        '''Check this and an other instance for type compatiblility.

           @param other: The other instance (any).
           @keyword name: Optional, name for other (string).

           @return: None.

           @raise TypeError: Mismatch of this and type(I{other}).
        '''
        if not (isinstance(self, other.__class__) or
                isinstance(other, self.__class__)):
            raise TypeError('type(%s) mismatch: %s vs %s' % (name,
                             self.classname(other), self.classname()))

    def toStr(self, **args):
        '''(INTERNAL) Must be overloaded.

           @param args: Optional, positional arguments.
        '''
        raise AssertionError('%s.toStr%r' % (self.__class__.__name__, args))

    def toStr2(self, **kwds):
        '''(INTERNAL) To be overloaded.

           @keyword kwds: Optional, keyword arguments.

           @return: toStr() plus keyword arguments (string).
        '''
        t = self.toStr(**kwds).lstrip('([{').rstrip('}])')
        return '%s(%s)' % (self.__class__.__name__, t)


class LatLonHeightBase(Base):
    '''(INTERNAL) Base class for I{LatLon} points on
       spherical or ellipsiodal earth models.
    '''
    _ab     = ()    #: (INTERNAL) Cache (lat, lon) radians (2-tuple)
    _height = 0     #: (INTERNAL) Height (meter)
    _lat    = 0     #: (INTERNAL) Latitude (degrees)
    _lon    = 0     #: (INTERNAL) Longitude (degrees)

    def __init__(self, lat, lon, height=0):
        '''New I{LatLon}.

           @param lat: Latitude (degrees or DMS string with N or S suffix).
           @param lon: Longitude (degrees or DMS string with E or W suffix).
           @keyword height: Optional height (meter above or below the earth surface).

           @return: New instance (I{LatLon}).

           @raise ValueError: Invalid I{lat} or I{lon}.

           @example:

           >>> p = LatLon(50.06632, -5.71475)
           >>> q = LatLon('50°03′59″N', """005°42'53"W""")
        '''
        self._lat, self._lon = parseDMS2(lat, lon)
        if height:  # elevation
            self._height = scalar(height, None, name='height')

    def __eq__(self, other):
        return self.equals(other)

    def __ne__(self, other):
        return not self.equals(other)

    def __str__(self):
        return self.toStr(form=F_D, prec=6)

    def _havg(self, other, f=0.5):
        '''(INTERNAL) Weighted, average height.

           @param other: An other point (I{LatLon}).
           @keyword f: Optional fraction (float).

           @return: Average, fractional height (float).
        '''
        return favg(self.height, other.height, f=f)

    def _update(self, updated):
        '''(INTERNAL) Reset caches if updated.
        '''
        if updated:  # reset caches
            self._ab = self._wm = None

    def antipode(self, height=None):
        '''Return the antipode, the point diametrically opposite
           to this point.

           @keyword height: Optional height of the antipode, height
                            of this point otherwise (meter).

           @return: The antipodal point (I{LatLon}).
        '''
        a, b = antipode(self.lat, self.lon)
        h = self.height if height is None else height
        return self.classof(a, b, height=h)

    def bounds(self, wide, high, radius=R_M):
        '''Return the SE and NW lat-/longitude of a great circle
           bounding box centered at this location.

           @param wide: Longitudinal box width (meter, like radius).
           @param high: Latitudinal box height (meter, like radius).
           @keyword radius: Optional, mean earth radius (meter).

           @return: 2-Tuple (LatLonSW, LatLonNE) of (I{LatLon}s).

           @see: U{http://www.movable-type.co.uk/scripts/latlong-db.html}
        '''
        a, _ = self.to2ab()
        ca = cos(a)

        if ca > EPS:
            w = abs(degrees(asin(wide * 0.5 / radius) / ca))
        else:
            w = 0  # XXX
        h = abs(degrees(high * 0.5 / radius))

        return self.classof(self.lat - h, self.lon - w, height=self.height), \
               self.classof(self.lat + h, self.lon + w, height=self.height)

    def copy(self):
        '''Copy this point.

           @return: A copy of this point (I{LatLon}).
        '''
        return self.classof(self.lat, self.lon, height=self.height)  # XXX

    def equals(self, other, eps=None):
        '''Compare this point with an other point.

           @param other: The other point (I{LatLon}).
           @keyword eps: Tolerance for equality (degrees).

           @return: True if both points are identical,
                    ignoring height (bool).

           @raise TypeError: The I{other} point is not I{LatLon}.

           @see: Method L{equals3}.

           @example:

           >>> p = LatLon(52.205, 0.119)
           >>> q = LatLon(52.205, 0.119)
           >>> e = p.equals(q)  # True
        '''
        self.others(other)

        if eps and eps > 0:
            return abs(self.lat - other.lat) < eps and \
                   abs(self.lon - other.lon) < eps
        else:
            return self.lat == other.lat and \
                   self.lon == other.lon

    def equals3(self, other, eps=None):
        '''Compare this point with an other point.

           @param other: The other point (I{LatLon}).
           @keyword eps: Tolerance for equality (degrees).

           @return: True if both points are identical,
                    I{including} height (bool).

           @raise TypeError: The I{other} point is not I{LatLon}.

           @see: Method L{equals}.

           @example:

           >>> p = LatLon(52.205, 0.119, 42)
           >>> q = LatLon(52.205, 0.119)
           >>> e = p.equals3(q)  # False
        '''
        return self.equals(other, eps=eps) and self.height == other.height

    @property
    def height(self):
        '''Get the height (meter).
        '''
        return self._height

    @height.setter  # PYCHOK setter!
    def height(self, height):
        '''Set the height.

           @param height: New height (meter).

           @raise TypeError: Invalid I{height}.

           @raise ValueError: Invalid I{height}.
        '''
        h = scalar(height, None, name='height')
        self._update(h != self._height)
        self._height = h

    def isantipode(self, other, eps=EPS):
        '''Check whether this and an other point are antipodal,
           on diametrically opposite sides of the earth.

           @param other: The other point (I{LatLon}).
           @keyword eps: Tolerance for near-equality (degrees).

           @return: True if points are antipodal within the given
                    tolerance, False otherwise.
        '''
        return isantipode(self.lat,  self.lon,
                         other.lat, other.lon, eps=eps)

    @property
    def lat(self):
        '''Get the latitude (degrees).
        '''
        return self._lat

    @lat.setter  # PYCHOK setter!
    def lat(self, lat):
        '''Set the latitude.

           @param lat: New latitude (degrees).

           @raise ValueError: Invalid I{lat}.
        '''
        lat = parseDMS(lat, suffix='NS', clip=90)
        self._update(lat != self._lat)
        self._lat = lat

    @property
    def latlon(self):
        '''Get the latitude and longitude (2-tuple of degrees).
        '''
        return self._lat, self._lon

    @latlon.setter  # PYCHOK setter!
    def latlon(self, latlonh):
        '''Set the lat- and longitude and optionally the height.

           @param latlonh: New lat-, longitude and height (2- or
                           3-tuple of degrees and meter).

           @raise TypeError: Height of I{latlonh} not scalar or
                             I{latlonh} not tuple or list.

           @raise ValueError: Invalid I{latlonh} or M{len(latlonh)}.

           @see: Function L{parse3llh} to parse a I{latlonh} string
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

    @property
    def lon(self):
        '''Get the longitude (degrees).
        '''
        return self._lon

    @lon.setter  # PYCHOK setter!
    def lon(self, lon):
        '''Set the longitude.

           @param lon: New longitude (degrees).

           @raise ValueError: Invalid I{lon}.
        '''
        lon = parseDMS(lon, suffix='EW', clip=180)
        self._update(lon != self._lon)
        self._lon = lon

    def points(self, points, closed=True):
        '''Check a polygon given as list, set or tuple of points.

           @param points: The points of the polygon (I{LatLon}[])
           @keyword closed: Optionally, treat polygon as closed (bool).

           @return: 2-Tuple (number, sequence) of points (int, sequence).

           @raise TypeError: Some I{points} are not I{LatLon}.

           @raise ValueError: Insufficient number of I{points}.
        '''
        return polygon(points, closed=closed, base=self)

    def to2ab(self):
        '''Return this point's lat-/longitude in radians.

           @return: 2-Tuple (lat, lon) in (radians, radians).
        '''
        if not self._ab:
            self._ab = map1(radians, self.lat, self.lon)
        return self._ab

    def to3llh(self):
        '''Return this point's lat-, longitude and height.

           @return: 3-Tuple (lat, lon, h) in (degrees, degrees, meter).
        '''
        return self.lat, self.lon, self.height

    def to3xyz(self):
        '''Convert this (geodetic) point to (n-)vector (normal
           to the earth's surface) x/y/z components, ignoring
           the height.

           @return: 3-Tuple (x, y, z) in (units, NOT meter).
        '''
        # Kenneth Gade eqn 3, but using right-handed
        # vector x -> 0°E,0°N, y -> 90°E,0°N, z -> 90°N
        a, b = self.to2ab()
        ca = cos(a)
        return ca * cos(b), ca * sin(b), sin(a)

    def toStr(self, form=F_DMS, prec=None, m='m', sep=', '):  # PYCHOK expected
        '''Convert this point to a "lat, lon [+/-height]" string,
           formatted in the given form.

           @keyword form: Optional format, F_D, F_DM, F_DMS for
                          deg°, deg°min′, deg°min′sec″ (string).
           @keyword prec: Optional number of decimal digits (0..8 or None).
           @keyword m: Optional unit of the height (string).
           @keyword sep: Optional separator to join (string).

           @return: Point in the specified form (string).

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


class Named(object):
    '''(INTERNAL) Named base class.
    '''
    _name = ''  #: (INTERNAL) name (string)

    def __init__(self, name):
        self._name = name

    @property
    def name(self):
        '''Gets the name (string).
        '''
        return self._name


VectorBase = Base  #: (INTERNAL) Used by vector3d.

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
