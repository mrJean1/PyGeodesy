
# -*- coding: utf-8 -*-

'''(INTERNAL) Common base classes and functions.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**,
see U{http://www.movable-type.co.uk/scripts/latlong.html}
and U{http://www.movable-type.co.uk/scripts/latlong-vectors.html}.

@newfield example: Example, Examples
'''

from dms import F_D, F_DMS, latDMS, lonDMS, parseDMS

from math import cos, radians, sin

# XXX the following classes are listed only to get
# Epydoc to include class and method documentation
__all__ = ('Base', 'LatLonHeightBase', 'Named', 'VectorBase')
__version__ = '17.02.09'


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

    def notImplemented(self, attr):
        '''Raise a NotImplementedError.

           @param attr: Attibute name (string).

           @raise NotImplementedError: No such L{attr}ibute.
        '''
        c = self.__class__.__name__
        return NotImplementedError('%s.%s' % (c, attr))

    def others(self, other, name='other'):
        '''Check this and an other instance for mutual compatiblility.

           @param other: The other instance (any).
           @keyword name: Other's name (string).

           @raise TypeError: Incompatible instances.
        '''
        if not (isinstance(self, other.__class__) or
                isinstance(other, self.__class__)):
            o = other.__class__.__name__
            s = self.__class__.__name__
            raise TypeError('%s %s mismatch: %s.%s vs %s.%s' % (name,
                            o, other.__module__, o, self.__module__, s))

    def topsub(self, *args, **kwds):
        '''New instance of this "top- or sub-most" class.

           @param args: Optional, positional arguments.
           @keyword kwds: Optional, keyword arguments.
        '''
        return self.__class__(*args, **kwds)

    def toStr(self, **args):
        '''(INTERNAL) Must be overloaded.
           @param args: Optional, positional arguments.
        '''
        raise AssertionError('%s.toStr%r' % (self.__class__.__name__, args))

    def toStr2(self, **kwds):
        '''(INTERNAL) To be overloaded.
           @keyword kwds: Optional, keyword arguments.
        '''
        t = self.toStr(**kwds).lstrip('([{').rstrip('}])')
        return '%s(%s)' % (self.__class__.__name__, t)


class LatLonHeightBase(Base):
    '''(INTERNAL) Base class for LatLon points on
       spherical or ellipsiodal earth models.
    '''
    _height = 0  #: (INTERNAL) Height (meter)
    _lat    = 0  #: (INTERNAL) Latitude (degrees)
    _lon    = 0  #: (INTERNAL) Longitude (degrees)

    def __init__(self, lat, lon, height=0):
        '''New LatLon.

           @param lat: Latitude (degrees or DMS string with N or S suffix).
           @param lon: Longitude (degrees or DMS string with E or W suffix).
           @keyword height: Optional height (meter above or below the earth surface).

           @return: New instance (LatLon).

           @raise ValueError: Invalid L{lat}- or L{lon}gitude.

           @example:

           >>> p = LatLon(50.06632, -5.71475)
           >>> q = LatLon('50°03′59″N', """005°42'53"W""")
        '''
        self._lat = parseDMS(lat, suffix='NS')
        self._lon = parseDMS(lon, suffix='EW')
        if height:  # elevation
            self._height = float(height)

    def __eq__(self, other):
        return self.equals(other)

    def __ne__(self, other):
        return not self.equals(other)

    def __str__(self):
        return self.toStr(form=F_D, prec=6)

    def _alter(self, other, f=0.5):
        '''(INTERNAL) Adjust elevations.
        '''
        return self.height + f * (other.height - self.height)

    def copy(self):
        '''Copy this point.

           @return: A copy of this point (LatLon).
        '''
        return self.topsub(self.lat, self.lon, height=self.height)  # XXX

    def equals(self, other, eps=None):
        '''Compare this to an other point.

           @param other: The other point (LatLon).

           @return: True if points are identical (bool).

           @example:

           >>> p = LatLon(52.205, 0.119)
           >>> q = LatLon(52.205, 0.119)
           >>> e = p.equals(q)  # True
        '''
        self.others(other)

        if eps and eps > 0:
            return max(abs(self.lat - other.lat),
                       abs(self.lon - other.lon)) < eps
        else:
            return self.lat == other.lat and \
                   self.lon == other.lon  # and \
#                  self.height == other.height

    @property
    def height(self):
        '''Get the height (meter).
        '''
        return self._height

    @height.setter  # PYCHOK setter!
    def height(self, height):
        '''Set the height.

           @param height: New height (meter).
        '''
        self._update(height != self._height)
        self._height = height

    @property
    def lat(self):
        '''Get the latitude (degrees).
        '''
        return self._lat

    @lat.setter  # PYCHOK setter!
    def lat(self, lat):
        '''Set the latitude.

           @param lat: New latitude (degrees).
        '''
        self._update(lat != self._lat)
        self._lat = lat

    @property
    def lon(self):
        '''Get the longitude (degrees).
        '''
        return self._lon

    @lon.setter  # PYCHOK setter!
    def lon(self, lon):
        '''Set the longitude.

           @param lon: New longitude (degrees).
        '''
        self._update(lon != self._lon)
        self._lon = lon

    def toradians(self):
        '''Return this point's lat-/longitude in radians.

           @return: 2-Tuple (lat, lon) in (radians, ...).
        '''
        return radians(self.lat), radians(self.lon)

    def toStr(self, form=F_DMS, prec=None, m='m', sep=', '):  # PYCHOK expected
        '''Convert this point to a "lat, lon [+/-height]" string,
           formatted in the given form.

           @keyword form: Use F_D, F_DM, F_DMS for deg°, deg°min', deg°min'sec" (string).
           @keyword prec: Number of decimal digits (0..8 or None).
           @keyword m: Unit of the height (string).
           @keyword sep: Separator to join (string).

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

    def to3xyz(self):
        '''Convert this (geodetic) point to n-vector (normal
           to the earth's surface) x/y/z components.

           @return: 3-Tuple (x, y, z) in (meter, ...).
        '''
        # Kenneth Gade eqn 3, but using right-handed
        # vector x -> 0°E,0°N, y -> 90°E,0°N, z -> 90°N
        a, b = self.toradians()
        ca = cos(a)
        return ca * cos(b), ca * sin(b), sin(a)


class Named(object):
    '''(INTERNAL) Named base class.
    '''
    _name = ''  #: (INTERNAL) name (string)

    def __init__(self, name):
        self._name = name

    @property
    def name(self):
        '''Get the name (string).
        '''
        return self._name


VectorBase = Base  #: (INTERNAL) Used by vector3d.

# **) MIT License
#
# Copyright (C) 2016-2017 -- mrJean1@Gmail.com
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
