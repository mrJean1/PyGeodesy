
# -*- coding: utf-8 -*-

# Common base classes and functions.

# After (C) Chris Veness 2011-2015 published under the same MIT Licence**,
# see <http://www.movable-type.co.uk/scripts/latlong.html>
# and <http://www.movable-type.co.uk/scripts/latlong-vectors.html>

from dms import F_D, F_DMS, latDMS, lonDMS, parseDMS
from math import cos, radians, sin

# all public contants, classes and functions
__all__ = ()  # none
__version__ = '16.11.11'


class _Base(object):

    def __repr__(self):
        return self.toStr2()

    def __str__(self):
        return self.toStr()

    def _update(self, unused):
        pass

    def notImplemented(self, attr):
        c = self.__class__.__name__
        return NotImplementedError('%s.%s' % (c, attr))

    def others(self, other, name='other'):
        '''Check mutually compatible classes.
        '''
        if not (isinstance(self, other.__class__) or
                isinstance(other, self.__class__)):
            o = other.__class__.__name__
            s = self.__class__.__name__
            raise TypeError('%s %s mismatch: %s.%s vs %s.%s' % (name,
                            o, other.__module__, o, self.__module__, s))

    def Top(self, *args, **kwds):
        '''Return a super.super... class instance.
        '''
        return self.__class__(*args, **kwds)

    def toStr(self, **unused):
        return ''

    def toStr2(self, **kwds):
        t = self.toStr(**kwds).lstrip('([{').rstrip('}])')
        return '%s(%s)' % (self.__class__.__name__, t)


_VectorBase = _Base  # used by vector3d


class _LatLonHeightBase(_Base):
    '''Base class for LatLon points on sphereical
       or ellipsiodal earth models.
    '''
    _height = 0

    def __init__(self, lat, lon, height=0):
        '''Create a new LatLon instance from the given lat-,
           longitude and height.

           @param {degrees|DMSstrNS} lat - Latitude in degrees
                           or as DMS string with N or S suffix.
           @param {degrees|DMSstrEW} lon - Longitude in degrees
                            or as DMS string with E or W suffix.
           @param {meter} [height=0] - Height above or below the
                                       earth surface in meter.

           @returns {LatLon} LatLon instance.

           @throws {ValueError} For invalid lat- or longitude
                                DMS strings or suffixes.

           @example
           p = LatLon(50.06632, -5.71475)
           q = LatLon('50°03′59″N', """005°42'53"W""")
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
        # adjust elevations
        return self.height + f * (other.height - self.height)

    def copy(self):
        '''Return a copy of this LatLon point.

           @returns {LatLon} A LatLon copy of this point.
        '''
        return self.Top(self.lat, self.lon, height=self.height)  # XXX

    def equals(self, other, eps=None):
        '''Check if this point is equal to an other point.

           @param {LatLon} other - The other point.

           @returns {bool} True if points are identical.

           @example
           p = LatLon(52.205, 0.119)
           q = LatLon(52.205, 0.119)
           e = p.equals(q)  # True
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
        '''Height in meter.
        '''
        return self._height

    @height.setter  # PYCHOK setter!
    def height(self, height):
        '''Set height in meter.
        '''
        self._update(height != self._height)
        self._height = height

    @property
    def lat(self):
        '''Latitude in degrees.
        '''
        return self._lat

    @lat.setter  # PYCHOK setter!
    def lat(self, lat):
        '''Set latitude in degrees.
        '''
        self._update(lat != self._lat)
        self._lat = lat

    @property
    def lon(self):
        '''Longitude in degrees.
        '''
        return self._lon

    @lon.setter  # PYCHOK setter!
    def lon(self, lon):
        '''Set longitude in degrees.
        '''
        self._update(lon != self._lon)
        self._lon = lon

    def toradians(self):
        '''Return this point's lat-/longitude in radians.

           @returns {(radLat, radLon)} 2-Tuple of lat and lon in radians.
        '''
        return radians(self.lat), radians(self.lon)

    def toStr(self, form=F_DMS, prec=None, m='m', sep=', '):  # PYCHOK expected
        '''Convert this point to a "lat, lon [+/-height]" string, formatted
           in the given form.

           @param {string} [form=D_DMS] - Use F_D, F_DM, F_DMS for deg°, deg°min', deg°min'sec".
           @param {number} [prec=0..8] - Number of decimal digits.
           @param {string} [m='m'] - Unit of the height, default meter
           @param {string} [sep=', '] - Separator to join.

           @returns {string} Point as string in the specified form.

           @example
           LatLon(51.4778, -0.0016).toStr()  # 51°28′40″N, 000°00′06″W
           LatLon(51.4778, -0.0016).toStr(F_D)  # 51.4778°N, 000.0016°W
           LatLon(51.4778, -0.0016, 42).toStr()  # 51°28′40″N, 000°00′06″W, +42.00m
        '''
        t = [latDMS(self.lat, form=form, prec=prec),
             lonDMS(self.lon, form=form, prec=prec)]
        if self.height:
            t += ['%+.2f%s' % (self.height, m)]
        return sep.join(t)

    def to3xyz(self):
        '''Convert this (geodetic) LatLon point to n-vector
           (normal to the earth's surface) x/y/z components.

           @returns {(meter, meter, meter)} 3-Tuple (x, y, z).
        '''
        # Kenneth Gade eqn 3, but using right-handed
        # vector x -> 0°E,0°N, y -> 90°E,0°N, z -> 90°N
        a, b = self.toradians()
        ca = cos(a)
        return ca * cos(b), ca * sin(b), sin(a)

# **) MIT License
#
# Copyright (c) 2016-2017 -- mrJean1@Gmail.com
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
