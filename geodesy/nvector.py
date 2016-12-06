
# -*- coding: utf-8 -*-

# Python implementation of n-vector-based geodesy tools for an
# ellipsoidal earth model.  Transcribed from JavaScript originals by
# (C) Chris Veness 2005-2016 and published under the same MIT Licence**,
# see <http://www.movable-type.co.uk/scripts/latlong-vectors.html>

from bases import _LatLonHeightBase
from utils import fsum, len2
from vector3d import Vector3d
# from math import cos, sin

# all public constants, classes and functions
__all__ = ('NorthPole', 'Nvector', 'SouthPole',  # constants
           'sumOf')  # functions
__version__ = '16.11.11'


class Nvector(Vector3d):  # XXX kept private
    '''Base class for ellipsoidal and spherical Nvector.
    '''
    _h = 0

    H = ''  # or '↑' XXX

    def __init__(self, x, y, z, h=0):
        '''Create an n-vector normal to the earth's surface.

           @param {number} x - X component.
           @param {number} y - Y component.
           @param {number} z - Z component.
           @param {number} [h=0] - Height above surface in meter.

           @example
           from ellipsoidalNvector import Nvector
           v = Nvector(0.5, 0.5, 0.7071, 1)
           v.toLatLon()  # 45.0°N, 045.0°E, +1.00m
        '''
        Vector3d.__init__(self, x, y, z)
        if h:
            self._h = float(h)

    def copy(self):
        '''Return a copy of this vector.

           @returns {Nvector} Copy of this vector.
        '''
        n = Vector3d.copy(self)
        if n.h != self.h:
            n.h = self.h
        return n

    @property
    def h(self):
        '''Height above surface in meter.
        '''
        return self._h

    @h.setter  # PYCHOK setter!
    def h(self, h):
        '''Set height above surface in meter.
        '''
        self._update(h != self._h)
        self._h = h

    def to3latlonheight(self):
        '''Convert this n-vector to (geodetic) lat-, longitude
           and height.

           @returns {(degrees90, degrees180, meter)} 3-Tuple of
                 (lat, lon, height) equivalent to this n-vector.
        '''
        return Vector3d.to2latlon(self) + (self.h,)

    def to4tuple(self):
        '''Return this n-vector as a 4-tuple.

           @returns {(x, y, z, h)} 4-Tuple with the components of
                                   this n-vector.
        '''
        return self.x, self.y, self.z, self.h

    def toStr(self, prec=5, fmt='(%s)', sep=', '):  # PYCHOK expected
        '''Return a string representation of this n-vector.

           Height component is only included if non-zero.

           @param {number} [prec=5] - Number of decimals, unstripped.
           @param {string} [fmt='[%s]'] - Enclosing backets format.
           @param {string} [sep=', '] - Separator between components.

           @returns {string} Comma-separated x, y, z [, h] values.

           @example
           Nvector(0.5, 0.5, 0.7071).toStr()  # (0.5, 0.5, 0.7071)
           Nvector(0.5, 0.5, 0.7071, 1).toStr(-3)  # (0.500, 0.500, 0.707, +1.00)
        '''
        t = Vector3d.toStr(self, prec=prec, fmt='%s', sep=sep)
        if self.h:
            t = '%s%s%s%+.2f' % (t, sep, self.H, self.h)
        return fmt % (t,)

    def unit(self):
        '''Normalize this vector to unit length.

           @returns {Nvector} Normalised vector.
        '''
        if self._united is None:
            u = Vector3d.unit(self).copy()
            if u.h != self.h:
                u.h = self.h
            self._united = u._united = u
        return self._united

NorthPole = Nvector(0, 0, +1)
SouthPole = Nvector(0, 0, -1)


class _LatLonNvectorBase(_LatLonHeightBase):
    '''Base class for n-vector-based ellipsoidal and spherical LatLon.
    '''

    def others(self, other, name='other'):
        '''Refine class comparison.
        '''
        try:
            _LatLonHeightBase.others(self, other, name=name)
        except TypeError:
            if not isinstance(other, Nvector):
                raise

    def to4xyzh(self):
        '''Convert this (geodetic) LatLon point to n-vector (normal
           to the earth's surface) x/y/z components and height.

           @returns {(meter, meter, meter, meter)} 4-Tuple (x, y, z, h).
        '''
        # Kenneth Gade eqn (3), but using right-handed
        # vector x -> 0°E,0°N, y -> 90°E,0°N, z -> 90°N
#       a, b = self.toradians()
#       ca = cos(a)
#       x, y, z = ca * cos(b), ca * sin(b), sin(a)
        return _LatLonHeightBase.to3xyz(self) + (self.height,)


def sumOf(nvectors):
    '''Return the vectorial sum of any number of n-vectors.

       @param {Nvector[]} nvectors - The n-vectors to be added.

       @returns {Nvector} New Nvector, vectorial sum.

       @throws {ValueError} No nvectors.
    '''
    n, nvectors = len2(nvectors)
    if n < 1:
        raise ValueError('no nvectors: %r' & (n,))
    return Nvector(fsum(n.x for n in nvectors),
                   fsum(n.y for n in nvectors),
                   fsum(n.z for n in nvectors),
                 h=fsum(n.h for n in nvectors))

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
