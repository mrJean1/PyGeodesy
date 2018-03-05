
# -*- coding: utf-8 -*-

u'''N-vector base class L{Nvector} and function L{sumOf}.

Pure Python implementation of I{n-vector}-based geodesy tools for
ellipsoidal earth models, transcribed from JavaScript originals by
I{(C) Chris Veness 2005-2016} and published under the same MIT Licence**,
see U{Vector-based geodesy
<http://www.movable-type.co.uk/scripts/latlong-vectors.html>}.

@newfield example: Example, Examples
'''

from bases import LatLonHeightBase
from fmath import fsum, len2, scalar
from vector3d import Vector3d, sumOf as _sumOf

# from math import cos, sin

# all public constants, classes and functions
__all__ = ('NorthPole', 'SouthPole',  # constants
           'Nvector',  # classes
           'sumOf')  # functions
__version__ = '18.03.04'


class Nvector(Vector3d):  # XXX kept private
    '''Base class for ellipsoidal and spherical L{Nvector}.
    '''
    _h = 0     #: (INTERNAL) Height (meter).

    H = ''  #: Heigth prefix (string), '↑' in JS version

    def __init__(self, x, y, z, h=0, ll=None):
        '''New n-vector normal to the earth's surface.

           @param x: X component (scalar).
           @param y: Y component (scalar).
           @param z: Z component (scalar).
           @keyword h: Optional height above surface (meter).
           @keyword ll: Optional, original latlon (I{LatLon}).

           @example:

           >>> from sphericalNvector import Nvector
           >>> v = Nvector(0.5, 0.5, 0.7071, 1)
           >>> v.toLatLon()  # 45.0°N, 045.0°E, +1.00m
        '''
        Vector3d.__init__(self, x, y, z, ll=ll)
        if h:
            self._h = scalar(h, None, name='h')

    def copy(self):
        '''Copy this vector.

           @return: Copy (L{Nvector}).
        '''
        n = Vector3d.copy(self)
        if n.h != self.h:
            n.h = self.h
        return n

    @property
    def h(self):
        '''Get the height above surface (meter).
        '''
        return self._h

    @h.setter  # PYCHOK setter!
    def h(self, h):
        '''Sets height above surface.

           @param h: New height (meter).

           @raise TypeError: If I{h} invalid.

           @raise ValueError: If I{h} invalid.
        '''
        h = scalar(h, None, name='h')
        self._update(h != self._h)
        self._h = h

    def to3llh(self):
        '''Convert this n-vector to (geodetic) lat-, longitude
           and height.

           @return: 3-Tuple (lat, lon, height) in (degrees90,
                    degrees180, meter).
        '''
        return Vector3d.to2ll(self) + (self.h,)

    def to4xyzh(self):
        '''Return this n-vector as a 4-tuple.

           @return: 4-Tuple (x, y, z, h) in (meter).
        '''
        return self.x, self.y, self.z, self.h

    def toStr(self, prec=5, fmt='(%s)', sep=', '):  # PYCHOK expected
        '''Return a string representation of this n-vector.

           Height component is only included if non-zero.

           @keyword prec: Optional number of decimals, unstripped (int).
           @keyword fmt: Optional enclosing backets format (string).
           @keyword sep: Optional separator between components (string).

           @return: Comma-separated "x, y, z [, h]" (string).

           @example:

           >>> Nvector(0.5, 0.5, 0.7071).toStr()  # (0.5, 0.5, 0.7071)
           >>> Nvector(0.5, 0.5, 0.7071, 1).toStr(-3)  # (0.500, 0.500, 0.707, +1.00)
        '''
        t = Vector3d.toStr(self, prec=prec, fmt='%s', sep=sep)
        if self.h:
            t = '%s%s%s%+.2f' % (t, sep, self.H, self.h)
        return fmt % (t,)

    def toVector3d(self):
        '''Convert this n-vector to a normalized 3-d vector,
           ignoring the height.

           @return: Normalized vector (L{Vector3d}).
        '''
        u = self.unit()
        return Vector3d(u.x, u.y, u.z)

    def unit(self):
        '''Normalize this vector to unit length.

           @return: Normalized vector (L{Nvector}).
        '''
        if self._united is None:
            u = Vector3d.unit(self)  # .copy()
            if u.h != self.h:
                u.h = self.h
            self._united = u._united = u
        return self._united


NorthPole = Nvector(0, 0, +1)  #: North pole (L{Nvector}).
SouthPole = Nvector(0, 0, -1)  #: South pole (L{Nvector}).


class LatLonNvectorBase(LatLonHeightBase):
    '''(INTERNAL) Base class for n-vector-based ellipsoidal
        and spherical LatLon.
    '''

    def others(self, other, name='other'):
        '''Refine the class comparison.

           @param other: The other point (L{LatLon}).
           @keyword name: Optional, other's name (string).

           @raise TypeError: This and type(I{other}) incompatible.
        '''
        try:
            LatLonHeightBase.others(self, other, name=name)
        except TypeError:
            if not isinstance(other, Nvector):
                raise

    def to4xyzh(self):
        '''Convert this (geodetic) point to n-vector (normal
           to the earth's surface) x/y/z components and height.

           @return: 4-Tuple (x, y, z, h) in (meter).
        '''
        # Kenneth Gade eqn (3), but using right-handed
        # vector x -> 0°E,0°N, y -> 90°E,0°N, z -> 90°N
#       a, b = self.to2ab()
#       ca = cos(a)
#       x, y, z = ca * cos(b), ca * sin(b), sin(a)
        # XXX don't use self.to3xyz() + ....
        return LatLonHeightBase.to3xyz(self) + (self.height,)


def sumOf(nvectors, Vector=Nvector, h=None, **kwds):
    '''Return the vectorial sum of any number of n-vectors.

       @param nvectors: Vectors to be added (L{Nvector}[]).
       @keyword Vector: Optional class for the vectorial sum (L{Nvector}).
       @keyword kwds: Optional, additional I{Vector} keyword argments.
       @keyword h: Optional height, overriding the mean height (meter).

       @return: Vectorial sum (I{Vector}).

       @raise ValueError: No I{nvectors}.
    '''
    n, nvectors = len2(nvectors)
    if n < 1:
        raise ValueError('no nvectors: %r' & (n,))
    if h is None:
        m = fsum(v.h for v in nvectors) / float(n)
    else:
        m = h
    return _sumOf(nvectors, Vector=Vector, h=m, **kwds)

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
