
# -*- coding: utf-8 -*-

u'''N-vector base class L{Nvector} and function L{sumOf}.

Pure Python implementation of I{n-vector}-based geodesy tools for
ellipsoidal earth models, transcribed from JavaScript originals by
I{(C) Chris Veness 2005-2016} and published under the same MIT Licence**,
see U{Vector-based geodesy
<http://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}.

@newfield example: Example, Examples
'''

from bases import LatLonHeightBase, _xattrs, _xnamed
from fmath import fsum, len2, scalar
from vector3d import Vector3d, sumOf as _sumOf

# from math import cos, sin

# all public constants, classes and functions
__all__ = ('NorthPole', 'SouthPole',  # constants
           'LatLonNvectorBase',  # for documentation
           'Nvector',  # classes
           'sumOf')  # functions
__version__ = '18.10.12'


class Nvector(Vector3d):  # XXX kept private
    '''Base class for ellipsoidal and spherical L{Nvector}.
    '''
    _h = 0     #: (INTERNAL) Height (C{meter}).

    H = ''  #: Heigth prefix (C{str}), '↑' in JS version

    def __init__(self, x, y, z, h=0, ll=None, name=''):
        '''New n-vector normal to the earth's surface.

           @param x: X component (C{scalar}).
           @param y: Y component (C{scalar}).
           @param z: Z component (C{scalar}).
           @keyword h: Optional height above surface (C{meter}).
           @keyword ll: Optional, original latlon (C{LatLon}).
           @keyword name: Optional name (C{str}).

           @example:

           >>> from sphericalNvector import Nvector
           >>> v = Nvector(0.5, 0.5, 0.7071, 1)
           >>> v.toLatLon()  # 45.0°N, 045.0°E, +1.00m
        '''
        Vector3d.__init__(self, x, y, z, ll=ll, name=name)
        if h:
            self._h = scalar(h, None, name='h')

    def _xcopy(self, *attrs):
        '''(INTERNAL) Make copy with add'l, subclass attributes.
        '''
        return Vector3d._xcopy(self, '_h', *attrs)

    def copy(self):
        '''Copy this vector.

           @return: The copy (L{Nvector} or subclass thereof).
        '''
        return self._xcopy()

    @property
    def h(self):
        '''Get the height above surface (C{meter}).
        '''
        return self._h

    @h.setter  # PYCHOK setter!
    def h(self, h):
        '''Sets height above surface.

           @param h: New height (C{meter}).

           @raise TypeError: If I{h} invalid.

           @raise ValueError: If I{h} invalid.
        '''
        h = scalar(h, None, name='h')
        self._update(h != self._h)
        self._h = h

    def to3abh(self):
        '''Convert this n-vector to (geodetic) lat-, longitude
           and height.

           @return: 3-Tuple (lat, lon, height) in (C{radians},
                    C{radians}, C{meter}).
        '''
        return Vector3d.to2ab(self) + (self.h,)

    def to3llh(self):
        '''Convert this n-vector to (geodetic) lat-, longitude
           and height.

           @return: 3-Tuple (lat, lon, height) in (C{degrees90},
                    C{degrees180}, C{meter}).
        '''
        return Vector3d.to2ll(self) + (self.h,)

    def _toLLh(self, LL, height, **kwds):
        '''(INTERNAL) Helper for I{subclass.toLatLon}.
        '''
        a, b = Vector3d.to2ll(self)
        h = self.h if height is None else height
        return (a, b, h) if LL is None else _xnamed(LL(
                a, b, height=h, **kwds), self.name)

    def to4xyzh(self):
        '''Return this n-vector as a 4-tuple.

           @return: 4-Tuple (x, y, z, h) in (C{meter}).
        '''
        return self.x, self.y, self.z, self.h

    def toStr(self, prec=5, fmt='(%s)', sep=', '):  # PYCHOK expected
        '''Return a string representation of this n-vector.

           Height component is only included if non-zero.

           @keyword prec: Optional number of decimals, unstripped (C{int}).
           @keyword fmt: Optional enclosing backets format (C{str}).
           @keyword sep: Optional separator between components (C{str}).

           @return: Comma-separated "x, y, z [, h]" (C{str}).

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
        return Vector3d(u.x, u.y, u.z, name=self.name)

    def unit(self):
        '''Normalize this vector to unit length.

           @return: Normalized vector (L{Nvector}).
        '''
        if self._united is None:
            u = Vector3d.unit(self)  # .copy()
            self._united = u._united = _xattrs(u, self, '_h')
        return self._united


NorthPole = Nvector(0, 0, +1, name='NorthPole')  #: North pole (L{Nvector}).
SouthPole = Nvector(0, 0, -1, name='SouthPole')  #: South pole (L{Nvector}).


class LatLonNvectorBase(LatLonHeightBase):
    '''(INTERNAL) Base class for n-vector-based ellipsoidal
        and spherical C{LatLon} classes.
    '''

    def others(self, other, name='other'):
        '''Refine the class comparison.

           @param other: The other point (C{LatLon}).
           @keyword name: Optional, other's name (C{str}).

           @raise TypeError: Incompatible I{other} C{type}.
        '''
        try:
            LatLonHeightBase.others(self, other, name=name)
        except TypeError:
            if not isinstance(other, Nvector):
                raise

    def to4xyzh(self):
        '''Convert this (geodetic) point to n-vector (normal
           to the earth's surface) x/y/z components and height.

           @return: 4-Tuple (x, y, z, h) in (C{meter}).
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
       @keyword kwds: Optional, additional I{Vector} keyword arguments.
       @keyword h: Optional height, overriding the mean height (C{meter}).

       @return: Vectorial sum (I{Vector}).

       @raise ValueError: No I{nvectors}.
    '''
    n, nvectors = len2(nvectors)
    if n < 1:
        raise ValueError('no nvectors: %r' & (n,))

    if h is None:
        h = fsum(v.h for v in nvectors) / float(n)
    return _sumOf(nvectors, Vector=Vector, h=h, **kwds)

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
