
# -*- coding: utf-8 -*-

u'''N-vector base class L{Nvector} and function L{sumOf}.

Pure Python implementation of C{n-vector}-based geodesy tools for
ellipsoidal earth models, transcribed from JavaScript originals by
I{(C) Chris Veness 2005-2016} and published under the same MIT Licence**,
see U{Vector-based geodesy
<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}.

@newfield example: Example, Examples
'''

from pygeodesy.bases import LatLonHeightBase
from pygeodesy.fmath import fsum, len2, scalar
from pygeodesy.lazily import _ALL_LAZY, _ALL_DOCS
from pygeodesy.named import Vector4Tuple, _xattrs
from pygeodesy.vector3d import Vector3d, VectorError, sumOf as _sumOf

# from math import cos, sin

# all public constants, classes and functions
__all__ = _ALL_LAZY.nvector + _ALL_DOCS('LatLonNvectorBase') + (
          'NorthPole', 'SouthPole',  # constants
          'Nvector',  # classes
          'sumOf')  # functions
__version__ = '19.07.12'


class Nvector(Vector3d):  # XXX kept private
    '''Base class for ellipsoidal and spherical L{Nvector}.
    '''
    _h = 0   #: (INTERNAL) Height (C{meter}).
    _H = ''  #: Heigth prefix (C{str}), '↑' in JS version

    def __init__(self, x, y, z, h=0, ll=None, name=''):
        '''New n-vector normal to the earth's surface.

           @param x: X component (C{scalar}).
           @param y: Y component (C{scalar}).
           @param z: Z component (C{scalar}).
           @keyword h: Optional height above surface (C{meter}).
           @keyword ll: Optional, original latlon (C{LatLon}).
           @keyword name: Optional name (C{str}).

           @example:

           >>> from pygeodesy.sphericalNvector import Nvector
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
        '''Set the height above surface.

           @param h: New height (C{meter}).

           @raise TypeError: If B{C{h}} invalid.

           @raise VectorError: If B{C{h}} invalid.
        '''
        h = scalar(h, None, name='h', Error=VectorError)
        self._update(h != self._h)
        self._h = h

    @property
    def H(self):
        '''Get the height prefix (C{str}).
        '''
        return self._H

    @H.setter  # PYCHOK setter!
    def H(self, H):
        '''Set the height prefix.

           @param H: New height prefix (C{str}).
        '''
        self._H = str(H) if H else ''

    def to3abh(self, height=None):
        '''Convert this n-vector to (geodetic) lat-, longitude
           and height.

           @keyword height: Optional height, overriding this
                            n-vector's height (C{meter}).

           @return: A L{PhiLam3Tuple}C{(phi, lambda, height)}.
        '''
        h = self.h if height is None else height
        return Vector3d.to2ab(self)._3Tuple(h)

    def to3llh(self, height=None):
        '''Convert this n-vector to (geodetic) lat-, longitude
           and height.

           @keyword height: Optional height, overriding this
                            n-vector's height (C{meter}).

           @return: A L{LatLon3Tuple}C{(lat, lon, height)}.
        '''
        return self._to3LLh(None, height)

    def _to3LLh(self, LL, height, **kwds):
        '''(INTERNAL) Helper for C{subclass.toLatLon} and C{.to3llh}.
        '''
        h = self.h if height is None else height
        r = Vector3d.to2ll(self)  # LatLon2Tuple
        if LL is None:
            r = r._3Tuple(h)  # already ._xnamed
        else:
            r = self._xnamed(LL(r.lat, r.lon, height=h, **kwds))
        return r

    def to4xyzh(self, h=None):
        '''Return this n-vector as a 4-tuple.

           @keyword h: Optional height, overriding this n-vector's
                       height (C{meter}).

           @return: A L{Vector4Tuple}C{(x, y, z, h)} in C{meter}.
        '''
        r = Vector4Tuple(self.x, self.y, self.z,
                         self.h if h is None else h)
        return self._xnamed(r)

    def toStr(self, prec=5, fmt='(%s)', sep=', '):  # PYCHOK expected
        '''Return a string representation of this n-vector.

           Height component is only included if non-zero.

           @keyword prec: Optional number of decimals, unstripped (C{int}).
           @keyword fmt: Optional enclosing backets format (C{str}).
           @keyword sep: Optional separator between components (C{str}).

           @return: Comma-separated C{"(x, y, z [, h])"} enclosed in
                    B{C{fmt}} brackets (C{str}).

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

    def unit(self, ll=None):
        '''Normalize this vector to unit length.

           @keyword ll: Optional, original latlon (C{LatLon}).

           @return: Normalized vector (L{Nvector}).
        '''
        if self._united is None:
            u = Vector3d.unit(self, ll=ll)  # .copy()
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

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        try:
            LatLonHeightBase.others(self, other, name=name)
        except TypeError:
            if not isinstance(other, Nvector):
                raise

    def to4xyzh(self, h=None):
        '''Convert this (geodetic) point to n-vector (normal
           to the earth's surface) x/y/z components and height.

           @keyword h: Optional height, overriding this point's
                       height (C{meter}).

           @return: A L{Vector4Tuple}C{(x, y, z, h)}, all in
                    (C{meter}).
        '''
        # Kenneth Gade eqn (3), but using right-handed
        # vector x -> 0°E,0°N, y -> 90°E,0°N, z -> 90°N
#       a, b = self.to2ab()
#       sa, ca, sb, cb = sincos2(a, b)
#       x, y, z = ca * cb, ca * sb, sa
        # XXX don't use self.to3xyz() + ....
        x, y, z = LatLonHeightBase.to3xyz(self)
        r = Vector4Tuple(x, y, z, self.height if h is None else h)
        return self._xnamed(r)


def sumOf(nvectors, Vector=Nvector, h=None, **kwds):
    '''Return the vectorial sum of two or more n-vectors.

       @param nvectors: Vectors to be added (L{Nvector}[]).
       @keyword Vector: Optional class for the vectorial sum (L{Nvector}).
       @keyword kwds: Optional, additional B{C{Vector}} keyword arguments.
       @keyword h: Optional height, overriding the mean height (C{meter}).

       @return: Vectorial sum (B{C{Vector}}).

       @raise VectorError: No B{C{nvectors}}.
    '''
    n, nvectors = len2(nvectors)
    if n < 1:
        raise VectorError('no nvectors: %r' & (n,))

    if h is None:
        h = fsum(v.h for v in nvectors) / float(n)
    return _sumOf(nvectors, Vector=Vector, h=h, **kwds)

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
