
# -*- coding: utf-8 -*-

u'''(INTERNAL) Base classes L{LatLonNvectorBase} and L{NvectorBase}
and function L{sumOf} for C{N-vectorial} ellipsoidal and spherical
C{Cartesian}s and C{LatLon}s.

Pure Python implementation of C{n-vector}-based geodesy tools for
ellipsoidal earth models, transcribed from JavaScript originals by
I{(C) Chris Veness 2005-2016} and published under the same MIT Licence**,
see U{Vector-based geodesy
<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}.

@newfield example: Example, Examples
'''

from pygeodesy.ecef import EcefVeness
from pygeodesy.fmath import fsum, hypot_, len2, scalar
from pygeodesy.latlonBase import LatLonBase
from pygeodesy.lazily import _ALL_DOCS, _2kwds
from pygeodesy.named import LatLon3Tuple, Vector3Tuple, \
                            Vector4Tuple, _xattrs
from pygeodesy.vector3d import Vector3d, VectorError, \
                               sumOf as _sumOf, _xyzhdn6
from pygeodesy.utily import property_RO

# from math import cos, sin

# all public constants, classes and functions
__all__ = _ALL_DOCS('LatLonNvectorBase') + (
          'NorthPole', 'SouthPole',  # constants
          'NvectorBase',  # classes
          'sumOf')  # functions
__version__ = '19.10.21'


class NvectorBase(Vector3d):  # XXX kept private
    '''Base class for ellipsoidal and spherical C{Nvector}s.
    '''
    _datum = None        #: (INTERNAL) L{Datum}, overriden.
    _Ecef  = EcefVeness  #: (INTERNAL) Preferred C{Ecef...} class, backward compatible.
    _e9t   = None        #: (INTERNAL) Cached toCartesian (L{Ecef9Tuple}).
    _h     = 0           #: (INTERNAL) Height (C{meter}).
    _H     = ''          #: Heigth prefix (C{str}), '↑' in JS version

    def __init__(self, x, y=None, z=None, h=0, ll=None, datum=None, name=''):
        '''New n-vector normal to the earth's surface.

           @param x: An C{Nvector}, L{Vector3Tuple}, L{Vector4Tuple} or
                     the C{X} coordinate (C{scalar}).
           @param y: The C{Y} coordinate (C{scalar}) if B{C{x}} C{scalar}.
           @param z: The C{Z} coordinate (C{scalar}) if B{C{x}} C{scalar}.
           @keyword h: Optional height above surface (C{meter}).
           @keyword ll: Optional, original latlon (C{LatLon}).
           @keyword datum: Optional, I{pass-thru} datum (C{Datum}).
           @keyword name: Optional name (C{str}).

           @raise TypeError: Non-scalar B{C{x}}, B{C{y}} or B{C{z}}
                             coordinate or B{C{x}} not an C{Nvector},
                             L{Vector3Tuple} or L{Vector4Tuple}.

           @example:

           >>> from pygeodesy.sphericalNvector import Nvector
           >>> v = Nvector(0.5, 0.5, 0.7071, 1)
           >>> v.toLatLon()  # 45.0°N, 045.0°E, +1.00m
        '''
        x, y, z, h, d, n = _xyzhdn6(x, y, z, h, datum, ll)
        Vector3d.__init__(self, x, y, z, ll=ll, name=name or n)
        if h:
            self.h = h
        if d:  # just pass-thru
            self._datum = d

    def _xcopy(self, *attrs):
        '''(INTERNAL) Make copy with add'l, subclass attributes.
        '''
        return Vector3d._xcopy(self, '_h', *attrs)

    def copy(self):
        '''Copy this vector.

           @return: The copy (C{Nvector} or subclass thereof).
        '''
        return self._xcopy()

    @property_RO
    def datum(self):
        '''Get the I{pass-thru} datum (C{Datum}) or C{None}.
        '''
        return self._datum

    @property_RO
    def Ecef(self):
        '''Get the ECEF I{class} (L{EcefKarney} or L{EcefVeness}).
        '''
        return self._Ecef

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
           in C{radians} and height.

           @keyword height: Optional height, overriding this
                            n-vector's height (C{meter}).

           @return: A L{PhiLam3Tuple}C{(phi, lam, height)}.
        '''
        h = self.h if height is None else height
        return Vector3d.to2ab(self)._3Tuple(h)

    def to3llh(self, height=None):
        '''Convert this n-vector to (geodetic) lat-, longitude
           in C{degrees} and height.

           @keyword height: Optional height, overriding this
                            n-vector's height (C{meter}).

           @return: A L{LatLon3Tuple}C{(lat, lon, height)}.
        '''
        r = self.toLatLon(height=height, LatLon=None)
        r = LatLon3Tuple(r.lat, r.lon, r.height)
        return self._xnamed(r)

    def to4xyzh(self, h=None):
        '''Return this n-vector's components as 4-tuple.

           @keyword h: Optional height, overriding this n-vector's
                       height (C{meter}).

           @return: A L{Vector4Tuple}C{(x, y, z, h)}.
        '''
        r = Vector4Tuple(self.x, self.y, self.z,
                         self.h if h is None else h)
        return self._xnamed(r)

    def toCartesian(self, h=None, Cartesian=None, datum=None, **kwds):
        '''Convert this n-vector to C{Nvector}-based cartesian (ECEF)
           coordinates.

           @keyword height: Optional height, overriding this n-vector's
                            height (C{meter}).
           @keyword Cartesian: Optional (sub-)class to return the
                               (ECEF)coordinates (L{Cartesian}).
           @keyword datum: Optional, spherical datum (C{Datum}).
           @keyword kwds: Optional, additional C{name=value} pairs
                          for B{C{Cartesian}} instance, provided
                          B{C{Cartesian}} is not C{None}.

           @return: Cartesian (ECEF) coordinates (B{C{Cartesian}}).

           @raise TypeError: Invalid B{C{Cartesian}}.

           @example:

           >>> v = Nvector(0.5, 0.5, 0.7071)
           >>> c = v.toCartesian()  # [3194434, 3194434, 4487327]
           >>> p = c.toLatLon()  # 45.0°N, 45.0°E
        '''
        x, y, z = self.x, self.y, self.z
        if h is None:
            h = self.h
        d = datum or self.datum

        E = d.ellipsoid
        # Kenneth Gade eqn (22)
        n = E.b / hypot_(x * E.a_b, y * E.a_b, z)
        r = h + n * E.a_b**2

        c = self.Ecef(d).reverse(x * r, y * r, z * (n + h), M=True)
        if Cartesian is not None:  # class or .classof
            c = Cartesian(c, **kwds)
        return self._xnamed(c)

    def toLatLon(self, height=None, LatLon=None, datum=None, **kwds):
        '''Convert this n-vector to an C{Nvector}-based geodetic point.

           @keyword height: Optional height, overriding this n-vector's
                            height (C{meter}).
           @keyword LatLon: Optional (sub-)class to return the
                            point (L{LatLon}) or C{None}.
           @keyword datum: Optional, spherical datum (C{Datum}).
           @keyword kwds: Optional, additional C{name=value} pairs
                          for B{C{LatLon}} instance, provided
                          B{C{LatLon}} is not C{None}.

           @return: The B{C{LatLon}} point (L{LatLon}) or if
                    C{B{LatLon}=None} or a L{LatLon3Tuple}C{(lat,
                    lon, height)} if B{C{LatLon}} is C{None}.

           @raise TypeError: Invalid B{C{LatLon}}.

           @example:

           >>> v = Nvector(0.5, 0.5, 0.7071)
           >>> p = v.toLatLon()  # 45.0°N, 45.0°E
        '''
        h = self.h if height is None else height
        d = datum or self.datum

        # XXX use self.Cartesian(Cartesian=None) if h == self.h
        # and d == self.datum, for better accuracy of the height
        r = self.Ecef(d).forward(Vector3d.to2ll(self), height=h, M=True)
        if LatLon is not None:  # class or .classof
            r = LatLon(r.lat, r.lon, r.height, datum=r.datum, **kwds)
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
           I{ignoring the height}.

           @return: Normalized vector (L{Vector3d}).
        '''
        u = self.unit()
        return Vector3d(u.x, u.y, u.z, name=self.name)

    def unit(self, ll=None):
        '''Normalize this vector to unit length.

           @keyword ll: Optional, original latlon (C{LatLon}).

           @return: Normalized vector (C{Nvector}).
        '''
        if self._united is None:
            u = Vector3d.unit(self, ll=ll)  # .copy()
            self._united = u._united = _xattrs(u, self, '_h')
        return self._united


NorthPole = NvectorBase(0, 0, +1, name='NorthPole')  #: North pole (C{Nvector}).
SouthPole = NvectorBase(0, 0, -1, name='SouthPole')  #: South pole (C{Nvector}).


class LatLonNvectorBase(LatLonBase):
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
            LatLonBase.others(self, other, name=name)
        except TypeError:
            if not isinstance(other, NvectorBase):
                raise

    def toNvector(self, **kwds):  # PYCHOK signature
        '''Convert this point to C{Nvector} components, I{including
           height}.

           @keyword kwds: Optional, additional B{C{Nvector}} keyword
                          arguments, ignored if C{B{Nvector}=None}.
                          Specify C{Nvector=...} to override this
                          C{Nvector} class or set C{B{Nvector}=None}.

           @return: The B{C{Nvector}} components (C{Nvector}) or a
                    L{Vector4Tuple}C{(x, y, z, h)} if C{B{Nvector}=None}.

           @raise TypeError: Invalid B{C{Nvector}} or B{C{kwds}}.
        '''
        kwds = _2kwds(kwds, Nvector=NvectorBase)
        return LatLonBase.toNvector(self, **kwds)


def sumOf(nvectors, Vector=None, h=None, **kwds):
    '''Return the vectorial sum of two or more n-vectors.

       @param nvectors: Vectors to be added (C{Nvector}[]).
       @keyword Vector: Optional class for the vectorial sum
                        (C{Nvector}) or C{None}.
       @keyword h: Optional height, overriding the mean height (C{meter}).
       @keyword kwds: Optional, additional B{C{Vector}} keyword arguments,
                      ignored if C{B{Vector}=None}.

       @return: Vectorial sum (B{C{Vector}}) or a L{Vector4Tuple}C{(x, y,
                z, h)} if C{B{Vector}=None}.

       @raise VectorError: No B{C{nvectors}}.
    '''
    n, nvectors = len2(nvectors)
    if n < 1:
        raise VectorError('no nvectors: %r' & (n,))

    if h is None:
        h = fsum(v.h for v in nvectors) / float(n)

    if Vector is None:
        r = _sumOf(nvectors, Vector=Vector3Tuple)._to4Tuple(h)
    else:
        r = _sumOf(nvectors, Vector=Vector, h=h, **kwds)
    return r

# **) MIT License
#
# Copyright (C) 2016-2020 -- mrJean1 at Gmail -- All Rights Reserved.
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
