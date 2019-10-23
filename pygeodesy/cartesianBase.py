
# -*- coding: utf-8 -*-

u'''(INTERNAL) Base classes for elliposiodal, spherical and N-/vectorial
C{Cartesian}s.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**,
see U{https://www.Movable-Type.co.UK/scripts/latlong.html},
U{https://www.Movable-Type.co.UK/scripts/latlong-vectors.html} and
U{https://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html}..

@newfield example: Example, Examples
'''

from pygeodesy.datum import Datum, Datums
from pygeodesy.ecef import EcefKarney
from pygeodesy.fmath import EPS, cbrt, fsum_, _IsNotError
from pygeodesy.lazily import _ALL_DOCS
from pygeodesy.named import LatLon4Tuple, Vector4Tuple
from pygeodesy.utily import property_RO, _TypeError
from pygeodesy.vector3d import Vector3d, _xyzhdn6

from math import hypot, sqrt

# XXX the following classes are listed only to get
# Epydoc to include class and method documentation
__all__ = _ALL_DOCS('CartesianBase')
__version__ = '19.10.21'


class CartesianBase(Vector3d):
    '''(INTERNAL) Base class for ellipsoidal and spherical C{Cartesian}.
    '''
    _datum = None        #: (INTERNAL) L{Datum}, to be overriden.
    _Ecef  = EcefKarney  #: (INTERNAL) Preferred C{Ecef...} class.
    _e9t   = None        #: (INTERNAL) Cached toEcef (L{Ecef9Tuple}).
    _v4t   = None        #: (INTERNAL) Cached toNvector (L{Vector4Tuple}).

    def __init__(self, xyz, y=None, z=None, datum=None, ll=None, name=''):
        '''New C{Cartesian...}.

           @param xyz: An L{Ecef9Tuple}, L{Vector3Tuple}, L{Vector4Tuple}
                       or the C{X} coordinate (C{scalar}).
           @param y: The C{Y} coordinate (C{scalar}) if B{C{xyz}} C{scalar}.
           @param z: The C{Z} coordinate (C{scalar}) if B{C{xyz}} C{scalar}.
           @keyword datum: Optional datum (L{Datum}).
           @keyword ll: Optional, original latlon (C{LatLon}).
           @keyword name: Optional name (C{str}).

           @raise TypeError: Non-scalar B{C{xyz}}, B{C{y}} or B{C{z}}
                             coordinate or B{C{xyz}} not an L{Ecef9Tuple},
                             L{Vector3Tuple} or L{Vector4Tuple}.
        '''
        x, y, z, _, d, n = _xyzhdn6(xyz, y, z, None, datum, ll)
        Vector3d.__init__(self, x, y, z, ll=ll, name=name or n)
        if d:
            self.datum = d

    def _update(self, updated):
        if updated:  # reset cached attrs
            self._e9t = self._v4t = None
            Vector3d._update(self, updated)

    def _xcopy(self, *attrs):
        '''(INTERNAL) Make copy with add'l, subclass attributes.
        '''
        return Vector3d._xcopy(self, '_datum', '_Ecef', '_e9t', '_v4t', *attrs)

    def _applyHelmert(self, transform, inverse=False, **datum):
        '''(INTERNAL) Return a new cartesian point by applying a
           Helmert transform to this point.

           @param transform: Transform to apply (L{Transform}).
           @keyword inverse: Optionally, apply the inverse
                             Helmert transform (C{bool}).

           @return: The transformed point (C{Cartesian}).
        '''
        xyz = transform.transform(self.x, self.y, self.z, inverse)
        return self._xnamed(self.classof(xyz, **datum))

    def convertDatum(self, datum2, datum=None):
        '''Convert this cartesian point from one to an other datum.

           @param datum2: Datum to convert I{to} (L{Datum}).
           @keyword datum: Datum to convert I{from} (L{Datum}).

           @return: The converted point (C{Cartesian}).

           @raise TypeError: B{C{datum2}} or B{C{datum}} not a
                             L{Datum}.
        '''
        _TypeError(Datum, datum2=datum2)

        if datum and self.datum != datum:
            c = self.convertDatum(datum)
        else:
            c = self

        i, d = False, c.datum
        if d == datum2:
            return c.copy() if c is self else c

        elif d == Datums.WGS84:
            d = datum2  # convert from WGS 84 to datum2

        elif datum2 == Datums.WGS84:
            i = True  # conver to WGS84, use inverse transform

        else:  # neither datum2 nor c.datum is WGS84, invert to WGS84 first
            c = c._applyHelmert(d.transform, True, datum=d)
            d = datum2

        return c._applyHelmert(d.transform, i, datum=d)

    @property
    def datum(self):
        '''Get this point's datum (L{Datum}).
        '''
        return self._datum

    @datum.setter  # PYCHOK setter!
    def datum(self, datum):
        '''Set this geocentric point's C{datum} I{without conversion}.

           @param datum: New datum (L{Datum}).

           @raise TypeError: The B{C{datum}} is not a L{Datum}.
        '''
        _TypeError(Datum, datum=datum)
        if self.datum.isEllipsoidal and not datum.isEllipsoidal:
            raise _IsNotError('ellipsoidal', datum=datum)
        elif self.datum.isSpherical and not datum.isSpherical:
            raise _IsNotError('spherical', datum=datum)
        self._update(datum != self._datum)
        self._datum = datum

    @property_RO
    def Ecef(self):
        '''Get the ECEF I{class} (L{EcefKarney} or L{EcefVeness}).
        '''
        return self._Ecef

    def toEcef(self):
        '''Convert this cartesian to geodetic coordinates.

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height,
                    C, M, datum)} with C{C} and C{M} if available.

           @raise EcefError: A C{.datum} or an ECEF issue.
        '''
        if self._e9t is None:
            r = self.Ecef(self.datum).reverse(self, M=True)
            self._e9t = self._xnamed(r)
        return self._e9t

    def to3llh(self, datum=None):
        '''DEPRECATED, use method C{toLatLon}.

           Convert this cartesian to geodetic lat-, longitude and
           height.

           @keyword datum: Optional datum to use (L{Datum}).

           @return: A L{LatLon4Tuple}C{(lat, lon, height, datum)}.

           @raise TypeError: Invalid B{C{datum}}.
        '''
        t = self.toLatLon(datum=datum, LatLon=None)
        r = LatLon4Tuple(t.lat, t.lon, t.height, t.datum)
        return self._xnamed(r)

#   def _to3LLh(self, datum, LL, **pairs):  # OBSOLETE
#       '''(INTERNAL) Helper for C{subclass.toLatLon} and C{.to3llh}.
#       '''
#       r = self.to3llh(datum)  # LatLon3Tuple
#       if LL is not None:
#           r = LL(r.lat, r.lon, height=r.height, datum=datum)
#           for n, v in pairs.items():
#               setattr(r, n, v)
#           r = self._xnamed(r)
#       return r

    def toLatLon(self, datum=None, LatLon=None, **kwds):
        '''Convert this cartesian point to a geodetic point.

           @keyword datum: Optional datum (L{Datum}) or C{None}.
           @keyword LatLon: Optional (sub-)class to return the
                            geodetic point (C{LatLon}) or C{None}.
           @keyword kwds: Optional, additional B{C{LatLon}} keyword
                          arguments, ignored if C{B{LatLon}=None}.

           @return: The B{C{LatLon}} point or if C{B{LatLon}=None},
                    an L{Ecef9Tuple}C{(x, y, z, lat, lon, height,
                    C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{datum}} or B{C{kwds}}.
        '''
        c = self
        if datum is not None:
            _TypeError(Datum, datum=datum)
            if datum != self.datum:
                c = self.convertDatum(datum)

        r = c.Ecef(c.datum).reverse(c, M=True)
        if LatLon is not None:  # class or .classof
            r = LatLon(r.lat, r.lon, height=r.height,
                                      datum=c.datum, **kwds)
        return self._xnamed(r)

    def toNvector(self, Nvector=None, datum=None, **kwds):  # PYCHOK Datums.WGS84
        '''Convert this cartesian to C{n-vector} components.

           @keyword Nvector: Optional (sub-)class to return the
                             C{n-vector} components (C{Nvector})
                             or C{None}.
           @keyword datum: Optional datum (L{Datum}) overriding this
                           cartesian's datum.
           @keyword kwds: Optional, additional B{C{Nvector}} keyword
                          arguments, ignored if C{B{Nvector}=None}.

           @return: Unit vector B{C{Nvector}} or a L{Vector4Tuple}C{(x,
                    y, z, h)} if B{C{Nvector}=None}.

           @raise ValueError: The B{C{Cartesian}} at origin.
        '''
        d = datum or self.datum
        r = self._v4t
        if r is None or self.datum != d:
            E = d.ellipsoid
            x, y, z = self.to3xyz()

            # Kenneth Gade eqn 23
            p = (x**2 + y**2) * E.a2_
            q = (z**2 * E.e12) * E.a2_
            r = fsum_(p, q, -E.e4) / 6
            s = (p * q * E.e4) / (4 * r**3)
            t = cbrt(fsum_(1, s, sqrt(s * (2 + s))))

            u = r * fsum_(1, t, 1 / t)
            v = sqrt(u**2 + E.e4 * q)
            w = E.e2 * fsum_(u, v, -q) / (2 * v)

            k = sqrt(fsum_(u, v, w**2)) - w
            if abs(k) < EPS:
                raise ValueError('%s: %r' % ('origin', self))
            e = k / (k + E.e2)

            t = hypot(e * hypot(x, y), z)
            if t < EPS:
                raise ValueError('%s: %r' % ('origin', self))
            h = fsum_(k, E.e2, -1) / k * t

            s = e / t
            r = Vector4Tuple(x * s, y * s, z / t, h)
            self._v4t = r if d == self.datum else None

        if Nvector is not None:
            r = Nvector(r.x, r.y, r.z, h=r.h, datum=d, **kwds)
        return self._xnamed(r)

    def toStr(self, prec=3, fmt='[%s]', sep=', '):  # PYCHOK expected
        '''Return the string representation of this cartesian.

           @keyword prec: Optional number of decimals, unstripped (C{int}).
           @keyword fmt: Optional enclosing backets format (string).
           @keyword sep: Optional separator to join (string).

           @return: Cartesian represented as "[x, y, z]" (string).
        '''
        return Vector3d.toStr(self, prec=prec, fmt=fmt, sep=sep)

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
