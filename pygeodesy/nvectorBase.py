
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

from pygeodesy.basics import len2, property_doc_, \
                             property_RO, _xattrs
from pygeodesy.ecef import EcefVeness
from pygeodesy.errors import _xkwds_pop
from pygeodesy.fmath import fsum, hypot_
from pygeodesy.formy import n_xyz2latlon, n_xyz2philam
from pygeodesy.interns import _COMMA_SPACE_, _h_, _Missing, NN, \
                              _NorthPole_, _other_, _PARENTH_, \
                              _SouthPole_, _sumOf_
from pygeodesy.latlonBase import LatLonBase
from pygeodesy.lazily import _ALL_DOCS
from pygeodesy.named import Vector3Tuple, Vector4Tuple
from pygeodesy.streprs import hstr
from pygeodesy.units import Height
from pygeodesy.vector3d import Vector3d, VectorError, \
                               sumOf as _sumOf, _xyzhdn6

# from math import atan2, cos, sin

__all__ = (_NorthPole_, _SouthPole_,  # constants
           _sumOf_)  # functions
__version__ = '20.07.24'


class NvectorBase(Vector3d):  # XXX kept private
    '''Base class for ellipsoidal and spherical C{Nvector}s.
    '''
    _datum  = None        #: (INTERNAL) L{Datum}, overriden.
    _Ecef   = EcefVeness  #: (INTERNAL) Preferred C{Ecef...} class, backward compatible.
    _h      = 0           #: (INTERNAL) Height (C{meter}).
    _H      = NN          #: Heigth prefix (C{str}), '↑' in JS version
    _latlon = None        #: (INTERNAL) Cached latlon (L{LatlLon2Tuple}).
    _philam = None        #: (INTERNAL) Cached philam (L{PhiLam2Tuple}).

    def __init__(self, x, y=None, z=None, h=0, ll=None, datum=None, name=NN):
        '''New n-vector normal to the earth's surface.

           @arg x: An C{Nvector}, L{Vector3Tuple}, L{Vector4Tuple} or
                     the C{X} coordinate (C{scalar}).
           @arg y: The C{Y} coordinate (C{scalar}) if B{C{x}} C{scalar}.
           @arg z: The C{Z} coordinate (C{scalar}) if B{C{x}} C{scalar}.
           @kwarg h: Optional height above surface (C{meter}).
           @kwarg ll: Optional, original latlon (C{LatLon}).
           @kwarg datum: Optional, I{pass-thru} datum (C{Datum}).
           @kwarg name: Optional name (C{str}).

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

    def _update(self, updated, *attrs):
        '''(INTERNAL) Zap cached attributes if updated.
        '''
        if updated:
            Vector3d._update(self, updated, '_latlon', '_philam', *attrs)

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

    @property_doc_(''' the height above surface (C{meter}).''')
    def h(self):
        '''Get the height above surface (C{meter}).
        '''
        return self._h

    @h.setter  # PYCHOK setter!
    def h(self, h):
        '''Set the height above surface.

           @arg h: New height (C{meter}).

           @raise TypeError: If B{C{h}} invalid.

           @raise VectorError: If B{C{h}} invalid.
        '''
        h = Height(h, name=_h_, Error=VectorError)
        self._update(h != self._h, '_latlon', '_philam')
        self._h = h

    @property_doc_(''' the height prefix (C{str}).''')
    def H(self):
        '''Get the height prefix (C{str}).
        '''
        return self._H

    @H.setter  # PYCHOK setter!
    def H(self, H):
        '''Set the height prefix.

           @arg H: New height prefix (C{str}).
        '''
        self._H = str(H) if H else NN

    def hStr(self, prec=-2, m=NN):
        '''Return a string for the height B{C{h}}.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg m: Optional unit of the height (C{str}).

           @see: Function L{hstr}.
        '''
        return self.H + hstr(self.h, prec=prec, m=m)

    @property_RO
    def isEllipsoidal(self):
        '''Check whether this n-vector is ellipsoidal (C{bool} or C{None} if unknown).
        '''
        return self.datum.isEllipsoidal if self._datum else None

    @property_RO
    def isSpherical(self):
        '''Check whether this n-vector is spherical (C{bool} or C{None} if unknown).
        '''
        return self.datum.isSpherical if self._datum else None

    @property_RO
    def lam(self):
        '''Get the (geodetic) longitude in C{radians} (C{float}).
        '''
        return self.philamheight.lam

    @property_RO
    def lat(self):
        '''Get the (geodetic) latitude in C{degrees} (C{float}).
        '''
        return self.latlonheight.lat

    @property_RO
    def latlon(self):
        '''Get the (geodetic) lat-, longitude in C{degrees} (L{LatLon2Tuple}C{(lat, lon)}).
        '''
        if self._latlon is None:
            self._latlon = n_xyz2latlon(self.x, self.y, self.z)
        return self._xnamed(self._latlon)

    @property_RO
    def latlonheight(self):
        '''Get the (geodetic) lat-, longitude in C{degrees} and height (L{LatLon3Tuple}C{(lat, lon, height)}).
        '''
        return self._xnamed(self.latlon.to3Tuple(self.h))

    @property_RO
    def lon(self):
        '''Get the (geodetic) longitude in C{degrees} (C{float}).
        '''
        return self.latlonheight.lon

    @property_RO
    def phi(self):
        '''Get the (geodetic) latitude in C{radians} (C{float}).
        '''
        return self.philamheight.phi

    @property_RO
    def philam(self):
        '''Get the (geodetic) lat-, longitude in C{radians} (L{PhiLam2Tuple}C{(phi, lam)}).
        '''
        if self._philam is None:
            self._philam = n_xyz2philam(self.x, self.y, self.z)
        return self._xnamed(self._philam)

    @property_RO
    def philamheight(self):
        '''Get the (geodetic) lat-, longitude in C{radians} and height (L{PhiLam3Tuple}C{(phi, lam, height)}).
        '''
        return self._xnamed(self.philam.to3Tuple(self.h))

    def to2ab(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{philam}.

           @return: A L{PhiLam2Tuple}C{(phi, lam)}.
        '''
        return self.philam

    def to3abh(self, height=None):  # PYCHOK no cover
        '''DEPRECATED, use method C{philamheight} or C{philam.to3Tuple}C{(}B{C{height}}C{)}.

           @kwarg height: Optional height, overriding this
                          n-vector's height (C{meter}).

           @return: A L{PhiLam3Tuple}C{(phi, lam, height)}.

           @raise ValueError: Invalid B{C{height}}.
        '''
        return self.philamheight if height in (None, self.h) else \
               self.philam.to3Tuple(height)

    def toCartesian(self, h=None, Cartesian=None, datum=None, **Cartesian_kwds):
        '''Convert this n-vector to C{Nvector}-based cartesian (ECEF)
           coordinates.

           @kwarg h: Optional height, overriding this n-vector's height (C{meter}).
           @kwarg Cartesian: Optional class to return the (ECEF)
                             coordinates (L{Cartesian}).
           @kwarg datum: Optional, spherical datum (C{Datum}).
           @kwarg Cartesian_kwds: Optional, additional B{C{Cartesian}}
                                  keyword arguments, ignored if
                                  B{C{Cartesian=None}}.

           @return: The cartesian (ECEF) coordinates (B{C{Cartesian}}) or
                    if B{C{Cartesian}} is C{None}, an L{Ecef9Tuple}C{(x, y,
                    z, lat, lon, height, C, M, datum)} with C{C} and C{M}
                    if available.

           @raise TypeError: Invalid B{C{Cartesian}}.

           @raise ValueError: Invalid B{C{h}}.

           @example:

           >>> v = Nvector(0.5, 0.5, 0.7071)
           >>> c = v.toCartesian()  # [3194434, 3194434, 4487327]
           >>> p = c.toLatLon()  # 45.0°N, 45.0°E
        '''
        x, y, z = self.x, self.y, self.z

        h = self.h if h is None else Height(h, name=_h_)
        d = datum or self.datum

        E = d.ellipsoid
        # Kenneth Gade eqn (22)
        n = E.b / hypot_(x * E.a_b, y * E.a_b, z)
        r = h + n * E.a_b**2

        c = self.Ecef(d).reverse(x * r, y * r, z * (n + h), M=True)
        if Cartesian is not None:  # class or .classof
            c = Cartesian(c, **Cartesian_kwds)
        return self._xnamed(c)

    def to2ll(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{latlon}.

           @return: A L{LatLon2Tuple}C{(lat, lon)}.
        '''
        return self.latlon

    def to3llh(self, height=None):  # PYCHOK no cover
        '''DEPRECATED, use property C{latlonheight} or C{latlon.to3Tuple}C{)}B{C{height}}C{)}.

           @kwarg height: Optional height, overriding this
                          n-vector's height (C{meter}).

           @return: A L{LatLon3Tuple}C{(lat, lon, height)}.

           @raise ValueError: Invalid B{C{height}}.
        '''
        return self.latlonheight if height in (None, self.h) else \
               self.latlon.to3Tuple(height)

    def toLatLon(self, height=None, LatLon=None, datum=None, **LatLon_kwds):
        '''Convert this n-vector to an C{Nvector}-based geodetic point.

           @kwarg height: Optional height, overriding this n-vector's
                          height (C{meter}).
           @kwarg LatLon: Optional class to return the geodetic point
                          (L{LatLon}) or C{None}.
           @kwarg datum: Optional, spherical datum (C{Datum}).
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                               arguments, ignored if B{C{LatLon=None}}.

           @return: The geodetic point (L{LatLon}) or if B{C{LatLon}} is
                    is C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon,
                    height, C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{LatLon}}.

           @raise ValueError: Invalid B{C{height}}.

           @example:

           >>> v = Nvector(0.5, 0.5, 0.7071)
           >>> p = v.toLatLon()  # 45.0°N, 45.0°E
        '''
        # use self.Cartesian(Cartesian=None) if h == self.h and
        # d == self.datum, for better accuracy of the height
        h = self.h if height is None else Height(height)
        r = self.Ecef(datum or self.datum).forward(self.latlon, height=h, M=True)
        if LatLon is not None:  # class or .classof
            r = LatLon(r.lat, r.lon, r.height, datum=r.datum, **LatLon_kwds)
        return self._xnamed(r)

    def toStr(self, prec=5, fmt=_PARENTH_, sep=_COMMA_SPACE_):  # PYCHOK expected
        '''Return a string representation of this n-vector.

           Height component is only included if non-zero.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg fmt: Optional enclosing backets format (C{str}).
           @kwarg sep: Optional separator between components (C{str}).

           @return: Comma-separated C{"(x, y, z [, h])"} enclosed in
                    B{C{fmt}} brackets (C{str}).

           @example:

           >>> Nvector(0.5, 0.5, 0.7071).toStr()  # (0.5, 0.5, 0.7071)
           >>> Nvector(0.5, 0.5, 0.7071, 1).toStr(-3)  # (0.500, 0.500, 0.707, +1.00)
        '''
        t = Vector3d.toStr(self, prec=prec, fmt=NN, sep=sep)
        if self.h:
            t = sep.join((t, self.hStr()))
        return t if not fmt else (fmt % (t,))

    def toVector3d(self, norm=True):
        '''Convert this n-vector to a 3-D vector, I{ignoring
           the height}.

           @kwarg norm: Normalize the 3-D vector (C{bool}).

           @return: The (normalized) vector (L{Vector3d}).
        '''
        u = self.unit()
        v = Vector3d(u.x, u.y, u.z, name=self.name)
        return v.unit() if norm else v

    def to4xyzh(self, h=None):  # PYCHOK no cover
        '''DEPRECATED, use property C{xyzh} or C{xyz.to4Tuple}C{(}B{C{h}}C{)}.
        '''
        return self.xyzh if h in (None, self.h) else \
               self._xnamed(Vector4Tuple(self.x, self.y, self.z, h))

    def unit(self, ll=None):
        '''Normalize this n-vector to unit length.

           @kwarg ll: Optional, original latlon (C{LatLon}).

           @return: Normalized vector (C{Nvector}).
        '''
        if self._united is None:
            u = Vector3d.unit(self, ll=ll)  # .copy()
            self._united = u._united = _xattrs(u, self, '_h')
        return self._united

    @property_RO
    def xyzh(self):
        '''Get this n-vector's components (L{Vector4Tuple}C{(x, y, z, h)})
        '''
        return self._xnamed(self.xyz.to4Tuple(self.h))


NorthPole = NvectorBase(0, 0, +1, name=_NorthPole_)  #: North pole (C{Nvector}).
SouthPole = NvectorBase(0, 0, -1, name=_SouthPole_)  #: South pole (C{Nvector}).


class _N_vector_(NvectorBase):
    '''(INTERNAL) Minimal, low-overhead C{n-vector}.
    '''
    def __init__(self, x, y, z, h=0):
        self._x, self._y, self._z = x, y, z
        if h:
            self._h = h


class LatLonNvectorBase(LatLonBase):
    '''(INTERNAL) Base class for n-vector-based ellipsoidal
       and spherical C{LatLon} classes.
    '''

    def _update(self, updated, *attrs, **kwds):  # PYCHOK _Nv=None
        '''(INTERNAL) Zap cached attributes if updated.

           @see: C{ellipsoidalNvector.LatLon} and C{sphericalNvector.LatLon}
                 for the special case of B{C{_Nv}}.
        '''
        if updated:
            _Nv = _xkwds_pop(kwds, _Nv=None)
            if _Nv is not None:
                if _Nv._fromll is not None:
                    _Nv._fromll = None
                self._Nv = None
            LatLonBase._update(self, updated, *attrs)

    def others(self, other, name=_other_, up=1):
        '''Refine the class comparison.

           @arg other: The other point (C{LatLon}).
           @kwarg name: Optional, other's name (C{str}).
           @kwarg up: Number of call stack frames up (C{int}).

           @return: The B{C{other}} if compatible.

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        try:
            LatLonBase.others(self, other, name=name, up=up + 1)
        except TypeError:
            if not isinstance(other, NvectorBase):
                raise
        return other

    def toNvector(self, Nvector=NvectorBase, **Nvector_kwds):  # PYCHOK signature
        '''Convert this point to C{Nvector} components, I{including
           height}.

           @kwarg Nvector_kwds: Optional, additional B{C{Nvector}} keyword
                                arguments, ignored if B{C{Nvector=None}}.

           @return: An B{C{Nvector}} or a L{Vector4Tuple}C{(x, y, z, h)} if
                    B{C{Nvector}} is C{None}.

           @raise TypeError: Invalid B{C{Nvector}} or B{C{Nvector_kwds}}.
        '''
        return LatLonBase.toNvector(self, Nvector=Nvector, **Nvector_kwds)


def sumOf(nvectors, Vector=None, h=None, **Vector_kwds):
    '''Return the vectorial sum of two or more n-vectors.

       @arg nvectors: Vectors to be added (C{Nvector}[]).
       @kwarg Vector: Optional class for the vectorial sum (C{Nvector})
                      or C{None}.
       @kwarg h: Optional height, overriding the mean height (C{meter}).
       @kwarg Vector_kwds: Optional, additional B{C{Vector}} keyword
                           arguments, ignored if B{C{Vector=None}}.

       @return: Vectorial sum (B{C{Vector}}) or a L{Vector4Tuple}C{(x, y,
                z, h)} if B{C{Vector}} is C{None}.

       @raise VectorError: No B{C{nvectors}}.
    '''
    n, nvectors = len2(nvectors)
    if n < 1:
        raise VectorError(nvectors=n, txt=_Missing)

    if h is None:
        h = fsum(v.h for v in nvectors) / float(n)

    if Vector is None:
        r = _sumOf(nvectors, Vector=Vector3Tuple).to4Tuple(h)
    else:
        r = _sumOf(nvectors, Vector=Vector, h=h, **Vector_kwds)
    return r


__all__ += _ALL_DOCS(LatLonNvectorBase, NvectorBase)  # classes

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
