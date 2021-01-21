
# -*- coding: utf-8 -*-

u'''(INTERNAL) Base classes for elliposiodal, spherical and N-/vectorial
C{Cartesian}s.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**,
see U{https://www.Movable-Type.co.UK/scripts/latlong.html},
U{https://www.Movable-Type.co.UK/scripts/latlong-vectors.html} and
U{https://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html}..

@newfield example: Example, Examples
'''

from pygeodesy.basics import _xinstanceof
from pygeodesy.datums import Datum, Datums, _ellipsoidal_datum
from pygeodesy.errors import _datum_datum, _IsnotError, _ValueError, _xkwds
from pygeodesy.fmath import cbrt, fsum_, hypot_, hypot2
from pygeodesy.interns import EPS0, NN, _COMMASPACE_, _ellipsoidal_, \
                             _1_0, _2_0, _4_0, _6_0
from pygeodesy.interns import _spherical_  # PYCHOK used!
from pygeodesy.lazily import _ALL_DOCS
from pygeodesy.namedTuples import LatLon4Tuple, Vector4Tuple
from pygeodesy.props import Property_RO, property_RO, property_doc_
from pygeodesy.streprs import Fmt
from pygeodesy.units import Height
from pygeodesy.vector3d import Vector3d, _xyzhdn6

from math import sqrt  # hypot

__all__ = ()
__version__ = '21.01.20'


class CartesianBase(Vector3d):
    '''(INTERNAL) Base class for ellipsoidal and spherical C{Cartesian}.
    '''
    _datum   = None   # L{Datum}, to be overriden
    _Ecef    = None   # preferred C{EcefKarney} class
    _height  = 0      # height (L{Height})

    def __init__(self, xyz, y=None, z=None, datum=None, ll=None, name=NN):
        '''New C{Cartesian...}.

           @arg xyz: An L{Ecef9Tuple}, L{Vector3Tuple}, L{Vector4Tuple}
                     or the C{X} coordinate (C{scalar}).
           @arg y: The C{Y} coordinate (C{scalar}) if B{C{xyz}} C{scalar}.
           @arg z: The C{Z} coordinate (C{scalar}) if B{C{xyz}} C{scalar}.
           @kwarg datum: Optional datum (L{Datum}, L{Ellipsoid}, L{Ellipsoid2}
                         or L{a_f2Tuple}).
           @kwarg ll: Optional, original latlon (C{LatLon}).
           @kwarg name: Optional name (C{str}).

           @raise TypeError: Non-scalar B{C{xyz}}, B{C{y}} or B{C{z}}
                             coordinate or B{C{xyz}} not an L{Ecef9Tuple},
                             L{Vector3Tuple} or L{Vector4Tuple}.
        '''
        x, y, z, h, d, n = _xyzhdn6(xyz, y, z, None, datum, ll)
        Vector3d.__init__(self, x, y, z, ll=ll, name=name or n)
        if h:
            self._height = h
        if d:
            self.datum = _ellipsoidal_datum(d, name=name)

    def _applyHelmert(self, transform, inverse=False, **datum):
        '''(INTERNAL) Return a new cartesian by applying a Helmert
           transform to this cartesian.

           @arg transform: Transform to apply (L{Transform}).
           @kwarg inverse: Optionally, apply the inverse
                           Helmert transform (C{bool}).
           @kwarg datum: Optional datum of the returned point,
                         (C{B{datum}=}L{Datum}).

           @return: The transformed point (C{Cartesian}).

           @note: For C{B{inverse}=True} keyword B{C{datum}} must
                  be C{B{datum}=}L{Datums.WGS84}.
        '''
        xyz = transform.transform(self.x, self.y, self.z, inverse)
        return self._xnamed(self.classof(xyz, **datum))

    @property_doc_(''' this cartesian's datum (L{Datum}).''')
    def datum(self):
        '''Get this cartesian's datum (L{Datum}).
        '''
        return self._datum

    @datum.setter  # PYCHOK setter!
    def datum(self, datum):
        '''Set this cartesian's C{datum} I{without conversion}.

           @arg datum: New datum (L{Datum}).

           @raise TypeError: The B{C{datum}} is not a L{Datum}.
        '''
        _xinstanceof(Datum, datum=datum)
        d = self.datum
        if d is not None:
            if d.isEllipsoidal and not datum.isEllipsoidal:
                raise _IsnotError(_ellipsoidal_, datum=datum)
            elif d.isSpherical and not datum.isSpherical:
                raise _IsnotError(_spherical_, datum=datum)
            self._update(datum != d)
        self._datum = datum

    @property_RO
    def Ecef(self):
        '''Get the ECEF I{class} (L{EcefKarney}).
        '''
        if CartesianBase._Ecef is None:
            from pygeodesy.ecef import EcefKarney
            CartesianBase._Ecef = EcefKarney  # default
        return CartesianBase._Ecef

    @Property_RO
    def _ecef9(self):
        '''(INTERNAL) Helper for L{toCartesian} and L{toEcef}.
        '''
        return self._xnamed(self.Ecef(self.datum).reverse(self, M=True))

    @property_doc_(''' the height (C{meter}).''')
    def height(self):
        '''Get the height (C{meter}).
        '''
        return self._height

    @height.setter  # PYCHOK setter!
    def height(self, height):
        '''Set the height.

           @arg height: New height (C{meter}).

           @raise TypeError: Invalid B{C{height}} C{type}.

           @raise ValueError: Invalid B{C{height}}.
        '''
        h = Height(height)
        self._update(h != self.height)
        self._height = h

    @Property_RO
    def isEllipsoidal(self):
        '''Check whether this cartesian is ellipsoidal (C{bool} or C{None} if unknown).
        '''
        return self.datum.isEllipsoidal if self._datum else None

    @Property_RO
    def isSpherical(self):
        '''Check whether this cartesian is spherical (C{bool} or C{None} if unknown).
        '''
        return self.datum.isSpherical if self._datum else None

    @Property_RO
    def latlon(self):
        '''Get this cartesian's (geodetic) lat- and longitude in C{degrees} (L{LatLon2Tuple}C{(lat, lon)}).
        '''
        return self.toEcef().latlon

    @Property_RO
    def latlonheight(self):
        '''Get this cartesian's (geodetic) lat-, longitude in C{degrees} with height (L{LatLon3Tuple}C{(lat, lon, height)}).
        '''
        return self.toEcef().latlonheight

    @Property_RO
    def latlonheightdatum(self):
        '''Get this cartesian's (geodetic) lat-, longitude in C{degrees} with height and datum (L{LatLon4Tuple}C{(lat, lon, height, datum)}).
        '''
        return self.toEcef().latlonheightdatum

    @Property_RO
    def _N_vector(self):
        '''(INTERNAL) Get the (C{nvectorBase._N_vector_}).
        '''
        from pygeodesy.nvectorBase import _N_vector_
        x, y, z, h = self._n_xyzh(self.datum)
        return _N_vector_(x, y, z, h=h, name=self.name)

    def _n_xyzh(self, datum):
        '''(INTERNAL) Get the n-vector components as L{Vector4Tuple}.
        '''
        # <https://www.Movable-Type.co.UK/scripts/geodesy/docs/
        #        latlon-nvector-ellipsoidal.js.html#line309>
        E = datum.ellipsoid
        x, y, z = self.xyz

        # Kenneth Gade eqn 23
        p = hypot2(x, y) * E.a2_
        q = (z**2 * E.e12) * E.a2_
        r = fsum_(p, q, -E.e4) / _6_0
        s = (p * q * E.e4) / (_4_0 * r**3)
        t = cbrt(fsum_(_1_0, s, sqrt(s * (_2_0 + s))))
        if abs(t) < EPS0:
            raise _ValueError(origin=self, txt=Fmt.EPS0(t))

        u = r * fsum_(_1_0, t, _1_0 / t)
        v = sqrt(u**2 + E.e4 * q)
        t = v * _2_0
        if t < EPS0:
            raise _ValueError(origin=self, txt=Fmt.EPS0(t))
        w = E.e2 * fsum_(u, v, -q) / t

        k = sqrt(fsum_(u, v, w**2)) - w
        if abs(k) < EPS0:
            raise _ValueError(origin=self, txt=Fmt.EPS0(k))
        t = k + E.e2
        if abs(t) < EPS0:
            raise _ValueError(origin=self, txt=Fmt.EPS0(t))
        e = k / t
#       d = e * hypot(x, y)

#       tmp = 1 / hypot(d, z) == 1 / hypot(e * hypot(x, y), z)
        t = hypot_(x * e, y * e, z)  # == 1 / tmp
        if t < EPS0:
            raise _ValueError(origin=self, txt=Fmt.EPS0(t))
        h = fsum_(k, E.e2, -_1_0) / k * t
        s = e / t  # == e * tmp
        return Vector4Tuple(x * s, y * s, z / t, h, name=self.name)

    @Property_RO
    def philam(self):
        '''Get this cartesian's (geodetic) lat- and longitude in C{radians} (L{PhiLam2Tuple}C{(phi, lam)}).
        '''
        return self.toEcef().philam

    @Property_RO
    def philamheight(self):
        '''Get this cartesian's (geodetic) lat-, longitude in C{radians} with height (L{PhiLam3Tuple}C{(phi, lam, height)}).
        '''
        return self.toEcef().philamheight

    @Property_RO
    def philamheightdatum(self):
        '''Get this cartesian's (geodetic) lat-, longitude in C{radians} with height and datum (L{PhiLam4Tuple}C{(phi, lam, height, datum)}).
        '''
        return self.toEcef().philamheightdatum

    def to3llh(self, datum=None):  # PYCHOK no cover
        '''DEPRECATED, use property C{latlonheightdatum} or property C{latlonheight}.

           @return: A L{LatLon4Tuple}C{(lat, lon, height, datum)}.

           @note: This method returns a B{C{-4Tuple}} I{and not a} C{-3Tuple}
                  as its name suggests.
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

    def toDatum(self, datum2, datum=None):
        '''Convert this cartesian from one datum to an other.

           @arg datum2: Datum to convert I{to} (L{Datum}).
           @kwarg datum: Datum to convert I{from} (L{Datum}).

           @return: The converted point (C{Cartesian}).

           @raise TypeError: B{C{datum2}} or B{C{datum}}
                             invalid.
        '''
        _xinstanceof(Datum, datum2=datum2)

        if datum not in (None, self.datum):
            c = self.toDatum(datum)
        else:
            c = self

        i, d = False, c.datum
        if d == datum2:
            return c.copy() if c is self else c

        elif d == Datums.WGS84:
            d = datum2  # convert from WGS84 to datum2

        elif datum2 == Datums.WGS84:
            i = True  # convert to WGS84 by inverse transform

        else:  # neither datum2 nor c.datum is WGS84, invert to WGS84 first
            c = c._applyHelmert(d.transform, True, datum=Datums.WGS84)
            d = datum2

        return c._applyHelmert(d.transform, i, datum=datum2)

    convertDatum = toDatum  # for backward compatibility

    def toEcef(self):
        '''Convert this cartesian to geodetic (lat-/longitude) coordinates.

           @return: An L{Ecef9Tuple}C{(x, y, z, lat, lon, height,
                    C, M, datum)} with C{C} and C{M} if available.

           @raise EcefError: A C{.datum} or an ECEF issue.
        '''
        return self._ecef9

    def toLatLon(self, datum=None, LatLon=None, **LatLon_kwds):  # see .ecef.Ecef9Tuple.toDatum
        '''Convert this cartesian to a geodetic (lat-/longitude) point.

           @kwarg datum: Optional datum (L{Datum}, L{Ellipsoid}, L{Ellipsoid2}
                         or L{a_f2Tuple}).
           @kwarg LatLon: Optional class to return the geodetic point
                          (C{LatLon}) or C{None}.
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}}
                               keyword arguments, ignored if
                               C{B{LatLon}=None}.

           @return: The geodetic point (B{C{LatLon}}) or if B{C{LatLon}}
                    is C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon,
                    height, C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{datum}} or B{C{LatLon_kwds}}.
        '''
        d = self.datum
        if datum in (None, d):
            r = self.toEcef()
        else:
            c = self.toDatum(datum)
            d = c.datum
            r = c.Ecef(d).reverse(c, M=True)

        if LatLon is not None:  # class or .classof
            kwds = _xkwds(LatLon_kwds, datum=r.datum, height=r.height, name=self.name)
            r = LatLon(r.lat, r.lon, **kwds)
        _datum_datum(r.datum, d)
        return r

    def toNvector(self, Nvector=None, datum=None, **Nvector_kwds):  # PYCHOK Datums.WGS84
        '''Convert this cartesian to C{n-vector} components.

           @kwarg Nvector: Optional class to return the C{n-vector}
                           components (C{Nvector}) or C{None}.
           @kwarg datum: Optional datum (L{Datum}, L{Ellipsoid}, L{Ellipsoid2}
                         or L{a_f2Tuple}) overriding this cartesian's datum.
           @kwarg Nvector_kwds: Optional, additional B{C{Nvector}} keyword
                                arguments, ignored if B{C{Nvector=None}}.

           @return: The C{unit, n-vector} components (B{C{Nvector}}) or a
                    L{Vector4Tuple}C{(x, y, z, h)} if B{C{Nvector}} is C{None}.

           @raise TypeError: Invalid B{C{datum}}.

           @raise ValueError: The B{C{Cartesian}} at origin.

           @example:

           >>> c = Cartesian(3980581, 97, 4966825)
           >>> n = c.toNvector()  # (x=0.622818, y=0.00002, z=0.782367, h=0.242887)
        '''
        d = _ellipsoidal_datum(datum or self.datum, name=self.name)
        r =  self._N_vector.xyzh if d == self.datum else self._n_xyzh(d)

        if Nvector is not None:
            kwds = _xkwds(Nvector_kwds, datum=d, h=r.h, name=self.name)
            r = Nvector(r.x, r.y, r.z, **kwds)
        return r

    def toStr(self, prec=3, fmt=Fmt.SQUARE, sep=_COMMASPACE_):  # PYCHOK expected
        '''Return the string representation of this cartesian.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg fmt: Optional enclosing backets format (string).
           @kwarg sep: Optional separator to join (string).

           @return: Cartesian represented as "[x, y, z]" (string).
        '''
        return Vector3d.toStr(self, prec=prec, fmt=fmt, sep=sep)

    def toVector(self, Vector=None, **Vector_kwds):
        '''Return this cartesian's components as vector.

           @kwarg Vector: Optional class to return the C{n-vector}
                          components (L{Vector3d}) or C{None}.
           @kwarg Vector_kwds: Optional, additional B{C{Vector}} keyword
                               arguments, ignored if C{B{Vector}=None}.

           @return: A B{C{Vector}} or an L{Vector3Tuple}C{(x, y, z)}
                    if B{C{Vector}} is C{None}.

           @raise TypeError: Invalid B{C{Vector}} or B{C{Vector_kwds}}.
        '''
        return self.xyz if Vector is None else \
              self._xnamed(Vector(self.x, self.y, self.z, **Vector_kwds))


__all__ += _ALL_DOCS(CartesianBase)

#   xyz = Vector3d.xyz
#   '''Get this cartesian's X, Y and Z components (L{Vector3Tuple}C{(x, y, z)}).
#   '''

# **) MIT License
#
# Copyright (C) 2016-2021 -- mrJean1 at Gmail -- All Rights Reserved.
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
