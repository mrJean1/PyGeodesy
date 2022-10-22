
# -*- coding: utf-8 -*-

u'''Class L{JacobiConformal}, Jacobi's conformal projection of a triaxial ellipsoid, trancoded
from I{Charles Karney}'s C++ class U{JacobiConformal<https://GeographicLib.SourceForge.io/C++/
doc/classGeographicLib_1_1JacobiConformal.html#details>} to pure Python, class L{Triaxial} and
classes L{BetaOmega3Tuple},L{JacobiError}, L{Jacobi2Tuple} and L{TriaxialError}.

@see: U{Geodesics on a triaxial ellipsoid<https://WikiPedia.org/wiki/Geodesics_on_an_ellipsoid#
      Geodesics_on_a_triaxial_ellipsoid>} and U{Triaxial coordinate systems and their geometrical
      interpretation<https://www.Topo.Auth.GR/wp-content/uploads/sites/111/2021/12/09_Panou.pdf>}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

# from pygeodesy.basics import isscalar  # from .fsums
from pygeodesy.constants import EPS, EPS0, EPS02, PI2, PI_3, PI4, _0_0, _0_5, \
                               _1_0, _3_0, _4_0, isfinite
from pygeodesy.elliptic import Elliptic, Property_RO
# from pygeodesy.errors import _ValueError  # from .fsums
from pygeodesy.fmath import fdot, Fmt, hypot, hypot_, hypot2, hypot2_, norm2
from pygeodesy.fsums import Fsum, fsum_, isscalar, _ValueError
from pygeodesy.interns import NN, _distant_, _height_, _inside_, _near_, _not_, _null_, _opposite_, \
                             _outside_, _spherical_, _too_, _x_, _y_
# from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS  # from .vector3d
from pygeodesy.named import _NamedBase, _NamedTuple, _Pass
from pygeodesy.namedTuples import LatLon3Tuple, Vector3Tuple, Vector4Tuple
# from pygeodesy.streprs import Fmt  # from .fmath
from pygeodesy.units import Degrees, Height_, Meter, Meter2, Meter3, Radians, Radius, Scalar
from pygeodesy.utily import asin1, atan2d, sincos2, sincos2d, sincos2d_
from pygeodesy.vector3d import _ALL_LAZY, _MODS, _otherV3d, Vector3d

from math import atan2, fabs, sqrt

__all__ = _ALL_LAZY.triaxials
__version__ = '22.10.22'

_not_ordered_ = _not_('ordered')
_TRIPS        =  1074  # Eberly root finder


class BetaOmega2Tuple(_NamedTuple):
    '''2-Tuple C{(beta, omega)} with I{ellipsoidal} lat- and
       longitude C{beta} and C{omega} both in C{Radians} (or
       C{Degrees}).
    '''
    _Names_ = ('beta', 'omega')
    _Units_ = (_Pass,  _Pass)

    def toDegrees(self, **toDMS_kwds):
        '''Convert this L{BetaOmega2Tuple} to C{Degrees} or C{toDMS}.

           @return: L{BetaOmega2Tuple}C{(beta, omega)} with
                    C{beta} and C{omega} both in C{Degrees}
                    or as an L{toDMS} string provided some
                    B{C{toDMS_kwds}} are supplied.
        '''
        return _toDegrees(self, *self, **toDMS_kwds)

    def toRadians(self):
        '''Convert this L{BetaOmega2Tuple} to C{Radians}.

           @return: L{BetaOmega2Tuple}C{(beta, omega)} with
                    C{beta} and C{omega} both in C{Radians}.
        '''
        return _toRadians(self, *self)


class BetaOmega3Tuple(_NamedTuple):
    '''3-Tuple C{(beta, omega, height)} with I{ellipsoidal} lat- and
       longitude C{beta} and C{omega} both in C{Radians} (or C{Degrees})
       and the C{height}, rather the (signed) I{distance} to the triaxial's
       surface (measured along the radial line to the triaxial's center)
       in C{meter}, conventionally.
    '''
    _Names_ = BetaOmega2Tuple._Names_ + (_height_,)
    _Units_ = BetaOmega2Tuple._Units_ + ( Meter,)

    def toDegrees(self, **toDMS_kwds):
        '''Convert this L{BetaOmega3Tuple} to C{Degrees} or C{toDMS}.

           @return: L{BetaOmega3Tuple}C{(beta, omega, height)} with
                    C{beta} and C{omega} both in C{Degrees} or as an
                    L{toDMS} string provided some B{C{toDMS_kwds}}
                    are supplied.
        '''
        return _toDegrees(self, *self, **toDMS_kwds)

    def toRadians(self):
        '''Convert this L{BetaOmega3Tuple} to C{Radians}.

           @return: L{BetaOmega3Tuple}C{(beta, omega, height)} with
                    C{beta} and C{omega} both in C{Radians}.
        '''
        return _toRadians(self, *self)

    def to2Tuple(self):
        '''Reduce this L{BetaOmega3Tuple} to a L{BetaOmega2Tuple}.
        '''
        return BetaOmega2Tuple(*self[:2])


class JacobiError(_ValueError):
    '''Raised for L{JacobiConformal} issues.
    '''
    pass  # ...


class Jacobi2Tuple(_NamedTuple):
    '''2-Tuple C{(x, y)} with a Jacobi Conformal C{x} and C{y}
       projection, both in C{Radians} (or C{Degrees}).
    '''
    _Names_ = (_x_,   _y_)
    _Units_ = (_Pass, _Pass)

    def toDegrees(self, **toDMS_kwds):
        '''Convert this L{Jacobi2Tuple} to C{Degrees} or C{toDMS}.

           @return: L{Jacobi2Tuple}C{(x, y)} with C{x} and C{y}
                    both in C{Degrees} or as an L{toDMS} string
                    provided some B{C{toDMS_kwds}} are supplied.
        '''
        return _toDegrees(self, *self, **toDMS_kwds)

    def toRadians(self):
        '''Convert this L{Jacobi2Tuple} to C{Radians}.

           @return: L{Jacobi2Tuple}C{(x, y)} with C{x}
                    and C{y} both in C{Radians}.
        '''
        return _toRadians(self, *self)


class TriaxialError(_ValueError):
    '''Raised for L{Triaxial} issues.
    '''
    pass  # ...


class Triaxial(_NamedBase):  # _NamedEnumItem
    '''Triaxial ellipsoid with semi-axes C{a}, C{b} and C{c}, oriented such that
       the large principal ellipse C{ab} is the equator I{Z}=0, I{beta}=0, while
       the small principal ellipse C{ac} is the prime meridian, I{Y}=0, I{omega}=0.
       The four umbilic points, C{abs}(I{omega}) = C{abs}(I{beta}) = C{PI/2}, lie
       on the middle principal ellipse C{bc} in plane I{X}=0, I{omega}=C{PI/2}.

       @note: Geodetic lat- and longitudes are in C{degrees}, I{ellipsoidal} lat-
              and longitude C{beta} and C{omega} are in C{Radians} by default
              (or C{Degrees}.
    '''
    _Error = TriaxialError

    def __init__(self, a, b, c, name=NN):
        '''New L{Triaxial}.

           @arg a: The largest semi-axis (C{meter}, conventionally).
           @arg b: The middle semi-axis (C{meter}, same units as B{C{a}}).
           @arg c: The smallest semi-axis (C{meter}, same units as B{C{a}}).
           @kwarg name: Optional name (C{str}).

           @note: The semi-axes must be ordered as C{B{a} >= B{b} >= B{c} > 0}
                  and be ellipsoidal, C{B{a} > B{c}}.

           @raise TriaxialError: Semi-axes not ordered, spherical or invalid.
        '''
        if name:
            self.name = name
        E = self._Error

        self._a = a = Radius(a=a, Error=E)
        self._b = b = Radius(b=b, Error=E)
        self._c = c = Radius(c=c, Error=E)
        if not (isfinite(a) and a >= b >= c > 0):
            raise E(a=a, b=b, c=c, txt=_not_ordered_)

        self._a2b2 = (a - b) * (a + b)  # == a**2 - b**2 == E_sub_e**2
        self._b2c2 = (b - c) * (b + c)  # == b**2 - c**2 == E_sub_y**2
        self._a2c2 = (a - c) * (a + c)  # == a**2 - c**2 == E_sub_x**2
        if not (a > c and self._a2c2 > 0 and self.e2ac > 0):
            raise E(a=a, c=c, e2ac=self.e2ac, txt=_spherical_)

        self._a2_b2 = _1_0 if a == b else (a / b)**2
        self._c2_b2 = _1_0 if c == b else (c / b)**2

    def __str__(self):
        return self.toStr()

    @Property_RO
    def a(self):
        '''Get the largest semi-axis (C{meter}, conventionally).
        '''
        return self._a

    @Property_RO
    def area(self):
        '''Get the surface area (C{meter} I{squared}).

           @see: U{Surface area<https://WikiPedia.org/wiki/Ellipsoid#Surface_area>}.
        '''
        kp2, k2 = self._k2_kp2  # swapped!
        aE = Elliptic(k2, _0_0, kp2, _1_0)
        c2 = self._1e2ac  # cos(phi)**2 == (c / a)**2
        s2 = self.  e2ac  # sin(phi)**2 ==  1 - c2
        s = sqrt(s2)
        r = asin1(s)  # phi == atan2(sqrt(c2), s)
        a = self.c**2 + (aE.fE(r) * s2 + aE.fF(r) * c2) / s * self.a * self.b
        return Meter2(area=PI2 * a)

    def area_p(self, p=1.6075):
        '''I{Approximate} the surface area (C{meter} I{squared}).

           @kwarg p: Exponent (C{scalar}), 1.6 for near-spherical or 1.5849625007
                     for "near-flat" triaxials

           @see: U{Surface area<https://WikiPedia.org/wiki/Ellipsoid#Approximate_formula>}.
        '''
        a, b, c, _p = self.a, self.b, self.c, pow
        a = _p(fsum_(_p(a * b, p), _p(a * c, p), _p(b * c, p)) / _3_0, _1_0 / p)
        return Meter2(area_p=PI4 * a)

    @Property_RO
    def b(self):
        '''Get the middle semi-axis (C{meter}, same units as B{C{a}}).
        '''
        return self._b

    @Property_RO
    def c(self):
        '''Get the smallest semi-axis (C{meter}, same units as B{C{a}}).
        '''
        return self._c

    @Property_RO
    def e2ab(self):
        '''Get the C{ab} ellipse' I{(1st) eccentricity squared} (C{scalar}), M{1 - (b / a)**2}.
        '''
        return Scalar(e2ab=(_1_0 - self._1e2ab) if self.b != self.a else _0_0)

    @Property_RO
    def _1e2ab(self):
        '''(INTERNAL) Get C{1 - e2ab}.
        '''
        return (_1_0 / self._a2_b2) if self.b != self.a else _1_0

    @Property_RO
    def e2bc(self):
        '''Get the C{bc} ellipse' I{(1st) eccentricity squared} (C{scalar}), M{1 - (c / b)**2}.
        '''
        return Scalar(e2bc=(_1_0 - self._1e2bc) if self.c != self.b else _0_0)

    @Property_RO
    def _1e2bc(self):
        '''(INTERNAL) Get C{1 - e2bc}.
        '''
        return self._c2_b2

    @Property_RO
    def e2ac(self):
        '''Get the C{ac} ellipse' I{(1st) eccentricity squared} (C{scalar}), M{1 - (c / a)**2}.
        '''
        return Scalar(e2ac=(_1_0 - self._1e2ac) if self.c != self.a else _0_0)

    @Property_RO
    def _1e2ac(self):
        '''(INTERNAL) Get C{1 - e2ac}.
        '''
        return (self.c / self.a)**2 if self.c != self.a else _1_0

    @Property_RO
    def _Exycr4(self):
        '''(INTERNAL) Helper for C{.forwardBetOmg}.
        '''
        return self._Exyur4(self.c)

    def _Exyur4(self, u):
        '''(INTERNAL) Helper for C{.forwardBetOmg}.
        '''
        if u > 0:
            u2 = u**2
            x  = self._sqrt(_1_0 + self._a2c2 / u2) * u
            y  = self._sqrt(_1_0 + self._b2c2 / u2) * u
        else:
            x  = y = u = _0_0
        return x, y, u, (self._a2b2 / self._a2c2)

    def forwardBetaOmega(self, beta, omega, height=0, name=NN):
        '''Convert I{ellipsoidal} lat- and longitude C{beta}, C{omega}
           and height to cartesian.

           @arg beta: Ellipsoidal latitude (C{radians} or L{Degrees}).
           @arg omega: Ellipsoidal longitude (C{radians} or L{Degrees}).
           @arg height: Height above or below the ellipsoid's surface
                        (C{meter}, same units as this triaxial's C{a},
                        C{b} and C{c} semi-axes).
           @kwarg name: Optional name (C{str}).

           @return: A L{Vector3Tuple}C{(x, y, z)}.

           @see: U{Expressions (23-25)<https://www.Topo.Auth.GR/wp-content/
                 uploads/sites/111/2021/12/09_Panou.pdf>}.
        '''
        sa, ca = _sincos2(beta)
        sb, cb = _sincos2(omega)

        if height:
            h = Height_(height=height, low=-self.c, Error=self._Error)
            x, y, z, r = self._Exyur4(h + self.c)
        else:
            x, y, z, r = self._Exycr4
        if z:  # and x and y
            x *= cb * self._sqrt(ca**2 + r * sa**2)
            y *= ca * sb
            z *= sa * self._sqrt(_1_0  - r * cb**2)
        return Vector3Tuple(x, y, z, name=name)

    def forwardBetaOmega_(self, sbeta, cbeta, somega, comega, name=NN):
        '''Convert I{ellipsoidal} lat- and longitude C{beta} and C{omega} to
           cartesian coordinates I{on the ellipsoid's surface}.

           @arg sbeta: Ellipsoidal latitude C{sin(beta)} (C{scalar}).
           @arg cbeta: Ellipsoidal latitude C{cos(beta)} (C{scalar}).
           @arg somega: Ellipsoidal longitude C{sin(omega)} (C{scalar}).
           @arg comega: Ellipsoidal longitude C{cos(omega)} (C{scalar}).
           @kwarg name: Optional name (C{str}).

           @return: A L{Vector3Tuple}C{(x, y, z)}.

           @see: U{Triaxial ellipsoid coordinate system<https://WikiPedia.org/wiki/
                 Geodesics_on_an_ellipsoid#Triaxial_ellipsoid_coordinate_system>}
                 and U{expressions (23-25)<https://www.Topo.Auth.GR/wp-content/
                 uploads/sites/111/2021/12/09_Panou.pdf>}.
        '''
        sa, ca  =  self._norm2(sbeta,  cbeta)
        sb, cb  =  self._norm2(somega, comega)

        b2_a2   =  self._1e2ab  # == (b / a)**2
        c2_a2   = -self._1e2ac  # == (c / a)**2
        a2c2_a2 =  self.  e2ac  # (a**2 - c**2) / a**2 == 1 - (c / a)**2

        x = Fsum(_1_0, -b2_a2 * sa**2, c2_a2 * ca**2).fover(a2c2_a2)
        x = self.a * cb * self._sqrt(x)
        y = self.b * ca * sb
        z = Fsum(c2_a2,         sb**2, b2_a2 * cb**2).fover(a2c2_a2)
        z = self.c * sa * self._sqrt(z)
        return Vector3Tuple(x, y, z, name=name)

    def forwardLatLon(self, lat, lon, height=0, name=NN):
        '''Convert I{geodetic} lat-, longitude and heigth to cartesian.

           @arg lat: Geodetic latitude (C{degrees}).
           @arg lon: Geodetic longitude (C{degrees}).
           @arg height: Height above the ellipsoid (C{meter}, same units
                        as this triaxial's C{a}, C{b} and C{c} axes).
           @kwarg name: Optional name (C{str}).

           @return: A L{Vector3Tuple}C{(x, y, z)}.

           @see: U{Expressions (9-11)<https://www.Topo.Auth.GR/wp-content/
                 uploads/sites/111/2021/12/09_Panou.pdf>}.
        '''
        return self._forward3(height, name, *sincos2d_(lat, lon))

    def forwardLatLon_(self, slat, clat, slon, clon, height=0, name=NN):
        '''Convert I{geodetic} lat-, longitude and heigth to cartesian.

           @arg slat: Geodetic latitude C{sin(lat)} (C{scalar}).
           @arg clat: Geodetic latitude C{cos(lat)} (C{scalar}).
           @arg slon: Geodetic longitude C{sin(lon)} (C{scalar}).
           @arg clon: Geodetic longitude C{cos(lon)} (C{scalar}).
           @arg height: Height above the ellipsoid (C{meter}, same units
                        as this triaxial's axes C{a}, C{b} and C{c}).
           @kwarg name: Optional name (C{str}).

           @return: A L{Vector3Tuple}C{(x, y, z)}.

           @see: U{Expressions (9-11)<https://www.Topo.Auth.GR/wp-content/
                 uploads/sites/111/2021/12/09_Panou.pdf>}.
        '''
        sa, ca = self._norm2(slat, clat)
        sb, cb = self._norm2(slon, clon)
        return self._forward3(height, name, sa, ca, sb, cb)

    def _forward3(self, h, name, sa, ca, sb, cb):
        '''(INTERNAL) Helper for C{.forwardLatLon} and C{.forwardLatLon_}.
        '''
        ca_sb = ca * sb
        # 1 - (1 - (c/a)**2) * sa**2 - (1 - (b / a)**2) * ca**2 * sb**2
        t = fsum_(_1_0, -self.e2ac * sa**2, -self.e2ab * ca_sb**2)
        N = self.a / self._sqrt(t)  # prime vertical
        x = (h + N)               * ca * cb
        y = (h + N * self._1e2ab) * ca_sb
        z = (h + N * self._1e2ac) * sa
        return Vector3Tuple(x, y, z, name=name)

    def hartzell(self, pov, los=None, name=NN):
        '''Compute the intersection of a Line-Of-Sight from a Point-Of-View in space
           with the surface of this triaxial.

           @arg pov: Point-Of-View outside this triaxial (C{Cartesian}, L{Ecef9Tuple}
                     or L{Vector3d}).
           @kwarg los: Line-Of-Sight, I{direction} to this triaxial (L{Vector3d}) or
                       C{None} to point to this triaxial's center.
           @kwarg name: Optional name (C{str}).

           @return: L{Vector4Tuple}C{(x, y, z, h)} on this triaxial's surface.

           @raise TriaxialError: Null B{C{pov}} or B{C{los}} vector, B{C{pov}} is
                                 inside this triaxial or B{C{los}} points outside
                                 this triaxial or points in an opposite direction.

           @raise TypeError: Invalid B{C{pov}} or B{C{los}}.

           @see: Function L{pygeodesy.tyr3d} for B{C{los}} and Hartzell, S. U{I{Satellite
                 Line-of-Sight Intersection with Earth}<https://StephenHartzell.Medium.com/
                 satellite-line-of-sight-intersection-with-earth-d786b4a6a9b6>}.
        '''
        def _Error(txt):
            return TriaxialError(pov=pov, los=los, triaxial=self, txt=txt)

        v, h = _hartzell3d2(pov, los, self.a, self.b, self.c, _Error)
        return Vector4Tuple(v.x, v.y, v.z, h, name=name or self.hartzell.__name__)

    def height4(self, x_xyz, y=None, z=None, normal=True, eps=EPS):
        '''Compute the projection on and the height of a cartesian above or below
           this triaxial's surface.

           @arg x_xyz: X component (C{scalar}) or a cartesian (C{Cartesian},
                       L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple} or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), needed if B{C{x_xyz}} if C{scalar}.
           @kwarg z: Z component (C{scalar}), needed if B{C{x_xyz}} if C{scalar}.
           @kwarg normal: If C{True} the projection is the nearest point on this
                          triaxial's surface, otherwise the intersection of the
                          radial line to the center and this triaxial's surface.
           @kwarg eps: Tolerance for root finding (C{scalar}).

           @return: L{Vector4Tuple}C{(x, y, z, h)} with the cartesian coordinates
                    C{x}, C{y} and C{z} of the projection on or the intersection
                    with and with the height C{h} above or below the triaxial's
                    surface in C{meter}, conventionally.

           @raise TriaxialError: Non-cartesian B{C{xyz}}, invalid B{C{eps}} or
                                 no convergence in root finding.

           @see: U{Eberly<https://www.GeometricTools.com/Documentation/DistancePointEllipseEllipsoid.pdf>}.
        '''
        v = Vector3d(x_xyz, y, z) if isscalar(x_xyz) else \
           _otherV3d(x_xyz=x_xyz)

        if normal:  # perpendicular to triaxial
            x, y, z, h, i = _normalTo5(v.x, v.y, v.z, self, eps=eps)
        else:  # radial to triaxial's center
            x, y, z = self.forwardBetaOmega_(v.z, hypot(v.x, v.y), v.y, v.x)
            h, i = v.minus_(x, y, z).length, None

        if h and hypot2_(v.x / self.a, v.y / self.b, v.z / self.c) < _1_0:
            h = -h  # below the surface
        return Vector4Tuple(x, y, z, h, iteration=i, name=self.height4.__name__)

    @Property_RO
    def _k2_kp2(self):
        '''(INTERNAL) Get C{k2} and C{kp2} for C{._xE}, C{._yE} and C{.area}.
        '''
        # k2  = a2b2 / a2c2 * c2_b2
        # kp2 = b2c2 / a2c2 * a2_b2
        # b2  = b**2
        # xE  = Elliptic(k2,  -a2b2 / b2, kp2, a2_b2)
        # yE  = Elliptic(kp2, +b2c2 / b2, k2,  c2_b2)
        # aE  = Elliptic(kp2,  0,         k2,  1)
        return (self._a2b2 / self._a2c2 * self._c2_b2,
                self._b2c2 / self._a2c2 * self._a2_b2)

    def _norm2(self, s, c, *a):
        '''(INTERNAL) Normalize C{s} and C{c} iff needed.
        '''
        if fabs(s) > _1_0 or fabs(c) > _1_0:
            s, c = norm2(s, c)
        if a:
            s, c = norm2(s * self.b, c * a[0])
        return (s or _0_0), (c or _0_0)

    def reverseBetaOmega(self, x, y, z, name=NN):
        '''Convert cartesian to I{ellipsoidal} lat- and longitude, C{beta}, C{omega}
           and height.

           @arg x: X coordinate along C{a}-axis (C{meter}, same units as the axes).
           @arg y: Y coordinate along C{b}-axis (C{meter}, same units as the axes).
           @arg z: Z coordinate along C{c}-axis (C{meter}, same units as the axes).
           @kwarg name: Optional name (C{str}).

           @return: A L{BetaOmega3Tuple}C{(beta, omega, height)} with C{beta} and
                    C{omega} in C{Radians} and C{height} in C{meter}, same units
                    as this triaxial's axes.

           @see: U{Expressions (21-22)<https://www.Topo.Auth.GR/wp-content/uploads/
                 sites/111/2021/12/09_Panou.pdf>}.
        '''
        v = Vector3d(x, y, z)
        a, b, h = self._reverse3(x, y, z, atan2, v, self.forwardBetaOmega_)
        return BetaOmega3Tuple(Radians(beta=a), Radians(omega=b), h, name=name)

    def reverseLatLon(self, x, y, z, name=NN):
        '''Convert cartesian to I{geodetic} lat-, longitude and height.

           @arg x: X coordinate along C{a}-axis (C{meter}, same units as the axes).
           @arg y: Y coordinate along C{b}-axis (C{meter}, same units as the axes).
           @arg z: Z coordinate along C{c}-axis (C{meter}, same units as the axes).
           @kwarg name: Optional name (C{str}).

           @return: A L{LatLon3Tuple}C{(lat, lon, height)} with C{lat} and C{lon}
                    in C{degrees} and C{height} in C{meter}, same units as this
                    triaxial's axes.

           @see: U{Expressions (4-5)<https://www.Topo.Auth.GR/wp-content/uploads/
                 sites/111/2021/12/09_Panou.pdf>}.
        '''
        v  = Vector3d(x, y, z)
        x *= self._1e2ac  # == 1 - e_sub_x**2
        y *= self._1e2bc  # == 1 - e_sub_y**2
        t  = self._reverse3(x, y, z, atan2d, v, self.forwardLatLon_)
        return LatLon3Tuple(*t, name=name)

    def _reverse3(self, x, y, z, atan2_, v, forward_):
        '''(INTERNAL) Helper for C{.reverseBetOmg} and C{.reverseLatLon}.
        '''
        d = hypot( x, y)
        a = atan2_(z, d)
        b = atan2_(y, x)
        h = v.minus_(*forward_(z, d, y, x)).length
        return a, b, h

    def _sqrt(self, x):
        '''(INTERNAL) Helper.
        '''
        if x < 0:
            raise self._Error(Fmt.PAREN(sqrt=x))
        return sqrt(x) if x > EPS02 else _0_0

    def toStr(self, prec=9, name=NN, **unused):  # PYCHOK signature
        '''Return this C{Triaxial} as a string.

           @kwarg prec: Precision, number of decimal digits (0..9).
           @kwarg name: Override name (C{str}) or C{None} to exclude
                        this ellipsoid's name.

           @return: This C{Triaxial}'s attributes (C{str}).
        '''
        T = self.__class__
        try:
            t = (T.xyQ2.name,)
        except AttributeError:
            t = ()
        return self._instr(name, prec, T.a.name, T.b.name, T.c.name,
                                       T.e2ab.name, T.e2bc.name, T.e2ac.name,
                                       T.area.name, T.volume.name, *t)

    @Property_RO
    def volume(self):
        '''Get the volume (C{meter**3}), M{4 / 3 * PI * a * b * c}.
        '''
        return Meter3(volume=self.a * self.b * self.c * PI_3 * _4_0)


class JacobiConformal(Triaxial):
    '''This is a conformal projection of the ellipsoid to a plane in which the grid
       lines are straight, see Jacobi, C. G. J. U{I{Vorlesungen über Dynamik}
       <https://Books.Google.com/books?id=ryEOAAAAQAAJ&pg=PA212>}, page 212ff.

       Ellipsoidal coordinates I{beta} and I{omega} are converted to Jacobi Conformal
       I{y} respectively I{x} separately.  Jacobi's coordinates have been multiplied
       by C{sqrt(B{a}**2 - B{c}**2) / (2 * B{b})} so that the customary results are
       returned in the case of an ellipsoid of revolution (or a sphere, I{currently
       not supported}).

       Copyright (C) U{Charles Karney<mailto:Charles@Karney.com>} (2014-2020) and
       licensed under the MIT/X11 License.

       @note: This constructor can not be used to specify a sphere.

       @see: L{Triaxial}, C++ class U{JacobiConformal<https://GeographicLib.SourceForge.io/
             C++/doc/classGeographicLib_1_1JacobiConformal.html#details>} and U{Jacobi's
             conformal projection<https://GeographicLib.SourceForge.io/C++/doc/jacobi.html>}.
    '''
    _Error = JacobiError

#   @Property_RO
#   def ab(self):
#       '''Get relative magnitude C{ab} (C{None} or C{meter}, same units as B{C{a}}).
#       '''
#       return self._ab

#   @Property_RO
#   def bc(self):
#       '''Get relative magnitude C{bc} (C{None} or C{meter}, same units as B{C{a}}).
#       '''
#       return self._bc

    @Property_RO
    def _xE(self):
        '''(INTERNAL) Get the x-elliptic function.
        '''
        k2, kp2 = self._k2_kp2
        # -a2b2 / b2 == (b2 - a2) / b2 == 1 - a2 / b2 == 1 - a2_b2
        return Elliptic(k2, _1_0 - self._a2_b2, kp2, self._a2_b2)

    def xR(self, omega):
        '''Compute a Jacobi Conformal C{x} projection.

           @arg omega: Ellipsoidal longitude (C{radians} or L{Degrees}).

           @return: The C{x} projection (C{Radians}).
        '''
        return self.xR_(*_sincos2(omega))

    def xR_(self, somega, comega):
        '''Compute a Jacobi Conformal C{x} projection.

           @arg somega: Ellipsoidal longitude C{sin(omega)} (C{scalar}).
           @arg comega: Ellipsoidal longitude C{cos(omega)} (C{scalar}).

           @return: The C{x} projection (C{Radians}).
        '''
        s, c = self._norm2(somega, comega, self.a)
        return Radians(x=self._xE.fPi(s, c) * self._a2_b2)

    def xyR2(self, beta, omega, name=NN):
        '''Compute a Jacobi Conformal C{x} and C{y} projection.

           @arg beta: Ellipsoidal latitude (C{radians} or L{Degrees}).
           @arg omega: Ellipsoidal longitude (C{radians} or L{Degrees}).
           @kwarg name: Optional name (C{str}).

           @return: A L{Jacobi2Tuple}C{(x, y)}.
        '''
        return self.xyR2_(*(_sincos2(beta) + _sincos2(omega)),
                          name=name or self.xyR2.__name__)

    def xyR2_(self, sbeta, cbeta, somega, comega, name=NN):
        '''Compute a Jacobi Conformal C{x} and C{y} projection.

           @arg sbeta: Ellipsoidal latitude C{sin(beta)} (C{scalar}).
           @arg cbeta: Ellipsoidal latitude C{cos(beta)} (C{scalar}).
           @arg somega: Ellipsoidal longitude C{sin(omega)} (C{scalar}).
           @arg comega: Ellipsoidal longitude C{cos(omega)} (C{scalar}).
           @kwarg name: Optional name (C{str}).

           @return: A L{Jacobi2Tuple}C{(x, y)}.
        '''
        return Jacobi2Tuple(self.xR_(somega, comega),
                            self.yR_(sbeta,  cbeta),
                            name=name or self.xyR2_.__name__)

    @Property_RO
    def xyQ2(self):
        '''Get the Jacobi Conformal quadrant size (L{Jacobi2Tuple}C{(x, y)}).
        '''
        return Jacobi2Tuple(Radians(x=self._a2_b2 * self._xE.cPi),
                            Radians(y=self._c2_b2 * self._yE.cPi),
                            name=JacobiConformal.xyQ2.name)

    @Property_RO
    def _yE(self):
        '''(INTERNAL) Get the x-elliptic function.
        '''
        kp2, k2 = self._k2_kp2  # swapped!
        # b2c2 / b2 == (b2 - c2) / b2 == 1 - c2 / b2 == e2bc
        return Elliptic(k2, self.e2bc, kp2, self._c2_b2)

    def yR(self, beta):
        '''Compute a Jacobi Conformal C{y} projection.

           @arg beta: Ellipsoidal latitude (C{radians} or L{Degrees}).

           @return: The C{y} projection (C{Radians}).
        '''
        return self.yR_(*_sincos2(beta))

    def yR_(self, sbeta, cbeta):
        '''Compute a Jacobi Conformal C{y} projection.

           @arg sbeta: Ellipsoidal latitude C{sin(beta)} (C{scalar}).
           @arg cbeta: Ellipsoidal latitude C{cos(beta)} (C{scalar}).

           @return: The C{y} projection (C{Radians}).
        '''
        s, c = self._norm2(sbeta, cbeta, self.c)
        return Radians(y=self._yE.fPi(s, c) * self._c2_b2)


def _hartzell3d2(pov, los, a, b, c, Error):
    '''(INTERNAL) Hartzell's "Satellite Lin-of-Sight Intersection ..."
    '''
    a2 = a**2
    if a == b:
        b2 = a2
        p2 = _1_0
    else:
        b2 = b**2
        p2 = b2 / a2
    c2 = c**2
    q2 = c2 / a2

    V3 = _MODS.vector3d._otherV3d
    p3 =  V3(pov=pov)
    u3 =  V3(los=los) if los else p3.negate()
    u3 =  u3.unit()  # unit vector, opposing signs

    x2, y2, z2 = p3.x2y2z2  # p3.times_(p3).xyz
    ux, vy, wz = u3.times_(p3).xyz
    u2, v2, w2 = u3.x2y2z2  # u3.times_(u3).xyz

    t = p2 * c2, c2, b2
    m = fdot(t, u2, v2, w2)  # a2 factored out
    if m < EPS0:  # zero or near-null LOS vector
        raise Error(_near_(_null_))

    r = fsum_(b2 * w2,       c2 * v2,      -v2 * z2,      vy * wz * 2,
             -w2 * y2,       b2 * u2 * q2, -u2 * z2 * p2, ux * wz * 2 * p2,
             -w2 * x2 * p2, -u2 * y2 * q2, -v2 * x2 * q2, ux * vy * 2 * q2, floats=True)
    if r > 0:  # a2 factored out
        r = sqrt(r) * b * c  # == a * a * b * c / a2
    elif r < 0:  # LOS pointing away from or missing the triaxial
        raise Error(_opposite_ if max(ux, vy, wz) > 0 else _outside_)

    n = fdot(t, ux, vy, wz)  # a2 factored out
    d = (n + r) / m  # (n - r) / m for antipode
    if d > 0:  # POV inside or LOS missing the triaxial
        raise Error(_inside_ if min(x2 - a2, y2 - b2, z2 - c2) < EPS else _outside_)
    elif fsum_(x2, y2, z2) < d**2:  # d past triaxial's center
        raise Error(_too_(_distant_))

    v = p3.minus(u3.times(d))  # Vector3d
    h = p3.minus(v).length
    return v, h


def _normalTo4(x, y, a, b, eps=EPS):
    '''(INTERNAL) Nearest point on and distance to a 2-D ellipse.

       @see: Function C{.ellipsoids._normalTo3} and U{Eberly<https://
             www.GeometricTools.com/Documentation/DistancePointEllipseEllipsoid.pdf>}.
    '''
    def _root2d(r, u, v, g, eps):
        # robust root finder
        u *=  r
        t0 =  v - _1_0
        t1 = _0_0 if g < 0 else (hypot(u, v) - _1_0)
        _f =  fabs
        _h =  hypot2
        for i in range(1, _TRIPS):
            t =   (t0 + t1) * _0_5
            e = _f(t0 - t1)
            if e < eps or t in (t0, t1):
                break
            g = _h(u / (t + r), v / (t + _1_0)) - _1_0
            if g > 0:
                t0 = t
            elif g < 0:
                t1 = t
            else:
                break
        else:  # PYCHOK no cover
            t = _root2d.__name__
            raise TriaxialError(Fmt.no_convergence(e, eps), txt=t)
        return t, i

    if not (a >= b > 0 and eps > 0):
        raise TriaxialError(a=a, b=b, eps=eps)

    i = None
    if y:
        if x:
            u = fabs(x / a)
            v = fabs(y / b)
            g = hypot2(u, v) - _1_0
            if g:  # on the ellipse
                r = (a / b)**2
                t, i = _root2d(r, u, v, g, eps)
                a =  x / (t / r + _1_0)
                b =  y / (t     + _1_0)
                d =  hypot(a - x, b - y)
            else:
                a, b, d = x, y, _0_0
        else:  # x == 0
            a, b, d = _0_0, y, fabs(b - y)

    else:  # PYCHOK no cover
        n =  a * x
        d = (a + b) * (a - b)
        if d > fabs(n):
            r  = n / d
            a *= r
            b *= sqrt(_1_0 - r**2)
            d  = hypot(a - x, b)
        else:
            a, b, d = x, _0_0, fabs(a - x)
    return a, b, d, i


def _normalTo5(x, y, z, T, eps=EPS):  # MCCABE 16
    '''(INTERNAL) Nearest point on and distance to a 3-D triaxial.

       @see: U{Eberly<https://www.GeometricTools.com/Documentation/
             DistancePointEllipseEllipsoid.pdf>}.
    '''
    def _root3d(r, s, u, v, w, g, eps):
        # robust root finder
        u *=  r
        v *=  s
        t0 =  w - _1_0
        t1 = _0_0 if g < 0 else (hypot_(u, v, z) - _1_0)
        _f =  fabs
        _h =  hypot2_
        for i in range(1, _TRIPS):
            t =   (t0 + t1) * _0_5
            e = _f(t0 - t1)
            if e < eps or t in (t0, t1):
                break
            g = _h(u / (t + r), v / (t + s), w / (t + _1_0)) - _1_0
            if g > 0:
                t0 = t
            elif g < 0:
                t1 = t
            else:
                break
        else:  # PYCHOK no cover
            t = _root3d.__name__
            raise TriaxialError(Fmt.no_convergence(e, eps), txt=t)
        return t, i

    a, b, c = T.a, T.b, T.c
    if not (a >= b >= c > 0 and eps > 0):
        raise TriaxialError(a=a, b=b, c=c, eps=eps)

    i = None
    if z:
        if y:
            if x:
                u = fabs(x / a)
                v = fabs(y / b)
                w = fabs(z / c)
                g = hypot2_(u, v, w) - _1_0
                if g:
                    r = (a / c)**2
                    s = (b / c)**2
                    t, i = _root3d(r, s, u, v, w, g, eps)
                    a =  x / (t / r + _1_0)
                    b =  y / (t / s + _1_0)
                    c =  z / (t     + _1_0)
                    d =  hypot_(a - x, b - y, c - z)
                else:  # on the ellipsoid
                    a, b, c, d = x, y, z, _0_0
            else:  # x == 0
                a          = _0_0
                b, c, d, i = _normalTo4(y, z, b, c, eps=eps)
        elif x:  # y == 0
            b          = _0_0
            a, c, d, i = _normalTo4(x, z, a, c, eps=eps)
        else:  # x == y == 0
            a, b, c, d = x, y, z, fabs(c - z)

    else:  # z == 0
        t =  False
        n =  a * x
        d = (a + c) * (a - c)
        if d > fabs(n):
            u =  n / d
            n =  b * y
            d = (b + c) * (b - c)
            if d > fabs(n):  # PYCHOK no cover
                v =  n / d
                d = _1_0 - hypot2(u, v)
                if d > 0:
                    a *= u
                    b *= v
                    c *= sqrt(d)
                    d  = hypot_(a - x, b - y, c)
                    t  = True
        if not t:
            c          = _0_0
            a, b, d, i = _normalTo4(x, y, a, b, eps=eps)

    e = hypot2_(a / T.a, b / T.b, c / T.c) - _1_0
    if fabs(e) > eps:
        raise TriaxialError(x=a, y=b, z=c, d=d, triaxial=T,
                            e=e, eps=eps, txt=_not_('on'))
    return a, b, c, d, i


def _sincos2(x):
    '''Get C{sin} and C{cos} of C{x} in C{radians} or {degrees}.
    '''
    return sincos2d(x) if isinstance(x, Degrees) else (
           sincos2(x)  if isinstance(x, Radians) else
           sincos2(float(x)))  # assume C{radians}


def _toDegrees(inst, a, b, *c, **toDMS_kwds):
    '''Helper for L{BetaOmega3Tuple} and L{Jacobi2Tuple}
    '''
    if toDMS_kwds:
        toDMS = _MODS.dms.toDMS
        a = toDMS(a.toDegrees(), **toDMS_kwds)
        b = toDMS(b.toDegrees(), **toDMS_kwds)
    elif isinstance(a, Degrees) and \
         isinstance(b, Degrees):
        return inst
    else:
        a, b = a.toDegrees(), b.toDegrees()
    return inst.classof(a, b, *c, name=inst.name)


def _toRadians(inst, a, b, *c):
    '''Helper for L{BetaOmega3Tuple} and L{Jacobi2Tuple}
    '''
    return inst if isinstance(a, Radians) and \
                   isinstance(b, Radians) else \
           inst.classof(a.toRadians(), b.toRadians(),
                       *c, name=inst.name)


if __name__ == '__main__':

    from pygeodesy.lazily import printf

    def _km2m(*abc):
        for m in abc:
            yield m * 1e3

    # <https://ArxIV.org/pdf/1909.06452.pdf> Table 1 Semi-axes in km    # Planet
    # <https://www.JPS.NASA.gov/education/images/pdf/ss-moons.pdf>
    # <https://link.Springer.com/article/10.1007/s00190-022-01650-9>
    for n, a, b, c in (('Amalthea',  125.0,        73.0,      64),      # Jupiter
                       ('Ariel',     581.1,       577.9,     577.7),    # Uranus
                       ('Earth',    6378.173435, 6378.1039, 6356.7544),
                       ('Enceladus', 256.6,       251.4,     248.3),    # Saturn
                       ('Europa',   1564.13,     1561.23,   1560.93),   # Jupiter
                       ('Io',       1829.4,      1819.3,    1815.7),    # Jupiter
                       ('Mars',     3394.6,      3393.3,    3376.3),
                       ('Mimas',     207.4,       196.8,     190.6),    # Saturn
                       ('Miranda',   240.4,       234.2,     232.9),    # Uranus
                       ('Moon',     1735.55,     1735.324,  1734.898),  # Earth
                       ('Tethys',    535.6,       528.2,     525.8)):   # Saturn
        t = Triaxial(*_km2m(a, b, c), name=n)
        printf('%r, area_p=%g', t, t.area_p())

# **) MIT License
#
# Copyright (C) 2022-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
