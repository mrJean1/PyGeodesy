
# -*- coding: utf-8 -*-

u'''Triaxal ellipsoid classes L{JacobiConformal}, Jacobi's conformal projection, trancoded
from I{Charles Karney}'s C++ class U{JacobiConformal<https://GeographicLib.SourceForge.io/C++/
doc/classGeographicLib_1_1JacobiConformal.html#details>} to pure Python, I{ordered} L{Triaxial}
and I{unordered} L{Triaxial_} and miscellaneous classes L{BetaOmega2Tuple}, L{BetaOmega3Tuple},
L{Jacobi2Tuple} and L{TriaxialError}.

@see: U{Geodesics on a triaxial ellipsoid<https://WikiPedia.org/wiki/Geodesics_on_an_ellipsoid#
      Geodesics_on_a_triaxial_ellipsoid>} and U{Triaxial coordinate systems and their geometrical
      interpretation<https://www.Topo.Auth.GR/wp-content/uploads/sites/111/2021/12/09_Panou.pdf>}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

# from pygeodesy.basics import isscalar, map1, _zip  # from .fsums, .namedTuples, .streprs
from pygeodesy.constants import EPS, EPS0, EPS02, EPS4, _EPS2e4, INT0, PI2, PI_3, PI4, \
                               _0_0, _0_5, _1_0, _N_1_0, isfinite, isnear1, \
                               _4_0  # PYCHOK used!
from pygeodesy.datums import Datum, Ellipsoid, _spherical_datum, _WGS84
# from pygeodesy.ellipsoids import Ellipsoid  # from .datums
# from pygeodesy.elliptic import Elliptic  # ._MODS
# from pygeodesy.errors import _ValueError  # from .streprs
from pygeodesy.fmath import Fdot, fdot, fmean_, hypot, hypot_, hypot2, hypot2_, norm2
from pygeodesy.fsums import Fsum, fsum, fsum_, isscalar, Property_RO
from pygeodesy.interns import NN, _a_, _b_, _c_, _distant_, _height_, _inside_, \
                             _invalid_, _near_, _not_, _NOTEQUAL_, _null_, _opposite_, \
                             _outside_, _SPACE_, _spherical_, _too_, _x_, _y_
# from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS  # from .vector3d
from pygeodesy.named import _NamedBase, _NamedTuple, _Pass
from pygeodesy.namedTuples import LatLon3Tuple, map1, Vector3Tuple, Vector4Tuple
# from pygeodesy.props import Property_RO  # from .fsums
from pygeodesy.streprs import Fmt, _ValueError, _zip
from pygeodesy.units import Degrees, Float, Height_, Meter, Meter2, Meter3, Radians, Radius
from pygeodesy.utily import asin1, atan2d, km2m, sincos2, sincos2d, sincos2d_
from pygeodesy.vector3d import _ALL_LAZY, _MODS, _otherV3d, Vector3d

from math import atan2, fabs, sqrt

__all__ = _ALL_LAZY.triaxials
__version__ = '22.11.02'

_not_ordered_ = _not_('ordered')
_TRIPS        =  256  # max 55, Eberly 1074?


class _ToNamedBase(_NamedBase):
    '''(INTERNAL) C{-.toDegrees}, C{-.toRadians} base.
    '''
    def _toDegrees(self, a, b, *c, **toDMS_kwds):
        if toDMS_kwds:
            toDMS = _MODS.dms.toDMS
            a = toDMS(a.toDegrees(), **toDMS_kwds)
            b = toDMS(b.toDegrees(), **toDMS_kwds)
        elif isinstance(a, Degrees) and \
             isinstance(b, Degrees):
            return self
        else:
            a, b = a.toDegrees(), b.toDegrees()
        return self.classof(a, b, *c, name=self.name)

    def _toRadians(self, a, b, *c):
        return self if isinstance(a, Radians) and \
                       isinstance(b, Radians) else \
               self.classof(a.toRadians(), b.toRadians(),
                           *c, name=self.name)


class BetaOmega2Tuple(_NamedTuple, _ToNamedBase):
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
        return _ToNamedBase._toDegrees(self, *self, **toDMS_kwds)

    def toRadians(self):
        '''Convert this L{BetaOmega2Tuple} to C{Radians}.

           @return: L{BetaOmega2Tuple}C{(beta, omega)} with
                    C{beta} and C{omega} both in C{Radians}.
        '''
        return _ToNamedBase._toRadians(self, *self)


class BetaOmega3Tuple(_NamedTuple, _ToNamedBase):
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
        return _ToNamedBase._toDegrees(self, *self, **toDMS_kwds)

    def toRadians(self):
        '''Convert this L{BetaOmega3Tuple} to C{Radians}.

           @return: L{BetaOmega3Tuple}C{(beta, omega, height)} with
                    C{beta} and C{omega} both in C{Radians}.
        '''
        return _ToNamedBase._toRadians(self, *self)

    def to2Tuple(self):
        '''Reduce this L{BetaOmega3Tuple} to a L{BetaOmega2Tuple}.
        '''
        return BetaOmega2Tuple(*self[:2])


class Jacobi2Tuple(_NamedTuple, _ToNamedBase):
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
        return _ToNamedBase._toDegrees(self, *self, **toDMS_kwds)

    def toRadians(self):
        '''Convert this L{Jacobi2Tuple} to C{Radians}.

           @return: L{Jacobi2Tuple}C{(x, y)} with C{x}
                    and C{y} both in C{Radians}.
        '''
        return _ToNamedBase._toRadians(self, *self)


class Triaxial_(_NamedBase):  # _NamedEnumItem
    '''I{Unordered} triaxial ellipsoid and base class.

       Triaxial ellipsoids with right-handed semi-axes C{a}, C{b} and C{c}, oriented
       such that the large principal ellipse C{ab} is the equator I{Z}=0, I{beta}=0,
       while the small principal ellipse C{ac} is the prime meridian, plane I{Y}=0,
       I{omega}=0.

       The four umbilic points, C{abs}(I{omega}) = C{abs}(I{beta}) = C{PI/2}, lie on
       the middle principal ellipse C{bc} in plane I{X}=0, I{omega}=C{PI/2}.

       @note: I{Geodetic} C{lat}- and C{lon}gitudes are in C{degrees}, I{geodetic}
              C{phi} and C{lam}bda are in C{radians}, but I{ellipsoidal} lat- and
              longitude C{beta} and C{omega} are in C{Radians} by default (or in
              C{Degrees} if converted).
    '''
    _unordered = True

    def __init__(self, a_triax, b=None, c=None, name=NN):
        '''New I{unordered} L{Triaxial_}.

           @arg a_triax: C{X} semi-axis (C{scalar}, conventionally in C{meter})
                         or an other L{Triaxial} or L{Triaxial_} instance.
           @kwarg b: C{Y} semi-axis (C{meter}, same units as B{C{a}}), required
                     if C{B{a_triax} is scalar}, ignored otherwise.
           @kwarg c: C{Z} semi-axis (C{meter}, same units as B{C{a}}), required
                     if C{B{a_triax} is scalar}, ignored otherwise.
           @kwarg name: Optional name (C{str}).

           @raise TriaxialError: Invalid semi-axis or -axes.
        '''
        try:
            a = a_triax
            t = a._abc3 if isinstance(a, Triaxial_) else map1(float, a, b, c)
        except (TypeError, ValueError) as x:
            raise TriaxialError(a=a, b=b, c=c, cause=x)
        if name:
            self.name = name

        a, b, c = self._abc3 = t
        if self._unordered:  # == not isinstance(self, Triaxial)
            s, _, t = sorted(t)
            if not (isfinite(t) and s > 0):
                raise TriaxialError(a=a, b=b, c=c)  # txt=_invalid_
        elif not (isfinite(a) and a >= b >= c > 0):
            raise TriaxialError(a=a, b=b, c=c, txt=_not_ordered_)
        elif not (a > c and self._a2c2 > 0 and self.e2ac > 0):
            raise TriaxialError(a=a, c=c, e2ac=self.e2ac, txt=_spherical_)

    def __str__(self):
        return self.toStr()

    @Property_RO
    def a(self):
        '''Get the (largest) C{x} semi-axis (C{meter}, conventionally).
        '''
        a, _, _ = self._abc3
        return Radius(a=a)

    @Property_RO
    def _a2b2(self):
        '''(INTERNAL) Get C{a**2 - b**2} == E_sub_e**2.
        '''
        a, b, _ = self._abc3
        return ((a - b) * (a + b)) if a != b else _0_0

    @Property_RO
    def _a2_b2(self):
        '''(INTERNAL) Get C{(a/b)**2}.
        '''
        a, b, _ = self._abc3
        return (a / b)**2 if a != b else _1_0

    @Property_RO
    def _a2c2(self):
        '''(INTERNAL) Get C{a**2 - c**2} == E_sub_x**2.
        '''
        a, _, c = self._abc3
        return ((a - c) * (a + c)) if a != c else _0_0

    @Property_RO
    def area(self):
        '''Get the surface area (C{meter} I{squared}).
        '''
        c, b, a = sorted(self._abc3)
        if a > c:
            a = Triaxial(a, b, c).area if a > b else \
                Ellipsoid(a, b=c).areax  # a == b
        else:  # a == c == b
            a = Meter2(area=a**2 * PI4)
        return a

    def area_p(self, p=1.6075):
        '''I{Approximate} the surface area (C{meter} I{squared}).

           @kwarg p: Exponent (C{scalar} > 0), 1.6 for near-spherical or 1.5849625007
                     for "near-flat" triaxials.

           @see: U{Surface area<https://WikiPedia.org/wiki/Ellipsoid#Approximate_formula>}.
        '''
        a, b, c = self._abc3
        if a == b == c:
            a *= a
        else:
            _p = pow
            a = _p(fmean_(_p(a * b, p), _p(a * c, p), _p(b * c, p)), _1_0 / p)
        return Meter2(area_p=a * PI4)

    @Property_RO
    def b(self):
        '''Get the (middle) C{y} semi-axis (C{meter}, same units as B{C{a}}).
        '''
        _, b, _ = self._abc3
        return Radius(b=b)

    @Property_RO
    def _b2c2(self):
        '''(INTERNAL) Get C{b**2 - c**2} == E_sub_y**2.
        '''
        _, b, c = self._abc3
        return ((b - c) * (b + c)) if b != c else _0_0

    @Property_RO
    def c(self):
        '''Get the (smallest) C{z} semi-axis (C{meter}, same units as B{C{a}}).
        '''
        _, _, c = self._abc3
        return Radius(c=c)

    @Property_RO
    def _c2_b2(self):
        '''(INTERNAL) Get C{(c/b)**2}.
        '''
        _, b, c = self._abc3
        return (c / b)**2 if b != c else _1_0

    @Property_RO
    def e2ab(self):
        '''Get the C{ab} ellipse' I{(1st) eccentricity squared} (C{scalar}), M{1 - (b/a)**2}.
        '''
        return Float(e2ab=(_1_0 - self._1e2ab) or _0_0)

    @Property_RO
    def _1e2ab(self):
        '''(INTERNAL) Get C{1 - e2ab} == C{(b/a)**2}.
        '''
        a, b, _ = self._abc3
        return (b / a)**2 if a != b else _1_0

    @Property_RO
    def e2bc(self):
        '''Get the C{bc} ellipse' I{(1st) eccentricity squared} (C{scalar}), M{1 - (c/b)**2}.
        '''
        return Float(e2bc=(_1_0 - self._1e2bc) or _0_0)

    _1e2bc = _c2_b2  # C{1 - e2bc} == C{(c/b)**2}

    @Property_RO
    def e2ac(self):
        '''Get the C{ac} ellipse' I{(1st) eccentricity squared} (C{scalar}), M{1 - (c/a)**2}.
        '''
        return Float(e2ac=(_1_0 - self._1e2ac) or _0_0)

    @Property_RO
    def _1e2ac(self):
        '''(INTERNAL) Get C{1 - e2ac} == C{(c/a)**2}.
        '''
        a, _, c = self._abc3
        return (c / a)**2 if a != c else _1_0

    @Property_RO
    def _Elliptic(self):
        '''(INTERNAL) Get class L{Elliptic} once.
        '''
        return _MODS.elliptic.Elliptic

    def hartzell4(self, pov, los=None, name=NN):
        '''Compute the intersection of this triaxial's surface with a Line-Of-Sight
           from a Point-Of-View in space.

           @see: Function L{pygeodesy.hartzell4} for further details.
        '''
        return hartzell4(pov, los=los, tri_biax=self, name=name)

    def height4(self, x_xyz, y=None, z=None, normal=True, eps=EPS):
        '''Compute the projection on and the height of a cartesian above or below
           this triaxial's surface.

           @arg x_xyz: X component (C{scalar}) or a cartesian (C{Cartesian},
                       L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple} or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} if C{scalar}.
           @kwarg z: Z component (C{scalar}), required if B{C{x_xyz}} if C{scalar}.
           @kwarg normal: If C{True} the projection is perpendicular to (the nearest
                          point on) this triaxial's surface, otherwise the C{radial}
                          line to this triaxial's center (C{bool}).
           @kwarg eps: Tolerance for root finding and validation (C{scalar}), use a
                       negative value to skip validation.

           @return: L{Vector4Tuple}C{(x, y, z, h)} with the cartesian coordinates
                    C{x}, C{y} and C{z} of the projection on or the intersection
                    with and with the height C{h} above or below the triaxial's
                    surface in C{meter}, conventionally.

           @raise TriaxialError: Non-cartesian B{C{xyz}}, invalid B{C{eps}}, no
                                 convergence in root finding or validation failed.

           @see: Method L{Ellipsoid.height4} and I{Eberly}'s U{Distance from a Point
                 to ... an Ellipsoid ...<https://www.GeometricTools.com/Documentation/
                 DistancePointEllipseEllipsoid.pdf>}.
        '''
        v = Vector3d(x_xyz, y, z) if isscalar(x_xyz) else _otherV3d(x_xyz=x_xyz)

        i, h = None, v.length
        if h < EPS0:  # EPS
            x = y = z = _0_0
            h -= min(self._abc3)
        elif self.isSpherical:
            r, _, _ = self._abc3
            x, y, z = v.times(r / h).xyz
            h -= r
        elif normal:  # perpendicular to triaxial
            a, b, c = self._abc3
            x, y, z = v.xyz
            x, y, z, h, i = _normalTo5( x, y, z, self,  eps=eps) if self.isOrdered else \
                            _normalTo5u(x, y, z, a,b,c, eps=eps)
        else:  # radially to triaxial's center
            x, y, z = v.xyz
            x, y, z = self._radialTo3(z, hypot(x, y), y, x)
            h, i = v.minus_(x, y, z).length, None

        if h > 0 and self.sideOf(v, eps=EPS0) < 0:
            h = -h  # below the surface
        return Vector4Tuple(x, y, z, h, iteration=i, name=self.height4.__name__)

    @Property_RO
    def isOrdered(self):
        '''Is this triaxial I{ordered} and not I{spherical} (C{bool})?
        '''
        a, b, c = self._abc3
        return bool(a >= b > c)  # b > c!

    @Property_RO
    def isSpherical(self):
        '''Is this triaxial I{spherical} (C{bool})?
        '''
        a, b, c = self._abc3
        return bool(a == b == c)

    def normal3d(self, x, y, z, length=_1_0):
        '''Get a 3-D vector at a cartesian on, perpendicular to this triaxial's surface.

           @arg x: X coordinate along C{a}-axis (C{meter}, same units as the axes).
           @arg y: Y coordinate along C{b}-axis (C{meter}, same units as B{C{x}}).
           @arg z: Z coordinate along C{c}-axis (C{meter}, same units as B{C{x}}).
           @kwarg length: Optional length and in-/outward direction (C{scalar}).

           @return: A C{Vector3d(x_, y_, z_)} normalized to B{C{length}}, pointing
                    in- or outward for neg- respectively positive B{C{length}}.

           @note: Cartesian location C{(B{x}, B{y}, B{z})} must be on this triaxial's
                  surface, use method L{Triaxial.sideOf} to validate.
        '''
        # n = 2 * (x / a2, y / b2, z / c2)
        #  == 2 * (x, y * a2 / b2, z * a2 / c2) / a2  # iff ordered
        #  == 2 * (x, y / _1e2ab, z / _1e2ac) / a2
        #  == unit(x, y / _1e2ab, z / _1e2ac).times(length)
        n = self._normal3d.times_(x, y, z)
        if n.length < EPS0:
            raise TriaxialError(x=x, y=y, z=z, txt=_null_)
        return n.times(length / n.length)

    @Property_RO
    def _normal3d(self):
        '''(INTERNAL) Get M{Vector3d((d/a)**2, (d/b)**2, (d/c)**2)}, M{d = max(a, b, c)}.
        '''
        d = max(self._abc3)
        t = tuple(((d / x)**2 if x != d else _1_0) for x in self._abc3)
        return Vector3d(*t)

    def _norm2(self, s, c, *a):
        '''(INTERNAL) Normalize C{s} and C{c} iff not already.
        '''
        if fabs(s) > _1_0 or fabs(c) > _1_0:  # or \
            # fabs(fsum(s**2, c**2, _N_1_0)) > EPS:
            s, c = norm2(s, c)
        if a:
            s, c = norm2(s * self.b, c * a[0])
        return (s or _0_0), (c or _0_0)

    def _radialTo3(self, sbeta, cbeta, somega, comega):
        '''(INTERNAL) I{Unordered} helper for C{.height4}.
        '''
        def _rphi(a, b, sphi, cphi):
            # <https://WikiPedia.org/wiki/Ellipse#Polar_form_relative_to_focus>
            # polar form: radius(phi) = a * b / hypot(a * sphi, b * cphi)
            return (b / hypot(sphi, b / a * cphi)) if a > b else (
                   (a / hypot(cphi, a / b * sphi)) if a < b else a)

        sa, ca = self._norm2(sbeta,  cbeta)
        sb, cb = self._norm2(somega, comega)

        a, b, c = self._abc3
        if a != b:
            a = _rphi(a, b, sb, cb)
        if a != c:
            c = _rphi(a, c, sa, ca)
        z, r = c * sa, c * ca
        x, y = r * cb, r * sb
        return x, y, z

    def sideOf(self, x_xyz, y=None, z=None, eps=EPS4):
        '''Is a cartesian above, below or on the surface of this triaxial?

           @arg x_xyz: X component (C{scalar}) or a cartesian (C{Cartesian},
                       L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple} or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} if C{scalar}.
           @kwarg z: Z component (C{scalar}), required if B{C{x_xyz}} if C{scalar}.
           @kwarg eps: Near surface tolerance(C{scalar}).

           @return: C{INT0} if C{(B{x}, B{y}, B{z})} is near this triaxial's surface
                    within tolerance B{C{eps}}, otherwise a neg- or positive C{float}
                    if in- respectively outside this triaxial.

           @see: Methods L{Triaxial.height4} and L{Triaxial.normal3d}.
        '''
        def _x2_a2_1(xyz, abc):
            for x, a in _zip(xyz, abc):  # strict=True
                yield (x / a)**2
            yield _N_1_0

        t = (x_xyz, y, z) if isscalar(x_xyz) else _otherV3d(x_xyz=x_xyz).xyz
        s = fsum(_x2_a2_1(t, self._abc3), floats=True)
        return s if s > eps or s < -eps else INT0

    def _sqrt(self, x):
        '''(INTERNAL) Helper.
        '''
        if x < 0:
            raise TriaxialError(Fmt.PAREN(sqrt=x))
        return _0_0 if x < EPS02 else sqrt(x)

    def toEllipsoid(self, name=NN):
        '''Convert this triaxial to an L{Ellipsoid}, provided C{a == b} or C{b == c}.

           @return: An L{Ellipsoid} with north along this C{Z} axis if C{a == b},
                    this C{Y} axis if C{a == c} or this C{X} axis if C{b == c}.

           @raise TriaxialError: This C{a} != C{b}, C{b} != C{c} and C{c} != C{a}.

           @see: Method L{Ellipsoid.toTriaxial}.
        '''
        a, b, c = self._abc3
        if a == b:  # N = Z
            b = c
        elif b == c:  # N = X
            a, b = b, a
        elif a != c:
            t = _SPACE_(_a_, _NOTEQUAL_, _b_, _NOTEQUAL_, _c_)
            raise TriaxialError(a=a, b=b, c=c, txt=t)
        return Ellipsoid(a, b=b, name=name or self.name)

    def toStr(self, prec=9, name=NN, **unused):  # PYCHOK signature
        '''Return this C{Triaxial} as a string.

           @kwarg prec: Precision, number of decimal digits (0..9).
           @kwarg name: Override name (C{str}) or C{None} to exclude
                        this triaxial's name.

           @return: This C{Triaxial}'s attributes (C{str}).
        '''
        T = Triaxial_
        t = T.a, T.b, T.c, T.e2ab, T.e2bc, T.e2ac
        if isinstance(self, JacobiConformal):
            t += JacobiConformal.xyQ2,
        t += T.volume, T.area
        return self._instr(name, prec, props=t, area_p=self.area_p())

    @Property_RO
    def volume(self):
        '''Get the volume (C{meter**3}), M{4 / 3 * PI * a * b * c}.
        '''
        a, b, c = self._abc3
        return Meter3(volume=a * b * c * PI_3 * _4_0)


class Triaxial(Triaxial_):
    '''I{Ordered} triaxial ellipsoid.

       @see: L{Triaxial_} for more information.
    '''
    _unordered = False

    def __init__(self, a_triax, b=None, c=None, name=NN):
        '''New I{ordered} L{Triaxial}.

           @arg a_triax: Largest semi-axis (C{scalar}, conventionally in C{meter})
                         or an other L{Triaxial} or L{Triaxial_} instance.
           @kwarg b: Middle semi-axis (C{meter}, same units as B{C{a}}), required
                     if C{B{a_triax} is scalar}, ignored otherwise.
           @kwarg c: Smallest semi-axis (C{meter}, same units as B{C{a}}), required
                     if C{B{a_triax} is scalar}, ignored otherwise.
           @kwarg name: Optional name (C{str}).

           @note: The semi-axes must be ordered as C{B{a} >= B{b} >= B{c} > 0} and
                  must be ellipsoidal, C{B{a} > B{c}}.

           @raise TriaxialError: Semi-axes not ordered, spherical or invalid.
        '''
        Triaxial_.__init__(self, a_triax, b=b, c=c, name=name)

    @Property_RO
    def _a2b2_a2c2(self):
        ''' @see: Method C{forwardBetaOmega}.
        '''
        return self._a2b2 / self._a2c2

    @Property_RO
    def area(self):
        '''Get the surface area (C{meter} I{squared}).

           @see: U{Surface area<https://WikiPedia.org/wiki/Ellipsoid#Surface_area>}.
        '''
        a, b, c = self._abc3
        if a != b:
            kp2, k2 = self._k2_kp2  # swapped!
            aE = self._Elliptic(k2, _0_0, kp2, _1_0)
            c2 = self._1e2ac  # cos(phi)**2 == (c/a)**2
            s2 = self.  e2ac  # sin(phi)**2 == 1 - c2
            s  = sqrt(s2)
            r  = asin1(s)  # phi == atan2(sqrt(c2), s)
            b *= fsum_(aE.fE(r) * s, aE.fF(r) * c2 / s, c / a * c / b, floats=True)
            a  = Meter2(area=a * b * PI2)
        else:  # a == b > c
            a  = Ellipsoid(a, b=c).areax
        return a

    def _exyz3(self, u):
        '''(INTERNAL) Helper for C{.forwardBetOmg}.
        '''
        if u > 0:
            u2 = u**2
            x  = self._sqrt(_1_0 + self._a2c2 / u2) * u
            y  = self._sqrt(_1_0 + self._b2c2 / u2) * u
        else:
            x  = y = u = _0_0
        return x, y, u

    def forwardBetaOmega(self, beta, omega, height=0, name=NN):
        '''Convert I{ellipsoidal} lat- and longitude C{beta}, C{omega}
           and height to cartesian.

           @arg beta: Ellipsoidal latitude (C{radians} or L{Degrees}).
           @arg omega: Ellipsoidal longitude (C{radians} or L{Degrees}).
           @arg height: Height above or below the ellipsoid's surface (C{meter}, same
                        units as this triaxial's C{a}, C{b} and C{c} semi-axes).
           @kwarg name: Optional name (C{str}).

           @return: A L{Vector3Tuple}C{(x, y, z)}.

           @see: Method L{Triaxial.reverseBetaOmega} and U{Expressions (23-25)<https://
                 www.Topo.Auth.GR/wp-content/uploads/sites/111/2021/12/09_Panou.pdf>}.
        '''
        if height:
            h = Height_(height=height, low=-self.c, Error=TriaxialError)
            x, y, z = self._exyz3(h + self.c)
        else:
            x, y, z = self._abc3  # == self._exyz3(self.c)
        if z:  # and x and y:
            sa, ca = _SinCos2(beta)
            sb, cb = _SinCos2(omega)

            r  =      self._a2b2_a2c2
            x *= cb * self._sqrt(ca**2 + r * sa**2)
            y *= ca * sb
            z *= sa * self._sqrt(_1_0  - r * cb**2)
        return Vector3Tuple(x, y, z, name=name)

    def forwardBetaOmega_(self, sbeta, cbeta, somega, comega, name=NN):
        '''Convert I{ellipsoidal} lat- and longitude C{beta} and C{omega}
           to cartesian coordinates I{on the triaxial's surface}.

           @arg sbeta: Ellipsoidal latitude C{sin(beta)} (C{scalar}).
           @arg cbeta: Ellipsoidal latitude C{cos(beta)} (C{scalar}).
           @arg somega: Ellipsoidal longitude C{sin(omega)} (C{scalar}).
           @arg comega: Ellipsoidal longitude C{cos(omega)} (C{scalar}).
           @kwarg name: Optional name (C{str}).

           @return: A L{Vector3Tuple}C{(x, y, z)} on the surface.

           @raise TriaxialError: This triaxial is near-spherical.

           @see: Method L{Triaxial.reverseBetaOmega}, U{Triaxial ellipsoid coordinate
                 system<https://WikiPedia.org/wiki/Geodesics_on_an_ellipsoid#
                 Triaxial_ellipsoid_coordinate_system>} and U{expressions (23-25)<https://
                 www.Topo.Auth.GR/wp-content/uploads/sites/111/2021/12/09_Panou.pdf>}.
        '''
        t = self._radialTo3(sbeta, cbeta, somega, comega)
        return Vector3Tuple(*t, name=name)

    def forwardCartesian(self, x, y, z, name=NN, **normal_eps):
        '''Project a cartesian on this triaxial.

           @arg x: X coordinate along C{a}-axis (C{meter}, same units as the axes).
           @arg y: Y coordinate along C{b}-axis (C{meter}, same units as B{C{x}}).
           @arg z: Z coordinate along C{c}-axis (C{meter}, same units as B{C{x}}).
           @kwarg name: Optional name (C{str}).
           @kwarg normal_eps: Optional keyword arguments C{B{normal}=True} and
                              C{B{eps}=EPS}, see method L{Triaxial.height4}.

           @see: Method L{Triaxial.height4} for further information and method
                 L{Triaxial.reverseCartesian} to reverse the projection.
        '''
        t = self.height4(x, y, z, **normal_eps)
        _ = t.rename(name)
        return t

    def forwardLatLon(self, lat, lon, height=0, name=NN):
        '''Convert I{geodetic} lat-, longitude and heigth to cartesian.

           @arg lat: Geodetic latitude (C{degrees}).
           @arg lon: Geodetic longitude (C{degrees}).
           @arg height: Height above the ellipsoid (C{meter}, same units
                        as this triaxial's C{a}, C{b} and C{c} axes).
           @kwarg name: Optional name (C{str}).

           @return: A L{Vector3Tuple}C{(x, y, z)}.

           @see: Method L{Triaxial.reverseLatLon} and U{Expressions (9-11)<https://
                 www.Topo.Auth.GR/wp-content/uploads/sites/111/2021/12/09_Panou.pdf>}.
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

           @see: Method L{Triaxial.reverseLatLon} and U{Expressions (9-11)<https://
                 www.Topo.Auth.GR/wp-content/uploads/sites/111/2021/12/09_Panou.pdf>}.
        '''
        sa, ca = self._norm2(slat, clat)
        sb, cb = self._norm2(slon, clon)
        return self._forward3(height, name, sa, ca, sb, cb)

    def _forward3(self, h, name, sa, ca, sb, cb):
        '''(INTERNAL) Helper for C{.forwardLatLon} and C{.forwardLatLon_}.
        '''
        ca_x_sb = ca * sb
        # 1 - (1 - (c/a)**2) * sa**2 - (1 - (b/a)**2) * ca**2 * sb**2
        t = fsum_(_1_0, -self.e2ac * sa**2, -self.e2ab * ca_x_sb**2, floats=True)
        N = self.a / self._sqrt(t)  # prime vertical
        x = (h + N)               * ca * cb
        y = (h + N * self._1e2ab) * ca_x_sb
        z = (h + N * self._1e2ac) * sa
        return Vector3Tuple(x, y, z, name=name)

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

    def _radialTo3(self, sbeta, cbeta, somega, comega):
        '''(INTERNAL) Convert I{ellipsoidal} lat- and longitude C{beta} and
           C{omega} to cartesian coordinates I{on the triaxial's surface},
           also I{ordered} helper for C{.height4}.
        '''
        sa, ca  =  self._norm2(sbeta,  cbeta)
        sb, cb  =  self._norm2(somega, comega)

        b2_a2   =  self._1e2ab  # ==  (b/a)**2
        c2_a2   = -self._1e2ac  # == -(c/a)**2
        a2c2_a2 =  self.  e2ac  # (a**2 - c**2) / a**2 == 1 - (c/a)**2

        x = Fsum(_1_0, -b2_a2 * sa**2, c2_a2 * ca**2).fover(a2c2_a2)
        z = Fsum(c2_a2,         sb**2, b2_a2 * cb**2).fover(a2c2_a2)

        x = self.a * cb * self._sqrt(x)
        y = self.b * ca * sb
        z = self.c * sa * self._sqrt(z)
        return x, y, z

    def reverseBetaOmega(self, x, y, z, name=NN):
        '''Convert cartesian to I{ellipsoidal} lat- and longitude, C{beta}, C{omega}
           and height.

           @arg x: X coordinate along C{a}-axis (C{meter}, same units as the axes).
           @arg y: Y coordinate along C{b}-axis (C{meter}, same units as B{C{x}}).
           @arg z: Z coordinate along C{c}-axis (C{meter}, same units as B{C{x}}).
           @kwarg name: Optional name (C{str}).

           @return: A L{BetaOmega3Tuple}C{(beta, omega, height)} with C{beta} and
                    C{omega} in C{Radians} and (radial) C{height} in C{meter}, same
                    units as this triaxial's axes.

           @see: Methods L{Triaxial.forwardBetaOmega} and L{Triaxial.forwardBetaOmega_}
                 and U{Expressions (21-22)<https://www.Topo.Auth.GR/wp-content/uploads/
                 sites/111/2021/12/09_Panou.pdf>}.
        '''
        v = Vector3d(x, y, z)
        a, b, h = self._reverse3(x, y, z, atan2, v, self.forwardBetaOmega_)
        return BetaOmega3Tuple(Radians(beta=a), Radians(omega=b), h, name=name)

    def reverseCartesian(self, x, y, z, h, normal=True, eps=_EPS2e4, name=NN):
        '''"Unproject" a cartesian on to a cartesion off this triaxial's surface.

           @arg x: X coordinate along C{a}-axis (C{meter}, same units as the axes).
           @arg y: Y coordinate along C{b}-axis (C{meter}, same units as B{C{x}}).
           @arg z: Z coordinate along C{c}-axis (C{meter}, same units as B{C{x}}).
           @arg h: Height above or below this triaxial's surface (C{meter}, same units
                   as B{C{x}}).
           @kwarg normal: If C{True} the height is C{normal} to the surface, otherwise
                          C{radially} to the center of this triaxial (C{bool}).
           @kwarg eps: Tolerance for surface test (C{scalar}).
           @kwarg name: Optional name (C{str}).

           @return: A L{Vector3Tuple}C{(x, y, z)}.

           @see: Methods L{Triaxial.forwardCartesian} and L{Triaxial.height4}.
        '''
        v = Vector3d(x, y, z, name=name)
        s = self.sideOf(v.x, v.y, v.z, eps=eps)
        if s:  # PYCHOK no cover
            t = _SPACE_((_inside_ if s < 0 else _outside_), self)
            raise TriaxialError(eps=eps, sideOf=s, x=x, y=y, z=z, txt=t)
        if h:
            if normal:
                v = v.plus(self.normal3d(v.x, v.y, v.z, length=h))
            elif v.length > EPS0:
                v = v.times(_1_0 + (h / v.length))
        return v.xyz  # Vector3Tuple

    def reverseLatLon(self, x, y, z, name=NN):
        '''Convert cartesian to I{geodetic} lat-, longitude and height.

           @arg x: X coordinate along C{a}-axis (C{meter}, same units as the axes).
           @arg y: Y coordinate along C{b}-axis (C{meter}, same units as B{C{x}}).
           @arg z: Z coordinate along C{c}-axis (C{meter}, same units as B{C{x}}).
           @kwarg name: Optional name (C{str}).

           @return: A L{LatLon3Tuple}C{(lat, lon, height)} with C{lat} and C{lon}
                    in C{degrees} and (radial) C{height} in C{meter}, same units
                    as this triaxial's axes.

           @see: Methods L{Triaxial.forwardLatLon} and L{Triaxial.forwardLatLon_}
                 and U{Expressions (4-5)<https://www.Topo.Auth.GR/wp-content/uploads/
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


class JacobiConformal(Triaxial):
    '''This is a conformal projection of a triaxial ellipsoid to a plane in which the
       C{X} and C{Y} grid lines are straight.

       Ellipsoidal coordinates I{beta} and I{omega} are converted to Jacobi Conformal
       I{y} respectively I{x} separately.  Jacobi's coordinates have been multiplied
       by C{sqrt(B{a}**2 - B{c}**2) / (2 * B{b})} so that the customary results are
       returned in the case of an ellipsoid of revolution (or a sphere, I{currently
       not supported}).

       Copyright (C) U{Charles Karney<mailto:Charles@Karney.com>} (2014-2020) and
       licensed under the MIT/X11 License.

       @note: This constructor can not be used to specify a sphere.

       @see: L{Triaxial}, C++ class U{JacobiConformal<https://GeographicLib.SourceForge.io/
             C++/doc/classGeographicLib_1_1JacobiConformal.html#details>}, U{Jacobi's conformal
             projection<https://GeographicLib.SourceForge.io/C++/doc/jacobi.html>} and Jacobi,
             C. G. J. I{U{Vorlesungen über Dynamik<https://Books.Google.com/books?
             id=ryEOAAAAQAAJ&pg=PA212>}}, page 212ff,
    '''
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
        return self._Elliptic(k2, _1_0 - self._a2_b2, kp2, self._a2_b2)

    def xR(self, omega):
        '''Compute a Jacobi Conformal C{x} projection.

           @arg omega: Ellipsoidal longitude (C{radians} or L{Degrees}).

           @return: The C{x} projection (C{Radians}).
        '''
        return self.xR_(*_SinCos2(omega))

    def xR_(self, somega, comega):
        '''Compute a Jacobi Conformal C{x} projection.

           @arg somega: Ellipsoidal longitude C{sin(omega)} (C{scalar}).
           @arg comega: Ellipsoidal longitude C{cos(omega)} (C{scalar}).

           @return: The C{x} projection (C{Radians}).
        '''
        s, c = self._norm2(somega, comega, self.a)
        return Radians(x=self._xE.fPi(s, c) * self._a2_b2)

    @Property_RO
    def xyQ2(self):
        '''Get the Jacobi Conformal quadrant size (L{Jacobi2Tuple}C{(x, y)}).
        '''
        return Jacobi2Tuple(Radians(x=self._a2_b2 * self._xE.cPi),
                            Radians(y=self._c2_b2 * self._yE.cPi),
                            name=JacobiConformal.xyQ2.name)

    def xyR2(self, beta, omega, name=NN):
        '''Compute a Jacobi Conformal C{x} and C{y} projection.

           @arg beta: Ellipsoidal latitude (C{radians} or L{Degrees}).
           @arg omega: Ellipsoidal longitude (C{radians} or L{Degrees}).
           @kwarg name: Optional name (C{str}).

           @return: A L{Jacobi2Tuple}C{(x, y)}.
        '''
        return self.xyR2_(*(_SinCos2(beta) + _SinCos2(omega)),
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
    def _yE(self):
        '''(INTERNAL) Get the x-elliptic function.
        '''
        kp2, k2 = self._k2_kp2  # swapped!
        # b2c2 / b2 == (b2 - c2) / b2 == 1 - c2 / b2 == e2bc
        return self._Elliptic(k2, self.e2bc, kp2, self._c2_b2)

    def yR(self, beta):
        '''Compute a Jacobi Conformal C{y} projection.

           @arg beta: Ellipsoidal latitude (C{radians} or L{Degrees}).

           @return: The C{y} projection (C{Radians}).
        '''
        return self.yR_(*_SinCos2(beta))

    def yR_(self, sbeta, cbeta):
        '''Compute a Jacobi Conformal C{y} projection.

           @arg sbeta: Ellipsoidal latitude C{sin(beta)} (C{scalar}).
           @arg cbeta: Ellipsoidal latitude C{cos(beta)} (C{scalar}).

           @return: The C{y} projection (C{Radians}).
        '''
        s, c = self._norm2(sbeta, cbeta, self.c)
        return Radians(y=self._yE.fPi(s, c) * self._c2_b2)


class TriaxialError(_ValueError):
    '''Raised for L{Triaxial} issues.
    '''
    pass  # ...


def _hartzell3d2(pov, los, abc3, Error=TriaxialError):  # MCCABE 13 in .formy.hartzell
    '''(INTERNAL) Hartzell's "Satellite Lin-of-Sight Intersection ...", I{unordered}.
    '''
    def _order2(a, b, c):
        # @return: 2-Tuple ((a, b, c), un) ordered a >= b >= c
        if a < b or b < c:
            t = (a, 0), (b, 1), (c, 2)
            return _zip(*reversed(sorted(t)))  # strict=True
        else:  # a >= b >= c
            return (a, b, c), None  # or ()

    def _order3d(un, v):
        # @return: Vector3d(x, y, z) un/ordered
        if un:
            _, xyz = _zip(*sorted(_zip(un, v.xyz)))  # strict=True
            v = Vector3d(*xyz)
        return v

    (a, b, c), un = _order2(*abc3)
    if not (isfinite(a) and c > 0):
        raise Error(_invalid_)

    a2     =   a**2  # largest, factored out
    b2, p2 =  (b**2, (b / a)**2) if a != b else (a2, _1_0)
    c2, q2 = ((c**2, (c / a)**2) if a != c else (a2, _1_0)) if b != c else (b2, _1_0)

    p3 = _order3d(un, _otherV3d(pov=pov))
    u3 = _order3d(un, _otherV3d(los=los)) if los else p3.negate()
    u3 =  u3.unit()  # unit vector, opposing signs

    x2, y2, z2 = p3.x2y2z2  # p3.times_(p3).xyz
    ux, vy, wz = u3.times_(p3).xyz
    u2, v2, w2 = u3.x2y2z2  # u3.times_(u3).xyz

    t = (p2 * c2),  c2, b2
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

    d = Fdot(t, ux, vy, wz).fadd_(r).fover(m)  # -r for antipode, a2 factored out
    if d > 0:  # POV inside or LOS missing, outside the triaxial
        raise Error(_inside_ if min(x2 - a2, y2 - b2, z2 - c2) < EPS else _outside_)
    elif fsum_(x2, y2, z2, floats=True) < d**2:  # d past triaxial's center
        raise Error(_too_(_distant_))

    v = p3.minus(u3.times(d))  # Vector3d
    h = p3.minus(v).length  # distance to triaxial
    return _order3d(un, v), h


E = _WGS84.ellipsoid
_WGS84 = Triaxial(E.a + 35, E.a - 35, E.b, name=E.name + '+/-35')
del E


def hartzell4(pov, los=None, tri_biax=_WGS84, name=NN):
    '''Compute the intersection of a tri-/biaxial ellipsoid and a Line-Of-Sight
       from a Point-Of-View in space.

       @arg pov: Point-Of-View outside the tri-/biaxial (C{Cartesian}, L{Ecef9Tuple}
                 or L{Vector3d}).
       @kwarg los: Line-Of-Sight, I{direction} to the tri-/biaxial (L{Vector3d}) or
                   C{None} to point to the tri-/biaxial's center.
       @kwarg tri_biax: A triaxial (L{Triaxial}, L{Triaxial_}, L{JacobiConformal})
                        or biaxial ellipsoid (L{Datum}, L{Ellipsoid}, L{Ellipsoid2},
                        L{a_f2Tuple} or C{scalar} radius in C{meter}).
       @kwarg name: Optional name (C{str}).

       @return: L{Vector4Tuple}C{(x, y, z, h)} on the tri-/biaxial's surface, with C{h}
                the distance from B{C{pov}} to C{(x, y, z)} along B{C{los}}.

       @raise TriaxialError: Null B{C{pov}} or B{C{los}} vector, B{C{pov}} is inside
                             the tri-/biaxial or B{C{los}} points outside the
                             tri-/biaxial or points in an opposite direction.

       @raise TypeError: Invalid B{C{pov}} or B{C{los}}.

       @see: Function L{pygeodesy.hartzell}, L{pygeodesy.tyr3d} for B{C{los}} and
             U{I{Satellite Line-of-Sight Intersection with Earth}<https://StephenHartzell.
             Medium.com/satellite-line-of-sight-intersection-with-earth-d786b4a6a9b6>}.
    '''
    def _Error(txt):
        return TriaxialError(pov=pov, los=los, tri_biax=tri_biax, txt=txt)

    if isinstance(tri_biax, Triaxial_):
        T = tri_biax
    else:
        D = tri_biax if isinstance(tri_biax, Datum) else \
                  _spherical_datum(tri_biax, name=hartzell4.__name__)
        E = D.ellipsoid
        T = Triaxial_(E.a, E.a, E.b, name=E.name)

    v, h = _hartzell3d2(pov, los, T._abc3, _Error)
    return Vector4Tuple(v.x, v.y, v.z, h, name=name or hartzell4.__name__)


def _normalTo4(x, y, a, b, eps=EPS):  # MCCABE 13
    '''(INTERNAL) Nearest point on and distance to an I{ordered} 2-D ellipse.

       @see: Function C{pygeodesy.ellipsoids._normalTo3} and I{Eberly}'s U{Distance
             from a Point to ... an Ellipsoid ...<https://www.GeometricTools.com/
             Documentation/DistancePointEllipseEllipsoid.pdf>}.
    '''
    def _root2d(r, u, v, g, eps):
        # robust root finder
        _1, __2 = _1_0, _0_5
        _a, _h2 = fabs, hypot2
        u *=  r
        t0 =  v - _1
        t1 = _0_0 if g < 0 else (hypot(u, v) - _1)
        for i in range(1, _TRIPS):
            t =   (t0 + t1) * __2
            e = _a(t0 - t1)
            if t in (t0, t1) or e < eps:
                break
            g = _h2(u / (t + r), v / (t + _1)) - _1
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

    if not (a >= b > 0):
        raise TriaxialError(a=a, b=b, txt=_not_ordered_)

    i = None
    if y:
        if x:
            u = fabs(x / a)
            v = fabs(y / b)
            g = hypot2(u, v) - _1_0
            if g:
                r = (a / b)**2
                t, i = _root2d(r, u, v, g, eps)
                a = x / (t / r + _1_0)
                b = y / (t     + _1_0)
                d = hypot(x - a, y - b)
            else:  # on the ellipse
                a, b, d = x, y, _0_0
        else:  # x == 0
            if y < 0:
                b = -b
            a, d = x, fabs(y - b)

    else:  # y == 0
        n =  a * x
        d = (a + b) * (a - b)
        if d > fabs(n):  # PYCHOK no cover
            r  = n / d
            a *= r
            b *= sqrt(_1_0 - r**2)
            d  = hypot(x - a, b)
        else:
            if x < 0:
                a = -a
            b, d = y, fabs(x - a)
    return a, b, d, i


def _normalTo4u(x, y, a, b, eps=EPS):  # PYCHOK no cover
    '''(INTERNAL) Nearest point on and distance to an I{unordered} 2-D ellipse.
    '''
    if a < b:
        b, a, d, i = _normalTo4(b, a, y, x, eps=eps)
    else:
        a, b, d, i = _normalTo4(a, b, x, y, eps=eps)
    return a, b, d, i


def _normalTo5(x, y, z, T, eps=EPS):  # MCCABE 23
    '''(INTERNAL) Nearest point on and distance to an I{ordered} triaxial.

       @see: I{Eberly}'s U{Distance from a Point to ... an Ellipsoid ...<https://
             www.GeometricTools.com/Documentation/DistancePointEllipseEllipsoid.pdf>}.
    '''
    def _root3d(r, s, u, v, w, g, eps):
        # robust root finder
        _1, __2 = _1_0, _0_5
        _a, _h2 = fabs, hypot2_
        u *=  r
        v *=  s
        t0 =  w - _1
        t1 = _0_0 if g < 0 else (hypot_(u, v, w) - _1)
        for i in range(1, _TRIPS):
            t =   (t0 + t1) * __2
            e = _a(t0 - t1)
            if t in (t0, t1) or e < eps:
                break
            g = _h2(u / (t + r), v / (t + s), w / (t + _1)) - _1
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

    a, b, c = T._abc3
    if not (a >= b >= c > 0):
        raise TriaxialError(a=a, b=b, c=c, triaxial=T, txt=_not_ordered_)

    if eps > 0:
        val = max(eps * 1e8, EPS)
    else:  # no validation
        val, eps = 0, -eps

    i = None
    if z:
        if y:
            if x:
                u = fabs(x / a)
                v = fabs(y / b)
                w = fabs(z / c)
                g = hypot2_(u, v, w) - _1_0
                if g:
                    r = T._1e2ac  # (c / a)**2
                    s = T._1e2bc  # (c / b)**2
                    t, i = _root3d(_1_0 / r, _1_0 / s, u, v, w, g, eps)
                    a = x / (t * r + _1_0)
                    b = y / (t * s + _1_0)
                    c = z / (t     + _1_0)
                    d = hypot_(x - a, y - b, z - c)
                else:  # on the ellipsoid
                    a, b, c, d = x, y, z, _0_0
            else:  # x == 0
                a          =  x  # 0
                b, c, d, i = _normalTo4(y, z, b, c, eps=eps)
        elif x:  # y == 0
            b          =  y  # 0
            a, c, d, i = _normalTo4(x, z, a, c, eps=eps)
        else:  # x == y == 0
            if z < 0:
                c = -c
            a, b, d = x, y, fabs(z - c)

    else:  # z == 0
        t = False
        n = a * x
        d = T._a2c2  # (a + c) * (a - c)
        if d > fabs(n):
            u = n / d
            n = b * y
            d = T._b2c2  # (b + c) * (b - c)
            if d > fabs(n):
                v =  n / d
                d = _1_0 - hypot2(u, v)
                if d > 0:
                    a *= u
                    b *= v
                    c *= sqrt(d)
                    d  = hypot_(x - a, y - b, c)
                    t  = True
        if not t:
            c          =  z  # 0
            a, b, d, i = _normalTo4(x, y, a, b, eps=eps)

    if val > 0:  # validate
        e = T.sideOf(a, b, c, eps=val)
        if e:  # not near the ellipsoid's surface
            raise TriaxialError(a=a, b=b, c=c, d=d,
                                x=x, y=y, z=z, eps=val,
                                sideOf=e, txt=T.toRepr())
        if d:  # angle of delta and normal vector
            m = Vector3d(x, y, z).minus_(a, b, c)
            if m.euclid > val:
                m = m.unit()
                n = T.normal3d(a, b, c)
                e = n.dot(m)  # n.negate().dot(m)
                if not isnear1(fabs(e), eps1=val):
                    raise TriaxialError(n=n, m=m, eps=val,
                                        dot=e, txt=T.toRepr())
    return a, b, c, d, i


def _normalTo5u(x, y, z, a, b, c, eps=EPS):
    '''(INTERNAL) Nearest point on and distance to an I{unordered} 3-D, triaxial ellipsoid.
    '''
    def _order3(a, b, c, x, y, z):
        # @return: 3-Tuple ((a, b, c), un, (x, y, z)), ordered a >= b >= c
        t = (a, 0, x), (b, 1, y), (c, 2, z)
        return _zip(*reversed(sorted(t)))  # strict=True

    def _unorder5(un, abcdi):
        # @return: 5-Tuple (a, b, c, d, i), unordered
        un += 3, 4  # keep d and i in place
        return tuple(_zip(*sorted(_zip(un, abcdi))))[1]  # strict=True

    abc, un, xyzT = _order3(a, b, c, x, y, z)
    xyzT += Triaxial(*abc),
    return _unorder5(un, _normalTo5(*xyzT, eps=eps))


def _SinCos2(x):
    '''Get C{sin} and C{cos} of C{x} in C{Degrees}, C{Radians} or {radians}.
    '''
    return sincos2d(x) if isinstance(x, Degrees) else (
           sincos2(x)  if isinstance(x, Radians) else
           sincos2(float(x)))  # assume C{radians}


if __name__ == '__main__':

    from pygeodesy import printf

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
        t = Triaxial(km2m(a), km2m(b), km2m(c), name=n)
        printf('# %r', t)
        if n == 'Earth':
            printf('# %r', JacobiConformal(t.a, t.b, t.c, name=n))
    printf('# %r', _WGS84)

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

# % python3 -m pygeodesy.triaxials

# Triaxial(name='Amalthea', a=125000, b=73000, c=64000, e2ab=0.658944, e2bc=0.231375493, e2ac=0.737856, volume=2446253479595252, area=93239507787.490371704, area_p=93212299402.670425415)
# Triaxial(name='Ariel', a=581100, b=577900, c=577700, e2ab=0.01098327, e2bc=0.000692042, e2ac=0.011667711, volume=812633172614203904, area=4211301462766.580078125, area_p=4211301574065.829589844)
# Triaxial(name='Earth', a=6378173.435, b=6378103.9, c=6356754.399999999, e2ab=0.000021804, e2bc=0.006683418, e2ac=0.006705077, volume=1083208241574987694080, area=510065911057441.0625, area_p=510065915922713.6875)
# JacobiConformal(name='Earth', a=6378173.435, b=6378103.9, c=6356754.399999999, e2ab=0.000021804, e2bc=0.006683418, e2ac=0.006705077, xyQ2=xyQ2(x=1.572084, y=4.249876), volume=1083208241574987694080, area=510065911057441.0625, area_p=510065915922713.6875)
# Triaxial(name='Enceladus', a=256600, b=251400, c=248300, e2ab=0.040119337, e2bc=0.024509841, e2ac=0.06364586, volume=67094551514082248, area=798618496278.596679688, area_p=798619018175.109863281)
# Triaxial(name='Europa', a=1564130, b=1561230, c=1560930, e2ab=0.003704694, e2bc=0.000384275, e2ac=0.004087546, volume=15966575194402123776, area=30663773697323.51953125, area_p=30663773794562.45703125)
# Triaxial(name='Io', a=1829400, b=1819300, c=1815700, e2ab=0.011011391, e2bc=0.003953651, e2ac=0.014921506, volume=25313121117889765376, area=41691875849096.7421875, area_p=41691877397441.2109375)
# Triaxial(name='Mars', a=3394600, b=3393300, c=3376300, e2ab=0.000765776, e2bc=0.009994646, e2ac=0.010752768, volume=162907283585817247744, area=144249140795107.4375, area_p=144249144150662.15625)
# Triaxial(name='Mimas', a=207400, b=196800, c=190600, e2ab=0.09960581, e2bc=0.062015624, e2ac=0.155444317, volume=32587072869017956, area=493855762247.691894531, area_p=493857714107.9375)
# Triaxial(name='Miranda', a=240400, b=234200, c=232900, e2ab=0.050915557, e2bc=0.011070811, e2ac=0.061422691, volume=54926187094835456, area=698880863325.756958008, area_p=698881306767.950317383)
# Triaxial(name='Moon', a=1735550, b=1735324, c=1734898, e2ab=0.000260419, e2bc=0.000490914, e2ac=0.000751206, volume=21886698675223740416, area=37838824729886.09375, area_p=37838824733332.2265625)
# Triaxial(name='Tethys', a=535600, b=528200, c=525800, e2ab=0.027441672, e2bc=0.009066821, e2ac=0.036259685, volume=623086233855821440, area=3528073490771.394042969, area_p=3528074261832.738769531)
# Triaxial(name='WGS84+/-35', a=6378172, b=6378102, c=6356752.314245179, e2ab=0.00002195, e2bc=0.006683478, e2ac=0.006705281, volume=1083207319768789942272, area=510065621722018.125, area_p=510065626587483.3125)
