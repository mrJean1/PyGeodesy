
# -*- coding: utf-8 -*-

u'''Triaxal ellipsoid classes I{ordered} L{Triaxial} and I{unordered} L{Triaxial_} and Jacobi
conformal projections L{JacobiConformal} and L{JacobiConformalSpherical}, transcoded from
I{Charles Karney}'s C++ class U{JacobiConformal<https://GeographicLib.SourceForge.io/C++/doc/
classGeographicLib_1_1JacobiConformal.html#details>} to pure Python and miscellaneous classes
L{BetaOmega2Tuple}, L{BetaOmega3Tuple}, L{Jacobi2Tuple} and L{TriaxialError}.

Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2008-2024).  For more information,
see the U{GeographicLib<https://GeographicLib.SourceForge.io>} documentation.

@see: U{Geodesics on a triaxial ellipsoid<https://WikiPedia.org/wiki/Geodesics_on_an_ellipsoid#
      Geodesics_on_a_triaxial_ellipsoid>} and U{Triaxial coordinate systems and their geometrical
      interpretation<https://www.Topo.Auth.GR/wp-content/uploads/sites/111/2021/12/09_Panou.pdf>}.

@var Triaxials.Amalthea: Triaxial(name='Amalthea', a=125000, b=73000, c=64000, e2ab=0.658944, e2bc=0.231375493, e2ac=0.737856, volume=2446253479595252, area=93239507787.490371704, area_p=93212299402.670425415)
@var Triaxials.Ariel: Triaxial(name='Ariel', a=581100, b=577900, c=577700, e2ab=0.01098327, e2bc=0.000692042, e2ac=0.011667711, volume=812633172614203904, area=4211301462766.580078125, area_p=4211301574065.829589844)
@var Triaxials.Earth: Triaxial(name='Earth', a=6378173.435, b=6378103.9, c=6356754.399999999, e2ab=0.000021804, e2bc=0.006683418, e2ac=0.006705077, volume=1083208241574987694080, area=510065911057441.0625, area_p=510065915922713.6875)
@var Triaxials.Enceladus: Triaxial(name='Enceladus', a=256600, b=251400, c=248300, e2ab=0.040119337, e2bc=0.024509841, e2ac=0.06364586, volume=67094551514082248, area=798618496278.596679688, area_p=798619018175.109863281)
@var Triaxials.Europa: Triaxial(name='Europa', a=1564130, b=1561230, c=1560930, e2ab=0.003704694, e2bc=0.000384275, e2ac=0.004087546, volume=15966575194402123776, area=30663773697323.51953125, area_p=30663773794562.45703125)
@var Triaxials.Io: Triaxial(name='Io', a=1829400, b=1819300, c=1815700, e2ab=0.011011391, e2bc=0.003953651, e2ac=0.014921506, volume=25313121117889765376, area=41691875849096.7421875, area_p=41691877397441.2109375)
@var Triaxials.Mars: Triaxial(name='Mars', a=3394600, b=3393300, c=3376300, e2ab=0.000765776, e2bc=0.009994646, e2ac=0.010752768, volume=162907283585817247744, area=144249140795107.4375, area_p=144249144150662.15625)
@var Triaxials.Mimas: Triaxial(name='Mimas', a=207400, b=196800, c=190600, e2ab=0.09960581, e2bc=0.062015624, e2ac=0.155444317, volume=32587072869017956, area=493855762247.691894531, area_p=493857714107.9375)
@var Triaxials.Miranda: Triaxial(name='Miranda', a=240400, b=234200, c=232900, e2ab=0.050915557, e2bc=0.011070811, e2ac=0.061422691, volume=54926187094835456, area=698880863325.756958008, area_p=698881306767.950317383)
@var Triaxials.Moon: Triaxial(name='Moon', a=1735550, b=1735324, c=1734898, e2ab=0.000260419, e2bc=0.000490914, e2ac=0.000751206, volume=21886698675223740416, area=37838824729886.09375, area_p=37838824733332.2265625)
@var Triaxials.Tethys: Triaxial(name='Tethys', a=535600, b=528200, c=525800, e2ab=0.027441672, e2bc=0.009066821, e2ac=0.036259685, volume=623086233855821440, area=3528073490771.394042969, area_p=3528074261832.738769531)
@var Triaxials.WGS84_35: Triaxial(name='WGS84_35', a=6378172, b=6378102, c=6356752.314245179, e2ab=0.00002195, e2bc=0.006683478, e2ac=0.006705281, volume=1083207319768789942272, area=510065621722018.125, area_p=510065626587483.3125)
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import isLatLon, isscalar
from pygeodesy.constants import EPS, EPS0, EPS02, EPS4, INT0, PI2, PI_3, PI4, \
                               _EPS2e4, _SQRT2_2, float0_, isfinite, isnear1, _over, \
                               _0_0, _0_5, _1_0, _N_1_0, _64_0,  _4_0  # PYCHOK used!
from pygeodesy.datums import Datum, _spherical_datum, _WGS84,  Ellipsoid, _EWGS84, Fmt
# from pygeodesy.ellipsoids import Ellipsoid, _EWGS84  # from .datums
# from pygeodesy.elliptic import Elliptic  # _MODS
from pygeodesy.errors import _AssertionError, _ValueError
from pygeodesy.fmath import Fdot, fdot, fmean_, hypot, hypot_, norm2, sqrt0
from pygeodesy.fsums import _Fsumf_, fsumf_, fsum1f_
from pygeodesy.interns import NN, _a_, _b_, _beta_, _c_, _distant_, _finite_, \
                             _height_, _inside_, _near_, _negative_, _not_, \
                             _NOTEQUAL_, _null_, _opposite_, _outside_, _SPACE_, \
                             _spherical_, _too_, _x_, _y_
# from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS  # from .vector3d
from pygeodesy.named import _lazyNamedEnumItem as _lazy, _name__, _NamedEnum, \
                                _NamedEnumItem, _Pass
from pygeodesy.namedTuples import LatLon3Tuple, _NamedTupleTo, Vector3Tuple, \
                                  Vector4Tuple
from pygeodesy.props import Property_RO, property_ROver
# from pygeodesy.streprs import Fmt  # from .datums
from pygeodesy.units import Degrees, Float, Height_, Meter, Meter2, Meter3, \
                            Radians, Radius, Scalar_
from pygeodesy.utily import asin1, atan2, atan2d, km2m, m2km, SinCos2, sincos2d_
from pygeodesy.vector3d import _otherV3d, Vector3d,  _ALL_LAZY, _MODS

from math import fabs, sqrt

__all__ = _ALL_LAZY.triaxials
__version__ = '24.11.24'

_not_ordered_ = _not_('ordered')
_omega_       = 'omega'
_TRIPS        =  359  # Eberly 1074?


class _NamedTupleToX(_NamedTupleTo):  # in .testNamedTuples
    '''(INTERNAL) Base class for L{BetaOmega2Tuple},
       L{BetaOmega3Tuple} and L{Jacobi2Tuple}.
    '''
    def _toDegrees(self, name, **toDMS_kwds):
        '''(INTERNAL) Convert C{self[0:2]} to L{Degrees} or C{toDMS}.
        '''
        return self._toX3U(_NamedTupleTo._Degrees3, Degrees, name, *self, **toDMS_kwds)

    def _toRadians(self, name):
        '''(INTERNAL) Convert C{self[0:2]} to L{Radians}.
        '''
        return self._toX3U(_NamedTupleTo._Radians3, Radians, name, *self)

    def _toX3U(self, _X3, U, name, a, b, *c, **kwds):
        a, b, s = _X3(self, a, b, **kwds)
        if s is None or name:
            n = self._name__(name)
            s = self.classof(a, b, *c, name=n).reUnit(U, U).toUnits()
        return s


class BetaOmega2Tuple(_NamedTupleToX):
    '''2-Tuple C{(beta, omega)} with I{ellipsoidal} lat- and
       longitude C{beta} and C{omega} both in L{Radians} (or
       L{Degrees}).
    '''
    _Names_ = (_beta_, _omega_)
    _Units_ = (_Pass,  _Pass)

    def toDegrees(self, name=NN, **toDMS_kwds):
        '''Convert this L{BetaOmega2Tuple} to L{Degrees} or C{toDMS}.

           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: L{BetaOmega2Tuple}C{(beta, omega)} with C{beta} and
                    C{omega} both in L{Degrees} or as L{toDMS} strings
                    provided some B{C{toDMS_kwds}} keyword arguments are
                    specified.
        '''
        return self._toDegrees(name, **toDMS_kwds)

    def toRadians(self, **name):
        '''Convert this L{BetaOmega2Tuple} to L{Radians}.

           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: L{BetaOmega2Tuple}C{(beta, omega)} with C{beta} and C{omega}
                    both in L{Radians}.
        '''
        return self._toRadians(name)


class BetaOmega3Tuple(_NamedTupleToX):
    '''3-Tuple C{(beta, omega, height)} with I{ellipsoidal} lat- and
       longitude C{beta} and C{omega} both in L{Radians} (or L{Degrees})
       and the C{height}, rather the (signed) I{distance} to the triaxial's
       surface (measured along the radial line to the triaxial's center)
       in C{meter}, conventionally.
    '''
    _Names_ = BetaOmega2Tuple._Names_ + (_height_,)
    _Units_ = BetaOmega2Tuple._Units_ + ( Meter,)

    def toDegrees(self, name=NN, **toDMS_kwds):
        '''Convert this L{BetaOmega3Tuple} to L{Degrees} or C{toDMS}.

           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: L{BetaOmega3Tuple}C{(beta, omega, height)} with
                    C{beta} and C{omega} both in L{Degrees} or as
                    L{toDMS} strings provided some B{C{toDMS_kwds}}
                    keyword arguments are specified.
        '''
        return self._toDegrees(name, **toDMS_kwds)

    def toRadians(self, **name):
        '''Convert this L{BetaOmega3Tuple} to L{Radians}.

           @kwarg name: Optional C{B{name}=NN} (C{str}), overriding this name.

           @return: L{BetaOmega3Tuple}C{(beta, omega, height)} with C{beta}
                    and C{omega} both in L{Radians}.
        '''
        return self._toRadians(name)

    def to2Tuple(self, **name):
        '''Reduce this L{BetaOmega3Tuple} to a L{BetaOmega2Tuple}.

           @kwarg name: Optional C{B{name}=NN} (C{str}), overriding this name.
        '''
        return BetaOmega2Tuple(*self[:2], name=self._name__(name))


class Jacobi2Tuple(_NamedTupleToX):
    '''2-Tuple C{(x, y)} with a Jacobi Conformal C{x} and C{y}
       projection, both in L{Radians} (or L{Degrees}).
    '''
    _Names_ = (_x_,   _y_)
    _Units_ = (_Pass, _Pass)

    def toDegrees(self, name=NN, **toDMS_kwds):
        '''Convert this L{Jacobi2Tuple} to L{Degrees} or C{toDMS}.

           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: L{Jacobi2Tuple}C{(x, y)} with C{x} and C{y} both
                    in L{Degrees} or as L{toDMS} strings provided some
                    B{C{toDMS_kwds}} keyword arguments are specified.
        '''
        return self._toDegrees(name, **toDMS_kwds)

    def toRadians(self, **name):
        '''Convert this L{Jacobi2Tuple} to L{Radians}.

           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: L{Jacobi2Tuple}C{(x, y)} with C{x} and C{y} both in L{Radians}.
        '''
        return self._toRadians(name)


class Triaxial_(_NamedEnumItem):
    '''I{Unordered} triaxial ellipsoid and base class.

       Triaxial ellipsoids with right-handed semi-axes C{a}, C{b} and C{c}, oriented
       such that the large principal ellipse C{ab} is the equator I{Z}=0, I{beta}=0,
       while the small principal ellipse C{ac} is the prime meridian, plane I{Y}=0,
       I{omega}=0.

       The four umbilic points, C{abs}(I{omega}) = C{abs}(I{beta}) = C{PI/2}, lie on
       the middle principal ellipse C{bc} in plane I{X}=0, I{omega}=C{PI/2}.

       @note: I{Geodetic} C{lat}- and C{lon}gitudes are in C{degrees}, I{geodetic}
              C{phi} and C{lam}bda are in C{radians}, but I{ellipsoidal} lat- and
              longitude C{beta} and C{omega} are in L{Radians} by default (or in
              L{Degrees} if converted).
    '''
    _ijk = _kji = None
    _unordered  = True

    def __init__(self, a_triaxial, b=None, c=None, **name):
        '''New I{unordered} L{Triaxial_}.

           @arg a_triaxial: Large, C{X} semi-axis (C{scalar}, conventionally in
                            C{meter}) or an other L{Triaxial} or L{Triaxial_} instance.
           @kwarg b: Middle, C{Y} semi-axis (C{meter}, same units as B{C{a}}), required
                     if C{B{a_triaxial} is scalar}, ignored otherwise.
           @kwarg c: Small, C{Z} semi-axis (C{meter}, like B{C{b}}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @raise TriaxialError: Invalid semi-axis or -axes.
        '''
        try:
            try:
                a = a_triaxial
                t = a._abc3
            except AttributeError:
                t = Radius(a=a), Radius(b=b), Radius(c=c)
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
        '''Get the C{largest, x} semi-axis (C{meter}, conventionally).
        '''
        a, _, _ = self._abc3
        return a

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
        '''Get the C{middle, y} semi-axis (C{meter}, same units as B{C{a}}).
        '''
        _, b, _ = self._abc3
        return b

    @Property_RO
    def _b2c2(self):
        '''(INTERNAL) Get C{b**2 - c**2} == E_sub_y**2.
        '''
        _, b, c = self._abc3
        return ((b - c) * (b + c)) if b != c else _0_0

    @Property_RO
    def c(self):
        '''Get the C{smallest, z} semi-axis (C{meter}, same units as B{C{a}}).
        '''
        _, _, c = self._abc3
        return c

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
    def e2bc(self):
        '''Get the C{bc} ellipse' I{(1st) eccentricity squared} (C{scalar}), M{1 - (c/b)**2}.
        '''
        return Float(e2bc=(_1_0 - self._1e2bc) or _0_0)

    _1e2bc = _c2_b2  # C{1 - e2bc} == C{(c/b)**2}

    @property_ROver
    def _Elliptic(self):
        '''(INTERNAL) Get class L{Elliptic}, I{once}.
        '''
        return _MODS.elliptic.Elliptic  # overwrite property_ROver

    def hartzell4(self, pov, los=False, **name):
        '''Compute the intersection of this triaxial's surface with a Line-Of-Sight
           from a Point-Of-View in space.

           @see: Function L{hartzell4<triaxials.hartzell4>} for further details.
        '''
        return hartzell4(pov, los=los, tri_biax=self, **name)

    def height4(self, x_xyz, y=None, z=None, normal=True, eps=EPS, **name):
        '''Compute the projection on and the height above or below this triaxial's surface.

           @arg x_xyz: X component (C{scalar}) or a cartesian (C{Cartesian}, L{Ecef9Tuple},
                       L{Vector3d}, L{Vector3Tuple} or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} if C{scalar}, ignored
                     otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg normal: If C{True}, the projection is the I{perpendicular, plumb} to the
                          triaxial's surface, otherwise the C{radial} line to the center of
                          this triaxial (C{bool}).
           @kwarg eps: Tolerance for root finding and validation (C{scalar}), use a negative
                       value to skip validation.
           @kwarg name: Optional C{B{name}="heigh4"} (C{str}).

           @return: L{Vector4Tuple}C{(x, y, z, h)} with the cartesian coordinates C{x}, C{y}
                    and C{z} of the projection on or the intersection with and with the height
                    C{h} above or below the triaxial's surface in C{meter}, conventionally.

           @raise TriaxialError: Non-cartesian B{C{xyz}}, invalid B{C{eps}}, no convergence in
                  root finding or validation failed.

           @see: Methods L{Triaxial.normal3d} and L{Ellipsoid.height4}, I{Eberly}'s U{Distance from a Point to ...
                 <https://www.GeometricTools.com/Documentation/DistancePointEllipseEllipsoid.pdf>} and I{Bektas}'
                 U{Shortest Distance from a Point to Triaxial Ellipsoid<https://www.ResearchGate.net/publication/
                 272149005_SHORTEST_DISTANCE_FROM_A_POINT_TO_TRIAXIAL_ELLIPSOID>}.
        '''
        v, r = _otherV3d_(x_xyz, y, z), self.isSpherical

        i, h = None, v.length
        if h < EPS0:  # EPS
            x = y = z = _0_0
            h -= min(self._abc3)  # nearest
        elif r:  # .isSpherical
            x, y, z = v.times(r / h).xyz3
            h -= r
        else:
            x, y, z = v.xyz3
            try:
                if normal:  # plumb to surface
                    x, y, z, h, i = _plumbTo5(x, y, z, self, eps=eps)
                else:  # radial to center
                    x, y, z = self._radialTo3(z, hypot(x, y), y, x)
                    h = v.minus_(x, y, z).length
            except Exception as e:
                raise TriaxialError(x=x, y=y, z=z, cause=e)
            if h > 0 and self.sideOf(v, eps=EPS0) < 0:
                h = -h  # inside
        n = _name__(name, name__=self.height4)  # _DUNDER_nameof
        return Vector4Tuple(x, y, z, h, iteration=i, name=n)

    @Property_RO
    def isOrdered(self):
        '''Is this triaxial I{ordered} and I{not spherical} (C{bool})?
        '''
        a, b, c = self._abc3
        return bool(a >= b > c)  # b > c!

    @Property_RO
    def isSpherical(self):
        '''Is this triaxial I{spherical} (C{Radius} or INT0)?
        '''
        a, b, c = self._abc3
        return a if a == b == c else INT0

    def _norm2(self, s, c, *a):
        '''(INTERNAL) Normalize C{s} and C{c} iff not already.
        '''
        if fabs(_hypot2_1(s, c)) > EPS02:
            s, c = norm2(s, c)
        if a:
            s, c = norm2(s * self.b, c * a[0])
        return float0_(s, c)

    def normal3d(self, x_xyz, y=None, z=None, length=_1_0):
        '''Get a 3-D vector at a cartesian I{on and perpendicular to} this triaxial's surface.

           @arg x_xyz: X component (C{scalar}) or a cartesian (C{Cartesian},
                       L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple} or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} if C{scalar}, ignored
                     otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg length: Optional, signed length and in-/outward direction (C{scalar}).

           @return: A C{Vector3d(x_, y_, z_)} normalized to B{C{length}}, pointing in- or
                    outward for neg- respectively positive B{C{length}}.

           @raise TriaxialError: Zero length cartesian or vector.

           @note: Cartesian C{(B{x}, B{y}, B{z})} I{must be on} this triaxial's surface, use
                  method L{Triaxial.sideOf} to validate.

           @see: Methods L{Triaxial.height4} and L{Triaxial.sideOf}.
        '''
        # n = 2 * (x / a2, y / b2, z / c2)
        #  == 2 * (x, y * a2 / b2, z * a2 / c2) / a2  # iff ordered
        #  == 2 * (x, y / _1e2ab, z / _1e2ac) / a2
        #  == unit(x, y / _1e2ab, z / _1e2ac).times(length)
        x, y, z = _otherV3d_(x_xyz, y, z).xyz3
        n = Vector3d(x, y / self._1e2ab,
                        z / self._1e2ac, name__=self.normal3d)
        u = n.length
        if u < EPS0:
            raise TriaxialError(x=x_xyz, y=y, z=z, txt=_null_)
        return n.times(length / u)

    def _order3(self, *abc, **reverse):  # reverse=False
        '''(INTERNAL) Un-/Order C{a}, C{b} and C{c}.

           @return: 3-Tuple C{(a, b, c)} ordered by or un-ordered
                    (reverse-ordered) C{ijk} if C{B{reverse}=True}.
        '''
        ijk = self._order_ijk(**reverse)
        return _getitems(abc, *ijk) if ijk else abc

    def _order3d(self, v, **reverse):  # reverse=False
        '''(INTERNAL) Un-/Order a C{Vector3d}.

           @return: Vector3d(x, y, z) un-/ordered.
        '''
        ijk = self._order_ijk(**reverse)
        return v.classof(*_getitems(v.xyz3, *ijk)) if ijk else v

    @Property_RO
    def _ordered4(self):
        '''(INTERNAL) Helper for C{_hartzell3} and C{_plumbTo5}.
        '''
        def _order2(reverse, a, b, c):
            '''(INTERNAL) Un-Order C{a}, C{b} and C{c}.

               @return: 2-Tuple C{((a, b, c), ijk)} with C{a} >= C{b} >= C{c}
                        and C{ijk} a 3-tuple with the initial indices.
            '''
            i, j, k = 0, 1, 2  # range(3)
            if a < b:
                a, b, i, j = b, a, j, i
            if a < c:
                a, c, i, k = c, a, k, i
            if b < c:
                b, c, j, k = c, b, k, j
            # reverse (k, j, i) since (a, b, c) is reversed-sorted
            ijk = (k, j, i) if reverse else (None if i < j < k else (i, j, k))
            return (a, b, c), ijk

        abc, T = self._abc3, self
        if not self.isOrdered:
            abc, ijk = _order2(False, *abc)
            if ijk:
                _, kji = _order2(True, *ijk)
                T = Triaxial_(*abc)
                T._ijk, T._kji = ijk, kji
        return abc + (T,)

    def _order_ijk(self, reverse=False):
        '''(INTERNAL) Get the un-/order indices.
        '''
        return self._kji if reverse else self._ijk

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
        t = c * ca
        return (t * cb), (t * sb), (c * sa)

    def sideOf(self, x_xyz, y=None, z=None, eps=EPS4):
        '''Is a cartesian on, above or below the surface of this triaxial?

           @arg x_xyz: X component (C{scalar}) or a cartesian (C{Cartesian},
                       L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple} or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} if C{scalar},
                     ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg eps: On-surface tolerance (C{scalar}, distance I{squared}).

           @return: C{INT0} if C{(B{x}, B{y}, B{z})} is near this triaxial's surface
                    within tolerance B{C{eps}}, otherwise the signed, radial distance
                    I{squared} (C{float}), negative for in- or positive for outside
                    this triaxial.

           @see: Methods L{Triaxial.height4} and L{Triaxial.normal3d}.
        '''
        v = _otherV3d_(x_xyz, y, z)
        s =  fsumf_(_N_1_0, *map(_over02, v.xyz3, self._abc3))
        return INT0 if fabs(s) < eps else s

    def toEllipsoid(self, **name):
        '''Convert this triaxial to an L{Ellipsoid}, provided 2 axes match.

           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: An L{Ellipsoid} with north along this C{Z} axis if C{a == b},
                    this C{Y} axis if C{a == c} or this C{X} axis if C{b == c}.

           @raise TriaxialError: This C{a != b}, C{b != c} and C{c != a}.

           @see: Method L{Ellipsoid.toTriaxial}.
        '''
        a, b, c = self._abc3
        if a == b:
            b = c  # N = c-Z
        elif b == c:  # N = a-X
            a, b = b, a
        elif a != c:  # N = b-Y
            t = _SPACE_(_a_, _NOTEQUAL_, _b_, _NOTEQUAL_, _c_)
            raise TriaxialError(a=a, b=b, c=c, txt=t)
        return Ellipsoid(a, b=b, name=self._name__(name))

    def toStr(self, prec=9, **name):  # PYCHOK signature
        '''Return this C{Triaxial} as a string.

           @kwarg prec: Precision, number of decimal digits (0..9).
           @kwarg name: Optional name (C{str}), to override or C{None}
                        to exclude this triaxial's name.

           @return: This C{Triaxial}'s attributes (C{str}).
        '''
        T  = Triaxial_
        t  = T.a,
        J  = JacobiConformalSpherical
        t += (J.ab, J.bc) if isinstance(self, J) else (T.b, T.c)
        t += T.e2ab, T.e2bc, T.e2ac
        J  = JacobiConformal
        if isinstance(self, J):
            t += J.xyQ2,
        t += T.volume, T.area
        return self._instr(area_p=self.area_p(), prec=prec, props=t, **name)

    @Property_RO
    def unOrdered(self):
        '''Is this triaxial I{un-ordered} and I{not spherical} (C{bool})?
        '''
        return not (self.isOrdered or bool(self.isSpherical))

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

    def __init__(self, a_triaxial, b=None, c=None, **name):
        '''New I{ordered} L{Triaxial}.

           @arg a_triaxial: Largest semi-axis (C{scalar}, conventionally in C{meter})
                            or an other L{Triaxial} or L{Triaxial_} instance.
           @kwarg b: Middle semi-axis (C{meter}, same units as B{C{a}}), required
                     if C{B{a_triaxial} is scalar}, ignored otherwise.
           @kwarg c: Smallest semi-axis (C{meter}, like B{C{b}}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @note: The semi-axes must be ordered as C{B{a} >= B{b} >= B{c} > 0} and
                  must be ellipsoidal, C{B{a} > B{c}}.

           @raise TriaxialError: Semi-axes unordered, spherical or invalid.
        '''
        Triaxial_.__init__(self, a_triaxial, b=b, c=c, **name)

    @Property_RO
    def _a2b2_a2c2(self):
        '''@see: Methods C{.forwardBetaOmega} and C{._k2_kp2}.
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
            c2 = self._1e2ac      # cos(phi)**2 = (c/a)**2
            s  = sqrt(self.e2ac)  # sin(phi)**2 =  1 - c2
            r  = asin1(s)  # phi = atan2(sqrt(c2), s)
            b *= fsum1f_(aE.fE(r) * s, (c / a) * (c / b),
                         aE.fF(r) * c2 / s)
            a  = Meter2(area=a * b * PI2)
        else:  # a == b > c
            a  = Ellipsoid(a, b=c).areax
        return a

    def forwardBetaOmega(self, beta, omega, height=0, **name):
        '''Convert I{ellipsoidal} lat- C{beta}, longitude C{omega} and C{height}
           to cartesian.

           @arg beta: Ellipsoidal latitude (C{radians} or L{Degrees}).
           @arg omega: Ellipsoidal longitude (C{radians} or L{Degrees}).
           @arg height: Height above or below the ellipsoid's surface (C{meter}, same
                        units as this triaxial's C{a}, C{b} and C{c} semi-axes).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: A L{Vector3Tuple}C{(x, y, z)}.

           @see: Method L{Triaxial.reverseBetaOmega} and U{expressions (23-25)<https://
                 www.Topo.Auth.GR/wp-content/uploads/sites/111/2021/12/09_Panou.pdf>}.
        '''
        if height:
            z = self._Height(height) + self.c
            if z > 0:
                z2 = z**2
                x  = z * _sqrt0(_1_0 + self._a2c2 / z2)
                y  = z * _sqrt0(_1_0 + self._b2c2 / z2)
            else:
                x  = y = z = _0_0
        else:
            x, y, z = self._abc3
        if z:  # and x and y:
            sa, ca = SinCos2(beta)
            sb, cb = SinCos2(omega)

            r  = self._a2b2_a2c2
            x *= cb * _sqrt0(ca**2 + sa**2 * r)
            y *= ca *  sb
            z *= sa * _sqrt0(_1_0  - cb**2 * r)
        return Vector3Tuple(x, y, z, **name)

    def forwardBetaOmega_(self, sbeta, cbeta, somega, comega, **name):
        '''Convert I{ellipsoidal} lat- and longitude C{beta} and C{omega}
           to cartesian coordinates I{on the triaxial's surface}.

           @arg sbeta: Ellipsoidal latitude C{sin(beta)} (C{scalar}).
           @arg cbeta: Ellipsoidal latitude C{cos(beta)} (C{scalar}).
           @arg somega: Ellipsoidal longitude C{sin(omega)} (C{scalar}).
           @arg comega: Ellipsoidal longitude C{cos(omega)} (C{scalar}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: A L{Vector3Tuple}C{(x, y, z)} on the surface.

           @raise TriaxialError: This triaxial is near-spherical.

           @see: Method L{Triaxial.reverseBetaOmega}, U{Triaxial ellipsoid coordinate
                 system<https://WikiPedia.org/wiki/Geodesics_on_an_ellipsoid#
                 Triaxial_ellipsoid_coordinate_system>} and U{expressions (23-25)<https://
                 www.Topo.Auth.GR/wp-content/uploads/sites/111/2021/12/09_Panou.pdf>}.
        '''
        t = self._radialTo3(sbeta, cbeta, somega, comega)
        return Vector3Tuple(*t, **name)

    def forwardCartesian(self, x_xyz, y=None, z=None, **normal_eps_name):
        '''Project a cartesian on this triaxial.

           @arg x_xyz: X component (C{scalar}) or a cartesian (C{Cartesian},
                       L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple} or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} if C{scalar},
                     ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg normal_eps_name: Optional keyword arguments C{B{normal}=True},
                         C{B{eps}=EPS} and overriding C{B{name}="height4"} (C{str}),
                         see method L{Triaxial.height4}.

           @see: Method L{Triaxial.height4} for further information and method
                 L{Triaxial.reverseCartesian} to reverse the projection.
        '''
        return self.height4(x_xyz, y, z, **normal_eps_name)

    def forwardLatLon(self, lat, lon, height=0, **name):
        '''Convert I{geodetic} lat-, longitude and height to cartesian.

           @arg lat: Geodetic latitude (C{degrees}).
           @arg lon: Geodetic longitude (C{degrees}).
           @arg height: Height above the ellipsoid (C{meter}, same units
                        as this triaxial's C{a}, C{b} and C{c} axes).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: A L{Vector3Tuple}C{(x, y, z)}.

           @see: Method L{Triaxial.reverseLatLon} and U{expressions (9-11)<https://
                 www.Topo.Auth.GR/wp-content/uploads/sites/111/2021/12/09_Panou.pdf>}.
        '''
        return self._forwardLatLon3(height, name, *sincos2d_(lat, lon))

    def forwardLatLon_(self, slat, clat, slon, clon, height=0, **name):
        '''Convert I{geodetic} lat-, longitude and height to cartesian.

           @arg slat: Geodetic latitude C{sin(lat)} (C{scalar}).
           @arg clat: Geodetic latitude C{cos(lat)} (C{scalar}).
           @arg slon: Geodetic longitude C{sin(lon)} (C{scalar}).
           @arg clon: Geodetic longitude C{cos(lon)} (C{scalar}).
           @arg height: Height above the ellipsoid (C{meter}, same units
                        as this triaxial's axes C{a}, C{b} and C{c}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: A L{Vector3Tuple}C{(x, y, z)}.

           @see: Method L{Triaxial.reverseLatLon} and U{expressions (9-11)<https://
                 www.Topo.Auth.GR/wp-content/uploads/sites/111/2021/12/09_Panou.pdf>}.
        '''
        sa, ca = self._norm2(slat, clat)
        sb, cb = self._norm2(slon, clon)
        return self._forwardLatLon3(height, name, sa, ca, sb, cb)

    def _forwardLatLon3(self, height, name, sa, ca, sb, cb):  # name always **name
        '''(INTERNAL) Helper for C{.forwardLatLon} and C{.forwardLatLon_}.
        '''
        ca_x_sb = ca * sb
        h = self._Height(height)
        # 1 - (1 - (c/a)**2) * sa**2 - (1 - (b/a)**2) * ca**2 * sb**2
        t = fsumf_(_1_0, -self.e2ac * sa**2, -self.e2ab * ca_x_sb**2)
        n = self.a / _sqrt0(t)  # prime vertical
        x = (h + n)               * ca * cb
        y = (h + n * self._1e2ab) * ca_x_sb
        z = (h + n * self._1e2ac) * sa
        return Vector3Tuple(x, y, z, **name)

    def _Height(self, height):
        '''(INTERNAL) Validate a C{height}.
        '''
        return Height_(height=height, low=-self.c, Error=TriaxialError)

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
        return (self._a2b2_a2c2         * self._c2_b2,
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

        x2 = _Fsumf_(_1_0, -b2_a2 * sa**2, c2_a2 * ca**2).fover(a2c2_a2)
        z2 = _Fsumf_(c2_a2,         sb**2, b2_a2 * cb**2).fover(a2c2_a2)

        x, y, z =  self._abc3
        x *= cb * _sqrt0(x2)
        y *= ca *  sb
        z *= sa * _sqrt0(z2)
        return x, y, z

    def reverseBetaOmega(self, x_xyz, y=None, z=None, **name):
        '''Convert cartesian to I{ellipsoidal} lat- and longitude, C{beta}, C{omega}
           and height.

           @arg x_xyz: X component (C{scalar}) or a cartesian (C{Cartesian},
                       L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple} or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} if C{scalar},
                     ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: A L{BetaOmega3Tuple}C{(beta, omega, height)} with C{beta} and
                    C{omega} in L{Radians} and (radial) C{height} in C{meter}, same
                    units as this triaxial's axes.

           @see: Methods L{Triaxial.forwardBetaOmega} and L{Triaxial.forwardBetaOmega_}
                 and U{expressions (21-22)<https://www.Topo.Auth.GR/wp-content/uploads/
                 sites/111/2021/12/09_Panou.pdf>}.
        '''
        v = _otherV3d_(x_xyz, y, z)
        a, b, h = self._reverseLatLon3(v, atan2, v, self.forwardBetaOmega_)
        return BetaOmega3Tuple(Radians(beta=a), Radians(omega=b), h, **name)

    def reverseCartesian(self, x_xyz, y=None, z=None, h=0, normal=True, eps=_EPS2e4, **name):
        '''Unproject" a cartesian I{on} this triaxial's surface to a cartesion I{off}
           this triaxial's surface.

           @arg x_xyz: X component (C{scalar}) or a cartesian (C{Cartesian},
                       L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple} or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} if C{scalar},
                       ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @arg h: Height above or below this triaxial's surface (C{meter}, same units
                   as this triaxial's axes).
           @kwarg normal: If C{True}, the height is C{normal} to the surface, otherwise
                          C{radially} to the center of this triaxial (C{bool}).
           @kwarg eps: Tolerance for on-surface test (C{scalar}), see method L{Triaxial.sideOf}.
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: A L{Vector3Tuple}C{(x, y, z)}.

           @raise TrialError: Cartesian C{(x, y, z)} not on this triaxial's surface.

           @see: Methods L{Triaxial.forwardCartesian} and L{Triaxial.height4}.
        '''
        v = _otherV3d_(x_xyz, y, z, **name)
        s =  self.sideOf(v.xyz, eps=eps)
        if s:  # PYCHOK no cover
            t = _SPACE_((_inside_ if s < 0 else _outside_), self.toRepr())
            raise TriaxialError(eps=eps, sideOf=s, x=v.x, y=v.y, z=v.z, txt=t)

        if h:
            if normal:
                v = v.plus(self.normal3d(*v.xyz, length=h))
            elif v.length > EPS0:
                v = v.times(_1_0 + (h / v.length))
        return v.xyz  # Vector3Tuple

    def reverseLatLon(self, x_xyz, y=None, z=None, **name):
        '''Convert cartesian to I{geodetic} lat-, longitude and height.

           @arg x_xyz: X component (C{scalar}) or a cartesian (C{Cartesian},
                       L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple} or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} if C{scalar},
                     ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: A L{LatLon3Tuple}C{(lat, lon, height)} with C{lat} and C{lon}
                    in C{degrees} and (radial) C{height} in C{meter}, same units
                    as this triaxial's axes.

           @see: Methods L{Triaxial.forwardLatLon} and L{Triaxial.forwardLatLon_}
                 and U{expressions (4-5)<https://www.Topo.Auth.GR/wp-content/uploads/
                 sites/111/2021/12/09_Panou.pdf>}.
        '''
        v = _otherV3d_(x_xyz, y, z)
        s =  v.times_(self._1e2ac,  # == 1 - e_sub_x**2
                      self._1e2bc,  # == 1 - e_sub_y**2
                     _1_0)
        a, b, h = self._reverseLatLon3(s, atan2d, v, self.forwardLatLon_)
        return LatLon3Tuple(Degrees(lat=a), Degrees(lon=b), h, **name)

    def _reverseLatLon3(self, s, atan2_, v, forward_):
        '''(INTERNAL) Helper for C{.reverseBetOmg} and C{.reverseLatLon}.
        '''
        x, y, z = s.xyz3
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
       returned in the case of an ellipsoid of revolution.

       Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2014-2024) and
       licensed under the MIT/X11 License.

       @note: This constructor can I{not be used to specify a sphere}, see alternate
              L{JacobiConformalSpherical}.

       @see: L{Triaxial}, C++ class U{JacobiConformal<https://GeographicLib.SourceForge.io/
             C++/doc/classGeographicLib_1_1JacobiConformal.html#details>}, U{Jacobi's conformal
             projection<https://GeographicLib.SourceForge.io/C++/doc/jacobi.html>} and Jacobi,
             C. G. J. I{U{Vorlesungen über Dynamik<https://Books.Google.com/books?
             id=ryEOAAAAQAAJ&pg=PA212>}}, page 212ff.
    '''

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

           @return: The C{x} projection (L{Radians}).
        '''
        return self.xR_(*SinCos2(omega))

    def xR_(self, somega, comega):
        '''Compute a Jacobi Conformal C{x} projection.

           @arg somega: Ellipsoidal longitude C{sin(omega)} (C{scalar}).
           @arg comega: Ellipsoidal longitude C{cos(omega)} (C{scalar}).

           @return: The C{x} projection (L{Radians}).
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

    def xyR2(self, beta, omega, **name):
        '''Compute a Jacobi Conformal C{x} and C{y} projection.

           @arg beta: Ellipsoidal latitude (C{radians} or L{Degrees}).
           @arg omega: Ellipsoidal longitude (C{radians} or L{Degrees}).
           @kwarg name: Optional name (C{str}), overriding C{B{name}="xyR2"}.

           @return: A L{Jacobi2Tuple}C{(x, y)}.
        '''
        return self.xyR2_(*(SinCos2(beta) + SinCos2(omega)),
                          name=_name__(name, name__=self.xyR2))

    def xyR2_(self, sbeta, cbeta, somega, comega, **name):
        '''Compute a Jacobi Conformal C{x} and C{y} projection.

           @arg sbeta: Ellipsoidal latitude C{sin(beta)} (C{scalar}).
           @arg cbeta: Ellipsoidal latitude C{cos(beta)} (C{scalar}).
           @arg somega: Ellipsoidal longitude C{sin(omega)} (C{scalar}).
           @arg comega: Ellipsoidal longitude C{cos(omega)} (C{scalar}).
           @kwarg name: Optional name (C{str}), overriding C{B{name}="xyR2_"}.

           @return: A L{Jacobi2Tuple}C{(x, y)}.
        '''
        return Jacobi2Tuple(self.xR_(somega, comega),
                            self.yR_(sbeta,  cbeta),
                            name=_name__(name, name__=self.xyR2_))

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

           @return: The C{y} projection (L{Radians}).
        '''
        return self.yR_(*SinCos2(beta))

    def yR_(self, sbeta, cbeta):
        '''Compute a Jacobi Conformal C{y} projection.

           @arg sbeta: Ellipsoidal latitude C{sin(beta)} (C{scalar}).
           @arg cbeta: Ellipsoidal latitude C{cos(beta)} (C{scalar}).

           @return: The C{y} projection (L{Radians}).
        '''
        s, c = self._norm2(sbeta, cbeta, self.c)
        return Radians(y=self._yE.fPi(s, c) * self._c2_b2)


class JacobiConformalSpherical(JacobiConformal):
    '''An alternate, I{spherical} L{JacobiConformal} projection.

       @see: L{JacobiConformal} for other and more details.
    '''
    _ab = _bc = 0

    def __init__(self, radius_triaxial, ab=0, bc=0, **name):
        '''New L{JacobiConformalSpherical}.

           @arg radius_triaxial: Radius (C{scalar}, conventionally in
                       C{meter}) or an other L{JacobiConformalSpherical},
                       L{JacobiConformal} or ordered L{Triaxial}.
           @kwarg ab: Relative magnitude of C{B{a} - B{b}} (C{meter},
                      same units as C{scalar B{radius}}.
           @kwarg bc: Relative magnitude of C{B{b} - B{c}} (C{meter},
                      same units as C{scalar B{radius}}.
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @raise TriaxialError: Invalid B{C{radius_triaxial}}, negative
                                 B{C{ab}}, negative B{C{bc}} or C{(B{ab}
                                 + B{bc})} not positive.

           @note: If B{C{radius_triaxial}} is a L{JacobiConformalSpherical}
                  and if B{C{ab}} and B{C{bc}} are both zero or C{None},
                  the B{C{radius_triaxial}}'s C{ab}, C{bc}, C{a}, C{b}
                  and C{c} are copied.
        '''
        try:
            r = radius_triaxial
            if isinstance(r, Triaxial):  # ordered only
                t = r._abc3
                j = isinstance(r, JacobiConformalSpherical) and not bool(ab or bc)
            else:
                t = (Radius(radius=r),) * 3
                j =  False
            self._ab = r.ab if j else Scalar_(ab=ab)  # low=0
            self._bc = r.bc if j else Scalar_(bc=bc)  # low=0
            if (self.ab + self.bc) <= 0:
                raise ValueError(_negative_)
            a, _,  c = self._abc3 = t
            if not (a >= c and isfinite(self._a2b2)
                           and isfinite(self._a2c2)):
                raise ValueError(_not_(_finite_))
        except (TypeError, ValueError) as x:
            raise TriaxialError(radius_triaxial=r, ab=ab, bc=bc, cause=x)
        if name:
            self.name = name

    @Property_RO
    def ab(self):
        '''Get relative magnitude C{a - b} (C{meter}, same units as B{C{a}}).
        '''
        return self._ab

    @Property_RO
    def _a2b2(self):
        '''(INTERNAL) Get C{a**2 - b**2} == ab * (a + b).
        '''
        a, b, _ = self._abc3
        return self.ab * (a + b)

    @Property_RO
    def _a2c2(self):
        '''(INTERNAL) Get C{a**2 - c**2} == a2b2 + b2c2.
        '''
        return self._a2b2 + self._b2c2

    @Property_RO
    def bc(self):
        '''Get relative magnitude C{b - c} (C{meter}, same units as B{C{a}}).
        '''
        return self._bc

    @Property_RO
    def _b2c2(self):
        '''(INTERNAL) Get C{b**2 - c**2} == bc * (b + c).
        '''
        _, b, c = self._abc3
        return self.bc * (b + c)

    @Property_RO
    def radius(self):
        '''Get radius (C{meter}, conventionally).
        '''
        return self.a


class TriaxialError(_ValueError):
    '''Raised for L{Triaxial} issues.
    '''
    pass  # ...


class Triaxials(_NamedEnum):
    '''(INTERNAL) L{Triaxial} registry, I{must} be a sub-class
       to accommodate the L{_LazyNamedEnumItem} properties.
    '''
    def _Lazy(self, *abc, **name):
        '''(INTERNAL) Instantiate the C{Triaxial}.
        '''
        a, b, c = map(km2m, abc)
        return Triaxial(a, b, c, **name)

Triaxials = Triaxials(Triaxial, Triaxial_)  # PYCHOK singleton
'''Some pre-defined L{Triaxial}s, all I{lazily} instantiated.'''
# <https://ArxIV.org/pdf/1909.06452.pdf> Table 1 Semi-axes in Km
# <https://www.JPS.NASA.gov/education/images/pdf/ss-moons.pdf>
# <https://link.Springer.com/article/10.1007/s00190-022-01650-9>
_abc84_35 = (_EWGS84.a + 35), (_EWGS84.a - 35), _EWGS84.b
Triaxials._assert(                 # a (Km)       b (Km)     c (Km)     planet
    Amalthea  = _lazy('Amalthea',  125.0,        73.0,     _64_0),    # Jupiter
    Ariel     = _lazy('Ariel',     581.1,       577.9,     577.7),    # Uranus
    Earth     = _lazy('Earth',    6378.173435, 6378.1039, 6356.7544),
    Enceladus = _lazy('Enceladus', 256.6,       251.4,     248.3),    # Saturn
    Europa    = _lazy('Europa',   1564.13,     1561.23,   1560.93),   # Jupiter
    Io        = _lazy('Io',       1829.4,      1819.3,    1815.7),    # Jupiter
    Mars      = _lazy('Mars',     3394.6,      3393.3,    3376.3),
    Mimas     = _lazy('Mimas',     207.4,       196.8,     190.6),    # Saturn
    Miranda   = _lazy('Miranda',   240.4,       234.2,     232.9),    # Uranus
    Moon      = _lazy('Moon',     1735.55,     1735.324,  1734.898),  # Earth
    Tethys    = _lazy('Tethys',    535.6,       528.2,     525.8),    # Saturn
    WGS84_35  = _lazy('WGS84_35', *map(m2km, _abc84_35)))
del _abc84_35, _EWGS84


def _getitems(items, *indices):
    '''(INTERNAL) Get the C{items} at the given I{indices}.

       @return: C{Type(items[i] for i in indices)} with
                C{Type = type(items)}, any C{type} having
                the special method C{__getitem__}.
    '''
    return type(items)(map(items.__getitem__, indices))


def _hartzell3(pov, los, Tun):  # in .Ellipsoid.hartzell4, .formy.hartzell
    '''(INTERNAL) Hartzell's "Satellite Line-of-Sight Intersection ...",
       formula from a Point-Of-View to an I{un-/ordered} Triaxial.
    '''
    def _toUvwV3d(los, pov):
        try:  # pov must be LatLon or Cartesian if los is a Los
            v = los.toUvw(pov)
        except (AttributeError, TypeError):
            v = _otherV3d(los=los)
        return v

    p3 = _otherV3d(pov=pov.toCartesian() if isLatLon(pov) else pov)
    if los is True:  # normal
        a, b, c, d, i = _plumbTo5(p3.x, p3.y, p3.z, Tun)
        return type(p3)(a, b, c), d, i

    u3 = p3.negate() if los is False or los is None else _toUvwV3d(los, pov)

    a, b, c, T = Tun._ordered4

    a2     =  a**2  # largest, factored out
    b2, p2 = (b**2, T._1e2ab) if b != a else (a2, _1_0)
    c2, q2 = (c**2, T._1e2ac) if c != a else (a2, _1_0)

    p3 = T._order3d(p3)
    u3 = T._order3d(u3).unit()  # unit vector, opposing signs

    x2, y2, z2 = p3.x2y2z2  # p3.times_(p3).xyz3
    ux, vy, wz = u3.times_(p3).xyz3
    u2, v2, w2 = u3.x2y2z2  # u3.times_(u3).xyz3

    t = (p2 * c2),  c2, b2
    m = fdot(t, u2, v2, w2)  # a2 factored out
    if m < EPS0:  # zero or near-null LOS vector
        raise _ValueError(_near_(_null_))

    r = fsumf_(b2 * w2,      c2 * v2,      -v2 * z2,      vy * wz * 2,
              -w2 * y2,     -u2 * y2 * q2, -u2 * z2 * p2, ux * wz * 2 * p2,
              -w2 * x2 * p2, b2 * u2 * q2, -v2 * x2 * q2, ux * vy * 2 * q2)
    if r > 0:  # a2 factored out
        r = sqrt(r) * b * c  # == a * a * b * c / a2
    elif r < 0:  # LOS pointing away from or missing the triaxial
        raise _ValueError(_opposite_ if max(ux, vy, wz) > 0 else _outside_)

    d = Fdot(t, ux, vy, wz).fadd_(r).fover(m)  # -r for antipode, a2 factored out
    if d > 0:  # POV inside or LOS outside or missing the triaxial
        s = fsumf_(_N_1_0, _over(x2, a2), _over(y2, b2), _over(z2, c2))  # like _sideOf
        raise _ValueError(_outside_ if s > 0 else _inside_)
    elif fsum1f_(x2, y2, z2, -d*d) < 0:  # d past triaxial's center
        raise _ValueError(_too_(_distant_))

    v = p3.minus(u3.times(d))  # cartesian type(pov) or Vector3d
    h = p3.minus(v).length  # distance to pov == -d
    return T._order3d(v, reverse=True), h, None


def hartzell4(pov, los=False, tri_biax=_WGS84, **name):
    '''Compute the intersection of a tri-/biaxial ellipsoid and a Line-Of-Sight from
       a Point-Of-View outside.

       @arg pov: Point-Of-View outside the tri-/biaxial (C{Cartesian}, L{Ecef9Tuple}
                 C{LatLon} or L{Vector3d}).
       @kwarg los: Line-Of-Sight, I{direction} to the tri-/biaxial (L{Los}, L{Vector3d}),
                   or C{True} for the I{normal, perpendicular, plumb} to the surface of
                   the tri-/biaxial or C{False} or C{None} to point to its center.
       @kwarg tri_biax: A triaxial (L{Triaxial}, L{Triaxial_}, L{JacobiConformal} or
                        L{JacobiConformalSpherical}) or biaxial ellipsoid (L{Datum},
                        L{Ellipsoid}, L{Ellipsoid2}, L{a_f2Tuple}) or spherical earth
                        radius (C{scalar}, conventionally in C{meter}).
       @kwarg name: Optional name (C{str}), overriding C{B{name}="hartzell4"}.

       @return: L{Vector4Tuple}C{(x, y, z, h)} on the tri-/biaxial's surface, with C{h}
                the distance from B{C{pov}} to C{(x, y, z)} I{along the} B{C{los}}, all
                in C{meter}, conventionally.

       @raise TriaxialError: Invalid B{C{pov}} or B{C{pov}} inside the tri-/biaxial or
                             invalid B{C{los}} or B{C{los}} points outside or away from
                             the tri-/biaxial.

       @raise TypeError: Invalid B{C{tri_biax}}, C{ellipsoid} or C{datum}.

       @see: Class L{pygeodesy3.Los}, functions L{pygeodesy.tyr3d} and L{pygeodesy.hartzell}
             and U{lookAtSpheroid<https://PyPI.org/project/pymap3d>} and U{"Satellite
             Line-of-Sight Intersection with Earth"<https://StephenHartzell.Medium.com/
             satellite-line-of-sight-intersection-with-earth-d786b4a6a9b6>}.
    '''
    if isinstance(tri_biax, Triaxial_):
        T = tri_biax
    else:
        D = tri_biax if isinstance(tri_biax, Datum) else \
                  _spherical_datum(tri_biax, name__=hartzell4)  # _DUNDER_nameof
        T = D.ellipsoid._triaxial
    try:
        v, h, i = _hartzell3(pov, los, T)
    except Exception as x:
        raise TriaxialError(pov=pov, los=los, tri_biax=tri_biax, cause=x)
    n = _name__(name, name__=hartzell4)  # _DUNDER_nameof
    return Vector4Tuple(v.x, v.y, v.z, h, iteration=i, name=n)


def _hypot2_1(x, y, z=0):
    '''(INTERNAL) Compute M{x**2 + y**2 + z**2 - 1} with C{max(fabs(x), fabs(y),
       fabs(z))} rarely greater than 1.0.
    '''
    return fsumf_(_N_1_0, x*x, y*y, z*z)


def _otherV3d_(x_xyz, y, z, **name):
    '''(INTERNAL) Get a Vector3d from C{x_xyz}, C{y} and C{z}.
    '''
    return Vector3d(x_xyz, y, z, **name) if isscalar(x_xyz) else \
          _otherV3d(x_xyz=x_xyz, **name)


def _over0(p, q):
    '''(INTERNAL) Return C{p / q} or C{0}.
    '''
    return (p / q) if q > fabs(p) else _0_0


def _over02(p, q):
    '''(INTERNAL) Return C{(p / q)**2} or C{0}.
    '''
    return (p / q)**2 if p and q else _0_0


def _plumbTo3(px, py, E, eps=EPS):  # in .ellipsoids.Ellipsoid.height4
    '''(INTERNAL) Nearest point on a 2-D ellipse in 1st quadrant.
    '''
    a, b = E.a, E.b
    if min(px, py, a, b) < EPS0:
        raise _AssertionError(px=px, py=py, a=a, b=b, E=E)

    a2 = a - b * E.b_a
    b2 = b - a * E.a_b
    tx = ty = _SQRT2_2
    for i in range(16):  # max 5
        ex = tx**3 * a2
        ey = ty**3 * b2

        qx = px - ex
        qy = py - ey
        q  = hypot(qx, qy)
        if q < EPS0:  # PYCHOK no cover
            break
        r = hypot(ex - tx * a,
                  ey - ty * b) / q

        sx, tx = tx, min(_1_0, max(0, (ex + qx * r) / a))
        sy, ty = ty, min(_1_0, max(0, (ey + qy * r) / b))
        t = hypot(ty, tx)
        if t < EPS0:  # PYCHOK no cover
            break
        tx = tx / t  # /= chokes PyChecker
        ty = ty / t
        if fabs(sx - tx) < eps and \
           fabs(sy - ty) < eps:
            break

    tx *= a / px
    ty *= b / py
    return tx, ty, i  # x and y as fractions


def _plumbTo4(x, y, a, b, eps=EPS):
    '''(INTERNAL) Nearest point on and distance to a 2-D ellipse, I{unordered}.

       @see: Function C{_plumbTo3} and I{Eberly}'s U{Distance from a Point to ... an Ellipse ...
             <https://www.GeometricTools.com/Documentation/DistancePointEllipseEllipsoid.pdf>}.
    '''
    if b > a:
        b, a, d, i = _plumbTo4(y, x, b, a, eps=eps)
        return a, b, d, i

    if not (b > 0 and isfinite(a)):
        raise _ValueError(a=a, b=b)

    i = None
    if y:
        if x:
            u =  fabs(x / a)
            w =  fabs(y / b)
            g = _hypot2_1(u, w)
            if fabs(g) > EPS02:
                r = (b / a)**2
                t, i = _rootNd(_1_0 / r, 0, u, 0, w, g)  # eps
                a = _over(x, t * r + _1_0)
                b = _over(y, t     + _1_0)
                d =  hypot(x - a, y - b)
            else:  # on the ellipse
                a, b, d = x, y, _0_0
        else:  # x == 0
            if y < 0:
                b = -b
            a = x  # signed-0
            d = fabs(y - b)

    elif x:  # y == 0
        d, r = None, _over0(a * x, (a + b) * (a - b))
        if r:
            a *= r
            r = _1_0 - r**2
            if r > EPS02:
                b *= sqrt(r)
                d  = hypot(x - a, y - b)
        elif x < 0:
            a = -a
        if d is None:
            b = y  # signed-0
            d = fabs(x - a)

    else:  # x == y == 0
        a = x  # signed-0
        d = b

    return a, b, d, i


def _plumbTo5(x, y, z, Tun, eps=EPS):  # in .testTriaxials
    '''(INTERNAL) Nearest point on and distance to an I{un-/ordered} triaxial.

       @see: I{Eberly}'s U{Distance from a Point to ... an Ellipsoid ...<https://
             www.GeometricTools.com/Documentation/DistancePointEllipseEllipsoid.pdf>}.
    '''
    a, b, c, T = Tun._ordered4
    if Tun is not T:  # T is ordered, Tun isn't
        t = T._order3(x, y, z) + (T,)
        a, b, c, d, i = _plumbTo5(*t, eps=eps)
        return T._order3(a, b, c, reverse=True) + (d, i)

    if not (c > 0 and isfinite(a)):
        raise _ValueError(a=a, b=b, c=c)

    if eps > 0:
        val = max(eps * 1e8, EPS)
    else:  # no validation
        val, eps = 0, max(EPS0, -eps)

    i = None
    if z:
        if y:
            if x:
                u =  fabs(x / a)
                v =  fabs(y / b)
                w =  fabs(z / c)
                g = _hypot2_1(u, v, w)
                if fabs(g) > EPS02:
                    r = T._1e2ac  # (c / a)**2
                    s = T._1e2bc  # (c / b)**2
                    t, i = _rootNd(_1_0 / r, _1_0 / s, u, v, w, g)  # eps
                    a = _over(x, t * r + _1_0)
                    b = _over(y, t * s + _1_0)
                    c = _over(z, t     + _1_0)
                    d =  hypot_(x - a, y - b, z - c)
                else:  # on the ellipsoid
                    a, b, c, d = x, y, z, _0_0
            else:  # x == 0
                a          =  x  # 0
                b, c, d, i = _plumbTo4(y, z, b, c, eps=eps)
        elif x:  # y == 0
            b          =  y  # 0
            a, c, d, i = _plumbTo4(x, z, a, c, eps=eps)
        else:  # x == y == 0
            if z < 0:
                c = -c
            a, b, d = x, y, fabs(z - c)

    else:  # z == 0
        u = _over0(a * x, T._a2c2)  # (a + c) * (a - c)
        v = _over0(b * y, T._b2c2)  # (b + c) * (b - c)
        s = _hypot2_1(u, v)
        if u and v and s < 0:
            a *= u
            b *= v
            c *= sqrt(-s)
            d  = hypot_(x - a, y - b, c)
        else:
            c          =  z  # signed-0
            a, b, d, i = _plumbTo4(x, y, a, b, eps=eps)

    if val > 0:
        _validate(a, b, c, d, T, x, y, z, val)
    return a, b, c, d, i


def _rootNd(r, s, u, v, w, g, eps=EPS0):
    '''(INTERNAL) Robust 2-D or 3-D root finder: 2-D if C{s == v == 0} else 3-D root.

       @see: I{Eberly}'s U{Robust Root Finders ... and Listing 4<https://
             www.GeometricTools.com/Documentation/DistancePointEllipseEllipsoid.pdf>}.
    '''
    u *=  r
    v *=  s  # 0 for 2-D root
    t0 =  w + _N_1_0
    t1 = _0_0 if g < 0 else (hypot_(u, w, v) + _N_1_0)
    # assert t0 <= t1
    _e = -eps
    for i in range(1, _TRIPS):  # 47..65
        t = (t1 + t0) * _0_5
        e =  t1 - t0
        if _e < e < eps or t in (t0, t1):
            break
        g = fsumf_(_N_1_0,  # ~= _hypot2_1
                   _over02(u, t +  r),
                   _over02(w, t - _N_1_0),
                   _over02(v, t +  s) if v else _0_0)
        if g > 0:
            t0 = t
        elif g < 0:
            t1 = t
        else:
            break
    else:  # PYCHOK no cover
        t = Fmt.no_convergence(e, eps)
        raise _ValueError(t, txt__=_rootNd)
    return t, i


def _sqrt0(x):
    '''(INTERNAL) C{sqrt0} with C{TriaxialError}.
    '''
    return sqrt0(x, Error=TriaxialError)


def _validate(a, b, c, d, T, x, y, z, val):
    '''(INTERNAL) Validate an C{_plumTo5} result.
    '''
    e = T.sideOf(a, b, c, eps=val)
    if e:  # not near the ellipsoid's surface
        raise _ValueError(a=a, b=b, c=c, d=d,
                          sideOf=e, eps=val)
    if d:  # angle of delta and normal vector
        m = Vector3d(x, y, z).minus_(a, b, c)
        if m.euclid > val:
            m = m.unit()
            n = T.normal3d(a, b, c)
            e = n.dot(m)  # n.negate().dot(m)
            if not isnear1(fabs(e), eps1=val):
                raise _ValueError(n=n, m=m,
                                dot=e, eps=val)


if __name__ == '__main__':

    from pygeodesy import printf
    from pygeodesy.interns import _COMMA_, _NL_, _NLATvar_

    t = Triaxial_(6378388.0, 6378318.0, 6356911.9461)
    t = t.height4(3909863.9271, 3909778.123, 3170932.5016)
    printf('# Bektas: %r', t)

    # __doc__ of this file, force all into registery
    t = [NN] + Triaxials.toRepr(all=True, asorted=True).split(_NL_)
    printf(_NLATvar_.join(i.strip(_COMMA_) for i in t))

# **) MIT License
#
# Copyright (C) 2022-2025 -- mrJean1 at Gmail -- All Rights Reserved.
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
