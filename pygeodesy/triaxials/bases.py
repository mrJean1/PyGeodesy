
# -*- coding: utf-8 -*-

u'''(INTERNAL) Base classes for I{ordered} triaxial ellipsoid classes L{Conformal}, L{Conformal3},
L{Triaxial}, L{Triaxial3} and I{unordered} L{Triaxial_}.

Transcoded to pure Python from I{Karney}'s GeographicLib 2.7 C++ classes U{Ellipsoid3<https://
GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Triaxial_1_1Ellipsoid3.html>},
U{Cartesian3<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Triaxial_1_1Cartesian3.html>} and
U{Conformal3<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Triaxial_1_1Conformal3.html>}.

GeographicLib 2.5.2 C++ class U{JacobiConformal<https://GeographicLib.SourceForge.io/C++/doc/
classGeographicLib_1_1JacobiConformal.html#details>}.

Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2008-2024, 2025) and licensed under the MIT/X11 License.
For more information, see the U{GeographicLib 2.5.2 and 2.7<https://GeographicLib.SourceForge.io/>} documentation.

Enum-like C{Lat-/Longitude Kinds (LLK)}, see I{Karney}'s U{coord<https://GeographicLib.SourceForge.io/
C++/doc/classGeographicLib_1_1Triaxial_1_1Cartesian3.html>}:

@var LLK.CONFORMAL: Jacobi conformal X and Y projection
@var LLK.ELLIPSOIDAL: Ellipsoidal lat-, longitude and heading C{bet}, C{omg}, C{alp} (L{Ang})
@var LLK.GEOCENTRIC: Geocentric lat-, longitude and heading C{phi}", C{lam}" and C{zet} (L{Ang})
@var LLK.GEOCENTRIC_X: Geocentric with pole along major X axis
@var LLK.GEODETIC: Geodetic lat-, longitude and heading C{phi}, C{lam} and C{zet} (L{Ang})
@var LLK.GEODETIC_X: Geodetic with pole along major X axis
@var LLK.GEODETIC_LON0: Geodetic lat-, longitude I{- lon0} and heading C{phi}, C{lam} and C{zet} (L{Ang})
@var LLK.GEOGRAPHIC = LLK.GEODETIC
@var LLK.PARAMETRIC: Parametric lat-, longitude and heading C{phi}', C{lam}' and C{zet} (L{Ang})
@var LLK.PARAMETRIC_X: Parametric with pole along major X axis
@var LLK.PLANETODETIC = LLK.GEODETIC
@var LLK.PLANETOCENTRIC = LLK.GEOCENTRIC
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # noqa: E702 ;

# from pygeodesy.angles import Ang, isAng  # _MODS
from pygeodesy.basics import map1, isscalar
from pygeodesy.constants import EPS, EPS0, EPS02, EPS4, _EPS2e4, INT0, NAN, PI2, PI_3, PI4, \
                               _isfinite, float0_, _0_0, _1_0, _N_1_0,  _4_0  # PYCHOK used!
# from pygeodesy.ellipsoids import Ellipsoid  # _MODS
# from pygeodesy.elliptic import Elliptic  # _MODS
# from pygeodesy.errors import _ValueError, _xkwds  # from .formy
from pygeodesy.fmath import fmean_, hypot, norm2, sqrt0,  fabs, sqrt
from pygeodesy.formy import elliperim,  _ValueError, _xkwds
from pygeodesy.fsums import _Fsumf_, fsumf_, fsum1f_
from pygeodesy.interns import _a_, _b_, _c_, _inside_, _not_, _NOTEQUAL_, _null_, \
                              _outside_, _scale_, _SPACE_, _spherical_, _x_, _y_, _z_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import _NamedEnumItem, _NamedTuple, _Pass
from pygeodesy.namedTuples import Vector4Tuple
from pygeodesy.props import Property_RO, property_doc_, property_RO, property_ROver
from pygeodesy.units import Degrees, Easting, Float, Height, Height_, Meter2, Meter3, \
                            Northing, Radius_, Scalar
from pygeodesy.utily import asin1
from pygeodesy.vector3d import _otherV3d, Vector3d

# from math import fabs, sqrt  # from .fmath

__all__ = _ALL_LAZY.triaxials_bases
__version__ = '25.12.12'

_bet_         = 'bet'  # PYCHOK shared
_llk_         = 'llk'  # PYCHOK shared
_MAXIT        =  33  # 20  # PYCHOK shared
_not_ordered_ = _not_('ordered')
_omg_         = 'omg'  # PYCHOK shared


class Conformal5Tuple(_NamedTuple):  # see .Forward4Tuple
    '''5-Tuple C{(x, y, z, scale, llk)} with the easting C{x} and
       northing C{y} projection, C{scale} or C{NAN} I{but with}
       C{z=INT0} I{and kind} C{llk=LLK.CONFORMAL} I{always}.
    '''
    _Names_ = (_x_,     _y_,       _z_,  _scale_, _llk_)
    _Units_ = ( Easting, Northing, _Pass, Scalar, _Pass)

    def __new__(cls, x, y, z=INT0, scale=NAN, llk=None, **kwds):  # **iteration_name
        args = x, y, (z or INT0), scale, (llk or LLK.CONFORMAL)
        return _NamedTuple.__new__(cls, args, **kwds)


class _LLK(str):
    '''(INTERNAL) Lat-/Longitude Kind.
    '''
    def __init__(self, llk):  # aka C++ alt
        self._X = bool(llk.endswith('_X'))
        str.__init__(llk)


class LLK(object):
    '''Enum-like C{Lat-/Longitude Kinds (LLK)}, see U{coord<https://GeographicLib.
       SourceForge.io/C++/doc/classGeographicLib_1_1Triaxial_1_1Cartesian3.html>}.
    '''
    CONFORMAL        = _LLK('CONFORMAL')

    ELLIPSOIDAL      = _LLK('ELLIPSOIDAL')    # bet, omg, alp
    GEOCENTRIC       = _LLK('GEOCENTRIC')     # phi2p, lam2p, zet
    GEOCENTRIC_X     = _LLK('GEOCENTRIC_X')
    GEODETIC         = _LLK('GEODETIC')       # phi, lam, zet
    GEODETIC_LON0    = _LLK('GEODETIC_LON0')
    GEODETIC_X       = _LLK('GEODETIC_X')
    GEOGRAPHIC       =  GEODETIC
    PARAMETRIC       = _LLK('PARAMETRIC')     # phi1p, lam1p, zet
    PARAMETRIC_X     = _LLK('PARAMETRIC_X')
    PLANETODETIC     =  GEODETIC
    PLANETOCENTRIC   =  GEOCENTRIC

    _CENTRICS = (GEOCENTRIC, GEOCENTRIC_X, PLANETOCENTRIC)
    _DETICS   = (GEODETIC, GEODETIC_X, GEODETIC_LON0, GEOGRAPHIC, PLANETODETIC)
    _METRICS  = (PARAMETRIC, PARAMETRIC_X)
    _NOIDAL   = (None, ELLIPSOIDAL)
#   _XCLUDE   = (CONFORMAL, GEOGRAPHIC, PLANETOCENTRIC, PLANETODETIC)

    def items(self):
        for n, llk in LLK.__class__.__dict__.items():
            if isinstance(llk, _LLK):
                yield n, llk

    def keys(self):
        for n, _ in self.items():
            yield n

    def values(self):
        for _, llk in self.items():
            yield llk

LLK = LLK()  # PYCHOK singleton
# del _LLK


def _HeightINT0(h):
    return h if h is INT0 else Height(h=h)


class _UnOrderedTriaxialBase(_NamedEnumItem):
    '''(INTERNAL) Base class for all I{unordered} triaxial classes.
    '''
    _ijk = _kji = None
    _unordered  = True

    def __init__(self, a_triaxial, b=None, c=None, **name):
        '''New I{unordered} C{Triaxial_}.

           @arg a_triaxial: Large, C{X} semi-axis (C{scalar}, conventionally in
                            C{meter}) or an other L{Triaxial}, L{Triaxial_} or
                            L{TriaxialB} instance.
           @kwarg b: Middle, C{Y} semi-axis (C{meter}, same units as B{C{a}}),
                     required if C{B{a_triaxial} is scalar}, ignored otherwise.
           @kwarg c: Small, C{Z} semi-axis (C{meter}, like B{C{b}}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @raise TriaxialError: Invalid semi-axis or -axes.
        '''
        try:
            try:
                a = a_triaxial
                t = a._abc3
                name = _xkwds(name, name=a.name)
            except AttributeError:
                t = Radius_(a=a), Radius_(b=b), Radius_(c=c)
        except (TypeError, ValueError) as x:
            raise TriaxialError(a=a, b=b, c=c, cause=x)
        if name:
            self.name = name

        a, b, c = self._abc3 = t
        if self._unordered:  # == not isinstance(self, Triaxial)
            s, _, t = sorted(t)
            if not (_isfinite(t) and _isfinite(s) and s > 0):
                raise TriaxialError(a=a, b=b, c=c)  # txt=_invalid_
        elif not (_isfinite(a) and a >= b >= c > 0):  # see TriaxialB
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
    def a2(self):
        '''Get C{a**2}.
        '''
        return self.a**2

    @Property_RO
    def _a2b2(self):
        '''(INTERNAL) Get C{a**2 - b**2} == E_sub_e**2.
        '''
        a, b, _ = self._abc3
        d = a - b
        return (d * (a + b)) if d else _0_0

    @Property_RO
    def _a2_b2(self):
        '''(INTERNAL) Get C{(a / b)**2}.
        '''
        a, b, _ = self._abc3
        return (a / b)**2 if a != b else _1_0

    @Property_RO
    def _a2b2c23(self):
        '''(INTERNAL) Get 3-tuple C{(a**2, b**2, c**2)}.
        '''
        a, b, c = self._abc3
        return self.a2, self.b2, self.c2

    @Property_RO
    def _ab_elliperim(self):
        '''(INTERNAL) Get C{ab} ellipse' perimeter.
        '''
        a, b, _ = self._abc3
        return elliperim(a, b)

    @Property_RO
    def _a2c2(self):
        '''(INTERNAL) Get C{a**2 - c**2} == E_sub_x**2.
        '''
        a, _, c = self._abc3
        d = a - c
        return (d * (a + c)) if d else _0_0

    @Property_RO
    def area(self):
        '''Get the surface area (C{meter} I{squared}).
        '''
        c, b, a = sorted(self._abc3)
        return _OrderedTriaxialBase(a, b, c).area if a > c else \
                Meter2(area=self.a2 * PI4)  # a == c == b

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
    def b2(self):
        '''Get C{b**2}.
        '''
        return self.b**2

    @Property_RO
    def _b2_a2(self):
        '''(INTERNAL) Get C{(b / a)**2}.
        '''
        a, b, _ = self._abc3
        return (b / a)**2 if a != b else _1_0

    @Property_RO
    def _b2c2(self):
        '''(INTERNAL) Get C{b**2 - c**2} == E_sub_y**2.
        '''
        _, b, c = self._abc3
        d = b - c
        return (d * (b + c)) if d else _0_0

    @Property_RO
    def _bc_elliperim(self):
        '''(INTERNAL) Get C{bc} ellipse' perimeter.
        '''
        _, b, c = self._abc3
        return elliperim(b, c)

    @Property_RO
    def c(self):
        '''Get the C{smallest, z} semi-axis (C{meter}, same units as B{C{a}}).
        '''
        _, _, c = self._abc3
        return c

    @Property_RO
    def c2(self):
        '''Get C{c**2}.
        '''
        return self.c**2

    @Property_RO
    def _c2_a2(self):
        '''(INTERNAL) Get C{(c / a)**2}.
        '''
        a, _, c = self._abc3
        return (c / a)**2 if a != c else _1_0

    @Property_RO
    def _c2_b2(self):
        '''(INTERNAL) Get C{(c / b)**2}.
        '''
        _, b, c = self._abc3
        return (c / b)**2 if b != c else _1_0

    @Property_RO
    def e2ab(self):
        '''Get the C{ab} ellipse' I{(1st) eccentricity squared} (C{scalar}), M{1 - (b/a)**2}.
        '''
        return Float(e2ab=(_1_0 - self._b2_a2) or _0_0)

#   _1e2ab = _b2_a2  # == C{1 - e2ab} == C{(b/a)**2}

    @Property_RO
    def e2ac(self):
        '''Get the C{ac} ellipse' I{(1st) eccentricity squared} (C{scalar}), M{1 - (c/a)**2}.
        '''
        return Float(e2ac=(_1_0 - self._c2_a2) or _0_0)

#   _1e2ac = _c2_a2  # == C{1 - e2ac} == C{(c/a)**2}

    @Property_RO
    def e2bc(self):
        '''Get the C{bc} ellipse' I{(1st) eccentricity squared} (C{scalar}), M{1 - (c/b)**2}.
        '''
        return Float(e2bc=(_1_0 - self._c2_b2) or _0_0)

#   _1e2bc = _c2_b2  # == C{1 - e2bc} == C{(c/b)**2}

    @property_ROver
    def _Ellipsoid(self):
        '''(INTERNAL) Get class L{Ellipsoid}, I{once}.
        '''
        return _MODS.ellipsoids.Ellipsoid  # overwrite property_ROver

    @property_ROver
    def _Elliptic(self):
        '''(INTERNAL) Get class L{Elliptic}, I{once}.
        '''
        return _MODS.elliptic.Elliptic  # overwrite property_ROver

    def hartzell4(self, pov, los=False, **name):
        '''Compute the intersection of this triaxial's surface with a Line-Of-Sight
           from a Point-Of-View in space.

           @see: Function L{hartzell4<triaxials.triaxial5.hartzell4>} for further details.
        '''
        return self._triaxials_triaxial5.hartzell4(pov, los=los, tri_biax=self, **name)

    def height4(self, x_xyz, y=None, z=None, normal=True, eps=EPS, **name):
        '''Compute the projection on and the height above or below this triaxial's surface.

           @see: Function L{height4<triaxials.triaxial5.height4>} for further details.
        '''
        m = self._triaxials_triaxial5
        return m.height4(x_xyz, y=y, z=z, tri_biax=self, normal=normal, eps=eps, **name)

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
        '''Get a 3-D vector I{on and perpendicular to} this triaxial's surface.

           @arg x_xyz: X component (C{scalar}) or a cartesian (C{Cartesian},
                       L{Ecef9Tuple}, L{Vector3d}, L{Vector3Tuple} or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} if C{scalar}, ignored
                     otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg length: Optional, signed length in out-/inward direction (C{scalar}).

           @return: A C{Vector3d(x_, y_, z_)} normalized to B{C{length}}, pointing out-
                    or inward for postive respectively negative B{C{length}}.

           @raise TriaxialError: Zero length cartesian or vector.

           @note: Cartesian C{(B{x}, B{y}, B{z})} I{must be on} this triaxial's surface,
                  use method L{Triaxial.sideOf} to validate.

           @see: Methods L{Triaxial.height4} and L{Triaxial.sideOf}.
        '''
        # n = 2 * (x / a2, y / b2, z / c2)
        #  == 2 * (x, y * a2 / b2, z * a2 / c2) / a2  # iff ordered
        #  == 2 * (x, y / _b2_a2, z / _c2_a2) / a2
        #  == unit(x, y / _b2_a2, z / _c2_a2).times(length)
        x, y, z = _otherV3d_(x_xyz, y, z).xyz3
        n = Vector3d(x, y / self._b2_a2,
                        z / self._c2_a2, name__=self.normal3d)
        u = n.length
        if u < EPS0:
            raise TriaxialError(x=x_xyz, y=y, z=z, txt=_null_)
        return n.times(length / u)

    def normal4(self, x_xyz, y=None, z=None, height=0, normal=True):
        '''Compute a cartesian at a height above or below this triaxial's surface.

           @arg x_xyz: X component (C{scalar}) or a cartesian (C{Cartesian}, L{Ecef9Tuple},
                       L{Vector3d}, L{Vector3Tuple} or L{Vector4Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} if C{scalar}, ignored
                     otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg height: The signed height (C{scalar}, same units as the triaxial axes).
           @kwarg normal: If C{True}, the B{C{height}} is I{perpendicular, plumb} to the
                          triaxial's surface, otherwise C{radially} to the center of this
                          triaxial (C{bool}).
           @kwarg name: Optional C{B{name}="normal4"} (C{str}).

           @return: L{Vector4Tuple}C{(x, y, z, h)} with the cartesian coordinates C{x},
                    C{y} and C{z} and C{h} the I{signed, normal distance} to the triaxial's
                    surface in C{meter}, conventionally.  Positive C{h} indicates, the
                    cartesian is outside the triaxial, negative C{h} means inside.

           @raise TriaxialError: Zero length cartesian or vector.

           @note: Cartesian C{(B{x}, B{y}, B{z})} I{must be on} this triaxial's surface,
                  use method L{Triaxial.sideOf} to validate.

           @see: Methods L{Triaxial.normal3d} and L{Triaxial.height4}.
        '''
        v, h = _otherV3d(x_xyz, y, z), Height_(height, low=None)
        if h:
            if v.length < EPS0:
                raise TriaxialError(x=x_xyz, y=y, z=z, txt=_null_)
            if normal:
                n  = self.normal3d(v, length=h)
                h  = n.length
                n += v
            else:
                h = h / v.length
                n = v.times(h + _1_0)
        else:
            n = v
        return Vector4Tuple(n.x, n.y, n.z, h, name__=self.normal4)

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
            i, j, k = range(3)
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
                T = _UnOrderedTriaxialBase(*abc)
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
           @kwarg y: Y component (C{scalar}), required if B{C{x_xyz}} is C{scalar},
                     ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg eps: On-surface tolerance (C{scalar}, distance I{squared}).

           @return: C{INT0} if C{(B{x}, B{y}, B{z})} is near this triaxial's surface
                    within tolerance B{C{eps}}, otherwise the signed, radial distance
                    I{squared} (C{float}), nega-/positive for in- respectively outside
                    this triaxial.

           @see: Methods L{Triaxial.height4} and L{Triaxial.normal3d}.
        '''
        v = _otherV3d_(x_xyz, y, z)
        s =  fsumf_(_N_1_0, *map(_over02, v.xyz3, self._abc3))
        return INT0 if fabs(s) < eps else s

    def _sideOn(self, v, eps=_EPS2e4):
        s = self.sideOf(v.xyz, eps=eps)
        if s:  # PYCHOK no cover
            t = _SPACE_((_inside_ if s < 0 else _outside_), self.toRepr())
            raise TriaxialError(eps=eps, sideOf=s, x=v.x, y=v.y, z=v.z, txt=t)
        return s

    def toEllipsoid(self, **name):
        '''Convert this triaxial to a I{biaxial} L{Ellipsoid}, provided 2 axes match.

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
        return self._Ellipsoid(a, b=b, name=self._name__(name))

    toBiaxial = toEllipsoid

    def toStr(self, prec=9, **name):  # PYCHOK signature
        '''Return this C{Triaxial} as a string.

           @kwarg prec: Precision, number of decimal digits (0..9).
           @kwarg name: Optional name (C{str}), to override or C{None}
                        to exclude this triaxial's name.

           @return: This C{Triaxial}'s attributes (C{str}).
        '''
        T = _UnOrderedTriaxialBase
        C =  self._triaxials_triaxial3.Triaxial3B
        if isinstance(self, C):
            t  =  T.b, C.e2, C.k2, C.kp2
        else:
            t  =  T.a,  # props
            C  =  self._triaxials_triaxial5.ConformalSphere
            t += (C.ab, C.bc) if isinstance(self, C) else (T.b, T.c)
            C  = _Triaxial3Base
            t += (C.k2, C.kp2) if isinstance(self, C) else \
                 (T.e2ab, T.e2bc, T.e2ac)
        for C in (self._triaxials_triaxial5.Conformal,
                  self._triaxials_conformal3.Conformal3):
            if isinstance(self, C):
                t += C.xyQ2,
                break
        t += T.volume, T.area
        return self._instr(area_p=self.area_p(), prec=prec, props=t, **name)

    @property_ROver
    def _triaxials_conformal3(self):
        '''(INTERNAL) Get module L{pygeodesy.triaxials.conformal3}, I{once}.
        '''
        return _MODS.triaxials.conformal3  # overwrite property_ROver

    @property_ROver
    def _triaxials_triaxial3(self):
        '''(INTERNAL) Get module L{pygeodesy.triaxials.triaxial3}, I{once}.
        '''
        return _MODS.triaxials.triaxial3  # overwrite property_ROver

    @property_ROver
    def _triaxials_triaxial5(self):
        '''(INTERNAL) Get module L{pygeodesy.triaxials.triaxial5}, I{once}.
        '''
        return _MODS.triaxials.triaxial5  # overwrite property_ROver

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


class _OrderedTriaxialBase(_UnOrderedTriaxialBase):
    '''(INTERNAL) Base class for all I{ordered} triaxial classes.
    '''
    _unordered = False

    def __init__(self, a_triaxial, b=None, c=None, **name):
        '''New I{ordered} L{Triaxial}, L{Triaxial3}, L{Conformal} or L{Conformal3}.

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
        _UnOrderedTriaxialBase.__init__(self, a_triaxial, b=b, c=c, **name)

    @Property_RO
    def _a2b2_a2c2(self):
        '''@see: Methods C{.forwardBetaOmega} and property C{._k2E_kp2E}.
        '''
        s = self._a2c2
        if s:
            s = self._a2b2 / s
        return s or _0_0

    @Property_RO
    def area(self):
        '''Get the surface area (C{meter} I{squared}).

           @see: U{Surface area<https://WikiPedia.org/wiki/Ellipsoid#Surface_area>}.
        '''
        a, b, c = self._abc3
        if a != b:
            kp2, k2 = self._k2E_kp2E  # swapped!
            aE = self._Elliptic(k2, _0_0, kp2, _1_0)
            c2 = self._c2_a2      # cos(phi)**2 = (c/a)**2
            s  = sqrt(self.e2ac)  # sin(phi)**2 =  1 - c2
            r  = asin1(s)  # phi = atan2(sqrt(c2), s)
            b *= fsum1f_(aE.fE(r) * s, (c / a) * (c / b),
                         aE.fF(r) * c2 / s)
            a  = Meter2(area=a * b * PI2)
        else:  # a == b > c
            a  = self._Ellipsoid(a, b=c).areax
        return a

    @Property_RO
    def _k2E_kp2E(self):
        '''(INTERNAL) Get elliptic C{k2} and C{kp2} for C{._xE}, C{._yE} and C{.area}.
        '''
        # k2  = a2b2 / a2c2 * c2_b2
        # kp2 = b2c2 / a2c2 * a2_b2
        # b2  = b**2
        # xE  = Elliptic(k2,  -a2b2 / b2, kp2, a2_b2)
        # yE  = Elliptic(kp2, +b2c2 / b2, k2,  c2_b2)
        # aE  = Elliptic(kp2,  0,         k2,  1)
        k2  = (self._c2_b2 * self._a2b2_a2c2) or _0_0
        kp2 = (self._a2_b2 * self._b2c2 / self._a2c2) if k2 else _1_0
        return k2, kp2

    def _radialTo3(self, sbeta, cbeta, somega, comega):
        '''(INTERNAL) Convert I{ellipsoidal} lat- C{beta} and longitude
           C{omega} to a cartesian I{on this triaxial's surface}, also
           I{ordered} helper for C{.height4 with normal=False}.
        '''
        sa, ca  =  self._norm2(sbeta,  cbeta)
        sb, cb  =  self._norm2(somega, comega)

        b2_a2   =  self._b2_a2  # ==  (b/a)**2
        c2_a2   = -self._c2_a2  # == -(c/a)**2
        a2c2_a2 =  self.  e2ac  # (a**2 - c**2) / a**2 == 1 - (c/a)**2

        x2 = _Fsumf_(_1_0, -b2_a2 * sa**2, c2_a2 * ca**2).fover(a2c2_a2)
        z2 = _Fsumf_(c2_a2,         sb**2, b2_a2 * cb**2).fover(a2c2_a2)

        x, y, z =  self._abc3
        x *= cb * _sqrt0(x2)
        y *= ca *  sb
        z *= sa * _sqrt0(z2)
        return x, y, z


class _Triaxial3Base(_OrderedTriaxialBase):
    '''(INTERNAL) Base class for I{unordered} triaxialC{3} classes.
    '''
    _e2_k2_kp2 = None
    _Lon0      = None

    @Property_RO
    def e2(self):
        '''Get the I{squared eccentricity} (C{scalar}), M{(a**2 - c**2) / b**2}.
        '''
        if self._e2_k2_kp2:
            e2, _, _ = self._e2_k2_kp2
        else:
            e2 = self._a2c2 / self.b2
        return Float(e2=e2)

    def _init_abc3_e2_k2_kp2(self, b, e2, k2, kp2, **name):
        '''(INTERNAL) C{Triaxial3B.__init__}.
        '''
        if name:
            self.name = name
        s = k2 + kp2
        if s > 0 and s != _1_0:
            k2  = k2 / s  # /= chokes PyChecker
            kp2 = kp2 / s
        if min(e2, k2, kp2) < 0 or not s > 0:
            raise TriaxialError(e2=e2, k2=k2, kp2=kp2)
        if e2:
            a = Radius_(a=_sqrt0(_1_0 + e2 * kp2) * b) if kp2 else b
            c = Radius_(c=_sqrt0(_1_0 - e2 * k2)  * b) if k2  else b
        else:  # spherical
            a = c = b
        if not (_isfinite(b) and a >= b >= c > 0):
            raise TriaxialError(b=b, a=a, c=c, e2=e2,
                                k2=k2, kp2=kp2, txt=_not_ordered_)
        self._abc3 = a, b, c
        self._e2_k2_kp2 = e2, k2, kp2

    @property_RO
    def isBiaxial(self):
        '''Is this triaxial I{biaxial} (C{bool}), C{a} == C{b} or C{b} == C{c}?
        '''
        return self.isOblate or self.isProlate

    @property_RO
    def isOblate(self):
        '''Is this triaxial I{oblate} (C{bool}), C{a} == C{b}?
        '''
        return bool(self.kp2 == 0)

    @property_RO
    def isProlate(self):
        '''Is this triaxial I{prolate} (C{bool}), C{b} == C{c}?
        '''
        return bool(self.k2 == 0)

    @Property_RO
    def _k_kp(self):
        '''(INTERNAL) Get the oblate C{k} and prolate C{kp} parameters.
        '''
        return map1(_sqrt0, *self._k2_kp2)

    @Property_RO
    def k2(self):
        '''(INTERNAL) Get the oblate C{k2} parameter I{squared}.
        '''
        k2, _ = self._k2_kp2
        return k2

    @Property_RO
    def _k2_kp2(self):
        '''(INTERNAL) Get the oblate C{k2} and prolate C{kp2} parameters I{squared}.
        '''
        if self._e2_k2_kp2:
            _, k2, kp2 = self._e2_k2_kp2
        else:
            s   =  self._a2c2
            k2  = (self._b2c2 / s) if s else _1_0
            kp2 = (self._a2b2 / s) if s else _0_0
        return k2, kp2

    @Property_RO
    def kp2(self):
        '''(INTERNAL) Get the prolate C{kp2} parameter I{squared}.
        '''
        _, kp2 = self._k2_kp2
        return kp2

    @Property_RO
    def _lcc23(self):
        return self._a2c2, self._b2c2, _0_0

    @property_doc_(" longitude of the I{earth}'s major semi-axis C{a}, (L{Ang}), Karney's C{Triaxial_Earth_lon0}.")
    def Lon0(self):
        if self._Lon0 is None:
            self.Lon0 = -(1493 / 100) if self.name.startswith('WGS84_3') else 0
        return self._Lon0

    @Lon0.setter  # PYCHOK setter!
    def Lon0(self, lon0):
        m = _MODS.angles
        d =  lon0.degrees if m.isAng(lon0) else Degrees(lon0)
        n = _Triaxial3Base.Lon0.name
        self._Lon0 = m.Ang.fromDegrees(d, name=n)

    @Property_RO
    def _xE(self):
        '''(INTERNAL) Get the x-elliptic function.
        '''
        return self._xyE(self.e2, self.k2, self.kp2)

    def _xyE(self, e2, k2, kp2):
        '''(INTERNAL) Helper for C{._xE} and C{._yE}.
        '''
        if e2:
            a2   = -kp2 * e2
            ap2  = _1_0 - a2
            kp2 *= _1_0 - k2 * e2
            k2  *=  ap2
        else:
            a2, ap2 = _0_0, _1_0
        return self._Elliptic(kp2, a2, k2, ap2)

    @Property_RO
    def _yE(self):
        '''(INTERNAL) Get the y-elliptic function.
        '''
        return self._xyE(-self.e2, self.kp2, self.k2)


class TriaxialError(_ValueError):
    '''Raised for any cartesian or conformal triaxial issues.
    '''
    pass  # ...


def _getitems(items, *indices):
    '''(INTERNAL) Get the C{items} at the given I{indices}.

       @return: C{Type(items[i] for i in indices)} with
                C{Type = type(items)}, any C{type} having
                the special method C{__getitem__}.
    '''
    return type(items)(map(items.__getitem__, indices))


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


def _sqrt0(x):
    '''(INTERNAL) C{sqrt0} with C{TriaxialError}.
    '''
    return sqrt0(x, Error=TriaxialError)


__all__ += _ALL_DOCS(_OrderedTriaxialBase, _Triaxial3Base, _UnOrderedTriaxialBase)

# **) MIT License
#
# Copyright (C) 2025-2026 -- mrJean1 at Gmail -- All Rights Reserved.
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
