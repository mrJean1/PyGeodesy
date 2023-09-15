
# -*- coding: utf-8 -*-

u'''3-Point resection functions L{cassini}, L{collins5}, L{pierlot}, L{pierlotx} and
L{tienstra7}, survey functions L{snellius3} and L{wildberger3} and triangle functions
L{triAngle}, L{triAngle4}, L{triSide}, L{triSide2} and L{triSide4}.

@note: Functions L{pierlot} and L{pierlotx} are transcoded to Python with permission from
       U{triangulationPierlot<http://www.Telecom.ULg.ac.BE/triangulation/doc/total_8c.html>} and
       U{Pierlot<http://www.Telecom.ULg.ac.BE/publi/publications/pierlot/Pierlot2014ANewThree>}.
'''
# make sure int/int division yields float quotient
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import map1, _ALL_LAZY
from pygeodesy.constants import EPS, EPS0, EPS02, INT0, NEG0, PI, PI2, PI_2, PI_4, \
                               _0_0, _0_5, _1_0, _N_1_0, _2_0, _N_2_0, _4_0, _16_0, \
                               _180_0, _360_0, isnear0, _over, _umod_360
from pygeodesy.errors import _and, _or, TriangleError, _ValueError, _xkwds
from pygeodesy.fmath import favg, Fdot, fidw, fmean, hypot, hypot2_
from pygeodesy.fsums import Fsum, fsumf_, fsum1, fsum1f_
from pygeodesy.interns import _a_, _A_, _area_, _b_, _B_, _c_, _C_, _coincident_, \
                              _colinear_, _d_, _eps_, _invalid_, _negative_, \
                              _not_, _rIn_, _SPACE_
# from pygeodesy.lazily import _ALL_LAZY  # from .basics
from pygeodesy.named import _NamedTuple, _Pass,  Fmt
# from pygeodesy.streprs import Fmt  # from .named
from pygeodesy.units import Degrees, Distance, Radians
from pygeodesy.utily import acos1, asin1, sincos2, sincos2_, sincos2d, sincos2d_
from pygeodesy.vector3d import _otherV3d, Vector3d

from math import cos, atan2, degrees, fabs, radians, sin, sqrt

__all__ = _ALL_LAZY.resections
__version__ = '23.09.12'

_concyclic_ = 'concyclic'
_PA_        = 'PA'
_PB_        = 'PB'
_PC_        = 'PC'
_pointH_    = 'pointH'
_pointP_    = 'pointP'
_positive_  = 'positive'
_R3__       = 'R3 '
_radA_      = 'radA'
_radB_      = 'radB'
_radC_      = 'radC'


class Collins5Tuple(_NamedTuple):
    '''5-Tuple C{(pointP, pointH, a, b, c)} with survey C{pointP}, auxiliary
       C{pointH}, each an instance of B{C{pointA}}'s (sub-)class and triangle
       sides C{a}, C{b} and C{c} in C{meter}, conventionally.
    '''
    _Names_ = (_pointP_, _pointH_, _a_,      _b_,      _c_)
    _Units_ = (_Pass,    _Pass,     Distance, Distance, Distance)


class ResectionError(_ValueError):
    '''Error raised for issues in L{pygeodesy.resections}.
    '''
    pass


class Survey3Tuple(_NamedTuple):
    '''3-Tuple C{(PA, PB, PC)} with distance from survey point C{P} to each of
       the triangle corners C{A}, C{B} and C{C} in C{meter}, conventionally.
    '''
    _Names_ = (_PA_,     _PB_,     _PC_)
    _Units_ = ( Distance, Distance, Distance)


class Tienstra7Tuple(_NamedTuple):
    '''7-Tuple C{(pointP, A, B, C, a, b, c)} with survey C{pointP}, interior
       triangle angles C{A}, C{B} and C{C} in C{degrees} and triangle sides
       C{a}, C{b} and C{c} in C{meter}, conventionally.
    '''
    _Names_ = (_pointP_, _A_,     _B_,     _C_,     _a_,      _b_,      _c_)
    _Units_ = (_Pass,     Degrees, Degrees, Degrees, Distance, Distance, Distance)


class TriAngle5Tuple(_NamedTuple):
    '''5-Tuple C{(radA, radB, radC, rIn, area)} with the interior angles at
       triangle corners C{A}, C{B} and C{C} in C{radians}, the C{InCircle}
       radius C{rIn} aka C{inradius} in C{meter} and the triangle C{area}
       in C{meter} I{squared}, conventionally.
    '''
    _Names_ = (_radA_,  _radB_,  _radC_,  _rIn_,     _area_)
    _Units_ = ( Radians, Radians, Radians, Distance, _Pass)


class TriSide2Tuple(_NamedTuple):
    '''2-Tuple C{(a, radA)} with triangle side C{a} in C{meter}, conventionally
       and angle C{radA} at the opposite triangle corner in C{radians}.
    '''
    _Names_ = (_a_,      _radA_)
    _Units_ = ( Distance, Radians)


class TriSide4Tuple(_NamedTuple):
    '''4-Tuple C{(a, b, radC, d)} with interior angle C{radC} at triangle corner
       C{C} in C{radians}with the length of triangle sides C{a} and C{b} and
       with triangle height C{d} perpendicular to triangle side C{c}, in the
       same units as triangle sides C{a} and C{b}.
    '''
    _Names_ = (_a_,      _b_,      _radC_,  _d_)
    _Units_ = ( Distance, Distance, Radians, Distance)


def _ABC3(useZ, pointA, pointB, pointC):
    '''(INTERNAL) Helper for L{cassini} and L{tienstra7}.
    '''
    return (_otherV3d(useZ=useZ, pointA=pointA),
            _otherV3d(useZ=useZ, pointB=pointB),
            _otherV3d(useZ=useZ, pointC=pointC))


def _B3(useZ, point1, point2, point3):
    '''(INTERNAL) Helper for L{pierlot} and L{pierlotx}.
    '''
    return (_otherV3d(useZ=useZ, point1=point1),
            _otherV3d(useZ=useZ, point2=point2),
            _otherV3d(useZ=useZ, point3=point3))


def cassini(pointA, pointB, pointC, alpha, beta, useZ=False, Clas=None, **Clas_kwds):
    '''3-Point resection using U{Cassini<https://NL.WikiPedia.org/wiki/Achterwaartse_insnijding>}'s method.

       @arg pointA: First point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg pointB: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg pointC: Center point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg alpha: Angle subtended by triangle side B{C{pointA}} to B{C{pointC}}
                   (C{degrees}, non-negative).
       @arg beta: Angle subtended by triangle side B{C{pointB}} to B{C{pointC}}
                  (C{degrees}, non-negative).
       @kwarg useZ: If C{True}, use and interpolate the Z component, otherwise
                    force C{z=INT0} (C{bool}).
       @kwarg Clas: Optional class to return the survey point or C{None} for
                    B{C{pointA}}'s (sub-)class.
       @kwarg Clas_kwds: Optional, additional keyword argument for the survey
                         point instance.

       @note: Typically, B{C{pointC}} is between B{C{pointA}} and B{C{pointB}}.

       @return: The survey point, an instance of B{C{Clas}} or if C{B{Clas} is
                None} an instance of B{C{pointA}}'s (sub-)class.

       @raise ResectionError: Near-coincident, -colinear or -concyclic points
                              or negative or invalid B{C{alpha}} or B{C{beta}}.

       @raise TypeError: Invalid B{C{pointA}}, B{C{pointB}} or B{C{pointM}}.

       @see: U{Three Point Resection Problem<https://Dokumen.tips/documents/
             three-point-resection-problem-introduction-kaestner-burkhardt-method.html>}
             and functions L{pygeodesy.collins5}, L{pygeodesy.pierlot} and
             L{pygeodesy.tienstra7}.
    '''

    def _H(A, C, sa):
        s, c = sincos2d(sa)
        if isnear0(s):
            raise ValueError(_or(_coincident_, _colinear_))
        t = s, c, c
        x = Fdot(t, A.x,  C.y, -A.y).fover(s)
        y = Fdot(t, A.y, -C.x,  A.x).fover(s)
        return x, y

    A, B, C = _ABC3(useZ, pointA, pointB, pointC)
    try:
        sa, sb = map1(float, alpha, beta)
        if min(sa, sb) < 0:
            raise ValueError(_negative_)
        if fsumf_(_360_0, -sa, -sb) < EPS0:
            raise ValueError()

        x1, y1 = _H(A, C,  sa)
        x2, y2 = _H(B, C, -sb)

        x = x1 - x2
        y = y1 - y2
        if isnear0(x) or isnear0(y):
            raise ValueError(_SPACE_(_concyclic_, (x, y)))

        m = y / x
        n = x / y
        N = n + m
        if isnear0(N):
            raise ValueError(_SPACE_(_concyclic_, (m, n, N)))

        t =  n, m, _1_0, _N_1_0
        x =  Fdot(t, C.x, x1, C.y, y1).fover(N)
        y =  Fdot(t, y1, C.y, C.x, x1).fover(N)
        z = _zidw(x, y, useZ, A, B, C)

        clas = Clas or pointA.classof
        return clas(x, y, z, **_xkwds(Clas_kwds, name=cassini.__name__))

    except (TypeError, ValueError) as x:
        raise ResectionError(pointA=pointA, pointB=pointB, pointC=pointC,
                             alpha=alpha, beta=beta, cause=x)


def collins5(pointA, pointB, pointC, alpha, beta, useZ=False, Clas=None, **Clas_kwds):
    '''3-Point resection using U{Collins<https://Dokumen.tips/documents/
       three-point-resection-problem-introduction-kaestner-burkhardt-method.html>}' method.

       @arg pointA: First point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg pointB: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg pointC: Center point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg alpha: Angle subtended by triangle side C{b} from B{C{pointA}} to
                   B{C{pointC}} (C{degrees}, non-negative).
       @arg beta: Angle subtended by triangle side C{a} from B{C{pointB}} to
                  B{C{pointC}} (C{degrees}, non-negative).
       @kwarg useZ: If C{True}, use and interpolate the Z component, otherwise
                    force C{z=INT0} (C{bool}).
       @kwarg Clas: Optional class to return the survey and auxiliary point
                    or C{None} for B{C{pointA}}'s (sub-)class.
       @kwarg Clas_kwds: Optional, additional keyword argument for the survey
                         and auxiliary point instance.

       @note: Typically, B{C{pointC}} is between B{C{pointA}} and B{C{pointB}}.

       @return: L{Collins5Tuple}C{(pointP, pointH, a, b, c)} with survey C{pointP},
                auxiliary C{pointH}, each an instance of B{C{Clas}} or if C{B{Clas}
                is None} an instance of B{C{pointA}}'s (sub-)class and triangle
                sides C{a}, C{b} and C{c} in C{meter}, conventionally.

       @raise ResectionError: Near-coincident, -colinear or -concyclic points
                              or negative or invalid B{C{alpha}} or B{C{beta}}.

       @raise TypeError: Invalid B{C{pointA}}, B{C{pointB}} or B{C{pointM}}.

       @see: U{Collins' methode<https://NL.WikiPedia.org/wiki/Achterwaartse_insnijding>}
             and functions L{pygeodesy.cassini}, L{pygeodesy.pierlot} and
             L{pygeodesy.tienstra7}.
    '''

    def _azi_len2(A, B, pi2):
        v = B.minus(A)
        r = atan2(v.x, v.y)
        if pi2 and r < 0:
            r += pi2
        return r, v.length

    def _cV3(d, r, A, B, C, useZ, V3, **kwds):
        s, c = sincos2(r)
        x =  A.x + d * s
        y =  A.y + d * c
        z = _zidw(x, y, useZ, A, B, C)
        return V3(x, y, z, **kwds)

    A, B, C = _ABC3(useZ, pointA, pointB, pointC)
    try:
        ra, rb = radians(alpha), radians(beta)
        if min(ra, rb) < 0:
            raise ValueError(_negative_)

        sra, srH = sin(ra), sin(ra + rb - PI)  # rH = PI - ((PI - ra) + (PI - rb))
        if isnear0(sra) or isnear0(srH):
            raise ValueError(_or(_coincident_, _colinear_, _concyclic_))

        clas =  Clas or pointA.classof
        kwds = _xkwds(Clas_kwds, name=collins5.__name__)

#       za, a = _azi_len2(C, B, PI2)
        zb, b = _azi_len2(C, A, PI2)
        zc, c = _azi_len2(A, B, 0)

#       d =  c * sin(PI - rb) / srH  # B.minus(H).length
        d =  c * sin(PI - ra) / srH  # A.minus(H).length
        r =  zc + PI - rb  # zh = zc + (PI - rb)
        H = _cV3(d, r, A, B, C, useZ, Vector3d)

        zh, _ = _azi_len2(C, H, PI2)

#       d =  a * sin(za - zh) / sin(rb)  # B.minus(P).length
        d =  b * sin(zb - zh) / sra  # A.minus(P).length
        r =  zh - ra  # zb - PI + (PI - ra - (zb - zh))
        P = _cV3(d, r, A, B, C, useZ, clas, **kwds)

        H =  clas(H.x, H.y, H.z, **kwds)
        a =  B.minus(C).length

        return Collins5Tuple(P, H, a, b, c, name=collins5.__name__)

    except (TypeError, ValueError) as x:
        raise ResectionError(pointA=pointA, pointB=pointB, pointC=pointC,
                             alpha=alpha, beta=beta, cause=x)


def pierlot(point1, point2, point3, alpha12, alpha23, useZ=False, eps=EPS,
                                                      Clas=None, **Clas_kwds):
    '''3-Point resection using U{Pierlot<http://www.Telecom.ULg.ac.BE/publi/publications/
       pierlot/Pierlot2014ANewThree>}'s method C{ToTal} with I{approximate} limits for
       the (pseudo-)singularities.

       @arg point1: First point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg point2: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg point3: Third point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg alpha12: Angle subtended from B{C{point1}} to B{C{point2}} or
                     B{C{alpha2 - alpha1}} (C{degrees}).
       @arg alpha23: Angle subtended from B{C{point2}} to B{C{point3}} or
                     B{C{alpha3 - alpha2}}(C{degrees}).
       @kwarg useZ: If C{True}, interpolate the survey point's Z component,
                    otherwise use C{z=INT0} (C{bool}).
       @kwarg eps: Tolerance for C{cot} (pseudo-)singularities (C{float}).
       @kwarg Clas: Optional class to return the survey point, if C{None} use
                    B{C{point1}}'s (sub-)class.
       @kwarg Clas_kwds: Optional, additional keyword arguments for the survey
                         point instance.

       @note: Typically, B{C{point1}}, B{C{point2}} and B{C{point3}} are ordered
              by angle, modulo 360, counter-clockwise.

       @return: The survey (or robot) point, an instance of B{C{Clas}} or if
                C{B{Clas} is None} an instance of B{C{point1}}'s (sub-)class.

       @raise ResectionError: Near-coincident, -colinear or -concyclic points
                              or invalid B{C{alpha12}} or B{C{alpha23}} or
                              non-positive B{C{eps}}.

       @raise TypeError: Invalid B{C{point1}}, B{C{point2}} or B{C{point3}}.

       @see: I{Pierlot's} C function U{triangulationPierlot<http://www.Telecom.ULg.ac.BE/
             triangulation/doc/total_8c_source.html>}, U{V. Pierlot, M. Van Droogenbroeck,
             "A New Three Object Triangulation Algorithm for Mobile Robot Positioning"
             <https://ORBi.ULiege.BE/bitstream/2268/157469/1/Pierlot2014ANewThree.pdf>},
             U{Vincent Pierlot, Marc Van Droogenbroeck, "18 Triangulation Algorithms for 2D
             Positioning (also known as the Resection Problem)"<http://www.Telecom.ULg.ac.BE/
             triangulation>} and functions L{pygeodesy.pierlotx}, L{pygeodesy.cassini},
             L{pygeodesy.collins5} and L{pygeodesy.tienstra7}.
    '''

    def _cot(s, c):  # -eps < I{approximate} cotangent < eps
        if eps > 0:
            return c / (min(s, -eps) if s < 0 else max(s, eps))
        raise ValueError(_SPACE_(_eps_, _not_, _positive_))

    B1, B2, B3 = _B3(useZ, point1, point2, point3)
    try:
        x, y, z = _pierlot3(B1, B2, B3, alpha12, alpha23, useZ, _cot)
        clas = Clas or point1.classof
        return clas(x, y, z, **_xkwds(Clas_kwds, name=pierlot.__name__))

    except (TypeError, ValueError) as x:
        raise ResectionError(point1=point1, point2=point2, point3=point3,
                             alpha12=alpha12, alpha23=alpha23, eps=eps, cause=x)


def _pierlot3(B1, B2, B3, a12, a23, useZ, cot):
    '''(INTERNAL) Shared L{pierlot} and L{pierlotx}.
    '''
    x1_, y1_, _ = B1.minus(B2).xyz
    x3_, y3_, _ = B3.minus(B2).xyz

    s12, c12, s23, c23 = sincos2d_(a12, a23)
    # cot31 = (1 - cot12 * cot23) / (cot12 + cot32)
    #       = (1 - c12 / s12 * c23 / s23) / (c12 / s12 + c23 / s23)
    #       = (1 - (c12 * c23) / (s12 * s23)) / (c12 * s23 + s12 * c23) / (s12 * s23)
    #       = (s12 * s23 - c12 * c23) / (c12 * s23 + s12 * c23)
    #       =           c31           /           s31
    cot31 = cot(fsum1f_(c12 * s23,  s12 * c23),  # s31
                fsum1f_(s12 * s23, -c12 * c23))  # c31

    K = Fsum(x3_ * x1_,  cot31 * (y3_ * x1_),
             y3_ * y1_, -cot31 * (x3_ * y1_))
    if K:
        cot12 = cot(s12, c12)
        cot23 = cot(s23, c23)

        # x12 = x1_ + cot12 * y1_
        # y12 = y1_ - cot12 * x1_

        # x23 = x3_ - cot23 * y3_
        # y23 = y3_ + cot23 * x3_

        # x31 = x3_ + x1_ + cot31 * (y3_ - y1_)
        # y31 = y3_ + y1_ - cot31 * (x3_ - x1_)

        # x12 - x23 = x1_ + cot12 * y1_ - x3_ + cot23 * y3_
        X12_23 = Fsum(x1_,  cot12 * y1_, -x3_,  cot23 * y3_)
        # y12 - y23 = y1_ - cot12 * x1_ - y3_ - cot23 * x3_
        Y12_23 = Fsum(y1_, -cot12 * x1_, -y3_, -cot23 * x3_)

        # x31 - x23 = x3_ + x1_ + cot31 * (y3_ - y1_) - x3_ + cot23 * y3_
        #           = x1_ + cot31 * y3_ - cot31 * y1_ + cot23 * y3_
        X31_23 = Fsum(x1_, -cot31 * y1_,  cot31 * y3_,  cot23 * y3_)
        # y31 - y23 = y3_ + y1_ - cot31 * (x3_ - x1_) - y3_ - cot23 * x3_
        #           = y1_ - cot31 * x3_ + cot31 * x1_ - cot23 * x3_
        Y31_23 = Fsum(y1_,  cot31 * x1_, -cot31 * x3_, -cot23 * x3_)

        # d = (x12 - x23) * (y23 - y31) + (x31 - x23) * (y12 - y23)
        #   = (x31 - x23) * (y12 - y23) - (x12 - x23) * (y12 - y23)
        d = float(X31_23 * Y12_23 - X12_23 * Y31_23)
        if isnear0(d):
            raise ValueError(_or(_coincident_, _colinear_, _concyclic_))

        x = (B2.x * d + K * Y12_23).fover(d)
        y = (B2.y * d - K * X12_23).fover(d)
    else:
        x, y, _ = B2.xyz
    return x, y, _zidw(x, y, useZ, B1, B2, B3)


def _pierlotx3(a_z_Bs, useZ, cot, C):
    '''(INTERNAL) Core of L{pierlotx}.
    '''
    (a12, z12, B1), \
    (a23, z23, B2), \
    (a31, z31, B3) = a_z_Bs
    if z12 and not z23:
        C(1)
    elif z23 and not z31:
        C(2)
        a23, B1, B2, B3 = a31, B2, B3, B1
    elif z31 and not z12:
        C(3)
        a23,     B2, B3 = a12,     B3, B2
    else:
        C(4)
        return _pierlot3(B1, B2, B3, a12, a23, useZ, cot)

    x1_, y1_, _ = B1.minus(B3).xyz
    x2_, y2_, _ = B2.minus(B3).xyz

    K = Fsum(_1_0, y1_ * x2_, -x1_ * y2_, _N_1_0)  # 1-primed
    if K:
        cot23 = cot(*sincos2d(a23))

        # x23 = x2_ + cot23 * y2_
        # y23 = y2_ - cot23 * x2_

        # x31 = x1_ + cot23 * y1_
        # y31 = y1_ - cot23 * x1_

        # x31 - x23 = x1_ + cot23 * y1_ - x2_ - cot23 * y2_
        X31_23 = Fsum(x1_,  cot23 * y1_, -x2_, -cot23 * y2_)
        # y31 - y23 = y1_ - cot23 * x1_ - y2_ + cot23 * x2_
        Y31_23 = Fsum(y1_, -cot23 * x1_, -y2_,  cot23 * x2_)

        # d = (x31 - x23) * (x2_ - x1_) - (y31 - y23) * (y1_ - y2_)
        d = float(X31_23 * x2_ - X31_23 * x1_ +
                  Y31_23 * y2_ - Y31_23 * y1_)
        if isnear0(d):
            raise ValueError(_or(_coincident_, _colinear_, _concyclic_))

        x = (B3.x * d - K * Y31_23).fover(d)
        y = (B3.y * d + K * X31_23).fover(d)
    else:
        x, y, _ = B3.xyz
    return x, y, _zidw(x, y, useZ, B1, B2, B3)


def pierlotx(point1, point2, point3, alpha1, alpha2, alpha3, useZ=False,
                                                     Clas=None, **Clas_kwds):
    '''3-Point resection using U{Pierlot<http://www.Telecom.ULg.ac.BE/publi/
       publications/pierlot/Pierlot2014ANewThree>}'s method C{ToTal} with
       I{exact} limits for the (pseudo-)singularities.

       @arg point1: First point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg point2: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg point3: Third point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg alpha1: Angle at B{C{point1}} (C{degrees}, counter-clockwise).
       @arg alpha2: Angle at B{C{point2}} (C{degrees}, counter-clockwise).
       @arg alpha3: Angle at B{C{point3}} (C{degrees}, counter-clockwise).
       @kwarg useZ: If C{True}, interpolate the survey point's Z component,
                    otherwise use C{z=INT0} (C{bool}).
       @kwarg Clas: Optional class to return the survey point, if C{None} use
                    B{C{point1}}'s (sub-)class.
       @kwarg Clas_kwds: Optional, additional keyword arguments for the survey
                         point instance.

       @return: The survey (or robot) point, an instance of B{C{Clas}} or if
                C{B{Clas} is None} an instance of B{C{point1}}'s (sub-)class.

       @raise ResectionError: Near-coincident, -colinear or -concyclic points or
                              invalid B{C{alpha1}}, B{C{alpha2}} or B{C{alpha3}}.

       @raise TypeError: Invalid B{C{point1}}, B{C{point2}} or B{C{point3}}.

       @see: I{Pierlot's} C function U{triangulationPierlot2<http://www.Telecom.ULg.ac.BE/
             triangulation/doc/total_8c_source.html>} and function L{pygeodesy.pierlot}.
    '''

    def _a_z_Bs(Bs, *alphas):
        a3 = map(_umod_360, alphas)  # 0 <= alphas < 360
        (a1, a2, a3), (B1, B2, B3) = zip(*sorted(zip(a3, Bs)))
        for a, d, B in ((a1, a2, B1), (a2, a3, B2), (a3, a1, B3)):
            d -= a  # a12 = a2 - a1, ...
            z  = isnear0(fabs(d) % _180_0)
            yield d, z, B

    def _cot(s, c):  # I{exact} cotangent
        try:
            return (c / s) if c else (NEG0 if s < 0 else _0_0)
        except ZeroDivisionError:
            raise ValueError(_or(_coincident_, _colinear_))

    Bs = _B3(useZ, point1, point2, point3)
    try:
        C = [0]  # pseudo-global, passing the exception Case
        x, y, z = _pierlotx3(_a_z_Bs(Bs, alpha1, alpha2, alpha3),
                              useZ, _cot, C.append)
        clas = Clas or point1.classof
        return clas(x, y, z, **_xkwds(Clas_kwds, name=pierlotx.__name__))

    except (TypeError, ValueError) as x:
        raise ResectionError(point1=point1, point2=point2, point3=point3, C=C.pop(),
                             alpha1=alpha1, alpha2=alpha2, alpha3=alpha3, cause=x)


def snellius3(a, b, degC, alpha, beta):
    '''Snellius' surveying using U{Snellius Pothenot<https://WikiPedia.org/wiki/Snellius–Pothenot_problem>}.

       @arg a: Length of the triangle side between corners C{B} and C{C} and opposite of
               triangle corner C{A} (C{scalar}, non-negative C{meter}, conventionally).
       @arg b: Length of the triangle side between corners C{C} and C{A} and opposite of
               triangle corner C{B} (C{scalar}, non-negative C{meter}, conventionally).
       @arg degC: Angle at triangle corner C{C}, opposite triangle side C{c} (non-negative C{degrees}).
       @arg alpha: Angle subtended by triangle side B{C{b}} (non-negative C{degrees}).
       @arg beta: Angle subtended by triangle side B{C{a}} (non-negative C{degrees}).

       @return: L{Survey3Tuple}C{(PA, PB, PC)} with distance from survey point C{P} to
                each of the triangle corners C{A}, C{B} and C{C}, same units as triangle
                sides B{C{a}}, B{C{b}} and B{C{c}}.

       @raise TriangleError: Invalid B{C{a}}, B{C{b}} or B{C{degC}} or negative B{C{alpha}}
                             or B{C{beta}}.

       @see: Function L{pygeodesy.wildberger3}.
    '''
    try:
        a, b, degC, alpha, beta = t = map1(float, a, b, degC, alpha, beta)
        if min(t) < 0:
            raise ValueError(_negative_)
        ra, rb, rC = map1(radians, alpha, beta, degC)

        r = fsumf_(ra, rb, rC) * _0_5
        k = PI - r
        if min(k, r) < 0:
            raise ValueError(_or(_coincident_, _colinear_))

        sa, sb = sin(ra), sin(rb)
        p = atan2(a * sa, b * sb)
        sp, cp, sr, cr = sincos2_(PI_4 - p, r)
        w = atan2(sp * sr, cp * cr)
        x = k + w
        y = k - w

        s = fabs(sa)
        if fabs(sb) > s:
            pc = fabs(a * sin(y) / sb)
        elif s:
            pc = fabs(b * sin(x) / sa)
        else:
            raise ValueError(_or(_colinear_, _coincident_))

        pa = _triSide(b, pc, fsumf_(PI, -ra, -x))
        pb = _triSide(a, pc, fsumf_(PI, -rb, -y))
        return Survey3Tuple(pa, pb, pc, name=snellius3.__name__)

    except (TypeError, ValueError) as x:
        raise TriangleError(a=a, b=b, degC=degC, alpha=alpha, beta=beta, cause=x)


def tienstra7(pointA, pointB, pointC, alpha, beta=None, gamma=None,
                                             useZ=False, Clas=None, **Clas_kwds):
    '''3-Point resection using U{Tienstra<https://WikiPedia.org/wiki/Tienstra_formula>}'s formula.

       @arg pointA: First point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple}, C{Vector4Tuple} or
                    C{Vector2Tuple} if C{B{useZ}=False}).
       @arg pointB: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple}, C{Vector4Tuple} or
                    C{Vector2Tuple} if C{B{useZ}=False}).
       @arg pointC: Third point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple}, C{Vector4Tuple} or
                    C{Vector2Tuple} if C{B{useZ}=False}).
       @arg alpha: Angle subtended by triangle side C{a} from B{C{pointB}} to B{C{pointC}}
                   (C{degrees}, non-negative).
       @kwarg beta: Angle subtended by triangle side C{b} from B{C{pointA}} to B{C{pointC}}
                    (C{degrees}, non-negative) or C{None} if C{B{gamma} is not None}.
       @kwarg gamma: Angle subtended by triangle side C{c} from B{C{pointA}} to B{C{pointB}}
                     (C{degrees}, non-negative) or C{None} if C{B{beta} is not None}.
       @kwarg useZ: If C{True}, use and interpolate the Z component, otherwise force C{z=INT0}
                    (C{bool}).
       @kwarg Clas: Optional class to return the survey point or C{None} for B{C{pointA}}'s
                    (sub-)class.
       @kwarg Clas_kwds: Optional, additional keyword arguments for the survey point instance.

       @note: Points B{C{pointA}}, B{C{pointB}} and B{C{pointC}} are ordered clockwise.

       @return: L{Tienstra7Tuple}C{(pointP, A, B, C, a, b, c)} with survey C{pointP}, an
                instance of B{C{Clas}} or if C{B{Clas} is None} an instance of B{C{pointA}}'s
                (sub-)class, with triangle angles C{A} at B{C{pointA}}, C{B} at B{C{pointB}}
                and C{C} at B{C{pointC}} in C{degrees} and with triangle sides C{a}, C{b} and
                C{c} in C{meter}, conventionally.

       @raise ResectionError: Near-coincident, -colinear or -concyclic points or sum of
                              B{C{alpha}}, B{C{beta}} and B{C{gamma}} not C{360} or negative
                              B{C{alpha}}, B{C{beta}} or B{C{gamma}}.

       @raise TypeError: Invalid B{C{pointA}}, B{C{pointB}} or B{C{pointC}}.

       @see: U{3-Point Resection Solver<http://MesaMike.org/geocache/GC1B0Q9/tienstra/>},
             U{V. Pierlot, M. Van Droogenbroeck, "A New Three Object Triangulation..."
             <http://www.Telecom.ULg.ac.BE/publi/publications/pierlot/Pierlot2014ANewThree/>},
             U{18 Triangulation Algorithms...<http://www.Telecom.ULg.ac.BE/triangulation/>} and
             functions L{pygeodesy.cassini}, L{pygeodesy.collins5} and L{pygeodesy.pierlot}.
    '''

    def _deg_ks(r, s, ks, N):
        if isnear0(fsumf_(PI, r, -s)):  # r + (PI2 - s) == PI
            raise ValueError(Fmt.PARENSPACED(concyclic=N))
        # k = 1 / (cot(r) - cot(s))
        #   = 1 / (cos(r) / sin(r) - cos(s) / sin(s))
        #   = 1 / (cos(r) * sin(s) - cos(s) * sin(r)) / (sin(r) * sin(s))
        #   = sin(r) * sin(s) / (cos(r) * sin(s) - cos(s) * sin(r))
        sr, cr, ss, cs = sincos2_(r, s)
        c = fsum1f_(cr * ss, -cs * sr)
        if isnear0(c):
            raise ValueError(Fmt.PARENSPACED(cotan=N))
        ks.append(sr * ss / c)
        return Degrees(degrees(r), name=N)  # C degrees

    A, B, C = _ABC3(useZ, pointA, pointB, pointC)
    try:
        sa, sb, sc = map1(radians, alpha, (beta or 0), (gamma or 0))
        if beta is None:
            if gamma is None:
                raise ValueError(_and(Fmt.EQUAL(beta=beta), Fmt.EQUAL(gamma=gamma)))
            sb = fsumf_(PI2, -sa, -sc)
        elif gamma is None:
            sc = fsumf_(PI2, -sa, -sb)
        else:  # subtended angles must add to 360 degrees
            r = fsum1f_(sa, sb, sc)
            if fabs(r - PI2) > EPS:
                raise ValueError(Fmt.EQUAL(sum=degrees(r)))
        if min(sa, sb, sc) < 0:
            raise ValueError(_negative_)

        # triangle sides
        a = B.minus(C).length
        b = A.minus(C).length
        c = A.minus(B).length

        ks = []  # 3 Ks and triangle angles
        dA = _deg_ks(_triAngle(b, c, a), sa, ks, _A_)
        dB = _deg_ks(_triAngle(a, c, b), sb, ks, _B_)
        dC = _deg_ks(_triAngle(a, b, c), sc, ks, _C_)

        k = fsum1(ks, floats=True)
        if isnear0(k):
            raise ValueError(Fmt.EQUAL(K=k))
        x =  Fdot(ks, A.x, B.x, C.x).fover(k)
        y =  Fdot(ks, A.y, B.y, C.y).fover(k)
        z = _zidw(x, y, useZ, A, B, C)

        n    = tienstra7.__name__
        clas = Clas or pointA.classof
        P    = clas(x, y, z, **_xkwds(Clas_kwds, name=n))
        return Tienstra7Tuple(P, dA, dB, dC, a, b, c, name=n)

    except (TypeError, ValueError) as x:
        raise ResectionError(pointA=pointA, pointB=pointB, pointC=pointC,
                             alpha=alpha, beta=beta, gamma=gamma, cause=x)


def triAngle(a, b, c):
    '''Compute one angle of a triangle.

       @arg a: Adjacent triangle side length (C{scalar}, non-negative
               C{meter}, conventionally).
       @arg b: Adjacent triangle side length (C{scalar}, non-negative
               C{meter}, conventionally).
       @arg c: Opposite triangle side length (C{scalar}, non-negative
               C{meter}, conventionally).

       @return: Angle in C{radians} at triangle corner C{C}, opposite
                triangle side B{C{c}}.

       @raise TriangleError: Invalid or negative B{C{a}}, B{C{b}} or B{C{c}}.

       @see: Functions L{pygeodesy.triAngle4} and L{pygeodesy.triSide}.
    '''
    try:
        return _triAngle(a, b, c)
    except (TypeError, ValueError) as x:
        raise TriangleError(a=a, b=b, c=c, cause=x)


def _triAngle(a, b, c):
    # (INTERNAL) To allow callers to embellish errors
    a, b, c = map1(float, a, b, c)
    if a < b:
        a, b = b, a
    if b < 0 or c < 0:
        raise ValueError(_negative_)
    if a < EPS0:
        raise ValueError(_coincident_)
    b_a = b / a
    if b_a < EPS0:
        raise ValueError(_coincident_)
    t = fsumf_(_1_0, b_a**2, -(c / a)**2) / (b_a * _2_0)
    return acos1(t)


def triAngle5(a, b, c):
    '''Compute the angles of a triangle.

       @arg a: Length of the triangle side opposite of triangle corner C{A}
               (C{scalar}, non-negative C{meter}, conventionally).
       @arg b: Length of the triangle side opposite of triangle corner C{B}
               (C{scalar}, non-negative C{meter}, conventionally).
       @arg c: Length of the triangle side opposite of triangle corner C{C}
               (C{scalar}, non-negative C{meter}, conventionally).

       @return: L{TriAngle5Tuple}C{(radA, radB, radC, rIn, area)} with angles
                C{radA}, C{radB} and C{radC} at triangle corners C{A}, C{B}
                and C{C}, all in C{radians}, the C{InCircle} radius C{rIn}
                aka C{inradius}, same units as triangle sides B{C{a}},
                B{C{b}} and B{C{c}} and the triangle C{area} in those same
                units I{squared}.

       @raise TriangleError: Invalid or negative B{C{a}}, B{C{b}} or B{C{c}}.

       @see: Functions L{pygeodesy.triAngle} and L{pygeodesy.triArea}.
    '''
    try:
        x, y, z = map1(float, a, b, c)
        ab = x < y
        if ab:
            x, y = y, x
        bc = y < z
        if bc:
            y, z = z, y

        if z > EPS0:  # z = min(a, b, c)
            s  = fsumf_(z, y, x) * _0_5
            sa, sb, r = (s - x), (s - y), (s - z)
            r *= _over(sa * sb, s)
            if r < EPS02:
                raise ValueError(_coincident_)
            r  = sqrt(r)
            rA = atan2(r, sa) * _2_0
            rB = atan2(r, sb) * _2_0
            rC = fsumf_(PI, -rA, -rB)
            if min(rA, rB, rC) < 0:
                raise ValueError(_colinear_)
            s *= r  # Heron's area
        elif z < 0:
            raise ValueError(_negative_)
        else:  # 0 <= c <= EPS0
            rA = rB = PI_2
            rC = r = s = _0_0

        if bc:
            rB, rC = rC, rB
        if ab:
            rA, rB = rB, rA
        return TriAngle5Tuple(rA, rB, rC, r, s, name=triAngle5.__name__)

    except (TypeError, ValueError) as x:
        raise TriangleError(a=a, b=b, c=c, cause=x)


def triArea(a, b, c):
    '''Compute the area of a triangle using U{Heron's<https://
       WikiPedia.org/wiki/Heron%27s_formula>} C{stable} formula.

       @arg a: Length of the triangle side opposite of triangle corner C{A}
               (C{scalar}, non-negative C{meter}, conventionally).
       @arg b: Length of the triangle side opposite of triangle corner C{B}
               (C{scalar}, non-negative C{meter}, conventionally).
       @arg c: Length of the triangle side opposite of triangle corner C{C}
               (C{scalar}, non-negative C{meter}, conventionally).

       @return: The triangle area (C{float}, conventionally C{meter} or
                same units as B{C{a}}, B{C{b}} and B{C{c}} I{squared}).

       @raise TriangleError: Invalid or negative B{C{a}}, B{C{b}} or B{C{c}}.
    '''
    try:
        r, y, x = sorted(map1(float, a, b, c))
        if r > 0:  # r = min(a, b, c)
            ab = x - y
            bc = y - r
            y += r
            r = (x + y) * (r - ab) * (r + ab) * (x + bc)
            if r:
                r = sqrt(r / _16_0)
        elif r < 0:
            raise ValueError(_negative_)
        return r

    except (TypeError, ValueError) as x:
        raise TriangleError(a=a, b=b, c=c, cause=x)


def triSide(a, b, radC):
    '''Compute one side of a triangle.

       @arg a: Adjacent triangle side length (C{scalar},
               non-negative C{meter}, conventionally).
       @arg b: Adjacent triangle side length (C{scalar},
               non-negative C{meter}, conventionally).
       @arg radC: Angle included by sides B{C{a}} and B{C{b}},
                  opposite triangle side C{c} (C{radians}).

       @return: Length of triangle side C{c}, opposite triangle
                corner C{C} and angle B{C{radC}}, same units as
                B{C{a}} and B{C{b}}.

       @raise TriangleError: Invalid B{C{a}}, B{C{b}} or B{C{radC}}.

       @see: Functions L{pygeodesy.sqrt_a}, L{pygeodesy.triAngle},
             L{pygeodesy.triSide2} and L{pygeodesy.triSide4}.
    '''
    try:
        return _triSide(a, b, radC)
    except (TypeError, ValueError) as x:
        raise TriangleError(a=a, b=b, radC=radC, cause=x)


def _triSide(a, b, radC):
    # (INTERNAL) To allow callers to embellish errors
    a, b, r = t = map1(float, a, b, radC)
    if min(t) < 0:
        raise ValueError(_negative_)

    if a < b:
        a, b = b, a
    if a > EPS0:
        ba = b / a
        c2 = fsumf_(_1_0, ba**2, _N_2_0 * ba * cos(r))
        if c2 > EPS02:
            return a * sqrt(c2)
        elif c2 < 0:
            raise ValueError(_invalid_)
    return hypot(a, b)


def triSide2(b, c, radB):
    '''Compute a side and its opposite angle of a triangle.

       @arg b: Adjacent triangle side length (C{scalar},
               non-negative C{meter}, conventionally).
       @arg c: Adjacent triangle side length (C{scalar},
               non-negative C{meter}, conventionally).
       @arg radB: Angle included by sides B{C{a}} and B{C{c}},
                  opposite triangle side C{b} (C{radians}).

       @return: L{TriSide2Tuple}C{(a, radA)} with triangle angle
                C{radA} in C{radians} and length of the opposite
                triangle side C{a}, same units as B{C{b}} and B{C{c}}.

       @raise TriangleError: Invalid B{C{b}} or B{C{c}} or either
                             B{C{b}} or B{C{radB}} near zero.

       @see: Functions L{pygeodesy.sqrt_a}, L{pygeodesy.triSide}
             and L{pygeodesy.triSide4}.
    '''
    try:
        return _triSide2(b, c, radB)
    except (TypeError, ValueError) as x:
        raise TriangleError(b=b, c=c, radB=radB, cause=x)


def _triSide2(b, c, radB):
    # (INTERNAL) To allow callers to embellish errors
    b, c, rB = map1(float, b, c, radB)
    if min(b, c, rB) < 0:
        raise ValueError(_negative_)
    sB, cB = sincos2(rB)
    if isnear0(sB):
        if not isnear0(b):
            raise ValueError(_invalid_)
        if cB < 0:
            a, rA = (b + c), PI
        else:
            a, rA = fabs(b - c), _0_0
    elif isnear0(b):
        raise ValueError(_invalid_)
    else:
        rA = fsumf_(PI, -rB, -asin1(c * sB / b))
        a = sin(rA) * b / sB
    return TriSide2Tuple(a, rA, name=triSide2.__name__)


def triSide4(radA, radB, c):
    '''Compute two sides and the height of a triangle.

       @arg radA: Angle at triangle corner C{A}, opposite triangle side C{a}
                  (non-negative C{radians}).
       @arg radB: Angle at triangle corner C{B}, opposite triangle side C{b}
                  (non-negative C{radians}).
       @arg c: Length of triangle side between triangle corners C{A} and C{B},
               (C{scalar}, non-negative C{meter}, conventionally).

       @return: L{TriSide4Tuple}C{(a, b, radC, d)} with triangle sides C{a} and
                C{b} and triangle height C{d} perpendicular to triangle side
                B{C{c}}, all in the same units as B{C{c}} and interior angle
                C{radC} in C{radians} at triangle corner C{C}, opposite
                triangle side B{C{c}}.

       @raise TriangleError: Invalid or negative B{C{radA}}, B{C{radB}} or B{C{c}}.

       @see: U{Triangulation, Surveying<https://WikiPedia.org/wiki/Triangulation_(surveying)>}
             and functions L{pygeodesy.sqrt_a}, L{pygeodesy.triSide} and L{pygeodesy.triSide2}.
    '''
    try:
        rA, rB, c = map1(float, radA, radB, c)
        rC = fsumf_(PI, -rA, -rB)
        if min(rC, rA, rB, c) < 0:
            raise ValueError(_negative_)
        sa, ca, sb, cb = sincos2_(rA, rB)
        sc = fsum1f_(sa * cb, sb * ca)
        if sc < EPS0 or min(sa, sb) < 0:
            raise ValueError(_invalid_)
        sc = c / sc
        return TriSide4Tuple((sa * sc), (sb * sc), rC, (sa * sb * sc),
                             name=triSide4.__name__)

    except (TypeError, ValueError) as x:
        raise TriangleError(radA=radA, radB=radB, c=c, cause=x)


def wildberger3(a, b, c, alpha, beta, R3=min):
    '''Snellius' surveying using U{Rational Trigonometry
       <https://WikiPedia.org/wiki/Snellius–Pothenot_problem>}.

       @arg a: Length of the triangle side between corners C{B} and C{C} and opposite of
               triangle corner C{A} (C{scalar}, non-negative C{meter}, conventionally).
       @arg b: Length of the triangle side between corners C{C} and C{A} and opposite of
               triangle corner C{B} (C{scalar}, non-negative C{meter}, conventionally).
       @arg c: Length of the triangle side between corners C{A} and C{B} and opposite of
               triangle corner C{C} (C{scalar}, non-negative C{meter}, conventionally).
       @arg alpha: Angle subtended by triangle side B{C{b}} (C{degrees}, non-negative).
       @arg beta: Angle subtended by triangle side B{C{a}} (C{degrees}, non-negative).
       @kwarg R3: Callable to determine C{R3} from C{(R3 - C)**2 = D}, typically standard
                  Python function C{min} or C{max}, invoked with 2 arguments.

       @return: L{Survey3Tuple}C{(PA, PB, PC)} with distance from survey point C{P} to
                each of the triangle corners C{A}, C{B} and C{C}, same units as B{C{a}},
                B{C{b}} and B{C{c}}.

       @raise TriangleError: Invalid B{C{a}}, B{C{b}} or B{C{c}} or negative B{C{alpha}} or
                             B{C{beta}} or B{C{R3}} not C{callable}.

       @see: U{Wildberger, Norman J.<https://Math.Sc.Chula.ac.TH/cjm/content/
             survey-article-greek-geometry-rational-trigonometry-and-snellius-–-pothenot-surveying>},
             U{Devine Proportions, page 252<http://www.MS.LT/derlius/WildbergerDivineProportions.pdf>}
             and function L{pygeodesy.snellius3}.
    '''
    def _s(x):
        return sin(x)**2

    def _vpa(r1, r3, q2, q3, s3):
        r = r1 * r3 * _4_0
        n = (r - Fsum(r1, r3, -q2).fpow(2)).fover(s3)
        if n < 0 or isnear0(r):
            raise ValueError(_coincident_)
        return sqrt((n / r) * q3) if n else _0_0

    try:
        a, b, c, da, db = t = map1(float, a, b, c, alpha, beta)
        if min(t) < 0:
            raise ValueError(_negative_)

        ra, rb = radians(da), radians(db)
        s1, s2, s3 = s = map1(_s, rb, ra, ra + rb)  # rb, ra!
        if min(s) < EPS02:
            raise ValueError(_or(_coincident_, _colinear_))

        q1, q2, q3 = q = a**2, b**2, c**2
        if min(q) < EPS02:
            raise ValueError(_coincident_)

        r1 = s2 * q3 / s3  # s2!
        r2 = s1 * q3 / s3  # s1!
        Qs = Fsum(*q)  # == hypot2_(a, b, c)
        Ss = Fsum(*s)  # == fsum1(s, floats=True)
        s += Qs * _0_5,  # tuple!
        C0 = Fdot(s, q1, q2, q3, -Ss)
        r3 = C0.fover(-s3)
        d0 = Qs.fpow(2).fsub_(hypot2_(*q) * _2_0).fmul(s1 * s2).fover(s3)
        if d0 > EPS02:  # > c0
            d0 = sqrt(d0)
            if not callable(R3):
                raise ValueError(_R3__ + _not_(callable.__name__))
            r3 = R3(float(C0 + d0), float(C0 - d0))  # XXX min or max
        elif d0 < 0:
            raise ValueError(_negative_)

        pa = _vpa(r1, r3, q2, q3, s3)
        pb = _vpa(r2, r3, q1, q3, s3)
        pc =  favg(_triSide2(b, pa, ra).a,
                   _triSide2(a, pb, rb).a)
        return Survey3Tuple(pa, pb, pc, name=wildberger3.__name__)

    except (TypeError, ValueError) as x:
        raise TriangleError(a=a, b=b, c=c, alpha=alpha, beta=beta, R3=R3, cause=x)


def _zidw(x, y, useZ, *ABC):
    if useZ:  # interpolate z or coplanar with A, B and C?
        t = tuple(_.z for _ in ABC)
        v = Vector3d(x, y, fmean(t))
        z = fidw(t, (v.minus(T).length for T in ABC))
    else:
        z = INT0
    return z

# **) MIT License
#
# Copyright (C) 2016-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
