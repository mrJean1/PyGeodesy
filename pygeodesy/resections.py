
# -*- coding: utf-8 -*-

u'''3-Point resection functions L{cassini}, L{collins}, L{pierlot} and L{tienstra},
survey functions L{snellius3} and L{wildberger3} and triangle functions L{triAngle},
L{triAngle4}, L{triSide}, L{triSide2} and L{triSide4}.

@note: Function L{pierlot} transcoded with permission from U{triangulationPierlot
       <http://www.Telecom.ULg.ac.Be/triangulation/doc/total_8c.html>} and U{Pierlot
       <http://www.Telecom.ULg.ac.Be/publi/publications/pierlot/Pierlot2014ANewThree>}.
'''

from pygeodesy.basics import isnear0, map1
from pygeodesy.errors import _and, _or, ResectionError, TriangleError, _xkwds
from pygeodesy.fmath import favg, fdot, fidw, fmean, fsum, fsum_, fsum1, fsum1_, \
                            hypot, hypot2_
from pygeodesy.interns import EPS, EPS0, EPS02, PI, PI2, PI_2, PI_4, _a_, _A_, \
                             _b_, _B_, _c_, _C_, _coincident_, _colinear_, _d_, \
                             _invalid_, _negative_, _not_, _rIn_, _SPACE_, \
                             _0_0, _0_5, _1_0, _N_1_0, _2_0, _N_2_0, _4_0, _360_0
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import Fmt, _NamedTuple, _Pass
# from pygeodesy.streprs import Fmt  # from .named
from pygeodesy.units import Degrees, Distance, Radians
from pygeodesy.utily import acos1, asin1, sincos2, sincos2_, sincos2d, sincos2d_
from pygeodesy.vector3d import _otherV3d, Vector3d

from math import cos, atan2, degrees, radians, sin, sqrt

__all__ = _ALL_LAZY.resections
__version__ = '21.11.30'

_concyclic_ = 'concyclic'
_PA_        = 'PA'
_PB_        = 'PB'
_PC_        = 'PC'
_pointH_    = 'pointH'
_pointP_    = 'pointP'
_R3_        = 'R3'
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


class Survey3Tuple(_NamedTuple):
    '''3-Tuple C{(PA, PB, PC)} with distance from survey point C{P} to each
       of the triangle corners C{A}, C{B} and C{C}.
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


class TriAngle4Tuple(_NamedTuple):
    '''4-Tuple C{(radA, radB, radC, rIn)} with the interior angles at triangle
       corner C{A}, C{B} and C{C} and the C{InCircle} radius C{rIn} aka C{inradius}.
    '''
    _Names_ = (_radA_,  _radB_,  _radC_,  _rIn_)
    _Units_ = ( Radians, Radians, Radians, Distance)


class TriSide2Tuple(_NamedTuple):
    '''2-Tuple C{(a, radA)} with triangle side C{a} (C{meter} conventionally)
       and the opposite triangle angle C{radA} (C{radians}).
    '''
    _Names_ = (_a_,      _radA_)
    _Units_ = ( Distance, Radians)


class TriSide4Tuple(_NamedTuple):
    '''4-Tuple C{(a, b, radC, d)} with the length of triangle sides C{a} and
       C{b}, the interior angle C{radC} at triangle corner C{radC} (C{radians})
       and triangle height C{d}, perpendicular to triangle side C{c}.
    '''
    _Names_ = (_a_,      _b_,      _radC_,  _d_)
    _Units_ = ( Distance, Distance, Radians, Distance)


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
                    force C{z=0} (C{bool}).
       @kwarg Clas: Optional class to return the survey and auxiliary point
                    or C{None} for B{C{pointA}}'s (sub-)class.
       @kwarg Clas_kwds: Optional additional keyword argument for the survey
                         and auxiliary point instance.

       @note: Typically, B{C{pointC}} is between B{C{pointA}} and B{C{pointB}}.

       @return: The survey point, an instance of B{C{Clas}} or if C{B{Clas} is
                None} of B{C{pointA}}'s (sub-)class.

       @raise ResectionError: Near-coincident, -colinear or -concyclic points
                              or negative or invalid B{C{alpha}} or B{C{beta}}.

       @raise TypeError: Invalid B{C{pointA}}, B{C{pointB}} or B{C{pointM}}.

       @see: U{Three Point Resection Problem<https://Dokumen.tips/documents/
             three-point-resection-problem-introduction-kaestner-burkhardt-method.html>}
             and functions L{pygeodesy.collins}, L{pygeodesy.pierlot} and
             L{pygeodesy.tienstra}.
    '''

    def _H(A, C, sa):
        s, c = sincos2d(sa)
        if isnear0(s):
            raise ValueError(_or(_coincident_, _colinear_))
        t = s, c, c
        x = fdot(t, A.x,  C.y, -A.y) / s
        y = fdot(t, A.y, -C.x,  A.x) / s
        return Vector3d(x, y, 0)

    A = _otherV3d(useZ=useZ, pointA=pointA)
    B = _otherV3d(useZ=useZ, pointB=pointB)
    C = _otherV3d(useZ=useZ, pointC=pointC)

    try:
        sa, sb = map1(float, alpha, beta)
        if min(sa, sb) < 0:
            raise ValueError(_negative_)
        if fsum_(_360_0, -sa, -sb) < EPS0:
            raise ValueError(_invalid_)

        H1 = _H(A, C,  sa)
        H2 = _H(B, C, -sb)

        x = H1.x - H2.x
        y = H1.y - H2.y
        if isnear0(x) or isnear0(y):
            raise ValueError(_SPACE_(_concyclic_, (x, y)))

        m = y / x
        n = x / y
        N = n + m
        if isnear0(N):
            raise ValueError(_SPACE_(_concyclic_, (m, n, N)))

        t = n, m, _1_0, _N_1_0
        x = fdot(t,  C.x, H1.x, C.y, H1.y) / N
        y = fdot(t, H1.y,  C.y, C.x, H1.x) / N
        z = _zidw(A, B, C, x, y) if useZ else 0

        clas = Clas or pointA.classof
        return clas(x, y, z, **_xkwds(Clas_kwds, name=cassini.__name__))

    except (TypeError, ValueError) as x:
        raise ResectionError(pointA=pointA, pointB=pointB, pointC=pointC,
                             alpha=alpha, beta=beta, txt=str(x))


def collins(pointA, pointB, pointC, alpha, beta, useZ=False, Clas=None, **Clas_kwds):
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
                    force C{z=0} (C{bool}).
       @kwarg Clas: Optional class to return the survey point or C{None} for
                    B{C{pointA}}'s (sub-)class.
       @kwarg Clas_kwds: Optional additional keyword argument for the survey
                         point instance.

       @note: Typically, B{C{pointC}} is between B{C{pointA}} and B{C{pointB}}.

       @return: L{Collins5Tuple}C{(pointP, pointH, a, b, c)} with survey C{pointP},
                auxiliary C{pointH}, each an instance of B{C{Clas}} or if C{B{Clas}
                is None} of B{C{pointA}}'s (sub-)class and triangle sides C{a},
                C{b} and C{c}.

       @raise ResectionError: Near-coincident, -colinear or -concyclic points
                              or negative or invalid B{C{alpha}} or B{C{beta}}.

       @raise TypeError: Invalid B{C{pointA}}, B{C{pointB}} or B{C{pointM}}.

       @see: U{Collins' methode<https://NL.WikiPedia.org/wiki/Achterwaartse_insnijding>}
             and functions L{pygeodesy.cassini}, L{pygeodesy.pierlot} and
             L{pygeodesy.tienstra}.
    '''

    def _azi_len2(A, B, pi2):
        v = B.minus(A)
        r = atan2(v.x, v.y)
        if pi2 and r < 0:
            r += pi2
        return r, v.length

    def _cV3(d, r, A, B, C, useZ, V3, **kwds):
        s, c = sincos2(r)
        x = A.x + d * s
        y = A.y + d * c
        z = _zidw(A, B, C, x, y) if useZ else 0
        return V3(x, y, z, **kwds)

    A = _otherV3d(useZ=useZ, pointA=pointA)
    B = _otherV3d(useZ=useZ, pointB=pointB)
    C = _otherV3d(useZ=useZ, pointC=pointC)

    try:
        ra, rb = radians(alpha), radians(beta)
        if min(ra, rb) < 0:
            raise ValueError(_negative_)

        sra, srH = sin(ra), sin(ra + rb - PI)  # rH = PI - ((PI - ra) + (PI - rb))
        if isnear0(sra) or isnear0(srH):
            raise ValueError(_or(_coincident_, _colinear_, _concyclic_))

        clas =  Clas or pointA.classof
        kwds = _xkwds(Clas_kwds, name=collins.__name__)

#       za, a = _azi_len2(C, B, PI2)
        zb, b = _azi_len2(C, A, PI2)
        zc, c = _azi_len2(A, B, 0)

#       d = c * sin(PI - rb) / srH  # B.minus(H).length
        d = c * sin(PI - ra) / srH  # A.minus(H).length
        r = zc + PI - rb  # zh = zc + (PI - rb)
        H = _cV3(d, r, A, B, C, useZ, Vector3d)

        zh, _ = _azi_len2(C, H, PI2)

#       d = a * sin(za - zh) / sin(rb)  # B.minus(P).length
        d = b * sin(zb - zh) / sra  # A.minus(P).length
        r = zh - ra  # zb - PI + (PI - ra - (zb - zh))
        P = _cV3(d, r, A, B, C, useZ, clas, **kwds)

        H = clas(H.x, H.y, H.z, **kwds)
        a = B.minus(C).length
        return Collins5Tuple(P, H, a, b, c)

    except (TypeError, ValueError) as x:
        raise ResectionError(pointA=pointA, pointB=pointB, pointC=pointC,
                             alpha=alpha, beta=beta, txt=str(x))


def pierlot(point1, point2, point3, alpha12, alpha23, useZ=False, Clas=None, **Clas_kwds):
    '''3-Point resection using U{Pierlot<http://www.Telecom.ULg.ac.Be/publi/publications/
       pierlot/Pierlot2014ANewThree>}'s method C{ToTal}.

       @arg point1: First point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg point2: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg point3: Third point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg alpha12: Angle subtended from B{C{point1}} to B{C{point2}} (C{degrees}).
       @arg alpha23: Angle subtended from B{C{point2}} to B{C{point3}} (C{degrees}).
       @kwarg useZ: If C{True}, interpolate the Z component, otherwise use C{z=0}
                    (C{bool}).
       @kwarg Clas: Optional class to return the survey point or C{None} for
                    B{C{point1}}'s (sub-)class.
       @kwarg Clas_kwds: Optional additional keyword arguments for the survey
                         point instance.

       @note: Points B{C{point1}}, B{C{point2}} and B{C{point3}} are ordered
              counter-clockwise, typically.

       @return: The survey (or robot) point, an instance of B{C{Clas}} or if
                C{B{Clas} is None} of B{C{point1}}'s (sub-)class.

       @raise ResectionError: Near-coincident, -colinear or -concyclic points
                              or invalid B{C{alpha12}} or B{C{alpha23}}.

       @raise TypeError: Invalid B{C{point1}}, B{C{point2}} or B{C{point3}}.

       @see: U{V. Pierlot, M. Van Droogenbroeck, "A New Three Object Triangulation
             Algorithm for Mobile Robot Positioning"<https://ORBi.ULiege.Be/
             bitstream/2268/157469/1/Pierlot2014ANewThree.pdf>}, U{Vincent Pierlot,
             Marc Van Droogenbroeck, "18 Triangulation Algorithms for 2D Positioning
             (also known as the Resection Problem)"<http://www.Telecom.ULg.ac.Be/
             triangulation>} and functions L{pygeodesy.cassini}, L{pygeodesy.collins}
             and L{pygeodesy.tienstra}.
    '''
    B1 = _otherV3d(useZ=useZ, point1=point1)
    B2 = _otherV3d(useZ=useZ, point2=point2)
    B3 = _otherV3d(useZ=useZ, point3=point3)

    try:  # (INTERNAL) Raises error for (pseudo-)singularities
        s12, c12, s23, c23 = sincos2d_(alpha12, alpha23)
        if isnear0(s12) or isnear0(s23):
            raise ValueError(_or(_coincident_, _colinear_))
        cot12 = c12 / s12
        cot23 = c23 / s23
#       cot31 = (1 - cot12 * cot23) / (cot12 + cot32)
        d = fsum1_(c12 * s23, s12 * c23)
        if isnear0(d):
            raise ValueError(_or(_coincident_, _colinear_))
        cot31 = fsum1_(s12 * s23, -c12 * c23) / d

        x1_, y1_, _ = B1.minus(B2).xyz
        x3_, y3_, _ = B3.minus(B2).xyz

        x23 = x3_ - cot23 * y3_
        y23 = y3_ + cot23 * x3_

        x12_23 = fsum_(x1_,  cot12 * y1_, -x23)
        y12_23 = fsum_(y1_, -cot12 * x1_, -y23)

        x31_23 = fsum_(x1_, -cot31 * y1_, x3_,  cot31 * y3_, -x23)
        y31_23 = fsum_(y1_,  cot31 * x1_, y3_, -cot31 * x3_, -y23)

        d = x31_23 * y12_23 - x12_23 * y31_23
        if isnear0(d):
            raise ValueError(_or(_coincident_, _colinear_, _concyclic_))
        k = fsum_(x3_ * x1_,  cot31 * (y3_ * x1_),
                  y3_ * y1_, -cot31 * (x3_ * y1_)) / d

        x = B2.x + k * y12_23
        y = B2.y - k * x12_23
        z = _zidw(B1, B2, B3, x, y) if useZ else 0

        clas = Clas or point1.classof
        return clas(x, y, z, **_xkwds(Clas_kwds, name=pierlot.__name__))

    except (TypeError, ValueError) as x:
        raise ResectionError(point1=point1, point2=point2, point3=point3,
                             alpha12=alpha12, alpha23=alpha23, txt=str(x))


def snellius3(a, b, degC, alpha, beta):
    '''Snellius' surveying using U{Snellius Pothenot<https://WikiPedia.org/wiki/Snellius–Pothenot_problem>}.

       @arg a: Length of triangle side C{BC}, opposite of triangle corner C{A} (C{scalar},
               non-negative C{meter} conventionally).
       @arg b: Length of triangle side C{AC}, opposite of triangle corner C{B} (C{scalar},
               non-negative C{meter} conventionally).
       @arg degC: Angle at triangle corner C{C}, opposite triangle side C{c} (C{degrees},
                  non-negative).
       @arg alpha: Angle subtended by triangle side B{C{b}} (C{degrees}, non-negative).
       @kwarg beta: Angle subtended by triangle side B{C{a}} (C{degrees}, non-negative).

       @return: L{Survey3Tuple}C{(PA, PB, PC)} with distance from survey point C{P} to
                each of the triangle corners C{A}, C{B} and C{C} (same units as B{C{a}},
                B{C{b}} and B{C{c}}).

       @raise TriangleError: Invalid B{C{a}}, B{C{b}} or B{C{degC}} or negative B{C{alpha}}
                             or B{C{beta}}.

       @see: Function L{pygeodesy.wildberger3}.
    '''
    try:
        a, b, degC, alpha, beta = map1(float, a, b, degC, alpha, beta)
        ra, rb, rC = map1(radians, alpha, beta, degC)
        if min(ra, rb, rC, a, b) < 0:
            raise ValueError(_negative_)

        r = fsum_(ra, rb, rC) * _0_5
        k = PI - r
        if min(k, r) < 0:
            raise ValueError(_or(_coincident_, _colinear_))

        sa, _, sb, _ = sincos2_(ra, rb)
        p = atan2(a * sa, b * sb)
        sp, cp, sr, cr = sincos2_(PI_4 - p, r)
        w = atan2(sp * sr, cp * cr)
        x = k + w
        y = k - w

        pc = abs((a * sin(y) / sb) if abs(sb) > abs(sa) else
                 (b * sin(x) / sa))
        pa = _triSide(b, pc, fsum_(PI, -ra, -x))
        pb = _triSide(a, pc, fsum_(PI, -rb, -y))
        return Survey3Tuple(pa, pb, pc)

    except (TypeError, ValueError) as x:
        raise TriangleError(a=a, b=b, degC=degC, alpha=alpha, beta=beta, txt=str(x))


def tienstra(pointA, pointB, pointC, alpha, beta=None, gamma=None,
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
       @kwarg useZ: If C{True}, use and interpolate the Z component, otherwise force C{z=0}
                    (C{bool}).
       @kwarg Clas: Optional class to return the survey point or C{None} for B{C{pointA}}'s
                    (sub-)class.
       @kwarg Clas_kwds: Optional additional keyword arguments for the survey point instance.

       @note: Points B{C{pointA}}, B{C{pointB}} and B{C{pointC}} are ordered clockwise.

       @return: L{Tienstra7Tuple}C{(pointP, A, B, C, a, b, c)} with survey C{pointP},
                an instance of B{C{Clas}} or if C{B{Clas} is None} of B{C{pointA}}'s
                (sub-)class and triangle angle C{A} at B{C{pointA}}, C{B} at B{C{pointB}}
                and C{C} at B{C{pointC}} in C{degrees} and triangle sides C{a}, C{b} and C{c}.

       @raise ResectionError: Near-coincident, -colinear or -concyclic points or sum of
                              B{C{alpha}}, B{C{beta}} and B{C{gamma}} not C{360} or
                              negative B{C{alpha}}, B{C{beta}} or B{C{gamma}}.

       @raise TypeError: Invalid B{C{pointA}}, B{C{pointB}} or B{C{pointC}}.

       @see: U{3-Point Resection Solver<http://MesaMike.org/geocache/GC1B0Q9/tienstra/>},
             U{V. Pierlot, M. Van Droogenbroeck, "A New Three Object Triangulation..."
             <http://www.Telecom.ULg.ac.Be/publi/publications/pierlot/Pierlot2014ANewThree/>},
             U{18 Triangulation Algorithms...<http://www.Telecom.ULg.ac.Be/triangulation/>} and
             functions L{pygeodesy.cassini}, L{pygeodesy.collins} and L{pygeodesy.pierlot}.
    '''

    def _deg_ks(r, s, ks, N):
        if isnear0(fsum1_(PI, r, -s)):  # r + (PI2 - s) == PI
            raise ValueError(Fmt.PARENSPACED(concyclic=N))
        # k = 1 / (cot(r) - cot(s))
        #   = 1 / (cos(r) / sin(r) - cos(s) / sin(s))
        #   = 1 / (cos(r) * sin(s) - cos(s) * sin(r)) / (sin(r) * sin(s))
        #   = sin(r) * sin(s) / (cos(r) * sin(s) - cos(s) * sin(r))
        sr, cr, ss, cs = sincos2_(r, s)
        c = cr * ss - cs * sr
        if isnear0(c):
            raise ValueError(Fmt.PARENSPACED(cotan=N))
        ks.append(sr * ss / c)
        return Degrees(degrees(r), name=N)  # C degrees

    A = _otherV3d(useZ=useZ, pointA=pointA)
    B = _otherV3d(useZ=useZ, pointB=pointB)
    C = _otherV3d(useZ=useZ, pointC=pointC)

    try:
        sa, sb, sc = map1(radians, alpha, (beta or 0), (gamma or 0))
        if beta is None:
            if gamma is None:
                raise ValueError(_and(Fmt.EQUAL(beta=beta), Fmt.EQUAL(gamma=gamma)))
            sb = fsum1_(PI2, -sa, -sc)
        elif gamma is None:
            sc = fsum1_(PI2, -sa, -sb)
        else:  # subtended angles must add to 360 degrees
            r = fsum1_(sa, sb, sc)
            if abs(r - PI2) > EPS:
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

        k = fsum1(ks)
        if isnear0(k):
            raise ValueError(Fmt.EQUAL(K=k))
        x = fdot(ks, A.x, B.x, C.x) / k
        y = fdot(ks, A.y, B.y, C.y) / k
        z = _zidw(A, B, C, x, y) if useZ else 0

        clas = Clas or pointA.classof
        P = clas(x, y, z, **_xkwds(Clas_kwds, name=tienstra.__name__))
        return Tienstra7Tuple(P, dA, dB, dC, a, b, c)

    except (TypeError, ValueError) as x:
        raise ResectionError(pointA=pointA, pointB=pointB, pointC=pointC,
                             alpha=alpha, beta=beta, gamma=gamma, txt=str(x))


def triAngle(a, b, c):
    '''Compute an interior angle of a triangle.

       @arg a: Adjacent triangle side length (C{scalar}, non-negative
               C{meter} conventionally).
       @arg b: Adjacent triangle side length (C{scalar}, non-negative
               C{meter} conventionally).
       @arg c: Opposite triangle side length (C{scalar}, non-negative
               C{meter} conventionally).

       @return: Angle at triangle corner C{C}, opposite triangle side
                B{C{c}} (C{radians}).

       @raise TriangleError: Invalid B{C{a}}, B{C{b}} or B{C{c}}.

       @see: Function L{pygeodesy.triSide}.
    '''
    try:
        return _triAngle(a, b, c)
    except (TypeError, ValueError) as x:
        raise TriangleError(a=a, b=b, c=c, tx=str(x))


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
    return acos1(fsum_(_1_0, b_a**2, -(c / a)**2) / (b_a * _2_0))


def triAngle4(a, b, c):
    '''Compute the angles of a triangle side.

       @arg a: Length of triangle side C{BC}, opposite of triangle corner C{A} (C{scalar},
               non-negative C{meter} conventionally).
       @arg b: Length of triangle side C{AC}, opposite of triangle corner C{B} (C{scalar},
               non-negative C{meter} conventionally).
       @arg c: Length of triangle side C{AB}, opposite of triangle corner C{C} (C{scalar},
               non-negative C{meter} conventionally).

       @return: L{TriAngle4Tuple}C{(radA, radB, radC, rIn)} with the triangle angles at
                corner C{A}, C{B} and C{C} (each in radians) and the C{InCircle} radius
                C{rIn} aka C{inradius} (same units and triangle sides B{C{a}}, B{C{b}}
                and B{C{c}}).

       @raise TriangleError: Invalid B{C{a}}, B{C{b}} or B{C{b}}.
    '''
    try:
        a, b, c = map1(float, a, b, c)
        ab = a < b
        if ab:
            a, b = b, a
        bc = b < c
        if bc:
            b, c = c, b

        if c > EPS0:  # c = min(a, b, c)
            s = fsum_(a, b, c) * _0_5
            if s < EPS0:
                raise ValueError(_negative_)
            sa, sb, sc = (s - a), (s - b), (s - c)
            r = sa * sb * sc / s
            if r < EPS02:
                raise ValueError(_coincident_)
            r = sqrt(r)
            rA = atan2(r, sa) * _2_0
            rB = atan2(r, sb) * _2_0
            rC = fsum_(PI, -rA, -rB)
            if min(rA, rB, rC) < 0:
                raise ValueError(_colinear_)
        elif c < 0:
            raise ValueError(_negative_)
        else:  # 0 <= c <= EPS0
            rA = rB = PI_2
            rC = r  = _0_0

        if bc:
            rB, rC = rC, rB
        if ab:
            rA, rB = rB, rA
        return TriAngle4Tuple(rA, rB, rC, r)

    except (TypeError, ValueError) as x:
        raise TriangleError(a=a, b=b, c=c, txt=str(x))


def triSide(a, b, radC):
    '''Compute the length of a triangle side.

       @arg a: Adjacent triangle side length (C{scalar},
               non-negative C{meter} conventionally).
       @arg b: Adjacent triangle side length (C{scalar},
               non-negative C{meter} conventionally).
       @arg radC: Angle included by sides B{C{a}} and B{C{b}},
                  opposite triangle side C{c} (C{radians}).

       @return: Length of triangle side C{c}, opposite angle B{C{rC}}
                (same units as B{C{a}} and B{C{b}}).

       @raise TriangleError: Invalid B{C{a}}, B{C{b}} or B{C{radC}}.
    '''
    try:
        return _triSide(a, b, radC)
    except (TypeError, ValueError) as x:
        raise TriangleError(a=a, b=b, radC=radC, txt=str(x))


def _triSide(a, b, radC):
    # (INTERNAL) To allow callers to embellish errors
    a, b, r = map1(float, a, b, radC)
    if min(a, b, r) < 0:
        raise ValueError(_negative_)

    if isnear0(r):
        return abs(a - b)
    ab = _N_2_0 * a * b * cos(r)
    if abs(ab) < EPS02:
        return hypot(a, b)

    c2 = fsum_(a**2, b**2, ab)
    if abs(c2) < EPS02:
        c = _0_0
    elif c2 < 0:
        raise ValueError(_invalid_)
    else:
        c = sqrt(c2)
    return c


def triSide2(b, c, radB):
    '''Compute the length of a triangle side and the angle.

       @arg b: Adjacent triangle side length (C{scalar},
               non-negative C{meter} conventionally).
       @arg c: Adjacent triangle side length (C{scalar},
               non-negative C{meter} conventionally).
       @arg radB: Angle included by sides B{C{a}} and B{C{c}},
                  opposite triangle side C{b} (C{radians}).

       @return: L{TriSide2Tuple}C{(a, radA)} with triangle angle
                C{radA} (radians) and the length of the opposite
                triangle side C{a} (same units as B{C{b}} and B{C{c}}).

       @raise TriangleError: Invalid B{C{b}} or B{C{c}} or either
                             B{C{b}} or B{C{radB}} near zero and
                             not both.
    '''
    try:
        return _triSide2(b, c, radB)
    except (TypeError, ValueError) as x:
        raise TriangleError(b=b, c=c, radB=radB, txt=str(x))


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
            a, rA = abs(b - c), _0_0
    elif isnear0(b):
        raise ValueError(_invalid_)
    else:
        rA = fsum1_(PI, -rB, -asin1(c * sB / b))
        a = sin(rA) * b / sB
    return TriSide2Tuple(a, rA)


def triSide4(radA, radB, c):
    '''Compute the length of two triangle sides and the triangle height.

       @arg radA: Angle at triangle corner C{A}, opposite triangle side C{a}
                  (C{scalar}, non-negative).
       @arg radB: Angle at triangle corner C{B}, opposite triangle side C{b}
                  (C{scalar}, non-negative).
       @arg c: Length of triangle side between corners C{A} and C{B}, (C{scalar},
               non-negative C{meter} conventionally).

       @return: L{TriSide4Tuple}C{(a, b, radC, d)} with triangle sides C{a} and
                C{b} and triangle height C{d} perpendicular to triangle side
                B{C{c}} (all in same units as B{C{c}}) and the interior angle
                at triangle corner C{C} (C{radians}), opposite of triangle
                side B{C{c}}.

       @raise TriangleError: Invalid B{C{a}} or B{C{b}}.

       @see: Function L{pygeodesy.triSide} and U{Triangulation, Surveying
             <https://WikiPedia.org/wiki/Triangulation_(surveying)>}.
    '''
    try:
        rA, rB, c = map1(float, radA, radB, c)
        rC = fsum_(PI, -rA, -rB)
        if min(rC, rA, rB, c) < 0:
            raise ValueError(_negative_)
        sa, ca, sb, cb = sincos2_(rA, rB)
        sc = fsum1_(sa * cb, sb * ca)
        if sc < EPS0 or min(sa, sb) < 0:
            raise ValueError(_invalid_)
        sc = c / sc
        return TriSide4Tuple((sa * sc), (sb * sc), rC, (sa * sb * sc))

    except (TypeError, ValueError) as x:
        raise TriangleError(radA=radA, radB=radB, c=c, txt=str(x))


def wildberger3(a, b, c, alpha, beta, R3=min):
    '''Snellius' surveying using U{Rational Trigonometry<https://WikiPedia.org/wiki/Snellius–Pothenot_problem>}.

       @arg a: Length of triangle side C{BC}, opposite of triangle corner C{A} (C{scalar},
               non-negative C{meter} conventionally).
       @arg b: Length of triangle side C{AC}, opposite of triangle corner C{B} (C{scalar},
               non-negative C{meter} conventionally).
       @arg c: Length of triangle side C{AB}, opposite of triangle corner C{C} (C{scalar},
               non-negative C{meter} conventionally).
       @arg alpha: Angle subtended by triangle side B{C{b}} (C{degrees}, non-negative).
       @arg beta: Angle subtended by triangle side B{C{a}} (C{degrees}, non-negative).
       @kwarg R3: Callable to determine C{R3} from C{(R3 - C)**2 = D}, typically standard
                  function C{min} or C{max}, invoked with 2 arguments.

       @return: L{Survey3Tuple}C{(PA, PB, PC)} with distance from survey point C{P} to
                each of the triangle corners C{A}, C{B} and C{C} (same units as B{C{a}},
                B{C{b}} and B{C{c}}).

       @raise TriangleError: Invalid B{C{a}}, B{C{b}} or B{C{c}} or negative B{C{alpha}} or
                             B{C{beta}} or B{C{R3}} not C{callable}.

       @see: U{Wildberger, Norman J.<https://math.sc.chula.ac.th/cjm/content/
             survey-article-greek-geometry-rational-trigonometry-and-snellius-–-pothenot-surveying>},
             U{Devine Proportions, page 252<http://www.ms.LT/derlius/WildbergerDivineProportions.pdf>}
             and function L{pygeodesy.snellius3}.
    '''
    def _s(x):
        return sin(x)**2

    def _vpa(r1, r3, q2, q3_r3):
        r = r1 * r3 * _4_0
        n = r - fsum_(r1, r3, -q2)**2
        if n < 0 or isnear0(r):
            raise ValueError(_coincident_)
        return sqrt((n / r) * q3_r3)

    try:
        a, b, c, da, db = map1(float, a, b, c, alpha, beta)
        if min(a, b, c, da, db) < 0:
            raise ValueError(_negative_)

        ra, rb = radians(da), radians(db)
        s1, s2, s3 = s = map1(_s, rb, ra, ra + rb)  # rb, ra!
        if min(s) < EPS02:
            raise ValueError(_or(_coincident_, _colinear_))

        q1, q2, q3 = q = a**2, b**2, c**2
        if min(q) < EPS02:
            raise ValueError(_coincident_)
        q3_s3 = q3 / s3

        r1 = s2 * q3_s3  # s2!
        r2 = s1 * q3_s3  # s1!
        qs = fsum(q)  # == hypot2_(a, b, c)
        ss = fsum1(s)
        s += (qs * _0_5),
        c0 = -fdot(s, q1, q2, q3, -ss) / s3
        d0 = (qs**2 - hypot2_(*q) * _2_0) * s1 * s2 / s3
        if d0 > EPS02:  # > c0
            d0 = sqrt(d0)
            if not callable(R3):
                raise ValueError(_SPACE_(_R3_, _not_(callable.__name__)))
            r3 = R3(c0 + d0, c0 - d0)  # XXX min or max
        elif d0 < 0:
            raise ValueError(_negative_)
        else:  # isnear0(sqrt(d0))
            r3 = c0

        pa = _vpa(r1, r3, q2, q3_s3)
        pb = _vpa(r2, r3, q1, q3_s3)
        pc = favg(_triSide2(b, pa, ra).a,
                  _triSide2(a, pb, rb).a)
        return Survey3Tuple(pa, pb, pc)

    except (TypeError, ValueError) as x:
        raise TriangleError(a=a, b=b, c=c, alpha=alpha, beta=beta, R3=R3, txt=str(x))


def _zidw(A, B, C, x, y):
    # interpolate z or coplanar with A, B and C?
    t = A.z, B.z, C.z
    v = Vector3d(x, y, fmean(t))
    return fidw(t, (v.minus(A).length, v.minus(B).length, v.minus(C).length))

# **) MIT License
#
# Copyright (C) 2016-2022 -- mrJean1 at Gmail -- All Rights Reserved.
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
