
# -*- coding: utf-8 -*-

u'''#-Point resection functions L{ciassini}, L{collins}, L{pierot} and L{tienstra}.
'''

from pygeodesy.basics import isnear0, map1
from pygeodesy.errors import _and, _or, _ValueError
from pygeodesy.fmath import fdot, fidw, fmean, fsum_, fsum1, fsum1_
from pygeodesy.formy import triAngle
from pygeodesy.interns import EPS, EPS0, PI, PI2, _a_, \
                             _A_, _B_, _C_, _coincident_, _colinear_, \
                             _invalid_, _negative_, _SPACE_, \
                             _1_0, _N_1_0, _360_0
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import Fmt, _NamedTuple, _Pass
# from pygeodesy.streprs import Fmt  # from .named
from pygeodesy.units import Degrees, Distance
from pygeodesy.utily import sincos2, sincos2_, sincos2d
from pygeodesy.vector3d import _otherV3d, Vector3d  # PYCHOK unused

from math import atan2, degrees, radians, sin

__all__ = _ALL_LAZY.resections
__version__ = '21.10.23'

_b_         = 'b'
_c_         = 'c'
_concyclic_ = 'concyclic'
_pointH_    = 'pointH'
_pointP_    = 'pointP'


class Collins5Tuple(_NamedTuple):
    '''5-Tuple C{(pointP, pointH, a, b, c)} with survey C{pointP}, auxiliary
       C{pointH}, each an instance of B{C{pointA}}'s (sub-)class and triangle
       sides C{a}, C{b} and C{c} in C{meter}, conventionally.
    '''
    _Names_ = (_pointP_, _pointH_, _a_,      _b_,      _c_)
    _Units_ = (_Pass,    _Pass,     Distance, Distance, Distance)


class ResectionError(_ValueError):
    '''Error raised for resection issues.
    '''
    pass


class Tienstra7Tuple(_NamedTuple):
    '''7-Tuple C{(pointP, A, B, C, a, b, c)} with survey C{pointP}, interior
       triangle angles C{A}, C{B} and C{C} in C{degrees} and triangle sides
       C{a}, C{b} and C{c} in C{meter}, conventionally.
    '''
    _Names_ = (_pointP_, _A_,     _B_,     _C_,     _a_,      _b_,      _c_)
    _Units_ = (_Pass,     Degrees, Degrees, Degrees, Distance, Distance, Distance)


def cassini(pointA, pointB, pointC, alpha, beta, useZ=False):
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

       @note: B{C{PointC}} is between B{C{pointA}} and B{C{pointB}}, typically.

       @return: Survey point, an instance of B{C{pointA}}'s (sub-)class.

       @raise ResectionError: Near-coincident, -colinear or -concyclic points
                              or negative or invalid B{C{alpha}} or B{C{beta}}.

       @raise TypeError: Invalid B{C{pointA}}, B{C{pointB}} or B{C{pointM}}.

       @see: U{Three Point Resection Problem<https://Dokumen.tips/documents/
             three-point-resection-problem-introduction-kaestner-burkhardt-method.html>}
             and functions L{pygeodesy.collins} and L{pygeodesy.tienstra}.
    '''
    try:
        return _cassini(_otherV3d(useZ=useZ, pointA=pointA),
                        _otherV3d(useZ=useZ, pointB=pointB),
                        _otherV3d(useZ=useZ, pointC=pointC),
                         alpha, beta, useZ=useZ, clas=pointA.classof)
    except (TypeError, ValueError) as x:
        raise ResectionError(pointA=pointA, pointB=pointB, pointC=pointC,
                             alpha=alpha, beta=beta, txt=str(x))


def _cassini(A, B, C, alpha, beta, useZ=False, clas=Vector3d):
    # (INTERNAL) Cassini's 3-point resection, note pointC == -M in reference

    def _H(A, C, sa):
        s, c = sincos2d(sa)
        if isnear0(s):
            raise ValueError(_or(_coincident_, _colinear_))
        t = s, c, c
        x = fdot(t, A.x,  C.y, -A.y) / s
        y = fdot(t, A.y, -C.x,  A.x) / s
        return Vector3d(x, y, 0)

    sa, sb = map1(float, alpha, beta)
    if min(sa, sb) < 0:
        raise ValueError(_negative_)
    if fsum_(_360_0, -sa, -sb) < EPS0:
        raise ValueError(_invalid_)

    h1 = _H(A, C,  sa)
    h2 = _H(B, C, -sb)

    x = h1.x - h2.x
    y = h1.y - h2.y
    if isnear0(x) or isnear0(y):
        raise ValueError(_SPACE_(_concyclic_, (x, y)))

    n = x / y
    m = y / x
    N = n + m
    if isnear0(N):
        raise ValueError(_SPACE_(_concyclic_, (N, n, m)))

    t = n, m, _1_0, _N_1_0
    x = fdot(t,  C.x, h1.x, C.y, h1.y) / N
    y = fdot(t, h1.y,  C.y, C.x, h1.x) / N
    z = _zidw(A, B, C, x, y) if useZ else 0
    return clas(x, y, z, name=cassini.__name__)


def collins(pointA, pointB, pointC, alpha, beta, useZ=False):
    '''3-Point resection using U{Collins<https://Dokumen.tips/documents/
       three-point-resection-problem-introduction-kaestner-burkhardt-method.html>}' method.

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

       @note: B{C{PointC}} is between B{C{pointA}} and B{C{pointB}}, typically.

       @return: L{Collins5Tuple}C{(pointP, pointH, a, b, c)} with survey C{pointP},
                auxiliary C{pointH}, each an instance of B{C{pointA}}'s (sub-)class
                and triangle sides C{a}, C{b} and C{c}.

       @raise ResectionError: Near-coincident, -colinear or -concyclic points
                              or negative or invalid B{C{alpha}} or B{C{beta}}.

       @raise TypeError: Invalid B{C{pointA}}, B{C{pointB}} or B{C{pointM}}.

       @see: U{Collins' methode<https://NL.WikiPedia.org/wiki/Achterwaartse_insnijding>}
             and functions L{pygeodesy.cassini} and L{pygeodesy.tienstra}.
    '''
    try:
        return _collins(_otherV3d(useZ=useZ, pointA=pointA),
                        _otherV3d(useZ=useZ, pointB=pointB),
                        _otherV3d(useZ=useZ, pointC=pointC),
                         alpha, beta, useZ=useZ, clas=pointA.classof)
    except (TypeError, ValueError) as x:
        raise ResectionError(pointA=pointA, pointB=pointB, pointC=pointC,
                             alpha=alpha, beta=beta, txt=str(x))


def _collins(A, B, C, alpha, beta, useZ=False, clas=Vector3d):
    # (INTERNAL) Collins' 3-point resection, note C{pointC} center

    def _azi_len2(A, B, pi2):
        v = B.minus(A)
        r = atan2(v.x, v.y)
        if pi2 and r < 0:
            r += pi2
        return r, v.length

    def _cV3(d, r, A, B, C, useZ, V3, **name):
        s, c = sincos2(r)
        x = A.x + d * s
        y = A.y + d * c
        z = _zidw(A, B, C, x, y) if useZ else 0
        return V3(x, y, z, **name)

    ra, rb = radians(alpha), radians(beta)
    if min(ra, rb) < 0:
        raise ValueError(_negative_)

    sra, srH = sin(ra), sin(ra + rb - PI)  # rH = PI - ((PI - ra) + (PI - rb))
    if isnear0(sra) or isnear0(srH):
        raise ValueError(_or(_coincident_, _colinear_, _concyclic_))

#   za, a = _azi_len2(C, B, PI2)
    zb, b = _azi_len2(C, A, PI2)
    zc, c = _azi_len2(A, B, 0)

#   d = c * sin(PI - rb) / srH  # B.minus(H).length
    d = c * sin(PI - ra) / srH  # A.minus(H).length
    r = zc + PI - rb  # zh = zc + (PI - rb)
    H = _cV3(d, r, A, B, C, useZ, Vector3d)

    zh, _ = _azi_len2(C, H, PI2)

#   d = a * sin(za - zh) / sin(rb)  # B.minus(P).length
    d = b * sin(zb - zh) / sra  # A.minus(P).length
    r = zh - ra  # zb - PI + (PI - ra - (zb - zh))
    P = _cV3(d, r, A, B, C, useZ, clas, name=collins.__name__)

    H = clas(H.x, H.y, H.z, name=collins.__name__)
    a = B.minus(C).length
    return Collins5Tuple(P, H, a, b, c)


def tienstra(pointA, pointB, pointC, alpha, beta=None, gamma=None, useZ=False):
    '''3-Point resection using U{Tienstra<https://WikiPedia.org/wiki/Tienstra_formula>}'s formula.

       @arg pointA: First point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg pointB: Second point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg pointC: Third point (C{Cartesian}, L{Vector3d}, C{Vector3Tuple},
                    C{Vector4Tuple} or C{Vector2Tuple} if C{B{useZ}=False}).
       @arg alpha: Angle subtended by triangle side B{C{pointB}} to B{C{pointC}}
                   (C{degrees}, non-negative).
       @kwarg beta: Angle subtended by triangle side B{C{pointA}} to B{C{pointC}}
                    (C{degrees}, non-negative) or C{None} if C{B{gamma} is not None}.
       @kwarg gamma: Angle subtended by triangle side B{C{pointA}} to B{C{pointB}}
                     (C{degrees}, non-negative) or C{None} if C{B{beta} is not None}.
       @kwarg useZ: If C{True}, use and interpolate the Z component, otherwise
                    force C{z=0} (C{bool}).

       @note: Points B{C{pointA}}, B{C{pointB}} and B{C{pointC}} are ordered
              clockwise.

       @return: L{Tienstra7Tuple}C{(pointP, A, B, C, a, b, c)} with survey C{pointP},
                an instance of B{C{pointA}}'s (sub-)class and triangle angles C{A},
                C{B} and C{C} in C{degrees} and triangle sides C{a}, C{b} and C{c}.

       @raise ResectionError: Near-coincident, -colinear or -concyclic points or sum
                              of B{C{alpha}}, B{C{beta}} and B{C{gamma}} not C{360}
                              or negative B{C{alpha}}, B{C{beta}} or B{C{gamma}}.

       @raise TypeError: Invalid B{C{pointA}}, B{C{pointB}} or B{C{pointC}}.

       @see: U{3-Point Resection Solver<http://MesaMike.org/geocache/GC1B0Q9/tienstra/>},
             U{V. Pierlot, M. Van Droogenbroeck, "A New Three Object Triangulation..."
             <http://Telecom.ULG.ac.BE/publi/publications/pierlot/Pierlot2014ANewThree/>},
             U{18 Triangulation Algorithms...<http://Telecom.ULG.ac.BE/triangulation/>} and
             functions L{pygeodesy.cassini} and L{pygeodesy.collins}.
    '''
    try:
        return _tienstra(_otherV3d(useZ=useZ, pointA=pointA),
                         _otherV3d(useZ=useZ, pointB=pointB),
                         _otherV3d(useZ=useZ, pointC=pointC),
                          alpha, beta, gamma, useZ=useZ, clas=pointA.classof)
    except (TypeError, ValueError) as x:
        raise ResectionError(pointA=pointA, pointB=pointB, pointC=pointC,
                             alpha=alpha, beta=beta, gamma=gamma, txt=str(x))


def _tienstra(A, B, C, alpha, beta, gamma, useZ=False, clas=Vector3d):
    # (INTERNAL) Tienstra's 3-point resection

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
    aA = _deg_ks(triAngle(b, c, a), sa, ks, _A_)
    aB = _deg_ks(triAngle(a, c, b), sb, ks, _B_)
    aC = _deg_ks(triAngle(a, b, c), sc, ks, _C_)

    k = fsum1(ks)
    if isnear0(k):
        raise ValueError(Fmt.EQUAL(K=k))
    x = fdot(ks, A.x, B.x, C.x) / k
    y = fdot(ks, A.y, B.y, C.y) / k
    z = _zidw(A, B, C, x, y) if useZ else 0
    P = clas(x, y, z, name=tienstra.__name__)
    return Tienstra7Tuple(P, aA, aB, aC, a, b, c)


def _zidw(A, B, C, x, y):
    # interpolate z or coplanar with A, B and C?
    t = A.z, B.z, C.z
    v = Vector3d(x, y, fmean(t))
    return fidw(t, (v.minus(A).length, v.minus(B).length, v.minus(C).length))

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
