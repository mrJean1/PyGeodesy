
# -*- coding: utf-8 -*-

u'''Class L{Intersector}, a pure Python version of parts of I{Karney}'s C++ class U{Intersect
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Intersect.html>} to intersect
geodesic lines.

Only C++ member functions C{All}, C{Closest} and C{All} have been transcoded into Python as methods
L{Intersector.All}, L{Intersector.Closest} and L{Intersector.Next} producing 4-item L{XDist}s.

Adjacent methods L{Intersector.All4}, L{Intersector.Closest4}, L{Intersector.Next4} and
L{Intersector.Next4s} return or yield L{Intersector4Tuple}s with the lat-, longitude, azimuth of
each intersection as a C{Position} L{GDict} on each geodesic line.

For more details, see the C++ U{GeographicLib<https://GeographicLib.SourceForge.io/C++/doc/index.html>}
documentation, I{Charles F.F. Karney}'s paper U{Geodesics intersections<https://arxiv.org/abs/2308.00495>}
and I{S. Baselga Moreno & J.C. Martinez-Llario}'s U{Intersection and point-to-line solutions for geodesics
on the ellipsoid<https://riunet.UPV.ES/bitstream/handle/10251/122902/Revised_Manuscript.pdf>}.
'''
# make sure int/int division yields float quotient
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import _copy, _enumereverse, _xinstanceof
from pygeodesy.constants import EPS, INF, INT0, PI, PI2, PI_4, _0_0, \
                               _0_5, _1_0, _1_5, _2_0, _3_0, _90_0
from pygeodesy.ellipsoids import _EWGS84
from pygeodesy.errors import GeodesicError, IntersectionError, \
                            _xgeodesics, _xkwds_get
from pygeodesy.fmath import euclid, favg, fdot
from pygeodesy.fsums import Fsum, fsum1_,  _ceil
from pygeodesy.interns import _A_, _B_, _c_, _too_
from pygeodesy.karney import Caps, _diff182, _sincos2de
from pygeodesy.lazily import _ALL_LAZY  # _ALL_MODS as _MODS
from pygeodesy.named import ADict, _NamedBase, _NamedTuple,  Fmt
from pygeodesy.namedTuples import Int, Meter, _Pass
from pygeodesy.props import Property_RO, property_RO
# from pygeodesy.streprs import Fmt  # from .named
# from pygeodesy.units import Int, Meter  # from .namedTuples
from pygeodesy.utily import sincos2,  atan2, fabs

# from math import atan2, ceil as _ceil, fabs  # .fsums, .utily

__all__ = _ALL_LAZY.geodesici
__version__ = '24.06.24'

_0t     =  0,  # int
_1_1t   = -1, +1
_1_0_1t = -1, 0, +1
_EPS3   =  EPS * _3_0
_EPSr5  =  pow(EPS, 0.2)  # PYCHOK used!
_TRIPS  =  128
_CWGS84 =  max(_EWGS84.a, _EWGS84.b) * PI2  # circumference


def _L1(a, b):
    return fabs(a) + fabs(b)


class XDist(ADict):
    '''4-Item result from L{Intersector.All}, L{Intersector.Closest} and
       L{Intersector.Next} with the intersection offsets C{sA}, C{sB} and
       C{sX0} in C{meter} and the coincidence indicator C{c}, an C{int},
       +1 for parallel, -1 for anti-parallel, 0 otherwise.
    '''
#   Delta = _EWGS84.R2 * PI * _EPSr5  # equality margin, ~15 Km

    def __init__(self, sA=0, sB=0, c=0, sX0=INT0):
        '''New L{XDist}.

          @kwarg sA: Offset on geodesic line A (C{meter}).
          @kwarg sB: Offset on geodesic line B (C{meter}).
          @kwarg c: Coincidence indicator (C{int}, +1 for parallel
                    -1 for anti-parallel, 0 otherwise.
          @kwarg sX0: Offset to C{X0} ({Cmeter}) or L{INT0}.
        '''
        self.set_(sA=sA, sB=sB, c=c, sX0=sX0)

    def __add__(self, other):
        X  = _copy(self)
        X +=  other
        return X

    def __eq__(self, other):
        # _xinstanceof(XDist, other=other)
        return self is other or self.L1(other) < EPS  # <= self.Delta

    def __iadd__(self, other):
        if isinstance(other, tuple) and len(other) == 2:
            a, b = other
        else:
            _xinstanceof(XDist, other=other)
            a = other.sA  # PYCHOK sA
            b = other.sB  # PYCHOK sB
            if other.c:   # PYCHOK c
                self.c = other.c
        self.sA += a  # PYCHOK sA
        self.sB += b  # PYCHOK sB
        return self

    def __le__(self, other):
        # _xinstanceof(XDist, other=other)
        return self == other or self < other

    def __lt__(self, other):
        # _xinstanceof(XDist, other=other)
        return (self.sA < other.sA or (self.sA == other.sA and  # PYCHOK sA
                self.sB < other.sB) and self != other)  # PYCHOK sB

    def __ne__(self, other):
        return not self.__eq__(other)

    def _fixCoincident(self, X, *c0):
        # modify C{X} to the mid-point if X is anti-/parallel
        c = c0[0] if c0 else X.c
        if c:
            X  = _copy(X)
            s  = (self.sA - X.sA  +  # PYCHOK sA
                 (self.sB - X.sB) * c) * _0_5  # PYCHOK sB
            X +=  s, (s * c)
        return X

    def L1(self, other=None):
        '''Return the C{L1} distance.
        '''
        a, b = self.sA, self.sB  # PYCHOK sA, sB
        if other is not None:
            # _xinstanceof(XDist, other=other)
            a -= other.sA
            b -= other.sB
        return _L1(a, b)

    def _nD1(self, n, D1):
        C = XDist
        return self + (C._1A[n] * D1,
                       C._1B[n] * D1)

    _1A = (0, 1, -1, 0,  0)
    _1B = (0, 0,  0, 1, -1)

    def _nD2(self, n, D2):
        C = XDist
        return C(C._2A[n] * D2,
                 C._2B[n] * D2)

    _2A = (-1, -1,  1, 1, -2, 0, 2,  0)
    _2B = (-1,  1, -1, 1,  0, 2, 0, -2)

    def _nmD3(self, n, m, D3):
        for i in range(n, m, 2):
            for j in range(n, m, 2):
                if i or j:
                    yield self + ((i + j) * D3,
                                  (i - j) * D3)

_X000 = XDist()  # PYCHOK origin
_XINF = XDist(INF)


class Intersector(_NamedBase):
    '''Finder of intersections between two goedesic lines, each an instance
       of L{GeodesicLineExact<pygeodesy.geodesicx.GeodesicLineExact>},
       wrapped L{GeodesicLine<pygeodesy.geodesicw.GeodesicLine>} or
       L{GeodesicLineSolve<pygeodesy.geodsolve.GeodesicLineSolve>}.

       @see: I{Karney}'s C++ class U{Intersect<https://GeographicLib.sourceforge.io/
             C++/doc/classGeographicLib_1_1Intersect.html#details>} for more details.
    '''
    # _D1 =  0
    # _D2 =  0
    # _g  =  None
    # _T1 =  0
    # +T5 =  0

    def __init__(self, geodesic, **name):
        '''New L{Intersector}.

           @arg geodesic: The geodesic (L{GeodesicExact<pygeodesy.geodesicx.GeodesicExact>},
                                        wrapped L{Geodesic<pygeodesy.geodesicw.Geodesic>} or
                                        L{GeodesicSolve<pygeodesy.geodsolve.GeodesicSolve>}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @raise GeodesicError: The eccentricity of the B{C{geodesic}}'s ellipsoid is too
                                 large or no initial convergence.

           @see: The B{Note} at I{Karney}'s C++ U{Intersect<https://GeographicLib.sourceforge.io/
                 C++/doc/classGeographicLib_1_1Intersect.html#ae41f54c9a44836f6c8f140f6994930cf>}.
        '''
        _xinstanceof(*_EWGS84._Geodesics, geodesic=geodesic)
        self._g = geodesic
        if name:
            self.name = name
        E = self.ellipsoid

        t1 = E.b * PI  # min distance between intersects
        t2 = self._polarDist2(_90_0)[0] * _2_0  # furthest closest intersect
        t5 = self._Invers12(_90_0) * _2_0  # longest shortest geodesic
        if self.f > 0:
            t3 = self._obliqDist4()[0]
            t4 = t1
        else:  # PYCHOK no cover
            t1, t2, t3 = t2, t1, t5
            t4 = self._polarB3()[0]
        d1 = t2 * _0_5
        d2 = t3 / _1_5
        d3 = t4 - self.Delta
        if not (d1 < d3 and d2 < d3 and d2 < (t1 * _2_0)):
            t = Fmt.PARENSPACED(_too_('eccentric'), E.e)
            raise GeodesicError(ellipsoid=E.toStr(terse=2), txt=t)
        self._D1 = d1  # tile spacing for Closest
        self._D2 = d2  # tile spacing for Next
        self._D3 = d3  # tile spacing for All
        self._T1 = t1  # min distance between intersects
#       self._T5 = t5

    @Property_RO
    def a(self):
        '''Get the I{equatorial} radius, semi-axis (C{meter}).
        '''
        return self.ellipsoid.a

    equatoradius = a  # = Requatorial

    def All(self, glA, glB, X0=_X000, sMax=_CWGS84):
        '''Find all intersection of two geodesic lines up to a limit.

           @arg glA: A geodesic line (L{Line<Intersector.Line>}).
           @arg glB: An other geodesic line (L{Line<Intersector.Line>}).
           @kwarg X0: Optional offsets along the geodesic lines (L{XDist}).
           @kwarg sMax: Upper limit for the distance (C{meter}).

           @return: Yield an L{XDist} for each intersection found.

           @raise GeodesicError: Geodesic line B{C{glA}} or B{C{glB}}
                                 invalid, incompatible or ill-configured.

           @raise IntersectionError: No convergence.
        '''
        self._xLines(glA, glB)
        D, _D = self.Delta, self._RPI
        xMax  = sMax + D
        m     = int(_ceil(xMax / self._D3))  # m x m tiles
        d3    = xMax / m
        T1d3D = self._T1d3Delta(d3)
        _X0fx = X0._fixCoincident

        c0 =  0
        c_ = _List(D)  # closest coincident
        x_ = _List(D)  # intersections found
        s_ =  list(X0._nmD3(1 - m, m, d3 * _0_5))
        # assert len(s_) + 1 == m * m + (m - 1) % 2
        while s_:
            Q, i = self._Basic2(glA, glB, s_.pop(0))
            if Q in x_ or (c0 and _X0fx(Q) in c_):
                continue
            # assert Q.c == c0 or not c0
            a, c0 = len(x_), Q.c
            if c0:  # coincident intersection
                Q, _ = c_.add(_X0fx(Q))
                # elimate all existing intersections
                # on this line (which didn't set c0)
                for j, X in _enumereverse(x_):
                    if _X0fx(X, c0).L1(Q) <= D:  # X' == Q
                        x_.pop(j)

                a, s0 = len(x_), Q.sA
                m_M_M = self._m12_M12_M21_3(glA, s0)
                _cjD  = self._conjDist
                for s in (-_D, _D):
                    s += s0
                    sa = 0
                    while True:
                        sa = _cjD(glA, s + sa, *m_M_M) - s0
                        X  = Q + (sa, (sa * c0))
                        i += 1
                        if x_.add(X, X0.L1(X), i) > xMax:
                            break

            x_.add(Q, X0.L1(Q), i + 1)
            for X in x_[a:]:  # added Xs
                for j, S in _enumereverse(s_):
                    if X.L1(S) < T1d3D:
                        s_.pop(j)

        return x_.sortrim(X0, sMax)  # generator!

    def All4(self, glA, glB, X0=_X000, aMax=0, sMax=_CWGS84):
        '''Find all intersection of two geodesic lines up to a limit.

           @kwarg aMax: Upper limit for the angular distance (C{degrees}).

           @return: Yield an L{Intersector4Tuple}C{(A, B, sAB, c)} for
                    each intersection found.

           @see: Methods L{All} for further details.
        '''
        aA = aB = _0_0
        for X in self.All(glA, glB, X0=X0, sMax=sMax):
            X = self._In4T(glA, glB, X, X)
            yield X
            if aMax:
                aA += X.A.a12
                aB += X.B.a12
                if fabs(aA) > aMax or fabs(aB) > aMax:
                    break

    def _Basic2(self, glA, glB, X, i=0):
        '''(INTERNAL) Get a basic solution.
        '''
        _S, _T = self._Spherical, self._Tol
        for x in range(_TRIPS):
            S = _S(glA, glB, X)
            X += S
            i += 1
            if X.c or S.L1() <= _T:  # or isnan
                return X.set_(iteration=x), i
        raise IntersectionError(Fmt.no_convergence(S.L1(), _T))

    def Closest(self, glA, glB, X0=_X000):
        '''Find the closest intersection of two geodesic lines.

           @arg glA: A geodesic line (L{Line<Intersector.Line>}).
           @arg glB: An other geodesic line (L{Line<Intersector.Line>}).
           @kwarg X0: Optional offsets along the geodesic lines (L{XDist}).

           @return: The intersection (L{XDist}) or C{None} if none found.

           @raise GeodesicError: Geodesic line B{C{glA}} or B{C{glB}}
                                 invalid, incompatible or ill-configured.

           @raise IntersectionError: No convergence.
        '''
        self._xLines(glA, glB)
        i,  d,  Q     = 1, INF, X0  # best so far
        D1, T1, T1D2D = self._D1, self._T1, self._T1D2Delta
        sk, K         = set(), len(XDist._1A)
        for n in (n for n in range(K) if n not in sk):
            X    = X0._nD1(n, D1)
            X, i = self._Basic2(glA, glB, X, i)
            X    = X0._fixCoincident(X)
            if X.L1(Q) > self.Delta:  # X != Q
                d0 = X.L1(X0)
                if d0 < T1:
                    Q, d = X, d0
                    break
                if d0 < d or not n:
                    Q, d = X, d0
                    i += 1
                for m in range(n + 1, K):
                    if m not in sk and \
                       X.L1(X0._nD1(m, D1)) < T1D2D:
                        sk.add(m)
        return None if Q is X0 else Q.set_(sX0=d0, iteration=i)

    def Closest4(self, glA, glB, X0=_X000):
        '''Find the closest intersection of two geodesic lines.

           @return: An L{Intersector4Tuple}C{(A, B, sAB, c)} or
                    C{None} if none found.

           @see: Method L{Closest} for further details.
        '''
        X = self.Closest(glA, glB, X0=X0)
        return X if X is None else self._In4T(glA, glB, X, X)

    def _conjDist(self, gl, s, m12=0, M12=1, M21=1, semi=False):
        # Find semi-/conjugate point relative to s0 which is close to s1.
        # if semi:
        #     solve for M23 = 0 using dM23 / ds3 = - (1 - M23 * M32) / m23
        # else:
        #     solve for m23 = 0 using dm23 / ds3 = M32
        _1,    _T = _1_0,  self._Tol
        _abs, _S2 =  fabs, Fsum(s).fsum2_
        for _ in range(_TRIPS):
            m13, M13, M31 = self._m12_M12_M21_3(gl, s)
            # see "Algorithms for geodesics", eqs. 31, 32, 33.
            m23 = m13 * M12
            M32 = M31 * M12
            if m12:
                m23 -= m12 * M13
                if m13:
                    M32 += (_1 - M13 * M31) * m12 / m13
            if semi:
                M23 = M13 * M21
                # when m12 -> eps, (1 - M12 * M21) -> eps^2, I suppose.
                if m12 and m13:
                    M23 += (_1 - M12 * M21) * m13 / m12
                d =  m23 * M23 / (_1 - M23 * M32)
            else:
                d = -m23 / M32
            s, d = _S2(d)
            if _abs(d) <= _T:
                break
        return s

    def _conjDist5(self, azi):
        gl   = self._Line(azi1=azi)
        s    = self._conjDist(gl, self._RPI)
        X, _ = self._Basic2(gl, gl, XDist(s * _0_5, -s * _1_5))
        return s, (X.L1() - s * _2_0), azi, X.sA, X.sB

    @Property_RO
    def Delta(self):
        '''Get the equality and tiling margin (C{meter}).
        '''
        return self._RPI * _EPSr5  # ~15 Km WGS84

    @Property_RO
    def ellipsoid(self):
        '''Get the C{geodesic}'s ellipsoid (C{Ellipsoid}).
        '''
        return self.geodesic.datum.ellipsoid

    @Property_RO
    def _EPS3R(self):
        return _EPS3 * self.R

    @Property_RO
    def f(self):
        '''Get the I{flattening} (C{scalar}), C{0} for spherical, negative for prolate.
        '''
        return self.ellipsoid.f

    flattening = f

    @Property_RO
    def _faPI_2(self):
        return (self.f + _2_0) * self.a * PI_4

    @property_RO
    def geodesic(self):
        '''Get the C{geodesic} (C{Geodesic...}).
        '''
        return self._g

    @Property_RO
    def _GeodesicLines(self):
        '''(INTERNAL) Get the C{Geodesic...Line} class(es).
        '''
        return type(self._Line()),

    def _In4T(self, glA, glB, S, X):
        # Return an intersection as C{Intersector4Tuple}.
        A = self._Position(glA, S.sA)
        B = self._Position(glB, S.sB)
        s = self._Invers12(A, B)
        r = Intersector4Tuple(A, B, s, X.c, iteration=X.iteration)
        return r

    def _Invers12(self, A, B=None):
        lls = (0, 0, A, 0) if B is None else (A.lat2, A.lon2,
                                              B.lat2, B.lon2)
        return self._g.Inverse(*lls, outmask=Caps.DISTANCE).s12

    def _Inverse(self, A, B):  # Caps.STANDARD
        return self._g.Inverse(A.lat2, A.lon2, B.lat2, B.lon2)

    def Line(self, lat1, lon1, azi_lat2, *lon2, **name):
        '''Return a geodesic line from this C{Intersector}'s geodesic, specified by
           two (goedetic) points or a (goedetic) point and an (initial) azimuth.

           @arg lat1: Latitude of the first point (C{degrees}).
           @arg lon1: Longitude of the first point (C{degrees}).
           @arg azi_lat2: Azimuth at the first point (compass C{degrees}) if no
                          B{C{lon2}} argument is given, otherwise the latitude of
                          the second point (C{degrees}).
           @arg lon2: If given, the longitude of the second point (C{degrees}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: A line (from L{geodesic<Intersector.geodesic>}C{.Line} or
                    C{-.InverseLine}), properly configured for L{Intersector}.
        '''
        args = (lat1, lon1, azi_lat2) + lon2
        gl   = self._g.InverseLine(*args, caps=Caps.LINE_CAPS) if lon2 else \
               self._g.Line(       *args, caps=Caps.LINE_CAPS)
        if name:
            gl.name= name
        return gl

    def _Line(self, lat1=0, lon1=0, azi1=0):
        return self._g.Line(lat1, lon1, azi1, caps=Caps.LINE_CAPS)

    def _m12_M12_M21_3(self, gl, s):
        P = gl.Position(s, outmask=Caps._REDUCEDLENGTH_GEODESICSCALE)
        return P.m12, P.M12, P.M21

    def Next(self, glA, glB, eps=_EPS3):
        '''Yield the next intersection of two I{intersecting} geodesic lines.

           @arg glA: A geodesic line (L{Line<Intersector.Line>}).
           @arg glB: An other geodesic line (L{Line<Intersector.Line>}).
           @kwarg eps: Equality margin (C{degrees}).

           @return: The intersection (L{XDist}) or C{None} if none found.

           @raise GeodesicError: Geodesic line B{C{glA}} or B{C{glB}}
                                 invalid, incompatible, ill-configured or
                                 do not have near-equal C{(lat1, lon1)}.

           @raise IntersectionError: No convergence.

           @note: Offset C{X0} is implicit, zeros.
        '''
        self._xLines(glA, glB)
        a = glA.lat1 - glB.lat1
        b = glA.lon1 - glB.lon1
        if fabs(a) > eps or fabs(b) > eps:
            raise GeodesicError(lat1=a, lon1=b, eps=eps)
        return self._Next(glA, glB)

    def Next4(self, glA, glB, eps=_EPS3):
        '''Yield the next intersection of two I{intersecting} geodesic lines.

           @return: An L{Intersector4Tuple}C{(A, B, sAB, c)} or C{None} if
                    none found.

           @see: Method L{Next} for further details.
        '''
        X = self.Next(glA, glB, eps=eps)
        return X if X is None else self._In4T(glA, glB, X, X)

    def Next4s(self, glA, glB, X0=_X000, aMax=1801, sMax=0, avg=False, **delta):
        '''Yield C{Next} intersections up to a maximal (angular) distance.

           @kwarg aMax: Upper limit for the angular distance (C{degrees}).
           @kwarg sMax: Upper limit for the distance (C{meter}).
           @kwarg avg: If C{True}, set the next intersection lat- and longitude
                       to the mid-point of the previous ones (C{bool}).
           @kwarg delta: Optional, threshold C{B{delta}=Delta} for the
                         incremental distance (C{meter}).

           @return: Yield an L{Intersector4Tuple}C{(A, B, sAB, c)} for
                    each intersection found.

           @see: Methods L{Next4} for further details.
        '''
        X = self.Closest(glA, glB, X0=X0)
        if X is not None:
            d = _xkwds_get(delta, delta=self.Delta) * _2_0
            S,  _L, _abs = X, self._Line, fabs
            while True:
                A, B, _, _ = r = self._In4T(glA, glB, S, X)
                yield r
                if (aMax and (_abs(A.a12) > aMax or _abs(B.a12) > aMax)) or \
                   (sMax and (_abs(A.s12) > sMax or _abs(B.s12) > sMax)):
                    break
                latA, lonA = A.lat2, A.lon2
                latB, lonB = B.lat2, B.lon2
                if avg:
                    latA = latB = favg(latA, latB)
                    lonA = lonB = favg(lonA, lonB)
                X = self._Next(_L(latA, lonA, A.azi2),
                               _L(latB, lonB, B.azi2))
                if X is None or X.L1() < d:
                    break
                S += X.sA, X.sB

    def _Next(self, glA, glB):
        '''(INTERNAL) Find the next intersection.
        '''
        i, X0, Q, d  = 1, _X000, _XINF, INF
        D, D2, T1D2D = self._RPI, self._D2, self._T1D2Delta
        sk, K        = set(), len(XDist._2A)
        for n in (n for n in range(K) if n not in sk):
            X    = X0._nD2(n, D2)
            X, i = self._Basic2(glA, glB, X, i)
            X    = X0._fixCoincident(X)
            c, z = X.c, (X.L1(X0) <= self.Delta)  # X == X0
            if z:
                if not c:  # next n
                    continue
                for s in (-D, 0, D):
                    s =  self._conjDist(glA, s, semi=False)
                    t = _L1(s, s * c)
                    if t < d:
                        Q, d = XDist(s, s * c, c), t
                        i += 1
            else:
                t = X.L1()
                if t < d:
                    Q, d = X, t
                    i += 1
            n += 1
            for s in ((_1_1t if z else _1_0_1t)
                             if c else _0t):
                if c and s:
                    s *= D2
                    T  = XDist(s, s * c)
                    T += X
                else:
                    T  = X
                for m in range(n, K):
                    if m not in sk and \
                       T.L1(X0._nD2(m, D2)) < T1D2D:
                        sk.add(m)
        return None if Q is _XINF else Q.set_(sX0=d, iteration=i)

    def _obliqDist4(self):
        zx = 45.0
        if self.f:
            _abs, _cjD5 = fabs, self._conjDist5

            _,  ds0, z0, _,   _   = _cjD5(zx + _1_0)
            s1, ds1, z1, sAx, sBx = _cjD5(zx - _1_0)
            sx, dsx, zx = s1, _abs(ds1), z1
            # find ds(azi) = 0 by secant method
            for _ in range(16):
                if ds1 == ds0:
                    break
                z = (z0 * ds1 - z1 * ds0) / (ds1 - ds0)
                _,  ds0, z0 = s1, ds1, z1
                s1, ds1, z1, a, b = _cjD5(z)
                if _abs(ds1) < dsx:
                    sx, dsx, zx, sAx, sBx = s1, _abs(ds1), z, a, b
                    if not dsx:
                        break
        else:
            sx, sAx, sBx = self._RPI, _0_5, -_1_5
        return sx, zx, sAx, sBx

    def _polarB3(self, lats=False):  # PYCHOK no cover
        latx =  64.0
        lat  = _90_0 - latx
        if self.f:
            _d, _pD2 =  fdot, self._polarDist2

            s0, lat0 = _pD2(latx - _1_0)
            s1, lat1 = _pD2(latx + _1_0)
            s2, lat2 = \
            sx, latx = _pD2(latx)
            prolate  =  self.f < 0
            # solve for ds(lat) / dlat = 0 with a quadratic fit
            for _ in range(_TRIPS):
                v = (lat1 - lat0), (lat0 - lat2), (lat2 - lat1)
                d = _d(v, s2, s1, s0) * _2_0
                if not d:  # or isnan(d)
                    break
                lat = _d(v, (lat1 + lat0) * s2,
                            (lat0 + lat2) * s1,
                            (lat2 + lat1) * s0) / d
                s0, lat0 =  s1, lat1
                s1, lat1 =  s2, lat2
                s2, lat2 = _pD2(lat)
                if (s2 < sx) if prolate else (s2 > sx):
                    sx, latx = s2, lat2
            if lats:
                _, lat = _pD2(latx, lat2=True)
            sx += sx
        else:
            sx  = self._RPI
        return sx, latx, lat

    def _polarDist2(self, lat1, lat2=False):
        gl = self._Line(lat1=lat1)
        s  = self._conjDist(gl, self._faPI_2, semi=True)
        if lat2:
            lat1 = gl.Position(s, outmask=Caps.LATITUDE).lat2
        return s, lat1

    def _Position(self, gl, s):
        return gl.Position(s, outmask=Caps._STD_LINE)

    @Property_RO
    def R(self):
        '''Get the I{authalic} earth radius (C{meter}).
        '''
        return self.ellipsoid.R2

    @Property_RO
    def _RPI(self):  # normalizer, semi-perimeter, C++ _d
        return self.R * PI  # ~20K Km WGS84

    def _Spherical(self, glA, glB, T):
        '''(INTERNAL) Get solution based from a spherical triangle.
        '''
        # threshold for coincident geodesics and intersections;
        # this corresponds to about 4.3 nm on WGS84.
        A = self._Position(glA, T.sA)
        B = self._Position(glB, T.sB)
        D = self._Inverse(A, B)

        # a = interior angle at A, b = exterior angle at B
        a, da = _diff182(A.azi2, D.azi1)
        b, db = _diff182(B.azi2, D.azi2)
        c, dc = _diff182(a, b)
        if fsum1_(dc, db, -da, c) < 0:  # inverted triangle
            a, da = -a, -da
            b, db = -b, -db
        sa, ca = _sincos2de(a, da)
        sb, cb = _sincos2de(b, db)

        e, z, _abs = _EPS3, D.s12, fabs
        if _abs(z) <= self._EPS3R:
            sA = sB = 0  # at intersection
            c  = 1 if _abs(sa - sb) <= e and _abs(ca - cb) <= e else (
                -1 if _abs(sa + sb) <= e and _abs(ca + cb) <= e else 0)
        elif _abs(sa) <= e and _abs(sb) <= e:  # coincident
            sA =  ca * z * _0_5  # choose midpoint
            sB = -cb * z * _0_5
            c  =  1 if (ca * cb) > 0 else -1
            # alt1: A =  ca * z; B = 0
            # alt2: B = -cb * z; A = 0
        else:  # general case
            sz, cz = sincos2(z / self.R)
            # [SKIP: Divide args by |sz| to avoid possible underflow
            # in {sa, sb} * sz; this is probably not necessary].
            # Definitely need to treat sz < 0 (z > pi*R) correctly in
            # order to avoid some convergence failures in _Basic2.
            sA = atan2(sb * sz,  sb * ca * cz - cb * sa) * self.R
            sB = atan2(sa * sz, -sa * cb * cz + ca * sb) * self.R
            c  = 0
        return XDist(sA, sB, c)

    @Property_RO
    def _T1D2Delta(self):
        return self._T1 * _2_0 - self._D2 - self.Delta

    def _T1d3Delta(self, d3):
        return self._T1 * _2_0 - d3 - self.Delta

    @Property_RO
    def _Tol(self):  # convergence tolerance
        return self._RPI * pow(EPS, 0.75)  # _0_75

    def toStr(self, **prec_sep_name):  # PYCHOK signature
        '''Return this C{Intersector} as string.

           @see: L{Ellipsoid.toStr<pygeodesy.ellipsoids.Ellipsoid.toStr>}
                 for further details.

           @return: C{Intersector} (C{str}).
        '''
        return self._instr(props=(Intersector.geodesic,), **prec_sep_name)

    def _xLines(self, glA, glB):
        # check two geodesic lines vs this geodesic
        _xinstanceof(*self._GeodesicLines, glA=glA, glB=glB)
        C, g = Caps.LINE_CAPS, self.geodesic
        for gl in (glA, glB):
            _xgeodesics(gl.geodesic, g, Error=GeodesicError)
            c = gl.caps & C
            if c != C:
                raise GeodesicError(caps=c, LINE_CAPS=C)  # _invalid_


class Intersector4Tuple(_NamedTuple):
    '''4-Tuple C{(A, B, sAB, c)} with C{A} and C{B} the C{Position}
       of the intersection on each geodesic line, the distance C{sAB}
       between C{A} and C{B} in C{meter} and the coincidence indicator
       C{c} (C{int}), see L{XDist}.

       @note: C{A} and C{B} are each a C{GeodesicLine...Position} for
              C{outmask=Caps.STANDARD} with the intersection location
              in C{lat2}, C{lon2}, azimuth in C{azi2}, the distance
              C{s12} in C{meter} and the angular distance C{a12} in
              C{degrees}.
    '''
    _Names_ = (_A_,   _B_,  'sAB', _c_)
    _Units_ = (_Pass, _Pass, Meter, Int)


class _List(list):

    Delta = 0  # equality margin

    def __init__(self, Delta):
        self.Delta = Delta
#       list.__init__(self)

    def __contains__(self, other):
        # handle C{if X in this: ...}
        a, b   = other.sA, other.sB
        D, _D1 = self.Delta, _L1
        for X in self:
            if _D1(X.sA - a, X.sB - b) <= D:
                return True
        return False

    def add(self, X, *L1_i):
        # append an item, updated
        if L1_i:
            d, i = L1_i
            X.set_(sX0=d, iteration=i)
        self.append(X)
        return X.sX0

    def sortrim(self, X0, sMax):
        # trim and sort the X items

        def _key(Xk):
            _, k = Xk
            return k  # rank of X

        for X, _ in sorted(self.trim(X0, sMax), key=_key):
            yield X  # de-tuple (X, k)

    def trim(self, X0, sMax):
        # trim and yield 2-tuple (X, rank)
        a, b, _eu = X0.sA, X0.sB, euclid

        for X in self:
            k = X.sX0
            if k <= sMax:
                k += _eu(X.sA - a, X.sB - b)
                yield X, k  # rank of X


if __name__ == '__main__':

    from pygeodesy import GeodesicExact, printf

    I = Intersector(GeodesicExact(), name='Test')  # PYCHOK I

    # <https://GeographicLib.sourceforge.io/C++/doc/classGeographicLib_1_1Intersect.html>
    a = I.Line( 0,  0,  45)
    b = I.Line(45, 10, 135)
    printf('Closest: %r', I.Closest(a, b))
    printf('Closest4: %r', I.Closest4(a, b), nt=1)

    for i, t in enumerate(I.Next4s(a, b)):
        printf('Next4s %s: %r (%s)', i, t, t.iteration)

    # <https://GeographicLib.sourceforge.io/C++/doc/IntersectTool.1.html>
    a = I.Line(50, -4, -147.7)
    b = I.Line( 0,  0,   90)
    printf('Closest: %r', I.Closest(a, b), nl=1)
    printf('Closest4: %r', I.Closest4(a, b), nt=1)

    a = I.Line(50,  -4, -147.7)
    b = I.Line( 0, 180,    0)
    for i, X in enumerate(I.All(a, b)):
        printf('All %s: %r (%s)', i, X, X.iteration)
        if i > 9:
            break
    printf('')
    for i, t in enumerate(I.All4(a, b)):
        printf('All4 %s: %r (%s)', i, t, t.iteration)
        if i > 9:
            break

    # % echo 50N 4W 147.7W 0 0 90 | IntersectTool -e 6371e3 0 -c -p 0 -C
    # 6077191 -3318019 0
    # -0.00000 -29.83966 -0.00000 -29.83966 0
    I = Intersector(GeodesicExact(6371e3, 0), name='Test')  # PYCHOK I
    a = I.Line(50, -4, -147.7)
    b = I.Line( 0,  0,   90)
    printf('Closest: %r', I.Closest(a, b), nl=1)
    printf('Closest4: %r', I.Closest4(a, b))

# **) MIT License
#
# Copyright (C) 2024-2024 -- mrJean1 at Gmail -- All Rights Reserved.
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
