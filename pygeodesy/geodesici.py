
# -*- coding: utf-8 -*-

u'''Classes L{Intersectool} and L{Intersector} to find the intersections of two geodesic lines or line segments.

Class L{Intersector} is a pure Python version of I{Karney}'s C++ class U{Intersect
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Intersect.html>}.

Class L{Intersectool} is a wrapper to invoke I{Karney}'s U{IntersectTool
<https://GeographicLib.SourceForge.io/C++/doc/IntersectTool.1.html>} utility, but intended I{for testing purposes only}.

Set env variable C{PYGEODESY_INTERSECTTOOL} to the (fully qualified) path of the C{IntersectTool} executable.  For usage
and some examples run C{"env PYGEODESY_INTERSECTTOOL=<IntersectTool-path> python3 -m pygeodesy.geodesici --help"}.

Both L{Intersectool} and L{Intersector} provide methods C{All}, C{Closest}, C{Next} and C{Segment} and produce
L{XDict} instances with 4 or more items.  Adjacent methods C{All5}, C{Closest5}, C{Next5} and C{Segment} return
or yield L{Intersectool5Tuple} or L{Intersector5Tuple}s with the lat-, longitude and azimuth of each intersection
as an extended, geodesic C{Position}-like L{GDict} instance.

For more details, see the C++ U{GeographicLib<https://GeographicLib.SourceForge.io/C++/doc/index.html>}
documentation, I{Charles F.F. Karney}'s paper U{Geodesics intersections<https://arxiv.org/abs/2308.00495>}
and I{S. Baselga Moreno & J.C. Martinez-Llario}'s U{Intersection and point-to-line solutions for geodesics
on the ellipsoid<https://riunet.UPV.ES/bitstream/handle/10251/122902/Revised_Manuscript.pdf>}.
'''
# make sure int/int division yields float quotient
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import _copy, _enumereverse, map1, \
                             _xinstanceof, _xor
from pygeodesy.constants import EPS, INF, INT0, PI, PI2, PI_4, \
                               _0_0, _0_5, _1_0, _1_5, _2_0, \
                               _3_0, _64_0, _90_0, isfinite, \
                               _EPSjam  # PYCHOK used!
from pygeodesy.ellipsoids import _EWGS84,  Fmt, unstr
from pygeodesy.errors import GeodesicError, IntersectionError, _an, \
                            _xgeodesics, _xkwds_get, _xkwds_kwds, \
                            _xkwds_pop2
# from pygeodesy.errors import exception_chaining  # _MODS
from pygeodesy.fmath import euclid, fdot
from pygeodesy.fsums import Fsum, fsum1_,  _ceil
from pygeodesy.interns import NN, _A_, _B_, _c_, _COMMASPACE_, _HASH_, \
                             _M_, _not_, _SPACE_, _too_
from pygeodesy.karney import Caps, _diff182, GDict, _sincos2de, _Xables
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import ADict, _NamedBase, _NamedTuple, _Pass
# from pygeodesy.namedTuples import _LL4Tuple  # _MODS
from pygeodesy.props import deprecated_method, Property, \
                            Property_RO, property_RO, property_ROver
from pygeodesy.solveBase import _SolveCapsBase,  pairs
# from pygeodesy.streprs import pairs  # from .solveBase
# from pygeodesy.streprs import Fmt, unstr  # from .ellipsoids
from pygeodesy.units import Azimuth as Azi, Degrees, Float, Int, \
                           _isDegrees,  Lat, Lon, Meter, Meter_
from pygeodesy.utily import atan2, sincos2,  fabs, radians

# from math import ceil as _ceil, fabs, radians  # .fsums, .utily

__all__ = _ALL_LAZY.geodesici
__version__ = '24.11.24'

_0t     =  0,  # int
_1_1t   = -1, +1
_1_0_1t = -1, 0, +1
_aAB_   = 'aAB'
_c__    = '-c'  # PYCHOK used!
_cWGS84 = _EWGS84.a * PI2  # outer circumference
_EPS3   =  EPS * _3_0
_EPSr5  =  pow(EPS, 0.2)  # PYCHOK used! 7.4e-4 or ~3"
_i__    = '-i'  # PYCHOK used!
_latA_  = 'latA'
_lonA_  = 'lonA'
_n__    = '-n'  # PYCHOK used!
_o__    = '-o'  # PYCHOK used!
_R__    = '-R'
_sAB_   = 'sAB'
_sX0_   = 'sX0'
_TRIPS  =  128


class XDict(ADict):
    '''4+Item result from L{Intersectool} and L{Intersector} methods
       C{All}, C{Closest}, C{Next} and C{Segment} with the intersection
       offsets C{sA}, C{sB} and C{sX0} in C{meter} and the coincidence
       indicator C{c}, an C{int}, +1 for parallel, -1 for anti-parallel
       or 0 otherwise.

       Offsets C{sA} and C{sB} are distances measured I{along} geodesic
       line C{glA} respectively C{glB}, but C{sX0} is the I{L1-distance}
       between the intersection and the I{origin} C{X0}.

       If present, distance C{sAB} and angular distance C{aAB} represent
       the difference between the intersection point on geodesic lines
       C{glA} and C{glB} in C{meter} respectively C{degrees}, typically
       below C{5e-9 meter} or C{5 nm} and C{5e-14 degrees} or C{1 n"}.

       For segments, indicators C{kA} and C{kB} are C{0} if the segments
       intersect or C{-1} or C{+1} if the intersection is I{before} the
       start, respectively I{after} the end of the segment, similar to
       L{Intersection3Tuple<Intersection3Tuple>}.  Segment indicator
       C{k} is I{Karney}'s C{segmode}, equal C{kA * 3 + kB}.
    '''
    _Delta = EPS  # default margin, see C{Intersector._Delto}

    def __add__(self, other):
        X  = _copy(self)
        X +=  other
        return X

    def __eq__(self, other):
        return not self.__ne__(other)

    def __iadd__(self, other):
        if isinstance(other, tuple):  # and len(other) == 2:
            a, b = other
        else:
            # _xinstanceof(XDict, other=other)
            a = other.sA
            b = other.sB
            if other.c:
                self.c = other.c
        self.sA += a  # PYCHOK sA
        self.sB += b  # PYCHOK sB
        return self

    def __le__(self, other):
        # _xinstanceof(XDict, other=other)
        return self == other or self < other

    def __lt__(self, other):
        # _xinstanceof(XDict, other=other)
        return (self.sA < other.sA or (self.sA == other.sA and  # PYCHOK sA
                self.sB < other.sB) and self != other)  # PYCHOK sB

    def __ne__(self, other):
        # _xinstanceof(XDict, other=other)
        return self is not other and self.L1(other) > self._Delta

    def _corners(self, sA, sB, T2):
        # yield all corners further than C{T2}
        a, b = self.sA, self.sB  # PYCHOK sA, sB
        for x in (0, sA):
            for y in (0, sB):
                if _L1(x - a, y - b) >= T2:
                    yield XDict_(x, y)

    def _fixCoincident(self, X, c0=0):
        # return the mid-point if C{X} is anti-/parallel
        c = c0 or X.c
        if c:
            s = (self.sA - X.sA  +  # PYCHOK sA
                (self.sB - X.sB) * c) * _0_5  # PYCHOK sB
            X = X + (s, s * c)  # NOT +=
        return X

    def _fixSegment(self, sA, sB):  # PYCHOK no cover
        # modify this anti-/parallel C{XDict}
        a, b, c = self.sA, self.sB, self.c  # PYCHOK sA, sB, c

        def _g():  # intersection in smallest gap
            if c > 0:  # distance to [A, B] is |(a - b) - (A - B)|
                t = a - b  # consider corners [0, sB] and [sA, 0]
                t = fabs(t + sB) < fabs(t - sA)
                s = a + b
            else:  # distance to [A, B] is |(a + b) - (A + B)|
                t = a + b  # consider corner [0, 0] and [sA, sB]
                t = fabs(t) < fabs(t - (sA + sB))
                s = sB + (a - b)
            return (sB if t else sA) - s

        ta = -a
        tb = sA -  a
        tc = -c *  b
        td = -c * (b - sB)

        ga = 0 <= (b + c * ta) <= sB
        gb = 0 <= (b + c * tb) <= sB
        gc = 0 <= (a +     tc) <= sA
        gd = 0 <= (a +     td) <= sA

        # test opposite rectangle sides first
        s = ((ta + tb) if ga and gb else (
             (tc + td) if gc and gd else (
             (ta + tc) if ga and gc else (
             (ta + td) if ga and gd else (
             (tb + tc) if gb and gc else (
             (tb + td) if gb and gd else _g())))))) * _0_5
        self += s, s * c

    @property_RO
    def _is00(self):
        return not (self.sA or self.sB)  # PYCHOK sA, sB

    def L1(self, other=None):
        '''Return the C{L1} distance.
        '''
        a, b = self.sA, self.sB  # PYCHOK sA, sB
        if other is not None:
            # _xinstanceof(XDict, other=other)
            a -= other.sA
            b -= other.sB
        return _L1(a, b)

    def _nD1(self, D1):
        # yield the C{Closest} starts
        D_ = 0, D1, -D1
        for a, b in zip((0, 1, -1, 0,  0),
                        (0, 0,  0, 1, -1)):
            yield self + (D_[a], D_[b])

    def _nD2(self, D2):
        # yield the C{Next} starts
        D22 = D2 * _2_0
        D_  = 0, D2, D22, -D22, -D2
        for a, b in zip((-1, -1,  1, 1, -2, 0, 2,  0),
                        (-1,  1, -1, 1,  0, 2, 0, -2)):
            yield self + (D_[a], D_[b])

    def _nmD3(self, n, m, D3):  # d3 / 2
        # yield the C{All} starts
        yield self
        for i in range(n, m, 2):
            for j in range(n, m, 2):
                if i or j:  # skip self
                    yield self + ((i + j) * D3,
                                  (i - j) * D3)

    def _outSide(self, sA, sB):
        # is this C{Xdist} outside one or both segments?
        a, b = self.sA, self.sB  # PYCHOK sA, sB
        kA = -1 if a < 0 else (+1 if a > sA else INT0)
        kB = -1 if b < 0 else (+1 if b > sB else INT0)
        self.set_(kA=kA, kB=kB, k=(kA * 3 + kB) or INT0)
        return bool(kA or kB)

    def _skip(self, S_, T1_Delta):
        # remove starts from list C{S_} near this C{XDict}
        for j, S in _enumereverse(S_):
            if S.L1(self) < T1_Delta:
                S_.pop(j)


def XDict_(sA=_0_0, sB=_0_0, c=INT0, sX0=_0_0):
    '''(INTERNAL) New L{XDict} from positionals.
    '''
    return XDict(sA=sA, sB=sB, c=c, sX0=sX0)

_X000 = XDict_()  # PYCHOK origin
_XINF = XDict_(INF)


class _IntersectBase(_NamedBase):
    '''(INTERNAL) Base class for L{Intersectool} and L{Intersector}.
    '''
    # _g  = None

    def __init__(self, geodesic, **name):
        _xinstanceof(*_EWGS84._Geodesics, geodesic=geodesic)
        self._g = geodesic
        if name:
            self.name = name

    @Property_RO
    def a(self):
        '''Get the I{equatorial} radius, semi-axis (C{meter}).
        '''
        return self.ellipsoid.a

    equatoradius = a  # = Requatorial

    def All(self, glA, glB, **kwds):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}.'''
        self._notOverloaded(glA, glB, **kwds)

    @Property_RO
    def _cHalf(self):  # normalizer, semi-circumference
        return self.R * PI  # ~20K Km WGS84

    @Property_RO
    def _cMax(self):  # outer circumference
        return max(self.a, self.ellipsoid.b, self.R) * PI2

    @property_RO
    def datum(self):
        '''Get the geodesic's datum (C{Datum}).
        '''
        return self.geodesic.datum

    @Property_RO
    def ellipsoid(self):
        '''Get the C{geodesic}'s ellipsoid (C{Ellipsoid}).
        '''
        return self.geodesic.datum.ellipsoid

    @Property_RO
    def f(self):
        '''Get the I{flattening} (C{scalar}), C{0} for spherical, negative for prolate.
        '''
        return self.ellipsoid.f

    flattening = f

    @property_RO
    def geodesic(self):
        '''Get the C{geodesic} (C{Geodesic...}).
        '''
        return self._g

    def _illz2G(self, G, il):
        '''(INTERNAL) Set C{InverseLine} 1-/2-attrs into C{G}, a C{GDict}.
        '''
        try:
            G.set_(lat1=il.lat1, lon1=il.lon1, azi1=il.azi1, a12=il.a13,  # .Arc()
                   lat2=il.lat2, lon2=il.lon2, azi2=il.azi2, s12=il.s13)  # .Distance()
        except AttributeError:
            r = il.Position(il.s13, outmask=Caps._STD_LINE)  # isfinite(il.s13)
            G.set_(**r)
#           for n, v in r.items():
#               if not hasattr(il, n):
#                   setattr(il, n, v)
        return G

    def intersect7(self, start1, end1, start2, end2, X0=_X000, aMaX0=0, sMaX0=_cWGS84,
                                                               **LatLon_and_kwds):
        '''Yield the intersection points of two lines, each defined by two (ellipsoidal)
           points or by an (ellipsoidal) start point and an azimuth from North.

           @arg start1: Start point of the first line (C{LatLon}).
           @arg end1: End point of the first line (C{LatLon}) or the azimuth at the
                      B{C{start1}} point (compass C{degrees360}).
           @arg start2: Start point of the second line (C{LatLon}).
           @arg end2: End point of the second line (C{LatLon}) or the azimuth at the
                      B{C{start2}} point (compass C{degrees360}).
           @kwarg X0: Optional I{origin} for I{L1-distances} (L{XDict}) or C{None} for
                      the L{Middle<Intersector.Middle>}, otherwise C{XDiff_(0, 0)}.
           @kwarg aMaX0: Upper limit for the I{angular L1-distance}
                         (C{degrees}) or C{None} or C{0} for unlimited.
           @kwarg sMaX0_C: Optional, upper limit C{B{sMaX0}=2*PI*R} for the
                        I{L1-distance} to B{C{X0}} (C{meter}).
           @kwarg LatLon_and_kwds: Optional class C{B{LatLon}=None} to return intersection
                         points and optional, additional B{C{LatLon}} keyword arguments.

           @note: The C{lat} and C{lon} attr of B{C{start1}}, B{C{end1}}, B{C{start2}} and
                  B{C{end2}} are used I{verbatim}, ignoring C{datum} or C{ellipsoid}.

           @return: Yield an L{Intersect7Tuple}C{(A, B, sAB, aAB, c, kA, kB)} for every
                    intersection found, with C{A} and C{B} each a B{C{LatLon}} or if
                    C{B{LatLon} is None} or not specified, a L{LatLon4Tuple}C{(lat, lon,
                    height, datum)} with C{height 0} and this C{datum}.

           @raise GeodesicError: Invalid B{C{start1}}, B{C{end1}}, B{C{start2}} or
                                 B{C{end2}} or B{C{end1}} and B{C{end2}} differ in type.

           @raise IntersectionError: No convergence.
        '''

        def _args(s, e):
            t = (e,) if _isDegrees(e) else (e.lat, e.lon)
            return (s.lat, s.lon) + t

        try:
            glA = self.Line(*_args(start1, end1))
            glB = self.Line(*_args(start2, end2))
        except Exception as x:
            raise GeodesicError(start1=start1, end1=end1, start2=start2, end2=end2, cause=x)

        LL, kwds = _xkwds_pop2(LatLon_and_kwds, LatLon=None)
        d,  kwds = _xkwds_pop2(kwds, datum=self.datum)
        h,  kwds = _xkwds_pop2(kwds, height=0)

        _LL4T = _MODS.namedTuples._LL4Tuple
        for X in self.All(glA, glB, X0=X0, aMaX0=aMaX0, sMaX0=sMaX0, _C=True):
            A = B = _LL4T(X.latA, X.lonA, h, d, LL, kwds, iteration=X.iteration)
            if X.sAB or X.latA != X.latB or X.lonA != X.lonB:
                B = _LL4T(X.latB, X.lonB, h, d, LL, kwds, iteration=X.iteration)
            yield Intersect7Tuple(A, B, X.sAB, X.aAB, X.c, _xkwds_get(X, kA=0),
                                                           _xkwds_get(X, kB=0))

    def _Inversa12(self, A, B=None):
        lls = (0, 0, A, 0) if B is None else (A.lat2, A.lon2,
                                              B.lat2, B.lon2)
        r = self._g.Inverse(*lls, outmask=Caps.DISTANCE)
        return r.s12, r.a12  # .a12 always in r

    def k2kAkB(self, k):
        '''Unravel C{k} into C{kA} and C{kB}.

           @arg k: Segment indicator C{kA * 3 + kB} (C{int}).

           @return: An C{ADict(k=k, kA=kA, kB=kB)}.

           @raise GeodesicError: Invalid B{C{k}}.
        '''
        for kA in range(-1, 2):
            for kB in range(-1, 2):
                if (kA * 3 + kB) == k:
                    return ADict(k=k, kA=kA, kB=kB)
        raise GeodesicError(k=k)

#   def k2kAkB(self, k):
#       # unravel C{k} into C{kA} and C{kB}.
#       kA, kB = divmod(k, 3)
#       if kB > 1:
#           kA += 1
#           kB -= 3
#       return kA, kB

    def Line(self, lat1, lon1, azi1_lat2, *lon2, **name):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}.'''
        self._notOverloaded(lat1, lon1, azi1_lat2, *lon2, **name)

    def _ll3z4ll(self, lat1, lon1, azi1_lat2, *lon2):
        t = Lat(lat1=lat1), Lon(lon1=lon1)
        if lon2:  # get azis for All, keep lat-/lons
            t += Lat(lat2=azi1_lat2), Lon(lon2=lon2[0])
        else:
            t += Azi(azi1=azi1_lat2),
        return t

    @deprecated_method
    def Next5s(self, glA, glB, X0=_X000, aMax=1801, sMax=0, **unused):  # PYCHOK no cover
        '''DEPRECATED on 2024.07.02, use method C{All5}.'''
        return self.All5(glA, glB, X0=X0, aMaX0=aMax, sMaX0=sMax)  # PYCHOK attr

    @Property_RO
    def R(self):
        '''Get the I{authalic} earth radius (C{meter}).
        '''
        return self.ellipsoid.R2

    def _sMaX0_C2(self, aMaX0=0, **sMaX0_C):
        _g = _xkwds_get
        s  = _g(sMaX0_C, sMaX0=self._cMax)
        s  = _g(sMaX0_C, sMax=s)  # for backward ...
        a  = _g(sMaX0_C, aMax=aMaX0)  # ... compatibility
        if a:  # degrees to meter, approx.
            s = min(s, self.R * radians(a))  # ellipsoid.degrees2m(a)
        s  = _g(sMaX0_C, _R=s)
        if s < _EPS3:
            s = _EPS3  # raise GeodesicError(sMaX0=s)
        return s, _g(sMaX0_C, _C=False)

    def _xNext(self, glA, glB, eps1, **eps_C):  # PYCHOK no cover
        eps1 = _xkwds_get(eps_C, eps=eps1)  # eps for backward compatibility
        if eps1 is not None:
            a = glA.lat1 - glB.lat1
            b = glA.lon1 - glB.lon1
            if euclid(a, b) > eps1:
                raise GeodesicError(lat_=a, lon_=b, eps1=eps1)
        return _xkwds_kwds(eps_C, _C=False)


class Intersectool(_IntersectBase, _SolveCapsBase):
    '''Wrapper to invoke I{Karney}'s utility U{IntersectTool
       <https://GeographicLib.SourceForge.io/C++/doc/IntersectTool.1.html>}
       similar to class L{Intersector<geodesici.Intersector>}.

       @note: Use property C{IntersectTool} or env variable C{PYGEODESY_INTERSECTTOOL}
              to specify the (fully qualified) path to the C{IntersectTool} executable.

       @note: This C{Intersectool} is intended I{for testing purposes only}, it invokes
              the C{IntersectTool} executable for I{every} method call.
    '''
    _c_alt       = _c__,  # Closest latA lonA aziA  latB lonB aziB
    _C_option    = '-C',
    _Error       =  GeodesicError
    _i_alt       = _i__,  # Segment latA1 lonA1 latA2 lonA2  latB1 lonB1 latB2 lonB2
    _linelimit   =  1200  # line printer width X 10
    _n_alt       = _n__,  # Next latA lonA aziA  aziB
    _Names_ABs   = _latA_, _lonA_, 'latB', 'lonB', _sAB_  # -C to stderr
    _Names_XDict = 'sA', 'sB', _c_  # plus 'k' from -i or 'sX0' from -R
    _o_alt       = _o__,  # Offset latA lonA aziA  latB lonB aziB  x0 y0
    _Xable_name  = _Xables.IntersectTool.__name__
    _Xable_path  = _Xables.IntersectTool()

    def __init__(self, a_geodesic=None, f=None, **name):
        '''New L{IntersectTool}.

           @arg a_geodesic: Earth' equatorial axis (C{meter}) or a geodesic
                            (L{GeodesicExact<pygeodesy.geodesicx.GeodesicExact>},
                            wrapped L{Geodesic<pygeodesy.geodesicw.Geodesic>} or
                            L{GeodesicSolve<pygeodesy.geodsolve.GeodesicSolve>}).
           @kwarg f: Earth' flattening (C{scalar}), required if B{C{a_geodesic}}
                     is in C{meter}, ignored otherwise.
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @raise GeodesicError: The eccentricity of the B{C{geodesic}}'s ellipsoid is too
                                 large or no initial convergence.

           @see: The B{Note} at I{Karney}'s C++ U{Intersect<https://GeographicLib.sourceforge.io/
                 C++/doc/classGeographicLib_1_1Intersect.html#ae41f54c9a44836f6c8f140f6994930cf>}.
        '''
        g = self._GeodesicExact() if a_geodesic is None else (a_geodesic if f is None else
            self._GeodesicExact(a_geodesic, f))
        _IntersectBase.__init__(self, g, **name)

    def All(self, glA, glB, X0=_X000, eps1=_0_0, aMaX0=0, **sMaX0_C):  # PYCHOK signature
        '''Yield all intersection of two geodesic lines up to a limit.

           @kwarg eps1: Optional margin for the L{euclid<pygeodesy.euclid>}ean distance
                        (C{degrees}) between the C{(lat1, lon1)} points of both lines for
                        using the L{IntersectTool<Intersectool.IntersectTool>}'s C{"-n"}
                        option, unless C{B{eps1}=None}.

           @return: An L{XDict} for each intersection.
        '''
        for X, _ in self._All2(glA, glB, X0, eps1, aMaX0=aMaX0, **sMaX0_C):
            yield X

    def _All2(self, glA, glB, X0, eps1, **aMaX0_sMaX0_C):  # MCCABE 13
        '''(INTERNAL) Helper for methods C{.All} and C{.All5}.
        '''
        def _xz2(**gl):
            try:
                n, gl = gl.popitem()  # _xkwds_item2(gl)
                try:
                    return self._c_alt, (gl.azi1,)
                except (AttributeError, KeyError):
                    return self._i_alt, (gl.lat2, gl.lon2)
            except Exception as x:
                raise GeodesicError(n, gl, cause=x)

        _t, a = _xz2(glA=glA)
        _x, b = _xz2(glB=glB)
        if _x is not _t:
            raise GeodesicError(glA=glA, glB=glB)

        A = glA.lat1, glA.lon1
        B = glB.lat1, glB.lon1
        if _x is self._c_alt:
            if X0 is _X000 or X0._is00:
                if eps1 is not None and \
                   euclid(glA.lat1 - glB.lat1,
                          glA.lon1 - glB.lon1) <= eps1:
                    _x, B = self._n_alt, ()
            else:  # non-zero offset
                _x = self._o_alt
                b += X0.sA, X0.sB

        sMaX0, _C = self._sMaX0_C2(**aMaX0_sMaX0_C)
        for X in self._XDictInvoke(_x, _sX0_, (A + a + B + b),
                                   _C=_C, _R=sMaX0):
            if _C:
                T = self._In5T(glA, glB, X, X)
                if _aAB_ not in X:
                    X.set_(sAB=T.sAB, aAB=T.aAB)
            else:
                T = None
            yield X.set_(c=int(X.c)), T

    def All5(self, glA, glB, X0=_X000, **aMaX0_sMaX0):
        '''Yield all intersection of two geodesic lines up to a limit.

           @return: An L{Intersectool5Tuple} for each intersection.
        '''
        for _, T in self._All2(glA, glB, X0, _0_0, _C=True, **aMaX0_sMaX0):
            yield T

    @Property_RO
    def _cmdBasic(self):
        '''(INTERNAL) Get the basic C{IntersectTool} cmd (C{tuple}).
        '''
        return (self.IntersectTool,) + (self._e_option +
                                        self._E_option +
                                        self._p_option)

    def Closest(self, glA, glB, X0=_X000, _C=False):
        '''Find the closest intersection of two geodesic lines.

           @kwarg _C: Use C{B{_C}=True} to include the C{"-C"} results (C{bool}).

           @return: An L{XDict}.
        '''
        args = glA.lat1, glA.lon1, glA.azi1, \
               glB.lat1, glB.lon1, glB.azi1
        if X0 is _X000 or X0._is000:
            _x = self._c_alt
        else:
            _x = self._o_alt
            args += X0.sA, X0.sB
        return self._XDictInvoke(_x, NN, args, _C=_C)  # _R=None)

    def Closest5(self, glA, glB, **unused):
        '''Find the closest intersection of two geodesic lines.

           @return: An L{Intersectool5Tuple}.
        '''
        X = self.Closest(glA, glB, _C=True)
        return self._In5T(glA, glB, X, X)

    @property_ROver
    def _GeodesicExact(self):
        '''Get the I{class} L{GeodesicExact}, I{once}.
        '''
        return _MODS.geodesicx.GeodesicExact  # overwrite property_ROver

    def _In5T(self, glA, glB, S, X, k2=False, **_2X):
        A = GDict(glA).set_(lat2=X.latA, lon2=X.lonA, s12=S.sA)
        B = GDict(glB).set_(lat2=X.latB, lon2=X.lonB, s12=S.sB)
        if k2:
            A.set_(k2=X.kA)
            B.set_(k2=X.kB)
        s, a =  self._Inversa12(A, B)
        sAB  = _xkwds_get(X, sAB=s)
        if a and s and s != sAB:
            a *= sAB / s  # adjust a
        return Intersectool5Tuple(A._2X(glA, **_2X),
                                  B._2X(glB, **_2X), sAB, a, X.c)

    @Property
    def IntersectTool(self):
        '''Get the U{IntersectTool<https://GeographicLib.SourceForge.io/C++/doc/IntersectTool.1.html>}
           executable (C{filename}).
        '''
        return self._Xable_path

    @IntersectTool.setter  # PYCHOK setter!
    def IntersectTool(self, path):
        '''Set the U{IntersectTool<https://GeographicLib.SourceForge.io/C++/doc/IntersectTool.1.html>}
           executable (C{filename}), the (fully qualified) path to the C{IntersectTool} executable.

           @raise GeodesicError: Invalid B{C{path}}, B{C{path}} doesn't exist or isn't the
                                 C{IntersectTool} executable.
        '''
        self._setXable(path)

    def Line(self, lat1, lon1, azi1_lat2, *lon2, **name):
        '''Return a geodesic line from this C{Intersector}'s geodesic, specified by
           two (goedetic) points or a (goedetic) point and an (forward) azimuth.

           @return: A 3- or 6-item, named L{GDict}.
        '''
        args = self._ll3z4ll(lat1, lon1, azi1_lat2, *lon2)
        gl   = GDict((u.name, u) for u in args)
#       if lon2:  # get azis for All, use lat-/lons as given
#           r = self._g.Inverse(outmask=Caps.AZIMUTH, *args)
#           gl.set_(azi1=Azi(azi1=r.azi1), azi2=Azi(azi2=r.azi2))
        if name:
            gl.name= name
        return gl

    def Middle(self, glA, glB, **_C):
        '''Get the mid-points on two geodesic line segments.

           @kwarg _C: Use C{B{_C}=True} to include the C{"-C"} results (C{bool}).

           @return: An L{XDict}.
        '''
        X, _, _, _, _ = self._middle5(glA, glB, **_C)
        return X

    def _middle5(self, glA, glB, _C=False, **unused):
        # return intersections C{A} and C{B} and the
        # center C{X0} of rectangle [sA, sB]

        def _smi4(**gl):
            try:
                n, gl = gl.popitem()
                il = self._g.InverseLine(gl.lat1, gl.lon1, gl.lat2, gl.lon2)
            except Exception as x:
                raise GeodesicError(n, gl, cause=x)
            s = il.s13
            m = s * _0_5
            return s, m, il, (il.Position(m, outmask=Caps._STD_LINE) if _C else None)

        sA, mA, iA, A = _smi4(glA=glA)
        sB, mB, iB, B = _smi4(glB=glB)
        X = XDict_(mA, mB)  # centers
        _ = X._outSide(sA, sB)
        if _C:  # _Names_ABs
            s, a = self._Inversa12(A, B)
            X.set_(latA=A.lat2, lonA=A.lon2, aMM=a,  # assert sA == A.s12
                   latB=B.lat2, lonB=B.lon2, sMM=s)  # assert sB == B.s12
        return X, A, iA, B, iB

    def Middle5(self, glA, glB, **unused):
        '''Get the mid-points on two geodesic line segments and their distance.

           @return: A L{Middle5Tuple}.
        '''
        X, A, iA, B, iB = self._middle5(glA, glB, _C=True)
        A, B,  s, a, c  = self._In5T(A, B, X, X, _2X=_M_)
        return Middle5Tuple(self._illz2G(A, iA),
                            self._illz2G(B, iB), s, a, c)

    def Next(self, glA, glB, eps1=None, **_C):  # PYCHOK no cover
        '''Find the next intersection of two I{intersecting} geodesic lines.

           @kwarg _C: Use C{B{_C}=True} to include the option C{"-C"} results (C{bool}).

           @return: An L{XDict}.
        '''
        if eps1 or _C:
            _C = self._xNext(glA, glB, eps1, **_C)
        return self._XDictInvoke(self._n_alt, NN,
                                      (glA.lat1, glA.lon1, glA.azi1, glB.azi1),
                                    **_C)  # _R=None

    def Next5(self, glA, glB, **eps1):  # PYCHOK no cover
        '''Find the next intersection of two I{intersecting} geodesic lines.

           @return: An L{Intersectool5Tuple}.
        '''
        X = self.Next(glA, glB, _C=True, **eps1)
        return self._In5T(glA, glB, X, X)

    def _R_option(self, _R=None):
        '''(INTERNAL) Get the C{-R maxdist} option.
        '''
        return () if _R is None else (_R__, str(_R))  # -R maxdist

    def Segment(self, glA, glB, **_C_unused):
        '''Find the intersection between two geodesic line segments.

           @kwarg _C: Use C{B{_C}=True} to include the option C{"-C"} results (C{bool}).

           @return: An L{XDict}.
        '''
        X = self._XDictInvoke(self._i_alt, 'k',
                                   (glA.lat1, glA.lon1, glA.lat2, glA.lon2,
                                    glB.lat1, glB.lon1, glB.lat2, glB.lon2),
                                   _C=_xkwds_get(_C_unused, _C=False))  # _R=None
        try:
            ks = self.k2kAkB(int(X.k))
        except Exception as x:
            raise GeodesicError(glA=glA, glB=glB, X=str(X), cause=x)
        return X.set_(**ks)

    def Segment5(self, glA, glB, **unused):
        '''Find the next intersection of two I{intersecting} geodesic lines.

           @return: An L{Intersectool5Tuple}.
        '''
        X = self.Segment(glA, glB, _C=True)
        return self._In5T(glA, glB, X, X, k2=True)

    def toStr(self, prec=6, sep=_COMMASPACE_, **unused):  # PYCHOK signature
        '''Return this C{Intersectool} as string.

           @kwarg prec_sep: Keyword argumens C{B{prec}=6} and C{B{sep}=", "}
                      for the C{float} C{prec}ision, number of decimal digits
                      (0..9) and the C{sep}arator string to join.  Trailing
                      zero decimals are stripped for B{C{prec}} values of 1
                      and above, but kept for negative B{C{prec}} values.

           @return: Intersectool items (C{str}).
        '''
        d = dict(geodesic=self.geodesic, invokation=self.invokation,
                                         status=self.status,
                                         IntersectTool=self.IntersectTool)
        return sep.join(pairs(d, prec=prec))

    def _XDictInvoke(self, alt, _k_sX0, args, _C=False, **_R):
        '''(INTERNAL) Invoke C{IntersectTool}, return results as C{XDict} or
           a C{generator} if keyword argument C{B{_R}=sMaX0} is specified.
        '''
        # assert len(args) == {self._c_alt: 6,
        #                      self._i_alt: 8,
        #                      self._n_alt: 4,
        #                      self._o_alt: 8}.get(alt, len(args))
        cmd   = self._cmdBasic
        Names = self._Names_XDict  # has _c_ always
        if _k_sX0:
            Names += _k_sX0,
        if _C:
            cmd   += self._C_option
            Names += self._Names_ABs
        if _R:
            cmd   += self._R_option(**_R)
        X, _R = self._DictInvoke2(cmd + alt, args, Names, XDict, **_R)
        return X if _R else X.set_(c=int(X.c))  # generator or XDict


class Intersector(_IntersectBase):
    '''Finder of intersections between two goedesic lines, each an instance
       of L{GeodesicLineExact<pygeodesy.geodesicx.GeodesicLineExact>},
       wrapped L{GeodesicLine<pygeodesy.geodesicw.GeodesicLine>} or
       L{GeodesicLineSolve<pygeodesy.geodsolve.GeodesicLineSolve>}.

       @see: I{Karney}'s C++ class U{Intersect<https://GeographicLib.sourceforge.io/
             C++/doc/classGeographicLib_1_1Intersect.html#details>} for more details.
    '''

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
        _IntersectBase.__init__(self, geodesic, **name)
        E  = self.ellipsoid
        t1 = E.b * PI  # min distance between intersects
        t2 = self._polarDist2(_90_0)[0] * _2_0  # furthest, closest intersect
        t5 = self._Inversa12( _90_0)[0] * _2_0  # longest, shortest geodesic
        if self.f > 0:
            t3 = self._obliqDist4()[0]
            t4 = t1
        else:  # PYCHOK no cover
            t1, t2, t3 = t2, t1, t5
            t4,  _,  _ = self._polarB3()

        self._D1 = d1 = t2 * _0_5  # ~E.L tile spacing for Closest
        self._D2 = d2 = t3 / _1_5  # tile spacing for Next
        self._D3 = d3 = t4 - self.Delta  # tile spacing for All
        self._T1 = t1  # min distance between intersects
        self._T2 = t2 = t1 * _2_0
#       self._T5 = t5  # not used
        if not (d1 < d3 and d2 < d3 and d2 < t2):
            t = Fmt.PARENSPACED(_too_('eccentric'), E.e)
            raise GeodesicError(ellipsoid=E.toStr(terse=2), txt=t)

    def All(self, glA, glB, X0=None, aMaX0=0, **sMaX0_C):  # MCCABE 13
        '''Yield all intersection of two geodesic lines up to a limit.

           @arg glA: A geodesic line (L{Line<Intersector.Line>}).
           @arg glB: An other geodesic line (L{Line<Intersector.Line>}).
           @kwarg X0: Optional I{origin} for I{L1-distances} (L{XDict}) or
                      C{None} for the L{Middle<Intersector.Middle>} of both
                      lines if both are a 4-C{args} L{Line<Intersector.Line>}
                      or C{InverseLine}, otherwise C{XDiff_(0, 0)}.
           @kwarg aMaX0: Upper limit for the I{angular L1-distance}
                         (C{degrees}) or C{None} or C{0} for unlimited.
           @kwarg sMaX0_C: Optional, upper limit C{B{sMaX0}=2*PI*R} for the
                        I{L1-distance} to B{C{X0}} (C{meter}) and option
                        C{B{_C}=False} to include the intersection lat-/
                        longitudes C{latA}, C{lonA}, C{latB}, C{lonB} and
                        distances C{sAB} and C{aSB}.

           @return: Yield an L{XDict} for each intersection found.

           @raise GeodesicError: Geodesic line B{C{glA}} or B{C{glB}}
                                 invalid, incompatible or ill-configured.

           @raise IntersectionError: No convergence.
        '''
        self._xLines(glA, glB)
        if X0 is None:
            try:  # determine X0
                X0, _, _ = self._middle3(glA, glB, True)
            except GeodesicError:  # no .Distance
                X0 = _X000
        sMaX0, _C = self._sMaX0_C2(aMaX0, **sMaX0_C)

        D, _D = self.Delta, self._cHalf  # C++ _d
        xMaX0 = sMaX0 + D
        m     = int(_ceil(xMaX0 / self._D3))  # m x m tiles
        d3    = xMaX0 / m
        T2d3D = self._T2d3Delta(d3)

        C_ = _List(D)  # closest coincident
        X_ = _List(D)  # intersections found
        c0 =  0
        S_ =  list(X0._nmD3(1 - m, m, d3 * _0_5))
        # assert len(S_) == m * m + (m - 1) % 2
        while S_:
            Q, i = self._Basic2(glA, glB, S_.pop(0))
            if Q in X_:
                continue
            if Q.c:  # coincident intersection # PYCHOK no cover
                _X0fx =  X0._fixCoincident
                Q     = _X0fx(Q)  # Q = Q'
                if c0 and Q in C_:
                    continue
                C_.addend(Q)
                # elimate all existing intersections
                # on this line (which didn't set c0)
                c0 = Q.c
                for j, X in _enumereverse(X_):
                    if _X0fx(X, c0).L1(Q) <= D:  # X' == Q
                        X_.pop(j)

                a, s0 = len(X_), Q.sA
                args  = self._m12_M12_M21(glA, s0)
                _cjD  = self._conjDist
                for s in (-_D, _D):
                    s += s0
                    sa = 0
                    while True:
                        i += 1
                        sa = _cjD(glA, s + sa, *args) - s0
                        X  = Q + (sa, sa * c0)
                        if X_.addend(X, X0.L1(X), i) > xMaX0:
                            break

            elif c0 and Q in C_:  # Q.c == 0
                continue
            else:
                a = len(X_)

            X_.addend(Q, X0.L1(Q), i + 1)
            for X in X_[a:]:  # addended Xs
                X._skip(S_, T2d3D)

        return X_.sorter(sMaX0, self._C, glA, glB, _C=_C)  # generator

    def All5(self, glA, glB, X0=_X000, **aMaX0_sMaX0_C):
        '''Yield all intersection of two geodesic lines up to a limit.

           @return: Yield an L{Intersector5Tuple}C{(A, B, sAB, aAB, c)}
                    for each intersection found.

           @see: Methods L{All} for further details.
        '''
        for X in self.All(glA, glB, X0=X0, **aMaX0_sMaX0_C):
            yield self._In5T(glA, glB, X, X)

    def _Basic2(self, glA, glB, S, i=0):
        '''(INTERNAL) Get a basic solution.
        '''
        X = _copy(S)
        for _ in range(_TRIPS):
            S  = self._Spherical(glA, glB, X)
            X += S
            i += 1
            if X.c or S.L1() <= self._Tol:  # or isnan
                return self._Delto(X), i

        raise IntersectionError(Fmt.no_convergence(S.L1(), self._Tol))

    def _C(self, X, glA, glB, _C=False, _MM=False):
        # add the C{_C} items to C{X}, if requested.
        if _C:
            A = self._Position(glA, X.sA)
            B = self._Position(glB, X.sB)
            s, a = self._Inversa12(A, B)
            X.set_(latA=A.lat2, lonA=A.lon2,
                   latB=B.lat2, lonB=B.lon2)
            if _MM:  # in .Middle5
                X.set_(sMM=s, aMM=a)
            else:
                X.set_(sAB=s, aAB=a)
        return X

    def Closest(self, glA, glB, X0=_X000, **_C):
        '''Find the closest intersection of two geodesic lines.

           @arg glA: A geodesic line (L{Line<Intersector.Line>}).
           @arg glB: An other geodesic line (L{Line<Intersector.Line>}).
           @kwarg X0: Optional I{origin} for I{L1-closeness} (L{XDict}).
           @kwarg _C: If C{True}, include the lat-/longitudes C{latA},
                      C{lonA}, C{latB}, C{lonB} oon and distances C{sAB}
                      and C{aSB} between the intersections.

           @return: The intersection (L{XDict}) or C{None} if none found.

           @raise GeodesicError: Geodesic line B{C{glA}} or B{C{glB}}
                                 invalid, incompatible or ill-configured.

           @raise IntersectionError: No convergence.
        '''
        self._xLines(glA, glB)
        Q, d, S_, i = X0, INF, list(X0._nD1(self._D1)), 0
        while S_:
            X, i = self._Basic2(glA, glB, S_.pop(0), i)
            X    = X0._fixCoincident(X)
            if X.L1(Q) > self.Delta:  # X != Q
                d0 = X.L1(X0)
                if d0 < self._T1:
                    Q, d, q = X, d0, i
                    break
                if d0 < d or Q is X0:
                    Q, d, q = X, d0, i
                X._skip(S_, self._T2D1Delta)

        return None if Q is X0 else self._C(Q, glA, glB, **_C).set_(sX0=d, iteration=q)

    def Closest5(self, glA, glB, X0=_X000):
        '''Find the closest intersection of two geodesic lines.

           @return: An L{Intersector5Tuple}C{(A, B, sAB, aAB, c)}
                    or C{None} if none found.

           @see: Method L{Closest} for further details.
        '''
        X = self.Closest(glA, glB, X0=X0)
        return X if X is None else self._In5T(glA, glB, X, X)

    def _conjDist(self, gl, s, m12=0, M12=1, M21=1, semi=False):
        # Find semi-/conjugate point relative to s0 which is close to s1.
        # if semi:
        #     solve for M23 = 0 using dM23 / ds3 = - (1 - M23 * M32) / m23
        # else:
        #     solve for m23 = 0 using dm23 / ds3 = M32
        _S2, _abs, _1 = Fsum(s).fsum2_, fabs, _1_0
        for _ in range(_TRIPS):
            m13, M13, M31 = self._m12_M12_M21(gl, s)
            # see "Algorithms for geodesics", eqs. 31, 32, 33.
            m23 = m13 * M12
            M32 = M31 * M12
            if m12:  # PYCHOK no cover
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
            if _abs(d) <= self._Tol:
                break
        return s

    _gl3 = None

    @Property
    def _conjDist3s(self):
        gl, self._gl3, _D = self._gl3, None, self._cHalf
        return tuple(self._conjDist(gl, s) for s in (-_D, 0, _D))

    @_conjDist3s.setter  # PYCHOK setter!
    def _conjDist3(self, gl):
        # _XLines(gl, gl)
        self._gl3 = gl

    def _conjDist3Tt_(self, c, X0=_X000):
        for s in self._conjDist3s:
            T = XDict_(s, s * c, c)
            yield self._Delto(T), T.L1(X0)

    def _conjDist5(self, azi):
        gl   = self._Line(azi1=azi)
        s    = self._conjDist(gl, self._cHalf)
        X, _ = self._Basic2(gl, gl, XDict_(s * _0_5, -s * _1_5))
        return s, (X.L1() - s * _2_0), azi, X.sA, X.sB

    @Property_RO
    def Delta(self):
        '''Get the equality and tiling margin (C{meter}).
        '''
        return self._cHalf * _EPSr5  # ~15 Km WGS84

    def _Delto(self, X):
        # copy Delta into X, overriding X's default
        X._Delta = self.Delta  # NOT X.set_(self.Delta)
        return X

    @Property_RO
    def _EPS3R(self):
        return _EPS3 * self.R

    @Property_RO
    def _faPI_4(self):
        return (self.f + _2_0) * self.a * PI_4

    @Property_RO
    def _GeodesicLines(self):
        '''(INTERNAL) Get the C{Geodesic...Line} class(es).
        '''
        return type(self._Line()),

    def _In5T(self, glA, glB, S, X, k2=False, **_2X):
        # Return an intersection as C{Intersector5Tuple}.
        A = self._Position(glA, S.sA)
        B = self._Position(glB, S.sB)
        if k2:
            A.set_(k2=X.kA)
            B.set_(k2=X.kB)
        s, a = self._Inversa12(A, B)
        return Intersector5Tuple(A._2X(glA, **_2X),
                                 B._2X(glB, **_2X), s, a, X.c, iteration=X.iteration)

    def _Inverse(self, A, B):  # caps=Caps.STANDARD
        return self._g.Inverse(A.lat2, A.lon2, B.lat2, B.lon2)

    def Line(self, lat1, lon1, azi1_lat2, *lon2, **name):
        '''Return a geodesic line from this C{Intersector}'s geodesic, specified by
           two (goedetic) points or a (goedetic) point and an (initial) azimuth.

           @arg lat1: Latitude of the first point (C{degrees}).
           @arg lon1: Longitude of the first point (C{degrees}).
           @arg azi1_lat2: Azimuth at the first point (compass C{degrees}) if no
                          B{C{lon2}} argument is given, otherwise the latitude of
                          the second point (C{degrees}).
           @arg lon2: If given, the longitude of the second point (C{degrees}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: A line (from L{geodesic<Intersector.geodesic>}C{.Line} or
                    C{.InverseLine} method) with C{LINE_CAPS}.
        '''
        args = self._ll3z4ll(lat1, lon1, azi1_lat2, *lon2)
        gl   = self._g.InverseLine(*args, caps=Caps.LINE_CAPS) if lon2 else \
               self._g.Line(       *args, caps=Caps.LINE_CAPS)
        if name:
            gl.name= name
        return gl

    def _Line(self, lat1=0, lon1=0, azi1=0):
        return self._g.Line(lat1, lon1, azi1, caps=Caps.LINE_CAPS)

    def Middle(self, glA, glB, raiser=True, **_C):
        '''Get the mid-points on two geodesic line segments.

           @arg glA: A geodesic line (L{Line<Intersector.Line>}, 4-C{args}).
           @arg glB: An other geodesic line (L{Line<Intersector.Line>}, 4-C{args}).
           @kwarg raiser: If C{True}, check that B{C{glA}} and B{C{glB}} are a
                          4-C{args} L{Line<Intersector.Line>} or C{InverseLine}
                          (C{bool}).
           @kwarg _C: If C{True}, include the lat-/longitudes C{latA}, C{lonA},
                      C{latB}, C{lonB} of the mid-points and half-lengths C{sA}
                      and C{sB} in C{meter} of the respective line segments.

           @return: The mid-point and half-length of each segment (L{XDict}),
                    B{C{_C}} above.

           @raise GeodesicError: Geodesic line B{C{glA}} or B{C{glB}} invalid,
                                 incompatible, ill-configured or not a 4-C{args
                                 Line} or other C{InverseLine}.
        '''
        M, _, _ = self._middle3(glA, glB, raiser)
        return self._C(M, glA, glB, **_C) if _C else M

    def _middle3(self, glA, glB, raiser):  # in .All, .Segment
        # return segment length C{sA} and C{sB} and the
        # center C{X0} of rectangle [sA, sB]
        self._xLines(glA, glB, s13=raiser)  # need .Arc, .Distance
        sA = glA.Distance()
        sB = glB.Distance()
        X  = XDict_(sA * _0_5, sB * _0_5)
        # _ = X._outSide(sA, sB)
        return self._Delto(X), sA, sB

    def Middle5(self, glA, glB, raiser=True):
        '''Get the mid-points of two geodesic line segments and distances.

           @return: A L{Middle5Tuple}C{(A, B, sMM, aMM, c)}.

           @see: Method L{Middle} for further details.
        '''
        M, _, _ = self._middle3(glA, glB, raiser)
        M  = self._C(M, glA, glB, _C=True, _MM=True)
        A, B, s, a, c = self._In5T(glA, glB, M, M, _2X=_M_)
        return Middle5Tuple(self._illz2G(A, glA),
                            self._illz2G(B, glB), s, a, c)

    def _m12_M12_M21(self, gl, s):
        P = gl.Position(s, outmask=Caps._REDUCEDLENGTH_GEODESICSCALE)
        return P.m12, P.M12, P.M21

    def Next(self, glA, glB, eps1=None, **_C):  # PYCHOK no cover
        '''Yield the next intersection of two I{intersecting} geodesic lines.

           @arg glA: A geodesic line (L{Line<Intersector.Line>}).
           @arg glB: An other geodesic line (L{Line<Intersector.Line>}).
           @kwarg eps1: Optional margin for the L{euclid<pygeodesy.euclid>}ean
                        distance (C{degrees}) between the C{(lat1, lon1)} points
                        of both lines or C{None} for unchecked.
           @kwarg _C: If C{True}, include the lat-/longitudes C{latA}, C{lonA},
                      C{latB}, C{lonB} of and distances C{sAB} and C{aSB}
                      between the intersections.

           @return: The intersection (L{XDict}) or C{None} if none found.

           @raise GeodesicError: Geodesic line B{C{glA}} or B{C{glB}} invalid,
                                 incompatible, ill-configured or C{(lat1, lon1)}
                                 not B{C{eps1}}-equal.

           @raise IntersectionError: No convergence.

           @note: Offset C{X0} is implicit, zeros.
        '''
        self._xLines(glA, glB)
        if eps1 or _C:  # eps
            _C = self._xNext(glA, glB, eps1, **_C)

        X0, self._conjDist3s = _X000, glA  # reset Property
        Q, d, S_, i = _XINF, INF, list(X0._nD2(self._D2)), 0
        while S_:
            X, i = self._Basic2(glA, glB, S_.pop(0), i)
            X    = X0._fixCoincident(X)
            t    = X.L1(X0)  # == X.L1()
            c, z = X.c, (t <= self.Delta)  # X == X0
            if z:
                if not c:
                    continue
                Tt_ = self._conjDist3Tt_(c, X0)
            else:
                Tt_ = (X, t),

            for T, t in Tt_:
                if t < d or Q is _XINF:
                    Q, d, q = T, t, i
                i += 1

            for s in ((_1_1t if z else _1_0_1t)
                             if c else _0t):
                T = X
                if s and c:
                    s *= self._D2
                    T  = X + (s, s * c)  # NOT +=
                T._skip(S_, self._T2D2Delta)

        return None if Q is _XINF else self._C(Q, glA, glB, **_C).set_(sX0=d, iteration=q)

    def Next5(self, glA, glB, **eps1):  # PYCHOK no cover
        '''Yield the next intersection of two I{intersecting} geodesic lines.

           @return: An L{Intersector5Tuple}C{(A, B, sAB, aAB, c)} or C{None}
                    if none found.

           @see: Method L{Next} for further details.
        '''
        X = self.Next(glA, glB, **eps1)
        return X if X is None else self._In5T(glA, glB, X, X)

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
            sx, sAx, sBx = self._cHalf, _0_5, -_1_5
        return sx, zx, sAx, sBx

    def _polarB3(self, lats=False):  # PYCHOK no cover
        latx = _64_0
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
                t = (lat1 - lat0), (lat0 - lat2), (lat2 - lat1)
                d = _d(t, s2, s1, s0) * _2_0
                if not d:  # or isnan(d)
                    break
                lat = _d(t, (lat1 + lat0) * s2,
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
            sx  = self._cHalf
        return sx, latx, lat

    def _polarDist2(self, lat1, lat2=False):
        gl = self._Line(lat1=lat1)
        s  = self._conjDist(gl, self._faPI_4, semi=True)
        if lat2:
            lat1 = gl.Position(s, outmask=Caps.LATITUDE).lat2
        return s, lat1

    def _Position(self, gl, s):
        return gl.Position(s, outmask=Caps._STD_LINE)

    def Segment(self, glA, glB, proven=None, raiser=True, **_C):
        '''Find the intersection between two geodesic line segments.

           @kwarg proven: Conjecture is that whenever two geodesic line
                          segments intersect, the intersection is the
                          one closest to the mid-points of segments.
                          If so, use C{B{proven}=True}, otherwise find
                          intersections on the segments and specify
                          C{B{proven}=None} to return the first or
                          C{B{proven}=False} the closest (C{bool}).
           @kwarg raiser: If C{True}, check that B{C{glA}} and B{C{glB}}
                          are a 4-C{args} L{Line<Intersector.Line>} or
                          C{InverseLine} (C{bool}).
           @kwarg _C: If C{True}, include the lat-/longitudes C{latA},
                      C{lonA}, C{latB}, C{lonB} of and distances C{sAB}
                      and C{aSB} between the intersections.

           @return: The intersection of the segments (L{XDict}) with
                    indicators C{kA}, C{kB} and C{k} set or if no
                    intersection is found, C{None}.

           @raise GeodesicError: Geodesic line B{C{glA}} or B{C{glB}}
                                 invalid, incompatible, ill-configured or
                                 not an C{InverseLine} or 4-C{args Line}.

           @raise IntersectionError: No convergence.

           @see: Method L{Middle<Intersector.Middle>} for further details.
        '''
        X0, sA, sB = self._middle3(glA, glB, raiser)
        Q = self.Closest(glA, glB, X0)  # to X0
        if Q is not None:
            if Q.c:  # anti-/parallel
                Q._fixSegment(sA, sB)
            # are rectangle [sA, sB] corners further from X0 than Q?
            d0 = X0.L1(Q)
            if Q._outSide(sA, sB) and d0 <= X0.L1() and not proven:
                i = Q.iteration
                for T in Q._corners(sA, sB, self._T2):
                    X, i = self._Basic2(glA, glB, T, i)
                    X    = T._fixCoincident(X)
                    if not X._outSide(sA, sB):
                        d = X0.L1(X)
                        if d < d0 or proven is None:
                            Q, d0 = X, d
                            if proven is None:
                                break
                Q.set_(iteration=i)

            Q = self._C(Q, glA, glB, **_C).set_(sX0=d0)
        return Q

    def Segment5(self, glA, glB, **proven_raiser):
        '''Find the intersection between two geodesic line segments.

           @return: An L{Intersector5Tuple}C{(A, B, sAB, aAB, c)}
                    or C{None} if none found.

           @see: Method L{Segment} for further details.
        '''
        X = self.Segment(glA, glB, **proven_raiser)
        return X if X is None else self._In5T(glA, glB, X, X, k2=True)

    def _Spherical(self, glA, glB, S):
        '''(INTERNAL) Get solution based from a spherical triangle.
        '''
        # threshold for coincident geodesics/intersections ~4.3 nm WGS84.
        A = self._Position(glA, S.sA)
        B = self._Position(glB, S.sB)
        D = self._Inverse(A, B)

        a, da = _diff182(A.azi2, D.azi1)  # interior angle at A
        b, db = _diff182(B.azi2, D.azi2)  # exterior angle at B
        c, dc = _diff182(a, b)
        if fsum1_(dc, db, -da, c) < 0:  # inverted triangle
            a, da = -a, -da
            b, db = -b, -db
        sa, ca = _sincos2de(a, da)
        sb, cb = _sincos2de(b, db)

        e, z, _abs = _EPS3, D.s12, fabs
        if _abs(z) <= self._EPS3R:  # XXX z <= ...
            sA = sB = 0  # at intersection
            c  = 1 if _abs(sa - sb) <= e and _abs(ca - cb) <= e else (
                -1 if _abs(sa + sb) <= e and _abs(ca + cb) <= e else 0)
        elif _abs(sa) <= e and _abs(sb) <= e:  # coincident
            sA =  ca * z * _0_5  # choose mid-point
            sB = -cb * z * _0_5
            c  =  1 if (ca * cb) > 0 else -1
            # alt1: sA =  ca * z; sB = 0
            # alt2: sB = -cb * z; sA = 0
        else:  # general case
            sz, cz = sincos2(z / self.R)
            # [SKIP: Divide args by |sz| to avoid possible underflow
            # in {sa, sb} * sz; this is probably not necessary].
            # Definitely need to treat sz < 0 (z > PI*R) correctly in
            # order to avoid some convergence failures in _Basic2.
            sA = atan2(sb * sz,  sb * ca * cz - sa * cb) * self.R
            sB = atan2(sa * sz, -sa * cb * cz + sb * ca) * self.R
            c  = 0
        return XDict_(sA, sB, c)  # no ._Delto

    @Property_RO
    def _T2D1Delta(self):
        return self._T2d3Delta(self._D1)

    @Property_RO
    def _T2D2Delta(self):
        return self._T2d3Delta(self._D2)

    def _T2d3Delta(self, d3):
        return self._T2 - d3 - self.Delta

    @Property_RO
    def _Tol(self):  # convergence tolerance
        return self._cHalf * _EPSjam

    def toStr(self, **prec_sep_name):  # PYCHOK signature
        '''Return this C{Intersector} as string.

           @see: L{Ellipsoid.toStr<pygeodesy.ellipsoids.Ellipsoid.toStr>}
                 for further details.

           @return: C{Intersector} (C{str}).
        '''
        return self._instr(props=(Intersector.geodesic,), **prec_sep_name)

    def _xLines(self, glA, glB, s13=False):
        # check two geodesic lines vs this geodesic
        C, gls = Caps.LINE_CAPS, dict(glA=glA, glB=glB)
        _xinstanceof(*self._GeodesicLines, **gls)
        for n, gl in gls.items():
            try:
                _xgeodesics(gl.geodesic, self.geodesic)
                if s13 and not isfinite(gl.s13):  # or not gl.caps & Caps.DISTANCE_IN
                    t = gl.geodesic.InverseLine.__name__
                    raise TypeError(_not_(_an(t)))
                c = gl.caps & C
                if c != C:  # not gl.caps_(C)
                    c, C, x = map1(bin, c, C, _xor(c, C))
                    x = _SPACE_(_xor.__name__, repr(x))[1:]
                    raise GeodesicError(caps=c, LINE_CAPS=C, txt=x)
            except Exception as x:
                raise GeodesicError(n, gl, cause=x)


class Intersect7Tuple(_NamedTuple):
    '''7-Tuple C{(A, B, sAB, aAB, c, kA, kB)} with C{A} and C{B} each
       a C{LatLon} or L{LatLon4Tuple}C{(lat, lon, height, datum)} of
       the intersection on each geodesic line, the distance C{sAB} in
       in C{meter} and angular distance C{aAB} in C{degrees} between
       C{A} and C{B}, coincidence indicator C{c} and segment indicators
       C{kA} and C{kB} all C{int}, see L{XDict} and method U{intersect7
       <_IntersectBase.intersect7>}.
    '''
    _Names_ = (_A_,   _B_,  _sAB_, _aAB_,   _c_, 'kA', 'kB')
    _Units_ = (_Pass, _Pass, Meter, Degrees, Int, Int,  Int)


class Intersectool5Tuple(_NamedTuple):
    '''5-Tuple C{(A, B, sAB, aAB, c)} with C{A} and C{B} the C{Position}
       of the intersection on each geodesic line, the distance C{sAB}
       between C{A} and C{B} in C{meter}, the angular distance C{aAB} in
       C{degrees} and coincidence indicator C{c} (C{int}), see L{XDict}.

       @note: C{A} and C{B} are each a C{GDict} with C{lat1}, C{lon1} and
              C{azi1} or C{lat2}, C{lon2} from the geodesic line C{glA}
              respectively C{glB} and the intersection location in C{latX},
              C{lonX}, distance C{s1X} in C{meter} and angular distance
              C{a1M} in C{degrees} and the segment indicator C{kX}.  See
              L{XDict} for more details.
    '''
    _Names_ = Intersect7Tuple._Names_[:5]
    _Units_ = Intersect7Tuple._Units_[:5]


class Intersector5Tuple(Intersectool5Tuple):
    '''5-Tuple C{(A, B, sAB, aAB, c)} with C{A} and C{B} the C{Position}
       of the intersection on each geodesic line, the distance C{sAB}
       between C{A} and C{B} in C{meter}, angular distance C{aAB} in
       C{degrees} and coincidence indicator C{c} (C{int}), see L{XDict}.

       @note: C{A} and C{B} are each a C{GeodesicLine...Position} for
              C{outmask=Caps.STANDARD} with the intersection location in
              C{latX}, C{lonX}, azimuth in C{aziX}, distance C{s1X} in
              C{meter} and angular distance C{a1X} in C{degrees} and the
              segment indicator C{kX}.  See L{XDict} for more details.
    '''
    _Names_ = Intersectool5Tuple._Names_


class Middle5Tuple(Intersectool5Tuple):
    '''5-Tuple C{(A, B, sMM, aMM, c)} with C{A} and C{B} the I{line segments}
       including the mid-point location in C{latM}, C{lonM}, distance C{s1M}
       in C{meter} and angular distance C{a1M} in C{degrees}, the distance
       between both mid-points C{sMM} in C{meter} and angular distance C{aMM}
       in C{degrees} and coincidence indicator C{c} (C{int}).  See L{XDict}
       for more details.
    '''
    _Names_ = (_A_, _B_, 'sMM', 'aMM', _c_)


class _List(list):

    _Delta = 0  # equality margin

    def __init__(self, Delta):
        self._Delta = Delta
#       list.__init__(self)

    def __contains__(self, other):
        # handle C{if X in this: ...}
        a,  b  = other.sA, other.sB
        D, _D1 = self._Delta, _L1
        for X in self:
            if _D1(X.sA - a, X.sB - b) <= D:
                return True
        return False

    def addend(self, X, *d0_i):
        # append an item, updated
        if d0_i:
            d0, i = d0_i
            X.set_(sX0=d0, iteration=i)
        self.append(X)
        return X.sX0

    def sorter(self, sMaX0, dot_C, glA, glB, **_C):
        # trim and sort the X items

        def _key(X):
            return X.sX0  # rank of X

        t = (X for X in self if X.sX0 <= sMaX0)
        for X in sorted(t, key=_key):
            yield dot_C(X, glA, glB, **_C) if _C else X


def _L1(a, b):
    '''(INTERNAL) Return the I{L1} distance.
    '''
    return fabs(a) + fabs(b)


__all__ += _ALL_DOCS(_IntersectBase)

if __name__ == '__main__':  # MCCABE 14

    from pygeodesy import printf
    __help_ = '--help'

    def _main(args):

        from pygeodesy import GeodesicExact
        from pygeodesy.internals import _plural, _usage
        from pygeodesy.interns import _COLONSPACE_, _DOT_, _EQUAL_, \
                                      _i_, _m_, _n_, _version_, _X_
        import re

        class XY0(Float):
            pass

        def _opts(_h):  # for _usage()
            ll4  = ' latA1 lonA1'
            ll4 += ll4.replace('1', '2')
            ll4 += ll4.replace(_A_, _B_)
            llz  = _SPACE_(NN, _latA_, _lonA_, 'aziA')
            llz2 =  llz + llz.replace(_A_, _B_)
            return dict(opts='-Verbose|V--version|v--help|h--Tool|T--Check|C-R <meter>-',
                        alts=((_c_ + llz2),
                              (_i_ + ll4),
                              (_m_ + ll4),
                              (_n_ + llz  + ' aziB'),
                              ('o' + llz2 + ' x0 y0')),
                        help=_h if isinstance(_h, str) else NN)

        def _starts(Opt, arg):
            return arg == Opt[1:3] or (len(arg) > 2 and Opt.startswith(arg))

        _isopt = re.compile('^[-]+[a-z]*$', flags=re.IGNORECASE).match

        I  =  Intersector(GeodesicExact())  # PYCHOK I
        M  =  m = _R =  None
        _T = _V = _h = _C = False

        while args and _isopt(args[0]):
            arg = args.pop(0)
            if arg == _c__:
                M, m = I.Closest, 6  # latA lonA aziA  latB lonB aziB
            elif _starts('--Check', arg):
                _C = True
            elif _starts(__help_, arg):
                _h = args[0] if args and _isopt(args[0]) else True
            elif arg == _i__:
                M, m = I.Segment, 8  # latA1 lonA1  latA2 lonA2  latB1 lonB1  latB2 lonB2
            elif arg == '-m':
                M, m = I.Middle,  8  # latA1 lonA1  latA2 lonA2  latB1 lonB1  latB2 lonB2
                _R = None  # zap -R
            elif arg == _n__:
                M, m = I.Next, 4  # latA lonA aziA  aziB
            elif arg == _o__:
                M, m = I.Closest, 8  # latA lonA aziA  latB lonB aziB  x0 y0
            elif arg == _R__ and args:
                _R = args.pop(0)
            elif _starts('--Tool', arg):
                I = Intersectool()  # PYCHOK I
                if _V:
                    I.verbose = True
                if not _Xables.X_OK(I.IntersectTool):
                    I.IntersectTool = _Xables.IntersectTool(_Xables.bin_)
                elif _V:
                    _ = I.version
                M, _T = None, True
            elif _starts('--Verbose', arg):
                _V = True
                if _T:
                    I.verbose = True
            elif _starts('--version', arg):
                printf(_COLONSPACE_(*((_version_, I.version) if _T else
                                      (__version__, repr(I)))))
            else:
                raise ValueError('invalid option %r' % (arg,))

        if _h or M is None:
            printf(_usage(__file__, **_opts(_h)), nl=1)
        else:
            n = len(args)
            if n < m:
                n = _plural('only %s arg' % (n,), n) if n else 'no args'
                raise ValueError('%s, need %s' % (n, m))
            args[:] = args[:m]

            kwds = dict(_C=True) if _C else {}
            if M == I.Next:  # -n
                # get latA lonA aziA latA lonA aziB
                args[3:] = args[:2] + args[3:4]
            elif M == I.Closest and m > 6:  # -o
                y0 = Meter(y0=args.pop())
                x0 = Meter(x0=args.pop())
                kwds.update(X0=XDict_(x0, y0))
            if _R:
                m = Meter_(_R, name=_R__, low=0)
                kwds.update(sMaX0=m)
                M = I.All

            n   = len(args) // 2
            glA = I.Line(*args[:n])
            glB = I.Line(*args[n:])

            m = _DOT_(I.__class__.__name__, M.__name__)
            if _V:
                X = _SPACE_(_X_, _EQUAL_, m)
                printf(unstr(X, glA, glB, **kwds))

            X = M(glA, glB, **kwds)
            if X is None or isinstance(X, (XDict, tuple)):
                printf(_COLONSPACE_(m, repr(X)))
            else:
                for i, X in enumerate(X):
                    printf(_COLONSPACE_(Fmt.INDEX(m, i), repr(X)))

    def _examples():

        from pygeodesy.internals import _usage_argv

        s = _SPACE_(*_usage_argv(__file__))
        for t in ('-h', '-h -n',
                  '-c 0 0 45  40 10 135',
                  '-C -c 0 0 45  40 10 135',
                  '-T -R 2.6e7 -c 0 0 45  40 10 135',
                  '-c 50 -4 -147.7  0 0 90',
                  '-C -c 50 -4 -147.7  0 0 90',
                  '# % echo 0 0  10 10  50 -4  50S 4W | IntersectTool -i  -p 0  -C',
                  '# -631414 5988887 0 -3',
                  '# -4.05187 -4.00000 -4.05187 -4.00000 0',
                  '-m 0 0  10 10  50 -4  50S 4W',
                  '-C -m 0 0  10 10  50 -4  50S 4W',
                  '-i 0 0  10 10  50 -4  50S 4W',
                  '-T -i 0 0  10 10  50 -4  50S 4W',
                  '-C -i 0 0  10 10  50 -4  50S 4W',
                  '-T -C -i 0 0  10 10  50 -4  50S 4W',
                  '-V -T -i 0 0  10 10  50 -4  -50 -4',
                  '-C -R 4e7 -c 50 -4 -147.7  0 0 90',
                  '-T -C -R 4e7 -c 50 -4 -147.7  0 0 90',
                  '-R 4e7 -i 0 0  10 10  50 -4  -50 -4',
                  '-T -R 4e7 -i 0 0  10 10  50 -4  -50 -4'):
            if t.startswith(_HASH_):
                printf(t, nl=int(t[2] == '%'))
            else:
                printf(_SPACE_(_HASH_, s, t), nl=1)
                argv[1:] = t = t.split()
                _main(t)

    from sys import argv, stderr
    try:
        if len(argv) == 2 and argv[1] == __help_:
            _examples()
        else:
            _main(argv[1:])

    except Exception as x:
        x = _SPACE_(x, NN, _HASH_, *argv)
        printf(x, file=stderr, nl=1)
        if '-V' in x or _MODS.errors.exception_chaining():
            raise
        exit(1)

# % env PYGEODESY_INTERSECTTOOL=... python3 -m pygeodesy.geodesici --help

# % python3 -m pygeodesy.geodesici -h
#
# usage: python3 -m ....pygeodesy.geodesici [--Verbose | -V] [--version | -v] [--help | -h] [--Tool | -T] [--Check | -C] [-R meter]
#                                           [-c latA lonA aziA latB lonB aziB |
#                                            -i latA1 lonA1 latA2 lonA2 latB1 lonB1 latB2 lonB2 |
#                                            -m latA1 lonA1 latA2 lonA2 latB1 lonB1 latB2 lonB2 |
#                                            -n latA lonA aziA aziB |
#                                            -o latA lonA aziA latB lonB aziB x0 y0]

# % python3 -m ....pygeodesy.geodesici -h -n
#
# usage: python3 -m ....pygeodesy.geodesici -n latA lonA aziA aziB

# % python3 -m ....pygeodesy.geodesici -c 0 0 45  40 10 135
# Intersector.Closest: XDict(c=0, sA=3862290.547855, sB=2339969.547699, sX0=6202260.095554)

# % python3 -m ....pygeodesy.geodesici -C -c 0 0 45  40 10 135
# Intersector.Closest: XDict(aAB=0.0, c=0, latA=23.875306, latB=23.875306, lonA=26.094096, lonB=26.094096, sA=3862290.547855, sAB=0.0, sB=2339969.547699, sX0=6202260.095554)

# % env PYGEODESY_INTERSECTTOOL=...python3 -m ....pygeodesy.geodesici -T -R 2.6e7 -c 0 0 45  40 10 135
# Intersectool.All[0]: XDict(c=0, sA=3862290.547855, sB=2339969.547699, sX0=6202260.095554)

# % python3 -m ....pygeodesy.geodesici -c 50 -4 -147.7  0 0 90
# Intersector.Closest: XDict(c=0, sA=6058048.653081, sB=-3311252.995823, sX0=9369301.648903)

# % python3 -m ....pygeodesy.geodesici -C -c 50 -4 -147.7  0 0 90
# Intersector.Closest: XDict(aAB=0.0, c=0, latA=0.0, latB=-0.0, lonA=-29.745492, lonB=-29.745492, sA=6058048.653081, sAB=0.0, sB=-3311252.995823, sX0=9369301.648903)

# % echo 0 0  10 10  50 -4  50S 4W | IntersectTool -i  -p 0  -C
# -631414 5988887 0 -3
# -4.05187 -4.00000 -4.05187 -4.00000 0

# % python3 -m ....pygeodesy.geodesici -m 0 0  10 10  50 -4  50S 4W
# Intersector.Middle: XDict(c=0, sA=782554.549609, sB=5536835.161499, sX0=0.0)

# % python3 -m ....pygeodesy.geodesici -C -m 0 0  10 10  50 -4  50S 4W
# Intersector.Middle: XDict(aAB=10.262308, c=0, latA=5.019509, latB=0.036282, lonA=4.961883, lonB=-4.0, sA=782554.549609, sAB=1138574.546746, sB=5536835.161499, sX0=0.0)

# % python3 -m ....pygeodesy.geodesici -i 0 0  10 10  50 -4  50S 4W
# Intersector.Segment: XDict(c=0, k=-3, kA=-1, kB=0, sA=-631414.26877, sB=5988887.278435, sX0=1866020.935315)

# % env PYGEODESY_INTERSECTTOOL=... python3 -m ....pygeodesy.geodesici -T -i 0 0  10 10  50 -4  50S 4W
# Intersectool.Segment: XDict(c=0, k=-3, kA=-1, kB=0, sA=-631414.26877, sB=5988887.278435)

# % python3 -m ....pygeodesy.geodesici -C -i 0 0  10 10  50 -4  50S 4W
# Intersector.Segment: XDict(aAB=0.0, c=0, k=-3, kA=-1, kB=0, latA=-4.051871, latB=-4.051871, lonA=-4.0, lonB=-4.0, sA=-631414.26877, sAB=0.0, sB=5988887.278435, sX0=1866020.935315)

# % env PYGEODESY_INTERSECTTOOL=... python3 -m ....pygeodesy.geodesici -T -C -i 0 0  10 10  50 -4  50S 4W
# Intersectool.Segment: XDict(c=0, k=-3, kA=-1, kB=0, latA=-4.051871, latB=-4.051871, lonA=-4.0, lonB=-4.0, sA=-631414.26877, sAB=0.0, sB=5988887.278435)

# % env PYGEODESY_INTERSECTTOOL=... python3 -m ....pygeodesy.geodesici -V -T -i 0 0  10 10  50 -4  -50 -4
# Intersectool@1: /opt/local/bin/IntersectTool --version (invoke)
# Intersectool@1: '/opt/local/bin/IntersectTool: GeographicLib version 2.3' (0)
# Intersectool@1: /opt/local/bin/IntersectTool: GeographicLib version 2.3 (0)
# X = Intersectool.Segment(GDict(lat1=0.0, lat2=10.0, lon1=0.0, lon2=10.0), GDict(lat1=50.0, lat2=-50.0, lon1=-4.0, lon2=-4.0))
# Intersectool@2: /opt/local/bin/IntersectTool -E -p 10 -i \ 0.0 0.0 10.0 10.0 50.0 -4.0 -50.0 -4.0 (Segment)
# Intersectool@2: '-631414.2687702414 5988887.2784352796 0 -3' (0)
# Intersectool@2: sA=-631414.2687702414, sB=5988887.2784352796, c=0, k=-3 (0)
# Intersectool.Segment: XDict(c=0, k=-3, kA=-1, kB=0, sA=-631414.26877, sB=5988887.278435)

# % python3 -m ....pygeodesy.geodesici -C -R 4e7 -c 50 -4 -147.7  0 0 90
# Intersector.All[0]: XDict(aAB=0.0, c=0, latA=0.0, latB=-0.0, lonA=-29.745492, lonB=-29.745492, sA=6058048.653081, sAB=0.0, sB=-3311252.995823, sX0=9369301.648903)
# Intersector.All[1]: XDict(aAB=0.0, c=0, latA=0.0, latB=0.0, lonA=150.046964, lonB=150.046964, sA=-13941907.021445, sAB=0.0, sB=16703151.659744, sX0=30645058.681189)
# Intersector.All[2]: XDict(aAB=0.0, c=0, latA=-0.0, latB=-0.0, lonA=-30.16058, lonB=-30.16058, sA=-33941862.69597, sAB=0.0, sB=-3357460.370268, sX0=37299323.066238)
# Intersector.All[3]: XDict(aAB=0.0, c=0, latA=-0.0, latB=0.0, lonA=150.046964, lonB=150.046964, sA=-13941907.021445, sAB=0.0, sB=-23371865.025835, sX0=37313772.047279)

# % env PYGEODESY_INTERSECTTOOL=... python3 -m ....pygeodesy.geodesici -T -C -R 4e7 -c 50 -4 -147.7  0 0 90
# Intersectool.All[0]: XDict(c=0, latA=-0.0, latB=-0.0, lonA=-29.745492, lonB=-29.745492, sA=6058048.653081, sAB=0.0, sB=-3311252.995823, sX0=9369301.648903)
# Intersectool.All[1]: XDict(c=0, latA=0.0, latB=0.0, lonA=150.046964, lonB=150.046964, sA=-13941907.021445, sAB=0.0, sB=16703151.659744, sX0=30645058.681189)
# Intersectool.All[2]: XDict(c=0, latA=-0.0, latB=-0.0, lonA=-30.16058, lonB=-30.16058, sA=-33941862.69597, sAB=0.0, sB=-3357460.370268, sX0=37299323.066238)
# Intersectool.All[3]: XDict(c=0, latA=-0.0, latB=0.0, lonA=150.046964, lonB=150.046964, sA=-13941907.021445, sAB=0.0, sB=-23371865.025835, sX0=37313772.047279)

# % python3 -m ....pygeodesy.geodesici -R 4e7 -i 0 0  10 10  50 -4  -50 -4
# Intersector.All[0]: XDict(c=0, sA=-631414.26877, sB=5988887.278435, sX0=1866020.935315)
# Intersector.All[1]: XDict(c=0, sA=19422725.117572, sB=-14062417.105648, sX0=38239422.83511)
# Intersector.All[2]: XDict(c=0, sA=19422725.117572, sB=25945445.811603, sX0=39048781.218067)
# Intersector.All[3]: XDict(c=0, sA=39476927.464575, sB=5894074.699478, sX0=39051612.452944)

# % env PYGEODESY_INTERSECTTOOL=... python3 -m ....pygeodesy.geodesici -T -R 4e7 -i 0 0  10 10  50 -4  -50 -4
# Intersectool.All[0]: XDict(c=0, sA=-631414.26877, sB=5988887.278435, sX0=1862009.05513)
# Intersectool.All[1]: XDict(c=0, sA=19422725.117572, sB=-14062417.105648, sX0=38243434.715295)
# Intersectool.All[2]: XDict(c=0, sA=19422725.117572, sB=25945445.811603, sX0=39044769.337882)
# Intersectool.All[3]: XDict(c=0, sA=39476927.464575, sB=5894074.699478, sX0=39047600.57276)


# **) MIT License
#
# Copyright (C) 2024-2025 -- mrJean1 at Gmail -- All Rights Reserved.
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
