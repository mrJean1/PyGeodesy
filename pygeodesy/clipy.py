
# -*- coding: utf-8 -*-

u'''Clip a path or polygon of C{LatLon} points against a rectangular box or
an arbitrary (convex) region.

Box clip functions L{clipCS4} I{Cohen-Sutherland} and L{clipLB6} I{Liang-Barsky},
region clip functions L{clipFHP4} I{Foster-Hormann-Popa}, L{clipGH4}
I{Greiner-Hormann} and L{clipSH} and L{clipSH3} I{Sutherland-Hodgeman}.
.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

# from pygeodesy.basics import len2  # from .fmath
from pygeodesy.constants import EPS, _0_0, _1_0
from pygeodesy.errors import _AssertionError, ClipError, PointsError
from pygeodesy.fmath import fabs, len2
from pygeodesy.fsums import fsumf_, Property_RO
from pygeodesy.interns import NN, _clipid_, _convex_, _DOT_, _end_, _few_, \
                             _fi_, _height_, _i_, _invalid_, _j_, _lat_, \
                             _lon_, _near_, _not_, _points_, _start_, _too_
from pygeodesy.iters import _imdex2, points2
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import _Named, _NamedTuple, _Pass
from pygeodesy.points import areaOf, boundsOf, isconvex_, LatLon_
# from pygeodesy.props import Property_RO  # from .fsums
from pygeodesy.units import Bool, FIx, HeightX, Lat, Lon, Number_

# from math import fabs  # from .fmath

__all__ = _ALL_LAZY.clipy
__version__ = '23.05.15'

_fj_       = 'fj'
_original_ = 'original'


def _box4(lowerleft, upperight, name):
    '''(INTERNAL) Get the clip box edges.

       @see: Class C{_Box} in .ellipsoidalBaseDI.py.
    '''
    try:
        yb, yt = lowerleft.lat, upperight.lat
        xl, xr = lowerleft.lon, upperight.lon
        if xl > xr or yb > yt:
            raise ValueError(_invalid_)
    except (AttributeError, TypeError, ValueError) as x:
        raise ClipError(name, 2, (lowerleft, upperight), cause=x)
    return xl, yb, xr, yt


def _4corners(corners):
    '''(INTERNAL) Clip region or box.
    '''
    n, cs = len2(corners)
    if n == 2:  # make a box
        yb, xl, yt, xr = boundsOf(cs, wrap=False)
        cs = (LatLon_(yb, xl), LatLon_(yt, xl),
              LatLon_(yt, xr), LatLon_(yb, xr))
    return cs


def _eq(p1, p2, eps=EPS):
    '''(INTERNAL) Check for near-equal points.
    '''
    return not _neq(p1, p2, eps)


def _neq(p1, p2, eps=EPS):
    '''(INTERNAL) Check for not near-equal points.
    '''
    return fabs(p1.lat - p2.lat) > eps or \
           fabs(p1.lon - p2.lon) > eps


def _pts2(points, closed, inull):
    '''(INTERNAL) Get the points to clip as a list.
    '''
    if closed and inull:
        n, pts = len2(points)
        # only remove the final, closing point
        if n > 1 and _eq(pts[n-1], pts[0]):
            n -= 1
            pts = pts[:n]
        if n < 2:
            raise PointsError(points=n, txt=_too_(_few_))
    else:
        n, pts = points2(points, closed=closed)
    return n, list(pts)


class _CS(_Named):
    '''(INTERNAL) Cohen-Sutherland line clipping.
    '''
    # single-bit clip codes
    _IN =  0  # inside clip box
    _XR =  1  # right of upperight.lon
    _XL =  2  # left of lowerleft.lon
    _YT =  4  # above upperight.lat
    _YB =  8  # below lowerleft.lat

    _dx = _0_0  # pts edge delta lon
    _dy = _0_0  # pts edge delta lat
    _x1 = _0_0  # pts corner
    _y1 = _0_0  # pts corner

    _xr = _0_0  # clip box upperight.lon
    _xl = _0_0  # clip box lowerleft.lon
    _yt = _0_0  # clip box upperight.lat
    _yb = _0_0  # clip box lowerleft.lat

    def __init__(self, lowerleft, upperight, name=__name__):
        self._xl, self._yb, \
        self._xr, self._yt = _box4(lowerleft, upperight, name)
        self.name = name

#   def clip4(self, p, c):  # clip point p for code c
#       if c & _CS._YB:
#           return self.lon4(p, self._yb)
#       elif c & _CS._YT:
#           return self.lon4(p, self._yt)
#       elif c & _CS._XL:
#           return self.lat4(p, self._xl)
#       elif c & _CS._XR:
#           return self.lat4(p, self._xr)
#       # should never get here
#       raise _AssertionError(self._DOT_(self.clip4.__name__))

    def code4(self, p):  # compute code for point p
        if p.lat < self._yb:
            c, m, b = _CS._YB, self.lon4, self._yb
        elif p.lat > self._yt:
            c, m, b = _CS._YT, self.lon4, self._yt
        else:
            c, m, b = _CS._IN, self.nop4, None
        if p.lon < self._xl:
            c |= _CS._XL
            m, b = self.lat4, self._xl
        elif p.lon > self._xr:
            c |= _CS._XR
            m, b = self.lat4, self._xr
        return c, m, b, p

    def edge(self, p1, p2):  # set edge p1 to p2
        self._y1, self._dy = p1.lat, float(p2.lat - p1.lat)
        self._x1, self._dx = p1.lon, float(p2.lon - p1.lon)
        return fabs(self._dx) > EPS or fabs(self._dy) > EPS

    def lat4(self, x, p):  # new lat and code at lon x
        y = self._y1 + self._dy * float(x - self._x1) / self._dx
        if y < self._yb:  # still outside
            return _CS._YB, self.lon4, self._yb, p
        elif y > self._yt:  # still outside
            return _CS._YT, self.lon4, self._yt, p
        else:  # inside
            return _CS._IN, self.nop4, None, p.classof(y, x)

    def lon4(self, y, p):  # new lon and code at lat y
        x = self._x1 + self._dx * float(y - self._y1) / self._dy
        if x < self._xl:  # still outside
            return _CS._XL, self.lat4, self._xl, p
        elif x > self._xr:  # still outside
            return _CS._XR, self.lat4, self._xr, p
        else:  # inside
            return _CS._IN, self.nop4, None, p.classof(y, x)

    def nop4(self, b, p):  # PYCHOK no cover
        if p:  # should never get here
            raise _AssertionError(self._DOT_(self.nop4.__name__))
        return _CS._IN, self.nop4, b, p


class ClipCS4Tuple(_NamedTuple):
    '''4-Tuple C{(start, end, i, j)} for each edge of a I{clipped}
       path with the C{start} and C{end} points (C{LatLon}) of the
       portion of the edge inside or on the clip box and the indices
       C{i} and C{j} (C{int}) of the edge start and end points in
       the original path.
    '''
    _Names_ = (_start_, _end_, _i_,     _j_)
    _Units_ = (_Pass,   _Pass,  Number_, Number_)


def clipCS4(points, lowerleft, upperight, closed=False, inull=False):
    '''Clip a path against a rectangular clip box using the U{Cohen-Sutherland
       <https://WikiPedia.org/wiki/Cohen-Sutherland_algorithm>} algorithm.

       @arg points: The points (C{LatLon}[]).
       @arg lowerleft: Bottom-left corner of the clip box (C{LatLon}).
       @arg upperight: Top-right corner of the clip box (C{LatLon}).
       @kwarg closed: Optionally, close the path (C{bool}).
       @kwarg inull: Optionally, retain null edges if inside (C{bool}).

       @return: Yield a L{ClipCS4Tuple}C{(start, end, i, j)} for each
                edge of the I{clipped} path.

       @raise ClipError: The B{C{lowerleft}} and B{C{upperight}} corners
                         specify an invalid clip box.

       @raise PointsError: Insufficient number of B{C{points}}.
    '''
    T4 =  ClipCS4Tuple
    cs = _CS(lowerleft, upperight, name=clipCS4.__name__)
    n, pts = _pts2(points, closed, inull)

    i, m = _imdex2(closed, n)
    cmbp =  cs.code4(pts[i])
    for j in range(m, n):
        c1, m1, b1, p1 = cmbp
        c2, m2, b2, p2 = cmbp = cs.code4(pts[j])
        if c1 & c2:  # edge outside
            pass
        elif cs.edge(p1, p2):
            for _ in range(5):
                if c1:  # clip p1
                    c1, m1, b1, p1 = m1(b1, p1)
                elif c2:  # clip p2
                    c2, m2, b2, p2 = m2(b2, p2)
                else:  # inside
                    if inull or _neq(p1, p2):
                        yield T4(p1, p2, i, j)
                    break
                if c1 & c2:  # edge outside
                    break
            else:  # PYCHOK no cover
                raise _AssertionError(_DOT_(cs.name, 'for_else'))

        elif inull and not c1:  # null edge
            yield T4(p1, p1, i, j)
        elif inull and not c2:
            yield T4(p2, p2, i, j)

        i = j


class ClipFHP4Tuple(_NamedTuple):
    '''4-Tuple C{(lat, lon, height, clipid)} for each point of the
       L{clipFHP4} result with the C{lat}-, C{lon}gitude, C{height}
       and C{clipid} of the polygon or clip.

       @note: The C{height} is a L{HeightX} instance if this point is
              an intersection, otherwise a L{Height} or C{int(0)}.
    '''
    _Names_ = (_lat_, _lon_, _height_, _clipid_)
    _Units_ = ( Lat,   Lon,  _Pass,     Number_)

    @Property_RO
    def isintersection(self):
        '''Is this an intersection?
        '''
        return isinstance(self.height, HeightX)

    @Property_RO
    def ispoint(self):
        '''Is this an original (polygon) point?
        '''
        return not self.isintersection


def clipFHP4(points, corners, closed=False, inull=False, raiser=False, eps=EPS):
    '''Clip one or more polygons against a clip region or box using U{Forster-Hormann-Popa
       <https://www.ScienceDirect.com/science/article/pii/S259014861930007X>}'s C++
       implementation transcoded to pure Python.

       @arg points: The polygon points and clips (C{LatLon}[]).
       @arg corners: Three or more points defining the clip regions (C{LatLon}[])
                     or two points to specify a single, rectangular clip box.
       @kwarg closed: If C{True}, close each result clip (C{bool}).
       @kwarg inull: If C{True}, retain null edges in result clips (C{bool}).
       @kwarg raiser: If C{True}, throw L{ClipError} exceptions (C{bool}).
       @kwarg esp: Tolerance for eliminating null edges (C{degrees}, same units
                   as the B{C{points}} and B{C{corners}} coordinates).

       @return: Yield a L{ClipFHP4Tuple}C{(lat, lon, height, clipid)} for each
                clipped point.  The result may consist of several clips, each
                a (closed) polygon with a unique C{clipid}.

       @raise ClipError: Insufficient B{C{points}} or B{C{corners}} or an open clip.

       @see: U{Forster, Hormann and Popa<https://www.ScienceDirect.com/science/
             article/pii/S259014861930007X>}, class L{BooleanFHP} and function
             L{clipGH4}.
    '''
    P = _MODS.booleans._CompositeFHP(points, kind=_points_, raiser=raiser,
                                             name=clipFHP4.__name__, eps=eps)
    Q = _4corners(corners)
    return P._clip(Q, Union=False, Clas=ClipFHP4Tuple, closed=closed,
                                   inull=inull, raiser=P._raiser, eps=eps)


class ClipGH4Tuple(ClipFHP4Tuple):
    '''4-Tuple C{(lat, lon, height, clipid)} for each point of the
       L{clipGH4} result with the C{lat}-, C{lon}gitude, C{height}
       and C{clipid} of the polygon or clip.

       @note: The C{height} is a L{HeightX} instance if this is
              an intersection, otherwise a L{Height} or C{int(0)}.
    '''
    _Names_ = ClipFHP4Tuple._Names_
    _Units_ = ClipFHP4Tuple._Units_


def clipGH4(points, corners, closed=False, inull=False, raiser=True, xtend=False, eps=EPS):
    '''Clip one or more polygons against a clip region or box using the U{Greiner-Hormann
       <http://www.Inf.USI.CH/hormann/papers/Greiner.1998.ECO.pdf>} algorithm, U{Correia
       <https://GitHub.com/helderco/univ-polyclip>}'s implementation modified and extended.

       @arg points: The polygon points and clips (C{LatLon}[]).
       @arg corners: Three or more points defining the clip regions (C{LatLon}[])
                     or two points to specify a single, rectangular clip box.
       @kwarg closed: If C{True}, close each result clip (C{bool}).
       @kwarg inull: If C{True}, retain null edges in result clips (C{bool}).
       @kwarg raiser: If C{True}, throw L{ClipError} exceptions (C{bool}).
       @kwarg xtend: If C{True}, extend edges of I{degenerate cases}, an attempt
                     to handle the latter (C{bool}).
       @kwarg esp: Tolerance for eliminating null edges (C{degrees}, same units
                   as the B{C{points}} and B{C{corners}} coordinates).

       @return: Yield a L{ClipGH4Tuple}C{(lat, lon, height, clipid)} for each
                clipped point.  The result may consist of several clips, each
                a (closed) polygon with a unique C{clipid}.

       @raise ClipError: Insufficient B{C{points}} or B{C{corners}}, an open clip,
                         a I{degenerate case} or I{unhandled} intersection.

       @note: To handle I{degenerate cases} like C{point-edge} and C{point-point}
              intersections, use function L{clipFHP4}.

       @see: U{Greiner-Hormann<https://WikiPedia.org/wiki/Greiner–Hormann_clipping_algorithm>},
             U{Ionel Daniel Stroe<https://Davis.WPI.edu/~matt/courses/clipping/>}, I{Correia}'s
             U{univ-polyclip<https://GitHub.com/helderco/univ-polyclip>}, class L{BooleanGH}
             and function L{clipFHP4}.
    '''
    S = _MODS.booleans._CompositeGH(points, raiser=raiser, xtend=xtend, eps=eps,
                                            name=clipGH4.__name__, kind=_points_)
    C = _4corners(corners)
    return S._clip(C, False, False, Clas=ClipGH4Tuple, closed=closed, inull=inull,
                                    raiser=S._raiser, xtend=S._xtend, eps=eps)


def _LBtrim(p, q, t):
    # Liang-Barsky trim t[0] or t[1]
    if p < 0:
        r = q / p
        if r > t[1]:
            return False  # too far above
        elif r > t[0]:
            t[0] = r
    elif p > 0:
        r = q / p
        if r < t[0]:
            return False  # too far below
        elif r < t[1]:
            t[1] = r
    elif q < 0:  # vertical or horizontal
        return False  # ... outside
    return True


class ClipLB6Tuple(_NamedTuple):
    '''6-Tuple C{(start, end, i, fi, fj, j)} for each edge of the
       I{clipped} path with the C{start} and C{end} points (C{LatLon})
       of the portion of the edge inside or on the clip box, indices
       C{i} and C{j} (both C{int}) of the original path edge start
       and end points and I{fractional} indices C{fi} and C{fj}
       (both L{FIx}) of the C{start} and C{end} points along the
       edge of the original path.

       @see: Class L{FIx} and function L{pygeodesy.fractional}.
    '''
    _Names_ = (_start_, _end_, _i_,      _fi_,  _fj_, _j_)
    _Units_ = (_Pass,   _Pass,  Number_, _Pass, _Pass, Number_)


def clipLB6(points, lowerleft, upperight, closed=False, inull=False):
    '''Clip a path against a rectangular clip box using the U{Liang-Barsky
       <https://www.CSE.UNT.edu/~renka/4230/LineClipping.pdf>} algorithm.

       @arg points: The points (C{LatLon}[]).
       @arg lowerleft: Bottom-left corner of the clip box (C{LatLon}).
       @arg upperight: Top-right corner of the clip box (C{LatLon}).
       @kwarg closed: Optionally, close the path (C{bool}).
       @kwarg inull: Optionally, retain null edges if inside (C{bool}).

       @return: Yield a L{ClipLB6Tuple}C{(start, end, i, fi, fj, j)} for
                each edge of the I{clipped} path.

       @raise ClipError: The B{C{lowerleft}} and B{C{upperight}} corners
                         specify an invalid clip box.

       @raise PointsError: Insufficient number of B{C{points}}.

       @see: U{Liang-Barsky Line Clipping<https://www.CS.Helsinki.FI/group/goa/
             viewing/leikkaus/intro.html>}, U{Liang-Barsky line clipping algorithm
             <https://www.Skytopia.com/project/articles/compsci/clipping.html>} and
             U{Liang-Barsky algorithm<https://WikiPedia.org/wiki/Liang-Barsky_algorithm>}.
    '''
    xl, yb, \
    xr, yt = _box4(lowerleft, upperight, clipLB6.__name__)
    n, pts = _pts2(points, closed, inull)

    T6  =  ClipLB6Tuple
    fin =  n if closed else None  # wrapping fi [n] to [0]
    _LB = _LBtrim

    i, m = _imdex2(closed, n)
    for j in range(m, n):
        p1 = pts[i]
        y1 = p1.lat
        x1 = p1.lon

        p2 = pts[j]
        dy = float(p2.lat - y1)
        dx = float(p2.lon - x1)
        if fabs(dx) > EPS or fabs(dy) > EPS:
            # non-null edge pts[i]...pts[j]
            t = [_0_0, _1_0]
            if _LB(-dx, -xl + x1, t) and \
               _LB( dx,  xr - x1, t) and \
               _LB(-dy, -yb + y1, t) and \
               _LB( dy,  yt - y1, t):
                # clip edge pts[i]...pts[j]
                # at fractions t[0] to t[1]
                f, t = t
                if f > _0_0:  # EPS
                    p1 = p1.classof(y1 + f * dy,
                                    x1 + f * dx)
                    fi = i + f
                else:
                    fi = i

                if (t - f) > EPS:  # EPS0
                    if t < _1_0:  # EPS1
                        p2 = p2.classof(y1 + t * dy,
                                        x1 + t * dx)
                        fj = i + t
                    else:
                        fj = j
                    fi = FIx(fi, fin=fin)
                    fj = FIx(fj, fin=fin)
                    yield T6(p1, p2, i, fi, fj, j)

                elif inull:
                    fi = FIx(fi, fin=fin)
                    yield T6(p1, p1, i, fi, fi, j)
#           else:  # outside
#               pass
        elif inull:  # null edge
            yield T6(p1, p2, i, FIx(i, fin=fin),
                                FIx(j, fin=fin), j)
        i = j


class _SH(_Named):
    '''(INTERNAL) Sutherland-Hodgman polyon clipping.
    '''
    _cs  = ()    # clip corners
    _cw  =  0    # counter-/clockwise
    _ccw =  0    # clock-/counterwise
    _dx  = _0_0  # clip edge[e] delta lon
    _dy  = _0_0  # clip edge[e] delta lat
    _nc  =  0    # len(._cs)
    _x1  = _0_0  # clip edge[e] lon origin
    _xy  = _0_0  # see .clipedges
    _y1  = _0_0  # clip edge[e] lat origin

    def __init__(self, corners, name=__name__):
        n, cs = 0, corners
        try:  # check the clip box/region
            cs = _4corners(cs)
            n, cs = len2(cs)
            n, cs = points2(cs, closed=True)
            self._cs = cs = cs[:n]
            self._nc = n
            self._cw = isconvex_(cs, adjust=False, wrap=False)
            if not self._cw:
                raise ValueError(_not_(_convex_))
            if areaOf(cs, adjust=True, radius=1, wrap=True) < EPS:
                raise ValueError(NN(_near_, 'zero area'))
            self._ccw = -self._cw
        except (PointsError, TypeError, ValueError) as x:
            raise ClipError(name, n, cs, cause=x)
        self.name = name

    def clip2(self, points, closed, inull):  # MCCABE 13, clip points
        np, pts = _pts2(points, closed, inull)
        pcs = _SHlist(inull)  # clipped points
        _ap = pcs.append
        _d2 = self.dot2
        _in = self.intersect

        ne = 0  # number of non-null clip edges
        for e in self.clipedges():
            ne += 1  # non-null clip edge

            # clip points, closed always
            d1, p1 = _d2(pts[np - 1])
            for i in range(np):
                d2, p2 = _d2(pts[i])
                if d1 < 0:  # p1 inside, p2 ...
                    # _ap(p1)
                    _ap(p2 if d2 < 0 else  # ... in-
                        _in(p1, p2, e))  # ... outside
                elif d2 < 0:  # p1 out-, p2 inside
                    _ap(_in(p1, p2, e))
                    _ap(p2)
#                   elif d1 > 0:  # both outside
#                       pass
                d1, p1 = d2, p2

            # replace points, in-place
            pts[:] = pcs
            pcs[:] = []
            np = len(pts)
            if not np:  # all outside
                break
        else:
            if ne < 3:
                raise ClipError(self.name, ne, self._cs, txt=_too_(_few_))

        if np > 1:
            p = pts[0]
            if closed:  # close clipped pts
                if _neq(pts[np - 1], p):
                    pts.append(p)
                    np += 1
            elif not inull:  # open clipped pts
                while np > 0 and _eq(pts[np - 1], p):
                    pts.pop()
                    np -= 1
        # assert len(pts) == np
        return np, pts

    def clipedges(self):  # yield clip edge index
        # and set self._x1, ._y1, ._dx, ._dy and
        # ._xy for each non-null clip edge
        nc = self._nc
        cs = self._cs
        c  = cs[nc - 1]
        for e in range(nc):
            y, x, c = c.lat, c.lon, cs[e]
            dy = float(c.lat - y)
            dx = float(c.lon - x)
            if fabs(dx) > EPS or fabs(dy) > EPS:
                self._y1, self._dy = y, dy
                self._x1, self._dx = x, dx
                self._xy = y * dx - x * dy
                yield e + 1

    def clipped2(self, p):  # return (clipped point [i], edge)
        if isinstance(p, _SHlli):  # intersection point
            return p.classof(p.lat, p.lon), p.edge
        else:  # original point
            return p, 0

    def dot2(self, p):  # dot product of point p to clip
        # corner c1 and clip edge c1 to c2, indicating where
        # point p is located: to the right, to the left or
        # on top of the (extended) clip edge from c1 to c2
        d = float(p.lat - self._y1) * self._dx - \
            float(p.lon - self._x1) * self._dy
        # clockwise corners, +1 means point p is to the right
        # of, -1 means on the left of, 0 means on edge c1 to c2
        d = self._ccw if d < 0 else (self._cw if d > 0 else 0)
        return d, p

    def intersect(self, p1, p2, edge):  # compute intersection
        # of polygon edge p1 to p2 and the current clip edge,
        # where p1 and p2 are known to NOT be located on the
        # same side of or on the current, non-null clip edge
        # <https://StackOverflow.com/questions/563198/
        #        how-do-you-detect-where-two-line-segments-intersect>
        y, dy = p1.lat, self._dy
        x, dx = p1.lon, self._dx
        fy = float(p2.lat - y)
        fx = float(p2.lon - x)
        d  = fy * dx - fx * dy
        if fabs(d) < EPS:  # PYCHOK no cover
            raise _AssertionError(self._DOT_(self.intersect.__name__))
        d  = fsumf_(self._xy, -y * dx, x * dy) / d
        y += d * fy
        x += d * fx
        return _SHlli(y, x, p1.classof, edge)


class _SHlist(list):
    '''(INTERNAL) List of _SH clipped points.
    '''
    _inull = False

    def __init__(self, inull):
        self._inull = inull
        list.__init__(self)

    def append(self, p):
        if (not self) or self._inull or _neq(p, self[-1]):
            list.append(self, p)


class _SHlli(LatLon_):
    '''(INTERNAL) LatLon_ for _SH intersections.
    '''
    # __slots__ are no longer space savers, see
    # the comments at the class .points.LatLon_
    # __slots__ = _lat_, _lon_, 'classof', 'edge', _name_

    def __init__(self, lat, lon, classof, edge):
        self.lat     = lat
        self.lon     = lon
        self.classof = classof
        self.edge    = edge  # clip edge
        self.name    = NN


class ClipSH3Tuple(_NamedTuple):
    '''3-Tuple C{(start, end, original)} for each edge of a I{clipped}
       polygon, the C{start} and C{end} points (C{LatLon}) of the
       portion of the edge inside or on the clip region and C{original}
       indicates whether the edge is part of the original polygon or
       part of the clip region (C{bool}).
    '''
    _Names_ = (_start_, _end_, _original_)
    _Units_ = (_Pass,   _Pass,  Bool)


def clipSH(points, corners, closed=False, inull=False):
    '''Clip a polygon against a clip region or box using the U{Sutherland-Hodgman
       <https://WikiPedia.org/wiki/Sutherland-Hodgman_algorithm>} algorithm.

       @arg points: The polygon points (C{LatLon}[]).
       @arg corners: Three or more points defining a convex clip
                     region (C{LatLon}[]) or two points to specify
                     a rectangular clip box.
       @kwarg closed: Close the clipped points (C{bool}).
       @kwarg inull: Optionally, include null edges (C{bool}).

       @return: Yield the clipped points (C{LatLon}[]).

       @raise ClipError: The B{C{corners}} specify a polar, zero-area,
                         non-convex or otherwise invalid clip box or
                         region.

       @raise PointsError: Insufficient number of B{C{points}}.
    '''
    sh = _SH(corners, name=clipSH.__name__)
    n, pts = sh.clip2(points, closed, inull)
    for i in range(n):
        p, _ = sh.clipped2(pts[i])
        yield p


def clipSH3(points, corners, closed=False, inull=False):
    '''Clip a polygon against a clip region or box using the U{Sutherland-Hodgman
       <https://WikiPedia.org/wiki/Sutherland-Hodgman_algorithm>} algorithm.

       @arg points: The polygon points (C{LatLon}[]).
       @arg corners: Three or more points defining a convex clip
                     region (C{LatLon}[]) or two points to specify
                     a rectangular clip box.
       @kwarg closed: Close the clipped points (C{bool}).
       @kwarg inull: Optionally, include null edges (C{bool}).

       @return: Yield a L{ClipSH3Tuple}C{(start, end, original)} for
                each edge of the I{clipped} polygon.

       @raise ClipError: The B{C{corners}} specify a polar, zero-area,
                         non-convex or otherwise invalid clip box or
                         region.

       @raise PointsError: Insufficient number of B{C{points}} or B{C{corners}}.
    '''
    sh = _SH(corners, name=clipSH3.__name__)
    n, pts = sh.clip2(points, closed, inull)
    if n > 1:
        T3 = ClipSH3Tuple
        p1, e1 = sh.clipped2(pts[0])
        for i in range(1, n):
            p2, e2 = sh.clipped2(pts[i])
            yield T3(p1, p2, not bool(e1 and e2 and e1 == e2))
            p1, e1 = p2, e2

# **) MIT License
#
# Copyright (C) 2018-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
