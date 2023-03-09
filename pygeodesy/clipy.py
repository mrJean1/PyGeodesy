
# -*- coding: utf-8 -*-

u'''Clip a path or polygon of C{LatLon} points against a rectangular box or a
(convex or arbitrary) clip region.

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
from pygeodesy.fsums import fsum_, Property_RO
from pygeodesy.interns import NN, _convex_, _DOT_, _end_, _few_, _fi_, \
                             _height_, _i_, _j_, _lat_, _lon_, _near_, \
                             _not_, _points_, _start_, _too_
from pygeodesy.iters import _imdex2, points2
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import _Named, _NamedTuple, _Pass
from pygeodesy.points import areaOf, boundsOf, isconvex_, LatLon_
# from pygeodesy.props import Property_RO  # from .fsums
from pygeodesy.units import Bool, FIx, HeightX, Lat, Lon, Number_

# from math import fabs  # from .fmath

__all__ = _ALL_LAZY.clipy
__version__ = '23.03.09'

_clipid_ = 'clipid'
_fj_     = 'fj'


def _box4(lowerleft, upperright, name):
    '''(INTERNAL) Get the clip box edges.

       @see: Class C{_Box} in .ellipsoidalBaseDI.py.
    '''
    try:
        ymin, ymax = lowerleft.lat, upperright.lat
        xmin, xmax = lowerleft.lon, upperright.lon
        if xmin > xmax or ymin > ymax:
            raise ValueError
    except (AttributeError, TypeError, ValueError) as x:
        raise ClipError(name, 2, (lowerleft, upperright), cause=x)
    return xmin, ymin, xmax, ymax


def _4corners(corners):
    '''(INTERNAL) Clip region or box.
    '''
    n, cs = len2(corners)
    if n == 2:  # make a box
        b, l, t, r = boundsOf(cs, wrap=False)
        cs = (LatLon_(b, l), LatLon_(t, l),
              LatLon_(t, r), LatLon_(b, r))
    return cs


def _eq(p1, p2):
    '''(INTERNAL) Check for near-equal points.
    '''
    return not _neq(p1, p2)


def _neq(p1, p2):
    '''(INTERNAL) Check for not near-equal points.
    '''
    return fabs(p1.lat - p2.lat) > EPS or \
           fabs(p1.lon - p2.lon) > EPS


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
    _IN   =  0  # inside clip box
    _XMAX =  1  # right of upperright.lon
    _XMIN =  2  # left of lowerleft.lon
    _YMAX =  4  # above upperright.lat
    _YMIN =  8  # below lowerleft.lat

    _dx   = _0_0  # pts edge delta lon
    _dy   = _0_0  # pts edge delta lat
    _x1   = _0_0  # pts corner
    _y1   = _0_0  # pts corner

    _xmax = _0_0  # clip box upperright.lon
    _xmin = _0_0  # clip box lowerleft.lon
    _ymax = _0_0  # clip box upperright.lat
    _ymin = _0_0  # clip box lowerleft.lat

    def __init__(self, lowerleft, upperright, name=__name__):
        self._xmin, self._ymin, \
        self._xmax, self._ymax = _box4(lowerleft, upperright, name)
        self.name = name

#   def clip4(self, p, c):  # clip point p for code c
#       if c & _CS._YMIN:
#           return self.lon4(p, self._ymin)
#       elif c & _CS._YMAX:
#           return self.lon4(p, self._ymax)
#       elif c & _CS._XMIN:
#           return self.lat4(p, self._xmin)
#       elif c & _CS._XMAX:
#           return self.lat4(p, self._xmax)
#       # should never get here
#       raise _AssertionError(self._DOT_(self.clip4.__name__))

    def code4(self, p):  # compute code for point p
        if p.lat < self._ymin:
            c, m, b = _CS._YMIN, self.lon4, self._ymin
        elif p.lat > self._ymax:
            c, m, b = _CS._YMAX, self.lon4, self._ymax
        else:
            c, m, b = _CS._IN, self.nop4, None
        if p.lon < self._xmin:
            c |= _CS._XMIN
            m, b = self.lat4, self._xmin
        elif p.lon > self._xmax:
            c |= _CS._XMAX
            m, b = self.lat4, self._xmax
        return c, m, b, p

    def edge(self, p1, p2):  # set edge p1 to p2
        self._y1, self._dy = p1.lat, float(p2.lat - p1.lat)
        self._x1, self._dx = p1.lon, float(p2.lon - p1.lon)
        return fabs(self._dx) > EPS or fabs(self._dy) > EPS

    def lat4(self, x, p):  # new lat and code at lon x
        y = self._y1 + self._dy * float(x - self._x1) / self._dx
        if y < self._ymin:  # still outside
            return _CS._YMIN, self.lon4, self._ymin, p
        elif y > self._ymax:  # still outside
            return _CS._YMAX, self.lon4, self._ymax, p
        else:  # inside
            return _CS._IN, self.nop4, None, p.classof(y, x)

    def lon4(self, y, p):  # new lon and code at lat y
        x = self._x1 + self._dx * float(y - self._y1) / self._dy
        if x < self._xmin:  # still outside
            return _CS._XMIN, self.lat4, self._xmin, p
        elif x > self._xmax:  # still outside
            return _CS._XMAX, self.lat4, self._xmax, p
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


def clipCS4(points, lowerleft, upperright, closed=False, inull=False):
    '''Clip a path against a rectangular clip box using the U{Cohen-Sutherland
       <https://WikiPedia.org/wiki/Cohen-Sutherland_algorithm>} algorithm.

       @arg points: The points (C{LatLon}[]).
       @arg lowerleft: Bottom-left corner of the clip box (C{LatLon}).
       @arg upperright: Top-right corner of the clip box (C{LatLon}).
       @kwarg closed: Optionally, close the path (C{bool}).
       @kwarg inull: Optionally, retain null edges if inside (C{bool}).

       @return: Yield a L{ClipCS4Tuple}C{(start, end, i, j)} for each
                edge of the I{clipped} path.

       @raise ClipError: The B{C{lowerleft}} and B{C{upperright}} corners
                         specify an invalid clip box.

       @raise PointsError: Insufficient number of B{C{points}}.
    '''
    cs = _CS(lowerleft, upperright, name=clipCS4.__name__)
    n, pts = _pts2(points, closed, inull)

    i, m = _imdex2(closed, n)
    cmbp = cs.code4(pts[i])
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
                        yield ClipCS4Tuple(p1, p2, i, j)
                    break
                if c1 & c2:  # edge outside
                    break
            else:  # PYCHOK no cover
                raise _AssertionError(_DOT_(cs.name, 'for_else'))

        elif inull and not c1:  # null edge
            yield ClipCS4Tuple(p1, p1, i, j)
        elif inull and not c2:
            yield ClipCS4Tuple(p2, p2, i, j)

        i = j


class ClipFHP4Tuple(_NamedTuple):
    '''4-Tuple C{(lat, lon, height, clipid)} for each point of the
       L{clipFHP4} result with the C{lat}-, C{lon}gitude, C{height}
       and C{clipid} of the polygon or polygonal clip.

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
    '''Clip a polygon against a clip region or box using U{Forster-Hormann-Popa
       <https://www.ScienceDirect.com/science/article/pii/S259014861930007X>}'s
       C++ implementation transcoded to pure Python.

       @arg points: The polygon points and clips (C{LatLon}[]).
       @arg corners: Three or more points defining a clip region (C{LatLon}[])
                     or two points to specify a rectangular clip box.
       @kwarg closed: Optionally, close each result clip (C{bool}).
       @kwarg inull: Optionally, retain null edges in result clips (C{bool}).
       @kwarg raiser: If C{True} throw L{ClipError} exceptions.
       @kwarg esp: Tolerance for eliminating null edges (C{scalar}).

       @return: Yield a L{ClipFHP4Tuple}C{(lat, lon, height, clipid)} for each
                clipped point.  The result may consist of several clips, each
                a (closed) polygon with a unique C{clipid}.

       @raise ClipError: No B{C{points}}, no B{C{corners}} or an open clip.

       @see: U{Forster, Hormann and Popa<https://www.ScienceDirect.com/science/
             article/pii/S259014861930007X>}, class L{BooleanFHP} and function
             L{clipGH4}.
    '''
    P = _MODS.booleans._CompositeFHP(points, kind=_points_, raiser=raiser,
                                             name=clipFHP4.__name__, eps=eps)
    Q = _4corners(corners)
    return P._clip(Q, Union=False, Clas=ClipFHP4Tuple, closed=closed,
                                   inull=inull, raiser=P._raiser)


class ClipGH4Tuple(ClipFHP4Tuple):
    '''4-Tuple C{(lat, lon, height, clipid)} for each point of the
       L{clipGH4} result with the C{lat}-, C{lon}gitude, C{height}
       and C{clipid} of the polygon or polygonal clip.

       @note: The C{height} is a L{HeightX} instance if this is
              an intersection, otherwise a L{Height} or C{int(0)}.
    '''
    _Names_ = ClipFHP4Tuple._Names_
    _Units_ = ClipFHP4Tuple._Units_


def clipGH4(points, corners, closed=False, inull=False, raiser=False, xtend=False):
    '''Clip a polygon against a clip region or box using the U{Greiner-Hormann
       <http://www.inf.USI.CH/hormann/papers/Greiner.1998.ECO.pdf>} algorithm,
       U{Correia<https://GitHub.com/helderco/univ-polyclip>}'s implementation
       modified and extended.

       @arg points: The polygon points and clips (C{LatLon}[]).
       @arg corners: Three or more points defining a clip region (C{LatLon}[])
                     or two points to specify a rectangular clip box.
       @kwarg closed: Optionally, close each result clip (C{bool}).
       @kwarg inull: Optionally, retain null edges in result clips (C{bool}).
       @kwarg raiser: If C{True} throw L{ClipError} exceptions.
       @kwarg xtend: Cover degenerate cases (C{bool}), but not (yet) based on
                     U{Forster-Hormann-Popa<https://www.ScienceDirect.com/
                     science/article/pii/S259014861930007X>}.

       @return: Yield a L{ClipGH4Tuple}C{(lat, lon, height, clipid)} for each
                clipped point.  The result may consist of several clips, each
                a (closed) polygon with a unique C{clipid}.

       @raise ClipError: No B{C{points}}, no B{C{corners}}, an open clip, a
                         I{degenerate case} or I{unhandled} intersection.

       @note: To handle I{degenerate cases} like C{point-edge} and C{point-point}
              intersections, use function L{clipFHP4}.

       @see: U{Greiner-Hormann<https://WikiPedia.org/wiki/Greiner–Hormann_clipping_algorithm>},
             U{Ionel Daniel Stroe<https://Davis.WPI.edu/~matt/courses/clipping/>}, U{Forster,
             Hormann and Popa<https://www.ScienceDirect.com/science/article/pii/S259014861930007X>},
             I{Correia}'s U{univ-polyclip<https://GitHub.com/helderco/univ-polyclip>}, class
             L{BooleanGH} and function L{clipFHP4}.
    '''
    S = _MODS.booleans._CompositeGH(points, raiser=raiser, xtend=xtend,
                                            name=clipGH4.__name__, kind=_points_)
    C = _4corners(corners)
    return S._clip(C, False, False, Clas=ClipGH4Tuple, closed=closed, inull=inull,
                                    raiser=S._raiser, xtend=S._xtend)


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


def clipLB6(points, lowerleft, upperright, closed=False, inull=False):
    '''Clip a path against a rectangular clip box using the U{Liang-Barsky
       <https://www.CSE.UNT.edu/~renka/4230/LineClipping.pdf>} algorithm.

       @arg points: The points (C{LatLon}[]).
       @arg lowerleft: Bottom-left corner of the clip box (C{LatLon}).
       @arg upperright: Top-right corner of the clip box (C{LatLon}).
       @kwarg closed: Optionally, close the path (C{bool}).
       @kwarg inull: Optionally, retain null edges if inside (C{bool}).

       @return: Yield a L{ClipLB6Tuple}C{(start, end, i, fi, fj, j)} for
                each edge of the I{clipped} path.

       @raise ClipError: The B{C{lowerleft}} and B{C{upperright}} corners
                         specify an invalid clip box.

       @raise PointsError: Insufficient number of B{C{points}}.

       @see: U{Liang-Barsky Line Clipping<https://www.CS.Helsinki.FI/group/goa/
             viewing/leikkaus/intro.html>}, U{Liang-Barsky line clipping algorithm
             <https://www.Skytopia.com/project/articles/compsci/clipping.html>}
             and U{Liang-Barsky algorithm<https://WikiPedia.org/wiki/Liang-Barsky_algorithm>}.
    '''
    xmin, ymin, \
    xmax, ymax = _box4(lowerleft, upperright, clipLB6.__name__)

    n, pts = _pts2(points, closed, inull)

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
            if _LB(-dx, -xmin + x1, t) and \
               _LB( dx,  xmax - x1, t) and \
               _LB(-dy, -ymin + y1, t) and \
               _LB( dy,  ymax - y1, t):
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
                    yield ClipLB6Tuple(p1, p2, i, fi, fj, j)

                elif inull:
                    fi = FIx(fi, fin=fin)
                    yield ClipLB6Tuple(p1, p1, i, fi, fi, j)
#           else:  # outside
#               pass
        elif inull:  # null edge
            yield ClipLB6Tuple(p1, p2, i, FIx(i, fin=fin),
                                          FIx(j, fin=fin), j)
        i = j


class _SH(_Named):
    '''(INTERNAL) Sutherland-Hodgman polyon clipping.
    '''
    _cs = ()    # clip corners
    _cw =  0    # counter-/clockwise
    _dx = _0_0  # clip edge[e] delta lon
    _dy = _0_0  # clip edge[e] delta lat
    _nc =  0    # len(._cs)
    _x1 = _0_0  # clip edge[e] lon origin
    _xy = _0_0  # see .clipedges
    _y1 = _0_0  # clip edge[e] lat origin

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
        except (PointsError, TypeError, ValueError) as x:
            raise ClipError(name, n, cs, cause=x)
        self.name = name

    def clip2(self, points, closed, inull):  # MCCABE 13, clip points
        np, pts = _pts2(points, closed, inull)
        pcs = _SHlist(inull)  # clipped points
        pca = pcs.append

        ne = 0  # number of non-null clip edges
        for e in self.clipedges():
            ne += 1  # non-null clip edge

            # clip points, closed always
            d2, p2 = self.dot2(pts[np - 1])
            for i in range(np):
                d1, p1 = d2, p2
                d2, p2 = self.dot2(pts[i])
                if d1 < 0:  # p1 inside, ...
                    # pca(p1)
                    if d2 < 0:  # ... p2 inside
                        p = p2
                    else:  # ... p2 outside
                        p = self.intersect(p1, p2, e)
                    pca(p)
                elif d2 < 0:  # p1 out-, p2 inside
                    p = self.intersect(p1, p2, e)
                    pca(p)
                    pca(p2)
#               elif d1 > 0:  # both outside
#                   pass
            # replace points, in-place
            pts[:] = pcs
            pcs[:] = []
            np = len(pts)
            if not np:  # all out
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
        c2 = self._cs[self._nc - 1]
        for e in range(self._nc):
            c1, c2 = c2, self._cs[e]
            self._y1, self._dy = c1.lat, float(c2.lat - c1.lat)
            self._x1, self._dx = c1.lon, float(c2.lon - c1.lon)
            if fabs(self._dx) > EPS or fabs(self._dy) > EPS:
                self._xy = self._y1 * self._dx - self._x1 * self._dy
                yield e + 1

    def clipped2(self, p):  # return (clipped point [i], edge)
        if isinstance(p, _SHlli):  # intersection point
            return p.classof(p.lat, p.lon), p.edge
        else:  # original point
            return p, 0

    def dot2(self, p):  # dot product of point p to the current
        # clip corner c1 and clip edge c1 to c2, indicating whether
        # point p is located to the right, to the left or on top
        # of the (extended) clip edge from c1 to c2
        d = self._dx * float(p.lat - self._y1) - \
            self._dy * float(p.lon - self._x1)
        # clockwise corners, +1 means point p is to the right
        # of, -1 means on the left of, 0 means on edge c1 to c2
        d = (-self._cw) if d < 0 else (self._cw if d > 0 else 0)
        return d, p

    def intersect(self, p1, p2, edge):  # compute intersection
        # of polygon edge p1 to p2 and the current clip edge,
        # where p1 and p2 are known to NOT be located on the
        # same side or on top of the current clip edge
        # <https://StackOverflow.com/questions/563198/
        #        how-do-you-detect-where-two-line-segments-intersect>
        y, x = p1.lat, p1.lon
        fy = float(p2.lat - y)
        fx = float(p2.lon - x)
        fp = fy * self._dx - fx * self._dy
        if fabs(fp) < EPS:  # PYCHOK no cover
            raise _AssertionError(self._DOT_(self.intersect.__name__))
        r = fsum_(self._xy, -y * self._dx, x * self._dy) / fp
        y += r * fy
        x += r * fx
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
    _Names_ = (_start_, _end_, 'original')
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
        yield sh.clipped2(pts[i])[0]


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
    if n > 0:
        p2, e2 = sh.clipped2(pts[0])
        for i in range(1, n):
            p1, e1 = p2, e2
            p2, e2 = sh.clipped2(pts[i])
            yield ClipSH3Tuple(p1, p2, not bool(e1 and e2 and e1 == e2))

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

# ***) GNU GPL 3
#
# Copyright (C) 2011-2012 Helder Correia <Helder.MC@Gmail.com>
#
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.GNU.org/licenses/>.
#
# You should have received the README file along with this program.
# If not, see <https://GitHub.com/helderco/univ-polyclip>.
