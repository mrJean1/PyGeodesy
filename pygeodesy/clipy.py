
# -*- coding: utf-8 -*-

u'''Functions to I{clip} a path or polygon of C{LatLon} points
against a rectangular box or a (convex) clip region.

@newfield example: Example, Examples
'''

from pygeodesy.basics import len2
from pygeodesy.errors import _AssertionError, PointsError, _ValueError
from pygeodesy.fmath import fsum_
from pygeodesy.formy import points2
from pygeodesy.interns import EPS, NN, _convex_, _DOT_, _end_, _few_, \
                             _i_, _j_, _not_, _SPACE_, _start_, _too_, \
                             _0_0, _1_0
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import _Named, _NamedTuple, _Pass
from pygeodesy.points import areaOf, _imdex2, boundsOf, isconvex_, \
                             LatLon_
from pygeodesy.units import Bool, FIx, Number_

__all__ = _ALL_LAZY.clipy
__version__ = '21.01.11'

_fi_ = 'fi'
_fj_ = 'fj'


class ClipError(_ValueError):
    '''Clip box or clip region issue.
    '''
    def __init__(self, *name_n_corners, **txt):
        '''New L{ClipError}.

           @arg name_n_corners: Either just a name (C{str}) or
                                name, number, corners (C{str},
                                C{int}, C{tuple}).
           @kwarg txt: Optional explanation of the error (C{str}).
        '''
        if len(name_n_corners) == 3:
            t, n, v = name_n_corners
            n = _SPACE_(t, 'clip', 'box' if n == 2 else 'region')
            name_n_corners = n, v
        _ValueError.__init__(self, *name_n_corners, **txt)


def _box4(lowerleft, upperright, name):
    '''(INTERNAL) Get the clip box edges.
    '''
    try:
        ymin, ymax = lowerleft.lat, upperright.lat
        xmin, xmax = lowerleft.lon, upperright.lon
        if xmin > xmax or ymin > ymax:
            raise ValueError
    except (AttributeError, TypeError, ValueError) as x:
        raise ClipError(name, 2, (lowerleft, upperright), txt=str(x))
    return xmin, ymin, xmax, ymax


def _eq(p1, p2):
    '''(INTERNAL) Check for near-equal points.
    '''
    return not _neq(p1, p2)


def _neq(p1, p2):
    '''(INTERNAL) Check for not near-equal points.
    '''
    return abs(p1.lat - p2.lat) > EPS or \
           abs(p1.lon - p2.lon) > EPS


def _pts2(points, closed, inull):
    '''(INTERNAL) Get the points to clip.
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
        return abs(self._dx) > EPS or abs(self._dy) > EPS

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


class ClipLB6Tuple(_NamedTuple):
    '''6-Tuple C{(start, end, i, fi, fj, j)} for each edge of
       the I{clipped} path with the C{start} and C{end} points
       (C{LatLon}) of the portion of the edge inside or on the
       clip box, indices C{i} and C{j} (C{int}) of the original
       path edge start and end points and I{fractional} indices
       C{fi} and C{fj} (L{FIx}) of the C{start} and C{end} points
       along the edge of the original path.

       @see: Class L{FIx} and function L{fractional}.
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
             and U{Liang–Barsky algorithm<https://WikiPedia.org/wiki/Liang–Barsky_algorithm>}.
    '''
    xmin, ymin, \
    xmax, ymax = _box4(lowerleft, upperright, clipLB6.__name__)

    n, pts = _pts2(points, closed, inull)

    fin =  n if closed else None  # wrapping fi [n] to [0]
    _LB = _LBtrim

    i, m = _imdex2(closed, n)
    pi = pts[i]
    for j in range(m, n):
        p1 = pi
        y1 = p1.lat
        x1 = p1.lon

        p2 = pi = pts[j]
        dy = float(p2.lat - y1)
        dx = float(p2.lon - x1)
        if abs(dx) > EPS or abs(dy) > EPS:
            # non-null edge pts[i]...pts[j]
            t = [_0_0, _1_0]
            if _LB(-dx, -xmin + x1, t) and \
               _LB( dx,  xmax - x1, t) and \
               _LB(-dy, -ymin + y1, t) and \
               _LB( dy,  ymax - y1, t):
                # clip edge pts[i]...pts[j]
                # at fractions t[0] to t[1]
                if t[0] > _0_0:  # EPS
                    p1 = p1.classof(y1 + t[0] * dy,
                                    x1 + t[0] * dx)
                    fi = i + t[0]
                else:
                    fi = i
                if t[0] < t[1]:
                    if t[1] < _1_0:  # EPS1
                        p2 = p2.classof(y1 + t[1] * dy,
                                        x1 + t[1] * dx)
                        fj = i + t[1]
                    else:
                        fj = j
                    fi = FIx(fi, fin=fin)
                    fj = FIx(fj, fin=fin)
                    yield ClipLB6Tuple(p1, p2, i, fi, fj, j)

                elif inull and (t[0] + EPS) > t[1]:
                    fi = FIx(fi, fin=fin)
                    yield ClipLB6Tuple(p1, p1, i, fi, fi, j)
#           else:  # outside
#               pass
        elif inull:  # null edge
            yield ClipLB6Tuple(p1, p2, i, FIx(i, fin=fin),
                                          FIx(j, fin=fin), j)
        i = j


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


class _List(list):
    '''(INTERNAL) List of clipped points.
    '''
    _inull = False

    def __init__(self, inull):
        self._inull = inull
        list.__init__(self)

    def append(self, p):
        if (not self) or self._inull or _neq(p, self[-1]):
            list.append(self, p)


class _LLi_(LatLon_):
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
            n, cs = len2(cs)
            if n == 2:  # make a box
                b, l, t, r = boundsOf(cs, wrap=False)
                cs = (LatLon_(b, l), LatLon_(t, l),
                      LatLon_(t, r), LatLon_(b, r))
            n, cs = points2(cs, closed=True)
            self._cs = cs = cs[:n]
            self._nc = n
            self._cw = isconvex_(cs, adjust=False, wrap=False)
            if not self._cw:
                raise ValueError(_not_(_convex_))
            if areaOf(cs, adjust=True, radius=1, wrap=True) < EPS:
                raise ValueError('near-zero area')
        except (PointsError, TypeError, ValueError) as x:
            raise ClipError(name, n, cs, txt=str(x))
        self.name = name

    def clip2(self, points, closed, inull):  # MCCABE 17, clip points
        np, pts = _pts2(points, closed, inull)
        pcs = _List(inull)  # clipped points

        ne = 0  # number of non-null clip edges
        no = ni = True  # all out- and inside?
        for e in self.clipedges():
            ne += 1  # non-null clip edge

            # clip points, closed always
            d2, p2 = self.dot2(pts[np - 1])
            for i in range(np):
                d1, p1 = d2, p2
                d2, p2 = self.dot2(pts[i])
                if d1 < 0:  # p1 inside, ...
                    # pcs.append(p1)
                    if d2 < 0:  # ... p2 inside
                        p = p2
                    else:  # ... p2 outside
                        p = self.intersect(p1, p2, e)
                        if d2 > 0:
                            no = False
                    pcs.append(p)
                elif d2 < 0:  # p1 out-, p2 inside
                    p = self.intersect(p1, p2, e)
                    pcs.append(p)
                    pcs.append(p2)
                    if d1 > 0:
                        no = ni = False
                elif d1 > 0:  # both out
                    ni = False

            if pcs:  # replace points
                pts[:] = pcs
                pcs[:] = []
                np = len(pts)

        if ne < 3:
            raise ClipError(self.name, ne, self._cs, txt=_too_(_few_))

        # ni is True iff all points are on or on the
        # right side (i.e. inside) of all clip edges,
        # no is True iff all points are on or at one
        # side (left or right) of each clip edge: if
        # none are inside (ni is False) and if all
        # are on the same side (no is True), then all
        # must be outside
        if no and not ni:
            np, pts = 0, []
        elif np > 1:
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
            if abs(self._dx) > EPS or abs(self._dy) > EPS:
                self._xy = self._y1 * self._dx - self._x1 * self._dy
                yield e + 1

    def clipped2(self, p):  # return (clipped point [i], edge)
        if isinstance(p, _LLi_):  # intersection point
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
        fy = float(p2.lat - p1.lat)
        fx = float(p2.lon - p1.lon)
        fp = fy * self._dx - fx * self._dy
        if abs(fp) < EPS:  # PYCHOK no cover
            raise _AssertionError(self._DOT_(self.intersect.__name__))
        r = fsum_(self._xy, -p1.lat * self._dx, p1.lon * self._dy) / fp
        y = p1.lat + r * fy
        x = p1.lon + r * fx
        return _LLi_(y, x, p1.classof, edge)


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
    '''Clip a polygon against a clip region or box using the
       U{Sutherland-Hodgman
       <https://WikiPedia.org/wiki/Sutherland_Hodgman_algorithm>} algorithm.

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
    '''Clip a polygon against a clip region or box using the
       U{Sutherland-Hodgman
       <https://WikiPedia.org/wiki/Sutherland_Hodgman_algorithm>} algorithm.

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

       @raise PointsError: Insufficient number of B{C{points}}.
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
# Copyright (C) 2018-2021 -- mrJean1 at Gmail -- All Rights Reserved.
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
