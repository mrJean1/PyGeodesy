
# -*- coding: utf-8 -*-

u'''Functions to I{clip} a path or polygon of C{LatLon} points
against a rectangular box or clip region.

@newfield example: Example, Examples
'''

from pygeodesy.basics import EPS, len2
from pygeodesy.errors import _AssertionError, PointsError, _ValueError
from pygeodesy.fmath import fsum_
from pygeodesy.formy import points2
from pygeodesy.interns import _dot_, _end_, _lat_, _lon_, _name_, \
                               NN, _not_convex_, _start_, _too_few_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import _Named, _NamedTuple
from pygeodesy.points import areaOf, _imdex2, boundsOf, isconvex_, \
                             LatLon_ as LL_

__all__ = _ALL_LAZY.clipy
__version__ = '20.07.08'


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
            n = 'box' if n == 2 else 'region'
            n = '%s clip %s' % (t, n)
            name_n_corners = n, v
        _ValueError.__init__(self, *name_n_corners, **txt)


def _eq(p1, p2):
    '''(INTERNAL) Check for near-equal points.
    '''
    return not _neq(p1, p2)


def _neq(p1, p2):
    '''(INTERNAL) Check for not near-equal points.
    '''
    return abs(p1.lat - p2.lat) > EPS or \
           abs(p1.lon - p2.lon) > EPS


def _points2(points, closed, inull):
    '''(INTERNAL) Get the points to clip.
    '''
    if closed and inull:
        n, pts = len2(points)
        # only remove the final, closing point
        if n > 1 and _eq(pts[n-1], pts[0]):
            n -= 1
            pts = pts[:n]
        if n < 3:
            raise PointsError(points=n, txt=_too_few_)
    else:
        n, pts = points2(points, closed=closed)
    return n, list(pts)


class _CS(_Named):
    '''(INTERNAL) Cohen-Sutherland line clipping.
    '''
    # single-bit clip codes
    _IN   = 0  # inside clip box
    _YMIN = 1  # below lowerleft.lat
    _YMAX = 2  # above upperright.lat
    _XMIN = 4  # left of lowerleft.lon
    _XMAX = 8  # right of upperright.lon

    _dx   = 0  # pts edge delta lon
    _dy   = 0  # pts edge delta lat
    _x1   = 0  # pts corner
    _y1   = 0  # pts corner

    _xmax = 0  # clip box upperright.lon
    _xmin = 0  # clip box lowerleft.lon
    _ymax = 0  # clip box upperright.lat
    _ymin = 0  # clip box lowerleft.lat

    def __init__(self, lowerleft, upperright, name=__name__):
        try:
            self._ymin, self._ymax = lowerleft.lat, upperright.lat
            self._xmin, self._xmax = lowerleft.lon, upperright.lon
            if self._xmin > self._xmax or self._ymin > self._ymax:
                raise ValueError
        except (AttributeError, TypeError, ValueError) as x:
            raise ClipError(name, 2, (lowerleft, upperright), txt=str(x))
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
#       raise _AssertionError(_dot_(self._name, self.clip4.__name__))

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
            raise _AssertionError(_dot_(self.name, self.nop4.__name__))
        return _CS._IN, self.nop4, b, p


class ClipCS3Tuple(_NamedTuple):
    '''3-Tuple C{(start, end, index)} for each edge of a I{clipped}
       path with the C{start} and C{end} points (C{LatLon}) of the
       portion of the edge inside or on the clip box and the C{index}
       (C{int}) of the edge in the original path.
    '''
    _Names_ = (_start_, _end_, 'index')


def clipCS3(points, lowerleft, upperright, closed=False, inull=False):
    '''Clip a path against a rectangular clip box using the
       U{Cohen-Sutherland
       <https://WikiPedia.org/wiki/Cohen-Sutherland_algorithm>} algorithm.

       @arg points: The points (C{LatLon}[]).
       @arg lowerleft: Bottom-left corner of the clip box (C{LatLon}).
       @arg upperright: Top-right corner of the clip box (C{LatLon}).
       @kwarg closed: Optionally, close the path (C{bool}).
       @kwarg inull: Optionally, include null edges if inside (C{bool}).

       @return: Yield a L{ClipCS3Tuple}C{(start, end, index)} for each
                edge of the clipped path.

       @raise ClipError: The B{C{lowerleft}} and B{C{upperright}} corners
                         specify an invalid clip box.

       @raise PointsError: Insufficient number of B{C{points}}.
    '''

    cs = _CS(lowerleft, upperright, name=clipCS3.__name__)
    n, pts = _points2(points, closed, inull)

    i, m = _imdex2(closed, n)
    cmbp = cs.code4(pts[i])
    for i in range(m, n):
        c1, m1, b1, p1 = cmbp
        c2, m2, b2, p2 = cmbp = cs.code4(pts[i])
        if c1 & c2:  # edge outside
            continue

        if not cs.edge(p1, p2):
            if inull:  # null edge
                if not c1:
                    yield ClipCS3Tuple(p1, p1, i)
                elif not c2:
                    yield ClipCS3Tuple(p2, p2, i)
            continue

        for _ in range(5):
            if c1:  # clip p1
                c1, m1, b1, p1 = m1(b1, p1)
            elif c2:  # clip p2
                c2, m2, b2, p2 = m2(b2, p2)
            else:  # inside
                if inull or _neq(p1, p2):
                    yield ClipCS3Tuple(p1, p2, i)
                break
            if c1 & c2:  # edge outside
                break
        else:  # PYCHOK no cover
            raise _AssertionError(_dot_(cs.name, 'for_else'))


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


class _LLi_(LL_):
    '''(INTERNAL) LatLon_ for _SH intersections.
    '''
    __slots__ = _lat_, _lon_, 'classof', 'edge', _name_

    def __init__(self, lat, lon, classof, edge):
        self.lat = lat
        self.lon = lon
        self.classof = classof
        self.edge = edge  # clip edge
        self.name = NN


class _SH(_Named):
    '''(INTERNAL) Sutherland-Hodgman polyon clipping.
    '''
    _cs = ()  # clip corners
    _cw = 0   # counter-/clockwise
    _dx = 0   # clip edge[e] delta lon
    _dy = 0   # clip edge[e] delta lat
    _nc = 0   # len(._cs)
    _x1 = 0   # clip edge[e] lon origin
    _xy = 0   # see .clipedges
    _y1 = 0   # clip edge[e] lat origin

    def __init__(self, corners, name=__name__):
        n, cs = 0, corners
        try:  # check the clip box/region
            n, cs = len2(cs)
            if n == 2:  # make a box
                b, l, t, r = boundsOf(cs, wrap=False)
                cs = LL_(b, l), LL_(t, l), LL_(t, r), LL_(b, r)
            n, cs = points2(cs, closed=True)
            self._cs = cs = cs[:n]
            self._nc = n
            self._cw = isconvex_(cs, adjust=False, wrap=False)
            if not self._cw:
                raise ValueError(_not_convex_)
            if areaOf(cs, adjust=True, radius=1, wrap=True) < EPS:
                raise ValueError('near-zero area')
        except (PointsError, TypeError, ValueError) as x:
            raise ClipError(name, n, cs, txt=str(x))
        self.name = name

    def clip2(self, points, closed, inull):  # MCCABE 17, clip points
        np, pts = _points2(points, closed, inull)
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
            raise ClipError(self.name, ne, self._cs, txt=_too_few_)

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
        # points[i] is located to the right, to the left or on top
        # of the (extended) clip edge from c1 to c2
        d = self._dx * float(p.lat - self._y1) - \
            self._dy * float(p.lon - self._x1)
        # clockwise corners, +1 means points[i] is to the right
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
            raise _AssertionError(_dot_(self.name, self.intersect.__name__))
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
                each edge of the clipped polygon.

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


__all__ += _ALL_DOCS(ClipCS3Tuple, ClipSH3Tuple)

# **) MIT License
#
# Copyright (C) 2018-2020 -- mrJean1 at Gmail -- All Rights Reserved.
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
