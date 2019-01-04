
# -*- coding: utf-8 -*-

u'''Functions to I{clip} a path or polygon of C{LatLon} points.

@newfield example: Example, Examples
'''

from fmath  import EPS, fsum_, len2
from lazily import _ALL_LAZY
from points import _imdex2, boundsOf, isclockwise, isconvex_, \
                   LatLon_ as LL_
from utily  import points2

__all__ = _ALL_LAZY.clipy
__version__ = '19.01.02'


def _eq(p1, p2):  # near-equal points
    return not _neq(p1, p2)


def _neq(p1, p2):  # not near-equal points
    return abs(p1.lat - p2.lat) > EPS or \
           abs(p1.lon - p2.lon) > EPS


class _CS(object):
    '''(INTERNAL) Cohen-Sutherland line clipping.
    '''
    # single-bit codes
    _IN   = 0  # inside clip box
    _YMIN = 1  # below lowerleft.lat
    _YMAX = 2  # above upperright.lat
    _XMIN = 4  # left of lowerleft.lon
    _XMAX = 8  # right of upperright.lon

    def __init__(self, lowerleft, upperright):
        self._ymin, self._ymax = lowerleft.lat, upperright.lat
        self._xmin, self._xmax = lowerleft.lon, upperright.lon
        if self._xmin > self._xmax or self._ymin > self._ymax:
            raise ValueError('invalid %s: %r to %r' % (
                             'clip box', lowerleft, upperright))
        self._y1 = self._dy = 0
        self._x1 = self._dx = 0

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
#       raise AssertionError('clipCS3.clip2')

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

    def nop4(self, b, p):
        if p:  # should never get here
            raise AssertionError('clipCS3.nop4')
        return _CS._IN, self.nop4, b, p


def clipCS3(points, lowerleft, upperright, closed=False, inull=False):  # MCCABE 25
    '''Clip a path against a rectangular clip box using the
       U{Cohen-Sutherland
       <http://WikiPedia.org/wiki/Cohen-Sutherland_algorithm>} algorithm.

       @param points: The points (C{LatLon}[]).
       @param lowerleft: Bottom-left corner of the clip box (C{LatLon}).
       @param upperright: Top-right corner of the clip box (C{LatLon}).
       @keyword closed: Optionally, close the path (C{bool}).
       @keyword inull: Optionally, include null edges if inside (C{bool}).

       @return: Yield a 3-tuple (start, end, index) for each edge
                of the clipped path with the start and end points
                C{LatLon} of the portion of the edge inside or on the
                clip box and the index C{int} of the edge in the
                original path.

       @raise ValueError: The I{lowerleft} corner is not below and/or
                          not to the left of the I{upperright} corner.
    '''

    cs = _CS(lowerleft, upperright)
    n, points = points2(points, closed=closed)

    i, m = _imdex2(closed, n)
    cmbp = cs.code4(points[i])
    for i in range(m, n):
        c1, m1, b1, p1 = cmbp
        c2, m2, b2, p2 = cmbp = cs.code4(points[i])
        if c1 & c2:  # edge outside
            continue

        if not cs.edge(p1, p2):
            if inull:  # null edge
                if not c1:
                    yield p1, p1, i
                elif not c2:
                    yield p2, p2, i
            continue

        for _ in range(5):
            if c1:  # clip p1
                c1, m1, b1, p1 = m1(b1, p1)
            elif c2:  # clip p2
                c2, m2, b2, p2 = m2(b2, p2)
            else:  # inside
                if inull or _neq(p1, p2):
                    yield p1, p2, i
                break
            if c1 & c2:  # edge outside
                break
        else:  # should never get here
            raise AssertionError('clipCS3.for_')


class _LLi_(LL_):
    '''(INTERNAL) LatLon_ for _SH intersections.
    '''
    __slots__ = 'lat', 'lon', 'classof', 'edge', 'name'

    def __init__(self, lat, lon, classof, edge):
        self.lat = lat
        self.lon = lon
        self.classof = classof
        self.edge = edge  # clip edge
        self.name = ''


class _SH(object):
    '''(INTERNAL) Sutherland-Hodgman polyon clipping.
    '''
    def __init__(self, corners):
        n = ''
        try:
            n, cs = len2(corners)
            if n == 2:  # make a box
                b, l, t, r = boundsOf(cs, wrap=False)
                cs = LL_(b, l), LL_(t, l), LL_(t, r), LL_(b, r)
            n, cs = points2(cs, closed=True)
            self._corners = cs = cs[:n]
            self._nc = n
            self._cw = 1 if isclockwise(cs, adjust=False, wrap=False) else -1
            if self._cw != isconvex_(cs, adjust=False, wrap=False):
                raise ValueError
        except ValueError:
            raise ValueError('invalid %s[%s]: %r' % ('corners', n, corners))
        self._clipped = self._points = []

    def append(self, p, inull):  # save a clipped point
        if inull or (not self._clipped) or _neq(p, self._clipped[-1]):
            self._clipped.append(p)

    def clip(self, points, inull, closed):  # MCCABE 19, clip points
        np, self._points = len2(points)
        if np > 0:  # initial points, opened
            p = self._points[0]
            while np > 1 and _eq(self._points[np - 1], p):
                np -= 1
        if np < 3:
            raise ValueError('too few %s: %s' % ('points', np))

        no = ni = True  # all out- or inside?
        for e in self.clips():
            # clip the points, closed
            d2, p2 = self.dot2(np - 1)
            for i in range(np):
                d1, p1 = d2, p2
                d2, p2 = self.dot2(i)
                if d1 < 0:  # p1 inside, ...
                    # self.append(p1, False, e)
                    if d2 < 0:  # ... p2 inside
                        self.append(p2, inull)
                    else:  # ... p2 outside
                        p = self.intersect(p1, p2, e)
                        self.append(p, inull)
                        if d2 > 0:
                            no = False
                elif d2 < 0:  # p1 out-, p2 inside
                    p = self.intersect(p1, p2, e)
                    self.append(p, inull)
                    self.append(p2, inull)
                    if d1 > 0:
                        no = ni = False
                elif d1 > 0:
                    ni = False
            if self._clipped:  # replace points
                self._points = self._clipped
                self._clipped = []
                np = len(self._points)

        # no is True iff all points are on or at one
        # side (left or right) of each clip edge,
        # ni is True iff all points are on or on the
        # right side (i.e. inside) of all clip edges
        if no and not ni:
            self._points = ()
            np = 0
        elif np > 1:
            p = self._points[0]
            if closed:  # close clipped polygon
                if _neq(self._points[np - 1], p):
                    self._points.append(p)
                    np += 1
            elif not inull:  # open clipped polygon
                while np > 0 and _eq(self._points[np - 1], p):
                    np -= 1
        return np

    def clipped2(self, i):  # return (clipped point[i], edge)
        p = self._points[i]
        if isinstance(p, _LLi_):  # intersection point
            return p.classof(p.lat, p.lon), p.edge
        else:  # original point
            return p, 0

    def clips(self):  # yield clip edge index
        c2 = self._corners[self._nc - 1]
        for i in range(self._nc):
            c1, c2 = c2, self._corners[i]
            self._y1, self._dy = c1.lat, float(c2.lat - c1.lat)
            self._x1, self._dx = c1.lon, float(c2.lon - c1.lon)
            if abs(self._dx) > EPS or abs(self._dy) > EPS:
                self._xy = self._y1 * self._dx - self._x1 * self._dy
                yield i + 1

    def dot2(self, i):  # dot product of points[i] to the current
        # clip corner c1 and clip edge c1 to c2, indicating whether
        # points[i] is located to the right, to the left or on the
        # (extended) clip edge from c1 to c2
        p = self._points[i]
        d = self._dx * float(p.lat - self._y1) - \
            self._dy * float(p.lon - self._x1)
        # clockwise corners, +1 means points[i] is to the right
        # of, -1 means on the left of, 0 means on edge c1 to c2
        d = self._cw if d > 0 else (-self._cw if d < 0 else 0)
        return d, p

    def intersect(self, p1, p2, edge):  # compute intersection
        # of polygon edge p1 to p2 and the current clip edge,
        # where p1 and p2 are known to NOT be located on the
        # same side of or on the current clip edge
        # <http://StackOverflow.com/questions/563198/
        #       how-do-you-detect-where-two-line-segments-intersect>
        fy = float(p2.lat - p1.lat)
        fx = float(p2.lon - p1.lon)
        fp = fy * self._dx - fx * self._dy
        if abs(fp) < EPS:
            raise AssertionError('clipSH.intersect')
        h = fsum_(self._xy, -p1.lat * self._dx, p1.lon * self._dy) / fp
        y = p1.lat + h * fy
        x = p1.lon + h * fx
        return _LLi_(y, x, p1.classof, edge)


def clipSH(points, corners, inull=False, closed=False):
    '''Clip a polygon against a clip region or box using the
       U{Sutherland_Hodgman
       <http://WikiPedia.org/wiki/Sutherland_Hodgman_algorithm>} algorithm.

       @param points: The polygon points (C{LatLon}[]).
       @param corners: Three or more points defining a convex clip
                       region (C{LatLon}[]) or two points to specify
                       a rectangular clip box.
       @keyword inull: Optionally, include null edges (C{bool}).
       @keyword closed: Close the clipped points (C{bool}).

       @return: Yield the clipped points (C{LatLon}[]).

       @raise ValueError: Insufficient number of I{points} or the
                          I{corners} specify a polar, zero-area,
                          non-convex or otherwise invalid clip region.
    '''
    sh = _SH(corners)
    n = sh.clip(points, inull, closed)
    for i in range(n):
        yield sh.clipped2(i)[0]


def clipSH3(points, corners, inull=False, closed=False):
    '''Clip a polygon against a clip region or box using the
       U{Sutherland_Hodgman
       <http://WikiPedia.org/wiki/Sutherland_Hodgman_algorithm>} algorithm.

       @param points: The polygon points (C{LatLon}[]).
       @param corners: Three or more points defining a convex clip
                       region (C{LatLon}[]) or two points to specify
                       a rectangular clip box.
       @keyword inull: Optionally, include null edges (C{bool}).
       @keyword closed: Close the clipped points (C{bool}).

       @return: Yield a 3-tuple (start, end, original) for each edge
                of the clipped polygon.  The start and end points
                C{LatLon} of the portion of the edge inside or on
                the clip region.  The original C{bool} indicates
                whether the edge is part of the original polygon or
                part of the clip region.

       @raise ValueError: Insufficient number of I{points} or the
                          I{corners} specify a polar, zero-area,
                          non-convex or otherwise invalid clip region.
    '''
    sh = _SH(corners)
    n = sh.clip(points, inull, closed)
    if n > 0:
        p2, e2 = sh.clipped2(0)
        for i in range(1, n):
            p1, e1 = p2, e2
            p2, e2 = sh.clipped2(i)
            yield p1, p2, not bool(e1 == e2 and e1 and e2)


# **) MIT License
#
# Copyright (C) 2018-2019 -- mrJean1 at Gmail dot com
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
