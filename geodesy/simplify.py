
# -*- coding: utf-8 -*-

'''This module offers 6 functions to simplify, linearize a path
given as a list, sequence or tuple of LatLon points.  Note, each
function produces different simplified results for the same path.

Function L{simplify1} eliminates points based on edge length.  Function
L{simplify2} slides a pipe over each edge, removing subsequent points
up to the first point outside the pipe.

The functions L{simplifyRDP} and L{simplifyRDPm} use the original,
respectively modified Ramer-Douglas-Peucker (RDP) algorithm, recursively
finding the points farthest from each path edge.  The difference is that
function L{simplifyRDP} exhaustively searches the single, most distant
point in each iteration, while function L{simplifyRDPm} stops at the
first point exceeding the distance tolerance.

Functions L{simplifyVW} and L{simplifyVWm} are based on the original,
respectively modified Visvalingam-Wyatt (VW) method using the area of
the triangle formed by three neigboring points.  The original L{simplifyVW}
method removes only one point per iteration, one with the smallest
triangle area.  Function L{simplifyVWm} removes all points with small
triangle area in every iteration.

Functions L{simplify2}, L{simplifyRDP} and L{simplifyRDPm} provide
keyword I{shortest} to choose between the shortest distance to path
edges/ends or the perpendicular distance to (extended) path edges.

For all functions, keyword I{adjust} scales the longitudinal distance
between two points by the cosine of the mean latitudes.

See:
 - U{https://en.m.wikipedia.org/wiki/Ramer–Douglas–Peucker_algorithm}
 - U{https://hydra.hull.ac.uk/assets/hull:8338/content}
 - U{https://www.cs.ubc.ca/cgi-bin/tr/1992/TR-92-07.pdf}
 - U{http://web.cs.sunyit.edu/~poissad/projects/Curve/about_project.php}
 - U{http://www.bdcc.co.uk/Gmaps/GDouglasPeuker.js}
 - U{https://github.com/mourner/simplify-js/}
 - U{https://github.com/omarestrella/simplify.py/}
 - U{https://bost.ocks.org/mike/simplify/}
 - U{https://pypi.python.org/pypi/visvalingam}
 - U{https://pypi.python.org/pypi/simplification/}
 - U{https://news.ycombinator.com/item?id=4055445}

Tested with 64-bit Python 2.7.13 and 3.6.0 on macOS X 10.12.4 Sierra.

@newfield example: Example, Examples
'''

from datum import R_M
from utils import EPS, len2, radiansPI, wrap180

from math  import cos, degrees, radians

__all__ = ('simplify1', 'simplify2',
           'simplifyRDP', 'simplifyRDPm',
           'simplifyVW', 'simplifyVWm')
__version__ = '17.04.06'


# try:
#     from collections import namedtuple
#     _T2 = namedtuple('_T2', 'ix, h2')
# except ImportError:
#    class _T2(object):
#        ...
# namedtuple not used because (a) values can not
# be updated and (b) it produces PyChecker warning
# "<string>:28: self is not first method argument"
# which can not be suppressed by option --stdlib
class _T2(object):
    '''(INTERNAL) VW 2-tuple (index, area).
    '''
    __slots__ = 'ix', 'h2'

    def __init__(self, ix, h2):
        self.ix = ix
        self.h2 = h2


class _Sy(object):
    '''(INTERNAL) State class.
    '''
    adjust = False
    d2i    = None  # d2i1 or d2i2
    d2xyse = ()
    eps    = EPS
    n      = 0
    pts    = []
    radius = R_M
    r      = {}  # indices or 2-tuples
    s2     = EPS

    def __init__(self, points, tolerance, radius, adjust, shortest):
        n, self.pts = len2(points)
        if n > 0:
            self.n = n
            self.r = {0: True, n-1: True}  # dict to avoid duplicates

        if adjust:
            self.adjust = True

        self.radius = radius

        # tolerance in degrees squared
        s2  = degrees(float(tolerance) / radius)
        s2 *= s2
        self.s2 = s2

        if s2 > EPS and shortest:
            self.eps = s2

        self.d2i = self.d2i2 if shortest else self.d2i1  # PYCHOK false

    def d21(self, s, e):
        # initial edge
        d21, x21, y21 = self.d2xy(s, e)
        self.d2xyse = d21, x21, y21, s, e
        return d21 > self.eps

    def d2i1(self, j, n):
        # consider perpendicular distance
        d21, x21, y21, s, _ = self.d2xyse
        eps, d2xy = self.eps, self.d2xy
        for i in range(j, n):
            d2, x01, y01 = d2xy(s, i)
            if d2 > eps:
                d2  = x01 * y21 + y01 * x21
                d2 *= d2 / d21
                yield d2, i

    def d2i2(self, j, n):
        # consider shortest distance
        d21, x21, y21, s, e = self.d2xyse
        eps, d2xy = self.eps, self.d2xy
        for i in range(j, n):
            d2, x01, y01 = d2xy(s, i)
            if d2 > eps:
                t = x01 * x21 - y01 * y21
                if t > 0:  # triangle sides ...
                    if (t * t) > d21:
                        d2, _, _ = d2xy(e, i)
                    else:  # ... and height
                        d2  = x01 * y21 + y01 * x21
                        d2 *= d2 / d21
                yield d2, i

    def d2xy(self, i, j):
        # edge [i, j]
        p1 = self.pts[i]
        p2 = self.pts[j]

        dx = wrap180(p2.lon - p1.lon)
        if self.adjust:  # lon
            dx *= cos(radiansPI(p1.lat + p2.lat) * 0.5)

        dy = p2.lat - p1.lat
        d2 = dx * dx + dy * dy  # squared!

        return d2, dx, dy

    def h2t(self, i1, i0, i2):
        # compute triangle area (double)
        d21, x21, y21 = self.d2xy(i1, i2)
        if d21 > 0:  # self.eps
            d01, x01, y01 = self.d2xy(i1, i0)
            if d01 > 0:  # self.eps
                h2 = x01 * y21 + y01 * x21
                # triangle height h = sqrt(h2 * h2 / d21)
                # and area = h * sqrt(d21) / 2 == h2 / 2
                return abs(h2)
        return 0

    def points(self, r):
        # simplified points
        return [self.pts[i] for i in sorted(r.keys())]

    def rm(self, m, eps):
        h2t, r = self.h2t, self.r
        # remove one point and redo the triangle area
        # of both neighboring points and remove the
        # latter if the area doesn't exceed epsilon
        r.pop(m)
        for i in (m, m - 1):
            if 0 < i < (len(r) - 1):
                h2 = h2t(r[i-1].ix, r[i].ix, r[i+1].ix)
                if h2 > eps:
                    r[i].h2 = h2
                else:
                    self.rm(i, eps)

    def vw(self):
        # VW initialization
        n, h2t, s2 = self.n, self.h2t, self.s2
        # make list of 2-tuples (ix, h2)
        if n > 2:
            s2 *= 2
            r = [_T2(0, s2 + 1)]
            r.extend(_T2(i, h2t(i-1, i, i+1)) for i in range(1, n-1))
            r.append(_T2(n-1, s2 + 1))
        elif n > 0:
            r = [_T2(i, 0) for i in range(0, n)]
        else:
            r = []

        self.r, self.s2 = r, s2
        return n, r

    def vwr(self, attr):
        # VW results
        pts, r, radius = self.pts, self.r, self.radius

        if attr:  # save triangle area (double)
            for t2 in r:  # convert back to meter
                setattr(pts[t2.ix], attr, radians(t2.h2) * radius)

        # double check for duplicates
        n = len(r)
        r = dict((t2.ix, True) for t2 in r)
        assert len(r) == n
        return r  # as dict


def simplify1(points, distance, radius=R_M, adjust=True):
    '''Basic simplification of a path: eliminate any points closer
       together than the given distance.

       @param points: Path points (LatLons).
       @param distance: Tolerance (meter, same units a radius).
       @keyword radius: Earth radius (meter).
       @keyword adjust: Adjust longitudes by cos(latitutes) (bool).

       @return: Simplified points (list of LatLons).
    '''
    S = _Sy(points, distance, radius, adjust, True)

    n, r = S.n, S.r
    if n > 1:
        s2, d2xy = S.s2, S.d2xy

        i = 0
        for j in range(1, n):
            d2, _, _ = d2xy(i, j)
            if d2 > s2:
                r[j] = True
                i = j

    return S.points(r)


def simplify2(points, band2, radius=R_M, adjust=True, shortest=False):
    '''Pipe simplification of a path: eliminate any points too close
       together or within the given band along an edge.

       @param points: Path points (LatLons).
       @param band2: Half band width (meter, same units a radius).
       @keyword radius: Earth radius (meter).
       @keyword adjust: Adjust longitudes by cos(latitutes) (bool).
       @keyword shortest: Shortest or perpendicular distance (bool).

       @return: Simplified points (list of LatLons).
    '''
    S = _Sy(points, band2, radius, adjust, shortest)

    n, r = S.n, S.r
    if n > 1:
        s2, d21, d2i = S.s2, S.d21, S.d2i

        s, e = 0, 1
        while s < e < n:
            if d21(s, e):
                for d2, i in d2i(e+1, n):
                    if d2 > s2:
                        r[s] = r[i] = True
                        s, e = i, i + 1
                        break  # for loop
                else:
                    r[s] = True  # r[n-1] = True
                    break  # while loop
            else:  # drop points[e]
                e += 1

    return S.points(r)


def simplifyRDP(points, distance, radius=R_M, adjust=True, shortest=False):
    '''Ramer-Douglas-Peucker (RDP) simplification of a path: eliminate
       any points too close together or within the given distance of an
       edge.

       This RDP method exhaustively searches for the point with the
       greatest distance, resulting in worst-case complexity O(n**2)
       where n is the number of points.

       @param points: Path points (LatLons).
       @param distance: Tolerance (meter, same units a radius).
       @keyword radius: Earth radius (meter).
       @keyword adjust: Adjust longitudes by cos(latitutes) (bool).
       @keyword shortest: Shortest or perpendicular distance (bool).

       @return: Simplified points (list of LatLons).
    '''
    S = _Sy(points, distance, radius, adjust, shortest)

    n, r = S.n, S.r
    if n > 1:
        s2, d21, d2i = S.s2, S.d21, S.d2i

        se = [(0, n-1)]
        while se:
            s, e = se.pop()
            if (e - s) > 1:
                if d21(s, e):
                    h2 = h = -1
                    for d2, i in d2i(s+1, e):
                        if d2 > h2:  # nearest so far
                            h2, h = d2, i
                    if h2 > s2:  # split at nearest
                        r[s] = r[h] = True
                        se.append((h, e))
                        se.append((s, h))
                    else:  # all near
                        r[s] = True
                else:  # split halfway
                    i = (e + s) // 2
                    se.append((i, e))
                    se.append((s, i))

    return S.points(r)


def simplifyRDPm(points, distance, radius=R_M, adjust=True, shortest=False):
    '''Modified Ramer-Douglas-Peucker (RDP) simplification of a path:
       eliminate any points too close together or within the given
       distance of an edge.

       This RDP method stops at the first point farther than the given
       distance, reducing the wort-case run time (but producing results
       which usually differ from the original RDP method).

       @param points: Path points (LatLons).
       @param distance: Tolerance (meter, same units a radius).
       @keyword radius: Earth radius (meter).
       @keyword adjust: Adjust longitudes by cos(latitutes) (bool).
       @keyword shortest: Shortest or perpendicular distance (bool).

       @return: Simplified points (list of LatLons).
    '''
    S = _Sy(points, distance, radius, adjust, shortest)

    n, r = S.n, S.r
    if n > 1:
        s2, d21, d2i = S.s2, S.d21, S.d2i

        se = [(0, n-1)]
        while se:
            s, e = se.pop()
            if (e - s) > 1:
                if d21(s, e):
                    for d2, i in d2i(s+1, e):
                        if d2 > s2:
                            r[s] = r[i] = True
                            se.append((i, e))
                            break
                    else:
                        r[s] = True
                else:  # split halfway
                    i = (e + s) // 2
                    se.append((i, e))
                    se.append((s, i))

    return S.points(r)


def simplifyVW(points, area2, radius=R_M, adjust=True, attr=None):
    '''Visvalingam-Wyatt (VW) simplification of a path: eliminate
       any points too close together or with a triangle area below
       the given area tolerance (squared).

       This VW method exhaustively searches for the point with the
       smallest triangle area, resulting in worst-case complexity
       O(n**2) where n is the number of points.

       @param points: Path points (LatLons).
       @param area2: Tolerance (meter, same units a radius).
       @keyword radius: Earth radius (meter).
       @keyword adjust: Adjust longitudes by cos(latitutes) (bool).
       @keyword attr: Points attribute save area value (string).

       @return: Simplified points (list of LatLons).
    '''
    S = _Sy(points, area2, radius, adjust, False)

    n, r = S.vw()
    if n > 2:
        rm, s2 = S.rm, S.s2

        # remove any points too close
        # together or with zero area
        for i in range(len(r)-2, 0, -1):
            if not r[i].h2:
                rm(i, 0)

        # keep removing the point with the smallest
        # area until latter exceeds the tolerance
        while len(r) > 2:
            m, m2 = 0, s2 + 1
            for i in range(1, len(r)-1):
                h2 = r[i].h2
                if h2 < m2:
                    m, m2 = i, h2
            if m2 > s2:
                break
            rm(m, 0)

    return S.points(S.vwr(attr))


def simplifyVWm(points, area2, radius=R_M, adjust=True, attr=None):
    '''Modified Visvalingam-Wyatt (VW) simplification of a path:
       eliminate any points too close together or with a triangle
       area below the given area tolerance (squared).

       This VW method repeatedly removes all points with a triangle
       area below the tolerance, significantly reducing the run time
       (but producing results which usually differ from the original
       VW method).

       @param points: Path points (LatLons).
       @param area2: Tolerance (meter, same units a radius).
       @keyword radius: Earth radius (meter).
       @keyword adjust: Adjust longitudes by cos(latitutes) (bool).
       @keyword attr: Points attribute save area value (string).

       @return: Simplified points (list of LatLons).
    '''
    S = _Sy(points, area2, radius, adjust, False)

    n, r = S.vw()
    if n > 2:
        rm, s2 = S.rm, S.s2

        # keep removing all points with
        # an area below the tolerance
        m = n + 1
        while m > n > 2:
            m = n = len(r)
            for i in range(m-2, 0, -1):
                if i < n and r[i].h2 <= s2:
                    rm(i, s2)
                    n = len(r)

    return S.points(S.vwr(attr))

# **) MIT License
#
# Copyright (C) 2016-2017 -- mrJean1 at Gmail dot com
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
