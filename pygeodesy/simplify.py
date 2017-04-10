
# -*- coding: utf-8 -*-

'''Several functions to simplify or linearize a path given as a list,
sequence or tuple of LatLon points.

Each of the simplify functions is based on a different algorithm and
produces different simplified results in (very) different run times
for the same path of LatLon points.

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
method removes only a single point per iteration, while the modified
L{simplifyVWm} removes all points with areas not exceeding the
tolerance in each iteration.

Functions L{simplify2}, L{simplifyRDP} and L{simplifyRDPm} provide
keyword I{shortest} to select the computation of the distance between
a point and a path edge.  If True, use the shortest distance to path
edge or path end points, False use the perpendicular distance to the
extended path edge.

For all functions, keyword I{adjust} scales the longitudinal distance
between two points by the cosine of the mean of the latitudes.

See:
 - U{https://en.m.wikipedia.org/wiki/Ramer-Douglas-Peucker_algorithm}
 - U{https://hydra.hull.ac.uk/assets/hull:8338}
 - U{https://www.cs.ubc.ca/cgi-bin/tr/1992/TR-92-07.pdf}
 - U{http://web.cs.sunyit.edu/~poissad/projects/Curve/about_project.php}
 - U{http://www.bdcc.co.uk/Gmaps/GDouglasPeuker.js}
 - U{https://github.com/mourner/simplify-js/}
 - U{https://github.com/omarestrella/simplify.py/}
 - U{https://bost.ocks.org/mike/simplify/}
 - U{https://pypi.python.org/pypi/visvalingam}
 - U{https://pypi.python.org/pypi/simplification/}
 - U{https://news.ycombinator.com/item?id=4055445}

Tested with 64-bit Python 2.6.9, 2.7.13, 3.5.3 and 3.6.0 on macOS
10.12.3 and 10.12.4 Sierra.  On macOS, Python 3 runs the simplify
tests nearly 2x faster than Python 2.

@newfield example: Example, Examples
'''

from datum import R_M
from utils import EPS, len2, radiansPI, wrap180

from math  import cos, degrees, radians

__all__ = ('simplify1', 'simplify2',
           'simplifyRDP', 'simplifyRDPm',
           'simplifyVW', 'simplifyVWm')
__version__ = '17.04.09'


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
    '''(INTERNAL) Simplify state.
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
        '''New state.
        '''
        n, self.pts = len2(points)
        if n > 0:
            self.n = n
            self.r = {0: True, n-1: True}  # dict to avoid duplicates

        if adjust:
            self.adjust = True

        self.radius = radius

        # tolerance converted to degrees squared
        s2  = degrees(float(tolerance) / radius)
        s2 *= s2
        self.s2 = s2

        if s2 > EPS and shortest:
            self.eps = s2

        self.d2i = self.d2i2 if shortest else self.d2i1  # PYCHOK false

    def d21(self, s, e):
        '''Sets path edge or line thru points[s] to [e].
        '''
        d21, x21, y21 = self.d2xy(s, e)
        self.d2xyse = d21, x21, y21, s, e
        return d21 > self.eps

    def d2i1(self, j, n):
        '''Yields perpendicular distance for all points[j..n]
           to the path edge or line thru points[s] to [e].
        '''
        d21, x21, y21, s, _ = self.d2xyse
        eps, d2xy = self.eps, self.d2xy
        for i in range(j, n):
            d2, x01, y01 = d2xy(s, i)
            if d2 > eps:
                d2  = x01 * y21 + y01 * x21
                d2 *= d2 / d21
                yield d2, i

    def d2i2(self, j, n):
        '''Yields shortest distance for all points[j..n]
           to the path edge or line thru points[s] to [e].
        '''
        d21, x21, y21, s, e = self.d2xyse
        eps, d2xy = self.eps, self.d2xy
        for i in range(j, n):
            # distance points[i] to [s]
            d2, x01, y01 = d2xy(s, i)
            if d2 > eps:
                t = x01 * x21 - y01 * y21
                if t > 0:
                    if (t * t) > d21:
                        # distance points[i] to [e]
                        d2, _, _ = d2xy(e, i)
                    else:  # perpendicular distance
                        d2  = x01 * y21 + y01 * x21
                        d2 *= d2 / d21
                yield d2, i

    def d2xy(self, i, j):
        '''Returns points[i] to [j] deltas.
        '''
        p1 = self.pts[i]
        p2 = self.pts[j]

        dx = wrap180(p2.lon - p1.lon)
        dy = wrap180(p2.lat - p1.lat)

        if self.adjust:  # scale lon
            dx *= cos(radiansPI(p1.lat + p2.lat) * 0.5)

        d2 = dx * dx + dy * dy  # squared!
        return d2, dx, dy

    def h2t(self, i1, i0, i2):
        '''Computes the Visvalingam-Wyatt triangular area,
           points[i1] to [i2] form the base and points[i0]
           is the top of the triangle.
        '''
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
        '''Returns the list of simplified points.
        '''
        return [self.pts[i] for i in sorted(r.keys())]

    def rm1(self, m, eps):
        '''Eliminates one Visvalingam-Wyatt point and recompute
           the trangular area of each neighboring point, but
           remove any of those too until its recomputed area
           exceeds the tolerance.
        '''
        h2t, r = self.h2t, self.r

        r.pop(m)
        for n in (m, m - 1):
            while 0 < n < (len(r) - 1):
                h2 = h2t(r[n-1].ix, r[n].ix, r[n+1].ix)
                if h2 > eps:
                    r[n].h2 = h2
                    break  # while
                else:
                    r.pop(n)

    def rm2(self, eps):
        '''Eliminates all Visvalingam-Wyatt points with a
           triangular area not exceeding the tollerance.
        '''
        r, rm1 = self.r, self.rm1

        i = len(r) - 1
        while i > 1:
            i -= 1
            if r[i].h2 <= eps:
                rm1(i, eps)
                i = min(i, len(r) - 1)

    def vw(self):
        '''Initializes Visvalingam-Wyatt as list of 2-tuples
           (ix, h2) where ix is the points[] index and h2
           twice the triangular area of that points[ix].
        '''
        n, h2t, s2 = self.n, self.h2t, self.s2

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
        return len(r), r

    def vwr(self, attr):
        '''Returns Visvalingam-Wyatt results as dict,
           optionally including the triangular area
           (in meters) for each simplified point.
        '''
        pts, r, radius, s2 = self.pts, self.r, self.radius, self.s2

        # double check the minimal triangular area
        assert min(t2.h2 for t2 in r) > s2

        if attr:  # return triangular area (double)
            r[0].h2 = r[-1].h2 = 0
            for t2 in r:  # convert back to meter
                setattr(pts[t2.ix], attr, radians(t2.h2) * radius)

        # double check for duplicates
        n = len(r)
        r = dict((t2.ix, True) for t2 in r)
        assert len(r) == n
        return r  # as dict


def simplify1(points, distance, radius=R_M, adjust=True):
    '''Basic simplification of a path of LatLon points.

       Eliminate any points closer together than the given
       distance tolerance.

       @param points: Path points (LatLons).
       @param distance: Tolerance (meter, same units a radius).
       @keyword radius: Earth radius (meter).
       @keyword adjust: Adjust longitudes (bool).

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
    '''Pipe simplification of a path of LatLon points.

       Eliminate any points too close together or within the given
       band tolerance along an edge.

       @param points: Path points (LatLons).
       @param band2: Half band width (meter, same units a radius).
       @keyword radius: Earth radius (meter).
       @keyword adjust: Adjust longitudes (bool).
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
    '''Ramer-Douglas-Peucker (RDP) simplification of a path of
       LatLon points.

       Eliminate any points too close together or closer to an
       edge than the given distance tolerance.

       This RDP method exhaustively searches for the point with
       the largest distance, resulting in worst-case complexity
       O(n**2) where n is the number of points.

       @param points: Path points (LatLons).
       @param distance: Tolerance (meter, same units a radius).
       @keyword radius: Earth radius (meter).
       @keyword adjust: Adjust longitudes (bool).
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
    '''Modified Ramer-Douglas-Peucker (RDP) simplification of a path
       of LatLon points.

       Eliminate any points too close together or closer to an edge
       then the given distance tolerance.

       This RDP method stops at the first point farther than the
       given distance tolerance, significantly reducing the run time
       (but producing results different from the original RDP method).

       @param points: Path points (LatLons).
       @param distance: Tolerance (meter, same units a radius).
       @keyword radius: Earth radius (meter).
       @keyword adjust: Adjust longitudes (bool).
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
    '''Visvalingam-Wyatt (VW) simplification of a path of LatLon
       points.

       Eliminate any points too close together or with a triangular
       area not exceeding the given area tolerance squared.

       This VW method exhaustively searches for the single point
       with the smallest triangular area, resulting in worst-case
       complexity O(n**2) where n is the number of points.

       @param points: Path points (LatLons).
       @param area2: Tolerance (meter, same units a radius).
       @keyword radius: Earth radius (meter).
       @keyword adjust: Adjust longitudes (bool).
       @keyword attr: Points attribute save area value (string).

       @return: Simplified points (list of LatLons).
    '''
    S = _Sy(points, area2, radius, adjust, False)

    n, r = S.vw()
    if n > 2:
        # remove any points too close or
        # with a zero triangular area
        S.rm2(0)

        # keep removing the point with the smallest
        # area until latter exceeds the tolerance
        while len(r) > 2:
            m, m2 = 0, S.s2 + 1
            for i in range(1, len(r)-1):
                h2 = r[i].h2
                if h2 < m2:
                    m, m2 = i, h2
            if m2 > S.s2:
                break
            S.rm1(m, 0)

    return S.points(S.vwr(attr))


def simplifyVWm(points, area2, radius=R_M, adjust=True, attr=None):
    '''Modified Visvalingam-Wyatt (VW) simplification of a path of
       LatLon points.

       Eliminate any points too close together or with a triangular
       area not exceeding the given area tolerance squared.

       This VW method removes all points with a triangular area
       below the tolerance per iteration, significantly reducing the
       run time (but producing results different from the original
       VW method).

       @param points: Path points (LatLons).
       @param area2: Tolerance (meter, same units a radius).
       @keyword radius: Earth radius (meter).
       @keyword adjust: Adjust longitudes (bool).
       @keyword attr: Attribute to save the area value (string).

       @return: Simplified points (list of LatLons).
    '''
    S = _Sy(points, area2, radius, adjust, False)

    n, _ = S.vw()
    if n > 2:
        # remove all points with an area
        # not exceeding the tolerance
        S.rm2(S.s2)

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
