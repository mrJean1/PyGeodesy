
# -*- coding: utf-8 -*-

u'''Simplify or linearize a path.

Each of the I{simplify} functions is based on a different algorithm and
produces different, simplified results in (very) different run times for
the same path of C{LatLon} points.

Function L{simplify1} eliminates points based on edge lengths shorter
than a given tolerance.

The functions L{simplifyRDP} and L{simplifyRDPm} use the original,
respectively modified I{Ramer-Douglas-Peucker} (RDP) algorithm, iteratively
finding the point farthest from each path edge.  The difference is that
function L{simplifyRDP} exhaustively searches the most distant point in
each iteration, while modified L{simplifyRDPm} stops at the first point
exceeding the distance tolerance.

Function L{simplifyRW} use the I{Reumann-Witkam} method, sliding a "pipe"
over each path edge, removing all subsequent points within, closer than
the pipe radius up to the first point outside the pipe.

Functions L{simplifyVW} and L{simplifyVWm} are based on the original,
respectively modified I{Visvalingam-Whyatt} (VW) method using the area of
the triangle formed by three neigboring points.  Function L{simplifyVW}
removes only a single point per iteration, while modified L{simplifyVWm}
eliminates in each iteration all points with a triangular area not
exceeding the tolerance.

Functions L{simplifyRDP}, L{simplifyRDPm} and L{simplifyRW} provide
keyword argument I{shortest} to specify of the distance between a point
and a path edge.  If C{True}, use the I{shortest} distance to the path
edge or edge points, otherwise use the I{perpendicular} distance to
the extended line through both path edge points.

Keyword argument B{C{radius}} of all fuctions is set to the mean earth
radius in C{meter}.  Other units can be used, provided that the radius
and tolerance are always specified in the same units.

Use keyword argument C{B{indices}=True} in any function to return a
list of simplified point I{indices} instead of the simplified points.
The first and last index are always the first and last original index.

Finally, any additional keyword arguments B{C{options}} to all functions
are passed thru to function L{pygeodesy.equirectangular_} to specify the
distance approximation.

To process C{NumPy} arrays containing rows of lat-, longitude and
possibly other values, use class L{Numpy2LatLon} to wrap the C{NumPy}
array into I{on-the-fly-LatLon} points.  Pass the L{Numpy2LatLon}
instance to any I{simplify} function and the returned result will be
a C{NumPy} array containing the simplified subset, a partial copy of
the original C{NumPy} array.  Use keyword argument C{B{indices}=True}
to return a list of array row indices inlieu of the simplified array
subset.

See:
 - U{https://Bost.Ocks.org/mike/simplify}
 - U{https://WikiPedia.org/wiki/Ramer-Douglas-Peucker_algorithm}
 - U{https://www.ScienceDirect.com/science/article/pii/S0098300402000092}
 - U{https://hydra.Hull.ac.UK/resources/hull:8338}
 - U{https://psimpl.SourceForge.net/reumann-witkam.html}
 - U{https://www.CS.UBC.Ca/cgi-bin/tr/1992/TR-92-07.pdf}
 - U{https://GitHub.com/FlorianWilhelm/gps_data_with_python}
 - U{https://www.BDCC.co.UK/Gmaps/GDouglasPeuker.js}
 - U{https://GitHub.com/mourner/simplify-js}
 - U{https://GitHub.com/OmarEstrella/simplify.py}
 - U{https://PyPI.org/project/rdp}
 - U{https://PyPI.org/project/visvalingam}
 - U{https://PyPI.org/project/simplification}
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

# from pygeodesy.basics import len2  # from .fmath
from pygeodesy.constants import EPS, R_M, _1_0
from pygeodesy.errors import _AttributeError, _ValueError
from pygeodesy.fmath import len2, sqrt0
from pygeodesy.formy import equirectangular_
from pygeodesy.interns import _small_, _too_
from pygeodesy.iters import isNumpy2, isTuple2
# from pygeodesy.lazily import _ALL_LAZY  # from .units
from pygeodesy.units import _ALL_LAZY, _1mm, Radius_

from math import degrees, fabs, radians

__all__ = _ALL_LAZY.simplify
__version__ = '23.03.22'


# try:
#     from collections import namedtuple
#     _T2 = namedtuple('_T2', 'ix, h2')
# except ImportError:
#     class _T2(object):
#         ...
# namedtuple (and .named._NamedTuple) can not be
# used because (a) values can not be updated and
# (b) it produces PyChecker warning "<string>:28:
# self is not first method argument" which can't
# be suppressed with command line option --stdlib
class _T2(object):
    '''(INTERNAL) VW 2-tuple (index, area).
    '''
    # __slots__ are no longer space savers, see
    # the comments at the class .points.LatLon_
    # __slots__ = ('ix', 'h2')

    def __init__(self, ix, h2):
        self.ix = ix
        self.h2 = h2


class _Sy(object):
    '''(INTERNAL) Simplify state.
    '''
    d2i2    = None  # d2iP2 or d2iS2
    d2yxse5 = ()    # 5-tuple
    eps     = EPS   # system epsilon
    indices = False
    n       = 0
    options = {}
    pts     = []
    radius  = R_M   # mean earth radius
    r       = {}    # RDP indices or VW 2-tuples
    s2      = EPS   # tolerance squared
    s2e     = EPS   # sentinel
    subset  = None  # isNumpy2 or isTuple2

    def __init__(self, points, tolerance, radius, shortest,
                                          indices, **options):
        '''New C{Simplify} state.
        '''
        n, self.pts = len2(points)
        if n > 0:
            self.n = n
            self.r = {0: True, n-1: True}  # dict to avoid duplicates

        if isNumpy2(points) or isTuple2(points):  # NOT self.pts
            self.subset = points.subset

        if indices:
            self.indices = True

        if radius:
            self.radius = Radius_(radius, low=self.eps)
        elif self.radius < self.eps:
            raise _ValueError(radius=radius, txt=_too_(_small_))

        if options:
            self.options = options

        # tolerance converted to degrees squared
        self.s2 = degrees(tolerance / self.radius)**2
        if min(self.s2, tolerance) < self.eps:
            raise _ValueError(tolerance=tolerance, txt=_too_(_small_))
        self.s2e = self.s2 + 1  # sentinel

        # compute either the shortest or perpendicular distance
        self.d2i2 = self.d2iS2 if shortest else self.d2iP2  # PYCHOK attr

    def d21(self, s, e):
        '''Set path edge or line thru points[s] to -[e].
        '''
        d21, y21, x21, _ = self.d2yxu4(s, e)
        self.d2yxse5 = d21, y21, x21, s, e
        return d21 > self.eps

    def d2ih2(self, n, m, brk):
        '''Find the tallest distance among all points[n..m]
           to points[s] to -[e] exceeding the tolerance.
        '''
        _, _, _, s, _ = self.d2yxse5
        eps, _d2yxu4 = self.eps, self.d2yxu4
        t2, t = self.s2, 0  # tallest
        for i in range(n, m):
            d2, _, _, _ = _d2yxu4(s, i)
            if d2 > t2:
                t2, t = d2, i
                if brk and d2 > eps:
                    break
        return t2, t

    def d2iP2(self, n, m, brk):
        '''Find the tallest I{perpendicular} distance among all
           points[n..m] to the path edge or line thru points[s]
           to -[e] exceeding the tolerance.
        '''
        d21, y21, x21, s, _ = self.d2yxse5
        eps, _d2yxu4 = self.eps, self.d2yxu4
        t2, t = self.s2, 0  # tallest
        for i in range(n, m):
            d2, y01, x01, _ = _d2yxu4(s, i)
            if d2 > eps:
                # perpendicular distance squared
                d2 = (y01 * x21 - x01 * y21)**2 / d21
                if d2 > t2:
                    t2, t = d2, i
                    if brk:
                        break
        return t2, t

    def d2iS2(self, n, m, brk):
        '''Find the tallest I{shortest} distance among all
           points[n..m] to the path edge or line thru points[s]
           to -[e] exceeding the tolerance.
        '''
        # point (x, y) on axis rotated by angle a ccw:
        #   x' = y * sin(a) + x * cos(a)
        #   y' = y * cos(a) - x * sin(a)
        #
        # distance (w) along and perpendicular (h) to
        # a line thru point (dx, dy) and the origin:
        #   w = (y * dy + x * dx) / hypot(dx, dy)
        #   h = (y * dx - x * dy) / hypot(dx, dy)

        d21, y21, x21, s, e = self.d2yxse5
        eps, _d2yxu4 = self.eps, self.d2yxu4
        t2, t = self.s2, 0  # tallest
        for i in range(n, m):
            # distance points[i] to -[s]
            d2, y01, x01, _ = _d2yxu4(s, i)
            if d2 > eps:
                w = y01 * y21 + x01 * x21
                if w > 0:
                    if w < d21:
                        # perpendicular distance squared
                        d2 = (y01 * x21 - x01 * y21)**2 / d21
                    else:  # distance points[i] to -[e]
                        d2, _, _, _ = _d2yxu4(e, i)
                if d2 > t2:
                    t2, t = d2, i
                    if brk:
                        break
        return t2, t

    def d2yxu4(self, i, j):
        '''Return distance I{squared}, points[i] to -[j] deltas
           and the (longitudinal) unrollment.
        '''
        p1, p2= self.pts[i], self.pts[j]
        return equirectangular_(p1.lat, p1.lon,
                                p2.lat, p2.lon, **self.options)

    def h2t(self, i1, i0, i2):
        '''Compute the Visvalingam-Whyatt triangular area,
           points[i0] is the top and points[i1] to -[i2]
           form the base of the triangle.
        '''
        d21, y21, x21 , _= self.d2yxu4(i1, i2)
        if d21 > self.eps:
            d01, y01, x01, _ = self.d2yxu4(i1, i0)
            if d01 > self.eps:
                h2 = fabs(y01 * x21 - x01 * y21)
                # triangle height h = h2 / sqrt(d21) and
                # the area = h * sqrt(d21) / 2 == h2 / 2
                return h2  # double triangle area
        return 0

    def points(self, r):
        '''Return the list of simplified points or indices.
        '''
        r = sorted(r.keys())
        if self.indices:
            return list(r)
        elif self.subset:
            return self.subset(r)
        else:
            return [self.pts[i] for i in r]

    def rdp(self, modified):
        '''Ramer-Douglas-Peucker (RDP) simplification of a
           path of C{LatLon} points.

           @arg modified: Use modified RDP (C{bool}).
        '''
        n, r = self.n, self.r
        if n > 1:
            s2,    _d21   = self.s2,   self.d21
            _d2i2, _d2ih2 = self.d2i2, self.d2ih2

            se = [(0, n-1)]
            _a = se.append
            _p = se.pop
            while se:
                s, e = _p()
                s1 = s + 1
                if e > s1:
                    if _d21(s, e):  # points[] to edge [s, e]
                        d2, i = _d2i2(s1, e, modified)
                    else:  # points[] to point [s]
                        d2, i = _d2ih2(s1, e, modified)
                    if d2 > s2 and i > 0:
                        r[i] = True
                        _a((i, e))
                        if not modified:
                            _a((s, i))
                    r[s] = True

        return self.points(r)

    def rm1(self, m, tol):
        '''Eliminate one Visvalingam-Whyatt point and recomputes
           the trangular area of both neighboring points, but
           removes those too unless the recomputed area exceeds
           the tolerance.
        '''
        _h2t, r = self.h2t, self.r

        r.pop(m)
        for n in (m, m - 1):
            while 0 < n < (len(r) - 1):
                h2 = _h2t(r[n-1].ix, r[n].ix, r[n+1].ix)
                if h2 > tol:
                    r[n].h2 = h2
                    break  # while
                else:
                    r.pop(n)

    def rm2(self, tol):
        '''Eliminate all Visvalingam-Whyatt points with a
           triangular area not exceeding the tolerance.
        '''
        r, _rm1 = self.r, self.rm1

        i = len(r) - 1
        while i > 1:
            i -= 1
            if r[i].h2 <= tol:
                _rm1(i, tol)
                i = min(i, len(r) - 1)

    def vwn(self):
        '''Initialize Visvalingam-Whyatt as list of 2-tuples
           _T2(ix, h2) where ix is the points[] index and h2
           is the triangular area I{(times 2)} of that point.
        '''
        n, _h2t, s2e, T2 = self.n, self.h2t, self.s2e, _T2

        self.r = r = []
        if n > 2:
            r[:] = [T2(i, _h2t(i-1, i, i+1)) for i in range(1, n-1)]
        if n > 1:
            r.append(T2(n-1, s2e))
        if n > 0:
            r.insert(0, T2(0, s2e))
        return len(r)

    def vwr(self, attr):
        '''Return the Visvalingam-Whyatt results, optionally
           including the triangular area (in meters) as
           attribute attr to each simplified point.
        '''
        pts, r = self.pts, self.r

        # double check the minimal triangular area
        assert min(t2.h2 for t2 in r) > self.s2 > 0

        if attr:  # return the trangular area (actually
            # the sqrt of double the triangular area)
            # converted back from degrees to meter
            if isNumpy2(pts):
                raise _AttributeError(attr=attr)
            m = radians(_1_0) * self.radius
            r[0].h2 = r[-1].h2 = 0  # zap sentinels
            for t2 in r:  # convert back to meter
                setattr(pts[t2.ix], attr, sqrt0(t2.h2) * m)

        # double check for duplicates
        n = len(r)
        r = dict((t2.ix, True) for t2 in r)
        assert len(r) == n
        return self.points(r)


def simplify1(points, distance=_1mm, radius=R_M, indices=False, **options):
    '''Basic simplification of a path of C{LatLon} points.

       Eliminates any points closer together than the given I{distance}
       tolerance.

       @arg points: Path points (C{LatLon}[]).
       @kwarg distance: Tolerance (C{meter}, same units as B{C{radius}}).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg indices: If C{True} return the simplified point indices
                       instead of the simplified points (C{bool}).
       @kwarg options: Optional keyword arguments passed thru to
                       function L{pygeodesy.equirectangular_}.

       @return: Simplified points (C{LatLon}[]).

       @raise LimitError: Lat- and/or longitudinal delta exceeds the B{C{limit}},
                          see function L{pygeodesy.equirectangular_}.

       @raise ValueError: Tolerance B{C{distance}} or B{C{radius}} too small.
    '''
    S = _Sy(points, distance, radius, True, indices, **options)

    r, n = S.r, S.n
    if n > 1:
        s2, _d2yxu4 = S.s2, S.d2yxu4

        i = 0
        for j in range(1, n):
            d2, _, _, _= _d2yxu4(i, j)
            if d2 > s2:
                r[j] = True
                i = j

    return S.points(r)


def simplifyRDP(points, distance=_1mm, radius=R_M, shortest=False,
                                       indices=False, **options):
    '''I{Ramer-Douglas-Peucker} (RDP) simplification of a path of
       C{LatLon} points.

       Eliminates any points too close together or closer to an
       edge than the given I{distance} tolerance.

       This C{RDP} method exhaustively searches for the point with
       the largest distance, resulting in worst-case complexity
       M{O(n**2)} where M{n} is the number of points.

       @arg points: Path points (C{LatLon}[]).
       @kwarg distance: Tolerance (C{meter}, same units as B{C{radius}}).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg shortest: If C{True} use the I{shortest} otherwise the
                        I{perpendicular} distance (C{bool}).
       @kwarg indices: If C{True} return the simplified point indices
                       instead of the simplified points (C{bool}).
       @kwarg options: Optional keyword arguments passed thru to
                       function L{pygeodesy.equirectangular_}.

       @return: Simplified points (C{LatLon}[]).

       @raise LimitError: Lat- and/or longitudinal delta exceeds the B{C{limit}},
                          see function L{pygeodesy.equirectangular_}.

       @raise ValueError: Tolerance B{C{distance}} or B{C{radius}} too small.
    '''
    S = _Sy(points, distance, radius, shortest, indices, **options)

    return S.rdp(False)


def simplifyRDPm(points, distance=_1mm, radius=R_M, shortest=False,
                                        indices=False, **options):
    '''Modified I{Ramer-Douglas-Peucker} (RDPm) simplification of a
       path of C{LatLon} points.

       Eliminates any points too close together or closer to an edge
       than the given I{distance} tolerance.

       This C{RDPm} method stops at the first point farther than the
       given distance tolerance, significantly reducing the run time
       (but producing results different from the original C{RDP} method).

       @arg points: Path points (C{LatLon}[]).
       @kwarg distance: Tolerance (C{meter}, same units as B{C{radius}}).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg shortest: If C{True} use the I{shortest} otherwise the
                        I{perpendicular} distance (C{bool}).
       @kwarg indices: If C{True} return the simplified point indices
                       instead of the simplified points (C{bool}).
       @kwarg options: Optional keyword arguments passed thru to
                       function L{pygeodesy.equirectangular_}.

       @return: Simplified points (C{LatLon}[]).

       @raise LimitError: Lat- and/or longitudinal delta exceeds the B{C{limit}},
                          see function L{pygeodesy.equirectangular_}.

       @raise ValueError: Tolerance B{C{distance}} or B{C{radius}} too small.
    '''
    S = _Sy(points, distance, radius, shortest, indices, **options)

    return S.rdp(True)


def simplifyRW(points, pipe=_1mm, radius=R_M, shortest=False,
                                  indices=False, **options):
    '''I{Reumann-Witkam} (RW) simplification of a path of C{LatLon}
       points.

       Eliminates any points too close together or within the given
       I{pipe} tolerance along an edge.

       @arg points: Path points (C{LatLon}[]).
       @kwarg pipe: Pipe radius, half-width (C{meter}, same units as
                    B{C{radius}}).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg shortest: If C{True} use the I{shortest} otherwise the
                        I{perpendicular} distance (C{bool}).
       @kwarg indices: If C{True} return the simplified point indices
                       instead of the simplified points (C{bool}).
       @kwarg options: Optional keyword arguments passed thru to
                       function L{pygeodesy.equirectangular_}.

       @return: Simplified points (C{LatLon}[]).

       @raise LimitError: Lat- and/or longitudinal delta exceeds the B{C{limit}},
                          see function L{pygeodesy.equirectangular_}.

       @raise ValueError: Tolerance B{C{pipe}} or B{C{radius}} too small.
    '''
    S = _Sy(points, pipe, radius, shortest, indices, **options)

    n, r = S.n, S.r
    if n > 1:
        s2, _d21, _d2i2 = S.s2, S.d21, S.d2i2

        s, e = 0, 1
        while s < e < n:
            if _d21(s, e):
                d2, i = _d2i2(e + 1, n, True)
                if d2 > s2 and i > 0:
                    r[s] = r[i] = True
                    s, e = i, i + 1
                else:
                    r[s] = True  # r[n-1] = True
                    break  # while loop
            else:  # drop points[e]
                e += 1

    return S.points(r)


def simplifyVW(points, area=_1mm, radius=R_M, attr=None,
                                  indices=False, **options):
    '''I{Visvalingam-Whyatt} (VW) simplification of a path of
       C{LatLon} points.

       Eliminates any points too close together or with a triangular
       area not exceeding the given I{area} tolerance I{squared}.

       This C{VW} method exhaustively searches for the single point
       with the smallest triangular area, resulting in worst-case
       complexity M{O(n**2)} where M{n} is the number of points.

       @arg points: Path points (C{LatLon}[]).
       @kwarg area: Tolerance (C{meter}, same units as B{C{radius}}).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg attr: Optional, points attribute to save the area value
                    (C{str}).
       @kwarg indices: If C{True} return the simplified point indices
                       instead of the simplified points (C{bool}).
       @kwarg options: Optional keyword arguments passed thru to
                       function L{pygeodesy.equirectangular_}.

       @return: Simplified points (C{LatLon}[]).

       @raise AttributeError: If an B{C{attr}} is specified for I{Numpy2}
                              B{C{points}}.

       @raise LimitError: Lat- and/or longitudinal delta exceeds the B{C{limit}},
                          see function L{pygeodesy.equirectangular_}.

       @raise ValueError: Tolerance B{C{area}} or B{C{radius}} too small.
    '''
    S = _Sy(points, area, radius, False, indices, **options)

    if S.vwn() > 2:
        # remove any points too close or
        # with a zero triangular area
        S.rm2(0)

        r, s2, s2e = S.r, S.s2, S.s2e
        # keep removing the point with the smallest
        # area until latter exceeds the tolerance
        while len(r) > 2:
            m, m2 = 0, s2e
            for i in range(1, len(r) - 1):
                h2 = r[i].h2
                if h2 < m2:
                    m, m2 = i, h2
            if m2 > s2:
                break
            S.rm1(m, 0)

    return S.vwr(attr)


def simplifyVWm(points, area=_1mm, radius=R_M, attr=None,
                                   indices=False, **options):
    '''I{Modified Visvalingam-Whyatt} (VWm) simplification of a path
       of C{LatLon} points.

       Eliminates any points too close together or with a triangular
       area not exceeding the given area tolerance I{squared}.

       This C{VWm} method removes all points with a triangular area
       below the tolerance in each iteration, significantly reducing
       the run time (but producing results different from the
       original C{VW} method).

       @arg points: Path points (C{LatLon}[]).
       @kwarg area: Tolerance (C{meter}, same units as B{C{radius}}).
       @kwarg radius: Mean earth radius (C{meter}).
       @kwarg attr: Optional, points attribute to save the area value
                    (C{str}).
       @kwarg indices: If C{True} return the simplified point indices
                       instead of the simplified points (C{bool}).
       @kwarg options: Optional keyword arguments passed thru to
                       function L{pygeodesy.equirectangular_}.

       @return: Simplified points (C{LatLon}[]).

       @raise AttributeError: If an B{C{attr}} is specified for I{Numpy2}
                              B{C{points}}.

       @raise LimitError: Lat- and/or longitudinal delta exceeds the B{C{limit}},
                          see function L{pygeodesy.equirectangular_}.

       @raise ValueError: Tolerance B{C{area}} or B{C{radius}} too small.
    '''
    S = _Sy(points, area, radius, False, indices, **options)

    if S.vwn() > 2:
        # remove all points with an area
        # not exceeding the tolerance
        S.rm2(S.s2)

    return S.vwr(attr)

# **) MIT License
#
# Copyright (C) 2016-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
