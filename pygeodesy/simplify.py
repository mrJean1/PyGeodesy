
# -*- coding: utf-8 -*-

u'''Simplify or linearize a path of C{LatLon} points.

Each of the 4 I{simplify} functions is based on a different algorithm and
produces different, simplified results in (very) different run times for the
same path:

 - Function L{simplify1} eliminates points with edge lengths shorter than
   the given tolerance.

 - Function L{simplifyRDP} implements the I{Ramer-Douglas-Peucker} (RDP)
   algorithm, iteratively finding the point farthest from each path edge.
   Original RDP exhaustively searches the most distant point in each iteration,
   I{modified} RDP stops at the first point exceeding the distance tolerance.

 - Function L{simplifyRW} uses the I{Reumann-Witkam} (RW) method, sliding a
   "pipe" over each path edge, removing all subsequent points within the pipe
   radius, up to the first point outside the pipe.

 - Function L{simplifyVW} provides the I{Visvalingam-Whyatt} (VW) method
   using the area of the triangle formed by three neigboring points.  Original
   VW removes only a single point per iteration, I{modified} VW eliminates all
   points with a triangular area not exceeding the tolerance in each iteration.

Keyword argument I{shortest} of functions L{simplifyRDP} and L{simplifyRW}
specifies of the distance between a point and a path edge.  If C{True}, use
the I{shortest} distance to the path edge or edge points, otherwise use the
I{perpendicular} distance to the (extended) edge through both points.

Keyword argument B{C{radius}} of all fuctions is set to the mean earth
radius in C{meter}, conventionally.  Other units may be used, provided
that radius and tolerance are specified in the same units.

Use keyword argument C{B{indices}=True} in any function to return a list
of I{indices} of simplified point instead of the simplified points with
the first and last index are always the first and last original index.

Finally, any additional keyword arguments B{C{options}} to all functions
are passed thru to function L{pygeodesy.equirectangular4} to specify the
distance approximation.

To process C{NumPy} arrays containing rows of lat-, longitude and possibly
other values, use class L{Numpy2LatLon} to wrap the C{NumPy} array into
I{on-the-fly-LatLon} points.  Pass the L{Numpy2LatLon} instance to any
I{simplify} function and the returned result will be a C{NumPy} array
containing the simplified subset, a partial copy of the original C{NumPy}
array.  Use keyword argument C{B{indices}=True} to return a list of array
row indices inlieu of the simplified array subset.

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
from pygeodesy.errors import _AttributeError, _ValueError, _xkwds_pop2
from pygeodesy.fmath import fdot_, len2, sqrt0
from pygeodesy.formy import equirectangular4
from pygeodesy.interns import _small_, _too_
from pygeodesy.iters import isNumpy2, isTuple2
# from pygeodesy.lazily import _ALL_LAZY  # from .units
from pygeodesy.units import _ALL_LAZY, _1mm, Radius_

from math import degrees, fabs, radians

__all__ = _ALL_LAZY.simplify
__version__ = '24.12.02'


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
    '''(INTERNAL) VW 2-tuple (index, area2).
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
    d2yxse5  = ()     # 5-tuple
    eps      = EPS    # system epsilon
    indices  = False
    ixs      = set()  # set(indices)
    n        = 0
    options  = {}
    pts      = []
    radius   = R_M    # mean earth radius
    s2       = EPS    # tolerance squared
    s2e      = EPS    # VW sentinel
    shortest = False  # i.e. perpendicular
    subset   = None   # isNumpy2 or isTuple2
    t2s      = []     # list(_T2s)

    def __init__(self, points, tolerance, radius, shortest,
                                          indices, **options):
        '''New C{Simplify} state.
        '''
        n, self.pts = len2(points)
        if n > 0:
            self.n   = n
            self.ixs = set((0, n-1))

        if radius is not R_M:
            self.radius = Radius_(radius, low=self.eps)
        # tolerance converted to degrees squared
        self.s2 = s2 = degrees(tolerance / self.radius)**2
        if min(s2, tolerance) < self.eps:
            raise _ValueError(tolerance=tolerance, txt=_too_(_small_))
        self.s2e = s2 + _1_0  # VW sentinel
        # assert self.s2e > s2

        if indices:
            self.indices = True
        if options:
            _, self.options = _xkwds_pop2(options, modified=None)
        if shortest:
            self.shortest = True
        if isNumpy2(points) or isTuple2(points):  # NOT self.pts
            self.subset = points.subset

    def d21(self, s, e):
        '''Set path edge or line thru (points[s], -[e]).
        '''
        d21, y21, x21, _ = self.d2yxu4(s, e)
        self.d2yxse5 = d21, y21, x21, s, e
        return d21 > self.eps

    def d2i2(self, m, n, modified):
        '''Find the tallest distance among all points[m..n]
           to (points[s], -[e]) exceeding the tolerance.
        '''
        _, _, _, s, _ = self.d2yxse5
        t2, t = self.s2, 0  # tallest
        for i in range(m, n):
            d2, _, _, _ = self.d2yxu4(s, i)
            if d2 > t2:
                t2, t = d2, i
                if modified and d2 > self.eps:
                    break
        return t2, t

    def d2ix2(self, m, n, modified):
        '''Find the tallest I{perpendicular B{or} shortest} distance
           among all points[m..n] to the path edge or line through
           (points[s], -[e]) exceeding the tolerance.
        '''
        h = not self.shortest
        # point (x, y) on axis rotated by angle a ccw:
        #   x' = y * sin(a) + x * cos(a)
        #   y' = y * cos(a) - x * sin(a)
        #
        # distance along (w) and perpendicular to (h)
        # a line from the origin to point (dx, dy):
        #   w = (y * dy + x * dx) / hypot(dx, dy)
        #   h = (y * dx - x * dy) / hypot(dx, dy)
        d21, y21, x21, s, e = self.d2yxse5
        t2, t = self.s2, 0  # tallest
        for i in range(m, n):
            # distance points[s] to -[i], ...
            d2, y01, x01, _ = self.d2yxu4(s, i)
            if d2 > self.eps:
                if h:  # perpendicular distance
                    d2 = fdot_(y01, x21, -x01, y21)**2 / d21
                else:
                    w  = fdot_(y01, y21,  x01, x21)
                    if w > 0:
                        if w < d21:  # ... perpendicular ...
                            d2 = fdot_(y01, x21, -x01, y21)**2 / d21
                        else:  # ... or points[e] to -[i]
                            d2, _, _, _ = self.d2yxu4(e, i)
                if d2 > t2:
                    t2, t = d2, i
                    if modified:
                        break
        return t2, t

    def d2yxu4(self, i, j):
        '''Return the distance I{squared}, the deltas and the
           (longitudinal) unrollment between (points[i], -[j]).
        '''
        p1, p2 = self.pts[i], self.pts[j]
        return equirectangular4(p1.lat, p1.lon,
                                p2.lat, p2.lon, **self.options)

    def h2t(self, i1, i2, i3):
        '''Compute (double) the triangle area, points[i2] is
           the top and edge (points[i1], -[i3]) is the base
           of the triangle.
        '''
        d21, y21, x21 , _ = self.d2yxu4(i1, i3)
        if d21 > self.eps:
            d01, y01, x01, _ = self.d2yxu4(i1, i2)
            if d01 > self.eps:
                h2 = fabs(fdot_(y01, x21, -x01, y21))
                # triangle height h = h2 / sqrt(d21) and
                # the area = h * sqrt(d21) / 2 == h2 / 2
                return h2  # double triangle area
        return 0

    def rdp(self, modified):
        '''Ramer-Douglas-Peucker (RDP) simplification of a
           path of C{LatLon} points.

           @arg modified: Use I{modified} RDP (C{bool}).
        '''
        r, n = self.ixs, self.n
        if n > 1:
            s2, se = self.s2, [(0, n-1)]
            while se:
                s, e = se.pop()
                s1 = s + 1
                if e > s1:
                    if self.d21(s, e):  # points[] to edge [s, e]
                        d2, i = self.d2ix2(s1, e, modified)
                    else:  # points[] to point [s]
                        d2, i = self.d2i2( s1, e, modified)
                    if d2 > s2 and i > 0:
                        se.append((i, e))
                        if not modified:
                            se.append((s, i))
                        r.add(i)
                    r.add(s)
        return self.result(r)

    def result(self, r):
        '''Return the simplified points or indices.
        '''
        r = sorted(r)
        if self.indices:
            return list(r)
        elif self.subset:
            return self.subset(r)
        else:
            return list(self.pts[i] for i in r)

    def rw(self):
        '''Reumann-Witkam simplification.
        '''
        r, n = self.ixs, self.n
        if n > 1:
            s, e, s2 = 0, 1, self.s2
            while s < e < n:
                if self.d21(s, e):
                    d2, i = self.d2ix2(e + 1, n, True)
                    r.add(s)
                    if d2 > s2 and i > 0:
                        r.add(i)
                        s = e = i
                    else:
                        # r.add(n - 1)
                        break
                e += 1
        return self.result(r)

    def sy1(self):
        '''Basic simplification.
        '''
        r, n = self.ixs, self.n
        if n > 1:
            s2, i = self.s2, 0
            for j in range(1, n):
                d2, _, _, _ = self.d2yxu4(i, j)
                if d2 > s2:
                    r.add(j)
                    i = j
        return self.result(r)

    def vwn(self):
        '''Initialize VW as list of 2-tuples _T2(ix, h2) where
           ix is the points[] index and h2 is the triangular
           area I{(times 2)} of that point.
        '''
        self.t2s = t = []
        n, T2 = self.n, _T2
        if n > 2:
            _h2t = self.h2t
            t[:] = [T2(i, _h2t(i-1, i, i+1)) for i in range(1, n - 1)]
        if n > 1:
            t.append(T2(n - 1, self.s2e))
        if n > 0:
            t.insert(0, T2(0, self.s2e))
        return len(t)

    def vwr(self, attr):
        '''Return the VW results, optionally including the
           triangular area (in C{meter}) as attribute C{attr}
           to each simplified point.
        '''
        pts, t = self.pts, self.t2s

        # double check the minimal triangular area
        assert min(t2.h2 for t2 in t) > self.s2 > 0

        if attr:  # return each trangular area (actually
            # the sqrt of double the triangular area)
            # converted back from degrees to meter
            if isNumpy2(pts):
                raise _AttributeError(attr=attr)
            t[0].h2 = t[-1].h2 = 0  # zap sentinels
            m = radians(_1_0) * self.radius  # meter
            for t2 in t:  # convert back to meter
                setattr(pts[t2.ix], attr, sqrt0(t2.h2) * m)

        n = len(t)  # double check for duplicates
        r = set(t2.ix for t2 in t)
        assert len(r) == n
        return self.result(r)

    def vwrm(self):
        '''Keep removing the VW point with the smallest triangular
           area until that area exceeds the tolerance.
        '''
        s2, t = self.s2, self.t2s
        while len(t) > 2:
            m2, m = t[1].h2, 1
            for i in range(2, len(t) - 1):
                h2 = t[i].h2
                if h2 < m2:
                    m2, m = h2, i
            if m2 > s2:
                break
            self.vwrm1(m, 0)

    def vwrm1(self, m, tol):
        '''Eliminate VW point[m], keep recomputing the trangular
           area of both neighboring points and removing those
           too until the recomputed area exceeds C{tol}.
        '''
        t, _h2t = self.t2s, self.h2t
        t.pop(m)
        for n in (m, m - 1):  # neighbors
            while 0 < n < (len(t) - 1):
                h2 = _h2t(t[n-1].ix, t[n].ix, t[n+1].ix)
                if h2 > tol:
                    t[n].h2 = h2
                    break  # while
                t.pop(n)

    def vwrm2(self, tol):
        '''Eliminate all VW points with a triangular area not
           exceeding C{tol}.
        '''
        t = self.t2s
        m = len(t) - 1
        while m > 1:
            m -= 1
            if t[m].h2 <= tol:
                self.vwrm1(m, tol)
                m = min(m, len(t) - 1)


def simplify1(points, distance=_1mm, radius=R_M, indices=False, **options):
    '''Basic simplification of a path of C{LatLon} points by eliminating
       any points closer together than the given I{distance} tolerance.

       @arg points: Iterable with the path points (C{LatLon}[]).
       @kwarg distance: Tolerance (C{meter}, same units as B{C{radius}}).
       @kwarg radius: Mean earth radius (C{meter}, conventionally).
       @kwarg indices: If C{True}, return B{C{points}} indices instead
                       of the simplified points (C{bool}).
       @kwarg options: Optional keyword arguments passed thru to function
                       L{pygeodesy.equirectangular4}.

       @return: Simplified points (C{LatLon}[]) or B{C{points}} indices.

       @raise LimitError: Lat- and/or longitudinal delta exceeds the B{C{limit}},
                          see function L{pygeodesy.equirectangular4}.

       @raise ValueError: Tolerance B{C{distance}} or B{C{radius}} too small.
    '''
    S = _Sy(points, distance, radius, True, indices, **options)
    return S.sy1()


def simplifyRDP(points, distance=_1mm, radius=R_M, shortest=False,
                        indices=False, modified=False, **options):
    '''I{Ramer-Douglas-Peucker} (RDP) simplification of a path of C{LatLon}
       points by eliminating any points too close together or closer to an
       edge than the given I{distance} tolerance.

       @arg points: Iterable with the path points (C{LatLon}[]).
       @kwarg distance: Tolerance (C{meter}, same units as B{C{radius}}).
       @kwarg radius: Mean earth radius (C{meter}, conventionally).
       @kwarg shortest: If C{True}, use the I{shortest} otherwise the
                        I{perpendicular} distance (C{bool}).
       @kwarg indices: If C{True}, return B{C{points}} indices instead
                       of the simplified points (C{bool}).
       @kwarg modified: If C{True}, use the C{modified RDP} method (C{bool}),
                        see the B{note}.
       @kwarg options: Optional keyword arguments passed thru to function
                       L{pygeodesy.equirectangular4}.

       @return: Simplified points (C{LatLon}[]) or B{C{points}} indices.

       @raise LimitError: Lat- and/or longitudinal delta exceeds the B{C{limit}},
                          see function L{pygeodesy.equirectangular4}.

       @raise ValueError: Tolerance B{C{distance}} or B{C{radius}} too small.

       @note: The original C{RDP} method exhaustively searches for the point
              with the largest distance (resulting in complexity M{O(n**2)}
              with M{n} is the number of points).  The B{C{modified}} C{RDP}
              method stops at the first point farther than the B{C{distance}}
              tolerance, significantly reducing the run time (but producing
              results different from the original C{RDP} method).
    '''
    S = _Sy(points, distance, radius, shortest, indices, **options)
    return S.rdp(bool(modified))


def simplifyRW(points, pipe=_1mm, radius=R_M, shortest=False,
                                  indices=False, **options):
    '''I{Reumann-Witkam} (RW) simplification of a path of C{LatLon} points
       by eliminating any points too close together or within the given
       I{pipe} tolerance along an edge.

       @arg points: Iterable with the path points (C{LatLon}[]).
       @kwarg pipe: Pipe radius, half-width (C{meter}, same units as
                    B{C{radius}}).
       @kwarg radius: Mean earth radius (C{meter}, conventionally).
       @kwarg shortest: If C{True}, use the I{shortest} otherwise the
                        I{perpendicular} distance (C{bool}).
       @kwarg indices: If C{True}, return B{C{points}} indices instead
                       of the simplified points (C{bool}).
       @kwarg options: Optional keyword arguments passed thru to function
                       L{pygeodesy.equirectangular4}.

       @return: Simplified points (C{LatLon}[]) or B{C{points}} indices.

       @raise LimitError: Lat- and/or longitudinal delta exceeds the B{C{limit}},
                          see function L{pygeodesy.equirectangular4}.

       @raise ValueError: Tolerance B{C{pipe}} or B{C{radius}} too small.
    '''
    S = _Sy(points, pipe, radius, shortest, indices, **options)
    return S.rw()


def simplifyVW(points, area=_1mm, radius=R_M, indices=False,
                       attr=None, modified=False, **options):
    '''I{Visvalingam-Whyatt} (VW) simplification of a path of C{LatLon}
       points by eliminating any points too close or with a triangular
       area not exceeding the given I{area} tolerance I{squared}.

       @arg points: Iterable with the path points (C{LatLon}[]).
       @kwarg area: Tolerance (C{meter}, same units as B{C{radius}}).
       @kwarg radius: Mean earth radius (C{meter}, conventionally).
       @kwarg indices: If C{True}, return B{C{points}} indices instead
                       of the simplified points (C{bool}).
       @kwarg attr: Optional, B{C{points}} attribute to save the area
                    value (C{str}).
       @kwarg modified: If C{True}, use the C{modified VW} method (C{bool}),
                        see the B{note}.
       @kwarg options: Optional keyword arguments passed thru to function
                       L{pygeodesy.equirectangular4}.

       @return: Simplified points (C{LatLon}[]) or B{C{points}} indices.

       @raise AttributeError: An B{C{attr}} isinvalid for I{Numpy2} B{C{points}}.

       @raise LimitError: Lat- and/or longitudinal delta exceeds the B{C{limit}},
                          see function L{pygeodesy.equirectangular4}.

       @raise ValueError: Tolerance B{C{area}} or B{C{radius}} too small.

       @note: The original C{VW} method exhaustively searches for the point
              with the smallest triangular I{area} (resulting in complexity
              M{O(n**2)} with M{n} the number of points).  The B{C{modified}}
              C{VW} method removes I{all} points with a triangular I{area}
              below the tolerance in each iteration, significantly reducing
              the run time (but producing results different from the original
              C{VW} method).
    '''
    S = _Sy(points, area, radius, False, indices, **options)
    if S.vwn() > 2:
        if modified:
            S.vwrm2(S.s2)
        else:
            S.vwrm2(0)
            S.vwrm()
    return S.vwr(attr)

# **) MIT License
#
# Copyright (C) 2016-2025 -- mrJean1 at Gmail -- All Rights Reserved.
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
