
# -*- coding: utf-8 -*-

u'''Clip a path or polygon of C{LatLon} points against a rectangular box or a
(convex or arbitrary) clip region.

Box clip functions L{clipCS4} I{Cohen-Sutherland} and L{clipLB6} I{Liang-Barsky},
region clip functions L{clipGH4} I{Greiner-Hormann} and L{clipSH} and L{clipSH3}
I{Sutherland-Hodgeman}.

Class L{BooleanGH} for I{boolean} polygon operations based on L{clipGH4}
I{Greiner-Hormann}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

# from pygeodesy.basics import len2, _xcopy  # from .iters, .named
from pygeodesy.constants import EPS, MAX, _0_0, _1_0
from pygeodesy.errors import _AssertionError, PointsError, _ValueError
from pygeodesy.fmath import fabs, favg, hypot, hypot2
from pygeodesy.fsums import fsum_, Property_RO, property_RO
from pygeodesy.interns import NN, _COMMASPACE_, _convex_, _DOT_, _ELLIPSIS_, \
                             _end_, _few_, _fi_, _height_, _i_, _j_, _lat_, \
                             _lon_, _near_, _no_, _not_, _SPACE_, _STAR_, \
                             _start_, _too_
from pygeodesy.iters import _imdex2, len2, points2
from pygeodesy.lazily import _ALL_LAZY, _ALL_DOCS, _FOR_DOCS
from pygeodesy.named import _Named, _NamedTuple, _Pass, _xcopy
from pygeodesy.points import areaOf, boundsOf, isconvex_, LatLon_
# from pygeodesy.props import Property_RO, property_RO  # from .fsums
from pygeodesy.units import Bool, FIx, Height, Lat, Lon, Number_

# from math import fabs  # from .fmath

__all__ = _ALL_LAZY.clipy
__version__ = '23.02.08'

_0_EPS =  EPS  # near-zero, positive
_EPS_0 = -EPS  # near-zero, negative
_1_EPS = _1_0 + EPS  # over-one
_EPS_1 = _1_0 - EPS  # near-one

_box_    = 'box'
_clip_   = 'clip'
_clipid_ = 'clipid'
_fj_     = 'fj'
_open_   = 'open'
_region_ = 'region'


class ClipError(_ValueError):
    '''Clip box or clip region issue.
    '''
    def __init__(self, *name_n_corners, **txt_cause):
        '''New L{ClipError}.

           @arg name_n_corners: Either just a name (C{str}) or
                                name, number, corners (C{str},
                                C{int}, C{tuple}).
           @kwarg txt_cause: Optional C{B{txt}=str} explanation
                             of the error and C{B{cause}=None}
                             for exception chaining.
        '''
        if len(name_n_corners) == 3:
            t, n, v = name_n_corners
            n = _SPACE_(t, _clip_, _box_ if n == 2 else _region_)
            name_n_corners = n, v
        _ValueError.__init__(self, *name_n_corners, **txt_cause)


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


def _min_max_eps(*xs):
    '''(INTERNAL) Return 2-tuple C{(min, max)}, oversized.
    '''
    lo, hi = min(xs), max(xs)
    lo *= _1_EPS if lo < 0 else _EPS_1
    hi *= _EPS_1 if hi < 0 else _1_EPS
    return (lo or _EPS_0), (hi or _0_EPS)


def _neq(p1, p2):
    '''(INTERNAL) Check for not near-equal points.
    '''
    return fabs(p1.lat - p2.lat) > EPS or \
           fabs(p1.lon - p2.lon) > EPS


def _outside(x1, x2, lo, hi):
    '''(INTERNAL) Is C{(x1, x2)} outside C{(lo, hi)}?
    '''
    return max(x1, x2) < lo or min(x1, x2) > hi


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


class _GHcircular(_Named):
    '''(INTERNAL) A I{circular, doubly-linked} list of L{LatLonGH}s,
       representing the original points of a polygon plus -if any-
       the intersections with another.
    '''

    _clipstr = 'corners '
    _first   =  None  # first latlon
    _last    =  None  # last-appended
    _len     =  0
    _raiser  =  False
    _xtend   =  False

    def __init__(self, lls, name=NN, raiser=False, xtend=False):
        '''New L{BooleanGH}.

           @arg lls: The polygon points (iterable of 2 or more
                     C{LatLon}s, L{LatLonGH}s or L{ClipGH4Tuple}s).
           @kwarg name: Optional name (C{str}).
           @kwarg raiser: If C{True} throw C{ClipError} exceptions
                          for errors or I{degenerate cases}.
           @kwarg xtend: If C{True} handle I{degenerate cases}.
        '''
        if name:
            self.name = name
        _a, LL = self._append, LatLonGH
        for ll in lls:
            _a(LL(ll))

        if xtend:
            self._xtend = True
        elif raiser:
            self._raiser = True

#   def __contains__(self, lalon):
#       # Is C{latlon} in this polygon?
#       for v in self:
#           if v is vertex:
#               return True
#       return False

    def __iter__(self):
        # Yield all points and intersections.
        v = f = self._first
        while True:
            yield v
            v = v._next
            if v is f:
                break

    def __len__(self):
        # Return the number of latlons.
        return self._len

    def __repr__(self):
        # String C{repr} of this C{GHcdl}.
        return '%s[%d](%s)' % (self.__class__.__name__, len(self), self)

    def __str__(self):
        '''String C{str} of this polygon.
        '''
        return ', '.join(str(v) for v in self)

    def _append(self, x_ll, *y_clipid):
        # Append a point given as C{x}, C{y} and C{clipid} args
        # or as a L{LatLonGH} instance.
        self._last = v = LatLonGH(x_ll, *y_clipid) if y_clipid else x_ll
        self._len += 1

        v._next = n = self._first
        if n:
            v._prev = p = n._prev
        else:
            self._first = p = n = v
        p._next = n._prev = v
        return v

#   def _appendedup(self, x_v, *y_clipid):
#       # Append a point if not a duplicate of the previous one.
#       p = self._last
#       return None if p and p == LatLonGH(x_v, *y_clipid) else \
#                             self._append(x_v, *y_clipid)

    @property_RO
    def _bottom_top_eps(self):
        # Get the bottom and top bounds of C{y}, oversized.
        return _min_max_eps(min(v.y for v in self),
                            max(v.y for v in self))

    def _clip(self, corners, s_entry, c_entry, **raiser_xtend):
        # Clip this polygon with another one, C{corners}.

        # Core of Greiner/Hormann's algorithm, enhanced U{Correia's
        # <https://GitHub.com/helderco/univ-polyclip>} implementation***
        # and extended to optionally handle so-called "degenerate cases"
        S, C = self, _GHcircular(corners or ())
        if C:
            lr = C._left_right_eps
            bt = C._bottom_top_eps

            # 1. find intersections
            for s1, s2 in S._edges2(corners=False, **raiser_xtend):
                if not (_outside(s1.x, s2.x, *lr) or
                        _outside(s1.y, s2.y, *bt)):
                    se = _GHedge(s1, s2, **raiser_xtend)
                    _, sx, _, sy = se._x_sx_y_sy
                    if sx or sy:
                        _i3 = se._intersections3
                        for c1, c2 in C._edges2(**raiser_xtend):
                            for xyc, sa, ca in _i3(c1, c2):
                                S._insert(xyc, s1, s2, sa,
                                C._insert(xyc, c1, c2, ca))

            # 2. identify entry/exit intersections
            if S._first:
                s_entry ^= S._first._inside(C, *bt)
                for v in S._intersections():
                    v._entry = s_entry = not s_entry

            if C._first:
                c_entry ^= C._first._inside(S)
                for v in C._intersections():
                    v._entry = c_entry = not c_entry

            # 3. yield the result(s)
            return self._results()

        raise ClipError(_no_(self._clipstr).strip(), txt=repr(corners))

    def clipids(self):
        '''Return a tuple with all C{clipid}s.
        '''
        return tuple(set(v.clipid for v in self))

    def _edges2(self, corners=True, raiser=False, **unused):
        # Yield each edge as a pair of points, 2-tuple C{(LatLonGH, LatLonGH)}.
        def _OpenClipError(s, e):
            t = self._clipstr if corners else NN
            t = NN(_open_, _SPACE_, t, _clip_)
            return ClipError(t, txt=NN(s, _ELLIPSIS_(_COMMASPACE_, e)))

        s = p1 = f = self._first
        assert not p1._linked
        while p1:
            p2 = p1._next
            while p2._linked:
                p2 = p2._next
            if p1.clipid != p2.clipid:  # next clip
                if p1 != s:  # close old one iff open
                    if raiser and s is not f:
                        raise _OpenClipError(s, p1)
                    yield p1, s
                s = p2  # start new
            else:
                yield p1, p2
            if p2 is f:
                if raiser and s is not f and p2 != s:
                    raise _OpenClipError(s, p2)
                break
            p1 = p2

    def _insert(self, xyc, start, end, alpha, link=None):
        # Create and insert an intersection between points
        # C{start} and C{end}, ordered by C{alpha}.
        v = LatLonGH(*xyc)
        v._alpha = alpha
#       v._entry = False
        if link:
            v._linked = k = link
            k._linked = v
        h = favg(start.height, end.height, f=alpha)
        v.height = HeightGH(h)  # i.e. intersection
        self._len += 1

#       assert start is not end
        n = start  # .next
        while n is not end and n._alpha < alpha:
            n = n._next
        v._next = n
        v._prev = p = n._prev
        p._next = n._prev = v
        return v

    def _intersections(self):
        # Yield the intersections.
        for v in self:
            if v._linked:
                yield v

    @property_RO
    def _left_right_eps(self):
        # Get the left and right bounds of C{x}, oversized.
        return _min_max_eps(min(v.x for v in self),
                            max(v.x for v in self))

    def _points(self):
        # Yield the points in original order.
        for v in self:
            if not v._linked:
                yield v

    def _results(self):
        # Yield each clipped point as a L{ClipGH4Tuple}.
        for c, v in enumerate(self._unchecked()):
            p = s = v._to4Tuple(c)  # new clip
            while not v._checked:
                v._check()
                t = v._to4Tuple(c)
                if t != p or p is s:  # dedup
                    yield t
                    p = t
                e = v._entry
                while True:
                    v = v._next if e else v._prev
                    t = v._to4Tuple(c)
                    if t != p:  # dedup
                        yield t
                        p = t
                    if v._linked:
                        break
                v = v._linked  # switch
            if p != s:  # close clip
                yield s

    def _unchecked(self):
        # Yield the (first) unchecked, intersection(s).
        while True:
            for v in self._intersections():
                if not v._checked:
                    yield v
                    break  # for
            else:
                break  # while


class _GHedge(object):
    # An edge between two L{LatLonGH} points.

    _raiser = False
    _xtend  = False

    def __init__(self, s1, s2, raiser=False, xtend=False):
        # New edge between points C{s1} and C{s2}, each a L{LatLonGH}.
        self._s1, self._s2 = s1, s2
        self._x_sx_y_sy = (s1.x, s2.x - s1.x,
                           s1.y, s2.y - s1.y)

        self._left_right_eps = _min_max_eps(s1.x, s2.x)
        self._bottom_top_eps = _min_max_eps(s1.y, s2.y)

        if xtend:
            self._xtend = True
        elif raiser:
            self._raiser = True

    def __str__(self):
        return 'edge(%s, %s)' % (self._s1, self._s2)

    def _alpha1(self, alpha):
        a = max(_0_0, min(alpha, _1_0))
        if fabs(a - alpha) < _0_EPS:
            return a
        raise ClipError('alpha outside [0..1]: %.3f' % (alpha,))

    def _alpha2(self, x, y, dx, dy):
        # Return C{(alpha)}, see .points.nearestOn5
        d = (y * dx - x * dy) / self._hypot0
        a = (y * dy + x * dx) / self._hypot2
        return a, fabs(d)

    def _error(self, n, c1, c2):
        t = _SPACE_(self, 'intersects', _GHedge(c1, c2))
        raise ClipError(t, txt='unhandled case %s' % (n,))

    @Property_RO
    def _hypot0(self):
        _, sx, _, sy = self._x_sx_y_sy
        return hypot( sx, sy) * _0_EPS

    @Property_RO
    def _hypot2(self):
        _, sx, _, sy = self._x_sx_y_sy
        return hypot2(sx, sy)

    def _intersections3(self, c1, c2, parallel=True):  # MCCABE 14
        # Yield the intersections this and another edge.

        # @return: None, 1 or 2 intersections, each a 3-Tuple
        #          (xyc, s_alpha, c_alpha) with intersection
        #          point C{(x, y, clipid)} and both alphas.

        # @raise ClipError: Intersection unhandled.

        # @see: U{Intersection point of two line segments
        #       <http://PaulBourke.net/geometry/pointlineplane/>}.
        c1_x, c1_y = c1.x, c1.y
        if not (_outside(c1_x, c2.x, *self._left_right_eps) and
                _outside(c1_y, c2.y, *self._bottom_top_eps)):
            x, sx, \
            y, sy = self._x_sx_y_sy

            cx = c2.x - c1_x
            cy = c2.y - c1_y
            d  = cy * sx - cx * sy

            if fabs(d) > _0_EPS:
                dx = x - c1_x
                dy = y - c1_y
                ca = (sx * dy - sy * dx) / d
                if _0_EPS < ca < _EPS_1 or (self._xtend and
                   _EPS_0 < ca < _1_EPS):
                    sa = (cx * dy - cy * dx) / d
                    if _0_EPS < sa < _EPS_1 or (self._xtend and
                       _EPS_0 < sa < _1_EPS):
                        x += sa * sx
                        y += sa * sy
                        yield (x, y, 0), sa, ca

                    # unhandled, "degenerate" cases 1, 2 or 3
                    elif self._raiser and not (sa < _EPS_0 or sa > _1_EPS):
                        self._error(1, c1, c2)  # intersection at s1 or s2

                elif self._raiser and not (ca < _EPS_0 or ca > _1_EPS):
                    # intersection at c1 or c2 or at c1 or c2 and s1 or s2
                    sa = (cx * dy - cy * dx) / d
                    e = 2 if sa < _EPS_0 or sa > _1_EPS else 3
                    self._error(e, c1, c2)

            elif parallel and (sx or sy) and (cx or cy):
                # non-null, parallel or colinear edges
                sa1, d1 = self._alpha2(c1_x - x, c1_y - y, sx, sy)
                sa2, d2 = self._alpha2(c2.x - x, c2.y - y, sx, sy)
                if max(d1, d2) < _0_EPS:
                    if self._xtend and not _outside(sa1, sa2, _EPS_0, _1_EPS):
                        if sa1 > sa2:  # anti-parallel
                            sa1, sa2 =  sa2,  sa1
                            ca1, ca2 = _1_0, _0_0
                        else:  # parallel
                            ca1, ca2 = _0_0, _1_0
                        ca = fabs((sx / cx) if cx else (sy / cy))
                        #  = hypot(sx, sy) / hypot(cx, cy)
                        if sa1 < 0:  # s1 is between c1 and c2
                            ca *= ca1 + sa1
                            yield (x, y, 0), ca1, self._alpha1(ca)
                        else:  # c1 is between s1 and s2
                            t = (x + sa1 * sx), (y + sa1 * sy), 0
                            yield t, sa1, ca1
                        if sa2 > 1:  # s2 is between c1 and c2
                            t = (x + sx), (y + sy), 0  # s2
                            ca *= sa2 - _1_0
                            yield t, ca2, self._alpha1(ca2 - ca)
                        else:  # c2 is between s1 and s2
                            t = (x + sa2 * sx), (y + sa2 * sy), 0
                            yield t, sa2, ca2
                    elif self._raiser and not _outside(sa1, sa2, _0_0, _1_EPS):
                        self._error(4, c1, c2)

#   def _intersect3(self, c1, c2):
#       # Intersect this and another edge.
#
#       # @return: 3-Tuple (xyc, s_alpha, c_alpha) with the
#       #          intersection point C{xyc} and both alphas.
#
#       # @raise ClipError: Intersection unhandled.
#
#       # @see: U{Intersection point of two line segments
#       #       <http://PaulBourke.net/geometry/pointlineplane/>}.
#
#       if _outside(c1.x, c2.x, *self._left_right_eps) or \
#          _outside(c1.y, c2.y, *self._bottom_top_eps):
#           sa = ca = None
#       else:
#           x, sx, \
#           y, sy = self._x_sx_y_sy
#
#           cx = c2.x - c1.x
#           cy = c2.y - c1.y
#           d  = cy * sx - cx * sy
#
#           if fabs(d) > _0_EPS:
#               dx = x - c1.x
#               dy = y - c1.y
#               sa = (cx * dy - cy * dx) / d
#               ca = (sx * dy - sy * dx) / d
#               if _0_EPS < ca < _EPS_1:
#                   if _0_EPS < sa < _EPS_1:
#                       x += sa * sx
#                       y += sa * sy
#                       return (x, y, 0), sa, ca
#
#                   # unhandled, "degenerate" cases 1, 2 or 3
#                   if self.raiser and not (sa < _EPS_0 or sa > _1_EPS):
#                       self._error(1, c1, c2)  # insection at s1 or s2
#
#               elif self.raiser and not (ca < _EPS_0 or ca > _1_EPS):
#                   # intersection at c1 or c2 or at c1 or c2 and s1 or s2
#                   self._error(2 if sa < _EPS_0 or sa > _1_EPS else 3, c1, c2)
#
#           else:  # null, parallel or colinear edges
#               if self.raiser and (sx or sy) and (cx or cy):
#                   h = hypot(sx, sy) * _0_EPS
#                   if min(_perpendicular(c1.x - x, c1.y - y, sx, sy),
#                          _perpendicular(c2.x - x, c2.y - y, sx, sy)) < h:
#                       self._error(4, c1, c2)  # colinear, overlapping
#               sa = ca = None
#
#       return None, sa, ca  # no intersection


class HeightGH(Height):
    '''Like L{Height} but to distnguish the interpolated
       height at an intersection from an original L{Height}.
    '''
    pass


class LatLonGH(object):
    '''A point or intersection in L{BooleanGH} polygon.
    '''

    _alpha   = 0      # relative length iff intersection
    _checked = False  # checked in phase 3 iff intersection
    _linked  = None   # link to neighbor iff intersection
    _entry   = None   # entry or exit iff intersection
    _next    = None   # link to the next vertex
    _prev    = None   # link to the previous vertex

    clipid   = 0      # polygonal clip identifier, number
    height   = 0      # interpolated height, usually meter

    def __init__(self, x_ll, *y_clipid):
        '''New C{LatLonGH} from separate C{x}, C{y} and C{clipid},
           or from a previous L{LatLonGH} or some other C{LatLon}
           instance.

           @arg x_ll: X coordinate (C{scalar}) or a lat/longitude
                      (L{LatLonGH}, C{LatLon} or L{ClipGH4Tuple}).
           @arg y_clipid: Scalar Y coordinate and C{clipid} iff
                          B{C{x_ll}} is scalar, ignored otherwise.
        '''
        if y_clipid:
            self.x = x_ll
            self.y, c = y_clipid
        else:
            self.x, self.y = x_ll.lon, x_ll.lat
            c = getattr(x_ll, _clipid_, 0)
            h = getattr(x_ll, _height_, 0)
            if self.height != h:
                self.height = h

        if self.clipid != c:
            self.clipid = c

    def __repr__(self):
        '''String C{repr} of this lat-/longitude.
        '''
        s = str(self)
        if self._linked:
            s += 'e' if self._entry else 'x'
            if not self._checked:
                s += '!'
        return _SPACE_(self._prev, '<->', s, '<->', self._next)

    def __str__(self):
        '''String C{str} of this lat-/longitude.
        '''
        t = tuple('%s=%.8g' % t for t in ((_lat_, self.lat),
                                          (_lon_, self.lon)))
        if self.height:
            x  = _STAR_ if self.isintersection else NN
            t += ('%s=%.6g%s' % (_height_, self.height, x)),
        if self.clipid:
            t += ('%s=%s' % (_clipid_, self.clipid)),
        return '(%s)' % (_COMMASPACE_.join(t),)

    def _check(self):
        # Check-mark this vertex and it's linked one.
        self._checked = True
        b = self._linked
        if b and not b._checked:
            b._check()

    def copy(self, deep=False, **unused):  # see .named._Named.copy
        '''Make a shallow or deep copy of this instance.

           @kwarg deep: If C{True} make a deep, otherwise
                        a shallow copy (C{bool}).

           @return: The copy (C{This class} or sub-class thereof).
        '''
        return _xcopy(self, deep=deep)

    def _inside(self, polygon, *bottom_top):
        # Is this vertex inside a polygon?  I{Odd-even rule}.

        # The I{odd-even} rule counts the number of edges
        # intersecting a ray emitted East-bound from this
        # point to infinity.  When I{odd} this point lies
        # inside, if I{even} outside.
        r, y = False, self.y
        if bottom_top and not _outside(y, y, *bottom_top):
            _i3 = _GHedge(self, LatLonGH(MAX, y, 0))._intersections3
            for p1, p2 in polygon._edges2():
                for xyc, _, _ in _i3(p1, p2, False):
                    if xyc:
                        r = not r
        return r

    def isinside(self, points):
        '''Is this point inside a polygon based on C{odd-even rule}?

           @arg points: Iterable of the polygon points (L{ClipGH4Tuple},
                        L{LatLonGH}, L{LatLon_}, etc.
        '''
        self._inside(BooleanGH(points, name=self.isinside.__name__))

    @property_RO
    def isintersection(self):
        '''Is this an intersection?
        '''
        return bool(self._linked)

    @property_RO
    def ispoint(self):
        '''Is this an original (boolean) point?
        '''
        return not bool(self._linked)

    @property_RO
    def lat(self):
        return self.y

    @property_RO
    def lon(self):
        return self.x

    def _to4Tuple(self, clipid=0):
        # Return this vertex as a L{ClipGH4Tuple}.
        h = self.height
        if self.isintersection:
            h = HeightGH(h)
        elif h:
            h = Height(h)
        else:
            h = 0
        return ClipGH4Tuple(self.y, self.x, h, clipid)


class BooleanGH(_GHcircular):
    '''Polygon class providing I{boolean} operations between
       two polygons based on L{clipGH4}.  The supported
       operations between polygon A and B are:

        -  C = A | B  or  A |= B,  union of A and B

        -  C = A & B  or  A &= B,  intersection of A and B

        -  C = A - B  or  A -= B,  difference A less B

        -  C = B - A  or  B -= A,  difference B less A
    '''
    _clipstr = 'boolean '
    _raiser  =  False
    _xtend   =  True

    if _FOR_DOCS:
        __init__ = _GHcircular.__init__

    def _boolean(self, other, s_entry, c_entry, op):
        s = BooleanGH(self, name=self.name)
        c = self._other(other)
        r = s._clip(c, s_entry, c_entry, raiser=self._raiser, xtend=self._xtend)
        return BooleanGH(r, name=op.__name__)

    def _inplace(self, r):
        # Replace this with a L{BooleanGH} result.
        self._first, self._last, self._len = \
           r._first,    r._last,    r._len
        r._first = r._last = r._len = None
        return self

    def _other(self, other):
        if isinstance(other, BooleanGH):
            return other
        raise TypeError('not %s: %r' % (BooleanGH.__name__, other))

    def __and__(self, other):
        '''Intersection: C{this & other}.
        '''
        return self._boolean(other, False, False, self.__and__)

    def __iand__(self, other):
        '''In-place intersection: C{this &= other}.
        '''
        return self._inplace(self.__and__(other))

    def __ior__(self, other):
        '''In-place union: C{this |= other}.
        '''
        return self._inplace(self.__or__(other))

    def __isub__(self, other):
        '''In-place difference: C{this -= other}.
        '''
        return self._inplace(self.__sub__(other))

    def __or__(self, other):
        '''Union: C{this | other}.
        '''
        return self._boolean(other, True, True, self.__or__)

    def __rand__(self, other):
        ''' Reverse intersection: C{other & this}
        '''
        return self._other(other).__and__(self)

    def __ror__(self, other):
        ''' Reverse union: C{other | this}
        '''
        return self._other(other).__or__(self)

    def __rsub__(self, other):
        ''' Reverse difference: C{other - this}
        '''
        return self._other(other).__sub__(self)

    def __sub__(self, other):
        '''Difference: C{this - other}.
        '''
        return self._boolean(other, True, False, self.__sub__)


class ClipGH4Tuple(_NamedTuple):
    '''4-Tuple C{(lat, lon, height, clipid)} for each point of the
       L{clipGH4} result with the C{lat}-, C{lon}gitude, C{height}
       and C{clipid} of the polygon or polygonal clip.

       @note: The C{height} is a L{HeightGH} instance if this is
              an intersection, otherwise a L{Height} or C{int(0)}.
    '''
    _Names_ = (_lat_, _lon_, _height_, _clipid_)
    _Units_ = ( Lat,   Lon,  _Pass,     Number_)

    @property_RO
    def isintersection(self):
        '''Is this an intersection?
        '''
        return isinstance(self.height, HeightGH)  # PYCHOK attr

    @property_RO
    def ispoint(self):
        '''Is this an original (polygon) point?
        '''
        return not isinstance(self.height, HeightGH)  # PYCHOK attr


def clipGH4(points, corners, raiser=False, xtend=False):
    '''Clip a polygon against a clip region or box using the U{Greiner-Hormann
       <http://www.inf.USI.CH/hormann/papers/Greiner.1998.ECO.pdf>} algorithm,
       U{Correia<https://GitHub.com/helderco/univ-polyclip>}'s implementation
       modified and extended.

       @arg points: The polygon points (C{LatLon}[]).
       @arg corners: Three or more points defining a clip region (C{LatLon}[])
                     or two points to specify a rectangular clip box.
       @kwarg raiser: If C{True} throw L{ClipError} exceptions.
       @kwarg xtend: Cover degenerate cases (C{bool}), but not (yet) based on
                     U{Forster-Hormann-Popa<https://www.ScienceDirect.com/
                     science/article/pii/S259014861930007X>}.

       @return: Yield a L{ClipGH4Tuple}C{(lat, lon, height, clipid)} for each
                clipped point.  The result may consist of several clips, each
                a closed polygon with a unique C{clipid}.

       @raise ClipError: Null C{clipper}, an open polygon clip or an I{unhandled}
                         intersection.

       @see: U{Greiner-Hormann<https://WikiPedia.org/wiki/Greiner–Hormann_clipping_algorithm>},
             U{Ionel Daniel Stroe<https://Davis.WPI.edu/~matt/courses/clipping/>}, U{Forster,
             Hormann and Popa<https://www.ScienceDirect.com/science/article/pii/S259014861930007X>}
             and I{Correia}'s U{univ-polyclip<https://GitHub.com/helderco/univ-polyclip>}.
    '''
    s = _GHcircular(points, xtend=xtend)
    c = _4corners(corners)
    return s._clip(c, False, False, raiser=raiser)


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


__all__ += _ALL_DOCS(_GHcircular)

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
