
# -*- coding: utf-8 -*-

u'''I{Boolean} operations on I{composite} polygons and I{clip}s.

Classes L{BooleanFHP} and L{BooleanGH} are I{composites} and provide
I{boolean} operations C{intersection}, C{difference},
C{reverse-difference} and C{union}.

@note: A I{clip} is defined as a single, usually closed polygon, a
       I{composite} is a collection of one or more I{clip}s.

@see: U{Forster-Hormann-Popa<https://www.ScienceDirect.com/science/
      article/pii/S259014861930007X>} and U{Greiner-Hormann
      <http://www.inf.USI.CH/hormann/papers/Greiner.1998.ECO.pdf>}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import isscalar, isodd
from pygeodesy.constants import EPS, EPS2, MAX, _0_0, _0_5, _1_0
from pygeodesy.errors import ClipError, _IsnotError
from pygeodesy.fmath import fabs, favg, hypot, hypot2
from pygeodesy.fsums import Property_RO, property_RO
from pygeodesy.interns import NN, _BANG_, _clip_, _COMMASPACE_, _DOT_, \
                             _e_, _ELLIPSIS_, _few_, _height_, _lat_, \
                             _lon_, _name_, _scalar_, _SPACE_, _too_, \
                             _X_, _x_,  _B_, _d_, _R_  # PYCHOK used!
from pygeodesy.lazily import _ALL_LAZY, _ALL_DOCS, _FOR_DOCS
from pygeodesy.named import Fmt, _Named, pairs
# from pygeodesy.props import Property_RO, property_RO  # from .fsums
# from pygeodesy.streprs import Fmt, pairs  # from .named
from pygeodesy.units import Height, HeightX

# from math import fabs  # from .fmath

__all__ = _ALL_LAZY.booleans
__version__ = '23.03.09'

_0_EPS =  EPS    # near-zero, positive
_EPS_0 = -EPS    # near-zero, negative
_1_EPS = _1_0 + EPS  # over-one
_EPS_1 = _1_0 - EPS  # near-one
_10EPS =  EPS * 10   # see ._2Abs

_alpha_   = 'alpha'
_boolean_ = 'boolean'
_clipid_  = 'clipid'
_corners_ = 'corners'
_open_    = 'open'


def _Enum(txt, enum):  # PYCHOK unused
    return txt  # NN(txt, _TILDE_, str(enum))


class _L(object):
    # Intersection labels
    CROSSING      = _Enum(_X_, 1)  # C++ enum
    CROSSING_D    = _Enum(_X_ + _d_, 8)
    CROSSINGs     = (CROSSING, CROSSING_D)
    BOUNCING      = _Enum(_B_, 2)
    BOUNCING_D    = _Enum(_B_ + _d_, 9)
    BOUNCINGs     = (BOUNCING, BOUNCING_D) + CROSSINGs
    LEFT_ON       = _Enum('Lo', 3)
    ON_ON         = _Enum('oo', 5)
    ON_LEFT       = _Enum('oL', 6)
    ON_RIGHT      = _Enum('oR', 7)
    RIGHT_ON      = _Enum('Ro', 4)
    RIGHT_LEFT_ON = (RIGHT_ON, LEFT_ON)
    # Entry/Exit flags
    ENTRY         = _Enum(_e_, 1)
    EXIT          = _Enum(_x_, 0)
    Toggle        = {ENTRY: EXIT,
                     EXIT: ENTRY,
                     None: None}

_L = _L()  # PYCHOK singleton


class _RP(object):
    # RelativePosition types
    IS_Pm = _Enum('Pm', 2)  # C++ enum
    IS_Pp = _Enum('Pp', 3)
    LEFT  = _Enum('L',  0)
    RIGHT = _Enum(_R_,  1)

_RP = _RP()  # PYCHOK singleton

_RP2L = {(_RP.LEFT,  _RP.RIGHT): _L.CROSSING,
         (_RP.RIGHT, _RP.LEFT):  _L.CROSSING,
         (_RP.LEFT,  _RP.LEFT):  _L.BOUNCING,
         (_RP.RIGHT, _RP.RIGHT): _L.BOUNCING,
         # overlapping cases
         (_RP.RIGHT, _RP.IS_Pp): _L.LEFT_ON,
         (_RP.IS_Pp, _RP.RIGHT): _L.LEFT_ON,
         (_RP.LEFT,  _RP.IS_Pp): _L.RIGHT_ON,
         (_RP.IS_Pp, _RP.LEFT):  _L.RIGHT_ON,
         (_RP.IS_Pm, _RP.IS_Pp): _L.ON_ON,
         (_RP.IS_Pp, _RP.IS_Pm): _L.ON_ON,
         (_RP.IS_Pm, _RP.RIGHT): _L.ON_LEFT,
         (_RP.RIGHT, _RP.IS_Pm): _L.ON_LEFT,
         (_RP.LEFT,  _RP.IS_Pm): _L.ON_RIGHT,
         (_RP.IS_Pm, _RP.LEFT):  _L.ON_RIGHT}


def _alpha4(a):
    # Return 4-tuple (alpha, abs(alpha) near 0, 0 < alpha < 1, not 0 < alpha < 1)
    return (a, False, True,  False) if _0_EPS < a < _EPS_1 else (
           (a, False, False, True)  if _0_EPS < fabs(a) else
           (a, True,  False, False))


def _min_max_eps(*xs):
    '''(INTERNAL) Return 2-tuple C{(min, max)}, oversized.
    '''
    lo, hi = min(xs), max(xs)
    lo *= _1_EPS if lo < 0 else _EPS_1
    hi *= _EPS_1 if hi < 0 else _1_EPS
    return (lo or _EPS_0), (hi or _0_EPS)


def _outside(x1, x2, lo, hi):
    '''(INTERNAL) Is C{(x1, x2)} outside C{(lo, hi)}?
    '''
    return max(x1, x2) < lo or min(x1, x2) > hi


class _LatLonBase(_Named):
    '''(INTERNAL) Base class for L{LatLonFHP} and L{LatLonGH}.
    '''
    _alpha   = None       # point AND intersection else length
    _checked = False      # checked in phase 3 iff intersection
    _clipid  = 0          # polygonal clip identifier, number
    _dupof   = None       # original of a duplicate
#   _e_x_str = NN         # shut up PyChecker
    _height  = Height(0)  # interpolated height, usually meter
    _linked  = None       # link to neighbor iff intersection
    _next    = None       # link to the next vertex
    _prev    = None       # link to the previous vertex

    def __init__(self, lat_ll, lon=None, height=0, clipid=0, **name):
        '''New C{LatLon[FHP|GH]} from separate C{lat}, C{lon}, C{height}
           and C{clipid} scalars or from a previous C{LatLon[FHP|GH]},
           a C{Clip[FHP|GH]4Tuple} or some other C{LatLon} instance.

           @arg lat_ll: Latitude (C{scalar}) or a lat/longitude
                        (C{LatLon[FHP|GH]}, aC{Clip[FHP|GH]4Tuple}
                        or some other C{LatLon}).
           @kwarg lon: Longitude (C{scalar}), iff B{C{lat_ll}} is
                       scalar, ignored otherwise.
           @kwarg height: Height (C{scalar}), conventionally C{meter}.
           @kwarg clipid: Clip identifier (C{int}).
           @kwarg name: Optional name (C{str}).
        '''
        if lon is None:
            self.y, self.x = lat_ll.lat, lat_ll.lon
            h = getattr(lat_ll, _height_, height)
            c = getattr(lat_ll, _clipid_, clipid)
        else:
            self.y, self.x = lat_ll, lon
            h, c = height, clipid
        # don't override defaults
        if self._height != h:
            self._height = h
        if self._clipid != c:
            self._clipid = c
        if name:
            self.name = name.get(_name_, NN)

    def __abs__(self):
        return max(fabs(self.x), fabs(self.y))

    def __eq__(self, other):
        return other is self or (other and other.x == self.x
                                       and other.y == self.y)

    def __repr__(self):
        '''String C{repr} of this lat-/longitude.
        '''
        t = _ELLIPSIS_(self._prev, self._next)
        return _SPACE_(self, Fmt.ANGLE(t))

    def __str__(self):
        '''String C{str} of this lat-/longitude.
        '''
        t = (_lat_, self.lat), (_lon_, self.lon)
        if self._height:
            X  = _X_ if self.isintersection else NN
            t += (_height_ + X, self._height),
        if self._clipid:
            t += (_clipid_, self._clipid),
        if self._alpha is not None:
            t += (_alpha_, self._alpha),
#       if self._dupof:  # recursion risk
#           t += (_dupof_, self._dupof.name),
        t = pairs(t, prec=8, fmt=Fmt.g, ints=True)
        t = Fmt.PAREN(_COMMASPACE_.join(t))
        if self._linked:
            k = _DOT_ if self._checked else _BANG_
            t += self._e_x_str(k)  # PYCHOK expected
        return self.name + t

    def __sub__(self, other):
        return self.__class__(self.y - other.y, self.x - other.x)

    def _2A(self, p2, p3):
        # I{Signed} area of a triangle, I{doubled}.
        x, y = self.x, self.y
        return (p2.x - x) * (p3.y - y) - \
               (p3.x - x) * (p2.y - y)

    def _2Abs(self, p2, p3, eps=_10EPS):
        # I{Unsigned} area of a triangle, I{doubled}
        # or 0 if below th given threshold C{eps}.
        a = fabs(self._2A(p2, p3))
        return 0 if a < eps else a

    @property_RO
    def clipid(self):
        '''Get the I{clipid} (C{int} or C{0}).
        '''
        return self._clipid

    @property_RO
    def height(self):
        '''Get the I{height} (C{Height} or C{int}).
        '''
        return self._height

    @property_RO
    def isintersection(self):
        '''Is this an intersection?
        '''
        return bool(self._linked)

    @property_RO
    def ispoint(self):
        '''Is this an original (boolean) point?
        '''
        return self._alpha is None

    @property_RO
    def lat(self):
        return self.y

    def _link(self, other):
        # Make this and an other point are neighbors.
        # assert isinstance(other, self.__class__)
        self._linked = other
        other._linked = self

    @property_RO
    def lon(self):
        return self.x

    def _toClas4(self, Clas, cid):
        # Return this vertex as a L{[Clip|LatLon][FHP|GH]4Tuple}.
        h = self._height
        if self.isintersection:
            h = HeightX(h)
        elif h:
            h = Height(h)
        return Clas(self.lat, self.lon, h, cid)


class LatLonFHP(_LatLonBase):
    '''A point or intersection in a L{BooleanFHP} clip.
    '''
    _en_ex  = None
    _label  = None
#   _prep2  = None, None  # shup up PyChecker
    _2split = None  # or C{._Clip}
    _2xing  = False

    def __init__(self, lat_ll, *lon_h_clipid, **name):
        '''New C{LatLonFHP} from separate C{lat}, C{lon}, C{h}eight
           and C{clipid} scalars, or from a previous L{LatLonFHP},
           a L{ClipFHP4Tuple} or some other C{LatLon} instance.

           @arg lat_ll: Latitude (C{scalar}) or a lat/longitude
                        (L{LatLonFHP}, C{LatLon} or L{ClipFHP4Tuple}).
           @arg lon_h_clipid: Longitude (C{scalar}), C{h}eight and
                              C{clipid} iff B{C{lat_ll}} is scalar,
                              ignored otherwise.
           @kwarg name: Optional name (C{str}).
        '''
        _LatLonBase.__init__(self, lat_ll, *lon_h_clipid, **name)

    def __add__(self, other):
        return self.__class__(self.y + other.y, self.x + other.x)

    def __mod__(self, other):  # cross product
        return self.x * other.y - self.y * other.x

    def __mul__(self, other):  # dot product
        if not isinstance(other, self.__class__):
            raise _IsnotError(self.__class__.__name__, other=other)
        return self.x * other.x + self.y * other.y

    def __rmul__(self, other):  # scalar product
        if not isscalar(other):
            raise _IsnotError(_scalar_, other=other)
        return self.__class__(self.y * other, self.x * other)

    def _e_x_str(self, t):
        if self._label:
            t = NN(self._label, t)
        if self._en_ex:
            t = NN(t, self._en_ex)
        return t

    @property_RO
    def _isduplicate(self):
        # Is this point a I{duplicate} intersection?
        p = self._dupof
        return bool(p and self._linked
                      and p is not self
                      and p == self
#                     and p._alpha in (None, self._alpha)
                      and self._alpha in (_0_0, p._alpha))

#   @property_RO
#   def _isduplicated(self):
#       # Return the number of I{duplicates}?
#       d, v = 0, self
#       while v:
#           if v._dupof is self:
#               d += 1
#           v = v._next
#           if v is self:
#              break
#       return d

    def _isinside(self, composite, *excluded):
        # Is this point inside a composite,
        # excluding certain C{_Clip}s.
        x, y, i = self.x, self.y, False
        for c in composite._clips:
            if c not in excluded:
                w = 0
                for p1, p2 in c._edges2():
                    # edge [p1,p2] must straddle y
                    if (p1.y < y) is not (p2.y < y):
                        r = p2.x > x
                        s = p2.y > p1.y
                        if p1.x < x:
                            b = r and (s is (p1._2A(p2, self) > 0))
                        else:
                            b = r or  (s is (p1._2A(p2, self) > 0))
                        if b:
                            w += 1 if s else -1
                if isodd(w):
                    i = not i
        return i

    def isinside(self, *composites):
        '''Is this point inside I{combined} composites based on C{winding number}?

           @arg composites: One or more iterables or composites of clips and points
                            (L{ClipFHP4Tuple}, L{ClipGH4Tuple}, L{LatLonFHP},
                            L{LatLonGH}, L{LatLon_} or any other C{LatLon}).

           @see: U{Algorithm 6<https://www.ScienceDirect.com/science/article/pii/
                 S0925772101000128>}.
        '''
        self._isinside(_CompositeEdges(self.__class__, composites))

    @property_RO
    def isintersection(self):
        '''Is this an intersection?  May be C{ispoint} too!
        '''
        return bool(self._linked)

    @property_RO
    def ispoint(self):
        '''Is this an I{original} point?  May be C{isintersection} too!
        '''
        return self._alpha is None

    @property_RO
    def _prev_next2(self):
        # Adjust 2-tuple (._prev, ._next) iff a I{duplicate} intersection
        p, n = self._prev, self._next
        if self._isduplicate:
            p = self._dupof
            while p._isduplicate:
                p = p._dupof
            p = p._prev
            while n._isduplicate:
                n = n._next
        return p, n

#   def _edge2(self):
#       # Return the start and end point of the
#       # edge containing I{intersection} C{v}.
#       n = p = self
#       while p.isintersection:
#           p = p._prev
#           if p is self:
#               break
#       while n.isintersection:
#           n = n._next
#           if n is self:
#               break
#       # assert p == self or not p._2Abs(self, n)
#       return p, n

    def _RPoracle(self, p1, p2, p3):
        # Relative Position oracle
        if p1._linked is self:  # or p1._linked2(self):
            T = _RP.IS_Pm
        elif p3._linked is self:  # or p3._linked2(self):
            T = _RP.IS_Pp
        elif p1._2A(p2, p3) > 0:  # left turn
            T = _RP.LEFT if self._2A(p1, p2) > 0 and \
                            self._2A(p2, p3) > 0 else \
                _RP.RIGHT  # PYCHOK indent
        else:  # right turn (or straight)
            T = _RP.RIGHT if self._2A(p1, p2) < 0 and \
                             self._2A(p2, p3) < 0 else \
                _RP.LEFT  # PYCHOK indent
        return T


class LatLonGH(_LatLonBase):
    '''A point or intersection in a L{BooleanGH} clip.
    '''
    _entry  = None   # entry or exit iff intersection
    _extend = False

    def __init__(self, lat_ll, *lon_h_clipid, **name):
        '''New C{LatLonGH} from separate C{lat}, C{lon}, C{h}eight
           and C{clipid} scalars, or from a previous L{LatLonGH},
           a L{ClipGH4Tuple} or some other C{LatLon} instance.

           @arg lat_ll: Latitude (C{scalar}) or a lat/longitude
                        (L{LatLonGH}, C{LatLon} or L{ClipGH4Tuple}).
           @arg lon_h_clipid: Longitude (C{scalar}), C{h}eight and
                              C{clipid} iff B{C{lat_ll}} is scalar,
                              ignored otherwise.
           @kwarg name: Optional name (C{str}).
        '''
        _LatLonBase.__init__(self, lat_ll, *lon_h_clipid, **name)

    def _check(self):
        # Check-mark this vertex and its link.
        self._checked = True
        b = self._linked
        if b and not b._checked:
            b._check()

    def _e_x_str(self, t):
        return NN(t, (NN if self._entry is None else (
                     _e_ if self._entry else _x_)))

    def _isinside(self, composite, *bottom_top):
        # Is this vertex inside the composite? I{Odd-even rule}.

        # The I{odd-even} rule counts the number of edges
        # intersecting a ray emitted East-bound from this
        # point to infinity.  When I{odd} this point lies
        # inside, if I{even} outside.
        r, y = False, self.y
        if not (bottom_top and _outside(y, y, *bottom_top)):
            e   =  self.__class__(MAX, y, clipid=self.clipid)
            _i4 = _EdgeGH(self, e)._intersect4
            for p1, p2, _ in composite._edges3():
                for _ in _i4(p1, p2, False):
                    r = not r
        return r

    def isinside(self, *composites):
        '''Is this point inside I{combined} composites based on C{odd-even rule}?

           @arg composites: One or more iterables or composites of clips and points
                            (L{ClipFHP4Tuple}, L{ClipGH4Tuple}, L{LatLonFHP},
                            L{LatLonGH}, L{LatLon_} or any other C{LatLon}).
        '''
        self._isinside(_CompositeEdges(self.__class__, composites))


class _Clip(_Named):
    '''(INTERNAL) A I{doubly-linked} list representing a I{closed}
       polygon of C{LatLon} points and intersections with other
       polygons in the same C{_Composite}.
    '''
    _composite = None
    _dups      = 0
    _first     = None
    _id        = 0
    _identical = False
    _noInters  = False
    _last      = None
    _LL        = None
    _len       = 0
    _pushback  = False
    _Xs        = 0  # number of intersections

    def __init__(self, composite, clipid=0):
        # assert isinstance(composite, _CompositeBase)
        self._composite  = composite
        self._id         = clipid
        self._LL         = composite._LL
        composite._clips = composite._clips + (self,)

    def __contains__(self, p):
        # Is C{p} one of this clip's vertices?
        for v in self:
            if v is p:  # or ==?
                return True
        return False

    def __iter__(self):
        # Yield all points and intersections.
        v = f = self._first
        while v:
            yield v
            v = v._next
            if v is f:
                break

    def __len__(self):
        # Return the number of latlons.
        return self._len

    _all = __iter__

    @property_RO
    def _all_ON_ON(self):
        # Check whether all verices are ON_ON.
        L_ON_ON = _L.ON_ON
        for v in self:
            if v._label is not L_ON_ON:
                return False
        return True

    def _append(self, y_v, *x_h_clipid):
        # Append a point given as C{y}, C{x}, C{h}eight and C{clipid}
        # args or as a C{LatLon[FHP|GH]} or C{Clip[FHP|GH}4Tuple}.
        self._last = v = self._LL(y_v, *x_h_clipid) if x_h_clipid else y_v
        self._len += 1
        # assert v._clipid == self._id

        v._next = n = self._first
        if n is None:  # set _first
            self._first = p = n = v
        else:  # insert before _first
            v._prev = p = n._prev
        p._next = n._prev = v
        return v

#   def _appendedup(self, v, cid=0):
#       # Like C{._append}, but only append C{v} if not a
#       # duplicate of the one previously append[edup]'ed.
#       y, x, p = v.y, v.x, self._last
#       if p is None or y != p.y or x != p.x or cid != p._clipid:
#           p = self._append(y, x, v.height, cid)
#           if v._linked:
#               p._linked = True  # to force errors
#       return p

    def _closed(self, unused):  # raiser
        # End a clip, close it and check length.
        p, f = self._last, self._first
        if f and f._prev is p and p == f and \
                 p._next is f and p is not f:
            # "un-close" the clip
            f._prev = p = p._prev
            p._next = f
            self._len -= 1
#       elif f and raiser:
#           raise self._OpenClipError(p, f)
        if len(self) < 3:
            raise self._Error(_too_(_few_))

#   def _dup(self, q, *alpha):
#       # Duplicate vertex, an intersection iff C{alpha}
#       # is non-negative C{float} and less than _1_0.
#       return self._insert(q.x, q.y, q, q._next, *alpha)

    def _dup2(self, q):
        # Duplicate a point (or intersection) as intersection.
        v = self._insert(q.x, q.y, q, q._next)
        v._alpha = q._alpha or _0_0  # _0_0 replaces None
        v._dupof = q._dupof or q
        # assert v._prev is q
        # assert q._next is v
        return v

    def _edges2(self, raiser=False, **unused):
        # Yield each I{original} edge as a 2-tuple
        # C{(LatLon[FHP|GH], LatLon[FHP|GH])}.
        p1 = p2 = f = self._first
        while p2:
            p2 = p2._next
            if p2.ispoint:
                yield p1, p2
                p1 = p2
            if p2 is f:
                break
        if raiser and p2 is not f:
            raise self._OpenClipError(p2, f)

    def _Error(self, txt):  # PYCHOK no cover
        # Build a C{ClipError} instance
        kwds = dict(len=len(self), txt=txt)
        if self._dups:
            kwds.update(dups=self._dups)
        cp = self._composite
        if self._id:
            try:
                i = cp._clips.index(self)
                if i != self._id:
                    kwds[_clip_] = i
            except ValueError:
                pass
            kwds[_clipid_] = self._id
        return ClipError(cp.name, cp._kind, **kwds)

    def _insert(self, x, y, start, end, *alpha):
        # insertVertex between points C{start} and
        # C{end}, ordered by C{alpha} iff given.
        v = self._LL(y, x, start._height, start._clipid)
        n = start._next
        if alpha:
            v._alpha  = alpha = alpha[0]
            v._height = favg(v.height, end.height, f=alpha)
            # assert start is not end
            while n is not end and n._alpha < alpha:
                n = n._next
        v._next = n
        v._prev = p = n._prev
        p._next = n._prev = v
        self._len += 1
        return v

    def _intersection(self, unused, q, *p1_p2_alpha):
        # insert an intersection or duplicate
        if p1_p2_alpha:  # intersection on edge
            v = self._insert(q.x, q.y, *p1_p2_alpha)
        else:  # intersection at point
            v = q
            # assert not v._linked
            # assert v._alpha is None
        return v

    def _intersections(self):
        # Yield all intersections.
        for v in self:
            if v.isintersection:
                yield v

    @Property_RO
    def _ishole(self):
        # Is this clip a hole inside its composite?
        v = self._first
        return v._isinside(self._composite, self) if v else False

    def _noXings(self, Union):
        # Are all intersections non-CROSSINGs?
        Ls = _L.BOUNCINGs if Union else _L.CROSSINGs
        for v in self._intersections():
            if v._label in Ls:
                return False
        return True

    def _OpenClipError(self, s, e):  # PYCHOK no cover
        # Return a C{CloseError} instance
        t = NN(s, _ELLIPSIS_(_COMMASPACE_, e))
        return self._Error(_SPACE_(_open_, t))

    def _point2(self, insert):
        # getNonIntersectionPoint and -Vertex
        if not (insert and self._noInters):
            for p in self._points():
                if not p.isintersection:  # or p._isduplicated?
                    return p, None
            for n in self._intersections():
                p, _ = n._prev_next2
                k = p._linked
                if k:
                    if n._linked not in k._prev_next2:
                        # create a pseudo-point
                        k = _0_5 * (p + n)
                        if insert:
                            k = self._insert(k.x, k.y, n._prev, n)
                            r = k
                        else:  # no ._prev, ._next
                            k._clipid = n._clipid
                            r = None
                        return k, r
        return None, None

    def _points(self):
        # Yield all points in original order.
        for v in self:
            if v.ispoint:
                yield v

    def _remove2(self, v):
        # Remove C{v}.
        # n = v._isduplicated
        # if n:
        #     raise NotImplementedError(dups=n, v=v)
        if self._len > 1:
            p = v._prev
            p._next = n = v._next
            n._prev = p
        else:
            p = n = None
        if self._first is v:
            self._first = n
        if self._last is v:
            self._last = n
        self._len -= 1
        return p, n

    def _Xings(self):
        # Yield all I{un-checked} CROSSING intersections.
        CROSSING = _L.CROSSING
        for v in self._intersections():
            if v._label is CROSSING and not v._checked:
                yield v


class _ClipEdges(_Clip):  # PYCHOK no cover
    # Wrapper yielding edges, eliminating
    # null edges and adding closing edge.

    def __init__(self, cps, LL):  # PYCHOK signature
        self._cps = cps  # _CompositeEdges
        self._LL  = LL

    def _edges2(self, *unused):  # PYCHOK signature
        # Yield C{LatLon[FHP|GH]} point pairs.
        for cp in self._cps:
            self._composite = cp
            p2 = s = None
            for ll in cp:
                p1, p2 = p2, self._LL(ll)
                if s is None:
                    s = p1 = p2
                elif p1._clipid != p2._clipid:
                    if p1 != s:
                        yield p1, s
                    s = p1 = p2
                elif p1 != p2:
                    yield p1, p2
                self._id = p2._clipid
            if p1 != s:
                yield p1, s


class _CompositeBase(_Named):
    '''(INTERNAL) A list of C{_Clips} representing a
       composite polygon of C{LatLon} points and
       mutual intersections.
    '''
    _clips  =  ()   # tuple of C{_Clips}
    _eps    =  EPS  # to pass to Q
    _kind   = _corners_
    _LL     = _LatLonBase  # shut up PyChecker
    _raiser =  False
    _xtend  =  False

    def __init__(self, lls, name=NN, kind=NN, eps=EPS):
        # New L{BooleanFHP} or L{BooleanGH}.
        n = name or getattr(lls, _name_, NN)
        if n:
            self.name = n
        if kind:
            self._kind = kind
        if self._eps != eps:
            self._eps = eps

        c = _Clip(self)
        lp = None
        for ll in lls:
            ll = self._LL(ll)
            if lp is None:
                c._id = ll._clipid  # keep clipid
                lp = c._append(ll)
            elif ll._clipid != lp._clipid:  # new clip
                c._closed(self._raiser)
                c = _Clip(self, ll._clipid)
                lp = c._append(ll)
            elif abs(ll - lp) > eps:  # PYCHOK lp
                lp = c._append(ll)
            else:
                c._dups += 1
        c._closed(self._raiser)

    def __contains__(self, v):  # PYCHOK no cover
        # Is C{p} one of the clips' points?
        for c in self._clips:
            if v in c:
                return True
        return False

    def __iter__(self):
        # Yield all points and intersections.
        for c in self._clips:
            for v in c:
                yield v

    def __len__(self):
        # Return the total number of latlons.
        return sum(map(len, self._clips)) if self._clips else 0

    def __repr__(self):
        '''String C{repr} of this composite.
        '''
        c = len(self._clips)
        c = Fmt.SQUARE(c) if c > 1 else NN
        n = Fmt.SQUARE(len(self))
        t = Fmt.PAREN(self)
        return NN(self.__class__.__name__, c, n, t)

    def __str__(self):
        '''String C{str} of this composite.
        '''
        return _COMMASPACE_.join(map(str, self))

    @property_RO
    def _bottom_top_eps(self):
        # Get the bottom and top C{y} bounds, oversized.
        return _min_max_eps(min(v.y for v in self),
                            max(v.y for v in self))

    def _class(self, corners, kwds, **dflts):
        # Return a new instance
        _g = kwds.get
        kwds = dict((n, _g(n, v)) for n, v in dflts.items())
        return self.__class__(corners or (), **kwds)

    @property_RO
    def _clipids(self):  # PYCHOK no cover
        for c in self._clips:
            yield c._id

    def clipids(self):
        '''Return a tuple with all C{clipid}s, I{ordered}.
        '''
        return tuple(self._clipids)

    def _edges3(self, **raiser):
        # Yield each I{original} edge as a 3-tuple
        # C{(LatLon[FHP|GH], LatLon[FHP|GH], _Clip)}.
        for c in self._clips:
            for p1, p2 in c._edges2(**raiser):
                yield p1, p2, c

    def _intersections(self):
        # Yield all intersections.
        for c in self._clips:
            for v in c._intersections():
                yield v

    @property_RO
    def _left_right_eps(self):
        # Get the left and right C{x} bounds, oversized.
        return _min_max_eps(min(v.x for v in self),
                            max(v.x for v in self))

    def _points(self):  # PYCHOK no cover
        # Yield all I{original} points,
        # some may be intersection too.
        for c in self._clips:
            for v in c._points():
                yield v

    def _results(self, _presults, Clas, closed=False, inull=False, **unused):
        # Yield the dedup'd results, as L{ClipFHP4Tuple}s
        C = self._LL if Clas is None else Clas
        for cid, ns in enumerate(_presults):
            f = p = v = None
            for n in ns:
                if f is None:
                    yield n._toClas4(C, cid)
                    f = p = n
                elif v is None:
                    v = n  # got f, p, v
                elif inull or p._2Abs(v, n):
                    yield v._toClas4(C, cid)
                    p, v = v, n
                else:  # null, colinear, ... skipped
                    v = n
            if v and (inull or p._2Abs(v, f)):
                yield v._toClas4(C, cid)
                p = v
            if f and p != f and closed:  # close clip
                yield f._toClas4(C, cid)


class _CompositeEdges(_CompositeBase):  # PYCHOK no cover
    # Polygon wrapper yielding clips and edges,
    # eliminating duplicates and closing open clips.

    def __init__(self, LL, composites):  # PYCHOK signature
        self._cps = composites
        self._LL  = LL

    @property_RO
    def _clips(self):
        # Use one C{_Clip}.
        return (_ClipEdges(self._cps, self._LL),)


class _CompositeFHP(_CompositeBase):
    '''(INTERNAL) A list of C{_Clips} representing a
       composite polygon of L{LatLonFHP} points and
       mutual intersections.
    '''
    _LL    = LatLonFHP
    _Union = False

    def __init__(self, lls, raiser=False, **name_kind_eps):
        '''New L{BooleanFHP}.

           @arg lls: The polygon points (iterable of 2 or more
                     C{LatLon}s, L{LatLonFHP}s or L{ClipFHP4Tuple}s).
           @kwarg raiser: If C{True} throw L{ClipError} exceptions
                          for errors.
           @kwarg name_kind_eps: Optional C{B{name}=NN} (C{str})
                       and/or optional C{B{kind}='corners'} of
                       B{C{lls}} points (C{str}) and C{B{eps}=EPS}
                       tolerance for duplicate removal (C{str}).
        '''
        if raiser:
            self._raiser = True
        _CompositeBase.__init__(self, lls, **name_kind_eps)

    def _classify(self):
        # 2) Classify intersection chains.
        L = _L
        for v in self._intersections():
            n, b = v, v._label
            if b in L.RIGHT_LEFT_ON:  # next chain
                while True:
                    # n.__dict__.pop('_label')
                    n._label = None
                    n = n._next
                    if n is v or n._label is not L.ON_ON:  # n._label and ...
                        break
                a = L.LEFT_ON if n._label is L.ON_LEFT else L.RIGHT_ON
                v._label = n._label = L.BOUNCING_D if a is b else L.CROSSING_D

        # 3) Copy labels
        for v in self._intersections():
            v._linked._label = v._label

    def _clip(self, corners, Union=False, Clas=None,
                           **closed_inull_raiser):
        # Clip this polygon with another one, C{corners},
        # using the Foster/Hormann/Popa's algorithm.
        P = self
        Q = self._class(corners, closed_inull_raiser,
                                 eps=P._eps, raiser=False)
        P._reset(Union, name='P')
        Q._reset(Union, name='Q')

        bt = Q._bottom_top_eps
        lr = Q._left_right_eps
        # compute and insert intersections
        for p1, p2, Pc in P._edges3(**closed_inull_raiser):
            if not (_outside(p1.x, p2.x, *lr) or
                    _outside(p1.y, p2.y, *bt)):
                e = _EdgeFHP(p1, p2)
                if e._dp2 > EPS2:  # non-null edge
                    for q1, q2, Qc in Q._edges3(**closed_inull_raiser):
                        for T, p, q in e._intersect3(q1, q2):
                            p = Pc._intersection(T, *p)
                            q = Qc._intersection(T, *q)
                            # assert not p._linked
                            # assert not q._linked
                            p._link(q)

        # label and classify intersections
        P._labelize()
        P._classify()

        # check for special cases
        P._special_cases(Q)
        Q._special_cases(P)
        # handle identicals
        P._special_identicals(Q)

        # set Entry/Exit flags
        P._set_entry_exits(Q)
        Q._set_entry_exits(P)

        # handle splits and crossings
        P._splits_xings(Q)

        # yield the results
        return P._results(self._presults(Q), Clas, **closed_inull_raiser)

    @property_RO
    def _identicals(self):
        # Yield all clips marked C{._identical}.
        for c in self._clips:
            if c._identical:
                yield c

    def _labelize(self):
        # 1) Intersections classification
        for p in self._intersections():
            q = p._linked
            # determine local configuration at this intersection
            # and positions of Q- and Q+ relative to (P-, I, P+)
            p1, p3 = p._prev_next2
            q1, q3 = q._prev_next2
            t = (q1._RPoracle(p1, p, p3),
                 q3._RPoracle(p1, p, p3))
            # check intersecting and overlapping cases
            p._label = _RP2L.get(t, None)

    def _presults(self, other):
        # Yield the result clips, each as
        # a generator of L{_LatLonFHP}s
        for cp in (self, other):
            for c in cp._clips:
                if c._pushback:
                    yield c._all()
        for c in self._clips:
            for X in c._Xings():
                yield self._resultX(X)

    def _reset(self, Union=False, **unused):
        if Union:
            self._Union = True

    def _resultX(self, X):
        # Yield the result from an unchecked CROSSING.
        L, U, v = _L, self._Union, X
        while v:  # and not v._checked:
            # yield v
            v._checked = True
            r = v  # in P or Q
            s = L.Toggle[v._en_ex]
            e = (s is L.EXIT) ^ U
            while True:
                v = v._next if e else v._prev
                yield v
                v._checked = True
                if v._en_ex is s or v is X:
                    break
                if v is r:  # full circle
                    raise ClipError(full_circle=v, clipid=v._clipid)
            if v is not X:
                v = v._linked
                v._checked = True
            if v is X:
                break

    def _set_entry_exits(self, other):  # MCCABE 17
        # 4) Set entry/exit flags
        L, U = _L, self._Union
        for c in self._clips:
            f, k = c._point2(True)
            if f:
                n = f
                s = L.EXIT if n._isinside(other) else L.ENTRY
                t = L.EXIT  # first_chain_vertex = True
                while True:
                    if n.isintersection:
                        b = n._label
                        if b is L.CROSSING:
                            n._en_ex = s
                            s = L.Toggle[s]
                        elif b is L.BOUNCING and ((s is L.EXIT) ^ U):
                            n._2split = c  # see ._splits_xings
                        elif b is L.CROSSING_D:
                            n._en_ex = s
                            if (s is t) ^ U:
                                n._label = L.CROSSING
                            t = L.Toggle[t]
                            if t is L.EXIT:
                                s = L.Toggle[s]
                        elif b is L.BOUNCING_D:
                            n._en_ex = s
                            if (s is t) ^ U:
                                n._2xing = c  # see ._splits_xings
                            s = L.Toggle[s]
                            t = L.Toggle[t]
                    n = n._next  # _, n = n._prev_next2
                    if n is f:
                        break  # PYCHOK attr?
                if k:
                    c._remove2(k)

    def _special_cases(self, other):
        # 3.5) Check special cases
        U = self._Union
        for c in self._clips:
            if c._noXings(U):
                c._noInters = True
                if c._all_ON_ON:
                    c._identical = True
                else:
                    p, _ = c._point2(False)
                    if p and (p._isinside(other) ^ U):
                        c._pushback = True

    def _special_identicals(self, other):
        # 3.5) Handle identicals
        cds = dict((c._id, c) for c in other._identicals)
        if cds:
            for c in self._identicals:
                for v in c._intersections():
                    d = cds.get(v._linked._clipid, None)
                    if d and d._ishole is c._ishole:
                        c._pushback = True
                        break  # next clip

    @property_RO
    def _2splits(self):
        # Yield all intersections marked C{._2split}
        for p in self._intersections():
            if p._2split:
                # assert isinstance(p._2split, _Clip)
                yield p

    def _splits_xings(self, other):  # MCCABE 15
        # 5) Handle split pairs and 6) crossing candidates

        def _2A_dup2(p, P):
            p1, p2 = p._prev_next2
            a2 = p1._2A(p, p2)
            Pc = p._2split
            # assert Pc in P._clips
            return a2, Pc._dup2(p)

        def _links2(ps, qs):  # PYCHOK P unused?
            # Yield each link as a 2-tuple(p, q)
            id_qs = set(map(id, qs))
            if id_qs:
                for p in ps:
                    q = p._linked
                    if id(q) in id_qs:
                        yield p, q

        L = _L
        E =  L.ENTRY if self._Union else L.EXIT
        X =  L.Toggle[E]
        for p, q in _links2(self._2splits, other._2splits):
            ap, pp = _2A_dup2(p, self)
            aq, qq = _2A_dup2(q, other)
            if (ap * aq) > 0:
                # overwrites p-q link
                p._link(qq)
                q._link(pp)
            else:
                pp._link(qq)
            p._en_ex  = q._en_ex  = E
            pp._en_ex = qq._en_ex = X
            p._label  = pp._label = \
            q._label  = qq._label = L.CROSSING

        for p, q in _links2(self._2xings, other._2xings):
            p._label = q._label = L.CROSSING

    @property_RO
    def _2xings(self):
        # Yield all intersections marked C{._2xing}
        for p in self._intersections():
            if p._2xing:
                yield p


class _CompositeGH(_CompositeBase):
    '''(INTERNAL) A I{circular, doubly-linked} list of L{LatLonGH}s,
       representing the original points of a polygon plus -if any-
       the intersections with another.
    '''
    _LL    = LatLonGH
    _xtend = False

    def __init__(self, lls, raiser=False, xtend=False, **name_kind_eps):
        '''New L{BooleanGH}.

           @arg lls: The polygon points (iterable of 2 or more
                     C{LatLon}s, L{LatLonGH}s or L{ClipGH4Tuple}s).
           @kwarg name: Optional name (C{str}).
           @kwarg raiser: If C{True} throw C{ClipError} exceptions
                          for errors or I{degenerate cases}.
           @kwarg xtend: If C{True} handle I{degenerate cases}.
           @kwarg name_kind_eps: Optional C{B{name}=NN} (C{str})
                       and/or optional C{B{kind}='corners'} of
                       B{C{lls}} points (C{str}) and C{B{eps}=EPS}
                       tolerance for duplicate removal (C{float}).
        '''
        if xtend:
            self._xtend = True
        elif raiser:
            self._raiser = True
        _CompositeBase.__init__(self, lls, **name_kind_eps)

    def _clip(self, corners, s_entry, c_entry, Clas=None,
                           **closed_inull_raiser_xtend):
        # Clip this polygon with another one, C{corners}.

        # Core of Greiner/Hormann's algorithm, enhanced U{Correia's
        # <https://GitHub.com/helderco/univ-polyclip>} implementation***
        # and extended to optionally handle so-called "degenerate cases"
        S = self
        C = self._class(corners, closed_inull_raiser_xtend,
                                 raiser=False, xtend=False)
        lr = C._left_right_eps
        bt = C._bottom_top_eps
        # 1. find intersections
        for s1, s2, Sc in S._edges3(**closed_inull_raiser_xtend):
            if not (_outside(s1.x, s2.x, *lr) or
                    _outside(s1.y, s2.y, *bt)):
                e = _EdgeGH(s1, s2, **closed_inull_raiser_xtend)
                _, sx, _, sy = e._x_sx_y_sy
                if sx or sy:
                    for c1, c2, Cc in C._edges3(**closed_inull_raiser_xtend):
                        for x, y, sa, ca in e._intersect4(c1, c2):
                            s = Sc._insert(x, y, s1, s2, sa)
                            c = Cc._insert(x, y, c1, c2, ca)
                            s._link(c)

        # 2. identify entry/exit intersections
        if S._first:
            s_entry ^= S._first._isinside(C, *bt)
            for v in S._intersections():
                v._entry = s_entry = not s_entry

        if C._first:
            c_entry ^= C._first._isinside(S)
            for v in C._intersections():
                v._entry = c_entry = not c_entry

        # 3. yield the result(s)
        return self._results(self._presults(), Clas, **closed_inull_raiser_xtend)

    @property_RO
    def _first(self):
        # Get the very first vertex
        for v in self:
            return v
        return None

    def _presults(self):
        # Yield the (first) unchecked, intersection(s).
        for c in self._clips:
            for v in c._intersections():
                if not v._checked:
                    yield self._resultU(v)

    def _resultU(self, v):
        # Yield the result from an un-checked intersection.
        while v and not v._checked:
            v._check()
            yield v
            r = v
            e = v._entry
            while True:
                v = v._next if e else v._prev
                yield v
                if v._linked:
                    break
                if v is r:
                    raise ClipError(full_circle=v, clipid=v._clipid)
            v = v._linked  # switch


class _EdgeFHP(object):
    # An edge between two L{LatLonFHP} points.

    X_INTERSECT = _Enum('Xi', 1)  # C++ enum
    X_OVERLAP   = _Enum('Xo', 5)
    P_INTERSECT = _Enum('Pi', 3)
    P_OVERLAP   = _Enum('Po', 7)
    Ps          = (P_INTERSECT, P_OVERLAP, X_OVERLAP)
    Q_INTERSECT = _Enum('Qi', 2)
    Q_OVERLAP   = _Enum('Qo', 6)
    Qs          = (Q_INTERSECT, Q_OVERLAP, X_OVERLAP)
    V_INTERSECT = _Enum('Vi', 4)
    V_OVERLAP   = _Enum('Vo', 8)
    Vs          = (V_INTERSECT, V_OVERLAP)

    def __init__(self, p1, p2, **unused):
        # New edge between points C{p1} and C{p2}, each a L{LatLonFHP}.
        self._p1_p2 = p1, p2
        self._dp    = dp = p2 - p1
        self._dp2   = dp * dp  # dot product, hypot2

        self._left_right_eps = _min_max_eps(p1.x, p2.x)
        self._bottom_top_eps = _min_max_eps(p1.y, p2.y)

    def _intersect3(self, q1, q2):
        # Return intersection Type or C{None}
        if not (_outside(q1.x, q2.x, *self._left_right_eps) or
                _outside(q1.y, q2.y, *self._bottom_top_eps)):
            T, E = None, _EdgeFHP  # self.__class__
            dq  = q2 - q1
            dq2 = dq * dq  # dot product, hypot2
            if dq2 > EPS2:  # like ._clip
                p1, p2 = self._p1_p2
                ap1   = p1._2A(q1, q2)
                ap2_1 = p2._2A(q1, q2) - ap1
                if fabs(ap2_1) > _0_EPS:  # non-parallel edges
                    aq1   = q1._2A(p1, p2)
                    aq2_1 = q2._2A(p1, p2) - aq1
                    if fabs(aq2_1) > _0_EPS:
                        # compute and classify alpha and beta
                        a, a_0, a_0_1, _ = _alpha4(-ap1 / ap2_1)
                        b, b_0, b_0_1, _ = _alpha4(-aq1 / aq2_1)
                        # distinguish intersection types
                        T = E.X_INTERSECT if a_0_1 and b_0_1 else (
                            E.P_INTERSECT if a_0_1 and b_0   else (
                            E.Q_INTERSECT if a_0   and b_0_1 else (
                            E.V_INTERSECT if a_0   and b_0   else None)))

                elif fabs(ap1) < _0_EPS:  # parallel or colinear edges
                    dp = self._dp
                    d1 = q1 - p1
                    # compute and classify alpha and beta
                    a, a_0, a_0_1, _a_0_1 = _alpha4((d1 * dp) / self._dp2)
                    b, b_0, b_0_1, _b_0_1 = _alpha4((d1 * dq) /     (-dq2))
                    # distinguish overlap type
                    T = E.X_OVERLAP if  a_0_1 and  b_0_1 else (
                        E.P_OVERLAP if  a_0_1 and _b_0_1 else (
                        E.Q_OVERLAP if _a_0_1 and  b_0_1 else (
                        E.V_OVERLAP if  a_0   and  b_0   else None)))

            if T:
                if T is E.X_INTERSECT:
                    v = p1 + a * self._dp
                    yield T, (v, p1, p2, a), (v, q1, q2, b)
                elif T in E.Vs:
                    yield T, (p1,), (q1,)
                else:
                    if T in E.Qs:
                        yield T, (p1,), (p1, q1, q2, b)
                    if T in E.Ps:
                        yield T, (q1, p1, p2, a), (q1,)


class _EdgeGH(object):
    # An edge between two L{LatLonGH} points.

    _raiser = False
    _xtend  = False

    def __init__(self, s1, s2, raiser=False, xtend=False, **unused):
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
        a = (y * dy + x * dx) / self._hypot2
        d = (y * dx - x * dy) / self._hypot0
        return a, fabs(d)

    def _error(self, n, c1, c2):
        t = _SPACE_(self._intersect4.__name__, _EdgeGH(c1, c2))
        raise ClipError('unhandled case %s' % (n,), txt=t)

    @Property_RO
    def _hypot0(self):
        _, sx, _, sy = self._x_sx_y_sy
        return hypot(sx, sy) * _0_EPS

    @Property_RO
    def _hypot2(self):
        _, sx, _, sy = self._x_sx_y_sy
        return hypot2(sx, sy)

    def _intersect4(self, c1, c2, parallel=True):  # MCCABE 14
        # Yield the intersections of this and another edge.

        # @return: None, 1 or 2 intersections, each a 4-Tuple
        #          (x, y, s_alpha, c_alpha) with intersection
        #          x, y and both alphas.

        # @raise ClipError: Intersection unhandled.

        # @see: U{Intersection point of two line segments
        #       <http://PaulBourke.net/geometry/pointlineplane/>}.
        c1_x, c1_y = c1.x, c1.y
        if not (_outside(c1_x, c2.x, *self._left_right_eps) or
                _outside(c1_y, c2.y, *self._bottom_top_eps)):
            x, sx, \
            y, sy = self._x_sx_y_sy

            cx = c2.x - c1_x
            cy = c2.y - c1_y
            d  = cy * sx - cx * sy

            if fabs(d) > _0_EPS:  # non-parallel edges
                dx = x - c1_x
                dy = y - c1_y
                ca = (sx * dy - sy * dx) / d
                if _0_EPS < ca < _EPS_1 or (self._xtend and
                   _EPS_0 < ca < _1_EPS):
                    sa = (cx * dy - cy * dx) / d
                    if _0_EPS < sa < _EPS_1 or (self._xtend and
                       _EPS_0 < sa < _1_EPS):
                        yield (x + sa * sx), (y + sa * sy), sa, ca

                    # unhandled, "degenerate" cases 1, 2 or 3
                    elif self._raiser and not (sa < _EPS_0 or sa > _1_EPS):
                        self._error(1, c1, c2)  # intersection at s1 or s2

                elif self._raiser and not (ca < _EPS_0 or ca > _1_EPS):
                    # intersection at c1 or c2 or at c1 or c2 and s1 or s2
                    sa = (cx * dy - cy * dx) / d
                    e = 2 if sa < _EPS_0 or sa > _1_EPS else 3
                    self._error(e, c1, c2)

            elif parallel and (sx or sy) and (cx or cy):  # PYCHOK no cover
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
                            yield x, y, ca1, self._alpha1(ca)
                        else:  # c1 is between s1 and s2
                            yield (x + sa1 * sx), (y + sa1 * sy), sa1, ca1
                        if sa2 > 1:  # s2 is between c1 and c2
                            ca *= sa2 - _1_0
                            yield (x + sx), (y + sy), ca2, self._alpha1(ca2 - ca)
                        else:  # c2 is between s1 and s2
                            yield (x + sa2 * sx), (y + sa2 * sy), sa2, ca2
                    elif self._raiser and not _outside(sa1, sa2, _0_0, _1_EPS):
                        self._error(4, c1, c2)

#   def _intersect4(self, c1, c2, **unused):
#       # Intersect this and another edge.
#
#       # @return: 4-Tuple (x, y, s_alpha, c_alpha) with the
#       #          intersection point x, y and both alphas.
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
#           if fabs(d) > _0_EPS:  # non-parallel edges
#               dx = x - c1.x
#               dy = y - c1.y
#               sa = (cx * dy - cy * dx) / d
#               ca = (sx * dy - sy * dx) / d
#               if _0_EPS < ca < _EPS_1:
#                   if _0_EPS < sa < _EPS_1:
#                       yield (x + sa * sx), (y + sa * sy), sa, ca
#
#                   # unhandled, "degenerate" cases 1, 2 or 3
#                   elif self.raiser and not (sa < _EPS_0 or sa > _1_EPS):
#                       self._error(1, c1, c2)  # insection at s1 or s2
#
#               elif self.raiser and not (ca < _EPS_0 or ca > _1_EPS):
#                   # intersection at c1 or c2 or at c1 or c2 and s1 or s2
#                   self._error(2 if sa < _EPS_0 or sa > _1_EPS else 3, c1, c2)
#
#           elif self.raiser and (sx or sy) and (cx or cy):
#                   # null, parallel or colinear edges
#                   h = hypot(sx, sy) * _0_EPS
#                   if min(_perpendicular(c1.x - x, c1.y - y, sx, sy),
#                          _perpendicular(c2.x - x, c2.y - y, sx, sy)) < h:
#                       self._error(4, c1, c2)  # colinear, overlapping


class _BooleanBase(object):
    # Shared C{Boolean[FHP|GH]} methods.

    def __and__(self, other):
        '''Intersection: C{this & other}.
        '''
        return self._boolean(other, False, False, self.__and__)  # PYCHOK OK

    def __iand__(self, other):
        '''In-place intersection: C{this &= other}.
        '''
        return self._inplace(self.__and__(other))

    def __ior__(self, other):
        '''In-place union: C{this |= other}.
        '''
        return self._inplace(self.__or__(other))

    def __or__(self, other):
        '''Union: C{this | other}.
        '''
        return self._boolean(other, True, True, self.__or__)  # PYCHOK OK

    def __rand__(self, other):
        ''' Reverse intersection: C{other & this}
        '''
        return self._other(other).__and__(self)

    def __ror__(self, other):
        ''' Reverse union: C{other | this}
        '''
        return self._other(other).__or__(self)

    def _inplace(self, r):
        # Replace this with a L{BooleanFHP} result.
        self._clips, r._clips = r._clips, None
        return self

    def _other(self, other):
        B = self.__class__
        if isinstance(other, B):
            return other
        raise TypeError('not %s: %r' % (B.__name__, other))


class BooleanFHP(_CompositeFHP, _BooleanBase):
    '''I{Composite} polygon class providing I{boolean} operations
       between two C{composite} polygons using U{Forster-Hormann-Popa
       <https://www.ScienceDirect.com/science/article/pii/S259014861930007X>}'s
       C++ implementation transcoded to pure Python.

       The supported operations between (composite) polygon A and B are:

        -  C = A | B  or  A |= B,  union of A and B

        -  C = A & B  or  A &= B,  intersection of A and B

       @see: Function L{clipFHP4} and class L{BooleanGH}.
    '''
    _kind = _boolean_

    if _FOR_DOCS:
        __init__ = _CompositeFHP.__init__

    def _boolean(self, other, Union, unused, op):
        P = BooleanFHP(self, raiser=self._raiser, name=self.name)
        Q = self._other(other)
        r = P._clip(Q, Union=Union, raiser=P._raiser)
        return BooleanFHP(r, raiser=P._raiser, name=op.__name__)

    def _ErrorOp(self, op):
        return NotImplementedError(_DOT_(self.named2, op.__name__))

    def __isub__(self, unused):
        '''In-place difference: C{this -= other}.
        '''
        raise self._ErrorOp(self.__isub__)

    def __rsub__(self, unused):
        ''' Reverse difference: C{other - this}
        '''
        raise self._ErrorOp(self.__rsub__)

    def __sub__(self, unused):
        '''Difference: C{this - other}.
        '''
        raise self._ErrorOp(self.__sub__)


class BooleanGH(_CompositeGH, _BooleanBase):
    '''I{Composite} polygon class providing I{boolean} operations
       between two C{composite} polygons using the U{Greiner-Hormann
       <http://www.inf.USI.CH/hormann/papers/Greiner.1998.ECO.pdf>}
       algorithm, U{Correia<https://GitHub.com/helderco/univ-polyclip>}'s
       implementation modified and extended.

       The supported operations between (composite) polygon A and B are:

        -  C = A | B  or  A |= B,  union of A and B

        -  C = A & B  or  A &= B,  intersection of A and B

        -  C = A - B  or  A -= B,  difference A less B

        -  C = B - A  or  B -= A,  difference B less A

       @note: To handle I{degenerate cases} like C{point-edge} and
              C{point-point} intersections, use class L{BooleanFHP}.

       @see: Function L{clipGH4} and class L{BooleanFHP}.
    '''
    _kind  = _boolean_
    _xtend =  True

    if _FOR_DOCS:
        __init__ = _CompositeGH.__init__

    def _boolean(self, other, s_entry, c_entry, op):
        S = BooleanGH(self, raiser=self._raiser, name=self.name)
        C = self._other(other)
        r = S._clip(C, s_entry, c_entry, raiser=S._raiser, xtend=S._xtend)
        return BooleanGH(r, raiser=S._raiser, name=op.__name__)

    def __isub__(self, other):
        '''In-place difference: C{this -= other}.
        '''
        return self._inplace(self.__sub__(other))

    def __rsub__(self, other):
        ''' Reverse difference: C{other - this}
        '''
        return self._other(other).__sub__(self)

    def __sub__(self, other):
        '''Difference: C{this - other}.
        '''
        return self._boolean(other, True, False, self.__sub__)


__all__ += _ALL_DOCS(_CompositeFHP, _CompositeGH, _LatLonBase)

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
