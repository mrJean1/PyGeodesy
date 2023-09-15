
# -*- coding: utf-8 -*-

u'''I{Boolean} operations on I{composite} polygons and I{clip}s.

Classes L{BooleanFHP} and L{BooleanGH} are I{composites} and
provide I{boolean} operations C{intersection}, C{difference},
C{reverse-difference}, C{sum} and C{union}.

@note: A I{clip} is defined as a single, usually closed polygon,
       a I{composite} is a collection of one or more I{clip}s.

@see: U{Forster-Hormann-Popa<https://www.ScienceDirect.com/science/
      article/pii/S259014861930007X>} and U{Greiner-Hormann
      <http://www.Inf.USI.CH/hormann/papers/Greiner.1998.ECO.pdf>}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import isodd, isscalar, issubclassof, map2
from pygeodesy.constants import EPS, EPS2, INT0, _0_0, _0_5, _1_0
from pygeodesy.errors import ClipError, _IsnotError, _TypeError, \
                            _ValueError, _xattr, _xkwds_get
from pygeodesy.fmath import favg, hypot, hypot2
# from pygeodesy.fsums import fsum1  # _MODS
from pygeodesy.interns import NN, _BANG_, _clip_, _clipid_, _COMMASPACE_, \
                             _composite_, _DOT_, _e_, _ELLIPSIS_, _few_, \
                             _height_, _lat_,_LatLon_, _lon_, _not_, \
                             _points_, _scalar_,_SPACE_, _too_, _X_, _x_, \
                             _B_, _d_, _R_  # PYCHOK used!
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.latlonBase import LatLonBase, \
                                 LatLon2Tuple, Property_RO, property_RO
from pygeodesy.named import Fmt, _Named, _NotImplemented, pairs, unstr
# from pygeodesy.namedTuples import LatLon2Tupe  # from .latlonBase
# from pygeodesy.points import boundsOf  # _MODS
# from pygeodesy.props import Property_RO, property_RO  # from .latlonBase
# from pygeodesy.streprs import Fmt, pairs, unstr  # from .named
from pygeodesy.units import Height, HeightX
from pygeodesy.utily import fabs, _unrollon, _Wrap

# from math import fabs  # from .utily

__all__ = _ALL_LAZY.booleans
__version__ = '23.09.11'

_0_EPS =  EPS  # near-zero, positive
_EPS_0 = -EPS  # near-zero, negative
_1_EPS = _1_0 + EPS   # near-one, over
_EPS_1 = _1_0 - EPS   # near-one, under
_10EPS =  EPS * 10    # see ._2Abs, ._10eps

_alpha_     = 'alpha'
_boolean_   = 'boolean'
_case_      = 'case'
_corners_   = 'corners'
_duplicate_ = 'duplicate'
_open_      = 'open'


def _Enum(txt, enum):  # PYCHOK unused
    return txt  # NN(txt, _TILDE_, enum)


class _L(object):  # Intersection labels
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


class _RP(object):  # RelativePositions
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


class _LatLonBool(_Named):
    '''(INTERNAL) Base class for L{LatLonFHP} and L{LatLonGH}.
    '''
    _alpha   = None       # point AND intersection else length
    _checked = False      # checked in phase 3 iff intersection
    _clipid  = INT0       # (polygonal) clip identifier, number
    _dupof   = None       # original of a duplicate
#   _e_x_str = NN         # shut up PyChecker
    _height  = Height(0)  # interpolated height, usually meter
    _linked  = None       # link to neighbor iff intersection
    _next    = None       # link to the next vertex
    _prev    = None       # link to the previous vertex

    def __init__(self, lat_ll, lon=None, height=0, clipid=INT0,
                                           wrap=False, name=NN):
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
           @kwarg wrap: If C{True}, wrap or I{normalize} B{C{lat}}
                        and B{C{lon}} (C{bool}).
           @kwarg name: Optional name (C{str}).
        '''
        if lon is None:
            y, x = lat_ll.lat, lat_ll.lon
            h = _xattr(lat_ll, height=height)
            c = _xattr(lat_ll, clipid=clipid)
        else:
            y, x = lat_ll, lon
            h, c = height, clipid
        self.y, self.x = _Wrap.latlon(y, x) if wrap else (y, x)
        # don't duplicate defaults
        if self._height != h:
            self._height = h
        if self._clipid != c:
            self._clipid = c
        if name:
            self.name = name

    def __abs__(self):
        return max(fabs(self.x), fabs(self.y))

    def __eq__(self, other):
        return other is self or bool(_other(self, other) and
                                      other.x == self.x  and
                                      other.y == self.y)

    def __ne__(self, other):  # required for Python 2
        return not self.__eq__(other)

    def __repr__(self):
        '''String C{repr} of this lat-/longitude.
        '''
        if self._prev or self._next:
            t = _ELLIPSIS_(self._prev, self._next)
            t = _SPACE_(self, Fmt.ANGLE(t))
        else:
            t =  str(self)
        return t

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
            t =  NN(t, self._e_x_str(k))  # PYCHOK expected
        return NN(self.name, t)

    def __sub__(self, other):
        _other(self, other)
        return self.__class__(self.y - other.y,  # classof
                              self.x - other.x)

    def _2A(self, p2, p3):
        # I{Signed} area of a triangle, I{doubled}.
        x, y = self.x, self.y
        return (p2.x - x) * (p3.y - y) - \
               (p3.x - x) * (p2.y - y)

    def _2Abs(self, p2, p3, eps=_10EPS):
        # I{Unsigned} area of a triangle, I{doubled}
        # or 0 if below the given threshold C{eps}.
        a = fabs(self._2A(p2, p3))
        return 0 if a < eps else a

    @property_RO
    def clipid(self):
        '''Get the I{clipid} (C{int} or C{0}).
        '''
        return self._clipid

    def _equi(self, llb, eps):
        # Is this LLB I{equivalent} to B{C{llb}} within
        # the given I{non-negative} tolerance B{C{eps}}?
        return not (fabs(llb.lon - self.x) > eps or
                    fabs(llb.lat - self.y) > eps)

    @property_RO
    def height(self):
        '''Get the I{height} (C{Height} or C{int}).
        '''
        h = self._height
        return HeightX(h) if self.isintersection else (
               Height(h)  if h else _LatLonBool._height)

    def isequalTo(self, other, eps=None):
        '''Is this point equal to an B{C{other}} within a given,
           I{non-negative} tolerance, ignoring C{height}?

           @arg other: The other point (C{LatLon}).
           @kwarg eps: Tolerance for equality (C{degrees} or C{None}).

           @return: C{True} if equivalent, C{False} otherwise (C{bool}).

           @raise TypeError: Invalid B{C{other}}.
        '''
        try:
            return self._equi(other, _eps0(eps))
        except (AttributeError, TypeError, ValueError):
            raise _IsnotError(_LatLon_, other=other)

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
    def lat(self):
        '''Get the latitude (C{scalar}).
        '''
        return self.y

    @property_RO
    def latlon(self):
        '''Get the lat- and longitude (L{LatLon2Tuple}).
        '''
        return LatLon2Tuple(self.y, self.x)

    def _link(self, other):
        # Make this and an other point neighbors.
        # assert _other(self, other)
        self._linked = other
        other._linked = self

    @property_RO
    def lon(self):
        '''Get the longitude (C{scalar}).
        '''
        return self.x

    def _toClas(self, Clas, clipid):
        # Return this vertex as a C{Clas} instance
        # (L{Clip[FHP|GH]4Tuple} or L{LatLon[FHP|GH]}).
        return Clas(self.lat, self.lon, self.height, clipid)


class LatLonFHP(_LatLonBool):
    '''A point or intersection in a L{BooleanFHP} clip or composite.
    '''
    _en_ex  = None
    _label  = None
    _2split = None  # or C{._Clip}
    _2xing  = False

    def __init__(self, lat_ll, *lon_h_clipid, **wrap_name):
        '''New C{LatLonFHP} from separate C{lat}, C{lon}, C{h}eight
           and C{clipid} scalars, or from a previous L{LatLonFHP},
           a L{ClipFHP4Tuple} or some other C{LatLon} instance.

           @arg lat_ll: Latitude (C{scalar}) or a lat/longitude
                        (L{LatLonFHP}, C{LatLon} or L{ClipFHP4Tuple}).
           @arg lon_h_clipid: Longitude (C{scalar}), C{h}eight and
                              C{clipid} iff B{C{lat_ll}} is scalar,
                              ignored otherwise.
           @kwarg wrap_name: Keyword arguments C{B{wrap}=False} and
                       C{B{name}=NN}.  If C{B{wrap} is True}, wrap
                       or I{normalize} the lat- and longitude
                       (C{bool}).  Optional B{C{name}} (C{str}).
        '''
        _LatLonBool.__init__(self, lat_ll, *lon_h_clipid, **wrap_name)

    def __add__(self, other):
        _other(self, other)
        return self.__class__(self.y + other.y, self.x + other.x)

    def __mod__(self, other):  # cross product
        _other(self, other)
        return self.x * other.y - self.y * other.x

    def __mul__(self, other):  # dot product
        _other(self, other)
        return self.x * other.x + self.y * other.y

    def __rmul__(self, other):  # scalar product
        if not isscalar(other):
            raise _IsnotError(_scalar_, other=other)
        return self.__class__(self.y * other, self.x * other)

    def _e_x_str(self, t):  # PYCHOK no cover
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

    def isenclosedBy(self, *composites_points, **wrap):
        '''Is this point inside one or more composites or polygons based
           the U{winding number<https://www.ScienceDirect.com/science/
           article/pii/S0925772101000128>}?

           @arg composites_points: Composites and/or iterables of points
                           (L{ClipFHP4Tuple}, L{ClipGH4Tuple}, L{LatLonFHP},
                           L{LatLonGH} or any C{LatLon}).
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        C{points} (C{bool}).

           @raise ValueError: Some C{points} invalid.

           @see: U{Algorithm 6<https://www.ScienceDirect.com/science/
                 article/pii/S0925772101000128>}.
        '''
        class _Pseudo(object):
            # Pseudo-_CompositeBase._clips tuple

            @property_RO
            def _clips(self):
                for cp in _Cps(_CompositeFHP, composites_points,
                                   LatLonFHP.isenclosedBy):  # PYCHOK yield
                    for c in cp._clips:
                        yield c

        return self._isinside(_Pseudo(), **wrap)

    def _isinside(self, composite, *excludes, **wrap):
        # Is this point inside a composite, excluding
        # certain C{_Clip}s?  I{winding number}?
        x, y, i = self.x, self.y, False
        for c in composite._clips:
            if c not in excludes:
                w = 0
                for p1, p2 in c._edges2(**wrap):
                    # edge [p1,p2] must straddle y
                    if (p1.y < y) is not (p2.y < y):  # or ^
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

    @property_RO
    def _prev_next2(self):
        # Adjust 2-tuple (._prev, ._next) iff a I{duplicate} intersection
        p, n = self, self._next
        if self._isduplicate:
            p = self._dupof
            while p._isduplicate:
                p = p._dupof
            while n._isduplicate:
                n = n._next
        return p._prev, n

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


class LatLonGH(_LatLonBool):
    '''A point or intersection in a L{BooleanGH} clip or composite.
    '''
    _entry = None   # entry or exit iff intersection

    def __init__(self, lat_ll, *lon_h_clipid, **wrap_name):
        '''New C{LatLonGH} from separate C{lat}, C{lon}, C{h}eight
           and C{clipid} scalars, or from a previous L{LatLonGH},
           L{ClipGH4Tuple} or some other C{LatLon} instance.

           @arg lat_ll: Latitude (C{scalar}) or a lat/longitude
                        (L{LatLonGH}, C{LatLon} or L{ClipGH4Tuple}).
           @arg lon_h_clipid: Longitude (C{scalar}), C{h}eight and
                              C{clipid} iff B{C{lat_ll}} is scalar,
                              ignored otherwise.
           @kwarg wrap_name: Keyword arguments C{B{wrap}=False} and
                       C{B{name}=NN}.  If C{B{wrap} is True}, wrap
                       or I{normalize} the lat- and longitude
                       (C{bool}).  Optional B{C{name}} (C{str}).
        '''
        _LatLonBool.__init__(self, lat_ll, *lon_h_clipid, **wrap_name)

    def _check(self):
        # Check-mark this vertex and its link.
        self._checked = True
        k = self._linked
        if k and not k._checked:
            k._checked = True

    def _e_x_str(self, t):  # PYCHOK no cover
        return t  if self._entry is None else NN(t,
             (_e_ if self._entry else _x_))

    def isenclosedBy(self, *composites_points, **wrap):
        '''Is this point inside one or more composites or polygons based
           on the U{even-odd-rule<https://www.ScienceDirect.com/science/
           article/pii/S0925772101000128>}?

           @arg composites_points: Composites and/or iterables of points
                           (L{ClipFHP4Tuple}, L{ClipGH4Tuple}, L{LatLonFHP},
                           L{LatLonGH} or any C{LatLon}).
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        C{points} (C{bool}).

           @raise ValueError: Some B{C{points}} invalid.
        '''
        class _Pseudo(object):
            # Pseudo-_CompositeBase._edges3 method

            def _edges3(self, **kwds):
                for cp in _Cps(_CompositeGH, composites_points,
                                   LatLonGH.isenclosedBy):  # PYCHOK yield
                    for e in cp._edges3(**kwds):
                        yield e

        return self._isinside(_Pseudo(), **wrap)

    def _isinside(self, composite, *bottom_top, **wrap):
        # Is this vertex inside the composite? I{even-odd rule}?

        def _x(y, p1, p2):
            # return C{x} at given C{y} on edge [p1,p2]
            return (y - p1.y) / (p2.y - p1.y) * (p2.x - p1.x)

        # The I{even-odd} rule counts the number of edges
        # intersecting a ray emitted from this point to
        # east-bound infinity.  When I{odd} this point lies
        # inside, if I{even} outside.
        y, i = self.y, False
        if not (bottom_top and _outside(y, y, *bottom_top)):
            x = self.x
            for p1, p2, _ in composite._edges3(**wrap):
                if (p1.y < y) is not (p2.y < y):  # or ^
                    r = p2.x > x
                    if p1.x < x:
                        b = r and (_x(y, p1, p2) > x)
                    else:
                        b = r or  (_x(y, p1, p2) > x)
                    if b:
                        i = not i
        return i


class _Clip(_Named):
    '''(INTERNAL) A I{doubly-linked} list representing a I{closed}
       polygon of L{LatLonFHP} or L{LatLonGH} points, duplicates
       and intersections with other clips.
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

    def __init__(self, composite, clipid=INT0):
        '''(INTERNAL) New C{_Clip}.
        '''
        # assert isinstance(composite, _CompositeBase)
        if clipid in composite._clipids:
            raise ClipError(clipid=clipid, txt=_duplicate_)
        self._composite  = composite
        self._id         = clipid
        self._LL         = composite._LL
        composite._clips = composite._clips + (self,)

    def __contains__(self, point):  # PYCHOK no cover
        '''Is the B{C{point}} in this clip?
        '''
        for v in self:
            if v is point:  # or ==?
                return True
        return False

    def __eq__(self, other):
        '''Is this clip I{equivalent} to an B{C{other}} clip,
           do both have the same C{len}, the same points, in
           the same order, possibly rotated?
        '''
        return self._equi(_other(self, other), 0)

    def __ge__(self, other):
        '''See method C{__lt__}.
        '''
        return not self.__lt__(other)

    def __gt__(self, other):
        '''Is this clip C{"above"} an B{C{other}} clip,
           located or stretched farther North or East?
        '''
        return self._bltr4 > _other(self, other)._bltr4

    def __hash__(self):  # PYCHOK no over
        return hash(self._bltr4)

    def __iter__(self):
        '''Yield the points, duplicates and intersections.
        '''
        v = f = self._first
        while v:
            yield v
            v = v._next
            if v is f:
                break

    def __le__(self, other):
        '''See method C{__gt__}.
        '''
        return not self.__gt__(other)

    def __len__(self):
        '''Return the number of points, duplicates and
           intersections in this clip.
        '''
        return self._len

    def __lt__(self, other):
        '''Is this clip C{"below"} an B{C{other}} clip,
           located or stretched farther South or West?
        '''
        return self._bltr4 < _other(self, other)._bltr4

    def __ne__(self, other):  # required for Python 2
        '''See method C{__eq__}.
        '''
        return not self.__eq__(other)

    _all = __iter__

    @property_RO
    def _all_ON_ON(self):
        # Check whether all vertices are ON_ON.
        L_ON_ON = _L.ON_ON
        return all(v._label is L_ON_ON for v in self)

    def _append(self, y_v, *x_h_clipid):
        # Append a point given as C{y}, C{x}, C{h}eight and
        # C{clipid} args or as a C{LatLon[FHP|GH]}.
        self._last = v = self._LL(y_v, *x_h_clipid) if x_h_clipid else y_v
        self._len += 1
        # assert v._clipid == self._id

        v._next = n = self._first
        if n is None:  # set ._first
            self._first = p = n = v
        else:  # insert before ._first
            v._prev = p = n._prev
        p._next = n._prev = v
        return v

#   def _appendedup(self, v, clipid=0):
#       # Like C{._append}, but only append C{v} if not a
#       # duplicate of the one previously append[edup]'ed.
#       y, x, p = v.y, v.x, self._last
#       if p is None or y != p.y or x != p.x or clipid != p._clipid:
#           p = self._append(y, x, v._height, clipid)
#           if v._linked:
#               p._linked = True  # to force errors
#       return p

    @Property_RO
    def _bltr4(self):
        # Get the bounds as 4-tuple C{(bottom, left, top, right)}.
        return map2(float, _MODS.points.boundsOf(self, wrap=False))

    def _bltr4eps(self, eps):
        # Get the ._bltr4 bounds tuple, oversized.
        if eps > 0:  # > EPS
            yb, xl, yt, xr = self._bltr4
            yb, yt = _low_high_eps2(yb, yt, eps)
            xl, xr = _low_high_eps2(xl, xr, eps)
            t = yb, xl, yt, xr
        else:
            t = self._bltr4
        return t

    def _closed(self, raiser):  # PYCHOK unused
        # End a clip, un-close it and check C{len}.
        p, f = self._last, self._first
        if f and f._prev is p and p is not f and \
                 p._next is f and p == f:  # PYCHOK no cover
            # un-close the clip
            f._prev = p = p._prev
            p._next = f
            self._len -= 1
#       elif f and raiser:
#           raise self._OpenClipError(p, f)
        if len(self) < 3:
            raise self._Error(_too_(_few_))

    def _dup(self, q):
        # Duplicate a point (or intersection) as intersection.
        v = self._insert(q.y, q.x, q)
        v._alpha = q._alpha or _0_0  # _0_0 replaces None
        v._dupof = q._dupof or q
        # assert v._prev is q
        # assert q._next is v
        return v

    def _edges2(self, wrap=False, **unused):
        # Yield each I{original} edge as a 2-tuple
        # (p1, p2), a pair of C{LatLon[FHP|GH])}s.
        p1 = p = f = self._first
        while p:
            p2 = p = p._next
            if p.ispoint:
                if wrap and p is not f:
                    p2 = _unrollon(p1, p)
                yield p1, p2
                p1 = p2
            if p is f:
                break

    def _equi(self, clip, eps):
        # Is this clip I{equivalent} to B{C{clip}} within
        # the given I{non-negative} tolerance B{C{eps}}?
        r, f = len(self), self._first
        if f and r == len(clip) and self._bltr4eps(eps) \
                                 == clip._bltr4eps(eps):
            _equi = _LatLonBool._equi
            for v in clip:
                if _equi(f, v, eps):
                    s, n = f, v
                    for _ in range(r):
                        s, n = s._next, n._next
                        if not _equi(s, n, eps):
                            break  # next v
                    else:  # equivalent
                        return True
        return False

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
        return ClipError(cp._kind, cp.name, **kwds)

    def _index(self, clips, eps):
        # see _CompositeBase._equi
        for i, c in enumerate(clips):
            if c._equi(self, eps):
                return i
        raise ValueError(NN)  # like clips.index(self)

    def _insert(self, y, x, start, *end_alpha):
        # insertVertex between points C{start} and
        # C{end}, ordered by C{alpha} iff given.
        v = self._LL(y, x, start._height, start._clipid)
        n = start._next
        if end_alpha:
            end, alpha = end_alpha
            v._alpha   = alpha
            v._height  = favg(v._height, end._height, f=alpha)
            # assert start is not end
            while n is not end and n._alpha < alpha:
                n = n._next
        v._next = n
        v._prev = p = n._prev
        p._next = n._prev = v
        self._len += 1
#       _Clip._bltr4._update(self)
#       _Clip._ishole._update(self)
        return v

    def _intersection(self, unused, q, *p1_p2_alpha):
        # insert an intersection or make a point both
        if p1_p2_alpha:  # intersection on edge
            v = self._insert(q.y, q.x, *p1_p2_alpha)
        else:  # intersection at point
            v = q
            # assert not v._linked
            # assert v._alpha is None
        return v

    def _intersections(self):
        # Yield all intersections, some may be points too.
        for v in self:
            if v.isintersection:
                yield v

    @Property_RO
    def _ishole(self):  # PYCHOK no cover
        # Is this clip a hole inside its composite?
        v = self._first
        return v._isinside(self._composite, self) if v else False

    @property_RO
    def _nodups(self):
        # Yield all non-duplicates.
        for v in self:
            if not v._dupof:
                yield v

    def _noXings(self, Union):
        # Are all intersections non-CROSSINGs, -BOUNCINGs?
        Ls = _L.BOUNCINGs if Union else _L.CROSSINGs
        return all(v._label not in Ls for v in self._intersections())

    def _OpenClipError(self, s, e):  # PYCHOK no cover
        # Return a C{CloseError} instance
        t = NN(s, _ELLIPSIS_(_COMMASPACE_, e))
        return self._Error(_SPACE_(_open_, t))

    def _point2(self, insert):
        # getNonIntersectionPoint and -Vertex
        if not (insert and self._noInters):
            for p in self._points(may_be=False):  # not p._isduplicated?
                return p, None
            for n in self._intersections():
                p, _ = n._prev_next2
                k = p._linked
                if k:
                    if n._linked not in k._prev_next2:
                        # create a pseudo-point
                        k = _0_5 * (p + n)
                        if insert:
                            k = self._insert(k.y, k.x, n._prev)
                            r = k  # to remove later
                        else:  # no ._prev, ._next
                            k._clipid = n._clipid
                            r = None
                        return k, r
        return None, None

    def _points(self, may_be=True):
        # Yield all points I{in original order}, which may be intersections too.
        for v in self:
            if v.ispoint and (may_be or not v.isintersection):
                yield v

    def _remove2(self, v):
        # Remove vertex C{v}.
        # assert not v._isduplicated
        if len(self) > 1:
            p = v._prev
            p._next = n = v._next
            n._prev = p
            if self._first is v:
                self._first = n
            if self._last is v:
                self._last = p
            self._len -= 1
        else:
            n = self._last = \
            p = self._first = None
            self._len = 0
        return p, n

    def _update_all(self):  # PYCHOK no cover
        # Zap the I{cached} properties.
        _Clip._bltr4._update( self)
        _Clip._ishole._update(self)
        return self  # for _special_identicals

    def _Xings(self):
        # Yield all I{un-checked} CROSSING intersections.
        CROSSING = _L.CROSSING
        for v in self._intersections():
            if v._label is CROSSING and not v._checked:
                yield v


class _CompositeBase(_Named):
    '''(INTERNAL) Base class for L{BooleanFHP} and L{BooleanGH}
       (C{_CompositeFHP} and C{_CompositeGH}).
    '''
    _clips  =  ()   # tuple of C{_Clips}
    _eps    =  EPS  # null edges
    _kind   = _corners_
    _LL     = _LatLonBool  # shut up PyChecker
    _raiser =  False
    _xtend  =  False

    def __init__(self, lls, name=NN, kind=NN, eps=EPS):
        '''(INTERNAL) See L{BooleanFHP} and L{BooleanGH}.
        '''
        n = name or _xattr(lls, name=NN)
        if n:
            self.name = n
        if kind:
            self._kind = kind
        if self._eps < eps:
            self._eps = eps

        c = _Clip(self)
        lp = None
        for ll in lls:
            ll = self._LL(ll)
            if lp is None:
                c._id = ll._clipid  # keep clipid
                lp = c._append(ll)
            elif ll._clipid != lp._clipid:  # new clip
                c._closed(self.raiser)
                c = _Clip(self, ll._clipid)
                lp = c._append(ll)
            elif abs(ll - lp) > eps:  # PYCHOK lp
                lp = c._append(ll)
            else:
                c._dups += 1
        c._closed(self.raiser)

    def __contains__(self, point):  # PYCHOK no cover
        '''Is the B{C{point}} in one of the clips?
        '''
        for c in self._clips:
            if point in c:
                return True
        return False

    def __eq__(self, other):
        '''Is this I{composite} equivalent to an B{C{other}}, i.e.
           do both contain I{equivalent} clips in the same or in a
           different order?  Two clips are considered I{equivalent}
           if both have the same points etc. in the same order,
           possibly rotated.
        '''
        return self._equi(_other(self, other), 0)

    def __iter__(self):
        '''Yield all points, duplicates and intersections.
        '''
        for c in self._clips:
            for v in c:
                yield v

    def __ne__(self, other):  # required for Python 2
        '''See method C{__eq__}.
        '''
        return not self.__eq__(other)

    def __len__(self):
        '''Return the I{total} number of points.
        '''
        return sum(map(len, self._clips)) if self._clips else 0

    def __repr__(self):
        '''String C{repr} of this composite.
        '''
        c = len(self._clips)
        c = Fmt.SQUARE(c) if c > 1 else NN
        n = Fmt.SQUARE(len(self))
        t = Fmt.PAREN(self)  # XXX not unstr
        return NN(self.__class__.__name__, c, n, t)

    def __str__(self):
        '''String C{str} of this composite.
        '''
        return _COMMASPACE_.join(map(str, self))

    @property_RO
    def _bottom_top_eps2(self):
        # Get the bottom and top C{y} bounds, oversized.
        return _min_max_eps2(min(v.y for v in self),
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

#   def _clipidups(self, other):
#       # Number common C{clipid}s between this and an C{other} composite
#       return len(set(self._clipids).intersection(set(other._clipids)))

    def _edges3(self, **raiser_wrap):
        # Yield each I{original} edge as a 3-tuple
        # C{(LatLon[FHP|GH], LatLon[FHP|GH], _Clip)}.
        for c in self._clips:
            for p1, p2 in c._edges2(**raiser_wrap):
                yield p1, p2, c

    def _encloses(self, lat, lon, **wrap):
        # see function .points.isenclosedBy
        return self._LL(lat, lon).isenclosedBy(self, **wrap)

    @property
    def eps(self):
        '''Get the null edges tolerance (C{degrees}, usually).
        '''
        return self._eps

    @eps.setter  # PYCHOK setter!
    def eps(self, eps):
        '''Set the null edges tolerance (C{degrees}, usually).
        '''
        self._eps = eps

    def _10eps(self, **eps):
        # Get eps for _LatLonBool._2Abs
        e = _xkwds_get(eps, eps=self._eps)
        if e != EPS:
            e *= _10EPS / EPS
        else:
            e  = _10EPS
        return e

    def _equi(self, other, eps):
        # Is this composite I{equivalent} to an B{C{other}} within
        # the given, I{non-negative} tolerance B{C{eps}}?
        cs, co = self._clips, other._clips
        if cs and len(cs) == len(co):
            if eps > 0:
                _index = _Clip._index
            else:
                def _index(c, cs, unused):
                    return cs.index(c)
            try:
                cs = list(sorted(cs))
                for c in sorted(co):
                    cs.pop(_index(c, cs, eps))
            except ValueError:  # from ._index
                pass
            return False if cs else True
        else:  # both null?
            return False if cs or co else True

    def _intersections(self):
        # Yield all intersections.
        for c in self._clips:
            for v in c._intersections():
                yield v

    def isequalTo(self, other, eps=None):
        '''Is this boolean/composite equal to an B{C{other}} within
           a given, I{non-negative} tolerance?

           @arg other: The other boolean/composite (C{Boolean[FHP|GB]}).
           @kwarg eps: Tolerance for equality (C{degrees} or C{None}).

           @return: C{True} if equivalent, C{False} otherwise (C{bool}).

           @raise TypeError: Invalid B{C{other}}.

           @see: Method C{__eq__}.
        '''
        if isinstance(other, _CompositeBase):
            return self._equi(other, _eps0(eps))
        raise _IsnotError(_boolean_, _composite_, other=other)

    def _kwds(self, op, **more):
        # Get all keyword arguments as C{dict}.
        kwds = dict(raiser=self.raiser, eps=self.eps,
                      name=self.name or op.__name__)
        kwds.update(more)
        return kwds

    @property_RO
    def _left_right_eps2(self):
        # Get the left and right C{x} bounds, oversized.
        return _min_max_eps2(min(v.x for v in self),
                             max(v.x for v in self))

    def _points(self, may_be=True):  # PYCHOK no cover
        # Yield all I{original} points, which may be intersections too.
        for c in self._clips:
            for v in c._points(may_be=may_be):
                yield v

    @property
    def raiser(self):
        '''Get the option to throw L{ClipError} exceptions (C{bool}).
        '''
        return self._raiser

    @raiser.setter  # PYCHOK setter!
    def raiser(self, throw):
        '''Set the option to throw L{ClipError} exceptions (C{bool}).
        '''
        self._raiser = bool(throw)

    def _results(self, _presults, Clas, closed=False, inull=False, **eps):
        # Yield the dedup'd results, as L{ClipFHP4Tuple}s
        C = self._LL if Clas is None else Clas
        e = self._10eps(**eps)
        for clipid, ns in enumerate(_presults):
            f = p = v = None
            for n in ns:
                if f is None:
                    yield n._toClas(C, clipid)
                    f = p = n
                elif v is None:
                    v = n  # got f, p, v
                elif inull or p._2Abs(v, n, eps=e):
                    yield v._toClas(C, clipid)
                    p, v = v, n
                else:  # null, colinear, ... skipped
                    v = n
            if v and (inull or p._2Abs(v, f, eps=e)):
                yield v._toClas(C, clipid)
                p = v
            if f and p != f and closed:  # close clip
                yield f._toClas(C, clipid)

    def _sum(self, other, op):
        # Combine this and an C{other} composite
        LL = self._LL
        sp = self.copy(name=self.name or op.__name__)
        sp._clips, sid = (), INT0  # new clips
        for cp in (self, other):
            for c in cp._clips:
                _ap = _Clip(sp, sid)._append
                for v in c._nodups:
                    _ap(LL(v.y, v.x, v.height, sid))
                sid += 1
        return sp

    def _sum1(self, _a_p, *args, **kwds):  # in .karney, .points
        # Sum the area or perimeter of all clips
        return _MODS.fsums.fsum1((_a_p(c, *args, **kwds) for c in self._clips), floats=True)

    def _sum2(self, LL, _a_p, *args, **kwds):  # in .sphericalNvector, -Trigonometry
        # Sum the area or perimeter of all clips

        def _lls(clip):  # convert clip to LLs
            _LL = LL
            for v in clip:
                yield _LL(v.lat, v.lon)  # datum=Sphere

        return _MODS.fsums.fsum1((_a_p(_lls(c), *args, **kwds) for c in self._clips), floats=True)

    def toLatLon(self, LatLon=None, closed=False, **LatLon_kwds):
        '''Yield all (non-duplicate) points and intersections
           as an instance of B{C{LatLon}}.

           @kwarg LatLon: Class to use (C{LatLon}) or if C{None},
                          L{LatLonFHP} or L{LatLonGH}.
           @kwarg closed: If C{True}, close each clip (C{bool}).
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}}
                               keyword arguments, ignore if
                               C{B{LatLon} is None}.

           @raise TypeError: Invalid B{C{LatLon}}.

           @note: For intersections, C{height} is an instance
                  of L{HeightX}, otherwise of L{Height}.
        '''
        if LatLon is None:
            LL, kwds = self._LL, {}
        elif issubclassof(LatLon, _LatLonBool, LatLonBase):
            LL, kwds = LatLon, LatLon_kwds
        else:
            raise _TypeError(LatLon=LatLon)

        for c in self._clips:
            lf, cid = None, c._id
            for v in c._nodups:
                ll = LL(v.y, v.x, **kwds)
                ll._height = v.height
                if ll._clipid != cid:
                    ll._clipid = cid
                yield ll
                if lf is None:
                    lf = ll
            if closed and lf:
                yield lf


class _CompositeFHP(_CompositeBase):
    '''(INTERNAL) A list of clips representing a I{composite}
       of L{LatLonFHP} points, duplicates and intersections
       with an other I{composite}.
    '''
    _LL    = LatLonFHP
    _Union = False

    def __init__(self, lls, raiser=False, **name_kind_eps):
        # New L{_CompositeFHP}.
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
                    n._label = None  # n.__dict__.pop('_label')
                    n = n._next
                    if n is v or n._label is not L.ON_ON:  # n._label and ...
                        break
                a = L.LEFT_ON if n._label is L.ON_LEFT else L.RIGHT_ON
                v._label = n._label = L.BOUNCING_D if a is b else L.CROSSING_D

        # 3) Copy labels
        for v in self._intersections():
            v._linked._label = v._label

    def _clip(self, corners, Union=False, Clas=None,
                           **closed_inull_raiser_eps):
        # Clip this composite with another one, C{corners},
        # using Foster-Hormann-Popa's algorithm.
        P = self
        Q = self._class(corners, closed_inull_raiser_eps,
                                 eps=P._eps, raiser=False)
        if Union:
            P._Union = Q._Union = True

        bt = Q._bottom_top_eps2
        lr = Q._left_right_eps2
        # compute and insert intersections
        for p1, p2, Pc in P._edges3(**closed_inull_raiser_eps):
            if not (_outside(p1.x, p2.x, *lr) or
                    _outside(p1.y, p2.y, *bt)):
                e = _EdgeFHP(p1, p2)
                if e._dp2 > EPS2:  # non-null edge
                    for q1, q2, Qc in Q._edges3(**closed_inull_raiser_eps):
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
        return P._results(P._presults(Q), Clas, **closed_inull_raiser_eps)

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
        # Yield the result clips, each as a generator
        # of the L{_LatLonFHP}s in that clip
        for cp in (self, other):
            for c in cp._clips:
                if c._pushback:
                    yield c._all()
        for c in self._clips:
            for X in c._Xings():
                yield self._resultX(X)

    def _resultX(self, X):
        # Yield the results from CROSSING C{X}.
        L, U, v = _L, self._Union, X
        while v:
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
            if v is X:
                break

    def _set_entry_exits(self, other):  # MCCABE 14
        # 4) Set entry/exit flags
        L, U = _L, self._Union
        for c in self._clips:
            n, k = c._point2(True)
            if n:
                f = n
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
                            if t is L.EXIT:  # first_chain_vertex == True
                                s = L.Toggle[s]
                        elif b is L.BOUNCING_D:
                            n._en_ex = s
                            if (s is t) ^ U:
                                n._2xing = True  # see ._splits_xings
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
        _u = _Clip._update_all
        cds = dict((c._id, _u(c)) for c in other._identicals)
        # assert len(cds) == len(other._identicals)
        if cds:  # PYCHOK no cover
            for c in self._identicals:
                c._update_all()
                for v in c._intersections():
                    d = cds.get(v._linked._clipid, None)
                    if d and d._ishole is c._ishole:
                        c._pushback = True
                        break  # next c

    @property_RO
    def _2splits(self):
        # Yield all intersections marked C{._2split}
        for p in self._intersections():
            if p._2split:
                # assert isinstance(p._2split, _Clip)
                yield p

    def _splits_xings(self, other):  # MCCABE 15
        # 5) Handle split pairs and 6) crossing candidates

        def _2A_dup2(p, P):  # PYCHOK unused
            p1, p2 = p._prev_next2
            ap = p1._2A(p, p2)
            Pc = p._2split
            # assert Pc in P._clips
            # assert p in Pc
            return ap, Pc._dup(p)

        def _links2(ps, qs):  # PYCHOK P unused?
            # Yield each link as a 2-tuple(p, q)
            id_qs = set(map(id, qs))
            if id_qs:
                for p in ps:
                    q = p._linked
                    if q and id(q) in id_qs:
                        yield p, q

        L = _L
        E =  L.ENTRY if self._Union else L.EXIT
        X =  L.Toggle[E]
        for p, q in _links2(self._2splits, other._2splits):
            ap, pp = _2A_dup2(p, self)
            aq, qq = _2A_dup2(q, other)
            if (ap * aq) > 0:  # PYCHOK no cover
                p._link(qq)  # overwrites ...
                q._link(pp)  # ... p-q link
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
    '''(INTERNAL) A list of clips representing a I{composite}
       of L{LatLonGH} points, duplicates and intersections
       with an other I{composite}.
    '''
    _LL    = LatLonGH
    _xtend = False

    def __init__(self, lls, raiser=False, xtend=False, **name_kind_eps):
        # New L{_CompositeGH}.
        if xtend:
            self._xtend = True
        elif raiser:
            self._raiser = True
        _CompositeBase.__init__(self, lls, **name_kind_eps)

    def _clip(self, corners, s_entry, c_entry, Clas=None,
                           **closed_inull_raiser_xtend_eps):
        # Clip this polygon with another one, C{corners}.

        # Core of Greiner/Hormann's algorithm, enhanced U{Correia's
        # <https://GitHub.com/helderco/univ-polyclip>} implementation***
        # and extended to optionally handle so-called "degenerate cases"
        S = self
        C = self._class(corners, closed_inull_raiser_xtend_eps,
                                 raiser=False, xtend=False)
        bt = C._bottom_top_eps2
        lr = C._left_right_eps2
        # 1. find intersections
        for s1, s2, Sc in S._edges3(**closed_inull_raiser_xtend_eps):
            if not (_outside(s1.x, s2.x, *lr) or
                    _outside(s1.y, s2.y, *bt)):
                e = _EdgeGH(s1, s2, **closed_inull_raiser_xtend_eps)
                if e._hypot2 > EPS2:  # non-null edge
                    for c1, c2, Cc in C._edges3(**closed_inull_raiser_xtend_eps):
                        for y, x, sa, ca in e._intersect4(c1, c2):
                            s = Sc._insert(y, x, s1, s2, sa)
                            c = Cc._insert(y, x, c1, c2, ca)
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
        return S._results(S._presults(), Clas, **closed_inull_raiser_xtend_eps)

    @property_RO
    def _first(self):
        # Get the very first vertex of the first clip
        for v in self:
            return v
        return None  # PYCHOK no cover

    def _kwds(self, op, **more):
        # Get the kwds C{dict}.
        return _CompositeBase._kwds(self, op, xtend=self.xtend, **more)

    def _presults(self):
        # Yield the unchecked intersection(s).
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

    @property
    def xtend(self):
        '''Get the option to handle I{degenerate cases} (C{bool}).
        '''
        return self._xtend

    @xtend.setter  # PYCHOK setter!
    def xtend(self, xtend):
        '''Set the option to handle I{degenerate cases} (C{bool}).
        '''
        self._xtend = bool(xtend)


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

        self._lr, \
        self._bt  = _left_right_bottom_top_eps2(p1, p2)

    def _intersect3(self, q1, q2):
        # Yield intersection(s) Type or C{None}
        if not (_outside(q1.x, q2.x, *self._lr) or
                _outside(q1.y, q2.y, *self._bt)):
            dq  = q2 - q1
            dq2 = dq * dq  # dot product, hypot2
            if dq2 > EPS2:  # like ._clip
                T,  E  = None, _EdgeFHP  # self.__class__
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
        self._lr, \
        self._bt  = _left_right_bottom_top_eps2(s1, s2)

        if xtend:
            self._xtend = True
        elif raiser:
            self._raiser = True

    def _alpha2(self, x, y, dx, dy):
        # Return C{(alpha)}, see .points.nearestOn5
        a = (y * dy + x * dx) / self._hypot2
        d = (y * dx - x * dy) / self._hypot0
        return a, fabs(d)

    def _Error(self, n, *args, **kwds):  # PYCHOK no cover
        t =  unstr(_EdgeGH.__name__, self._s1, self._s2)
        t = _DOT_(t, _EdgeGH._intersect4.__name__)
        t =  unstr(t, *args, **kwds)
        return ClipError(_case_, n, txt=t)

    @Property_RO
    def _hypot0(self):
        _, sx, _, sy = self._x_sx_y_sy
        return hypot(sx, sy) * _0_EPS

    @Property_RO
    def _hypot2(self):
        _, sx, _, sy = self._x_sx_y_sy
        return hypot2(sx, sy)

    def _intersect4(self, c1, c2, parallel=True):  # MCCABE 14
        # Yield the intersection(s) of this and another edge.

        # @return: None, 1 or 2 intersections, each a 4-Tuple
        #          (y, x, s_alpha, c_alpha) with intersection
        #          coordinates x and y and both alphas.

        # @raise ClipError: Intersection unhandled.

        # @see: U{Intersection point of two line segments
        #       <http://PaulBourke.net/geometry/pointlineplane/>}.
        c1_x, c1_y = c1.x, c1.y
        if not (_outside(c1_x, c2.x, *self._lr) or
                _outside(c1_y, c2.y, *self._bt)):
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
                        yield (y + sa * sy), (x + sa * sx), sa, ca

                    # unhandled, "degenerate" cases 1, 2 or 3
                    elif self._raiser and not (sa < _EPS_0 or sa > _1_EPS):  # PYCHOK no cover
                        raise self._Error(1, c1, c2, sa=sa)  # intersection at s1 or s2

                elif self._raiser and not (ca < _EPS_0 or ca > _1_EPS):  # PYCHOK no cover
                    # intersection at c1 or c2 or at c1 or c2 and s1 or s2
                    sa = (cx * dy - cy * dx) / d
                    e = 2 if sa < _EPS_0 or sa > _1_EPS else 3
                    raise self._Error(e, c1, c2, ca=ca)

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
                            yield y, x, ca1, _alpha1(ca)
                        else:  # c1 is between s1 and s2
                            yield (y + sa1 * sy), (x + sa1 * sx), sa1, ca1
                        if sa2 > 1:  # s2 is between c1 and c2
                            ca *= sa2 - _1_0
                            yield (y + sy), (x + sx), ca2, _alpha1(ca2 - ca)
                        else:  # c2 is between s1 and s2
                            yield (y + sa2 * sy), (x + sa2 * sx), sa2, ca2
                    elif self._raiser and not _outside(sa1, sa2, _0_0, _1_EPS):
                        raise self._Error(4, c1, c2, d1=d1, d2=d2)


class _BooleanBase(object):
    # Shared C{Boolean[FHP|GH]} methods.

    def __add__(self, other):
        '''Sum: C{this + other} clips.
        '''
        return self._sum(_other(self, other), self.__add__)  # PYCHOK OK

    def __and__(self, other):
        '''Intersection: C{this & other}.
        '''
        return self._boolean(other, False, False, self.__and__)  # PYCHOK OK

    def __iadd__(self, other):
        '''In-place sum: C{this += other} clips.
        '''
        return self._inplace(self.__add__(other))

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

    def __radd__(self, other):
        '''Reverse sum: C{other + this} clips.
        '''
        return _other(self, other)._sum(self, self.__radd__)

    def __rand__(self, other):
        '''Reverse intersection: C{other & this}
        '''
        return _other(self, other).__and__(self)

    def __ror__(self, other):
        '''Reverse union: C{other | this}
        '''
        return _other(self, other).__or__(self)

    def _boolean4(self, other, op):
        # Set up a new C{Boolean[FHP|GH]}.
        C = self.__class__
        kwds = C._kwds(self, op)
        a =  C(self, **kwds)
        b = _other(self, other)
        return a, b, C, kwds

    def _inplace(self, r):
        # Replace this with a L{Boolean*} result.
        self._clips, r._clips = r._clips, None
#       if self._raiser != r._raiser:
#           self._raiser = r._raiser
#       if self._xtend != r._xtend:
#           self._xtend = r._xtend
#       if self._eps != r._eps:
#           self._eps = r._eps
        return self


class BooleanFHP(_CompositeFHP, _BooleanBase):
    '''I{Composite} class providing I{boolean} operations between two
       I{composites} using U{Forster-Hormann-Popa<https://www.ScienceDirect.com/
       science/article/pii/S259014861930007X>}'s C++ implementation, transcoded
       to pure Python.

       The supported operations between (composite) polygon A and B are:

        -  C = A & B  or  A &= B,  intersection of A and B

        -  C = A + B  or  A += B,  sum of A and B clips

        -  C = A | B  or  A |= B,  union of A and B

        -  A == B     or  A != B,  equivalent A and B clips

        -  A.isequalTo(B, eps), equivalent within tolerance

       @see: Methods C{__eq__} and C{isequalTo}, function L{clipFHP4}
             and class L{BooleanGH}.
    '''
    _kind = _boolean_

    def __init__(self, lls, raiser=False, eps=EPS, name=NN):
        '''New L{BooleanFHP} operand for I{boolean} operation.

           @arg lls: The polygon points and clips (iterable of L{LatLonFHP}s,
                      L{ClipFHP4Tuple}s or other C{LatLon}s).
           @kwarg raiser: If C{True}, throw L{ClipError} exceptions (C{bool}).
           @kwarg esp: Tolerance for eliminating null edges (C{degrees}, same
                       units as the B{C{lls}} coordinates).
           @kwarg name: Optional name (C{str}).
        '''
        _CompositeFHP.__init__(self, lls, raiser=raiser,
                                          eps=eps, name=name)

    def __isub__(self, other):
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __rsub__(self, other):
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __sub__(self, other):
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def _boolean(self, other, Union, unused, op):
        # One C{BooleanFHP} operation.
        p, q, C, kwds = self._boolean4(other, op)
        r = p._clip(q, Union=Union, **kwds)
        return C(r, **kwds)


class BooleanGH(_CompositeGH, _BooleanBase):
    '''I{Composite} class providing I{boolean} operations between two
       I{composites} using the U{Greiner-Hormann<http://www.Inf.USI.CH/
       hormann/papers/Greiner.1998.ECO.pdf>} algorithm and U{Correia
       <https://GitHub.com/helderco/univ-polyclip>}'s implementation,
       modified and extended.

       The supported operations between (composite) polygon A and B are:

        -  C = A - B  or  A -= B,  difference A less B

        -  C = B - A  or  B -= A,  difference B less B

        -  C = A & B  or  A &= B,  intersection of A and B

        -  C = A + B  or  A += B,  sum of A and B clips

        -  C = A | B  or  A |= B,  union of A and B

        -  A == B     or  A != B,  equivalent A and B clips

        -  A.isequalTo(B, eps), equivalent within tolerance

       @note: To handle I{degenerate cases} like C{point-edge} and
              C{point-point} intersections, use class L{BooleanFHP}.

       @see: Methods C{__eq__} and C{isequalTo}, function L{clipGH4}
             and class L{BooleanFHP}.
    '''
    _kind = _boolean_

    def __init__(self, lls, raiser=True, xtend=False, eps=EPS, name=NN):
        '''New L{BooleanFHP} operand for I{boolean} operation.

           @arg lls: The polygon points and clips (iterable of L{LatLonGH}s,
                      L{ClipGH4Tuple}s or other C{LatLon}s).
           @kwarg raiser: If C{True}, throw L{ClipError} exceptions (C{bool}).
           @kwarg xtend: If C{True}, extend edges of I{degenerate cases}, an
                         attempt to handle the latter (C{bool}).
           @kwarg esp: Tolerance for eliminating null edges (C{degrees}, same
                       units as the B{C{lls}} coordinates).
           @kwarg name: Optional name (C{str}).
        '''
        _CompositeGH.__init__(self, lls, raiser=raiser, xtend=xtend,
                                         eps=eps, name=name)

    def _boolean(self, other, s_entry, c_entry, op):
        # One C{BooleanGH} operation.
        s, c, C, kwds = self._boolean4(other, op)
        r = s._clip(c, s_entry, c_entry, **kwds)
        return C(r, **kwds)

    def __isub__(self, other):
        '''In-place difference: C{this -= other}.
        '''
        return self._inplace(self.__sub__(other))

    def __rsub__(self, other):
        ''' Reverse difference: C{other - this}
        '''
        return _other(self, other).__sub__(self)

    def __sub__(self, other):
        '''Difference: C{this - other}.
        '''
        return self._boolean(other, True, False, self.__sub__)


def _alpha1(alpha):
    # Return C{alpha} in C{[0..1]} range
    if _EPS_0 < alpha < _1_EPS:
        return max(_0_0, min(alpha, _1_0))
    t = _not_(Fmt.SQUARE(_ELLIPSIS_(0, 1)))
    raise ClipError(_alpha_, alpha, txt=t)


def _alpha4(a):
    # Return 4-tuple (alpha, -EPS < alpha < EPS,
    #                           0 < alpha < 1,
    #                       not 0 < alpha < 1)
    return (a, False, True,  False) if _0_EPS < a < _EPS_1 else (
           (a, False, False, True)  if _0_EPS < fabs(a) else
           (a, True,  False, False))


def _Cps(Cp, composites_points, where):
    # Yield composites and points as a C{Cp} composite.
    try:
        kwds = dict(kind=_points_, name=where.__name__)
        for cp in composites_points:
            yield cp if isBoolean(cp) else Cp(cp, **kwds)
    except (AttributeError, ClipError, TypeError, ValueError) as x:
        raise _ValueError(points=cp, cause=x)


def _eps0(eps):
    # Adjust C{eps} or C{None}.
    return eps if eps and eps > EPS else 0


def isBoolean(obj):
    '''Check for C{Boolean} composites.

       @arg obj: The object (any C{type}).

       @return: C{True} if B{C{obj}} is L{BooleanFHP},
                L{BooleanGH} oe some other composite,
                C{False} otherwise.
    '''
    return isinstance(obj, _CompositeBase)


def _left_right_bottom_top_eps2(p1, p2):
    '''(INTERNAL) Return 2-tuple C{(left, right), (bottom, top)}, oversized.
    '''
    return (_min_max_eps2(p1.x, p2.x),
            _min_max_eps2(p1.y, p2.y))


def _low_high_eps2(lo, hi, eps):
    '''(INTERNAL) Return 2-tuple C{(lo, hi)}, oversized.
    '''
    lo *= (_1_0 + eps) if lo < 0 else (_1_0 - eps)
    hi *= (_1_0 - eps) if hi < 0 else (_1_0 + eps)
    return (lo or -eps), (hi or eps)


def _min_max_eps2(*xs):
    '''(INTERNAL) Return 2-tuple C{(min, max)}, oversized.
    '''
    lo, hi = min(xs), max(xs)
    lo *= _1_EPS if lo < 0 else _EPS_1
    hi *= _EPS_1 if hi < 0 else _1_EPS
    return (lo or _EPS_0), (hi or _0_EPS)


def _other(this, other):
    '''(INTERNAL) Check for compatible C{type}s.
    '''
    C = this.__class__
    if isinstance(other, C):
        return other
    raise _IsnotError(C.__name__, other=other)


def _outside(x1, x2, lo, hi):
    '''(INTERNAL) Is C{(x1, x2)} outside C{(lo, hi)}?
    '''
    return max(x1, x2) < lo or min(x1, x2) > hi


__all__ += _ALL_DOCS(_BooleanBase, _Clip,
                     _CompositeBase, _CompositeFHP, _CompositeGH,
                     _LatLonBool)

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
