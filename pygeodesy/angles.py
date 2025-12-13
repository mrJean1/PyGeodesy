
# -*- coding: utf-8 -*-

u'''Classes L{Ang}, L{Deg}, L{Rad} and L{Lambertian} accurately representing an angle
as a 3-tuple C{(sine, cosine, turns)}, with C{turns} the number of full turns.

Transcoded to pure Python from I{Karney}'s GeographicLib 2.7 C++ class U{AngleT
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AngleT.html>}.

Copyright (C) U{Charles Karney <mailto:Karney@Alum.MIT.edu>} (2024-2025) and licensed
under the MIT/X11 License.  For more information, see the U{GeographicLib 2.7
<https://GeographicLib.SourceForge.io/>} documentation.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # noqa: E702 ;

from pygeodesy.basics import _copysign, map1, signBit, _signOf
from pygeodesy.constants import EPS, EPS0, NAN, PI2, _0_0, _N_0_0, \
                               _0_25, _1_0, _N_1_0, _4_0, _360_0, \
                               _copysign_0_0, _copysign_1_0, \
                               _flipsign, float_, _isfinite, \
                               _over, _pos_self, remainder
from pygeodesy.errors import _xkwds, _xkwds_get, _xkwds_pop2
from pygeodesy.fmath import hypot,  _ALL_LAZY, _MODS
# from pygeodesy.interns import _COMMASPACE_  # from .streprs
# from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS  # from .fmath
from pygeodesy.named import _Named, _NamedTuple, _Pass
from pygeodesy.props import Property_RO, property_doc_, property_RO, \
                           _allPropertiesOf_n, _update_all
from pygeodesy.streprs import Fmt, fstr, unstr,  _COMMASPACE_
from pygeodesy.units import Degrees, _isDegrees, _isRadians, Radians
from pygeodesy.utily import atan2, atan2d, sincos2, sincos2d, SinCos2

from math import asinh, ceil as _ceil, fabs, floor as _floor, \
                 isinf, isnan, sinh

__all__ = _ALL_LAZY.angles
__version__ = '25.12.02'

_EPS03 = EPS / (1 << 20)
# _HD = _180_0
# _QD =  _90_0
# _TD = _360_0
# _DM = _SM = _60_0
# _DS = _3600_0
_ZRND = _1_0 / 1024

_CARDINAL2 = {-2: (_N_0_0, _N_1_0),
              -1: (_N_1_0,   _0_0),
               1: (  _1_0,   _0_0),
               2: (  _0_0, _N_1_0)}.get


def _fint(f):
    # float of C{int(f)} preserving signed C{0}.
    i = int(f)
    return float_(i) if i else _copysign_0_0(f)


def _ncardinal(s, c, n):
    if n:
        n *= _4_0
    i = (1 if (-c) < fabs(s) else 2) if signBit(c) else \
        (1 if   c  < fabs(s) else 0)
    if i:
        n += _copysign(i, s)
    return n


def _normalize2(s, c):
    h = hypot(s, c)
    if _isfinite(h):
        sc = ((s / h), (c / h)) if h else (
             # If y is +/-0 and x = -0, +/-pi is returned,
             # or y is +/-0 and x = +0, +/-0 is returned,
             # so, retain the sign of s = +/-0
             _orthogonal2(False, s, c))
    elif isnan(h) or (isinf(s) and isinf(c)):
        sc =  NAN, NAN
    else:
        sc = _orthogonal2(isinf(s), s, c)
    return sc


def _other(x, unit=Radians, **unused):
    # get C{x} as C{Ang} from C{Degrees}, C{Radians} or C{Lambertian}
    return Ang.fromLambertian(x) if                       unit is Lambertian else (
           Ang.fromRadians(x)    if _isRadians(x, iscalar=unit is Radians)   else (
           Ang.fromDegrees(x)    if _isDegrees(x, iscalar=unit is Degrees)   else
          _raiseError(unit, x)))  # PYCHOK indent


def _orthogonal2(pred, s, c):
    return (_copysign_1_0(s), _copysign_0_0(c)) if pred else \
           (_copysign_0_0(s), _copysign_1_0(c))


def _raiseError(unit, arg, **kwds):
    raise TypeError(unstr(unit, arg, **kwds))


def _rnd(x):
    w = _ZRND - fabs(x)
    if w > 0:
        x = _copysign(_ZRND - w, x)
    return x


def _scnu4(s, c, n, unit=Radians, **unused):  # unit=Ang._unit
    s, c, n = map1(float, s, c, n)
    return _normalize2(s, c) + (n, unit)


class Ang(_Named):
    '''An accurate representation of angles, as 3-tuple C{(s, c, n)}.

       This class represents an angle via its sine C{s}, cosine C{c} and
       the number of full turns C{n}.  The angle is then C{atan2(s, c) +
       n * PI2}.  This representation offers several advantages:

         - cardinal directions (multiples of 90 degrees) are exactly represented
           (a benefit shared by representing angles as degrees)

         - angles very close to any cardinal direction are accurately represented

         - there's no loss of precision with large angles (outside the "normal"
           range [-180, +180])

         - various operations, such as adding a multiple of 90 degrees to an
           angle are performed exactly.

       @note: B{C{n}} is a C{float}, this allows it to be NAN, INF or NINF.
    '''
    _unit = Radians  # see _scnu4

    def __init__(self, s_ang=0, c=None, n=0, normal=True, **unit_name):
        '''New L{Ang}.

           @kwarg s_ang: A previous L{Ang}, C{Degrees}, C{Radians} if C{B{c}
                         is None}, otherwise the sine component (C{float}).
           @kwarg c: The cosine component (C{float}) iff C{not None}.
           @kwarg n: The number of L{PI2} turns (C{float}).
           @kwarg normal: If C{True}, B{C{s}} and B{C{c}} are normalized, i.e.
                          on the unit circle (C{boo}).
           @kwarg unit_name: Type C{B{unit}=}L{Radians} or L{Degrees} of scalar
                       scalar values (L{Degrees} or L{Radians}).

           @note: Either B{C{s}} or B{C{c}} can be INF or NINF, but not both.

           @note: By default, the point B{C{(s, c)}} is scaled to lie on the
                  unit circle.
        '''
        s, c, n, u = s_ang.scnu4 if isAng(s_ang) else (
                    _other(s_ang, **unit_name).scnu4 if c is None else
                    _scnu4(s_ang, c, n, **unit_name))
        if unit_name:
            u, name = _xkwds_pop2(unit_name, unit=u)
            if name:
                self.name = name
        self._n = _fint(n)
        self._s, self._c = (s, c) if normal else _normalize2(s, c)
        self.unit = u

    def __abs__(self):
        s, _ = self._float2()
        return self._float1(fabs(s))

    def __add__(self, other):
        return self.copy().__iadd__(other)

    def __bool__(self):  # PYCHOK Python 3+
        s, c, n = self.scn3
        return bool(s or c or n)

#   def __call__(self, *args, **kwds):  # PYCHOK no cover
#       return self._NotImplemented(*args, **kwds)

    def __ceil__(self):  # PYCHOK not special in Python 2-
        s, _ = self._float2()
        return self._float1(_ceil(s))

    def __cmp__(self, other):  # PYCHOK no cover
        s, r = self._float2(other)
        return _signOf(s, r)  # -1, 0, +1

    def __divmod__(self, other):
        s, r = self._float2(other)
        q, r = divmod(s, r)
        return q, self._float1(r)

    def __eq__(self, other):
        s, r = self._float2(other)
        return fabs(s - r) < EPS0

    def __float__(self):
        u = self.unit
        return self.radians    if u is Radians else (
               self.degrees    if u is Degrees else (
               self.lambertian if u is Lambertian else
              _raiseError(float, u)))  # PYCHOK indent

    def __floor__(self):  # PYCHOK not special in Python 2-
        s, _ = self._float2()
        return self._float1(_floor(s))

    def __floordiv__(self, other):
        return self.copy().__ifloordiv__(other)

#   def __format__(self, *other):  # PYCHOK no cover
#       return self._NotImplemented(self, *other)

    def __ge__(self, other):
        s, r = self._float2(other)
        return s >= r

    def __gt__(self, other):
        s, r = self._float2(other)
        return s > r

    def __hash__(self):  # PYCHOK no cover
        # @see: U{Notes for type implementors<https://docs.Python.org/
        #       3/library/numbers.html#numbers.Rational>}
        return hash(self.scn3)  # tuple.__hash__()

    def __iadd__(self, other):
        p = self._other(other)
        q = p.ncardinal + self.ncardinal
        s, c, n = self.scn3
        s, c = _normalize2(s * p.c + c * p.s,
                           c * p.c - s * p.s)
        q -= _ncardinal(s, c, n)
        n  = _fint(q * _0_25) + p.n
        if n:
            self._n += n
            self._s  = s
            self._c  = c
            _update_all(self)
        return self._update(s, c)

    def __ifloordiv__(self, other):
        s, r = self._float2(other)
        return self._ifloat(s // r)

    def __imatmul__(self, other):  # PYCHOK no cover
        return self._notImplemented()

    def __imod__(self, other):
        s, r = self._float2(other)
        return self._ifloat(s % r)

    def __imul__(self, other):
        s, r = self._float2(other)
        return self._ifloat(s * r)

    def __int__(self):
        s, _ = self._float2(0)
        return int(s)

    def __invert__(self):  # PYCHOK no cover
        # Luciano Ramalho, "Fluent Python", O'Reilly, 2nd Ed, 2022 p. 567
        return self._notImplemented()

    def __ipow__(self, other, *mod):  # PYCHOK 2 vs 3 args
        s, r = self._float2(other)
        return self._ifloat(pow(s, r, *mod))

    def __isub__(self, other):
        return self.__iadd__(-other)

#   def __iter__(self):
#       '''
#       return self._NotImplemented()

    def __itruediv__(self, other):
        s, r = self._float2(other)
        return self._ifloat(s / r)

    def __le__(self, other):
        s, r = self._float2(other)
        return s <= r

    def __lt__(self, other):
        s, r = self._float2(other)
        return s < r

    def __matmul__(self, other):  # PYCHOK no cover
        return self._notImplemented(other)

    def __mod__(self, other):
        s, r = self._float2(other)
        return self._float1(s % r)

    def __mul__(self, other):
        return self.copy().__imul__(other)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __neg__(self):
        s, c, n =  self.scn3
        s,    n = _flipsign(s), _flipsign(n)
        return self._Ang(s, c, n)  # normal=True

    def __pos__(self):
        return self if _pos_self else self.copy()

    def __pow__(self, other, *mod):  # PYCHOK 2 vs 3 args
        return self.copy().__ipow__(other, *mod)

    def __radd__(self, other):
        return self._other(other) + self

    def __rdivmod__(self, other):
        return divmod(self._other(other), self)

    def __repr__(self):
        return self.toRepr()

    def __rfloordiv__(self, other):
        return self._other(other) // self

    def __rmatmul__(self, other):  # PYCHOK no cover
        return self._notImplemented(self, other)

    def __rmod__(self, other):
        return self._other(other) % self

    def __rmul__(self, other):
        return self._other(other) * self

    def __round__(self, *ndigits):  # PYCHOK Python 3+
        return self.round(*ndigits)

    def __rpow__(self, other, *mod):
        return pow(self._other(other), self, *mod)

    def __rsub__(self, other):
        return self._other(other) - self

    def __rtruediv__(self, other):
        return self._other(other) / self

    def __str__(self):
        return self.toStr(0)  # ignore turns

    def __sub__(self, other):
        return self.copy().__isub__(other)

    def __truediv__(self, other):
        return self.copy().__itruediv__(other)

    __trunc__ = __int__

    if _MODS.sys_version_info2 < (3, 0):  # PYCHOK no cover
        # <https://docs.Python.org/2/library/operator.html#mapping-operators-to-functions>
        __div__     = __truediv__
        __idiv__    = __itruediv__
        __long__    = __int__
        __nonzero__ = __bool__
        __rdiv__    = __rtruediv__

    def _Ang(self, s, *cn, **normal_unit_name):
        # return an C{Ang} like C{self}
        return Ang(s, *cn, **self._kwds(normal_unit_name))

    def base(self, *center):
        '''Return this C{Angle}'s base, optionally centered.
        '''
        r = self.copy()
        if center:
            c  = self._other(center[0])
            b  = self - c
            b  = b.base()
            b += c
            r.n0 = b.n0
        else:
            r.n = 0
        return r

    @property_RO
    def c(self):
        '''Get the cosine of this C{Angle} (C{float}).
        '''
        return self._c

    @staticmethod
    def cardinal(q=0, **unit_name):
        '''A cardinal direction.

           @kwarg q: The number of I{quarter} turns (C{scalar}).

           @return: An C{Ang} equivalent to B{C{q}} quarter turns.

           @note: B{C{q}} is truncated to an integer and signed
                  C{0} is distinguished.  C{Ang.NAN} is returned
                  if B{C{q}} is not finite.
        '''
        if _isfinite(q):
            if q:
                q = _fint(q)
                i =  int(remainder(q, _4_0))  # i is in [-2, 2]
                n = _fint((q - i) * _0_25)
                s, c = _CARDINAL2(i, ((_0_0 if q else q), _1_0))
                t = s is not q
            else:
                s, c, n, t = _copysign_0_0(q), 1, 0, True
            r = Ang(s, c, n, normal=t, **unit_name)
        else:
            r = Ang.NAN(**unit_name)
        return r

    def copy(self, **unit_name):  # PYCHOK signature
        '''Return a copy of this C{Ang}.
        '''
        return self._Ang(self, **self._kwds(unit_name))

    @Property_RO
    def degrees(self):
        '''Get this C{Ang} in C{degrees}.
        '''
        d = self.degrees0
        if self.n:
            d += self.n * _360_0
        return d  # XXX Degrees(d, self.name)

    @Property_RO
    def degrees0(self):
        '''Get this C{Ang} in C{degrees} ignoring the turns.
        '''
        return atan2d(*self.sc2)  # XXX Degrees(d, self.name)

    divmod = __divmod__

    @staticmethod
    def EPS0(**unit_name):
        '''Get a tiny C{Ang}.

           @note: This allows angles extremely close to the cardinal
                  directions to be generated.  The C{.round} method
                  will flush this angle to C{0}.
        '''
        return Ang(_EPS03, 1, **unit_name)

    @staticmethod
    def _flip(bet, omg, alp=None):
        '''(INTERNAL) Reflect C{bet}, C{omg} and C{alp} inplace.
        '''  # Ellipsoid3.Flip
        bet.reflect(flipc=True)
        omg.reflect(flips=True)
        if alp:
            alp.reflect(flips=True, flipc=True)

    def flipsign(self, mul=-1, **name):
        '''Copy this C{Ang} with sign flipped.
        '''
        r = (-self) if signBit(mul) else self
        return self._Ang(r, **name) if name else r

    def _float1(self, f, **name):
        # return C{f} as C{Ang} in this C{unit}
        return _Ang_from[self.unit](f, **name)

    def _float2(self, other=None):
        # get self and C{other} as floats
        r = other if other is None or isinstance(other, int) else \
            float(_Ang_from[self.unit](other))
        return float(self), r

    def _ifloat(self, f):  # PYCHOK expected
        # set self to C{f} degrees or radians
        scn = self._float1(f).scn3
        self._s, self._c, self._n = scn
        return self._update()

    @staticmethod
    def fromDegrees(deg, **unit_name):
        '''Get an C{Ang} from degrees.
        '''
        if isAng(deg):
            s, c, n = deg.scn3
            d       = deg.degrees0
        elif _isDegrees(deg, iscalar=True):
            s, c = sincos2d(deg)
            d    = atan2d(s, c)
            n    = round((deg - d) / _360_0)
        else:
            _raiseError(Ang.fromDegrees, deg, **unit_name)
        a = Ang(s, c, n, **_xkwds(unit_name, unit=Degrees))
        a.__dict__.update(degrees0=d)  # Property_RO
        return a

    @staticmethod
    def fromLambertian(psi, **unit_name):
        '''Get an C{Ang} from C{lamberterian} radians.
        '''
        s = psi.lambertian if isAng(psi) else sinh(psi)
        return Ang(s, 1, normal=False, **_xkwds(unit_name, unit=Lambertian))

    @staticmethod
    def fromRadians(rad, **unit_name):
        '''Get an C{Ang} from radians.
        '''
        if isAng(rad):
            s, c, n = rad.scn3
            r       = rad.radians0
        elif _isRadians(rad, iscalar=True):
            s, c = sincos2(rad)
            r = atan2(s, c)
            n = round((rad - r) / PI2)
        else:
            _raiseError(Ang.fromRadians, rad, **unit_name)
        a = Ang(s, c, n, **_xkwds(unit_name, unit=Radians))
        a.__dict__.update(radians0=r)  # Property_RO
        return a

    @staticmethod
    def fromScalar(ang, **unit_name):
        '''Get an C{Ang} from C{Degrees}, C{Radians} or another C{Ang}.
        '''
        if isAng(ang):
            r = Ang(ang, **_xkwds(unit_name, unit=ang.unit))
        else:
            u = _xkwds_get(unit_name, unit=None)
            if u is Lambertian:
                r = Ang.fromLambertian(ang, **unit_name)
            elif _isDegrees(ang, iscalar=u is Degrees):
                r = Ang.fromDegrees(ang, **unit_name)
            elif _isRadians(ang, iscalar=u is Radians):
                r = Ang.fromRadians(ang, **unit_name)
            else:
                _raiseError(Ang.fromScalar, ang, **unit_name)
        return r

    def is_integer(self, *n):
        '''Is this C{Ang}'s degrees C{integer}? (C{bool}).
        '''
        return self.toDegrees(*n).is_integer()

    def isnear0(self, eps0=EPS0):  # aka zerop
        '''Is this C{Ang} near C{0} within a tolerance?
        '''
        s, c, n = self.scn3
        return bool(n == 0 and c > 0 and fabs(s) <= eps0)

    def _kwds(self, kwds, **dflt):
        return _xkwds(kwds, **_xkwds(dflt, unit=self.unit,
                                           name=self.name))

    @Property_RO
    def lambertian(self):
        '''Get this C{Ang}'s Lambertian, C{asinh(tan(radians))}.
        '''
        return asinh(self.t)  # XXX Lambertian(self.t)

    def mod(self, mul=_1_0, **unit_name):
        '''Return the I{reduced latitude} C{atan(B{mul} *
           tan(B{this}))} as an C{Ang}.

           @arg mul: Factor (C{scalar}, positive).

           @note: The quadrant of the result tracks that of
                  this C{Ang} through multiples turns.
        '''
        kwds = self._kwds(unit_name)
        if signBit(mul):
            r = self._Ang(Ang.NAN(), **kwds)
        else:
            s, c, n = self.scn3
            if mul > 1:
                c = c / mul  # /= chokes PyChecker
            else:  # mul <= 1
                s *= mul
            r = self._Ang(s, c, n, normal=False, **kwds)
        return r

    @staticmethod
    def N(**unit_name):
        '''Get North C{Ang}.
        '''
        return Ang(0, 1, **unit_name)

    @property
    def n(self):
        '''Return the number of turns (C{float}) or C{0.0}.
        '''
        return self._n or _0_0

    @n.setter  # PYCHOK setter!
    def n(self, n):
        self._n_0(_fint(n))

    def _n_0(self, n):
        '''(INTERNAL) Set C{n} or C{n0}.
        '''
        if self._n != n:
            self._n, n = n, self._n
            self._update()
        return n

    @property
    def n0(self):
        '''Return the number of turns, treating C{-180} as C{180 - 1 turn} (C{float}).
        '''
        return (self.n - self._n01) or _0_0

    @n0.setter  # PYCHOK setter!
    def n0(self, n):
        self._n_0(_fint(n) + self._n01)

    @Property_RO
    def _n01(self):
        s, c = self.sc2
        return int(c < 0 and s == 0 and signBit(s))

    @staticmethod
    def NAN(**unit_name):
        '''Get an invalid C{Ang}.
        '''
        return Ang(NAN, NAN, **unit_name)

    @Property_RO
    def ncardinal(self):
        '''Get the nearest cardinal direction (C{float_int}).

           @note: This is the reverse of C{cardinal}.
        '''
        return _ncardinal(*self.scn3)

    def nearest(self, ind=0, **name):
        '''Return the closest cardinal direction (C{Ang}).

           @arg ind: An indicator, if C{B{ind}=0} the closest cardinal
                     direction, otherwise, if B{C{ind}} is even, the
                     closest even (N/S) cardinal direction or if B{C{ind}}
                     is odd, the closest odd (E/W) cardinal direction.
        '''
        s, c, n = self.scn3
        p = (ind == 0 and fabs(s) > fabs(c)) or (ind & 1)
        s, c = _orthogonal2(p, s, c)
        return self._Ang(s, c, n, **self._kwds(name))

    @staticmethod
    def _norm(bet, omg, alp=None, alt=False):
        '''(INTERNAL) Put C{bet}, C{ong} and C{alp} in range.
        '''  # Ellipsoid3.AngNorm
        flip = signBit(omg.s if alt else bet.c)
        if flip:
            Ang._flip(bet, omg, alp)
        return flip

    def normalize(self, *n):
        '''Re-normalize this C{Ang}, optionally replacing turns.
        '''
        sc = _normalize2(*self.sc2)
        if n:
            self.n, n = n[0], self.n
            if self.n != n:  # updated
                self._s, self._c = sc
                return self
        return self._update(*sc)

    def _other(self, other):
        # get C{other} as C{Ang} from C{unit}
        return other if isAng(other) else _other(other, self.unit)

    pow = __pow__

    @Property_RO
    def _quadrant(self):
        s, c = map(int, map(signBit, self.sc2))
        return s + s + (c ^ s)

    @property_doc_("this C{Ang}'s quadrant (C{int} 0..3)")
    def quadrant(self):
        return self._quadrant

    @quadrant.setter  # PYCHOK setter!
    def quadrant(self, quadrant):
        s, c = map(fabs, self.sc2)
        q = int(quadrant)
        if (q & 2):
            s = -s  # _copysign(self.s, -1 if (q & 2) else 1)
        if (((q >> 1) ^ q) & 1):
            c = -c  # _copysign(self.c, -1 if (((q >> 1) ^ q) & 1) else 1)
        self._update(s, c)

    @Property_RO
    def radians(self):
        '''Get this C{Ang} in C{radians}.
        '''
        r = self.radians0
        if self.n:
            r += self.n * PI2
        return r  # XXX Radians(r, self.name)

    @Property_RO
    def radians0(self):
        '''Get this C{Ang} in C{radians} ignoring the turns.
        '''
        return atan2(*self.sc2)  # XXX Radians(r, self.name)

    def reflect(self, flips=False, flipc=False, swapsc=False):
        '''Reflect this C{Ang} in various ways.

           @kwarg flips: Flip the sign of C{s}.
           @kwarg flipc: Flip the sign of C{c}.
           @kwarg swapsc: Swap C{s} and C{c}.

           @note: The operations are carried out in the order
                  of the arguments.
        '''
        s, c = self.sc2
        if flips:
            s = -s
        if flipc:
            c = -c
        if swapsc:
            s, c = c, s
        return self._update(s, c)

    def round(self, *ndigits, **name):
        '''Return this C{Ang}, optionally rounded to C{ndigits} (C{Ang}).
        '''
        s, c, n = self.scn3
        if ndigits:
            s = round(s, *ndigits)
            c = round(c, *ndigits)
        else:
            s, c = map1(_rnd, s, c)
        return self._Ang(s, c, n, **self._kwds(name))

    @property_RO
    def s(self):
        '''Get the sine of this C{Ang} (C{float}).
        '''
        return self._s

    @property_RO
    def sc2(self):
        '''Get the 2-tuple C{(s, c)}.
        '''
        return self.s, self.c

    @Property_RO
    def scn3(self):
        '''Get the 3-tuple C{(s, c, n)}.
        '''
        return self.s, self.c, self.n

    @property_RO
    def scnu4(self):
        '''Get the 4-tuple C{(s, c, n, unit)}.
        '''
        return self.s, self.c, self.n, self.unit

    def shift(self, q=0, **unit_name):
        '''Shift this C{Ang} by C{q} I{quarter} turns (C{scalar}).
        '''
        kwds = self._kwds(unit_name)
        if _isfinite(q):
            s = self.copy(**kwds)
            if q:
                s -= Ang.cardinal(q)
        else:
            s = Ang.NAN(**kwds)
        return s

    def signOf(self, *n):
        '''Determine this C{Ang}'s sign, optionally replacing the turns.

           @return: The sign (C{int}, -1, 0 or +1).
        '''
        return _signOf(self.toDegrees(*n), 0)

    @Property_RO
    def t(self):
        '''Get the tangent of this C{Ang} (C{float}).
        '''
        return _over(*self.sc2)

    def toDegrees(self, *n):
        '''Return this C{Ang} as C{Degrees}, optionally replacing the turns.
        '''
        if n:
            d = self.degrees0
            n = float(n[0])
            if n:
                d += n * _360_0
        else:
            d = self.degrees
        return Degrees(d, self.name)

    def toLambertian(self, **name):
        '''Return this C{Ang} as L{Lambertian}.
        '''
        name = _xkwds(name, name=self.name)
        return Lambertian(self.lambertian, **name)

    def toRadians(self, *n):
        '''Return this C{Ang} as C{Radians}, optionally replacing the turns.
        '''
        if n:
            r = self.radians0
            n = float(n[0])
            if n:
                r += n * PI2
        else:
            r = self.radians
        return Radians(r, self.name)

    def toRepr(self, *n, **prec_fmt):  # PYCHOK signature
        '''Return this C{Ang} as C{"<name>(<value>)"} with/out turns (C{str}).
        '''
        return self.toUnit(*n).toRepr(**prec_fmt)

    def toStr(self, *n, **prec_fmt):  # PYCHOK signature
        '''Return this C{Ang} as C{"<value>"} with/out turns (C{str}).
        '''
        return self.toUnit(*n).toStr(**prec_fmt)

    def toTuple(self, **prec_fmt_sep):
        '''Return string C{"(s, c, n)"} or tuple C{('s', 'c', 'n')} if C{sep is None}.
        '''
        return fstr(self.scn3, **prec_fmt_sep)

    def toUnit(self, *n):
        '''Return this C{Ang} as C{self.unit}s, optionally replacing the turns.
        '''
        u = self.unit
        return self.toRadians(*n)  if u is Radians else (
               self.toDegrees(*n)  if u is Degrees else (
               self.toLambertian() if u is Lambertian else
              _raiseError(self.toUnit, u)))  # PYCHOK indent

    @property_doc_(' the scalar unit to L{Degrees} or L{Radians}')
    def unit(self):
        return self._unit

    @unit.setter  # PYCHOK setter!
    def unit(self, unit):
        if unit not in _Ang_types:  # PYCHOK no cover
            _raiseError(Ang.unit, unit)
        if self._unit != unit:
            self._unit = unit

    def _update(self, *sc):
        if sc:
            if sc == self.sc2:
                return self
            self._s, self._c = sc
        _update_all(self)
        return self

_allPropertiesOf_n(14, Ang)  # PYCHOK assert


class _Ang3Tuple(_NamedTuple):
    '''(INTERNAL) Methods C{.toDegrees}, C{.toLambertian}, C{.toRadians} and C{.toUnit}.
    '''
    _Names_ = (Ang.__name__,) * 3  # needed for ...
    _Units_ =  Ang, Ang, _Pass  # ...testNamedTuples

    def toDegrees(self, *n, **fmt_prec_sep):
        '''Change any C{Ang} to C{unit Degrees} or to C{Degrees.toStr} if any B{C{fmt_prec_sep}}.
        '''
        t = self.toUnit(Degrees, *n)
        if fmt_prec_sep:  # see C{Degrees.toStr}
            sep, fmt_prec = _xkwds_pop2(fmt_prec_sep, sep=_COMMASPACE_)
            s = self.toStr(sep=None) if sep else self
            t = (a.toStr(**fmt_prec) if isAng(a) else s for a, s in zip(t, s))
            t = Fmt.PAREN(sep.join(t)) if sep else tuple(t)
        return t

    def toLambertian(self):
        '''Change any C{Ang} to C{unit Lambertian}.
        '''
        return self.toUnit(Lambertian)

    def toRadians(self, *n):
        '''Change any C{Ang} to C{unit Radians}.
        '''
        return self.toUnit(Radians, *n)

    def toUnit(self, unit, *n):
        '''Change any C{Ang} to C{unit}, .
        '''
        for a in self:
            if isAng(a):  # and a.unit is not unit:
                a.unit = unit
                if n:
                    a.n = n[0]
        return self


class Lambertian(Radians):
    '''A C{Lambertian} in C{radians}.
    '''
    def __new__(cls, *args, **kwds):
        return Radians.__new__(cls, *args, **_xkwds(kwds, name='psi'))


_Ang_from = {Radians:    Ang.fromRadians,
             Degrees:    Ang.fromDegrees,
             Lambertian: Ang.fromLambertian}
_Ang_types = tuple(_Ang_from.keys())  # PYCHOK used!


def Ang_(s, c=None, n=1, **unit_name):
    '''(INTERNAL) New, non-normal C{Ang}.
    '''
    return Ang(s, c, n, **_xkwds(unit_name, normal=False))


def Deg(deg, **name):
    '''Return an L{Ang} from C{deg} degrees or an other L{Ang}.
    '''
    return Ang(deg, unit=Degrees, **name)


def isAng(ang):
    '''Is C{ang} an L{Ang} instance?
    '''
    return isinstance(ang, Ang)


def Rad(rad, **name):
    '''Return an L{Ang} from C{rad} radians or an other L{Ang}.
    '''
    return Ang(rad, unit=Radians, **name)


def _SinCos2(ang, *unit):
    '''Get C{sin} and C{cos} of an L{Ang}, any I{typed} C{ang}le
       or C{unit} if C{ang}le is scalar.

       @see: Function L{SinCos2<pygeodesy.utily.SinCos2>}.
    '''
    return ang.sc2 if isAng(ang) else SinCos2(ang, *unit)

# **) MIT License
#
# Copyright (C) 2025-2026 -- mrJean1 at Gmail -- All Rights Reserved.
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
