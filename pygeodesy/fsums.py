
# -*- coding: utf-8 -*-

u'''Class L{Fsum} for running, precision floating point summation.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import _isfinite, isint, isscalar, neg
from pygeodesy.errors import _NotImplementedError, _OverflowError, \
                             _TypeError, _ValueError, _xkwds_get
from pygeodesy.interns import NN, _COMMASPACE_, _finite_, _iadd_, _not_, \
                             _SPACE_, _supported_, _0_0, _1_0, _N_1_0
from pygeodesy.lazily import _ALL_LAZY, _sys_version_info2
from pygeodesy.named import _Named, _NamedTuple, _NotImplemented
# from pygeodesy.props import property_RO  # from .units
from pygeodesy.streprs import Fmt, pairs, unstr
from pygeodesy.units import Float, property_RO

from math import ceil as _ceil, floor as _floor

__all__ = _ALL_LAZY.fsums
__version__ = '22.01.21'

_idiv_ = '/='
_imul_ = '*='
_ipow_ = '**='
_isub_ = '-='

_non_zero_ = 'non-zero'
_residual_ = 'residual'


def _2even(s, r, p):
    '''(INTERNAL) Half-even rounding.
    '''
    if (r > 0 and p > 0) or \
       (r < 0 and p < 0):  # signs match
        t, p = _2sum(s, p * 2)
        if not p:
            s = t
    return s


def _2float(index=None, **name_value):
    '''(INTERNAL) Raise C{TypeError} or C{ValueError} if not scalar or infinite.
    '''
    n, v = name_value.popitem()
    try:
        if _isfinite(v):
            return v if isinstance(v, float) else float(v)
        X, t = _ValueError, _not_(_finite_)
    except TypeError as x:
        X, t = _TypeError, str(x)
    except ValueError as x:
        X, t = _ValueError, str(x)
    except Exception as x:
        X, t = _NotImplementedError, repr(x)
    if index is not None:
        n = Fmt.SQUARE(n, index)
    raise X(n, v, txt=t)


def _2Fsum(x, name=NN):
    '''(INTERNAL) Return B{C{x}} as L{Fsum} instance.
    '''
    return x.fcopy(name=name) if isinstance(x, Fsum) else \
           Fsum(x, name=name)


def _2sum(a, b):  # by .testFmath
    '''(INTERNAL) Precision C{2sum} of M{a + b} as 2-tuple (sum, residual).
    '''
    s = a + b
    if not _isfinite(s):
        raise _OverflowError(unstr(_2sum.__name__, a, b), txt=str(s))
    if abs(a) < abs(b):
        a, b = b, a
    return s, (b - (s - a))  # abs(b) <= abs(a)


class _ResidualError(Exception):
    pass


class Fsum(_Named):
    '''Precision summation similar to standard Python function C{math.fsum}.

       Unlike C{math.fsum}, this class accumulates the values and provides
       intermediate, precision running sums.  Accumulation may continue
       after intermediate summations.

       @note: Handling of exceptions, C{inf}, C{INF}, C{nan} and C{NAN}
              values is different from C{math.fsum}.

       @see: U{Hettinger<https://GitHub.com/ActiveState/code/blob/master/recipes/Python/
             393090_Binary_floating_point_summatiaccurate_full/recipe-393090.py>},
             U{Kahan<https://WikiPedia.org/wiki/Kahan_summation_algorithm>},
             U{Klein<https://Link.Springer.com/article/10.1007/s00607-005-0139-x>},
             Python 2.6+ file I{Modules/mathmodule.c} and the issue log
             U{Full precision summation<https://Bugs.Python.org/issue2819>}.
    '''
    _fsum2_ = None
    _n      = 0
    _ps     = []  # partials

    def __init__(self, *starts, **name):
        '''Initialize a new accumulator with one or more start values.

           @arg starts: No, one or more start values (C{scalar}s).
           @kwarg name: Optional name (C{str}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{starts}} value.

           @raise ValueError: Invalid or non-finite B{C{starts}} value.
        '''
#       self._n  = 0
        self._ps = []
        if name:
            self.name = _xkwds_get(name, name=NN)
        if starts:
            self.fadd(starts)

    def __abs__(self):
        '''Return absolute value of this instance.
        '''
        return self.__neg__() if self < 0 else self

    def __add__(self, other):
        '''Sum of this and a scalar or an other instance.

           @arg other: L{Fsum} instance or C{scalar}.

           @return: The sum (L{Fsum}).

           @see: Method L{Fsum.__iadd__}.
        '''
        f  = self.fcopy(name=self.__add__.__name__)
        f += other
        return f

    def __bool__(self):  # PYCHOK PyChecker
        '''Is this instance non-zero?.
        '''
        return self.__ne__(0)

    def __ceil__(self):  # PYCHOK not special in Python 2-
        '''Return the C{ceil} of this instance as C{float}.

           @see: Methods L{__floor__} and L{__int__}.
        '''
        s, p = self._cmp2()
        return _ceil(max(p, neg(s)))

    def __divmod__(self, other):
        '''Return C{divmod(this_instance, B{other})} as 2-tuple
           C{(quotient, remainder)} both C{float}.

           @see: Method L{__itruediv__}.
        '''
        f  = self.fcopy(name=self.__divmod__.__name__)
        i  = int(f / other)  # no // __floordiv__
        f -= other * i
        return float(i), f

    def __eq__(self, other):
        '''Compare this and an other instance or scalar.
        '''
        s, p = self._cmp2(other)
        return p == s

    def __float__(self):
        '''Convert this instance to C{float} as C{float(self.fsum())}.
        '''
        s = self._fsum2_
        return self.fsum() if s is None else s

    def __floor__(self):  # PYCHOK not special in Python 2-
        '''Return the C{floor} of this instance as C{float}.

           @see: Methods L{__ceil__} and L{__int__}.
        '''
        s, p = self._cmp2()
        return _floor(min(p, neg(s)))

    def __floordiv__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __format__(self, *other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, *other)

    def __ge__(self, other):
        '''Compare this and an other instance or scalar.
        '''
        s, p = self._cmp2(other)
        return p >= s

    def __gt__(self, other):
        '''Compare this and an other instance or scalar.
        '''
        s, p = self._cmp2(other)
        return p > s

    def __hash__(self):  # PYCHOK no cover
        '''Return the C{hash} of this instance.
        '''
        return hash(self._ps)

    def __iadd__(self, other):
        '''Add a scalar or an other instance to this instance.

           @arg other: L{Fsum} instance or C{scalar}.

           @return: This instance, updated (L{Fsum}).

           @raise TypeError: Invalid B{C{other}} type.

           @see: Method L{Fsum.fadd}.
        '''
        if isscalar(other):
            if other:
                self.fadd_(other)
        elif isinstance(other, Fsum):
            if other is self:  # or self.__eq__(other):
                self.fmul(2)
            else:
                self.fadd(other._ps)
        else:
            raise self._Error(_iadd_, other)
        return self

    def __imatmul__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __imul__(self, other):
        '''Multiply this instance by a scalar or an other instance.

           @arg other: L{Fsum} instance or C{scalar}.

           @return: This instance, updated (L{Fsum}).

           @raise TypeError: Invalid B{C{other}} type.

           @see: Method L{Fsum.fmul}.
        '''
        if isscalar(other):
            self.fmul(other)
        elif isinstance(other, Fsum):
            ps = list(other._ps)  # copy
            p = ps.pop() if ps else 0
            if p:
                s = self.fcopy()
                self.fmul(p)
                while ps:  # self += s * ps.pop()
                    p = ps.pop()
                    if p:
                        p = s.fcopy().fmul(p)
                        self.fadd(p._ps)
            else:  # PYCHOK no cover
                self._iset(_0_0)
        else:
            raise self._Error(_imul_, other)
        return self

    def __int__(self):
        '''Convert this instance to C{int} as C{int(self.fsum() + partials)}.

           @see: Methods L{__ceil__} and L{__floor__}.
        '''
        s, p = self._cmp2()
        s = neg(s)
        return int(p if (s > 0 and p < s) or
                        (s < 0 and p > s) else s)

    def __ipow__(self, other, *mod):  # PYCHOK 2 vs 3 args
        '''Raise this instance to power B{C{other}}.

           @arg other: The exponent (C{scalar}).
           @arg mod: Not implemented (C{scalar}).

           @return: This instance, updated (L{Fsum}).

           @raise NotImplementedError: Argument B{C{mod}} used.

           @raise TypeError: Non-scalar B{C{other}}.

           @raise ValueError: Fractional B{C{other}} and this
                              instance is negative or has a
                              non-zero C{residual} or negative
                              B{C{other}} and this instance C{0}.

           @see: CPython function U{float_pow<https://GitHub.com/
                 python/cpython/blob/main/Objects/floatobject.c>}.
        '''
        if mod:
            return _NotImplemented(self, other, *mod)
        x = other
        if not isscalar(x):
            raise self._Error(_ipow_, other)

        s, r = self.fsum2()
        if isint(x, both=True):
            x = int(x)  # integer exponent
            if not -2 < x < 2:
                if r:
                    s = self.fcopy()
                    for _ in range(1, abs(x)):
                        s *= self
                    if x < 0:
                        s, r = s.fsum2()
                        if r:
                            raise self._Erres(_ipow_, other)
                        # use **= -1 for the CPython float_pow
                        # error if s is zero, and not s = 1 / s
                        x = -1
                    else:
                        x = None
            elif x < 0:  # x == -1
                if r:
                    raise self._Erres(_ipow_, other)
            else:  # x == 0 or x == 1
                s = self if x else _1_0
                x = None
        elif r:  # fractional exponent
            raise self._Erres(_ipow_, other)
        elif s < 0:  # negative**fractional yields complex
            raise self._Error(_ipow_, other, Error=_ValueError,
                              txt=_not_(_supported_))
        if x is not None:
            try:
                s **= x
            except Exception as X:
                raise self._Error(_ipow_, other, Error=_ValueError,
                                  txt=str(X))
        return self._iset(s)  # s is float or an Fsum, perhaps self

    def __isub__(self, other):
        '''Subtract a scalar or an other instance from this instance.

           @arg other: L{Fsum} instance or C{scalar}.

           @return: This instance, updated (L{Fsum}).

           @raise TypeError: Invalid B{C{other}} type.

           @see: Method L{Fsum.fadd}.
        '''
        if isscalar(other):
            if other:
                self.fadd_(neg(other))
        elif isinstance(other, Fsum):
            if other is self:  # or self.__eq__(other):
                self._iset(_0_0)
            else:
                self.fadd(map(neg, other._ps))
        else:
            raise self._Error(_isub_, other)
        return self

    def __itruediv__(self, other):
        '''Devide this instance by a I{scalar} divisor only.

           @arg other: The denominator (C{scalar}).

           @raise TypeError: Non-scalar B{C{other}}.

           @return: This instance, updated (L{Fsum}).

           @raise ValueError: Zero, invalid or non-finite B{C{other}}
                              or an L{Fsum} B{C{other}} with non-zero
                              C{partials}.

           @raise ZerDivisionError: Zero B{C{other}}.
        '''
        if isscalar(other):
            d = other
        elif isinstance(other, Fsum):
            if other is self or self.__eq__(other):
                self._iset(_1_0)
                d = 1
            else:
                f = other.fcopy()
                d, r = f.fsum2()
                if r:
                    raise _ValueError(_SPACE_(self, _idiv_, f.toRepr()),
                                      txt=_SPACE_(_non_zero_, _residual_))
                # d == 0 throws ZeroDivisionError
        else:
            raise self._Error(_idiv_, other)
        if d not in (1, _1_0):
            self.fdiv(d)
        return self

    def __le__(self, other):
        '''Compare this and an other instance or scalar.
        '''
        s, p = self._cmp2(other)
        return p <= s

    def __len__(self):
        '''Return the I{total} number of accumulated values (C{int}).
        '''
        return self._n

    def __lt__(self, other):
        '''Compare this and an other instance or scalar.
        '''
        s, p = self._cmp2(other)
        return p < s

    def __matmul__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __mod__(self, other):
        '''Return C{this_instance % B{other}}.

           @see: Method L{__divmod__}.
        '''
        _, m = divmod(self, other)
        return m

    def __mul__(self, other):
        '''Product of this and an other instance or a scalar.

           @arg other: L{Fsum} instance or C{scalar}.

           @return: The product (L{Fsum}).

           @see: Method L{Fsum.__imul__}.
        '''
        f  = self.fcopy(name=self.__mul__.__name__)
        f *= other
        return f

    def __ne__(self, other):
        '''Compare this and an other instance or scalar.
        '''
        s, p = self._cmp2(other)
        return p != s

    def __neg__(self):
        '''Return a copy of this instance, negated.
        '''
        f = self.fcopy(name=self.__neg__.__name__)
        f._ps[:]  = map(neg, f._ps)
        f._fsum2_ = None
        return f

    def __pos__(self):  # PYCHOK no cover
        '''Return this instance, I{as-is}.
        '''
        return self

    def __pow__(self, other, *mod):  # PYCHOK 2 vs 3 args
        '''Return C{self ** other} as L{Fsum}, see L{Fsum.__ipow__}.'''
        f = self.fcopy()  # deep=False
        f.__ipow__(other, *mod)
        return f

    __radd__ = __add__

    def __rdivmod__ (self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __rfloordiv__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __rmatmul__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __rmod__(self, other):  # PYCHOK no cover
        '''Return C{other % self} as L{Fsum}.'''
        f = _2Fsum(other, name=self.__rmod__.__name__)
        _, m = divmod(f, self)  # %= chokes PyChecker
        return m

    __rmul__ = __mul__

    def __round__(self, ndigits=None):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, ndigits=ndigits)

    def __rpow__(self, other, *mod):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other, *mod)

    def __rsub__(self, other):
        '''Return C{other - this} as L{Fsum}.'''
        f  = _2Fsum(other, name=self.__rsub__.__name__)
        f -=  self
        return f

    def __rtruediv__(self, other):  # PYCHOK no cover
        '''Return C{other / self} as L{Fsum}.'''
        f = _2Fsum(other, name=self.__rtruediv__.__name__)
        f.__itruediv__(self)  # /= chokes PyChecker
        return f

    def __sizeof__(self):  # PYCHOK no cover
        '''Size of this instance in C{bytes}.
        '''
        return sum(p.__sizeof__() for p in (self._ps + [self._ps]))

    def __str__(self):
        '''Default C{str(self)}.
        '''
        return self.toStr()

    def __sub__(self, other):
        '''Difference of this and an other instance or a scalar.

           @arg other: L{Fsum} instance or C{scalar}.

           @return: The difference (L{Fsum}).

           @see: Method L{Fsum.__isub__}.
        '''
        f  = self.fcopy(name=self.__sub__.__name__)
        f -= other
        return f

    def __truediv__(self, other):
        '''Quotient of this instance and a I{scalar}.

           @arg other: The denominator (C{scalar}).

           @raise TypeError: Non-scalar B{C{other}}.

           @raise ValueError: Zero, invalid or non-finite B{C{other}}
                              or an L{Fsum} B{C{other}} with too many
                              C{partials}.

           @return: The quotient (L{Fsum}).
        '''
        f = self.fcopy(name=self.__truediv__.__name__)
        f.__itruediv__(other)  # /= chokes PyChecker
        return f

    __trunc__ = __int__

    if _sys_version_info2 < (3, 0):  # PYCHOK no cover
        # <https://docs.Python.org/2/library/operator.html#mapping-operators-to-functions>
        __div__     = __truediv__
        __idiv__    = __itruediv__
        __long__    = __int__
        __nonzero__ = __bool__

    def _cmp2(self, *other):
        '''(INTERNAL) Diff this and another instance or C{0}.
        '''
        f = self.fcopy()
        if other:
            f = f.__isub__(*other)
        s = neg(f.fsum())  # negative sum!
        f = f._ps
        p = fsum1_(fsum(p for p in f if p > 0),
                   fsum(p for p in f if p < 0)) if f else _0_0
        return s, p

    def _Erres(self, op, other):
        '''(INTERNAL) Return an residual/partials C{ValueError}.
        '''
        return self._Error(op, other, Error=_ValueError,
                           txt=_SPACE_(_non_zero_, _residual_))

    def _Error(self, op, other, Error=_TypeError, **txt):
        '''(INTERNAL) Return an operation B{C{Error}}.
        '''
        return Error(_SPACE_(self.toRepr(), op, repr(other)), **txt)

    def _iset(self, other):
        '''(INTERNAL) Override this instance with an other.
        '''
        if isscalar(other):
            self._fsum2_ = f = _2float(other=other)
            self._n     =  1
            self._ps[:] = [f] if f else []
        elif isinstance(other, Fsum):  # PYCHOK no cover
            if other is not self:
                self._fsum2_ = other._fsum2_
                self._n      = other._n
                self._ps[:]  = other._ps
        else:
            raise self._Error(_isub_, other)
        return self

    def fadd(self, xs):
        '''Accumulate more scalar values from an iterable.

           @arg xs: Iterable, list, tuple, etc. (C{scalar}s).

           @return: This instance (L{Fsum}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{xs}} value.

           @raise ValueError: Invalid or non-finite B{C{xs}} value.
        '''
        if isscalar(xs):  # for backward compatibility
            xs = (xs,)  # PYCHOK no cover

        ps, n = self._ps, -1
        for n, x in enumerate(xs):  # _iter()
            if x:
                x = _2float(xs=x, index=n)
                i = 0
                for p in ps:
                    x, p = _2sum(x, p)
                    if p:
                        ps[i] = p
                        i += 1
                ps[i:] = [x]
        # assert self._ps is ps
        self._n += n + 1
        self._fsum2_ = None
        return self

    def fadd_(self, *xs):
        '''Accumulate more I{scalar} values from positional arguments.

           @arg xs: Values to add (C{scalar}s), all positional.

           @return: This instance (L{Fsum}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{xs}} value.

           @raise ValueError: Invalid or non-finite B{C{xs}} value.
        '''
        return self.fadd(xs)

    def fcopy(self, deep=False, name=NN):
        '''Copy this instance, C{shallow} or B{C{deep}}.

           @return: The copy (L{Fsum}).
         '''
        f = _Named.copy(self, deep=deep, name=name)  # see .__neg__
#       f._fsum2_  = self._fsum2_
        f._n       = self._n if deep else 1
        f._ps = list(self._ps)  # list copy
        return f

    copy = fcopy

    def fdiv(self, divisor):
        '''Devide this instance by a I{scalar}.

           @arg divisor: The denominator (C{scalar}).

           @raise TypeError: Non-scalar B{C{divisor}}.

           @raise ValueError: Zero, invalid or non-finite B{C{divisor}}.

           @return: This instance (L{Fsum}).
        '''
        try:
            self.fmul(_1_0 / _2float(divisor=divisor))
        except (TypeError, ValueError, ZeroDivisionError) as x:
            raise self._Error(_idiv_, divisor, Error=_ValueError, txt=str(x))
        return self

    def fmul(self, factor):
        '''Multiple this instance by a I{scalar}.

           @arg factor: The multiplier (C{scalar}).

           @raise TypeError: Non-scalar B{C{factor}}.

           @raise ValueError: Invalid or non-finite B{C{factor}}.

           @return: This instance (L{Fsum}).

           @see: Method L{Fsum.fadd}.
        '''
        f, ps = _2float(factor=factor), self._ps
        if ps:
            if abs(f) != 1:
                # multiply and adjust partial sums
                ps[:] = [p * f for p in ps]
                self.fadd_(ps.pop())
                self._n -= 1
            elif f < 0:  # == -1
                ps[:] = map(neg, ps)
                self._fsum2_ = None
        # assert self._ps is ps
        return self

    def fsub(self, xs):
        '''Subtract several I{scalar} values.

           @arg xs: Iterable, list, tuple. etc. (C{scalar}s).

           @return: This instance (L{Fsum}).

           @see: Method L{Fsum.fadd}.
        '''
        if xs:
            self.fadd(map(neg, xs))
        return self

    def fsub_(self, *xs):
        '''Subtract all I{scalar} positional values.

           @arg xs: Values to subtract (C{scalar}s), all positional.

           @return: This instance (L{Fsum}).

           @see: Method L{Fsum.fadd}.
        '''
        return self.fsub(xs)

    def fsum(self, xs=None):
        '''Accumulate more I{scalar} values and sum all.

           @kwarg xs: Iterable, list, tuple, etc. (C{scalar}s).

           @return: Accurate, running sum (C{float}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{xs}} value.

           @raise ValueError: Invalid or non-finite B{C{xs}} value.

           @note: Accumulation can continue after summation.
        '''
        if xs:
            self.fadd(xs)

        s = self._fsum2_
        if s is None:
            ps = self._ps
            i = len(ps) - 1
            if i < 0:
                s = _0_0
            else:
                s = ps[i]
                while i > 0:
                    i -= 1
                    s, p = _2sum(s, ps[i])
                    ps[i:] = [s]
                    if p:  # sum(ps) became inexact
                        ps.append(p)
                        if i > 0:  # half-even round if signs match
                            s = _2even(s, ps[i-1], p)
                        break
                # assert self._ps is ps
            self._fsum2_ = s
        return s

    def fsum_(self, *xs):
        '''Accumulate all I{scalar} positional values and sum all.

           @arg xs: Values to add (C{scalar}s), all positional.

           @return: Accurate, running sum (C{float}).

           @see: Method L{Fsum.fsum}.

           @note: Accumulation can continue after summation.
        '''
        return self.fsum(xs)

    def fsum2(self, xs=None):
        '''Accumulate more I{scalar} values and return the
           sum and residual.

           @kwarg xs: Iterable, list, tuple, etc. (C{scalar}s).

           @return: L{Fsum2Tuple}C{(fsum, residual)} with the
                    accurate, running C{fsum} and C{residual}
                    the precision sum of the remaining partials.

           @see: Methods L{Fsum.fsum} and L{Fsum.fsum2_}
        '''
        r, s = _0_0, self._fsum2_
        if s is None or xs:
            s = self.fsum(xs)
        if len(self._ps) > 1:
            r = fsum1_(neg(s), *self._ps)
        return Fsum2Tuple(s, r) if s else Fsum2Tuple(r, _0_0)

    def fsum2_(self, *xs):
        '''Accumulate all I{scalar} positional values and provide
           the sum and delta.

           @arg xs: Values to add (C{scalar}s), all positional.

           @return: 2-Tuple C{(sum, delta)} with the accurate,
                    running C{sum} and the C{delta} with the
                    previous running C{sum} (C{float}s).

           @see: Method L{Fsum.fsum_}.

           @note: Accumulation can continue after summation.
        '''
        p = self._fsum2_
        if p is None:
            p = self.fsum()
        s = self.fsum(xs)  # if xs else p
        return s, (s - p)

    @property_RO
    def imag(self):
        '''Return the imaginary part of this instance.
        '''
        return _0_0

    def is_integer(self):
        '''Return C{True} if this instance is an integer.
        '''
        s, p = self._cmp2()
        return s.is_integer() and p == neg(s)

    @property_RO
    def real(self):
        '''Return the real part of this instance.
        '''
        return float(self)

    def toRepr(self, prec=8, sep=_COMMASPACE_, fmt=Fmt.g, **unused):  # PYCHOK signature
        '''Return this C{Fsum} instance as representation.

           @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                        Trailing zero decimals are stripped for B{C{prec}} values
                        of 1 and above, but kept for negative B{C{prec}} values.
           @kwarg sep: Optional separator to join (C{str}).
           @kwarg fmt: Optional, C{float} format (C{str}).

           @return: This instance (C{str}).
        '''
        t = sep.join(pairs(self.fsum2().items(), prec=prec, fmt=fmt))
        return _SPACE_(Fmt.SQUARE(self.named3, len(self)), Fmt.PAREN(t))

    def toStr(self, prec=8, sep=_COMMASPACE_, fmt=Fmt.g, **unused):  # PYCHOK signature
        '''Return this C{Fsum} instance as string.

           @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                        Trailing zero decimals are stripped for B{C{prec}} values
                        of 1 and above, but kept for negative B{C{prec}} values.
           @kwarg sep: Optional separator to join (C{str}).
           @kwarg fmt: Optional, C{float} format (C{str}).

           @return: This instance (C{repr}).
        '''
        t = self.fsum2().toStr(prec=prec, sep=sep, fmt=fmt)
        return _SPACE_(Fmt.SQUARE(self.named3, len(self)), t)


class Fsum2Tuple(_NamedTuple):
    '''2-Tuple C{(fsum, residual)} with the accurate,
       running C{fsum} and C{residual} the precision
       sum of the remaining partials, both C{float}.
    '''
    _Names_ = (Fsum.fsum.__name__, _residual_)
    _Units_ = (Float,               Float)


try:
    from math import fsum  # precision IEEE-754 sum, Python 2.6+

    # make sure fsum works as expected (XXX check
    # float.__getformat__('float')[:4] == 'IEEE'?)
    if fsum((1, 1e101, 1, -1e101)) != 2:  # PYCHOK no cover
        del fsum  # nope, remove fsum ...
        raise ImportError  # ... use fsum below

except ImportError:

    def fsum(xs):
        '''Precision summation similar to standard Python function C{math.fsum}.

           Exception and I{non-finite} handling differs from C{math.fsum}.

           @arg xs: Iterable, list, tuple, etc. of values (C{scalar}s).

           @return: Accurate C{sum} (C{float}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{xs}} value.

           @raise ValueError: Invalid or non-finite B{C{xs}} value.

           @see: Class L{Fsum} and method L{Fsum.fsum}.
        '''
        return Fsum(name=fsum.__name__).fsum(xs) if xs else _0_0


def fsum_(*xs):
    '''Precision summation of all positional arguments.

       @arg xs: Values to be added (C{scalar}s).

       @return: Accurate L{fsum} (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise TypeError: Non-scalar B{C{xs}} value.

       @raise ValueError: Invalid or non-finite B{C{xs}} value.
    '''
    return fsum(map(float, xs))


def fsum1(xs):
    '''Precision summation, primed with C{1.0}.

       @arg xs: Iterable, list, tuple, etc. of values (C{scalar}s).

       @return: Accurate L{fsum} (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise TypeError: Non-scalar B{C{xs}} value.

       @raise ValueError: Invalid or non-finite B{C{xs}} value.
    '''
    def _xs(xs):
        yield _1_0
        for x in xs:
            yield x
        yield _N_1_0

    return fsum(_xs(xs)) if xs else _0_0


def fsum1_(*xs):
    '''Precision summation of a few arguments, primed with C{1.0}.

       @arg xs: Values to be added (C{scalar}s), all positional.

       @return: Accurate L{fsum} (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise TypeError: Non-scalar B{C{xs}} value.

       @raise ValueError: Invalid or non-finite B{C{xs}} value.
    '''
    return fsum1(xs)


# **) MIT License
#
# Copyright (C) 2016-2022 -- mrJean1 at Gmail -- All Rights Reserved.
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
