
# -*- coding: utf-8 -*-

u'''Precision floating point summation and utilities.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import copysign0, _isfinite, isint, isnear0, \
                             isscalar, len2, neg
from pygeodesy.errors import _IsnotError, LenError, _NotImplementedError, \
                             _OverflowError, _TypeError, _ValueError, \
                             _xkwds_get
from pygeodesy.interns import EPS0, EPS02, EPS1, MISSING, NN, PI, PI_2, PI_4, \
                             _finite_, _few_, _h_, _iadd_, _negative_, _not_, \
                             _singular_, _SPACE_, _too_, _0_0, _0_5, \
                             _1_0, _N_1_0, _1_5, _2_0, _3_0
from pygeodesy.lazily import _ALL_LAZY, _sys_version_info2
from pygeodesy.named import _Named, _NotImplemented
# from pygeodesy.props import property_RO  # from .units
from pygeodesy.streprs import Fmt, unstr
from pygeodesy.units import Int_, property_RO

from math import ceil as _ceil, floor as _floor, sqrt  # pow
from operator import mul as _mul

__all__ = _ALL_LAZY.fmath
__version__ = '22.01.16'

# sqrt(2) <https://WikiPedia.org/wiki/Square_root_of_2>
_0_4142 =  0.414213562373095  # sqrt(_2_0) - _1_0
_1_3rd  = _1_0 / _3_0
_2_3rd  = _2_0 / _3_0

_idiv_     = '/='
_imul_     = '*='
_isub_     = '-='
_partials_ = 'partials'


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
            return float(v)
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


def _2sum(a, b):  # by .testFmath
    '''(INTERNAL) Precision C{2sum} of M{a + b} as 2-tuple (sum, residual).
    '''
    s = a + b
    if not _isfinite(s):
        raise _OverflowError(unstr(_2sum.__name__, a, b), txt=str(s))
    if abs(a) < abs(b):
        a, b = b, a
    return s, (b - (s - a))  # abs(b) <= abs(a)


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
        self._n = 0
        self._ps = []
        if starts:
            self.fadd(starts)
        if name:
            self.name = _xkwds_get(name, name=NN)

    def __abs__(self):
        '''Return absolute value of this instance.
        '''
        return self.__neg__() if self.__lt__(0) else self

    def __add__(self, other):
        '''Sum of this and a scalar or an other instance.

           @arg other: L{Fsum} instance or C{scalar}.

           @return: The sum (L{Fsum}).

           @see: Method L{Fsum.__iadd__}.
        '''
        f = self.fcopy(name=self.__add__.__name__)
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
        return _ceil(max(p, -s))

    def __divmod__(self, other):
        '''Return the 2-tuple C{divmod(this_instance, B{other})}.

           @see: Method L{__itruediv__}.
        '''
        f = self.fcopy(name=self.__divmod__.__name__)
        i = f.__truediv__(other).__int__()
        return float(i), f.__isub__(other * i)

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
        return _floor(min(p, -s))

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
            if other not in (_1_0, 1):
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
                self._ps = []  # zero
                self._fsum2_ = _0_0
        else:
            raise self._Error(_imul_, other)
        return self

    def __int__(self):
        '''Convert this instance to C{int} as C{int(self.fsum() + partials)}.

           @see: Methods L{__ceil__} and L{__floor__}.
        '''
        s, p = self._cmp2()
        s = -s
        return int(p if (s > 0 and p < s) or
                        (s < 0 and p > s) else s)

    def __iset__(self, other):  # PYCHOK not special
        '''Override this instance with an other.

           @arg other: L{Fsum} instance or C{scalar}.

           @return: This instance, overridden (L{Fsum}).

           @raise TypeError: Invalid B{C{other}} type.
        '''
        if isscalar(other):
            self._n  =     1
            self._ps =    [float(other)]
            self._fsum2_ = float(other)
        elif isinstance(other, Fsum):
            if other is not self:
                self._n  =      other._n
                self._ps = list(other._ps)
                self._fsum2_ =  other._fsum2_
        else:
            raise self._Error(_isub_, other)
        return self

    def __isub__(self, other):
        '''Subtract a scalar or an other instance from this instance.

           @arg other: L{Fsum} instance or C{scalar}.

           @return: This instance, updated (L{Fsum}).

           @raise TypeError: Invalid B{C{other}} type.

           @see: Method L{Fsum.fadd}.
        '''
        if isscalar(other):
            if other:
                self.fadd_(-other)
        elif isinstance(other, Fsum):
            if other is self:  # or self.__eq__(other):
                self._ps = []  # zero
                self._fsum2_ = _0_0
            else:
                self.fadd(-p for p in other._ps)  # neg_(*other._ps)
        else:
            raise self._Error(_isub_, other)
        return self

    def __itruediv__(self, other):
        '''Devide this instance by a I{scalar} divisor only.

           @arg other: The denominator (C{scalar}).

           @raise TypeError: Non-scalar B{C{other}}.

           @raise ValueError: Zero, invalid or non-finite B{C{other}}
                              or an L{Fsum} B{C{other}} with too many
                              C{partials}.

           @return: This instance (L{Fsum}).
        '''
        if isscalar(other):
            d = other
        elif isinstance(other, Fsum):
            if other.__eq__(0):
                d = 0  # throw ZeroDivisionError
            elif other is self or self.__eq__(other):
                self._ps[:] = [_1_0]
                self._fsum2_ = _1_0
                return self
            else:
                f = other.fcopy()
                d = f.fsum()
                p = len(f._ps) - 1  # NOT len(f)!
                if p > 0:
                    raise self._Error(_idiv_, other, Error=_ValueError,
                                txt=_SPACE_(_partials_, Fmt.PAREN(p)))
        else:
            raise self._Error(_idiv_, other)
        if d not in (_1_0, 1):
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
        _, m = self.__divmod__(other)
        return m

    def __mul__(self, other):
        '''Product of this and an other instance or a scalar.

           @arg other: L{Fsum} instance or C{scalar}.

           @return: The product (L{Fsum}).

           @see: Method L{Fsum.__imul__}.
        '''
        f = self.fcopy(name=self.__mul__.__name__)
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
        f *= -1  # negates
        return f

    def __pos__(self):  # PYCHOK no cover
        '''Return this instance, I{as-is}.
        '''
        return self

    def __pow__(self, other, *mod):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other, *mod)

    def __radd__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

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
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __rmul__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __round__(self, ndigits=None):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, ndigits=ndigits)

    def __rpow__(self, other, *mod):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other, *mod)

    def __rsub__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __rtruediv__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __sizeof__(self):  # PYCHOK no cover
        '''Size of this instance in C{bytes}.
        '''
        return sum(p.__sizeof__() for p in (self._ps + [self._ps]))

    def __str__(self):
        return Fmt.SQUARE(self.named3, len(self))

    def __sub__(self, other):
        '''Difference of this and an other instance or a scalar.

           @arg other: L{Fsum} instance or C{scalar}.

           @return: The difference (L{Fsum}).

           @see: Method L{Fsum.__isub__}.
        '''
        f = self.fcopy(name=self.__sub__.__name__)
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
        p = (fsum(p for p in f._ps if p > 0) +
             fsum(p for p in f._ps if p < 0)) if f._ps else _0_0
        return s, p

    def _Error(self, op, other, Error=_TypeError, **txt):  # PYCHOK no cover
        '''(INTERNAL) Return an operation error.
        '''
        return Error(_SPACE_(self, op, repr(other)), **txt)

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
        # f._fsum2_ = self._fsum2_
        # f._n = self._n
        f._ps = list(self._ps)  # separate copy
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
                ps[:] = [-p for p in ps]
                if self._fsum2_:
                    self._fsum2_ = -self._fsum2_
        # assert self._ps is ps
        return self

    def fsub(self, xs):
        '''Subtract more I{scalar} values from an iterable.

           @arg xs: Iterable, list, tuple (C{scalar}s).

           @return: This instance (L{Fsum}).

           @see: Method L{Fsum.fadd}.
        '''
        if xs:
            self.fadd(-s for s in xs)
        return self

    def fsub_(self, *xs):
        '''Subtract more I{scalar} values from positional arguments.

           @arg xs: Values to subtract (C{scalar}s), all positional.

           @return: This instance (L{Fsum}).

           @see: Method L{Fsum.fadd}.
        '''
        return self.fsub(xs)

    def fsum(self, xs=None):
        '''Accumulate more I{scalar} values from an iterable and sum all.

           @kwarg xs: Iterable, list, tuple (C{scalar}s) or a single
                      C{scalar} of L{Fsum} instance.

           @return: Accurate, running sum (C{float}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{xs}} value.

           @raise ValueError: Invalid or non-finite B{C{xs}} value.

           @note: Accumulation can continue after summation.
        '''
        if isinstance(xs, _Scalar):
            self += xs
        elif xs:
            self.fadd(xs)

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
        '''Accumulate more I{scalar} values from positional arguments and sum all.

           @arg xs: Values to add (C{scalar}s), all positional.

           @return: Accurate, running sum (C{float}).

           @see: Method L{Fsum.fsum}.

           @note: Accumulation can continue after summation.
        '''
        return self.fsum(xs)

    def fsum2_(self, *xs):
        '''Accumulate more I{scalar} values from positional arguments,
           sum all and provide the sum and delta.

           @arg xs: Values to add (C{scalar}s), all positional.

           @return: 2-Tuple C{(sum, delta)} with the accurate,
                    running C{sum} and the C{delta} with the
                    previous running C{sum}, both (C{float}).

           @see: Method L{Fsum.fsum_}.

           @note: Accumulation can continue after summation.
        '''
        p = self._fsum2_
        if p is None:
            p = self.fsum()
        s = self.fsum(xs)  # if xs else self._fsum2_
        return s, s - p

    @property_RO
    def imag(self):
        '''Return the imaginary part of this instance.
        '''
        return _0_0

    def is_integer(self):
        '''Return C{True} if this instance is an integer.
        '''
        s, p = self._cmp2()
        return s.is_integer() and p == -s

    @property_RO
    def real(self):
        '''Return the real part of this instance.
        '''
        return float(self)


_Float = Fsum, float  # in .fstats
try:  # XXX basics._Ints is ABCMeta
    _Scalar = _Float + (int, long)
except NameError:  # Python 3+
    _Scalar = _Float + (int,)


class Fdot(Fsum):
    '''Precision dot product.
    '''
    def __init__(self, a, *b):
        '''New L{Fdot} precision dot product M{sum(a[i] * b[i]
           for i=0..len(a))}.

           @arg a: List, sequence, tuple, etc. (C{scalar}s).
           @arg b: All positional arguments (C{scalar}s).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise LenError: Unequal C{len(B{a})} and C{len(B{b})}.

           @see: Function L{fdot} and method L{Fsum.fadd}.
        '''
        Fsum.__init__(self)
        self.fadd(_map_a_x_b(a, b))


class Fhorner(Fsum):
    '''Precision polynomial evaluation using the Horner form.
    '''
    def __init__(self, x, *cs):
        '''New L{Fhorner} evaluation of the polynomial
           M{sum(cs[i] * x**i for i=0..len(cs))}.

           @arg x: Polynomial argument (C{scalar}).
           @arg cs: Polynomial coeffients (C{scalar}[]).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{x}}.

           @raise ValueError: Non-finite B{C{x}}.

           @see: Function L{fhorner} and methods L{Fsum.fadd} and L{Fsum.fmul}.
        '''
        Fsum.__init__(self, *cs[-1:])
        if len(cs) > 1:
            x  = _2float(x=x)
            a_ =  self.fadd_
            ps =  self._ps
            for c in reversed(cs[:-1]):  # multiply-accumulate
                ps[:] = [p * x for p in ps]
                a_(c)
            # assert self._ps is ps


class Fpolynomial(Fsum):
    '''Precision polynomial evaluation.
    '''
    def __init__(self, x, *cs):
        '''New L{Fpolynomial} evaluation of the polynomial
           M{sum(cs[i] * x**i for i=0..len(cs))}.

           @arg x: Polynomial argument (C{scalar}).
           @arg cs: Polynomial coeffients (C{scalar}[]).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{x}}.

           @raise ValueError: Non-finite B{C{x}}.

           @see: Function L{fpolynomial} and method L{Fsum.fadd}.
        '''
        Fsum.__init__(self, *cs[:1])
        n = len(cs) - 1
        if n > 0:
            self.fadd(_map_a_x_b(cs[1:], fpowers(x, n)))


def cbrt(x3):
    '''Compute the cube root M{x3**(1/3)}.

       @arg x3: Value (C{scalar}).

       @return: Cubic root (C{float}).

       @see: Functions L{cbrt2} and L{sqrt3}.
    '''
    # <https://archive.lib.MSU.edu/crcmath/math/math/r/r021.htm>
    # simpler and more accurate than Ken Turkowski's CubeRoot, see
    # <https://People.FreeBSD.org/~lstewart/references/apple_tr_kt32_cuberoot.pdf>
    return copysign0(pow(abs(x3), _1_3rd), x3)


def cbrt2(x3):
    '''Compute the cube root I{squared} M{x3**(2/3)}.

       @arg x3: Value (C{scalar}).

       @return: Cube root I{squared} (C{float}).

       @see: Functions L{cbrt} and L{sqrt3}.
    '''
    return pow(abs(x3), _2_3rd)  # XXX pow(abs(x3), _1_3rd)**2


def euclid(x, y):
    '''I{Appoximate} the norm M{sqrt(x**2 + y**2)} by
       M{max(abs(x), abs(y)) + min(abs(x), abs(y)) * 0.4142...}.

       @arg x: X component (C{scalar}).
       @arg y: Y component (C{scalar}).

       @return: Appoximate norm (C{float}).

       @see: Function L{euclid_}.
    '''
    x, y = abs(x), abs(y)
    if y > x:
        x, y = y, x
    return x + y * _0_4142  # XXX * _0_5 before 20.10.02


def euclid_(*xs):
    '''I{Appoximate} the norm M{sqrt(sum(x**2 for x in xs))}
       by cascaded L{euclid}.

       @arg xs: X arguments, positional (C{scalar}[]).

       @return: Appoximate norm (C{float}).

       @see: Function L{euclid}.
    '''
    e = _0_0
    for x in sorted(map(abs, xs)):  # XXX not reverse=True
        # e = euclid(x, e)
        if x > e:
            e, x = x, e
        if x:
            e += x * _0_4142
    return e


def facos1(x):
    '''Fast approximation of L{pygeodesy.acos1}C{(B{x})}.

       @see: U{ShaderFastLibs.h<https://GitHub.com/michaldrobot/
             ShaderFastLibs/blob/master/ShaderFastMathLib.h>}.
    '''
    a = abs(x)
    if a < EPS0:
        r = PI_2
    elif a < EPS1:
        H = Fhorner(-a, 1.5707288, 0.2121144, 0.0742610, 0.0187293)
        r = H.fmul(sqrt(_1_0 - a)).fsum()
        if x < 0:
            r = PI - r
    else:
        r = PI if x < 0 else _0_0
    return r


def fasin1(x):  # PYCHOK no cover
    '''Fast approximation of L{pygeodesy.asin1}C{(B{x})}.

       @see: L{facos1}.
    '''
    return PI_2 - facos1(x)


def fatan(x):
    '''Fast approximation of C{atan(B{x})}.
    '''
    a = abs(x)
    if a < _1_0:
        r = fatan1(a) if a else _0_0
    elif a > _1_0:
        r = PI_2 - fatan1(_1_0 / a)  # == fatan2(a, _1_0)
    else:
        r = PI_4
    return (-r) if x < 0 else r  # copysign0(r, x)


def fatan1(x):
    '''Fast approximation of C{atan(B{x})} for C{0 <= B{x} <= 1}, I{unchecked}.

       @see: U{ShaderFastLibs.h<https://GitHub.com/michaldrobot/
             ShaderFastLibs/blob/master/ShaderFastMathLib.h>} and
             U{Efficient approximations for the arctangent function
             <http://www-Labs.IRO.UMontreal.Ca/~mignotte/IFT2425/Documents/
             EfficientApproximationArctgFunction.pdf>}, IEEE Signal
             Processing Magazine, 111, May 2006.
    '''
    # Eq (9): PI_4 * x - x * (x - 1) * (0.2447 + 0.0663 * x**2)
    # == x * (1.0300982 + x * (-0.2447 + x * 0.0663 * (1 - x)))
    return x * fhorner(x, 1.0300982, -0.2447, 0.0663, -0.0663)


def fatan2(y, x):
    '''Fast approximation of C{atan2(B{y}, B{x})}.

       @see: U{fastApproximateAtan(x, y)<https://GitHub.com/CesiumGS/cesium/blob/
             master/Source/Shaders/Builtin/Functions/fastApproximateAtan.glsl>}
             and L{fatan1}.
    '''
    b, a = abs(y), abs(x)
    if a < b:
        r = (PI_2 - fatan1(a / b)) if a else PI_2
    elif b < a:
        r = fatan1(b / a) if b else _0_0
    elif a:  # == b != 0
        r = PI_4
    else:  # a == b == 0
        return _0_0
    if x < 0:
        r = PI - r
    return (-r) if y < 0 else r  # copysign0(r, y)


def favg(v1, v2, f=_0_5):
    '''Return the average of two values.

       @arg v1: One value (C{scalar}).
       @arg v2: Other value (C{scalar}).
       @kwarg f: Optional fraction (C{float}).

       @return: M{v1 + f * (v2 - v1)} (C{float}).
    '''
#      @raise ValueError: Fraction out of range.
#   '''
#   if not 0 <= f <= 1:  # XXX restrict fraction?
#       raise _ValueError(fraction=f)
    # v1 + f * (v2 - v1) == v1 * (1 - f) + v2 * f
    return fsum1_(v1, -f * v1, f * v2)


def fdot(a, *b):
    '''Return the precision dot product M{sum(a[i] * b[i] for
       i=0..len(a))}.

       @arg a: List, sequence, tuple, etc. (C{scalar}s).
       @arg b: All positional arguments (C{scalar}s).

       @return: Dot product (C{float}).

       @raise LenError: Unequal C{len(B{a})} and C{len(B{b})}.

       @see: Class L{Fdot} and U{Algorithm 5.10 B{DotK}
             <https://www.TUHH.De/ti3/paper/rump/OgRuOi05.pdf>}.
    '''
    return fsum(_map_a_x_b(a, b))


def fdot3(a, b, c, start=0):
    '''Return the precision dot product M{start +
       sum(a[i] * b[i] * c[i] for i=0..len(a))}.

       @arg a: List, sequence, tuple, etc. (C{scalar}[]).
       @arg b: List, sequence, tuple, etc. (C{scalar}[]).
       @arg c: List, sequence, tuple, etc. (C{scalar}[]).
       @kwarg start: Optional bias (C{scalar}).

       @return: Dot product (C{float}).

       @raise LenError: Unequal C{len(B{a})}, C{len(B{b})}
                        and/or C{len(B{c})}.

       @raise OverflowError: Partial C{2sum} overflow.
    '''
    def _mul3(a, b, c):  # map function
        return a * b * c

    def _muly(a, b, c, start):
        yield start
        for abc in map(_mul3, a, b, c):
            yield abc

    if not len(a) == len(b) == len(c):
        raise LenError(fdot3, a=len(a), b=len(b), c=len(c))

    return fsum(_muly(a, b, c, start) if start else map(_mul3, a, b, c))


def fhorner(x, *cs):
    '''Evaluate the polynomial M{sum(cs[i] * x**i for
       i=0..len(cs))} using the Horner form.

       @arg x: Polynomial argument (C{scalar}).
       @arg cs: Polynomial coeffients (C{scalar}[]).

       @return: Horner value (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise TypeError: Non-scalar B{C{x}}.

       @raise ValueError: No B{C{cs}} coefficients or B{C{x}} is not finite.

       @see: Function L{fpolynomial} and class L{Fhorner}.
    '''
    h = Fhorner(x, *cs)
    return h.fsum()


def fidw(xs, ds, beta=2):
    '''Interpolate using using U{Inverse Distance Weighting
       <https://WikiPedia.org/wiki/Inverse_distance_weighting>} (IDW).

       @arg xs: Known values (C{scalar}[]).
       @arg ds: Non-negative distances (C{scalar}[]).
       @kwarg beta: Inverse distance power (C{int}, 0, 1, 2, or 3).

       @return: Interpolated value C{x} (C{float}).

       @raise LenError: Unequal or zero C{len(B{ds})} and C{len(B{xs})}.

       @raise ValueError: Invalid B{C{beta}}, negative B{C{ds}} value,
                          weighted B{C{ds}} below L{EPS}.

       @note: Using C{B{beta}=0} returns the mean of B{C{xs}}.
    '''
    n, xs = len2(xs)
    d, ds = len2(ds)
    if n != d or n < 1:
        raise LenError(fidw, xs=n, ds=d)

    d, x = min(zip(ds, xs))
    if d > EPS0 and n > 1:
        b = -Int_(beta=beta, low=0, high=3)
        if b < 0:
            ds = tuple(d**b for d in ds)
            d  = fsum(ds)
            if isnear0(d):  # PYCHOK no cover
                n = Fmt.PAREN(fsum='ds')
                raise _ValueError(n, d, txt=_singular_)
            x = fdot(xs, *ds) / d
        else:  # b == 0
            x = fsum(xs) / n  # fmean(xs)
    elif d < 0:  # PYCHOK no cover
        n = Fmt.SQUARE(ds=ds.index(d))
        raise _ValueError(n, d, txt=_negative_)
    return x


def fmean(xs):
    '''Compute the accurate mean M{sum(xs[i] for
       i=0..len(xs)) / len(xs)}.

       @arg xs: Values (C{scalar}s).

       @return: Mean value (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise ValueError: No B{C{xs}} values.
    '''
    n, xs = len2(xs)
    if n > 0:
        return fsum(xs) / n
    raise _ValueError(xs=xs)


def fmean_(*xs):
    '''Compute the accurate mean M{sum(xs[i] for
       i=0..len(xs)) / len(xs)}.

       @arg xs: Values (C{scalar}s).

       @return: Mean value (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise ValueError: No B{C{xs}} values.
    '''
    return fmean(xs)


def fpolynomial(x, *cs):
    '''Evaluate the polynomial M{sum(cs[i] * x**i for
       i=0..len(cs))}.

       @arg x: Polynomial argument (C{scalar}).
       @arg cs: Polynomial coeffients (C{scalar}[]).

       @return: Polynomial value (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise TypeError: Non-scalar B{C{x}}.

       @raise ValueError: No B{C{cs}} coefficients or B{C{x}} is not finite.

       @see: Function L{fhorner} and class L{Fpolynomial}.
    '''
    p = Fpolynomial(x, *cs)
    return p.fsum()


def fpowers(x, n, alts=0):
    '''Return a series of powers M{[x**i for i=1..n]}.

       @arg x: Value (C{scalar}).
       @arg n: Highest exponent (C{int}).
       @kwarg alts: Only alternating powers, starting with
                    this exponent (C{int}).

       @return: Powers of B{C{x}} (C{float}s or C{int}s).

       @raise TypeError: Non-scalar B{C{x}} or B{C{n}} not C{int}.

       @raise ValueError: Non-finite B{C{x}} or non-positive B{C{n}}.
    '''
    if not isint(n):
        raise _IsnotError(int.__name__, n=n)
    elif n < 1:
        raise _ValueError(n=n)

    p  = t = x if isint(x) else _2float(x=x)
    ps = [p]
    a_ = ps.append
    for _ in range(1, n):
        p *= t
        a_(p)

    if alts > 0:  # x**2, x**4, ...
        # ps[alts-1::2] chokes PyChecker
        ps = ps[slice(alts-1, None, 2)]

    return ps


try:
    from math import prod as fprod  # Python 3.8
except ImportError:

    def fprod(iterable, start=_1_0):
        '''Iterable product, like C{math.prod} or C{numpy.prod}.

           @arg iterable: Terms to be multiplied (C{scalar}[]).
           @kwarg start: Initial term, also the value returned
                         for an empty iterable (C{scalar}).

           @return: The product (C{float}).

           @see: U{NumPy.prod<https://docs.SciPy.org/doc/
                 numpy/reference/generated/numpy.prod.html>}.
        '''
        return freduce(_mul, iterable, start)


def frange(start, number, step=1):
    '''Generate a range of C{float}s.

       @arg start: First value (C{float}).
       @arg number: The number of C{float}s to generate (C{int}).
       @kwarg step: Increment value (C{float}).

       @return: A generator (C{float}s).

       @see: U{NumPy.prod<https://docs.SciPy.org/doc/
             numpy/reference/generated/numpy.arange.html>}.
    '''
    if not isint(number):
        raise _IsnotError(int.__name__, number=number)
    for i in range(number):
        yield start + i * step


try:
    from functools import reduce as freduce
except ImportError:
    try:
        freduce = reduce  # PYCHOK expected
    except NameError:  # Python 3+

        def freduce(f, iterable, *start):
            '''For missing C{functools.reduce}.
            '''
            if start:
                r = v = start[0]
            else:
                r, v = 0, MISSING
            for v in iterable:
                r = f(r, v)
            if v is MISSING:
                raise _TypeError(iterable=(), start=MISSING)
            return r

try:
    from math import fsum  # precision IEEE-754 sum, Python 2.6+

    # make sure fsum works as expected (XXX check
    # float.__getformat__('float')[:4] == 'IEEE'?)
    if fsum((1, 1e101, 1, -1e101)) != 2:  # PYCHOK no cover
        del fsum  # nope, remove fsum ...
        raise ImportError  # ... use fsum below

except ImportError:

    def fsum(iterable):
        '''Precision summation similar to standard Python function C{math.fsum}.

           Exception and I{non-finite} handling differs from C{math.fsum}.

           @arg iterable: Values to be added (C{scalar}[]).

           @return: Accurate C{sum} (C{float}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{iterable}} value.

           @raise ValueError: Invalid or non-finite B{C{iterable}} value.

           @see: Class L{Fsum}.
        '''
        f = Fsum(name=fsum.__name__)
        return f.fsum(iterable)


def fsum1(iterable):
    '''Precision summation, primed with C{1.0}.

       @arg iterable: Values to be added (C{scalar}[]).

       @return: Accurate L{fsum} (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise TypeError: Non-scalar B{C{iterable}} value.

       @raise ValueError: Invalid or non-finite B{C{iterable}} value.
    '''
    def _xs(iterable):
        yield _1_0
        for x in iterable:
            yield x
        yield _N_1_0

    return fsum(_xs(iterable))


def fsum_(*xs):
    '''Precision summation of all positional arguments.

       @arg xs: Values to be added (C{scalar}s).

       @return: Accurate L{fsum} (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise TypeError: Non-scalar B{C{xs}} value.

       @raise ValueError: Invalid or non-finite B{C{xs}} value.
    '''
    return fsum(map(float, xs))


def fsum1_(*xs):
    '''Precision summation of a few arguments, primed with C{1.0}.

       @arg xs: Values to be added (C{scalar}s).

       @return: Accurate L{fsum} (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise TypeError: Non-scalar B{C{xs}} value.

       @raise ValueError: Invalid or non-finite B{C{xs}} value.
    '''
    return fsum1(xs)


if _sys_version_info2 < (3, 8):

    from math import hypot  # OK in Python 3.7-

    def hypot_(*xs):
        '''Compute the norm M{sqrt(sum(x**2 for x in xs))}.

           Similar to Python 3.8+ n-dimension U{math.hypot
           <https://docs.Python.org/3.8/library/math.html#math.hypot>},
           but exceptions, C{nan} and C{infinite} values are
           handled differently.

           @arg xs: X arguments, positional (C{scalar}[]).

           @return: Norm (C{float}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise ValueError: Invalid or no B{C{xs}} values.

           @note: The Python 3.8+ Euclidian distance U{math.dist
                  <https://docs.Python.org/3.8/library/math.html#math.dist>}
                  between 2 I{n}-dimensional points I{p1} and I{p2} can be
                  computed as M{hypot_(*((c1 - c2) for c1, c2 in zip(p1, p2)))},
                  provided I{p1} and I{p2} have the same, non-zero length I{n}.
        '''
        h, x2 = _h_x2(xs)
        return (h * sqrt(x2)) if x2 else _0_0

elif _sys_version_info2 < (3, 10):

    # In Python 3.8 and 3.9 C{math.hypot} is inaccurate, see
    # agdhruv <https://GitHub.com/geopy/geopy/issues/466>,
    # cffk <https://Bugs.Python.org/issue43088> and module
    # geomath.py <https://PyPI.org/project/geographiclib/1.52>

    def hypot(x, y):
        '''Compute the norm M{sqrt(x**2 + y**2)}.

           @arg x: Argument (C{scalar}).
           @arg y: Argument (C{scalar}).

           @return: C{sqrt(B{x}**2 + B{y}**2)} (C{float}).
        '''
        if x:
            h = sqrt(fsum1_(x**2, y**2)) if y else abs(x)
        elif y:
            h = abs(y)
        else:
            h = _0_0
        return h

    from math import hypot as hypot_  # PYCHOK in Python 3.8 and 3.9
else:
    from math import hypot  # PYCHOK in Python 3.10+
    hypot_ = hypot


def _h_x2(xs):
    '''(INTERNAL) Helper for L{hypot_} and L{hypot2_}.
    '''
    def _x2s(xs, h):
        yield _1_0
        for x in xs:
            if x:
                yield (x / h)**2
        yield _N_1_0

    if xs:
        n, xs = len2(xs)
        if n > 0:
            h  = float(max(map(abs, xs)))
            x2 = fsum(_x2s(xs, h)) if h else _0_0
            return h, x2

    raise _ValueError(xs=xs, txt=_too_(_few_))


def hypot1(x):
    '''Compute the norm M{sqrt(1 + x**2)}.

       @arg x: Argument (C{scalar}).

       @return: Norm (C{float}).
    '''
    return hypot(_1_0, x) if x else _1_0


def hypot2(x, y):
    '''Compute the I{squared} norm M{x**2 + y**2}.

       @arg x: Argument (C{scalar}).
       @arg y: Argument (C{scalar}).

       @return: C{B{x}**2 + B{y}**2} (C{float}).
    '''
    if x:
        x *= x
        h2 = fsum1_(x, y**2) if y else x
    elif y:
        h2 = y**2
    else:
        h2 = _0_0
    return h2


def hypot2_(*xs):
    '''Compute the I{squared} norm C{sum(x**2 for x in B{xs})}.

       @arg xs: X arguments, positional (C{scalar}[]).

       @return: Squared norm (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise ValueError: Invalid or no B{C{xs}} value.

       @see: Function L{hypot_}.
    '''
    h, x2 = _h_x2(xs)
    return (h**2 * x2) if x2 else _0_0


def _map_a_x_b(a, b):
    '''(INTERNAL) Yield B{C{a * b}}.
    '''
    n = len(b)
    if len(a) != n:
        raise LenError(fdot, a=len(a), b=n)
    return map(_mul, a, b) if n > 3 else _map_a_x_b1(a, b)


def _map_a_x_b1(a, b):
    '''(INTERNAL) Yield B{C{a * b}}, primed with C{1.0}.
    '''
    yield _1_0
    for ab in map(_mul, a, b):
        if ab:
            yield ab
    yield _N_1_0


def norm2(x, y):
    '''Normalize a 2-dimensional vector.

       @arg x: X component (C{scalar}).
       @arg y: Y component (C{scalar}).

       @return: 2-Tuple C{(x, y)}, normalized.

       @raise ValueError: Invalid B{C{x}} or B{C{y}}
              or zero norm.
    '''
    h = hypot(x, y)
    try:
        return x / h, y / h
    except (TypeError, ValueError) as X:
        raise _ValueError(x=x, y=y, h=h, txt=str(X))


def norm_(*xs):
    '''Normalize all n-dimensional vector components.

       @arg xs: The component (C{scalar}[]).

       @return: Yield each component, normalized.

       @raise ValueError: Invalid or insufficent B{C{xs}}
              or zero norm.
    '''
    h = hypot_(*xs)
    try:
        for i, x in enumerate(xs):
            yield x / h
    except (TypeError, ValueError) as X:
        raise _ValueError(Fmt.SQUARE(xs=i), x, _h_, h, txt=str(X))


def sqrt0(x2):
    '''Compute the square root iff C{B{x2} >} L{EPS02}.

       @arg x2: Value (C{scalar}).

       @return: Square root (C{float}) or C{0.0}.

       @note: Any C{B{x2} <} L{EPS02} I{including} C{B{x2} < 0}
              returns C{0.0}.
    '''
    return sqrt(x2) if x2 > EPS02 else (_0_0 if x2 < EPS02 else EPS0)


def sqrt3(x2):
    '''Compute the square root, I{cubed} M{sqrt(x)**3} or M{sqrt(x**3)}.

       @arg x2: Value (C{scalar}).

       @return: Cubed square root (C{float}).

       @raise ValueError: Negative B{C{x2}}.

       @see: Functions L{cbrt} and L{cbrt2}.
    '''
    if x2 < 0:
        raise _ValueError(x2=x2, txt=_negative_)
    return pow(x2, _1_5) if x2 else _0_0

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
