
# -*- coding: utf-8 -*-

u'''Class L{Fsum} for precision floating point summation and I{running}
summation, based on respectively similar to Python's C{math.fsum}.

Generally, an L{Fsum} instance is considered a C{float} plus a small or
zero C{residual} value, see property L{Fsum.residual}.  However, there
are several C{integer} L{Fsum} cases, for example the result of C{ceil},
C{floor}, C{Fsum.__floordiv__} and methods L{Fsum.fint} and L{Fsum.fint2}.

Also, L{Fsum} methods L{Fsum.pow}, L{Fsum.__ipow__}, L{Fsum.__pow__} and
L{Fsum.__rpow__} return a (very long) C{int} if invoked with optional
argument C{mod} set to C{None}.  The C{residual} of an C{integer} L{Fsum}
may be anywhere between C{-1.0} and C{+1.0}, including C{INT0} if considered
to be I{exact}.

Set env variable C{PYGEODESY_FSUM_RESIDUAL} to any non-empty string to throw
a L{ResidualError} for division or exponention by an L{Fsum} instance with
a non-zero C{residual}, see methods L{Fsum.RESIDUAL}, L{Fsum.pow},
L{Fsum.__ipow__} and L{Fsum.__itruediv__}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import iscomplex, _isfinite, isinf, isint, isnan, \
                             isscalar, map1, neg, signOf, _signOf
from pygeodesy.errors import _OverflowError, _TypeError, _ValueError, \
                             _xError2, _xkwds_get, _ZeroDivisionError
from pygeodesy.interns import INT0, NN, _arg_, _COMMASPACE_, _DASH_, _EQUAL_, \
                             _from_, _iadd_, _not_finite_, _not_scalar_, \
                             _PERCENT_, _PLUS_, _SLASH_, _SPACE_, _STAR_, \
                             _0_0, _1_0, _N_1_0
from pygeodesy.lazily import _ALL_LAZY, _getenv, _sys_version_info2
from pygeodesy.named import _Named, _NamedTuple, _NotImplemented
from pygeodesy.props import deprecated_property_RO, Property_RO, property_RO
from pygeodesy.streprs import Fmt, pairs, unstr
from pygeodesy.units import Float, Int

from math import ceil as _ceil, floor as _floor  # PYCHOK used!

__all__ = _ALL_LAZY.fsums
__version__ = '22.06.07'

_eq_ = _EQUAL_ * 2  # _DEQUAL_
_ge_ = '>='
_gt_ = '>'  # _RANGLE_
_le_ = '<='
_lt_ = '<'  # _LANGLE_
_ne_ = '!='

_add_      = _PLUS_
_floordiv_ = _SLASH_ * 2
_fset_     = _EQUAL_
_integer_  = 'integer'
_mod_      = _PERCENT_
_divmod_   = _floordiv_ + _mod_
_mul_      = _STAR_
_non_zero_ = 'non-zero'
_pow_      = _STAR_ * 2  # _DSTAR_
_residual_ = 'residual'
_sub_      = _DASH_
_truediv_  = _SLASH_

_pos_self  = _1_0.__pos__() is _1_0


def _2even(s, r, p):
    '''(INTERNAL) Half-even rounding.
    '''
    if (r > 0 and p > 0) or \
       (r < 0 and p < 0):  # signs match
        r, p = _2sum(s, p * 2)
        if not p:
            s = r
    return s


def _2float(index=None, **name_value):
    '''(INTERNAL) Raise C{TypeError} or C{ValueError} if not scalar or infinite.
    '''
    n, v = name_value.popitem()  # _xkwds_popitem(name_value)
    try:
        v = float(v)
        if _isfinite(v):
            return v
        E, t = _ValueError, _not_finite_
    except Exception as x:
        E, t = _xError2(x)
    if index is not None:
        n = Fmt.SQUARE(n, index)
    raise E(n, v, txt=t)


def _floats(floats=False, **unused):
    '''(INERNAL) Unravel the optional C{floats} keyword argument for
       functions C{fsum}, C{fsum_}, C{fsum1}, and C{fsum1_} below.
    '''
    return floats


def _2floats(xs, origin=0, primed=False, sub=False, floats=False):
    '''(INTERNAL) Yield all B{C{xs}} as C{float}s.
    '''
    _2f = _2float
    if primed:
        yield _1_0
    i = origin
    for x in xs:
        if floats:
            yield x
        elif isinstance(x, Fsum):
            ps = x._ps
            if ps:
                if sub:
                    ps = map(neg, ps)
                for x in ps:
                    yield x
        else:
            x = _2f(index=i, xs=x)
            if x:
                yield (-x) if sub else x
        i += 1
    if primed:
        yield _N_1_0


def _2Fsum(other, name=NN):
    '''(INTERNAL) Return B{C{other}} as an L{Fsum} instance.
    '''
    return other.copy(name=name) if isinstance(other, Fsum) else \
           Fsum(name=name)._fset(_2float(other=other), asis=True)


def _2scalar(other, raiser=False):
    '''(INTERNAL) Return B{C{other}} as C{int}, C{float} or C{as-is}.
    '''
    if isinstance(other, Fsum):
        s, r = other._fint2
        if r:
            s, r = other._fprs2
            if r:  # PYCHOK no cover
                if raiser:
                    raise ValueError(_2stresidual(_non_zero_, r))
                s = other  # L{Fsum} as-is
    else:
        s = other  # C{type} as-is
        if isint(s, both=True):
            s = int(s)
    return s


def _2strcomplex(s, *args):
    '''(INTERNAL) C{Complex} 2- or 3-arg C{pow} error C{str}.
    '''
    c =  iscomplex.__name__[2:]
    n = _DASH_(len(args), _arg_)
    t = _SPACE_(c, s, _from_, n, pow.__name__)
    return unstr(t, *args)


def _2stresidual(prefix, residual, **name_values):
    '''(INTERNAL) Residual error C{str}.
    '''
    def _fmt(s):
        return str(s) if isint(s) and abs(s) < 1e9 else (
               Fmt.g(s, prec=9) if isscalar(s) else repr(s))

    p = _SPACE_(prefix, _residual_)
    t =  Fmt.PARENSPACED(p, _fmt(residual))
    for n, v in name_values.items():
        n =  Fmt.PARENSPACED(n, _fmt(v))
        t = _COMMASPACE_(t, n)
    return t


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
    '''Precision floating point I{running} summation similar to standard Python's C{math.fsum}.

       Unlike C{math.fsum}, this class accumulates values and provides I{intermediate}
       precision floating point summation.  Accumulation may continue after I{intermediate}
       summuation, aka I{running} summation.

       @note: Accumulated values may be other L{Fsum} or C{scalar} instances with C{scalar}
              meaning type C{float}, C{int} or any C{type} convertible into a single-instance
              C{float}.

       @note: Handling of exceptions and of values C{inf}, C{INF}, C{nan} and C{NAN} differs
              from standard Python's C{math.fsum}.

       @see: U{Hettinger<https://GitHub.com/ActiveState/code/blob/master/recipes/Python/
             393090_Binary_floating_point_summatiaccurate_full/recipe-393090.py>}, U{Kahan
             <https://WikiPedia.org/wiki/Kahan_summation_algorithm>}, U{Klein
             <https://Link.Springer.com/article/10.1007/s00607-005-0139-x>}, Python 2.6+
             file I{Modules/mathmodule.c} and the issue log U{Full precision summation
             <https://Bugs.Python.org/issue2819>}.
    '''
    _math_fsum = None
    _n         = 0
#   _ps        = []  # partial sums
    _RESIDUAL  = bool(_getenv('PYGEODESY_FSUM_RESIDUAL', NN))

    def __init__(self, *xs, **name_NN):
        '''New L{Fsum} for precision floating point I{running} summation.

           @arg xs: No, one or more initial values (C{scalar} or
                    L{Fsum} instances).
           @kwarg name_NN: Optional name (C{str}).

           @see: Method L{Fsum.fadd}.
        '''
#       self._n  = 0
        self._ps = []  # [_0_0], see L{Fsum._fprs}
        if name_NN:
            self.name = _xkwds_get(name_NN, name=NN)
        if len(xs) > 1:
            self._facc(_2floats(xs, origin=1))
        elif xs:  # len(xs) == 1
            self._ps.append(_2float(xs=xs[0]))
            self._n = 1

    def __abs__(self):
        '''Return this instance' absolute value as an L{Fsum}.
        '''
        f = self.copy(name=self.__abs__.__name__)
        return f._fneg() if f < 0 else f

    def __add__(self, other):
        '''Return the sum C{B{self} + B{other}} as an L{Fsum}.

           @arg other: An L{Fsum} or C{scalar}.

           @return: The sum (L{Fsum}).

           @see: Method L{Fsum.__iadd__}.
        '''
        f = self.copy(name=self.__add__.__name__)
        return f._fadd(other, _add_)

    def __bool__(self):  # PYCHOK not special in Python 2-
        '''Return C{True} if this instance is non-zero.
        '''
        s, r = self._fprs2
        return bool(s or r)

    def __ceil__(self):  # PYCHOK not special in Python 2-
        '''Return this instance' C{math.ceil} as C{int} or C{float}.

           @return: An C{int} in Python 3+, but C{float} in Python 2-.

           @see: Methods L{Fsum.__floor__} and property L{Fsum.ceil}.
        '''
        return self.ceil

    def __divmod__(self, other):
        '''Return C{divmod(B{self}, B{other})} as 2-tuple C{(quotient,
           remainder)}, an C{int} in Python 3+ or C{float} in Python 2-
           and an L{Fsum}.

           @arg other: An L{Fsum} or C{scalar} modulus.

           @see: Method L{Fsum.__itruediv__}.
        '''
        f = self.copy(name=self.__divmod__.__name__)
        return f._fdivmod(other, _divmod_)

    def __eq__(self, other):
        '''Compare this with an other instance or scalar.
        '''
        s, r = self._cmp2(other, _eq_)
        return not bool(s or r)

    def __float__(self):
        '''Return this instance' current precision running sum as C{float}.

           @see: Methods L{Fsum.fsum} and L{Fsum.int_float}.
        '''
        return float(self._fprs)

    def __floor__(self):  # PYCHOK not special in Python 2-
        '''Return this instance' C{math.floor} as C{int} or C{float}.

           @return: An C{int} in Python 3+, but C{float} in Python 2-.

           @see: Methods L{Fsum.__ceil__} and property L{Fsum.floor}.
        '''
        return self.floor

    def __floordiv__(self, other):
        '''Return C{B{self} // B{other}} as an L{Fsum}.

           @arg other: An L{Fsum} or C{scalar} divisor.

           @return: The C{floor} quotient (L{Fsum}).

           @see: Methods L{Fsum.__ifloordiv__}.
        '''
        f = self.copy(name=self.__floordiv__.__name__)
        return f._floordiv(other, _floordiv_)

    def __format__(self, *other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, *other)

    def __ge__(self, other):
        '''Compare this with an other instance or scalar.
        '''
        s, r = self._cmp2(other, _ge_)
        return r >= -s

    def __gt__(self, other):
        '''Compare this with an other instance or scalar.
        '''
        s, r = self._cmp2(other, _gt_)
        return r > -s

    def __hash__(self):  # PYCHOK no cover
        '''Return this instance' C{hash}.
        '''
        return hash(self._ps)  # XXX id(self)?

    def __iadd__(self, other):
        '''Apply C{B{self} += B{other}} to this instance.

           @arg other: An L{Fsum} or C{scalar} instance.

           @return: This instance, updated (L{Fsum}).

           @raise TypeError: Invalid B{C{other}}, not
                             C{scalar} nor L{Fsum}.

           @see: Methods L{Fsum.fadd} and L{Fsum.fadd_}.
        '''
        return self._fadd(other, _iadd_)

    def __ifloordiv__(self, other):
        '''Apply C{B{self} //= B{other}} to this instance.

           @arg other: An L{Fsum} or C{scalar} divisor.

           @return: This instance, updated (L{Fsum}).

           @raise ResidualError: Non-zero residual in B{C{other}}.

           @raise TypeError: Invalid B{C{other}} type.

           @raise ValueError: Invalid or non-finite B{C{other}}.

           @raise ZeroDivisionError: Zero B{C{other}}.

           @see: Methods L{Fsum.__itruediv__}.
        '''
        return self._floordiv(other, _floordiv_ + _fset_)

    def __imatmul__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __imod__(self, other):
        '''Apply C{B{self} %= B{other}} to this instance.

           @arg other: An L{Fsum} or C{scalar} modulus.

           @return: This instance, updated (L{Fsum}).

           @see: Method L{Fsum.__divmod__}.
        '''
        self._fdivmod(other, _mod_ + _fset_)
        return self

    def __imul__(self, other):
        '''Apply C{B{self} *= B{other}} to this instance.

           @arg other: An L{Fsum} or C{scalar} factor.

           @return: This instance, updated (L{Fsum}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Invalid B{C{other}} type.

           @raise ValueError: Invalid or non-finite B{C{other}}.
        '''
        return self._fmul(other, _mul_ + _fset_)

    def __int__(self):
        '''Return this instance as an C{int}.

           @see: Methods L{Fsum.int_float}, L{Fsum.__ceil__}
                 and L{Fsum.__floor__} and properties
                 L{Fsum.ceil} and L{Fsum.floor}.
        '''
        return int(self._fprs)

    def __ipow__(self, other, *mod):  # PYCHOK 2 vs 3 args
        '''Apply C{B{self} **= B{other}} to this instance.

           @arg other: The exponent (L{Fsum} or C{scalar}).
           @arg mod: Optional modulus (C{int} or C{None}) for the
                     3-argument C{pow(B{self}, B{other}, B{mod})}
                     version.

           @return: This instance, updated (L{Fsum}).

           @note: If B{C{mod}} is given, the result will be an C{integer}
                  L{Fsum} in Python 3+ if this instance C{is_integer} or
                  set to C{as_integer} if B{C{mod}} given as C{None}.

           @raise OverflowError: Partial C{2sum} overflow.

           @raise ResidualError: Non-zero residual in B{C{other}} and
                                 env var C{PYGEODESY_FSUM_RESIDUAL}
                                 set or this instance has a non-zero
                                 residual and either B{C{mod}} is
                                 given and non-C{None} or B{C{other}}
                                 is a negative or fractional C{scalar}.

           @raise TypeError: Invalid B{C{other}} type or 3-argument
                             C{pow} invocation failed.

           @raise ValueError: If B{C{other}} is a negative C{scalar}
                              and this instance is C{0} or B{C{other}}
                              is a fractional C{scalar} and this
                              instance is negative or has a non-zero
                              residual or B{C{mod}} is given and C{0}.

           @see: CPython function U{float_pow<https://GitHub.com/
                 python/cpython/blob/main/Objects/floatobject.c>}.
        '''
        return self._fpow(other, _pow_ + _fset_, *mod)

    def __isub__(self, other):
        '''Apply C{B{self} -= B{other}} to this instance.

           @arg other: An L{Fsum} or C{scalar}.

           @return: This instance, updated (L{Fsum}).

           @raise TypeError: Invalid B{C{other}} type.

           @see: Method L{Fsum.fadd}.
        '''
        return self._fsub(other, _sub_ + _fset_)

    def __iter__(self):
        '''Return an C{iter}ator over a C{partials} duplicate.
        '''
        return iter(self.partials)

    def __itruediv__(self, other):
        '''Apply C{B{self} /= B{other}} to this instance.

           @arg other: An L{Fsum} or C{scalar} divisor.

           @return: This instance, updated (L{Fsum}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise ResidualError: Non-zero residual in B{C{other}} and
                                 env var C{PYGEODESY_FSUM_RESIDUAL} set.

           @raise TypeError: Invalid B{C{other}} type.

           @raise ValueError: Invalid or non-finite B{C{other}}.

           @raise ZeroDivisionError: Zero B{C{other}}.

           @see: Method L{Fsum.__ifloordiv__}.
        '''
        return self._ftruediv(other, _truediv_ + _fset_)

    def __le__(self, other):
        '''Compare this with an other instance or scalar.
        '''
        s, r = self._cmp2(other, _le_)
        return r <= -s

    def __len__(self):
        '''Return the I{total} number of values accumulated (C{int}).
        '''
        return self._n

    def __lt__(self, other):
        '''Compare this with an other instance or scalar.
        '''
        s, r = self._cmp2(other, _lt_)
        return r < -s

    def __matmul__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __mod__(self, other):
        '''Return C{B{self} % B{other}} as an L{Fsum}.

           @see: Method L{Fsum.__imod__}.
        '''
        f = self.copy(name=self.__mod__.__name__)
        return f._fdivmod(other, _mod_)[1]

    def __mul__(self, other):
        '''Return C{B{self} * B{other}} as an L{Fsum}.

           @see: Method L{Fsum.__imul__}.
        '''
        f = self.copy(name=self.__mul__.__name__)
        return f._fmul(other, _mul_)

    def __ne__(self, other):
        '''Compare this with an other instance or scalar.
        '''
        s, r = self._cmp2(other, _ne_)
        return bool(s or r)

    def __neg__(self):
        '''Return I{a copy of} this instance, negated.
        '''
        f = self.copy(name=self.__neg__.__name__)
        return f._fneg()

    def __pos__(self):
        '''Return this instance I{as-is}, like C{float.__pos__()}.
        '''
        return self if _pos_self else self.copy(name=self.__pos__.__name__)

    def __pow__(self, other, *mod):  # PYCHOK 2 vs 3 args
        '''Return C{B{self}**B{other}} as an L{Fsum}.

           @see: Method L{Fsum.__ipow__}.
        '''
        f = self.copy(name=self.__pow__.__name__)
        return f._fpow(other, _pow_, *mod)

    def __radd__(self, other):
        '''Return C{B{other} + B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__iadd__}.
        '''
        f = _2Fsum(other, name=self.__radd__.__name__)
        return f._fadd(self, _add_)

    def __rdivmod__(self, other):
        '''Return C{divmod(B{other}, B{self})} as 2-tuple C{(quotient,
           remainder)}.

           @see: Method L{Fsum.__divmod__}.
        '''
        f = _2Fsum(other, name=self.__rdivmod__.__name__)
        return f._fdivmod(self, _divmod_)

#   def __repr__(self):
#       '''Return the default C{repr(this)}.
#       '''
#       return self.toRepr()

    def __rfloordiv__(self, other):
        '''Return C{B{other} // B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__ifloordiv__}.
        '''
        f = _2Fsum(other, name=self.__rfloordiv__.__name__)
        return f._floordiv(self, _floordiv_)

    def __rmatmul__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __rmod__(self, other):
        '''Return C{B{other} % B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__imod__}.
        '''
        f = _2Fsum(other, name=self.__rmod__.__name__)
        return f._fdivmod(self, _mod_)[1]

    def __rmul__(self, other):
        '''Return C{B{other} * B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__imul__}.
        '''
        f = _2Fsum(other, name=self.__rmul__.__name__)
        return f._fmul(self, _mul_)

    def __round__(self, ndigits=None):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, ndigits=ndigits)

    def __rpow__(self, other, *mod):
        '''Return C{B{other}**B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__ipow__}.
        '''
        f = _2Fsum(other, name=self.__rpow__.__name__)
        return f._fpow(self, _pow_, *mod)

    def __rsub__(self, other):
        '''Return C{B{other} - B{self}} as L{Fsum}.

           @see: Method L{Fsum.__isub__}.
        '''
        f = _2Fsum(other, name=self.__rsub__.__name__)
        return f._fsub(self, _sub_)

    def __rtruediv__(self, other):
        '''Return C{B{other} / B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__itruediv__}.
        '''
        f = _2Fsum(other, name=self.__rtruediv__.__name__)
        return f._ftruediv(self, _truediv_)

    def __sizeof__(self):  # PYCHOK not special in Python 2-
        '''Return the current size of this instance in C{bytes}.
        '''
        from sys import getsizeof
        return sum(map1(getsizeof, self._fint2,
                                   self._fint2[0],
                                   self._fint2[1],
                                   self._fprs,
                                   self._fprs2,
                                   self._fprs2.fsum,
                                   self._fprs2.residual,
                                   self._n,
                                   self._ps, *self._ps))

    def __str__(self):
        '''Return the default C{str(self)}.
        '''
        return self.toStr()

    def __sub__(self, other):
        '''Return C{B{self} - B{other}} as an L{Fsum}.

           @arg other: An L{Fsum} or C{scalar}.

           @return: The difference (L{Fsum}).

           @see: Method L{Fsum.__isub__}.
        '''
        f = self.copy(name=self.__sub__.__name__)
        return f._fsub(other, _sub_)

    def __truediv__(self, other):
        '''Return C{B{self} / B{other}} as an L{Fsum}.

           @arg other: An L{Fsum} or C{scalar} divisor.

           @return: The quotient (L{Fsum}).

           @see: Method L{Fsum.__itruediv__}.
        '''
        f = self.copy(name=self.__truediv__.__name__)
        return f._ftruediv(other, _truediv_)

    __trunc__ = __int__

    if _sys_version_info2 < (3, 0):  # PYCHOK no cover
        # <https://docs.Python.org/2/library/operator.html#mapping-operators-to-functions>
        __div__     = __truediv__
        __idiv__    = __itruediv__
        __long__    = __int__
        __nonzero__ = __bool__
        __rdiv__    = __rtruediv__

    def as_integer_ratio(self):
        '''Return this instance as the integer ratio.

           @return: 2-Tuple C{(numerator, denominator)} both
                    C{int} and with positive C{denominator}.

           @see: Standard C{float.as_integer_ratio} in Python 3+.
        '''
        n, r = self._fint2
        if r:
            i, d = r.as_integer_ratio()
            n = (n * d) + i
        else:  # PYCHOK no cover
            d = 1
        return n, d

    @property_RO
    def ceil(self):
        '''Get this instance' C{ceil} value (C{int} in Python 3+, but
           C{float} in Python 2-).

           @note: The C{ceil} takes the C{residual} into account.

           @see: Method L{Fsum.int_float} and properties L{Fsum.floor},
                 L{Fsum.imag} and L{Fsum.real}.
        '''
        s, r = self._fprs2
        c = _ceil(s) + int(r) - 1
        while r > (c - s):  # (s + r) > c
            c += 1
        return c

    def _cmp2(self, other, op):
        '''(INTERNAL) Subtract an B{C{other}} instance or scalar and
           return an L{Fsum2Tuple}C{(fsum, residual)} for comparison
           operator B{C{op}}.
        '''
        if isinstance(other, Fsum):
            f = self._copy()._facc(map(neg, other._ps))
        elif not isscalar(other):
            raise self._TypeError(op, other, txt=_not_scalar_)
        elif self._finite(other, op):
            f = self._copy()._facc_(-other)
        else:
            f = self
        return f._fprs2

    def copy(self, deep=False, name=NN):
        '''Copy this instance, C{shallow} or B{C{deep}}.

           @return: The copy (L{Fsum}).
         '''
        f = _Named.copy(self, deep=deep, name=name)
#       f._update(other=self)
        f._n  = self._n if deep else 1
        f._ps = list(self._ps)  # separate list
        return f

    def _copy(self, _fprs2=False):
        '''(INTERNAL) Fast, un-named copy.
        '''
        f = Fsum()  # NN
        f._n     = self._n
        f._ps[:] = self._ps
        if _fprs2:  # only the ._fprs2 2-tuple
            Fsum._fprs2._update_from(f, self)
        return f

    def divmod(self, other):
        '''Return C{divmod(B{self}, B{other})} as 2-tuple C{(quotient,
           remainder)}.

           @arg other: An L{Fsum} or C{scalar} divisor.

           @return: 2-Tuple C{(quotient, remainder)}, with the C{quotient}
                    an C{int} in Python 3+ or a C{float} in Python 2- and
                    the C{remainder} an L{Fsum} instance.

           @see: Method L{Fsum.__itruediv__}.
        '''
        return self.copy(name=self.divmod.__name__)._fdivmod(other, _divmod_)

    def _Error(self, op, other, Error, **txt):
        '''(INTERNAL) Format an B{C{Error}} for C{{self} B{op} B{other}}.
        '''
        return Error(_SPACE_(self.toRepr(), op, repr(other)), **txt)

    def _facc(self, xs):
        '''(INTERNAL) Accumulate more C{scalar}s.
        '''
        ps, n, _2s = self._ps, 0, _2sum
        for x in xs:  # _iter()
            # assert isscalar(x) and isfinite(x)
            i = 0
            for p in ps:
                x, p = _2s(x, p)
                if p:
                    ps[i] = p
                    i += 1
            ps[i:] = [x]
            n += 1
        # assert self._ps is ps
        if n:
            self._n += n
            self._update()
        return self

    def _facc_(self, *xs):
        '''(INTERNAL) Accumulate all positional C{scalar}s.
        '''
        return self._facc(xs)

    def _facc_up(self):
        '''(INTERNAL) Update the C{partials}, by removing and
           re-adding the final C{partial}.
        '''
        while len(self._ps) > 1:
            p = self._ps.pop()
            if p:
                self._facc_(p)
                self._n -= 1
                break
        else:  # force zap
            self._update()
        return self  # ._fpsqz()

    def fadd(self, xs=()):
        '''Add an iterable of C{scalar} or L{Fsum} instances
           to this instance.

           @arg xs: Iterable, list, tuple, etc. (C{scalar} or
                    L{Fsum} instances).

           @return: This instance (L{Fsum}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: An invalid B{C{xs}} type, not C{scalar}
                             nor L{Fsum}.

           @raise ValueError: Invalid or non-finite B{C{xs}} value.
        '''
        if isinstance(xs, Fsum):
            self._facc(xs._ps)
        elif isscalar(xs):  # for backward compatibility
            self._facc_(_2float(x=xs))  # PYCHOK no cover
        elif xs:
            self._facc(_2floats(xs))
        return self

    def fadd_(self, *xs):
        '''Add all positional C{scalar} or L{Fsum} instances
           to this instance.

           @arg xs: Values to add (C{scalar} or L{Fsum} instances),
                    all positional.

           @return: This instance (L{Fsum}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: An invalid B{C{xs}} type, not C{scalar}
                             nor L{Fsum}.

           @raise ValueError: Invalid or non-finite B{C{xs}} value.
        '''
        return self._facc(_2floats(xs, origin=1))

    def _fadd(self, other, op):
        '''(INTERNAL) Apply C{B{self} += B{other}}.
        '''
        if isinstance(other, Fsum):
            if other._ps:
                if other is self:  # self *= 2
                    self._n    += len(self._ps)  # like ._mul_scalar
                    self._ps[:] = self._ps_x(op, 2)
                    self._facc_up()
                else:
                    self._facc(other._ps)
        elif not isscalar(other):
            raise self._TypeError(op, other, txt=_not_scalar_)
        elif other:
            self._facc_(other)
        return self

    fcopy   =   copy        # for backward compatibility
    fdiv    = __itruediv__  # for backward compatibility
    fdivmod = __divmod__    # for backward compatibility

    def _fdivmod(self, other, op):
        '''(INTERNAL) C{divmod(B{self}, B{other})} as 2-tuple
           (C{int} or C{float}, remainder C{self}).
        '''
        # result mostly follows CPython function U{float_divmod
        # <https://GitHub.com/python/cpython/blob/main/Objects/floatobject.c>},
        # but at least divmod(-3, 2) equals Cpython's result (-2, 1).
        q = self._copy(_fprs2=True)._ftruediv(other, op).floor
        if q:  # == float // other == floor(float / other)
            self -= other * q

        s = signOf(other)  # make signOf(self) == signOf(other)
        if s and self.signOf() == -s:  # PYCHOK no cover
            self += other
            q -= 1

#           t = self.signOf()
#           if t and t != s:
#               from pygeodesy.errors import _AssertionError
#               raise self._Error(op, other, _AssertionError, txt=signOf.__name__)
        return q, self  # q is C{int} in Python 3+, but C{float} in Python 2-

    def _finite(self, other, op=None):
        '''(INTERNAL) Return B{C{other}} if C{finite}.
        '''
        if _isfinite(other):
            return other
        raise ValueError(_not_finite_) if not op else \
              self._ValueError(op, other, txt=_not_finite_)

    def fint(self, raiser=True, name=NN):
        '''Return this instance' current running sum as C{integer}.

           @kwarg raiser: If C{True} throw a L{ResidualError} if the
                          I{integer} residual is non-zero.
           @kwarg name: Optional name (C{str}), overriding C{"fint"}.

           @return: The C{integer} (L{Fsum}).

           @raise ResidualError: Non-zero I{integer} residual.

           @see: Methods L{Fsum.int_float} and L{Fsum.is_integer}.
        '''
        i, r = self._fint2
        if r and raiser:
            t = _2stresidual(_integer_, r)
            raise ResidualError(_integer_, i, txt=t)
        n = name or self.fint.__name__
        return Fsum(name=n)._fset(i, asis=True)

    def fint2(self):
        '''Return this instance' current running sum as C{int} and
           the I{integer} residual.

           @return: L{Fsum2Tuple}C{(fsum, residual)} with C{fsum} an
                    C{int} and I{integer} C{residual} a C{float} or
                    C{INT0} if the C{fsum} is considered I{exact}.
        '''
        return Fsum2Tuple(*self._fint2)

    @Property_RO
    def _fint2(self):
        '''(INTERNAL) Get 2-tuple (C{int}, I{integer} residual).
        '''
        i = int(self._fprs)  # int(self)
        # assert len(self._ps) > 0
        r = _fsum(self._ps1(i)) if len(self._ps) > 1 else (
                 (self._ps[0] - i) or INT0)
        return i, r

    @deprecated_property_RO
    def float_int(self):  # PYCHOK no cover
        '''DEPRECATED, use method C{Fsum.int_float}.'''
        return self.int_float()  # raiser=False

    @property_RO
    def floor(self):
        '''Get this instance' C{floor} (C{int} in Python 3+, but
           C{float} in Python 2-).

           @note: The C{floor} takes the C{residual} into account.

           @see: Method L{Fsum.int_float} and properties L{Fsum.ceil},
                 L{Fsum.imag} and L{Fsum.real}.
        '''
        s, r = self._fprs2
        f = _floor(s) + _floor(r) + 1
        while r < (f - s):  # (s + r) < f
            f -= 1
        return f

#   floordiv = __floordiv__  # for naming consistency

    def _floordiv(self, other, op):  # rather _ffloordiv?
        '''Apply C{B{self} //= B{other}}.
        '''
        q = self._ftruediv(other, op)  # == self
        return self._fset(q.floor, asis=True)  # floor(q)

    fmul = __imul__  # for backward compatibility

    def _fmul(self, other, op):
        '''(INTERNAL) Apply C{B{self} *= B{other}}.
        '''
        a = False
        if isinstance(other, Fsum):
            if len(self._ps) != 1:
                f = self._mul_Fsum(other, op)
            elif len(other._ps) != 1:
                f = other._mul_scalar(self._ps[0], op)
                f._n = self._n
            else:
                f = self._finite(other._ps[0] * self._ps[0], op)
                a = isint(f)
        elif not isscalar(other):
            raise self._TypeError(op, other, txt=_not_scalar_)
        else:
            f = self._mul_scalar(other, op)
        return self._fset(f, asis=a)

    def _fneg(self):
        '''(INTERNAL) Negate this instance.
        '''
        if self._ps:
            self._ps[:] = map(neg, self._ps)
            self._facc_up()
        return self

    def fover(self, over):
        '''Apply C{B{self} /= B{over}} and summate.

           @arg over: An L{Fsum} or C{scalar} denominator.

           @return: Precision running sum (C{float}).

           @see: Methods L{Fsum.fsum} and L{Fsum.__itruediv__}.
        '''
        return float(self.fdiv(over)._fprs)

    fpow = __ipow__  # for backward compatibility

    def _fpow(self, other, op, *mod):
        '''Apply C{B{self} **= B{other}}, optional B{C{mod}} or C{None}.
        '''
        if mod and mod[0] is not None:  # == 3-arg C{pow}
            s = self._pow_3(other, mod[0], op)
        elif mod and mod[0] is None and self.is_integer():
            # return an exact C{int} for C{int}**C{int}
            i, _ = self._fint2   # assert _ == 0
            x = _2scalar(other)  # C{int}, C{float} or other
            if isscalar(x):
                s = self._Fsum(i)._pow_2(x, other, op)
            else:
                s = self._Fsum(i)._fpow(other, op)
        else:  # pow(self, other) == pow(self, other, None)
            p = None
            if isinstance(other, Fsum):
                x, r = other._fprs2
                if r:
                    if self._RESIDUAL:
                        raise self._ResidualError(op, other, r)
                    p =  self._pow_scalar(r, other, op)
                    p = _2scalar(p)  # raiser=False
            elif not isscalar(other):
                raise self._TypeError(op, other, txt=_not_scalar_)
            else:
                x = self._finite(other, op)
            s = self._pow_scalar(x, other, op)
            if p not in (None, _1_0, 1):
                s *= p  # each C{scalar} or L{Fsum}
        return self._fset(s, asis=isint(s))

    @Property_RO
    def _fprs(self):
        '''(INTERNAL) Get and cache this instance' precision
           running sum (C{float} or C{int}), ignoring C{residual}.

           @note: The precision running C{fsum} after a C{//=} or
                  C{//} C{floor} division is C{int} in Python 3+.
        '''
        ps = self._ps
        i = len(ps) - 1
        if i < 0:
            s = _0_0
            ps[:] = [s]
        else:
            s, _2s = ps[i], _2sum
            while i > 0:
                i -= 1
                s, p = _2s(s, ps[i])
                ps[i:] = [s]
                if p:  # sum(ps) became inexact
                    if s:
                        ps.append(p)
                        if i > 0:  # half-even round if signs match
                            s = _2even(s, ps[i-1], p)
                        break
                    else:  # PYCHOK no cover
                        ps[i] = s = p  # replace and continue
        # assert self._ps is ps
        # assert Fsum._fprs2.name not in self.__dict__
        return s

    @Property_RO
    def _fprs2(self):
        '''(INTERNAL) Get and cache this instance' precision
           running sum and residual (L{Fsum2Tuple}).
        '''
        s =  self._fprs
        r = _fsum(self._ps1(s)) if len(self._ps) > 1 else INT0
        return Fsum2Tuple(s, r)

#   def _fpsqz(self):
#       '''(INTERNAL) Compress, squeeze the C{partials}.
#       '''
#       if len(self._ps) > 2:
#           _ = self._fprs
#       return self

    def _fset(self, other, asis=False, n=1):
        '''(INTERNAL) Overwrite this instance with an other or a C{scalar}.
        '''
        if other is self:
            pass  # from ._pow_scalar
        elif isinstance(other, Fsum):
            self._n     = other._n
            self._ps[:] = other._ps
            self._update(other=other)
        elif isscalar(other):
            s = other if asis else float(other)
            self._n     =  n
            self._ps[:] = [s]
            self._update(_fprs=s, _fprs2=Fsum2Tuple(s, INT0))
        else:  # PYCHOK no cover
            raise self._TypeError(_fset_, other)  # txt=_invalid_
        return self

    def fsub(self, xs=()):
        '''Subtract an iterable of C{scalar} or L{Fsum} instances
           from this instance.

           @arg xs: Iterable, list, tuple. etc. (C{scalar}
                    or L{Fsum} instances).

           @return: This instance, updated (L{Fsum}).

           @see: Method L{Fsum.fadd}.
        '''
        return self._facc(_2floats(xs, sub=True)) if xs else self

    def fsub_(self, *xs):
        '''Subtract all positional C{scalar} or L{Fsum} instances
           from this instance.

           @arg xs: Values to subtract (C{scalar} or
                    L{Fsum} instances), all positional.

           @return: This instance, updated (L{Fsum}).

           @see: Method L{Fsum.fadd}.
        '''
        return self._facc(_2floats(xs, origin=1, sub=True)) if xs else self

    def _fsub(self, other, op):
        '''(INTERNAL) Apply C{B{self} -= B{other}}.
        '''
        if isinstance(other, Fsum):
            if other is self:  # or other._fprs2 == self._fprs2:
                self._fset(_0_0, asis=True, n=len(self) * 2)  # self -= self
            elif other._ps:
                self._facc(map(neg, other._ps))
        elif not isscalar(other):
            raise self._TypeError(op, other, txt=_not_scalar_)
        elif other:
            self._facc_(-other)
        return self

    def _Fsum(self, x):
        '''(INTERNAL) Fast, single-scalar C{_Fsum}, named the same.
        '''
        f = _Fsum(x)  # C{int} or C{scalar}
        f._name = self.name  # .rename calls _update_all
        return f

    def fsum(self, xs=()):
        '''Add more C{scalar} or L{Fsum} instances and summate.

           @kwarg xs: Iterable, list, tuple, etc. (C{scalar} or
                      L{Fsum} instances).

           @return: Precision running sum (C{float} or C{int}).

           @see: Method L{Fsum.fadd}.

           @note: Accumulation can continue after summation.
        '''
        f = self._facc(_2floats(xs)) if xs else self
        return f._fprs

    def fsum_(self, *xs):
        '''Add all positional C{scalar} or L{Fsum} instances and summate.

           @arg xs: Values to add (C{scalar} or L{Fsum} instances),
                    all positional.

           @return: Precision running sum (C{float} or C{int}).

           @see: Method L{Fsum.fsum}.
        '''
        f = self._facc(_2floats(xs, origin=1)) if xs else self
        return f._fprs

    def fsum2(self, xs=()):
        '''Add more C{scalar} or L{Fsum} instances and return the
           current precision running sum and the C{residual}.

           @kwarg xs: Iterable, list, tuple, etc. (C{scalar} or
                      L{Fsum} instances).

           @return: L{Fsum2Tuple}C{(fsum, residual)} with C{fsum} the
                    current precision running sum and C{residual}, the
                    (precision) sum of the remaining C{partials}.  The
                    C{residual is INT0} if the C{fsum} is considered
                    to be I{exact}.

           @see: Methods L{Fsum.fint2}, L{Fsum.fsum} and L{Fsum.fsum2_}
        '''
        f = self._facc(_2floats(xs)) if xs else self
        return f._fprs2

    def fsum2_(self, *xs):
        '''Add any positional C{scalar} or L{Fsum} instances and return
           the precision running sum and the C{differential}.

           @arg xs: Values to add (C{scalar} or L{Fsum} instances),
                    all positional.

           @return: 2-Tuple C{(fsum, delta)} with the current precision
                    running C{fsum} and C{delta}, the difference with
                    the previous running C{fsum} (C{float}s).

           @see: Methods L{Fsum.fsum_} and L{Fsum.fsum}.
        '''
        p, r = self._fprs2
        s, t = self._facc(_2floats(xs, origin=1))._fprs2 if xs else (p, r)
        d = (s - p) + (t - r)  # == _fsum((s, -p, t, -r))
        return s, d

#   ftruediv = __itruediv__   # for naming consistency

    def _ftruediv(self, other, op):
        '''(INTERNAL) Apply C{B{self} /= B{other}}.
        '''
        if isinstance(other, Fsum):
            if other is self or other._fprs2 == self._fprs2:
                return self._fset(_1_0, asis=True)
            d, r = other._fprs2
            if r:
                if self._RESIDUAL:
                    raise self._ResidualError(op, other, r)
                if d:
                    # self / (d + r) == self * n / d
                    # n = d / (d + r) = 1 / (1 + r / d)
                    # d' = d / n = d * (1 + r / d), but
                    # may be moot if (1 + r / d) == 1
                    d *= r / d + _1_0
                else:  # PYCHOK no cover
                    d  = r
        elif isscalar(other):
            d = other
        else:  # PYCHOK no cover
            raise self._TypeError(op, other)  # txt=_invalid_
        try:
            s = 0 if isinf(d) else (
                d if isnan(d) else self._finite(_1_0 / d))
        except Exception as x:
            E, t = _xError2(x)
            raise self._Error(op, d, E, txt=t)
        f = self._mul_scalar(s, _mul_)  # handles 0, NAN, etc.
        return self._fset(f)

    @property_RO
    def imag(self):
        '''Get the C{imaginary} part of this instance (C{0.0}, always).

           @see: Properties L{Fsum.ceil}, L{Fsum.floor} and L{Fsum.real}.
        '''
        return _0_0

    def int_float(self, raiser=False):
        '''Return this instance' current running sum as C{int} or C{float}.

           @kwarg raiser: If C{True} throw a L{ResidualError} if the
                          residual is non-zero.

           @return: This C{integer} sum if this instance C{is_integer} or
                    this C{float} sum if the residual is zero or ignored.

           @raise ResidualError: Non-zero residual.

           @see: Methods L{Fsum.fint} and L{Fsum.fint2}.
        '''
        s, r = self._fint2
        if r:
            s, r = self._fprs2
            if r and raiser:  # PYCHOK no cover
                t = _2stresidual(_non_zero_, r)
                raise ResidualError(int_float=s, txt=t)
            s = float(s)  # redundant
        return s

    def is_exact(self):
        '''Is this instance' current running C{fsum} considered to be exact? (C{bool}).
        '''
        return self.residual is INT0

    def is_integer(self):
        '''Is this instance' current running sum C{integer}? (C{bool}).

           @see: Methods L{Fsum.fint} and L{Fsum.fint2}.
        '''
        _, r = self._fint2
        return not r

    def is_math_fsum(self):
        '''Return C{True} if functions L{fsum}, L{fsum_}, L{fsum1} and L{fsum1_}
           are all based on Python's C{math.fsum}, C{False} otherwise.
        '''
        return bool(Fsum._math_fsum)

    def _mul_Fsum(self, other, op):
        '''(INTERNAL) Return C{B{self} * Fsum B{other}} as L{Fsum}.
        '''
        # assert isinstance(other, Fsum)
        return self._Fsum(_0_0)._facc(self._ps_x(op, *other._ps))  # ._fpsqz()

    def _mul_scalar(self, factor, op):
        '''(INTERNAL) Return C{B{self} * scalar B{factor}} as L{Fsum} or C{0}.
        '''
        # assert isscalar(factor)
        if self._finite(factor, op) and self._ps:
            if factor != _1_0:
                f = _Fsum(0)
                f._n     = self._n  # also in ._fadd
                f._ps[:] = self._ps_x(op, factor)
                f._facc_up()
            else:
                f = self
        else:  # to allow setting f._n
            f = _Fsum(_0_0)
        return f

    @property_RO
    def partials(self):
        '''Get this instance' current partial sums (C{tuple} of C{float}s and/or C{int}s).
        '''
        return tuple(self._ps)

    def pow(self, x, *mod):
        '''Return C{B{self}**B{x}} as L{Fsum}.

           @arg x: The exponent (L{Fsum} or C{scalar}).
           @arg mod: Optional modulus (C{int} or C{None}) for the 3-argument
                     C{pow(B{self}, B{other}, B{mod})} version.

           @return: The C{pow(self, B{x})} or C{pow(self, B{x}, *B{mod})}
                    result (L{Fsum}).

           @note: If B{C{mod}} is given as C{None}, the result will be an
                  C{integer} L{Fsum} provided this instance C{is_integer}
                  or set C{integer} with L{Fsum.fint}.

           @see: Methods L{Fsum.__ipow__}, L{Fsum.fint} and L{Fsum.is_integer}.
        '''
        f = self.copy(name=self.pow.__name__)
        if f and isint(x) and x >= 0 and not mod:
            f._pow_int(x, x, _pow_)  # f **= x
        else:
            f._fpow(x, _pow_, *mod)  # f = pow(f, x, *mod)
        return f

    def _pow_0_1(self, x):
        '''(INTERNAL) Return B{C{self}**1} or C{B{self}**0 == 1.0}.
        '''
        return self if x else _1_0

    def _pow_2(self, x, other, op):
        '''(INTERNAL) 2-arg C{pow(B{self}, scalar B{x})} embellishing errors.
        '''
        # assert len(self._ps) == 1 and isscalar(x)
        b = self._ps[0]  # assert isscalar(b)
        try:  # type(s) == type(x) if x in (_1_0, 1)
            s = pow(b, x)  # -1**2.3 == -(1**2.2)
            if not iscomplex(s):
                return self._finite(s)  # 0**INF == 0.0, 1**INF==1.0
            # neg**frac == complex in Python 3+, but ValueError in 2-
            E, t = _ValueError, _2strcomplex(s, b, x)  # PYCHOK no cover
        except Exception as x:
            E, t = _xError2(x)
        raise self._Error(op, other, E, txt=t)

    def _pow_3(self, other, mod, op):
        '''(INTERNAL) 3-arg C{pow(B{self}, B{other}, int B{mod} or C{None})}.
        '''
        b, r = self._fprs2 if mod is None else self._fint2
        if r:  # and self._RESIDUAL:
            t = _non_zero_ if mod is None else _integer_
            E, t = ResidualError, _2stresidual(t, r, mod=mod)
        else:
            try:  # b, other, mod all C{int}, unless C{mod} is C{None}
                x = _2scalar(other, raiser=self._RESIDUAL)
                s =  pow(b, x, mod)
                if not iscomplex(s):
                    return self._finite(s)
                # neg**frac == complex in Python 3+, but ValueError in 2-
                E, t = _ValueError, _2strcomplex(s, b, x, mod)  # PYCHOK no cover
            except Exception as x:
                E, t = _xError2(x)
                t = _COMMASPACE_(Fmt.PARENSPACED(mod=mod), t)
        raise self._Error(op, other, E, txt=t)

    def _pow_int(self, x, other, op):
        '''(INTERNAL) Return C{B{self} **= B{x}} for C{int B{x} >= 0}.
        '''
        # assert isint(x) and x >= 0
        if len(self._ps) > 1:
            if x > 2:
                p = self._copy()
                m = 1  # single-bit mask
                if x & m:
                    x -= m  # x ^= m
                    f  = p._copy()
                else:
                    f = _Fsum(_1_0)
                while x:
                    p  = p._mul_Fsum(p, op)  # p **= 2
                    m += m  # m <<= 1
                    if x & m:
                        x -= m  # x ^= m
                        f  = f._mul_Fsum(p, op)  # f *= p
            elif x > 1:  # self**2
                f = self._mul_Fsum(self, op)
            else:  # self**1 or self**0
                f = self._pow_0_1(x)
        elif self._ps:  # self._ps[0]**x
            f = self._pow_2(x, other, op)
        else:  # PYCHOK no cover
            # 0**pos_int == 0, but 0**0 == 1
            f = 0 if x else 1  # like ._fprs
        return self._fset(f, asis=isint(f))

    def _pow_scalar(self, x, other, op):
        '''(INTERNAL) Return C{self**B{x}} for C{scalar B{x}}.
        '''
        s, r = self._fprs2
        if isint(x, both=True):
            x = int(x)  # Fsum**int
            y = abs(x)
            if y > 1:
                if r:
                    s = self._copy()._pow_int(y, other, op)
                    if x < 0:
                        _, r = s._fprs2
                        if r:
                            s, x = self._Fsum(_1_0)._ftruediv(s, op), None
                        else:
                            # use **= -1 for the CPython float_pow
                            # error if s is zero, and not s = 1 / s
                            s, x = s._fprs, -1
                    else:
                        x = None
            elif x < 0:  # self**-1 == 1 / self
                if r:
                    s, x = self._Fsum(_1_0)._ftruediv(self, op), None
            else:  # self**1 or self**0
                s, x = self._pow_0_1(x), None
        elif not isscalar(x):  # assert ...
            raise self._TypeError(op, other, txt=_not_scalar_)
        elif r:  # non-zero residual**fractional
            raise self._ResidualError(op, other, r, fractional=x)
        if x is not None:
            # assert isscalar(s) and isscalar(x)
            s = self._Fsum(s)._pow_2(x, other, op)
        return s  # C{int}, C{scalar}, C{self} or an L{Fsum}

    def _ps1(self, less):
        '''(INTERNAL) Yield partials, pseudo-sorted, 1-primed
            minus C{less} if non-zero.
        '''
        yield _1_0
        if less:
            yield -less
        for p in self._ps:
            if p > 0:
                yield p
        for p in self._ps:
            if p < 0:
                yield p
        yield _N_1_0

    def _ps_x(self, op, *factors):  # see .fmath.Fhorner
        '''(INTERNAL) Yield all C{partials} times each B{C{factor}},
           in total C{len(partials) * len(factors)} items.
        '''
        ps = self._ps
        if len(ps) < len(factors):
            ps, factors = factors, ps
        _f = _isfinite
        for f in factors:
            for p in ps:
                p *= f
                if _f(p):
                    yield p
                else:  # PYCHOK no cover
                    self._finite(p, op)  # throw ValueError

    @property_RO
    def real(self):
        '''Get the C{real} part of this instance (C{float}).

           @see: Methods L{Fsum.__float__} and L{Fsum.fsum}
                 and properties L{Fsum.ceil}, L{Fsum.floor},
                 L{Fsum.imag} and L{Fsum.residual}.
        '''
        return float(self._fprs)

    @property_RO
    def residual(self):
        '''Get this instance' residual (C{float} or C{int}), the sum of
           the C{partials} less the precision running sum C{fsum}.

           @note: If the C{residual is INT0}, the precision running
                  C{fsum} is considered to be I{exact}.

           @see: Methods L{Fsum.fsum}, L{Fsum.fsum2} and L{Fsum.is_exact}.
        '''
        return self._fprs2.residual

    def RESIDUAL(self, *raiser):
        '''Do or don't raise L{ResidualError} exceptions for this instance,
           overriding the default from env var C{PYGEODESY_FSUM_RESIDUAL}.

           @arg raiser: If C{True} throw L{ResidualError}s for division
                        and exponention, if C{False} don't, if C{None}
                        restore the default setting (C{bool}) and if
                        omitted, maintain the current setting.

           @return: The previous C{RESIDUAL} setting (C{bool}).
        '''
        r = self._RESIDUAL
        if raiser:
            s = raiser[0]
            self._RESIDUAL = Fsum._RESIDUAL if s is None else bool(s)
        return r

    def _ResidualError(self, op, other, residual, **name_values):
        '''(INTERNAL) Non-zero B{C{residual}} and C{name=value, ...} error.
        '''
        t = _2stresidual(_non_zero_, residual, **name_values)
        return self._Error(op, other, ResidualError, txt=t)

    def signOf(self, res=True):
        '''Determine the sign of this instance.

           @kwarg res: If C{True} consider, otherwise
                       ignore the residual (C{bool}).

           @return: The sign (C{int}, -1, 0 or +1).
        '''
        s, r = self._fprs2 if res else (self._fprs, 0)
        return _signOf(r, -s)

    def toRepr(self, prec=6, sep=_COMMASPACE_, fmt=Fmt.g, **unused):  # PYCHOK signature
        '''Return this C{Fsum} instance as representation.

           @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                        Trailing zero decimals are stripped for B{C{prec}} values
                        of 1 and above, but kept for negative B{C{prec}} values.
           @kwarg sep: Separator to join (C{str}).
           @kwarg fmt: Optional C{float} format (C{str}).

           @return: This instance (C{str}).
        '''
        t = sep.join(pairs(self._fprs2.items(), prec=prec, fmt=fmt))
        return _SPACE_(Fmt.SQUARE(self.named3, len(self)), Fmt.PAREN(t))

    def toStr(self, prec=6, sep=_COMMASPACE_, fmt=Fmt.g, **unused):  # PYCHOK signature
        '''Return this C{Fsum} instance as string.

           @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                        Trailing zero decimals are stripped for B{C{prec}} values
                        of 1 and above, but kept for negative B{C{prec}} values.
           @kwarg sep: Separator to join (C{str}).
           @kwarg fmt: Optional C{float} format (C{str}).

           @return: This instance (C{repr}).
        '''
        t = self._fprs2.toStr(prec=prec, sep=sep, fmt=fmt)
        return _SPACE_(Fmt.SQUARE(self.named3, len(self)), t)

    def _TypeError(self, op, other, **txt):  # PYCHOK no cover
        '''(INTERNAL) Return a C{TypeError}.
        '''
        return self._Error(op, other, _TypeError, **txt)

    def _update(self, other=None, **setters):
        '''(INTERNAL) Copy, set or zap all cached C{Property_RO} values.
        '''
        if other is None:  # zap if present
            Fsum._fint2._update(self)
            Fsum._fprs ._update(self)
            Fsum._fprs2._update(self)
        else:  # dup if present, otherwise zap
            Fsum._fint2._update_from(self, other)
            Fsum._fprs ._update_from(self, other)
            Fsum._fprs2._update_from(self, other)
        if setters:
            # Property_RO _fint2, _fprs and _fprs2 can't be a Property:
            # Property's _fset zaps the value just set by the @setter
            self.__dict__.update(setters)
        return self

    def _ValueError(self, op, other, **txt):  # PYCHOK no cover
        '''(INTERNAL) Return a C{ValueError}.
        '''
        return self._Error(op, other, _ValueError, **txt)

    def _ZeroDivisionError(self, op, other, **txt):  # PYCHOK no cover
        '''(INTERNAL) Return a C{ZeroDivisionError}.
        '''
        return self._Error(op, other, _ZeroDivisionError, **txt)


class _Fsum(Fsum):
    '''(INTERNAL) Fast, single-scalar L{Fsum}.
    '''
    _n = 1

    def __init__(self, x):
        self._ps = [x]  # assert isscalar(x) and isfinite(x)

_Fsum._name = _Fsum.__name__  # PYCHOK set once


def _Float_Int(arg, **name_Error):
    '''(INTERNAL) Unit of L{Fsum2Tuple} items.
    '''
    U = Int if isint(arg) else Float
    return U(arg, **name_Error)


class Fsum2Tuple(_NamedTuple):
    '''2-Tuple C{(fsum, residual)} with the precision running C{fsum}
       and the C{residual}, the sum of the remaining partials.  Each
       item is either C{float} or C{int}.

       @note: If the C{residual is INT0}, the C{fsum} is considered
              to be I{exact}.
    '''
    _Names_ = ( Fsum.fsum.__name__, _residual_)
    _Units_ = (_Float_Int,          _Float_Int)


class ResidualError(_ValueError):
    '''Error raised for an operation involving a L{pygeodesy.sums.Fsum}
       instance with a non-zero C{residual}, I{integer} or otherwise.

       @see: Module L{pygeodesy.fsums} and method L{Fsum.RESIDUAL}.
    '''
    pass


try:
    from math import fsum as _fsum  # precision IEEE-754 sum, Python 2.6+

    # make sure _fsum works as expected (XXX check
    # float.__getformat__('float')[:4] == 'IEEE'?)
    if _fsum((1, 1e101, 1, -1e101)) != 2:  # PYCHOK no cover
        del _fsum  # nope, remove _fsum ...
        raise ImportError  # ... use _fsum below

    Fsum._math_fsum = _fsum

except ImportError:

    def _fsum(xs):
        '''(INTERNAL) Precision summation, Python 2.5-.
        '''
        return Fsum(name=_fsum.__name__)._facc(xs)._fprs


def fsum(xs, floats=False):
    '''Precision floating point summation based on or like Python's C{math.fsum}.

       @arg xs: Iterable, list, tuple, etc. of values (C{scalar} or
                L{Fsum} instances).
       @kwarg floats: Optionally, set C{B{floats}=True} iff I{all}
                      B{C{xs}} are known to be C{float}.

       @return: Precision C{fsum} (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise TypeError: Non-scalar B{C{xs}} value.

       @raise ValueError: Invalid or non-finite B{C{xs}} value.

       @note: Exceptions and I{non-finite} handling may differ if not
              based on Python's C{math.fsum}.

       @see: Class L{Fsum} and methods L{Fsum.fsum} and L{Fsum.fadd}.
    '''
    return _fsum(_2floats(xs, floats=floats)) if xs else _0_0


def fsum_(*xs, **floats):
    '''Precision floating point summation of all positional arguments.

       @arg xs: Values to be added (C{scalar} or L{Fsum} instances),
                all positional.
       @kwarg floats: Optionally, set C{B{floats}=True} iff I{all}
                      B{C{xs}} are known to be C{float}.

       @return: Precision C{fsum} (C{float}).

       @see: Function C{fsum}.
    '''
    return _fsum(_2floats(xs, origin=1, **_floats(floats))) if xs else _0_0


def fsum1(xs, floats=False):
    '''Precision floating point summation of a few values, 1-primed.

       @arg xs: Iterable, list, tuple, etc. of values (C{scalar} or
                L{Fsum} instances).
       @kwarg floats: Optionally, set C{B{floats}=True} iff I{all}
                      B{C{xs}} are known to be C{float}.

       @return: Precision C{fsum} (C{float}).

       @see: Function C{fsum}.
    '''
    return _fsum(_2floats(xs, primed=True, floats=floats)) if xs else _0_0


def fsum1_(*xs, **floats):
    '''Precision floating point summation of a few arguments, 1-primed.

       @arg xs: Values to be added (C{scalar} or L{Fsum} instances),
                all positional.
       @kwarg floats: Optionally, set C{B{floats}=True} iff I{all}
                      B{C{xs}} are known to be C{float}.

       @return: Precision C{fsum} (C{float}).

       @see: Function C{fsum}
    '''
    return _fsum(_2floats(xs, origin=1, primed=True, **_floats(floats))) if xs else _0_0


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
