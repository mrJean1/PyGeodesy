
# -*- coding: utf-8 -*-

u'''Class L{Fsum} for precision I{running} floating point summation.

Set env variable C{PYGEODESY_FSUM_RESIDUAL} to "std" to throw an
exception for division or exponention by an L{Fsum} instance with
a non-zero C{residual}, see method L{Fsum.fsum2}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import _isfinite, isint, isscalar, map1, neg
from pygeodesy.errors import _NotImplementedError, _OverflowError, \
                             _TypeError, _ValueError, _xkwds_get
from pygeodesy.interns import NN, _COMMASPACE_, _EQUAL_, _finite_, \
                             _iadd_, _not_, _SPACE_, _supported_, \
                             _0_0, _1_0, _N_1_0
from pygeodesy.lazily import _ALL_LAZY, _getenv, _sys_version_info2
from pygeodesy.named import _Named, _NamedTuple, _NotImplemented
# from pygeodesy.props import Property_RO, property_RO  # from .units
from pygeodesy.streprs import Fmt, fstr, pairs, unstr
from pygeodesy.units import Float, Property_RO, property_RO

from math import ceil as _ceil, floor as _floor

__all__ = _ALL_LAZY.fsums
__version__ = '22.02.01'

_eq_   = _EQUAL_ * 2
_ge_   = '>='
_gt_   = '>'
_le_   = '<='
_lt_   = '<'
_ne_   = '!='
_idiv_ = '/='
_iflr_ = '//='
_imul_ = '*='
_ipow_ = '**='
_iset_ = _EQUAL_
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


def _2Fsum(other, name=NN):
    '''(INTERNAL) Return B{C{other}} as an L{Fsum}.
    '''
    if not isinstance(other, Fsum):
        other = _2float(other=other)
    return Fsum(name=name)._iset(other)


def _2sum(a, b):  # by .testFmath
    '''(INTERNAL) Precision C{2sum} of M{a + b} as 2-tuple (sum, residual).
    '''
    s = a + b
    if not _isfinite(s):
        raise _OverflowError(unstr(_2sum.__name__, a, b), txt=str(s))
    if abs(a) < abs(b):
        a, b = b, a
    return s, (b - (s - a))  # abs(b) <= abs(a)


class ResidualError(_ValueError):
    '''Error raised for an operation involving an L{Fsum} with a non-zero residual.
    '''
    pass


class Fsum(_Named):
    '''Precision I{running} floating point summation similar to standard Python's C{math.fsum}.

       Unlike C{math.fsum}, this class accumulates values and provides intermediate,
       I{running} precision floating point summation.  Accumulation may continue
       after intermediate, I{running} summations.

       @note: Handling of exceptions, C{inf}, C{INF}, C{nan} and C{NAN} values is
              different from Python's C{math.fsum}.

       @see: U{Hettinger<https://GitHub.com/ActiveState/code/blob/master/recipes/Python/
             393090_Binary_floating_point_summatiaccurate_full/recipe-393090.py>},
             U{Kahan<https://WikiPedia.org/wiki/Kahan_summation_algorithm>},
             U{Klein<https://Link.Springer.com/article/10.1007/s00607-005-0139-x>},
             Python 2.6+ file I{Modules/mathmodule.c} and the issue log
             U{Full precision summation<https://Bugs.Python.org/issue2819>}.
    '''
    _n  = 0
    _ps = []  # partials
    _Rx = bool(_getenv('PYGEODESY_FSUM_RESIDUAL', NN))

    def __init__(self, *starts, **name_NN):
        '''New L{Fsum} for I{running} precision floating point summation.

           @arg starts: No, one or more initial values (C{scalar}s).
           @kwarg name_NN: Optional name (C{str}).

           @see: Method L{Fsum.fadd}.
        '''
#       self._n  = 0
        self._ps = []
        if name_NN:
            self.name = _xkwds_get(name_NN, name=NN)
        if starts:
            self.fadd(starts)

    def __abs__(self):
        '''Return this instance' absolute value as an L{Fsum}.
        '''
        f = self.fcopy(name=self.__abs__.__name__)
        return f._ineg() if f < 0 else f

    def __add__(self, other):
        '''Return the sum C{B{self} + B{other}} as an L{Fsum}.

           @arg other: An L{Fsum} or C{scalar}.

           @return: The sum (L{Fsum}).

           @see: Method L{Fsum.__iadd__}.
        '''
        f  = self.fcopy(name=self.__add__.__name__)
        f += other
        return f

    def __bool__(self):  # PYCHOK not special in Python 2-
        '''Return C{True} if this instance is non-zero.
        '''
        s, r = self._fsum2
        return bool(s or r)

    def __ceil__(self):  # PYCHOK not special in Python 2-
        '''Return this instance' C{math.ceil} as C{float}.

           @see: Methods L{Fsum.__floor__} and L{Fsum.__int__}.
        '''
        return _ceil(max(self._fsums2))

    def __divmod__(self, other):
        '''Return C{divmod(self, B{other})} as 2-tuple C{(quotient,
           remainder)} as an integer C{float} and an L{Fsum}.

           @arg other: An L{Fsum} or C{scalar} divisor.

           @see: Method L{Fsum.__itruediv__}.
        '''
        return self._divmod(other, Fsum.__divmod__.__name__)

    def __eq__(self, other):
        '''Compare this with an other instance or scalar.
        '''
        s, r = self._cmp2(other, _eq_)
        return not bool(s or r)

    def __float__(self):
        '''Return this instance as C{float}, without residual.

           @see: Method L{Fsum.fsum2}.
        '''
        return self._fsum

    def __floor__(self):  # PYCHOK not special in Python 2-
        '''Return this instance' C{math.floor} as C{float}.

           @see: Methods L{Fsum.__ceil__} and L{Fsum.__int__}.
        '''
        return _floor(min(self._fsums2))

    def __floordiv__(self, other):
        '''Return C{B{self} // B{other}} as an L{Fsum}.

           @arg other: An L{Fsum} or C{scalar} divisor.

           @return: The C{floor} quotient (L{Fsum}).

           @see: Methods L{Fsum.__ifloordiv__}.
        '''
        f = self.fcopy(name=self.__floordiv__.__name__)
        return f.__ifloordiv__(other)  # //= chokes PyChecker

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

           @arg other: An L{Fsum} or C{scalar}.

           @return: This instance, updated (L{Fsum}).

           @raise TypeError: Invalid B{C{other}} type.

           @see: Method L{Fsum.fadd}.
        '''
        if isinstance(other, Fsum):
            if other is self:  # or self.__eq__(other):
                self._fmul(2)
            elif other._ps:
                self._fadd(other._ps)
        elif not isscalar(other):
            raise self._TypeError(_iadd_, other)
        elif other:
            self._fadd_(other)
        return self

    def __ifloordiv__(self, other):
        '''Apply C{B{self} //= B{other}} to this instance.

           @arg other: An L{Fsum} or C{scalar} divisor.

           @return: This instance, updated (L{Fsum}).

           @raise ResidualError: Non-zero residual in B{C{other}}.

           @raise TypeError: Invalid B{C{other}} type.

           @raise ValueError: Zero, invalid or non-finite B{C{other}}.

           @see: Methods L{Fsum.__itruediv__}.
        '''
        q = float(self._idiv(other, _iflr_))
        return self._iset(_floor(q))

    def __imatmul__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __imod__(self, other):
        '''Apply C{B{self} %= B{other}} to this instance.

           @arg other: An L{Fsum} or C{scalar} divisor.

           @return: This instance, updated (L{Fsum}).

           @raise TypeError: Invalid B{C{other}} type.

           @see: Method L{Fsum.__divmod__}.
        '''
        _, m = self.__divmod__(other)
        return self._iset(m)

    def __imul__(self, other):
        '''Apply C{B{self} *= B{other}} to this instance.

           @arg other: An L{Fsum} or C{scalar} factor.

           @return: This instance, updated (L{Fsum}).

           @raise TypeError: Invalid B{C{other}} type.

           @see: Method L{Fsum.fmul}.
        '''
        if isinstance(other, Fsum):
            ps = list(other._ps)
            if ps and self._ps:
                s = self.fcopy()._fmul(ps.pop())
                while ps:  # s += self * ps.pop()
                    p = ps.pop()
                    if p:
                        s._fadd(self._ps_x(p), update=False)
                s._update()  # XXX s._fadd_1()?
            else:
                s = _0_0
            self._iset(s)
        elif isscalar(other):
            self._fmul(other)
        else:
            raise self._TypeError(_imul_, other)
        return self

    def __int__(self):
        '''Return this instance as an C{int}.

           @see: Methods L{Fsum.__ceil__} and L{Fsum.__floor__}.
        '''
        return int(self._fsum)

    def __ipow__(self, other, *mod):  # PYCHOK 2 vs 3 args
        '''Apply C{B{self} **= B{other}} to this instance.

           @arg other: An L{Fsum} or C{scalar} exponent.
           @arg mod: Not implemented (C{scalar}).

           @return: This instance, updated (L{Fsum}).

           @raise NotImplementedError: Argument B{C{mod}} used.

           @raise ResidualError: This residual non-zero and negative
                                 or fractional B{C{other}} or non-zero
                                 residual in B{C{other}}.

           @raise TypeError: Invalid B{C{other}} type.

           @raise ValueError: If B{C{other}} is a negative C{scalar}
                              and this instance is C{0} or B{C{other}}
                              is a fractional C{scalar} and this
                              instance is negative or has a non-zero
                              residual or B{C{other}} is an L{Fsum}
                              with a non-zero residual.

           @see: CPython function U{float_pow<https://GitHub.com/
                 python/cpython/blob/main/Objects/floatobject.c>}.
        '''
        if mod:
            return _NotImplemented(self, other, *mod)

        x, p = other, _1_0
        if isinstance(x, Fsum):
            x, r = x._fsum2
            if r:
                if Fsum._Rx:
                    raise self._ResidualError(_ipow_, other, r)
                p = self._pow(r, other)  # scalar or Fsum
        s = self._pow(x, other)  # scalar or Fsum
        if p != _1_0:
            s *= p
        return self._iset(s)

    def __isub__(self, other):
        '''Apply C{B{self} -= B{other}} to this instance.

           @arg other: An L{Fsum} or C{scalar}.

           @return: This instance, updated (L{Fsum}).

           @raise TypeError: Invalid B{C{other}} type.

           @see: Method L{Fsum.fadd}.
        '''
        if isinstance(other, Fsum):
            if other is self:  # or self.__eq__(other):
                self._iset(_0_0)
            elif other._ps:
                self._fadd(map(neg, other._ps))
        elif not isscalar(other):
            raise self._TypeError(_isub_, other)
        elif other:
            self._fadd_(neg(other))
        return self

    def __itruediv__(self, other):
        '''Apply C{B{self} /= B{other}} to this instance.

           @arg other: An L{Fsum} or C{scalar} divisor.

           @return: This instance, updated (L{Fsum}).

           @raise TypeError: Invalid B{C{other}} type.

           @raise ValueError: Zero, invalid or non-finite B{C{other}}
                              or an B{C{other}} is an L{Fsum} with a
                              non-zero residual.

           @see: Method L{Fsum.__ifloordiv__}.
        '''
        return self._idiv(other, _idiv_)

    def __le__(self, other):
        '''Compare this with an other instance or scalar.
        '''
        s, r = self._cmp2(other, _le_)
        return r <= -s

    def __len__(self):
        '''Return the I{total} number of value accumulated (C{int}).
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
        f = self.fcopy(name=self.__mod__.__name__)
        return f.__imod__(other)

    def __mul__(self, other):
        '''Return C{B{self} * B{other}} as an L{Fsum}.

           @see: Method L{Fsum.__imul__}.
        '''
        f  = self.fcopy(name=self.__mul__.__name__)
        f *= other
        return f

    def __ne__(self, other):
        '''Compare this with an other instance or scalar.
        '''
        s, r = self._cmp2(other, _ne_)
        return bool(s or r)

    def __neg__(self):
        '''Return I{a copy of} this instance, negated.
        '''
        f = self.fcopy(name=self.__neg__.__name__)
        return f._ineg()

    def __pos__(self):  # PYCHOK no cover
        '''Return this instance, I{as-is}.
        '''
        return self.fcopy(name=self.__pos__.__name__)

    def __pow__(self, other, *mod):  # PYCHOK 2 vs 3 args
        '''Return C{B{self} ** B{other}} as L{Fsum}.

           @see: Method L{Fsum.__ipow__}.
        '''
        f = self.fcopy(name=self.__pow__.__name__)
        return f.__ipow__(other, *mod)  # to pass C{*mod}

    def __radd__(self, other):
        '''Return C{B{other} + B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__iadd__}.
        '''
        f = _2Fsum(other, name=self.__radd__.__name__)
        return f._fadd(self._ps)

    def __rdivmod__(self, other):
        '''Return C{divmod(B{other}, B{self})} as 2-tuple C{(quotient,
           remainder)} of an integer C{float} and an L{Fsum}.

           @see: Method L{Fsum.__divmod__}.
        '''
        f = _2Fsum(other)
        return f._divmod(self, Fsum.__rdivmod__.__name__)

#   def __repr__(self):
#       '''Return the default C{repr(this)}.
#       '''
#       return self.toRepr()

    def __rfloordiv__(self, other):
        '''Return C{B{other} // B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__ifloordiv__}.
        '''
        f = _2Fsum(other, name=self.__rfloordiv__.__name__)
        return f.__ifloordiv__(self)  # //= chokes PyChecker

    def __rmatmul__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __rmod__(self, other):
        '''Return C{B{other} % B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__imod__}.
        '''
        f  = _2Fsum(other, name=self.__rmod__.__name__)
        f %=  self
        return f

    def __rmul__(self, other):
        '''Return C{B{other} * B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__imul__}.
        '''
        f  = _2Fsum(other, name=self.__rmul__.__name__)
        f *=  self
        return f

    def __round__(self, ndigits=None):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, ndigits=ndigits)

    def __rpow__(self, other, *mod):
        '''Return C{B{other} ** B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__ipow__}.
        '''
        f = _2Fsum(other, name=self.__rpow__.__name__)
        return f.__ipow__(self, *mod)  # to pass *mod

    def __rsub__(self, other):
        '''Return C{B{other} - B{self}} as L{Fsum}.

           @see: Method L{Fsum.__isub__}.
        '''
        f  = _2Fsum(other, name=self.__rsub__.__name__)
        f -=  self
        return f

    def __rtruediv__(self, other):
        '''Return C{B{other} / B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__itruediv__}.
        '''
        f = _2Fsum(other, name=self.__rtruediv__.__name__)
        return f.__itruediv__(self)  # /= chokes PyChecker

    def __sizeof__(self):  # PYCHOK not special in Python 2-
        '''Return the size of this instance in C{bytes}.
        '''
        from sys import getsizeof
        return sum(map1(getsizeof, self._fsum, self._fsum2, self._n, self._ps, *self._ps))

    def __str__(self):
        '''Return the default C{str(this)}.
        '''
        return self.toStr()

    def __sub__(self, other):
        '''Return C{B{self} - B{other}} as an L{Fsum}.

           @arg other: An L{Fsum} or C{scalar}.

           @return: The difference (L{Fsum}).

           @see: Method L{Fsum.__isub__}.
        '''
        f  = self.fcopy(name=self.__sub__.__name__)
        f -= other
        return f

    def __truediv__(self, other):
        '''Return C{B{self} / B{other}} as an L{Fsum}.

           @arg other: An L{Fsum} or C{scalar} divisor.

           @return: The quotient (L{Fsum}).

           @see: Method L{Fsum.__itruediv__}.
        '''
        f = self.fcopy(name=self.__truediv__.__name__)
        return f.__itruediv__(other)  # /= chockes PyChecker

    __trunc__ = __int__

    if _sys_version_info2 < (3, 0):  # PYCHOK no cover
        # <https://docs.Python.org/2/library/operator.html#mapping-operators-to-functions>
        __div__     = __truediv__
        __idiv__    = __itruediv__
        __long__    = __int__
        __nonzero__ = __bool__
        __rdiv__    = __rtruediv__

    def _cmp2(self, other, cop):
        '''(INTERNAL) Subtract an B{C{other}} instance or scalar
           and return 2-tuple C{(fsum, residual)}, both C{float}
           for comparison operator B{C{cop}}.
        '''
        if isinstance(other, Fsum):
            f = Fsum()  # name=self._cmp2.__name__
            f._ps[:] = self._ps  # fast fcopy
            f._fadd(map(neg, other._ps))  # fast f -= other
            s = f.fsum()  # fast _fsums2  XXX = fsum(f._ps)
            r = fsum(f._ps1(s)) if len(f._ps) > 1 else 0
        elif isscalar(other):
            s, r = self._fsum2
            s -= other
        else:
            raise self._TypeError(cop, other)
        return s, r

    def _divmod(self, other, name):
        '''(INTERNAL) C{divmod(B{self}, B{other})} as 2-tuple (C{float}, L{Fsum}).
        '''
        # mostly like like CPython function U{float_divmod
        # <https://GitHub.com/python/cpython/blob/main/Objects/floatobject.c>},
        # but at least _divmod(-3, 2) equals Cpython's result (-2, 1).
        s, r = self._fsum2
        if r or isinstance(other, Fsum):
            f = self.fcopy(name=name)
            q = float(f // other)
            if q:
                f -= other * q
            if (f < 0 and other > 0) or \
               (f > 0 and other < 0):  # make sign(f) == sign(other)
                f += other
                q -= _1_0
        else:
            q, f = divmod(s, _2float(other=other))
            q, f = float(q), Fsum(name=name)._iset(f)
        return q, f

    def _Error(self, op, other, Error, **txt):
        '''(INTERNAL) Format an B{C{Error}} for C{{self} B{op} B{other}}.
        '''
        return Error(_SPACE_(self.toRepr(), op, repr(other)), **txt)

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

        def _xs(xs):
            _2f = _2float
            for i, x in enumerate(xs):
                x = _2f(index=i, xs=x)
                if x:
                    yield x

        return self._fadd(_xs(xs))

    def fadd_(self, *xs):
        '''Accumulate more I{scalar} values from positional arguments.

           @arg xs: Values to add (C{scalar}s), all positional.

           @return: This instance (L{Fsum}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{xs}} value.

           @raise ValueError: Invalid or non-finite B{C{xs}} value.
        '''
        return self.fadd(xs)

    def _fadd(self, xs, update=True):
        '''(INTERNAL) Accumulate more C{scalar}s.
        '''
        ps, n, _2s = self._ps, 0, _2sum
        for x in xs:  # _iter()
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
            if update:
                self._update()
        return self

    def _fadd_(self, *xs):
        '''(INTERNAL) Add all positional C{scalar}s.
        '''
        return self._fadd(xs)

    def _fadd_1(self):
        '''(INTERNAL) Adjust the C{partials}, by re-adding
           the final C{partial}.
        '''
        while self._ps:
            p = self._ps.pop()
            if p:
                self._fadd_(p)
                self._n -= 1
                break
        else:
            self._iset(_0_0)
        return self

    def fcopy(self, deep=False, name=NN):
        '''Copy this instance, C{shallow} or B{C{deep}}.

           @return: The copy (L{Fsum}).
         '''
        f = _Named.copy(self, deep=deep, name=name)  # see .__neg__
#       f._update(other=self)
        f._n  = self._n if deep else 1
        f._ps = list(self._ps)  # separate list
        return f

    copy = fcopy

    def fdiv(self, divisor):
        '''Devide this instance by a I{scalar}.

           @arg divisor: The denominator (C{scalar}).

           @raise TypeError: Non-scalar B{C{divisor}}.

           @raise ValueError: Zero, invalid or non-finite B{C{divisor}}.

           @return: This instance (L{Fsum}).
        '''
        return self._fdiv(_2float(divisor=divisor))

    def _fdiv(self, divisor, op=_idiv_):
        '''(INTERNAL) Divide this instance by a C{scalar} only.
        '''
        if divisor != _1_0:
            try:
                self._fmul(_1_0 / divisor)
            except (TypeError, ValueError, ZeroDivisionError) as x:
                raise self._ValueError(op, divisor, txt=str(x))
        return self

    def fmul(self, factor):
        '''Multiple this instance by a I{scalar}.

           @arg factor: The multiplier (C{scalar}).

           @return: This instance, updated (L{Fsum}).

           @raise TypeError: Non-scalar B{C{factor}}.

           @raise ValueError: Invalid or non-finite B{C{factor}}.

           @see: Method L{Fsum.fadd}.
        '''
        return self._fmul(_2float(factor=factor))

    def _fmul(self, factor):
        '''(INTERNAL) Multiply this instance by a C{scalar} only.
        '''
        if factor and self._ps:
            if factor != _1_0:
                # multiply and adjust partials
                self._ps[:] = self._ps_x(factor)
                self._fadd_1()
        else:
            self._iset(_0_0)
        # assert self._ps is ps
        return self

    def fpow(self, x):
        '''Return I{a copy of} this instance raised to power B{C{x}}.

           @arg x: The exponent (C{scalar}).

           @return: The C{pow(self, B{x})} (L{Fsum}).

           @raise ResidualError: This residual non-zero and negative
                                 or fractional B{C{x}}.

           @raise TypeError: Non-scalar B{C{x}}.

           @raise ValueError: Invalid or non-finite B{C{factor}}.

           @see: Method L{Fsum.__ipow__}.
        '''
        s = Fsum(name=self.fpow.__name__)
        if isint(x) and x >= 0:
            m = 1
            p = self.fcopy()
            s._iset(_1_0)
            while True:
                if x & m:
                    s *= p
                    x -= m
                if x:
                    m += m
                    p *= p
                else:
                    break
        else:
            s._iset(self._pow(_2float(x=x), x))
        return s

    def fsub(self, xs):
        '''Subtract several I{scalar} values.

           @arg xs: Iterable, list, tuple. etc. (C{scalar}s).

           @return: This instance, updated (L{Fsum}).

           @see: Method L{Fsum.fadd}.
        '''
        if xs:
            self.fadd(map(neg, xs))
        return self

    def fsub_(self, *xs):
        '''Subtract any I{scalar} positional value.

           @arg xs: Values to subtract (C{scalar}s), all positional.

           @return: This instance, updated (L{Fsum}).

           @see: Method L{Fsum.fadd}.
        '''
        return self.fsub(xs)

    def fsum(self, xs=None):
        '''Add more I{scalar} values and sum all.

           @kwarg xs: Iterable, list, tuple, etc. (C{scalar}s).

           @return: Precision running sum (C{float}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{xs}} value.

           @raise ValueError: Invalid or non-finite B{C{xs}} value.

           @note: Accumulation can continue after summation.
        '''
        if xs:
            self.fadd(xs)

        ps = self._ps
        i = len(ps) - 1
        if i < 0:
            s = _0_0
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
                        ps[i] = s = p  # swap s and p, continue
            # assert self._ps is ps
        self._update(_fsum=s)
        return s

    @Property_RO
    def _fsum(self):
        '''(INTERNAL) Cache the L{Fsum.fsum} result.
        '''
        return self.fsum()

    def fsum_(self, *xs):
        '''Add any I{scalar} positional value and sum all.

           @arg xs: Values to add (C{scalar}s), all positional.

           @return: Precision running sum (C{float}).

           @see: Method L{Fsum.fsum}.

           @note: Accumulation can continue after summation.
        '''
        return self.fsum(xs)

    def fsum2(self, xs=None):
        '''Add more I{scalar} values and return the precision
           running sum and the residual.

           @kwarg xs: Iterable, list, tuple, etc. (C{scalar}s).

           @return: L{Fsum2Tuple}C{(fsum, residual)} with the
                    running C{fsum} and C{residual}, the sum
                    of the remaining C{partials}.

           @see: Methods L{Fsum.fsum} and L{Fsum.fsum2_}
        '''
        if xs:
            self.fsum(xs)
        return self._fsum2

    @Property_RO
    def _fsum2(self):
        '''(INTERNAL) Cache the L{Fsum.fsum2} result.
        '''
        s = self._fsum
        r = fsum(self._ps1(s)) if len(self._ps) > 1 else 0
        t = Fsum2Tuple(s, r) if s else Fsum2Tuple((r if r else _0_0), 0)
        return t

    def fsum2_(self, *xs):
        '''Add any I{scalar} positional value and return the
           precision running sum and the difference.

           @arg xs: Values to add (C{scalar}s), all positional.

           @return: 2-Tuple C{(sum, delta)} with the precision
                    running C{sum} and C{delta} with the prior
                    running C{sum} (C{float}s).

           @see: Method L{Fsum.fsum_}.

           @note: Accumulation can continue after summation.
        '''
        p = self._fsum
        s = self.fsum(xs)  # if xs else p
        return s, (s - p)

    @Property_RO
    def _fsums2(self):
        '''(INTERNAL) This instance as 2-tuple C{(L{Fsum.fsum},
           L{fsum1}(partials))}, both C{float}s.
        '''
        s = self._fsum
        p = fsum(self._ps1(0)) if len(self._ps) > 1 else s
        return s, p

    def _idiv(self, other, op):
        '''(INTERNAL) Apply C{B{self} /= B{other}}.
        '''
        if isinstance(other, Fsum):
            if other is self or self.__eq__(other):
                self._iset(_1_0)
                return self
            d, r = other._fsum2
            if r:
                if Fsum._Rx:
                    raise self._ResidualError(op, other, r)
                if d:
                    # self / (d + r) == self * n / d
                    # n = d / (d + r) = 1 / (1 + r / d)
                    # d' = d / n = d * (1 + r / d), may
                    # be pointless if (r / d + 1) == 1
                    r = d * (r / d + _1_0)
                d = r
        elif isscalar(other):
            d = other
        else:
            raise self._TypeError(op, other)
        return self._fdiv(d)

    @property_RO
    def imag(self):
        '''Return the C{imaginary} part of this instance (C{0.0}, always).
        '''
        return _0_0

    def _ineg(self):
        '''(INTERNAL) Negate this instance.
        '''
        if self._ps:
            self._ps[:] = map(neg, self._ps)
            self._fadd_1()
        return self

    def _iset(self, other):
        '''(INTERNAL) Overwrite this instance with an other or a C{scalar}.
        '''
        if other is self:
            pass  # from ._pow
        elif isinstance(other, Fsum):
            self._n     = other._n
            self._ps[:] = other._ps
            self._update(other=other)
        elif isscalar(other):
            f = float(other)
            t = Fsum2Tuple(f, 0)
            self._n     =  1
            self._ps[:] = [f] if f else []
            self._update(_fsum=f, _fsum2=t, _fsums2=(f, f))
        else:  # PYCHOK no cover
            raise self._TypeError(_iset_, other)
        return self

    def is_integer(self):
        '''Return C{True} if this instance is an integer, C{False} otherwise.
        '''
        s, r = self._fsum2
        return (not r) and s.is_integer()

    def _pow(self, x, other):  # MCCABE 14
        '''(INTERNAL) Return C{self ** B{x}} for C{scalar B{x}} only.
        '''
        s, r = self._fsum2
        if isint(x, both=True):
            x = int(x)  # Fsum ** integer
            y = abs(x)
            if y > 1:
                if r:
                    s = self.fpow(y)
                    if x < 0:
                        s, r = s._fsum2
                        if r:
                            raise self._ResidualError(_ipow_, other, r)
                        # use **= -1 for the CPython float_pow
                        # error if s is zero, and not s = 1 / s
                        x = -1
                    else:
                        x = None
            elif x < 0:  # x == -1
                if r:
                    raise self._ResidualError(_ipow_, other, r)
            else:  # x == 1 or x == 0
                s = self if x else _1_0
                x = None
        elif not isscalar(x):
            raise self._TypeError(_ipow_, other)
        elif r:  # residual Fsum ** fractional
            raise self._ResidualError(_ipow_, other, r)
        elif s < 0:  # negative**fractional yields complex
            raise self._ValueError(_ipow_, other, txt=_not_(_supported_))
        if x not in (None, _1_0):
            try:
                s **= x
            except Exception as e:
                raise self._ValueError(_ipow_, other, txt=str(e))
        return s  # C{scalar} or an L{Fsum}

    def _ps1(self, less):
        '''(INTERNAL) Yield partials, pseudo-sorted, 1-primed.
        '''
        yield _1_0
        for p in self._ps:
            if p > 0:
                yield p
        if less:
            yield -less
        for p in self._ps:
            if p < 0:
                yield p
        yield _N_1_0

    def _ps_x(self, factor):
        '''(INTERNAL) Yield partials times B{C{factor}}.
        '''
        for p in self._ps:
            yield p * factor

    @property_RO
    def real(self):
        '''Return the C{real} part of this instance (C{float}).

           @see: Method L{Fsum.__float__}.
        '''
        return self._fsum

    def _ResidualError(self, op, other, residual):
        '''(INTERNAL) Non-zero residual C{ValueError}.
        '''
        t = _SPACE_(_non_zero_, _residual_)
        t =  Fmt.PARENSPACED(t, fstr(residual, fmt=Fmt.e, prec=8))
        return self._Error(op, other, ResidualError, txt=t)

    def signOf(self, res=True):
        '''Determine the sign of this instance.

           @kwarg res: If C{True} include the residual,
                       otherwise ignore (C{bool}).

           @return: The sign (C{int}, -1, 0 or +1).
        '''
        s, r = self._fsum2 if res else (self._fsum, 0)
        return 1 if r > -s else (-1 if r < -s else 0)

    def toRepr(self, prec=6, sep=_COMMASPACE_, fmt=Fmt.g, **unused):  # PYCHOK signature
        '''Return this C{Fsum} instance as representation.

           @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                        Trailing zero decimals are stripped for B{C{prec}} values
                        of 1 and above, but kept for negative B{C{prec}} values.
           @kwarg sep: Optional separator to join (C{str}).
           @kwarg fmt: Optional, C{float} format (C{str}).

           @return: This instance (C{str}).
        '''
        t = sep.join(pairs(self._fsum2.items(), prec=prec, fmt=fmt))
        return _SPACE_(Fmt.SQUARE(self.named3, len(self)), Fmt.PAREN(t))

    def toStr(self, prec=6, sep=_COMMASPACE_, fmt=Fmt.g, **unused):  # PYCHOK signature
        '''Return this C{Fsum} instance as string.

           @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                        Trailing zero decimals are stripped for B{C{prec}} values
                        of 1 and above, but kept for negative B{C{prec}} values.
           @kwarg sep: Optional separator to join (C{str}).
           @kwarg fmt: Optional, C{float} format (C{str}).

           @return: This instance (C{repr}).
        '''
        t = self._fsum2.toStr(prec=prec, sep=sep, fmt=fmt)
        return _SPACE_(Fmt.SQUARE(self.named3, len(self)), t)

    def _TypeError(self, op, other, **txt):
        '''(INTERNAL) Operand C{TypeError}.
        '''
        return self._Error(op, other, _TypeError, **txt)

    def _update(self, other=None, **setters):
        '''(INTERNAL) Copy, set or zap all cached C{Property_RO} values.
        '''
        if other is None:  # zap all
            Fsum._fsum._update(self)
            Fsum._fsum2._update(self)
            Fsum._fsums2._update(self)
        else:  # copy if present, otherwise zap
            Fsum._fsum._update_from(self, other)
            Fsum._fsum2._update_from(self, other)
            Fsum._fsums2._update_from(self, other)
        if setters:
            # Property_RO ._fsum, _fsum2 and ._fsums2 can not be a Property
            # since Property's _fset zaps the value just set by the @setter
            # plus the Property and the value attribute are named the same
            self.__dict__.update(setters)
        return self

    def _ValueError(self, op, other, **txt):  # PYCHOK no cover
        '''(INTERNAL) Return a C{ValueError}.
        '''
        return self._Error(op, other, _ValueError, **txt)


class Fsum2Tuple(_NamedTuple):
    '''2-Tuple C{(fsum, residual)} with the precision running C{fsum} and
       C{residual} the sum of the remaining partials, both C{float}.
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
        '''Precision floating point summation similar to Python's C{math.fsum}.

           Exception and I{non-finite} handling differs from C{math.fsum}.

           @arg xs: Iterable, list, tuple, etc. of values (C{scalar}s).

           @return: Precision C{sum} (C{float}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{xs}} value.

           @raise ValueError: Invalid or non-finite B{C{xs}} value.

           @see: Class L{Fsum} and methods L{Fsum.fsum} and L{Fsum.fadd}.
        '''
        return Fsum(name=fsum.__name__).fsum(xs) if xs else _0_0


def fsum_(*xs):
    '''Precision floating point summation of all positional arguments.

       @arg xs: Values to be added (C{scalar}s).

       @return: Precision L{fsum} (C{float}).

       @see: Function C{fsum}.
    '''
    return fsum(xs)


def fsum1(xs):
    '''Precision floating point summation, 1-primed.

       @arg xs: Iterable, list, tuple, etc. of values (C{scalar}s).

       @return: Precision L{fsum} (C{float}).

       @see: Function C{fsum}.
    '''
    def _xs1(xs):
        yield _1_0
        for x in xs:
            yield x
        yield _N_1_0

    return fsum(_xs1(xs)) if xs else _0_0


def fsum1_(*xs):
    '''Precision floating point summation of a few arguments, 1-primed.

       @arg xs: Values to be added (C{scalar}s), all positional.

       @return: Precision L{fsum} (C{float}).

       @see: Function C{fsum}
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
