
# -*- coding: utf-8 -*-

u'''Class L{Fsum} for precision floating point summation and I{running}
summation based on, respectively similar to Python's C{math.fsum}.

Generally, an L{Fsum} instance is considered a C{float} plus a small or zero
C{residual} value, see property L{Fsum.residual}.  However, there are several
C{integer} L{Fsum} cases, for example the result of C{ceil}, C{floor},
C{Fsum.__floordiv__} and methods L{Fsum.fint} and L{Fsum.fint2}.

Also, L{Fsum} methods L{Fsum.pow}, L{Fsum.__ipow__}, L{Fsum.__pow__} and
L{Fsum.__rpow__} return a (very long) C{int} if invoked with optional argument
C{mod} set to C{None}.  The C{residual} of an C{integer} L{Fsum} may be between
C{-1.0} and C{+1.0}, including C{INT0} if considered to be I{exact}.

Set env variable C{PYGEODESY_FSUM_PARTIALS} to an empty string (or anything
other than C{"fsum"}) for backward compatible summation of L{Fsum} partials.

Set env variable C{PYGEODESY_FSUM_RESIDUAL} to a C{float} string greater
than C{"0.0"} as the threshold to throw a L{ResidualError} in division or
exponention of an L{Fsum} instance with a I{relative} C{residual} exceeding
the threshold, see methods L{Fsum.RESIDUAL}, L{Fsum.pow}, L{Fsum.__ipow__}
and L{Fsum.__itruediv__}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import iscomplex, isint, isscalar, signOf, _signOf
from pygeodesy.constants import INT0, _isfinite, isinf, isnan, _pos_self, \
                               _0_0, _1_0, _N_1_0,  Float, Int
from pygeodesy.errors import itemsorted, _OverflowError, _TypeError, \
                            _ValueError, _xError2, _xkwds_get, _xkwds_get_, \
                            _ZeroDivisionError
from pygeodesy.interns import NN, _arg_, _COMMASPACE_, _DASH_, _EQUAL_, \
                             _exceeds_, _from_, _iadd_op_, _LANGLE_, \
                             _negative_, _NOTEQUAL_, _not_finite_, \
                             _not_scalar_, _PERCENT_, _PLUS_, _R_, _RANGLE_, \
                             _SLASH_, _SPACE_, _STAR_, _UNDER_
from pygeodesy.lazily import _ALL_LAZY, _getenv, _sys_version_info2
from pygeodesy.named import _Named, _NamedTuple, _NotImplemented, Fmt, unstr
from pygeodesy.props import _allPropertiesOf_n, deprecated_property_RO, \
                             Property_RO, property_RO
# from pygeodesy.streprs import Fmt, unstr  # from .named
# from pygeodesy.units import Float, Int  # from .constants

from math import ceil as _ceil, fabs, floor as _floor  # PYCHOK used! .ltp

__all__ = _ALL_LAZY.fsums
__version__ = '23.08.27'

_add_op_       = _PLUS_  # in .auxilats.auxAngle
_eq_op_        = _EQUAL_ * 2  # _DEQUAL_
_COMMASPACE_R_ = _COMMASPACE_ + _R_
_exceeds_R_    = _SPACE_ + _exceeds_(_R_)
_floordiv_op_  = _SLASH_ * 2  # _DSLASH_
_fset_op_      = _EQUAL_
_ge_op_        = _RANGLE_ + _EQUAL_
_gt_op_        = _RANGLE_
_integer_      = 'integer'
_le_op_        = _LANGLE_ + _EQUAL_
_lt_op_        = _LANGLE_
_mod_op_       = _PERCENT_
_mul_op_       = _STAR_
_ne_op_        = _NOTEQUAL_
_non_zero_     = 'non-zero'
_pow_op_       = _STAR_ * 2  # _DSTAR_, in .fmath
_sub_op_       = _DASH_      # in .auxilats.auxAngle, .fsums
_truediv_op_   = _SLASH_
_divmod_op_    = _floordiv_op_ + _mod_op_
_isub_op_      = _sub_op_ + _fset_op_  # in .auxilats.auxAngle, .fsums


def _2float(index=None, **name_value):  # in .fmath, .fstats
    '''(INTERNAL) Raise C{TypeError} or C{ValueError} if not scalar or infinite.
    '''
    n, v = name_value.popitem()  # _xkwds_popitem(name_value)
    try:
        v = float(v)
        if _isfinite(v):
            return v
        E, t = _ValueError, _not_finite_
    except Exception as e:
        E, t = _xError2(e)
    if index is not None:
        n = Fmt.SQUARE(n, index)
    raise E(n, v, txt=t)


def _2floats(xs, origin=0, sub=False):
    '''(INTERNAL) Yield each B{C{xs}} as a C{float}.
    '''
    try:
        i, x  =  origin, None
        _fin  = _isfinite
        _Fsum =  Fsum
        for x in xs:
            if isinstance(x, _Fsum):
                for p in x._ps:
                    yield (-p) if sub else p
            else:
                f = float(x)
                if not _fin(f):
                    raise ValueError(_not_finite_)
                if f:
                    yield (-f) if sub else f
            i += 1
    except Exception as e:
        E, t = _xError2(e)
        n = Fmt.SQUARE(xs=i)
        raise E(n, x, txt=t)


def _Powers(power, xs, origin=1):  # in .fmath
    '''(INTERNAL) Yield each C{xs} as C{float(x**power)}.
    '''
    if not isscalar(power):
        raise _TypeError(power=power, txt=_not_scalar_)
    try:
        i, x  =  origin, None
        _fin  = _isfinite
        _Fsum =  Fsum
        _pow  =  pow  # XXX math.pow
        for x in xs:
            if isinstance(x, _Fsum):
                P = x.pow(power)
                for p in P._ps:
                    yield p
            else:
                p = _pow(float(x), power)
                if not _fin(p):
                    raise ValueError(_not_finite_)
                yield p
            i += 1
    except Exception as e:
        E, t = _xError2(e)
        n = Fmt.SQUARE(xs=i)
        raise E(n, x, txt=t)


def _1primed(xs):
    '''(INTERNAL) 1-Prime the summation of C{xs}
       arguments I{known} to be C{finite float}.
    '''
    yield _1_0
    for x in xs:
        if x:
            yield x
    yield _N_1_0


def _psum(ps):  # PYCHOK used!
    '''(INTERNAL) Partials summation updating C{ps}, I{overridden below}.
    '''
    i   =  len(ps) - 1  # len(ps) > 2
    s   =  ps[i]
    _2s = _2sum
    while i > 0:
        i -= 1
        s, r = _2s(s, ps[i])
        if r:  # sum(ps) became inexact
            ps[i:] = [s, r] if s else [r]
            if i > 0:
                p = ps[i-1]  # round half-even
                if (p > 0 and r > 0) or \
                   (p < 0 and r < 0):  # signs match
                    r *= 2
                    t  = s + r
                    if r == (t - s):
                        s = t
            break
        ps[i:] = [s]
    return s


def _2scalar(other, _raiser=None):
    '''(INTERNAL) Return B{C{other}} as C{int}, C{float} or C{as-is}.
    '''
    if isinstance(other, Fsum):
        s, r = other._fint2
        if r:
            s, r = other._fprs2
            if r:  # PYCHOK no cover
                if _raiser and _raiser(r, s):
                    raise ValueError(_stresidual(_non_zero_, r))
                s = other  # L{Fsum} as-is
    else:
        s = other  # C{type} as-is
        if isint(s, both=True):
            s = int(s)
    return s


def _strcomplex(s, *args):
    '''(INTERNAL) C{Complex} 2- or 3-arg C{pow} error C{str}.
    '''
    c =  iscomplex.__name__[2:]
    n = _DASH_(len(args), _arg_)
    t = _SPACE_(c, s, _from_, n, pow.__name__)
    return unstr(t, *args)


def _stresidual(prefix, residual, **name_values):
    '''(INTERNAL) Residual error C{str}.
    '''
    p = _SPACE_(prefix, Fsum.residual.name)
    t =  Fmt.PARENSPACED(p, Fmt(residual))
    for n, v in itemsorted(name_values):
        n =  n.replace(_UNDER_, _SPACE_)
        p =  Fmt.PARENSPACED(n, Fmt(v))
        t = _COMMASPACE_(t, p)
    return t


def _2sum(a, b):  # by .testFmath
    '''(INTERNAL) Return C{a + b} as 2-tuple (sum, residual).
    '''
    s = a + b
    if not _isfinite(s):
        u = unstr(_2sum.__name__, a, b)
        t = Fmt.PARENSPACED(_not_finite_, s)
        raise _OverflowError(u, txt=t)
    if fabs(a) < fabs(b):
        a, b = b, a
    return s, (b - (s - a))


class Fsum(_Named):  # sync __methods__ with .vector3dBase.Vector3dBase
    '''Precision floating point I{running} summation.

       Unlike Python's C{math.fsum}, this class accumulates values and provides intermediate,
       I{running} precision floating point summation.  Accumulation may continue after
       intermediate, I{running} summuation.

       @note: Accumulated values may be L{Fsum} or C{scalar} instances with C{scalar} meaning
              type C{float}, C{int} or any C{type} convertible to a single C{float}, having
              method C{__float__}.

       @note: Handling of exceptions and C{inf}, C{INF}, C{nan} and C{NAN} differs from
              Python's C{math.fsum}.

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
#   _px        = 0
    _ratio     = None
    _RESIDUAL  = max(float(_getenv('PYGEODESY_FSUM_RESIDUAL', _0_0)), _0_0)

    def __init__(self, *xs, **name_RESIDUAL):
        '''New L{Fsum} for precision floating point I{running} summation.

           @arg xs: No, one or more initial values (each C{scalar} or an
                    L{Fsum} instance).
           @kwarg name_RESIDUAL: Optional C{B{name}=NN} for this L{Fsum}
                                 (C{str}) and C{B{RESIDUAL}=None} for the
                                 L{ResidualError} threshold.

           @see: Methods L{Fsum.fadd} and L{Fsum.RESIDUAL}.
        '''
        if name_RESIDUAL:
            n, r = _xkwds_get_(name_RESIDUAL, name=NN, RESIDUAL=None)
            if n:  # set name ...
                self.name = n
            if r is not None:
                self.RESIDUAL(r)  # ... for ResidualError
#       self._n  = 0
        self._ps = []  # [_0_0], see L{Fsum._fprs}
        if len(xs) > 1:
            self._facc(_2floats(xs, origin=1), up=False)  # PYCHOK yield
        elif xs:  # len(xs) == 1
            self._ps = [_2float(x=xs[0])]
            self._n  = 1

    def __abs__(self):
        '''Return this instance' absolute value as an L{Fsum}.
        '''
        s = _fsum(self._ps_1())  # == self._cmp_0(0, ...)
        return self._copy_n(self.__abs__) if s < 0 else \
               self._copy_2(self.__abs__)

    def __add__(self, other):
        '''Return the C{Fsum(B{self}, B{other})}.

           @arg other: An L{Fsum} or C{scalar}.

           @return: The sum (L{Fsum}).

           @see: Method L{Fsum.__iadd__}.
        '''
        f = self._copy_2(self.__add__)
        return f._fadd(other, _add_op_)

    def __bool__(self):  # PYCHOK not special in Python 2-
        '''Return C{True} if this instance is I{exactly} non-zero.
        '''
        s, r = self._fprs2
        return bool(s or r) and s != -r  # == self != 0

    def __ceil__(self):  # PYCHOK not special in Python 2-
        '''Return this instance' C{math.ceil} as C{int} or C{float}.

           @return: An C{int} in Python 3+, but C{float} in Python 2-.

           @see: Methods L{Fsum.__floor__} and property L{Fsum.ceil}.
        '''
        return self.ceil

    def __cmp__(self, other):  # Python 2-
        '''Compare this with an other instance or C{scalar}.

           @return: -1, 0 or +1 (C{int}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        s = self._cmp_0(other, self.cmp.__name__)
        return _signOf(s, 0)

    cmp = __cmp__

    def __divmod__(self, other):
        '''Return C{divmod(B{self}, B{other})} as 2-tuple C{(quotient,
           remainder)}, an C{int} in Python 3+ or C{float} in Python 2-
           and an L{Fsum}.

           @arg other: An L{Fsum} or C{scalar} modulus.

           @see: Method L{Fsum.__itruediv__}.
        '''
        f = self._copy_2(self.__divmod__)
        return f._fdivmod2(other, _divmod_op_)

    def __eq__(self, other):
        '''Compare this with an other instance or C{scalar}.
        '''
        return self._cmp_0(other, _eq_op_) == 0

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
        f = self._copy_2(self.__floordiv__)
        return f._floordiv(other, _floordiv_op_)

    def __format__(self, *other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, *other)

    def __ge__(self, other):
        '''Compare this with an other instance or C{scalar}.
        '''
        return self._cmp_0(other, _ge_op_) >= 0

    def __gt__(self, other):
        '''Compare this with an other instance or C{scalar}.
        '''
        return self._cmp_0(other, _gt_op_) > 0

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
        return self._fadd(other, _iadd_op_)

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
        return self._floordiv(other, _floordiv_op_ + _fset_op_)

    def __imatmul__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __imod__(self, other):
        '''Apply C{B{self} %= B{other}} to this instance.

           @arg other: An L{Fsum} or C{scalar} modulus.

           @return: This instance, updated (L{Fsum}).

           @see: Method L{Fsum.__divmod__}.
        '''
        self._fdivmod2(other, _mod_op_ + _fset_op_)
        return self

    def __imul__(self, other):
        '''Apply C{B{self} *= B{other}} to this instance.

           @arg other: An L{Fsum} or C{scalar} factor.

           @return: This instance, updated (L{Fsum}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Invalid B{C{other}} type.

           @raise ValueError: Invalid or non-finite B{C{other}}.
        '''
        return self._fmul(other, _mul_op_ + _fset_op_)

    def __int__(self):
        '''Return this instance as an C{int}.

           @see: Methods L{Fsum.int_float}, L{Fsum.__ceil__}
                 and L{Fsum.__floor__} and properties
                 L{Fsum.ceil} and L{Fsum.floor}.
        '''
        i, _ = self._fint2
        return i

    def __invert__(self):  # PYCHOK no cover
        '''Not implemented.'''
        # Luciano Ramalho, "Fluent Python", O'Reilly, 2nd Ed, 2022 p. 567
        return _NotImplemented(self)

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
        return self._fpow(other, _pow_op_ + _fset_op_, *mod)

    def __isub__(self, other):
        '''Apply C{B{self} -= B{other}} to this instance.

           @arg other: An L{Fsum} or C{scalar}.

           @return: This instance, updated (L{Fsum}).

           @raise TypeError: Invalid B{C{other}} type.

           @see: Method L{Fsum.fadd}.
        '''
        return self._fsub(other, _isub_op_)

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
        return self._ftruediv(other, _truediv_op_ + _fset_op_)

    def __le__(self, other):
        '''Compare this with an other instance or C{scalar}.
        '''
        return self._cmp_0(other, _le_op_) <= 0

    def __len__(self):
        '''Return the number of values accumulated (C{int}).
        '''
        return self._n

    def __lt__(self, other):
        '''Compare this with an other instance or C{scalar}.
        '''
        return self._cmp_0(other, _lt_op_) < 0

    def __matmul__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __mod__(self, other):
        '''Return C{B{self} % B{other}} as an L{Fsum}.

           @see: Method L{Fsum.__imod__}.
        '''
        f = self._copy_2(self.__mod__)
        return f._fdivmod2(other, _mod_op_)[1]

    def __mul__(self, other):
        '''Return C{B{self} * B{other}} as an L{Fsum}.

           @see: Method L{Fsum.__imul__}.
        '''
        f = self._copy_2(self.__mul__)
        return f._fmul(other, _mul_op_)

    def __ne__(self, other):
        '''Compare this with an other instance or C{scalar}.
        '''
        return self._cmp_0(other, _ne_op_) != 0

    def __neg__(self):
        '''Return I{a copy of} this instance, negated.
        '''
        return self._copy_n(self.__neg__)

    def __pos__(self):
        '''Return this instance I{as-is}, like C{float.__pos__()}.
        '''
        return self if _pos_self else self._copy_2(self.__pos__)

    def __pow__(self, other, *mod):  # PYCHOK 2 vs 3 args
        '''Return C{B{self}**B{other}} as an L{Fsum}.

           @see: Method L{Fsum.__ipow__}.
        '''
        f = self._copy_2(self.__pow__)
        return f._fpow(other, _pow_op_, *mod)

    def __radd__(self, other):
        '''Return C{B{other} + B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__iadd__}.
        '''
        f = self._copy_r2(other, self.__radd__)
        return f._fadd(self, _add_op_)

    def __rdivmod__(self, other):
        '''Return C{divmod(B{other}, B{self})} as 2-tuple C{(quotient,
           remainder)}.

           @see: Method L{Fsum.__divmod__}.
        '''
        f = self._copy_r2(other, self.__rdivmod__)
        return f._fdivmod2(self, _divmod_op_)

#   def __repr__(self):
#       '''Return the default C{repr(this)}.
#       '''
#       return self.toRepr(lenc=True)

    def __rfloordiv__(self, other):
        '''Return C{B{other} // B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__ifloordiv__}.
        '''
        f = self._copy_r2(other, self.__rfloordiv__)
        return f._floordiv(self, _floordiv_op_)

    def __rmatmul__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __rmod__(self, other):
        '''Return C{B{other} % B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__imod__}.
        '''
        f = self._copy_r2(other, self.__rmod__)
        return f._fdivmod2(self, _mod_op_)[1]

    def __rmul__(self, other):
        '''Return C{B{other} * B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__imul__}.
        '''
        f = self._copy_r2(other, self.__rmul__)
        return f._fmul(self, _mul_op_)

    def __round__(self, ndigits=None):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, ndigits=ndigits)

    def __rpow__(self, other, *mod):
        '''Return C{B{other}**B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__ipow__}.
        '''
        f = self._copy_r2(other, self.__rpow__)
        return f._fpow(self, _pow_op_, *mod)

    def __rsub__(self, other):
        '''Return C{B{other} - B{self}} as L{Fsum}.

           @see: Method L{Fsum.__isub__}.
        '''
        f = self._copy_r2(other, self.__rsub__)
        return f._fsub(self, _sub_op_)

    def __rtruediv__(self, other):
        '''Return C{B{other} / B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__itruediv__}.
        '''
        f = self._copy_r2(other, self.__rtruediv__)
        return f._ftruediv(self, _truediv_op_)

    def __str__(self):
        '''Return the default C{str(self)}.
        '''
        return self.toStr(lenc=True)

    def __sub__(self, other):
        '''Return C{B{self} - B{other}} as an L{Fsum}.

           @arg other: An L{Fsum} or C{scalar}.

           @return: The difference (L{Fsum}).

           @see: Method L{Fsum.__isub__}.
        '''
        f = self._copy_2(self.__sub__)
        return f._fsub(other, _sub_op_)

    def __truediv__(self, other):
        '''Return C{B{self} / B{other}} as an L{Fsum}.

           @arg other: An L{Fsum} or C{scalar} divisor.

           @return: The quotient (L{Fsum}).

           @see: Method L{Fsum.__itruediv__}.
        '''
        f = self._copy_2(self.__truediv__)
        return f._ftruediv(other, _truediv_op_)

    __trunc__ = __int__

    if _sys_version_info2 < (3, 0):  # PYCHOK no cover
        # <https://docs.Python.org/2/library/operator.html#mapping-operators-to-functions>
        __div__     = __truediv__
        __idiv__    = __itruediv__
        __long__    = __int__
        __nonzero__ = __bool__
        __rdiv__    = __rtruediv__

    def as_integer_ratio(self):
        '''Return this instance as the ratio of 2 integers.

           @return: 2-Tuple C{(numerator, denominator)} both
                    C{int} and with positive C{denominator}.

           @see: Standard C{float.as_integer_ratio} in Python 3+.
        '''
        n, r = self._fint2
        if r:
            i, d = r.as_integer_ratio()
            n *= d
            n += i
        else:  # PYCHOK no cover
            d = 1
        return n, d

    @property_RO
    def ceil(self):
        '''Get this instance' C{ceil} value (C{int} in Python 3+,
           but C{float} in Python 2-).

           @note: The C{ceil} takes the C{residual} into account.

           @see: Method L{Fsum.int_float} and properties L{Fsum.floor},
                 L{Fsum.imag} and L{Fsum.real}.
        '''
        s, r = self._fprs2
        c = _ceil(s) + int(r) - 1
        while r > (c - s):  # (s + r) > c
            c += 1
        return c

    def _cmp_0(self, other, op):
        '''(INTERNAL) Return C{scalar(self - B{other})} for 0-comparison.
        '''
        if isscalar(other):
            if other:
                s = _fsum(self._ps_1(other))
            else:
                s, r = self._fprs2
                s = _signOf(s, -r)
        elif isinstance(other, Fsum):
            s = _fsum(self._ps_1(*other._ps))
        else:
            raise self._TypeError(op, other)  # txt=_invalid_
        return s

    def copy(self, deep=False, name=NN):
        '''Copy this instance, C{shallow} or B{C{deep}}.

           @return: The copy (L{Fsum}).
         '''
        f = _Named.copy(self, deep=deep, name=name)
        f._n  = self._n if deep else 1
        f._ps = list(self._ps)  # separate list
        return f

    def _copy_0(self, *xs):
        '''(INTERNAL) Copy with/-out overriding C{partials}.
        '''
        # for x in xs:
        #     assert isscalar(x)
        f = self._Fsum(self._n + len(xs), *xs)
        if self.name:
            f._name = self.name  # .rename calls _update_attrs
        return f

    def _copy_2(self, which):
        '''(INTERNAL) Copy for I{dyadic} operators.
        '''
        # NOT .classof due to .Fdot(a, *b) args, etc.
        f = _Named.copy(self, deep=False, name=which.__name__)
        # assert f._n == self._n
        f._ps = list(self._ps)  # separate list
        return f

    def _copy_n(self, which):
        '''(INTERNAL) Negated copy for I{monadic} C{__abs__} and C{__neg__}.
        '''
        if self._ps:
            f = self._Fsum(self._n)
            f._ps[:] = self._ps_n()
#           f._facc_up(up=False)
        else:
            f = self._Fsum(self._n, _0_0)
        f._name = which.__name__  # .rename calls _update_attrs
        return f

    def _copy_r2(self, other, which):
        '''(INTERNAL) Copy for I{reverse-dyadic} operators.
        '''
        return other._copy_2(which) if isinstance(other, Fsum) else \
               Fsum(other, name=which.__name__)  # see ._copy_2

    def _copy_RESIDUAL(self, other):
        '''(INTERNAL) Copy C{other._RESIDUAL}.
        '''
        R = other._RESIDUAL
        if R is not Fsum._RESIDUAL:
            self._RESIDUAL = R

    def _copy_up(self, _fprs2=False):
        '''(INTERNAL) Minimal, anonymous copy.
        '''
        f = self._Fsum(self._n, *self._ps)
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
        f = self._copy_2(self.divmod)
        return f._fdivmod2(other, _divmod_op_)

    def _Error(self, op, other, Error, **txt):
        '''(INTERNAL) Format an B{C{Error}} for C{{self} B{op} B{other}}.
        '''
        return Error(_SPACE_(self.toRepr(), op, repr(other)), **txt)

    def _ErrorX(self, X, xs, **kwds):  # in .fmath
        '''(INTERNAL) Format a caught exception.
        '''
        E, t = _xError2(X)
        n = unstr(self.named3, *xs[:3], _ELLIPSIS=len(xs) > 3, **kwds)
        return E(n, txt=t, cause=X)

    def _facc(self, xs, up=True):  # from .elliptic._Defer.Fsum
        '''(INTERNAL) Accumulate more known C{scalar}s.
        '''
        n, ps, _2s = 0, self._ps, _2sum
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
            # Fsum._px = max(Fsum._px, len(ps))
            if up:
                self._update()
        return self

    def _facc_(self, *xs, **up):
        '''(INTERNAL) Accumulate all positional C{scalar}s.
        '''
        return self._facc(xs, **up) if xs else self

#   def _facc_up(self, up=True):
#       '''(INTERNAL) Update the C{partials}, by removing
#          and re-accumulating the final C{partial}.
#       '''
#       while len(self._ps) > 1:
#           p = self._ps.pop()
#           if p:
#               n = self._n
#               self._facc_(p, up=False)
#               self._n = n
#               break
#       return self._update() if up else self  # ._fpsqz()

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
            self._facc(_2floats(xs))  # PYCHOK yield
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
        return self._facc(_2floats(xs, origin=1))  # PYCHOK yield

    def _fadd(self, other, op):  # in .fmath.Fhorner
        '''(INTERNAL) Apply C{B{self} += B{other}}.
        '''
        if isinstance(other, Fsum):
            if other is self:
                self._facc_(*other._ps)  # == ._facc(tuple(other._ps))
            elif other._ps:
                self._facc(other._ps)
        elif not isscalar(other):
            raise self._TypeError(op, other)  # txt=_invalid_
        elif other:
            self._facc_(other)
        return self

    fcopy   =   copy        # for backward compatibility
    fdiv    = __itruediv__  # for backward compatibility
    fdivmod = __divmod__    # for backward compatibility

    def _fdivmod2(self, other, op):
        '''(INTERNAL) C{divmod(B{self}, B{other})} as 2-tuple
           (C{int} or C{float}, remainder C{self}).
        '''
        # result mostly follows CPython function U{float_divmod
        # <https://GitHub.com/python/cpython/blob/main/Objects/floatobject.c>},
        # but at least divmod(-3, 2) equals Cpython's result (-2, 1).
        q = self._copy_up(_fprs2=True)._ftruediv(other, op).floor
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
            t = _stresidual(_integer_, r)
            raise ResidualError(_integer_, i, txt=t)
        n = name or self.fint.__name__
        return Fsum(name=n)._fset(i, asis=True)

    def fint2(self, **name):
        '''Return this instance' current running sum as C{int} and
           the I{integer} residual.

           @kwarg name: Optional name (C{str}).

           @return: An L{Fsum2Tuple}C{(fsum, residual)} with C{fsum}
                    an C{int} and I{integer} C{residual} a C{float} or
                    C{INT0} if the C{fsum} is considered to be I{exact}.
        '''
        return Fsum2Tuple(*self._fint2, **name)

    @Property_RO
    def _fint2(self):  # see ._fset
        '''(INTERNAL) Get 2-tuple (C{int}, I{integer} residual).
        '''
        i =  int(self._fprs)  # int(self)
        r = _fsum(self._ps_1(i)) if len(self._ps) > 1 else (
               (self._ps[0] - i) if self._ps else -i)
        return i, (r or INT0)

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
        if isscalar(other):
            f = self._mul_scalar(other, op)
        elif not isinstance(other, Fsum):
            raise self._TypeError(op, other)  # txt=_invalid_
        elif len(self._ps) != 1:
            f = self._mul_Fsum(other, op)
        elif len(other._ps) != 1:  # len(self._ps) == 1
            f = other._copy_up()._mul_scalar(self._ps[0], op)
        else:  # len(other._ps) == len(self._ps) == 1
            s = self._finite(self._ps[0] * other._ps[0])
            return self._fset(s, asis=True, n=len(self) + 1)
        return self._fset(f)

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
            i =  self._copy_0(self._fint2[0])  # assert _fint2[1] == 0
            x = _2scalar(other)  # C{int}, C{float} or other
            s =  i._pow_2(x, other, op) if isscalar(x) else i._fpow(x, op)
        else:  # pow(self, other) == pow(self, other, None)
            p = None
            if isinstance(other, Fsum):
                x, r = other._fprs2
                if r:
                    if self._raiser(r, x):
                        raise self._ResidualError(op, other, r)
                    p =  self._pow_scalar(r, other, op)
#                   p = _2scalar(p)  # _raiser = None
            elif not isscalar(other):
                raise self._TypeError(op, other)  # txt=_invalid_
            else:
                x = self._finite(other, op)
            s = self._pow_scalar(x, other, op)
            if p is not None:
                s *= p
        return self._fset(s, asis=isint(s), n=max(len(self), 1))

    @Property_RO
    def _fprs(self):
        '''(INTERNAL) Get and cache this instance' precision
           running sum (C{float} or C{int}), ignoring C{residual}.

           @note: The precision running C{fsum} after a C{//=} or
                  C{//} C{floor} division is C{int} in Python 3+.
        '''
        ps = self._ps
        n = len(ps) - 1
        if n > 1:
            s = _psum(ps)
        elif n > 0:  # len(ps) == 2
            s, p = _2sum(*ps) if ps[1] else ps
            ps[:] = ([p, s] if s else [p]) if p else [s]
        elif n < 0:  # see L{Fsum.__init__}
            s = _0_0
            ps[:] = [s]
        else:  # len(ps) == 1
            s = ps[0]
        # assert self._ps is ps
        # assert Fsum._fprs2.name not in self.__dict__
        return s

    @Property_RO
    def _fprs2(self):
        '''(INTERNAL) Get and cache this instance' precision
           running sum and residual (L{Fsum2Tuple}).
        '''
        s =  self._fprs
        r = _fsum(self._ps_1(s)) if len(self._ps) > 1 else INT0
        return Fsum2Tuple(s, r)  # name=Fsum.fsum2.__name__

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
            pass  # from ._fmul, ._ftruediv and ._pow_scalar
        elif isinstance(other, Fsum):
            self._n     = other._n
            self._ps[:] = other._ps
            self._copy_RESIDUAL(other)
            # use or zap the C{Property_RO} values
            Fsum._fint2._update_from(self, other)
            Fsum._fprs ._update_from(self, other)
            Fsum._fprs2._update_from(self, other)
        elif isscalar(other):
            s = other if asis else float(other)
            i = int(s)  # see ._fint2
            t = i, ((s - i) or INT0)
            self._n     =  n
            self._ps[:] = [s]
            # Property_RO _fint2, _fprs and _fprs2 can't be a Property:
            # Property's _fset zaps the value just set by the @setter
            self.__dict__.update(_fint2=t, _fprs=s, _fprs2=Fsum2Tuple(s, INT0))
        else:  # PYCHOK no cover
            raise self._TypeError(_fset_op_, other)  # txt=_invalid_
        return self

    def fsub(self, xs=()):
        '''Subtract an iterable of C{scalar} or L{Fsum} instances
           from this instance.

           @arg xs: Iterable, list, tuple. etc. (C{scalar}
                    or L{Fsum} instances).

           @return: This instance, updated (L{Fsum}).

           @see: Method L{Fsum.fadd}.
        '''
        return self._facc(_2floats(xs, sub=True)) if xs else self  # PYCHOK yield

    def fsub_(self, *xs):
        '''Subtract all positional C{scalar} or L{Fsum} instances
           from this instance.

           @arg xs: Values to subtract (C{scalar} or
                    L{Fsum} instances), all positional.

           @return: This instance, updated (L{Fsum}).

           @see: Method L{Fsum.fadd}.
        '''
        return self._facc(_2floats(xs, origin=1, sub=True)) if xs else self  # PYCHOK yield

    def _fsub(self, other, op):
        '''(INTERNAL) Apply C{B{self} -= B{other}}.
        '''
        if isinstance(other, Fsum):
            if other is self:  # or other._fprs2 == self._fprs2:
                self._fset(_0_0, asis=True, n=len(self) * 2)  # self -= self
            elif other._ps:
                self._facc(other._ps_n())
        elif not isscalar(other):
            raise self._TypeError(op, other)  # txt=_invalid_
        elif self._finite(other, op):
            self._facc_(-other)
        return self

    def _Fsum(self, n, *ps):
        '''(INTERNAL) New L{Fsum} instance.
        '''
        f = Fsum()
        f._n = n
        if ps:
            f._ps[:] = ps
        f._copy_RESIDUAL(self)
        return f

    def fsum(self, xs=()):
        '''Add more C{scalar} or L{Fsum} instances and summate.

           @kwarg xs: Iterable, list, tuple, etc. (C{scalar} or
                      L{Fsum} instances).

           @return: Precision running sum (C{float} or C{int}).

           @see: Method L{Fsum.fadd}.

           @note: Accumulation can continue after summation.
        '''
        f = self._facc(_2floats(xs)) if xs else self  # PYCHOK yield
        return f._fprs

    def fsum_(self, *xs):
        '''Add all positional C{scalar} or L{Fsum} instances and summate.

           @arg xs: Values to add (C{scalar} or L{Fsum} instances),
                    all positional.

           @return: Precision running sum (C{float} or C{int}).

           @see: Method L{Fsum.fsum}.
        '''
        f = self._facc(_2floats(xs, origin=1)) if xs else self  # PYCHOK yield
        return f._fprs

    def fsum2(self, xs=(), **name):
        '''Add more C{scalar} or L{Fsum} instances and return the
           current precision running sum and the C{residual}.

           @kwarg xs: Iterable, list, tuple, etc. (C{scalar} or
                      L{Fsum} instances).
           @kwarg name: Optional name (C{str}).

           @return: L{Fsum2Tuple}C{(fsum, residual)} with C{fsum} the
                    current precision running sum and C{residual}, the
                    (precision) sum of the remaining C{partials}.  The
                    C{residual is INT0} if the C{fsum} is considered
                    to be I{exact}.

           @see: Methods L{Fsum.fint2}, L{Fsum.fsum} and L{Fsum.fsum2_}
        '''
        f = self._facc(_2floats(xs)) if xs else self  # PYCHOK yield
        t = f._fprs2
        if name:
            t = t.dup(name=_xkwds_get(name, name=NN))
        return t

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
        if xs:
            s, t = self._facc(_2floats(xs, origin=1))._fprs2  # PYCHOK yield
            return s, _fsum((s, -p, r, -t))  # ((s - p) + (r - t))
        else:  # PYCHOK no cover
            return p, _0_0

#   ftruediv = __itruediv__   # for naming consistency

    def _ftruediv(self, other, op):
        '''(INTERNAL) Apply C{B{self} /= B{other}}.
        '''
        n = _1_0
        if isinstance(other, Fsum):
            if other is self or other._fprs2 == self._fprs2:
                return self._fset(_1_0, asis=True, n=len(self))
            d, r = other._fprs2
            if r:
                if not d:  # PYCHOK no cover
                    d = r
                elif self._raiser(r, d):
                    raise self._ResidualError(op, other, r)
                else:
                    d, n = other.as_integer_ratio()
        elif isscalar(other):
            d = other
        else:  # PYCHOK no cover
            raise self._TypeError(op, other)  # txt=_invalid_
        try:
            s = 0 if isinf(d) else (
                d if isnan(d) else self._finite(n / d))
        except Exception as x:
            E, t = _xError2(x)
            raise self._Error(op, other, E, txt=t)
        f = self._mul_scalar(s, _mul_op_)  # handles 0, NAN, etc.
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

           @return: This C{integer} sum if this instance C{is_integer},
                    otherwise return the C{float} sum if the residual
                    is zero or if C{B{raiser}=False}.

           @raise ResidualError: Non-zero residual and C{B{raiser}=True}.

           @see: Methods L{Fsum.fint} and L{Fsum.fint2}.
        '''
        s, r = self._fint2
        if r:
            s, r = self._fprs2
            if r and raiser:  # PYCHOK no cover
                t = _stresidual(_non_zero_, r)
                raise ResidualError(int_float=s, txt=t)
            s = float(s)  # redundant
        return s

    def is_exact(self):
        '''Is this instance' current running C{fsum} considered to
           be exact? (C{bool}).
        '''
        return self.residual is INT0

    def is_integer(self):
        '''Is this instance' current running sum C{integer}? (C{bool}).

           @see: Methods L{Fsum.fint} and L{Fsum.fint2}.
        '''
        _, r = self._fint2
        return False if r else True

    def is_math_fsum(self):
        '''Return whether functions L{fsum}, L{fsum_}, L{fsum1}
           and L{fsum1_} plus partials summation are based on
           Python's C{math.fsum} or not.

           @return: C{2} if all functions and partials summation
                    are based on C{math.fsum}, C{True} if only
                    the functions are based on C{math.fsum} (and
                    partials summation is not) or C{False} if
                    none are.
        '''
        f = Fsum._math_fsum
        return 2 if _psum is f else bool(f)

    def _mul_Fsum(self, other, op=_mul_op_):
        '''(INTERNAL) Return C{B{self} * Fsum B{other}} as L{Fsum}.
        '''
        # assert isinstance(other, Fsum)
        return self._copy_0()._facc(self._ps_x(op, *other._ps), up=False)

    def _mul_scalar(self, factor, op):
        '''(INTERNAL) Return C{B{self} * scalar B{factor}} as L{Fsum} or C{0}.
        '''
        # assert isscalar(factor)
        if self._finite(factor, op) and self._ps:
            if factor == _1_0:
                return self
            f = self._copy_0()._facc(self._ps_x(op, factor), up=False)
        else:
            f = self._copy_0(_0_0)
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
        f = self._copy_2(self.pow)
        if f and isint(x) and x >= 0 and not mod:
            f._pow_int(x, x, _pow_op_)  # f **= x
        else:
            f._fpow(x, _pow_op_, *mod)  # f = pow(f, x, *mod)
        return f

    def _pow_0_1(self, x, other):
        '''(INTERNAL) Return B{C{self}**1} or C{B{self}**0 == 1.0}.
        '''
        return self if x else (1 if self.is_integer() and isint(other) else _1_0)

    def _pow_2(self, x, other, op):
        '''(INTERNAL) 2-arg C{pow(B{self}, scalar B{x})} embellishing errors.
        '''
        # assert len(self._ps) == 1 and isscalar(x)
        b = self._ps[0]  # assert isscalar(b)
        try:  # type(s) == type(x) if x in (_1_0, 1)
            s = pow(b, x)  # -1**2.3 == -(1**2.3)
            if not iscomplex(s):
                return self._finite(s)  # 0**INF == 0.0, 1**INF==1.0
            # neg**frac == complex in Python 3+, but ValueError in 2-
            E, t = _ValueError, _strcomplex(s, b, x)  # PYCHOK no cover
        except Exception as x:
            E, t = _xError2(x)
        raise self._Error(op, other, E, txt=t)

    def _pow_3(self, other, mod, op):
        '''(INTERNAL) 3-arg C{pow(B{self}, B{other}, int B{mod} or C{None})}.
        '''
        b, r = self._fprs2 if mod is None else self._fint2
        if r and self._raiser(r, b):
            t = _non_zero_ if mod is None else _integer_
            E, t = ResidualError, _stresidual(t, r, mod=mod)
        else:
            try:  # b, other, mod all C{int}, unless C{mod} is C{None}
                x = _2scalar(other, _raiser=self._raiser)
                s =  pow(b, x, mod)
                if not iscomplex(s):
                    return self._finite(s)
                # neg**frac == complex in Python 3+, but ValueError in 2-
                E, t = _ValueError, _strcomplex(s, b, x, mod)  # PYCHOK no cover
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
                p = self._copy_up()
                m = 1  # single-bit mask
                if x & m:
                    x -= m  # x ^= m
                    f  = p._copy_up()
                else:
                    f = self._copy_0(_1_0)
                while x:
                    p  = p._mul_Fsum(p, op)  # p **= 2
                    m += m  # m <<= 1
                    if x & m:
                        x -= m  # x ^= m
                        f  = f._mul_Fsum(p, op)  # f *= p
            elif x > 1:  # self**2
                f = self._mul_Fsum(self, op)
            else:  # self**1 or self**0
                f = self._pow_0_1(x, other)
        elif self._ps:  # self._ps[0]**x
            f = self._pow_2(x, other, op)
        else:  # PYCHOK no cover
            # 0**pos_int == 0, but 0**0 == 1
            f = 0 if x else 1  # like ._fprs
        return self._fset(f, asis=isint(f), n=len(self))

    def _pow_scalar(self, x, other, op):
        '''(INTERNAL) Return C{self**B{x}} for C{scalar B{x}}.
        '''
        s, r = self._fprs2
        if isint(x, both=True):
            x = int(x)  # Fsum**int
            y = abs(x)
            if y > 1:
                if r:
                    f = self._copy_up()._pow_int(y, other, op)
                    if x > 0:  # > 1
                        return f
                    # assert x < 0  # < -1
                    s, r = f._fprs2
                    if r:
                        return self._copy_0(_1_0)._ftruediv(f, op)
                    # use **= -1 for the CPython float_pow
                    # error if s is zero, and not s = 1 / s
                    x = -1
#           elif y > 1:  # self**2 or self**-2
#               f = self._mul_Fsum(self, op)
#               if x < 0:
#                   f = f._copy_0(_1_0)._ftruediv(f, op)
#               return f
            elif x < 0:  # self**-1 == 1 / self
                if r:
                    return self._copy_0(_1_0)._ftruediv(self, op)
            else:  # self**1 or self**0
                return self._pow_0_1(x, other)  # self or 0.0
        elif not isscalar(x):  # assert ...
            raise self._TypeError(op, other, txt=_not_scalar_)
        elif r and self._raiser(r, s):  # non-zero residual**fractional
            # raise self._ResidualError(op, other, r, fractional_power=x)
            t = _stresidual(_non_zero_, r, fractional_power=x)
            raise self._Error(op, other, ResidualError, txt=t)
        # assert isscalar(s) and isscalar(x)
        return self._copy_0(s)._pow_2(x, other, op)

    def _ps_1(self, *less):
        '''(INTERNAL) Yield partials, 1-primed and subtract any C{less}.
        '''
        yield _1_0
        for p in self._ps:
            if p:
                yield p
        for p in less:
            if p:
                yield -p
        yield _N_1_0

    def _ps_n(self):
        '''(INTERNAL) Yield partials, negated.
        '''
        for p in self._ps:
            if p:
                yield -p

    def _ps_x(self, op, *factors):  # see .fmath.Fhorner
        '''(INTERNAL) Yield all C{partials} times each B{C{factor}},
           in total, up to C{len(partials) * len(factors)} items.
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
        '''Get this instance' residual (C{float} or C{int}), the
           C{sum(partials)} less the precision running sum C{fsum}.

           @note: If the C{residual is INT0}, the precision running
                  C{fsum} is considered to be I{exact}.

           @see: Methods L{Fsum.fsum}, L{Fsum.fsum2} and L{Fsum.is_exact}.
        '''
        return self._fprs2.residual

    def _raiser(self, r, s):
        '''(INTERNAL) Does the ratio C{r / s} exceed threshold?
        '''
        self._ratio = t = fabs((r / s) if s else r)
        return t > self._RESIDUAL

    def RESIDUAL(self, *threshold):
        '''Get and set this instance' I{ratio} for raising L{ResidualError}s,
           overriding the default from env variable C{PYGEODESY_FSUM_RESIDUAL}.

           @arg threshold: If C{scalar}, the I{ratio} to exceed for raising
                           L{ResidualError}s in division and exponention, if
                           C{None} restore the default set with env variable
                           C{PYGEODESY_FSUM_RESIDUAL} or if omitted, keep the
                           current setting.

           @return: The previous C{RESIDUAL} setting (C{float}).

           @raise ValueError: Negative B{C{threshold}}.

           @note: A L{ResidualError} is thrown if the non-zero I{ratio}
                  C{residual} / C{fsum} exceeds the B{C{threshold}}.
        '''
        r = self._RESIDUAL
        if threshold:
            t = threshold[0]
            t = Fsum._RESIDUAL if t is None else (
                      float(t) if isscalar(t) else (  # for backward ...
                          _0_0 if bool(t) else _1_0))  # ... compatibility
            if t < 0:
                u = self._unstr(self.RESIDUAL, *threshold)
                raise _ValueError(u, RESIDUAL=t, txt=_negative_)
            self._RESIDUAL = t
        return r

    def _ResidualError(self, op, other, residual):
        '''(INTERNAL) Non-zero B{C{residual}} etc.
        '''
        t = _stresidual(_non_zero_, residual, ratio=self._ratio,
                                            RESIDUAL=self._RESIDUAL)
        t =  t.replace(_COMMASPACE_R_, _exceeds_R_)
        return self._Error(op, other, ResidualError, txt=t)

    def signOf(self, res=True):
        '''Determine the sign of this instance.

           @kwarg res: If C{True} consider, otherwise
                       ignore the residual (C{bool}).

           @return: The sign (C{int}, -1, 0 or +1).
        '''
        s, r = self._fprs2 if res else (self._fprs, 0)
        return _signOf(s, -r)

    def toRepr(self, **prec_sep_fmt_lenc):  # PYCHOK signature
        '''Return this C{Fsum} instance as representation.

           @kwarg prec_sep_fmt_lenc: Optional keyword arguments for
                       method L{Fsum2Tuple.toRepr} plus C{B{lenc}=True}
                       (C{bool}) to in-/exclude the current C{[len]}
                       of this L{Fsum} enclosed in I{[brackets]}.

           @return: This instance (C{repr}).
        '''
        return self._toT(self._fprs2.toRepr, **prec_sep_fmt_lenc)

    def toStr(self, **prec_sep_fmt_lenc):  # PYCHOK signature
        '''Return this C{Fsum} instance as string.

           @kwarg prec_sep_fmt_lenc: Optional keyword arguments for
                       method L{Fsum2Tuple.toStr} plus C{B{lenc}=True}
                       (C{bool}) to in-/exclude the current C{[len]}
                       of this L{Fsum} enclosed in I{[brackets]}.

           @return: This instance (C{str}).
        '''
        return self._toT(self._fprs2.toStr, **prec_sep_fmt_lenc)

    def _toT(self, toT, fmt=Fmt.g, lenc=True, **kwds):
        '''(INTERNAL) Helper for C{toRepr} and C{toStr}.
        '''
        n = self.named3
        if lenc:
            n = Fmt.SQUARE(n, len(self))
        return _SPACE_(n, toT(fmt=fmt, **kwds))

    def _TypeError(self, op, other, **txt):  # PYCHOK no cover
        '''(INTERNAL) Return a C{TypeError}.
        '''
        return self._Error(op, other, _TypeError, **txt)

    def _update(self):  # see ._fset
        '''(INTERNAL) Zap all cached C{Property_RO} values.
        '''
        Fsum._fint2._update(self)
        Fsum._fprs ._update(self)
        Fsum._fprs2._update(self)
        return self

    def _ValueError(self, op, other, **txt):  # PYCHOK no cover
        '''(INTERNAL) Return a C{ValueError}.
        '''
        return self._Error(op, other, _ValueError, **txt)

    def _ZeroDivisionError(self, op, other, **txt):  # PYCHOK no cover
        '''(INTERNAL) Return a C{ZeroDivisionError}.
        '''
        return self._Error(op, other, _ZeroDivisionError, **txt)

_allPropertiesOf_n(3, Fsum, Property_RO)  # PYCHOK assert, see Fsum._fset, -._update


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
              to be I{exact}, see method L{Fsum2Tuple.is_exact}.
    '''
    _Names_ = ( Fsum.fsum.__name__, Fsum.residual.name)
    _Units_ = (_Float_Int,         _Float_Int)

    @Property_RO
    def Fsum(self):
        '''Get this L{Fsum2Tuple} as an L{Fsum}.
        '''
        f = Fsum(name=self.name)
        return f._copy_0(*(s for s in reversed(self) if s))

    def is_exact(self):
        '''Is this L{Fsum2Tuple} considered to be exact? (C{bool}).
        '''
        return self.Fsum.is_exact()

    def is_integer(self):
        '''Is this L{Fsum2Tuple} C{integer}? (C{bool}).
        '''
        return self.Fsum.is_integer()


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

    Fsum._math_fsum = _sum = _fsum  # PYCHOK exported

    if _getenv('PYGEODESY_FSUM_PARTIALS', _fsum.__name__) == _fsum.__name__:
        _psum = _fsum  # PYCHOK redef

except ImportError:
    _sum = sum  # Fsum(NAN) exception fall-back

    def _fsum(xs):
        '''(INTERNAL) Precision summation, Python 2.5-.
        '''
        return Fsum(name=_fsum.__name__)._facc(xs, up=False)._fprs


def fsum(xs, floats=False):
    '''Precision floating point summation based on or like Python's C{math.fsum}.

       @arg xs: Iterable, list, tuple, etc. of values (C{scalar} or
                L{Fsum} instances).
       @kwarg floats: Optionally, use C{B{floats}=True} iff I{all}
                      B{C{xs}} are known to be C{float}.

       @return: Precision C{fsum} (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise TypeError: Non-scalar B{C{xs}} value.

       @raise ValueError: Invalid or non-finite B{C{xs}} value.

       @note: Exceptions and I{non-finite} handling may differ if not
              based on Python's C{math.fsum}.

       @see: Class L{Fsum} and methods L{Fsum.fsum} and L{Fsum.fadd}.
    '''
    return _fsum(xs if floats else _2floats(xs)) if xs else _0_0  # PYCHOK yield


def fsum_(*xs, **floats):
    '''Precision floating point summation of all positional arguments.

       @arg xs: Values to be added (C{scalar} or L{Fsum} instances),
                all positional.
       @kwarg floats: Optionally, use C{B{floats}=True} iff I{all}
                      B{C{xs}} are known to be C{float}.

       @return: Precision C{fsum} (C{float}).

       @see: Function C{fsum}.
    '''
    return _fsum(xs if _xkwds_get(floats, floats=False) else
                _2floats(xs, origin=1)) if xs else _0_0  # PYCHOK yield


def fsumf_(*xs):
    '''Precision floating point summation L{fsum_}C{(*xs, floats=True)}.
    '''
    return _fsum(xs) if xs else _0_0


def fsum1(xs, floats=False):
    '''Precision floating point summation of a few arguments, 1-primed.

       @arg xs: Iterable, list, tuple, etc. of values (C{scalar} or
                L{Fsum} instances).
       @kwarg floats: Optionally, use C{B{floats}=True} iff I{all}
                      B{C{xs}} are known to be C{float}.

       @return: Precision C{fsum} (C{float}).

       @see: Function C{fsum}.
    '''
    return _fsum(_1primed(xs if floats else _2floats(xs))) if xs else _0_0  # PYCHOK yield


def fsum1_(*xs, **floats):
    '''Precision floating point summation of a few arguments, 1-primed.

       @arg xs: Values to be added (C{scalar} or L{Fsum} instances),
                all positional.
       @kwarg floats: Optionally, use C{B{floats}=True} iff I{all}
                      B{C{xs}} are known to be C{float}.

       @return: Precision C{fsum} (C{float}).

       @see: Function C{fsum}
    '''
    return _fsum(_1primed(xs if _xkwds_get(floats, floats=False) else
                         _2floats(xs, origin=1))) if xs else _0_0  # PYCHOK yield


def fsum1f_(*xs):
    '''Precision floating point summation L{fsum1_}C{(*xs, floats=True)}.
    '''
    return _fsum(_1primed(xs)) if xs else _0_0


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
