
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

Set env variable C{PYGEODESY_FSUM_PARTIALS} to string C{"fsum"}) for summation
of L{Fsum} partials by Python function C{math.fsum}.

Set env variable C{PYGEODESY_FSUM_RESIDUAL} to a C{float} string greater than
C{"0.0"} as the threshold to throw a L{ResidualError} in division or exponention
of an L{Fsum} instance with a I{relative} C{residual} exceeding the threshold,
see methods L{Fsum.RESIDUAL}, L{Fsum.pow}, L{Fsum.__ipow__} and L{Fsum.__itruediv__}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import iscomplex, isint, isscalar, itemsorted, \
                             signOf, _signOf
from pygeodesy.constants import INT0, _isfinite, isinf, isnan, NEG0, _pos_self, \
                               _0_0, _1_0, _N_1_0,  Float, Int
from pygeodesy.errors import _OverflowError, _TypeError, _ValueError, _xError, \
                             _xError2, _xkwds_get, _ZeroDivisionError
from pygeodesy.interns import NN, _arg_, _COMMASPACE_, _DASH_, _DOT_, _EQUAL_, \
                             _exceeds_, _from_, _iadd_op_, _LANGLE_, _negative_, \
                             _NOTEQUAL_, _not_finite_, _not_scalar_, _PERCENT_, \
                             _PLUS_, _R_, _RANGLE_, _SLASH_, _SPACE_, _STAR_, _UNDER_
from pygeodesy.lazily import _ALL_LAZY, _getenv, _sys_version_info2
from pygeodesy.named import _Named, _NamedTuple, _NotImplemented,  Fmt, unstr
from pygeodesy.props import _allPropertiesOf_n, deprecated_property_RO, \
                             Property_RO, property_RO
# from pygeodesy.streprs import Fmt, unstr  # from .named
# from pygeodesy.units import Float, Int  # from .constants

from math import ceil as _ceil, fabs, floor as _floor  # PYCHOK used! .ltp

__all__ = _ALL_LAZY.fsums
__version__ = '24.04.09'

_add_op_       = _PLUS_  # in .auxilats.auxAngle
_eq_op_        = _EQUAL_ * 2  # _DEQUAL_
_COMMASPACE_R_ = _COMMASPACE_ + _R_
_div_          = 'div'
_exceeds_R_    = _SPACE_ + _exceeds_(_R_)
_floordiv_op_  = _SLASH_ * 2  # _DSLASH_
_fset_op_      = _EQUAL_
_ge_op_        = _RANGLE_ + _EQUAL_
_gt_op_        = _RANGLE_
_integer_      = 'integer'
_le_op_        = _LANGLE_ + _EQUAL_
_lt_op_        = _LANGLE_
_mod_          = 'mod'
_mod_op_       = _PERCENT_
_mul_op_       = _STAR_
_ne_op_        = _NOTEQUAL_
_non_zero_     = 'non-zero'
_pow_op_       = _STAR_ * 2  # _DSTAR_, in .fmath
_sub_op_       = _DASH_      # in .auxilats.auxAngle, .fsums
_truediv_op_   = _SLASH_
_divmod_op_    = _floordiv_op_ + _mod_op_
_isub_op_      = _sub_op_ + _fset_op_  # in .auxilats.auxAngle, .fsums


def _2delta(*ab):
    '''(INTERNAL) Helper for C{Fsum.fsum2f_}.
    '''
    try:
        a, b = _2sum(*ab)
    except _OverflowError:
        a, b =  ab
    return float(a if fabs(a) > fabs(b) else b)


def _2error(unused):
    '''(INTERNAL) Throw a C{not finite} exception.
    '''
    raise ValueError(_not_finite_)


def _2float(index=None, **name_value):  # in .fmath, .fstats
    '''(INTERNAL) Raise C{TypeError} or C{ValueError} if not scalar or infinite.
    '''
    n, v = name_value.popitem()  # _xkwds_item2(name_value)
    try:
        v = float(v)
        return v if _isfinite(v) else _2error(v)
    except Exception as X:
        raise _xError(X, Fmt.INDEX(n, index), v)


def _X_ps(X):  # for _2floats only
    return X._ps


def _2floats(xs, origin=0, _X=_X_ps, _x=float):
    '''(INTERNAL) Yield each B{C{xs}} as a C{float}.
    '''
    try:
        i, x =  origin, None
        _fin = _isfinite
        _Fs  =  Fsum
        for x in xs:
            if isinstance(x, _Fs):
                for p in _X(x):
                    yield p
            else:
                f = _x(x)
                yield f if _fin(f) else _2error(f)
            i += 1
    except Exception as X:
        raise _xError(X, Fmt.INDEX(xs=i), x)


def _2halfeven(s, r, p):
    '''(INTERNAL) Round half-even.
    '''
    if (p > 0 and r > 0) or \
       (p < 0 and r < 0):  # signs match
        r *= 2
        t  = s + r
        if r == (t - s):
            s = t
    return s


def _1primed(xs):  # in .fmath
    '''(INTERNAL) 1-Primed summation of iterable C{xs}
       items, all I{known} to be C{finite float}.
    '''
    yield _1_0
    for x in xs:
        yield x
    yield _N_1_0


def _2ps(s, r):
    '''(INTERNAL) Return a C{s} and C{r} pair, I{ps-ordered}.
    '''
    return (s, r) if fabs(s) < fabs(r) else (r, s)


def _psum(ps):  # PYCHOK used!
    '''(INTERNAL) Partials sum, updating C{ps}, I{overridden below}.
    '''
    # assert isinstance(ps, list)
    i   =  len(ps) - 1
    s   = _0_0 if i < 0 else ps[i]
    _2s = _2sum
    while i > 0:
        i -= 1
        s, r = _2s(s, ps[i])
        if r:  # sum(ps) became inexact
            if s:
                ps[i:] = r, s
                if i > 0:
                    s = _2halfeven(s, r, ps[i-1])
                break  # return s
            s = r  # PYCHOK no cover
        ps[i:] = s,
    return s


def _Psum(ps, **name):
    '''(INTERNAL) Return an C{Fsum} from I{ordered} partials C{ps}.
    '''
    f = Fsum(**name) if name else Fsum()
    if ps:
        f._ps[:] = ps
        f._n = len(f._ps)
    return f


def _Psum_1(p=_1_0, **name):
    '''(INTERNAL) Return an C{Fsum} from a single partial C{p}.
    '''
    f = Fsum(**name) if name else Fsum()
    f._ps[:] = p,
    f._n = 1  # len(f._ps)
    return f


def _2scalar(other, _raiser=None, **mod):
    '''(INTERNAL) Return B{C{other}} as C{int}, C{float} or C{as-is}.
    '''
    if isinstance(other, Fsum):
        s, r = other._fint2
        if r:
            s, r = other._fprs2
            if r:  # PYCHOK no cover
                if _raiser and _raiser(r, s):
                    t = _stresidual(_non_zero_, r, **mod)
                    raise ResidualError(t, txt=None)
                s = other  # L{Fsum} as-is
    else:
        s = other  # C{type} as-is
        if isint(s, both=True):
            s = int(s)
    return s


def _strcomplex(s, *args):
    '''(INTERNAL) C{Complex} 2- or 3-arg C{pow} error as C{str}.
    '''
    c = _strcomplex.__name__[4:]
    n = _DASH_(len(args), _arg_)
    t =  unstr(pow, *args)
    return _SPACE_(c, s, _from_, n, t)


def _stresidual(prefix, residual, **name_values):
    '''(INTERNAL) Residual error as C{str}.
    '''
    p = _stresidual.__name__[3:]
    t =  Fmt.PARENSPACED(p, Fmt(residual))
    for n, v in itemsorted(name_values):
        n =  n.replace(_UNDER_, _SPACE_)
        p =  Fmt.PARENSPACED(n, Fmt(v))
        t = _COMMASPACE_(t, p)
    return _SPACE_(prefix, t)


def _2sum(a, b):  # by .testFmath
    '''(INTERNAL) Return C{a + b} as 2-tuple (sum, residual).
    '''
    s = a + b
    if not _isfinite(s):
        u = unstr(_2sum, a, b)
        t = Fmt.PARENSPACED(_not_finite_, s)
        raise _OverflowError(u, txt=t)
    if fabs(a) < fabs(b):
        a, b = b, a
    return s, (b - (s - a))


class Fsum(_Named):  # sync __methods__ with .vector3dBase.Vector3dBase
    '''Precision floating point summation and I{running} summation.

       Unlike Python's C{math.fsum}, this class accumulates values and provides intermediate,
       I{running}, precision floating point summations.  Accumulation may continue after any
       intermediate, I{running} summuation.

       @note: Accumulated values may be L{Fsum} or C{scalar} instances, any C{type} having
              method C{__float__} to convert the C{scalar} to a single C{float}.

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
#   _ps_max    = 0   # max(Fsum._ps_max, len(Fsum._ps))
    _ratio     = None
    _recursive = bool(_getenv('PYGEODESY_FSUM_RECURSIVE', NN))
    _RESIDUAL  = max(float(_getenv('PYGEODESY_FSUM_RESIDUAL', _0_0)), _0_0)

    def __init__(self, *xs, **name_RESIDUAL):
        '''New L{Fsum} for I{running} precision floating point summation.

           @arg xs: No, one or more initial values, all positional (each C{scalar}
                    or an L{Fsum} instance).
           @kwarg name_RESIDUAL: Optional C{B{name}=NN} for this L{Fsum} and
                       C{B{RESIDUAL}=None} for the L{ResidualError} threshold.

           @see: Methods L{Fsum.fadd} and L{Fsum.RESIDUAL}.
        '''
        if name_RESIDUAL:
            n = _xkwds_get(name_RESIDUAL, name=NN)
            if n:  # set name before ...
                self.name = n
            r = _xkwds_get(name_RESIDUAL, RESIDUAL=None)
            if r is not None:
                self.RESIDUAL(r)  # ... ResidualError
        self._ps = []  # [_0_0], see L{Fsum._fprs}
        if xs:
            self._facc_any(xs, origin=1, up=False)

    def __abs__(self):
        '''Return this instance' absolute value as an L{Fsum}.
        '''
        s = _fsum(self._ps_1())  # == self._cmp_0(0, ...)
        return (-self) if s < 0 else self._copy_2(self.__abs__)

    def __add__(self, other):
        '''Return C{B{self} + B{other}} as an L{Fsum}.

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
        '''Return C{divmod(B{self}, B{other})} as a L{DivMod2Tuple}
           with quotient C{div} an C{int} in Python 3+ or C{float}
           in Python 2- and remainder C{mod} an L{Fsum}.

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
        '''Return this instance' current, precision running sum as C{float}.

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
        return self._fdivmod2(other, _mod_op_ + _fset_op_).mod

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
        return f._fdivmod2(other, _mod_op_).mod

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
        '''Return I{a copy of} this instance, I{negated}.
        '''
        f = self._copy_2(self.__neg__)
        return f._fset(self._neg)

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
        f = self._copy_2r(other, self.__radd__)
        return f._fadd(self, _add_op_)

    def __rdivmod__(self, other):
        '''Return C{divmod(B{other}, B{self})} as 2-tuple C{(quotient,
           remainder)}.

           @see: Method L{Fsum.__divmod__}.
        '''
        f = self._copy_2r(other, self.__rdivmod__)
        return f._fdivmod2(self, _divmod_op_)

#   def __repr__(self):
#       '''Return the default C{repr(this)}.
#       '''
#       return self.toRepr(lenc=True)

    def __rfloordiv__(self, other):
        '''Return C{B{other} // B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__ifloordiv__}.
        '''
        f = self._copy_2r(other, self.__rfloordiv__)
        return f._floordiv(self, _floordiv_op_)

    def __rmatmul__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __rmod__(self, other):
        '''Return C{B{other} % B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__imod__}.
        '''
        f = self._copy_2r(other, self.__rmod__)
        return f._fdivmod2(self, _mod_op_).mod

    def __rmul__(self, other):
        '''Return C{B{other} * B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__imul__}.
        '''
        f = self._copy_2r(other, self.__rmul__)
        return f._fmul(self, _mul_op_)

    def __round__(self, *ndigits):  # PYCHOK no cover
        '''Return C{round(B{self}, *B{ndigits}} as an L{Fsum}.

           @arg ndigits: Optional number of digits (C{int}).
        '''
        # <https://docs.Python.org/3.12/reference/datamodel.html?#object.__round__>
        return _Psum_1(round(float(self), *ndigits),  # can be C{int}
                       name=self.__round__.__name__)

    def __rpow__(self, other, *mod):
        '''Return C{B{other}**B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__ipow__}.
        '''
        f = self._copy_2r(other, self.__rpow__)
        return f._fpow(self, _pow_op_, *mod)

    def __rsub__(self, other):
        '''Return C{B{other} - B{self}} as L{Fsum}.

           @see: Method L{Fsum.__isub__}.
        '''
        f = self._copy_2r(other, self.__rsub__)
        return f._fsub(self, _sub_op_)

    def __rtruediv__(self, other):
        '''Return C{B{other} / B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__itruediv__}.
        '''
        f = self._copy_2r(other, self.__rtruediv__)
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
        if isinstance(other, Fsum):
            s = _fsum(self._ps_1(*other._ps))
        elif isscalar(other):
            if other:
                s = _fsum(self._ps_1(other))
            else:
                s, r = self._fprs2
                s = _signOf(s, -r)
        else:
            raise self._TypeError(op, other)  # txt=_invalid_
        return s

    def copy(self, deep=False, name=NN):
        '''Copy this instance, C{shallow} or B{C{deep}}.

           @return: The copy (L{Fsum}).
         '''
        f = _Named.copy(self, deep=deep, name=name)
        f._ps = list(self._ps)  # separate list
        f._n  = self._n if deep else 1
        return f

    def _copy_2(self, which, name=NN):
        '''(INTERNAL) Copy for I{dyadic} operators.
        '''
        n =  name or which.__name__
        # NOT .classof due to .Fdot(a, *b) args, etc.
        f = _Named.copy(self, deep=False, name=n)
        # assert f._n == self._n
        f._ps = list(self._ps)  # separate list
        return f

    def _copy_2r(self, other, which):
        '''(INTERNAL) Copy for I{reverse-dyadic} operators.
        '''
        return other._copy_2(which) if isinstance(other, Fsum) else \
               Fsum(other, name=which.__name__)

#   def _copy_RESIDUAL(self, other):
#       '''(INTERNAL) Copy C{other._RESIDUAL}.
#       '''
#       R = other._RESIDUAL
#       if R is not Fsum._RESIDUAL:
#           self._RESIDUAL = R

    def divmod(self, other):
        '''Return C{divmod(B{self}, B{other})} as 2-tuple C{(quotient,
           remainder)}.

           @arg other: An L{Fsum} or C{scalar} divisor.

           @return: A L{DivMod2Tuple}C{(div, mod)}, with quotient C{div}
                    an C{int} in Python 3+ or C{float} in Python 2- and
                    remainder C{mod} an L{Fsum} instance.

           @see: Method L{Fsum.__itruediv__}.
        '''
        f = self._copy_2(self.divmod)
        return f._fdivmod2(other, _divmod_op_)

    def _Error(self, op, other, Error, **txt_cause):
        '''(INTERNAL) Format an B{C{Error}} for C{{self} B{op} B{other}}.
        '''
        return Error(_SPACE_(self.toStr(), op, other), **txt_cause)

    def _ErrorX(self, X, op, other, *mod):
        '''(INTERNAL) Format the caught exception C{X}.
        '''
        E, t = _xError2(X)
        if mod:
            t = _COMMASPACE_(Fmt.PARENSPACED(mod=mod[0]), t)
        return self._Error(op, other, E, txt=t, cause=X)

    def _ErrorXs(self, X, xs, **kwds):  # in .fmath
        '''(INTERNAL) Format the caught exception C{X}.
        '''
        E, t = _xError2(X)
        n = unstr(self.named3, *xs[:3], _ELLIPSIS=len(xs) > 3, **kwds)
        return E(n, txt=t, cause=X)

    def _facc(self, xs, **up):
        '''(INTERNAL) Accumulate all C{xs}, known to be scalar.
        '''
        self._ps_acc(self._ps, xs, **up)
        return self

    def _facc_(self, *xs, **up):
        '''(INTERNAL) Accumulate all positional C{xs}, known to be scalar.
        '''
        if xs:
            self._ps_acc(self._ps, xs, **up)
        return self

    def _facc_any(self, xs, up=True, **origin_X_x):
        '''(INTERNAL) Accumulate more C{scalars} or L{Fsum}s.
        '''
        self._ps[:] = self._ps_acc(list(self._ps),
                                  _2floats(xs, **origin_X_x), up=up)  # PYCHOK yield
        return self

    def _facc_any_neg(self, xs, up=True, **origin):
        '''(INTERNAL) Accumulate more C{scalars} or L{Fsum}s, negated.
        '''
        def _neg(x):
            return -x

        self._ps[:] = self._ps_acc(list(self._ps), map(_neg,
                                  _2floats(xs, **origin)), up=up)  # PYCHOK yield
        return self

    def _facc_power(self, power, xs, which):  # in .fmath
        '''(INTERNAL) Add each C{xs} as C{float(x**power)}.
        '''
        p = power
        if isinstance(p, Fsum):
            if p.is_exact:
                return self._facc_power(p._fprs, xs, which)
            _Pow = Fsum._pow_any
        elif isint(p, both=True) and p >= 0:
            _Pow, p = Fsum._pow_int, int(p)
        else:
            _Pow, p = Fsum._pow_scalar, _2float(power=p)

        if p:
            from math import pow as _pow
            op  = which.__name__
            _Fs = Fsum

            def _X(X):
                f = _Pow(X, p, power, op)
                return f._ps if isinstance(f, _Fs) else (f,)

            def _x(x):
                return _pow(float(x), p)

            f = self._facc_any(xs, origin=1, _X=_X, _x=_x)
        else:
            f = self._facc_(float(len(xs)))  # x**0 == 1
        return f

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
            self._facc(xs._ps)  # tuple
        elif isscalar(xs):  # for backward compatibility
            self._facc_(_2float(x=xs))  # PYCHOK no cover
        elif xs:
            self._facc_any(xs)
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
        return self._facc_any(xs, origin=1)

    def _fadd(self, other, op, **up):  # in .fmath.Fhorner
        '''(INTERNAL) Apply C{B{self} += B{other}}.
        '''
        if isinstance(other, Fsum):
            self._facc(other._ps, **up)  # tuple
        elif not isscalar(other):
            raise self._TypeError(op, other)  # txt=_invalid_
        elif other:
            self._facc_(other, **up)
        return self

    fcopy   =   copy        # for backward compatibility
    fdiv    = __itruediv__  # for backward compatibility
    fdivmod = __divmod__    # for backward compatibility

    def _fdivmod2(self, other, op):
        '''(INTERNAL) Apply C{B{self} %= B{other}} and return a L{DivMod2Tuple}.
        '''
        # result mostly follows CPython function U{float_divmod
        # <https://GitHub.com/python/cpython/blob/main/Objects/floatobject.c>},
        # but at least divmod(-3, 2) equals Cpython's result (-2, 1).
        q = self._copy_2(self._fdivmod2)._ftruediv(other, op).floor
        if q:  # == float // other == floor(float / other)
            self -= other * q

        s = signOf(other)  # make signOf(self) == signOf(other)
        if s and self.signOf() == -s:  # PYCHOK no cover
            self += other
            q -= 1
#           t  = self.signOf()
#           if t and t != s:
#               raise self._Error(op, other, _AssertionError, txt=signOf.__name__)
        return DivMod2Tuple(q, self)  # q is C{int} in Python 3+, but C{float} in Python 2-

    def _finite(self, other, op=None):
        '''(INTERNAL) Return B{C{other}} if C{finite}.
        '''
        if _isfinite(other):
            return other
        raise ValueError(_not_finite_) if op is None else \
              self._ValueError(op, other, txt=_not_finite_)

    def fint(self, raiser=True, **name):
        '''Return this instance' current running sum as C{integer}.

           @kwarg raiser: If C{True} throw a L{ResidualError} if the
                          I{integer} residual is non-zero (C{bool}).
           @kwarg name: Optional name (C{str}), overriding C{"fint"}.

           @return: The C{integer} (L{Fsum}).

           @raise ResidualError: Non-zero I{integer} residual.

           @see: Methods L{Fsum.int_float} and L{Fsum.is_integer}.
        '''
        i, r = self._fint2
        if r and raiser:
            t = _stresidual(_integer_, r)
            raise ResidualError(_integer_, i, txt=t)
        f = self._copy_2(self.fint, **name)
        return f._fset(i)

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
        s, r = self._fprs2
        i =  int(s)
        r = _fsum(self._ps_1(i)) if r else float(s - i)
        return i, (r or INT0)  # Fsum2Tuple?

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
        while (f - s) > r:  # f > (s + r)
            f -= 1
        return f

#   floordiv = __floordiv__  # for naming consistency

    def _floordiv(self, other, op):  # rather _ffloordiv?
        '''Apply C{B{self} //= B{other}}.
        '''
        q = self._ftruediv(other, op)  # == self
        return self._fset(q.floor)  # floor(q)

    fmul = __imul__  # for backward compatibility

    def _fmul(self, other, op):
        '''(INTERNAL) Apply C{B{self} *= B{other}}.
        '''
        if isinstance(other, Fsum):
            if len(self._ps) != 1:
                f = self._mul_Fsum(other, op)
            elif len(other._ps) != 1:  # and len(self._ps) == 1
                f = other._mul_scalar(self._ps[0], op)
            else:  # len(other._ps) == len(self._ps) == 1
                f = self._finite(self._ps[0] * other._ps[0])
        elif isscalar(other):
            f = self._mul_scalar(other, op) if other != _1_0 else self
        else:
            raise self._TypeError(op, other)  # txt=_invalid_
        return self._fset(f)  # n=len(self) + 1

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
        if mod:
            if mod[0] is not None:  # == 3-arg C{pow}
                f = self._pow_2_3(self, other, other, op, *mod)
            elif self.is_integer():
                # return an exact C{int} for C{int}**C{int}
                i, _ = self._fint2  # assert _ == 0
                x = _2scalar(other)  # C{int}, C{float} or other
                f =  self._pow_2_3(i, x, other, op) if isscalar(x) else \
                    _Psum_1(i)._pow_any(x, other, op)
            else:  # mod[0] is None, power(self, other)
                f =  self._pow_any(other, other, op)
        else:  # pow(self, other) == pow(self, other, None)
            f = self._pow_any(other, other, op)
        return self._fset(f, asis=isint(f))  # n=max(len(self), 1)

    @Property_RO
    def _fprs(self):
        '''(INTERNAL) Get and cache this instance' precision
           running sum (C{float} or C{int}), ignoring C{residual}.

           @note: The precision running C{fsum} after a C{//=} or
                  C{//} C{floor} division is C{int} in Python 3+.
        '''
        return self._fprs2.fsum

    @Property_RO
    def _fprs2(self):
        '''(INTERNAL) Get and cache this instance' precision
           running sum and residual (L{Fsum2Tuple}).
        '''
        ps = self._ps
        n  = len(ps) - 2
        if n > 0:  # len(ps) > 2
            s = _psum(ps)
            n =  len(ps) - 2
            if n > 0:
                r = _fsum(self._ps_1(s)) or INT0
                return Fsum2Tuple(s, r)
        if n == 0:  # len(ps) == 2
            ps[:] = _2ps(*_2sum(*ps))
            r, s  = (INT0, ps[0]) if len(ps) != 2 else ps
        elif ps:  # len(ps) == 1
            s, r = ps[0], INT0
        else:  # len(ps) == 0
            s, r = _0_0, INT0
            ps[:] = s,
        # assert self._ps is ps
        return Fsum2Tuple(s, r)

#   def _fpsqz(self):
#       '''(INTERNAL) Compress, squeeze the C{partials}.
#       '''
#       if len(self._ps) > 2:
#           _ = self._fprs2
#       return self

    def fset_(self, *xs):
        '''Replace this instance' value with C{xs}.

           @arg xs: Optional, new values (C{scalar} or L{Fsum}
                    instances), all positional.

           @return: This instance (C{Fsum}).

           @see: Method L{Fsum.fadd} for further details.
        '''
        self._ps[:] = 0,
        self._n = 0
        return self.fadd(xs) if xs else self._update()

    def _fset(self, other, asis=True, n=0):
        '''(INTERNAL) Overwrite this instance with an other or a C{scalar}.
        '''
        if other is self:
            pass  # from ._fmul, ._ftruediv and ._pow_scalar
        elif isinstance(other, Fsum):
            self._ps[:] = other._ps
            self._n     = n or other._n
#           self._copy_RESIDUAL(other)
            # use or zap the C{Property_RO} values
            Fsum._fint2._update_from(self, other)
            Fsum._fprs ._update_from(self, other)
            Fsum._fprs2._update_from(self, other)
        elif isscalar(other):
            s = other if asis else float(other)
            i = int(s)  # see ._fint2
            t = i, ((s - i) or INT0)
            self._ps[:] = s,
            self._n     = n or 1
            # Property_ROs _fint2, _fprs and _fprs2 can't be a Property:
            # Property's _fset zaps the value just set by the @setter
            self.__dict__.update(_fint2=t, _fprs=s, _fprs2=Fsum2Tuple(s, INT0))
        else:  # PYCHOK no cover
            raise self._TypeError(_fset_op_, other)  # txt=_invalid_
        return self

    def _fset_ps(self, other, n=0):  # in .fmath
        '''(INTERNAL) Set partials from a known C{Fsum} or C{scalar}.
        '''
        if isinstance(other, Fsum):
            self._ps[:] = other._ps
            self._n     = n or other._n
        else:  # assert isscalar(other)
            self._ps[:] = other,
            self._n     = n or 1
        return self

    def fsub(self, xs=()):
        '''Subtract an iterable of C{scalar} or L{Fsum} instances from
           this instance.

           @arg xs: Iterable, list, tuple. etc. (C{scalar} or L{Fsum}
                    instances).

           @return: This instance, updated (L{Fsum}).

           @see: Method L{Fsum.fadd}.
        '''
        return self._facc_any_neg(xs) if xs else self

    def fsub_(self, *xs):
        '''Subtract all positional C{scalar} or L{Fsum} instances from
           this instance.

           @arg xs: Values to subtract (C{scalar} or L{Fsum} instances),
                    all positional.

           @return: This instance, updated (L{Fsum}).

           @see: Method L{Fsum.fadd}.
        '''
        return self._facc_any_neg(xs, origin=1) if xs else self

    def _fsub(self, other, op):
        '''(INTERNAL) Apply C{B{self} -= B{other}}.
        '''
        if isinstance(other, Fsum):
            if other is self:  # or other._fprs2 == self._fprs2:
                self._fset(_0_0)  # n=len(self) * 2, self -= self
            elif other._ps:
                self._facc(other._ps_neg)
        elif not isscalar(other):
            raise self._TypeError(op, other)  # txt=_invalid_
        elif self._finite(other, op):
            self._facc_(-other)
        return self

    def fsum(self, xs=()):
        '''Add more C{scalar} or L{Fsum} instances and summate.

           @kwarg xs: Iterable, list, tuple, etc. (C{scalar} or
                      L{Fsum} instances).

           @return: Precision running sum (C{float} or C{int}).

           @see: Method L{Fsum.fadd}.

           @note: Accumulation can continue after summation.
        '''
        f = self._facc_any(xs) if xs else self
        return f._fprs

    def fsum_(self, *xs):
        '''Add all positional C{scalar} or L{Fsum} instances and summate.

           @arg xs: Values to add (C{scalar} or L{Fsum} instances), all
                    positional.

           @return: Precision running sum (C{float} or C{int}).

           @see: Methods L{Fsum.fsum}, L{Fsum.Fsum_} and L{Fsum.fsumf_}.
        '''
        f = self._facc_any(xs, origin=1) if xs else self
        return f._fprs

    def Fsum_(self, *xs):
        '''Like method L{Fsum.fsum_} but returning an L{Fsum}.

           @return: Current, precision running sum (L{Fsum}).
        '''
        return self._facc_any(xs, origin=1)._copy_2(self.Fsum_)

    def fsum2(self, xs=(), name=NN):
        '''Add more C{scalar} or L{Fsum} instances and return the
           current precision running sum and the C{residual}.

           @kwarg xs: Iterable, list, tuple, etc. (C{scalar} or L{Fsum}
                      instances).
           @kwarg name: Optional name (C{str}).

           @return: L{Fsum2Tuple}C{(fsum, residual)} with C{fsum} the
                    current precision running sum and C{residual}, the
                    (precision) sum of the remaining C{partials}.  The
                    C{residual is INT0} if the C{fsum} is considered
                    to be I{exact}.

           @see: Methods L{Fsum.fint2}, L{Fsum.fsum} and L{Fsum.fsum2_}
        '''
        f = self._facc_any(xs) if xs else self
        t = f._fprs2
        if name:
            t = t.dup(name=name)
        return t

    def fsum2_(self, *xs):
        '''Add any positional C{scalar} or L{Fsum} instances and return
           the precision running sum and the C{differential}.

           @arg xs: Values to add (C{scalar} or L{Fsum} instances), all
                    positional.

           @return: 2Tuple C{(fsum, delta)} with the current, precision
                    running C{fsum} like method L{Fsum.fsum} and C{delta},
                    the difference with previous running C{fsum}, C{float}.

           @see: Methods L{Fsum.fsum_} and L{Fsum.fsum}.
        '''
        return self._fsum2f_any(xs, self._facc_any, origin=1)

    def fsumf_(self, *xs):
        '''Like method L{Fsum.fsum_} but only for I{known} C{float B{xs}}.
        '''
        f = self._facc(xs) if xs else self
        return f._fprs

    def Fsumf_(self, *xs):
        '''Like method L{Fsum.Fsum_} but only for I{known} C{float B{xs}}.
        '''
        return self._facc(xs)._copy_2(self.Fsumf_)

    def fsum2f_(self, *xs):
        '''Like method L{Fsum.fsum2_} but only for I{known} C{float B{xs}}.
        '''
        return self._fsum2f_any(xs, self._facc)

    def _fsum2f_any(self, xs, _facc, **origin):
        '''(INTERNAL) Helper for L{Fsum.fsum2_} and L{Fsum.fsum2f_}.
        '''
        p, q = self._fprs2
        if xs:
            s, r = _facc(xs, **origin)._fprs2
            return s, _2delta(s - p, r - q)  # _fsum(_1primed((s, -p, r, -q))
        else:
            return p, _0_0

#   ftruediv = __itruediv__   # for naming consistency?

    def _ftruediv(self, other, op):
        '''(INTERNAL) Apply C{B{self} /= B{other}}.
        '''
        n = _1_0
        if isinstance(other, Fsum):
            if other is self or other == self:
                return self._fset(_1_0)  # n=len(self)
            d, r = other._fprs2
            if r:
                if d:
                    if self._raiser(r, d):
                        raise self._ResidualError(op, other, r)
                    d, n = other.as_integer_ratio()
                else:  # PYCHOK no cover
                    d = r
        elif isscalar(other):
            d = other
        else:  # PYCHOK no cover
            raise self._TypeError(op, other)  # txt=_invalid_
        try:
            s = 0 if isinf(d) else (
                d if isnan(d) else self._finite(n / d))
        except Exception as X:
            raise self._ErrorX(X, op, other)
        f = self._mul_scalar(s, _mul_op_)  # handles 0, NAN, etc.
        return self._fset(f, asis=False)

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
        '''Is this instance' running C{fsum} considered to be exact? (C{bool}).
        '''
        return self.residual is INT0

    def is_integer(self):
        '''Is this instance' running sum C{integer}? (C{bool}).

           @see: Methods L{Fsum.fint} and L{Fsum.fint2}.
        '''
        _, r = self._fint2
        return False if r else True

    def is_math_fsum(self):
        '''Return whether functions L{fsum}, L{fsum_}, L{fsum1} and
           L{fsum1_} plus partials summation are based on Python's
           C{math.fsum} or not.

           @return: C{2} if all functions and partials summation
                    are based on C{math.fsum}, C{True} if only
                    the functions are based on C{math.fsum} (and
                    partials summation is not) or C{False} if
                    none are.
        '''
        f = Fsum._math_fsum
        return 2 if _psum is f else bool(f)

    def _mul_Fsum(self, other, op=_mul_op_):  # in .fmath.Fhorner
        '''(INTERNAL) Return C{B{self} * Fsum B{other}} as L{Fsum} or C{0}.
        '''
        # assert isinstance(other, Fsum)
        if self._ps and other._ps:
            f =  self._ps_mul(op, *other._ps)  # NO ._2scalar
        else:
            f = _0_0
        return f

    def _mul_scalar(self, factor, op):  # in .fmath.Fhorner
        '''(INTERNAL) Return C{B{self} * scalar B{factor}} as L{Fsum}, C{0} or C{self}.
        '''
        # assert isscalar(factor)
        if self._ps and self._finite(factor, op):
            f =  self      if factor == _1_0 else (
                 self._neg if factor == _N_1_0 else
                 self._ps_mul(op, factor)._2scalar)
        else:
            f = _0_0
        return f

    @property_RO
    def _neg(self):
        '''(INTERNAL) Return C{Fsum(-self)} or scalar C{NEG0}.
        '''
        return _Psum(self._ps_neg) if self._ps else NEG0

    @property_RO
    def partials(self):
        '''Get this instance' current, partial sums (C{tuple} of C{float}s).
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
                  or set to C{integer} by an L{Fsum.fint} call.

           @see: Methods L{Fsum.__ipow__}, L{Fsum.fint} and L{Fsum.is_integer}.
        '''
        f = self._copy_2(self.pow)
        return f._fpow(x, _pow_op_, *mod)  # f = pow(f, x, *mod)

    def _pow_0_1(self, x, other):
        '''(INTERNAL) Return B{C{self}**1} or C{B{self}**0 == 1.0}.
        '''
        return self if x else (1 if isint(other) and self.is_integer() else _1_0)

    def _pow_2_3(self, b, x, other, op, *mod):
        '''(INTERNAL) 2-arg C{pow(B{b}, scalar B{x})} and 3-arg C{pow(B{b},
           B{x}, int B{mod} or C{None})}, embellishing errors.
        '''
        try:
            if mod:  # b, x, mod all C{int}, unless C{mod} is C{None}
                m = mod[0]
                b, r = b._fprs2 if m is None else b._fint2
                if r and self._raiser(r, b):
                    t = _non_zero_ if m is None else _integer_
                    raise ResidualError(_stresidual(t, r, mod=m), txt=None)
                x = _2scalar(x, _raiser=self._raiser, mod=m)
            # 0**INF == 0.0, 1**INF == 1.0, -1**2.3 == -(1**2.3)
            s = pow(b, x, *mod)
            if iscomplex(s):
                # neg**frac == complex in Python 3+, but ValueError in 2-
                raise ValueError(_strcomplex(s, b, x, *mod))
            return self._finite(s)
        except Exception as X:
            raise self._ErrorX(X, op, other, *mod)

    def _pow_any(self, other, unused, op):
        '''Return C{B{self} ** B{other}}.
        '''
        if isinstance(other, Fsum):
            x, r = other._fprs2
            if r and self._raiser(r, x):
                raise self._ResidualError(op, other, r)
            f = self._pow_scalar(x, other, op)
            if r:
                f *= self._pow_scalar(r, other, op)
        elif isscalar(other):
            x = self._finite(other, op)
            f = self._pow_scalar(x, other, op)
        else:
            raise self._TypeError(op, other)  # txt=_invalid_
        return f

    def _pow_int(self, x, other, op):
        '''(INTERNAL) Return C{B{self} **= B{x}} for C{int B{x} >= 0}.
        '''
        # assert isint(x) and x >= 0
        ps = self._ps
        if len(ps) > 1:
            _mul_Fsum = Fsum._mul_Fsum
            if x > 4:
                p = self
                f = self if (x & 1) else _Psum_1()
                m = x >> 1  # // 2
                while m:
                    p = _mul_Fsum(p, p, op)  # p **= 2
                    if (m & 1):
                        f = _mul_Fsum(f, p, op)  # f *= p
                    m >>= 1  # //= 2
            elif x > 1:  # self**2, 3 or 4
                f = _mul_Fsum(self, self, op)
                if x > 2:  # self**3 or 4
                    p =  self if x < 4 else f
                    f = _mul_Fsum(f, p, op)._2scalar
            else:  # self**1 or self**0 == 1 or _1_0
                f = self._pow_0_1(x, other)
        elif ps:  # self._ps[0]**x
            f = self._pow_2_3(ps[0], x, other, op)
        else:  # PYCHOK no cover
            # 0**pos_int == 0, but 0**0 == 1
            f = 0 if x else 1
        return f

    def _pow_scalar(self, x, other, op):
        '''(INTERNAL) Return C{self**B{x}} for C{scalar B{x}}.
        '''
        s, r = self._fprs2
        if isint(x, both=True):
            x = int(x)  # Fsum**int
            y = abs(x)
            if y > 1:
                if r:
                    f = self._pow_int(y, other, op)
                    if x > 0:  # > 1
                        return f
                    # assert x < 0  # < -1
                    s, r = f._fprs2 if isinstance(f, Fsum) else (f, 0)
                    if r:
                        return _Psum_1()._ftruediv(f, op)
                    # use **= -1 for the CPython float_pow
                    # error if s is zero, and not s = 1 / s
                    x = -1
            elif x < 0:  # == -1: self**(-1) == 1 / self
                if r:
                    return _Psum_1()._ftruediv(self, op)
            else:  # self**1 or self**0
                return self._pow_0_1(x, other)  # self, 1 or 1.0
        elif not isscalar(x):  # assert ...
            raise self._TypeError(op, other, txt=_not_scalar_)
        elif r and self._raiser(r, s):  # non-zero residual**fractional
            # raise self._ResidualError(op, other, r, fractional_power=x)
            t = _stresidual(_non_zero_, r, fractional_power=x)
            raise self._Error(op, other, ResidualError, txt=t)
        # assert isscalar(s) and isscalar(x)
        return self._pow_2_3(s, x, other, op)

    def _ps_1(self, *less):
        '''(INTERNAL) Yield partials, 1-primed and subtract any C{less}.
        '''
        yield _1_0
        for p in self._ps:
            yield  p
        for p in less:
            yield -p
        yield _N_1_0

    def _ps_acc(self, ps, xs, up=True):
        '''(INTERNAL) Accumulate all scalar C{xs} into C{ps}.
        '''
        n   =  0
        _2s = _2sum
        for x in (tuple(xs) if xs is ps else xs):
            # assert isscalar(x) and _isfinite(x)
            if x:
                i = 0
                for p in ps:
                    x, p = _2s(x, p)
                    if p:
                        ps[i] = p
                        i += 1
                ps[i:] = (x,) if x else ()
                n += 1
        if n:
            self._n += n
            # Fsum._ps_max = max(Fsum._ps_max, len(ps))
            if up:
                self._update()
        return ps

    def _ps_mul(self, op, *factors):
        '''(INTERNAL) Multiply this instance' C{partials} with
           each of the B{C{factors}}, all known to be scalar.
        '''
        def _pfs(ps, fs):
            if len(ps) < len(fs):
                ps, fs = fs, ps
            _fin = _isfinite
            for f in fs:
                for p in ps:
                    p *= f
                    yield p if _fin(p) else self._finite(p, op)

        return _Psum(self._ps_acc([], _pfs(self._ps, factors)))

    @property_RO
    def _ps_neg(self):
        '''(INTERNAL) Yield the partials, I{negated}.
        '''
        for p in self._ps:
            yield -p

    def _raiser(self, r, s):
        '''(INTERNAL) Does ratio C{r / s} exceed threshold?
        '''
        self._ratio = t = fabs((r / s) if s else r)
        return t > self._RESIDUAL

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

    def RESIDUAL(self, *threshold):
        '''Get and set this instance' I{ratio} for raising L{ResidualError}s,
           overriding the default from env variable C{PYGEODESY_FSUM_RESIDUAL}.

           @arg threshold: If C{scalar}, the I{ratio} to exceed for raising
                           L{ResidualError}s in division and exponention, if
                           C{None} restore the default set with env variable
                           C{PYGEODESY_FSUM_RESIDUAL} or if omitted, keep the
                           current setting.

           @return: The previous C{RESIDUAL} setting (C{float}), default C{0}.

           @raise ValueError: Negative B{C{threshold}}.

           @note: L{ResidualError}s will be thrown if the non-zero I{ratio}
                  C{residual / fsum} exceeds the B{C{threshold}}.
        '''
        r = self._RESIDUAL
        if threshold:
            t = threshold[0]
            t = Fsum._RESIDUAL if t is None else (
                      float(t) if isscalar(t) else (  # for backward ...
                          _0_0 if bool(t) else _1_0))  # ... compatibility
            if t < 0:
                u = _DOT_(self, unstr(self.RESIDUAL, *threshold))
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

    @property_RO
    def _2scalar(self):
        '''(INTERNAL) Get this instance as C{scalar} or C{as-is}.
        '''
        s, r = self._fprs2
        return self if r else s

    def signOf(self, res=True):
        '''Determine the sign of this instance.

           @kwarg res: If C{True} consider, otherwise
                       ignore the residual (C{bool}).

           @return: The sign (C{int}, -1, 0 or +1).
        '''
        s, r = self._fprs2
        return _signOf(s, (-r) if res else 0)

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

    def _update(self, updated=True):  # see ._fset
        '''(INTERNAL) Zap all cached C{Property_RO} values.
        '''
        if updated:
            _pop = self.__dict__.pop
            for p in _ROs:
                _ = _pop(p, None)
#           Fsum._fint2._update(self)
#           Fsum._fprs ._update(self)
#           Fsum._fprs2._update(self)
        return self  # for .fset_

    def _ValueError(self, op, other, **txt):  # PYCHOK no cover
        '''(INTERNAL) Return a C{ValueError}.
        '''
        return self._Error(op, other, _ValueError, **txt)

    def _ZeroDivisionError(self, op, other, **txt):  # PYCHOK no cover
        '''(INTERNAL) Return a C{ZeroDivisionError}.
        '''
        return self._Error(op, other, _ZeroDivisionError, **txt)

_ROs = _allPropertiesOf_n(3, Fsum, Property_RO)  # PYCHOK assert, see Fsum._fset, -._update


def _Float_Int(arg, **name_Error):
    '''(INTERNAL) Unit of L{Fsum2Tuple} items.
    '''
    U = Int if isint(arg) else Float
    return U(arg, **name_Error)


class DivMod2Tuple(_NamedTuple):
    '''2-Tuple C{(div, mod)} with the quotient C{div} and remainder
       C{mod} results of a C{divmod} operation.

       @note: Quotient C{div} an C{int} in Python 3+ or a C{float} in
              Python 2-.  Remainder C{mod} an L{Fsum} instance.
    '''
    _Names_ = (_div_,     _mod_)
    _Units_ = (_Float_Int, Fsum)


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
    def _Fsum(self):
        '''(INTERNAL) Get this L{Fsum2Tuple} as an L{Fsum}.
        '''
        s, r = map(float, self)
        return _Psum(_2ps(s, r), name=self.name)

    def is_exact(self):
        '''Is this L{Fsum2Tuple} considered to be exact? (C{bool}).
        '''
        return self._Fsum.is_exact()

    def is_integer(self):
        '''Is this L{Fsum2Tuple} C{integer}? (C{bool}).
        '''
        return self._Fsum.is_integer()


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

    if _getenv('PYGEODESY_FSUM_PARTIALS', NN) == _fsum.__name__:
        _psum = _fsum  # PYCHOK re-def

except ImportError:
    _sum = sum  # Fsum(NAN) exception fall-back, in .elliptic

    def _fsum(xs):
        '''(INTERNAL) Precision summation, Python 2.5-.
        '''
        f = Fsum()
        f.name = _fsum.__name__
        return f.fsum(xs)


def fsum(xs, floats=False):
    '''Precision floating point summation based on or like Python's C{math.fsum}.

       @arg xs: Iterable, list, tuple, etc. of values (C{scalar} or L{Fsum}
                instances).
       @kwarg floats: Use C{B{floats}=True} iff I{all} B{C{xs}} are known
                      to be C{float} scalars (C{bool}).

       @return: Precision C{fsum} (C{float}).

       @raise OverflowError: Partial C{2sum} overflow.

       @raise TypeError: Non-scalar B{C{xs}} value.

       @raise ValueError: Invalid or non-finite B{C{xs}} value.

       @note: Exception and I{non-finite} handling may differ if not based
              on Python's C{math.fsum}.

       @see: Class L{Fsum} and methods L{Fsum.fsum} and L{Fsum.fadd}.
    '''
    return _fsum(xs if floats else _2floats(xs)) if xs else _0_0  # PYCHOK yield


def fsum_(*xs, **floats):
    '''Precision floating point summation of all positional arguments.

       @arg xs: Values to be added (C{scalar} or L{Fsum} instances), all
                positional.
       @kwarg floats: Use C{B{floats}=True} iff I{all} B{C{xs}} are known
                      to be C{float} scalars (C{bool}).

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
    '''Precision floating point summation, 1-primed.

       @arg xs: Iterable, list, tuple, etc. of values (C{scalar} or L{Fsum}
                instances).
       @kwarg floats: Use C{B{floats}=True} iff I{all} B{C{xs}} are known
                      to be C{float}.

       @return: Precision C{fsum} (C{float}).

       @see: Function C{fsum}.
    '''
    return _fsum(_1primed(xs if floats else _2floats(xs))) if xs else _0_0  # PYCHOK yield


def fsum1_(*xs, **floats):
    '''Precision floating point summation, 1-primed.

       @arg xs: Values to be added (C{scalar} or L{Fsum} instances), all
                positional.
       @kwarg floats: Use C{B{floats}=True} iff I{all} B{C{xs}} are known
                      to be C{float} scalars (C{bool}).

       @return: Precision C{fsum} (C{float}).

       @see: Function C{fsum}
    '''
    return _fsum(_1primed(xs if _xkwds_get(floats, floats=False) else
                         _2floats(xs, origin=1))) if xs else _0_0  # PYCHOK yield


def fsum1f_(*xs):
    '''Precision floating point summation, L{fsum1_}C{(*xs, floats=True)}.
    '''
    return _fsum(_1primed(xs)) if xs else _0_0


if __name__ == '__main__':

    # usage: [env PYGEODESY_FSUM_PARTIALS=fsum] python3 -m pygeodesy.fsums

    def _test(n):
        # copied from Hettinger, see L{Fsum} reference
        from pygeodesy import printf
        from random import gauss, random, shuffle

        printf(_fsum.__name__, end=_COMMASPACE_)
        printf(_psum.__name__, end=_COMMASPACE_)

        F = Fsum()
        if F.is_math_fsum():
            c = (7, 1e100, -7, -1e100, -9e-20, 8e-20) * 10
            for _ in range(n):
                t = list(c)
                s = 0
                for _ in range(n * 8):
                    v = gauss(0, random())**7 - s
                    t.append(v)
                    s += v
                shuffle(t)
                assert float(F.fset_(*t)) == _fsum(t)
                printf(_DOT_, end=NN)
        printf(NN)

    _test(128)

# **) MIT License
#
# Copyright (C) 2016-2024 -- mrJean1 at Gmail -- All Rights Reserved.
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
