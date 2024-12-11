
# -*- coding: utf-8 -*-

u'''Class L{Fsum} for precision floating point summation similar to
Python's C{math.fsum} enhanced with I{running} summation and as an
option, accurate I{TwoProduct} multiplication.

Accurate multiplication is based on the C{math.fma} function for
Python 3.13 and newer or one of two equivalent C{fma} implementations
for Python 3.12 and older.  To enable accurate multiplication, set
env variable C{PYGEODESY_FSUM_F2PRODUCT} to C{"std"} or any non-empty
string or invoke function C{pygeodesy.f2product(True)} or set.  With
C{"std"} the C{fma} implemention follows the C{math.fma} function,
otherwise the C{PyGeodesy 24.09.09} release.

Generally, an L{Fsum} instance is considered a C{float} plus a small or
zero C{residue} aka C{residual} value, see property L{Fsum.residual}.

Set env variable C{PYGEODESY_FSUM_RESIDUAL} to a C{float} string greater
than C{"0.0"} as the threshold to throw a L{ResidualError} for a division,
power or root operation of an L{Fsum} with a C{residual} I{ratio} exceeding
the threshold.  See methods L{Fsum.RESIDUAL}, L{Fsum.pow}, L{Fsum.__ipow__}
and L{Fsum.__itruediv__}.

There are several C{integer} L{Fsum} cases, for example the result from
functions C{ceil}, C{floor}, C{Fsum.__floordiv__} and methods L{Fsum.fint},
L{Fsum.fint2} and L{Fsum.is_integer}.  Also, L{Fsum} methods L{Fsum.pow},
L{Fsum.__ipow__}, L{Fsum.__pow__} and L{Fsum.__rpow__} return a (very long)
C{int} if invoked with optional argument C{mod} set to C{None}.  The
C{residual} of an C{integer} L{Fsum} is between C{-1.0} and C{+1.0} and
will be C{INT0} if that is considered to be I{exact}.

Set env variable C{PYGEODESY_FSUM_NONFINITES} to C{"std"} or use function
C{pygeodesy.nonfiniterrors(False)} to allow I{non-finite} C{float}s like
C{inf}, C{INF}, C{NINF}, C{nan} and C{NAN} and to ignore C{OverflowError}
respectively C{ValueError} exceptions.  However, in that case I{non-finite}
results may differ from Python's C{math.fsum} results.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import isbool, iscomplex, isint, isscalar, \
                            _signOf, itemsorted, signOf, _xiterable
from pygeodesy.constants import INF, INT0, MANT_DIG, NEG0, NINF, _0_0, \
                               _1_0, _N_1_0, _isfinite, _pos_self, \
                               Float, Int
from pygeodesy.errors import _AssertionError, _OverflowError, _TypeError, \
                             _ValueError, _xError, _xError2, _xkwds, \
                             _xkwds_get, _xkwds_get1, _xkwds_not, \
                             _xkwds_pop, _xsError
from pygeodesy.internals import _enquote, _getPYGEODESY, _MODS, _passarg
from pygeodesy.interns import NN, _arg_, _COMMASPACE_, _DOT_, _from_, \
                             _not_finite_, _SPACE_, _std_, _UNDER_
# from pygeodesy.lazily import _ALL_LAZY  # from .named
from pygeodesy.named import _name__, _name2__, _Named, _NamedTuple, \
                            _NotImplemented,  _ALL_LAZY
from pygeodesy.props import _allPropertiesOf_n, deprecated_method, \
                             deprecated_property_RO, Property, \
                             Property_RO, property_RO
from pygeodesy.streprs import Fmt, fstr, unstr
# from pygeodesy.units import Float, Int  # from .constants

from math import fabs, isinf, isnan, \
                 ceil as _ceil, floor as _floor  # PYCHOK used! .ltp

__all__ = _ALL_LAZY.fsums
__version__ = '24.12.02'

from pygeodesy.interns import (
  _PLUS_     as _add_op_,  # in .auxilats.auxAngle
  _EQUAL_    as _fset_op_,
  _RANGLE_   as _gt_op_,
  _LANGLE_   as _lt_op_,
  _PERCENT_  as _mod_op_,
  _STAR_     as _mul_op_,
  _NOTEQUAL_ as _ne_op_,
  _DASH_     as _sub_op_,  # in .auxilats.auxAngle
  _SLASH_    as _truediv_op_
)
_floordiv_op_ = _truediv_op_ * 2  # _DSLASH_
_divmod_op_   = _floordiv_op_ + _mod_op_
_F2PRODUCT    = _getPYGEODESY('FSUM_F2PRODUCT')
_iadd_op_     = _add_op_ + _fset_op_  # in .auxilats.auxAngle, .fstats
_integer_     = 'integer'
_isub_op_     = _sub_op_ + _fset_op_  # in .auxilats.auxAngle
_NONFINITEr   = _0_0  # NOT INT0!
_NONFINITES   = _getPYGEODESY('FSUM_NONFINITES')
_non_zero_    = 'non-zero'
_pow_op_      = _mul_op_ * 2  # _DSTAR_
_RESIDUAL_0_0 = _getPYGEODESY('FSUM_RESIDUAL', _0_0)
_significant_ = 'significant'
_threshold_   = 'threshold'


def _2finite(x, _isfine=_isfinite):  # in .fstats
    '''(INTERNAL) return C{float(x)} if finite.
    '''
    return (float(x) if _isfine(x)  # and isscalar(x)
            else _nfError(x))


def _2float(index=None, _isfine=_isfinite, **name_x):  # in .fmath, .fstats
    '''(INTERNAL) Raise C{TypeError} or C{Overflow-/ValueError} if C{x} not finite.
    '''
    n, x = name_x.popitem()  # _xkwds_item2(name_x)
    try:
        f = float(x)
        return f if _isfine(f) else _nfError(x)
    except Exception as X:
        raise _xError(X, Fmt.INDEX(n, index), x)


try:  # MCCABE 26
    from math import fma as _fma

    def _2products(x, ys, *zs):
        # yield(x * y for y in ys) + yield(z in zs)
        # TwoProductFMA U{Algorithm 3.5
        # <https://www.TUHH.De/ti3/paper/rump/OgRuOi05.pdf>}
        for y in ys:
            f = x * y
            yield f
            if _isfinite(f):
                yield _fma(x, y, -f)
        for z in zs:
            yield z

#   _2split3  = \
    _2split3s = _passarg  # in Fsum.is_math_fma

except  ImportError:  # PYCHOK DSPACE! Python 3.12-

    if _F2PRODUCT and _F2PRODUCT != _std_:
        # backward to PyGeodesy 24.09.09, with _fmaX

        def _fma(*a_b_c):  # PYCHOK no cover
            # mimick C{math.fma} from Python 3.13+,
            # the same accuracy, but ~14x slower
            (na, da), (nb, db), (nc, dc) = map(_2n_d, a_b_c)
            n  = na * nb * dc
            n += da * db * nc
            d  = da * db * dc
            try:
                n, d = _n_d2(n, d)
                r    = float(n / d)
            except OverflowError:  # "integer division result too large ..."
                r = NINF if (_signOf(n, 0) * _signOf(d, 0)) < 0 else INF
            return r if _isfinite(r) else _fmaX(r, *a_b_c)  # "overflow in fma"

        def _2n_d(x):  # PYCHOK no cover
            try:  # int.as_integer_ratio in 3.8+
                return  x.as_integer_ratio()
            except (AttributeError, OverflowError, TypeError, ValueError):
                return (x if isint(x) else float(x)), 1
    else:

        def _fma(a, b, c):  # PYCHOK redef
            # mimick C{math.fma} from Python 3.13+,
            # the same accuracy, but ~13x slower
            b3s = _2split3(b),  # 1-tuple of 3-tuple
            r   = _fsum(_2products(a, b3s, c))
            return r if _isfinite(r) else _fmaX(r, a, b, c)

        _2n_d = None  # redef

    def _fmaX(r, *a_b_c):  # PYCHOK no cover
        # handle non-finite as Python 3.13+ C-function U{math_fma_impl<https://
        # GitHub.com/python/cpython/blob/main/Modules/mathmodule.c#L2305>}:
        #  raise a ValueError for a NAN result from non-NAN C{a_b_c}s or an
        #  OverflowError for a non-NAN non-finite from all finite C{a_b_c}s.
        if isnan(r):
            def _x(x):
                return not isnan(x)
        else:  # non-NAN non-finite
            _x = _isfinite
        if all(map(_x, a_b_c)):
            raise _nfError(r, unstr(_fma, *a_b_c))
        return r

    def _2products(x, y3s, *zs):  # PYCHOK in _fma, ...
        # yield(x * y3 for y3 in y3s) + yield(z in zs)
        # TwoProduct U{Algorithm 3.3
        # <https://www.TUHH.De/ti3/paper/rump/OgRuOi05.pdf>}
        # also in Python 3.13+ C{Modules/mathmodule.c} under
        # #ifndef UNRELIABLE_FMA ... #else ... #endif
        _, a, b = _2split3(x)
        for y, c, d in y3s:
            y *= x
            yield y
            if _isfinite(y):
                # yield b * d - (((y - a * c) - b * c) - a * d)
                #     = b * d + (a * d - ((y - a * c) - b * c))
                #     = b * d + (a * d + (b * c - (y - a * c)))
                #     = b * d + (a * d + (b * c + (a * c - y)))
                yield a * c - y
                yield b * c
                if d:
                    yield a * d
                    yield b * d
        for z in zs:
            yield z

    _2FACTOR = pow(2, (MANT_DIG + 1) // 2) + _1_0  # 134217729 if MANT_DIG == 53

    def _2split3(x):
        # Split U{Algorithm 3.2
        # <https://www.TUHH.De/ti3/paper/rump/OgRuOi05.pdf>}
        a  = c = x * _2FACTOR
        a -= c - x
        b  = x - a
        return x, a, b

    def _2split3s(xs):  # in Fsum.is_math_fma
        return map(_2split3, xs)


def f2product(two=None):
    '''Turn accurate I{TwoProduct} multiplication on or off.

       @kwarg two: If C{True}, turn I{TwoProduct} on, if C{False} off or
                   if C{None} or omitted, keep the current setting.

       @return: The previous setting (C{bool}).

       @see: I{TwoProduct} multiplication is based on the I{TwoProductFMA}
             U{Algorithm 3.5 <https://www.TUHH.De/ti3/paper/rump/OgRuOi05.pdf>}
             using function C{math.fma} from Python 3.13 and later or an
             equivalent, slower implementation when not available.
    '''
    t = Fsum._f2product
    if two is not None:
        Fsum._f2product = bool(two)
    return t


def _Fsumf_(*xs):  # in .auxLat, ...
    '''(INTERNAL) An C{Fsum(xs)}, all C{scalar}, an L{Fsum} or L{Fsum2Tuple}.
    '''
    return Fsum()._facc_scalarf(xs, up=False)


def _Fsum1f_(*xs):  # in .albers
    '''(INTERNAL) An C{Fsum(xs)}, all C{scalar}, an L{Fsum} or L{Fsum2Tuple}, 1-primed.
    '''
    return Fsum()._facc_scalarf(_1primed(xs), origin=-1, up=False)


def _halfeven(s, r, p):
    '''(INTERNAL) Round half-even.
    '''
    if (p > 0 and r > 0) or \
       (p < 0 and r < 0):  # signs match
        r *= 2
        t  = s + r
        if r == (t - s):
            s = t
    return s


def _isFsum(x):  # in .fmath
    '''(INTERNAL) Is C{x} an C{Fsum} instance?
    '''
    return isinstance(x, Fsum)


def _isFsum_2Tuple(x):  # in .basics, .constants, .fmath, .fstats
    '''(INTERNAL) Is C{x} an C{Fsum} or C{Fsum2Tuple} instance?
    '''
    return isinstance(x, _Fsum_2Tuple_types)


def _isOK(unused):
    '''(INTERNAL) Helper for C{Fsum._fsum2} and C{Fsum.nonfinites}.
    '''
    return True


def _isOK_or_finite(x, _isfine=_isfinite):
    '''(INTERNAL) Is C{x} finite or is I{non-finite} OK?
    '''
    # assert _isfine in (_isOK, _isfinite)
    return _isfine(x)  # C{bool}


try:
    from math import gcd as _gcd

    def _n_d2(n, d):
        '''(INTERNAL) Reduce C{n} and C{d} by C{gcd}.
        '''
        if n and d:
            try:
                c = _gcd(n, d)
                if c > 1:
                    n, d = (n // c), (d // c)
            except TypeError:  # non-int float
                pass
        return n, d

except ImportError:  # 3.4-

    def _n_d2(*n_d):  # PYCHOK redef
        return n_d


def _nfError(x, *args):
    '''(INTERNAL) Throw a C{not-finite} exception.
    '''
    E = _NonfiniteError(x)
    t =  Fmt.PARENSPACED(_not_finite_, x)
    if args:  # in _fmaX, _2sum
        return E(txt=t, *args)
    raise E(t, txt=None)


def _NonfiniteError(x):
    '''(INTERNAL) Return the Error class for C{x}, I{non-finite}.
    '''
    return _OverflowError if isinf(x) else (
           _ValueError    if isnan(x) else _AssertionError)


def nonfiniterrors(raiser=None):
    '''Throw C{OverflowError} and C{ValueError} exceptions for or
       handle I{non-finite} C{float}s as C{inf}, C{INF}, C{NINF},
       C{nan} and C{NAN} in summations and multiplications.

       @kwarg raiser: If C{True}, throw exceptions, if C{False} handle
                      I{non-finites} or if C{None} or omitted, leave
                      the setting unchanged.

       @return: Previous setting (C{bool}).

       @note: C{inf}, C{INF} and C{NINF} throw an C{OverflowError},
              C{nan} and C{NAN} a C{ValueError}.
    '''
    d = Fsum._isfine
    if raiser is not None:
        Fsum._isfine = {} if bool(raiser) else Fsum._nonfinites_isfine_kwds[True]
    return (False if d is Fsum._nonfinites_isfine_kwds[True] else
         _xkwds_get1(d, _isfine=_isfinite) is _isfinite) if d else True


def _1primed(xs):  # in .fmath
    '''(INTERNAL) 1-Primed summation of iterable C{xs}
       items, all I{known} to be C{scalar}.
    '''
    yield _1_0
    for x in xs:
        yield x
    yield _N_1_0


def _psum(ps, **_isfine):  # PYCHOK used!
    '''(INTERNAL) Partials summation, updating C{ps}.
    '''
    # assert isinstance(ps, list)
    i =  len(ps) - 1
    s = _0_0 if i < 0 else ps[i]
    while i > 0:
        i -= 1
        s, r = _2sum(s, ps[i], **_isfine)
        if r:  # sum(ps) became inexact
            if s:
                ps[i:] = r, s
                if i > 0:
                    s = _halfeven(s, r, ps[i-1])
                break  # return s
            s = r  # PYCHOK no cover
        elif not _isfinite(s):  # non-finite OK
            i = 0  # collapse ps
            if ps:
                s += sum(ps)
        ps[i:] = s,
    return s


def _Psum(ps, **name_f2product_nonfinites_RESIDUAL):
    '''(INTERNAL) Return an C{Fsum} from I{ordered} partials C{ps}.
    '''
    F = Fsum(**name_f2product_nonfinites_RESIDUAL)
    if ps:
        F._ps[:] = ps
        F._n = len(F._ps)
    return F


def _Psum_(*ps, **name_f2product_nonfinites_RESIDUAL):  # in .fmath
    '''(INTERNAL) Return an C{Fsum} from I{known scalar} C{ps}.
    '''
    return _Psum(ps, **name_f2product_nonfinites_RESIDUAL)


def _residue(other):
    '''(INTERNAL) Return the C{residual} or C{None} for C{scalar}.
    '''
    try:
        r = other.residual
    except AttributeError:
        r = None  # float, int, other
    return r


def _s_r2(s, r):
    '''(INTERNAL) Return C{(s, r)}, I{ordered}.
    '''
    if _isfinite(s):
        if r:
            if fabs(s) < fabs(r):
                s, r = r, (s or INT0)
        else:
            r = INT0
    else:
        r = _NONFINITEr
    return s, r


def _strcomplex(s, *args):
    '''(INTERNAL) C{Complex} 2- or 3-arg C{pow} error as C{str}.
    '''
    c = _strcomplex.__name__[4:]
    n = _sub_op_(len(args), _arg_)
    t =  unstr(pow, *args)
    return _SPACE_(c, s, _from_, n, t)


def _stresidual(prefix, residual, R=0, **mod_ratio):
    '''(INTERNAL) Residual error txt C{str}.
    '''
    p = _stresidual.__name__[3:]
    t =  Fmt.PARENSPACED(p, Fmt(residual))
    for n, v in itemsorted(mod_ratio):
        p =  Fmt.PARENSPACED(n, Fmt(v))
        t = _COMMASPACE_(t, p)
    return _SPACE_(prefix, t, Fmt.exceeds_R(R), _threshold_)


def _2sum(a, b, _isfine=_isfinite):  # in .testFmath
    '''(INTERNAL) Return C{a + b} as 2-tuple C{(sum, residual)} with finite C{sum},
       otherwise as 2-tuple C{(nonfinite, 0)} iff I{non-finites} are OK.
    '''
    # FastTwoSum U{Algorithm 1.1<https://www.TUHH.De/ti3/paper/rump/OgRuOi05.pdf>}

    # Neumaier, A. U{Rundungsfehleranalyse einiger Verfahren zur Summation endlicher
    # Summen<https://OnlineLibrary.Wiley.com/doi/epdf/10.1002/zamm.19740540106>},
    # 1974, Zeitschrift für Angewandte Mathmatik und Mechanik, vol 51, nr 1, p 39-51
    # <https://StackOverflow.com/questions/78633770/can-neumaier-summation-be-sped-up>
    s = a + b
    if _isfinite(s):
        if fabs(a) < fabs(b):
            r = (b - s) + a
        else:
            r = (a - s) + b
    elif _isfine(s):
        r = _NONFINITEr
    else:  # non-finite and not OK
        t =  unstr(_2sum, a, b)
        raise _nfError(s, t)
    return s, r


def _threshold(threshold=_0_0, **kwds):
    '''(INTERNAL) Get the L{ResidualError}s threshold,
       optionally from single kwds C{B{RESIDUAL}=scalar}.
    '''
    if kwds:
        threshold = _xkwds_get1(kwds, RESIDUAL=threshold)
    try:
        return _2finite(threshold)  # PYCHOK None
    except Exception as x:
        raise ResidualError(threshold=threshold, cause=x)


def _2tuple2(other):
    '''(INTERNAL) Return 2-tuple C{(other, r)} with C{other} as C{int},
       C{float} or C{as-is} and C{r} the residual of C{as-is} or 0.
    '''
    if _isFsum_2Tuple(other):
        s, r = other._fint2
        if r:
            s, r = other._nfprs2
            if r:  # PYCHOK no cover
                s = other  # L{Fsum} as-is
    else:
        r = 0
        s = other  # C{type} as-is
        if isint(s, both=True):
            s = int(s)
    return s, r


class Fsum(_Named):  # sync __methods__ with .vector3dBase.Vector3dBase, .fstats, ...
    '''Precision floating point summation, I{running} summation and accurate multiplication.

       Unlike Python's C{math.fsum}, this class accumulates values and provides intermediate,
       I{running}, precision floating point summations.  Accumulation may continue after any
       intermediate, I{running} summuation.

       @note: Values may be L{Fsum}, L{Fsum2Tuple}, C{int}, C{float} or C{scalar} instances,
              i.e. any C{type} having method C{__float__}.

       @note: Handling of I{non-finites} as C{inf}, C{INF}, C{NINF}, C{nan} and C{NAN} is
              determined by function L{nonfiniterrors<fsums.nonfiniterrors>} for the default
              and by method L{nonfinites<Fsum.nonfinites>} for individual C{Fsum} instances,
              overruling the default.  For backward compatibility, I{non-finites} raise
              exceptions by default.

       @see: U{Hettinger<https://GitHub.com/ActiveState/code/tree/master/recipes/Python/
             393090_Binary_floating_point_summatiaccurate_full/recipe-393090.py>},
             U{Kahan<https://WikiPedia.org/wiki/Kahan_summation_algorithm>}, U{Klein
             <https://Link.Springer.com/article/10.1007/s00607-005-0139-x>}, Python 2.6+
             file I{Modules/mathmodule.c} and the issue log U{Full precision summation
             <https://Bugs.Python.org/issue2819>}.

       @see: Method L{f2product<Fsum.f2product>} for details about accurate I{TwoProduct}
             multiplication.

       @see: Module L{fsums<pygeodesy.fsums>} for env variables C{PYGEODESY_FSUM_F2PRODUCT},
             C{PYGEODESY_FSUM_NONFINITES} and C{PYGEODESY_FSUM_RESIDUAL}.
    '''
    _f2product = _MODS.sys_version_info2 > (3, 12) or bool(_F2PRODUCT)
    _isfine    = {}  # == _isfinite, see nonfiniterrors()
    _n         =  0
#   _ps        = []  # partial sums
#   _ps_max    =  0  # max(Fsum._ps_max, len(Fsum._ps))  # 41
    _RESIDUAL  = _threshold(_RESIDUAL_0_0)

    def __init__(self, *xs, **name_f2product_nonfinites_RESIDUAL):
        '''New L{Fsum}.

           @arg xs: No, one or more initial items to accumulate (each C{scalar}, an
                    L{Fsum} or L{Fsum2Tuple}), all positional.
           @kwarg name_f2product_nonfinites_RESIDUAL: Optional C{B{name}=NN} (C{str})
                       and settings C{B{f2product}=None} (C{bool}), C{B{nonfinites}=None}
                       (C{bool}) and C{B{RESIDUAL}=0.0} threshold (C{scalar}) for this
                       L{Fsum}.

           @see: Methods L{Fsum.f2product}, L{Fsum.nonfinites}, L{Fsum.RESIDUAL},
                 L{Fsum.fadd} and L{Fsum.fadd_}.
        '''
        if name_f2product_nonfinites_RESIDUAL:
            self._optionals(**name_f2product_nonfinites_RESIDUAL)
        self._ps = []  # [_0_0], see L{Fsum._fprs}
        if xs:
            self._facc_args(xs, up=False)

    def __abs__(self):
        '''Return C{abs(self)} as an L{Fsum}.
        '''
        s = self.signOf()  # == self._cmp_0(0)
        return (-self) if s < 0 else self._copy_2(self.__abs__)

    def __add__(self, other):
        '''Return C{B{self} + B{other}} as an L{Fsum}.

           @arg other: An L{Fsum}, L{Fsum2Tuple} or C{scalar}.

           @return: The sum (L{Fsum}).

           @see: Methods L{Fsum.fadd_} and L{Fsum.fadd}.
        '''
        f = self._copy_2(self.__add__)
        return f._fadd(other)

    def __bool__(self):  # PYCHOK Python 3+
        '''Return C{bool(B{self})}, C{True} iff C{residual} is zero.
        '''
        s, r = self._nfprs2
        return bool(s or r) and s != -r  # == self != 0

    def __call__(self, other, **up):  # in .fmath
        '''Reset this C{Fsum} to C{other}, default C{B{up}=True}.
        '''
        self._ps[:] = 0,  # clear for errors
        self._fset(other, op=_fset_op_, **up)
        return self


    def __ceil__(self):  # PYCHOK not special in Python 2-
        '''Return this instance' C{math.ceil} as C{int} or C{float}.

           @return: An C{int} in Python 3+, but C{float} in Python 2-.

           @see: Methods L{Fsum.__floor__} and property L{Fsum.ceil}.
        '''
        return self.ceil

    def __cmp__(self, other):  # PYCHOK no cover
        '''Compare this with an other instance or C{scalar}, Python 2-.

           @return: -1, 0 or +1 (C{int}).

           @raise TypeError: Incompatible B{C{other}} C{type}.
        '''
        s = self._cmp_0(other, self.cmp.__name__)
        return _signOf(s, 0)

    def __divmod__(self, other, **raiser_RESIDUAL):
        '''Return C{divmod(B{self}, B{other})} as a L{DivMod2Tuple}
           with quotient C{div} an C{int} in Python 3+ or C{float}
           in Python 2- and remainder C{mod} an L{Fsum} instance.

           @arg other: An L{Fsum}, L{Fsum2Tuple} or C{scalar} modulus.
           @kwarg raiser_RESIDUAL: Use C{B{raiser}=False} to ignore
                         L{ResidualError}s (C{bool}) and C{B{RESIDUAL}=scalar}
                         to override the current L{RESIDUAL<Fsum.RESIDUAL>}.

           @raise ResidualError: Non-zero, significant residual or invalid
                                 B{C{RESIDUAL}}.

           @see: Method L{Fsum.fdiv}.
        '''
        f = self._copy_2(self.__divmod__)
        return f._fdivmod2(other, _divmod_op_, **raiser_RESIDUAL)

    def __eq__(self, other):
        '''Return C{(B{self} == B{other})} as C{bool} where B{C{other}}
           is C{scalar}, an other L{Fsum} or L{Fsum2Tuple}.
        '''
        return self._cmp_0(other, _fset_op_ + _fset_op_) == 0

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

           @arg other: An L{Fsum}, L{Fsum2Tuple} or C{scalar} divisor.

           @return: The C{floor} quotient (L{Fsum}).

           @see: Methods L{Fsum.__ifloordiv__}.
        '''
        f = self._copy_2(self.__floordiv__)
        return f._floordiv(other, _floordiv_op_)

    def __ge__(self, other):
        '''Return C{(B{self} >= B{other})}, see C{__eq__}.
        '''
        return self._cmp_0(other, _gt_op_ + _fset_op_) >= 0

    def __gt__(self, other):
        '''Return C{(B{self} > B{other})}, see C{__eq__}.
        '''
        return self._cmp_0(other, _gt_op_) > 0

    def __hash__(self):  # PYCHOK no cover
        '''Return C{hash(B{self})} as C{float}.
        '''
        # @see: U{Notes for type implementors<https://docs.Python.org/
        #       3/library/numbers.html#numbers.Rational>}
        return hash(self.partials)  # tuple.__hash__()

    def __iadd__(self, other):
        '''Apply C{B{self} += B{other}} to this instance.

           @arg other: An L{Fsum}, L{Fsum2Tuple} or C{scalar} value or
                       an iterable of several of the former.

           @return: This instance, updated (L{Fsum}).

           @raise TypeError: Invalid B{C{other}}, not
                             C{scalar} nor L{Fsum}.

           @see: Methods L{Fsum.fadd_} and L{Fsum.fadd}.
        '''
        try:
            return self._fadd(other, op=_iadd_op_)
        except TypeError:
            pass
        _xiterable(other)
        return self._facc(other)

    def __ifloordiv__(self, other):
        '''Apply C{B{self} //= B{other}} to this instance.

           @arg other: An L{Fsum}, L{Fsum2Tuple} or C{scalar} divisor.

           @return: This instance, updated (L{Fsum}).

           @raise ResidualError: Non-zero, significant residual
                                 in B{C{other}}.

           @raise TypeError: Invalid B{C{other}} type.

           @raise ValueError: Invalid or I{non-finite} B{C{other}}.

           @raise ZeroDivisionError: Zero B{C{other}}.

           @see: Methods L{Fsum.__itruediv__}.
        '''
        return self._floordiv(other, _floordiv_op_ + _fset_op_)

    def __imatmul__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __imod__(self, other):
        '''Apply C{B{self} %= B{other}} to this instance.

           @arg other: An L{Fsum}, L{Fsum2Tuple} or C{scalar} modulus.

           @return: This instance, updated (L{Fsum}).

           @see: Method L{Fsum.__divmod__}.
        '''
        return self._fdivmod2(other, _mod_op_ + _fset_op_).mod

    def __imul__(self, other):
        '''Apply C{B{self} *= B{other}} to this instance.

           @arg other: An L{Fsum}, L{Fsum2Tuple} or C{scalar} factor.

           @return: This instance, updated (L{Fsum}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Invalid B{C{other}} type.

           @raise ValueError: Invalid or I{non-finite} B{C{other}}.
        '''
        return self._fmul(other, _mul_op_ + _fset_op_)

    def __int__(self):
        '''Return this instance as an C{int}.

        @see: Method L{Fsum.int_float} and properties L{Fsum.ceil}
              and L{Fsum.floor}.
        '''
        i, _ = self._fint2
        return i

    def __invert__(self):  # PYCHOK no cover
        '''Not implemented.'''
        # Luciano Ramalho, "Fluent Python", O'Reilly, 2nd Ed, 2022 p. 567
        return _NotImplemented(self)

    def __ipow__(self, other, *mod, **raiser_RESIDUAL):  # PYCHOK 2 vs 3 args
        '''Apply C{B{self} **= B{other}} to this instance.

           @arg other: The exponent (C{scalar}, an L{Fsum} or L{Fsum2Tuple}).
           @arg mod: Optional modulus (C{int} or C{None}) for the 3-argument
                     C{pow(B{self}, B{other}, B{mod})} version.
           @kwarg raiser_RESIDUAL: Use C{B{raiser}=False} to ignore
                         L{ResidualError}s (C{bool}) and C{B{RESIDUAL}=scalar}
                         to override the current L{RESIDUAL<Fsum.RESIDUAL>}.

           @return: This instance, updated (L{Fsum}).

           @note: If B{C{mod}} is given, the result will be an C{integer}
                  L{Fsum} in Python 3+ if this instance C{is_integer} or
                  set to C{as_integer} and B{C{mod}} is given and C{None}.

           @raise OverflowError: Partial C{2sum} overflow.

           @raise ResidualError: Invalid B{C{RESIDUAL}} or the residual
                                 is non-zero and significant and either
                                 B{C{other}} is a fractional or negative
                                 C{scalar} or B{C{mod}} is given and not
                                 C{None}.

           @raise TypeError: Invalid B{C{other}} type or 3-argument C{pow}
                             invocation failed.

           @raise ValueError: If B{C{other}} is a negative C{scalar} and this
                              instance is C{0} or B{C{other}} is a fractional
                              C{scalar} and this instance is negative or has a
                              non-zero and significant residual or B{C{mod}}
                              is given as C{0}.

           @see: CPython function U{float_pow<https://GitHub.com/
                 python/cpython/blob/main/Objects/floatobject.c>}.
        '''
        return self._fpow(other, _pow_op_ + _fset_op_, *mod, **raiser_RESIDUAL)

    def __isub__(self, other):
        '''Apply C{B{self} -= B{other}} to this instance.

           @arg other: An L{Fsum}, L{Fsum2Tuple} or C{scalar} value or
                       an iterable of several of the former.

           @return: This instance, updated (L{Fsum}).

           @raise TypeError: Invalid B{C{other}} type.

           @see: Methods L{Fsum.fsub_} and L{Fsum.fsub}.
        '''
        try:
            return self._fsub(other, _isub_op_)
        except TypeError:
            pass
        _xiterable(other)
        return self._facc_neg(other)

    def __iter__(self):
        '''Return an C{iter}ator over a C{partials} duplicate.
        '''
        return iter(self.partials)

    def __itruediv__(self, other, **raiser_RESIDUAL):
        '''Apply C{B{self} /= B{other}} to this instance.

           @arg other: An L{Fsum}, L{Fsum2Tuple} or C{scalar} divisor.
           @kwarg raiser_RESIDUAL: Use C{B{raiser}=False} to ignore
                         L{ResidualError}s (C{bool}) and C{B{RESIDUAL}=scalar}
                         to override the current L{RESIDUAL<Fsum.RESIDUAL>}.

           @return: This instance, updated (L{Fsum}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise ResidualError: Non-zero, significant residual or invalid
                                 B{C{RESIDUAL}}.

           @raise TypeError: Invalid B{C{other}} type.

           @raise ValueError: Invalid or I{non-finite} B{C{other}}.

           @raise ZeroDivisionError: Zero B{C{other}}.

           @see: Method L{Fsum.__ifloordiv__}.
        '''
        return self._ftruediv(other, _truediv_op_ + _fset_op_, **raiser_RESIDUAL)

    def __le__(self, other):
        '''Return C{(B{self} <= B{other})}, see C{__eq__}.
        '''
        return self._cmp_0(other, _lt_op_ + _fset_op_) <= 0

    def __len__(self):
        '''Return the number of values accumulated (C{int}).
        '''
        return self._n

    def __lt__(self, other):
        '''Return C{(B{self} < B{other})}, see C{__eq__}.
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
        '''Return C{(B{self} != B{other})}, see C{__eq__}.
        '''
        return self._cmp_0(other, _ne_op_) != 0

    def __neg__(self):
        '''Return C{copy(B{self})}, I{negated}.
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
        return f._fadd(self)

    def __rdivmod__(self, other):
        '''Return C{divmod(B{other}, B{self})} as 2-tuple
           C{(quotient, remainder)}.

           @see: Method L{Fsum.__divmod__}.
        '''
        f = self._copy_2r(other, self.__rdivmod__)
        return f._fdivmod2(self, _divmod_op_)

#       turned off, called by _deepcopy and _copy
#   def __reduce__(self):  # Python 3.8+
#       ''' Pickle, like std C{fractions.Fraction}, see U{__reduce__
#           <https://docs.Python.org/3/library/pickle.html#object.__reduce__>}
#       '''
#       dict_ = self._Fsum_as().__dict__  # no __setstate__
#       return (self.__class__, self.partials, dict_)

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

    def __rmatmul__(self, other):  # PYCHOK no coveS
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

    def __round__(self, *ndigits):  # PYCHOK Python 3+
        '''Return C{round(B{self}, *B{ndigits}} as an L{Fsum}.

           @arg ndigits: Optional number of digits (C{int}).
        '''
        f = self._copy_2(self.__round__)
        # <https://docs.Python.org/3.12/reference/datamodel.html?#object.__round__>
        return f._fset(round(float(self), *ndigits))  # can be C{int}

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

    def __rtruediv__(self, other, **raiser_RESIDUAL):
        '''Return C{B{other} / B{self}} as an L{Fsum}.

           @see: Method L{Fsum.__itruediv__}.
        '''
        f = self._copy_2r(other, self.__rtruediv__)
        return f._ftruediv(self, _truediv_op_, **raiser_RESIDUAL)

    def __str__(self):
        '''Return the default C{str(self)}.
        '''
        return self.toStr(lenc=True)

    def __sub__(self, other):
        '''Return C{B{self} - B{other}} as an L{Fsum}.

           @arg other: An L{Fsum}, L{Fsum2Tuple} or C{scalar}.

           @return: The difference (L{Fsum}).

           @see: Method L{Fsum.__isub__}.
        '''
        f = self._copy_2(self.__sub__)
        return f._fsub(other, _sub_op_)

    def __truediv__(self, other, **raiser_RESIDUAL):
        '''Return C{B{self} / B{other}} as an L{Fsum}.

           @arg other: An L{Fsum}, L{Fsum2Tuple} or C{scalar} divisor.
           @kwarg raiser_RESIDUAL: Use C{B{raiser}=False} to ignore
                         L{ResidualError}s (C{bool}) and C{B{RESIDUAL}=scalar}
                         to override the current L{RESIDUAL<Fsum.RESIDUAL>}.

           @return: The quotient (L{Fsum}).

           @raise ResidualError: Non-zero, significant residual or invalid
                                 B{C{RESIDUAL}}.

           @see: Method L{Fsum.__itruediv__}.
        '''
        return self._truediv(other, _truediv_op_, **raiser_RESIDUAL)

    __trunc__ = __int__

    if _MODS.sys_version_info2 < (3, 0):  # PYCHOK no cover
        # <https://docs.Python.org/2/library/operator.html#mapping-operators-to-functions>
        __div__     = __truediv__
        __idiv__    = __itruediv__
        __long__    = __int__
        __nonzero__ = __bool__
        __rdiv__    = __rtruediv__

    def as_integer_ratio(self):
        '''Return this instance as the ratio of 2 integers.

           @return: 2-Tuple C{(numerator, denominator)} both C{int} with
                    C{numerator} signed and C{denominator} non-zero and
                    positive.  The C{numerator} is I{non-finite} if this
                    instance is.

           @see: Method L{Fsum.fint2} and C{float.as_integer_ratio} in
                 Python 2.7+.
        '''
        n, r = self._fint2
        if r:
            i, d =  float(r).as_integer_ratio()
            n, d = _n_d2(n * d + i, d)
        else:  # PYCHOK no cover
            d = 1
        return n, d

    @property_RO
    def as_iscalar(self):
        '''Get this instance I{as-is} (L{Fsum} with C{non-zero residual},
           C{scalar} or I{non-finite}).
        '''
        s, r = self._nfprs2
        return self if r else s

    @property_RO
    def ceil(self):
        '''Get this instance' C{ceil} value (C{int} in Python 3+, but
           C{float} in Python 2-).

           @note: This C{ceil} takes the C{residual} into account.

           @see: Method L{Fsum.int_float} and properties L{Fsum.floor},
                 L{Fsum.imag} and L{Fsum.real}.
        '''
        s, r = self._fprs2
        c = _ceil(s) + int(r) - 1
        while r > (c - s):  # (s + r) > c
            c += 1
        return c  # _ceil(self._n_d)

    cmp = __cmp__

    def _cmp_0(self, other, op):
        '''(INTERNAL) Return C{scalar(self - B{other})} for 0-comparison.
        '''
        if _isFsum_2Tuple(other):
            s = self._ps_1sum(*other._ps)
        elif self._scalar(other, op):
            s = self._ps_1sum(other)
        else:
            s = self.signOf()  # res=True
        return s

    def copy(self, deep=False, **name):
        '''Copy this instance, C{shallow} or B{C{deep}}.

           @kwarg name: Optional, overriding C{B{name}='"copy"} (C{str}).

           @return: The copy (L{Fsum}).
         '''
        n = _name__(name, name__=self.copy)
        f = _Named.copy(self, deep=deep, name=n)
        if f._ps is self._ps:
            f._ps = list(self._ps)  # separate list
        if not deep:
            f._n = 1
        # assert f._f2product == self._f2product
        # assert f._Fsum is f
        # assert f._isfine is self._isfine
        # assert f._RESIDUAL is self._RESIDUAL
        return f

    def _copy_2(self, which, name=NN):
        '''(INTERNAL) Copy for I{dyadic} operators.
        '''
        n =  name or which.__name__  # _DUNDER_nameof
        # NOT .classof due to .Fdot(a, *b) args, etc.
        f = _Named.copy(self, deep=False, name=n)
        f._ps = list(self._ps)  # separate list
        # assert f._n == self._n
        # assert f._f2product == self._f2product
        # assert f._Fsum is f
        # assert f._isfine is self._isfine
        # assert f._RESIDUAL is self._RESIDUAL
        return f

    def _copy_2r(self, other, which):
        '''(INTERNAL) Copy for I{reverse-dyadic} operators.
        '''
        return other._copy_2(which) if _isFsum(other) else \
                self._copy_2(which)._fset(other)

    divmod = __divmod__

    def _Error(self, op, other, Error, **txt_cause):
        '''(INTERNAL) Format an B{C{Error}} for C{{self} B{op} B{other}}.
        '''
        # self.as_iscalar causes RecursionError for ._fprs2 errors
        s = _Psum(self._ps, nonfinites=True, name=self.name)
        return Error(_SPACE_(s.as_iscalar, op, other), **txt_cause)

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
        u = unstr(self.named3, *xs, _ELLIPSIS=4, **kwds)
        return E(u, txt=t, cause=X)

    def _facc(self, xs, up=True, **_X_x_origin):
        '''(INTERNAL) Accumulate more C{scalar}s or L{Fsum}s.
        '''
        if xs:
            kwds = self._isfine
            if _X_x_origin:
                kwds = _xkwds(_X_x_origin, **kwds)
            fs    = _xs(xs, **kwds)  # PYCHOK yield
            ps    =  self._ps
            ps[:] =  self._ps_acc(list(ps), fs, up=up)
#           if len(ps) > 16:
#               _ = _psum(ps, **self._isfine)
        return self

    def _facc_args(self, xs, **up):
        '''(INTERNAL) Accumulate 0, 1 or more C{xs}, all positional
           arguments in the caller of this method.
        '''
        return self._fadd(xs[0], **up) if len(xs) == 1 else \
               self._facc(xs, **up)  # origin=1?

    def _facc_dot(self, n, xs, ys, **kwds):  # in .fmath
        '''(INTERNAL) Accumulate C{fdot(B{xs}, *B{ys})}.
        '''
        if n > 0:
            _f = Fsum(**kwds)
            self._facc(_f(x).fmul(y) for x, y in zip(xs, ys))  # PYCHOK attr?
        return self

    def _facc_neg(self, xs, **up_origin):
        '''(INTERNAL) Accumulate more C{xs}, negated.
        '''
        def _N(X):
            return X._ps_neg

        def _n(x):
            return -float(x)

        return self._facc(xs, _X=_N, _x=_n, **up_origin)

    def _facc_power(self, power, xs, which, **raiser_RESIDUAL):  # in .fmath
        '''(INTERNAL) Add each C{xs} as C{float(x**power)}.
        '''
        def _Pow4(p):
            r = 0
            if _isFsum_2Tuple(p):
                s, r = p._fprs2
                if r:
                    m = Fsum._pow
                else:  # scalar
                    return _Pow4(s)
            elif isint(p, both=True) and int(p) >= 0:
                p = s = int(p)
                m = Fsum._pow_int
            else:
                p = s = _2float(power=p, **self._isfine)
                m = Fsum._pow_scalar
            return m, p, s, r

        _Pow, p, s, r = _Pow4(power)
        if p:  # and xs:
            op   =  which.__name__
            _FsT = _Fsum_2Tuple_types
            _pow =  self._pow_2_3

            def _P(X):
                f = _Pow(X, p, power, op, **raiser_RESIDUAL)
                return f._ps if isinstance(f, _FsT) else (f,)

            def _p(x):
                x =  float(x)
                f = _pow(x, s, power, op, **raiser_RESIDUAL)
                if f and r:
                    f *= _pow(x, r, power, op, **raiser_RESIDUAL)
                return f

            f = self._facc(xs, _X=_P, _x=_p)  # origin=1?
        else:
            f = self._facc_scalar_(float(len(xs)))  # x**0 == 1
        return f

    def _facc_scalar(self, xs, **up):
        '''(INTERNAL) Accumulate all C{xs}, each C{scalar}.
        '''
        if xs:
            ps    = self._ps
            ps[:] = self._ps_acc(list(ps), xs, **up)
        return self

    def _facc_scalar_(self, *xs, **up):
        '''(INTERNAL) Accumulate all positional C{xs}, each C{scalar}.
        '''
        return self._facc_scalar(xs, **up)

    def _facc_scalarf(self, xs, up=True, **origin_which):
        '''(INTERNAL) Accumulate all C{xs}, each C{scalar}, an L{Fsum} or
           L{Fsum2Tuple}, like function C{_xsum}.
        '''
        _C =  self.__class__
        fs = _xs(xs, **_x_isfine(self.nonfinitesOK, _Cdot=_C,
                               **origin_which))  # PYCHOK yield
        return self._facc_scalar(fs, up=up)

#   def _facc_up(self, up=True):
#       '''(INTERNAL) Update the C{partials}, by removing
#          and re-accumulating the final C{partial}.
#       '''
#       ps = self._ps
#       while len(ps) > 1:
#           p = ps.pop()
#           if p:
#               n = self._n
#               _ = self._ps_acc(ps, (p,), up=False)
#               self._n = n
#               break
#       return self._update() if up else self

    def fadd(self, xs=()):
        '''Add an iterable's items to this instance.

           @arg xs: Iterable of items to add (each C{scalar},
                    an L{Fsum} or L{Fsum2Tuple}).

           @return: This instance (L{Fsum}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: An invalid B{C{xs}} item.

           @raise ValueError: Invalid or I{non-finite} B{C{xs}} value.
        '''
        if _isFsum_2Tuple(xs):
            self._facc_scalar(xs._ps)
        elif isscalar(xs):    # for backward compatibility # PYCHOK no cover
            x = _2float(x=xs, **self._isfine)
            self._facc_scalar_(x)
        elif xs:  # _xiterable(xs)
            self._facc(xs)
        return self

    def fadd_(self, *xs):
        '''Add all positional items to this instance.

           @arg xs: Values to add (each C{scalar}, an L{Fsum}
                    or L{Fsum2Tuple}), all positional.

           @see: Method L{Fsum.fadd} for further details.
        '''
        return self._facc_args(xs)

    def _fadd(self, other, op=_add_op_, **up):
        '''(INTERNAL) Apply C{B{self} += B{other}}.
        '''
        if _isFsum_2Tuple(other):
            self._facc_scalar(other._ps, **up)
        elif self._scalar(other, op):
            self._facc_scalar_(other, **up)
        return self

    fcopy   =   copy  # for backward compatibility
    fdiv    = __itruediv__
    fdivmod = __divmod__

    def _fdivmod2(self, other, op, **raiser_RESIDUAL):
        '''(INTERNAL) Apply C{B{self} %= B{other}} and return a L{DivMod2Tuple}.
        '''
        # result mostly follows CPython function U{float_divmod
        # <https://GitHub.com/python/cpython/blob/main/Objects/floatobject.c>},
        # but at least divmod(-3, 2) equals Cpython's result (-2, 1).
        q = self._truediv(other, op, **raiser_RESIDUAL).floor
        if q:  # == float // other == floor(float / other)
            self -= self._Fsum_as(q) * other  # NOT other * q!

        s = signOf(other)  # make signOf(self) == signOf(other)
        if s and self.signOf() == -s:  # PYCHOK no cover
            self += other
            q -= 1
#           t  = self.signOf()
#           if t and t != s:
#               raise self._Error(op, other, _AssertionError, txt__=signOf)
        return DivMod2Tuple(q, self)  # q is C{int} in Python 3+, but C{float} in Python 2-

    def _fhorner(self, x, cs, where, incx=True):  # in .fmath
        '''(INTERNAL) Add an L{Fhorner} evaluation of polynomial
           C{sum(cs[i] * B{x}**i for i=0..len(cs)-1) if B{incx}
           else sum(... i=len(cs)-1..0)}.
        '''
        # assert _xiterablen(cs)
        try:
            n  = len(cs)
            H  = self._Fsum_as(name__=self._fhorner)
            _m = H._mul_Fsum if _isFsum_2Tuple(x) else \
                 H._mul_scalar
            if _2finite(x, **self._isfine) and n > 1:
                for c in (reversed(cs) if incx else cs):
                    H._fset(_m(x, _mul_op_), up=False)
                    H._fadd(c, up=False)
            else:  # x == 0
                H = cs[0] if n else 0
            self._fadd(H)
        except Exception as X:
            t = unstr(where, x, *cs, _ELLIPSIS=4, incx=incx)
            raise self._ErrorX(X, _add_op_, t)
        return self

    def _finite(self, other, op=None):
        '''(INTERNAL) Return B{C{other}} if C{finite}.
        '''
        if _isOK_or_finite(other, **self._isfine):
            return other
        E = _NonfiniteError(other)
        raise self._Error(op, other, E, txt=_not_finite_)

    def fint(self, name=NN, **raiser_RESIDUAL):
        '''Return this instance' current running sum as C{integer}.

           @kwarg name: Optional, overriding C{B{name}="fint"} (C{str}).
           @kwarg raiser_RESIDUAL: Use C{B{raiser}=False} to ignore
                         L{ResidualError}s (C{bool}) and C{B{RESIDUAL}=scalar}
                         to override the current L{RESIDUAL<Fsum.RESIDUAL>}.

           @return: The C{integer} sum (L{Fsum}) if this instance C{is_integer}
                    with a zero or insignificant I{integer} residual.

           @raise ResidualError: Non-zero, significant residual or invalid
                                 B{C{RESIDUAL}}.

           @see: Methods L{Fsum.fint2}, L{Fsum.int_float} and L{Fsum.is_integer}.
        '''
        i, r = self._fint2
        if r:
            R = self._raiser(r, i, **raiser_RESIDUAL)
            if R:
                t = _stresidual(_integer_, r, **R)
                raise ResidualError(_integer_, i, txt=t)
        return self._Fsum_as(i, name=_name__(name, name__=self.fint))

    def fint2(self, **name):
        '''Return this instance' current running sum as C{int} and the
           I{integer} residual.

           @kwarg name: Optional name (C{str}).

           @return: An L{Fsum2Tuple}C{(fsum, residual)} with C{fsum}
                    an C{int} and I{integer} C{residual} a C{float} or
                    C{INT0} if the C{fsum} is considered to be I{exact}.
                    The C{fsum} is I{non-finite} if this instance is.
        '''
        return Fsum2Tuple(*self._fint2, **name)

    @Property
    def _fint2(self):  # see ._fset
        '''(INTERNAL) Get 2-tuple (C{int}, I{integer} residual).
        '''
        s, r = self._nfprs2
        if _isfinite(s):
            i =  int(s)
            r = (self._ps_1sum(i) if len(self._ps) > 1 else
                 float(s - i)) or INT0
        else:  # INF, NAN, NINF
            i =  float(s)
#           r = _NONFINITEr
        return i, r  # Fsum2Tuple?

    @_fint2.setter_  # PYCHOK setter_UNDERscore!
    def _fint2(self, s):  # in _fset
        '''(INTERNAL) Replace the C{_fint2} value.
        '''
        if _isfinite(s):
            i =  int(s)
            r = (s - i) or INT0
        else:  # INF, NAN, NINF
            i =  float(s)
            r = _NONFINITEr
        return i, r  # like _fint2.getter

    @deprecated_property_RO
    def float_int(self):  # PYCHOK no cover
        '''DEPRECATED, use method C{Fsum.int_float}.'''
        return self.int_float()  # raiser=False

    @property_RO
    def floor(self):
        '''Get this instance' C{floor} (C{int} in Python 3+, but
           C{float} in Python 2-).

           @note: This C{floor} takes the C{residual} into account.

           @see: Method L{Fsum.int_float} and properties L{Fsum.ceil},
                 L{Fsum.imag} and L{Fsum.real}.
        '''
        s, r = self._fprs2
        f = _floor(s) + _floor(r) + 1
        while (f - s) > r:  # f > (s + r)
            f -= 1
        return f  # _floor(self._n_d)

#   ffloordiv = __ifloordiv__  # for naming consistency?
#   floordiv  = __floordiv__   # for naming consistency?

    def _floordiv(self, other, op, **raiser_RESIDUAL):  # rather _ffloordiv?
        '''Apply C{B{self} //= B{other}}.
        '''
        q = self._ftruediv(other, op, **raiser_RESIDUAL)  # == self
        return self._fset(q.floor)  # floor(q)

    def fma(self, other1, other2, **nonfinites):  # in .fmath.fma
        '''Fused-multiply-add C{self *= B{other1}; self += B{other2}}.

           @arg other1: Multiplier (C{scalar}, an L{Fsum} or L{Fsum2Tuple}).
           @arg other2: Addend (C{scalar}, an L{Fsum} or L{Fsum2Tuple}).
           @kwarg nonfinites: Use C{B{nonfinites}=True} or C{False}, to
                              override L{nonfinites<Fsum.nonfinites>} and
                              L{nonfiniterrors} default (C{bool}).
        '''
        op  = self.fma.__name__
        _fs = self._ps_other
        try:
            s, r = self._fprs2
            if r:
                f  =  self._f2mul(self.fma, (other1,), **nonfinites)
                f +=  other2
            elif _residue(other1) or _residue(other2):
                fs = _2split3s(_fs(op, other1))
                fs = _2products(s, fs, *_fs(op, other2))
                f  = _Psum(self._ps_acc([], fs, up=False), name=op)
            else:
                f  = _fma(s, other1, other2)
                f  = _2finite(f, **self._isfine)
        except TypeError as X:
            raise self._ErrorX(X, op, (other1, other2))
        except (OverflowError, ValueError) as X:  # from math.fma
            f  = self._mul_reduce(s, other1)  # INF, NAN, NINF
            f += sum(_fs(op, other2))
            f  = self._nonfiniteX(X, op, f, **nonfinites)
        return self._fset(f)

    fmul = __imul__

    def _fmul(self, other, op):
        '''(INTERNAL) Apply C{B{self} *= B{other}}.
        '''
        if _isFsum_2Tuple(other):
            if len(self._ps) != 1:
                f = self._mul_Fsum(other, op)
            elif len(other._ps) != 1:  # and len(self._ps) == 1
                f = self._ps_mul(op, *other._ps) if other._ps else _0_0
            elif self._f2product:  # len(other._ps) == 1
                f = self._mul_scalar(other._ps[0], op)
            else:  # len(other._ps) == len(self._ps) == 1
                f = self._finite(self._ps[0] * other._ps[0], op=op)
        else:
            s = self._scalar(other, op)
            f = self._mul_scalar(s, op)
        return self._fset(f)  # n=len(self) + 1

    @deprecated_method
    def f2mul(self, *others, **raiser):
        '''DEPRECATED on 2024.09.13, use method L{f2mul_<Fsum.f2mul_>}.'''
        return self._fset(self.f2mul_(others, **raiser))

    def f2mul_(self, *others, **f2product_nonfinites):  # in .fmath.f2mul
        '''Return C{B{self} * B{other} * B{other} ...} for all B{C{others}} using cascaded,
           accurate multiplication like with L{f2product<Fsum.f2product>}C{(B{True})}.

           @arg others: Multipliers (each C{scalar}, an L{Fsum} or L{Fsum2Tuple}), all
                        positional.
           @kwarg f2product_nonfinites: Use C{B{f2product=False}} to override the default
                            C{True} and C{B{nonfinites}=True} or C{False}, to override
                            settings L{nonfinites<Fsum.nonfinites>} and L{nonfiniterrors}.

           @return: The cascaded I{TwoProduct} (L{Fsum} or C{float}).

           @see: U{Equations 2.3<https://www.TUHH.De/ti3/paper/rump/OzOgRuOi06.pdf>}
        '''
        return self._f2mul(self.f2mul_, others, **f2product_nonfinites)

    def _f2mul(self, where, others, f2product=True, **nonfinites_raiser):
        '''(INTERNAL) See methods C{fma} and C{f2mul_}.
        '''
        f = _Psum(self._ps, f2product=f2product, name=where.__name__)
        if others and f:
            if f.f2product():
                def _pfs(f, ps):
                    return _2products(f, _2split3s(ps))
            else:
                def _pfs(f, ps):  # PYCHOK redef
                    return (f * p for p in ps)

            op, ps = where.__name__, f._ps
            try:  # as if self.f2product(True)
                for other in others:  # to pinpoint errors
                    for p in self._ps_other(op, other):
                        ps[:] = f._ps_acc([], _pfs(p, ps), update=False)
                f._update()
            except TypeError as X:
                raise self._ErrorX(X, op, other)
            except (OverflowError, ValueError) as X:
                r = self._mul_reduce(sum(ps), other)  # INF, NAN, NINF
                r = self._nonfiniteX(X, op, r, **nonfinites_raiser)
                f._fset(r)
        return f

    def fover(self, over, **raiser_RESIDUAL):
        '''Apply C{B{self} /= B{over}} and summate.

           @arg over: An L{Fsum} or C{scalar} denominator.
           @kwarg raiser_RESIDUAL: Use C{B{raiser}=False} to ignore
                         L{ResidualError}s (C{bool}) and C{B{RESIDUAL}=scalar}
                         to override the current L{RESIDUAL<Fsum.RESIDUAL>}.

           @return: Precision running sum (C{float}).

           @raise ResidualError: Non-zero, significant residual or invalid
                                 B{C{RESIDUAL}}.

           @see: Methods L{Fsum.fsum} and L{Fsum.__itruediv__}.
        '''
        return float(self.fdiv(over, **raiser_RESIDUAL)._fprs)

    fpow = __ipow__

    def _fpow(self, other, op, *mod, **raiser_RESIDUAL):
        '''Apply C{B{self} **= B{other}}, optional B{C{mod}} or C{None}.
        '''
        if mod:
            if mod[0] is not None:  # == 3-arg C{pow}
                f = self._pow_2_3(self, other, other, op, *mod, **raiser_RESIDUAL)
            elif self.is_integer():
                # return an exact C{int} for C{int}**C{int}
                i, _ =  self._fint2  # assert _ == 0
                x, r = _2tuple2(other)  # C{int}, C{float} or other
                f = self._Fsum_as(i)._pow_Fsum(other, op, **raiser_RESIDUAL) if r else \
                    self._pow_2_3(i, x,        other, op, **raiser_RESIDUAL)
            else:  # mod[0] is None, power(self, other)
                f = self._pow(other, other, op, **raiser_RESIDUAL)
        else:  # pow(self, other)
            f = self._pow(other, other, op, **raiser_RESIDUAL)
        return self._fset(f)  # n=max(len(self), 1)

    def f2product(self, *two):
        '''Get and set accurate I{TwoProduct} multiplication for this
           L{Fsum}, overriding the L{f2product} default.

           @arg two: If omitted, leave the override unchanged, if C{True},
                     turn I{TwoProduct} on, if C{False} off, if C{None}e
                     remove th override (C{bool} or C{None}).

           @return: The previous setting (C{bool} or C{None} if not set).

           @see: Function L{f2product<fsums.f2product>}.

           @note: Use C{f.f2product() or f2product()} to determine whether
                  multiplication is accurate for L{Fsum} C{f}.
        '''
        if two:  # delattrof(self, _f2product=None)
            t = _xkwds_pop(self.__dict__, _f2product=None)
            if two[0] is not None:
                self._f2product = bool(two[0])
        else:  # getattrof(self, _f2product=None)
            t = _xkwds_get(self.__dict__, _f2product=None)
        return t

    @Property
    def _fprs(self):
        '''(INTERNAL) Get and cache this instance' precision
           running sum (C{float} or C{int}), ignoring C{residual}.

           @note: The precision running C{fsum} after a C{//=} or
                  C{//} C{floor} division is C{int} in Python 3+.
        '''
        s, _ = self._fprs2
        return s  # ._fprs2.fsum

    @_fprs.setter_  # PYCHOK setter_UNDERscore!
    def _fprs(self, s):
        '''(INTERNAL) Replace the C{_fprs} value.
        '''
        return s

    @Property
    def _fprs2(self):
        '''(INTERNAL) Get and cache this instance' precision
           running sum and residual (L{Fsum2Tuple}).
        '''
        ps = self._ps
        n  = len(ps)
        try:
            if n > 2:
                s = _psum(ps, **self._isfine)
                if not _isfinite(s):
                    ps[:] = s,  # collapse ps
                    return Fsum2Tuple(s, _NONFINITEr)
                n = len(ps)
#               Fsum._ps_max = max(Fsum._ps_max, n)
                if n > 2:
                    r = self._ps_1sum(s)
                    return Fsum2Tuple(*_s_r2(s, r))
            if n > 1:  # len(ps) == 2
                s, r  = _s_r2(*_2sum(*ps, **self._isfine))
                ps[:] = (r, s) if r else (s,)
            elif ps:  # len(ps) == 1
                s = ps[0]
                r = INT0 if _isfinite(s) else _NONFINITEr
            else:  # len(ps) == 0
                s = _0_0
                r = INT0 if _isfinite(s) else _NONFINITEr
                ps[:] = s,
        except (OverflowError, ValueError) as X:
            op = _fset_op_  # INF, NAN, NINF
            ps[:] = sum(ps),  # collapse ps
            s =  self._nonfiniteX(X, op, ps[0])
            r = _NONFINITEr
        # assert self._ps is ps
        return Fsum2Tuple(s, r)

    @_fprs2.setter_  # PYCHOK setter_UNDERscore!
    def _fprs2(self, s_r):
        '''(INTERNAL) Replace the C{_fprs2} value.
        '''
        return Fsum2Tuple(s_r)

    def fset_(self, *xs):
        '''Apply C{B{self}.partials = Fsum(*B{xs}).partials}.

           @arg xs: Optional, new values (each C{scalar} or an L{Fsum}
                    or L{Fsum2Tuple} instance), all positional.

           @return: This instance, replaced (C{Fsum}).

           @see: Method L{Fsum.fadd} for further details.
        '''
        f = (xs[0] if xs else _0_0) if len(xs) < 2 else \
            Fsum(*xs, nonfinites=self.nonfinites())  # self._Fsum_as(*xs)
        return self._fset(f, op=_fset_op_)

    def _fset(self, other, n=0, up=True, **op):
        '''(INTERNAL) Overwrite this instance with an other or a C{scalar}.
        '''
        if other is self:
            pass  # from ._fmul, ._ftruediv and ._pow_0_1
        elif _isFsum_2Tuple(other):
            if op:  # and not self.nonfinitesOK:
                self._finite(other._fprs, **op)
            self._ps[:] = other._ps
            self._n     = n or other._n
            if up:  # use or zap the C{Property_RO} values
                Fsum._fint2._update_from(self, other)
                Fsum._fprs ._update_from(self, other)
                Fsum._fprs2._update_from(self, other)
        elif isscalar(other):
            s = float(self._finite(other, **op)) if op else other
            self._ps[:] = s,
            self._n     = n or 1
            if up:  # Property _fint2, _fprs and _fprs2 all have
                # @.setter_underscore and NOT @.setter because the
                # latter's _fset zaps the value set by @.setter
                self._fint2 = s
                self._fprs  = s
                self._fprs2 = s, INT0
                # assert self._fprs is s
        else:
            op = _xkwds_get1(op, op=_fset_op_)
            raise self._Error(op, other, _TypeError)
        return self

    def fsub(self, xs=()):
        '''Subtract an iterable's items from this instance.

           @see: Method L{Fsum.fadd} for further details.
        '''
        return self._facc_neg(xs)

    def fsub_(self, *xs):
        '''Subtract all positional items from this instance.

           @see: Method L{Fsum.fadd_} for further details.
        '''
        return self._fsub(xs[0], _sub_op_) if len(xs) == 1 else \
               self._facc_neg(xs)  # origin=1?

    def _fsub(self, other, op):
        '''(INTERNAL) Apply C{B{self} -= B{other}}.
        '''
        if _isFsum_2Tuple(other):
            if other is self:  # or other._fprs2 == self._fprs2:
                self._fset(_0_0, n=len(self) * 2)
            elif other._ps:
                self._facc_scalar(other._ps_neg)
        elif self._scalar(other, op):
            self._facc_scalar_(-other)
        return self

    def fsum(self, xs=()):
        '''Add an iterable's items, summate and return the current
           precision running sum.

           @arg xs: Iterable of items to add (each item C{scalar},
                    an L{Fsum} or L{Fsum2Tuple}).

           @return: Precision running sum (C{float} or C{int}).

           @see: Method L{Fsum.fadd}.

           @note: Accumulation can continue after summation.
        '''
        return self._facc(xs)._fprs

    def fsum_(self, *xs):
        '''Add any positional items, summate and return the current
           precision running sum.

           @arg xs: Items to add (each C{scalar}, an L{Fsum} or
                    L{Fsum2Tuple}), all positional.

           @return: Precision running sum (C{float} or C{int}).

           @see: Methods L{Fsum.fsum}, L{Fsum.Fsum_} and L{Fsum.fsumf_}.
        '''
        return self._facc_args(xs)._fprs

    def Fsum_(self, *xs, **name):
        '''Like method L{Fsum.fsum_} but returning a named L{Fsum}.

           @kwarg name: Optional name (C{str}).

           @return: Copy of this updated instance (L{Fsum}).
        '''
        return self._facc_args(xs)._copy_2(self.Fsum_, **name)

    def Fsum2Tuple_(self, *xs, **name):
        '''Like method L{Fsum.fsum_} but returning a named L{Fsum2Tuple}.

           @kwarg name: Optional name (C{str}).

           @return: Precision running sum (L{Fsum2Tuple}).
        '''
        return Fsum2Tuple(self._facc_args(xs)._nfprs2, **name)

    @property_RO
    def _Fsum(self):  # like L{Fsum2Tuple._Fsum}, in .fstats
        return self  # NOT @Property_RO, see .copy and ._copy_2

    def _Fsum_as(self, *xs, **name_f2product_nonfinites_RESIDUAL):
        '''(INTERNAL) Return an C{Fsum} with this C{Fsum}'s C{.f2product},
            C{.nonfinites} and C{.RESIDUAL} setting, optionally
            overridden with C{name_f2product_nonfinites_RESIDUAL} and
            with any C{xs} accumulated.
        '''
        kwds = _xkwds_not(None, Fsum._RESIDUAL, f2product =self.f2product(),
                                                nonfinites=self.nonfinites(),
                                                RESIDUAL  =self.RESIDUAL())
        if name_f2product_nonfinites_RESIDUAL:  # overwrites
            kwds.update(name_f2product_nonfinites_RESIDUAL)
        f = Fsum(**kwds)
        # assert all(v == self.__dict__[n] for n, v in f.__dict__.items())
        return f._fset(xs[0], op=_fset_op_) if len(xs) == 1 else (
               f._facc(xs, up=False) if xs else f)

    def fsum2(self, xs=(), **name):
        '''Add an iterable's items, summate and return the
           current precision running sum I{and} the C{residual}.

           @arg xs: Iterable of items to add (each item C{scalar},
                    an L{Fsum} or L{Fsum2Tuple}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: L{Fsum2Tuple}C{(fsum, residual)} with C{fsum} the
                    current precision running sum and C{residual}, the
                    (precision) sum of the remaining C{partials}.  The
                    C{residual is INT0} if the C{fsum} is considered
                    to be I{exact}.

           @see: Methods L{Fsum.fint2}, L{Fsum.fsum} and L{Fsum.fsum2_}
        '''
        t = self._facc(xs)._fprs2
        return t.dup(name=name) if name else t

    def fsum2_(self, *xs):
        '''Add any positional items, summate and return the current
           precision running sum and the I{differential}.

           @arg xs: Values to add (each C{scalar}, an L{Fsum} or
                    L{Fsum2Tuple}), all positional.

           @return: 2Tuple C{(fsum, delta)} with the current, precision
                    running C{fsum} like method L{Fsum.fsum} and C{delta},
                    the difference with previous running C{fsum}, C{float}.

           @see: Methods L{Fsum.fsum_} and L{Fsum.fsum}.
        '''
        return self._fsum2(xs, self._facc_args)

    def _fsum2(self, xs, _facc, **facc_kwds):
        '''(INTERNAL) Helper for L{Fsum.fsum2_} and L{Fsum.fsum2f_}.
        '''
        p, q = self._fprs2
        if xs:
            s, r = _facc(xs, **facc_kwds)._fprs2
            if _isfinite(s):  # _fsum(_1primed((s, -p, r, -q))
                d, r = _2sum(s - p, r - q, _isfine=_isOK)
                r, _ = _s_r2(d, r)
            return s, (r if _isfinite(r) else _NONFINITEr)
        else:
            return p, _0_0

    def fsumf_(self, *xs):
        '''Like method L{Fsum.fsum_} iff I{all} C{B{xs}}, each I{known to be}
           C{scalar}, an L{Fsum} or L{Fsum2Tuple}.
        '''
        return self._facc_scalarf(xs, which=self.fsumf_)._fprs  # origin=1?

    def Fsumf_(self, *xs):
        '''Like method L{Fsum.Fsum_} iff I{all} C{B{xs}}, each I{known to be}
           C{scalar}, an L{Fsum} or L{Fsum2Tuple}.
        '''
        return self._facc_scalarf(xs, which=self.Fsumf_)._copy_2(self.Fsumf_)  # origin=1?

    def fsum2f_(self, *xs):
        '''Like method L{Fsum.fsum2_} iff I{all} C{B{xs}}, each I{known to be}
           C{scalar}, an L{Fsum} or L{Fsum2Tuple}.
        '''
        return self._fsum2(xs, self._facc_scalarf, which=self.fsum2f_)  # origin=1?

#   ftruediv = __itruediv__   # for naming consistency?

    def _ftruediv(self, other, op, **raiser_RESIDUAL):
        '''(INTERNAL) Apply C{B{self} /= B{other}}.
        '''
        n = _1_0
        if _isFsum_2Tuple(other):
            if other is self or self == other:
                return self._fset(n, n=len(self))
            d, r = other._fprs2
            if r:
                R = self._raiser(r, d, **raiser_RESIDUAL)
                if R:
                    raise self._ResidualError(op, other, r, **R)
                d, n = other.as_integer_ratio()
        else:
            d = self._scalar(other, op)
        try:
            s = n / d
        except Exception as X:
            raise self._ErrorX(X, op, other)
        f = self._mul_scalar(s, _mul_op_)  # handles 0, INF, NAN
        return self._fset(f)

    @property_RO
    def imag(self):
        '''Get the C{imaginary} part of this instance (C{0.0}, always).

           @see: Property L{Fsum.real}.
        '''
        return _0_0

    def int_float(self, **raiser_RESIDUAL):
        '''Return this instance' current running sum as C{int} or C{float}.

           @kwarg raiser_RESIDUAL: Use C{B{raiser}=False} to ignore
                         L{ResidualError}s (C{bool}) and C{B{RESIDUAL}=scalar}
                         to override the current L{RESIDUAL<Fsum.RESIDUAL>}.

           @return: This C{int} sum if this instance C{is_integer} and
                    I{finite}, otherwise the C{float} sum if the residual
                    is zero or not significant.

           @raise ResidualError: Non-zero, significant residual or invalid
                                 B{C{RESIDUAL}}.

           @see: Methods L{Fsum.fint}, L{Fsum.fint2}, L{Fsum.is_integer},
                 L{Fsum.RESIDUAL} and property L{Fsum.as_iscalar}.
        '''
        s, r = self._fint2
        if r:
            s, r = self._fprs2
            if r:  # PYCHOK no cover
                R = self._raiser(r, s, **raiser_RESIDUAL)
                if R:
                    t = _stresidual(_non_zero_, r, **R)
                    raise ResidualError(int_float=s, txt=t)
            s = float(s)
        return s

    def is_exact(self):
        '''Is this instance' running C{fsum} considered to be exact?
           (C{bool}), C{True} only if the C{residual is }L{INT0}.
        '''
        return self.residual is INT0

    def is_finite(self):  # in .constants
        '''Is this instance C{finite}? (C{bool}).

           @see: Function L{isfinite<pygeodesy.isfinite>}.
        '''
        return _isfinite(sum(self._ps))  # == sum(self)

    def is_integer(self):
        '''Is this instance' running sum C{integer}? (C{bool}).

           @see: Methods L{Fsum.fint}, L{Fsum.fint2} and L{Fsum.is_scalar}.
        '''
        s, r = self._fint2
        return False if r else (_isfinite(s) and isint(s))

    def is_math_fma(self):
        '''Is accurate L{f2product} multiplication based on Python's C{math.fma}?

           @return: C{True} if accurate multiplication uses C{math.fma}, C{False}
                    an C{fma} implementation as C{math.fma} or C{None}, a previous
                    C{PyGeodesy} implementation.
        '''
        return (_2split3s is _passarg) or (False if _2n_d is None else None)

    def is_math_fsum(self):
        '''Are the summation functions L{fsum}, L{fsum_}, L{fsumf_}, L{fsum1},
           L{fsum1_} and L{fsum1f_} based on Python's C{math.fsum}?

           @return: C{True} if summation functions use C{math.fsum}, C{False}
                    otherwise.
        '''
        return _sum is _fsum  # _fsum.__module__ is fabs.__module__

    def is_scalar(self, **raiser_RESIDUAL):
        '''Is this instance' running sum C{scalar} with C{0} residual or with
           a residual I{ratio} not exceeding the RESIDUAL threshold?

           @kwarg raiser_RESIDUAL: Use C{B{raiser}=False} to ignore
                         L{ResidualError}s (C{bool}) and C{B{RESIDUAL}=scalar}
                         to override the current L{RESIDUAL<Fsum.RESIDUAL>}.

           @return: C{True} if this instance' residual is C{0} or C{insignificant},
                    i.e. its residual C{ratio} doesn't exceed the L{RESIDUAL
                    <Fsum.RESIDUAL>} threshold (C{bool}).

           @raise ResidualError: Non-zero, significant residual or invalid
                                 B{C{RESIDUAL}}.

           @see: Methods L{Fsum.RESIDUAL} and L{Fsum.is_integer} and property
                 L{Fsum.as_iscalar}.
        '''
        s, r = self._fprs2
        return False if r and self._raiser(r, s, **raiser_RESIDUAL) else True

    def _mul_Fsum(self, other, op):
        '''(INTERNAL) Return C{B{self} * B{other}} as L{Fsum} or C{0}.
        '''
        # assert _isFsum_2Tuple(other)
        if self._ps and other._ps:
            try:
                f = self._ps_mul(op, *other._ps)  # NO .as_iscalar!
            except Exception as X:
                raise self._ErrorX(X, op, other)
        else:
            f = _0_0
        return f

    def _mul_reduce(self, *others):
        '''(INTERNAL) Like fmath.fprod for I{non-finite} C{other}s.
        '''
        r = _1_0
        for f in others:
            r *= sum(f._ps) if _isFsum_2Tuple(f) else float(f)
        return r

    def _mul_scalar(self, factor, op):
        '''(INTERNAL) Return C{B{self} * scalar B{factor}} as L{Fsum}, C{0.0} or C{self}.
        '''
        # assert isscalar(factor)
        if self._ps and self._finite(factor, op=op):
            f =  self      if factor == _1_0   else (
                 self._neg if factor == _N_1_0 else
                 self._ps_mul(op, factor).as_iscalar)
        else:
            f = _0_0
        return f

#   @property_RO
#   def _n_d(self):
#       n, d = self.as_integer_ratio()
#       return n / d

    @property_RO
    def _neg(self):
        '''(INTERNAL) Return C{Fsum(-self)} or scalar C{NEG0}.
        '''
        return _Psum(self._ps_neg) if self._ps else NEG0

    @property_RO
    def _nfprs2(self):
        '''(INTERNAL) Handle I{non-finite} C{_fprs2}.
        '''
        try:  # to handle nonfiniterrors, etc.
            t = self._fprs2
        except (OverflowError, ValueError):
            t = Fsum2Tuple(sum(self._ps), _NONFINITEr)
        return t

    def nonfinites(self, *OK):
        '''Handle I{non-finite} C{float}s as C{inf}, C{INF}, C{NINF}, C{nan}
           and C{NAN} for this L{Fsum} or throw C{OverflowError} respectively
           C{ValueError} exceptions, overriding the L{nonfiniterrors} default.

           @arg OK: If omitted, leave the override unchanged, if C{True},
                    I{non-finites} are C{OK}, if C{False} throw exceptions
                    or if C{None} remove the override (C{bool} or C{None}).

           @return: The previous setting (C{bool} or C{None} if not set).

           @see: Function L{nonfiniterrors<fsums.nonfiniterrors>}.

           @note: Use property L{nonfinitesOK<Fsum.nonfinitesOK>} to determine
                  whether I{non-finites} are C{OK} for this L{Fsum} and by the
                  L{nonfiniterrors} default.
        '''
        _ks = Fsum._nonfinites_isfine_kwds
        if OK:  # delattrof(self, _isfine=None)
            k = _xkwds_pop(self.__dict__, _isfine=None)
            if OK[0] is not None:
                self._isfine = _ks[bool(OK[0])]
            self._update()
        else:  # getattrof(self, _isfine=None)
            k = _xkwds_get(self.__dict__, _isfine=None)
        # dict(map(reversed, _ks.items())).get(k, None)
        # raises a TypeError: unhashable type: 'dict'
        return True  if k is _ks[True]  else (
               False if k is _ks[False] else None)

    _nonfinites_isfine_kwds = {True:  dict(_isfine=_isOK),
                               False: dict(_isfine=_isfinite)}

    @property_RO
    def nonfinitesOK(self):
        '''Are I{non-finites} C{OK} for this L{Fsum} or by default? (C{bool}).
        '''
#       nf = self.nonfinites()
#       if nf is None:
#           nf = not nonfiniterrors()
        return _isOK_or_finite(INF, **self._isfine)

    def _nonfiniteX(self, X, op, f, nonfinites=None, raiser=None):
        '''(INTERNAL) Handle a I{non-finite} exception.
        '''
        if nonfinites is None:
            nonfinites = _isOK_or_finite(f, **self._isfine) if raiser is None else (not raiser)
        if not nonfinites:
            raise self._ErrorX(X, op, f)
        return f

    def _optionals(self, f2product=None, nonfinites=None, **name_RESIDUAL):
        '''(INTERNAL) Re/set options from keyword arguments.
        '''
        if f2product is not None:
            self.f2product(f2product)
        if nonfinites is not None:
            self.nonfinites(nonfinites)
        if name_RESIDUAL:  # MUST be last
            n, kwds = _name2__(**name_RESIDUAL)
            if kwds:
                R =  Fsum._RESIDUAL
                t = _threshold(R, **kwds)
                if t != R:
                    self._RESIDUAL = t
            if n:
                self.name = n  # self.rename(n)

    def _1_Over(self, x, op, **raiser_RESIDUAL):  # vs _1_over
        '''(INTERNAL) Return C{Fsum(1) / B{x}}.
        '''
        return self._Fsum_as(_1_0)._ftruediv(x, op, **raiser_RESIDUAL)

    @property_RO
    def partials(self):
        '''Get this instance' current, partial sums (C{tuple} of C{float}s).
        '''
        return tuple(self._ps)

    def pow(self, x, *mod, **raiser_RESIDUAL):
        '''Return C{B{self}**B{x}} as L{Fsum}.

           @arg x: The exponent (C{scalar}, an L{Fsum} or L{Fsum2Tuple}).
           @arg mod: Optional modulus (C{int} or C{None}) for the 3-argument
                     C{pow(B{self}, B{other}, B{mod})} version.
           @kwarg raiser_RESIDUAL: Use C{B{raiser}=False} to ignore
                         L{ResidualError}s (C{bool}) and C{B{RESIDUAL}=scalar}
                         to override the current L{RESIDUAL<Fsum.RESIDUAL>}.

           @return: The C{pow(self, B{x})} or C{pow(self, B{x}, *B{mod})}
                    result (L{Fsum}).

           @raise ResidualError: Non-zero, significant residual or invalid
                                 B{C{RESIDUAL}}.

           @note: If B{C{mod}} is given and C{None}, the result will be an
                  C{integer} L{Fsum} provided this instance C{is_integer}
                  or set to C{integer} by an L{Fsum.fint} call.

           @see: Methods L{Fsum.__ipow__}, L{Fsum.fint}, L{Fsum.is_integer}
                 and L{Fsum.root}.
        '''
        f = self._copy_2(self.pow)
        return f._fpow(x, _pow_op_, *mod, **raiser_RESIDUAL)  # f = pow(f, x, *mod)

    def _pow(self, other, unused, op, **raiser_RESIDUAL):
        '''Return C{B{self} ** B{other}}.
        '''
        if _isFsum_2Tuple(other):
            f = self._pow_Fsum(other, op, **raiser_RESIDUAL)
        elif self._scalar(other, op):
            x = self._finite(other, op=op)
            f = self._pow_scalar(x, other, op, **raiser_RESIDUAL)
        else:
            f = self._pow_0_1(0, other)
        return f

    def _pow_0_1(self, x, other):
        '''(INTERNAL) Return B{C{self}**1} or C{B{self}**0 == 1.0}.
        '''
        return self if x else (1 if isint(other) and self.is_integer() else _1_0)

    def _pow_2_3(self, b, x, other, op, *mod, **raiser_RESIDUAL):
        '''(INTERNAL) 2-arg C{pow(B{b}, scalar B{x})} and 3-arg C{pow(B{b},
           B{x}, int B{mod} or C{None})}, embellishing errors.
        '''

        if mod:  # b, x, mod all C{int}, unless C{mod} is C{None}
            m = mod[0]
            # assert _isFsum_2Tuple(b)

            def _s(s, r):
                R = self._raiser(r, s, **raiser_RESIDUAL)
                if R:
                    raise self._ResidualError(op, other, r, mod=m, **R)
                return s

            b = _s(*(b._fprs2 if m is None else b._fint2))
            x = _s(*_2tuple2(x))

        try:
            # 0**INF == 0.0, 1**INF == 1.0, -1**2.3 == -(1**2.3)
            s = pow(b, x, *mod)
            if iscomplex(s):
                # neg**frac == complex in Python 3+, but ValueError in 2-
                raise ValueError(_strcomplex(s, b, x, *mod))
            _ = _2finite(s, **self._isfine)  # ignore float
            return s
        except Exception as X:
            raise self._ErrorX(X, op, other, *mod)

    def _pow_Fsum(self, other, op, **raiser_RESIDUAL):
        '''(INTERNAL) Return C{B{self} **= B{other}} for C{_isFsum_2Tuple(other)}.
        '''
        # assert _isFsum_2Tuple(other)
        x, r = other._fprs2
        f = self._pow_scalar(x, other, op, **raiser_RESIDUAL)
        if f and r:
            f *= self._pow_scalar(r, other, op, **raiser_RESIDUAL)
        return f

    def _pow_int(self, x, other, op, **raiser_RESIDUAL):
        '''(INTERNAL) Return C{B{self} **= B{x}} for C{int B{x} >= 0}.
        '''
        # assert isint(x) and x >= 0
        ps = self._ps
        if len(ps) > 1:
            _mul_Fsum = Fsum._mul_Fsum
            if x > 4:
                p = self
                f = self if (x & 1) else self._Fsum_as(_1_0)
                m = x >> 1  # // 2
                while m:
                    p = _mul_Fsum(p, p, op)  # p **= 2
                    if (m & 1):
                        f = _mul_Fsum(f, p, op)  # f *= p
                    m >>= 1  # //= 2
            elif x > 1:  # self**2, 3, or 4
                f = _mul_Fsum(self, self, op)
                if x > 2:  # self**3 or 4
                    p =  self if x < 4 else f
                    f = _mul_Fsum(f, p, op)
            else:  # self**1 or self**0 == 1 or _1_0
                f = self._pow_0_1(x, other)
        elif ps:  # self._ps[0]**x
            f = self._pow_2_3(ps[0], x, other, op, **raiser_RESIDUAL)
        else:  # PYCHOK no cover
            # 0**pos_int == 0, but 0**0 == 1
            f = 0 if x else 1
        return f

    def _pow_scalar(self, x, other, op, **raiser_RESIDUAL):
        '''(INTERNAL) Return C{self**B{x}} for C{scalar B{x}}.
        '''
        s, r = self._fprs2
        if r:
            # assert s != 0
            if isint(x, both=True):  # self**int
                x = int(x)
                y = abs(x)
                if y > 1:
                    f = self._pow_int(y, other, op, **raiser_RESIDUAL)
                    if x > 0:  # i.e. > 1
                        return f  # Fsum or scalar
                    # assert x < 0  # i.e. < -1
                    if _isFsum(f):
                        s, r = f._fprs2
                        if r:
                            return self._1_Over(f, op, **raiser_RESIDUAL)
                    else:  # scalar
                        s = f
                    # use s**(-1) to get the CPython
                    # float_pow error iff s is zero
                    x = -1
                elif x < 0:  # self**(-1)
                    return self._1_Over(self, op, **raiser_RESIDUAL)  # 1 / self
                else:  # self**1 or self**0
                    return self._pow_0_1(x, other)  # self, 1 or 1.0
            else:  # self**fractional
                R = self._raiser(r, s, **raiser_RESIDUAL)
                if R:
                    raise self._ResidualError(op, other, r, **R)
                n, d = self.as_integer_ratio()
                if abs(n) > abs(d):
                    n, d, x = d, n, (-x)
                s = n / d
        # assert isscalar(s) and isscalar(x)
        return self._pow_2_3(s, x, other, op, **raiser_RESIDUAL)

    def _ps_acc(self, ps, xs, up=True, **unused):
        '''(INTERNAL) Accumulate C{xs} known scalars into list C{ps}.
        '''
        n   =  0
        _2s = _2sum
        _fi =  self._isfine
        for x in (tuple(xs) if xs is ps else xs):
            # assert isscalar(x) and _isOK_or_finite(x, **self._isfine)
            if x:
                i = 0
                for p in ps:
                    x, p = _2s(x, p, **_fi)
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
#               x = sum(ps)
#               if not _isOK_or_finite(x, **fi):
#                   ps[:] = x,  # collapse ps
        return ps

    def _ps_mul(self, op, *factors):
        '''(INTERNAL) Multiply this instance' C{partials} with
           each scalar C{factor} and accumulate into an C{Fsum}.
        '''
        def _psfs(ps, fs, _isfine=_isfinite):
            if len(ps) < len(fs):
                ps, fs = fs, ps
            if self._f2product:
                fs, p = _2split3s(fs), fs
                if len(ps) > 1 and fs is not p:
                    fs = tuple(fs)  # several ps
                _pfs = _2products
            else:
                def _pfs(p, fs):
                    return (p * f for f in fs)

            for p in ps:
                for x in _pfs(p, fs):
                    yield x if _isfine(x) else _nfError(x)

        xs = _psfs(self._ps, factors, **self._isfine)
        f  = _Psum(self._ps_acc([], xs, up=False), name=op)
        return f

    @property_RO
    def _ps_neg(self):
        '''(INTERNAL) Yield the partials, I{negated}.
        '''
        for p in self._ps:
            yield -p

    def _ps_other(self, op, other):
        '''(INTERNAL) Yield C{other} as C{scalar}s.
        '''
        if _isFsum_2Tuple(other):
            for p in other._ps:
                yield p
        else:
            yield self._scalar(other, op)

    def _ps_1sum(self, *less):
        '''(INTERNAL) Return the partials sum, 1-primed C{less} some scalars.
        '''
        def _1psls(ps, ls):
            yield _1_0
            for p in ps:
                yield  p
            for p in ls:
                yield -p
            yield _N_1_0

        return _fsum(_1psls(self._ps, less))

    def _raiser(self, r, s, raiser=True, **RESIDUAL):
        '''(INTERNAL) Does ratio C{r / s} exceed the RESIDUAL threshold
           I{and} is residual C{r} I{non-zero} or I{significant} (for a
           negative respectively positive C{RESIDUAL} threshold)?
        '''
        if r and raiser:
            t = self._RESIDUAL
            if RESIDUAL:
                t = _threshold(t, **RESIDUAL)
            if t < 0 or (s + r) != s:
                q = (r / s) if s else s  # == 0.
                if fabs(q) > fabs(t):
                    return dict(ratio=q, R=t)
        return {}

    rdiv = __rtruediv__

    @property_RO
    def real(self):
        '''Get the C{real} part of this instance (C{float}).

           @see: Methods L{Fsum.__float__} and L{Fsum.fsum}
                 and properties L{Fsum.ceil}, L{Fsum.floor},
                 L{Fsum.imag} and L{Fsum.residual}.
        '''
        return float(self)

    @property_RO
    def residual(self):
        '''Get this instance' residual or residue (C{float} or C{int}):
           the C{sum(partials)} less the precision running sum C{fsum}.

           @note: The C{residual is INT0} iff the precision running
                  C{fsum} is considered to be I{exact}.

           @see: Methods L{Fsum.fsum}, L{Fsum.fsum2} and L{Fsum.is_exact}.
        '''
        return self._fprs2.residual

    def RESIDUAL(self, *threshold):
        '''Get and set this instance' I{ratio} for raising L{ResidualError}s,
           overriding the default from env variable C{PYGEODESY_FSUM_RESIDUAL}.

           @arg threshold: If C{scalar}, the I{ratio} to exceed for raising
                           L{ResidualError}s in division and exponention, if
                           C{None}, restore the default set with env variable
                           C{PYGEODESY_FSUM_RESIDUAL} or if omitted, keep the
                           current setting.

           @return: The previous C{RESIDUAL} setting (C{float}), default C{0.0}.

           @raise ResidualError: Invalid B{C{threshold}}.

           @note: L{ResidualError}s may be thrown if (1) the non-zero I{ratio}
                  C{residual / fsum} exceeds the given B{C{threshold}} and (2)
                  the C{residual} is non-zero and (3) is I{significant} vs the
                  C{fsum}, i.e. C{(fsum + residual) != fsum} and (4) optional
                  keyword argument C{raiser=False} is missing.  Specify a
                  negative B{C{threshold}} for only non-zero C{residual}
                  testing without the I{significant} case.
        '''
        r = self._RESIDUAL
        if threshold:
            t = threshold[0]
            self._RESIDUAL = Fsum._RESIDUAL if t is None else (  # for ...
                           (_0_0 if t else _1_0) if isbool(t) else
                           _threshold(t))  # ... backward compatibility
        return r

    def _ResidualError(self, op, other, residual, **mod_R):
        '''(INTERNAL) Non-zero B{C{residual}} etc.
        '''
        def _p(mod=None, R=0, **unused):  # ratio=0
            return (_non_zero_ if R < 0 else _significant_) \
                               if mod is None else _integer_

        t = _stresidual(_p(**mod_R), residual, **mod_R)
        return self._Error(op, other, ResidualError, txt=t)

    def root(self, root, **raiser_RESIDUAL):
        '''Return C{B{self}**(1 / B{root})} as L{Fsum}.

           @arg root: Non-zero order (C{scalar}, an L{Fsum} or L{Fsum2Tuple}).
           @kwarg raiser_RESIDUAL: Use C{B{raiser}=False} to ignore any
                         L{ResidualError}s (C{bool}) or C{B{RESIDUAL}=scalar}
                         to override the current L{RESIDUAL<Fsum.RESIDUAL>}.

           @return: The C{self ** (1 / B{root})} result (L{Fsum}).

           @raise ResidualError: Non-zero, significant residual or invalid
                                 B{C{RESIDUAL}}.

           @see: Method L{Fsum.pow}.
        '''
        x = self._1_Over(root, _truediv_op_, **raiser_RESIDUAL)
        f = self._copy_2(self.root)
        return f._fpow(x, f.name, **raiser_RESIDUAL)  # == pow(f, x)

    def _scalar(self, other, op, **txt):
        '''(INTERNAL) Return scalar C{other} or throw a C{TypeError}.
        '''
        if isscalar(other):
            return other
        raise self._Error(op, other, _TypeError, **txt)  # _invalid_

    def signOf(self, res=True):
        '''Determine the sign of this instance.

           @kwarg res: If C{True}, consider the residual,
                       otherwise ignore the latter (C{bool}).

           @return: The sign (C{int}, -1, 0 or +1).
        '''
        s, r = self._nfprs2
        r = (-r) if res else 0
        return _signOf(s, r)

    def toRepr(self, **lenc_prec_sep_fmt):  # PYCHOK signature
        '''Return this C{Fsum} instance as representation.

           @kwarg lenc_prec_sep_fmt: Optional keyword arguments
                       for method L{Fsum.toStr}.

           @return: This instance (C{repr}).
        '''
        return Fmt.repr_at(self, self.toStr(**lenc_prec_sep_fmt))

    def toStr(self, lenc=True, **prec_sep_fmt):  # PYCHOK signature
        '''Return this C{Fsum} instance as string.

           @kwarg lenc: If C{True}, include the current C{[len]} of this
                        L{Fsum} enclosed in I{[brackets]} (C{bool}).
           @kwarg prec_sep_fmt: Optional keyword arguments for method
                       L{Fsum2Tuple.toStr}.

           @return: This instance (C{str}).
        '''
        p = self.classname
        if lenc:
            p = Fmt.SQUARE(p, len(self))
        n = _enquote(self.name, white=_UNDER_)
        t =  self._nfprs2.toStr(**prec_sep_fmt)
        return NN(p, _SPACE_, n, t)

    def _truediv(self, other, op, **raiser_RESIDUAL):
        '''(INTERNAL) Return C{B{self} / B{other}} as an L{Fsum}.
        '''
        f = self._copy_2(self.__truediv__)
        return f._ftruediv(other, op, **raiser_RESIDUAL)

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

_ROs = _allPropertiesOf_n(3, Fsum, Property_RO)  # PYCHOK see Fsum._update

if _NONFINITES == _std_:  # PYCHOK no cover
    _ = nonfiniterrors(False)


def _Float_Int(arg, **name_Error):
    '''(INTERNAL) L{DivMod2Tuple}, L{Fsum2Tuple} Unit.
    '''
    U = Int if isint(arg) else Float
    return U(arg, **name_Error)


class DivMod2Tuple(_NamedTuple):
    '''2-Tuple C{(div, mod)} with the quotient C{div} and remainder
       C{mod} results of a C{divmod} operation.

       @note: Quotient C{div} an C{int} in Python 3+ but a C{float}
              in Python 2-.  Remainder C{mod} an L{Fsum} instance.
    '''
    _Names_ = ('div',     'mod')
    _Units_ = (_Float_Int, Fsum)


class Fsum2Tuple(_NamedTuple):  # in .fstats
    '''2-Tuple C{(fsum, residual)} with the precision running C{fsum}
       and the C{residual}, the sum of the remaining partials.  Each
       item is C{float} or C{int}.

       @note: If the C{residual is INT0}, the C{fsum} is considered
              to be I{exact}, see method L{Fsum2Tuple.is_exact}.
    '''
    _Names_ = ( Fsum.fsum.__name__, Fsum.residual.name)
    _Units_ = (_Float_Int,         _Float_Int)

    def __abs__(self):  # in .fmath
        return self._Fsum.__abs__()

    def __bool__(self):  # PYCHOK Python 3+
        return bool(self._Fsum)

    def __eq__(self, other):
        return self._other_op(other, self.__eq__)

    def __float__(self):
        return self._Fsum.__float__()

    def __ge__(self, other):
        return self._other_op(other, self.__ge__)

    def __gt__(self, other):
        return self._other_op(other, self.__gt__)

    def __le__(self, other):
        return self._other_op(other, self.__le__)

    def __lt__(self, other):
        return self._other_op(other, self.__lt__)

    def __int__(self):
        return self._Fsum.__int__()

    def __ne__(self, other):
        return self._other_op(other, self.__ne__)

    def __neg__(self):
        return self._Fsum.__neg__()

    __nonzero__ = __bool__  # Python 2-

    def __pos__(self):
        return self._Fsum.__pos__()

    def as_integer_ratio(self):
        '''Return this instance as the ratio of 2 integers.

           @see: Method L{Fsum.as_integer_ratio} for further details.
        '''
        return self._Fsum.as_integer_ratio()

    @property_RO
    def _fint2(self):
        return self._Fsum._fint2

    @property_RO
    def _fprs2(self):
        return self._Fsum._fprs2

    @Property_RO
    def _Fsum(self):  # this C{Fsum2Tuple} as L{Fsum}, in .fstats
        s, r = _s_r2(*self)
        ps = (r, s) if r else (s,)
        return _Psum(ps, name=self.name)

    def Fsum_(self, *xs, **name_f2product_nonfinites_RESIDUAL):
        '''Return this C{Fsum2Tuple} as an L{Fsum} plus some C{xs}.
        '''
        return Fsum(self, *xs, **name_f2product_nonfinites_RESIDUAL)

    def is_exact(self):
        '''Is this L{Fsum2Tuple} considered to be exact? (C{bool}).
        '''
        return self._Fsum.is_exact()

    def is_finite(self):  # in .constants
        '''Is this L{Fsum2Tuple} C{finite}? (C{bool}).

           @see: Function L{isfinite<pygeodesy.isfinite>}.
        '''
        return self._Fsum.is_finite()

    def is_integer(self):
        '''Is this L{Fsum2Tuple} C{integer}? (C{bool}).
        '''
        return self._Fsum.is_integer()

    def _mul_scalar(self, other, op):  # for Fsum._fmul
        return self._Fsum._mul_scalar(other, op)

    @property_RO
    def _n(self):
        return self._Fsum._n

    def _other_op(self, other, which):
        C, s = (tuple, self) if isinstance(other, tuple) else (Fsum, self._Fsum)
        return getattr(C, which.__name__)(s, other)

    @property_RO
    def _ps(self):
        return self._Fsum._ps

    @property_RO
    def _ps_neg(self):
        return self._Fsum._ps_neg

    def signOf(self, **res):
        '''Like method L{Fsum.signOf}.
        '''
        return self._Fsum.signOf(**res)

    def toStr(self, fmt=Fmt.g, **prec_sep):  # PYCHOK signature
        '''Return this L{Fsum2Tuple} as string (C{str}).

           @kwarg fmt: Optional C{float} format (C{letter}).
           @kwarg prec_sep: Optional keyword arguments for function
                       L{fstr<streprs.fstr>}.
        '''
        return Fmt.PAREN(fstr(self, fmt=fmt, strepr=str, force=False, **prec_sep))

_Fsum_2Tuple_types = Fsum, Fsum2Tuple  # PYCHOK lines


class ResidualError(_ValueError):
    '''Error raised for a division, power or root operation of
       an L{Fsum} instance with a C{residual} I{ratio} exceeding
       the L{RESIDUAL<Fsum.RESIDUAL>} threshold.

       @see: Module L{pygeodesy.fsums} and method L{Fsum.RESIDUAL}.
    '''
    pass


try:
    from math import fsum as _fsum  # precision IEEE-754 sum, Python 2.6+

    # make sure _fsum works as expected (XXX check
    # float.__getformat__('float')[:4] == 'IEEE'?)
    if _fsum((1, 1e101, 1, -1e101)) != 2:  # PYCHOK no cover
        del _fsum  # nope, remove _fsum ...
        raise ImportError()  # ... use _fsum below

    _sum = _fsum  # in .elliptic
except ImportError:
    _sum =  sum  # in .elliptic

    def _fsum(xs):
        '''(INTERNAL) Precision summation, Python 2.5-.
        '''
        F = Fsum(name=_fsum.name, f2product=False, nonfinites=True)
        return float(F._facc(xs, up=False))


def fsum(xs, nonfinites=None, **floats):
    '''Precision floating point summation from Python's C{math.fsum}.

       @arg xs: Iterable of items to add (each C{scalar}, an L{Fsum} or L{Fsum2Tuple}).
       @kwarg nonfinites: Use C{B{nonfinites}=True} if I{non-finites} are C{OK}, if
                          C{False} I{non-finites} raise an Overflow-/ValueError or if
                          C{None}, L{nonfiniterrors} applies (C{bool} or C{None}).
       @kwarg floats: DEPRECATED keyword argument C{B{floats}=False} (C{bool}), use
                      keyword argument C{B{nonfinites}=False} instead.

       @return: Precision C{fsum} (C{float}).

       @raise OverflowError: Infinite B{C{xs}} item or intermediate C{math.fsum} overflow.

       @raise TypeError: Invalid B{C{xs}} item.

       @raise ValueError: Invalid or C{NAN} B{C{xs}} item.

       @see: Function L{nonfiniterrors}, class L{Fsum} and methods L{Fsum.nonfinites},
             L{Fsum.fsum}, L{Fsum.fadd} and L{Fsum.fadd_}.
    '''
    return _xsum(fsum, xs, nonfinites=nonfinites, **floats) if xs else _0_0


def fsum_(*xs, **nonfinites):
    '''Precision floating point summation of all positional items.

       @arg xs: Items to add (each C{scalar}, an L{Fsum} or L{Fsum2Tuple}), all positional.
       @kwarg nonfinites: Use C{B{nonfinites}=True} if I{non-finites} are C{OK} (C{bool}).

       @see: Function L{fsum<fsums.fsum>} for further details.
    '''
    return _xsum(fsum_, xs, **nonfinites) if xs else _0_0  # origin=1?


def fsumf_(*xs):
    '''Precision floating point summation of all positional items with I{non-finites} C{OK}.

       @arg xs: Items to add (each C{scalar}, an L{Fsum} or L{Fsum2Tuple}),
                all positional.

       @see: Function L{fsum_<fsums.fsum_>} for further details.
    '''
    return _xsum(fsumf_, xs, nonfinites=True) if xs else _0_0  # origin=1?


def fsum1(xs, **nonfinites):
    '''Precision floating point summation, 1-primed.

       @arg xs: Iterable of items to add (each C{scalar}, an L{Fsum} or L{Fsum2Tuple}).
       @kwarg nonfinites: Use C{B{nonfinites}=True} if I{non-finites} are C{OK} (C{bool}).

       @see: Function L{fsum<fsums.fsum>} for further details.
    '''
    return _xsum(fsum1, xs, primed=1, **nonfinites) if xs else _0_0


def fsum1_(*xs, **nonfinites):
    '''Precision floating point summation of all positional items, 1-primed.

       @arg xs: Items to add (each C{scalar}, an L{Fsum} or L{Fsum2Tuple}), all positional.
       @kwarg nonfinites: Use C{B{nonfinites}=True} if I{non-finites} are C{OK} (C{bool}).

       @see: Function L{fsum_<fsums.fsum_>} for further details.
    '''
    return _xsum(fsum1_, xs, primed=1, **nonfinites) if xs else _0_0  # origin=1?


def fsum1f_(*xs):
    '''Precision floating point summation of all positional items, 1-primed and
       with I{non-finites} C{OK}.

       @see: Function L{fsum_<fsums.fsum_>} for further details.
    '''
    return _xsum(fsum1f_, xs, nonfinites=True, primed=1) if xs else _0_0


def _x_isfine(nfOK, **kwds):  # get the C{_x} and C{_isfine} handlers.
    _x_kwds = dict(_x=     (_passarg if nfOK else _2finite),
                   _isfine=(_isOK    if nfOK else _isfinite))  # PYCHOK kwds
    _x_kwds.update(kwds)
    return _x_kwds


def _X_ps(X):  # default C{_X} handler
    return X._ps  # lambda X: X._ps


def _xs(xs, _X=_X_ps, _x=float, _isfine=_isfinite,  # defaults for Fsum._facc
                        origin=0, which=None, **_Cdot):
    '''(INTERNAL) Yield each C{xs} item as 1 or more C{float}s.
    '''
    i, x = 0, xs
    try:
        for i, x in enumerate(_xiterable(xs)):
            if _isFsum_2Tuple(x):
                for p in _X(x):
                    yield p if _isfine(p) else _nfError(p)
            else:
                f = _x(x)
                yield f if _isfine(f) else _nfError(f)

    except (OverflowError, TypeError, ValueError) as X:
        t = _xsError(X, xs, i + origin, x)
        if which:  # prefix invokation
            w =  unstr(which, *xs, _ELLIPSIS=4, **_Cdot)
            t = _COMMASPACE_(w, t)
        raise _xError(X, t, txt=None)


def _xsum(which, xs, nonfinites=None, primed=0, **floats):  # origin=0
    '''(INTERNAL) Precision summation of C{xs} with conditions.
    '''
    if floats:  # for backward compatibility
        nonfinites = _xkwds_get1(floats, floats=nonfinites)
    elif nonfinites is None:
        nonfinites =  not nonfiniterrors()
    fs = _xs(xs, **_x_isfine(nonfinites, which=which))  # PYCHOK yield
    return _fsum(_1primed(fs) if primed else fs)


# delete all decorators, etc.
del _allPropertiesOf_n, deprecated_method, deprecated_property_RO, \
     Property, Property_RO, property_RO, _ALL_LAZY, _F2PRODUCT, \
     MANT_DIG, _NONFINITES, _RESIDUAL_0_0, _getPYGEODESY, _std_

if __name__ == '__main__':

    # usage: python3 -m pygeodesy.fsums

    def _test(n):
        # copied from Hettinger, see L{Fsum} reference
        from pygeodesy import frandoms, printf

        printf(_fsum.__name__, end=_COMMASPACE_)
        printf(_psum.__name__, end=_COMMASPACE_)

        F = Fsum()
        if F.is_math_fsum():
            for t in frandoms(n, seeded=True):
                assert float(F.fset_(*t)) == _fsum(t)
                printf(_DOT_, end=NN)
        printf(NN)

    _test(128)

# **) MIT License
#
# Copyright (C) 2016-2025 -- mrJean1 at Gmail -- All Rights Reserved.
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
