# -*- coding: utf-8 -*-

u'''Single-instance C{float} and C{int} constants across C{pygeodesy}
modules and related functions L{pygeodesy.float_}, L{pygeodesy.isclose},
L{pygeodesy.isfinite}, L{pygeodesy.isinf}, L{pygeodesy.isint0},
L{pygeodesy.isnan}, L{pygeodesy.isnear0}, L{pygeodesy.isnear1},
L{pygeodesy.isnear90}, L{pygeodesy.isneg0}, L{pygeodesy.isninf},
L{pygeodesy.isnon0} and L{pygeodesy.remainder}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import _0_0, _copysign, isbool, iscomplex, isint
from pygeodesy.errors import _xError, _xError2, _xkwds_get
# from pygeodesy.fmath import fprod  # from _MODS
from pygeodesy.interns import _INF_, _NAN_, _UNDER_
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS
# from pygeodesy.streprs import Fmt  # from .unitsBase
from pygeodesy.unitsBase import Float, Int, Radius,  Fmt

from math import fabs, isinf, isnan, pi as _PI, sqrt
try:
    from math import inf as _inf, nan as _nan  # PYCHOK Python 3+
except ImportError:  # Python 2-
    _inf, _nan = float(_INF_), float(_NAN_)
try:
    from math import log2 as _log2
except ImportError:  # Python 3.3-
    from math import log as _log

    def _log2(x):  # in .rhumbaux, .auxilats.auxLat
        return _log(x, 2)

__all__ = _ALL_LAZY.constants
__version__ = '23.08.30'


def _copysign_0_0(y):
    '''(INTERNAL) copysign(0.0, y), only C{float}.
    '''
    return _N_0_0 if y < 0 else _0_0


def _copysign_1_0(y):
    '''(INTERNAL) copysign(1.0, y), only C{float}.
    '''
    return _N_1_0 if y < 0 else _1_0


def _copysignINF(y):
    '''(INTERNAL) copysign(INF, y), only C{float}.
    '''
    return NINF if y < 0 else INF


def _Float(**name_arg):
    '''(INTERNAL) New named, cached C{Float}.
    '''
    n, arg = name_arg.popitem()
    return Float(_float(arg), name=n)


def _Radius(**name_arg):
    '''(INTERNAL) New named, cached C{Radius}.
    '''
    n, arg = name_arg.popitem()
    return Radius(_float(arg), name=n)


def float_(*fs, **sets):  # sets=False
    '''Get scalars as C{float} or I{intern}'ed C{float}.

       @arg fs: One more values (C{scalar}), all positional.
       @kwarg sets: Use C{B{sets}=True} to C{intern} each
                    B{C{fs}}, otherwise don't C{intern}.

       @return: A single C{float} if only one B{C{fs}} is
                given, otherwise a tuple of C{float}s.

       @raise TypeError: Some B{C{fs}} is not C{scalar}.
    '''
    fl = []
    _a = fl.append
    _f = _floats.setdefault if _xkwds_get(sets, sets=False) else \
         _floats.get
    try:
        for i, f in enumerate(fs):
            f = float(f)
            _a(_f(f, f))
    except Exception as x:
        E, t = _xError2(x)
        fs_i =  Fmt.SQUARE(fs=i)
        raise E(fs_i, f, txt=t)
    return fl[0] if len(fl) == 1 else tuple(fl)


def _float(f):  # in .datums, .ellipsoids, ...
    '''(INTERNAL) Cache initial C{float}s.
    '''
    f = float(f)
    return _floats.setdefault(f, f)  # PYCHOK del _floats


def float0_(*xs):
    '''Yield C{B{x}s} as a non-NEG0 C{float}.
    '''
    for x in xs:
        yield float(x) if x else _0_0


def _float0(f):  # in .resections, .vector3dBase, ...
    '''(INTERNAL) Return C{float(B{f})} or C{INT0}.
    '''
    if f:
        f =  float(f)
        f = _floats.get(f, f)
    elif f is not INT0:
        f = _0_0
    return f


def _floatuple(*fs):
    '''(INTERNAL) Cache a tuple of C{float}s.
    '''
    return tuple(map(_float, fs))


def _naninf(*xs):
    '''(INTERNAL) Return C{NAN}, C{NINF} or C{INF}.
    '''
    return NAN if NAN in xs else _copysignINF(_MODS.fmath.fprod(xs))


def _over(p, q):
    '''(INTERNAL) Return C{B{p} / B{q}} avoiding ZeroDivisionError exceptions.
    '''
    try:
        return  p / q
    except ZeroDivisionError:
        return _naninf(p)


def _1_over(x):
    '''(INTERNAL) Return reciprocal C{1 / B{x}} avoiding ZeroDivisionError exceptions.
    '''
    try:
        return _1_0 / float(x)
    except ZeroDivisionError:
        return  INF


_floats  = {}     # PYCHOK floats cache, in .__main__
# _float = float  # PYCHOK expected
# del _floats     # XXX zap floats cache never

_0_0     = _float(  _0_0)     # PYCHOK expected
_0_0_1T  = _0_0,              # PYCHOK 1-tuple
_0_0_9T  = _0_0_1T * 9        # PYCHOK 9-tuple
_0_0001  = _float(   0.0001)  # PYCHOK expected
_0_001   = _float(   0.001)   # PYCHOK expected
_0_01    = _float(   0.01)    # PYCHOK expected
_0_1     = _float(   0.1)     # PYCHOK expected
_0_125   = _float(   0.125)   # PYCHOK expected
_0_25    = _float(   0.25)    # PYCHOK expected
_0_26    = _float(   0.26)    # PYCHOK expected
_0_5     = _float(   0.5)     # PYCHOK expected
_1_0     = _float(   1)       # PYCHOK expected
_1_0_1T  = _1_0,              # PYCHOK 1-tuple
_1_5     = _float(   1.5)     # PYCHOK expected
_1_75    = _float(   1.75)    # PYCHOK expected
_2_0     = _float(   2)       # PYCHOK expected
_3_0     = _float(   3)       # PYCHOK expected
_4_0     = _float(   4)       # PYCHOK expected
_5_0     = _float(   5)       # PYCHOK expected
_6_0     = _float(   6)       # PYCHOK expected
_8_0     = _float(   8)       # PYCHOK expected
_9_0     = _float(   9)       # PYCHOK expected
_10_0    = _float(  10)       # PYCHOK expected
_16_0    = _float(  16)       # PYCHOK expected
_32_0    = _float(  32)       # PYCHOK expected
_60_0    = _float(  60)       # PYCHOK expected
_90_0    = _float(  90)       # PYCHOK expected
_100_0   = _float( 100)       # PYCHOK expected
_180_0   = _float( 180)       # PYCHOK expected
_270_0   = _float( 270)       # PYCHOK expected
_360_0   = _float( 360)       # PYCHOK expected
_400_0   = _float( 400)       # PYCHOK expected
_720_0   = _float( 720)       # PYCHOK expected
_1000_0  = _float(1000)       # PYCHOK expected
_3600_0  = _float(3600)       # PYCHOK expected

_N_0_0   =  float(  '-0.0')  # PYCHOK NOT _float!
_N_0_5   = _float(  -_0_5)   # PYCHOK expected
_N_1_0   = _float(  -_1_0)   # PYCHOK expected
_N_2_0   = _float(  -_2_0)   # PYCHOK expected
_N_90_0  = _float( -_90_0)   # PYCHOK expected
_N_180_0 = _float(-_180_0)   # PYCHOK expected

_M_KM =       _1000_0     # meter per Kilo meter, see .utily
_M_NM = _float(1852.0)    # meter per Nautical Mile
_M_SM = _float(1609.344)  # meter per Statute Mile

try:
    from sys import float_info as _f_i
    # @see: <https://NumPy.org/doc/stable/reference/generated/numpy.finfo.html>
    DIG      =  Int(  DIG     =_f_i.dig)       # PYCHOK system's float decimal digits
    EPS      = _Float(EPS     =_f_i.epsilon)   # PYCHOK system's EPSilon
    MANT_DIG =  Int(  MANT_DIG=_f_i.mant_dig)  # PYCHOK system's float mantissa bits
    MAX      = _Float(MAX     =_f_i.max)       # PYCHOK system's MAX float 1.7976931348623157e+308
    MAX_EXP  =  Int(  MAX_EXP =_f_i.max_exp)   # PYTHON system's max base 2 exponent
    MIN      = _Float(MIN     =_f_i.min)       # PYCHOK system's MIN float 2.2250738585072014e-308
    MIN_EXP  =  Int(  MIN_EXP =_f_i.min_exp)   # PYTHON system's min base 2 exponent
#   RADIX    =  Int(  RADIX   =_f_i.radix)     # PYTHON system's float base
    del _f_i
except ImportError:  # PYCHOK no cover
    DIG      =  Int(  DIG     =15)          # PYCHOK system's 64-bit float decimal digits
    EPS      = _Float(EPS     =2.220446049250313e-16)  # PYCHOK EPSilon 2**-52, M{EPS +/- 1 != 1}
    MANT_DIG =  Int(  MANT_DIG=53)          # PYCHOK float mantissa bits ≈ 53 (C{int})
    MAX      = _Float(MAX     =pow(_2_0,  1023) * (_2_0 - EPS))  # PYCHOK ≈ 10**308
    MAX_EXP  =  Int(  MAX_ESP =_log2(MAX))  # 308 base 10
    MIN      = _Float(MIN     =pow(_2_0, -1021))  # PYCHOK ≈ 10**-308
    MIN_EXP  =  Int(MIN_EXP   =_log2(MIN))  # -307 base 10
#   RADIX    =  Int(Radix     =2)           # base

EPS0     = _Float( EPS0   = EPS**2)         # PYCHOK near-/non-zero comparison 4.930381e-32, or EPS or EPS_2
EPS02    = _Float( EPS02  = EPS**4)         # PYCHOK near-zero-squared comparison 2.430865e-63
EPS_2    = _Float( EPS_2  = EPS / _2_0)     # PYCHOK ≈ 1.110223024625e-16
EPS1     = _Float( EPS1   =_1_0 - EPS)      # PYCHOK ≈ 0.9999999999999998
EPS2     = _Float( EPS2   = EPS * _2_0)     # PYCHOK ≈ 4.440892098501e-16
EPS4     = _Float( EPS4   = EPS * _4_0)     # PYCHOK ≈ 8.881784197001e-16
# _1EPS  = _Float(_1EPS   =_1_0 + EPS)      # PYCHOK ≈ 1.0000000000000002
_1_EPS   = _Float(_1_EPS  =_1_0 / EPS)      # PYCHOK = 4503599627370496.0
# _2_EPS = _Float(_2_EPS  =_2_0 / EPS)      # PYCHOK = 9007199254740992.0
_EPS2e4  = _Float(_EPS2e4 = EPS2 * 1.e4)    # PYCHOK ≈ 4.440892098501e-12
_EPS4e8  = _Float(_EPS4e8 = EPS4 * 1.e8)    # PYCHOK ≈ 8.881784197001e-08
_EPSmin  = _Float(_EPSmin = sqrt(MIN))      # PYCHOK = 1.49166814624e-154
_EPSqrt  = _Float(_EPSqrt = sqrt(EPS))      # PYCHOK = 1.49011611938e5-08
_EPStol  = _Float(_EPStol =_EPSqrt * _0_1)  # PYCHOK = 1.49011611938e5-09 == sqrt(EPS * _0_01)

_89_999_ = _Float(_89_999_= EPS1 * _90_0)   # just below 90.0
# <https://Numbers.Computation.Free.FR/Constants/Miscellaneous/digits.html>
_1__90   = _Float(_1__90  =_1_0 / _90_0)    # PYCHOK = 0.011_111_111_111_111_111_111_111_111_111_111_111_111_111_111_11111
_2__PI   = _Float(_2__PI  =_2_0 / _PI)      # PYCHOK = 0.636_619_772_367_581_343_075_535_053_490_057_448_137_838_582_96182

_1_16th  = _Float(_1_16th =_1_0 / _16_0)  # PYCHOK in .ellipsoids, .karney
_1_64th  = _Float(_1_64th =_1_0 /  64)    # PYCHOK in .elliptic, pow(2.0, -6)
_1_3rd   = _Float(_1_3rd  =_1_0 /  _3_0)  # PYCHOK in .fmath
_1_6th   = _Float(_1_6th  =_1_0 /  _6_0)  # PYCHOK in .fmath
_2_3rd   = _Float(_2_3rd  =_2_0 /  _3_0)  # PYCHOK in .fmath

_K0_UTM  = _Float(_K0_UTM = 0.9996)  # PYCHOK in .etm, .ktm, .utm, UTM scale at central meridian
# sqrt(2) <https://WikiPedia.org/wiki/Square_root_of_2>
# 1.414213562373095_048_801_688_724_209_698_078_569_671_875_376_948_073_176_679_737_99
# _1SQRT2= _Float(_1SQRT2 =sqrt(_2_0) + 1)
_SQRT2_2 = _Float(_SQRT2_2=sqrt(_0_5))  # PYCHOK = 0.707106781186547_6 == sqrt(2) / 2

INF   =  Float(INF =_inf)    # PYCHOK INFinity, see function L{isinf}, L{isfinite}, NOT _Float!
INT0  =  Int(  INT0= 0)      # PYCHOK unique int(0) instance, see .fsums, useZ=False
NAN   =  Float(NAN =_nan)    # PYCHOK Not-A-Number, see function L{isnan}, NOT _Float!
NEG0  =  Float(NEG0=_N_0_0)  # PYCHOK NEGative 0.0, see function L{isneg0}, NOT _Float!
NINF  =  Float(NINF=-INF)    # PYCHOK Negative INFinity, NOT _Float!

PI    = _Float(PI   =_PI)
PI2   = _Float(PI2  =_PI * _2_0)  # PYCHOK Two PI, M{PI * 2} aka I{Tau}
PI_2  = _Float(PI_2 =_PI / _2_0)  # PYCHOK Half PI, M{PI / 2}
PI3   = _Float(PI3  =_PI * _3_0)  # PYCHOK Three PI, M{PI * 3}
PI3_2 = _Float(PI3_2=_PI * _1_5)  # PYCHOK PI and a half, M{PI * 3 / 2}
PI_3  = _Float(PI_3 =_PI / _3_0)  # PYCHOK One third PI, M{PI / 3}
PI4   = _Float(PI4  =_PI * _4_0)  # PYCHOK Four PI, M{PI * 4}
PI_4  = _Float(PI_4 =_PI / _4_0)  # PYCHOK Quarter PI, M{PI / 4}

R_MA  = _Radius(R_MA=6378137.0)       # PYCHOK equatorial earth radius (C{meter}), WGS84, EPSG:3785
R_MB  = _Radius(R_MB=6356752.3)       # PYCHOK polar earth radius (C{meter}), WGS84, EPSG:3785
R_M   = _Radius(R_M =6371008.771415)  # PYCHOK mean, spherical earth radius (C{meter})
R_KM  = _Radius(R_KM=R_M / _M_KM)     # PYCHOK mean, spherical earth radius (C{KM}, Kilo meter)
R_NM  = _Radius(R_NM=R_M / _M_NM)     # PYCHOK mean, spherical earth radius (C{NM}, nautical miles)
R_SM  = _Radius(R_SM=R_M / _M_SM)     # PYCHOK mean, spherical earth radius (C{SM}, statute miles)
# See <https://www.EdWilliams.org/avform.htm>, <https://www.DTIC.mil/dtic/tr/fulltext/u2/a216843.pdf>
# and <https://GitHub.com/NASA/MultiDop/blob/master/src/share/man/man3/geog_lib.3> based on the
# International Standard Nautical Mile of 1,852 meter (1' latitude)
R_FM  = _Radius(R_FM=6371000.0)        # PYCHOK former FAI Sphere earth radius (C{meter})
R_GM  = _Radius(R_GM=6371230.0)        # PYCHOK avg. radius, distance to geoid surface (C{meter})
# <http://Wiki.GIS.com/wiki/index.php/Ellipsoidal_quadratic_mean_radius>
R_QM  = _Radius(R_QM=6372797.560856)   # PYCHOK earth' quadratic mean radius (C{meter})
# Rtri= _Radius(Rtri=6372797.5559594)  # PYCHOK Rtriaxial quadratic mean radius (C{meter}), WGS84
# Rbi = _Radius(Rbi =6367453.6345163)  # PYCHOK Rbiaxial quadratic mean radius (C{meter}), WGS84
R_VM  = _Radius(R_VM=6366707.0194937)  # PYCHOK aViation/naVigation earth radius (C{meter})
# R_AU=  Meter( R_AU=149597870700.0)   # PYCHOK <https://WikiPedia.org/wiki/Astronomical_unit>

_INF_NAN_NINF =  INF, NAN, NINF
_pos_self     = _1_0.__pos__() is _1_0  # PYCHOK in .fsums, .vector3dBase


def _0_0s(n):
    '''(INTERNAL) Return an C{B{n}-tuple} of C{_0_0} zeros.
    '''
    return _0_0_9T[:n] if 0 <= n <= len(_0_0_9T) else (_0_0_1T * n)


try:
    from math import isclose as _isclose
except ImportError:  # Python 3.4-

    def _isclose(a, b, rel_tol=1e-9, abs_tol=0):
        '''Mimick Python 3.5+ C{math.isclose}.
        '''
        t, d = abs_tol, fabs(a - b)
        if d > t:
            r = max(fabs(a), fabs(b)) * rel_tol
            t = max(r, t)
        return d <= t


def isclose(a, b, rel_tol=1e-12, abs_tol=EPS0):
    '''Like C{math.isclose}, but with defaults such
       that C{isclose(0, EPS0)} is C{True} by default.
    '''
    return _isclose(a, b, rel_tol=rel_tol, abs_tol=abs_tol)


try:
    from math import isfinite as _isfinite  # in .ellipsoids, .fsums, .karney
except ImportError:  # Python 3.1-

    def _isfinite(x):
        '''Mimick Python 3.2+ C{math.isfinite}.
        '''
        return not (isinf(x) or isnan(x))


def isfinite(obj):
    '''Check a finite C{scalar} or C{complex} value.

       @arg obj: Value (C{scalar} or C{complex}).

       @return: C{False} if B{C{obj}} is C{INF}, C{NINF}
                or C{NAN}, C{True} otherwise.

       @raise TypeError: Non-scalar and non-complex B{C{obj}}.
    '''
    try:
        return (obj not in _INF_NAN_NINF) and _isfinite(obj)
    except Exception as x:
        if iscomplex(obj):  # _isfinite(complex) thows TypeError
            return isfinite(obj.real) and isfinite(obj.imag)
        raise _xError(x, Fmt.PAREN(isfinite.__name__, obj))


def isint0(obj, both=False):
    '''Check for L{INT0} or C{int(0)} value.

       @arg obj: The object (any C{type}).
       @kwarg both: If C{true}, also check C{float(0)} (C{bool}).

       @return: C{True} if B{C{obj}} is L{INT0}, C{int(0)} or
                C{float(0)}, C{False} otherwise.
    '''
    return (obj is INT0 or obj is int(0) or bool(both and
       (not obj) and isint(obj, both=True))) and not isbool(obj)


def isnear0(x, eps0=EPS0):
    '''Is B{C{x}} near I{zero} within a tolerance?

       @arg x: Value (C{scalar}).
       @kwarg eps0: Near-I{zero} tolerance (C{EPS0}).

       @return: C{True} if C{abs(B{x}) < B{eps0}},
                C{False} otherwise.

       @see: Function L{isnon0}.
    '''
    return bool(eps0 > x > -eps0)


def isnear1(x, eps1=EPS0):
    '''Is B{C{x}} near I{one} within a tolerance?

       @arg x: Value (C{scalar}).
       @kwarg eps1: Near-I{one} tolerance (C{EPS0}).

       @return: C{isnear0(B{x} - 1, eps0=B{eps1})}.

       @see: Function L{isnear0}.
    '''
    return bool(eps1 > (x - _1_0) > -eps1)


def isnear90(x, eps90=EPS0):
    '''Is B{C{x}} near I{90} within a tolerance?

       @arg x: Value (C{scalar}).
       @kwarg eps90: Near-I{90} tolerance (C{EPS0}).

       @return: C{isnear0(B{x} - 90)}.

       @see: Function L{isnear0}.
    '''
    return bool(eps90 > (x - _90_0) > -eps90)


def isneg0(x):
    '''Check for L{NEG0}, negative C{0.0}.

       @arg x: Value (C{scalar}).

       @return: C{True} if B{C{x}} is C{NEG0} or C{-0.0},
                C{False} otherwise.
    '''
    return (not x) and _copysign(1, x) < 0
#                  and str(x).startswith(_MINUS_)


def isninf(x):
    '''Check for L{NINF}, negative C{INF}.

       @arg x: Value (C{scalar}).

       @return: C{True} if B{C{x}} is C{NINF} or C{-inf},
                C{False} otherwise.
    '''
    return x is NINF or (x < 0 and not isfinite(x))


def isnon0(x, eps0=EPS0):
    '''Is B{C{x}} non-zero with a tolerance?

       @arg x: Value (C{scalar}).
       @kwarg eps0: Non-zero tolerance (C{EPS0}).

       @return: C{True} if C{abs(B{x}) > B{eps0}},
                C{False} otherwise.

       @see: Function L{isnear0}.
    '''
    return not bool(eps0 > x > -eps0)  # not isnear0


def _off90(lat):
    '''(INTERNAL) Off 90.0 for .gars and .wgrs.
    '''
    return max(min(lat, _89_999_), -_89_999_)


try:
    from math import remainder
except ImportError:  # Python 3.6-
    from math import fmod as _fmod

    def remainder(x, y):
        '''Mimick Python 3.7+ C{math.remainder}.
        '''
        if isnan(y):
            x =  NAN
        elif x and not isnan(x):
            y =  fabs(y)
            x = _fmod(x, y)
            h = _0_5 * y
            if x >= h:
                x -= y
            elif x < -h:
                x += y
        return x  # keep signed 0.0


def _umod_360(deg):
    '''(INTERNAL) Non-negative C{deg} modulo 360, basic C{.utily.wrap360}.
    '''
    return (deg % _360_0) or _0_0


def _umod_PI2(rad):
    '''(INTERNAL) Non-negative C{rad} modulo PI2, basic C{.utily.wrapPI2}.
    '''
    return (rad % PI2) or _0_0


if __name__ == '__main__':

    from pygeodesy.errors import itemsorted
    from pygeodesy.lazily import printf

    t = n = v = []
    for n, v in itemsorted(locals()):
        if isinstance(v, (Float, Int, Radius)):
            printf('%9s: %r', n, v.toRepr(std=False))
            if v.name != n:
                raise AssertionError('%r != %r' % (n, v))
            if v.name is not n:
                raise AssertionError('%r is not %r' % (n, v))
            if not n.startswith(_UNDER_):
                t.append(n)
    t.append(float_.__name__)
    printf('__all__ = %r', tuple(t))

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
