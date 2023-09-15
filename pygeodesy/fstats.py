
# -*- coding: utf-8 -*-

u'''Classes for I{running} statistics and regressions based on
L{pygeodesy.Fsum}, precision floating point summation.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import isodd, islistuple, _xinstanceof, \
                            _xsubclassof, _zip
from pygeodesy.constants import _0_0, _2_0, _3_0, _4_0, _6_0,  _xError
# from pygeodesy.errors import _xError  # from .constants
from pygeodesy.fmath import hypot2,  sqrt
from pygeodesy.fsums import _2float, Fsum,  Fmt
from pygeodesy.interns import NN, _iadd_op_, _invalid_, _other_, _SPACE_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import _Named, _NotImplemented, notOverloaded, \
                             property_RO
# from pygeodesy.props import property_RO  # from .named
# from pygeodesy.streprs import Fmt  # from .fsums

# from math import sqrt  # from .fmath

__all__ = _ALL_LAZY.fstats
__version__ = '23.09.08'

_Floats =  Fsum, float
_Scalar = _Floats + (int,)  # XXX basics._Ints is ABCMeta in Python 2
try:
    _Scalar += (long,)
except NameError:  # Python 3+
    pass


def _2Floats(xs, ys=False):
    '''(INTERNAL) Yield each value as C{float} or L{Fsum}.
    '''
    if ys:
        def _2f(i, x):
            return _2float(index=i, ys=x)
    else:
        def _2f(i, x):  # PYCHOK redef
            return _2float(index=i, xs=x)

    for i, x in enumerate(xs):  # don't unravel Fsums
        yield x if isinstance(x, _Floats) else _2f(i, x)


def _sampled(n, sample):
    '''(INTERNAL) Return the sample or the entire count.
    '''
    return (n - 1) if sample and n > 0 else n


class _FstatsNamed(_Named):
    '''(INTERNAL) Base class.
    '''
    _n = 0

    def __add__(self, other):
        '''Sum of this and a scalar, an L{Fsum} or an other instance.
        '''
        f  = self.fcopy(name=self.__add__.__name__)  # PYCHOK expected
        f += other
        return f

    def __float__(self):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self)

    def __int__(self):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self)

    def __len__(self):
        '''Return the I{total} number of accumulated values (C{int}).
        '''
        return self._n

    def __neg__(self):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self)

    def __radd__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __str__(self):
        return Fmt.SQUARE(self.named3, len(self))

    def fcopy(self, deep=False, name=NN):
        '''Copy this instance, C{shallow} or B{C{deep}}.
        '''
        n =  name or self.fcopy.__name__
        f = _Named.copy(self, deep=deep, name=n)
        return self._copy(f, self)  # PYCHOK expected

    copy = fcopy


class _FstatsBase(_FstatsNamed):
    '''(INTERNAL) Base running stats class.
    '''
    _Ms = ()

    def _copy(self, c, s):
        '''(INTERNAL) Copy C{B{c} = B{s}}.
        '''
        _xinstanceof(self.__class__, c=c, s=s)
        c._Ms = tuple(M.fcopy() for M in s._Ms)  # deep=False
        c._n  =                          s._n
        return c

    def fadd(self, xs, sample=False):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}, see function C{notOverloaded}.
        '''
        notOverloaded(self, xs, sample=sample)

    def fadd_(self, *xs, **sample):
        '''Accumulate and return the current count.

           @see: Method C{fadd}.
        '''
        return self.fadd(xs, **sample)

    def fmean(self, xs=None):
        '''Accumulate and return the current mean.

           @kwarg xs: Iterable with additional values (C{Scalar}s).

           @return: Current, running mean (C{float}).

           @see: Method C{fadd}.
        '''
        if xs:
            self.fadd(xs)
        return self._M1.fsum()

    def fmean_(self, *xs):
        '''Accumulate and return the current mean.

           @see: Method C{fmean}.
        '''
        return self.fmean(xs)

    def fstdev(self, xs=None, sample=False):
        '''Accumulate and return the current standard deviation.

           @kwarg xs: Iterable with additional values (C{Scalar}).
           @kwarg sample: Return the I{sample} instead of the entire
                          I{population} standard deviation (C{bool}).

           @return: Current, running (sample) standard deviation (C{float}).

           @see: Method C{fadd}.
        '''
        v = self.fvariance(xs, sample=sample)
        return sqrt(v) if v > 0 else _0_0

    def fstdev_(self, *xs, **sample):
        '''Accumulate and return the current standard deviation.

           @see: Method C{fstdev}.
        '''
        return self.fstdev(xs, **sample)

    def fvariance(self, xs=None, sample=False):
        '''Accumulate and return the current variance.

           @kwarg xs: Iterable with additional values (C{Scalar}s).
           @kwarg sample: Return the I{sample} instead of the entire
                          I{population} variance (C{bool}).

           @return: Current, running (sample) variance (C{float}).

           @see: Method C{fadd}.
        '''
        n = self.fadd(xs, sample=sample)
        return float(self._M2 / float(n)) if n > 0 else _0_0

    def fvariance_(self, *xs, **sample):
        '''Accumulate and return the current variance.

           @see: Method C{fvariance}.
        '''
        return self.fvariance(xs, **sample)

    def _iadd_other(self, other):
        '''(INTERNAL) Add Scalar or Scalars.
        '''
        if isinstance(other, _Scalar):
            self.fadd_(other)
        else:
            try:
                if not islistuple(other):
                    raise TypeError(_SPACE_(_invalid_, _other_))
                self.fadd(other)
            except Exception as x:
                raise _xError(x, _SPACE_(self, _iadd_op_, repr(other)))

    @property_RO
    def _M1(self):
        '''(INTERNAL) get the 1st Moment accumulator.'''
        return self._Ms[0]

    @property_RO
    def _M2(self):
        '''(INTERNAL) get the 2nd Moment accumulator.'''
        return self._Ms[1]


class Fcook(_FstatsBase):
    '''U{Cook<https://www.JohnDCook.com/blog/skewness_kurtosis>}'s
       C{RunningStats} computing the running mean, median and
       (sample) kurtosis, skewness, variance, standard deviation
       and Jarque-Bera normality.

       @see: L{Fwelford} and U{Higher-order statistics<https://
             WikiPedia.org/wiki/Algorithms_for_calculating_variance>}.
    '''
    def __init__(self, xs=None, name=NN):
        '''New L{Fcook} stats accumulator.

           @kwarg xs: Iterable with initial values (C{Scalar}s).
           @kwarg name: Optional name (C{str}).

           @see: Method L{Fcook.fadd}.
        '''
        self._Ms = tuple(Fsum() for _ in range(4))  # 1st, 2nd ... Moment
        if name:
            self.name = name
        if xs:
            self.fadd(xs)

    def __iadd__(self, other):
        '''Add B{C{other}} to this L{Fcook} instance.

           @arg other: An L{Fcook} instance or C{Scalar}s, meaning
                       one or more C{scalar} or L{Fsum} instances.

           @return: This instance, updated (L{Fcook}).

           @raise TypeError: Invalid B{C{other}} type.

           @raise ValueError: Invalid B{C{other}}.

           @see: Method L{Fcook.fadd}.
        '''
        if isinstance(other, Fcook):
            nb = len(other)
            if nb > 0:
                na = len(self)
                if na > 0:
                    A1, A2, A3, A4 = self._Ms
                    B1, B2, B3, B4 = other._Ms

                    n   = na + nb
                    n_  = float(n)
                    D   = A1 - B1  # b1 - a1
                    Dn  = D / n_
                    Dn2 = Dn**2  # d**2 / n**2
                    nab = na * nb
                    Dn3 = Dn2 * (D * nab)

                    na2 =  na**2
                    nb2 =  nb**2
                    A4 +=  B4
                    A4 += (B3 * na  - (A3 * nb))  * (Dn  * _4_0)
                    A4 += (B2 * na2 + (A2 * nb2)) * (Dn2 * _6_0)
                    A4 += (Dn * Dn3) * (na2 - nab + nb2)  # d**4 / n**3

                    A3 +=  B3
                    A3 += (A2 *  na - (B2 * nb)) * (Dn * _3_0)
                    A3 += Dn3 * (na - nb)

                    A2 += B2
                    A2 += Dn2 * (nab / n_)

                    B1n = B1 * nb  # if other is self
                    A1 *= na
                    A1 += B1n
                    A1 *= 1 / n_  # /= chokes PyChecker

#                   self._Ms = A1, A2, A3, A4
                    self._n  = n
                else:
                    self._copy(self, other)
        else:
            self._iadd_other(other)
        return self

    def fadd(self, xs, sample=False):
        '''Accumulate and return the current count.

           @arg xs: Iterable with additional values (C{Scalar}s,
                    meaning C{scalar} or L{Fsum} instances).
           @kwarg sample: Return the I{sample} instead of the entire
                          I{population} count (C{bool}).

           @return: Current, running (sample) count (C{int}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{xs}} value.

           @raise ValueError: Invalid or non-finite B{C{xs}} value.

           @see: U{online_kurtosis<https://WikiPedia.org/wiki/
                 Algorithms_for_calculating_variance>}.
        '''
        n = self._n
        if xs:
            M1, M2, M3, M4 = self._Ms
            for x in _2Floats(xs):
                n1 = n
                n += 1
                D  = x - M1
                Dn = D / n
                if Dn:
                    Dn2 = Dn**2
                    if n1 > 1:
                        T1 = D  * (Dn  *  n1)
                        T2 = T1 * (Dn  * (n1 - 1))
                        T3 = T1 * (Dn2 * (n**2 - 3 * n1))
                    elif n1 > 0:  # n1 == 1, n == 2
                        T1 =  D * Dn
                        T2 = _0_0
                        T3 = T1 * Dn2
                    else:
                        T1 = T2 = T3 = _0_0
                    M4 += T3
                    M4 -= M3 * (Dn  * _4_0)
                    M4 += M2 * (Dn2 * _6_0)

                    M3 += T2
                    M3 -= M2 * (Dn  * _3_0)

                    M2 += T1
                    M1 += Dn
#           self._Ms = M1, M2, M3, M4
            self._n  = n
        return _sampled(n, sample)

    def fjb(self, xs=None, sample=True, excess=True):
        '''Accumulate and compute the current U{Jarque-Bera
           <https://WikiPedia.org/wiki/Jarque–Bera_test>} normality.

           @kwarg xs: Iterable with additional values (C{Scalar}s).
           @kwarg sample: Return the I{sample} normality (C{bool}), default.
           @kwarg excess: Return the I{excess} kurtosis (C{bool}), default.

           @return: Current, running (sample) Jarque-Bera normality (C{float}).

           @see: Method L{Fcook.fadd}.
        '''
        n = self.fadd(xs,  sample=sample)
        k = self.fkurtosis(sample=sample, excess=excess) / _2_0
        s = self.fskewness(sample=sample)
        return n * hypot2(k, s) / _6_0

    def fjb_(self, *xs, **sample_excess):
        '''Accumulate and compute the current U{Jarque-Bera
           <https://WikiPedia.org/wiki/Jarque–Bera_test>} normality.

           @see: Method L{Fcook.fjb}.
        '''
        return self.fjb(xs, **sample_excess)

    def fkurtosis(self, xs=None, sample=False, excess=True):
        '''Accumulate and return the current kurtosis.

           @kwarg xs: Iterable with additional values (C{Scalar}s).
           @kwarg sample: Return the I{sample} instead of the entire
                          I{population} kurtosis (C{bool}).
           @kwarg excess: Return the I{excess} kurtosis (C{bool}), default.

           @return: Current, running (sample) kurtosis or I{excess} kurtosis (C{float}).

           @see: U{Kurtosis Formula<https://www.Macroption.com/kurtosis-formula>}
                 and U{Mantalos<https://www.ResearchGate.net/publication/227440210>}.

           @see: Method L{Fcook.fadd}.
        '''
        k, n = _0_0, self.fadd(xs, sample=sample)
        if n > 0:
            _, M2, _, M4 = self._Ms
            m2 = float(M2 * M2)
            if m2:
                K, x = (M4 * (n / m2)), _3_0
                if sample and 2 < n < len(self):
                    d  = float((n - 1) * (n - 2))
                    K *=       (n + 1) * (n + 2) / d
                    x *= n**2 / d
                if excess:
                    K -= x
                k = K.fsum()
        return k

    def fkurtosis_(self, *xs, **sample_excess):
        '''Accumulate and return the current kurtosis.

           @see: Method L{Fcook.fkurtosis}.
        '''
        return self.fkurtosis(xs, **sample_excess)

    def fmedian(self, xs=None):
        '''Accumulate and return the current median.

           @kwarg xs: Iterable with additional values (C{Scalar}s).

           @return: Current, running median (C{float}).

           @see: U{Pearson's Skewness Coefficients<https://MathWorld.Wolfram.com/
                 PearsonsSkewnessCoefficients.html>}, U{Skewness & Kurtosis Simplified
                 https://TowardsDataScience.com/skewness-kurtosis-simplified-1338e094fc85>}
                 and method L{Fcook.fadd}.
        '''
        # skewness = 3 * (mean - median) / stdev, i.e.
        # median = mean - skewness * stdef / 3
        m = float(self._M1) if xs is None else self.fmean(xs)
        return m - self.fskewness() * self.fstdev() / _3_0

    def fmedian_(self, *xs):
        '''Accumulate and return the current median.

           @see: Method L{Fcook.fmedian}.
        '''
        return self.fmedian(xs)

    def fskewness(self, xs=None, sample=False):
        '''Accumulate and return the current skewness.

           @kwarg xs: Iterable with additional values (C{Scalar}s).
           @kwarg sample: Return the I{sample} instead of the entire
                          I{population} skewness (C{bool}).

           @return: Current, running (sample) skewness (C{float}).

           @see: U{Skewness Formula<https://www.Macroption.com/skewness-formula/>}
                 and U{Mantalos<https://www.ResearchGate.net/publication/227440210>}.

           @see: Method L{Fcook.fadd}.
        '''
        s, n = _0_0, self.fadd(xs, sample=sample)
        if n > 0:
            _, M2, M3, _ = self._Ms
            m = float(M2**3)
            if m > 0:
                S = M3 * sqrt(float(n) / m)
                if sample and 1 < n < len(self):
                    S *= (n + 1) / float(n - 1)
                s = S.fsum()
        return s

    def fskewness_(self, *xs, **sample):
        '''Accumulate and return the current skewness.

           @see: Method L{Fcook.fskewness}.
        '''
        return self.fskewness(xs, **sample)

    def toFwelford(self, name=NN):
        '''Return an L{Fwelford} equivalent.
        '''
        f = Fwelford(name=name or self.name)
        f._Ms = self._M1.fcopy(), self._M2.fcopy()  # deep=False
        f._n  = self._n
        return f


class Fwelford(_FstatsBase):
    '''U{Welford<https://WikiPedia.org/wiki/Algorithms_for_calculating_variance>}'s
       accumulator computing the running mean, (sample) variance and standard deviation.

       @see: U{Cook<https://www.JohnDCook.com/blog/standard_deviation/>} and L{Fcook}.
    '''
    def __init__(self, xs=None, name=NN):
        '''New L{Fwelford} stats accumulator.

           @kwarg xs: Iterable with initial values (C{Scalar}s).
           @kwarg name: Optional name (C{str}).

           @see: Method L{Fwelford.fadd}.
        '''
        self._Ms = Fsum(), Fsum()  # 1st and 2nd Moment
        if name:
            self.name = name
        if xs:
            self.fadd(xs)

    def __iadd__(self, other):
        '''Add B{C{other}} to this L{Fwelford} instance.

           @arg other: An L{Fwelford} or L{Fcook} instance or C{Scalar}s,
                       meaning one or more C{scalar} or L{Fsum} instances.

           @return: This instance, updated (L{Fwelford}).

           @raise TypeError: Invalid B{C{other}} type.

           @raise ValueError: Invalid B{C{other}}.

           @see: Method L{Fwelford.fadd} and U{Parallel algorithm<https//
                 WikiPedia.org/wiki/Algorithms_for_calculating_variance>}.
        '''
        if isinstance(other, Fwelford):
            nb = len(other)
            if nb > 0:
                na = len(self)
                if na > 0:
                    M,  S  = self._Ms
                    M_, S_ = other._Ms

                    n  = na + nb
                    n_ = float(n)

                    D  = M_ - M
                    D *= D  # D**2
                    D *= na * nb / n_
                    S += D
                    S += S_

                    Mn = M_ * nb  # if other is self
                    M *= na
                    M += Mn
                    M *= 1 / n_  # /= chokes PyChecker

#                   self._Ms = M, S
                    self._n  = n
                else:
                    self._copy(self, other)

        elif isinstance(other, Fcook):
            self += other.toFwelford()
        else:
            self._iadd_other(other)
        return self

    def fadd(self, xs, sample=False):
        '''Accumulate and return the current count.

           @arg xs: Iterable with additional values (C{Scalar}s,
                    meaning C{scalar} or L{Fsum} instances).
           @kwarg sample: Return the I{sample} instead of the entire
                          I{population} count (C{bool}).

           @return: Current, running (sample) count (C{int}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{xs}} value.

           @raise ValueError: Invalid or non-finite B{C{xs}} value.
        '''
        n = self._n
        if xs:
            M, S = self._Ms
            for x in _2Floats(xs):
                n += 1
                D  = x - M
                M += D / n
                D *= x - M
                S += D
#           self._Ms = M, S
            self._n = n
        return _sampled(n, sample)


class Flinear(_FstatsNamed):
    '''U{Cook<https://www.JohnDCook.com/blog/running_regression>}'s
       C{RunningRegression} computing the running slope, intercept
       and correlation of a linear regression.
    '''
    def __init__(self, xs=None, ys=None, Fstats=Fwelford, name=NN):
        '''New L{Flinear} regression accumulator.

           @kwarg xs: Iterable with initial C{x} values (C{Scalar}s).
           @kwarg ys: Iterable with initial C{y} values (C{Scalar}s).
           @kwarg Fstats: Stats class for C{x} and C{y} values (L{Fcook}
                          or L{Fwelford}).
           @kwarg name: Optional name (C{str}).

           @raise TypeError: Invalid B{C{Fs}}, not L{Fcook} or
                              L{Fwelford}.
           @see: Method L{Flinear.fadd}.
        '''
        _xsubclassof(Fcook, Fwelford, Fstats=Fstats)
        if name:
            self.name = name

        self._S = Fsum(name=name)
        self._X = Fstats(name=name)
        self._Y = Fstats(name=name)
        if xs and ys:
            self.fadd(xs, ys)

    def __iadd__(self, other):
        '''Add B{C{other}} to this instance.

           @arg other: An L{Flinear} instance or C{Scalar} pairs,
                       meaning C{scalar} or L{Fsum} instances.

           @return: This instance, updated (L{Flinear}).

           @raise TypeError: Invalid B{C{other}} or the B{C{other}}
                             and these C{x} and C{y} accumulators
                             are not compatible.

           @raise ValueError: Invalid or odd-length B{C{other}}.

           @see: Method L{Flinear.fadd_}.
        '''
        if isinstance(other, Flinear):
            if len(other) > 0:
                if len(self) > 0:
                    n = other._n
                    S = other._S
                    X = other._X
                    Y = other._Y
                    D = (X._M1 - self._X._M1) * \
                        (Y._M1 - self._Y._M1) * \
                        (n     * self._n / float(n + self._n))
                    self._n += n
                    self._S += S + D
                    self._X += X
                    self._Y += Y
                else:
                    self._copy(self, other)
        else:
            try:
                if not islistuple(other):
                    raise TypeError(_SPACE_(_invalid_, _other_))
                elif isodd(len(other)):
                    raise ValueError(Fmt.PAREN(isodd=Fmt.PAREN(len=_other_)))
                self.fadd_(*other)
            except Exception as x:
                raise _xError(x, _SPACE_(self, _iadd_op_, repr(other)))
        return self

    def _copy(self, c, s):
        '''(INTERNAL) Copy C{B{c} = B{s}}.
        '''
        _xinstanceof(Flinear, c=c, s=s)
        c._n = s._n
        c._S = s._S.fcopy(deep=False)
        c._X = s._X.fcopy(deep=False)
        c._Y = s._Y.fcopy(deep=False)
        return c

    def fadd(self, xs, ys, sample=False):
        '''Accumulate and return the current count.

           @arg xs: Iterable with additional C{x} values (C{Scalar}s),
                    meaning C{scalar} or L{Fsum} instances).
           @arg ys: Iterable with additional C{y} values (C{Scalar}s,
                    meaning C{scalar} or L{Fsum} instances).
           @kwarg sample: Return the I{sample} instead of the entire
                          I{population} count (C{bool}).

           @return: Current, running (sample) count (C{int}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{xs}} or B{C{ys}} value.

           @raise ValueError: Invalid or non-finite B{C{xs}} or B{C{ys}} value.
        '''
        n = self._n
        if xs and ys:
            S = self._S
            X = self._X
            Y = self._Y
            for x, y in _zip(_2Floats(xs), _2Floats(ys, ys=True)):  # strict=True
                n1 = n
                n += 1
                if n1 > 0:
                    S += (X._M1 - x) * (Y._M1 - y) * (n1 / float(n))
                X += x
                Y += y
            self._n  = n
        return _sampled(n, sample)

    def fadd_(self, *x_ys, **sample):
        '''Accumulate and return the current count.

           @arg x_ys: Individual, alternating C{x, y, x, y, ...}
                      positional values (C{Scalar}s).

           @see: Method C{Flinear.fadd}.
        '''
        return self.fadd(x_ys[0::2], x_ys[1::2], **sample)

    def fcorrelation(self, sample=False):
        '''Return the current, running (sample) correlation (C{float}).

           @kwarg sample: Return the I{sample} instead of the entire
                          I{population} correlation (C{bool}).
        '''
        return self._sampled(self.x.fstdev(sample=sample) *
                             self.y.fstdev(sample=sample), sample)

    def fintercept(self, sample=False):
        '''Return the current, running (sample) intercept (C{float}).

           @kwarg sample: Return the I{sample} instead of the entire
                          I{population} intercept (C{bool}).
        '''
        return float(self.y._M1 -
                    (self.x._M1 * self.fslope(sample=sample)))

    def fslope(self, sample=False):
        '''Return the current, running (sample) slope (C{float}).

           @kwarg sample: Return the I{sample} instead of the entire
                          I{population} slope (C{bool}).
        '''
        return self._sampled(self.x.fvariance(sample=sample), sample)

    def _sampled(self, t, sample):
        '''(INTERNAL) Compute the sampled or entire population result.
        '''
        t *= float(_sampled(self._n, sample))
        return float(self._S / t) if t else _0_0

    @property_RO
    def x(self):
        '''Get the C{x} accumulator (L{Fcook} or L{Fwelford}).
        '''
        return self._X

    @property_RO
    def y(self):
        '''Get the C{y} accumulator (L{Fcook} or L{Fwelford}).
        '''
        return self._Y


__all__ += _ALL_DOCS(_FstatsBase, _FstatsNamed)

# **) MIT License
#
# Copyright (C) 2021-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
