
# -*- coding: utf-8 -*-

u'''Classes for running statistics based on L{Fsum}, precision
floating point summation.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

# from pygeodesy.basics import isscalar  # from .fmath
# from pygeodesy.errors import _TypeError  # from .fmath
from pygeodesy.fmath import _2float, Fsum, hypot2, isscalar, sqrt, \
                            _TypeError
from pygeodesy.interns import NN, _iadd_, _SPACE_, _0_0, \
                             _1_5, _2_0, _3_0, _4_0, _6_0
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import _Named, _NotImplemented

# from math import sqrt  # pow from .fmath

__all__ = _ALL_LAZY.fstats
__version__ = '22.01.09'

_Float =  Fsum, float
_N_3_0 = -_3_0
_N_4_0 = -_4_0


class _StatsBase(_Named):
    _n = 0

    def __len__(self):
        '''Return the I{total} number of accumulated values (C{int}).
        '''
        return self._n

    def fadd_(self, *xs, **sample):
        '''Accumulate and return the current count.

           @see: Method C{fadd}.
        '''
        return self.fadd(xs, **sample)  # PYCHOK expected

    def _2Float(self, xs):
        '''(INTERNAL) Yield each value as C{float}.
        '''
        for i, x in enumerate(xs):
            yield x if isinstance(x, _Float) else _2float(index=i, xs=x)

    def fmean_(self, *xs):
        '''Accumulate and return the current mean.

           @see: Method C{fmean}.
        '''
        return self.fmean(xs)  # PYCHOK expected

    def fstdev_(self, *xs, **sample):
        '''Accumulate and return the current standard deviation.

           @see: Method C{fstdev}.
        '''
        return self.fstdev(xs, **sample)  # PYCHOK expected

    def fvariance_(self, *xs, **sample):
        '''Accumulate and return the current variance.

           @see: Method C{fvariance}.
        '''
        return self.fvariance(xs, **sample)  # PYCHOK expected


class Fcook(_StatsBase):
    '''U{Cook's<https://www.JohnDCook.com/blog/skewness_kurtosis>}
       C{RunningStats} computing the running kurtosis, mean, skewness,
       (sample) variance, standard deviation and Jarque-Bera normality.

       @see: L{Fwelford}.
    '''
    def __init__(self, xs=(), name=NN):
        '''New L{Fcook} accumulator.

           @kwarg xs: Iterable with initial values (C{Scalar}s).
           @kwarg name: Optional name (C{str}).

           @see: Method L{Fcook.fadd}.
        '''
        self._Ms = tuple(Fsum() for _ in range(4))  # 1st, 2nd ... Moment
        self._ms = (_0_0,) * len(self._Ms)  # Moment values
        if name:
            self.name = name
        if xs:
            self.fadd(xs)

    def __add__(self, other):
        '''Sum of this and a scalar, an L{Fsum} or an other instance.

           @arg other: L{Fcook} instance or a C{Scalar}.

           @return: The sum (L{Fcook}).

           @see: Methods L{Fcook.__iadd__} and L{Fcook.fadd}.
        '''
        f = self.fcopy(name=self.__add__.__name__)
        f += other
        return f

    def __iadd__(self, other):
        '''Add a scalar, an L{Fsum} or an other instance to this instance.

           @arg other: L{Fcook} instance or a C{Scalar}.

           @return: This instance, updated (L{Fcook}).

           @raise TypeError: Invalid B{C{other}} type.

           @see: Method L{Fcook.fadd}.
        '''
        if isinstance(other, Fcook):
            nb = len(other)
            if nb > 0:
                na = len(self)
                if na > 0:
                    A1, A2, A3, A4 = self._Ms
                    B1, B2, B3, B4 = other._dupMs() if \
                                     other is self else other._Ms

                    d  = other._ms[0] - self._ms[0]  # b1 - a1
                    p  = na * nb
                    n  = na + nb
                    nf = float(n)
                    q1 = d / nf
                    q2 = q1**2  # d**2 / n**2
                    qp = q2 * p

                    na2 =  na**2
                    nb2 =  nb**2
                    A4 +=  B4
                    A4 += (B3 * na  - (A3 * nb))  * (_4_0 * q1)
                    A4 += (B2 * na2 + (A2 * nb2)) * (_6_0 * q2)
                    A4 += (na2 - p + nb2) * qp * (d * q1)  # d**4 / n**3

                    A3 +=  B3
                    A3 += (A2 * na - (B2 * nb)) * (_3_0 * q1)
                    A3 += (na - nb) * qp * d

                    A2 += B2
                    A2 += nf * qp  # d**2 * p / n

                    A1 *= na
                    A1 += B1 * nb
                    A1 *= 1 / nf  # /= chokes PyChecker

#                   self._Ms = A1, A2, A3, A4
                    self._ms = tuple(M.fsum() for M in self._Ms)
                    self._n  = n
                else:
                    self._Ms = tuple(other._dupMs())  # deep=False
                    self._ms = tuple(other._ms)
                    self._n  = nb

        elif isinstance(other, _Float) or isscalar(other):
            self.fadd_(other)
        else:
            raise _TypeError(_SPACE_(self, _iadd_, repr(other)))
        return self

    def __radd__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def _dupMs(self, deep=False):
        for M in self._Ms:
            yield M.fcopy(deep=deep)

    def fadd(self, xs, sample=False):
        '''Accumulate and return the current count.

           @kwarg xs: Iterable with additional values (C{Scalar}s).
           @kwarg sample: Return the I{sample} instead of the full
                          I{population} value (C{bool}).

           @note: C{Scalar} means an L{Fsum} instance or C{scalar}.

           @return: Current, running count (C{float}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{xs}} value.

           @raise ValueError: Invalid or non-finite B{C{xs}} value.
        '''
        n = self._n
        if xs:
            m1, m2, m3, m4 = self._ms
            M1, M2, M3, M4 = (M.fsum_ for M in self._Ms)
            for x in self._2Float(xs):
                n1 = n
                n += 1
                d  = x - m1
                d1 = float(d)
                q1 = float(d / n)
                q2 = q1**2
                if n1 > 0 and q1 and d1:
                    t1 = d1 * q1 *  n1
                    t2 = t1 * q1 * (n - 2)
                    t3 = t1 * q2 * (n**2 - n1 * 3)
                else:
                    t1 = t2 = t3 = _0_0
                m4 = M4(t3, _N_4_0 * q1 * m3, _6_0 * q2 * m2)
                m3 = M3(t2, _N_3_0 * q1 * m2)
                m2 = M2(t1)
                m1 = M1(q1)
            self._ms = m1, m2, m3, m4
            self._n = n
        return n - (1 if sample and n else 0)

    def fcopy(self, deep=False, name=NN):
        '''Copy this instance, C{shallow} or B{C{deep}}.

           @return: The copy (L{Fwelford}).
         '''
        f = _Named.copy(self, deep=deep, name=name)
        f._Ms = tuple(self._dupMs())  # deep=False
        f._ms = tuple(self._ms)
        f._n  =       self._n
        return f

    copy = fcopy

    def fjb(self, xs=(), sample=True, excess=True):
        '''Accumulate and compute the current U{Jarque-Bera
           <https://WikiPedia.org/wiki/Jarque–Bera_test>} normality.

           @kwarg xs: Iterable with additional values (C{Scalar}s).
           @kwarg sample: Return the I{sample} value (C{bool}), default.
           @kwarg excess: Return the I{excess} kurtosis (C{bool}), default.

           @return: Current, running Jarque-Bera normality (C{float}).

           @see: Method L{Fcook.fadd}.
        '''
        n = self.fadd(xs,  sample=sample)
        k = self.fkurtosis(sample=sample, excess=excess) / _2_0
        s = self.fskewness(sample=sample)
        return n * hypot2(k, s) / _6_0

    def fjb_(self, *xs, **excess_sample):
        '''Accumulate and compute the current U{Jarque-Bera
           <https://WikiPedia.org/wiki/Jarque–Bera_test>} normality.

           @see: Method L{Fcook.fjb}.
        '''
        return self.fjb(xs, **excess_sample)

    def fkurtosis(self, xs=(), sample=False, excess=True):
        '''Accumulate and return the current kurtosis.

           @kwarg xs: Iterable with additional values (C{Scalar}s).
           @kwarg sample: Return the I{sample} instead of the full
                          I{population} value (C{bool}).
           @kwarg excess: Return the I{excess} kurtosis (C{bool}), default.

           @return: Current, running kurtosis or I{excess} kurtosis (C{float}).

           @see: U{Kurtosis Formula<https://www.Macroption.com/kurtosis-formula>}
                 and U{Mantalos<https://www.researchgate.net/publication/
                 227440210_Three_different_measures_of_sample_skewness_and_kurtosis_and_their_effects_on_the_JarqueBera_test_for_normality>}.

           @see: Method L{Fcook.fadd}.
        '''
        k, n = _0_0, self.fadd(xs, sample=sample)
        if n > 0:
            _, m2, _, m4 = self._ms
            m2 *= m2
            if m2:
                k, x = n * m4 / m2, _3_0
                if sample and 2 < n < len(self):
                    d  = (n - 1) * (n - 2)
                    k *= (n + 1) * (n + 2) / float(d)
                    x *=  n**2 / float(d)
                if excess:
                    k -= x
        return k

    def fkurtosis_(self, *xs, **excess_sample):
        '''Accumulate and return the current kurtosis.

           @see: Method L{Fcook.fkurtosis}.
        '''
        return self.fkurtosis(xs, **excess_sample)

    def fmean(self, xs=()):
        '''Accumulate and return the current mean.

           @kwarg xs: Iterable with additional values (C{Scalar}s).

           @return: Current, running mean (C{float}).

           @see: Method L{Fcook.fadd}.
        '''
        if xs:
            self.fadd(xs)
        return self._ms[0]

    def fmedian(self, xs=()):
        '''Accumulate and return the current median.

           @kwarg xs: Iterable with additional values (C{Scalar}s).

           @return: Current, running median (C{float}).

           @see: U{Pearson's Skewness Coefficients<https://MathWorld.Wolfram.com/
                 PearsonsSkewnessCoefficients.html>} and U{Skewness & Kurtosis Simplified
                 https://TowardsDataScience.com/skewness-kurtosis-simplified-1338e094fc85>}.

           @see: Method L{Fcook.fadd}.
        '''
        # skewness = 3 * (mean - median) / stdev, i.e.
        # median = mean - skewness * stdef / 3
        m1 = self.fmean(xs) if xs else self._ms[0]
        return m1 - self.fskewness() * self.fstdev() / _3_0

    def fmedian_(self, *xs):
        '''Accumulate and return the current median.

           @see: Method L{Fcook.fmedian}.
        '''
        return self.fmedian(xs)

    def fskewness(self, xs=(), sample=False):
        '''Accumulate and return the current skewness.

           @kwarg xs: Iterable with additional values (C{Scalar}s).
           @kwarg sample: Return the I{sample} instead of the full
                          I{population} value (C{bool}).

           @return: Current, running skewness (C{float}).

           @see: U{Skewness Formula<https://www.Macroption.com/skewness-formula/>}
                 and U{Mantalos<https://www.researchgate.net/publication/
                 227440210_Three_different_measures_of_sample_skewness_and_kurtosis_and_their_effects_on_the_JarqueBera_test_for_normality>}.

           @see: Method L{Fcook.fadd}.
        '''
        s, n = _0_0, self.fadd(xs, sample=sample)
        if n > 0:
            _, m2, m3, _ = self._ms
            m2 = pow(m2, _1_5)
            if m2:
                s = sqrt(float(n)) * m3 / m2
                if sample and 1 < n < len(self):
                    s *= (n + 1) / float(n - 1)
        return s

    def fskewness_(self, *xs, **sample):
        '''Accumulate and return the current skewness.

           @see: Method L{Fcook.fskewness}.
        '''
        return self.fskewness(xs, **sample)

    def fstdev(self, xs=(), sample=False):
        '''Accumulate and return the current standard deviation.

           @kwarg xs: Iterable with additional values (C{Scalar}s).
           @kwarg sample: Return the I{sample} instead of the full
                          I{population} value (C{bool}).

           @return: Current, running standard deviation (C{float}).

           @see: Method L{Fcook.fadd}.
        '''
        v = self.fvariance(xs, sample=sample)
        return sqrt(v) if v > 0 else _0_0

    def fvariance(self, xs=(), sample=False):
        '''Accumulate and return the current variance.

           @kwarg xs: Iterable with additional values (C{Scalar}s).
           @kwarg sample: Return the I{sample} instead of the full
                          I{population} value (C{bool}).

           @return: Current, running variance (C{float}).

           @see: Method L{Fcook.fadd}.
        '''
        n = self.fadd(xs, sample=sample)
        return (self._ms[1] / float(n)) if n > 0 else _0_0


class Fwelford(_StatsBase):
    '''U{Welford's<https://WikiPedia.org/wiki/Algorithms_for_calculating_variance>}
       accumulator computing the running mean, (sample) variance and standard deviation.

       @see: U{Cook<https://www.JohnDCook.com/blog/standard_deviation/>} and L{Fcook}.
    '''
    def __init__(self, xs=(), name=NN):
        '''New L{Fwelford} accumulator.

           @kwarg xs: Iterable with initial values (C{Scalar}s).
           @kwarg name: Optional name (C{str}).

           @see: Method L{Fwelford.fadd}.
        '''
        self._M = Fsum()  # 1st Moment
        self._m = _0_0    # mean
        self._S = Fsum()  # 2nd Moment
        if name:
            self.name = name
        if xs:
            self.fadd(xs)

    def __len__(self):
        '''Return the I{total} number of accumulated values (C{int}).
        '''
        return self._n

    def fadd(self, xs, sample=False):
        '''Accumulate and return the current count.

           @kwarg xs: Iterable with additional values (C{Scalar}s).
           @kwarg sample: Return the I{sample} instead of the full
                          I{population} value (C{bool}).

           @note: C{Scalar} means an L{Fsum} instance or C{scalar}.

           @return: Current, running count (C{float}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{xs}} value.

           @raise ValueError: Invalid or non-finite B{C{xs}} value.
        '''
        n = self._n
        if xs:
            M = self._M.fsum_
            m = self._m
            S = self._S.fadd_
            for x in self._2Float(xs):
                n += 1
                d  = x - m
                m  = M(float(d / n))
                S(float((x - m) * d))
            self._m = m
            self._n = n
        return n - (1 if sample and n else 0)

    def fcopy(self, deep=False, name=NN):
        '''Copy this instance, C{shallow} or B{C{deep}}.

           @return: The copy (L{Fwelford}).
         '''
        f = _Named.copy(self, deep=deep, name=name)
        f._M = self._M.fcopy(deep=False)
        f._m = self._m
        f._n = self._n
        f._S = self._S.fcopy(deep=False)
        return f

    copy = fcopy

    def fmean(self, xs=()):
        '''Accumulate and return the current mean.

           @kwarg xs: Iterable with additional values (C{Scalar}s).

           @return: Current, running mean (C{float}).

           @see: Method L{Fwelford.fadd}.
        '''
        if xs:
            self.fadd(xs)
        return self._m

    def fstdev(self, xs=(), sample=False):
        '''Accumulate and return the current standard deviation.

           @kwarg xs: Iterable with additional values (C{Scalar}s).
           @kwarg sample: Return the I{sample} instead of the full
                          I{population} value (C{bool}).

           @return: Current, running standard deviation (C{float}).

           @see: Method L{Fwelford.fadd}.
        '''
        v = self.fvariance(xs, sample=sample)
        return sqrt(v) if v > 0 else _0_0

    def fvariance(self, xs=(), sample=False):
        '''Accumulate and return the current variance.

           @kwarg xs: Iterable with additional values (C{Scalar}s).
           @kwarg sample: Return the I{sample} instead of the full
                          I{population} value (C{bool}).

           @return: Current, running variance (C{float}).

           @see: Method L{Fwelford.fadd}.
        '''
        n = self.fadd(xs, sample=sample)
        return (self._S.fsum() / float(n)) if n > 0 else _0_0

# **) MIT License
#
# Copyright (C) 2021-2022 -- mrJean1 at Gmail -- All Rights Reserved.
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
