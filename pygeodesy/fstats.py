
# -*- coding: utf-8 -*-

u'''Classes for running statistics and regreesions based on
L{pygeodesy.Fsum}, precision floating point summation.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import isodd, _xinstanceof, _xsubclassof
from pygeodesy.errors import _xError
from pygeodesy.fmath import _Float, _2float, Fmt, Fsum, hypot2, \
                            _Scalar, sqrt
from pygeodesy.interns import NN, _iadd_, _invalid_, _other_, _SPACE_, \
                             _0_0, _1_5, _2_0, _3_0, _4_0, _6_0
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import _Named, _NotImplemented, Property_RO
# from pygeodesy.proprs import Property_RO  # from .named
# from pygeodesy.streprs import Fmt  # from .fmath

# from math import sqrt  # pow from .fmath

__all__ = _ALL_LAZY.fstats
__version__ = '22.01.16'

_N_3_0 = -_3_0
_N_4_0 = -_4_0


def _2Floats(xs, ys=False):
    '''(INTERNAL) Yield each value as C{float} or L{Fsum}.
    '''
    for i, x in enumerate(xs):
        yield x if isinstance(x, _Float) else (_2float(index=i, ys=x)
                                   if ys else  _2float(index=i, xs=x))


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
        f = self.fcopy(name=self.__add__.__name__)  # PYCHOK expected
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

    def fadd_(self, *xs, **sample):
        '''Accumulate and return the current count.

           @see: Method C{fadd}.
        '''
        return self.fadd(xs, **sample)  # PYCHOK expected

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

    def _iadd_other(self, other):
        '''(INTERNAL) Add Scalar or Scalars.
        '''
        if isinstance(other, _Scalar):
            self.fadd_(other)
        else:
            try:
                if not isinstance(other, (list, tuple)):
                    raise TypeError(_SPACE_(_invalid_, _other_))
                self.fadd(other)  # PYCHOK expected
            except Exception as x:
                raise _xError(x, _SPACE_(self, _iadd_, repr(other)))


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
        self._ms = (_0_0,) * len(self._Ms)  # Moment running values
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

                    d   = other._ms[0] - self._ms[0]  # b1 - a1
                    n   = na + nb
                    n_  = float(n)
                    q1  = d / n_
                    q2  = q1**2  # d**2 / n**2
                    nab = na * nb
                    nq2 = q2 * nab

                    na2 =  na**2
                    nb2 =  nb**2
                    A4 +=  B4
                    A4 += (B3 * na  - (A3 * nb))  * (_4_0 * q1)
                    A4 += (B2 * na2 + (A2 * nb2)) * (_6_0 * q2)
                    A4 += (na2 - nab + nb2) * nq2 * (d * q1)  # d**4 / n**3

                    A3 +=  B3
                    A3 += (A2 * na - (B2 * nb)) * (_3_0 * q1)
                    A3 += (na - nb) * nq2 * d

                    A2 += B2
                    A2 += n_ * nq2  # d**2 * na * nb / n

                    B1n = B1 * nb  # if other is self
                    A1 *= na
                    A1 += B1n
                    A1 *= 1 / n_  # /= chokes PyChecker

#                   self._Ms = A1, A2, A3, A4
                    self._ms = tuple(M.fsum() for M in self._Ms)
                    self._n  = n
                else:
                    self._copy(self, other)
        else:
            self._iadd_other(other)
        return self

    def _copy(self, c, s):
        '''(INTERNAL) Copy C{B{c} = B{s}}.
        '''
        _xinstanceof(Fcook, c=c, s=s)
        c._Ms = tuple(M.fcopy() for M in s._Ms)  # deep=False
        c._ms = tuple(s._ms)
        c._n  =       s._n
        return c

    def fadd(self, xs, sample=False):
        '''Accumulate and return the current count.

           @arg xs: Iterable with additional values (C{Scalar}s,
                    meaning C{scalar} or L{Fsum} instances).
           @kwarg sample: Return the I{sample} instead of the entire
                          I{population} value (C{bool}).

           @return: Current, running (sample) count (C{int}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{xs}} value.

           @raise ValueError: Invalid or non-finite B{C{xs}} value.

           @see: U{online_kurtosis<https://WikiPedia.org/wiki/
                 Algorithms_for_calculating_variance>}.
        '''
        n = self._n
        if xs:
            M1, M2, M3, M4 = (M.fsum_ for M in self._Ms)
            m1, m2, m3, m4 = self._ms
            for x in _2Floats(xs):
                n1 = n
                n += 1
                d  = x - m1  # Fsum or float
                dn = float(d / n)  # == d.fsum()
                if dn:
                    dn2 = dn**2
                    if n1 > 0:
                        d *= dn * n1   # Fsum or float
                        t1 = float(d)  # == d.fsum()
                        t2 = float(d * (dn  * (n - 2)))
                        t3 = float(d * (dn2 * (n**2 - 3 * n1)))
                    else:
                        t1 = t2 = t3 = _0_0
                    m4 = M4(t3, _N_4_0 * dn * m3, _6_0 * dn2 * m2)
                    m3 = M3(t2, _N_3_0 * dn * m2)
                    m2 = M2(t1)
                    m1 = M1(dn)
            self._ms = m1, m2, m3, m4
            self._n  = n
        return _sampled(n, sample)

    def fjb(self, xs=None, sample=True, excess=True):
        '''Accumulate and compute the current U{Jarque-Bera
           <https://WikiPedia.org/wiki/Jarque–Bera_test>} normality.

           @kwarg xs: Iterable with additional values (C{Scalar}s).
           @kwarg sample: Return the I{sample} value (C{bool}), default.
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
                          I{population} value (C{bool}).
           @kwarg excess: Return the I{excess} kurtosis (C{bool}), default.

           @return: Current, running (sample) kurtosis or I{excess} kurtosis (C{float}).

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

    def fkurtosis_(self, *xs, **sample_excess):
        '''Accumulate and return the current kurtosis.

           @see: Method L{Fcook.fkurtosis}.
        '''
        return self.fkurtosis(xs, **sample_excess)

    def fmean(self, xs=None):
        '''Accumulate and return the current mean.

           @kwarg xs: Iterable with additional values (C{Scalar}s).

           @return: Current, running mean (C{float}).

           @see: Method L{Fcook.fadd}.
        '''
        if xs:
            self.fadd(xs)
        return self._ms[0]

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
        m1 = self._ms[0] if xs is None else self.fmean(xs)
        return m1 - self.fskewness() * self.fstdev() / _3_0

    def fmedian_(self, *xs):
        '''Accumulate and return the current median.

           @see: Method L{Fcook.fmedian}.
        '''
        return self.fmedian(xs)

    def fskewness(self, xs=None, sample=False):
        '''Accumulate and return the current skewness.

           @kwarg xs: Iterable with additional values (C{Scalar}s).
           @kwarg sample: Return the I{sample} instead of the entire
                          I{population} value (C{bool}).

           @return: Current, running (sample) skewness (C{float}).

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

    def fstdev(self, xs=None, sample=False):
        '''Accumulate and return the current standard deviation.

           @kwarg xs: Iterable with additional values (C{Scalar}s).
           @kwarg sample: Return the I{sample} instead of the entire
                          I{population} value (C{bool}).

           @return: Current, running (sample) standard deviation (C{float}).

           @see: Method L{Fcook.fadd}.
        '''
        v = self.fvariance(xs, sample=sample)
        return sqrt(v) if v > 0 else _0_0

    def fvariance(self, xs=None, sample=False):
        '''Accumulate and return the current variance.

           @kwarg xs: Iterable with additional values (C{Scalar}s).
           @kwarg sample: Return the I{sample} instead of the entire
                          I{population} value (C{bool}).

           @return: Current, running (sample) variance (C{float}).

           @see: Method L{Fcook.fadd}.
        '''
        n = self.fadd(xs, sample=sample)
        return (self._ms[1] / float(n)) if n > 0 else _0_0

    def toFwelford(self, name=NN):
        '''Return an L{Fwelford} equivalent.
        '''
        f = Fwelford(name=name or self.name)
        f._M = self._Ms[0].fcopy()
        f._m = self._ms[0]
        f._n = self._n
        f._S = self._Ms[1].fcopy()
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
        self._M =  Fsum()  # 1st Moment
        self._m = _0_0     # running mean
        self._S =  Fsum()  # 2nd Moment
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
                    M, S = self._M, self._S

                    n  = na + nb
                    n_ = float(n)

                    d2 = (other._m - self._m)**2
                    S +=  other._S
                    S +=  d2 * (na * nb / n_)

                    Mn = other._M * nb  # if other is self
                    M *= na
                    M += Mn
                    M *= 1 / n_  # /= chokes PyChecker

                    self._m = M.fsum()
                    self._n = n
                else:
                    self._copy(self, other)

        elif isinstance(other, Fcook):
            self += other.toFwelford()
        else:
            self._iadd_other(other)
        return self

    def _copy(self, c, s):
        '''(INTERNAL) Copy C{B{c} = B{s}}.
        '''
        _xinstanceof(Fwelford, c=c, s=s)
        c._M = s._M.fcopy(deep=False)
        c._m = s._m
        c._n = s._n
        c._S = s._S.fcopy(deep=False)
        return c

    def fadd(self, xs, sample=False):
        '''Accumulate and return the current count.

           @arg xs: Iterable with additional values (C{Scalar}s,
                    meaning C{scalar} or L{Fsum} instances).
           @kwarg sample: Return the I{sample} instead of the entire
                          I{population} value (C{bool}).

           @return: Current, running (sample) count (C{int}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{xs}} value.

           @raise ValueError: Invalid or non-finite B{C{xs}} value.
        '''
        n = self._n
        if xs:
            M = self._M.fsum  # allows single Scalar
            m = self._m
            S = self._S
            for x in _2Floats(xs):
                n += 1
                d  = x - m
                m  = M(d / n)
                S += (x - m) * d
            self._m = m
            self._n = n
        return _sampled(n, sample)

    def fmean(self, xs=None):
        '''Accumulate and return the current mean.

           @kwarg xs: Iterable with additional values (C{Scalar}s).

           @return: Current, running mean (C{float}).

           @see: Method L{Fwelford.fadd}.
        '''
        if xs:
            self.fadd(xs)
        return self._m

    def fstdev(self, xs=None, sample=False):
        '''Accumulate and return the current standard deviation.

           @kwarg xs: Iterable with additional values (C{Scalar}s).
           @kwarg sample: Return the I{sample} instead of the entire
                          I{population} value (C{bool}).

           @return: Current, running (sample) standard deviation (C{float}).

           @see: Method L{Fwelford.fadd}.
        '''
        v = self.fvariance(xs, sample=sample)
        return sqrt(v) if v > 0 else _0_0

    def fvariance(self, xs=None, sample=False):
        '''Accumulate and return the current variance.

           @kwarg xs: Iterable with additional values (C{Scalar}s).
           @kwarg sample: Return the I{sample} instead of the entire
                          I{population} value (C{bool}).

           @return: Current, running (sample) variance (C{float}).

           @see: Method L{Fwelford.fadd}.
        '''
        n = self.fadd(xs, sample=sample)
        return (self._S.fsum() / float(n)) if n > 0 else _0_0


class Flinear(_FstatsNamed):
    '''U{Cook<https://www.JohnDCook.com/blog/running_regression>}'s
       C{RunningRegression} computing the running slope, intercept
       and correlation of a linear regression.
    '''
    _mx = _0_0  # running x mean
    _my = _0_0  # running y mean

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
                    n =  other._n  + self._n
                    s = (other._mx - self._mx) * \
                        (other._my - self._my) * \
                        (other._n  * self._n) / float(n)
                    self._S += other._S + s
                    self._X += other._X
                    self._Y += other._Y
                    self._mx = self._X.fmean()
                    self._my = self._Y.fmean()
                    self._n  = n
                else:
                    self._copy(self, other)
        else:
            try:
                if not isinstance(other, (list, tuple)):
                    raise TypeError(_SPACE_(_invalid_, _other_))
                elif isodd(len(other)):
                    raise ValueError(Fmt.PAREN(isodd=Fmt.PAREN(len=_other_)))
                self.fadd_(*other)
            except Exception as x:
                raise _xError(x, _SPACE_(self, _iadd_, repr(other)))
        return self

    def _copy(self, c, s):
        '''(INTERNAL) Copy C{B{c} = B{s}}.
        '''
        _xinstanceof(Flinear, c=c, s=s)
        c._mx = s._mx
        c._my = s._my
        c._n  = s._n
        c._S  = s._S.fcopy(deep=False)
        c._X  = s._X.fcopy(deep=False)
        c._Y  = s._Y.fcopy(deep=False)
        return c

    def fadd(self, xs, ys, sample=False):
        '''Accumulate and return the current count.

           @arg xs: Iterable with additional C{x} values (C{Scalar}s),
                    meaning C{scalar} or L{Fsum} instances).
           @arg ys: Iterable with additional C{y} values (C{Scalar}s,
                    meaning C{scalar} or L{Fsum} instances).
           @kwarg sample: Return the I{sample} instead of the entire
                          I{population} value (C{bool}).

           @return: Current, running (sample) count (C{int}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Non-scalar B{C{xs}} or B{C{ys}} value.

           @raise ValueError: Invalid or non-finite B{C{xs}} or B{C{ys}} value.
        '''
        n = self._n
        if xs and ys:
            mx = self._mx
            my = self._my
            S  = self._S
            X  = self._X.fmean_
            Y  = self._Y.fmean_
            for x, y in zip(_2Floats(xs), _2Floats(ys, ys=True)):
                n1 = n
                n += 1
                if n1 > 0:
                    S += (x - mx) * (y - my) * (n1 / float(n))
                mx, my = X(x), Y(y)
            self._mx = mx
            self._my = my
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
                          I{population} value (C{bool}).
        '''
        return self._sampled(self.x.fstdev(sample=sample) *
                             self.y.fstdev(sample=sample), sample)

    def fintercept(self, sample=False):
        '''Return the current, running (sample) intercept (C{float}).

           @kwarg sample: Return the I{sample} instead of the entire
                          I{population} value (C{bool}).
        '''
        return self._my - self.fslope(sample=sample) * self._mx

    def fslope(self, sample=False):
        '''Return the current, running (sample) slope (C{float}).

           @kwarg sample: Return the I{sample} instead of the entire
                          I{population} value (C{bool}).
        '''
        return self._sampled(self.x.fvariance(sample=sample), sample)

    def _sampled(self, t, sample):
        '''(INTERNAL) Compute the sampled or entire population result.
        '''
        t *= float(_sampled(self._n, sample))
        return (self._S.fsum() / t) if t else _0_0

    @Property_RO
    def x(self):
        '''Get the C{x} accumulator (L{Fcook} or L{Fwelford}).
        '''
        return self._X

    @Property_RO
    def y(self):
        '''Get the C{y} accumulator (L{Fcook} or L{Fwelford}).
        '''
        return self._Y


__all__ += _ALL_DOCS(_FstatsBase, _FstatsNamed)

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
