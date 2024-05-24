
# -*- coding: utf-8 -*-

u'''Classes for I{running} statistics and regressions based on
L{pygeodesy.Fsum}, precision floating point summation.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import isscalar, isodd, _xinstanceof, \
                            _xiterable, _xsubclassof, _zip
from pygeodesy.constants import _0_0, _1_0, _2_0, _3_0, _4_0, _6_0
from pygeodesy.errors import _AssertionError, _ValueError, _xError
from pygeodesy.fmath import Fsqrt,  Fmt
from pygeodesy.fsums import _2finite, Fsum, _iadd_op_, _isFsumTuple
from pygeodesy.interns import _odd_, _SPACE_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import _name__, _Named, _NotImplemented, \
                             property_RO
# from pygeodesy.props import property_RO  # from .named
# from pygeodesy.streprs import Fmt  # from .fmath

__all__ = _ALL_LAZY.fstats
__version__ = '24.05.21'


def _2Floats(**xs):
    '''(INTERNAL) Yield each value as C{float} or L{Fsum}.
    '''
    try:
        name, xs = xs.popitem()
    except Exception as X:
        raise _AssertionError(xs=xs, cause=X)

    try:
        i = None
        for i, x in enumerate(_xiterable(xs)):  # don't unravel Fsums
            yield x._Fsum if _isFsumTuple(x) else _2finite(x)
    except Exception as X:
        raise _xError(X, name, xs) if i is None else \
              _xError(X, Fmt.INDEX(name, i), x)


def _sampled(n, sample):
    '''(INTERNAL) Return the sample or the entire count.
    '''
    return (n - 1) if sample and n > 0 else n


class _FstatsNamed(_Named):
    '''(INTERNAL) Base class.
    '''
    _n = 0

    def __add__(self, other):
        '''Sum of this and an other instance or a C{scalar} or an
           L{Fsum}, L{Fsum2Tuple} or
           .
        '''
        f  = self.copy(name=self.__add__.__name__)  # PYCHOK expected
        f += other
        return f

    def __float__(self):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self)

    def __int__(self):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self)

    def __len__(self):
        '''Return the I{total} number of accumulated C{Scalars} (C{int}).
        '''
        return self._n

    def __neg__(self):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self)

    def __radd__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)

    def __str__(self):
        n =  self.name
        n = _SPACE_(self.classname, n) if n else self.classname
        return Fmt.SQUARE(n, len(self))

    def copy(self, deep=False, **name):
        '''Copy this instance, C{shallow} or B{C{deep}}.

           @kwarg name: Optional, overriding C{B{name}="copy"} (C{str}).

           @return: The copy instance.
        '''
        n = _name__(name, name__=self.copy)
        f = _Named.copy(self, deep=deep, name=n)
        return self._copy(f, self)  # PYCHOK expected

    fcopy = copy  # for backward compatibility


class _FstatsBase(_FstatsNamed):
    '''(INTERNAL) Base running stats class.
    '''
    _Ms = ()

    def _copy(self, d, s):
        '''(INTERNAL) Copy C{B{c} = B{s}}.
        '''
        _xinstanceof(self.__class__, d=d, s=s)
        d._Ms = tuple(M.copy() for M in s._Ms)  # deep=False
        d._n  =                         s._n
        return d

    def fadd(self, xs, sample=False):  # PYCHOK no cover
        '''I{Must be overloaded}.'''
        self._notOverloaded(xs, sample=sample)

    def fadd_(self, *xs, **sample):
        '''Accumulate and return the current count.

           @see: Method C{fadd} for further details.
        '''
        return self.fadd(xs, **sample)

    def fmean(self, xs=None):
        '''Accumulate and return the current mean.

           @kwarg xs: Iterable of additional values (each C{scalar} or
                      an L{Fsum} or L{Fsum2Tuple} instance).

           @return: Current, running mean (C{float}).

           @see: Method C{fadd}.
        '''
        return float(self._Mean(xs))

    def fmean_(self, *xs):
        '''Accumulate and return the current mean.

           @see: Method C{fmean} for further details.
        '''
        return self.fmean(xs)

    def fstdev(self, xs=None, **sample):
        '''Accumulate and return the current standard deviation.

           @arg xs: Iterable of additional values (each C{scalar} or an
                    L{Fsum} or L{Fsum2Tuple} instance).
           @kwarg sample: Use C{B{sample}=True} for the I{sample} deviation
                          instead of the I{population} deviation (C{bool}).

           @return: Current, running (sample) standard deviation (C{float}).

           @see: Method C{fadd}.
        '''
        return float(self._Stdev(xs, **sample))

    def fstdev_(self, *xs, **sample):
        '''Accumulate and return the current standard deviation.

           @see: Method C{fstdev} for further details.
        '''
        return self.fstdev(xs, **sample)

    def fvariance(self, xs=None, **sample):
        '''Accumulate and return the current variance.

           @arg xs: Iterable of additional values (each C{scalar} or an
                    L{Fsum} or L{Fsum2Tuple} instance).
           @kwarg sample: Use C{B{sample}=True} for the I{sample} variance
                          instead of the I{population} variance (C{bool}).

           @return: Current, running (sample) variance (C{float}).

           @see: Method C{fadd}.
        '''
        return float(self._Variance(xs, **sample))

    def fvariance_(self, *xs, **sample):
        '''Accumulate and return the current variance.

           @see: Method C{fvariance} for further details.
        '''
        return self.fvariance(xs, **sample)

    def _iadd_other(self, other):
        '''(INTERNAL) Add one or several values.
        '''
        try:
            if _isFsumTuple(other):
                self.fadd_(other._Fsum)
            elif isscalar(other):
                self.fadd_(_2finite(other))
            elif _xiterable(other):
                self.fadd(other)
        except Exception as X:
            t = _SPACE_(self, _iadd_op_, repr(other))
            raise _xError(X, t)

    @property_RO
    def _M1(self):
        '''(INTERNAL) get the 1st Moment accumulator.'''
        return self._Ms[0]

    @property_RO
    def _M2(self):
        '''(INTERNAL) get the 2nd Moment accumulator.'''
        return self._Ms[1]

    def _Mean(self, xs=None):
        '''(INTERNAL) Return the current mean as L{Fsum}.
        '''
        if xs:
            self.fadd(xs)
        return self._M1  # .copy()

    def _Stdev(self, xs=None, **sample):
        '''(INTERNAL) Return the current (sample) standard deviation as L{Fsum}.
        '''
        V = self._Variance(xs, **sample)
        return Fsqrt(V) if V > 0 else _0_0

    def _Variance(self, xs=None, **sample):
        '''(INTERNAL) Return the current (sample) variance as L{Fsum}.
        '''
        n = self.fadd(xs, **sample)
        return (self._M2 / n) if n > 0 else _0_0


class Fcook(_FstatsBase):
    '''U{Cook<https://www.JohnDCook.com/blog/skewness_kurtosis>}'s
       C{RunningStats} computing the running mean, median and
       (sample) kurtosis, skewness, variance, standard deviation
       and Jarque-Bera normality.

       @see: L{Fwelford} and U{Higher-order statistics<https://
             WikiPedia.org/wiki/Algorithms_for_calculating_variance>}.
    '''
    def __init__(self, xs=None, **name):
        '''New L{Fcook} stats accumulator.

           @arg xs: Iterable of additional values (each C{scalar} or
                    an L{Fsum} or L{Fsum2Tuple} instance).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @see: Method L{Fcook.fadd}.
        '''
        self._Ms = tuple(Fsum() for _ in range(4))  # 1st, 2nd ... Moment
        if name:
            self.name = name
        if xs:
            self.fadd(xs)

    def __iadd__(self, other):
        '''Add B{C{other}} to this L{Fcook} instance.

           @arg other: An L{Fcook} instance or value or iterable
                       of values (each C{scalar} or an L{Fsum}
                       or L{Fsum2Tuple} instance).

           @return: This instance, updated (L{Fcook}).

           @raise TypeError: Invalid B{C{other}}.

           @raise ValueError: Invalid or non-finite B{C{other}}.

           @see: Method L{Fcook.fadd}.
        '''
        if isinstance(other, Fcook):
            nb = len(other)
            if nb > 0:
                na = len(self)
                if na > 0:
                    A1, A2, A3, A4 = self._Ms
                    B1, B2, B3, B4 = other._Ms

                    n   =  na + nb
                    _n  = _1_0 / n
                    D   =  A1 - B1  # b1 - a1
                    Dn  =  D * _n
                    Dn2 =  Dn**2  # d**2 / n**2
                    nab =  na * nb
                    Dn3 =  Dn2 * (D * nab)

                    na2 =  na**2
                    nb2 =  nb**2
                    A4 +=  B4
                    A4 += (B3 * na  - (A3 * nb))  * (Dn  * _4_0)
                    A4 += (B2 * na2 + (A2 * nb2)) * (Dn2 * _6_0)
                    A4 += (Dn * Dn3) * (na2 - nab + nb2)  # d**4 / n**3

                    A3 +=  B3
                    A3 += (A2 * na - (B2 * nb)) * (Dn * _3_0)
                    A3 +=  Dn3 * (na - nb)

                    A2 += B2
                    A2 += Dn2 * (nab * _n)

                    B1n = B1 * nb  # if other is self
                    A1 *= na
                    A1 += B1n
                    A1 *= _n

#                   self._Ms = A1, A2, A3, A4
                    self._n  = n
                else:
                    self._copy(self, other)
        else:
            self._iadd_other(other)
        return self

    def fadd(self, xs, sample=False):
        '''Accumulate and return the current count.

           @arg xs: Iterable of additional values (each C{scalar} or an
                    L{Fsum} or L{Fsum2Tuple} instance).
           @kwarg sample: Use C{B{sample}=True} for the I{sample} count
                          instead of the I{population} count (C{bool}).

           @return: Current, running (sample) count (C{int}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Invalid B{C{xs}}.

           @raise ValueError: Invalid or non-finite B{C{xs}}.

           @see: U{online_kurtosis<https://WikiPedia.org/wiki/
                 Algorithms_for_calculating_variance>}.
        '''
        n = self._n
        if xs:
            M1, M2, M3, M4 = self._Ms
            for x in _2Floats(xs=xs):  # PYCHOK yield
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

    def fjb(self, xs=None, excess=True, sample=True):
        '''Accumulate and compute the current U{Jarque-Bera
           <https://WikiPedia.org/wiki/Jarque–Bera_test>} normality.

           @kwarg xs: Iterable of additional values (each C{scalar} or an
                      L{Fsum} or L{Fsum2Tuple}).
           @kwarg excess: Apply the I{excess} kurtosis (C{bool}), default.
           @kwarg sample: Use C{B{sample}=False} for the I{population}
                          normality instead of the I{sample} one (C{bool}).

           @return: Current, running (sample) Jarque-Bera normality (C{float}).

           @see: Method L{Fcook.fadd}.
        '''
        return float(self._JarqueBera(xs, excess, sample=sample))

    def fjb_(self, *xs, **sample_excess):
        '''Accumulate and compute the current U{Jarque-Bera
           <https://WikiPedia.org/wiki/Jarque–Bera_test>} normality.

           @see: Method L{Fcook.fjb} for further details.
        '''
        return self.fjb(xs, **sample_excess)

    def fkurtosis(self, xs=None, excess=True, **sample):
        '''Accumulate and return the current kurtosis.

           @arg xs: Iterable of additional values (each C{scalar} or an
                    L{Fsum} or L{Fsum2Tuple} instance).
           @kwarg excess: Return the I{excess} kurtosis (C{bool}), default.
           @kwarg sample: Use C{B{sample}=True} for the I{sample} kurtosis
                          instead of the I{population} kurtosis (C{bool}).

           @return: Current, running (sample) kurtosis or I{excess} kurtosis (C{float}).

           @see: U{Kurtosis Formula<https://www.Macroption.com/kurtosis-formula>}
                 and U{Mantalos<https://www.ResearchGate.net/publication/227440210>}.

           @see: Method L{Fcook.fadd}.
        '''
        n = self.fadd(xs, **sample)
        return float(self._Kurtosis(n, excess, **sample))

    def fkurtosis_(self, *xs, **excess_sample):
        '''Accumulate and return the current kurtosis.

           @see: Method L{Fcook.fkurtosis} for further details.
        '''
        return self.fkurtosis(xs, **excess_sample)

    def fmedian(self, xs=None):
        '''Accumulate and return the current median.

           @arg xs: Iterable of additional values (each C{scalar} or an
                    L{Fsum} or L{Fsum2Tuple} instance).

           @return: Current, running median (C{float}).

           @see: U{Pearson's Skewness Coefficients<https://MathWorld.Wolfram.com/
                 PearsonsSkewnessCoefficients.html>}, U{Skewness & Kurtosis Simplified
                 https://TowardsDataScience.com/skewness-kurtosis-simplified-1338e094fc85>}
                 and method L{Fcook.fadd}.
        '''
        return float(self._Median(xs))

    def fmedian_(self, *xs):
        '''Accumulate and return the current median.

           @see: Method L{Fcook.fmedian} for further details.
        '''
        return self.fmedian(xs)

    def fskewness(self, xs=None, **sample):
        '''Accumulate and return the current skewness.

           @arg xs: Iterable of additional values (each C{scalar} or an
                    L{Fsum} or L{Fsum2Tuple} instance).
           @kwarg sample: Use C{B{sample}=True} for the I{sample} skewness
                          instead of the I{population} skewness (C{bool}).

           @return: Current, running (sample) skewness (C{float}).

           @see: U{Skewness Formula<https://www.Macroption.com/skewness-formula/>}
                 and U{Mantalos<https://www.ResearchGate.net/publication/227440210>}.

           @see: Method L{Fcook.fadd}.
        '''
        n = self.fadd(xs, **sample)
        return float(self._Skewness(n, **sample))

    def fskewness_(self, *xs, **sample):
        '''Accumulate and return the current skewness.

           @see: Method L{Fcook.fskewness} for further details.
        '''
        return self.fskewness(xs, **sample)

    def _JarqueBera(self, xs, excess, **sample):
        '''(INTERNAL) Return the (sample) Jarque-Bera normality as L{Fsum}.
        '''
        N, n = _0_0, self.fadd(xs, **sample)
        if n > 0:
            K = self._Kurtosis(n, excess, **sample) / _2_0
            S = self._Skewness(n,         **sample)
            N = (K**2 + S**2) * (n / _6_0)  # Fpowers(2, K, S) * ...
        return N

    def _Kurtosis(self, n, excess, sample=False):
        '''(INTERNAL) Return the (sample) kurtosis as L{Fsum} or C{0.0}.
        '''
        K = _0_0
        if n > 0:
            _, M2, _, M4 = self._Ms
            M = M2**2
            if M > 0:
                K, x = M.rdiv(M4 * n, raiser=False), _3_0
                if sample and 2 < n < len(self):
                    d  = (n - 1) * (n - 2)
                    K *= (n + 1) * (n + 2) / d
                    x *=  n**2 / d
                if excess:
                    K -= x
        return K

    def _Median(self, xs=None):
        '''(INTERNAL) Return the median as L{Fsum}.
        '''
        # skewness = 3 * (mean - median) / stdev, i.e.
        # median = mean - (skewness * stdef) / 3
        return self._Mean(xs) - (self._Skewness(self._n) *
                                 self._Stdev()) / _3_0

    def _Skewness(self, n, sample=False):
        '''(INTERNAL) Return the (sample) skewness as L{Fsum} or C{0.0}.
        '''
        S = _0_0
        if n > 0:
            _, M2, M3, _ = self._Ms
            M = M2**3
            if M > 0:
                M = M.rdiv(n, raiser=False)
                S = M3 * Fsqrt(M, raiser=False)
                if sample and 1 < n < len(self):
                    S *= (n + 1) / (n - 1)
        return S

    def toFwelford(self, **name):
        '''Return a L{Fwelford} equivalent.

           @kwarg name: Optional C{B{name}=NN} (C{str}).
        '''
        f = Fwelford(name=self._name__(name))
        f._Ms = self._M1.copy(), self._M2.copy()  # deep=False
        f._n  = self._n
        return f


class Fwelford(_FstatsBase):
    '''U{Welford<https://WikiPedia.org/wiki/Algorithms_for_calculating_variance>}'s
       accumulator computing the running mean, (sample) variance and standard deviation.

       @see: U{Cook<https://www.JohnDCook.com/blog/standard_deviation/>} and L{Fcook}.
    '''
    def __init__(self, xs=None, **name):
        '''New L{Fwelford} stats accumulator.

           @arg xs: Iterable of initial values (each C{scalar} or an
                    L{Fsum} or L{Fsum2Tuple} instance).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @see: Method L{Fwelford.fadd}.
        '''
        self._Ms = Fsum(), Fsum()  # 1st and 2nd Moment
        if name:
            self.name = name
        if xs:
            self.fadd(xs)

    def __iadd__(self, other):
        '''Add B{C{other}} to this L{Fwelford} instance.

           @arg other: An L{Fwelford} or L{Fcook} instance or value
                       or an iterable of values (each C{scalar} or
                       an L{Fsum} or L{Fsum2Tuple} instance).

           @return: This instance, updated (L{Fwelford}).

           @raise TypeError: Invalid B{C{other}}.

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

                    n  =  na + nb
                    _n = _1_0 / n

                    D  = M_ - M
                    D *= D  # D**2
                    D *= na * nb * _n
                    S += D
                    S += S_

                    Mn = M_ * nb  # if other is self
                    M *= na
                    M += Mn
                    M *= _n

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

           @arg xs: Iterable of additional values (each C{scalar} or an
                    L{Fsum} or L{Fsum2Tuple} instance).
           @kwarg sample: Use C{B{sample}=True} for the I{sample} count
                          instead of the I{population} count (C{bool}).

           @return: Current, running (sample) count (C{int}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Invalid B{C{xs}}.

           @raise ValueError: Invalid or non-finite B{C{xs}}.
        '''
        n = self._n
        if xs:
            M, S = self._Ms
            for x in _2Floats(xs=xs):  # PYCHOK yield
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
    def __init__(self, xs=None, ys=None, Fstats=Fwelford, **name):
        '''New L{Flinear} regression accumulator.

           @kwarg xs: Iterable of initial C{x} values (each C{scalar} or
                      an L{Fsum} or L{Fsum2Tuple} instance).
           @kwarg ys: Iterable of initial C{y} values (each C{scalar} or
                      an L{Fsum} or L{Fsum2Tuple} instance).
           @kwarg Fstats: Class for C{xs} and C{ys} values (L{Fcook} or
                          L{Fwelford}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @raise TypeError: B{C{Fstats}} not L{Fcook} or L{Fwelford}.

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

           @arg other: An L{Flinear} instance or an iterable of
                       C{x_ys} values, see method C{fadd_}.

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
                    n =  other._n
                    D = (other._X._M1 - self._X._M1) * \
                        (other._Y._M1 - self._Y._M1) * \
                        (n * self._n / (self._n + n))
                    self._S += other._S + D
                    self._X += other._X
                    self._Y += other._Y
                    self._n += n
                else:
                    self._copy(self, other)
        else:
            try:
                if _xiterable(other):
                    self.fadd_(*other)
            except Exception as X:
                op = _SPACE_(self, _iadd_op_, repr(other))
                raise _xError(X, op)
        return self

    def _copy(self, d, s):
        '''(INTERNAL) Copy C{B{d} = B{s}}.
        '''
        _xinstanceof(Flinear, d=d, s=s)
        d._S = s._S.copy(deep=False)
        d._X = s._X.copy(deep=False)
        d._Y = s._Y.copy(deep=False)
        d._n = s._n
        return d

    def _Correlation(self, **sample):
        '''(INTERNAL) Return the current (sample) correlation as L{Fsum}.
        '''
        return self._Sampled(self._X._Stdev(**sample) *
                             self._Y._Stdev(**sample), **sample)

    def fadd(self, xs, ys, sample=False):
        '''Accumulate and return the current count.

           @arg xs: Iterable of additional C{x} values (each C{scalar}
                    or an L{Fsum} or L{Fsum2Tuple} instance).
           @arg ys: Iterable of additional C{y} values (each C{scalar}
                    or an L{Fsum} or L{Fsum2Tuple} instance).
           @kwarg sample: Use C{B{sample}=True} for the I{sample} count
                          instead of the I{population} count (C{bool}).

           @return: Current, running (sample) count (C{int}).

           @raise OverflowError: Partial C{2sum} overflow.

           @raise TypeError: Invalid B{C{xs}} or B{C{ys}}.

           @raise ValueError: Invalid or non-finite B{C{xs}} or B{C{ys}}.
        '''
        n = self._n
        if xs and ys:
            S = self._S
            X = self._X
            Y = self._Y
            for x, y in _zip(_2Floats(xs=xs), _2Floats(ys=ys)):  # PYCHOK strict=True
                n1 = n
                n += 1
                if n1 > 0:
                    S += (X._M1 - x) * (Y._M1 - y) * (n1 / n)
                X += x
                Y += y
            self._n = n
        return _sampled(n, sample)

    def fadd_(self, *x_ys, **sample):
        '''Accumulate and return the current count.

           @arg x_ys: Individual, alternating C{x, y, x, y, ...} values
                      (each C{scalar} or an L{Fsum} or L{Fsum2Tuple}
                      instance).

           @see: Method C{Flinear.fadd} for further details.
        '''
        if isodd(len(x_ys)):
            t = _SPACE_(_odd_, len.__name__)
            raise _ValueError(t, len(x_ys))
        return self.fadd(x_ys[0::2], x_ys[1::2], **sample)

    def fcorrelation(self, **sample):
        '''Return the current, running (sample) correlation (C{float}).

           @kwarg sample: Use C{B{sample}=True} for the I{sample} correlation
                          instead of the I{population} correlation (C{bool}).
        '''
        return float(self._Correlation(**sample))

    def fintercept(self, **sample):
        '''Return the current, running (sample) intercept (C{float}).

           @kwarg sample: Use C{B{sample}=True} for the I{sample} intercept
                          instead of the I{population} intercept (C{bool}).
        '''
        return float(self._Intercept(**sample))

    def fslope(self, **sample):
        '''Return the current, running (sample) slope (C{float}).

           @kwarg sample: Use C{B{sample}=True} for the I{sample} slope
                          instead of the I{population} slope (C{bool}).
        '''
        return float(self._Slope(**sample))

    def _Intercept(self, **sample):
        '''(INTERNAL) Return the current (sample) intercept as L{Fsum}.
        '''
        return self._Y._M1 - self._X._M1 * self._Slope(**sample)

    def _Sampled(self, T, sample=False):
        '''(INTERNAL) Compute the sampled or entire population result.
        '''
        T *= _sampled(self._n, sample)
        return self._S.copy().fdiv(T, raiser=False) if T else T

    def _Slope(self, **sample):
        '''(INTERNAL) Return the current (sample) slope as L{Fsum}.
        '''
        return self._Sampled(self._X._Variance(**sample), **sample)

    @property_RO
    def x(self):
        '''Get the C{x} accumulator (L{Fcook} or L{Fwelford}).
        '''
        return self._X  # .copy()

    @property_RO
    def y(self):
        '''Get the C{y} accumulator (L{Fcook} or L{Fwelford}).
        '''
        return self._Y  # .copy()


__all__ += _ALL_DOCS(_FstatsBase, _FstatsNamed)

# **) MIT License
#
# Copyright (C) 2021-2024 -- mrJean1 at Gmail -- All Rights Reserved.
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
