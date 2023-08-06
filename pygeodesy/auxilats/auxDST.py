
# -*- coding: utf-8 -*-

u'''Discrete Sine Transforms (AuxDST) in Python, transcoded from I{Karney}'s C++ class
U{DST<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1DST.html>}
in I{GeographicLib version 2.2+}.

Copyright (C) U{Charles Karney<mailto:Charles@Karney.com>} (2022-2023) and licensed
under the MIT/X11 License.  For more information, see the U{GeographicLib
<https://GeographicLib.SourceForge.io>} documentation.

@note: Class L{AuxDST} requires package U{numpy<https://PyPI.org/project/numpy>} to be
       installed, version 1.16 or newer and needed for I{exact} area calculations.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.auxilats.auxily import _2cos2x
from pygeodesy.basics import isodd, map2, neg, _xnumpy
from pygeodesy.constants import PI_2, PI_4, _0_0, _0_5
from pygeodesy.fsums import Fsum,  Property_RO, property_RO
from pygeodesy.lazily import _ALL_DOCS
# from pygeodesy.props import Property_RO, property_RO  # from .fsums

__all__ = ()
__version__ = '23.08.06'


class AuxDST(object):
    '''Discrete Sine Transforms (DST) for I{Auxiliary Latitudes}.

       @see: I{Karney}'s C++ class U{DST
             <https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1DST.html>}.
    '''
    _N = 0

    def __init__(self, N):
        '''New L{AuxDST} instance.

           @arg N: Size, number of points (C{int}).
        '''
        if N > 0:
            self._N = int(N)
        # kissfft(N, False)  # size, inverse

    @staticmethod
    def evaluate(sinx, cosx, F, *N):
        '''Compute the Fourier sum given the sine and cosine of the angle,
           using I{Clenshaw} summation C{sum(B{F}[i] * sin((2*i+1) * x))}
           for C{i in range(min(len(B{F}), *B{N}))}.

           @arg sinx: The sin(I{sigma}) (C{float}).
           @arg cosx: The cos(I{sigma}) (C{float}).
           @arg F: The Fourier coefficients (C{float}[]).
           @arg N: Optional, (smaller) number of terms to evaluate (C{int}).

           @return: Precison I{Clenshaw} sum (C{float}).

           @see: Methods C{AuxDST.integral} and C{AuxDST.integral2}.
        '''
        n, Y1 = len(F), Fsum()
        if N:
            n = min(n, *N)
        if isodd(n):
            n -= 1
            Y0 = Fsum(-F[n])
        else:
            Y0 = Fsum()
        a = -_2cos2x(cosx, sinx)
        while n > 0:  # Y0, Y1 negated
            n -= 1; Y1 -= Y0 * a + F[n]  # PYCHOK semicolon
            n -= 1; Y0 -= Y1 * a + F[n]  # PYCHOK semicolon
        return -float((Y0 + Y1) * sinx)

    @Property_RO
    def _fft_real(self):
        '''(INTERNAL) Get NumPy's I{kiss}-like C{transform_real},
           taking C{float}[:N] and returning C{complex}[:N*2].
        '''
        # <https://GitHub.com/mborgerding/kissfft/blob/master/test/testkiss.py>
        return _xnumpy(AuxDST, 1, 16).fft.rfftn

    def _ffts(self, data, cIV):
        '''(INTERNAL) Compute the DST-III or DST-IV FFTransforms.

           @arg data: Elements DST-III[1:N+1] or DST-IV[0:N] (C{float}[0:N]).
           @arg cIV: If C{True} DST-IV, otherwise DST-III.

           @return: FFTransforms (C{float}[0:N]).
        '''
        N = self.N
        if N > 0:
            N2 = N * 2
            d  = list(data)
            # assert len(d) == N + (0 if cIV else 1)

            if cIV:  # DST-IV
                from cmath import exp as _cexp

                def _cF(c, j, d=-PI_4 / N):
                    return c * _cexp(complex(0, d * j))

                i = 0
            else:  # DST-III
                i = 1
                # assert d[0] == _0_0

                def _cF(c, unused):  # PYCHOK redef
                    return c

            d += reversed(d[i:N])  # i == len(d) - N
            d += map(neg, d[:N2])
            c  = self._fft_real(d)  # complex[0:N*2]
            n2 = float(-N2)
            d  = tuple(_cF(c[j], j).imag / n2 for j in range(1, N2, 2))
        else:
            d = ()
        return d

    def _ffts2(self, data, F):
        '''(INTERNAL) Doubled FFTransforms.

           @arg data: Grid centers (C{float}[N]).
           @arg F: The transforms (C{float}[N])

           @return: Doubled FFTransforms (C{float}[N*2]).
        '''
        def _dmF(d, F):
            return (d - F) * _0_5

        def _dpF(d, F):
            return (d + F) * _0_5

        N = self._N
        # Copy DST-IV order N transform to d[0:N]
        d = self._ffts(data, True)
        # assert len(d) >= N and len(F) >= N
        # (DST-IV order N - DST-III order N) / 2
        M = map2(_dmF, d[:N], F[:N])
        # (DST-IV order N + DST-III order N) / 2
        P = map2(_dpF, d[:N], F[:N])
        return P + tuple(reversed(M))

    @staticmethod
    def integral(sinx, cosx, F, *N):
        '''Compute the integral of Fourier sum given the sine and cosine
           of the angle using I{Clenshaw} summation C{-sum(B{F}[i] / (2*i+1) *
           cos((2*i+1) * x))} for C{i in range(min(len(B{F}), *B{N}))}.

           @arg sinx: The sin(I{sigma}) (C{float}).
           @arg cosx: The cos(I{sigma}) (C{float}).
           @arg F: The Fourier coefficients (C{float}[]).
           @arg N: Optional, (smaller) number of terms to evaluate (C{int}).

           @return: Precison I{Clenshaw} intergral (C{float}).

           @see: Methods C{AuxDST.evaluate} and C{AuxDST.integral2}.
        '''
        a = _2cos2x(cosx - sinx, cosx + sinx)
        Y0, Y1 = Fsum(), Fsum()
        for r in _reveN(F, *N):
            Y1 -= Y0 * a + r
            Y0, Y1 = -Y1, Y0
        return float((Y1 - Y0) * cosx)

#   @staticmethod
#   def integral2(sin1, cos1, sin2, cos2, F, *N):
#       '''Compute the integral of Fourier sum given the sine and cosine
#          of the angles at the end points using I{Clenshaw} summation
#          C{integral(siny, cosy, F) - integral(sinx, cosx, F)}.
#
#          @arg sin1: The sin(I{sigma1}) (C{float}).
#          @arg cos1: The cos(I{sigma1}) (C{float}).
#          @arg sin2: The sin(I{sigma2}) (C{float}).
#          @arg cos2: The cos(I{sigma2}) (C{float}).
#          @arg F: The Fourier coefficients (C{float}[]).
#          @arg N: Optional, (smaller) number of terms to evaluate (C{int}).
#
#          @return: Precison I{Clenshaw} intergral (C{float}).
#
#          @see: Methods C{AuxDST.evaluate} and C{AuxDST.integral}.
#       '''
#       #  2 * cos(y - x)*cos(y + x) -> 2 * cos(2 * x)
#       a =  _2cos2x(cos2 * cos1, sin2 * sin1)
#       # -2 * sin(y - x)*sin(y + x) -> 0
#       b = -_2cos2x(sin2 * cos1, cos2 * sin1)
#       Y0, Y1 = Fsum(), Fsum()
#       Z0, Z1 = Fsum(), Fsum()
#       for r in _reveN(F, *N):
#           Y1 -= Y0 * a + Z0 * b + r
#           Z1 -= Y0 * b + Z0 * a
#           Y0, Y1 = -Y1, Y0
#           Z0, Z1 = -Z1, Z0
#       return float((Y1 - Y0) * (cos2 - cos1) +
#                    (Z1 - Z0) * (cos2 + cos1))

    @property_RO
    def N(self):
        '''Get this DST's size, number of points (C{int}).
        '''
        return self._N

    def refine(self, f, F):
        '''Double the number of sampled points on a Fourier series.

           @arg f: Single-argument function (C{callable(sigma)} with
                   C{sigma = j * PI_4 / N for j in range(1, N*2, 2)}.
           @arg F: The initial Fourier series coefficients (C{float}[:N]).

           @return: Fourier series coefficients (C{float}[:N*2]).
        '''
        def _data(_f, N):  # [:N]
            if N > 0:
                d = PI_4 / N
                for j in range(1, N*2, 2):
                    yield _f(d * j)

        return self._ffts2(_data(f, self.N), F)

    def reset(self, N):
        '''Reset this DST.

           @arg N: Size, number of points (C{int}).

           @return: The new size (C{int}, non-negative).
        '''
        self._N = N = max(0, N)
        # kissfft.assign(N*2, False)  # "reset" size, inverse
        return N

    def transform(self, f):
        '''Compute C{N} terms in the Fourier series.

           @arg f: Single-argument function (C{callable(sigma)} with
                   C{sigma = i * PI_2 / N for i in range(1, N + 1)}.

           @return: Fourier series coefficients (C{float}[:N]).
        '''
        def _data(_f, N):  # [:N + 1]
            yield _0_0  # data[0] = 0
            if N > 0:
                d = PI_2 / N
                for i in range(1, N + 1):
                    yield _f(d * i)

        return self._ffts(_data(f, self.N), False)


def _reveN(F, *N):
    # Yield F reversed and scaled
    n = len(F)
    if N:
        n = min(n, *N)
    if n > 0:
        n2 = n * 2 + 1
        for r in reversed(F[:n]):
            n2 -= 2
            yield r / n2


__all__ += _ALL_DOCS(AuxDST)

# **) MIT License
#
# Copyright (C) 2023-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
