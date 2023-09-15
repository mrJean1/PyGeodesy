
# -*- coding: utf-8 -*-

u'''Discrete Sine Transforms (AuxDST) in Python, transcoded from I{Karney}'s C++ class
U{DST<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1DST.html>}
in I{GeographicLib version 2.2+}.

Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2022-2023) and licensed
under the MIT/X11 License.  For more information, see the U{GeographicLib
<https://GeographicLib.SourceForge.io>} documentation.

@note: Class L{AuxDST} requires package U{numpy<https://PyPI.org/project/numpy>} to be
       installed, version 1.16 or newer and needed for I{exact} area calculations.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.auxilats.auxily import _Dm
from pygeodesy.basics import isodd, map2, neg, _reverange, _xnumpy
from pygeodesy.constants import PI_2, PI_4, isfinite, _0_0, _0_5, _naninf
# from pygeodesy.fsums import Fsum   # from .karney
from pygeodesy.karney import _2cos2x,  _ALL_DOCS, Fsum, property_RO
# from pygeodesy.lazily import _ALL_DOCS  # from .karney
# from pygeodesy.props import property_RO  # from .karney

__all__ = ()
__version__ = '23.09.14'


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
        '''Evaluate the Fourier sum given the sine and cosine of the angle,
           using precision I{Clenshaw} summation.

           @arg sinx: The sin(I{sigma}) (C{float}).
           @arg cosx: The cos(I{sigma}) (C{float}).
           @arg F: The Fourier coefficients (C{float}[]).
           @arg N: Optional, (smaller) number of terms to evaluate (C{int}).

           @return: Precison I{Clenshaw} sum (C{float}).

           @see: Methods C{AuxDST.integral} and C{AuxDST.integral2}.
        '''
        a = -_2cos2x(cosx, sinx)
        if isfinite(a):
            Y0, Y1 = Fsum(), Fsum()
            n  = _len_N(F, *N)
            Fn =  list(F[:n])
            _F =  Fn.pop
            if isodd(n):
                Y0 -= _F()
            while Fn:  # Y0, Y1 negated
                Y1 -= Y0 * a + _F()
                Y0 -= Y1 * a + _F()
            r = float(_Dm(-Y0, Y1, sinx))
        else:
            r = _naninf(-a)
        return r

    @property_RO
    def _fft_numpy(self):
        '''(INTERNAL) Get the C{numpy.fft} module, I{once}.
        '''
        AuxDST._fft_numpy = fft = _xnumpy(AuxDST, 1, 16).fft  # overwrite property_RO
        return fft

    def _fft_real(self, data):
        '''(INTERNAL) NumPy's I{kissfft}-like C{transform_real} function,
           taking C{float}[:N] B{C{data}} and returning C{complex}[:N*2].
        '''
        # <https://GitHub.com/mborgerding/kissfft/blob/master/test/testkiss.py>
        return self._fft_numpy.rfftn(data)

    def _ffts(self, data, cIV):
        '''(INTERNAL) Compute the DST-III or DST-IV FFTransforms.

           @arg data: Elements DST-III[0:N+1] or DST-IV[0:N] (C{float}[])
                      with DST_III[0] = 0.
           @arg cIV: If C{True} DST-IV, otherwise DST-III.

           @return: FFTransforms (C{float}[0:N]).
        '''
        t, N = (), self.N
        if N > 0:
            N2 = N * 2
            d  = tuple(data)
            # assert len(d) == N + (0 if cIV else 1)

            if cIV:  # DST-IV
                from cmath import exp as _cexp

                def _cF(c, j, r=-PI_4 / N):
                    return c * _cexp(complex(0, r * j))

                i = 0
            else:  # DST-III
                i = 1
                # assert d[0] == _0_0

                def _cF(c, unused):  # PYCHOK redef
                    return c

            d += tuple(reversed(d[i:N]))  # i == len(d) - N
            d += tuple(map(neg, d[:N2]))
            c  = self._fft_real(d)  # complex[0:N*2]
            n2 = float(-N2)
            t  = tuple(_cF(c[j], j).imag / n2 for j in range(1, N2, 2))
        return t

    def _ffts2(self, data, F):
        '''(INTERNAL) Doubled FFTransforms.

           @arg data: Grid centers (C{float}[N]).
           @arg F: The transforms (C{float}[N])

           @return: Doubled FFTransforms (C{float}[N*2]).
        '''
        __2 = _0_5

        def _dmF_2(d, F):
            return (d - F) * __2

        def _dpF_2(d, F):
            return (d + F) * __2

        N = self._N
        # copy DST-IV order N transform to d[0:N]
        d = self._ffts(data, True)
        # assert len(d) >= N and len(F) >= N
        # (DST-IV order N - DST-III order N) / 2
        m = map2(_dmF_2, d[:N], F[:N])
        # (DST-IV order N + DST-III order N) / 2
        p = map2(_dpF_2, d[:N], F[:N])
        return p + tuple(reversed(m))

    @staticmethod
    def integral(sinx, cosx, F, *N):
        '''Evaluate the integral of Fourier sum given the sine and
           cosine of the angle, using precision I{Clenshaw} summation.

           @arg sinx: The sin(I{sigma}) (C{float}).
           @arg cosx: The cos(I{sigma}) (C{float}).
           @arg F: The Fourier coefficients (C{float}[]).
           @arg N: Optional, C{len(B{F})} or a (smaller) number of
                   terms to evaluate (C{int}).

           @return: Precison I{Clenshaw} intergral (C{float}).

           @see: Methods C{AuxDST.evaluate} and C{AuxDST.integral2}.
        '''
        a = _2cos2x(cosx - sinx, cosx + sinx)
        if isfinite(a):
            Y0, Y1 = Fsum(), Fsum()
            for r in _reverscaled(F, *N):
                Y1 -= Y0 * a + r
                Y1, Y0 = Y0, -Y1
            r = float(_Dm(Y1, Y0, cosx))
        else:
            r = _naninf(a)
        return r

    @staticmethod
    def integral2(sinx, cosx, siny, cosy, F, *N):  # PYCHOK no cover
        '''Compute the definite integral of Fourier sum given the
           sine and cosine of the angles at the end points, using
           precision I{Clenshaw} summation.

           @arg sinx: The sin(I{sigma1}) (C{float}).
           @arg cosx: The cos(I{sigma1}) (C{float}).
           @arg siny: The sin(I{sigma2}) (C{float}).
           @arg cosy: The cos(I{sigma2}) (C{float}).
           @arg F: The Fourier coefficients (C{float}[]).
           @arg N: Optional, C{len(B{F})} or a (smaller) number of
                   terms to evaluate (C{int}).

           @return: Precison I{Clenshaw} integral (C{float}).

           @see: Methods C{AuxDST.evaluate} and C{AuxDST.integral}.
        '''
        #  2 * cos(y - x) * cos(y + x) -> 2 * cos(2 * x)
        c =  _2cos2x(cosy * cosx, siny * sinx)
        # -2 * sin(y - x) * sin(y + x) -> 0
        s = -_2cos2x(siny * cosx, cosy * sinx)
        if isfinite(c) and isfinite(s):
            Y0, Y1 = Fsum(), Fsum()
            Z0, Z1 = Fsum(), Fsum()
            for r in _reverscaled(F, *N):
                Y1 -= Y0 * c + Z0 * s + r
                Z1 -= Y0 * s + Z0 * c
                Y1, Y0 = Y0, -Y1
                Z1, Z0 = Z0, -Z1
            r = float(_Dm(Y1, Y0, cosy - cosx) +
                      _Dm(Z1, Z0, cosy + cosx))
        else:
            r = _naninf(c, s)
        return r

    @property_RO
    def N(self):
        '''Get this DST's size, number of points (C{int}).
        '''
        return self._N

    def refine(self, f, F):
        '''Refine the Fourier series by doubling the sampled points.

           @arg f: Single-argument callable (C{B{f}(sigma)}).
           @arg F: The initial Fourier series coefficients (C{float}[:N]).

           @return: Fourier series coefficients (C{float}[:N*2]).
        '''
        def _data(_f, N):  # [:N]
            if N > 0:
                r = PI_4 / N
                for j in range(1, N*2, 2):
                    yield _f(r * j)

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
        '''Determine C{N + 1} terms in the Fourier series.

           @arg f: Single-argument callable (C{B{f}(sigma)}).

           @return: Fourier series coefficients (C{float}[:N+1],
                    leading 0).
        '''
        def _data(_f, N):  # [:N + 1]
            yield _0_0  # data[0] = 0
            if N > 0:
                r = PI_2 / N
                for i in range(1, N + 1):
                    yield _f(r * i)

        return self._ffts(_data(f, self.N), False)


def _len_N(F, *N):
    # Adjusted C{len(B{F})}.
    return min(len(F), *N) if N else len(F)


def _reverscaled(F, *N):
    # Yield F[:N], reversed and scaled
    for n in _reverange(_len_N(F, *N)):
        yield F[n] / float(n * 2 + 1)


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
