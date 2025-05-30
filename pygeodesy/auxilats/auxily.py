
# -*- coding: utf-8 -*-

u'''(INTERNAL) I{Auxiliary} latitudes' classes, constants and functions.

Class L{AuxAngle} transcoded to Python from I{Karney}'s C++ class U{AuxAngle
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AuxAngle.html>}
in I{GeographicLib version 2.2+}.

Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2022-2024) and licensed
under the MIT/X11 License.  For more information, see the U{GeographicLib
<https://GeographicLib.SourceForge.io>} documentation.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # noqa: E702 ;

# from pygeodesy import auxilats  # _MODS
from pygeodesy.auxilats._CX_Rs import _Rkey
from pygeodesy.basics import _isin,  typename
from pygeodesy.constants import INF, NAN, isinf, isnan, _0_0, _0_5, _1_0, \
                               _copysign_1_0, _over, _1_over
from pygeodesy.errors import AuxError
from pygeodesy.fmath import hypot1 as _sc, hypot2_
# from pygeodesy.internals import typename  # from .basics
from pygeodesy.interns import NN,  _DOT_, _UNDER_  # PYCHOK used!
from pygeodesy.lazily import _ALL_DOCS,  _ALL_MODS as _MODS  # PYCHOK used!
from pygeodesy.utily import atan1

from math import asinh, copysign

__all__ = ()
__version__ = '25.05.12'


class Aux(object):
    '''Enum-style Aux names.
    '''
    _coeffs    = {}
    GEOGRAPHIC = PHI   =  GEODETIC = _Rkey(0)
    PARAMETRIC = BETA  =  REDUCED  = _Rkey(1)
    GEOCENTRIC = THETA = _Rkey(2)  # all ...
    RECTIFYING = MU    = _Rkey(3)  # use n^2
    CONFORMAL  = CHI   = _Rkey(4)  # use n
    AUTHALIC   = XI    = _Rkey(5)  # use n
    N  = 6
    N2 = 36

    def __index__(self, aux):
        # throws KeyError, not IndexError
        return _Aux2Greek[aux]

    def __len__(self):
        return Aux.N

    def _CXcoeffs(self, aL):  # in .auxLat.AuxLat._CXcoeffs
        '''(INTERNAL) Get the C{_CX_4._coeffs_4}, C{_CX_6._coeffs_6}
           or C{_CS_8._coeffs_8} coefficients, once.
        '''
        try:
            _coeffs = Aux._coeffs[aL]
        except KeyError:
            try:  # from pygeodesy.auxilats._CX_x import _coeffs_x as _coeffs
                _CX_x   = _DOT_(typename(_MODS.auxilats), _UNDER_('_CX', aL))
                _coeffs = _MODS.getattr(_CX_x, _UNDER_('_coeffs', aL))
            except (AttributeError, ImportError, KeyError, TypeError) as x:
                raise AuxError(ALorder=aL, cause=x)
            Aux._coeffs[aL] = _coeffs._validate(aL, Aux.len(aL))
        return _coeffs

    def _1d(self, auxout, auxin):
        '''Get the 1-d index into N^2 coeffs.
        '''
        N = Aux.N
        if 0 <= auxout < N and 0 <= auxin < N:
            return N * auxout + auxin
        raise AuxError(auxout=auxout, auxin=auxin, N=N)

    def Greek(self, aux):
        '''Get an angle's name (C{str}).
        '''
        return _Aux2Greek.get(aux, NN)

    def len(self, ALorder):  # PYCHOK no cover
        aL = ALorder  # aka Lmax
        mu = Aux.MU * (Aux.MU + 1)
        nu = Aux.N2 -  Aux.N - mu
        return (mu * (aL * (aL + 3) - (aL // 2) * 2) // 4 +
                nu * (aL * (aL + 1)) // 2)

    def power(self, auxout, auxin):
        '''Get the C{convert} exponent (C{int} or C{None}).
        '''
        self._1d(auxout, auxin)  # validate
        return (auxout - auxin) if max(auxin, auxout) < Aux.MU else None

    def use_n2(self, aux):
        return not _isin(aux, Aux.CHI, Aux.XI)

Aux = Aux()  # PYCHOK singleton

_Aux2Greek = {Aux.AUTHALIC:   'Xi',
              Aux.CONFORMAL:  'Chi',
              Aux.GEOCENTRIC: 'Theta',
              Aux.GEODETIC:   'Phi',   # == .GEOGRAPHIC
              Aux.PARAMETRIC: 'Beta',  # == .REDUCED
              Aux.RECTIFYING: 'Mu'}
_Greek2Aux = dict(map(reversed, _Aux2Greek.items()))  # PYCHOK exported
# _Greek2Aux.update((_g.upper(), _x) for _g, _x in _Greek2Aux.items())


def _Dasinh(x, y):
    d = y - x
    if isinf(d):  # PYCHOK no cover
        r = _0_0
    elif isnan(d):  # PYCHOK no cover
        r =  NAN
    elif d:
        xy = x * y
        if xy > 0:
            hx, hy = _sc(x), _sc(y)
            if xy < 1:
                hx, hy = hy, hx
            else:
                x = _1_0 / x
                y = _1_0 / y
            r = _over(x + y, hx * x + hy * y)
            r =  asinh(r * d) / d
        else:
            r = (asinh(y) - asinh(x)) / d
    else:
        r = _1_over(_sc(x))
    return r


def _Datan(x, y):
    xy = x * y
    r = xy + _1_0
    if isnan(r):  # PYCHOK no cover
        pass
    elif x == y:
        r = _1_over(r)
    elif x > 0 and isinf(xy):  # PYCHOK no cover
        r = _0_0
    else:
        d = y - x
        if (r + xy) > 0:
            r =  atan1(d, r) / d  # atan(d / r) / d
        else:
            r = (atan1(y) - atan1(x)) / d
    return r


def _Dh(x, y):
    r = x + y
    if isnan(r):
        pass  # N.B. NAN for inf-inf
    elif isinf(x):  # PYCHOK no cover
        r = copysign(_0_5, x)
    elif isinf(y):  # PYCHOK no cover
        r = copysign(_0_5, y)
    else:
        snx, sny = _sn(x), _sn(y)
        dy = sny * y
        dx = snx * x
        d  = dy + dx
        if (d * _0_5):
            if (x * y) > 0:
                r *= hypot2_(snx / _sc(y), snx * sny,
                             sny / _sc(x)) / (d + d)
            else:
                r = _over((dy - dx) * _0_5, y - x)
        else:  # underflow and x == y == d == 0
            r *= _0_5  # PYCHOK no cover
    return r


def _Dlam(x, y):  # Chi1.tan, Chi2.tan
    # I{Divided difference} of the isometric latitude
    # with respect to the conformal latitude
    if isnan(x) or isnan(y):  # PYCHOK no cover
        r =  NAN
    elif isinf(x) or isinf(y):  # PYCHOK no cover
        r =  INF
    elif x == y:
        r = _sc(x)
    else:
        r = _over(_Dasinh(x, y), _Datan(x, y))
    return r


def _Dm(X, Y, s):  # in .auxDLat, .auxDST
    # Return M{(X - Y) * s}, inplace X
    X -= Y
    X *= s
    return X  # Fsum


def _Dp0Dpsi(x, y):  # Chi1.tan, Chi2.tan
    # I{Divided difference} of the spherical rhumb area
    # term with respect to the isometric latitude
    r = x + y
    if isnan(r):  # PYCHOK no cover
        pass  # NAN for inf-inf
    elif isinf(x):  # PYCHOK no cover
        r = _copysign_1_0(x)
    elif isinf(y):  # PYCHOK no cover
        r = _copysign_1_0(y)
    elif x == y:
        r = _sn(x)
    else:
        r = _Dasinh(_h(x), _h(y))
        r = _over(_Dh(x, y) * r, _Dasinh(x, y))
    return r


def _h(tx):
    '''(INTERNAL) M{h(tan(x)) = tan(x) * sin(x) / 2}
    '''
    if tx:
        tx *= _sn(tx) * _0_5
    return tx  # preserve signed-0


def _sn(tx):
    '''(INTERNAL) M{sin(x) = tan(x) / sqrt(tan(x)**2 + 1)}.
    '''
    if tx:
        tx = _copysign_1_0(tx) if isinf(tx) else (
                           NAN if isnan(tx) else (tx / _sc(tx)))
    return tx  # preserve signed-0


__all__ += _ALL_DOCS(Aux.__class__)
del _Rkey

# **) MIT License
#
# Copyright (C) 2023-2025 -- mrJean1 at Gmail -- All Rights Reserved.
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
