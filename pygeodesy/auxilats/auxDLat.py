# -*- coding: utf-8 -*-

u'''Class L{AuxDLat} transcoded to Python from I{Karney}'s C++ class U{DAuxLatitude
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1DAuxLatitude.html>}
in I{GeographicLib version 2.2+}.

Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2022-2023) and licensed
under the MIT/X11 License.  For more information, see the U{GeographicLib
<https://GeographicLib.SourceForge.io>} documentation.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.auxilats.auxily import Aux, _Datan, _Dasinh, _Dm, _sc, _sn, \
                                      AuxError
from pygeodesy.auxilats.auxLat import AuxLat,  _ALL_DOCS
from pygeodesy.basics import map1, _reverange
from pygeodesy.constants import INF, NAN, isfinite, isinf, isnan, _0_0, _0_5, \
                               _1_0, _2_0, _N_2_0, _naninf, _over, _1_over
from pygeodesy.elliptic import Elliptic as _Ef,  Fsum
# from pygeodesy.errors import AuxError  # from .auxilats.auxily
# from pygeodesy.fsums import Fsum  # from .elliptic
# from pygeodesy.lazily import _ALL_DOCS  # from .auxilats.auxLat

from math import atan, atan2, cos, sin, sqrt

__all__ = ()
__version__ = '23.09.14'


class AuxDLat(AuxLat):
    '''Class to compute C{Divided Differences} of I{Auxiliary}
       latitudes and other C{Divided Differences} needed for
       L{RhumbAux} and L{RhumbLineAux} calculations.
    '''

    def CParametric(self, Zeta1, Zeta2):
        '''Short for C{.Dconvert(Aux.BETA, B{Zeta1}, B{Zeta2})}.
        '''
        return self.Dconvert(Aux.BETA, Zeta1, Zeta2)

    def CRectifying(self, Zeta1, Zeta2):
        '''Short for C{.Dconvert(Aux.MU, B{Zeta1}, B{Zeta2})}.
        '''
        return self.Dconvert(Aux.MU, Zeta1, Zeta2)

    def _Datanhee(self, x, y):
        # atan(e*sn(tphi))/e:
        #  Datan(e*sn(x),e*sn(y))*Dsn(x,y)/Datan(x,y)
        # asinh(e1*sn(fm1*tphi)):
        #  Dasinh(e1*sn(fm1*x)), e1*sn(fm1*y)) *
        # e1 * Dsn(fm1*x, fm1*y) *fm1 / (e * Datan(x,y))
        # = Dasinh(e1*sn(fm1*x)), e1*sn(fm1*y)) *
        #  Dsn(fm1*x, fm1*y) / Datan(x,y)
        if self.f < 0:
            e  =  self._e
            r  = _Datan(e * _sn(x), e * _sn(y))
        else:
            x *=  self._fm1
            y *=  self._fm1
            e1 =  self._e1
            r  = _Dasinh(e1 * _sn(x), e1 * _sn(y))
        return _Dsn(x, y) * r

    def Dconvert(self, auxout, Zeta1, Zeta2):
        '''I{Divided Difference} of one auxiliary latitude wrt another.
        '''
        auxin = Zeta1._AUX
        # assert Zeta2._AUX == auxin
        try:
            if auxout != auxin:
                cs =  self._coeffs(auxout, auxin)
                # assert len(cs) == self.ALorder
                r  = _DClenshaw(True, Zeta1, Zeta2, cs, self.ALorder)
            else:
                r  = _1_0
        except AuxError:  # no _coeffs
            r = NAN
        return r

    def DE(self, X, Y):
        # We assume that X and Y are in [-90d, 90d] and
        # have the same sign. If not we would include
        #    if (Xn.y() * Yn.y() < 0)
        #      return d != 0 ? (E(X) - E(Y)) / d : 1
        # The general formula fails for x = y = 0d and
        # x = y = 90d.  Probably this is fixable (the
        # formula works for other x = y.  But let's
        # also stipulate that x != y.

        # Make both y positive, so we can do the swap a <-> b trick
        sx, cx, x =  X._yxr_normalized(True)
        sy, cy, y =  Y._yxr_normalized(True)
        k2,     d = -self._e12, (y - x)
        # Switch prolate to oblate, then use formulas for k2 < 0
        if self.f < 0:  # XXX and False?
            sx, cx = cx, sx
            sy, cy = cy, sy
            d,  k2 = -d, self._e2
        # See DLMF: Eqs (19.11.2) and (19.11.4) letting
        Dt = _Dsin(x, y) * (sx + sy)
        if Dt:
            t   = _sxk2y(sx, sy, k2) + _sxk2y(sy, sx, k2)
            Dt  = _over(Dt, t * (cx + cy))
            t   =  d * Dt
            t2  = _1_0 + t**2
            Dt *= _2_0 / t2
            sk2 = (d * Dt)**2 * k2
            d2  = _1_0 - sk2
            c2  = ((_1_0 - t) * (_1_0 + t) / t2)**2 if t else _1_0
            # E(z)/sin(z)
            Dt *= _Ef._RFRD(c2, d2, _1_0, sk2) - k2 * sx * sy
        return Dt

    def DIsometric(self, Phi1, Phi2):
        '''I{Divided Difference} of the isometric wrt the geographic latitude.
        '''
        tx, ty = Phi1.tan, Phi2.tan
        if isnan(ty) or isnan(tx):  # PYCHOK no cover
            r =  NAN
        elif isinf(ty) or isinf(tx):  # PYCHOK no cover
            r =  INF
        else:  # psi = asinh(tan(Phi)) - e^2 * atanhee(tan(Phi))
            r =  self._Datanhee(tx, ty) * self._e2
            r = _over(_Dasinh(tx, ty) - r, _Datan(tx, ty))
        return r

    def DParametric(self, Phi1, Phi2):
        '''I{Divided Difference} of the parametric wrt the geographic latitude.
        '''
        fm1, e2m1 = self._fm1, self._e2m1
        tx,  ty   = Phi1.tan,  Phi2.tan
        # DbetaDphi = Datan(fm1*tx, fm1*ty) * fm1 / Datan(tx, ty)
        # Datan(x, y) = 1 / (1 + x^2),                   if x == y
        #             = (atan(y) - atan(x)) / (y-x),     if x*y < 0
        #             = atan((y-x) / (1 + x*y)) / (y-x), if x*y > 0
        txy = tx * ty
        if txy < 0 or (isinf(ty) and not tx):
            _a =  atan
            r  = _over(_a(fm1 * ty) - _a(fm1 * tx), _a(ty) - _a(tx))
        elif tx == ty:  # includes tx = ty = inf
            if txy > 1:  # == tx**2
                txy = _1_over(txy)
                r   =  txy + e2m1
            else:
                r   =  txy * e2m1 + _1_0
            r = _over(fm1 * (txy + _1_0), r)
        else:
            if txy > 1:
                tx  = _1_over(tx)
                ty  = _1_over(ty)
                txy =  tx * ty
                t   =  txy + e2m1
            else:
                t   =  txy * e2m1 + _1_0
            r =  ty - tx
            r = _over(atan2(r * fm1, t), atan2(r, _1_0 + txy))
        return r

    def DRectifying(self, Phi1, Phi2):
        '''I{Divided Difference} of the rectifying wrt the geographic latitude.
        '''
        # Stipulate that Phi1 and Phi2 are in [-90d, 90d]
        x, y = Phi1.toRadians, Phi2.toRadians
        if y == x:  # isnear0
            Mu1 = self.Rectifying(Phi1, diff=True)
            tphi1, r = Phi1.tan, Mu1.diff
            if isfinite(tphi1):
                r *= _over(_sc(tphi1), _sc(Mu1.tan))**2
            else:  # PYCHOK no cover
                r  = _1_over(r)
        elif (x * y) < 0:
            r = _over(self.Rectifying(Phi2).toRadians -
                      self.Rectifying(Phi1).toRadians, y - x)
        else:
            r  = _over(self.b, self.RectifyingRadius(True))
            r *= self.DE(*map1(self.Parametric, Phi1, Phi2))
            r *= self.DParametric(Phi1, Phi2)
        return r  # or INF or NAN


def _DClenshaw(sinp, Zeta1, Zeta2, cs, K):
    '''(INTERNAL) I{Divided Difference} of L{AuxLat._Clenshaw}.

        @return: C{Fsum} if B{C{sinp}} otherwise a C{float}.
    '''
    s1, c1, r1 = Zeta1._yxr_normalized(False)
    s2, c2, r2 = Zeta2._yxr_normalized(False)
    Delta = r2 - r1
    # Evaluate (Clenshaw(sinp, szeta2, czeta2, cs, K) -
    #           Clenshaw(sinp, szeta1, czeta1, cs, K)) / Delta
    #        or f = sin if sinp else cos
    #           sum(cs[k] * (f((2*k+2) * Zeta2) -
    #                        f((2*k+2) * Zeta2))) / Delta
    #
    # Delta is EITHER 1, giving the plain difference OR (Zeta2 - Zeta1)
    # in radians, giving the I{Divided Difference}.  Other values will
    # produce nonsense.
    #
    # Suffices a and b denote [1,1], [2,1] elements of matrix/vector
    cp = cm = c2 * c1
    t       = s2 * s1
    cp -= t  # not +
    cm += t  # not -

    sp = s2 * c1
    t  = c2 * s1
    smd = ((sin(Delta) / Delta) if Delta != _1_0 else
           (sp - t))  if Delta  else _1_0
    sp += t

    xa = cp * cm  * _2_0
    xb = sp * smd * _N_2_0
    xD = xb * Delta**2

    if isfinite(xD) and isfinite(xb) and isfinite(xa):
        U0a, U1a = Fsum(), Fsum()
        U0b, U1b = Fsum(), Fsum()
        for k in _reverange(K):  # assert len(cs) == K
            # t = x . U0 - U1 + cs[k] * I
            U1a -= U0a * xa + U0b * xD + cs[k]
            U1b -= U0a * xb + U0b * xa
            U1a, U0a = U0a, -U1a
            U1b, U0b = U0b, -U1b
        # F0a  = (sp if sinp else  cp) * cm
        # F0b  = (cp if sinp else -sp) * smd
        # Fm1a =   0 if sinp else   1  # Fm1b = 0
        # return (U0b * F0a + U0a * F0b - U1b * Fm1a) * 2
        if sinp:
            U1b = _0_0
        else:
            sp, cp = cp, -sp
        U0b *=  sp * cm
        U0a *=  cp * smd
        U0a +=  U0b
        U0a  = _Dm(U0a, U1b, _2_0)
        r = float(U0a) if sinp else U0a  # Fsum
    else:
        r = _naninf(xD, xb, xa)
    return r


def _Dsin(x, y):  # see also .rhumbx._Dsin
    r = cos((x + y) * _0_5)
    d =     (x - y) * _0_5
    if d:
        r *= sin(d) / d
    return r


def _Dsn(x, y):
    # (sn(y) - sn(x)) / (y - x)
    if x != y:
        snx, sny = _sn(x), _sn(y)
        if (x * y) > 0:
            scx, scy = _sc(x), _sc(y)
            r = _over((snx / scy) + (sny / scx),
                      (snx + sny) *  scy * scx)
        else:
            r = (sny - snx) / (y - x)
    elif x:
        r = _1_over(_sc(x) * (x**2 + _1_0))  # == 1 / sqrt3(x**2 + 1)
    else:
        r = _1_0
    return r


def _sxk2y(sx, sy, k2):
    # .DE helper
    sy *= sy * k2
    if sy:
        try:
            sx *= sqrt(_1_0 - sy)
        except ValueError:  # domain error
            sx  = NAN
    return sx


__all__ += _ALL_DOCS(AuxDLat)

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
