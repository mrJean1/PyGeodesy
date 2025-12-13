
# -*- coding: utf-8 -*-

u'''I{Jacobi Conformal projection} classes L{Conformal3}, L{Conformal3B} and L{Conformal3Sphere} on
triaxial ellipsoids and spheres using L{Ang}, L{Deg}, L{Rad} lat-, longitude, heading and meridian
(convergence) angles.

Transcoded to pure Python from I{Karney}'s GeographicLib 2.7 C++ class U{Conformal3<https://
GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1Triaxial_1_1Conformal3.html>}.

Copyright (C) U{Charles Karney<malto:Karney@Alum.MIT.edu>} (2024-2025) and licensed under the MIT/X11
License.  For more information, see the U{GeographicLib 2.7<https://GeographicLib.SourceForge.io/>}
documentation.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # noqa: E702 ;

from pygeodesy.angles import Ang, _Ang3Tuple, isAng,  Property_RO
from pygeodesy.basics import _copysign, signBit
from pygeodesy.constants import EPS, INF, NAN, PI_2, \
                               _copysign_1_0, _over, _1_over, remainder, \
                               _0_0, _0_01, _0_5, _0_75, _1_0, _2_0, _4_0
from pygeodesy.constants import _16_0  # PYCHOK used!
from pygeodesy.fmath import hypot1,  _ALL_LAZY, Fmt
from pygeodesy.interns import _DMAIN_, _scale_
# from pygeodesy.lazily import _ALL_LAZY  # from .fmath
# from pygeodesy.named import _NamedTuple, _Pass  # from .namedTuples
from pygeodesy.namedTuples import Vector2Tuple,  _Pass
# from pygeodesy.props import Property_RO  # from .triaxials.angles
# from pygeodesy.streprs import Fmt  # from .fmath
from pygeodesy.triaxials.bases import _bet_, Conformal5Tuple, LLK, _llk_, \
                                      _MAXIT, _omg_, TriaxialError, \
                                      _Triaxial3Base,  Vector3d
from pygeodesy.triaxials.triaxial3 import Triaxial3B
from pygeodesy.units import Easting, Northing, Radians, Radius_, Scalar
from pygeodesy.utily import sincos2
# from pygeodesy.vector3d import Vector3d  # from .triaxials.bases

from math import atan, atanh, exp, fabs, log, sinh, sqrt

__all__ = _ALL_LAZY.triaxials_conformal3
__version__ = '25.11.29'

_gam_    = 'gam'
_LOG2MIN =  log(EPS) * _2_0
# _N     = (log(_4_0) - log(PI)) / (PI_2 - log(_4_0))
_NLOG2   = -log(_2_0)
# _B     =  exp(_N * PI_2) - pow(_4_0, _N)


class BetOmgGam5Tuple(_Ang3Tuple):
    '''5-Tuple C{(bet, omg, gam, scale, llk)} with I{ellipsoidal} lat-
       C{bet}, longitude C{omg} and meridian convergence C{gam} all
       L{Ang}les, C{scale} I{and kind} C{llk} I{set to} C{LLK.CONFORMAL}.
    '''
    _Names_ = (_bet_, _omg_, _gam_, _scale_, _llk_)
    _Units_ = ( Ang,   Ang,  _Pass,  Scalar, _Pass)

    def __new__(cls, bet, omg, gam, scale=NAN, llk=None, unit=None, **kwds):  # **iteration_name
        args = bet, omg, gam, scale, (llk or LLK.CONFORMAL)
        self = _Ang3Tuple.__new__(cls, args, **kwds)
        return self.toUnit(unit) if unit else self


class Conformal3(_Triaxial3Base):
    '''I{Jacobi Conformal} projection of triaxial ellipsoid using class L{Ang}
       lat- and longitudes.

       @see: L{Triaxial<triaxials.triaxial5.Triaxial>} for details.
    '''
    @Property_RO
    def _a_b(self):
        return self.a / self.b

    @Property_RO
    def _a2_b(self):
        return self.a * self._a_b  # == self.b * self._a2_b2

    @Property_RO
    def _b_a(self):
        return self.b / self.a  # = sqrt(self._b2_a2)

    @Property_RO
    def _b_c(self):
        return self.b / self.c

    @Property_RO
    def _c_b(self):
        return self.c / self.b  # = sqrt(self._c2_b2)

    @Property_RO
    def _c2_b(self):
        return self.c * self._c_b  # == self.b * self._c2_b2

    @Property_RO
    def _C3S(self):
        '''(INTERNAL) Cache the I{equivalent Conformal Sphere}.
        '''
        return self.equi3Sphere(*self.xyQ2, name=self.name)

    def equi3Sphere(self, x, y, **name):
        '''Get this projection's I{equivalent Conformal Sphere}.

           @arg x: Quadrant x length, easting (C{meter}).
           @arg y: Quadrant y length, northing (C{meter}).

           @return: The C{Comformal3Sphere} of this projection.

           @see: Classes L{Conformal3Sphere<triaxials.spheres.Conformal3Sphere>} and
                 L{ConformalSphere<triaxials.triaxial5.ConformalSphere>} and I{Karney}'s
                 GeographicLib 2.7 C++ U{Triaxial::Conformal3<https://GeographicLib.
                 SourceForge.io/C++/doc/classGeographicLib_1_1Triaxial_1_1Conformal3.html>}
                 for more information.
        '''
        x, y = self.xyQ2
        # Find C{b, k2, kp2} s.t. C{b * K(kp2) = x},
        # C{b * K(k2) = y} and C{x * K(k2) - y * K(kp2) = 0}.
        _EF = self._Elliptic
        _xy = x < y
        if _xy:
            x, y = y, x
        # Now x >= y, k2 <= 1/2
        if x != y:
            s  =  x + y
            ny = _over(y, s)
            if ny:
                nx = _over(x, s)
                # assert nx != ny
                # Find initial guess assume K(k2) = pi/2, so K(kp2) = nx/ny * pi/2.
                # Invert using approximate k(K) given in https://arxiv.org/abs/2505.17159v4
                KK = _over(nx, ny) * PI_2
                # k2 = _16_0 / pow(exp(_N * KK) - _B, _2_0 / _N)
                # Alternatively using KK = 1/2*log(16/kp) A+S 17.3.26
                k2 = min(_0_5, exp(-KK * _2_0) * _16_0)  # Make sure guess is sane
                logk2 = log(k2)
                if logk2 > _LOG2MIN:
                    # Solve for log(k2) to preserve relative accuracy for tiny k2.
                    def _k2s2(logk2):
                        k2  =  exp(logk2)
                        kp2 = _1_0 - k2
                        xS  = _EF(kp2, 0, k2)
                        yS  = _EF(k2,  0, kp2)
                        f   =  nx *  yS.cK - ny  * xS.cK
                        fp  = (nx * (yS.cE - kp2 * yS.cK) +
                               ny * (xS.cE - k2  * xS.cK)) / (kp2 * _2_0)
                        return f, fp

                    logk2, _, _, _ = _root4(_k2s2, 0, logk2, _LOG2MIN, _NLOG2)
                    k2 = exp(logk2)
                # otherwise accept the asymptotic result
                kp2 = _1_0 - k2
            else:
                k2, kp2 = _0_0, _1_0
        else:
            k2 = kp2 = _0_5
        # b * K(kp2) = x
        # b * K(k2)  = y

        K = _EF(k2, 0, kp2).cK
        b = _over(y, K)
        if _xy:
            k2, kp2 = kp2, k2
        return Conformal3Sphere(b, k2, kp2, **name)

    def forwardBetOmg(self, bet, omg, M=False, **unit_name):
        '''Compute the projection to this conformal triaxial.

           @arg bet: Ellipsoidal latitude (C{Ang} or B{C{unit}}).
           @arg omg: Ellipsoidal longitude (C{Ang} or B{C{unit}})..
           @kwarg M: If C{True}, compute and include the scale (C{bool}).
           @kwarg unit_name: Optional C{B{unit}=}L{Radians} and C{B{name}=NN} (C{str}).

           @return: A L{Conformal5Tuple}C{(x, y, z, scale, llk)} with C{z = INT0}
                    I{always} and C{scale} if C{B{M}=True}, otherwise C{scale = NAN}.
        '''
        bet, omg, name = _bet_omg_name(bet, omg, **unit_name)
        m = _1_over(self._invScale(bet, omg)) if M else NAN
        return Conformal5Tuple(self._x(omg), self._y(bet), scale=m, **name)

    def forwardOther(self, other, bet, omg, M=False, **unit_name):
        '''Compute the projection to an other conformal triaxial.

           @arg other: A conformal triaxial (C{Conformal3}).
           @arg bet: Ellipsoidal latitude (C{Ang} or B{C{unit}}).
           @arg omg: Ellipsoidal longitude (C{Ang} or B{C{unit}}).
           @kwarg M: If C{True}, compute and include the scale (C{bool}).
           @kwarg unit_name: Optional C{B{unit}=}L{Radians} and C{B{name}=NN} (C{str}).

           @return: A L{BetOmgGam5Tuple}C{(bet, omg, gam, scale, llk)} with C{scale}
                    if C{B{M}=True}, otherwise C{scale = NAN}.
        '''
        if not isinstance(other, Conformal3):
            raise TriaxialError(other=other)
        bet, omg, name = _bet_omg_name(bet, omg, **unit_name)
        m = other._C3S.b / self._C3S.b
        ct, v, ma = self.forwardSphere3(bet, omg, M=M)
        ct = Vector3d(ct).times(m)
        bet, omg, gam, mb, _ = other.reverseSphere(ct, dir3d=v, M=M)
        if M:
            m *= _over(ma, mb)
        else:
            m  =  NAN
        return BetOmgGam5Tuple(bet, omg, gam, m, **name)

    def forwardSphere3(self, bet, omg, M=False, **unit_name):
        '''Compute the projection to and direction on the equivalent I{Conformal
           Sphere}.

           @arg bet: Ellipsoidal latitude (C{Ang} or B{C{unit}}).
           @arg omg: Ellipsoidal longitude (C{Ang} or B{C{unit}}).
           @kwarg M: If C{True}, compute and include the scale (C{bool}).
           @kwarg unit_name: Optional C{B{unit}=}L{Radians} and C{B{name}=NN} (C{str}).

           @return: 3-Tuple C{(cartesian, direction, scale)} with a C{cartesian}
                    L{Cartesian5Tuple}C{(x, y, z, h, llk)} on and C{direction} a
                    C{Vector3d} due North and tangent to the I{Conformal Sphere}
                    and C{scale} if C{B{M}=True}, otherwise C{scale = NAN}.
        '''
        if M:
            bet, omg, _ = _bet_omg_name(bet, omg, **unit_name)
        x, y, _, _, _ = self.forwardBetOmg(bet, omg, M=False, **unit_name)
        S = self._C3S
        bets = _invF(S._yE, y / S.b)
        omgs = _invF(S._xE, x / S.b).shift(-1)
        ct, dir3d = S.forwardBetOmgAlp2(bets, omgs, Ang.N(), **unit_name)
        if M:
            ma =  self._invScale(bet, omg)
            mb = _invScale(S, bets, omgs)
            m  =  self._sphScale(ma, mb)
        else:
            m  =  NAN
        return ct, dir3d, m

    def _invScale(self, bet, omg):
        # scale helper
        return _invScale(self, bet, omg.shift(+1))

    def reverseBetOmg(self, x_cf, y=None, M=False, **unit_name):
        '''Reverse a projection from this conformal triaxial.

           @arg x_cf: Easting (C{scalar}) or a conformal projection (C{Conformal5Tuple}).
           @kwarg y: Northing (C{scalar}), required if B{C{x_cf}} is C{scalar},
                     ignored otherwise.
           @kwarg M: If C{True}, compute and include the scale (C{bool}).
           @kwarg unit_name: Optional C{B{unit}=}L{Radians} and C{B{name}=NN} (C{str}).

           @return: A L{BetOmgGam5Tuple}C{(bet, omg, gam, scale, llk)} with
                    C{gam} set to C{None} and C{scale} only if C{B{M}=True},
                    otherwise C{scale is NAN}.
        '''
        x, y = _cf2en(x_cf, y)
        bet  = _invPi(self._yE, y / self._c2_b).mod(self._c_b)
        omg  = _invPi(self._xE, x / self._a2_b).mod(self._a_b).shift(-1)
        m    = _1_over(self._invScale(bet, omg)) if M else NAN
        return BetOmgGam5Tuple(bet, omg, None, m, **unit_name)

    def reverseOther(self, other, beto, omgo, M=False, **unit_name):
        '''Reverse a projection from an other conformal triaxial.

           @arg other: A conformal triaxial (C{Conformal3}).
           @arg beto: Ellipsoidal latitude on the B{C{other}} triaxial
                      (C{Ang} or B{C{unit}}).
           @arg omgo: Ellipsoidal longitude on the B{C{other}} triaxial
                      (C{Ang} or B{C{unit}}).
           @kwarg M: If C{True}, compute and include the scale (C{bool}).
           @kwarg unit_name: Optional C{B{unit}=}L{Radians} and C{B{name}=NN} (C{str}).

           @return: A L{BetOmgGam5Tuple}C{(bet, omg, gam, scale, llk)} with
                    C{scale} if C{B{M}=True}, otherwise C{scale = NAN}.
        '''
        if not isinstance(other, Conformal3):
            raise TriaxialError(other=other)
        bet, omg, gam, m, _ = other.forwardOther(self, beto, omgo, M=M)
        m = _1_over(m) if M else NAN
        return BetOmgGam5Tuple(bet, omg, _Neg(gam), m, **unit_name)

    def reverseSphere(self, x_ct, y=None, z=None, dir3d=None, M=False, **unit_name):
        '''Reverse a projection from this C{Conformal Sphere} to this triaxial.

           @arg x_ct: X component (C{scalar}) of or a cartesian on the C{Conformal
                      Sphere} (C{Cartesian5Tuple}).
           @kwarg y: Y component (C{scalar}), required if B{C{x_ct}} is C{scalar},
                     ignored otherwise.
           @kwarg z: Z component (C{scalar}), like B{C{y}}.
           @kwarg dir3d: The direction (C{Vector3d} or C{None}), reference.
           @kwarg M: If C{True}, compute and include the scale (C{bool}).
           @kwarg unit_name: Optional C{B{unit}=}L{Radians} and C{B{name}=NN} (C{str}).

           @return: A L{BetOmgGam5Tuple}C{(bet, omg, gam, scale, llk)} with
                    C{scale} only if C{B{M}=True}, otherwise C{scale = NAN}.
        '''
        S = self._C3S
        bets, omgs, alp, _, _ = S.reverseBetOmgAlp(x_ct, y, z, dir3d=dir3d)
        _  = Ang._norm(bets, omgs, alp, alt=bool(S.k2 == 0))
        x = S.b * _F(S._xE, omgs.shift(+1))
        y = S.b * _F(S._yE, bets)
        bet, omg, _, _, _ = self.reverseBetOmg(x, y, M=False)
        if M:
            mb = _invScale(S, bets, omgs)
            ma = _invScale(S, bet,  omg)
            m  =  self._sphScale(ma, mb)
        else:
            m  =  NAN
        return BetOmgGam5Tuple(bet, omg, _Neg(alp), m, **unit_name)

    def _sphScale(self, ma, mb):
        '''(INTERNAL) Compute the C{scale}.
        '''
        return _over(mb, ma) if ma and mb else self._sphScale_m

    @Property_RO
    def _sphScale_m(self):
        '''(INTERNAL) Cache the constant C{scale}.
        '''
        e = sqrt(self.e2)

        def _h1sg(atan_):
            sg = sinh(e * atan_(e))
            return hypot1(sg) + sg

        k2, kp2 = self._k2_kp2
        if not kp2:  # oblate pole
            m =  self._c_b * _h1sg(atanh)
        elif not k2:  # prolate pole
            m = _over(self._a_b, _h1sg(atan))
        else:  # trixial umbilical
            S =  self._C3S
            s = _over(S.k2 * S.kp2, k2 * kp2)
            m = (sqrt(s) * (self.b / S.b)) if s > 0 else _0_0
        return m

    def _x(self, omg):  # degrees?
        '''(INTERNAL) Compute an C{x} projection easting.
        '''
        omg = omg.shift(+1)
        _, c, n = omg.scn3
        if (n or signBit(c)) and not self.k2:
            x =  NAN
        else:
            x = _Pi(self._xE, omg.mod(self._b_a))
            x *= self._a2_b
        return x

    @Property_RO
    def xyQ2(self):
        '''Get the quadrant length in x and y direction (Vector2Tuple{x, y}).
        '''
        return Vector2Tuple(self._a2_b * self._xE.cPi,
                            self._c2_b * self._yE.cPi,
                            name=Conformal3.xyQ2.name)

    def _y(self, bet):  # degrees?
        '''(INTERNAL) Compute a C{y} projection northing.
        '''
        _, c, n = bet.scn3
        if (n or signBit(c)) and not self.kp2:
            y =  NAN
        else:
            y = _Pi(self._yE, bet.mod(self._b_c))
            y *= self._c2_b
        return y


class Conformal3B(Conformal3):
    '''I{Jacobi Conformal projection} on a triaxial ellipsoid
       specified by its middle semi-axis and shape.

       @see: L{Conformal3} for details and more information.
    '''
    def __init__(self, b, e2=_0_0, k2=_1_0, kp2=_0_0, **name):
        '''New L{Conformal3B} triaxial.

          @note: Use C{B{b}=radius} and C{B{e2}=0} for a conformal
                 I{spherical} projection.

          @see: L{Conformal<triaxials.triaxial5.Conformal.__init__>}.
        '''
        self._init_abc3_e2_k2_kp2(Radius_(b=b), e2, k2, kp2, **name)


class Conformal3Sphere(Triaxial3B):  # note C{Triaxial3}!
    '''I{Jacobi Conformal projection} on a I{spherical} triaxial.

       @see: Method L{equiv3Sphere<triaxials.conformal3.Conformal3.equi3Sphere>}.
    '''
    def __init__(self, radius, k2=_1_0, kp2=_0_0, **name):
        '''New, L{Conformal3Sphere} instance.

           @see: L{Triaxial3<triaxials.triaxial3.Triaxial3>} for more information.
        '''
        self._init_abc3_e2_k2_kp2(Radius_(radius), 0, k2, kp2, **name)
        # self._xE = self._Elliptic(kp2, 0, k2,  **name)  # _Triaxial3Base._xE
        # self._yE = self._Elliptic(k2,  0, kp2, **name)  # _Triaxial3Base._yE


def _bet_omg_name(bet, omg, unit=Radians, **name):
    '''(INTERNAL) Get C{(bet, omg, name)}.
    '''
    return (Ang.fromScalar(bet, unit=unit),
            Ang.fromScalar(omg, unit=unit), name)


def _cf2en(x_cf, y):
    '''(INTERNAL) Get easting C{x} and notrthing C{y}.
    '''
    return x_cf[:2] if isinstance(x_cf, Conformal5Tuple) else (
           Easting(x_cf), Northing(y))


def _F(eF, phi):  # -> float
    '''(INTERNAL) Elliptic function C{F(phi)}.
    '''
    s, c, n = phi.scn3
    if (n or signBit(c)) and not eF.kp2:
        p = NAN
    else:
        p = eF.fF(s, c, eF.fDelta(s, c))
        if n:
            p += eF.cK * n * _4_0
    return p


def _invF(eF, x):  # -> Ang
    '''(INTERNAL) Inverse elliptic function C{F(x)}.
    '''
    r, y = _invRy2(eF, eF.cK, x)
    if y:  # solve eF.fF(phi) = y for phi
        def _fF2(phi):  # -> pair<real, real>
            s, c = sincos2(phi)
            f  = eF.fF(s, c, eF.fDelta(s, c))
            fp = sqrt(eF.kp2 + c**2 * eF.k2)
            return f, _1_over(fp)

        z = fabs(y)
        z, _, _, _ = _root4(_fF2, z, z * PI_2 / eF.cK)
        r += Ang.fromRadians(_copysign(z, y))
    return r


def _invPi(eF, x):  # -> Ang
    '''(INTERNAL) Inverse elliptic function C{Pi(x)}.
    '''
    r, y = _invRy2(eF, eF.cPi, x)
    if y:  # solve eF.fPi(phi) = y for phi
        def _fPi2(phi):  # -> pair<real, real>
            s, c = sincos2(phi)
            f  = eF.fPi(s, c, eF.fDelta(s, c))
            a2 = eF.alpha2
            fp = (_1_0 - c**2 * a2) if a2 < 0 else \
                 (eF.alphap2 + s**2 * a2)
            fp *= sqrt(eF.kp2 + c**2 * eF.k2)
            return f, _1_over(fp)

        z = fabs(y)
        z, _, _, _ = _root4(_fPi2, z, z * PI_2 / eF.cPi)
        r += Ang.fromRadians(_copysign(z, y))
    return r


def _invRy2(eF, eF_c, x):
    # helper for C{_invF} and C{_invPi}
    if eF.kp2:
        n = eF_c * _2_0
        y = remainder(x, n)
        n = round((x - y) / n) * _2_0
    else:  # eF_c == N-/INF
        y, n = x, _0_0
    if not y:
        y, n = 0, (n if n else y)  # signed 0
    elif fabs(y) == eF_c:
        y, n = 0, (_copysign_1_0(y) + n)
    return Ang.cardinal(n), y


def _invScale(triax, bet, omg):
    # helper for triaxial, sphere scale
    k2, kp2 = triax._k2_kp2
    return sqrt(k2 * bet.c**2 + kp2 * omg.c**2)


def _Neg(ang):
    return (-ang) if isAng(ang) else (ang or None)


def _Pi(eF, phi):  # -> float
    '''(INTERNAL) Elliptic function C{Pi(phi)}.
    '''
    s, c, n = phi.scn3
    if (n or signBit(c)) and not eF.kp2:
        p = NAN
    else:
        p = eF.fPi(s, c, eF.fDelta(s, c))
        if n:
            p += eF.cPi * n * _4_0
    return p


def _root4(_fp2, z, x, xa=_0_0, xb=PI_2, xscale=1, zscale=1, s=1, tol=EPS):  # int s
    '''(INTERNAL) Solve v = _fp2(x) - z = 0.
    '''
    k = b = C = 0
    if xa < xb and xa <= x <= xb:
        # p = PI_2 * 0  #???
        # def _fp2z(x):
        #     f, fp = _fp2(x)
        #     f -= z
        #     # "DAT ", x, f, fp
        #     return f
        # a, _, b = map1(_fp2z, xa, x, xb)
        # if (a * b) > 0:
        #     raise TriaxalError('"DATBAD")
        # tol = max(tol, EPS)  # tol if tol > 0 else EPS
        vtol = tol * zscale * _0_01
        xtol = pow(tol, _0_75) * xscale
        oldv = oldx = olddx = INF
        for k in range(1, _MAXIT):
            # TODO: 20 60 -90 180 127.4974 24.6254 2.4377
            v, vp = _fp2(x)
            v -= z
            va =  fabs(v)
            dx = _over(-v, vp)
            # "XX ", k, (xa - p), (x - p), (xb - p), dx, (x + dx - p), v, vp
            if not (va > (0 if k < 2 else vtol)):
                C = 1  # k, va
                break
            elif (v * s) > 0:
                xb = min(xb, x)
            else:
                xa = max(xa, x)
            x  += dx
            dxa = fabs(dx)
            if x < xa or x > xb or va  > oldv or \
                        (k > 2 and dxa > olddx):
                b += 1  # k, xa, x, xb
                x  = (xa + xb) * _0_5
                if x == oldx:
                    C = 3  # k, x, dx
                    break
            elif not dxa > xtol:
                C = 2  # k, dx, xtol
                break
            # "GAPS ", k, dx, (x - xa), (xb - x), oldx, x, (oldx - x)
            oldx  = x
            oldv  = va
            olddx = dxa * _0_5
        else:
            t = Fmt.no_convergence(dx, xtol)
            raise TriaxialError(x=x, xa=xa, xb=xb, txt=t)
    else:
        x = NAN
    return x, k, b, C


if __name__ == _DMAIN_:

    from pygeodesy import Degrees, printf
    from pygeodesy.triaxials import Triaxials

    # <https://GeographicLib.SourceForge.io/C++/doc/Cart3Convert.1.html>
    T = Conformal3(Triaxials.WGS84_3)
    printf(T)
    # name='WGS84_3', a=6378171.36, b=6378101.609999999, c=6356751.84, ...
    t = T.forwardBetOmg(Ang.fromDegrees(33.3), Ang.fromDegrees(44.4), M=True)
    printf((t.x, t.y, t.scale))
    # (-5077802.461853351, 3922572.0186951873, 1.197034384522207)
    #  -5077732.396        3922571.859         1.1970343759  C++
    t = T.reverseBetOmg(*t[:2], M=True)
    printf((t.bet.degrees0, t.omg.degrees0, t.scale))
    # (33.47654394192169, 44.39937131735643, 1.1994622456567812)
    #  33.30000000        44.40000000        1.1970343759  C++

    T = Conformal3(Triaxials.WGS84_3r)  # rounded
    printf(T)
    # name='WGS84_3r', a=6378172, b=6378102, c=6356752, ...
    t = T.forwardBetOmg(Degrees(33.3), Degrees(44.4), M=True)
    printf((t.x, t.y, t.scale))
    # (-5077802.439189989, 3922571.859124643, 1.197034375926918)
    #  -5077732.396        3922571.859        1.1970343759  C++
    t = T.reverseBetOmg(*t[:2], M=True)
    printf((t.bet.degrees0, t.omg.degrees0, t.scale))
    # (33.47654654826102, 44.39937131735643, 1.1994622731246583)
    #  33.30000000        44.40000000        1.1970343759  C++

    c = 6378137 * (1 - 1 / (298257223563 / 1000000000))
    T = Conformal3(6378172, 6378102, c)
    printf(T)
    # name='', a=6378172, b=6378102, c=6356752.314245179, ...
    t = T.forwardBetOmg(Degrees(33.3), Degrees(44.4), M=True)
    printf((t.x, t.y, t.scale))
    # (-5077802.461853351, 3922572.0186951873, 1.197034384522207)
    #  -5077732.396        3922571.859         1.1970343759  C++
    t = T.reverseBetOmg(*t[:2], M=True)
    printf((t.bet.degrees0, t.omg.degrees0, t.scale))
    # (33.47654394192169, 44.39937131735643, 1.1994622456567812)
    #  33.30000000        44.40000000        1.1970343759  C++

# **) MIT License
#
# Copyright (C) 2025-2026 -- mrJean1 at Gmail -- All Rights Reserved.
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
