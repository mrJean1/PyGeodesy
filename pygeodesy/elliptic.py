
# -*- coding: utf-8 -*-

u'''I{Karney}'s elliptic functions and integrals.

Class L{Elliptic} transcoded from I{Charles Karney}'s C++ class U{EllipticFunction
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1EllipticFunction.html>}
to pure Python, including symmetric integrals L{Elliptic.fRC}, L{Elliptic.fRD},
L{Elliptic.fRF}, L{Elliptic.fRG} and L{Elliptic.fRJ} as C{static methods}.

Python method names follow the C++ member functions, I{except}:

 - member functions I{without arguments} are mapped to Python properties
   prefixed with C{"c"}, for example C{E()} is property C{cE},

 - member functions with 1 or 3 arguments are renamed to Python methods
   starting with an C{"f"}, example C{E(psi)} to C{fE(psi)} and C{E(sn,
   cn, dn)} to C{fE(sn, cn, dn)},

 - other Python method names conventionally start with a lower-case
   letter or an underscore if private.

Following is a copy of I{Karney}'s U{EllipticFunction.hpp
<https://GeographicLib.SourceForge.io/C++/doc/EllipticFunction_8hpp_source.html>}
file C{Header}.

Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2008-2023)
and licensed under the MIT/X11 License.  For more information, see the
U{GeographicLib<https://GeographicLib.SourceForge.io>} documentation.

B{Elliptic integrals and functions.}

This provides the elliptic functions and integrals needed for
C{Ellipsoid}, C{GeodesicExact}, and C{TransverseMercatorExact}.  Two
categories of function are provided:

 - functions to compute U{symmetric elliptic integrals
   <https://DLMF.NIST.gov/19.16.i>}

 - methods to compute U{Legrendre's elliptic integrals
   <https://DLMF.NIST.gov/19.2.ii>} and U{Jacobi elliptic
   functions<https://DLMF.NIST.gov/22.2>}.

In the latter case, an object is constructed giving the modulus
C{k} (and optionally the parameter C{alpha}).  The modulus (and
parameter) are always passed as squares which allows C{k} to be
pure imaginary.  (Confusingly, Abramowitz and Stegun call C{m = k**2}
the "parameter" and C{n = alpha**2} the "characteristic".)

In geodesic applications, it is convenient to separate the incomplete
integrals into secular and periodic components, e.g.

I{C{E(phi, k) = (2 E(k) / pi) [ phi + delta E(phi, k) ]}}

where I{C{delta E(phi, k)}} is an odd periodic function with
period I{C{pi}}.

The computation of the elliptic integrals uses the algorithms given
in U{B. C. Carlson, Computation of real or complex elliptic integrals
<https://DOI.org/10.1007/BF02198293>} (also available U{here
<https://ArXiv.org/pdf/math/9409227.pdf>}), Numerical Algorithms 10,
13--26 (1995) with the additional optimizations given U{here
<https://DLMF.NIST.gov/19.36.i>}.

The computation of the Jacobi elliptic functions uses the algorithm
given in U{R. Bulirsch, Numerical Calculation of Elliptic Integrals
and Elliptic Functions<https://DOI.org/10.1007/BF01397975>},
Numerische Mathematik 7, 78--90 (1965).

The notation follows U{NIST Digital Library of Mathematical Functions
<https://DLMF.NIST.gov>} chapters U{19<https://DLMF.NIST.gov/19>} and
U{22<https://DLMF.NIST.gov/22>}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import copysign0, map2, neg, neg_
from pygeodesy.constants import EPS, INF, NAN, PI, PI_2, PI_4, \
                               _EPStol as _TolJAC, _0_0, _1_64th, \
                               _0_25, _0_5, _1_0, _2_0, _N_2_0, \
                               _3_0, _4_0, _6_0, _8_0, _180_0, \
                               _360_0, _over
from pygeodesy.errors import _ValueError, _xattr, _xkwds_pop
from pygeodesy.fmath import fdot, hypot1, zqrt
from pygeodesy.fsums import Fsum, _sum
from pygeodesy.interns import NN, _delta_, _DOT_, _f_, _invalid_, \
                             _invokation_, _negative_, _SPACE_
from pygeodesy.karney import _K_2_0, _norm180, _signBit, _sincos2, \
                             _ALL_LAZY
# from pygeodesy.lazily import _ALL_LAZY  # from .karney
from pygeodesy.named import _Named, _NamedTuple,  Fmt, unstr
from pygeodesy.props import _allPropertiesOf_n, Property_RO, _update_all
# from pygeodesy.streprs import Fmt, unstr  # from .named
from pygeodesy.units import Scalar, Scalar_
# from pygeodesy.utily import sincos2 as _sincos2  # from .karney

from math import asinh, atan, atan2, ceil, cosh, fabs, floor, \
                 radians, sin, sqrt, tanh

__all__ = _ALL_LAZY.elliptic
__version__ = '23.09.14'

_TolRD  =  zqrt(EPS * 0.002)
_TolRF  =  zqrt(EPS * 0.030)
_TolRG0 = _TolJAC   * 2.7
_TRIPS  =  21  # Max depth, 7 might be sufficient


class _Cs(object):
    '''(INTERAL) Complete integrals cache.
    '''
    def __init__(self, **kwds):
        self.__dict__ = kwds


class _D(list):
    '''(INTERNAL) Deferred C{Fsum}.
    '''
    def __call__(self, s):
        try:  # Fsum *= s
            return Fsum(*self).fmul(s)
        except ValueError:  # Fsum(NAN) exception
            return _sum(self) * s

    def __iadd__(self, x):
        list.append(self, x)
        return self


class Elliptic(_Named):
    '''Elliptic integrals and functions.

       @see: I{Karney}'s U{Detailed Description<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1EllipticFunction.html#details>}.
    '''
#   _alpha2  = 0
#   _alphap2 = 0
#   _eps     = EPS
#   _k2      = 0
#   _kp2     = 0

    def __init__(self, k2=0, alpha2=0, kp2=None, alphap2=None, name=NN):
        '''Constructor, specifying the C{modulus} and C{parameter}.

           @kwarg name: Optional name (C{str}).

           @see: Method L{Elliptic.reset} for further details.

           @note: If only elliptic integrals of the first and second kinds
                  are needed, use C{B{alpha2}=0}, the default value.  In
                  that case, we have C{Π(φ, 0, k) = F(φ, k), G(φ, 0, k) =
                  E(φ, k)} and C{H(φ, 0, k) = F(φ, k) - D(φ, k)}.
        '''
        self.reset(k2=k2, alpha2=alpha2, kp2=kp2, alphap2=alphap2)

        if name:
            self.name = name

    @Property_RO
    def alpha2(self):
        '''Get α^2, the square of the parameter (C{float}).
        '''
        return self._alpha2

    @Property_RO
    def alphap2(self):
        '''Get α'^2, the square of the complementary parameter (C{float}).
        '''
        return self._alphap2

    @Property_RO
    def cD(self):
        '''Get Jahnke's complete integral C{D(k)} (C{float}),
           U{defined<https://DLMF.NIST.gov/19.2.E6>}.
        '''
        return self._cDEKEeps.cD

    @Property_RO
    def _cDEKEeps(self):
        '''(INTERNAL) Get the complete integrals D, E, K and KE plus C{eps}.
        '''
        k2, kp2 = self.k2, self.kp2
        if k2:
            if kp2:
                try:
                    self._iteration = 0
                    # D(k) = (K(k) - E(k))/k2, Carlson eq.4.3
                    # <https://DLMF.NIST.gov/19.25.E1>
                    D   = _RD(self, _0_0, kp2, _1_0, _3_0)
                    cD  =  float(D)
                    # Complete elliptic integral E(k), Carlson eq. 4.2
                    # <https://DLMF.NIST.gov/19.25.E1>
                    cE  = _rG2(self, kp2, _1_0, PI_=PI_2)
                    # Complete elliptic integral K(k), Carlson eq. 4.1
                    # <https://DLMF.NIST.gov/19.25.E1>
                    cK  = _rF2(self, kp2, _1_0)
                    cKE =  float(D.fmul(k2))
                    eps =  k2 / (sqrt(kp2) + _1_0)**2

                except Exception as e:
                    raise _ellipticError(self.reset, k2=k2, kp2=kp2, cause=e)
            else:
                cD  =  cK = cKE = INF
                cE  = _1_0
                eps =  k2
        else:
            cD  =  PI_4
            cE  =  cK = PI_2
            cKE = _0_0  # k2 * cD
            eps =  EPS

        return _Cs(cD=cD, cE=cE, cK=cK, cKE=cKE, eps=eps)

    @Property_RO
    def cE(self):
        '''Get the complete integral of the second kind C{E(k)}
           (C{float}), U{defined<https://DLMF.NIST.gov/19.2.E5>}.
        '''
        return self._cDEKEeps.cE

    @Property_RO
    def cG(self):
        '''Get Legendre's complete geodesic longitude integral
           C{G(α^2, k)} (C{float}).
        '''
        return self._cGHPi.cG

    @Property_RO
    def _cGHPi(self):
        '''(INTERNAL) Get the complete integrals G, H and Pi.
        '''
        alpha2, alphap2, kp2 = self.alpha2, self.alphap2, self.kp2
        try:
            self._iteration = 0
            if alpha2:
                if alphap2:
                    if kp2:  # <https://DLMF.NIST.gov/19.25.E2>
                        cK =  self.cK
                        Rj = _RJ(self, _0_0, kp2, _1_0, alphap2, _3_0)
                        cG =  float(Rj * (alpha2 - self.k2) + cK)  # G(alpha2, k)
                        cH = -float(Rj *  alphap2 - cK)  # H(alpha2, k)
                        cPi = float(Rj *  alpha2  + cK)  # Pi(alpha2, k)
                    else:
                        cG  = cH = _rC(self, _1_0, alphap2)
                        cPi = INF  # XXX or NAN?
                else:
                    cG = cH = cPi = INF  # XXX or NAN?
            else:
                cG, cPi = self.cE, self.cK
                # H = K - D but this involves large cancellations if k2 is near 1.
                # So write (for alpha2 = 0)
                #   H = int(cos(phi)**2 / sqrt(1-k2 * sin(phi)**2), phi, 0, pi/2)
                #     = 1 / sqrt(1-k2) * int(sin(phi)**2 / sqrt(1-k2/kp2 * sin(phi)**2,...)
                #     = 1 / kp * D(i * k/kp)
                # and use D(k) = RD(0, kp2, 1) / 3, so
                #   H = 1/kp * RD(0, 1/kp2, 1) / 3
                #     = kp2 * RD(0, 1, kp2) / 3
                # using <https://DLMF.NIST.gov/19.20.E18>.  Equivalently
                #   RF(x, 1) - RD(0, x, 1) / 3 = x * RD(0, 1, x) / 3 for x > 0
                # For k2 = 1 and alpha2 = 0, we have
                #   H = int(cos(phi),...) = 1
                cH = float(_RD(self, _0_0, _1_0, kp2, _3_0 / kp2)) if kp2 else _1_0

        except Exception as e:
            raise _ellipticError(self.reset, kp2=kp2, alpha2 =alpha2,
                                                      alphap2=alphap2, cause=e)
        return _Cs(cG=cG, cH=cH, cPi=cPi)

    @Property_RO
    def cH(self):
        '''Get Cayley's complete geodesic longitude difference integral
           C{H(α^2, k)} (C{float}).
        '''
        return self._cGHPi.cH

    @Property_RO
    def cK(self):
        '''Get the complete integral of the first kind C{K(k)}
           (C{float}), U{defined<https://DLMF.NIST.gov/19.2.E4>}.
        '''
        return self._cDEKEeps.cK

    @Property_RO
    def cKE(self):
        '''Get the difference between the complete integrals of the
           first and second kinds, C{K(k) − E(k)} (C{float}).
        '''
        return self._cDEKEeps.cKE

    @Property_RO
    def cPi(self):
        '''Get the complete integral of the third kind C{Pi(α^2, k)}
           (C{float}), U{defined<https://DLMF.NIST.gov/19.2.E7>}.
        '''
        return self._cGHPi.cPi

    def deltaD(self, sn, cn, dn):
        '''Jahnke's periodic incomplete elliptic integral.

           @arg sn: sin(φ).
           @arg cn: cos(φ).
           @arg dn: sqrt(1 − k2 * sin(2φ)).

           @return: Periodic function π D(φ, k) / (2 D(k)) - φ (C{float}).

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        return _deltaX(sn, cn, dn, self.cD, self.fD)

    def deltaE(self, sn, cn, dn):
        '''The periodic incomplete integral of the second kind.

           @arg sn: sin(φ).
           @arg cn: cos(φ).
           @arg dn: sqrt(1 − k2 * sin(2φ)).

           @return: Periodic function π E(φ, k) / (2 E(k)) - φ (C{float}).

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        return _deltaX(sn, cn, dn, self.cE, self.fE)

    def deltaEinv(self, stau, ctau):
        '''The periodic inverse of the incomplete integral of the second kind.

           @arg stau: sin(τ)
           @arg ctau: cos(τ)

           @return: Periodic function E^−1(τ (2 E(k)/π), k) - τ (C{float}).

           @raise EllipticError: No convergence.
        '''
        try:
            if _signBit(ctau):  # pi periodic
                stau, ctau = neg_(stau, ctau)
            t = atan2(stau, ctau)
            return self._Einv(t * self.cE / PI_2) - t

        except Exception as e:
            raise _ellipticError(self.deltaEinv, stau, ctau, cause=e)

    def deltaF(self, sn, cn, dn):
        '''The periodic incomplete integral of the first kind.

           @arg sn: sin(φ).
           @arg cn: cos(φ).
           @arg dn: sqrt(1 − k2 * sin(2φ)).

           @return: Periodic function π F(φ, k) / (2 K(k)) - φ (C{float}).

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        return _deltaX(sn, cn, dn, self.cK, self.fF)

    def deltaG(self, sn, cn, dn):
        '''Legendre's periodic geodesic longitude integral.

           @arg sn: sin(φ).
           @arg cn: cos(φ).
           @arg dn: sqrt(1 − k2 * sin(2φ)).

           @return: Periodic function π G(φ, k) / (2 G(k)) - φ (C{float}).

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        return _deltaX(sn, cn, dn, self.cG, self.fG)

    def deltaH(self, sn, cn, dn):
        '''Cayley's periodic geodesic longitude difference integral.

           @arg sn: sin(φ).
           @arg cn: cos(φ).
           @arg dn: sqrt(1 − k2 * sin(2φ)).

           @return: Periodic function π H(φ, k) / (2 H(k)) - φ (C{float}).

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        return _deltaX(sn, cn, dn, self.cH, self.fH)

    def deltaPi(self, sn, cn, dn):
        '''The periodic incomplete integral of the third kind.

           @arg sn: sin(φ).
           @arg cn: cos(φ).
           @arg dn: sqrt(1 − k2 * sin(2φ)).

           @return: Periodic function π Π(φ, α2, k) / (2 Π(α2, k)) - φ
                    (C{float}).

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        return _deltaX(sn, cn, dn, self.cPi, self.fPi)

    def _Einv(self, x):
        '''(INTERNAL) Helper for C{.deltaEinv} and C{.fEinv}.
        '''
        E2 = self.cE * _2_0
        n  = floor(x / E2 + _0_5)
        r  = x - E2 * n  # r in [-cE, cE)
        # linear approximation
        phi = PI * r / E2  # phi in [-PI_2, PI_2)
        Phi = Fsum(phi)
        # first order correction
        phi = Phi.fsum_(self.eps * sin(phi * _2_0) / _N_2_0)
        # For kp2 close to zero use asin(r / cE) or J. P. Boyd,
        # Applied Math. and Computation 218, 7005-7013 (2012)
        # <https://DOI.org/10.1016/j.amc.2011.12.021>
        _Phi2, self._iteration = Phi.fsum2_, 0  # aggregate
        for i in range(1, _TRIPS):  # GEOGRAPHICLIB_PANIC
            sn, cn, dn = self._sncndn3(phi)
            if dn:
                sn = self.fE(sn, cn, dn)
                phi, d = _Phi2((r - sn) / dn)
            else:  # PYCHOK no cover
                d = _0_0  # XXX continue?
            if fabs(d) < _TolJAC:  # 3-4 trips
                _iterations(self, i)
                break
        else:  # PYCHOK no cover
            raise _convergenceError(d, _TolJAC)
        return Phi.fsum_(n * PI) if n else phi

    @Property_RO
    def eps(self):
        '''Get epsilon (C{float}).
        '''
        return self._cDEKEeps.eps

    def fD(self, phi_or_sn, cn=None, dn=None):
        '''Jahnke's incomplete elliptic integral in terms of
           Jacobi elliptic functions.

           @arg phi_or_sn: φ or sin(φ).
           @kwarg cn: C{None} or cos(φ).
           @kwarg dn: C{None} or sqrt(1 − k2 * sin(2φ)).

           @return: D(φ, k) as though φ ∈ (−π, π] (C{float}),
                    U{defined<https://DLMF.NIST.gov/19.2.E4>}.

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        def _fD(sn, cn, dn):
            r = fabs(sn)**3
            if r:
                r = float(_RD(self, cn**2, dn**2, _1_0, _3_0 / r))
            return r

        return self._fXf(phi_or_sn, cn, dn, self.cD,
                                            self.deltaD, _fD)

    def fDelta(self, sn, cn):
        '''The C{Delta} amplitude function.

           @arg sn: sin(φ).
           @arg cn: cos(φ).

           @return: sqrt(1 − k2 * sin(2φ)) (C{float}).
        '''
        try:
            k2 =  self.k2
            s  = (self.kp2 + cn**2 * k2) if k2 > 0 else (
                     (_1_0 - sn**2 * k2) if k2 < 0 else self.kp2)
            return sqrt(s) if s else _0_0

        except Exception as e:
            raise _ellipticError(self.fDelta, sn, cn, k2=k2, cause=e)

    def fE(self, phi_or_sn, cn=None, dn=None):
        '''The incomplete integral of the second kind in terms of
           Jacobi elliptic functions.

           @arg phi_or_sn: φ or sin(φ).
           @kwarg cn: C{None} or cos(φ).
           @kwarg dn: C{None} or sqrt(1 − k2 * sin(2φ)).

           @return: E(φ, k) as though φ ∈ (−π, π] (C{float}),
                    U{defined<https://DLMF.NIST.gov/19.2.E5>}.

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        def _fE(sn, cn, dn):
            '''(INTERNAL) Core of C{.fE}.
            '''
            if sn:
                sn2, cn2, dn2 = sn**2, cn**2, dn**2
                kp2, k2 = self.kp2, self.k2
                if k2 <= 0:  # Carlson, eq. 4.6, <https://DLMF.NIST.gov/19.25.E9>
                    Ei = _RF3(self, cn2, dn2, _1_0)
                    if k2:
                        Ei -= _RD(self, cn2, dn2, _1_0, _3over(k2, sn2))
                elif kp2 >= 0:  # k2 > 0, <https://DLMF.NIST.gov/19.25.E10>
                    Ei = _over(k2 * fabs(cn), dn)  # float
                    if kp2:
                        Ei += (_RD( self, cn2, _1_0,  dn2, _3over(k2, sn2)) +
                               _RF3(self, cn2,  dn2, _1_0)) * kp2
                else:  # kp2 < 0, <https://DLMF.NIST.gov/19.25.E11>
                    Ei  = _over(dn, fabs(cn))
                    Ei -= _RD(self, dn2, _1_0, cn2, _3over(kp2, sn2))
                Ei *= fabs(sn)
                ei  = float(Ei)
            else:  # PYCHOK no cover
                ei = _0_0
            return ei

        return self._fXf(phi_or_sn, cn, dn, self.cE,
                                            self.deltaE, _fE)

    def fEd(self, deg):
        '''The incomplete integral of the second kind with
           the argument given in C{degrees}.

           @arg deg: Angle (C{degrees}).

           @return: E(π B{C{deg}} / 180, k) (C{float}).

           @raise EllipticError: No convergence.
        '''
        if _K_2_0:
            e =  round((deg - _norm180(deg)) / _360_0)
        elif fabs(deg) < _180_0:
            e = _0_0
        else:
            e =  ceil(deg / _360_0 - _0_5)
            deg -= e * _360_0
        return self.fE(radians(deg)) + e * self.cE * _4_0

    def fEinv(self, x):
        '''The inverse of the incomplete integral of the second kind.

           @arg x: Argument (C{float}).

           @return: φ = 1 / E(B{C{x}}, k), such that E(φ, k) = B{C{x}}
                    (C{float}).

           @raise EllipticError: No convergence.
        '''
        try:
            return self._Einv(x)
        except Exception as e:
            raise _ellipticError(self.fEinv, x, cause=e)

    def fF(self, phi_or_sn, cn=None, dn=None):
        '''The incomplete integral of the first kind in terms of
           Jacobi elliptic functions.

           @arg phi_or_sn: φ or sin(φ).
           @kwarg cn: C{None} or cos(φ).
           @kwarg dn: C{None} or sqrt(1 − k2 * sin(2φ)).

           @return: F(φ, k) as though φ ∈ (−π, π] (C{float}),
                    U{defined<https://DLMF.NIST.gov/19.2.E4>}.

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        def _fF(sn, cn, dn):
            r = fabs(sn)
            if r:
                r = float(_RF3(self, cn**2, dn**2, _1_0).fmul(r))
            return r

        return self._fXf(phi_or_sn, cn, dn, self.cK,
                                            self.deltaF, _fF)

    def fG(self, phi_or_sn, cn=None, dn=None):
        '''Legendre's geodesic longitude integral in terms of
           Jacobi elliptic functions.

           @arg phi_or_sn: φ or sin(φ).
           @kwarg cn: C{None} or cos(φ).
           @kwarg dn: C{None} or sqrt(1 − k2 * sin(2φ)).

           @return: G(φ, k) as though φ ∈ (−π, π] (C{float}).

           @raise EllipticError: Invalid invokation or no convergence.

           @note: Legendre expresses the longitude of a point on the
                  geodesic in terms of this combination of elliptic
                  integrals in U{Exercices de Calcul Intégral, Vol 1
                  (1811), p 181<https://Books.Google.com/books?id=
                  riIOAAAAQAAJ&pg=PA181>}.

           @see: U{Geodesics in terms of elliptic integrals<https://
                 GeographicLib.SourceForge.io/html/geodesic.html#geodellip>}
                 for the expression for the longitude in terms of this function.
        '''
        return self._fXa(phi_or_sn, cn, dn, self.alpha2 - self.k2,
                                            self.cG, self.deltaG)

    def fH(self, phi_or_sn, cn=None, dn=None):
        '''Cayley's geodesic longitude difference integral in terms of
           Jacobi elliptic functions.

           @arg phi_or_sn: φ or sin(φ).
           @kwarg cn: C{None} or cos(φ).
           @kwarg dn: C{None} or sqrt(1 − k2 * sin(2φ)).

           @return: H(φ, k) as though φ ∈ (−π, π] (C{float}).

           @raise EllipticError: Invalid invokation or no convergence.

           @note: Cayley expresses the longitude difference of a point
                  on the geodesic in terms of this combination of
                  elliptic integrals in U{Phil. Mag. B{40} (1870), p 333
                  <https://Books.Google.com/books?id=Zk0wAAAAIAAJ&pg=PA333>}.

           @see: U{Geodesics in terms of elliptic integrals<https://
                 GeographicLib.SourceForge.io/html/geodesic.html#geodellip>}
                 for the expression for the longitude in terms of this function.
        '''
        return self._fXa(phi_or_sn, cn, dn, -self.alphap2,
                                             self.cH, self.deltaH)

    def fPi(self, phi_or_sn, cn=None, dn=None):
        '''The incomplete integral of the third kind in terms of
           Jacobi elliptic functions.

           @arg phi_or_sn: φ or sin(φ).
           @kwarg cn: C{None} or cos(φ).
           @kwarg dn: C{None} or sqrt(1 − k2 * sin(2φ)).

           @return: Π(φ, α2, k) as though φ ∈ (−π, π] (C{float}).

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        if dn is None and cn is not None:  # and isscalar(phi_or_sn)
            dn = self.fDelta(phi_or_sn, cn)  # in .triaxial
        return self._fXa(phi_or_sn, cn, dn, self.alpha2,
                                            self.cPi, self.deltaPi)

    def _fXa(self, phi_or_sn, cn, dn, aX, cX, deltaX):
        '''(INTERNAL) Helper for C{.fG}, C{.fH} and C{.fPi}.
        '''
        def _fX(sn, cn, dn):
            if sn:
                cn2, dn2 = cn**2, dn**2
                R = _RF3(self, cn2, dn2, _1_0)
                if aX:
                    sn2 = sn**2
                    p   = sn2 * self.alphap2 + cn2
                    R += _RJ(self, cn2, dn2, _1_0, p, _3over(aX, sn2))
                R *= fabs(sn)
                r  = float(R)
            else:  # PYCHOK no cover
                r = _0_0
            return r

        return self._fXf(phi_or_sn, cn, dn, cX, deltaX, _fX)

    def _fXf(self, phi_or_sn, cn, dn, cX, deltaX, fX):
        '''(INTERNAL) Helper for C{.fD}, C{.fE}, C{.fF} and C{._fXa}.
        '''
        self._iteration = 0  # aggregate
        phi = sn = phi_or_sn
        if cn is dn is None:  # fX(phi) call
            sn, cn, dn = self._sncndn3(phi)
            if fabs(phi) >= PI:
                return (deltaX(sn, cn, dn) + phi) * cX / PI_2
            # fall through
        elif cn is None or dn is None:
            n = NN(_f_, deltaX.__name__[5:])
            raise _ellipticError(n, sn, cn, dn)

        if _signBit(cn):  # enforce usual trig-like symmetries
            xi = cX * _2_0 - fX(sn, cn, dn)
        else:
            xi = fX(sn, cn, dn) if cn > 0 else cX
        return copysign0(xi, sn)

    @Property_RO
    def k2(self):
        '''Get k^2, the square of the modulus (C{float}).
        '''
        return self._k2

    @Property_RO
    def kp2(self):
        '''Get k'^2, the square of the complementary modulus (C{float}).
        '''
        return self._kp2

    def reset(self, k2=0, alpha2=0, kp2=None, alphap2=None):  # MCCABE 13
        '''Reset the modulus, parameter and the complementaries.

           @kwarg k2: Modulus squared (C{float}, NINF <= k^2 <= 1).
           @kwarg alpha2: Parameter squared (C{float}, NINF <= α^2 <= 1).
           @kwarg kp2: Complementary modulus squared (C{float}, k'^2 >= 0).
           @kwarg alphap2: Complementary parameter squared (C{float}, α'^2 >= 0).

           @raise EllipticError: Invalid B{C{k2}}, B{C{alpha2}}, B{C{kp2}}
                                 or B{C{alphap2}}.

           @note: The arguments must satisfy C{B{k2} + B{kp2} = 1} and
                  C{B{alpha2} + B{alphap2} = 1}.  No checking is done
                  that these conditions are met to enable accuracy to be
                  maintained, e.g., when C{k} is very close to unity.
        '''
        if self.__dict__:
            _update_all(self, _Named.iteration._uname, Base=Property_RO)

        self._k2  = Scalar_(k2=k2, Error=EllipticError, low=None, high=_1_0)
        self._kp2 = Scalar_(kp2=((_1_0 - k2) if kp2 is None else kp2), Error=EllipticError)

        self._alpha2  = Scalar_(alpha2=alpha2, Error=EllipticError, low=None, high=_1_0)
        self._alphap2 = Scalar_(alphap2=((_1_0 - alpha2) if alphap2 is None else alphap2),
                                Error=EllipticError)

        # Values of complete elliptic integrals for k = 0,1 and alpha = 0,1
        #         K     E     D
        # k = 0:  pi/2  pi/2  pi/4
        # k = 1:  inf   1     inf
        #                    Pi    G     H
        # k = 0, alpha = 0:  pi/2  pi/2  pi/4
        # k = 1, alpha = 0:  inf   1     1
        # k = 0, alpha = 1:  inf   inf   pi/2
        # k = 1, alpha = 1:  inf   inf   inf
        #
        # G(0, k) = Pi(0, k) = H(1, k) = E(k)
        # H(0, k) = K(k) - D(k)
        # Pi(alpha2, 0) = G(alpha2, 0) = pi / (2 * sqrt(1 - alpha2))
        # H( alpha2, 0) = pi / (2 * (sqrt(1 - alpha2) + 1))
        # Pi(alpha2, 1) = inf
        # G( alpha2, 1) = H(alpha2, 1) = RC(1, alphap2)

    def sncndn(self, x):
        '''The Jacobi elliptic function.

           @arg x: The argument (C{float}).

           @return: An L{Elliptic3Tuple}C{(sn, cn, dn)} with
                    C{*n(B{x}, k)}.

           @raise EllipticError: No convergence.
        '''
        self._iteration = 0  # reset
        try:  # Bulirsch's sncndn routine, p 89.
            if self.kp2:
                c, d, cd, mn = self._sncndn4
                dn = _1_0
                sn, cn = _sincos2(x * cd)
                if sn:
                    a  = cn / sn
                    c *= a
                    for m, n in reversed(mn):
                        a *= c
                        c *= dn
                        dn = (n + a) / (m + a)
                        a  =  c / m
                    a  = _1_0 / hypot1(c)
                    sn =  neg(a) if _signBit(sn) else a
                    cn =  c * sn
                    if d and _signBit(self.kp2):
                        cn, dn = dn, cn
                        sn = sn / d  # /= chokes PyChecker
            else:
                sn = tanh(x)
                cn = dn = _1_0 / cosh(x)

        except Exception as e:
            raise _ellipticError(self.sncndn, x, kp2=self.kp2, cause=e)

        return Elliptic3Tuple(sn, cn, dn, iteration=self._iteration)

    def _sncndn3(self, phi):
        '''(INTERNAL) Helper for C{.fEinv} and C{._fXf}.
        '''
        sn, cn = _sincos2(phi)
        return sn, cn, self.fDelta(sn, cn)

    @Property_RO
    def _sncndn4(self):
        '''(INTERNAL) Get Bulirsch' 4-tuple C{(c, d, cd, mn)}.
        '''
        # Bulirsch's sncndn routine, p 89.
        d, mc = 0, self.kp2
        if _signBit(mc):
            d  = _1_0 - mc
            mc =  neg(mc / d)
            d  =  sqrt(d)

        mn, a = [], _1_0
        for i in range(1, _TRIPS):  # GEOGRAPHICLIB_PANIC
            mc = sqrt(mc)
            mn.append((a, mc))
            c = (a + mc) * _0_5
            r =  fabs(mc -  a)
            t = _TolJAC  *  a
            if r <= t:  # 6 trips, quadratic
                _iterations(self, i)
                break
            mc *= a
            a   = c
        else:  # PYCHOK no cover
            raise _convergenceError(r, t)
        cd = (c * d) if d else c
        return c, d, cd, mn

    @staticmethod
    def fRC(x, y):
        '''Degenerate symmetric integral of the first kind C{RC(x, y)}.

           @return: C{RC(x, y)}, equivalent to C{RF(x, y, y)}.

           @see: U{C{RC} definition<https://DLMF.NIST.gov/19.2.E17>} and
                 U{Carlson<https://ArXiv.org/pdf/math/9409227.pdf>}.
        '''
        return _rC(None, x, y)

    @staticmethod
    def fRD(x, y, z, *over):
        '''Degenerate symmetric integral of the third kind C{RD(x, y, z)}.

           @return: C{RD(x, y, z) / over}, equivalent to C{RJ(x, y, z, z)
                    / over} with C{over} typically 3.

           @see: U{C{RD} definition<https://DLMF.NIST.gov/19.16.E5>} and
                 U{Carlson<https://ArXiv.org/pdf/math/9409227.pdf>}.
        '''
        try:
            return float(_RD(None, x, y, z, *over))
        except Exception as e:
            raise _ellipticError(Elliptic.fRD, x, y, z, *over, cause=e)

    @staticmethod
    def fRF(x, y, z=0):
        '''Symmetric or complete symmetric integral of the first kind
           C{RF(x, y, z)} respectively C{RF(x, y)}.

           @return: C{RF(x, y, z)} or C{RF(x, y)} for missing or zero B{C{z}}.

           @see: U{C{RF} definition<https://DLMF.NIST.gov/19.16.E1>} and
                 U{Carlson<https://ArXiv.org/pdf/math/9409227.pdf>}.
        '''
        try:
            return float(_RF3(None, x, y, z)) if z else _rF2(None, x, y)
        except Exception as e:
            raise _ellipticError(Elliptic.fRF, x, y, z, cause=e)

    @staticmethod
    def fRG(x, y, z=0):
        '''Symmetric or complete symmetric integral of the second kind
           C{RG(x, y, z)} respectively C{RG(x, y)}.

           @return: C{RG(x, y, z)} or C{RG(x, y)} for missing or zero B{C{z}}.

           @see: U{C{RG} definition<https://DLMF.NIST.gov/19.16.E3>},
                 U{Carlson<https://ArXiv.org/pdf/math/9409227.pdf>} and
                 U{RG<https://GeographicLib.SourceForge.io/C++/doc/
                 EllipticFunction_8cpp_source.html#l00096>} version 2.3.
        '''
        try:
            return _rG2(None, x, y) if z == 0 else (
                   _rG2(None, z, x) if y == 0 else (
                   _rG2(None, y, z) if x == 0 else _rG3(None, x, y, z)))
        except Exception as e:
            t = _negative_ if min(x, y, z) < 0 else NN
            raise _ellipticError(Elliptic.fRG, x, y, z, cause=e, txt=t)

    @staticmethod
    def fRJ(x, y, z, p):
        '''Symmetric integral of the third kind C{RJ(x, y, z, p)}.

           @return: C{RJ(x, y, z, p)}.

           @see: U{C{RJ} definition<https://DLMF.NIST.gov/19.16.E2>} and
                 U{Carlson<https://ArXiv.org/pdf/math/9409227.pdf>}.
        '''
        try:
            return float(_RJ(None, x, y, z, p))
        except Exception as e:
            raise _ellipticError(Elliptic.fRJ, x, y, z, p, cause=e)

    @staticmethod
    def _RFRD(x, y, z, m):
        # in .auxilats.AuxDLat.DE, .auxilats.AuxLat.Rectifying
        try:  # float(RF(x, y, z) - RD(x, y, z, 3 / m))
            R = _RF3(None, x, y, z)
            if m:
                R -= _RD(None, x, y, z, _3_0 / m)
        except Exception as e:
            raise _ellipticError(Elliptic._RFRD, x, y, z, m, cause=e)
        return float(R)

_allPropertiesOf_n(15, Elliptic)  #  # PYCHOK assert, see Elliptic.reset


class EllipticError(_ValueError):
    '''Elliptic function, integral, convergence or other L{Elliptic} issue.
    '''
    pass


class Elliptic3Tuple(_NamedTuple):
    '''3-Tuple C{(sn, cn, dn)} all C{scalar}.
    '''
    _Names_ = ('sn',   'cn',   'dn')
    _Units_ = ( Scalar, Scalar, Scalar)


class _L(list):
    '''(INTERNAL) Helper for C{_RD}, C{_RF3} and C{_RJ}.
    '''
    _a0   = None
#   _xyzp = ()

    def __init__(self, *xyzp):  # x, y, z [, p]
        list.__init__(self, xyzp)
        self._xyzp = xyzp

    def a0(self, n):
        '''Compute the initial C{a}.
        '''
        t = tuple(self)
        m = n - len(t)
        if m > 0:
            t += t[-1:] * m
        try:
            a = Fsum(*t).fover(n)
        except ValueError:  # Fsum(NAN) exception
            a = _sum(t) / n
        self._a0 = a
        return a

    def amrs4(self, inst, y, Tol):
        '''Yield Carlson 4-tuples C{(An, mul, lam, s)} plus sentinel, with
           C{lam = fdot(sqrt(x), ... (z))} and C{s = (sqrt(x), ... (p))}.
        '''
        L = self
        a = L.a0(5 if y else 3)
        m = 1
        t = max(fabs(a - _) for _ in L) / Tol
        for i in range(_TRIPS):
            d = fabs(a * m)
            if d > t:  # 3-6 trips
                _iterations(inst, i)
                break
            s = map2(sqrt, L)  # sqrt(x), srqt(y), sqrt(z) [, sqrt(p)]
            try:
                r = fdot(s[:3], s[1], s[2], s[0])  # sqrt(x) * sqrt(y) + ...
            except ValueError:  # Fsum(NAN) exception
                r = _sum(s[i] * s[(i + 1) % 3] for i in range(3))
            L[:] = [(r + _) * _0_25 for _ in L]
            a    =  (r + a) * _0_25
            if y:  # yield only if used
                yield a, m, r, s  # L[2] is next z
            m *= 4
        else:  # PYCHOK no cover
            raise _convergenceError(d, t, thresh=True)
        yield a, m, None, ()  # sentinel: same a, next m, no r and s

    def rescale(self, am, *xs):
        '''Rescale C{x}, C{y}, ...
        '''
        # assert am
        a0 = self._a0
        for x in xs:
            yield (a0 - x) / am


def _ab2(inst, x, y):
    '''(INTERNAL) Yield Carlson 2-tuples C{(xn, yn)}.
    '''
    a, b = sqrt(x), sqrt(y)
    if b > a:
        a, b = b, a
    yield a, b  # initial x0, y0
    for i in range(_TRIPS):
        d =  fabs(a - b)
        t = _TolRG0 * a
        if d <= t:  # 3-4 trips
            _iterations(inst, i)
            break
        a, b = ((a + b) * _0_5), sqrt(a * b)
        yield a, b  # xn, yn
    else:  # PYCHOK no cover
        raise _convergenceError(d, t)


def _convergenceError(d, tol, **thresh):
    '''(INTERNAL) Format a no-convergence Error.
    '''
    t = Fmt.no_convergence(d, tol, **thresh)
    return ValueError(t)  # txt only


def _deltaX(sn, cn, dn, cX, fX):
    '''(INTERNAL) Helper for C{Elliptic.deltaD} thru C{.deltaPi}.
    '''
    try:
        if cn is None or dn is None:
            raise ValueError(_invalid_)

        if _signBit(cn):
            sn, cn = neg_(sn, cn)
        r = fX(sn, cn, dn) * PI_2 / cX
        return r - atan2(sn, cn)

    except Exception as e:
        n = NN(_delta_, fX.__name__[1:])
        raise _ellipticError(n, sn, cn, dn, cause=e)


def _ellipticError(where, *args, **kwds_cause_txt):
    '''(INTERNAL) Format an L{EllipticError}.
    '''
    c = _xkwds_pop(kwds_cause_txt, cause=None)
    t = _xkwds_pop(kwds_cause_txt, txt=NN)
    n = _xattr(where, __name__=where)  # _dunder_nameof(where, where)
    n = _DOT_(Elliptic.__name__, n)
    n = _SPACE_(_invokation_, n)
    u =  unstr(n, *args, **kwds_cause_txt)
    return EllipticError(u, cause=c, txt=t)


def _Horner(S, e1, E2, E3, E4, E5, *over):
    '''(INTERNAL) Horner form for C{_RD} and C{_RJ} below.
    '''
    E22 = E2**2
    # Polynomial is <https://DLMF.NIST.gov/19.36.E2>
    # (1 - 3*E2/14 + E3/6 + 9*E2**2/88 - 3*E4/22 - 9*E2*E3/52
    #    + 3*E5/26 - E2**3/16 + 3*E3**2/40 + 3*E2*E4/20
    #    + 45*E2**2*E3/272 - 9*(E3*E4+E2*E5)/68)
    # converted to Horner-like form ...
    F  = Fsum
    e  = e1 * 4084080
    S *= e
    S += F(E2 * -540540,                               471240).fmul(E5)
    S += F(E2 *  612612,                E3 * -540540, -556920).fmul(E4)
    S += F(E2 * -706860, E22 *  675675, E3 *  306306,  680680).fmul(E3)
    S += F(E2 *  417690, E22 * -255255,               -875160).fmul(E2)
    S += 4084080
    return S.fdiv((over[0] * e) if over else e)  # Fsum


def _iterations(inst, i):
    '''(INTERNAL) Aggregate iterations B{C{i}}.
    '''
    if inst and i > 0:
        inst._iteration += i


def _3over(a, b):
    '''(INTERNAL) Return C{3 / (a * b)}.
    '''
    return _over(_3_0, a * b)


def _rC(unused, x, y):
    '''(INTERNAL) Defined only for C{y != 0} and C{x >= 0}.
    '''
    d = x - y
    if d < 0:  # catch NaN
        # <https://DLMF.NIST.gov/19.2.E18>
        d = -d
        r = atan(sqrt(d / x)) if x > 0 else PI_2
    elif d == 0:  # XXX d < EPS0? or EPS02 or _EPSmin
        d, r = y, _1_0
    elif y > 0:  # <https://DLMF.NIST.gov/19.2.E19>
        r = asinh(sqrt(d / y))  # atanh(sqrt((x - y) / x))
    elif y < 0:  # <https://DLMF.NIST.gov/19.2.E20>
        r = asinh(sqrt(-x / y))  # atanh(sqrt(x / (x - y)))
    else:  # PYCHOK no cover
        raise _ellipticError(Elliptic.fRC, x, y)
    return r / sqrt(d)  # float


def _RD(inst, x, y, z, *over):
    '''(INTERNAL) Carlson, eqs 2.28 - 2.34.
    '''
    L = _L(x, y, z)
    S = _D()
    for a, m, r, s in L.amrs4(inst, True, _TolRF):
        if s:
            S += _over(_3_0, (r + z) * s[2] * m)
            z  =  L[2]  # s[2] = sqrt(z)
    x, y = L.rescale(-a * m, x, y)
    xy =  x * y
    z  = (x + y) / _3_0
    z2 =  z**2
    return _Horner(S(_1_0), sqrt(a) * a * m,
                   xy        - _6_0 * z2,
                  (xy * _3_0 - _8_0 * z2) * z,
                  (xy - z2)  * _3_0 * z2,
                   xy * z2 * z, *over)  # Fsum


def _rF2(inst, x, y):  # 2-arg version, z=0
    '''(INTERNAL) Carlson, eqs 2.36 - 2.38.
    '''
    for a, b in _ab2(inst, x, y):  # PYCHOK yield
        pass
    return _over(PI, a + b)  # float


def _RF3(inst, x, y, z):  # 3-arg version
    '''(INTERNAL) Carlson, eqs 2.2 - 2.7.
    '''
    L = _L(x, y, z)
    for a, m, _, _ in L.amrs4(inst, False, _TolRF):
        pass
    x, y = L.rescale(a * m, x, y)
    z  = neg(x + y)
    xy = x  * y
    e2 = xy - z**2
    e3 = xy * z
    e4 = e2**2
    # Polynomial is <https://DLMF.NIST.gov/19.36.E1>
    # (1 - E2/10 + E3/14 + E2**2/24 - 3*E2*E3/44
    #    - 5*E2**3/208 + 3*E3**2/104 + E2**2*E3/16)
    # converted to Horner-like form ...
    S  = Fsum(e4 * 15015, e3 * 6930, e2 * -16380,  17160).fmul(e3)
    S += Fsum(e4 * -5775,            e2 *  10010, -24024).fmul(e2)
    S += 240240
    return S.fdiv(sqrt(a) * 240240)  # Fsum


def _rG2(inst, x, y, PI_=PI_4):  # 2-args
    '''(INTERNAL) Carlson, eqs 2.36 - 2.39.
    '''
    m = -1  # neg!
    S = _D()
    # assert not S
    for a, b in _ab2(inst, x, y):  # PYCHOK yield
        if S:
            S += (a - b)**2 *  m
            m *=  2
        else:  # initial
            S += (a + b)**2 * _0_5
    return S(PI_).fover(a + b)


def _rG3(inst, x, y, z):  # 3-arg version
    '''(INTERNAL) C{x}, C{y} and C{z} all non-zero, see C{.fRG}.
    '''
    R  = _RF3(inst, x, y, z) * z
    rd = (x - z) * (z - y)  # - (y - z)
    if rd:  # Carlson, eq 1.7
        R += _RD(inst, x, y, z, _3_0 / rd)
    R += sqrt(x * y / z)
    return R.fover(_2_0)


def _RJ(inst, x, y, z, p, *over):
    '''(INTERNAL) Carlson, eqs 2.17 - 2.25.
    '''
    def _xyzp(x, y, z, p):
        return (x + p) * (y + p) * (z + p)

    L = _L(x, y, z, p)
    n =  neg(_xyzp(x, y, z, -p))
    S = _D()
    for a, m, _, s in L.amrs4(inst, True, _TolRD):
        if s:
            d = _xyzp(*s)
            if d:
                if n:
                    rc = _rC(inst, _1_0, n / d**2 + _1_0)
                    n *= _1_64th  # /= chokes PyChecker
                else:
                    rc = _1_0  # == _rC(None, _1_0, _1_0)
                S += rc / (d * m)
            else:  # PYCHOK no cover
                return NAN
    x, y, z = L.rescale(a * m, x, y, z)
    p   = Fsum(x, y, z).fover(_N_2_0)
    p2  = p**2
    p3  = p2 * p
    E2  = Fsum(x * y, x * z, y * z, -p2 * _3_0)
    E2p = E2 * p
    xyz = x * y * z
    return _Horner(S(_6_0), sqrt(a) * a * m, E2,
                   Fsum(p3 * _4_0, xyz, E2p * _2_0),
                   Fsum(p3 * _3_0, E2p, xyz * _2_0).fmul(p),
                   xyz * p2, *over)  # Fsum

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
