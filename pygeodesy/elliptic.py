
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

Copyright (C) U{Charles Karney<mailto:Charles@Karney.com>} (2008-2022)
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

from pygeodesy.basics import copysign0, map2, neg
# from pygeodesy.errors import _ValueError  # from .fmath
from pygeodesy.fmath import fdot, Fsum, hypot1, _ValueError
# from pygeodesy.fsums import Fsum  # from .fmath
from pygeodesy.interns import EPS, INF, NN, PI, PI_2, PI_4, \
                             _DOT_, _EPStol as _TolJAC, _convergence_, \
                             _f_, _no_, _SPACE_, _0_0, _0_125, _0_25, \
                             _0_5, _1_0, _2_0, _N_2_0, _3_0, _4_0, \
                             _5_0, _6_0, _8_0, _180_0, _360_0
from pygeodesy.karney import _ALL_LAZY, _signBit
# from pygeodesy.lazily import _ALL_LAZY  # from .karney
from pygeodesy.named import _Named, _NamedTuple, unstr
from pygeodesy.props import Property_RO, _update_all
# from pygeodesy.streprs import unstr  # from .named
from pygeodesy.units import Scalar, Scalar_
from pygeodesy.utily import sincos2, sincos2d

from math import asinh, atan, atan2, ceil, cosh, floor, sin, \
                 sqrt, tanh

__all__ = _ALL_LAZY.elliptic
__version__ = '22.07.08'

_delta_      = 'delta'
_invokation_ = 'invokation'
_1_64th      = _1_0 / 64  # pow(2.0, -6)
_TolRD       =  pow(EPS * 0.002, _0_125)  # 8th root: quadquadratic?, octic?, ocrt?
_TolRF       =  pow(EPS * 0.030, _0_125)  # 4th root: biquadratic, quartic, qurt?
_TolRG0      = _TolJAC  * 2.7
_TRIPS       =  31  # Max depth, 7 might be sufficient


class _Complete(object):
    '''(INTERAL) Hold complete integrals.
    '''
    def __init__(self, **kwds):
        self.__dict__ = kwds


class Elliptic(_Named):
    '''Elliptic integrals and functions.

       @see: I{Karney}'s U{Detailed Description<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1EllipticFunction.html#details>}.
    '''
    _alpha2  = 0
    _alphap2 = 0
    _eps     = EPS
    _k2      = 0
    _kp2     = 0

    def __init__(self, k2=0, alpha2=0, kp2=None, alphap2=None):
        '''Constructor, specifying the C{modulus} and C{parameter}.

           @see: Method L{Elliptic.reset} for further details.

           @note: If only elliptic integrals of the first and second kinds
                  are needed, use C{B{alpha2}=0}, the default value.  In
                  that case, we have C{Π(φ, 0, k) = F(φ, k), G(φ, 0, k) =
                  E(φ, k)} and C{H(φ, 0, k) = F(φ, k) - D(φ, k)}.
        '''
        self.reset(k2=k2, alpha2=alpha2, kp2=kp2, alphap2=alphap2)

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
        return self._reset4integrals.cD

    @Property_RO
    def cE(self):
        '''Get the complete integral of the second kind C{E(k)}
           (C{float}), U{defined<https://DLMF.NIST.gov/19.2.E5>}.
        '''
        return self._reset4integrals.cE

    @Property_RO
    def cG(self):
        '''Get Legendre's complete geodesic longitude integral
           C{G(α^2, k)} (C{float}).
        '''
        return self._reset3integrals.cG

    @Property_RO
    def cH(self):
        '''Get Cayley's complete geodesic longitude difference integral
           C{H(α^2, k)} (C{float}).
        '''
        return self._reset3integrals.cH

    @Property_RO
    def cK(self):
        '''Get the complete integral of the first kind C{K(k)}
           (C{float}), U{defined<https://DLMF.NIST.gov/19.2.E4>}.
        '''
        return self._reset4integrals.cK

    @Property_RO
    def cKE(self):
        '''Get the difference between the complete integrals of the
           first and second kinds, C{K(k) − E(k)} (C{float}).
        '''
        return self._reset4integrals.cKE

    @Property_RO
    def cPi(self):
        '''Get the complete integral of the third kind C{Pi(α^2, k)}
           (C{float}), U{defined<https://DLMF.NIST.gov/19.2.E7>}.
        '''
        return self._reset3integrals.cPi

    def deltaD(self, sn, cn, dn):
        '''The periodic Jahnke's incomplete elliptic integral.

           @arg sn: sin(φ).
           @arg cn: cos(φ).
           @arg dn: sqrt(1 − k2 sin(2φ)).

           @return: Periodic function π D(φ, k) / (2 D(k)) - φ (C{float}).

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        return self._deltaX(sn, cn, dn, self.cD, self.fD)

    def deltaE(self, sn, cn, dn):
        '''The periodic incomplete integral of the second kind.

           @arg sn: sin(φ).
           @arg cn: cos(φ).
           @arg dn: sqrt(1 − k2 sin(2φ)).

           @return: Periodic function π E(φ, k) / (2 E(k)) - φ (C{float}).

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        return self._deltaX(sn, cn, dn, self.cE, self.fE)

    def deltaEinv(self, stau, ctau):
        '''The periodic inverse of the incomplete integral of the second kind.

           @arg stau: sin(τ)
           @arg ctau: cos(τ)

           @return: Periodic function E^−1(τ (2 E(k)/π), k) - τ (C{float}).

           @raise EllipticError: No convergence.
        '''
        # Function is periodic with period pi
        t = atan2(-stau, -ctau) if _signBit(ctau) else atan2(stau, ctau)
        return self.fEinv(t * self.cE / PI_2) - t

    def deltaF(self, sn, cn, dn):
        '''The periodic incomplete integral of the first kind.

           @arg sn: sin(φ).
           @arg cn: cos(φ).
           @arg dn: sqrt(1 − k2 sin(2φ)).

           @return: Periodic function π F(φ, k) / (2 K(k)) - φ (C{float}).

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        return self._deltaX(sn, cn, dn, self.cK, self.fF)

    def deltaG(self, sn, cn, dn):
        '''Legendre's periodic geodesic longitude integral.

           @arg sn: sin(φ).
           @arg cn: cos(φ).
           @arg dn: sqrt(1 − k2 sin(2φ)).

           @return: Periodic function π G(φ, k) / (2 G(k)) - φ (C{float}).

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        return self._deltaX(sn, cn, dn, self.cG, self.fG)

    def deltaH(self, sn, cn, dn):
        '''Cayley's periodic geodesic longitude difference integral.

           @arg sn: sin(φ).
           @arg cn: cos(φ).
           @arg dn: sqrt(1 − k2 sin(2φ)).

           @return: Periodic function π H(φ, k) / (2 H(k)) - φ (C{float}).

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        return self._deltaX(sn, cn, dn, self.cH, self.fH)

    def deltaPi(self, sn, cn, dn):
        '''The periodic incomplete integral of the third kind.

           @arg sn: sin(φ).
           @arg cn: cos(φ).
           @arg dn: sqrt(1 − k2 sin(2φ)).

           @return: Periodic function π Π(φ, α2, k) / (2 Π(α2, k)) - φ
                    (C{float}).

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        return self._deltaX(sn, cn, dn, self.cPi, self.fPi)

    def _deltaX(self, sn, cn, dn, cX, fX):
        '''(INTERNAL) Helper for C{.deltaD} thru C{.deltaPi}.
        '''
        if cn is None or dn is None:
            n = NN(_delta_, fX.__name__[1:])
            raise _callError(n, sn, cn, dn)

        if _signBit(cn):
            cn, sn = -cn, -sn
        return fX(sn, cn, dn) * PI_2 / cX - atan2(sn, cn)

    @Property_RO
    def eps(self):
        '''Get epsilon (C{float}).
        '''
        return self._reset4integrals.eps

    def fD(self, phi_or_sn, cn=None, dn=None):
        '''Jahnke's incomplete elliptic integral in terms of
           Jacobi elliptic functions.

           @arg phi_or_sn: φ or sin(φ).
           @kwarg cn: C{None} or cos(φ).
           @kwarg dn: C{None} or sqrt(1 − k2 sin(2φ)).

           @return: D(φ, k) as though φ ∈ (−π, π] (C{float}),
                    U{defined<https://DLMF.NIST.gov/19.2.E4>}.

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        def _fD(sn, cn, dn):
            r = abs(sn)**3
            if r:
                r = _RD(self, cn**2, dn**2, _1_0, _3_0 / r)
            return r

        return self._fXf(phi_or_sn, cn, dn, self.cD,
                                            self.deltaD, _fD)

    def fDelta(self, sn, cn):
        '''The C{Delta} amplitude function.

           @arg sn: sin(φ).
           @arg cn: cos(φ).

           @return: sqrt(1 − k2 sin(2φ)) (C{float}).
        '''
        k2, kp2 = self.k2, self.kp2
        s = (_1_0 - k2 * sn**2) if k2 < 0 else (
             (kp2 + k2 * cn**2) if k2 > 0 else kp2)
        return sqrt(s) if s else _0_0

    def fE(self, phi_or_sn, cn=None, dn=None):
        '''The incomplete integral of the second kind in terms of
           Jacobi elliptic functions.

           @arg phi_or_sn: φ or sin(φ).
           @kwarg cn: C{None} or cos(φ).
           @kwarg dn: C{None} or sqrt(1 − k2 sin(2φ)).

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
                    ei = _RF3(self, cn2, dn2, _1_0)
                    if k2:
                        ei -= _RD(self, cn2, dn2, _1_0, _3rd(k2, sn2))
                elif kp2 >= 0:  # <https://DLMF.NIST.gov/19.25.E10>
                    ei = k2 * abs(cn) / dn
                    if kp2:
                        ei += (_RD(self, cn2, _1_0, dn2, _3rd(k2, sn2)) +
                               _RF3(self, cn2, dn2, _1_0)) * kp2
                else:  # <https://DLMF.NIST.gov/19.25.E11>
                    ei = dn / abs(cn) - _RD(self, dn2, _1_0, cn2, _3rd(kp2, sn2))
                ei *= abs(sn)
            else:  # PYCHOK no cover
                ei = _0_0
            return ei

        return self._fXf(phi_or_sn, cn, dn, self.cE,
                                            self.deltaE, _fE)

    def fEd(self, deg):
        '''The incomplete integral of the second kind with
           the argument given in degrees.

           @arg deg: Angle (C{degrees}).

           @return: E(π B{C{deg}}/180, k) (C{float}).

           @raise EllipticError: No convergence.
        '''
        if abs(deg) < _180_0:
            e = _0_0
        else:  # PYCHOK no cover
            e    = ceil(deg / _360_0 - _0_5)
            deg -= e * _360_0
            e   *= self.cE * _4_0
        sn, cn = sincos2d(deg)
        return self.fE(sn, cn, self.fDelta(sn, cn)) + e

    def fEinv(self, x):
        '''The inverse of the incomplete integral of the second kind.

           @arg x: Argument (C{float}).

           @return: φ = 1 / E(B{C{x}}, k), such that E(φ, k) = B{C{x}}
                    (C{float}).

           @raise EllipticError: No convergence.
        '''
        E2 = self.cE * _2_0
        n  = floor(x / E2 + _0_5)
        y  = x - E2 * n  # y now in [-ec, ec)
        # linear approximation
        phi = PI * y / E2  # phi in [-pi/2, pi/2)
        Phi = Fsum(phi)
        # first order correction
        phi = Phi.fsum_(self.eps * sin(phi * _2_0) / _N_2_0)
        # For kp2 close to zero use asin(x/.cE) or J. P. Boyd,
        # Applied Math. and Computation 218, 7005-7013 (2012)
        # <https://DOI.org/10.1016/j.amc.2011.12.021>
        Phi_fsum2_ = Phi.fsum2_
        fE         = self.fE
        _sncndnPhi = self._sncndnPhi
        self._iteration = 0  # aggregate
        for i in range(1, _TRIPS):  # GEOGRAPHICLIB_PANIC
            sn, cn, dn = _sncndnPhi(phi)
            phi, e = Phi_fsum2_((y - fE(sn, cn, dn)) / dn)
            if abs(e) < _TolJAC:
                self._iteration += i
                break
        else:  # PYCHOK no cover
            raise _convergenceError(_TolJAC, self.fEinv, x)
        return Phi.fsum_(n * PI) if n else phi

    def fF(self, phi_or_sn, cn=None, dn=None):
        '''The incomplete integral of the first kind in terms of
           Jacobi elliptic functions.

           @arg phi_or_sn: φ or sin(φ).
           @kwarg cn: C{None} or cos(φ).
           @kwarg dn: C{None} or sqrt(1 − k2 sin(2φ)).

           @return: F(φ, k) as though φ ∈ (−π, π] (C{float}),
                    U{defined<https://DLMF.NIST.gov/19.2.E4>}.

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        def _fF(sn, cn, dn):
            r = abs(sn)
            if r:
                r *= _RF3(self, cn**2, dn**2, _1_0)
            return r

        return self._fXf(phi_or_sn, cn, dn, self.cK,
                                            self.deltaF, _fF)

    def fG(self, phi_or_sn, cn=None, dn=None):
        '''Legendre's geodesic longitude integral in terms of
           Jacobi elliptic functions.

           @arg phi_or_sn: φ or sin(φ).
           @kwarg cn: C{None} or cos(φ).
           @kwarg dn: C{None} or sqrt(1 − k2 sin(2φ)).

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
           @kwarg dn: C{None} or sqrt(1 − k2 sin(2φ)).

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
           @kwarg dn: C{None} or sqrt(1 − k2 sin(2φ)).

           @return: Π(φ, α2, k) as though φ ∈ (−π, π] (C{float}).

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        return self._fXa(phi_or_sn, cn, dn, self.alpha2,
                                            self.cPi, self.deltaPi)

    def _fXa(self, phi_or_sn, cn, dn, aX, cX, deltaX):
        '''(INTERNAL) Helper for C{.fG}, C{.fH} and C{.fPi}.
        '''
        def _fX(sn, cn, dn):
            if sn:
                cn2, sn2, dn2 = cn**2, sn**2, dn**2
                r = _RF3(self, cn2, dn2, _1_0)
                if aX:
                    z  =  cn2 + sn2 * self.alphap2
                    r += _RJ(self, cn2, dn2, _1_0, z, _3rd(aX, sn2))
                r *= abs(sn)
            else:  # PYCHOK no cover
                r = _0_0
            return r

        return self._fXf(phi_or_sn, cn, dn, cX, deltaX, _fX)

    def _fXf(self, phi_or_sn, cn, dn, cX, deltaX, fX):
        '''(INTERNAL) Helper for C{f.D}, C{.fE}, C{.fF} and C{._fXa}.
        '''
        self._iteration = 0  # aggregate
        phi = sn = phi_or_sn
        if cn is dn is None:  # fX(phi) call
            sn, cn, dn = self._sncndnPhi(phi)
            if abs(phi) >= PI:  # PYCHOK no cover
                return (deltaX(sn, cn, dn) + phi) * cX / PI_2
            # fall through
        elif cn is None or dn is None:
            n = NN(_f_, deltaX.__name__[5:])
            raise _callError(n, sn, cn, dn)

        if _signBit(cn):  # enforce usual trig-like symmetries
            xi = _2_0 * cX - fX(sn, cn, dn)
        elif cn > 0:
            xi = fX(sn, cn, dn)
        else:
            xi = cX
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
        _update_all(self, _Named.iteration._uname)

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

    @Property_RO
    def _reset3integrals(self):
        '''(INTERNAL) Get the complete integrals G, H and Pi.
        '''
        self._iteration = 0
        alpha2 = self.alpha2
        if alpha2:
            alphap2 = self.alphap2
            if alphap2:
                kp2 = self.kp2
                if kp2:  # <https://DLMF.NIST.gov/19.25.E2>
                    rj   = _RJ(self, _0_0, kp2, _1_0, alphap2, _3_0)
                    cPi  =  cH = cG = self.cK
                    cG  += (alpha2 - self.k2) * rj  # G(alpha2, k)
                    cH  -=  alphap2 * rj  # H(alpha2, k)
                    cPi +=  alpha2  * rj  # Pi(alpha2, k)
                else:  # PYCHOK no cover
                    cG  = cH = _RC(self, _1_0, alphap2)
                    cPi = INF  # XXX or NAN?
            else:  # PYCHOK no cover
                cG = cH = cPi = INF  # XXX or NAN?
        else:
            cG, cPi, kp2 = self.cE, self.cK, self.kp2
            # H = K - D but this involves large cancellations if k2 is near 1.
            # So write (for alpha2 = 0)
            #   H = int(cos(phi)**2/sqrt(1-k2*sin(phi)**2),phi,0,pi/2)
            #     = 1/sqrt(1-k2) * int(sin(phi)**2/sqrt(1-k2/kp2*sin(phi)**2,...)
            #     = 1/kp * D(i*k/kp)
            # and use D(k) = RD(0, kp2, 1) / 3
            # so H = 1/kp * RD(0, 1/kp2, 1) / 3
            #      = kp2 * RD(0, 1, kp2) / 3
            # using <https://DLMF.NIST.gov/19.20.E18>.  Equivalently
            #   RF(x, 1) - RD(0, x, 1)/3 = x * RD(0, 1, x)/3 for x > 0
            # For k2 = 1 and alpha2 = 0, we have
            #   H = int(cos(phi),...) = 1
            cH = _RD(self, _0_0, _1_0, kp2, _3_0 / kp2) if kp2 else _1_0

        return _Complete(cG=cG, cH=cH, cPi=cPi)

    @Property_RO
    def _reset4integrals(self):
        '''(INTERNAL) Get the complete integrals D, E, K and KE plus C{eps}.
        '''
        k2 = self.k2
        if k2:
            kp2 = self.kp2
            if kp2:
                self._iteration = 0
                # D(k) = (K(k) - E(k))/k2, Carlson eq.4.3
                # <https://DLMF.NIST.gov/19.25.E1>
                cD  = _RD(self, _0_0, kp2, _1_0, _3_0)
                # Complete elliptic integral E(k), Carlson eq. 4.2
                # <https://DLMF.NIST.gov/19.25.E1>
                cE  = _RG2(self, kp2, _1_0)
                # Complete elliptic integral K(k), Carlson eq. 4.1
                # <https://DLMF.NIST.gov/19.25.E1>
                cK  = _RF2(self, kp2, _1_0)
                cKE =  k2 * cD
                eps =  k2 / (sqrt(kp2) + _1_0)**2
            else:  # PYCHOK no cover
                cD  =  cK = cKE = INF
                cE  = _1_0
                eps =  k2
        else:  # PYCHOK no cover
            cD  =  PI_4
            cE  =  cK = PI_2
            cKE = _0_0  # k2 * cD
            eps =  EPS

        return _Complete(cD=cD, cE=cE, cK=cK, cKE=cKE, eps=eps)

    def sncndn(self, x):
        '''The Jacobi elliptic function.

           @arg x: The argument (C{float}).

           @return: An L{Elliptic3Tuple}C{(sn, cn, dn)} with
                    C{*n(B{x}, k)}.

           @raise EllipticError: No convergence.
        '''
        self._iteration = 0  # reset
        # Bulirsch's sncndn routine, p 89.
        if self.kp2:
            c, d, mn_ = self._sncndnBulirsch3
            dn = _1_0
            x *= (c * d) if d else c
            sn, cn = sincos2(x)
            if sn:
                a  = cn / sn
                c *= a
                for m, n in mn_:
                    a *= c
                    c *= dn
                    dn = (n + a) / (m + a)
                    a  = c / m
                sn = copysign0(_1_0 / hypot1(c), sn)  # _signBit(sn)
                cn = c * sn
                if d and _signBit(self.kp2):  # PYCHOK no cover
                    cn, dn = dn, cn
                    sn = sn / d  # /= d chokes PyChecker
        else:
            sn = tanh(x)
            cn = dn = _1_0 / cosh(x)

        return Elliptic3Tuple(sn, cn, dn, iteration=self._iteration)

    @Property_RO
    def _sncndnBulirsch3(self):
        '''(INTERNAL) Get and cache Bulirsch' 3-tuple C{(c, d, mn_)}.
        '''
        # Bulirsch's sncndn routine, p 89.
        d, mc = 0, self.kp2
        if _signBit(mc):  # PYCHOK no cover
            d  = _1_0 - mc
            mc =  neg(mc / d)
            d  =  sqrt(d)

        a, mn, t = _1_0, [], _TolJAC
        for i in range(1, _TRIPS):  # GEOGRAPHICLIB_PANIC
            # This converges quadratically, max 6 trips
            mc = sqrt(mc)
            mn.append((a, mc))
            c = (a + mc) * _0_5
            if abs(a - mc) <= t:
                self._iteration += i  # accumulate
                break
            mc *= a
            a   = c
            t   = a * _TolJAC
        else:  # PYCHOK no cover
            raise _convergenceError(t, None, kp=self.kp, kp2=self.kp2)
        return c, d, tuple(reversed(mn))  # mn reversed!

    def _sncndnPhi(self, phi):
        '''(INTERNAL) Helper for C{.fEinv} and C{._fXf}.
        '''
        sn, cn = sincos2(phi)
        return Elliptic3Tuple(sn, cn, self.fDelta(sn, cn))

    @staticmethod
    def fRC(x, y):
        '''Degenerate symmetric integral of the first kind C{RC(x, y)}.

           @return: C{RC(x, y)}, equivalent to C{RF(x, y, y)}.

           @see: U{C{RC} definition<https://DLMF.NIST.gov/19.2.E17>} and
                 U{Carlson<https://ArXiv.org/pdf/math/9409227.pdf>}.
        '''
        return _RC(None, x, y)

    @staticmethod
    def fRD(x, y, z):
        '''Degenerate symmetric integral of the third kind C{RD(x, y, z)}.

           @return: C{RD(x, y, z)}, equivalent to C{RJ(x, y, z, z)}.

           @see: U{C{RD} definition<https://DLMF.NIST.gov/19.16.E5>} and
                 U{Carlson<https://ArXiv.org/pdf/math/9409227.pdf>}.
        '''
        return _RD(None, x, y, z)

    @staticmethod
    def fRF(x, y, *z):
        '''Symmetric or complete symmetric integral of the first kind
           C{RF(x, y, z)} respectively C{RF(x, y)}.

           @return: C{RF(x, y, z)} or C{RF(x, y)} for missing or zero B{C{z}}.

           @see: U{C{RF} definition<https://DLMF.NIST.gov/19.16.E1>} and
                 U{Carlson<https://ArXiv.org/pdf/math/9409227.pdf>}.
        '''
        return _RF3(None, x, y, *z) if z and z[0] else _RF2(None, x, y)

    @staticmethod
    def fRG(x, y, *z):
        '''Symmetric or complete symmetric integral of the second kind
           C{RG(x, y, z)} respectively C{RG(x, y)}.

           @return: C{RG(x, y, z)} or C{RG(x, y)} for missing or zero B{C{z}}.

           @see: U{C{RG} definition<https://DLMF.NIST.gov/19.16.E3>} and
                 U{Carlson<https://ArXiv.org/pdf/math/9409227.pdf>}.
        '''
        return _RG3(None, x, y, *z) if z and z[0] else (_RG2(None,x, y) * _0_5)

    @staticmethod
    def fRJ(x, y, z, p):
        '''Symmetric integral of the third kind C{RJ(x, y, z, p)}.

           @return: C{RJ(x, y, z, p)}.

           @see: U{C{RJ} definition<https://DLMF.NIST.gov/19.16.E2>} and
                 U{Carlson<https://ArXiv.org/pdf/math/9409227.pdf>}.
        '''
        return _RJ(None, x, y, z, p)


class EllipticError(_ValueError):
    '''Elliptic integral, function, convergence or other L{Elliptic} issue.
    '''
    pass


class Elliptic3Tuple(_NamedTuple):
    '''3-Tuple C{(sn, cn, dn)} all C{scalar}.
    '''
    _Names_ = ('sn',   'cn',   'dn')
    _Units_ = ( Scalar, Scalar, Scalar)


def _callError(name, *args):  # PYCHOK no cover
    '''(INTERNAL) Return an L{EllipticError}.
    '''
    n = _DOT_(Elliptic.__name__, name)
    n = _SPACE_(_invokation_, n)
    return EllipticError(NN(n, repr(args)))  # unstr


def _convergenceError(tol, where, *args, **kwds):  # PYCHOK no cover
    '''(INTERNAL) Return an L{EllipticError}.
    '''
    n = Elliptic.__name__
    if where:
        n = _DOT_(n, where.__name__)
    t = unstr(n, *args, **kwds)
    return EllipticError(_no_(_convergence_), tol, txt=t)


def _horner(S, e1, E2, E3, E4, E5, *over):
    '''(INTERNAL) Horner form for C{_RD} and C{_RJ} below.
    '''
    E22 = E2**2
    # Polynomial is <https://DLMF.NIST.gov/19.36.E2>
    # (1 - 3*E2/14 + E3/6 + 9*E2**2/88 - 3*E4/22 - 9*E2*E3/52
    #    + 3*E5/26 - E2**3/16 + 3*E3**2/40 + 3*E2*E4/20
    #    + 45*E2**2*E3/272 - 9*(E3*E4+E2*E5)/68)
    # converted to Horner form ...
    e  = e1 * 4084080
    S *= e
    S += Fsum(                             E2 * -540540,  471240).fmul(E5)
    S += Fsum(               E3 * -540540, E2 *  612612, -556920).fmul(E4)
    S += Fsum(E22 *  675675, E3 *  306306, E2 * -706860,  680680).fmul(E3)
    S += Fsum(E22 * -255255,               E2 *  417690, -875160).fmul(E2)
    return S.fadd_(4084080).fover((e * over[0]) if over else e)


def _iterations(inst, i):
    '''(INTERNAL) Aggregate iterations B{C{i}}.
    '''
    if inst:
        inst._iteration += i


def _Q(_Tol, A0, *ts):
    '''(INTERNAL) Helper for C{_RD}, C{_RF3} and C{_RJ}.
    '''
    return max(abs(A0 - t) for t in ts) / _Tol


def _RC(unused, x, y):
    '''(INTERNAL) Defined only for y != 0 and x >= 0.
    '''
    d = x - y
    if d < 0:  # catch _NaN
        # <https://DLMF.NIST.gov/19.2.E18>
        d = -d
        r = atan(sqrt(d / x)) if x > 0 else PI_2
    elif d == _0_0:  # XXX d < EPS0? or EPS02 or _EPSmin
        d, r = y, _1_0
    elif y > 0:  # <https://DLMF.NIST.gov/19.2.E19>
        r = asinh(sqrt(d / y))  # atanh(sqrt((x - y) / x))
    elif y < 0:  # <https://DLMF.NIST.gov/19.2.E20>
        r = asinh(sqrt(-x / y))  # atanh(sqrt(x / (x - y)))
    else:
        raise _callError(Elliptic.fRC.__name__, x, y)
    return r / sqrt(d)


def _RD(inst, x, y, z, *over):
    '''(INTERNAL) Carlson, eqs 2.28 - 2.34.
    '''
    A0 =  Fsum(x, y, _3_0 * z).fover(_5_0)
    T  = (A0, x, y, z)
    Q  = _Q(_TolRF, *T)
    S  =  Fsum()
    m  =  1
    for i in range(_TRIPS):
        An = T[0]
        Am = An * m
        if Q < abs(Am):  # max 7 trips
            _iterations(inst, i)
            break
        t = T[3]  # z0...n
        T, s, r = _Tsr3(T)
        S += _1_0 / ((t + r) * s[2] * m)
        m *=  4
    else:  # PYCHOK no cover
        raise _convergenceError(Q, Elliptic.fRD, x, y, z)

    x, y = _xyz(A0, -Am, x, y)
    z  = (x + y) / _3_0
    z2 =  z**2
    xy =  x * y
    S *= _3_0
    return _horner(S, Am * sqrt(An),
                   xy - _6_0 * z2,
                  (xy * _3_0 - _8_0 * z2) * z,
                  (xy - z2) * _3_0 * z2,
                   xy * z2 * z, *over)


def _3rd(a, b):
    '''(INTERNAL) Return _horner C{over} value.
    '''
    return _3_0 / (a * b)


def _RF2(inst, x, y):  # 2-arg version, z=0
    '''(INTERNAL) Carlson, eqs 2.36 - 2.38.
    '''
    a, b = sqrt(x), sqrt(y)
    if a < b:
        a, b = b, a
    for i in range(_TRIPS):
        t = _TolRG0 * a
        if abs(a - b) <= t:  # max 4 trips
            _iterations(inst, i)
            return (PI / (a + b))
        b, a = sqrt(a * b), (a + b) * _0_5
    else:  # PYCHOK no cover
        raise _convergenceError(t, Elliptic.fRF, x, y)


def _RF3(inst, x, y, z):  # 3-arg version
    '''(INTERNAL) Carlson, eqs 2.2 - 2.7.
    '''
    A0 =  Fsum(x, y, z).fover(_3_0)
    T  = (A0, x, y, z)
    Q  = _Q(_TolRF, *T)
    m  =  1
    for i in range(_TRIPS):
        An = T[0]
        Am = An * m
        if Q < abs(Am):  # max 6 trips
            _iterations(inst, i)
            break
        T, _, _ = _Tsr3(T)
        m *= 4
    else:  # PYCHOK no cover
        raise _convergenceError(Q, Elliptic.fRF, x, y, z)

    x, y = _xyz(A0, Am, x, y)
    z  = neg(x + y)
    e2 = x * y - z**2
    e3 = x * y * z
    e4 = e2**2
    # Polynomial is <https://DLMF.NIST.gov/19.36.E1>
    # (1 - E2/10 + E3/14 + E2**2/24 - 3*E2*E3/44
    #    - 5*E2**3/208 + 3*E3**2/104 + E2**2*E3/16)
    # converted to Horner form ...
    S  = Fsum(e4 * 15015, e3 * 6930, e2 * -16380,  17160).fmul(e3)
    S += Fsum(e4 * -5775,            e2 *  10010, -24024).fmul(e2)
    return S.fadd_(240240).fover(sqrt(An) * 240240)


def _RG2(inst, x, y):  # 2-args and I{doubled}
    '''(INTERNAL) Carlson, eqs 2.36 - 2.39.
    '''
    a, b = sqrt(x), sqrt(y)
    if a < b:
        a, b = b, a
    ab =  a - b  # abs(a - b)
    S  =  Fsum(_0_5 * (a + b)**2)
    m  = -1
    for i in range(_TRIPS):  # max 4 trips
        t = _TolRG0 * a
        if ab <= t:
            _iterations(inst, i)
            return S.fover((a + b) / PI_2)
        a, b = ((a + b) * _0_5), sqrt(a * b)
        ab = abs(a - b)
        S += ab**2 * m
        m *= 2
    else:  # PYCHOK no cover
        raise _convergenceError(t, Elliptic.fRG, x, y)


def _RG3(inst, x, y, z):  # 3-arg version
    '''(INTERNAL) Never called with zero B{C{z}}, see C{.fRG}.
    '''
#   if not z:
#       y, z = z, y
    rd = (x - z) * (z - y)  # - (y - z)
    if rd:  # Carlson, eq 1.7
        rd = _RD(inst, x, y, z, _3_0 * z / rd)
    xyz = x * y
    if xyz:
        xyz = sqrt(xyz / z**3)
    return Fsum(_RF3(inst, x, y, z), rd, xyz).fover(_2_0 / z)


def _RJ(inst, x, y, z, p, *over):
    '''(INTERNAL) Carlson, eqs 2.17 - 2.25.
    '''
    def _xyzp(x, y, z, p):
        return (x + p) * (y + p) * (z + p)

    A0 =  Fsum(x, y, z, _2_0 * p).fover(_5_0)
    T  = (A0, x, y, z, p)
    Q  = _Q(_TolRD, *T)
    S  =  Fsum()
    m  =  1
    Dn =  neg(_xyzp(x, y, z, -p))
    for i in range(_TRIPS):
        An = T[0]
        Am = An * m
        if Q < abs(Am):  # max 7 trips
            _iterations(inst, i)
            break
        T, s, _ = _Tsr3(T)
        d = _xyzp(*s)
        if Dn:
            rc  = _RC(inst, _1_0, Dn / d**2 + _1_0)
            Dn *= _1_64th
        else:
            rc  = _1_0  # == _RC(None, _1_0, _1_0)
        S += rc / (d * m)
        m *= 4
    else:  # PYCHOK no cover
        raise _convergenceError(Q, Elliptic.fRJ, x, y, z, p)

    x, y, z = _xyz(A0, Am, x, y, z)
    xyz =  x * y * z
    p   = -Fsum(x, y, z).fover(_2_0)
    p2  =  p**2
    p3  =  p**3
    E2  =  Fsum(x * y, x * z, y * z, -p2 * _3_0)
    E2p =  E2 * p
    S  *= _6_0
    return _horner(S, Am * sqrt(An), E2,
                   Fsum(p3 * _4_0, xyz, E2p * _2_0),
                   Fsum(p3 * _3_0, E2p, xyz * _2_0).fmul(p),
                   xyz * p2, *over)


def _Tsr3(T):
    '''(INTERNAL) Helper for C{_RD}, C{_RF3} and C{_RJ}.
    '''
    s = map2(sqrt, T[1:])  # sqrt(x), srqt(y), sqrt(z) [, sqrt(p)]
    r = fdot(s[:3], s[1], s[2], s[0])  # sqrt(x) * sqrt(y) + ...
    T = tuple((t + r) * _0_25 for t in T)  # T[:] = ...
    return T, s, r


def _xyz(A0, Am, *xyz):
    '''(INTERNAL) Rescale any C{xys}.
    '''
    return tuple((A0 - x) / Am for x in xyz)

# **) MIT License
#
# Copyright (C) 2016-2022 -- mrJean1 at Gmail -- All Rights Reserved.
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
