
# -*- coding: utf-8 -*-

u'''Elliptic integrals and functions transcribed from I{Charles Karney}'s
C++ class U{EllipticFunction
<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1EllipticFunction.html>}
into pure Python class L{Elliptic}.

Python method names follow the C++ member functions, except:

 - member functions I{without arguments} are mapped to Python properties
   prefixed with C{"c"}, for example C{E()} is property C{cE},

 - member functions with 1 or 3 arguments are renamed to Python methods
   starting with an C{"f"}, example C{E(psi)} to C{fE(psi)} and C{E(sn,
   cn, dn)} to C{fE(sn, cn, dn)},

 - other Python method names conventionally start with a lower-case
   letter or an underscore if private.

Following is a copy of I{Karney}'s U{EllipticFunction.hpp
<https://GeographicLib.SourceForge.io/html/EllipticFunction_8hpp_source.html>}
file C{Header}.

Copyright (C) U{Charles Karney<mailto:Charles@Karney.com>} (2008-2017)
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
from __future__ import division

from pygeodesy.basics import copysign, map2, neg
from pygeodesy.errors import _ValueError
from pygeodesy.fmath import fdot, fmean_, Fsum, fsum_, hypot1
from pygeodesy.interns import EPS, INF, NN, PI, PI_2, PI_4, \
                             _EPStol as _TolJAC, _convergence_, \
                             _no_, _SPACE_, _0_0, _0_125, _0_25, \
                             _0_5, _1_0, _2_0, _3_0, _4_0, _5_0, \
                             _6_0, _8_0, _180_0, _360_0
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import _Named, _NamedTuple
from pygeodesy.props import Property_RO, property_RO, _update_all
# from pygeodesy.streprs import unstr
from pygeodesy.units import Scalar, Scalar_
from pygeodesy.utily import sincos2, sincos2d

from math import asinh, atan, atan2, ceil, cosh, floor, sin, \
                 sqrt, tanh

__all__ = _ALL_LAZY.elliptic
__version__ = '21.04.15'

_TolRD  =  pow(EPS * 0.002, _0_125)
_TolRF  =  pow(EPS * 0.030, _0_125)
_TolRG0 = _TolJAC  * 2.7
_TRIPS  =  15  # Max depth, 7 might be enough, by .etm


class EllipticError(_ValueError):
    '''Elliptic integral, function, convergence or other L{Elliptic} issue.
    '''
    pass


class Elliptic3Tuple(_NamedTuple):
    '''3-Tuple C{(sn, cn, dn)} all C{scalar}.
    '''
    _Names_ = ('sn',   'cn',   'dn')
    _Units_ = ( Scalar, Scalar, Scalar)


class Elliptic(_Named):
    '''Elliptic integrals and functions.

       @see: I{Karney}'s U{Detailed Description<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1EllipticFunction.html#details>}.
    '''
    _alpha2    = 0
    _alphap2   = 0
    _eps       = EPS
    _iteration = None  # only .fEinv and .sncndn
    _k2        = 0
    _kp2       = 0

    _cD  =  INF
    _cE  = _1_0
    _cG  =  INF
    _cH  =  INF
    _cK  =  INF
    _cKE =  INF
    _cPi =  INF

    def __init__(self, k2=0, alpha2=0, kp2=None, alphap2=None):
        '''Constructor, specifying the C{modulus} and C{parameter}.

           @kwarg k2: Modulus squared (C{float}, 0 <= k^2 <= 1).
           @kwarg alpha2: Parameter squared (C{float}, 0 <= α^2 <= 1).
           @kwarg kp2: Complementary modulus squared (C{float}, k'^2 >= 0).
           @kwarg alphap2: Complementary parameter squared (C{float}, α'^2 >= 0).

           @see: Method L{reset} for further details.

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
        return self._cD

    @Property_RO
    def cE(self):
        '''Get the complete integral of the second kind C{E(k)}
           (C{float}), U{defined<https://DLMF.NIST.gov/19.2.E5>}.
        '''
        return self._cE

    @Property_RO
    def cG(self):
        '''Get Legendre's complete geodesic longitude integral
           C{G(α^2, k)} (C{float}).
        '''
        return self._cG

    @Property_RO
    def cH(self):
        '''Get Cayley's complete geodesic longitude difference integral
           C{H(α^2, k)} (C{float}).
        '''
        return self._cH

    @Property_RO
    def cK(self):
        '''Get the complete integral of the first kind C{K(k)}
           (C{float}), U{defined<https://DLMF.NIST.gov/19.2.E4>}.
        '''
        return self._cK

    @Property_RO
    def cKE(self):
        '''Get the difference between the complete integrals of the
           first and second kinds, C{K(k) − E(k)} (C{float}).
        '''
        return self._cKE

    @Property_RO
    def cPi(self):
        '''Get the complete integral of the third kind C{Pi(α^2, k)}
           (C{float}), U{defined<https://DLMF.NIST.gov/19.2.E7>}.
        '''
        return self._cPi

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
        t = atan2(-stau, -ctau) if ctau < 0 else atan2(stau, ctau)
        return self.fEinv(t * self._cE / PI_2) - t

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
        if None in (cn, dn):
            t = self.classname + '.delta' + fX.__name__[1:]
            raise _invokationError(t, sn, cn, dn)

        if cn < 0:
            cn, sn = -cn, neg(sn)
        return fX(sn, cn, dn) * PI_2 / cX - atan2(sn, cn)

    @Property_RO
    def eps(self):
        '''Get epsilon (C{float}).
        '''
        return self._eps

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
            return abs(sn) * sn**2 * _RD_3(cn**2, dn**2, _1_0)

        return self._fXf(phi_or_sn, cn, dn, self.cD,
                                            self.deltaD, _fD)

    def fDelta(self, sn, cn):
        '''The C{Delta} amplitude function.

           @arg sn: sin(φ).
           @arg cn: cos(φ).

           @return: sqrt(1 − k2 sin(2φ)) (C{float}).
        '''
        k2 = self.k2
        return sqrt((_1_0 - k2 * sn**2) if k2 < 0 else
                (self.kp2 + k2 * cn**2))

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
            sn2, cn2, dn2 = sn**2, cn**2, dn**2
            kp2, k2 = self.kp2, self.k2
            if k2 <= 0:  # Carlson, eq. 4.6, <https://DLMF.NIST.gov/19.25.E9>
                ei = _RF(cn2, dn2, _1_0) - k2 * sn2 * _RD_3(cn2, dn2, _1_0)
            elif kp2 >= 0:  # <https://DLMF.NIST.gov/19.25.E10>
                ei = fsum_(kp2 * _RF(cn2, dn2, _1_0),
                           kp2 * k2 * sn2 * _RD_3(cn2, _1_0, dn2),
                           k2 * abs(cn) / dn)
            else:  # <https://DLMF.NIST.gov/19.25.E11>
                ei = dn / abs(cn) - kp2 * sn2 * _RD_3(dn2, _1_0, cn2)
            return ei * abs(sn)

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
        else:
            n    = ceil(deg / _360_0 - _0_5)
            deg -= n * _360_0
            e    = n * _4_0 * self.cE
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
        x -= E2 * n  # x now in [-ec, ec)
        # linear approximation
        phi = PI * x / E2  # phi in [-pi/2, pi/2)
        Phi = Fsum(phi)
        # first order correction
        phi = Phi.fsum_(-self.eps * sin(2 * phi) * _0_5)
        # For kp2 close to zero use asin(x/.cE) or J. P. Boyd,
        # Applied Math. and Computation 218, 7005-7013 (2012)
        # <https://DOI.org/10.1016/j.amc.2011.12.021>
        for self._iteration in range(1, _TRIPS):  # GEOGRAPHICLIB_PANIC
            sn, cn, dn = self._sncndn3(phi)
            phi, e = Phi.fsum2_((x - self.fE(sn, cn, dn)) / dn)
            if abs(e) < _TolJAC:
                if n:
                    phi = Phi.fsum_(n * PI)
                return phi
        raise _convergenceError(self.fEinv, x)

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
            return abs(sn) * _RF(cn**2, dn**2, _1_0)

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

           @note: Legendre expresses the longitude of a point on
                  the geodesic in terms of this combination of
                  elliptic integrals in U{Exercices de Calcul
                  Intégral, Vol 1 (1811), p 181<https://Books.
                  Google.com/books?id=riIOAAAAQAAJ&pg=PA181>}.

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
        return self._fXa(phi_or_sn, cn, dn, self.alphap2,
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
            cn2, sn2, dn2 = cn**2, sn**2, dn**2
            return abs(sn) * (_RF(cn2, dn2, _1_0) + aX * sn2 *
                            _RJ_3(cn2, dn2, _1_0, cn2 + self.alphap2 * sn2))

        return self._fXf(phi_or_sn, cn, dn, cX, deltaX, _fX)

    def _fXf(self, phi_or_sn, cn, dn, cX, deltaX, fX):
        '''(INTERNAL) Helper for C{f.D}, C{.fE}, C{.fF} and C{._fXa}.
        '''
        phi = sn = phi_or_sn
        if cn is dn is None:  # fX(phi) call
            sn, cn, dn = self._sncndn3(phi)
            if abs(phi) >= PI:
                return (deltaX(sn, cn, dn) + phi) * cX / PI_2
            # fall through
        elif None in (cn, dn):
            t = self.classname + '.f' + deltaX.__name__[5:]
            raise _invokationError(t, sn, cn, dn)

        if cn > 0:  # enforce usual trig-like symmetries
            xi = fX(sn, cn, dn)
        elif cn < 0:
            xi = _2_0 * cX - fX(sn, cn, dn)
        else:
            xi = cX
        return copysign(xi, sn)

    @property_RO
    def iteration(self):
        '''Get the most recent C{Elliptic.fEinv} or C{Elliptic.sncndn}
           iteration number (C{int}) or C{None} if not available/applicable.
        '''
        return self._iteration

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
        '''Reset the modulus and parameter.

           @kwarg k2: modulus squared (C{float}, -INF <= k^2<= 1).
           @kwarg alpha2: parameter (C{float}, -INF <= α^2 <= 1).
           @kwarg kp2: complementary modulus squared (C{float}, k'^2 >= 0).
           @kwarg alphap2: complementary parameter squared (C{float}, α'^2 >= 0).

           @raise EllipticError: Invalid B{C{k2}}, B{C{alpha2}}, B{C{kp2}}
                                 or B{C{alphap2}} or no convergence.

           @note: The arguments must satisfy C{B{k2} + B{kp2} = 1} and
                  C{B{alpha2} + B{alphap2} = 1}.  No checking is done
                  that these conditions are met to enable accuracy to
                  be maintained, e.g., when C{k} is very close to unity.
        '''
        _update_all(self)

        self._k2 = k2 = Scalar_(k2=k2, Error=EllipticError, low=None, high=_1_0)

        self._alpha2 = alpha2 = Scalar_(alpha2=alpha2, Error=EllipticError, low=None, high=_1_0)

        self._kp2 = kp2 = Scalar_(kp2=((_1_0 - k2) if kp2 is None else kp2), Error=EllipticError)

        self._alphap2 = alphap2 = Scalar_(alphap2=((_1_0 - alpha2) if alphap2 is None else alphap2),
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
        if k2:
            if kp2:
                self._eps = k2 / (sqrt(kp2) + _1_0)**2
                # D(k) = (K(k) - E(k))/k2, Carlson eq.4.3
                # <https://DLMF.NIST.gov/19.25.E1>
                self._cD = _RD_3(_0_0, kp2, _1_0)
                self._cKE = k2 * self.cD
                # Complete elliptic integral E(k), Carlson eq. 4.2
                # <https://DLMF.NIST.gov/19.25.E1>
                self._cE = _RG_(kp2, _1_0) * _2_0
                # Complete elliptic integral K(k), Carlson eq. 4.1
                # <https://DLMF.NIST.gov/19.25.E1>
                self._cK = _RF_(kp2, _1_0)
            else:
                self._eps =  k2
                self._cD  =  self._cK = self._cKE = INF
                self._cE  = _1_0
        else:
            self._eps =  EPS  # delattr(self, '_eps')
            self._cD  =  PI_4
            self._cE  =  self._cK = PI_2
            self._cKE = _0_0  # k2 * self._cD

        if alpha2:
            if alphap2:
                # <https://DLMF.NIST.gov/19.25.E2>
                if kp2:
                    rj_3 = _RJ_3(_0_0, kp2, _1_0, alphap2)
                    # G(alpha2, k)
                    self._cG = self.cK + (alpha2 - k2) * rj_3
                    # H(alpha2, k)
                    self._cH = self.cK - alphap2 * rj_3
                    # Pi(alpha2, k)
                    self._cPi = self.cK + alpha2 * rj_3
                else:
                    self._cG = self._cH = _RC(_1_0, alphap2)
                    self._cPi = INF  # XXX or NAN?
            else:
                self._cG = self._cH = self._cPi = INF  # XXX or NAN?
        else:
            self._cG  = self.cE
            self._cPi = self.cK
            # cH = cK - cD but this involves large cancellations
            # if k2 is close to 1.  So write (for alpha2 = 0)
            #   cH = int(cos(phi)**2/sqrt(1-k2*sin(phi)**2),phi,0,pi/2)
            #      = 1/sqrt(1-k2) * int(sin(phi)**2/sqrt(1-k2/kp2*sin(phi)**2,...)
            #      = 1/kp * D(i*k/kp)
            # and use D(k) = RD(0, kp2, 1) / 3
            # so cH = 1/kp * RD(0, 1/kp2, 1) / 3
            #       = kp2 * RD(0, 1, kp2) / 3
            # using <https://DLMF.NIST.gov/19.20.E18>
            # Equivalently
            #   RF(x, 1) - RD(0, x, 1)/3 = x * RD(0, 1, x)/3 for x > 0
            # For k2 = 1 and alpha2 = 0, we have
            #   cH = int(cos(phi),...) = 1
            self._cH = kp2 * _RD_3(_0_0, _1_0, kp2) if kp2 else _1_0

#       self._iteration = 0

    def sncndn(self, x):  # PYCHOK x used!
        '''The Jacobi elliptic function.

           @arg x: The argument (C{float}).

           @return: An L{Elliptic3Tuple}C{(sn, cn, dn)} with
                    C{*n(B{x}, k)}.

           @raise EllipticError: No convergence.
        '''
        # Bulirsch's sncndn routine, p 89.
        mc = self.kp2
        if mc:  # never negative ...
            if mc < 0:  # PYCHOK no cover
                d  = _1_0 - mc
                mc = neg(mc / d)  # /= -d chokes PyChecker
                d  = sqrt(d)
                x *= d
            else:
                d = 0
            a, mn = _1_0, []
            for self._iteration in range(1, _TRIPS):  # GEOGRAPHICLIB_PANIC
                # This converges quadratically, max 6 trips
                mc = sqrt(mc)  # XXX sqrt0?
                mn.append((a, mc))
                c = (a + mc) * _0_5
                if abs(a - mc) <= (_TolJAC * a):
                    break
                mc *= a
                a = c
            else:
                raise _convergenceError(self.sncndn, x)
            x *=  c
            dn = _1_0
            sn, cn = sincos2(x)
            if sn:
                a = cn / sn
                c *= a
                while mn:
                    m, n = mn.pop()
                    a *= c
                    c *= dn
                    dn = (n + a) / (m + a)
                    a  = c / m
                sn = copysign(_1_0 / hypot1(c), sn)
                cn = c * sn
                if d:  # PYCHOK no cover
                    cn, dn = dn, cn
                    sn = sn / d  # /= d chokes PyChecker
        else:
            self._iteration = 0
            sn = tanh(x)
            cn = dn = _1_0 / cosh(x)

        r = Elliptic3Tuple(sn, cn, dn)
        r._iteration = self._iteration
        return r

    def _sncndn3(self, phi):
        '''(INTERNAL) Helper for C{.fEinv} and C{._fXf}.
        '''
        sn, cn = sincos2(phi)
        return Elliptic3Tuple(sn, cn, self.fDelta(sn, cn))


def _convergenceError(where, *args):  # PYCHOK no cover
    '''(INTERNAL) Return an L{EllipticError}.
    '''
    t = _SPACE_(where.__name__, repr(args))
    return EllipticError(_no_(_convergence_), txt=t)  # unstr


def _horner(e0, e1, e2, e3, e4, e5):
    '''(INTERNAL) Horner form for C{_RD} and C{_RJ} below.
    '''
    # Polynomial is <https://DLMF.NIST.gov/19.36.E2>
    # (1 - 3*E2/14 + E3/6 + 9*E2**2/88 - 3*E4/22 - 9*E2*E3/52 + 3*E5/26
    #    - E2**3/16 + 3*E3**2/40 + 3*E2*E4/20 + 45*E2**2*E3/272
    #    - 9*(E3*E4+E2*E5)/68)
    # converted to Horner form ...
    H  = Fsum(471240,      -540540 * e2) * e5
    H += Fsum(612612 * e2, -540540 * e3,    -556920) * e4
    H += Fsum(306306 * e3,  675675 * e2**2, -706860  * e2, 680680) * e3
    H += Fsum(417690 * e2, -255255 * e2**2, -875160) * e2
    e  = 4084080 * e1
    return H.fsum_(4084080, e * e0) / e


def _invokationError(name, *args):  # PYCHOK no cover
    '''(INTERNAL) Return an L{EllipticError}.
    '''
    return EllipticError(_SPACE_('invokation', NN(name, repr(args))))  # unstr


def _Q(A, T, tol):
    '''(INTERNAL) Helper for C{_RD}, C{_RF} and C{_RJ}.
    '''
    return max(abs(A - t) for t in T[1:]) / tol


def _RC(x, y):  # used by testElliptic.py
    '''Degenerate symmetric integral of the first kind C{_RC(x, y)}.

       @return: C{_RC(x, y) = _RF(x, y, y)}.

       @see: U{C{_RC} definition<https://DLMF.NIST.gov/19.2.E17>}.
    '''
    # Defined only for y != 0 and x >= 0.
    d = x - y
    if d < 0:  # catch _NaN
        # <https://DLMF.NIST.gov/19.2.E18>
        d = -d
        r = atan(sqrt(d / x)) if x > 0 else PI_2
    elif d == _0_0:  # XXX d < EPS0?
        r, d = _1_0, y
    elif y > 0:  # <https://DLMF.NIST.gov/19.2.E19>
        r = asinh(sqrt(d / y))  # atanh(sqrt((x - y) / x))
    elif y < 0:  # <https://DLMF.NIST.gov/19.2.E20>
        r = asinh(sqrt(-x / y))  # atanh(sqrt(x / (x - y)))
    else:
        raise _invokationError(_RC.__name__, x, y)
    return r / sqrt(d)


def _RD_3(x, y, z):
    '''Degenerate symmetric integral of the third kind C{_RD(x, y, z) / 3}.
    '''
    return _RD(x, y, z) / _3_0


def _RD(x, y, z):  # used by testElliptic.py
    '''Degenerate symmetric integral of the third kind C{_RD(x, y, z)}.

       @return: C{_RD(x, y, z) = _RJ(x, y, z, z)}.

       @see: U{C{_RD} definition<https://DLMF.NIST.gov/19.16.E5>} and
             U{Carlson<https://ArXiv.org/pdf/math/9409227.pdf>}.
    '''
    # Carlson, eqs 2.28 - 2.34
    m = _1_0
    A =  fsum_(x, y, _3_0 * z) / _5_0
    T = (A, x, y, z)
    Q = _Q(A, T, _TolRD)
    S =  Fsum()
    for _ in range(_TRIPS):
        An = T[0]
        if Q < abs(m * An):  # max 7 trips
            break
        t = T[3]  # z0
        r, s, T = _rsT3(T)
        S += _1_0 / (m * s[2] * (t + r))
        m *= _4_0
    else:
        raise _convergenceError(_RD, x, y, z)

    m *= An
    x = (x - A) / m
    y = (y - A) / m
    z = (x + y) / _3_0
    z2 = z**2
    xy = x * y
    return _horner(3 * S.fsum(), m * sqrt(An),
                   xy - _6_0 * z2,
                  (xy * _3_0 - _8_0 * z2) * z,
                  (xy - z2) * _3_0 * z2,
                   xy * z2 * z)


def _RF_(x, y):
    '''Symmetric integral of the first kind C{_RF_(x, y)}.

       @return: C{_RF(x, y)}.

       @see: U{C{_RF} definition<https://DLMF.NIST.gov/19.16.E1>}.
    '''
    # Carlson, eqs 2.36 - 2.38
    a, b = sqrt(x), sqrt(y)
    if a < b:
        a, b = b, a
    for _ in range(_TRIPS):
        if abs(a - b) <= (_TolRG0 * a):  # max 4 trips
            return PI / (a + b)
        b, a = sqrt(a * b), (a + b) * _0_5

    raise _convergenceError(_RF_, x, y)


def _RF(x, y, z):  # used by testElliptic.py
    '''Symmetric integral of the first kind C{_RF(x, y, z)}.

       @return: C{_RF(x, y, z)}.

       @see: U{C{_RF} definition<https://DLMF.NIST.gov/19.16.E1>} and
             U{Carlson<https://ArXiv.org/pdf/math/9409227.pdf>}.
    '''
    # Carlson, eqs 2.2 - 2.7
    m = _1_0
    A =  fmean_(x, y, z)
    T = (A, x, y, z)
    Q = _Q(A, T, _TolRF)
    for _ in range(_TRIPS):
        An = T[0]
        if Q < abs(m * An):  # max 6 trips
            break
        _, _, T = _rsT3(T)
        m *= _4_0
    else:
        raise _convergenceError(_RF, x, y, z)

    m *= An
    x = (A - x) / m
    y = (A - y) / m
    z = neg(x + y)

    e2 = x * y - z**2
    e3 = x * y * z
    # Polynomial is <https://DLMF.NIST.gov/19.36.E1>
    # (1 - E2/10 + E3/14 + E2**2/24 - 3*E2*E3/44
    #    - 5*E2**3/208 + 3*E3**2/104 + E2**2*E3/16)
    # converted to Horner form ...
    H  = Fsum( 6930 * e3, 15015 * e2**2, -16380  * e2, 17160) * e3
    H += Fsum(10010 * e2, -5775 * e2**2, -24024) * e2
    return H.fsum_(240240) / (240240 * sqrt(An))


def _RG_(x, y):
    '''Symmetric integral of the second kind C{_RG_(x, y)}.

       @return: C{_RG(x, y)}.

       @see: U{C{_RG} definition<https://DLMF.NIST.gov/19.16.E3>} and
             U{Carlson<https://ArXiv.org/pdf/math/9409227.pdf>}.
    '''
    # Carlson, eqs 2.36 - 2.39
    a, b = sqrt(x), sqrt(y)
    if a < b:
        a, b = b, a
    m = _0_5
    S =  Fsum(_0_25 * (a + b)**2)
    for _ in range(_TRIPS):  # max 4 trips
        if abs(a - b) <= (_TolRG0 * a):
            S *= PI_2 / (a + b)
            return S.fsum()
        b, a = sqrt(a * b), (a + b) * _0_5
        S -=  m * (a - b)**2
        m *= _2_0

    raise _convergenceError(_RG_, x, y)


def _RG(x, y, z):  # used by testElliptic.py
    '''Symmetric integral of the second kind C{_RG(x, y, z)}.

       @return: C{_RG(x, y, z)}.

       @see: U{C{_RG} definition<https://DLMF.NIST.gov/19.16.E3>} and
             U{Carlson<https://ArXiv.org/pdf/math/9409227.pdf>}.
    '''
    if not z:
        y, z = z, y
    # Carlson, eq 1.7
    return fsum_(_RF(x, y, z) * z,
                 _RD_3(x, y, z) * (x - z) * (z - y),
                  sqrt(x * y / z)) * _0_5


def _RJ_3(x, y, z, p):
    '''Symmetric integral of the third kind C{_RJ(x, y, z, p) / 3}.
    '''
    return _RJ(x, y, z, p) / _3_0


def _RJ(x, y, z, p):  # used by testElliptic.py
    '''Symmetric integral of the third kind C{_RJ(x, y, z, p)}.

       @return: C{_RJ(x, y, z, p)}.

       @see: U{C{_RJ} definition<https://DLMF.NIST.gov/19.16.E2>} and
             U{Carlson<https://ArXiv.org/pdf/math/9409227.pdf>}.
    '''
    def _xyzp(x, y, z, p):
        return (x + p) * (y + p) * (z + p)

    # Carlson, eqs 2.17 - 2.25
    m =  m3 = _1_0
    D =  neg(_xyzp(x, y, z, -p))
    A =  fsum_(x, y, z, _2_0 * p) / _5_0
    T = (A, x, y, z, p)
    Q = _Q(A, T, _TolRD)
    S =  Fsum()
    for _ in range(_TRIPS):
        An = T[0]
        if Q < abs(m * An):  # max 7 trips
            break
        _, s, T = _rsT3(T)
        d  = _xyzp(*s)
        e  =  D / (m3 * d**2)
        S += _RC(_1_0, _1_0 + e) / (m * d)
        m *= _4_0
        m3 *= 64
    else:
        raise _convergenceError(_RJ, x, y, z, p)

    m *= An
    x = (A - x) / m
    y = (A - y) / m
    z = (A - z) / m
    xyz = x * y * z
    p  = neg(x + y + z) * _0_5
    p2 = p**2
    p3 = p * p2

    e2 = fsum_(x * y, x * z, y * z, -p2 * _3_0)
    return _horner(_6_0 * S.fsum(), m * sqrt(An), e2,
                   fsum_(xyz, _2_0 * p * e2, _4_0 * p3),
                   fsum_(xyz * _2_0, p * e2, _3_0 * p3) * p,
                   p2 * xyz)


def _rsT3(T):
    '''(INTERNAL) Helper for C{_RD}, C{_RF} and C{_RJ}.
    '''
    s = map2(sqrt, T[1:])
    r = fdot(s[:3], s[1], s[2], s[0])
    T = tuple((t + r) * _0_25 for t in T)
    return r, s, T

# **) MIT License
#
# Copyright (C) 2016-2021 -- mrJean1 at Gmail -- All Rights Reserved.
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
