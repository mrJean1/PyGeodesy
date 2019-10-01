
# -*- coding: utf-8 -*-

u'''Elliptic integrals and functions transcribed from I{Charles Karney's}
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

Following is a copy of Karney's U{EllipticFunction.hpp
<https://GeographicLib.SourceForge.io/html/EllipticFunction_8hpp_source.html>}
file C{Header}.

Copyright (C) U{Charles Karney<mailto:Charles@Karney.com>} (2008-2017)
and licensed under the MIT/X11 License.  For more information, see the
U{GeographicLib<https://GeographicLib.SourceForge.io/>} documentation.

B{Elliptic integrals and functions.}

This provides the elliptic functions and integrals needed for
C{Ellipsoid}, C{GeodesicExact}, and C{TransverseMercatorExact}.  Two
categories of function are provided:

 - functions to compute U{symmetric elliptic integrals
   <https://DLMF.NIST.gov/19.16.i>}

 - methods to compute U{Legrendre's elliptic integrals
   <https://DLMF.NIST.gov/19.2.ii>} and the U{Jacobi elliptic
   functions<https://DLMF.NIST.gov/22.2>}.

In the latter case, an object is constructed giving the modulus
C{k} (and optionally the parameter C{alpha2}.  The modulus is always
passed as its square which allows C{k} to be pure imaginary.
(Confusingly, Abramowitz and Stegun call C{m = k**2} the "parameter"
and C{n = alpha**2} the "characteristic".)

In geodesic applications, it is convenient to separate the incomplete
integrals into secular and periodic components, e.g.

I{C{E(phi, k) = (2 E(k) / pi) [ phi + delta E(phi, k) ]}}

where I{C{delta E(phi, k)}} is an odd periodic function with
period I{C{pi}}.

The computation of the elliptic integrals uses the algorithms given
in U{B. C. Carlson, Computation of real or complex elliptic integrals
<https://DOI.org/10.1007/BF02198293>}, Numerical Algorithms 10, 13--26
(1995) with the additional optimizations given U{here
<https://DLMF.NIST.gov/19.36.i>}.

The computation of the Jacobi elliptic functions uses the algorithm
given in U{R. Bulirsch, Numerical Calculation of Elliptic Integrals
and Elliptic Functions <https://DOI.org/10.1007/BF01397975>},
Numerische Mathematik 7, 78--90 (1965).

The notation follows U{NIST Digital Library of Mathematical Functions
<https://DLMF.NIST.gov>} chapters U{19<https://DLMF.NIST.gov/19>} and
U{22<https://DLMF.NIST.gov/22>}.
'''

from pygeodesy.fmath import EPS, fdot, Fsum, fsum_, hypot1, \
                            INF as _INF, map2
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import _Named
from pygeodesy.utily import PI, PI_2, PI_4, property_RO, \
                            sincos2, sincos2d

from math import asinh, atan, atan2, ceil, copysign, cosh, \
                 floor, sin, sqrt, tanh

__all__ = _ALL_LAZY.elliptic
__version__ = '19.07.14'

_tolJAC = sqrt(EPS * 0.01)
_tolRD  =  pow(EPS * 0.002, 0.125)
_tolRF  =  pow(EPS * 0.030, 0.125)
_tolRG0 = _tolJAC * 2.7
_TRIPS  =  13  # Max depth for sncndn, etc, 5-7 might be enough


class EllipticError(ValueError):
    '''Elliptic integral, function, convergence or other L{Elliptic} issue.
    '''
    pass


class Elliptic(_Named):
    '''Elliptic integrals and functions.

       @see: Karney's U{Detailed Description<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1EllipticFunction.html#details>}.
    '''
    _alpha2  = 0
    _alphap2 = 0
    _eps     = EPS
    _k2      = 0
    _kp2     = 0
    _trips_  = _TRIPS

    _Dc      = _INF
    _Ec      = 1.0
    _Gc      = _INF
    _Hc      = _INF
    _Kc      = _INF
    _KEc     = _INF
    _Pic     = _INF

    def __init__(self, k2=0, alpha2=0, kp2=None, alphap2=None):
        '''Constructor, specifying the C{modulus} and C{parameter}.

           @keyword k2: Modulus squared (C{scalar} k^2 <= 1).
           @keyword alpha2: Parameter squared (C{scalar} α^2 <= 1).
           @keyword kp2: Complementary modulus squared (C{float} k'^2 >= 0).
           @keyword alphap2: Complementary parameter α'2 (C{float} α'^2 >= 0).

           @note: If only elliptic integrals of the first and
                  second kinds are needed, then set B{C{alpha2}} = 0
                  (the default value); in this case, we have
                  Π(φ, 0, k) = F(φ, k), G(φ, 0, k) = E(φ, k),
                  and H(φ, 0, k) = F(φ, k) - D(φ, k).

           @see: Method L{reset}.
        '''
        self.reset(k2, alpha2=alpha2, kp2=kp2, alphap2=alphap2)

    @property_RO
    def alpha2(self):
        '''Get the parameter B{C{alpha2}} (C{float}).
        '''
        return self._alpha2

    @property_RO
    def alphap2(self):
        '''Get the complementary parameter B{C{alphap2}} (C{float}).
        '''
        return self._alphap2

    @property_RO
    def cD(self):
        '''Get Jahnke's complete integral C{D(k)},
           U{defined<https://DLMF.NIST.gov/19.2.E6>}.
        '''
        return self._Dc

    @property_RO
    def cE(self):
        '''Get the complete integral of the second kind C{E(k)},
           U{defined<https://DLMF.NIST.gov/19.2.E5>}.
        '''
        return self._Ec

    @property_RO
    def cG(self):
        '''Get Legendre's complete geodesic longitude integral C{G(alpha2, k)}.
        '''
        return self._Gc

    @property_RO
    def cH(self):
        '''Get Cayley's complete geodesic longitude difference integral C{H(alpha2, k)}.
        '''
        return self._Hc

    @property_RO
    def cK(self):
        '''Get the complete integral of the first kind C{K(k)},
           U{defined<https://DLMF.NIST.gov/19.2.E4>}.
        '''
        return self._Kc

    @property_RO
    def cKE(self):
        '''Get the difference between the complete integrals of the
           first and second kinds, C{K(k) − E(k)}.
        '''
        return self._KEc

    @property_RO
    def cPi(self):
        '''Get the complete integral of the third kind C{Pi(alpha2, k)},
           U{defined<https://DLMF.NIST.gov/19.2.E7>}.
        '''
        return self._Pic

    def deltaD(self, sn, cn, dn):
        '''The periodic Jahnke's incomplete elliptic integral.

           @param sn: sinφ.
           @param cn: cosφ.
           @param dn: sqrt(1 − k2 sin2φ).

           @return: Periodic function π D(φ, k) / (2 D(k)) - φ.
        '''
        return self._deltaX(sn, cn, dn, self.cD, self.fD)

    def deltaE(self, sn, cn, dn):
        '''The periodic incomplete integral of the second kind.

           @param sn: sinφ.
           @param cn: cosφ.
           @param dn: sqrt(1 − k2 sin2φ).

           @return: Periodic function 𝜋π E(φ, k) / (2 E(k)) - φ.
        '''
        return self._deltaX(sn, cn, dn, self.cE, self.fE)

    def deltaEinv(self, stau, ctau):
        '''The periodic inverse of the incomplete integral of the second kind.

           @param stau: sinτ
           @param ctau: cosτ

           @return: Periodic function E^−1(τ (2 E(k)/π), k) - τ.
        '''
        # Function is periodic with period pi
        t = atan2(-stau, -ctau) if ctau < 0 else atan2(stau, ctau)
        return self.fEinv(t * self._Ec / PI_2) - t

    def deltaF(self, sn, cn, dn):
        '''The periodic incomplete integral of the first kind.

           @param sn: sinφ.
           @param cn: cosφ.
           @param dn: sqrt(1 − k2 sin2φ).

           @return: Periodic function π F(φ, k) / (2 K(k)) - φ.
        '''
        return self._deltaX(sn, cn, dn, self.cK, self.fF)

    def deltaG(self, sn, cn, dn):
        '''Legendre's periodic geodesic longitude integral.

           @param sn: sinφ.
           @param cn: cosφ.
           @param dn: sqrt(1 − k2 sin2φ).

           @return: Periodic function π G(φ, k) / (2 G(k)) - φ.
        '''
        return self._deltaX(sn, cn, dn, self.cG, self.fG)

    def deltaH(self, sn, cn, dn):
        '''Cayley's periodic geodesic longitude difference integral.

           @param sn: sinφ.
           @param cn: cosφ.
           @param dn: sqrt(1 − k2 sin2φ).

           @return: Periodic function π H(φ, k) / (2 H(k)) - φ.
        '''
        return self._deltaX(sn, cn, dn, self.cH, self.fH)

    def deltaPi(self, sn, cn, dn):
        '''The periodic incomplete integral of the third kind.

           @param sn: sinφ.
           @param cn: cosφ.
           @param dn: sqrt(1 − k2 sin2φ).

           @return: Periodic function π Π(φ, α2, k) / (2 Π(α2, k)) - φ.
        '''
        return self._deltaX(sn, cn, dn, self.cPi, self.fPi)

    def _deltaX(self, sn, cn, dn, cX, fX):
        '''(INTERNAL) Helper for C{.deltaD} thru C{.deltaPi}.
        '''
        if cn < 0:
            cn, sn = -cn, -sn
        return fX(sn, cn, dn) * PI_2 / cX - atan2(sn, cn)

    def fD(self, phi_or_sn, cn=None, dn=None):
        '''Jahnke's incomplete elliptic integral in terms
           of Jacobi elliptic functions.

           @param phi_or_sn: φ or sinφ.
           @keyword cn: cosφ.
           @keyword dn: sqrt(1 − k2 sin2φ).

           @return: D(φ, k) as though φ ∈ (−π, π],
                    U{defined<https://DLMF.NIST.gov/19.2.E4>}.
        '''
        def _fD(sn, cn, dn):
            return abs(sn) * sn**2 * _RD_3(cn**2, dn**2, 1)

        return self._fXf(phi_or_sn, cn, dn, self.cD,
                                            self.deltaD, _fD)

    def fDelta(self, sn, cn):
        '''The C{Delta} amplitude function.

           @param sn: sinφ.
           @param cn: cosφ.

           @return: C{sqrt(1 − k2 sin2φ)}.
        '''
        k2 = self._k2
        return sqrt((1 - k2 * sn**2) if k2 < 0 else
                        (k2 * cn**2 + self._kp2))

    def fE(self, phi_or_sn, cn=None, dn=None):
        '''The incomplete integral of the second kind in terms
           of Jacobi elliptic functions.

           @param phi_or_sn: φ or sinφ.
           @keyword cn: cosφ.
           @keyword dn: sqrt(1 − k2 sin2φ).

           @return: E(φ, k) as though φ ∈ (−π, π],
                    U{defined<https://DLMF.NIST.gov/19.2.E5>}.
        '''
        def _fE(sn, cn, dn):
            sn2, cn2, dn2 = sn**2, cn**2, dn**2
            kp2, k2 = self._kp2, self._k2
            if k2 <= 0:
                # Carlson, eq. 4.6 and <https://DLMF.NIST.gov/19.25.E9>
                ei = _RF(cn2, dn2, 1) - k2 * sn2 * _RD_3(cn2, dn2, 1)
            elif kp2 >= 0:
                # <https://DLMF.NIST.gov/19.25.E10>
                ei = fsum_(kp2 * _RF(cn2, dn2, 1),
                           kp2 * k2 * sn2 * _RD_3(cn2, 1, dn2),
                           k2 * abs(cn) / dn)
            else:
                # <https://DLMF.NIST.gov/19.25.E11>
                ei = dn / abs(cn) - kp2 * sn2 * _RD_3(dn2, 1, cn2)
            return ei * abs(sn)

        return self._fXf(phi_or_sn, cn, dn, self.cE,
                                            self.deltaE, _fE)

    def fEd(self, deg):
        '''The incomplete integral of the second kind with
           the argument given in degrees.

           @param deg: Angle (C{degrees}).

           @return: E(π deg/180, k).
        '''
        n = ceil(deg / 360.0 - 0.5)
        sn, cn = sincos2d(deg - n * 360.0)
        return self.fE(sn, cn, self.fDelta(sn, cn)) + 4 * self.cE * n

    def fEinv(self, x):
        '''The inverse of the incomplete integral of the second kind.

           @param x: Argument (C{float}).

           @return: φ = 1 / E(B{C{x}}, k), such that E(φ, k) = B{C{x}}.
        '''
        E2 = self.cE * 2.0
        n  = floor(x / E2 + 0.5)
        x -= E2 * n  # x now in [-ec, ec)
        # linear approximation
        phi = PI * x / E2  # phi in [-pi/2, pi/2)
        Phi = Fsum(phi)
        # first order correction
        phi = Phi.fsum_(-self._eps * sin(2 * phi) * 0.5)
        # For kp2 close to zero use asin(x/.cE) or J. P. Boyd,
        # Applied Math. and Computation 218, 7005-7013 (2012)
        # <https://DOI.org/10.1016/j.amc.2011.12.021>
        for _ in range(self._trips_):  # GEOGRAPHICLIB_PANIC
            sn, cn, dn = self._sncndn3(phi)
            phi, e = Phi.fsum2_((x - self.fE(sn, cn, dn)) / dn)
            if abs(e) < _tolJAC:
                return n * PI + phi
        raise EllipticError('no %s convergence' % ('fEinv',))

    def fF(self, phi_or_sn, cn=None, dn=None):
        '''The incomplete integral of the first kind in terms
           of Jacobi elliptic functions.

           @param phi_or_sn: φ or sinφ.
           @keyword cn: cosφ.
           @keyword dn: sqrt(1 − k2 sin2φ).

           @return: F(φ, k) as though φ ∈ (−π, π],
                    U{defined<https://DLMF.NIST.gov/19.2.E4>}.
        '''
        def _fF(sn, cn, dn):
            return abs(sn) * _RF(cn**2, dn**2, 1)

        return self._fXf(phi_or_sn, cn, dn, self.cK,
                                            self.deltaF, _fF)

    def fG(self, phi_or_sn, cn=None, dn=None):
        '''Legendre's geodesic longitude integral in terms
           of Jacobi elliptic functions.

           @param phi_or_sn: φ or sinφ.
           @keyword cn: cosφ.
           @keyword dn: sqrt(1 − k2 sin2φ).

           @return: G(φ, k) as though φ ∈ (−π, π].

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
        '''Cayley's geodesic longitude difference integral in terms
           of Jacobi elliptic functions.

           @param phi_or_sn: φ or sinφ.
           @keyword cn: cosφ.
           @keyword dn: sqrt(1 − k2 sin2φ).

           @return: H(φ, k) as though φ ∈ (−π, π].

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
        '''The incomplete integral of the third kind in terms
           of Jacobi elliptic functions.

           @param phi_or_sn: φ or sinφ.
           @keyword cn: cosφ.
           @keyword dn: sqrt(1 − k2 sin2φ).

           @return: Π(φ, α2, k) as though φ ∈ (−π, π].
        '''
        return self._fXa(phi_or_sn, cn, dn, self.alpha2,
                                            self.cPi, self.deltaPi)

    def _fXa(self, phi_or_sn, cn, dn, aX, cX, deltaX):
        '''(INTERNAL) Helper for C{.fG}, C{.fH} and C{.fPi}.
        '''
        def _fX(sn, cn, dn):
            cn2, sn2, dn2 = cn**2, sn**2, dn**2
            return abs(sn) * (_RF(cn2, dn2, 1) + aX * sn2 *
                            _RJ_3(cn2, dn2, 1, cn2 + self.alphap2 * sn2))

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
            n = self.classname + '.f' + deltaX.__name__[5:]
            raise EllipticError('%s invalid: %s%r' % ('args', n, (sn, cn, dn)))

        if cn < 0:  # enforce usual trig-like symmetries
            xi = 2 * cX - fX(sn, cn, dn)
        elif cn > 0:
            xi = fX(sn, cn, dn)
        else:
            xi = cX
        return copysign(xi, sn)

    @property_RO
    def k2(self):
        '''Get the square of the modulus (C{float}).
        '''
        return self._k2

    @property_RO
    def kp2(self):
        '''Get the square of the complementary modulus (C{float}).
        '''
        return self._kp2

    def reset(self, k2=0, alpha2=0, kp2=None, alphap2=None):  # MCCABE 13
        '''Reset the modulus and parameter.

           @keyword k2: modulus squared (C{float} k2 <= 1).
           @keyword alpha2: parameter (C{float} α2 <= 1).
           @keyword kp2: complementary modulus squared (C{float} k'2 >= 0).
           @keyword alphap2: complementary parameter α'2 (C{float} α'2 >= 0).

           @note: The arguments must satisfy I{B{C{k2}} + B{C{k2}} = 1}
                  and I{B{C{alpha2}} + B{C{alphap2}} = 1}.  No checking
                  is done that these conditions are met to enable
                  accuracy to be maintained, e.g., when C{k} is very
                  close to unity.
        '''
        if k2 > 1:
            raise EllipticError('%s invalid: %r' % ('k2', k2))
        else:
            self._k2 = k2 = float(k2)

        if alpha2 > 1:
            raise EllipticError('%s invalid: %r' % ('alpha2', alpha2))
        else:
            self._alpha2 = alpha2 = float(alpha2)

        if kp2 is None:
            self._kp2 = kp2 = 1 - k2
        elif kp2 < 0:
            raise EllipticError('%s invalid: %r' % ('kp2', kp2))
        else:
            self._kp2 = kp2 = float(kp2)

        if alphap2 is None:
            self._alphap2 = alphap2 = 1 - alpha2
        elif alphap2 < 0:
            raise EllipticError('%s invalid: %r' % ('alphap2', alphap2))
        else:
            self._alphap2 = alphap2 = float(alphap2)

        self._eps = k2 / (sqrt(kp2) + 1)**2

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
        # G( 0, k) = E(k)
        # H( 0, k) = K(k) - D(k)
        # Pi(0, k) = K(k)
        # H( 1, k) = K(k)
        # Pi(alpha2, 0) = pi/(2*sqrt(1-alpha2))
        # G( alpha2, 0) = pi/(2*sqrt(1-alpha2))
        # H( alpha2, 0) = pi/(2*(1 + sqrt(1-alpha2)))
        # Pi(alpha2, 1) = inf
        # G( alpha2, 1) = H(alpha2, 1) = RC(1, alphap2)
        if k2:
            if kp2:
                # D(k) = (K(k) - E(k))/k^2, Carlson eq.4.3
                # <https://DLMF.NIST.gov/19.25.E1>
                self._Dc = _RD_3(0, kp2, 1)
                self._KEc = k2 * self._Dc
                # Complete elliptic integral E(k), Carlson eq. 4.2
                # <https://DLMF.NIST.gov/19.25.E1>
                self._Ec = _RG(kp2, 1) * 2
                # Complete elliptic integral K(k), Carlson eq. 4.1
                # https://DLMF.NIST.gov/19.25.E1
                self._Kc = _RF(kp2, 1)
        else:
            self._Dc  = PI_4
            self._Ec  = PI_2
            self._Kc  = PI_2
            self._KEc = 0.0

        if alpha2:
            if kp2:
                # <https://DLMF.NIST.gov/19.25.E2>
                if alphap2:
                    rj_3 = _RJ_3(0, kp2, 1, alphap2)
                    # G(alpha^2, k)
                    self._Gc = self._Kc + (alpha2 - k2) * rj_3
                    # H(alpha^2, k)
                    self._Hc = self._Kc - alphap2 * rj_3
                    # Pi(alpha^2, k)
                    self._Pic = self._Kc + alpha2 * rj_3
                # else:
                #   self._Gc = self._Hc = self._Pc = _INF  # XXX _NAN?
            elif alphap2:
                self._Gc = self._Hc = _RC(1, alphap2)
        else:
            self._Gc  = self._Ec
            self._Pic = self._Kc
            # Hc = Kc - Dc but this involves large cancellations
            # if k2 is close to 1.  So write (for alpha2 = 0)
            #   Hc = int(cos(phi)^2/sqrt(1-k2*sin(phi)^2),phi,0,pi/2)
            #      = 1/sqrt(1-k2) * int(sin(phi)^2/sqrt(1-k2/kp2*sin(phi)^2,...)
            #      = 1/kp * D(i*k/kp)
            # and use D(k) = RD(0, kp2, 1) / 3
            # so Hc = 1/kp * RD(0, 1/kp2, 1) / 3
            #       = kp2 * RD(0, 1, kp2) / 3
            # using https://DLMF.NIST.gov/19.20.E18
            # Equivalently
            #   RF(x, 1) - RD(0, x, 1)/3 = x * RD(0, 1, x)/3 for x > 0
            # For k2 = 1 and alpha2 = 0, we have
            #   Hc = int(cos(phi),...) = 1
            self._Hc = kp2 * _RD_3(0, 1, kp2) if kp2 else 1.0

    def sncndn(self, x):
        '''The Jacobi elliptic function.

           @param x: The argument (C{float}).

           @return: 3-Tuple C{(sn(x, k), cn(x, k), dn(x, k))}.
        '''
        # Bulirsch's sncndn routine, p 89.
        mc = self._kp2
        if mc:
            if mc < 0:
                d = 1.0 - mc
                mc /= -d
                d = sqrt(d)
                x *= d
            else:
                d = 0
            a, c, mn = 1.0, 0, []
            for _ in range(self._trips_):  # GEOGRAPHICLIB_PANIC
                # This converges quadratically, max 6 trips
                mc = sqrt(mc)
                mn.append((a, mc))
                c = (a + mc) * 0.5
                if not abs(a - mc) > (_tolJAC * a):
                    break
                mc *= a
                a = c
            else:
                raise EllipticError('no %s convergence' % ('sncndn',))
            x *= c
            dn = 1
            sn, cn = sincos2(x)
            if sn:
                a = cn / sn
                c *= a
                while mn:
                    m, n = mn.pop()
                    a *= c
                    c *= dn
                    dn = (n + a) / (m + a)
                    a = c / m
                sn = copysign(1.0 / hypot1(c), sn)
                cn = c * sn
                if d:
                    cn, dn = dn, cn
                    sn /= d
        else:
            sn = tanh(x)
            cn = dn = 1.0 / cosh(x)
        return sn, cn, dn

    def _sncndn3(self, phi):
        '''(INTERNAL) Helper for C{.fEinv} and C{._fXf}.
        '''
        sn, cn = sincos2(phi)
        return sn, cn, self.fDelta(sn, cn)


def _Hf(e0, e1, e2, e3, e4, e5):
    '''(INTERNAL) Horner form for C{_RD} and C{_RJ} below.
    '''
    # Polynomial is <https://DLMF.NIST.gov/19.36.E2>
    # (1 - 3*E2/14 + E3/6 + 9*E2^2/88 - 3*E4/22 - 9*E2*E3/52 + 3*E5/26
    #    - E2^3/16 + 3*E3^2/40 + 3*E2*E4/20 + 45*E2^2*E3/272
    #    - 9*(E3*E4+E2*E5)/68)
    return fsum_(fsum_(471240,      -540540 * e2) * e5,
                 fsum_(612612 * e2, -540540 * e3,    -556920) * e4,
                 fsum_(306306 * e3,  675675 * e2**2, -706860 * e2, 680680) * e3,
                 fsum_(417690 * e2, -255255 * e2**2, -875160) * e2,
                 4084080) / (4084080 * e1) + e0


def _Q(A, T, tol):
    '''(INTERNAL) Helper for C{_RD}, C{_RF} and C{_RJ}.
    '''
    return max(abs(A - t) for t in T[1:]) / tol


def _RC(x, y):
    '''Degenerate symmetric integral of the first kind C{_RC}.

       @return: C{_RC(x, y) = _RF(x, y, y)}.

       @see: U{C{_RC} definition<https://DLMF.NIST.gov/19.2.E17>}.
    '''
    # Defined only for y != 0 and x >= 0.
    d = x - y
    if d < 0:  # catch _NaN
        # <https://DLMF.NIST.gov/19.2.E18>
        d = -d
        r = atan(sqrt(d / x))
    elif d == 0:
        r, d = 1.0, y
    elif y > 0:  # <https://DLMF.NIST.gov/19.2.E19>
        # atanh(sqrt((x - y) / x))
        r = asinh(sqrt(d / y))
    else:  # <https://DLMF.NIST.gov/19.2.E20>
        # atanh(sqrt(x / (x - y)))
        r = asinh(sqrt(-x / y))
    return r / sqrt(d)


def _RD_3(x, y, z):
    '''Degenerate symmetric integral of the third kind C{_RD} / 3.
    '''
    return _RD(x, y, z) / 3.0


def _RD(x, y, z):
    '''Degenerate symmetric integral of the third kind C{_RD}.

       @return: C{_RD(x, y, z) = _RJ(x, y, z, z)}.

       @see: U{C{_RD} definition<https://DLMF.NIST.gov/19.16.E5>}.
    '''
    # Carlson, eqs 2.28 - 2.34
    m = 1.0
    S = Fsum()
    A = fsum_(x, y, z, z, z) * 0.2
    T = [A, x, y, z]
    Q = _Q(A, T, _tolRD)
    for _ in range(_TRIPS):
        if Q < abs(m * T[0]):  # max 7 trips
            break
        t = T[3]  # z0
        r, s, T = _rsT(T)
        S.fadd_(1.0 / (m * s[2] * (t + r)))
        m *= 4
    else:
        raise EllipticError('no %s convergence' % ('RD',))

    S *= 3
    m *= T[0]  # An
    x = (x - A) / m
    y = (y - A) / m
    z = (x + y) / 3.0
    z2 = z**2
    xy = x * y
    return _Hf(S.fsum(), m * sqrt(T[0]),
               xy - 6 * z2,
              (xy * 3 - 8 * z2) * z,
              (xy - z2) * 3 * z2,
               xy * z2 * z)


def _RF(x, y, z=None):
    '''Symmetric integral of the first kind C{_RF}.

       @return: C{_RF(x, y, z)}.

       @see: U{C{_RF} definition<https://DLMF.NIST.gov/19.16.E1>}.
    '''
    if z is None:
        # Carlson, eqs 2.36 - 2.38
        a, b = sqrt(x), sqrt(y)
        if a < b:
            a, b = b, a
        while abs(a - b) > (_tolRG0 * a):  # max 4 trips
            b, a = sqrt(a * b), (a + b) * 0.5
        return PI / (a + b)

    # Carlson, eqs 2.2 - 2.7
    m = 1.0
    A = fsum_(x, y, z) / 3.0
    T = [A, x, y, z]
    Q = _Q(A, T, _tolRF)
    for _ in range(_TRIPS):
        if Q < abs(m * T[0]):  # max 6 trips
            break
        _, _, T = _rsT(T)
        m *= 4
    else:
        raise EllipticError('no %s convergence' % ('RF',))

    m *= T[0]  # An
    x = (A - x) / m
    y = (A - y) / m
    z = -(x + y)

    e2 = x * y - z**2
    e3 = x * y * z
    # Polynomial is <https://DLMF.NIST.gov/19.36.E1>
    # (1 - E2/10 + E3/14 + E2^2/24 - 3*E2*E3/44
    #    - 5*E2^3/208 + 3*E3^2/104 + E2^2*E3/16)
    # convert to Horner form ...
    return fsum_(fsum_( 6930 * e3, 15015 * e2**2, -16380 * e2, 17160) * e3,
                 fsum_(10010 * e2, -5775 * e2**2, -24024) * e2,
                 240240) / (240240 * sqrt(T[0]))


def _RG(x, y, z=None):
    '''Symmetric integral of the second kind C{_RG}.

       @return: C{_RG(x, y, z)}.

       @see: U{C{_RG} definition<https://DLMF.NIST.gov/19.16.E3>}
             and in Carlson, eq. 1.5.
    '''
    if z is None:
        # Carlson, eqs 2.36 - 2.39
        a, b = sqrt(x), sqrt(y)
        if a < b:
            a, b = b, a
        S = Fsum(0.25 * (a + b)**2)
        m = -0.25  # note, negative
        while abs(a - b) > (_tolRG0 * a):  # max 4 trips
            b, a = sqrt(a * b), (a + b) * 0.5
            m *= 2
            S.fadd_(m * (a - b)**2)
        return S.fsum() * PI_2 / (a + b)

    if not z:
        y, z = z, y
    # Carlson, eq 1.7
    return fsum_(_RF(x, y, z) * z,
                 _RD_3(x, y, z) * (x - z) * (z - y),
                 sqrt(x * y / z)) * 0.5


def _RJ_3(x, y, z, p):
    '''Symmetric integral of the third kind C{_RJ} / 3.
    '''
    return _RJ(x, y, z, p) / 3.0


def _RJ(x, y, z, p):
    '''Symmetric integral of the third kind C{_RJ}.

       @return: C{_RJ(x, y, z, p)}.

       @see: U{C{_RJ} definition<https://DLMF.NIST.gov/19.16.E2>}.
    '''
    def _xyzp(x, y, z, p):
        return (x + p) * (y + p) * (z + p)

    # Carlson, eqs 2.17 - 2.25
    m = m3 = 1.0
    S = Fsum()
    D = -_xyzp(x, y, z, -p)
    A = fsum_(x, y, z, 2 * p) * 0.2
    T = [A, x, y, z, p]
    Q = _Q(A, T, _tolRD)
    for _ in range(_TRIPS):
        if Q < abs(m * T[0]):  # max 7 trips
            break
        _, s, T = _rsT(T)
        d = _xyzp(*s)
        e = D / (m3 * d**2)
        S.fadd_(_RC(1, 1 + e) / (m * d))
        m *= 4
        m3 *= 64
    else:
        raise EllipticError('no %s convergence' % ('RJ',))

    S *= 6
    m *= T[0]  # An
    x = (A - x) / m
    y = (A - y) / m
    z = (A - z) / m
    xyz = x * y * z
    p = -(x + y + z) * 0.5
    p2 = p**2

    e2 = fsum_(x * y, x * z, y * z, -3 * p2)
    return _Hf(S.fsum(), m * sqrt(T[0]),
               e2,
               fsum_(xyz, 2 * p * e2, 4 * p * p2),
               fsum_(xyz * 2, p * e2, 3 * p * p2) * p,
               p2 * xyz)


def _rsT(T):
    '''(INTERNAL) Helper for C{_RD}, C{_RF} and C{_RJ}.
    '''
    s = map2(sqrt, T[1:])
    r = fdot(s[:3], s[1], s[2], s[0])
    T[:] = [(t + r) * 0.25 for t in T]
    return r, s, T

# **) MIT License
#
# Copyright (C) 2016-2020 -- mrJean1 at Gmail -- All Rights Reserved.
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
