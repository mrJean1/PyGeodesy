
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
from pygeodesy.constants import EPS, _EPStol as _TolJAC, INF, PI, PI_2, PI_4, \
                               _0_0, _0_125, _0_25, _0_5, _1_0, _1_64th, _2_0, \
                               _N_2_0, _3_0, _4_0, _6_0, _8_0, _180_0, _360_0
from pygeodesy.errors import _ValueError, _xkwds_pop
from pygeodesy.fmath import fdot, Fsum, hypot1
# from pygeodesy.fsums import Fsum  # from .fmath
from pygeodesy.interns import NN, _DOT_, _f_, _SPACE_
from pygeodesy.karney import _ALL_LAZY, _signBit
# from pygeodesy.lazily import _ALL_LAZY  # from .karney
from pygeodesy.named import _Named, _NamedTuple, _NotImplemented
from pygeodesy.props import _allPropertiesOf_n, Property_RO, property_RO, \
                            _update_all
from pygeodesy.streprs import Fmt, unstr
from pygeodesy.units import Scalar, Scalar_
from pygeodesy.utily import sincos2, sincos2d

from math import asinh, atan, atan2, ceil, cosh, floor, sin, sqrt, tanh

__all__ = _ALL_LAZY.elliptic
__version__ = '22.09.12'

_delta_      = 'delta'
_invokation_ = 'invokation'
_TolRD       =  pow(EPS * 0.002, _0_125)  # 8th root: quadquadratic?, octic?, ocrt?
_TolRF       =  pow(EPS * 0.030, _0_125)  # 4th root: biquadratic, quartic, qurt?
_TolRG0      = _TolJAC  * 2.7
_TRIPS       =  31  # Max depth, 7 might be sufficient


class _Complete(object):
    '''(INTERAL) Hold complete integrals.
    '''
    def __init__(self, **kwds):
        self.__dict__ = kwds


class _Deferred(list):
    '''(INTERNAL) Collector for L{_Deferred_Fsum}.
    '''
    def __init__(self, *xs):
        if xs:
            list.__init__(self, xs)

    def __add__(self, other):  # PYCHOK no cover
        return _NotImplemented(self, other)

    def __iadd__(self, x):  # overide C{list} += C{list}
        # assert isscalar(x)
        self.append(float(x))
        return self

    def __imul__(self, other):  # PYCHOK no cover
        return _NotImplemented(self, other)

    def __isub__(self, x):  # PYCHOK no cover
        # assert isscalar(x)
        self.append(-float(x))
        return self

    def __mul__(self, other):  # PYCHOK no cover
        return _NotImplemented(self, other)

    def __radd__(self, other):  # PYCHOK no cover
        return _NotImplemented(self, other)

    def __rsub__(self, other):  # PYCHOK no cover
        return _NotImplemented(self, other)

    def __sub__(self, other):  # PYCHOK no cover
        return _NotImplemented(self, other)

    @property_RO
    def Fsum(self):
        # get a L{_Deferred_Fsum} instance, pre-named
        return _Deferred_Fsum()._facc(self)  # known C{float}s


class _Deferred_Fsum(Fsum):
    '''(INTERNAL) Deferred L{Fsum}.
    '''
    name = NN  # pre-named, overridden below

    def _update(self, **other):  # PYCHOK don't ...
        # ... waste time zapping non-existing Property/_ROs
        if other or len(self.__dict__) > 2:
            Fsum._update(self, **other)

_Deferred_Fsum.name = _Deferred_Fsum.__name__  # PYCHOK once


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
        return self._reset_cDcEcKcKE_eps.cD

    @Property_RO
    def cE(self):
        '''Get the complete integral of the second kind C{E(k)}
           (C{float}), U{defined<https://DLMF.NIST.gov/19.2.E5>}.
        '''
        return self._reset_cDcEcKcKE_eps.cE

    @Property_RO
    def cG(self):
        '''Get Legendre's complete geodesic longitude integral
           C{G(α^2, k)} (C{float}).
        '''
        return self._reset_cGcHcPi.cG

    @Property_RO
    def cH(self):
        '''Get Cayley's complete geodesic longitude difference integral
           C{H(α^2, k)} (C{float}).
        '''
        return self._reset_cGcHcPi.cH

    @Property_RO
    def cK(self):
        '''Get the complete integral of the first kind C{K(k)}
           (C{float}), U{defined<https://DLMF.NIST.gov/19.2.E4>}.
        '''
        return self._reset_cDcEcKcKE_eps.cK

    @Property_RO
    def cKE(self):
        '''Get the difference between the complete integrals of the
           first and second kinds, C{K(k) − E(k)} (C{float}).
        '''
        return self._reset_cDcEcKcKE_eps.cKE

    @Property_RO
    def cPi(self):
        '''Get the complete integral of the third kind C{Pi(α^2, k)}
           (C{float}), U{defined<https://DLMF.NIST.gov/19.2.E7>}.
        '''
        return self._reset_cGcHcPi.cPi

    def deltaD(self, sn, cn, dn):
        '''The periodic Jahnke's incomplete elliptic integral.

           @arg sn: sin(φ).
           @arg cn: cos(φ).
           @arg dn: sqrt(1 − k2 * sin(2φ)).

           @return: Periodic function π D(φ, k) / (2 D(k)) - φ (C{float}).

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        return self._deltaX(sn, cn, dn, self.cD, self.fD)

    def deltaE(self, sn, cn, dn):
        '''The periodic incomplete integral of the second kind.

           @arg sn: sin(φ).
           @arg cn: cos(φ).
           @arg dn: sqrt(1 − k2 * sin(2φ)).

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
           @arg dn: sqrt(1 − k2 * sin(2φ)).

           @return: Periodic function π F(φ, k) / (2 K(k)) - φ (C{float}).

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        return self._deltaX(sn, cn, dn, self.cK, self.fF)

    def deltaG(self, sn, cn, dn):
        '''Legendre's periodic geodesic longitude integral.

           @arg sn: sin(φ).
           @arg cn: cos(φ).
           @arg dn: sqrt(1 − k2 * sin(2φ)).

           @return: Periodic function π G(φ, k) / (2 G(k)) - φ (C{float}).

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        return self._deltaX(sn, cn, dn, self.cG, self.fG)

    def deltaH(self, sn, cn, dn):
        '''Cayley's periodic geodesic longitude difference integral.

           @arg sn: sin(φ).
           @arg cn: cos(φ).
           @arg dn: sqrt(1 − k2 * sin(2φ)).

           @return: Periodic function π H(φ, k) / (2 H(k)) - φ (C{float}).

           @raise EllipticError: Invalid invokation or no convergence.
        '''
        return self._deltaX(sn, cn, dn, self.cH, self.fH)

    def deltaPi(self, sn, cn, dn):
        '''The periodic incomplete integral of the third kind.

           @arg sn: sin(φ).
           @arg cn: cos(φ).
           @arg dn: sqrt(1 − k2 * sin(2φ)).

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
            raise _invokationError(n, sn, cn, dn)

        if _signBit(cn):
            cn, sn = -cn, -sn
        return fX(sn, cn, dn) * PI_2 / cX - atan2(sn, cn)

    @Property_RO
    def eps(self):
        '''Get epsilon (C{float}).
        '''
        return self._reset_cDcEcKcKE_eps.eps

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

           @return: sqrt(1 − k2 * sin(2φ)) (C{float}).
        '''
        k2 = self.k2
        s  = (_1_0 - k2 * sn**2) if k2 < 0 else (self.kp2
                 + ((k2 * cn**2) if k2 > 0 else _0_0))
        return sqrt(s) if s else _0_0

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
                    ei = _RF3(self, cn2, dn2, _1_0)
                    if k2:
                        ei -= _RD(self, cn2, dn2, _1_0, _3over(k2, sn2))
                elif kp2 >= 0:  # <https://DLMF.NIST.gov/19.25.E10>
                    ei = k2 * abs(cn) / dn
                    if kp2:
                        ei += (_RD( self, cn2, _1_0,  dn2, _3over(k2, sn2)) +
                               _RF3(self, cn2,  dn2, _1_0)) * kp2
                else:  # <https://DLMF.NIST.gov/19.25.E11>
                    ei  =  dn / abs(cn)
                    ei -= _RD(self, dn2, _1_0, cn2, _3over(kp2, sn2))
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
        fE         = self.fE
        _sncndnPhi = self._sncndnPhi
        _Phi2      = Phi.fsum2_
        self._iteration = 0  # aggregate
        for i in range(1, _TRIPS):  # GEOGRAPHICLIB_PANIC
            sn, cn, dn = _sncndnPhi(phi)
            phi, e = _Phi2((y - fE(sn, cn, dn)) / dn)
            if abs(e) < _TolJAC:
                self._iteration += i
                break
        else:  # PYCHOK no cover
            raise _no_convergenceError(e, _TolJAC, self.fEinv, x)
        return Phi.fsum_(n * PI) if n else phi

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
                    r += _RJ(self, cn2, dn2, _1_0, z, _3over(aX, sn2))
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
            raise _invokationError(n, sn, cn, dn)

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

    @Property_RO
    def _reset_cDcEcKcKE_eps(self):
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

    @Property_RO
    def _reset_cGcHcPi(self):
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
            c, d, cd, mn_ = self._sncndnBulirsch4
            dn = _1_0
            sn, cn = sincos2(x * cd)
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
                    sn = sn / d  # /= chokes PyChecker
        else:
            sn = tanh(x)
            cn = dn = _1_0 / cosh(x)

        return Elliptic3Tuple(sn, cn, dn, iteration=self._iteration)

    @Property_RO
    def _sncndnBulirsch4(self):
        '''(INTERNAL) Get Bulirsch' 4-tuple C{(c, d, cd, mn_)}.
        '''
        # Bulirsch's sncndn routine, p 89.
        d, mc = 0, self.kp2
        if _signBit(mc):  # PYCHOK no cover
            d  = _1_0 - mc
            mc =  neg(mc / d)
            d  =  sqrt(d)

        mn, a = [], _1_0
        for i in range(1, _TRIPS):  # GEOGRAPHICLIB_PANIC
            # This converges quadratically, max 6 trips
            mc = sqrt(mc)
            mn.append((a, mc))
            c = (a + mc) * _0_5
            t = _TolJAC  *  a
            if abs(a - mc) <= t:
                self._iteration += i  # accumulate
                break
            mc *= a
            a   = c
        else:  # PYCHOK no cover
            raise _no_convergenceError(a - mc, t, None, kp=self.kp, kp2=self.kp2)
        cd = (c * d) if d else c
        return c, d, cd, tuple(reversed(mn))  # mn reversed!

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
        return _RG3(None, x, y, *z) if z and z[0] else (
               _RG2(None, x, y) * _0_5)

    @staticmethod
    def fRJ(x, y, z, p):
        '''Symmetric integral of the third kind C{RJ(x, y, z, p)}.

           @return: C{RJ(x, y, z, p)}.

           @see: U{C{RJ} definition<https://DLMF.NIST.gov/19.16.E2>} and
                 U{Carlson<https://ArXiv.org/pdf/math/9409227.pdf>}.
        '''
        return _RJ(None, x, y, z, p)

_allPropertiesOf_n(15, Elliptic)  #  # PYCHOK assert, see Elliptic.reset


class EllipticError(_ValueError):
    '''Elliptic integral, function, convergence or other L{Elliptic} issue.
    '''
    pass


class Elliptic3Tuple(_NamedTuple):
    '''3-Tuple C{(sn, cn, dn)} all C{scalar}.
    '''
    _Names_ = ('sn',   'cn',   'dn')
    _Units_ = ( Scalar, Scalar, Scalar)


class _Lxyz(list):
    '''(INTERNAL) Helper for C{_RD}, C{_RF3} and C{_RJ}.
    '''
    _a  = None
    _a0 = None

    def __init__(self, *xyz_):  # x, y, z [, p]
        list.__init__(self, xyz_)

    def a0(self, n):
        '''Compute the initial C{a}.
        '''
        t = tuple(self)
        m = n - len(t)
        if m > 0:
            t += t[-1:] * m
        self._a0 = self._a = Fsum(*t).fover(n)
        return self._a0

    def asr3(self, a):
        '''Compute next C{a}, C{sqrt(xyz_)} and C{fdot(sqrt(xyz))}.
        '''
        L    = self
        # assert a is L._a
        s    = map2(sqrt, L)  # sqrt(x), srqt(y), sqrt(z) [, sqrt(p)]
        r    = fdot(s[:3], s[1], s[2], s[0])  # sqrt(x) * sqrt(y) + ...
        L[:] = [(x + r) * _0_25 for x in L]
        # assert L is self
        L._a = a = (a + r) * _0_25
        return a, s, r

    def rescale(self, am, *xy_):
        '''Rescale C{x}, C{y}, ...
        '''
        for x in xy_:
            yield (self._a0 - x) / am

    def thresh(self, Tol):
        '''Return the convergence threshold.
        '''
        return max(abs(self._a0 - x) for x in self) / Tol


def _horner(S, e1, E2, E3, E4, E5, *over):
    '''(INTERNAL) Horner form for C{_RD} and C{_RJ} below.
    '''
    E22 = E2**2
    # Polynomial is <https://DLMF.NIST.gov/19.36.E2>
    # (1 - 3*E2/14 + E3/6 + 9*E2**2/88 - 3*E4/22 - 9*E2*E3/52
    #    + 3*E5/26 - E2**3/16 + 3*E3**2/40 + 3*E2*E4/20
    #    + 45*E2**2*E3/272 - 9*(E3*E4+E2*E5)/68)
    # converted to Horner-like form ...
    e  = e1 * 4084080
    S *= e
    S += Fsum(                             E2 * -540540,  471240).fmul(E5)
    S += Fsum(               E3 * -540540, E2 *  612612, -556920).fmul(E4)
    S += Fsum(E22 *  675675, E3 *  306306, E2 * -706860,  680680).fmul(E3)
    S += Fsum(E22 * -255255,               E2 *  417690, -875160).fmul(E2)
    S += 4084080
    return S.fover((e * over[0]) if over else e)


def _invokationError(name, *args):  # PYCHOK no cover
    '''(INTERNAL) Return an L{EllipticError}.
    '''
    n = _DOT_(Elliptic.__name__, name)
    n = _SPACE_(_invokation_, n)
    return EllipticError(NN(n, repr(args)))  # unstr


def _iterations(inst, i):
    '''(INTERNAL) Aggregate iterations B{C{i}}.
    '''
    if inst:
        inst._iteration += i


def _no_convergenceError(d, tol, where, *args, **kwds_thresh):  # PYCHOK no cover
    '''(INTERNAL) Return an L{EllipticError}.
    '''
    n = Elliptic.__name__
    if where:
        n = _DOT_(n, where.__name__)
    if kwds_thresh:
        q = _xkwds_pop(kwds_thresh, thresh=False)
        t =  unstr(n, *args, **kwds_thresh)
    else:
        q = False
        t = unstr(n, *args)
    return EllipticError(Fmt.no_convergence(d, tol, thresh=q), txt=t)


def _3over(a, b):
    '''(INTERNAL) Return C{3 / (a * b)}.
    '''
    return _3_0 / (a * b)


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
        raise _invokationError(Elliptic.fRC.__name__, x, y)
    return r / sqrt(d)


def _RD(inst, x, y, z, *over):
    '''(INTERNAL) Carlson, eqs 2.28 - 2.34.
    '''
    L = _Lxyz(x, y, z)
    a =  L.a0(5)
    q =  L.thresh(_TolRF)
    S = _Deferred()
    am, m = a, 1
    for i in range(_TRIPS):
        if abs(am) > q:  # max 7 trips
            _iterations(inst, i)
            break
        t = L[2]  # z0...n
        a, s, r = L.asr3(a)
        S += _1_0 / ((t + r) * s[2] * m)
        m *=  4
        am =  a * m
    else:  # PYCHOK no cover
        raise _no_convergenceError(am, q, Elliptic.fRD, x, y, z, *over,
                                          thresh=True)
    x, y = L.rescale(-am, x, y)
    xy =  x * y
    z  = (x + y) / _3_0
    z2 =  z**2
    S  =  S.Fsum.fmul(_3_0)
    return _horner(S, am * sqrt(a),
                   xy        - _6_0 * z2,
                  (xy * _3_0 - _8_0 * z2) * z,
                  (xy - z2)  * _3_0 * z2,
                   xy * z2 * z, *over)


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
        a, b = ((a + b) * _0_5), sqrt(a * b)
    else:  # PYCHOK no cover
        raise _no_convergenceError(a - b, t, Elliptic.fRF, x, y)


def _RF3(inst, x, y, z):  # 3-arg version
    '''(INTERNAL) Carlson, eqs 2.2 - 2.7.
    '''
    L = _Lxyz(x, y, z)
    a =  L.a0(3)
    q =  L.thresh(_TolRF)
    am, m = a, 1
    for i in range(_TRIPS):
        if abs(am) > q:  # max 6 trips
            _iterations(inst, i)
            break
        a, _, _ = L.asr3(a)
        m *= 4
        am = a * m
    else:  # PYCHOK no cover
        raise _no_convergenceError(am, q, Elliptic.fRF, x, y, z,
                                          thresh=True)
    x, y = L.rescale(am, x, y)
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
    S += Fsum(240240)
    return S.fover(sqrt(a) * 240240)


def _RG2(inst, x, y):  # 2-args and I{doubled}
    '''(INTERNAL) Carlson, eqs 2.36 - 2.39.
    '''
    a, b = sqrt(x), sqrt(y)
    if a < b:
        a, b = b, a
    ab =  a - b  # abs(a - b)
    S  = _Deferred(_0_5 * (a + b)**2)
    m  = -1
    for i in range(_TRIPS):  # max 4 trips
        t = _TolRG0 * a
        if ab <= t:
            _iterations(inst, i)
            return S.Fsum.fover((a + b) / PI_2)
        a, b = ((a + b) * _0_5), sqrt(a * b)
        ab = abs(a - b)
        S += ab**2 * m
        m *= 2
    else:  # PYCHOK no cover
        raise _no_convergenceError(ab, t, Elliptic.fRG, x, y)


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

    L = _Lxyz(x, y, z, p)
    a =  L.a0(5)
    q =  L.thresh(_TolRD)
    S = _Deferred()
    n =  neg(_xyzp(x, y, z, -p))
    am, m = a, 1
    for i in range(_TRIPS):
        if abs(am) > q:  # max 7 trips
            _iterations(inst, i)
            break
        a, s, _ = L.asr3(a)
        d = _xyzp(*s)
        if n:
            rc = _RC(inst, _1_0, n / d**2 + _1_0)
            n *= _1_64th  # /= chokes PyChecker
        else:
            rc = _1_0  # == _RC(None, _1_0, _1_0)
        S += rc / (d * m)
        m *= 4
        am = a * m
    else:  # PYCHOK no cover
        raise _no_convergenceError(am, q, Elliptic.fRJ, x, y, z, p,
                                          thresh=True)
    x, y, z = L.rescale(am, x, y, z)
    xyz = x * y * z
    p   = Fsum(x, y, z).fover(_N_2_0)
    p2  = p**2
    p3  = p2*p
    E2  = Fsum(x * y, x * z, y * z, -p2 * _3_0)
    E2p = E2 * p
    S   = S.Fsum.fmul(_6_0)
    return _horner(S, am * sqrt(a), E2,
                   Fsum(p3 * _4_0, xyz, E2p * _2_0),
                   Fsum(p3 * _3_0, E2p, xyz * _2_0).fmul(p),
                   xyz * p2, *over)

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
