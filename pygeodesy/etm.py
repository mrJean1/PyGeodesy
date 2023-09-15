
# -*- coding: utf-8 -*-

u'''A pure Python version of I{Karney}'s C{Exact Transverse Mercator} (ETM) projection.

Classes L{Etm}, L{ETMError} and L{ExactTransverseMercator}, transcoded from I{Karney}'s
C++ class U{TransverseMercatorExact<https://GeographicLib.SourceForge.io/C++/doc/
classGeographicLib_1_1TransverseMercatorExact.html>}, abbreviated as C{TMExact} below.

Class L{ExactTransverseMercator} provides C{Exact Transverse Mercator} projections while
instances of class L{Etm} represent ETM C{(easting, northing)} locations.  See also
I{Karney}'s utility U{TransverseMercatorProj<https://GeographicLib.SourceForge.io/C++/doc/
TransverseMercatorProj.1.html>} and use C{"python[3] -m pygeodesy.etm ..."} to compare
the results.

Following is a copy of I{Karney}'s U{TransverseMercatorExact.hpp
<https://GeographicLib.SourceForge.io/C++/doc/TransverseMercatorExact_8hpp_source.html>}
file C{Header}.

Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2008-2023) and licensed
under the MIT/X11 License.  For more information, see the U{GeographicLib<https://
GeographicLib.SourceForge.io>} documentation.

The method entails using the U{Thompson Transverse Mercator<https://WikiPedia.org/
wiki/Transverse_Mercator_projection>} as an intermediate projection.  The projections
from the intermediate coordinates to C{phi, lam} and C{x, y} are given by elliptic
functions.  The inverse of these projections are found by Newton's method with a
suitable starting guess.

The relevant section of L.P. Lee's paper U{Conformal Projections Based On Jacobian
Elliptic Functions<https://DOI.org/10.3138/X687-1574-4325-WM62>} in part V, pp
67-101.  The C++ implementation and notation closely follow Lee, with the following
exceptions::

  Lee   here   Description

  x/a   xi     Northing (unit Earth)

  y/a   eta    Easting (unit Earth)

  s/a   sigma  xi + i * eta

  y     x      Easting

  x     y      Northing

  k     e      Eccentricity

  k^2   mu     Elliptic function parameter

  k'^2  mv     Elliptic function complementary parameter

  m     k      Scale

  zeta  zeta   Complex longitude = Mercator = chi in paper

  s     sigma  Complex GK = zeta in paper

Minor alterations have been made in some of Lee's expressions in an attempt to
control round-off.  For example, C{atanh(sin(phi))} is replaced by C{asinh(tan(phi))}
which maintains accuracy near C{phi = pi/2}.  Such changes are noted in the code.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import map1, neg, neg_, _xinstanceof
from pygeodesy.constants import EPS, EPS02, PI_2, PI_4, _K0_UTM, \
                               _1_EPS, _0_0, _0_1, _0_5, _1_0, _2_0, \
                               _3_0, _4_0, _90_0, isnear0, isnear90
from pygeodesy.datums import _ellipsoidal_datum, _WGS84
from pygeodesy.elliptic import _ALL_LAZY, Elliptic
# from pygeodesy.errors import _incompatible  # from .named
from pygeodesy.fmath import cbrt, hypot, hypot1, hypot2
from pygeodesy.fsums import Fsum, fsum1f_
from pygeodesy.interns import NN, _COMMASPACE_, _DASH_, _near_, _SPACE_, \
                             _spherical_, _usage
from pygeodesy.karney import _copyBit, _diff182, _fix90, _norm2, _norm180, \
                             _tand, _unsigned2
# from pygeodesy.lazily import _ALL_LAZY  # from .elliptic
from pygeodesy.named import callername, _incompatible, _NamedBase
from pygeodesy.namedTuples import Forward4Tuple, Reverse4Tuple
from pygeodesy.props import deprecated_method, deprecated_property_RO, \
                            Property_RO, property_RO, _update_all, \
                            property_doc_
from pygeodesy.streprs import Fmt, pairs, unstr
from pygeodesy.units import Degrees, Scalar_
from pygeodesy.utily import atand, atan2d, _loneg, sincos2
from pygeodesy.utm import _cmlon, _LLEB, _parseUTM5, _toBand, _toXtm8, \
                          _to7zBlldfn, Utm, UTMError

from math import asinh, atan2, degrees, radians, sinh, sqrt

__all__ = _ALL_LAZY.etm
__version__ = '23.09.07'

_OVERFLOW = _1_EPS**2  # about 2e+31
_TAYTOL   =  pow(EPS, 0.6)
_TAYTOL2  = _TAYTOL * _2_0
_TOL_10   =  EPS * _0_1
_TRIPS    =  21  # C++ 10


def _overflow(x):
    '''(INTERNAL) Like C{copysign0(OVERFLOW, B{x})}.
    '''
    return _copyBit(_OVERFLOW, x)


class ETMError(UTMError):
    '''Exact Transverse Mercator (ETM) parse, projection or other
       L{Etm} issue or L{ExactTransverseMercator} conversion failure.
    '''
    pass


class Etm(Utm):
    '''Exact Transverse Mercator (ETM) coordinate, a sub-class of L{Utm},
       a Universal Transverse Mercator (UTM) coordinate using the
       L{ExactTransverseMercator} projection for highest accuracy.

       @note: Conversion of (geodetic) lat- and longitudes to/from L{Etm}
              coordinates is 3-4 times slower than to/from L{Utm}.

       @see: Karney's U{Detailed Description<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1TransverseMercatorExact.html#details>}.
    '''
    _Error   = ETMError  # see utm.UTMError
    _exactTM = None

    __init__ = Utm.__init__
    '''New L{Etm} Exact Transverse Mercator coordinate, raising L{ETMError}s.

       @see: L{Utm.__init__} for more information.

       @example:

        >>> import pygeodesy
        >>> u = pygeodesy.Etm(31, 'N', 448251, 5411932)
    '''

    @property_doc_(''' the ETM projection (L{ExactTransverseMercator}).''')
    def exactTM(self):
        '''Get the ETM projection (L{ExactTransverseMercator}).
        '''
        if self._exactTM is None:
            self.exactTM = self.datum.exactTM  # ExactTransverseMercator(datum=self.datum)
        return self._exactTM

    @exactTM.setter  # PYCHOK setter!
    def exactTM(self, exactTM):
        '''Set the ETM projection (L{ExactTransverseMercator}).

           @raise ETMError: The B{C{exacTM}}'s datum incompatible
                            with this ETM coordinate's C{datum}.
        '''
        _xinstanceof(ExactTransverseMercator, exactTM=exactTM)

        E = self.datum.ellipsoid
        if E != exactTM.ellipsoid:  # may be None
            raise ETMError(repr(exactTM), txt=_incompatible(repr(E)))
        self._exactTM = exactTM
        self._scale0  = exactTM.k0

    def parse(self, strETM, name=NN):
        '''Parse a string to a similar L{Etm} instance.

           @arg strETM: The ETM coordinate (C{str}),
                        see function L{parseETM5}.
           @kwarg name: Optional instance name (C{str}),
                        overriding this name.

           @return: The instance (L{Etm}).

           @raise ETMError: Invalid B{C{strETM}}.

           @see: Function L{pygeodesy.parseUPS5}, L{pygeodesy.parseUTM5}
                 and L{pygeodesy.parseUTMUPS5}.
        '''
        return parseETM5(strETM, datum=self.datum, Etm=self.classof,
                                 name=name or self.name)

    @deprecated_method
    def parseETM(self, strETM):  # PYCHOK no cover
        '''DEPRECATED, use method L{Etm.parse}.
        '''
        return self.parse(strETM)

    def toLatLon(self, LatLon=None, unfalse=True, **unused):  # PYCHOK expected
        '''Convert this ETM coordinate to an (ellipsoidal) geodetic point.

           @kwarg LatLon: Optional, ellipsoidal class to return the geodetic
                          point (C{LatLon}) or C{None}.
           @kwarg unfalse: Unfalse B{C{easting}} and B{C{northing}} if
                           C{falsed} (C{bool}).

           @return: This ETM coordinate as (B{C{LatLon}}) or a
                    L{LatLonDatum5Tuple}C{(lat, lon, datum, gamma,
                    scale)} if B{C{LatLon}} is C{None}.

           @raise ETMError: This ETM coordinate's C{exacTM} and this C{datum}
                            incompatible or no convergence transforming to
                            lat- and longitude.

           @raise TypeError: Invalid or non-ellipsoidal B{C{LatLon}}.

           @example:

            >>> from pygeodesy import ellipsoidalVincenty as eV, Etm
            >>> u = Etm(31, 'N', 448251.795, 5411932.678)
            >>> ll = u.toLatLon(eV.LatLon)  # 48°51′29.52″N, 002°17′40.20″E
        '''
        if not self._latlon or self._latlon._toLLEB_args != (unfalse, self.exactTM):
            self._toLLEB(unfalse=unfalse)
        return self._latlon5(LatLon)

    def _toLLEB(self, unfalse=True, **unused):  # PYCHOK signature
        '''(INTERNAL) Compute (ellipsoidal) lat- and longitude.
        '''
        xTM, d = self.exactTM, self.datum
        # double check that this and exactTM's ellipsoid match
        if xTM._E != d.ellipsoid:  # PYCHOK no cover
            t = repr(d.ellipsoid)
            raise ETMError(repr(xTM._E), txt=_incompatible(t))

        e, n =  self.eastingnorthing2(falsed=not unfalse)
        lon0 = _cmlon(self.zone) if bool(unfalse) == self.falsed else None
        lat, lon, g, k = xTM.reverse(e, n, lon0=lon0)

        ll = _LLEB(lat, lon, datum=d, name=self.name)  # utm._LLEB
        ll._gamma = g
        ll._scale = k
        self._latlon5args(ll, _toBand, unfalse, xTM)

    def toUtm(self):  # PYCHOK signature
        '''Copy this ETM to a UTM coordinate.

           @return: The UTM coordinate (L{Utm}).
        '''
        return self._xcopy2(Utm)


class ExactTransverseMercator(_NamedBase):
    '''Pure Python version of Karney's C++ class U{TransverseMercatorExact
       <https://GeographicLib.SourceForge.io/C++/doc/TransverseMercatorExact_8cpp_source.html>},
       a numerically exact transverse Mercator projection, further referred to as C{TMExact}.
    '''
    _datum     =  None    # Datum
    _E         =  None    # Ellipsoid
    _extendp   =  False   # use extended domain
#   _iteration =  None    # ._sigmaInv2 and ._zetaInv2
    _k0        = _K0_UTM  # central scale factor
    _lon0      = _0_0     # central meridian
    _mu        = _0_0     # ._E.e2, 1st eccentricity squared
    _mv        = _1_0     # _1_0 - ._mu
    _raiser    =  False   # throw Error
    _sigmaC    =  None    # most recent _sigmaInv04 case C{int}
    _zetaC     =  None    # most recent _zetaInv04 case C{int}

    def __init__(self, datum=_WGS84, lon0=0, k0=_K0_UTM, extendp=False, name=NN, raiser=False):
        '''New L{ExactTransverseMercator} projection.

           @kwarg datum: The I{non-spherical} datum or ellipsoid (L{Datum},
                         L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg lon0: Central meridian, default (C{degrees180}).
           @kwarg k0: Central scale factor (C{float}).
           @kwarg extendp: Use the I{extended} domain (C{bool}), I{standard} otherwise.
           @kwarg name: Optional name for the projection (C{str}).
           @kwarg raiser: If C{True}, throw an L{ETMError} for convergence failures (C{bool}).

           @raise ETMError: Near-spherical B{C{datum}} or C{ellipsoid} or invalid B{C{lon0}}
                            or B{C{k0}}.

           @see: U{Constructor TransverseMercatorExact<https://GeographicLib.SourceForge.io/
                 html/classGeographicLib_1_1TransverseMercatorExact.html>} for more details,
                 especially on B{X{extendp}}.

           @note: For all 255.5K U{TMcoords.dat<https://Zenodo.org/record/32470>} tests (with
                  C{0 <= lat <= 84} and C{0 <= lon}) the maximum error is C{5.2e-08 .forward}
                  (or 52 nano-meter) easting and northing and C{3.8e-13 .reverse} (or 0.38
                  pico-degrees) lat- and longitude (with Python 3.7.3+, 2.7.16+, PyPy6 3.5.3
                  and PyPy6 2.7.13, all in 64-bit on macOS 10.13.6 High Sierra C{x86_64} and
                  12.2 Monterey C{arm64} and C{"arm64_x86_64"}).
        '''
        if extendp:
            self._extendp = True
        if name:
            self.name = name
        if raiser:
            self.raiser = True

        self.datum = datum  # invokes ._reset
        self.k0    = k0
        self.lon0  = lon0

    @property_doc_(''' the datum (L{Datum}).''')
    def datum(self):
        '''Get the datum (L{Datum}) or C{None}.
        '''
        return self._datum

    @datum.setter  # PYCHOK setter!
    def datum(self, datum):
        '''Set the datum and ellipsoid (L{Datum}, L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}).

           @raise ETMError: Near-spherical B{C{datum}} or C{ellipsoid}.
        '''
        d = _ellipsoidal_datum(datum, name=self.name)  # raiser=_datum_)
        self._reset(d)
        self._datum = d

    @Property_RO
    def _e(self):
        '''(INTERNAL) Get and cache C{_e}.
        '''
        return self._E.e

    @Property_RO
    def _1_e_90(self):  # PYCHOK no cover
        '''(INTERNAL) Get and cache C{(1 - _e) * 90}.
        '''
        return (_1_0 - self._e) * _90_0

    @property_RO
    def ellipsoid(self):
        '''Get the ellipsoid (L{Ellipsoid}).
        '''
        return self._E

    @Property_RO
    def _e_PI_2(self):
        '''(INTERNAL) Get and cache C{_e * PI / 2}.
        '''
        return self._e * PI_2

    @Property_RO
    def _e_PI_4_(self):
        '''(INTERNAL) Get and cache C{-_e * PI / 4}.
        '''
        return -self._e * PI_4

    @Property_RO
    def _1_e_PI_2(self):
        '''(INTERNAL) Get and cache C{(1 - _e) * PI / 2}.
        '''
        return (_1_0 - self._e) * PI_2

    @Property_RO
    def _1_2e_PI_2(self):
        '''(INTERNAL) Get and cache C{(1 - 2 * _e) * PI / 2}.
        '''
        return (_1_0 - self._e * _2_0) * PI_2

    @property_RO
    def equatoradius(self):
        '''Get the C{ellipsoid}'s equatorial radius, semi-axis (C{meter}).
        '''
        return self._E.a

    a = equatoradius

    @Property_RO
    def _e_TAYTOL(self):
        '''(INTERNAL) Get and cache C{e * TAYTOL}.
        '''
        return self._e * _TAYTOL

    @Property_RO
    def _Eu(self):
        '''(INTERNAL) Get and cache C{Elliptic(_mu)}.
        '''
        return Elliptic(self._mu)

    @Property_RO
    def _Eu_cE(self):
        '''(INTERNAL) Get and cache C{_Eu.cE}.
        '''
        return self._Eu.cE

    def _Eu_2cE_(self, xi):
        '''(INTERNAL) Return C{_Eu.cE * 2 - B{xi}}.
        '''
        return self._Eu_cE * _2_0 - xi

    @Property_RO
    def _Eu_cE_4(self):
        '''(INTERNAL) Get and cache C{_Eu.cE / 4}.
        '''
        return self._Eu_cE / _4_0

    @Property_RO
    def _Eu_cK(self):
        '''(INTERNAL) Get and cache C{_Eu.cK}.
        '''
        return self._Eu.cK

    @Property_RO
    def _Eu_cK_cE(self):
        '''(INTERNAL) Get and cache C{_Eu.cK / _Eu.cE}.
        '''
        return self._Eu_cK / self._Eu_cE

    @Property_RO
    def _Eu_2cK_PI(self):
        '''(INTERNAL) Get and cache C{_Eu.cK * 2 / PI}.
        '''
        return self._Eu_cK / PI_2

    @Property_RO
    def _Ev(self):
        '''(INTERNAL) Get and cache C{Elliptic(_mv)}.
        '''
        return Elliptic(self._mv)

    @Property_RO
    def _Ev_cK(self):
        '''(INTERNAL) Get and cache C{_Ev.cK}.
        '''
        return self._Ev.cK

    @Property_RO
    def _Ev_cKE(self):
        '''(INTERNAL) Get and cache C{_Ev.cKE}.
        '''
        return self._Ev.cKE

    @Property_RO
    def _Ev_3cKE_4(self):
        '''(INTERNAL) Get and cache C{_Ev.cKE * 3 / 4}.
        '''
        return self._Ev_cKE * 0.75  # _0_75

    @Property_RO
    def _Ev_5cKE_4(self):
        '''(INTERNAL) Get and cache C{_Ev.cKE * 5 / 4}.
        '''
        return self._Ev_cKE * 1.25  # _1_25

    @Property_RO
    def extendp(self):
        '''Get the domain (C{bool}), I{extended} or I{standard}.
        '''
        return self._extendp

    @property_RO
    def flattening(self):
        '''Get the C{ellipsoid}'s flattening (C{scalar}).
        '''
        return self._E.f

    f = flattening

    def forward(self, lat, lon, lon0=None, name=NN):  # MCCABE 13
        '''Forward projection, from geographic to transverse Mercator.

           @arg lat: Latitude of point (C{degrees}).
           @arg lon: Longitude of point (C{degrees}).
           @kwarg lon0: Central meridian (C{degrees180}), overriding
                        the default if not C{None}.
           @kwarg name: Optional name (C{str}).

           @return: L{Forward4Tuple}C{(easting, northing, gamma, scale)}.

           @see: C{void TMExact::Forward(real lon0, real lat, real lon,
                                         real &x, real &y,
                                         real &gamma, real &k)}.

           @raise ETMError: No convergence, thrown iff property
                            C{B{raiser}=True}.
        '''
        lat    = _fix90(lat)
        lon, _ = _diff182((self.lon0 if lon0 is None else lon0), lon)
        if self.extendp:
            backside = _lat = _lon = False
        else:  # enforce the parity
            lat, _lat = _unsigned2(lat)
            lon, _lon = _unsigned2(lon)
            backside  =  lon > 90
            if backside:  # PYCHOK no cover
                lon = _loneg(lon)
                if lat == 0:
                    _lat = True

        # u, v = coordinates for the Thompson TM, Lee 54
        if lat == 90:  # isnear90(lat)
            u = self._Eu_cK
            v = self._iteration = self._zetaC = 0
        elif lat == 0 and lon == self._1_e_90:  # PYCHOK no cover
            u = self._iteration = self._zetaC = 0
            v = self._Ev_cK
        else:  # tau = tan(phi), taup = sinh(psi)
            tau, lam = _tand(lat), radians(lon)
            u, v = self._zetaInv2(self._E.es_taupf(tau), lam)

        sncndn6 = self._sncndn6(u, v)
        y, x, _ = self._sigma3(v,  *sncndn6)
        g, k    = (lon, self.k0) if isnear90(lat) else \
                  self._zetaScaled(sncndn6, ll=False)

        if backside:
            y, g = self._Eu_2cE_(y), _loneg(g)
        y *= self._k0_a
        x *= self._k0_a
        if _lat:
            y, g = neg_(y, g)
        if _lon:
            x, g = neg_(x, g)
        return Forward4Tuple(x, y, g, k, iteration=self._iteration,
                                         name=name or self.name)

    def _Inv03(self, psi, dlam, _3_mv_e):  # (xi, deta, _3_mv)
        '''(INTERNAL) Partial C{_zetaInv04} or C{_sigmaInv04}, Case 2
        '''
        # atan2(dlam-psi, psi+dlam) + 45d gives arg(zeta - zeta0) in
        # range [-135, 225).  Subtracting 180 (multiplier is negative)
        # makes range [-315, 45).  Multiplying by 1/3 (for cube root)
        # gives range [-105, 15).  In particular the range [-90, 180]
        # in zeta space maps to [-90, 0] in w space as required.
        a = atan2(dlam - psi, psi + dlam) / _3_0 - PI_4
        s, c = sincos2(a)
        h = hypot(psi, dlam)
        r = cbrt(h * _3_mv_e)
        u = r * c
        v = r * s + self._Ev_cK
        # Error using this guess is about 0.068 * rad^(5/3)
        return u, v, h

    @property_RO
    def iteration(self):
        '''Get the most recent C{ExactTransverseMercator.forward}
           or C{ExactTransverseMercator.reverse} iteration number
           (C{int}) or C{None} if not available/applicable.
        '''
        return self._iteration

    @property_doc_(''' the central scale factor (C{float}).''')
    def k0(self):
        '''Get the central scale factor (C{float}), aka I{C{scale0}}.
        '''
        return self._k0  # aka scale0

    @k0.setter  # PYCHOK setter!
    def k0(self, k0):
        '''Set the central scale factor (C{float}), aka I{C{scale0}}.

           @raise ETMError: Invalid B{C{k0}}.
        '''
        k0 = Scalar_(k0=k0, Error=ETMError, low=_TOL_10, high=_1_0)
        if self._k0 != k0:
            ExactTransverseMercator._k0_a._update(self)  # redo ._k0_a
            self._k0 = k0

    @Property_RO
    def _k0_a(self):
        '''(INTERNAL) Get and cache C{k0 * equatoradius}.
        '''
        return self.k0 * self.equatoradius

    @property_doc_(''' the central meridian (C{degrees180}).''')
    def lon0(self):
        '''Get the central meridian (C{degrees180}).
        '''
        return self._lon0

    @lon0.setter  # PYCHOK setter!
    def lon0(self, lon0):
        '''Set the central meridian (C{degrees180}).

           @raise ETMError: Invalid B{C{lon0}}.
        '''
        self._lon0 = _norm180(Degrees(lon0=lon0, Error=ETMError))

    @deprecated_property_RO
    def majoradius(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{equatoradius}.'''
        return self.equatoradius

    @Property_RO
    def _1_mu_2(self):
        '''(INTERNAL) Get and cache C{_mu / 2 + 1}.
        '''
        return self._mu * _0_5 + _1_0

    @Property_RO
    def _3_mv(self):
        '''(INTERNAL) Get and cache C{3 / _mv}.
        '''
        return _3_0 / self._mv

    @Property_RO
    def _3_mv_e(self):
        '''(INTERNAL) Get and cache C{3 / (_mv * _e)}.
        '''
        return _3_0 / (self._mv * self._e)

    def _Newton2(self, taup, lam, u, v, C, *psi):  # or (xi, eta, u, v)
        '''(INTERNAL) Invert C{_zetaInv2} or C{_sigmaInv2} using Newton's method.

           @return: 2-Tuple C{(u, v)}.

           @raise ETMError: No convergence.
        '''
        sca1, tol2 = _1_0, _TOL_10
        if psi:  # _zetaInv2
            sca1 = sca1 / hypot1(taup)  # /= chokes PyChecker
            tol2 = tol2 / max(psi[0], _1_0)**2

            _zeta3    = self._zeta3
            _zetaDwd2 = self._zetaDwd2
        else:  # _sigmaInv2
            _zeta3    = self._sigma3
            _zetaDwd2 = self._sigmaDwd2

        d2, r = tol2, self.raiser
        _U_2  = Fsum(u).fsum2_
        _V_2  = Fsum(v).fsum2_
        # min iterations 2, max 6 or 7, mean 3.9 or 4.0
        for i in range(1, _TRIPS):  # GEOGRAPHICLIB_PANIC
            sncndn6 =  self._sncndn6(u, v)
            du, dv  = _zetaDwd2(*sncndn6)
            T, L, _ = _zeta3(v, *sncndn6)
            T  = (taup - T) * sca1
            L -=  lam
            u, dU = _U_2(T * du,  L * dv)
            v, dV = _V_2(T * dv, -L * du)
            if d2 < tol2:
                r = False
                break
            d2 = hypot2(dU, dV)

        self._iteration = i
        if r:  # PYCHOK no cover
            n = callername(up=2, underOK=True)
            t = unstr(n, taup, lam, u, v, C=C)
            raise ETMError(Fmt.no_convergence(d2, tol2), txt=t)
        return u, v

    @property_doc_(''' raise an L{ETMError} for convergence failures (C{bool}).''')
    def raiser(self):
        '''Get the error setting (C{bool}).
        '''
        return self._raiser

    @raiser.setter  # PYCHOK setter!
    def raiser(self, raiser):
        '''Set the error setting (C{bool}), if C{True} throw an L{ETMError}
           for convergence failures.
        '''
        self._raiser = bool(raiser)

    def _reset(self, datum):
        '''(INTERNAL) Set the ellipsoid and elliptic moduli.

           @arg datum: Ellipsoidal datum (C{Datum}).

           @raise ETMError: Near-spherical B{C{datum}} or C{ellipsoid}.
        '''
        E  = datum.ellipsoid
        mu = E.e2   # .eccentricity1st2
        mv = E.e21  # _1_0 - mu
        if isnear0(E.e) or isnear0(mu, eps0=EPS02) \
                        or isnear0(mv, eps0=EPS02):  # or sqrt(mu) != E.e
            raise ETMError(ellipsoid=E, txt=_near_(_spherical_))

        if self._datum or self._E:
            _i = ExactTransverseMercator.iteration._uname
            _update_all(self, _i, '_sigmaC', '_zetaC')  # _under

        self._E  = E
        self._mu = mu
        self._mv = mv

    def reverse(self, x, y, lon0=None, name=NN):
        '''Reverse projection, from Transverse Mercator to geographic.

           @arg x: Easting of point (C{meters}).
           @arg y: Northing of point (C{meters}).
           @kwarg lon0: Central meridian (C{degrees180}), overriding
                        the default if not C{None}.
           @kwarg name: Optional name (C{str}).

           @return: L{Reverse4Tuple}C{(lat, lon, gamma, scale)}.

           @see: C{void TMExact::Reverse(real lon0, real x, real y,
                                         real &lat, real &lon,
                                         real &gamma, real &k)}

           @raise ETMError: No convergence, thrown iff property
                            C{B{raiser}=True}.
        '''
        # undoes the steps in .forward.
        xi  = y / self._k0_a
        eta = x / self._k0_a
        if self.extendp:
            backside = _lat = _lon = False
        else:  # enforce the parity
            eta, _lon = _unsigned2(eta)
            xi,  _lat = _unsigned2(xi)
            backside  =  xi > self._Eu_cE
            if backside:  # PYCHOK no cover
                xi = self._Eu_2cE_(xi)

        # u, v = coordinates for the Thompson TM, Lee 54
        if xi or eta != self._Ev_cKE:
            u, v = self._sigmaInv2(xi, eta)
        else:  # PYCHOK no cover
            u = self._iteration = self._sigmaC = 0
            v = self._Ev_cK

        if v or u != self._Eu_cK:
            g, k, lat, lon = self._zetaScaled(self._sncndn6(u, v))
        else:  # PYCHOK no cover
            g, k, lat, lon = _0_0, self.k0, _90_0, _0_0

        if backside:  # PYCHOK no cover
            lon, g = _loneg(lon), _loneg(g)
        if _lat:
            lat, g = neg_(lat, g)
        if _lon:
            lon, g = neg_(lon, g)
        lon += self.lon0 if lon0 is None else _norm180(lon0)
        return Reverse4Tuple(lat, _norm180(lon), g, k,  # _norm180(lat)
                                   iteration=self._iteration,
                                   name=name or self.name)

    def _scaled2(self, tau, d2, snu, cnu, dnu, snv, cnv, dnv):
        '''(INTERNAL) C{scaled}.

           @note: Argument B{C{d2}} is C{_mu * cnu**2 + _mv * cnv**2}
                  from C{._zeta3}.

           @return: 2-Tuple C{(convergence, scale)}.

           @see: C{void TMExact::Scale(real tau, real /*lam*/,
                                       real snu, real cnu, real dnu,
                                       real snv, real cnv, real dnv,
                                       real &gamma, real &k)}.
        '''
        mu, mv = self._mu, self._mv
        cnudnv = cnu * dnv
        # Lee 55.12 -- negated for our sign convention.  g gives
        # the bearing (clockwise from true north) of grid north
        g = atan2d(mv * cnv * snv * snu, cnudnv * dnu)
        # Lee 55.13 with nu given by Lee 9.1 -- in sqrt change
        # the numerator from (1 - snu^2 * dnv^2) to (_mv * snv^2
        # + cnu^2 * dnv^2) to maintain accuracy near phi = 90
        # and change the denomintor from (dnu^2 + dnv^2 - 1) to
        # (_mu * cnu^2 + _mv * cnv^2) to maintain accuracy near
        # phi = 0, lam = 90 * (1 - e).  Similarly rewrite sqrt in
        # 9.1 as _mv + _mu * c^2 instead of 1 - _mu * sin(phi)^2
        if d2 > 0:
            # originally: sec2 = 1 + tau**2  # sec(phi)^2
            #             d2   = (mu * cnu**2 + mv * cnv**2)
            #             q2   = (mv * snv**2 + cnudnv**2) / d2
            # k = sqrt(mv + mu / sec2) * sqrt(sec2) * sqrt(q2)
            #   = sqrt(mv * sec2 + mu) * sqrt(q2)
            #   = sqrt(mv + mv * tau**2 + mu) * sqrt(q2)
            k, q2 = _0_0, (mv * snv**2 + cnudnv**2)
            if q2 > 0:
                k2 = fsum1f_(mu, mv, mv * tau**2)
                if k2 > 0:
                    k = sqrt(k2) * sqrt(q2 / d2) * self.k0
        else:
            k = _OVERFLOW
        return g, k

    def _sigma3(self, v, snu, cnu, dnu, snv, cnv, dnv):
        '''(INTERNAL) C{sigma}.

           @return: 3-Tuple C{(xi, eta, d2)}.

           @see: C{void TMExact::sigma(real /*u*/, real snu, real cnu, real dnu,
                                       real   v,   real snv, real cnv, real dnv,
                                       real &xi, real &eta)}.

           @raise ETMError: No convergence.
        '''
        mu = self._mu * cnu
        mv = self._mv * cnv
        # Lee 55.4 writing
        # dnu^2 + dnv^2 - 1 = _mu * cnu^2 + _mv * cnv^2
        d2  = cnu * mu + cnv * mv
        mu *= snu * dnu
        mv *= snv * dnv
        if d2 > 0:  # /= chokes PyChecker
            mu = mu / d2
            mv = mv / d2
        else:
            mu, mv = map1(_overflow, mu, mv)
        xi = self._Eu.fE(snu, cnu, dnu) - mu
        v -= self._Ev.fE(snv, cnv, dnv) - mv
        return xi, v, d2

    def _sigmaDwd2(self, snu, cnu, dnu, snv, cnv, dnv):
        '''(INTERNAL) C{sigmaDwd}.

           @return: 2-Tuple C{(du, dv)}.

           @see: C{void TMExact::dwdsigma(real /*u*/, real snu, real cnu, real dnu,
                                          real /*v*/, real snv, real cnv, real dnv,
                                          real &du, real &dv)}.
        '''
        snuv = snu * snv
        # Reciprocal of 55.9: dw / ds = dn(w)^2/_mv,
        # expanding complex dn(w) using A+S 16.21.4
        d = self._mv * (cnv**2 + self._mu * snuv**2)**2
        r = cnv * dnu * dnv
        i = cnu * snuv * self._mu
        du = (r**2 - i**2) / d
        dv = neg(r * i * _2_0 / d)
        return du, dv

    def _sigmaInv2(self, xi, eta):
        '''(INTERNAL) Invert C{sigma} using Newton's method.

           @return: 2-Tuple C{(u, v)}.

           @see: C{void TMExact::sigmainv(real xi, real eta,
                                          real &u, real &v)}.

           @raise ETMError: No convergence.
        '''
        u, v, t, self._sigmaC = self._sigmaInv04(xi, eta)
        if not t:
            u, v = self._Newton2(xi, eta, u, v, self._sigmaC)
        return u, v

    def _sigmaInv04(self, xi, eta):
        '''(INTERNAL) Starting point for C{sigmaInv}.

           @return: 4-Tuple C{(u, v, trip, Case)}.

           @see: C{bool TMExact::sigmainv0(real xi, real eta,
                                           real &u, real &v)}.
        '''
        t = False
        d = eta - self._Ev_cKE
        if eta > self._Ev_5cKE_4 or (xi < d and xi < -self._Eu_cE_4):
            # sigma as a simple pole at
            #  w = w0 = Eu.K() + i * Ev.K()
            # and sigma is approximated by
            #  sigma = (Eu.E() + i * Ev.KE()) + 1 / (w - w0)
            u, v = _norm2(xi - self._Eu_cE, -d)
            u += self._Eu_cK
            v += self._Ev_cK
            C  = 1

        elif (eta > self._Ev_3cKE_4 and xi < self._Eu_cE_4) or d > 0:
            # At w = w0 = i * Ev.K(), we have
            #  sigma  = sigma0  = i * Ev.KE()
            #  sigma' = sigma'' = 0
            # including the next term in the Taylor series gives:
            #  sigma = sigma0 - _mv / 3 * (w - w0)^3
            # When inverting this, we map arg(w - w0) = [-pi/2, -pi/6]
            # to arg(sigma - sigma0) = [-pi/2, pi/2] mapping arg =
            # [-pi/2, -pi/6] to [-pi/2, pi/2]
            u, v, h = self._Inv03(xi, d, self._3_mv)
            t = h < _TAYTOL2
            C = 2

        else:  # use w = sigma * Eu.K/Eu.E (correct in limit _e -> 0)
            u = v = self._Eu_cK_cE
            u *= xi
            v *= eta
            C  = 3

        return u, v, t, C

    def _sncndn6(self, u, v):
        '''(INTERNAL) Get 6-tuple C{(snu, cnu, dnu, snv, cnv, dnv)}.
        '''
        # snu, cnu, dnu = self._Eu.sncndn(u)
        # snv, cnv, dnv = self._Ev.sncndn(v)
        return self._Eu.sncndn(u) + self._Ev.sncndn(v)

    def toStr(self, joined=_COMMASPACE_, **kwds):  # PYCHOK signature
        '''Return a C{str} representation.

           @kwarg joined: Separator to join the attribute strings
                         (C{str} or C{None} or C{NN} for non-joined).
           @kwarg kwds: Optional, overriding keyword arguments.
        '''
        d = dict(datum=self.datum.name, lon0=self.lon0,
                 k0=self.k0, extendp=self.extendp)
        if self.name:
            d.update(name=self.name)
        t = pairs(d, **kwds)
        return joined.join(t) if joined else t

    def _zeta3(self, unused, snu, cnu, dnu, snv, cnv, dnv):  # _sigma3 signature
        '''(INTERNAL) C{zeta}.

           @return: 3-Tuple C{(taup, lambda, d2)}.

           @see: C{void TMExact::zeta(real /*u*/, real snu, real cnu, real dnu,
                                      real /*v*/, real snv, real cnv, real dnv,
                                      real &taup, real &lam)}
        '''
        e, cnu2, mv = self._e, cnu**2, self._mv
        # Overflow value like atan(overflow) = pi/2
        t1 = t2 = _overflow(snu)
        # Lee 54.17 but write
        # atanh(snu * dnv)      = asinh(snu * dnv / sqrt(cnu^2 + _mv * snu^2 * snv^2))
        # atanh(_e * snu / dnv) = asinh(_e * snu / sqrt(_mu * cnu^2 + _mv * cnv^2))
        d1 = cnu2 + mv * (snu * snv)**2
        if d1 > EPS02:  # _EPSmin
            t1 = snu * dnv / sqrt(d1)
        else:
            d1 = 0
        d2 = self._mu * cnu2 + mv * cnv**2
        if d2 > EPS02:  # _EPSmin
            t2 = sinh(e * asinh(e * snu / sqrt(d2)))
        else:
            d2 = 0
        # psi = asinh(t1) - asinh(t2)
        # taup = sinh(psi)
        taup = t1 * hypot1(t2) - t2 * hypot1(t1)
        lam  = (atan2(dnu * snv,     cnu * cnv) -
                atan2(cnu * snv * e, dnu * cnv) * e) if d1 and d2 else _0_0
        return taup, lam, d2

    def _zetaDwd2(self, snu, cnu, dnu, snv, cnv, dnv):
        '''(INTERNAL) C{zetaDwd}.

           @return: 2-Tuple C{(du, dv)}.

           @see: C{void TMExact::dwdzeta(real /*u*/, real snu, real cnu, real dnu,
                                         real /*v*/, real snv, real cnv, real dnv,
                                         real &du, real &dv)}.
        '''
        cnu2  = cnu**2 * self._mu
        cnv2  = cnv**2
        dnuv  = dnu * dnv
        dnuv2 = dnuv**2
        snuv  = snu * snv
        snuv2 = snuv**2 * self._mu
        # Lee 54.21 but write (see A+S 16.21.4)
        # (1 - dnu^2 * snv^2) = (cnv^2 + _mu * snu^2 * snv^2)
        d  = self._mv   * (cnv2 + snuv2)**2  # max(d, EPS02)?
        du = cnu * dnuv * (cnv2 - snuv2) / d
        dv = cnv * snuv * (cnu2 + dnuv2) / d
        return du, neg(dv)

    def _zetaInv2(self, taup, lam):
        '''(INTERNAL) Invert C{zeta} using Newton's method.

           @return: 2-Tuple C{(u, v)}.

           @see: C{void TMExact::zetainv(real taup, real lam,
                                         real &u, real &v)}.

           @raise ETMError: No convergence.
        '''
        psi = asinh(taup)
        u, v, t, self._zetaC = self._zetaInv04(psi, lam)
        if not t:
            u, v = self._Newton2(taup, lam, u, v, self._zetaC, psi)
        return u, v

    def _zetaInv04(self, psi, lam):
        '''(INTERNAL) Starting point for C{zetaInv}.

           @return: 4-Tuple C{(u, v, trip, Case)}.

           @see: C{bool TMExact::zetainv0(real psi, real lam,  # radians
                                          real &u, real &v)}.
        '''
        if lam > self._1_2e_PI_2:
            d = lam - self._1_e_PI_2
            if psi < d and psi < self._e_PI_4_:  # PYCHOK no cover
                # N.B. this branch is normally *not* taken because psi < 0
                # is converted psi > 0 by .forward.  There's a log singularity
                # at w = w0 = Eu.K() + i * Ev.K(), corresponding to the south
                # pole, where we have, approximately
                #  psi = _e + i * pi/2 - _e * atanh(cos(i * (w - w0)/(1 + _mu/2)))
                # Inverting this gives:
                e = self._e  # eccentricity
                s, c = sincos2((PI_2 - lam) / e)
                h, r = sinh(_1_0 - psi / e), self._1_mu_2
                u = self._Eu_cK - r * asinh(s / hypot(c, h))
                v = self._Ev_cK - r * atan2(c, h)
                return u, v, False, 1

            elif psi < self._e_PI_2:
                # At w = w0 = i * Ev.K(), we have
                #  zeta  = zeta0  = i * (1 - _e) * pi/2
                #  zeta' = zeta'' = 0
                # including the next term in the Taylor series gives:
                #  zeta = zeta0 - (_mv * _e) / 3 * (w - w0)^3
                # When inverting this, we map arg(w - w0) = [-90, 0]
                # to arg(zeta - zeta0) = [-90, 180]
                u, v, h = self._Inv03(psi, d, self._3_mv_e)
                return u, v, (h < self._e_TAYTOL), 2

        # Use spherical TM, Lee 12.6 -- writing C{atanh(sin(lam) /
        # cosh(psi)) = asinh(sin(lam) / hypot(cos(lam), sinh(psi)))}.
        # This takes care of the log singularity at C{zeta = Eu.K()},
        # corresponding to the north pole.
        s, c = sincos2(lam)
        h, r = sinh(psi), self._Eu_2cK_PI
        # But scale to put 90, 0 on the right place
        u = r * atan2(h, c)
        v = r * asinh(s / hypot(h, c))
        return u, v, False, 3

    def _zetaScaled(self, sncndn6, ll=True):
        '''(INTERNAL) Recompute (T, L) from (u, v) to improve accuracy of Scale.

           @arg sncndn6: 6-Tuple C{(snu, cnu, dnu, snv, cnv, dnv)}.

           @return: 2-Tuple C{(g, k)} if not C{B{ll}} else
                    4-tuple C{(g, k, lat, lon)}.
        '''
        t, lam, d2 = self._zeta3(None, *sncndn6)
        tau = self._E.es_tauf(t)
        g_k = self._scaled2(tau, d2, *sncndn6)
        if ll:
            g_k += atand(tau), degrees(lam)
        return g_k  # or (g, k, lat, lon)


def parseETM5(strUTM, datum=_WGS84, Etm=Etm, falsed=True, name=NN):
    '''Parse a string representing a UTM coordinate, consisting
       of C{"zone[band] hemisphere easting northing"}.

       @arg strUTM: A UTM coordinate (C{str}).
       @kwarg datum: Optional datum to use (L{Datum}, L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}).
       @kwarg Etm: Optional class to return the UTM coordinate
                   (L{Etm}) or C{None}.
       @kwarg falsed: Both easting and northing are C{falsed} (C{bool}).
       @kwarg name: Optional B{C{Etm}} name (C{str}).

       @return: The UTM coordinate (B{C{Etm}}) or if B{C{Etm}} is
                C{None}, a L{UtmUps5Tuple}C{(zone, hemipole, easting,
                northing, band)}.  The C{hemipole} is the hemisphere
                C{'N'|'S'}.

       @raise ETMError: Invalid B{C{strUTM}}.

       @raise TypeError: Invalid or near-spherical B{C{datum}}.

       @example:

        >>> u = parseETM5('31 N 448251 5411932')
        >>> u.toRepr()  # [Z:31, H:N, E:448251, N:5411932]
        >>> u = parseETM5('31 N 448251.8 5411932.7')
        >>> u.toStr()  # 31 N 448252 5411933
    '''
    r = _parseUTM5(strUTM, datum, Etm, falsed, Error=ETMError, name=name)
    return r


def toEtm8(latlon, lon=None, datum=None, Etm=Etm, falsed=True,
                                         name=NN, strict=True,
                                         zone=None, **cmoff):
    '''Convert a geodetic lat-/longitude to an ETM coordinate.

       @arg latlon: Latitude (C{degrees}) or an (ellipsoidal)
                    geodetic C{LatLon} instance.
       @kwarg lon: Optional longitude (C{degrees}) or C{None}.
       @kwarg datum: Optional datum for the ETM coordinate,
                     overriding B{C{latlon}}'s datum (L{Datum},
                     L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}).
       @kwarg Etm: Optional class to return the ETM coordinate
                   (L{Etm}) or C{None}.
       @kwarg falsed: False both easting and northing (C{bool}).
       @kwarg name: Optional B{C{Utm}} name (C{str}).
       @kwarg strict: Restrict B{C{lat}} to UTM ranges (C{bool}).
       @kwarg zone: Optional UTM zone to enforce (C{int} or C{str}).
       @kwarg cmoff: DEPRECATED, use B{C{falsed}}.  Offset longitude
                     from the zone's central meridian (C{bool}).

       @return: The ETM coordinate as an B{C{Etm}} instance or a
                L{UtmUps8Tuple}C{(zone, hemipole, easting, northing,
                band, datum, gamma, scale)} if B{C{Etm}} is C{None}
                or not B{C{falsed}}.  The C{hemipole} is the C{'N'|'S'}
                hemisphere.

       @raise ETMError: No convergence transforming to ETM easting
                        and northing.

       @raise ETMError: Invalid B{C{zone}} or near-spherical or
                        incompatible B{C{datum}} or C{ellipsoid}.

       @raise RangeError: If B{C{lat}} outside the valid UTM bands or
                          if B{C{lat}} or B{C{lon}} outside the valid
                          range and L{pygeodesy.rangerrors} set to C{True}.

       @raise TypeError: Invalid or near-spherical B{C{datum}} or
                         B{C{latlon}} not ellipsoidal.

       @raise ValueError: The B{C{lon}} value is missing or B{C{latlon}}
                          is invalid.
    '''
    z, B, lat, lon, d, f, name = _to7zBlldfn(latlon, lon, datum,
                                             falsed, name, zone,
                                             strict, ETMError, **cmoff)
    lon0 = _cmlon(z) if f else None
    x, y, g, k = d.exactTM.forward(lat, lon, lon0=lon0)

    return _toXtm8(Etm, z, lat, x, y, B, d, g, k, f,
                        name, latlon, d.exactTM, Error=ETMError)


if __name__ == '__main__':  # MCCABE 13

    from pygeodesy import fstr, KTransverseMercator, printf
    from sys import argv, exit as _exit

    # mimick some of I{Karney}'s utility C{TransverseMercatorProj}
    _f = _r = _s = _t = False
    _as = argv[1:]
    while _as and _as[0].startswith(_DASH_):
        _a = _as.pop(0)
        if len(_a) < 2:
            _exit('%s: option %r invalid' % (_usage(*argv), _a))
        elif '-forward'.startswith(_a):
            _f, _r = True, False
        elif '-reverse'.startswith(_a):
            _f, _r = False, True
        elif '-series'.startswith(_a):
            _s, _t = True, False
        elif _a == '-t':
            _s, _t = False, True
        elif '-help'.startswith(_a):
            _exit(_usage(argv[0], '[-s | -t]',
                                  '[-f[orward] <lat> <lon>',
                                  '| -r[everse] <easting> <northing>',
                                  '| <lat> <lon>]',
                                  '| -h[elp]'))
        else:
            _exit('%s: option %r not supported' % (_usage(*argv), _a))
    if len(_as) > 1:
        f2 = map1(float, *_as[:2])
    else:
        _exit('%s ...: incomplete' % (_usage(*argv),))

    if _s:  # -series
        tm = KTransverseMercator()
    else:
        tm = ExactTransverseMercator(extendp=_t)

    if _f:
        t = tm.forward(*f2)
    elif _r:
        t = tm.reverse(*f2)
    else:
        t = tm.forward(*f2)
        printf('%s: %s', tm.classname, fstr(t, sep=_SPACE_))
        t = tm.reverse(t.easting, t.northing)
    printf('%s: %s', tm.classname, fstr(t, sep=_SPACE_))


# % python3 -m pygeodesy.etm 33.33 44.44
# ExactTransverseMercator: 4276926.114804 4727193.767015 28.375537 1.233325
# ExactTransverseMercator: 33.33 44.44 28.375537 1.233325

# % python3 -m pygeodesy.etm -s 33.33 44.44
# KTransverseMercator: 4276926.114804 4727193.767015 28.375537 1.233325
# KTransverseMercator: 33.33 44.44 28.375537 1.233325

# % echo 33.33 44.44 | .../bin/TransverseMercatorProj
# 4276926.114804 4727193.767015 28.375536563148 1.233325101778

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
