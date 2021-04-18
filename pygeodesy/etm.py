
# -*- coding: utf-8 -*-

u'''Classes L{ETMError} and L{Etm}, a pure Python transcription of
I{Charles Karney}'s C++ class U{TransverseMercatorExact
<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1TransverseMercatorExact.html>},
abbreviated as C{TMExact} below.

Python class L{ExactTransverseMercator} implements the C{Exact Transverse
Mercator} (ETM) projection.  Instances of class L{Etm} represent ETM
C{easting, nothing} locations.

Following is a copy of I{Karney}'s U{TransverseMercatorExact.hpp
<https://GeographicLib.SourceForge.io/html/TransverseMercatorExact_8hpp_source.html>}
file C{Header}.

Copyright (C) U{Charles Karney<mailto:Charles@Karney.com>} (2008-2017)
and licensed under the MIT/X11 License.  For more information, see the
U{GeographicLib<https://GeographicLib.SourceForge.io>} documentation.

The method entails using the C{Thompson Transverse Mercator} as an
intermediate projection.  The projections from the intermediate
coordinates to C{phi, lam} and C{x, y} are given by elliptic functions.
The inverse of these projections are found by Newton's method with a
suitable starting guess.

The relevant section of L.P. Lee's paper U{Conformal Projections Based On
Jacobian Elliptic Functions<https://DOI.org/10.3138/X687-1574-4325-WM62>}
is part V, pp 67--101.  The C++ implementation and notation closely
follow Lee, with the following exceptions::

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

Minor alterations have been made in some of Lee's expressions in an
attempt to control round-off.  For example, C{atanh(sin(phi))} is
replaced by C{asinh(tan(phi))} which maintains accuracy near
C{phi = pi/2}.  Such changes are noted in the code.

@newfield example: Example, Examples
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division

from pygeodesy.basics import copysign, neg, neg_, _xinstanceof
from pygeodesy.datums import _ellipsoidal_datum, _WGS84
from pygeodesy.elliptic import Elliptic, EllipticError, _TRIPS
from pygeodesy.errors import _incompatible
from pygeodesy.fmath import cbrt, Fsum, fsum_, hypot, hypot1, hypot2
from pygeodesy.interns import EPS, _1_EPS, NN, PI_2, PI_4, \
                             _COMMASPACE_, _convergence_, _easting_, \
                             _EPS0__2, _lat_, _lon_, _no_, _northing_, \
                             _scale_, _0_0, _0_1, _0_25, _0_5, _1_0, \
                             _2_0, _3_0, _90_0, _180_0
from pygeodesy.interns import _lon0_  # PYCHOK used!
from pygeodesy.karney import _diff182, _fix90, _norm180
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import _NamedBase, _NamedTuple
from pygeodesy.props import deprecated_method, deprecated_Property_RO, \
                            Property_RO, property_RO, property_doc_
from pygeodesy.streprs import pairs, unstr
from pygeodesy.units import Degrees, Easting, Lat,Lon, Northing, \
                            Scalar, Scalar_
from pygeodesy.utily import atand, atan2d, sincos2
from pygeodesy.utm import _cmlon, _K0_UTM, _LLEB, _parseUTM5, \
                           Utm, UTMError, _toXtm8, _to7zBlldfn

from math import asinh, atan2, degrees, radians, sinh, sqrt, tan

__all__ = _ALL_LAZY.etm
__version__ = '21.04.15'

_OVERFLOW = _1_EPS**2  # about 2e+31
_TOL_10   = _0_1 * EPS
_TAYTOL   =  pow(EPS, 0.6)
_TAYTOL2  = _2_0 * _TAYTOL


class EasNorExact4Tuple(_NamedTuple):
    '''4-Tuple C{(easting, northing, convergence, scale)} in
       C{meter}, C{meter}, C{degrees} and C{scalar}.
    '''
    _Names_ = (_easting_, _northing_, _convergence_, _scale_)
    _Units_ = ( Easting,   Northing,   Degrees,       Scalar)


class ETMError(UTMError):
    '''Exact Transverse Mercator (ETM) parse, projection or other L{Etm} issue.
    '''
    pass


class Etm(Utm):
    '''Exact Transverse Mercator (ETM) coordinate, a sub-class of
       L{Utm}, a Universal Transverse Mercator (UTM) coordinate
       using the L{ExactTransverseMercator} projection for highest
       accuracy.

       @note: Conversion of L{Etm} coordinates to and from (geodetic)
              lat- and longitude is 3-4 times slower than L{Utm}.

       @see: Karney's U{Detailed Description<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1TransverseMercatorExact.html#details>}.
    '''
    _Error   = ETMError
    _exactTM = None

    def __init__(self, zone, hemisphere, easting, northing, band=NN,  # PYCHOK expected
                             datum=_WGS84, falsed=True,
                             convergence=None, scale=None, name=NN):
        '''New L{Etm} coordinate.

           @arg zone: Longitudinal UTM zone (C{int}, 1..60) or zone
                      with/-out (latitudinal) Band letter (C{str},
                      '01C'..'60X').
           @arg hemisphere: Northern or southern hemisphere (C{str},
                            C{'N[orth]'} or C{'S[outh]'}).
           @arg easting: Easting, see B{C{falsed}} (C{meter}).
           @arg northing: Northing, see B{C{falsed}} (C{meter}).
           @kwarg band: Optional, (latitudinal) band (C{str}, 'C'..'X').
           @kwarg datum: Optional, this coordinate's datum (L{Datum},
                         L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg falsed: Both B{C{easting}} and B{C{northing}} are
                          falsed (C{bool}).
           @kwarg convergence: Optional meridian convergence, bearing
                               off grid North, clockwise from true North
                               (C{degrees}) or C{None}.
           @kwarg scale: Optional grid scale factor (C{scalar}) or C{None}.
           @kwarg name: Optional name (C{str}).

           @raise ETMError: Invalid B{C{zone}}, B{C{hemishere}} or
                            B{C{band}}.

           @raise TypeError: Invalid B{C{datum}}.

           @example:

            >>> import pygeodesy
            >>> u = pygeodesy.Etm(31, 'N', 448251, 5411932)
        '''
        Utm.__init__(self, zone, hemisphere, easting, northing,
                                 band=band, datum=datum, falsed=falsed,
                                 convergence=convergence, scale=scale,
                                 name=name)
        self.exactTM = self.datum.exactTM  # ExactTransverseMercator(datum=self.datum)

    @property_doc_(''' the ETM projection (L{ExactTransverseMercator}).''')
    def exactTM(self):
        '''Get the ETM projection (L{ExactTransverseMercator}).
        '''
        return self._exactTM

    @exactTM.setter  # PYCHOK setter!
    def exactTM(self, exactTM):
        '''Set the ETM projection (L{ExactTransverseMercator}).
        '''
        _xinstanceof(ExactTransverseMercator, exactTM=exactTM)

        E = self.datum.ellipsoid
        if exactTM._E != E or exactTM.equatoradius != E.a \
                           or exactTM.flattening   != E.f:
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

           @see: Function L{parseUPS5}, L{parseUTM5} and L{parseUTMUPS5}.
        '''
        return parseETM5(strETM, datum=self.datum, Etm=self.classof,
                                 name=name or self.name)

    @deprecated_method
    def parseETM(self, strETM):
        '''DEPRECATED, use method L{Etm.parse}.
        '''
        return self.parse(strETM)

    def toLatLon(self, LatLon=None, unfalse=True, **unused):  # PYCHOK expected
        '''Convert this ETM coordinate to an (ellipsoidal) geodetic point.

           @kwarg LatLon: Optional, ellipsoidal class to return the
                          geodetic point (C{LatLon}) or C{None}.
           @kwarg unfalse: Unfalse B{C{easting}} and B{C{northing}} if
                           falsed (C{bool}).

           @return: This ETM coordinate as (B{C{LatLon}}) or a
                    L{LatLonDatum5Tuple}C{(lat, lon, datum,
                    convergence, scale)} if B{C{LatLon}} is C{None}.

           @raise EllipticError: No convergence.

           @raise TypeError: If B{C{LatLon}} is not ellipsoidal.

           @example:

            >>> from pygeodesy import ellipsoidalVincenty as eV, Etm
            >>> u = Etm(31, 'N', 448251.795, 5411932.678)
            >>> ll = u.toLatLon(eV.LatLon)  # 48°51′29.52″N, 002°17′40.20″E
        '''
        xTM, d = self.exactTM, self.datum
        # double check that this and exactTM's ellipsoids stil match
        if xTM._E != d.ellipsoid:
            t = repr(d.ellipsoid)
            raise ETMError(repr(xTM._E), txt=_incompatible(t))

        if self._latlon and self._latlon_args == (xTM, unfalse):
            return self._latlon5(LatLon)

        f = not unfalse
        e, n = self.eastingnorthing2(falsed=f)
        # f  = unfalse == self.falsed
        #   == unfalse and self.falsed or (not unfalse and not self.falsed)
        #   == unfalse if self.falsed else not unfalse
        #   == unfalse if self.falsed else f
        if self.falsed:
            f = unfalse
        lon0 = _cmlon(self.zone) if f else None
        lat, lon, g, k = xTM.reverse(e, n, lon0=lon0)

        ll = _LLEB(lat, lon, datum=d, name=self.name)
        ll._convergence = g
        ll._scale = k

        self._latlon_to(ll, xTM, unfalse)
        return self._latlon5(LatLon)

    def _latlon_to(self, ll, xTM, unfalse):
        '''(INTERNAL) See C{.toLatLon}, C{toEtm8}, C{_toXtm8}.
        '''
        self._latlon, self._latlon_args = ll, (xTM, unfalse)

    def toUtm(self):  # PYCHOK signature
        '''Copy this ETM to a UTM coordinate.

           @return: The UTM coordinate (L{Utm}).
        '''
        return self._xcopy2(Utm)


class ExactTransverseMercator(_NamedBase):
    '''A Python version of Karney's U{TransverseMercatorExact
       <https://GeographicLib.SourceForge.io/html/TransverseMercatorExact_8cpp_source.html>}
       C++ class, a numerically exact transverse mercator projection, here referred to as
       C{TMExact}.

       @see: C{U{TMExact(real a, real f, real k0, bool extendp)<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1TransverseMercatorExact.html#a72ffcc89eee6f30a6d1f4d061518a6f1>}}.
    '''
    _datum     =  None  # Datum
    _E         =  None  # Ellipsoid
    _extendp   =  True
    _iteration =  None  # ._sigmaInv and ._zetaInv
    _k0        = _1_0   # central scale factor
    _k0_a      = _0_0
    _lon0      = _0_0   # central meridian

#   see ._reset() below:

#   _e_PI_2     = _e *  PI_2
#   _e_PI_4     = _e *  PI_4
#   _e_TAYTOL   = _e * _TAYTOL

#   _1_e_90     = (_1_0 - _e) * _90_0
#   _1_e_PI_2   = (_1_0 - _e) *  PI_2
#   _1_e2_PI_2  = (_1_0 - _e * 2) * PI_2

#   _mu         =  _e**2
#   _mu_2_1     = (_e**2 + _2_0) * _0_5

#   _Eu         =  Elliptic(_mu)
#   _Eu_cE_1_4  = _Eu.cE * _0_25
#   _Eu_cK_cE   = _Eu.cK / _Eu.cE
#   _Eu_cK_PI_2 = _Eu.cK / PI_2

#   _mv         =  _1_0 - _mu
#   _3_mv       =  _3_0 / _mv
#   _3_mv_e     = _3_mv / _e

#   _Ev         =  Elliptic(_mv)
#   _Ev_cKE_3_4 = _Ev.cKE * 0.75
#   _Ev_cKE_5_4 = _Ev.cKE * 1.25

    def __init__(self, datum=_WGS84, lon0=0, k0=_K0_UTM, extendp=True, name=NN):
        '''New L{ExactTransverseMercator} projection.

           @kwarg datum: The datum, ellipsoid to use (L{Datum},
                         L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg lon0: The central meridian (C{degrees180}).
           @kwarg k0: The central scale factor (C{float}).
           @kwarg extendp: Use the extended domain (C{bool}).
           @kwarg name: Optional name for the projection (C{str}).

           @raise EllipticError: No convergence.

           @raise ETMError: Invalid B{C{k0}}.

           @raise TypeError: Invalid B{C{datum}}.

           @raise ValueError: Invalid B{C{lon0}} or B{C{k0}}.

           @note: The maximum error for all 255.5K U{TMcoords.dat
                  <https://Zenodo.org/record/32470>} tests (with
                  C{0 <= lat <= 84} and C{0 <= lon}) is C{5.2e-08
                  .forward} or 52 nano-meter easting and northing
                  and C{3.8e-13 .reverse} or 0.38 pico-degrees lat-
                  and longitude (with Python 3.7.3, 2.7.16, PyPy6
                  3.5.3 and PyPy6 2.7.13, all in 64-bit on macOS
                  10.13.6 High Sierra).
        '''
        if not extendp:
            self._extendp = False
        if name:
            self.name = name

        self.datum = datum
        self.lon0  = Lon(lon0=lon0)
        self.k0    = k0

    @property_doc_(''' the datum (L{Datum}).''')
    def datum(self):
        '''Get the datum (L{Datum}) or C{None}.
        '''
        return self._datum

    @datum.setter  # PYCHOK setter!
    def datum(self, datum):
        '''Set the datum and ellipsoid (L{Datum}, L{Ellipsoid},
           L{Ellipsoid2} or L{a_f2Tuple}).

           @raise EllipticError: No convergence.

           @raise TypeError: Invalid B{C{datum}}.
        '''
        d = _ellipsoidal_datum(datum, name=self.name, raiser=True)
        self._reset(d.ellipsoid)
        self._datum = d

    @Property_RO
    def equatoradius(self):
        '''Get the equatorial radius, semi-axis (C{meter}).
        '''
        return self._E.a

    @Property_RO
    def extendp(self):
        '''Get using the extended domain (C{bool}).
        '''
        return self._extendp

    @Property_RO
    def flattening(self):
        '''Get the flattening (C{float}).
        '''
        return self._E.f

    def forward(self, lat, lon, lon0=None):  # MCCABE 13
        '''Forward projection, from geographic to transverse Mercator.

           @arg lat: Latitude of point (C{degrees}).
           @arg lon: Longitude of point (C{degrees}).
           @kwarg lon0: Central meridian of the projection (C{degrees}).

           @return: L{EasNorExact4Tuple}C{(easting, northing,
                    convergence, scale)} in C{meter}, C{meter},
                    C{degrees} and C{scalar}.

           @see: C{void TMExact::Forward(real lon0, real lat, real lon,
                                         real &x, real &y,
                                         real &gamma, real &k)}.

           @raise EllipticError: No convergence.
        '''
        lat    = _fix90(lat)
        lon, _ = _diff182((self._lon0 if lon0 is None else lon0), lon)
        # Explicitly enforce the parity
        backside = _lat = _lon = False
        if not self.extendp:
            if lat < 0:
                _lat, lat = True, -lat
            if lon < 0:
                _lon, lon = True, -lon
            if lon > _90_0:
                backside = True
                if lat == 0:
                    _lat = True
                lon = _180_0 - lon

        # u,v = coordinates for the Thompson TM, Lee 54
        if lat == _90_0:
            u = self._Eu.cK
            v = self._iteration = 0
        elif lat == 0 and lon == self._1_e_90:
            u = self._iteration = 0
            v = self._Ev.cK
        else:  # tau = tan(phi), taup = sinh(psi)
            tau, lam = tan(radians(lat)), radians(lon)
            u, v = self._zetaInv(self._E.es_taupf(tau), lam)

        sncndn6 = self._sncndn6(u, v)
        xi, eta, _ = self._sigma3(v, *sncndn6)
        if backside:
            xi = 2 * self._Eu.cE - xi
        y = xi  * self._k0_a
        x = eta * self._k0_a

        if lat == _90_0:
            g, k = lon, self._k0
        else:
            g, k = self._zetaScaled(sncndn6, ll=False)

        if backside:
            g = _180_0 - g
        if _lat:
            y, g = neg_(y, g)
        if _lon:
            x, g = neg_(x, g)

        r = EasNorExact4Tuple(x, y, g, k)
        r._iteration = self._iteration
        return r

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
        self._k0 = Scalar_(k0=k0, Error=ETMError, low=_TOL_10, high=_1_0)
        # if not self._k0 > 0:
        #     raise Scalar_.Error_(Scalar_, k0, name=_k0_, Error=ETMError)
        self._k0_a = self._k0 * self.equatoradius  # see ._reset

    @property_doc_(''' the central meridian (C{degrees180}).''')
    def lon0(self):
        '''Get the central meridian (C{degrees180}).
        '''
        return self._lon0

    @lon0.setter  # PYCHOK setter!
    def lon0(self, lon0):
        '''Set the central meridian (C{degrees180}).

           @raise ValueError: Invalid B{C{lon0}}.
        '''
        self._lon0 = _norm180(Lon(lon0=lon0))

    @deprecated_Property_RO
    def majoradius(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{equatoradius}.'''
        return self.equatoradius

    def _reset(self, E):
        '''(INTERNAL) Get elliptic functions and pre-compute
           some frequently used values.

           @arg E: An ellipsoid (L{Ellipsoid}).

           @raise EllipticError: No convergence.
        '''
        # _update_all(self)
        e, e2 = E.e, E.e2
        # assert e2 == e**2 != _0_0

        self._e_PI_2     = e *  PI_2
        self._e_PI_4     = e *  PI_4
        self._e_TAYTOL   = e * _TAYTOL

        self._1_e_90     = (_1_0 - e) * _90_0
        self._1_e_PI_2   = (_1_0 - e) *  PI_2
        self._1_e2_PI_2  = (_1_0 - e * 2) * PI_2

        self._mu         =  e2
        self._mu_2_1     = (e2 + _2_0) * _0_5

        self._Eu         = Elliptic(self._mu)
        self._Eu_cE_1_4  = self._Eu.cE * _0_25
        self._Eu_cK_cE   = self._Eu.cK / self._Eu.cE
        self._Eu_cK_PI_2 = self._Eu.cK / PI_2

        self._mv         = _1_0 - e2
        self._3_mv       = _3_0 / self._mv
        self._3_mv_e     = self._3_mv / e

        self._Ev         = Elliptic(self._mv)
        self._Ev_cKE_3_4 = self._Ev.cKE * 0.75
        self._Ev_cKE_5_4 = self._Ev.cKE * 1.25

        self._iteration  = None

        self._E    = E
        self._k0_a = E.a * self.k0  # see .k0 setter

    def reverse(self, x, y, lon0=None):
        '''Reverse projection, from Transverse Mercator to geographic.

           @arg x: Easting of point (C{meters}).
           @arg y: Northing of point (C{meters}).
           @kwarg lon0: Central meridian of the projection (C{degrees}).

           @return: L{LatLonExact4Tuple}C{(lat, lon, convergence, scale)}
                    in C{degrees}, C{degrees180}, C{degrees} and C{scalar}.

           @see: C{void TMExact::Reverse(real lon0, real x, real y,
                                         real &lat, real &lon,
                                         real &gamma, real &k)}

           @raise EllipticError: No convergence.
        '''
        # undoes the steps in .forward.
        xi  = y / self._k0_a
        eta = x / self._k0_a
        backside = _lat = _lon = False
        if not self.extendp:  # enforce the parity
            if y < 0:
                _lat, xi  = True, -xi
            if x < 0:
                _lon, eta = True, -eta
            if xi > self._Eu.cE:
                xi = 2 * self._Eu.cE - xi
                backside = True

        # u,v = coordinates for the Thompson TM, Lee 54
        if xi != _0_0 or eta != self._Ev.cKE:
            u, v = self._sigmaInv(xi, eta)
        else:
            u = self._iteration = 0
            v = self._Ev.cK

        if v != _0_0 or u != self._Eu.cK:
            g, k, lat, lon = self._zetaScaled(self._sncndn6(u, v))
        else:
            g, k, lat, lon = _0_0, self._k0, _90_0, _0_0

        if backside:
            lon, g = (_180_0 - lon), (_180_0 - g)
        if _lat:
            lat, g = neg_(lat, g)
        if _lon:
            lon, g = neg_(lon, g)

        lon += self._lon0 if lon0 is None else _norm180(lon0)
        r = LatLonExact4Tuple(_norm180(lat), _norm180(lon), g, k)
        r._iteration = self._iteration
        return r

    def _scaled(self, tau, d2, snu, cnu, dnu, snv, cnv, dnv):
        '''(INTERNAL) C{scaled}.

           @note: Argument B{C{d2}} is C{_mu * cnu**2 + _mv * cnv**2}
                  from C{._sigma3} or C{._zeta3}.

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
        # the numerator from
        #
        #  (1 - snu^2 * dnv^2) to (_mv * snv^2 + cnu^2 * dnv^2)
        #
        # to maintain accuracy near phi = 90 and change the
        # denomintor from
        #  (dnu^2 + dnv^2 - 1) to (_mu * cnu^2 + _mv * cnv^2)
        #
        # to maintain accuracy near phi = 0, lam = 90 * (1 - e).
        # Similarly rewrite sqrt term in 9.1 as
        #
        #  _mv + _mu * c^2 instead of 1 - _mu * sin(phi)^2
        q2 = (mv * snv**2 + cnudnv**2) / d2
        # originally: sec2 = 1 + tau**2  # sec(phi)^2
        # k = sqrt(mv + mu / sec2) * sqrt(sec2) * sqrt(q2)
        #   = sqrt(mv + mv * tau**2 + mu) * sqrt(q2)
        k = sqrt(fsum_(mu, mv, mv * tau**2)) * sqrt(q2)
        return g, k * self._k0

    def _sigma3(self, v, snu, cnu, dnu, snv, cnv, dnv):  # PYCHOK unused
        '''(INTERNAL) C{sigma}.

           @return: 3-Tuple C{(xi, eta, d2)}.

           @see: C{void TMExact::sigma(real /*u*/, real snu, real cnu, real dnu,
                                       real   v,   real snv, real cnv, real dnv,
                                       real &xi, real &eta)}.

           @raise EllipticError: No convergence.
        '''
        # Lee 55.4 writing
        # dnu^2 + dnv^2 - 1 = _mu * cnu^2 + _mv * cnv^2
        d2 = self._mu * cnu**2 + self._mv * cnv**2
        xi =      self._Eu.fE(snu, cnu, dnu) - self._mu * snu * cnu * dnu / d2
        eta = v - self._Ev.fE(snv, cnv, dnv) + self._mv * snv * cnv * dnv / d2
        return xi, eta, d2

    def _sigmaDwd(self, snu, cnu, dnu, snv, cnv, dnv):
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
        dv = neg(_2_0 * i * r / d)
        return du, dv

    def _sigmaInv(self, xi, eta):
        '''(INTERNAL) Invert C{sigma} using Newton's method.

           @return: 2-Tuple C{(u, v)}.

           @see: C{void TMExact::sigmainv(real xi, real eta,
                                          real &u, real &v)}.

           @raise EllipticError: No convergence.
        '''
        u, v, trip = self._sigmaInv0(xi, eta)
        if trip:
            self._iteration = 0
        else:
            U, V = Fsum(u), Fsum(v)
            # min iterations = 2, max = 7, mean = 3.9
            for self._iteration in range(1, _TRIPS):  # GEOGRAPHICLIB_PANIC
                sncndn6 = self._sncndn6(u, v)
                X, E, _ = self._sigma3(v, *sncndn6)
                dw, dv  = self._sigmaDwd( *sncndn6)
                X  = xi - X
                E -= eta
                u, du = U.fsum2_(X * dw,  E * dv)
                v, dv = V.fsum2_(X * dv, -E * dw)
                if trip:
                    break
                trip = hypot2(du, dv) < _TOL_10
            else:
                t = unstr(self._sigmaInv.__name__, xi, eta)
                raise EllipticError(_no_(_convergence_), txt=t)
        return u, v

    def _sigmaInv0(self, xi, eta):
        '''(INTERNAL) Starting point for C{sigmaInv}.

           @return: 3-Tuple C{(u, v, trip)}.

           @see: C{bool TMExact::sigmainv0(real xi, real eta,
                                           real &u, real &v)}.
        '''
        trip = False
        if eta > self._Ev_cKE_5_4 or xi < min(- self._Eu_cE_1_4,
                                          eta - self._Ev.cKE):
            # sigma as a simple pole at
            #  w = w0 = Eu.K() + i * Ev.K()
            # and sigma is approximated by
            #  sigma = (Eu.E() + i * Ev.KE()) + 1 / (w - w0)
            x = xi  - self._Eu.cE
            y = eta - self._Ev.cKE
            d = hypot2(x, y)
            u = self._Eu.cK + x / d
            v = self._Ev.cK - y / d

        elif eta > self._Ev.cKE or (xi < self._Eu_cE_1_4 and
                                   eta > self._Ev_cKE_3_4):
            # At w = w0 = i * Ev.K(), we have
            #  sigma  = sigma0  = i * Ev.KE()
            #  sigma' = sigma'' = 0
            #
            # including the next term in the Taylor series gives:
            #  sigma = sigma0 - _mv / 3 * (w - w0)^3
            #
            # When inverting this, we map arg(w - w0) = [-pi/2, -pi/6]
            # to arg(sigma - sigma0) = [-pi/2, pi/2]
            # mapping arg = [-pi/2, -pi/6] to [-pi/2, pi/2]
            d = eta - self._Ev.cKE
            r = hypot(xi, d)
            # Error using this guess is about 0.068 * rad^(5/3)
            trip = r < _TAYTOL2
            # Map the range [-90, 180] in sigma space to [-90, 0] in
            # w space.  See discussion in zetainv0 on the cut for ang.
            r = cbrt(r * self._3_mv)
            a = atan2(d - xi, xi + d) / _3_0 - PI_4
            s, c = sincos2(a)
            u = r * c
            v = r * s + self._Ev.cK

        else:  # use w = sigma * Eu.K/Eu.E (correct in the limit _e -> 0)
            u = xi  * self._Eu_cK_cE
            v = eta * self._Eu_cK_cE

        return u, v, trip

    def _sncndn6(self, u, v):
        '''(INTERNAL) Get 6-tuple C{(snu, cnu, dnu, snv, cnv, dnv)}.
        '''
        # snu, cnu, dnu = self._Eu.sncndn(u)
        # snv, cnv, dnv = self._Ev.sncndn(v)
        return self._Eu.sncndn(u) + self._Ev.sncndn(v)

    def toStr(self, **kwds):
        '''Return a C{str} representation.

           @arg kwds: Optional, overriding keyword arguments.
        '''
        d = dict(name=self.name) if self.name else {}
        d = dict(datum=self.datum.name, lon0=self.lon0,
                 k0=self.k0, extendp=self.extendp, **d)
        return _COMMASPACE_.join(pairs(d, **kwds))

    def _zeta3(self, snu, cnu, dnu, snv, cnv, dnv):
        '''(INTERNAL) C{zeta}.

           @return: 3-Tuple C{(taup, lambda, d2)}.

           @see: C{void TMExact::zeta(real /*u*/, real snu, real cnu, real dnu,
                                      real /*v*/, real snv, real cnv, real dnv,
                                      real &taup, real &lam)}
        '''
        e = self._E.e
        # Lee 54.17 but write
        # atanh(snu * dnv)      = asinh(snu * dnv / sqrt(cnu^2 + _mv * snu^2 * snv^2))
        # atanh(_e * snu / dnv) = asinh(_e * snu / sqrt(_mu * cnu^2 + _mv * cnv^2))
        d1 = cnu**2 + self._mv * (snu * snv)**2
        d2 = self._mu * cnu**2 + self._mv * cnv**2
        # Overflow value s.t. atan(overflow) = pi/2
        t1 = t2 = copysign(_OVERFLOW, snu)
        if d1 > _EPS0__2:
            t1 = snu * dnv / sqrt(d1)
        lam = _0_0
        if d2 > _EPS0__2:
            t2 = sinh(e * asinh(e * snu / sqrt(d2)))
            if d1 > _EPS0__2:
                lam = atan2(dnu * snv    , cnu * cnv) - \
                      atan2(cnu * snv * e, dnu * cnv) * e
        # psi = asinh(t1) - asinh(t2)
        # taup = sinh(psi)
        taup = t1 * hypot1(t2) - t2 * hypot1(t1)
        return taup, lam, d2

    def _zetaDwd(self, snu, cnu, dnu, snv, cnv, dnv):
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
        # Lee 54.21 but write
        # (1 - dnu^2 * snv^2) = (cnv^2 + _mu * snu^2 * snv^2)
        # (see A+S 16.21.4)
        d  = self._mv   * (cnv2 + snuv2)**2
        du = cnu * dnuv * (cnv2 - snuv2) / d
        dv = cnv * snuv * (cnu2 + dnuv2) / d
        return du, neg(dv)

    def _zetaInv(self, taup, lam):
        '''(INTERNAL) Invert C{zeta} using Newton's method.

           @return: 2-Tuple C{(u, v)}.

           @see: C{void TMExact::zetainv(real taup, real lam,
                                         real &u, real &v)}.

           @raise EllipticError: No convergence.
        '''
        psi = asinh(taup)
        u, v, trip = self._zetaInv0(psi, lam)
        if trip:
            self._iteration = 0
        else:
            sca   = _1_0 / hypot1(taup)
            stol2 = _TOL_10 / max(psi**2, _1_0)
            U, V  =  Fsum(u), Fsum(v)
            # min iterations = 2, max = 6, mean = 4.0
            for self._iteration in range(1, _TRIPS):  # GEOGRAPHICLIB_PANIC
                sncndn6 = self._sncndn6(u, v)
                T, L, _ = self._zeta3(  *sncndn6)
                dw, dv  = self._zetaDwd(*sncndn6)
                T  = (taup - T) * sca
                L -= lam
                u, du = U.fsum2_(T * dw,  L * dv)
                v, dv = V.fsum2_(T * dv, -L * dw)
                if trip:
                    break
                trip = hypot2(du, dv) < stol2
            else:
                t = unstr(self._zetaInv.__name__, taup, lam)
                raise EllipticError(_no_(_convergence_), txt=t)
        return u, v

    def _zetaInv0(self, psi, lam):
        '''(INTERNAL) Starting point for C{zetaInv}.

           @return: 3-Tuple C{(u, v, trip)}.

           @see: C{bool TMExact::zetainv0(real psi, real lam,  # radians
                                          real &u, real &v)}.
        '''
        trip = False
        if psi < -self._e_PI_4 and lam >        self._1_e2_PI_2 \
                               and psi < (lam - self._1_e_PI_2):
            # N.B. this branch is normally not taken because psi < 0
            # is converted psi > 0 by Forward.
            #
            # There's a log singularity at w = w0 = Eu.K() + i * Ev.K(),
            # corresponding to the south pole, where we have, approximately
            #
            #  psi = _e + i * pi/2 - _e * atanh(cos(i * (w - w0)/(1 + _mu/2)))
            #
            # Inverting this gives:
            e = self._E.e
            h = sinh(_1_0 - psi / e)
            a = (PI_2 - lam) / e
            s, c = sincos2(a)
            u = self._Eu.cK - asinh(s / hypot(c, h)) * self._mu_2_1
            v = self._Ev.cK - atan2(c, h) * self._mu_2_1

        elif psi < self._e_PI_2 and lam > self._1_e2_PI_2:
            # At w = w0 = i * Ev.K(), we have
            #
            #  zeta  = zeta0  = i * (1 - _e) * pi/2
            #  zeta' = zeta'' = 0
            #
            # including the next term in the Taylor series gives:
            #
            #  zeta = zeta0 - (_mv * _e) / 3 * (w - w0)^3
            #
            # When inverting this, we map arg(w - w0) = [-90, 0] to
            # arg(zeta - zeta0) = [-90, 180]

            d = lam - self._1_e_PI_2
            r = hypot(psi, d)
            # Error using this guess is about 0.21 * (rad/e)^(5/3)
            trip = r < self._e_TAYTOL
            # atan2(dlam-psi, psi+dlam) + 45d gives arg(zeta - zeta0)
            # in range [-135, 225).  Subtracting 180 (since multiplier
            # is negative) makes range [-315, 45).  Multiplying by 1/3
            # (for cube root) gives range [-105, 15).  In particular
            # the range [-90, 180] in zeta space maps to [-90, 0] in
            # w space as required.
            r = cbrt(r * self._3_mv_e)
            a = atan2(d - psi, psi + d) / _3_0 - PI_4
            s, c = sincos2(a)
            u = r * c
            v = r * s + self._Ev.cK

        else:
            # Use spherical TM, Lee 12.6 -- writing C{atanh(sin(lam) /
            # cosh(psi)) = asinh(sin(lam) / hypot(cos(lam), sinh(psi)))}.
            # This takes care of the log singularity at C{zeta = Eu.K()},
            # corresponding to the north pole.
            s, c = sincos2(lam)
            h, r = sinh(psi), self._Eu_cK_PI_2
            # But scale to put 90, 0 on the right place
            u = r * atan2(h, c)
            v = r * asinh(s / hypot(c, h))

        return u, v, trip

    def _zetaScaled(self, sncndn6, ll=True):
        '''(INTERNAL) Recompute (T, L) from (u, v) to improve accuracy of Scale.

           @arg sncndn6: 6-Tuple C{(snu, cnu, dnu, snv, cnv, dnv)}.

           @return: 2-Tuple C{(g, k)} if B{C{ll}} is C{False} else
                    4-tuple C{(g, k, lat, lon)}.
        '''
        t, lam, d2 = self._zeta3( *sncndn6)
        tau = self._E.es_tauf(t)
        r = self._scaled(tau, d2, *sncndn6)
        if ll:
            r += atand(tau), degrees(lam)
        return r


class LatLonExact4Tuple(_NamedTuple):
    '''4-Tuple C{(lat, lon, convergence, scale)} in C{degrees180},
       C{degrees180}, C{degrees} and C{scalar}.
    '''
    _Names_ = (_lat_, _lon_, _convergence_, _scale_)
    _Units_ = ( Lat,   Lon,   Degrees,       Scalar)


def parseETM5(strUTM, datum=_WGS84, Etm=Etm, falsed=True, name=NN):
    '''Parse a string representing a UTM coordinate, consisting
       of C{"zone[band] hemisphere easting northing"}.

       @arg strUTM: A UTM coordinate (C{str}).
       @kwarg datum: Optional datum to use (L{Datum}, L{Ellipsoid},
                     L{Ellipsoid2} or L{a_f2Tuple}).
       @kwarg Etm: Optional class to return the UTM coordinate
                   (L{Etm}) or C{None}.
       @kwarg falsed: Both easting and northing are falsed (C{bool}).
       @kwarg name: Optional B{C{Etm}} name (C{str}).

       @return: The UTM coordinate (B{C{Etm}}) or if B{C{Etm}} is
                C{None}, a L{UtmUps5Tuple}C{(zone, hemipole, easting,
                northing, band)}.  The C{hemipole} is the hemisphere
                C{'N'|'S'}.

       @raise ETMError: Invalid B{C{strUTM}}.

       @raise TypeError: Invalid B{C{datum}}.

       @example:

        >>> u = parseETM5('31 N 448251 5411932')
        >>> u.toRepr()  # [Z:31, H:N, E:448251, N:5411932]
        >>> u = parseETM5('31 N 448251.8 5411932.7')
        >>> u.toStr()  # 31 N 448252 5411933
    '''
    r = _parseUTM5(strUTM, datum, Etm, falsed, Error=ETMError, name=name)
    return r


def toEtm8(latlon, lon=None, datum=None, Etm=Etm, falsed=True, name=NN,
                                         zone=None, **cmoff):
    '''Convert a lat-/longitude point to an ETM coordinate.

       @arg latlon: Latitude (C{degrees}) or an (ellipsoidal)
                    geodetic C{LatLon} point.
       @kwarg lon: Optional longitude (C{degrees}) or C{None}.
       @kwarg datum: Optional datum for this ETM coordinate,
                     overriding B{C{latlon}}'s datum (L{Datum},
                     L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}).
       @kwarg Etm: Optional class to return the ETM coordinate
                   (L{Etm}) or C{None}.
       @kwarg falsed: False both easting and northing (C{bool}).
       @kwarg name: Optional B{C{Utm}} name (C{str}).
       @kwarg zone: Optional UTM zone to enforce (C{int} or C{str}).
       @kwarg cmoff: DEPRECATED, use B{C{falsed}}.  Offset longitude
                     from the zone's central meridian (C{bool}).

       @return: The ETM coordinate (B{C{Etm}}) or a
                L{UtmUps8Tuple}C{(zone, hemipole, easting, northing,
                band, datum, convergence, scale)} if B{C{Etm}} is
                C{None} or not B{C{falsed}}.  The C{hemipole} is the
                C{'N'|'S'} hemisphere.

       @raise EllipticError: No convergence.

       @raise ETMError: Invalid B{C{zone}}.

       @raise RangeError: If B{C{lat}} outside the valid UTM bands or
                          if B{C{lat}} or B{C{lon}} outside the valid
                          range and L{rangerrors} set to C{True}.

       @raise TypeError: Invalid B{C{datum}} or B{C{latlon}} not
                         ellipsoidal.

       @raise ValueError: If B{C{lon}} value is missing or if
                          B{C{latlon}} is invalid.
    '''
    z, B, lat, lon, d, f, name = _to7zBlldfn(latlon, lon, datum,
                                             falsed, name, zone,
                                             ETMError, **cmoff)
    lon0 = _cmlon(z) if f else None
    x, y, g, k = d.exactTM.forward(lat, lon, lon0=lon0)

    return _toXtm8(Etm, z, lat, x, y, B, d, g, k, f,
                        name, latlon, d.exactTM, Error=ETMError)

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
