
# -*- coding: utf-8 -*-

u'''Classes L{ETMError} and L{Etm}, a pure Python implementation of
I{Charles Karney}'s C++ class U{TransverseMercatorExact
<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1TransverseMercatorExact.html>},
abbreviated as C{TMExact} below.

Python class L{ExactTransverseMercator} implements the C{Exact Transverse
Mercator} (ETM) projection.  Instances of class L{Etm} represent ETM
C{easting, nothing} locations.

Following is a copy of Karney's U{TransverseMercatorExact.hpp
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
'''
# make sure int/int division yields float quotient
from __future__ import division
division = 1 / 2  # double check int division, see .datum.py, .utily.py
if not division:
    raise ImportError('%s 1/2 == %d' % ('division', division))
del division

from pygeodesy.basics import EPS, PI_2, PI_4, property_doc_, \
                             property_RO, _xinstanceof
from pygeodesy.datum import Datum, Datums
from pygeodesy.elliptic import Elliptic, EllipticError, _TRIPS
from pygeodesy.errors import _incompatible
from pygeodesy.fmath import cbrt, Fsum, fsum_, hypot, hypot1, hypot2
from pygeodesy.interns import _COMMA_SPACE_, _convergence_, _easting_, \
                              _k0_, _lat_, _lon_, _lon0_, NN, _northing_, \
                              _no_convergence_, _scale_  # PYCHOK used!
from pygeodesy.karney import _diff182, _fix90, _norm180
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import _NamedBase, _NamedTuple, _xnamed
from pygeodesy.streprs import pairs, unstr
from pygeodesy.units import Lon, Scalar_
from pygeodesy.utily import sincos2
from pygeodesy.utm import _cmlon, _K0, _parseUTM5, Utm, UTMError, \
                          _toXtm8, _to7zBlldfn
from pygeodesy.utmupsBase import _LLEB

from math import asinh, atan, atan2, copysign, degrees, radians, \
                 sinh, sqrt, tan

__all__ = _ALL_LAZY.etm
__version__ = '20.08.04'

_OVERFLOW = 1.0 / EPS**2
_TOL      = EPS
_TOL_10   = 0.1 * _TOL
_TAYTOL   = pow(_TOL, 0.6)
_TAYTOL2  = 2.0 * _TAYTOL


class EasNorExact4Tuple(_NamedTuple):
    '''4-Tuple C{(easting, northing, convergence, scale)} in
       C{meter}, C{meter}, C{degrees} and C{scalar}.
    '''
    _Names_ = (_easting_, _northing_, _convergence_, _scale_)


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
                             datum=Datums.WGS84, falsed=True,
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
           @kwarg datum: Optional, this coordinate's datum (L{Datum}).
           @kwarg falsed: Both B{C{easting}} and B{C{northing}} are
                          falsed (C{bool}).
           @kwarg convergence: Optional meridian convergence, bearing
                               off grid North, clockwise from true North
                               (C{degrees}) or C{None}.
           @kwarg scale: Optional grid scale factor (C{scalar}) or C{None}.
           @kwarg name: Optional name (C{str}).

           @raise ETMError: Invalid B{C{zone}}, B{C{hemishere}} or
                            B{C{band}}.

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
        if exactTM._E != E or exactTM.majoradius != E.a \
                           or exactTM.flattening != E.f:
            raise ETMError(repr(exactTM), txt=_incompatible(repr(E)))
        self._exactTM = exactTM
        self._scale0  = exactTM.k0

    def parseETM(self, strETM):
        '''Parse a string to a ETM coordinate.

           @return: The coordinate (L{Etm}).

           @see: Function L{parseETM5} in this module L{etm}.
        '''
        return parseETM5(strETM, datum=self.datum, Etm=self.classof)

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
        e, n = self.to2en(falsed=f)
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
        '''Coopy this ETM to a UTM coordinate.

           @return: The UTM coordinate (L{Utm}).
        '''
        return self._xnamed(self._xcopy2(Utm))


class ExactTransverseMercator(_NamedBase):
    '''A Python version of Karney's U{TransverseMercatorExact
       <https://GeographicLib.SourceForge.io/html/TransverseMercatorExact_8cpp_source.html>}
       C++ class, a numerically exact transverse mercator projection,
       referred to as C{TMExact} here.

       @see: C{U{TMExact(real a, real f, real k0, bool extendp)<https://geographiclib.sourceforge.io/
             html/classGeographicLib_1_1TransverseMercatorExact.html#a72ffcc89eee6f30a6d1f4d061518a6f1>}}.
    '''
    _a       = 0     # major radius
    _datum   = None  # Datum
    _e       = 0     # eccentricity
    _E       = None  # Ellipsoid
    _extendp = True
    _f       = 0     # flattening
    _k0      = 1     # central scale factor
    _k0_a    = 0
    _lon0    = 0     # central meridian
    _trips_  = _TRIPS

#   see ._reset() below:

#   _e_PI_2     = _e *  PI_2
#   _e_PI_4     = _e *  PI_4
#   _e_TAYTOL   = _e * _TAYTOL

#   _1_e_90     = (1 - _e) * 90
#   _1_e_PI_2   = (1 - _e) * PI_2
#   _1_e2_PI_2  = (1 - _e * 2) * PI_2

#   _mu         =  _e**2
#   _mu_2_1     = (_e**2 + 2) * 0.5

#   _Eu         =  Elliptic(_mu)
#   _Eu_cE_1_4  = _Eu.cE * 0.25
#   _Eu_cK_cE   = _Eu.cK / _Eu.cE
#   _Eu_cK_PI_2 = _Eu.cK / PI_2

#   _mv         =  1   - _mu
#   _3_mv       =  3.0 / _mv
#   _3_mv_e     = _3_mv / _e

#   _Ev         =  Elliptic(_mv)
#   _Ev_cKE_3_4 = _Ev.cKE * 0.75
#   _Ev_cKE_5_4 = _Ev.cKE * 1.25

    def __init__(self, datum=Datums.WGS84, lon0=0, k0=_K0, extendp=True, name=NN):
        '''New L{ExactTransverseMercator} projection.

           @kwarg datum: The datum and ellipsoid to use (C{Datum}).
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
        self.lon0  = lon0
        self.k0    = k0

    @property_doc_(''' the datum (L{Datum}).''')
    def datum(self):
        '''Get the datum (L{Datum}) or C{None}.
        '''
        return self._datum

    @datum.setter  # PYCHOK setter!
    def datum(self, datum):
        '''Set the datum and ellipsoid (L{Datum}).

           @raise EllipticError: No convergence.

           @raise TypeError: Invalid B{C{datum}}.
        '''
        _xinstanceof(Datum, datum=datum)

        E = datum.ellipsoid
        self._reset(E.e, E.e2)
        self._a = E.a
        self._f = E.f  # flattening = (a - b) / a

        self._datum = datum
        self._E     = E

    @property_RO
    def equatoradius(self):
        '''Get the equatorial (major) radius, semi-axis (C{meter}).
        '''
        return self._a

    majoradius = equatoradius  # for backward compatibility

    @property_RO
    def extendp(self):
        '''Get using the extended domain (C{bool}).
        '''
        return self._extendp

    @property_RO
    def flattening(self):
        '''Get the flattening (C{float}).
        '''
        return self._f

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
        lat = _fix90(lat)
        lon, _ = _diff182((self._lon0 if lon0 is None else lon0), lon)
        # Explicitly enforce the parity
        backside = _lat = _lon = False
        if not self.extendp:
            if lat < 0:
                _lat, lat = True, -lat
            if lon < 0:
                _lon, lon = True, -lon
            if lon > 90:
                backside = True
                if lat == 0:
                    _lat = True
                lon = 180 - lon

        # u,v = coordinates for the Thompson TM, Lee 54
        if lat == 90:
            u, v = self._Eu.cK, 0
        elif lat == 0 and lon == self._1_e_90:
            u, v = 0, self._Ev.cK
        else:  # tau = tan(phi), taup = sinh(psi)
            tau, lam = tan(radians(lat)), radians(lon)
            u, v = self._zetaInv(self._E.es_taupf(tau), lam)

        sncndn6 = self._sncndn6(u, v)
        xi, eta, _ = self._sigma3(v, *sncndn6)
        if backside:
            xi = 2 * self._Eu.cE - xi
        y = xi  * self._k0_a
        x = eta * self._k0_a

        if lat == 90:
            g, k = lon, self._k0
        else:
            g, k = self._zetaScaled(sncndn6, ll=False)

        if backside:
            g = 180 - g
        if _lat:
            y, g = -y, -g
        if _lon:
            x, g = -x, -g
        return EasNorExact4Tuple(x, y, g, k)

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
        self._k0 = Scalar_(k0, name=_k0_, Error=ETMError, low=_TOL_10, high=1.0)
        # if not self._k0 > 0:
        #     raise Scalar_.Error_(Scalar_, k0, name=_k0_, Error=ETMError)
        self._k0_a = self._k0 * self._a

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
        self._lon0 = _norm180(Lon(lon0, name=_lon0_))

    def _reset(self, e, e2):
        '''(INTERNAL) Get elliptic functions and pre-compute some
           frequently used values.

           @arg e: Eccentricity (C{float}).
           @arg e2: Eccentricity squared (C{float}).

           @raise EllipticError: No convergence.
        '''
        # assert e2 == e**2
        self._e          = e
        self._e_PI_2     = e *  PI_2
        self._e_PI_4     = e *  PI_4
        self._e_TAYTOL   = e * _TAYTOL

        self._1_e_90     = (1 - e) * 90
        self._1_e_PI_2   = (1 - e) * PI_2
        self._1_e2_PI_2  = (1 - e * 2) * PI_2

        self._mu         =  e2
        self._mu_2_1     = (e2 + 2) * 0.5

        self._Eu         = Elliptic(self._mu)
        self._Eu_cE_1_4  = self._Eu.cE * 0.25
        self._Eu_cK_cE   = self._Eu.cK / self._Eu.cE
        self._Eu_cK_PI_2 = self._Eu.cK / PI_2

        self._mv         = 1 - e2
        self._3_mv       = 3.0 / self._mv
        self._3_mv_e     = self._3_mv / e

        self._Ev         = Elliptic(self._mv)
        self._Ev_cKE_3_4 = self._Ev.cKE * 0.75
        self._Ev_cKE_5_4 = self._Ev.cKE * 1.25

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
        if xi != 0 or eta != self._Ev.cKE:
            u, v = self._sigmaInv(xi, eta)
        else:
            u, v = 0, self._Ev.cK

        if v != 0 or u != self._Eu.cK:
            g, k, lat, lon = self._zetaScaled(self._sncndn6(u, v))
        else:
            g, k, lat, lon = 0, self._k0, 90, 0

        if backside:
            lon, g = (180 - lon), (180 - g)
        if _lat:
            lat, g = -lat, -g
        if _lon:
            lon, g = -lon, -g

        lon += self._lon0 if lon0 is None else _norm180(lon0)
        return LatLonExact4Tuple(_norm180(lat), _norm180(lon), g, k)

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
        g = atan2(mv * cnv * snv * snu, cnudnv * dnu)
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
        return degrees(g), k * self._k0

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
        r =  cnv * dnu * dnv
        i = -cnu * snuv * self._mu
        du = (r**2 - i**2) / d
        dv = 2 * r * i / d
        return du, dv

    def _sigmaInv(self, xi, eta):
        '''(INTERNAL) Invert C{sigma} using Newton's method.

           @return: 2-Tuple C{(u, v)}.

           @see: C{void TMExact::sigmainv(real xi, real eta,
                                          real &u, real &v)}.

           @raise EllipticError: No convergence.
        '''
        u, v, trip = self._sigmaInv0(xi, eta)
        if not trip:
            U, V = Fsum(u), Fsum(v)
            # min iterations = 2, max = 7, mean = 3.9
            for _ in range(self._trips_):  # GEOGRAPHICLIB_PANIC
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
                raise EllipticError(_no_convergence_, txt=t)
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
            a = atan2(d - xi, xi + d) / 3.0 - PI_4
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
        return _COMMA_SPACE_.join(pairs(d, **kwds))

    def _zeta3(self, snu, cnu, dnu, snv, cnv, dnv):
        '''(INTERNAL) C{zeta}.

           @return: 3-Tuple C{(taup, lambda, d2)}.

           @see: C{void TMExact::zeta(real /*u*/, real snu, real cnu, real dnu,
                                      real /*v*/, real snv, real cnv, real dnv,
                                      real &taup, real &lam)}
        '''
        e = self._e
        # Lee 54.17 but write
        # atanh(snu * dnv)      = asinh(snu * dnv / sqrt(cnu^2 + _mv * snu^2 * snv^2))
        # atanh(_e * snu / dnv) = asinh(_e * snu / sqrt(_mu * cnu^2 + _mv * cnv^2))
        d1 = cnu**2 + self._mv * (snu * snv)**2
        d2 = self._mu * cnu**2 + self._mv * cnv**2
        # Overflow value s.t. atan(overflow) = pi/2
        t1 = t2 = copysign(_OVERFLOW, snu)
        if d1 > 0:
            t1 = snu * dnv / sqrt(d1)
        lam = 0
        if d2 > 0:
            t2 = sinh(e * asinh(e * snu / sqrt(d2)))
            if d1 > 0:
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
        d  =  self._mv   * (cnv2 + snuv2)**2
        du =  cnu * dnuv * (cnv2 - snuv2) / d
        dv = -cnv * snuv * (cnu2 + dnuv2) / d
        return du, dv

    def _zetaInv(self, taup, lam):
        '''(INTERNAL) Invert C{zeta} using Newton's method.

           @return: 2-Tuple C{(u, v)}.

           @see: C{void TMExact::zetainv(real taup, real lam,
                                         real &u, real &v)}.

           @raise EllipticError: No convergence.
        '''
        psi = asinh(taup)
        sca = 1.0 / hypot1(taup)
        u, v, trip = self._zetaInv0(psi, lam)
        if not trip:
            stol2 = _TOL_10 / max(psi**2, 1.0)
            U, V = Fsum(u), Fsum(v)
            # min iterations = 2, max = 6, mean = 4.0
            for _ in range(self._trips_):  # GEOGRAPHICLIB_PANIC
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
                raise EllipticError(_no_convergence_, txt=t)
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
            h = sinh(1 - psi / self._e)
            a = (PI_2 - lam) / self._e
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
            a = atan2(d - psi, psi + d) / 3.0 - PI_4
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
            r += degrees(atan(tau)), degrees(lam)
        return r


class LatLonExact4Tuple(_NamedTuple):
    '''4-Tuple C{(lat, lon, convergence, scale)} in C{degrees180},
       C{degrees}, C{degrees} and C{scalar}.
    '''
    _Names_ = (_lat_, _lon_, _convergence_, _scale_)


def parseETM5(strUTM, datum=Datums.WGS84, Etm=Etm, falsed=True, name=NN):
    '''Parse a string representing a UTM coordinate, consisting
       of C{"zone[band] hemisphere easting northing"}.

       @arg strUTM: A UTM coordinate (C{str}).
       @kwarg datum: Optional datum to use (L{Datum}).
       @kwarg Etm: Optional class to return the UTM coordinate
                   (L{Etm}) or C{None}.
       @kwarg falsed: Both easting and northing are falsed (C{bool}).
       @kwarg name: Optional B{C{Etm}} name (C{str}).

       @return: The UTM coordinate (B{C{Etm}}) or if B{C{Etm}} is
                C{None}, a L{UtmUps5Tuple}C{(zone, hemipole, easting,
                northing, band)}.  The C{hemipole} is the hemisphere
                C{'N'|'S'}.

       @raise ETMError: Invalid B{C{strUTM}}.

       @example:

       >>> u = parseETM5('31 N 448251 5411932')
       >>> u.toStr2()  # [Z:31, H:N, E:448251, N:5411932]
       >>> u = parseETM5('31 N 448251.8 5411932.7')
       >>> u.toStr()  # 31 N 448252 5411933
    '''
    r = _parseUTM5(strUTM, datum, Etm, falsed, Error=ETMError)
    return _xnamed(r, name)


def toEtm8(latlon, lon=None, datum=None, Etm=Etm, falsed=True, name=NN,
                                         zone=None, **cmoff):
    '''Convert a lat-/longitude point to an ETM coordinate.

       @arg latlon: Latitude (C{degrees}) or an (ellipsoidal)
                    geodetic C{LatLon} point.
       @kwarg lon: Optional longitude (C{degrees}) or C{None}.
       @kwarg datum: Optional datum for this ETM coordinate,
                     overriding B{C{latlon}}'s datum (C{Datum}).
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

       @raise TypeError: If B{C{latlon}} is not ellipsoidal.

       @raise RangeError: If B{C{lat}} outside the valid UTM bands or
                          if B{C{lat}} or B{C{lon}} outside the valid
                          range and L{rangerrors} set to C{True}.

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


__all__ += _ALL_DOCS(EasNorExact4Tuple, LatLonExact4Tuple)

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
