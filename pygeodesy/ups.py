
# -*- coding: utf-8 -*-

u'''I{Karney}'s Universal Polar Stereographic (UPS) projection.

Classes L{Ups} and L{UPSError} and functions L{parseUPS5}, L{toUps8} and L{upsZoneBand5}.

A pure Python implementation, partially transcoded from C++ class U{PolarStereographic
<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1PolarStereographic.html>}
by I{Charles Karney}.

The U{UPS<https://WikiPedia.org/wiki/Universal_polar_stereographic_coordinate_system>}
system is used in conjuction with U{UTM
<https://WikiPedia.org/wiki/Universal_Transverse_Mercator_coordinate_system>}
for locations on the polar regions of the earth.  UPS covers areas south of 79.5°S
and north of 83.5°N (slightly overlapping the UTM range from 80°S to 84°N by 30' at
each end).
'''

# from pygeodesy.basics import neg  # from .dms
from pygeodesy.datums import _ellipsoidal_datum, _WGS84
from pygeodesy.dms import degDMS, neg, parseDMS2
from pygeodesy.errors import RangeError, _ValueError
from pygeodesy.fmath import hypot, hypot1, sqrt0
from pygeodesy.interns import EPS, EPS0, NN, _COMMASPACE_, _EPSmin as _Tol90, \
                             _N_, _inside_, _pole_, _range_, _S_, _SPACE_, \
                             _to_, _UTM_, _0_0, _0_5, _1_0, _2_0, _90_0
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import nameof, _xnamed
from pygeodesy.namedTuples import EasNor2Tuple, UtmUps5Tuple, \
                                  UtmUps8Tuple, UtmUpsLatLon5Tuple
from pygeodesy.props import deprecated_method, property_doc_, \
                            Property_RO, _update_all
# from pygeodesy.streprs import Fmt  # from .utmupsBase
from pygeodesy.units import Meter, Lat, Scalar, Scalar_
from pygeodesy.utily import degrees90, degrees180, sincos2d
from pygeodesy.utmupsBase import Fmt, _LLEB, _hemi, _parseUTMUPS5, _to4lldn, \
                                _to3zBhp, _to3zll, _UPS_BANDS as _Bands, \
                                _UPS_LAT_MAX, _UPS_LAT_MIN, _UPS_ZONE, \
                                _UPS_ZONE_STR, UtmUpsBase

from math import atan, atan2, radians, tan

__all__ = _ALL_LAZY.ups
__version__ = '21.11.11'

_Falsing = Meter(2000e3)  # false easting and northing (C{meter})
_K0_UPS  = Scalar(0.994)  # central UPS scale factor
_K1_UPS  = Scalar(_1_0)   # rescale point scale factor


def _scale(E, rho, tau):
    # compute the point scale factor, ala Karney
    t = hypot1(tau)
    return Scalar((rho / E.a) * t * sqrt0(E.e12 + E.e2 / t**2))


def _toBand(a_lat, b_lon):  # see utm._toBand
    '''(INTERNAL) Get the I{polar} Band letter from
       the signs of (phi, lam) or (lat, lon).
    '''
    return _Bands[(0 if a_lat < 0 else 2) + (0 if b_lon < 0 else 1)]


class UPSError(_ValueError):
    '''Universal Polar Stereographic (UPS) parse or other L{Ups} issue.
    '''
    pass


class Ups(UtmUpsBase):
    '''Universal Polar Stereographic (UPS) coordinate.
    '''
#   _band        =  NN        # polar band ('A', 'B', 'Y' or 'Z')
    _Bands       = _Bands     # polar Band letters (tuple)
    _Error       =  UPSError  # Error class
    _pole        =  NN        # UPS projection top/center ('N' or 'S')
#   _scale       =  None      # point scale factor (C{scalar})
    _scale0      = _K0_UPS    # central scale factor (C{scalar})

    def __init__(self, zone, pole, easting, northing, band=NN,  # PYCHOK expected
                                   datum=_WGS84, falsed=True,
                                   convergence=None, scale=None, name=NN):
        '''New L{Ups} UPS coordinate.

           @arg zone: UPS zone (C{int}, zero) or zone with/-out I{polar} Band
                      letter (C{str}, '00', '00A', '00B', '00Y' or '00Z').
           @arg pole: Top/center of (stereographic) projection (C{str},
                      C{'N[orth]'} or C{'S[outh]'}).
           @arg easting: Easting, see B{C{falsed}} (C{meter}).
           @arg northing: Northing, see B{C{falsed}} (C{meter}).
           @kwarg band: Optional, I{polar} Band (C{str}, 'A'|'B'|'Y'|'Z').
           @kwarg datum: Optional, this coordinate's datum (L{Datum},
                         L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg falsed: If C{True}, both B{C{easting}} and B{C{northing}}
                          are falsed (C{bool}).
           @kwarg convergence: Optional, meridian convergence gamma
                               to save (C{degrees}).
           @kwarg scale: Optional, computed scale factor k to save
                         (C{scalar}).
           @kwarg name: Optional name (C{str}).

           @raise TypeError: Invalid B{C{datum}}.

           @raise UPSError: Invalid B{C{zone}}, B{C{pole}}, B{C{easting}},
                            B{C{northing}}, B{C{band}}, B{C{convergence}}
                            or B{C{scale}}.
        '''
        if name:
            self.name = name

        try:
            z, B, p = _to3zBhp(zone, band, hemipole=pole)
            if z != _UPS_ZONE or (B and (B not in _Bands)):
                raise ValueError
        except (TypeError, ValueError) as x:
            raise UPSError(zone=zone, pole=pole, band=band, txt=str(x))
        self._pole = p
        UtmUpsBase.__init__(self, easting, northing, band=B, datum=datum, falsed=falsed,
                                                     convergence=convergence, scale=scale)

    def __eq__(self, other):
        return isinstance(other, Ups) and other.zone     == self.zone \
                                      and other.pole     == self.pole \
                                      and other.easting  == self.easting \
                                      and other.northing == self.northing \
                                      and other.band     == self.band \
                                      and other.datum    == self.datum

    @property_doc_(''' the I{polar} band.''')
    def band(self):
        '''Get the I{polar} band (C{'A'|'B'|'Y'|'Z'}).
        '''
        if not self._band:
            self._toLLEB()
        return self._band

    @band.setter  # PYCHOK setter!
    def band(self, band):
        '''Set or reset the I{polar} band.

           @arg band: The I{polar} band letter (C{'A'|'B'|'Y'|'Z'})
                      or C{None} or C{""} to reset.

           @raise TypeError: Invalid B{C{band}}.

           @raise ValueError: Invalid B{C{band}}.
        '''
        self._band1(band)

    @Property_RO
    def falsed2(self):
        '''Get the easting and northing falsing (L{EasNor2Tuple}C{(easting, northing)}).
        '''
        f = _Falsing if self.falsed else 0
        return EasNor2Tuple(f, f)

    @Property_RO
    def _mgrs(self):
        '''(INTERNAL) Cache for L{toMgrs}.
        '''
        return self.toUtm(None).toMgrs()

    def parse(self, strUPS, name=NN):
        '''Parse a string to a similar L{Ups} instance.

           @arg strUPS: The UPS coordinate (C{str}),
                        see function L{parseUPS5}.
           @kwarg name: Optional instance name (C{str}),
                        overriding this name.

           @return: The similar instance (L{Ups}).

           @raise UTMError: Invalid B{C{strUPS}}.

           @see: Functions L{parseUTM5} and L{pygeodesy.parseUTMUPS5}.
        '''
        return parseUPS5(strUPS, datum=self.datum, Ups=self.classof,
                                 name=name or self.name)

    @deprecated_method
    def parseUPS(self, strUPS):  # PYCHOK no cover
        '''DEPRECATED, use method L{parse}.'''
        return self.parse(strUPS)

    @Property_RO
    def pole(self):
        '''Get the top/center of (stereographic) projection (C{'N'|'S'} or C{""}).
        '''
        return self._pole

    def rescale0(self, lat, scale0=_K0_UPS):
        '''Set the central scale factor for this UPS projection.

           @arg lat: Northern latitude (C{degrees}).
           @arg scale0: UPS k0 scale at B{C{lat}} latitude (C{scalar}).

           @raise RangeError: If B{C{lat}} outside the valid range and
                              L{pygeodesy.rangerrors} set to C{True}.

           @raise UPSError: Invalid B{C{scale}}.
        '''
        s0 = Scalar_(scale0=scale0, Error=UPSError, low=EPS)  # <= 1.003 or 1.0016?
        u  = toUps8(abs(Lat(lat)), _0_0, datum=self.datum, Ups=_Ups_K1)
        k  = s0 / u.scale
        if self.scale0 != k:
            _update_all(self)
            self._band = NN  # force re-compute
            self._latlon = self._utm = None
            self._scale0 = Scalar(scale0=k)

    def toLatLon(self, LatLon=None, unfalse=True, **LatLon_kwds):
        '''Convert this UPS coordinate to an (ellipsoidal) geodetic point.

           @kwarg LatLon: Optional, ellipsoidal class to return the
                          geodetic point (C{LatLon}) or C{None}.
           @kwarg unfalse: Unfalse B{C{easting}} and B{C{northing}}
                           if falsed (C{bool}).
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                               arguments, ignored if C{B{LatLon} is None}.

           @return: This UPS coordinate (B{C{LatLon}}) or if B{C{LatLon}}
                    is C{None}, a L{LatLonDatum5Tuple}C{(lat, lon, datum,
                    convergence, scale)}.

           @raise TypeError: If B{C{LatLon}} is not ellipsoidal.

           @raise UPSError: Invalid meridional radius or H-value.
        '''
        if self._latlon and self._latlon._toLLEB_args == (unfalse,):
            return self._latlon5(LatLon)
        else:
            self._toLLEB(unfalse=unfalse)
            return self._latlon5(LatLon, **LatLon_kwds)

    def _toLLEB(self, unfalse=True):  # PYCHOK signature
        '''(INTERNAL) Compute (ellipsoidal) lat- and longitude.
        '''
        E = self.datum.ellipsoid  # XXX vs LatLon.datum.ellipsoid

        x, y = self.eastingnorthing2(falsed=not unfalse)

        r = hypot(x, y)
        t = (r / (_2_0 * self.scale0 * E.a / E.es_c)) if r > 0 else EPS0
        t = E.es_tauf((1 / t - t) * _0_5)
        a = atan(t)
        if self._pole == _N_:
            b, c = atan2(x, -y), 1
        else:
            a, b, c = neg(a), atan2(x, y), -1
        ll = _LLEB(degrees90(a), degrees180(b), datum=self._datum, name=self.name)

        ll._convergence = b * c  # gamma
        ll._scale = _scale(E, r, t) if r > 0 else self.scale0
        self._latlon5args(ll, _toBand, unfalse)

    def toMgrs(self):
        '''Convert this UPS coordinate to an MGRS grid reference.

           @return: The MGRS grid reference (L{Mgrs}).

           @see: Methods L{Ups.toUtm} and L{Utm.toMgrs}.

           @note: If not specified, the I{polar} C{band} is computed
                  from the (geodetic) latitude and the C{datum}.
        '''
        return self._mgrs

    def toRepr(self, prec=0, fmt=Fmt.SQUARE, sep=_COMMASPACE_, B=False, cs=False, **unused):  # PYCHOK expected
        '''Return a string representation of this UPS coordinate.

           Note that UPS coordinates are rounded, not truncated (unlike
           MGRS grid references).

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg fmt: Optional, enclosing backets format (C{str}).
           @kwarg sep: Optional separator between name:value pairs (C{str}).
           @kwarg B: Optionally, include polar band letter (C{bool}).
           @kwarg cs: Optionally, include gamma meridian convergence and
                      point scale factor (C{bool} or non-zero C{int} to
                      specify the precison like B{C{prec}}).

           @return: This UPS as a string with C{00[Band] pole, easting,
                    northing, [convergence, scale]} as C{"[Z:00[Band],
                    P:N|S, E:meter, N:meter]"} plus C{", C:DMS, S:float"}
                    if B{C{cs}} is C{True}, where C{[Band]} is present and
                    C{'A'|'B'|'Y'|'Z'} only if B{C{B}} is C{True} and
                    convergence C{DMS} is in I{either} degrees, minutes
                    I{or} seconds (C{str}).

           @note: Pseudo zone zero (C{"00"}) for UPS follows I{Karney}'s U{zone UPS
                  <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1UTMUPS.html>}.
        '''
        return self._toRepr(fmt, B, cs, prec, sep)

    def toStr(self, prec=0, sep=_SPACE_, B=False, cs=False):  # PYCHOK expected
        '''Return a string representation of this UPS coordinate.

           Note that UPS coordinates are rounded, not truncated (unlike
           MGRS grid references).

           @kwarg prec: Optional number of decimals, unstripped (C{int}).
           @kwarg sep: Optional separator to join (C{str}) or C{None}
                       to return an unjoined C{tuple} of C{str}s.
           @kwarg B: Optionally, include and polar band letter (C{bool}).
           @kwarg cs: Optionally, include gamma meridian convergence and
                      point scale factor (C{bool} or non-zero C{int} to
                      specify the precison like B{C{prec}}).

           @return: This UPS as a string with C{00[Band] pole, easting,
                    northing, [convergence, scale]} as C{"00[B] N|S
                    meter meter"} plus C{" DMS float"} if B{C{cs}} is C{True},
                    where C{[Band]} is present and C{'A'|'B'|'Y'|'Z'} only
                    if B{C{B}} is C{True} and convergence C{DMS} is in
                    I{either} degrees, minutes I{or} seconds (C{str}).

           @note: Zone zero (C{"00"}) for UPS follows I{Karney}'s U{zone UPS
                  <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1UTMUPS.html>}.
        '''
        return self._toStr(self.pole, B, cs, prec, sep)  # PYCHOK pole

    def toUps(self, pole=NN, **unused):
        '''Duplicate this UPS coordinate.

           @kwarg pole: Optional top/center of the UPS projection,
                        (C{str}, 'N[orth]'|'S[outh]').

           @return: A copy of this UPS coordinate (L{Ups}).

           @raise UPSError: Invalid B{C{pole}} or attempt to transfer
                            the projection top/center.
        '''
        if self.pole == pole or not pole:
            return self.copy()
        t = _SPACE_(_pole_, repr(self.pole), _to_, repr(pole))
        raise UPSError('no transfer', txt=t)

    def toUtm(self, zone, falsed=True, **unused):
        '''Convert this UPS coordinate to a UTM coordinate.

           @arg zone: The UTM zone (C{int}).
           @kwarg falsed: False both easting and northing (C{bool}).

           @return: The UTM coordinate (L{Utm}).
        '''
        u = self._utm
        if u is None or u.zone != zone or falsed != u.falsed:
            from pygeodesy.utm import toUtm8, Utm
            ll = self.toLatLon(LatLon=None, unfalse=True)
            self._utm = toUtm8(ll, Utm=Utm, falsed=falsed, name=self.name, zone=zone)
        return self._utm

    @Property_RO
    def zone(self):
        '''Get the polar pseudo zone (C{0}), like I{Karney}'s U{zone UPS<https://
           GeographicLib.SourceForge.io/html/classGeographicLib_1_1UTMUPS.html>}.
        '''
        return _UPS_ZONE


class _Ups_K1(Ups):
    '''(INTERNAL) For method L{Ups.rescale0}.
    '''
    _scale0 = _K1_UPS


def parseUPS5(strUPS, datum=_WGS84, Ups=Ups, falsed=True, name=NN):
    '''Parse a string representing a UPS coordinate, consisting of
       C{"[zone][band] pole easting northing"} where B{C{zone}} is
       pseudo zone C{"00"|"0"|""} and C{band} is C{'A'|'B'|'Y'|'Z'|''}.

       @arg strUPS: A UPS coordinate (C{str}).
       @kwarg datum: Optional datum to use (L{Datum}).
       @kwarg Ups: Optional class to return the UPS coordinate (L{Ups})
                   or C{None}.
       @kwarg falsed: Both B{C{easting}} and B{C{northing}} are falsed (C{bool}).
       @kwarg name: Optional B{C{Ups}} name (C{str}).

       @return: The UPS coordinate (B{C{Ups}}) or a
                L{UtmUps5Tuple}C{(zone, hemipole, easting, northing,
                band)} if B{C{Ups}} is C{None}.  The C{hemipole} is
                the C{'N'|'S'} pole, the UPS projection top/center.

       @raise UPSError: Invalid B{C{strUPS}}.
    '''
    z, p, e, n, B = _parseUTMUPS5(strUPS, _UPS_ZONE_STR, Error=UPSError)
    if z != _UPS_ZONE or (B and B not in _Bands):
        raise UPSError(strUPS=strUPS, zone=z, band=B)

    if Ups is None:
        r = UtmUps5Tuple(z, p, e, n, B, Error=UPSError, name=name)
    else:
        r = Ups(z, p, e, n, band=B, falsed=falsed, datum=datum)
        r = _xnamed(r, name)
    return r


def toUps8(latlon, lon=None, datum=None, Ups=Ups, pole=NN,
                             falsed=True, strict=True, name=NN):
    '''Convert a lat-/longitude point to a UPS coordinate.

       @arg latlon: Latitude (C{degrees}) or an (ellipsoidal)
                    geodetic C{LatLon} point.
       @kwarg lon: Optional longitude (C{degrees}) or C{None} if
                   B{C{latlon}} is a C{LatLon}.
       @kwarg datum: Optional datum for this UPS coordinate,
                     overriding B{C{latlon}}'s datum (C{Datum},
                     L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}).
       @kwarg Ups: Optional class to return the UPS coordinate
                   (L{Ups}) or C{None}.
       @kwarg pole: Optional top/center of (stereographic) projection
                    (C{str}, C{'N[orth]'} or C{'S[outh]'}).
       @kwarg falsed: False both easting and northing (C{bool}).
       @kwarg strict: Restrict B{C{lat}} to UPS ranges (C{bool}).
       @kwarg name: Optional B{C{Ups}} name (C{str}).

       @return: The UPS coordinate (B{C{Ups}}) or a
                L{UtmUps8Tuple}C{(zone, hemipole, easting, northing,
                band, datum, convergence, scale)} if B{C{Ups}} is
                C{None}.  The C{hemipole} is the C{'N'|'S'} pole,
                the UPS projection top/center.

       @raise RangeError: If B{C{strict}} and B{C{lat}} outside the valid
                          UPS bands or if B{C{lat}} or B{C{lon}} outside
                          the valid range and L{pygeodesy.rangerrors} set
                          to C{True}.

       @raise TypeError: If B{C{latlon}} is not ellipsoidal or
                         B{C{datum}} invalid.

       @raise ValueError: If B{C{lon}} value is missing or if B{C{latlon}}
                          is invalid.

       @see: I{Karney}'s C++ class U{UPS
             <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1UPS.html>}.
    '''
    lat, lon, d, name = _to4lldn(latlon, lon, datum, name)
    z, B, p, lat, lon = upsZoneBand5(lat, lon, strict=strict)  # PYCHOK UtmUpsLatLon5Tuple

    d = _ellipsoidal_datum(d, name=name)
    E = d.ellipsoid

    p = str(pole or p)[:1].upper()
    N = p == _N_  # is north

    a = lat if N else -lat
    A = abs(a - _90_0) < _Tol90  # at pole

    t = tan(radians(a))
    T = E.es_taupf(t)
    r = hypot1(T) + abs(T)
    if T >= _0_0:
        r = _0_0 if A else _1_0 / r

    k0 = getattr(Ups, '_scale0', _K0_UPS)  # Ups is class or None
    r *= 2 * k0 * E.a / E.es_c

    k = k0 if A else _scale(E, r, t)
    c = lon  # [-180, 180) from .upsZoneBand5
    x, y = sincos2d(c)
    x *= r
    y *= r
    if N:
        y = neg(y)
    else:
        c = neg(c)

    if falsed:
        x += _Falsing
        y += _Falsing

    n = name or nameof(latlon)
    if Ups is None:
        r = UtmUps8Tuple(z, p, x, y, B, d, c, k, Error=UPSError, name=n)
    else:
        if z != _UPS_ZONE and not strict:
            z = _UPS_ZONE  # ignore UTM zone
        r = _xnamed(Ups(z, p, x, y, band=B, datum=d, falsed=falsed,
                                            convergence=c, scale=k), n)
        if isinstance(latlon, _LLEB) and d is latlon.datum:  # see utm._toXtm8
            r._latlon5args(latlon, _toBand, falsed)  # XXX weakref(latlon)?
            latlon._convergence = c
            latlon._scale = k
        else:
            r._hemisphere = _hemi(lat)
            if not r._band:
                r._band = _toBand(lat, lon)
    return r


def upsZoneBand5(lat, lon, strict=True, name=NN):
    '''Return the UTM/UPS zone number, I{polar} Band letter, pole and
       clipped lat- and longitude for a given location.

       @arg lat: Latitude in degrees (C{scalar} or C{str}).
       @arg lon: Longitude in degrees (C{scalar} or C{str}).
       @kwarg strict: Restrict B{C{lat}} to UPS ranges (C{bool}).
       @kwarg name: Optional name (C{str}).

       @return: A L{UtmUpsLatLon5Tuple}C{(zone, band, hemipole,
                lat, lon)} where C{hemipole} is the C{'N'|'S'} pole,
                the UPS projection top/center.

       @raise RangeError: If B{C{strict}} and B{C{lat}} in the UTM
                          and not the UPS range or if B{C{lat}} or
                          B{C{lon}} outside the valid range and
                          L{pygeodesy.rangerrors} set to C{True}.

       @raise ValueError: Invalid B{C{lat}} or B{C{lon}}.
    '''
    z, lat, lon = _to3zll(*parseDMS2(lat, lon))
    if lat < _UPS_LAT_MIN:  # includes 30' overlap
        z, B, p = _UPS_ZONE, _toBand(lat, lon), _S_

    elif lat > _UPS_LAT_MAX:  # includes 30' overlap
        z, B, p = _UPS_ZONE, _toBand(lat, lon), _N_

    elif strict:
        r = _range_(_UPS_LAT_MIN, _UPS_LAT_MAX)
        t = _SPACE_(_inside_, _UTM_, _range_, r)
        raise RangeError(lat=degDMS(lat), txt=t)

    else:
        B, p = NN, _hemi(lat)
    return UtmUpsLatLon5Tuple(z, B, p, lat, lon, Error=UPSError, name=name)

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
