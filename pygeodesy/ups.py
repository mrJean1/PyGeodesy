
# -*- coding: utf-8 -*-

u'''Universal Polar Stereographic (UPS) classes L{Ups} and L{UPSError}
and functions L{parseUPS5}, L{toUps8} and L{upsZoneBand5}.

A pure Python implementation, partially transcribed from C++ class U{PolarStereographic
<http://GeographicLib.SourceForge.io/html/classGeographicLib_1_1PolarStereographic.html>}
by I{Charles Karney}.

The U{UPS<http://WikiPedia.org/wiki/Universal_polar_stereographic_coordinate_system>}
system is used in conjuction with U{UTM
<http://WikiPedia.org/wiki/Universal_Transverse_Mercator_coordinate_system>}
to location on the polar regions of the earth.  UPS covers areas south of 79.5°S
and north of 83.5°N (slightly overlapping the UTM range from 80°S to 84°N).

@newfield example: Example, Examples
'''

from bases import _Based, _xattrs, _xnamed
from datum import Datums, _TOL
from dms import clipDMS, degDMS, parseDMS2, _parseUTMUPS, RangeError
from ellipsoidalBase import LatLonEllipsoidalBase as _LLEB, _hemi, \
                           _to4lldn, _to3zBhp, _to3zll, _UPS_LAT_MAX, \
                           _UPS_LAT_MIN, _UPS_ZONE, _UPS_ZONE_STR
from fmath import EPS, fStr, hypot, hypot1
from lazily import _ALL_LAZY
from utily import degrees90, degrees180, property_RO, radians, sincos2d

from math import atan, atan2, sqrt, tan

# all public contants, classes and functions
__all__ = _ALL_LAZY.ups
__version__ = '19.04.26'

_Bands   = 'A', 'B', 'Y', 'Z'    #: (INTERNAL) Polar bands.
_Falsing = 2000e3  #: (INTERNAL) False easting and northing (C{meter}).
_K0      = 0.994   #: (INTERNAL) Central UPS scale factor.
_K1      = 1.0     #: (INTERNAL) Rescale point scale factor.


class UPSError(ValueError):
    '''UPS parse or other error.
    '''
    pass


def _Band(a, b):
    # determine the polar band letter
    return _Bands[(0 if a < 0 else 2) + (0 if b < 0 else 1)]


def _scale(E, rho, tau):
    # compute the point scale factor, ala Karney
    t = hypot1(tau)
    return (rho / E.a) * t * sqrt(E.e12 + E.e2 / t**2)


class Ups(_Based):
    '''Universal Polar Stereographic (UPS) coordinate.
    '''
    _band        = ''    #: (INTERNAL) Polar band ('A', 'B', 'Y' or 'Z').
    _convergence = None  #: (INTERNAL) Gamma meridian conversion (C{degrees}).
    _datum       = Datums.WGS84  #: (INTERNAL) L{Datum}.
    _easting     = 0     #: (INTERNAL) Easting from false easting (C{meter}).
    _epsg        = None  #: (INTERNAL) toEpsg cache (L{Epsg}).
    _falsed      = 0     #: (INTERANL) Falsed easting and northing (C{meter}).
    _hemisphere  = ''    #: (INTERNAL) Hemisphere ('N' or 'S'), different from pole.
    # _latlon also set by ellipsoidalBase.LatLonEllipsoidalBase.toUtm, .toUps and .toUtmUps
    _latlon      = None  #: (INTERNAL) toLatLon cache (C{LatLon}).
    _mgrs        = None  #: (INTERNAL) toMgrs cache (L{Mgrs}).
    _northing    = 0     #: (INTERNAL) Northing from false northing (C{meter}).
    _pole        = ''    #: (INTERNAL) Pojection top/center ('N' or 'S').
    _scale       = None  #: (INTERNAL) Point scale factor (C{scalar}).
    _scale0      = _K0   #: (INTERNAL) Central scale factor (C{scalar}).
    _utm         = None  #: (INTERNAL) toUtm cache (L{Utm}).

    def __init__(self, zone, pole, easting, northing, band='',  # PYCHOK expected
                                   datum=Datums.WGS84, falsed=True,
                                   convergence=None, scale=None, name=''):
        '''New UPS coordinate.

           @param zone: UPS zone (C{int}, zero) or zone with/-out Band
                        letter (C{str}, '00', '00A', '00B', '00Y' or '00Z').
           @param pole: Top/center of (stereographic) projection
                        (C{str}, C{'N[orth]'} or C{'S[outh]'}).
           @param easting: Easting (C{meter}).
           @param northing: Northing (C{meter}).
           @keyword band: Optional, polar Band (C{str}, 'A'|'B'|'Y'|'Z').
           @keyword datum: Optional, this coordinate's datum (L{Datum}).
           @keyword falsed: Both I{easting} and I{northing} are falsed (C{str}).
           @keyword convergence: Optionally, save gamma meridian convergence
                                 (C{degrees}).
           @keyword scale: Optionally, save computed k scale (C{scalar}).
           @keyword name: Optional name (C{str}).

           @raise UPSError: Invalid I{zone}, I{pole} or I{band}.
        '''
        if name:
            self.name = name

        try:
            z, B, p = _to3zBhp(zone, band=band, hemipole=pole)
            if z != _UPS_ZONE or (B and B not in _Bands):
                raise ValueError
        except ValueError:
            raise UPSError('%s, %s or %s invalid: %r' %
                           ('zone', 'pole', 'band', (zone, pole, band)))

        self._band        = B
        self._convergence = convergence
        self._easting     = float(easting)
        self._falsed      = _Falsing if falsed else 0
        self._northing    = float(northing)
        self._pole        = p
        self._datum       = datum
        self._scale       = scale

    def __repr__(self):
        return self.toStr2(B=True)

    def __str__(self):
        return self.toStr()

    def _xcopy(self, *attrs):
        '''(INTERNAL) Make copy with add'l, subclass attributes.
        '''
        return _xattrs(self.classof(self.pole, self.easting, self.northing,
                                    band=self.band, datum=self.datum,
                                    convergence=self.convergence,
                                    scale=self.scale),
                       self, *attrs)

    @property_RO
    def band(self):
        '''Get the polar band letter ('A', 'B', 'Y' or 'Z').
        '''
        if not self._band:
            self.toLatLon(unfalse=True)
        return self._band

    @property_RO
    def convergence(self):
        '''Get the gamma meridian convergence (C{degrees}) or C{None}).
        '''
        return self._convergence

    def copy(self):
        '''Copy this UPS coordinate.

           @return: The copy (L{Ups} or subclass thereof).
        '''
        return self._xcopy()

    @property_RO
    def datum(self):
        '''Get the datum (L{Datum}).
        '''
        return self._datum

    @property_RO
    def easting(self):
        '''Get the (falsed) easting (C{meter}).
        '''
        return self._easting

    @property_RO
    def falsed(self):
        '''Get the easting and northing falsing (C{meter}) or C{0}.
        '''
        return self._falsed

    @property_RO
    def hemisphere(self):
        '''Get the hemisphere (C{str}, 'N'|'S').
        '''
        if not self._hemisphere:
            self.toLatLon(unfalse=True)
        return self._hemisphere

    @property_RO
    def northing(self):
        '''Get the (falsed) northing (C{meter}).
        '''
        return self._northing

    def parseUPS(self, strUPS):
        '''Parse a string to a UPS coordinate.

           @return: The coordinate (L{Ups}).

           @see: Function L{parseUPS5} in this module L{ups}.
        '''
        return parseUPS5(strUPS, datum=self.datum, Ups=self.classof)

    @property_RO
    def pole(self):
        '''Get the center of (stereographic) projection (N|S).
        '''
        return self._pole

    def rescale0(self, lat, scale0=_K0):
        '''Set the central scale factor for this UPS projection.

           @param lat: Northern latitude (C{degrees}).
           @param scale0: UPS k0 scale at I{lat} latitude (C{scalar}).

           @raise RangeError: If I{lat} outside the valid range
                              and I{rangerrrors} set to C{True}.

           @raise ValueError: Invalid I{scale}.
        '''
        try:
            s0 = float(scale0)
            if not 0 < s0:  # <= 1.003 or 1.0016?
                raise ValueError
        except (TypeError, ValueError):
            raise ValueError('%s invalid: %r' % ('scale', scale0))

        lat = clipDMS(lat, 90)  # clip and force N
        u = toUps8(abs(lat), 0, datum=self.datum, Ups=_UpsK1)
        k = s0 / u.scale
        if self._scale0 != k:
            self._band = ''  # force re-compute
            self._latlon = self._mgrs = self._utm = None
            self._scale0 = k

    @property_RO
    def scale(self):
        '''Get the point scale factor (C{scalar} or C{None}).
        '''
        return self._scale

    @property_RO
    def scale0(self):
        '''Get the central scale factor (C{scalar} or C{None}).
        '''
        return self._scale0

    def toEpsg(self):
        '''Determine the I{EPSG (European Petroleum Survey Group)} code.

           @return: C{EPSG} code (C{int}).

           @raise EPSGError: See L{Epsg}.
        '''
        if self._epsg is None:
            from epsg import Epsg  # PYCHOK circular import
            self._epsg = Epsg(self)
        return self._epsg

    def toLatLon(self, LatLon=None, unfalse=True):
        '''Convert this UPS coordinate to an (ellipsoidal) geodetic point.

           @keyword LatLon: Optional, ellipsoidal (sub-)class to return
                            the point (C{LatLon}) or C{None}.
           @keyword unfalse: Unfalse I{easting} and I{northing} if falsed
                             (C{bool}).

           @return: This UPS coordinate as (I{LatLon}) or 5-tuple
                    (C{lat, lon, datum, convergence, scale}) if
                    I{LatLon} is C{None}.

           @raise TypeError: If I{LatLon} is not ellipsoidal.

           @raise UPSError: Invalid meridional radius or H-value.
        '''
        if self._latlon:
            return self._latlon5(LatLon)

        E = self.datum.ellipsoid  # XXX vs LatLon.datum.ellipsoid

        f = self.falsed if unfalse else 0
        x = self.easting  - f
        y = self.northing - f

        r = hypot(x, y)
        t = (r / (2 * self.scale0 * E.a / E.es_c)) if r > 0 else EPS**2
        t = E.es_tauf((1 / t - t) * 0.5)
        if self.pole == 'N':
            a, b, c = atan(t), atan2(x, -y), 1
        else:
            a, b, c = -atan(t), atan2(x, y), -1

        a, b = degrees90(a), degrees180(b)
        if not self._band:
            self._band = _Band(a, b)
        if not self._hemisphere:
            self._hemisphere = _hemi(a)

        ll = _LLEB(a, b, datum=self._datum, name=self.name)
        ll._convergence = b * c  # gamma
        ll._scale = _scale(E, r, t) if r > 0 else self.scale0

        self._latlon = ll
        return self._latlon5(LatLon)

    def _latlon5(self, LatLon):
        '''(INTERNAL) Convert cached LatLon
        '''
        ll = self._latlon
        if LatLon is None:
            return ll.lat, ll.lon, ll.datum, ll.convergence, ll.scale
        elif issubclass(LatLon, _LLEB):
            return _xnamed(_xattrs(LatLon(ll.lat, ll.lon, datum=ll.datum),
                                   ll, '_convergence', '_scale'), ll.name)
        raise TypeError('%s not ellipsoidal: %r' % ('LatLon', LatLon))

    def toMgrs(self):
        '''Convert this UPS coordinate to an MGRS grid reference.

           @return: The MGRS grid reference (L{Mgrs}).

           @see: Methods L{Ups.toUtm} and L{Utm.toMgrs}.
        '''
        if self._mgrs is None:
            self._mgrs = self.toUtm(None).toMgrs()
        return self._mgrs

    def toStr(self, prec=0, sep=' ', B=False, cs=False):  # PYCHOK expected
        '''Return a string representation of this UPS coordinate.

           Note that UPS coordinates are rounded, not truncated
           (unlike MGRS grid references).

           @keyword prec: Optional number of decimals, unstripped (C{int}).
           @keyword sep: Optional separator to join (C{str}).
           @keyword B: Optionally, include and polar band letter (C{bool}).
           @keyword cs: Optionally, include gamma meridian convergence
                        and point scale factor (C{bool}).

           @return: This UPS as a string with I{00[Band] pole, easting,
                    northing, [convergence, scale]} as C{"00[B] N|S
                    meter meter"} plus C{" DMS float"} if I{cs} is C{True},
                    where C{[Band]} is present and C{'A'|'B'|'Y'|'Z'} only
                    if I{B} is C{True} and convergence C{DMS} is in
                    I{either} degrees, minutes I{or} seconds (C{str}).

           @note: Zone zero (C{"00"}) for UPS follows Karney's U{zone UPS
                  <http://GeographicLib.SourceForge.io/html/classGeographicLib_1_1UTMUPS.html>}.
        '''
        z = _UPS_ZONE_STR + (self.band if B else '')
        t = (z, self.pole, fStr(self.easting,  prec=prec),
                           fStr(self.northing, prec=prec))
        if cs:
            t += ('n/a' if self.convergence is None else
                    degDMS(self.convergence, prec=8, pos='+'),
                  'n/a' if self.scale is None else
                      fStr(self.scale, prec=8))
        return sep.join(t)

    def toStr2(self, prec=0, fmt='[%s]', sep=', ', B=False, cs=False):  # PYCHOK expected
        '''Return a string representation of this UPS coordinate.

           Note that UPS coordinates are rounded, not truncated
           (unlike MGRS grid references).

           @keyword prec: Optional number of decimals, unstripped (C{int}).
           @keyword fmt: Optional, enclosing backets format (C{str}).
           @keyword sep: Optional separator between name:value pairs (C{str}).
           @keyword B: Optionally, include polar band letter (C{bool}).
           @keyword cs: Optionally, include gamma meridian convergence
                        and point scale factor (C{bool}).

           @return: This UPS as a string  with I{00[Band] pole, easting,
                    northing, [convergence, scale]} as C{"[Z:00[Band],
                    P:N|S, E:meter, N:meter]"} plus C{", C:DMS, S:float"}
                    if I{cs} is C{True}, where C{[Band]} is present and
                    C{'A'|'B'|'Y'|'Z'} only if I{B} is C{True} and
                    convergence C{DMS} is in I{either} degrees, minutes
                    I{or} seconds (C{str}).

           @note: Zone zero (C{"00"}) for UPS follows Karney's U{zone UPS
                  <http://GeographicLib.SourceForge.io/html/classGeographicLib_1_1UTMUPS.html>}.
        '''
        t = self.toStr(prec=prec, sep=' ', B=B, cs=cs).split()
        return fmt % (sep.join('%s:%s' % t for t in zip('ZPENCS', t)),)

    def toUps(self, pole='', **unused):
        '''Duplicate this UPS coordinate.

           @keyword pole: Optional top/center of the UPS projection,
                          (C{str}, 'N[orth]'|'S[outh]').

           @return: A copt of this UPS coordinate (L{Ups}).

           @raise UPSError: Invalid I{pole} or attempt to transfer
                            the projection top/center.
        '''
        if self.pole == pole or not pole:
            return self.copy()
        raise UPSError('%s transfer invalid: %r to %r' %
                       ('pole', self.pole, pole))

    def toUtm(self, zone, **unused):
        '''Convert this UPS coordinate to a UTM coordinate.

           @param zone: The UTM zone (C{int}).

           @return: The UTM coordinate (L{Utm}).
        '''
        u = self._utm
        if u is None or u.zone != zone:
            from utm import toUtm8  # PYCHOK recursive import
            ll = self.toLatLon(LatLon=None, unfalse=True)
            self._utm = toUtm8(ll, name=self.name, zone=zone)
        return self._utm

    @property_RO
    def zone(self):
        '''Get the polar zone (c{int}), like Karney's U{zone UPS<http://
           GeographicLib.SourceForge.io/html/classGeographicLib_1_1UTMUPS.html>}.
        '''
        return _UPS_ZONE


class _UpsK1(Ups):
    '''(INTEERNAL) For method L{Ups.rescale}.
    '''
    _scale0 = _K1


def parseUPS5(strUPS, datum=Datums.WGS84, Ups=Ups, falsed=True, name=''):
    '''Parse a string representing a UPS coordinate, consisting of
       I{"[zone][band] pole easting northing"} where I{zone} is pseudo
       zone I{"00"|"0"|""} and I{band} is I{'A'|'B'|'Y'|'Z'|''}.

       @param strUPS: A UPS coordinate (C{str}).
       @keyword datum: Optional datum to use (L{Datum}).
       @keyword Ups: Optional (sub-)class to return the UPS
                     coordinate (L{Ups}) or C{None}.
       @keyword falsed: Both I{easting} and I{northing} are falsed (C{bool}).
       @keyword name: Optional I{Ups} name (C{str}).

       @return: The UPS coordinate (L{Ups}) or 5-tuple (C{zone, pole,
                easting, northing, band}) if I{Ups} is C{None}.

       @raise UPSError: Invalid I{strUPS}.
    '''
    try:
        u = strUPS.lstrip()
        if not u.startswith(_UPS_ZONE_STR):
            raise ValueError

        z, p, e, n, B = _parseUTMUPS(u)
        if z != _UPS_ZONE or (B and B not in _Bands):
            raise ValueError
    except (AttributeError, TypeError, ValueError):
        raise UPSError('%s invalid: %r' % ('strUPS', strUPS))

    return (z, p, e, n, B) if Ups is None else _xnamed(Ups(
            z, p, e, n, band=B, falsed=falsed, datum=datum), name)


def toUps8(latlon, lon=None, datum=None, Ups=Ups, pole='',
                             falsed=True, strict=True, name=''):
    '''Convert a lat-/longitude point to a UPS coordinate.

       @param latlon: Latitude (C{degrees}) or an (ellipsoidal)
                      geodetic C{LatLon} point.
       @keyword lon: Optional longitude (C{degrees}) or C{None}
                     if I{latlon} is a C{LatLon}.
       @keyword datum: Optional datum for this UPS coordinate,
                       overriding I{latlon}'s datum (C{Datum}).
       @keyword Ups: Optional (sub-)class to return the UPS
                     coordinate (L{Ups}) or C{None}.
       @keyword pole: Optional top/center of (stereographic) projection
                      (C{str}, C{'N[orth]'} or C{'S[outh]'}).
       @keyword falsed: False both easting and northing (C{bool}).
       @keyword strict: Restrict I{lat} to UPS ranges (C{bool}).
       @keyword name: Optional I{Ups} name (C{str}).

       @return: The UPS coordinate (L{Ups}) or a 8-tuple (zone, hemisphere,
                easting, northing, band, datum, convergence, scale) if
                I{Ups} is C{None}.

       @raise RangeError: If I{strict} and I{lat} outside the valid UPS
                          bands or if I{lat} or I{lon} outside the valid
                          range and I{rangerrrors} set to C{True}.

       @raise TypeError: If I{latlon} is not ellipsoidal.

       @raise ValueError: If I{lon} value is missing or if I{latlon}
                          is invalid.

       @see: Karney's C++ class U{UPS
             <http://GeographicLib.SourceForge.io/html/classGeographicLib_1_1UPS.html>}.
    '''
    lat, lon, d, name = _to4lldn(latlon, lon, datum, name)
    z, B, p, lat, lon = upsZoneBand5(lat, lon, strict=strict)

    p = str(pole or p)[:1]
    N = p in 'Nn'
    E = d.ellipsoid
    A = abs(lat - 90) < _TOL

    t = tan(radians(lat if N else -lat))
    T = E.es_taupf(t)
    r = hypot1(T) + abs(T)
    if T >= 0:
        r = 0 if A else 1 / r

    k0 = getattr(Ups, '_scale0', _K0)  # Ups is class or None
    r *= 2 * k0 * E.a / E.es_c

    k = k0 if A else _scale(E, r, t)
    c = lon  # [-180, 180) from .upsZoneBand5
    x, y = sincos2d(c)
    x *= r
    y *= r
    if N:
        y = -y
    else:
        c = -c

    if falsed:
        x += _Falsing
        y += _Falsing

    if Ups is None:
        r = z, p, x, y, B, d, c, k
    else:
        if z != _UPS_ZONE and not strict:
            z = _UPS_ZONE  # ignore UTM zone
        r = _xnamed(Ups(z, p, x, y, band=B, datum=d,
                                    convergence=c, scale=k,
                                    falsed=falsed), name)
        if hasattr(Ups, '_hemisphere'):
            r._hemisphere = _hemi(lat)
    return r


def upsZoneBand5(lat, lon, strict=True):
    '''Return the UTM/UPS zone number, (polar) Band letter, pole and
       clipped lat- and longitude for a given location.

       @param lat: Latitude in degrees (C{scalar} or C{str}).
       @param lon: Longitude in degrees (C{scalar} or C{str}).
       @keyword strict: Restrict I{lat} to UPS ranges (C{bool}).

       @return: 5-Tuple (C{zone, Band, hemisphere, lat, lon}) as
                (C{int, str, 'N'|'S', degrees90, degrees180}) where
                C{zone} is always C{0} for UPS and polar C{Band} is
                C{""} or C{'A'|'B'|'Y'|'Z'}.

       @raise RangeError: If I{strict} and I{lat} in UTM and not UPS
                          range or if I{lat} or I{lon} outside the valid
                          range and I{rangerrrors} set to C{True}.

       @raise ValueError: Invalid I{lat} or I{lon}.
    '''
    z, lat, lon = _to3zll(*parseDMS2(lat, lon))
    if lat < _UPS_LAT_MIN:  # includes 30' overlap
        return _UPS_ZONE, _Band(lat, lon), 'S', lat, lon

    elif lat > _UPS_LAT_MAX:  # includes 30' overlap
        return _UPS_ZONE, _Band(lat, lon), 'N', lat, lon

    elif strict:
        x = '%s [%s, %s]' % ('range', _UPS_LAT_MIN, _UPS_LAT_MAX)
        raise RangeError('%s inside UTM %s: %s' % ('lat', x, degDMS(lat)))

    return z, '', _hemi(lat), lat, lon

# **) MIT License
#
# Copyright (C) 2016-2019 -- mrJean1 at Gmail dot com
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
