
# -*- coding: utf-8 -*-

u'''Functions L{parseUTMUPS5}, L{toUtmUps8},  L{UtmUps} and
L{utmupsZoneBand5} to handle both I{Universal Transverse Mercator
(U{UTM<http://WikiPedia.org/wiki/Universal_Transverse_Mercator_coordinate_system>})}
and I{Universal Polar Stereographic
(U{UPS<http://WikiPedia.org/wiki/Universal_polar_stereographic_coordinate_system>})}
coordinates.

A pure Python implementation, partially transcribed from C++ class U{UTMUPS
<http://GeographicLib.SourceForge.io/html/classGeographicLib_1_1UTMUPS.html>}
by I{Charles Karney}.
'''

from datum import Datums
from dms import RangeError
from ellipsoidalBase import _to4lldn, _to3zBhp, \
                            _UPS_ZONE, _UPS_ZONE_STR, \
                            _UTMUPS_ZONE_MIN, _UTMUPS_ZONE_MAX
from lazily import _ALL_LAZY
from utily import OK
from ups import parseUPS5, toUps8, Ups, UPSError, upsZoneBand5
from utm import parseUTM5, toUtm8, Utm, UTMError, utmZoneBand5

# all public contants, classes and functions
__all__ = _ALL_LAZY.utmups
__version__ = '19.04.26'

_MGRS_TILE = 100e3  # block size (C{meter})

_UPS_N_MAX = 27 * _MGRS_TILE
_UPS_N_MIN = 13 * _MGRS_TILE
_UPS_S_MAX = 32 * _MGRS_TILE
_UPS_S_MIN =  8 * _MGRS_TILE

_UTM_C_MAX =   9 * _MGRS_TILE
_UTM_C_MIN =   1 * _MGRS_TILE
_UTM_N_MAX =  95 * _MGRS_TILE
_UTM_N_MIN =   0 * _MGRS_TILE
_UTM_S_MAX = 100 * _MGRS_TILE
_UTM_S_MIN =  10 * _MGRS_TILE

_UTM_N_SHIFT = _UTM_S_MAX - _UTM_N_MIN  # South minus North UTM northing


class _UpsMinMax(object):
    # UPS ranges for South, North pole
    eMax = _UPS_S_MAX, _UPS_N_MAX
    eMin = _UPS_S_MIN, _UPS_N_MIN
    nMax = _UPS_S_MAX, _UPS_N_MAX
    nMin = _UPS_S_MIN, _UPS_N_MIN


class _UtmMinMax(object):
    # UTM ranges for South-, Northern hemisphere
    eMax =  _UTM_C_MAX, _UTM_C_MAX
    eMin =  _UTM_C_MIN, _UTM_C_MIN
    nMax = (_UTM_N_MAX + _UTM_N_SHIFT), _UTM_N_MAX
    nMin =  _UTM_S_MIN, (_UTM_S_MIN - _UTM_N_SHIFT)


class UTMUPSError(ValueError):
    '''UTM/UPS parse, validate or other error.
    '''
    pass


def parseUTMUPS5(strUTMUPS, datum=Datums.WGS84, Utm=Utm, Ups=Ups, name=''):
    '''Parse a string representing a UTM or UPS coordinate, consisting
       of I{"zone[band] hemisphere/pole easting northing"}.

       @param strUTMUPS: A UTM or UPS coordinate (C{str}).
       @keyword datum: Optional datum to use (L{Datum}).
       @keyword Utm: Optional (sub-)class to return the UTM coordinate
                     (L{Utm}) or C{None}.
       @keyword Ups: Optional (sub-)class to return the UPS coordinate
                     (L{Ups}) or C{None}.
       @keyword name: Optional name (C{str}).

       @return: The UTM or UPS coordinate (L{Utm} or L{Ups}) or 5-tuple
                (C{zone, hemisphere/pole, easting, northing, Band}) if
                I{Utm} respectively I{Ups} or both are C{None} as (C{int,
                'N'|'S', meter, meter, str}) where C{zone} is C{1..60}
                for UTM or C{0} for UPS and C{Band} is C{""} or
                C{'C'|'D'..'W'|'X'} for UTM or C{'A'|'B'|'Y'|'Z'} for
                UPS.

       @raise UTMUPSError: Invalid I{strUTMUPS}.

       @see: Functions L{parseUTM5} and L{parseUPS5}.
    '''
    try:
        try:
            u = parseUTM5(strUTMUPS, datum=datum, Utm=Utm, name=name)
        except UTMError:
            u = parseUPS5(strUTMUPS, datum=datum, Ups=Ups, name=name)
    except (UTMError, UPSError):
        raise UTMUPSError('%s invalid: %r' % ('strUTMUPS', strUTMUPS))
    return u


def toUtmUps8(latlon, lon=None, datum=None, Utm=Utm, Ups=Ups, pole='', name=''):
    '''Convert a lat-/longitude point to a UTM or UPS coordinate.

       @param latlon: Latitude (C{degrees}) or an (ellipsoidal)
                      geodetic C{LatLon} point.
       @keyword lon: Optional longitude (C{degrees}) or C{None}.
       @keyword datum: Optional datum to use this UTM coordinate,
                       overriding I{latlon}'s datum (C{Datum}).
       @keyword Utm: Optional (sub-)class to return the UTM coordinate
                     (L{Utm}) or C{None}.
       @keyword Ups: Optional (sub-)class to return the UPS coordinate
                     (L{Ups}) or C{None}.
       @keyword pole: Optional top/center of UPS (stereographic)
                      projection (C{str}, C{'N[orth]'} or C{'S[outh]'}).
       @keyword name: Optional name (C{str}).

       @return: The UTM or UPS coordinate (L{Utm} respectively L{Ups})
                or an 8-tuple (C{zone, hemisphere/pole, easting, northing,
                Band, datum, convergence, scale}) if I{Utm} respectively
                I{Ups} is C{None} or I{cmoff} is C{False} as (C{int,
                'N'|'S', meter, meter, str, degrees, scalar}) where C{zone}
                is C{1..60} for UTM or C{0} for UPS and C{Band} is C{""}
                or C{'C'|'D'..'W'|'X'} for UTM or C{'A'|'B'|'Y'|'Z'} for
                UPS.

       @raise RangeError: If I{lat} outside the valid UTM or UPS bands
                          or if I{lat} or I{lon} outside the valid range
                          and I{rangerrrors} set to C{True}.

       @raise TypeError: If I{latlon} is not ellipsoidal or I{lon}
                         value is missing.

       @raise UTMUPSError: UTM or UPS validation failed.

       @raise ValueError: Invalid I{lat} or I{lon}.

       @see: Functions L{toUtm8} and L{toUps8}.
    '''
    lat, lon, d, name = _to4lldn(latlon, lon, datum, name)
    z, B, p, lat, lon = utmupsZoneBand5(lat, lon)

    if z == _UPS_ZONE:
        u = toUps8(lat, lon, datum=d, Ups=Ups, pole=pole or p, falsed=True, name=name)
    else:
        u = toUtm8(lat, lon, datum=d, Utm=Utm, cmoff=True, name=name)
    return u


def UtmUps(zone, hemipole, easting, northing, band='', datum=Datums.WGS84,
                                              falsed=True, name=''):
    '''Class-like function to create a UTM/UPS coordinate.

       @keyword zone: The UTM (longitudinal) zone with/-out Band
                      letter for UTM or for UPS zone C{"00"} or
                      C{0} (C{str} or C{int}).
       @keyword hemipole: UTM hemisphere or UPS top/center of projection
                          (C{str}, C{'N[orth]'} or C{'S[outh]'}).
       @param easting: Easting, see I{falsed} (C{meter}).
       @param northing: Northing, see I{falsed} (C{meter}).
       @keyword band: Optional, UTM (latitudinal) Band letter
                      C{'C'|'D'..'W'|'X'} or UPS (polar) Band letter
                      C{'A'|'B'|'Y'|'Z'} (C{str}).
       @keyword datum: Optional, the coordinate's datum (L{Datum}).
       @keyword falsed: Both I{easting} and I{northing} are falsed (C{bool}).
       @keyword name: Optional name (C{str}).

       @return: New UTM or UPS instance (L{Utm} or L{Ups}).

       @raise UTMUPSError: UTM or UPS validation failed.

       @see: Classes L{Utm} and L{Ups} and Karney's U{UTMUPS
             <http://GeographicLib.SourceForge.io/html/classGeographicLib_1_1UTMUPS.html>}.
    '''
    z, B, hp = _to3zBhp(zone, band=band, hemipole=hemipole)
    U = Ups if z in (_UPS_ZONE, _UPS_ZONE_STR) else Utm
    return U(z, hp, easting, northing, band=B, datum=datum, falsed=falsed, name=name)


def utmupsValidate(coord, falsed=False, MGRS=False):
    '''Check a UTM or UPS coordinate.

       @param coord: The UTM or UPS coordinate (L{Utm}, L{Ups} or C{5+Tuple}).
       @keyword falsed: C{5+Tuple} easting and northing are falsed (C{bool}).
       @keyword MGRS: Increase easting and northing ranges (C{bool}).

       @return: C{None} if validation passed.

       @raise UTMUPSError: Validation failed.

       @see: Function L{utmupsValidateOK}.
    '''

    def _en(en, lo, hi, ename):
        try:
            if lo <= float(en) <= hi:
                return
        except (TypeError, ValueError):
            pass
        t = '%s range [%.0f, %.0f]' % (U, lo, hi)
        raise UTMUPSError('%s outside %s: %g' % (ename, t, en))

    if isinstance(coord, (Ups, Utm)):
        zone = coord.zone
        hemi = coord.hemisphere
        e, n = coord.easting, coord.northing
        band = coord.band
        enMM = coord.falsed
    elif isinstance(coord, tuple) and len(coord) > 4:
        zone, hemi, e, n, band = coord[:5]
        enMM = falsed
    else:
        raise UTMUPSError('%s invalid: %r' % ('coord', coord))

    z, B, h = _to3zBhp(zone, band, hemipole=hemi)

    if z == _UPS_ZONE:  # UPS
        import ups as u  # PYCHOK expected
        U, M = 'UPS', _UpsMinMax
    else:  # UTM
        import utm as u  # PYCHOK expected
        U, M = 'UTM', _UtmMinMax

    if MGRS:
        U, s = 'MGRS', _MGRS_TILE
    else:
        s = 0

    U = '%s %s%s %s' % (U, z,B, h)

    i = 'SN'.find(h)
    if i < 0 or z < _UTMUPS_ZONE_MIN \
             or z > _UTMUPS_ZONE_MAX \
             or B not in u._Bands:
        raise UTMUPSError('%s %s, %s or %s invalid: %r' % (U,
                          'zone', 'hemisphere', 'band', (zone, hemi, band)))

    if enMM:
        _en(e, M.eMin[i] - s, M.eMax[i] + s, 'easting')   # PYCHOK .eMax .eMin
        _en(n, M.nMin[i] - s, M.nMax[i] + s, 'northing')  # PYCHOK .nMax .nMin


def utmupsValidateOK(coord, falsed=False, ok=OK):
    '''Check a UTM or UPS coordinate.

       @param coord: The UTM or UPS coordinate (L{Utm}, L{Ups} or C{5+Tuple}).
       @keyword falsed: C{5+Tuple} easting and northing are falsed (C{bool}).
       @keyword ok: Result to return if validation passed (I{OK}).

       @return: I{ok} if validation passed, the L{UTMUPSError} otherwise.

       @see: Function L{utmupsValidate}.
    '''
    try:
        utmupsValidate(coord, falsed=falsed)
        return ok
    except UTMUPSError as x:
        return x


def utmupsZoneBand5(lat, lon, cmoff=False):
    '''Return the UTM/UPS zone number, Band letter, hemisphere/pole
       and clipped lat- and longitude for a given location.

       @param lat: Latitude in degrees (C{scalar} or C{str}).
       @param lon: Longitude in degrees (C{scalar} or C{str}).
       @keyword cmoff: Offset longitude from the zone's central
                       meridian, for UTM only (C{bool}).

       @return: 5-Tuple (C{zone, Band, hemisphere/pole, lat, lon}) as
                (C{int, str, 'N'|'S', degrees90, degrees180}) where
                C{zone} is C{1..60} for UTM or C{0} for UPS and
                C{Band} is C{""} or C{'C'|'D'..'W'|'X'} for UTM or
                C{'A'|'B'|'Y'|'Z'} for UPS.

       @raise RangeError: If I{lat} outside the valid UTM or UPS bands
                          or if I{lat} or I{lon} outside the valid range
                          and I{rangerrrors} set to C{True}.

       @raise ValueError: Invalid I{lat} or I{lon}.

       @see: Functions L{utmZoneBand5} and L{upsZoneBand5}.
    '''
    try:
        return utmZoneBand5(lat, lon, cmoff=cmoff)
    except RangeError:
        return upsZoneBand5(lat, lon)

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
