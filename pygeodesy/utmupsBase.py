
# -*- coding: utf-8 -*-

u'''(INTERNAL) Base class C{UtmUpsBase} and private functions
for the UTM, UPS, Mgrs and Epsg classes/modules.
'''

from pygeodesy.basics import _isnotError, isscalar, isstr, \
                              issubclassof, map1, property_RO, \
                              _TypeError, _xkwds
from pygeodesy.ellipsoidalBase import LatLonEllipsoidalBase as _LLEB
from pygeodesy.datum import Datums
from pygeodesy.dms import degDMS, parseDMS2
from pygeodesy.lazily import _ALL_DOCS
from pygeodesy.named import EasNor2Tuple, LatLonDatum5Tuple, \
                           _NamedBase, nameof, notOverloaded, \
                           _xattrs, _xnamed
from pygeodesy.streprs import fstr
from pygeodesy.utily import wrap90, wrap360

__all__ = _ALL_DOCS('UtmUpsBase')
__version__ = '20.03.23'

_MGRS_TILE = 100e3  # PYCHOK block size (C{meter})

_UTM_LAT_MAX      =  84  # PYCHOK for export (C{degrees})
_UTM_LAT_MIN      = -80  # PYCHOK for export (C{degrees})
_UTM_ZONE_MAX     =  60  # PYCHOK for export
_UTM_ZONE_MIN     =   1  # PYCHOK for export
_UTM_ZONE_OFF_MAX =  60  # PYCHOK max Central meridian offset (C{degrees})

_UPS_LAT_MAX  = _UTM_LAT_MAX - 0.5     # PYCHOK includes 30' UTM overlap
_UPS_LAT_MIN  = _UTM_LAT_MIN + 0.5     # PYCHOK includes 30' UTM overlap
_UPS_ZONE     = _UTM_ZONE_MIN - 1      # PYCHOK for export
_UPS_ZONE_STR = '%02d' % (_UPS_ZONE,)  # PYCHOK for export

_UTMUPS_ZONE_INVALID = -4             # PYCHOK for export too
_UTMUPS_ZONE_MIN     = _UPS_ZONE      # PYCHOK for export too
_UTMUPS_ZONE_MAX     = _UTM_ZONE_MAX  # PYCHOK for export too

# _MAX_PSEUDO_ZONE      = -1
# _MIN_PSEUDO_ZONE      = -4
# _UTMUPS_ZONE_MATCH    = -3
# _UTMUPS_ZONE_STANDARD = -1
# _UTM                  = -2


def _hemi(lat):  # imported by .ups, .utm
    '''Return the hemisphere letter.

       @arg lat: Latitude (C{degrees} or C{radians}).

       @return: C{'N'|'S'} for north-/southern hemisphere.
    '''
    return 'S' if lat < 0 else 'N'


def _to4lldn(latlon, lon, datum, name):
    '''(INTERNAL) Return 4-tuple (C{lat, lon, datum, name}).
    '''
    try:
        # if lon is not None:
        #     raise AttributeError
        lat, lon = map1(float, latlon.lat, latlon.lon)
        _TypeError(_LLEB, LatLonDatum5Tuple, latlon=latlon)
        d = datum or latlon.datum
    except AttributeError:
        lat, lon = parseDMS2(latlon, lon)
        d = datum or Datums.WGS84
    return lat, lon, d, (name or nameof(latlon))


def _to3zBhp(zone, band, hemipole=''):  # imported by .epsg, .ups, .utm, .utmups
    '''Parse UTM/UPS zone, Band letter and hemisphere/pole letter.

       @arg zone: Zone with/-out Band (C{scalar} or C{str}).
       @kwarg band: Optional (longitudinal/polar) Band letter (C{str}).
       @kwarg hemipole: Optional hemisphere/pole letter (C{str}).

       @return: 3-Tuple (C{zone, Band, hemisphere/pole}) as (C{int,
                str, 'N'|'S'}) where C{zone} is C{0} for UPS or
                C{1..60} for UTM and C{Band} is C{'A'..'Z'} I{NOT}
                checked for valid UTM/UPS bands.

       @raise ValueError: Invalid B{C{zone}}, B{C{band}} or B{C{hemipole}}.
    '''
    B = band
    try:
        z = _UTMUPS_ZONE_INVALID
        if isscalar(zone) or zone.isdigit():
            z = int(zone)
        elif zone and isstr(zone):
            if len(zone) > 1:
                B = zone[-1:]
                z = int(zone[:-1])
            elif zone in 'AaBbYyZz':  # single letter
                B = zone
                z = _UPS_ZONE

        if _UTMUPS_ZONE_MIN <= z <= _UTMUPS_ZONE_MAX:
            hp = hemipole[:1].upper()
            if hp in ('N', 'S') or not hp:
                B = B.upper()
                if B.isalpha():
                    return z, B, (hp or ('S' if B < 'N' else 'N'))
                elif not B:
                    return z, B, hp

    except (AttributeError, TypeError, ValueError):
        pass
    raise ValueError('%s, %s or %s invalid: %r' %
                     ('zone', 'band', 'hemipole', (zone, B, hemipole)))


def _to3zll(lat, lon):  # imported by .ups, .utm
    '''Wrap lat- and longitude and determine UTM zone.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).

       @return: 3-Tuple (C{zone, lat, lon}) as (C{int}, C{degrees90},
                C{degrees180}) where C{zone} is C{1..60} for UTM.
    '''
    x = wrap360(lon + 180)  # use wrap360 to get ...
    z = int(x) // 6 + 1  # ... longitudinal UTM zone [1, 60] and ...
    lon = x - 180.0  # ... lon [-180, 180) i.e. -180 <= lon < 180
    return z, wrap90(lat), lon


class UtmUpsBase(_NamedBase):
    '''(INTERNAL) Base class for L{Utm} and L{Ups} coordinates.
    '''
    _band        = ''    #: (INTERNAL) Latitude band letter ('A..Z').
    _convergence = None  #: (INTERNAL) Meridian conversion (C{degrees}).
    _datum       = Datums.WGS84  #: (INTERNAL) L{Datum}.
    _easting     = 0     #: (INTERNAL) Easting, see B{C{falsed}} (C{meter}).
    _epsg        = None  #: (INTERNAL) toEpsg cache (L{Epsg}).
    _falsed      = True  #: (INTERNAL) Falsed easting and northing (C{bool}).
    _hemisphere  = ''    #: (INTERNAL) Hemisphere ('N' or 'S'), different from pole.
    _latlon      = None  #: (INTERNAL) toLatLon cache (C{LatLon}).
    _latlon_args = None  #: (INTERNAL) toLatLon args (varies).
    _mgrs        = None  #: (INTERNAL) toMgrs cache (L{Mgrs}.
    _northing    = 0     #: (INTERNAL) Northing, see B{C{falsed}} (C{meter}).
    _scale       = None  #: (INTERNAL) Grid scale factor (C{scalar}) or C{None}.
#   _scale0      = _K0   #: (INTERNAL) Central scale factor (C{scalar}).
    _ups         = None  #: (INTERNAL) toUps cache (L{Ups}).
    _utm         = None  #: (INTERNAL) toUtm cache (L{Utm}).

    def __repr__(self):
        return self.toStr2(B=True)

    def __str__(self):
        return self.toStr()

    @property_RO
    def convergence(self):
        '''Get the meridian convergence (C{degrees}) or C{None}.
        '''
        return self._convergence

    @property_RO
    def datum(self):
        '''Get the datum (L{Datum}).
        '''
        return self._datum

    @property_RO
    def easting(self):
        '''Get the easting (C{meter}).
        '''
        return self._easting

    @property_RO
    def eastingnorthing(self):
        '''Get easting and northing (L{EasNor2Tuple}C{(easting, northing)}) in C{meter}s.
        '''
        return EasNor2Tuple(self.easting, self.northing)

    def eastingnorthing2(self, falsed=True):
        '''Return easting and northing, falsed or unfalsed.

           @kwarg falsed: Return easting and northing falsed
                          (C{bool}), otherwise unfalsed.

           @return: An L{EasNor2Tuple}C{(easting, northing)} in C{meter}s.
        '''
        e, n = self.falsed2
        if self.falsed and not falsed:
            e, n = -e, -n
        elif falsed and not self.falsed:
            pass
        else:
            e = n = 0
        return EasNor2Tuple(e + self.easting, n + self.northing)

    @property_RO
    def falsed(self):
        '''Get easting and northing falsed (C{bool}).
        '''
        return self._falsed

    @property_RO
    def falsed2(self):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}.
        '''
        notOverloaded(self, self.falsed2.__name__)

    @property_RO
    def hemisphere(self):
        '''Get the hemisphere (C{str}, 'N'|'S').
        '''
        return self._hemisphere

    def _latlon5(self, LatLon, **LatLon_kwds):
        '''(INTERNAL) Convert cached LatLon
        '''
        ll = self._latlon
        if LatLon is None:
            r = LatLonDatum5Tuple(ll.lat, ll.lon, ll.datum,
                                  ll.convergence, ll.scale)
        elif issubclassof(LatLon, _LLEB):
            kwds = _xkwds(LatLon_kwds, datum=ll.datum)
            r = _xattrs(LatLon(ll.lat, ll.lon, **kwds),
                               ll, '_convergence', '_scale')
        else:
            raise _isnotError(_LLEB.__name__, LatLon=LatLon)
        return _xnamed(r, ll.name)

    @property_RO
    def northing(self):
        '''Get the northing (C{meter}).
        '''
        return self._northing

    @property_RO
    def scale(self):
        '''Get the grid scale (C{float}) or C{None}.
        '''
        return self._scale

    @property_RO
    def scale0(self):
        '''Get the central scale factor (C{float}).
        '''
        return self._scale0

    def to2en(self, falsed=True):
        '''DEPRECATED, use method C{eastingnorthing2}.

           @return: An L{EasNor2Tuple}C{(easting, northing)}.
        '''
        return self.eastingnorthing2(falsed=falsed)

    def toEpsg(self):
        '''Determine the B{EPSG (European Petroleum Survey Group)} code.

           @return: C{EPSG} code (C{int}).

           @raise EPSGError: See L{Epsg}.
        '''
        if self._epsg is None:
            from pygeodesy.epsg import Epsg  # PYCHOK circular import
            self._epsg = Epsg(self)
        return self._epsg

    def _toStr4_6(self, hemipole, B, cs, prec, sep):
        '''(INTERNAL) Return a string representation of this UTM/UPS coordinate.
        '''
        z = '%02d%s' % (self.zone, (self.band if B else ''))  # PYCHOK band
        t = (z, hemipole, fstr(self.easting,  prec=prec),
                          fstr(self.northing, prec=prec))
        if cs:
            t += ('n/a' if self.convergence is None else
                    degDMS(self.convergence, prec=8, pos='+'),
                  'n/a' if self.scale is None else
                      fstr(self.scale, prec=8))
        return sep.join(t)

    def _toStr2(self, prec=0, fmt='[%s]', sep=', ', B=False, cs=False):  # PYCHOK expected
        '''(INTERNAL) Return a string representation of this UTM/UPS coordinate.
        '''
        t = self.toStr(prec=prec, sep=' ', B=B, cs=cs).split()
        return fmt % (sep.join('%s:%s' % t for t in zip('ZHENCS', t)),)

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
