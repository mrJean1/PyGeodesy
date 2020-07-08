
# -*- coding: utf-8 -*-

u'''(INTERNAL) Base class C{UtmUpsBase} and private functions
for the UTM, UPS, Mgrs and Epsg classes/modules.
'''

from pygeodesy.basics import isscalar, isstr, map1, property_RO, \
                            _xattrs, _xinstanceof, _xkwds, \
                            _xsubclassof, _xzipairs
from pygeodesy.ellipsoidalBase import LatLonEllipsoidalBase as _LLEB
from pygeodesy.errors import ParseError, _parseX, _ValueError
from pygeodesy.datum import Datum, Datums
from pygeodesy.dms import degDMS, parseDMS2
from pygeodesy.interns import _band_, _COMMA_, _COMMA_SPACE_, _convergence_, \
                              _datum_, _easting_, _hemipole_, _invalid_, \
                              _lat_, _lon_, _N_, _n_a_, NN, _northing_, _NS_, \
                              _PLUS_, _scale_, _SPACE_, _SQUARE_, _zone_
from pygeodesy.lazily import _ALL_DOCS
from pygeodesy.named import EasNor2Tuple,_NamedBase, _NamedTuple, nameof, \
                            notOverloaded, _xnamed
from pygeodesy.streprs import fstr
from pygeodesy.units import Band, Easting, Lat, Lon, Northing, \
                            Scalar, Zone
from pygeodesy.utily import wrap90, wrap360

__all__ = ()
__version__ = '20.07.08'

_MGRS_TILE =  100e3  # PYCHOK block size (C{meter})

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
_UTMUPS_ZONE_MAX     = _UTM_ZONE_MAX  # PYCHOK for export too, by .units.py
_UTMUPS_ZONE_MIN     = _UPS_ZONE      # PYCHOK for export too, by .units.py

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
    return _NS_[int(lat < 0)]


def _to4lldn(latlon, lon, datum, name):
    '''(INTERNAL) Return 4-tuple (C{lat, lon, datum, name}).
    '''
    try:
        # if lon is not None:
        #     raise AttributeError
        lat, lon = map1(float, latlon.lat, latlon.lon)
        _xinstanceof(_LLEB, LatLonDatum5Tuple, latlon=latlon)
        d = datum or latlon.datum
    except AttributeError:
        lat, lon = parseDMS2(latlon, lon)
        d = datum or Datums.WGS84
    return lat, lon, d, (name or nameof(latlon))


def _to3zBhp(zone, band, hemipole=NN, Error=_ValueError):  # imported by .epsg, .ups, .utm, .utmups
    '''Parse UTM/UPS zone, Band letter and hemisphere/pole letter.

       @arg zone: Zone with/-out Band (C{scalar} or C{str}).
       @kwarg band: Optional (longitudinal/polar) Band letter (C{str}).
       @kwarg hemipole: Optional hemisphere/pole letter (C{str}).
       @kwarg Error: Optional error to raise, overriding the default
                     C{ValueError}.

       @return: 3-Tuple (C{zone, Band, hemisphere/pole}) as (C{int,
                str, 'N'|'S'}) where C{zone} is C{0} for UPS or
                C{1..60} for UTM and C{Band} is C{'A'..'Z'} I{NOT}
                checked for valid UTM/UPS bands.

       @raise ValueError: Invalid B{C{zone}}, B{C{band}} or B{C{hemipole}}.
    '''
    try:
        B, z = band, _UTMUPS_ZONE_INVALID
        if isscalar(zone):
            z = int(zone)
        elif zone and isstr(zone):
            if zone.isdigit():
                z = int(zone)
            elif len(zone) > 1:
                B = zone[-1:]
                z = int(zone[:-1])
            elif zone in 'AaBbYyZz':  # single letter
                B = zone
                z = _UPS_ZONE

        if _UTMUPS_ZONE_MIN <= z <= _UTMUPS_ZONE_MAX:
            hp = hemipole[:1].upper()
            if hp in _NS_ or not hp:
                z = Zone(z)
                B = Band(B.upper())
                if B.isalpha():
                    return z, B, (hp or _NS_[B < _N_])
                elif not B:
                    return z, B, hp

        t = _invalid_
    except (AttributeError, IndexError, TypeError, ValueError) as x:
        t = str(x)  # no Python 3+ exception chaining
    raise Error(zone=zone, band=B, hemipole=hemipole, txt=t)


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
    return Zone(z), wrap90(lat), lon


class LatLonDatum5Tuple(_NamedTuple):
    '''5-Tuple C{(lat, lon, datum, convergence, scale)} in
       C{degrees90}, C{degrees180}, L{Datum}, C{degrees}
       and C{float}.
    '''
    _Names_ = (_lat_, _lon_, _datum_, _convergence_, _scale_)

    def __new__(cls, lat, lon, d, c, s):
        return _NamedTuple.__new__(cls, Lat(lat), Lon(lon), d,
                                        Scalar(c, name=_convergence_),
                                        Scalar(s, name=_scale_))


class UtmUpsBase(_NamedBase):
    '''(INTERNAL) Base class for L{Utm} and L{Ups} coordinates.
    '''
    _band        = NN    #: (INTERNAL) Latitude band letter ('A..Z').
    _convergence = None  #: (INTERNAL) Meridian conversion (C{degrees}).
    _datum       = Datums.WGS84  #: (INTERNAL) L{Datum}.
    _easting     = 0     #: (INTERNAL) Easting, see B{C{falsed}} (C{meter}).
    _Error       = None  #: (INTERNAL) I{Must be overloaded}.
    _epsg        = None  #: (INTERNAL) toEpsg cache (L{Epsg}).
    _falsed      = True  #: (INTERNAL) Falsed easting and northing (C{bool}).
    _hemisphere  = NN    #: (INTERNAL) Hemisphere ('N' or 'S'), different from pole.
    _latlon      = None  #: (INTERNAL) toLatLon cache (C{LatLon}).
    _latlon_args = None  #: (INTERNAL) toLatLon args (varies).
    _mgrs        = None  #: (INTERNAL) toMgrs cache (L{Mgrs}.
    _northing    = 0     #: (INTERNAL) Northing, see B{C{falsed}} (C{meter}).
    _scale       = None  #: (INTERNAL) Grid scale factor (C{scalar}) or C{None}.
#   _scale0      = _K0   #: (INTERNAL) Central scale factor (C{scalar}).
    _ups         = None  #: (INTERNAL) toUps cache (L{Ups}).
    _utm         = None  #: (INTERNAL) toUtm cache (L{Utm}).

    def __init__(self, easting, northing, band=NN, datum=None, falsed=True,
                                          convergence=None, scale=None):
        '''(INTERNAL) New L{UtmUpsBase}.
        '''
        E = self._Error
        if not E:
            notOverloaded(self, '_Error')

        self._easting  = Easting(easting,   Error=E)
        self._northing = Northing(northing, Error=E)

        if band:
            _xinstanceof(str, band=band)
            self._band = band

        if datum:
            _xinstanceof(Datum, datum=datum)
            if datum != self._datum:
                self._datum = datum

        if not falsed:
            self._falsed = False

        if convergence is not self._convergence:
            self._convergence = Scalar(convergence, name=_convergence_, Error=E)
        if scale is not self._scale:
            self._scale = Scalar(scale, name=_scale_, Error=E)

    def __repr__(self):
        return self.toRepr(B=True)

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
        return EasNor2Tuple(Easting( e + self.easting,  Error=self._Error),
                            Northing(n + self.northing, Error=self._Error))

    @property_RO
    def falsed(self):
        '''Get easting and northing falsed (C{bool}).
        '''
        return self._falsed

    @property_RO
    def falsed2(self):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}.
        '''
        notOverloaded(self, self.falsed2)

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
        else:
            _xsubclassof(_LLEB, LatLon=LatLon)
            kwds = _xkwds(LatLon_kwds, datum=ll.datum)
            r = _xattrs(LatLon(ll.lat, ll.lon, **kwds),
                               ll, '_convergence', '_scale')
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

    def to2en(self, falsed=True):  # PYCHOK no cover
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

    def _toRepr(self, prec=0, fmt=_SQUARE_, sep=_COMMA_SPACE_, B=False, cs=False, **unused):  # PYCHOK expected
        '''(INTERNAL) Return a string representation of this UTM/UPS coordinate.
        '''
        t = self.toStr(prec=prec, sep=None, B=B, cs=cs)
        return _xzipairs('ZHENCS', t, sep=sep, fmt=fmt)  # 'ZHENCS'[:len(t)]

    def _toStr(self, hemipole, B, cs, prec, sep):
        '''(INTERNAL) Return a string representation of this UTM/UPS coordinate.
        '''
        z = '%02d%s' % (self.zone, (self.band if B else NN))  # PYCHOK band
        t = (z, hemipole, fstr(self.easting,  prec=prec),
                          fstr(self.northing, prec=prec))
        if cs:
            t += (_n_a_ if self.convergence is None else
                    degDMS(self.convergence, prec=8, pos=_PLUS_),
                  _n_a_ if self.scale is None else
                      fstr(self.scale, prec=8))
        return t if sep is None else sep.join(t)


class UtmUps2Tuple(_NamedTuple):  # imported by .espg
    '''2-Tuple C{(zone, hemipole)} as C{int} and C{str}, where
       C{zone} is C{1..60} for UTM or C{0} for UPS and C{hemipole}
       C{'N'|'S'} is the UTM hemisphere or the UPS pole.
    '''
    _Names_ = (_zone_, _hemipole_)


class UtmUps5Tuple(_NamedTuple):  # imported by .mgrs
    '''5-Tuple C{(zone, hemipole, easting, northing, band)} as C{int},
       C{str}, C{meter}, C{meter} and C{band} letter, where C{zone}
       is C{1..60} for UTM or C{0} for UPS, C{hemipole} C{'N'|'S'} is
       the UTM hemisphere or the UPS pole and C{band} is C{""} or the
       (longitudinal) UTM band C{'C'|'D'..'W'|'X'} or the (polar) UPS
       band C{'A'|'B'|'Y'|'Z'}.
    '''
    _Names_ = (_zone_, _hemipole_, _easting_, _northing_, _band_)

    def __new__(cls, z, h, e, n, B, Error=None):
        if Error is not None:
            e = Easting( e, Error=Error)
            n = Northing(n, Error=Error)
        return _NamedTuple.__new__(cls, z, h, e, n, B)


class UtmUps8Tuple(_NamedTuple):
    '''8-Tuple C{(zone, hemipole, easting, northing, band, datum,
       convergence, scale)} as C{int}, C{str}, C{meter}, C{meter},
       C{band} letter, C{Datum}, C{degrees} and C{float}, where
       C{zone} is C{1..60} for UTM or C{0} for UPS, C{hemipole}
       C{'N'|'S'} is the UTM hemisphere or the UPS pole and C{band}
       is C{""} or the (longitudinal) UTM band C{'C'|'D'..'W'|'X'}
       or the (polar) UPS band C{'A'|'B'|'Y'|'Z'}.
    '''
    _Names_ = (_zone_, _hemipole_, _easting_, _northing_,
               _band_, _datum_, _convergence_, _scale_)

    def __new__(cls, z, h, e, n, B, d, c, s, Error=None):
        if Error is not None:
            e = Easting( e, Error=Error)
            n = Northing(n, Error=Error)
            c = Scalar(c, name=_convergence_, Error=Error)
            s = Scalar(s, name=_scale_, Error=Error)
        return _NamedTuple.__new__(cls, z, h, e, n, B, d, c, s)


class UtmUpsLatLon5Tuple(_NamedTuple):
    '''5-Tuple C{(zone, band, hemipole, lat, lon)} as C{int},
       C{str}, C{str}, C{degrees90} and C{degrees180}, where
       C{zone} is C{1..60} for UTM or C{0} for UPS, C{band} is
       C{""} or the (longitudinal) UTM band C{'C'|'D'..'W'|'X'}
       or (polar) UPS band C{'A'|'B'|'Y'|'Z'} and C{hemipole}
       C{'N'|'S'} is the UTM hemisphere or the UPS pole.
    '''
    _Names_ = (_zone_, _band_, _hemipole_, _lat_, _lon_)

    def __new__(cls, z, B, h, lat, lon, Error=None):
        if Error is not None:
            lat = Lat(lat, Error=Error)
            lon = Lon(lon, Error=Error)
        return _NamedTuple.__new__(cls, z, B, h, lat, lon)


def _parseUTMUPS5(strUTMUPS, UPS, Error=ParseError, band=NN, sep=_COMMA_):
    '''(INTERNAL) Parse a string representing a UTM or UPS coordinate
       consisting of C{"zone[band] hemisphere/pole easting northing"}.

       @arg strUTMUPS: A UTM or UPS coordinate (C{str}).
       @kwarg band: Optional, default Band letter (C{str}).
       @kwarg sep: Optional, separator to split (",").

       @return: 5-Tuple (C{zone, hemisphere/pole, easting, northing,
                band}).

       @raise ParseError: Invalid B{C{strUTMUPS}}.
    '''
    def _UTMUPS5_(strUTMUPS, UPS, band, sep):
        u = strUTMUPS.lstrip()
        if UPS and not u.startswith(_UPS_ZONE_STR):
            raise ValueError

        u = u.replace(sep, _SPACE_).strip().split()
        if len(u) < 4:
            raise ValueError

        z, h = u[:2]
        if h[:1].upper() not in _NS_:
            raise ValueError

        if z.isdigit():
            z, B = int(z), band
        else:
            for i in range(len(z)):
                if not z[i].isdigit():
                    # int('') raises ValueError
                    z, B = int(z[:i]), z[i:]
                    break
            else:
                raise ValueError

        e, n = map(float, u[2:4])
        return z, h.upper(), e, n, B.upper()

    return _parseX(_UTMUPS5_, strUTMUPS, UPS, band, sep,
                              strUTMUPS=strUTMUPS, Error=Error)


__all__ += _ALL_DOCS(UtmUpsBase, LatLonDatum5Tuple,
                     UtmUps2Tuple,
                     UtmUps5Tuple,
                     UtmUps8Tuple,
                     UtmUpsLatLon5Tuple)

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
