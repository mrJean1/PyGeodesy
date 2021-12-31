
# -*- coding: utf-8 -*-

u'''(INTERNAL) ETM, Epsg, Mgrs, UTM and UPS bases.

Base class C{UtmUpsBase} and private functions and constants.
'''

from pygeodesy.basics import isint, isscalar, isstr, map1, neg_, \
                            _xinstanceof, _xsubclassof
from pygeodesy.datums import _ellipsoidal_datum, _WGS84
from pygeodesy.dms import degDMS, parseDMS2
from pygeodesy.ellipsoidalBase import LatLonEllipsoidalBase as _LLEB
from pygeodesy.errors import _or, ParseError, _parseX, _ValueError, _xkwds
from pygeodesy.interns import NN, _A_, _B_, _COMMA_, _float, _invalid_, \
                             _N_, _n_a_, _not_, _NS_, _PLUS_, _SPACE_, \
                             _0_0, _0_5, _180_0
from pygeodesy.lazily import _ALL_DOCS, _ALL_MODS as _MODS
from pygeodesy.named import _NamedBase, nameof, notOverloaded, _xnamed
from pygeodesy.namedTuples import EasNor2Tuple, LatLonDatum5Tuple
from pygeodesy.props import deprecated_method, property_doc_, \
                            Property_RO, property_RO
from pygeodesy.streprs import Fmt, fstr, _fstrENH2, _xattrs, _xzipairs
from pygeodesy.units import Band, Easting, Northing, Scalar, Zone
from pygeodesy.utily import wrap90, wrap360

__all__ = ()
__version__ = '21.12.28'

_MGRS_TILE =  100e3  # PYCHOK block size (C{meter})
_UPS_BANDS = _A_, _B_, 'Y', 'Z'  # UPS polar bands

_UTM_LAT_MAX  = _float( 84)          # PYCHOK for export (C{degrees})
_UTM_LAT_MIN  = _float(-80)          # PYCHOK for export (C{degrees})

_UPS_LAT_MAX  = _UTM_LAT_MAX - _0_5  # PYCHOK includes 30' UTM overlap
_UPS_LAT_MIN  = _UTM_LAT_MIN + _0_5  # PYCHOK includes 30' UTM overlap

_UTM_ZONE_MAX        =  60  # PYCHOK for export
_UTM_ZONE_MIN        =   1  # PYCHOK for export
_UTM_ZONE_OFF_MAX    =  60  # PYCHOK max Central meridian offset (C{degrees})

_UPS_ZONE            = _UTM_ZONE_MIN - 1     # PYCHOK for export
_UPS_ZONE_STR        =  Fmt.zone(_UPS_ZONE)  # PYCHOK for export

_UTMUPS_ZONE_INVALID = -4             # PYCHOK for export too
_UTMUPS_ZONE_MAX     = _UTM_ZONE_MAX  # PYCHOK for export too, by .units.py
_UTMUPS_ZONE_MIN     = _UPS_ZONE      # PYCHOK for export too, by .units.py

# _MAX_PSEUDO_ZONE      = -1
# _MIN_PSEUDO_ZONE      = -4
# _UTMUPS_ZONE_MATCH    = -3
# _UTMUPS_ZONE_STANDARD = -1
# _UTM                  = -2


def _hemi(lat, N=0):  # imported by .ups, .utm
    '''Return the hemisphere letter.

       @arg lat: Latitude (C{degrees} or C{radians}).
       @kwarg N: Minimal North latitude, C{0} or C{_N_}.

       @return: C{'N'|'S'} for north-/southern hemisphere.
    '''
    return _NS_[int(lat < N)]


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
        d = datum or _WGS84
    return lat, lon, d, (name or nameof(latlon))


def _to3zBhp(zone, band, hemipole=NN, Error=_ValueError):  # imported by .epsg, .ups, .utm, .utmups
    '''Parse UTM/UPS zone, Band letter and hemisphere/pole letter.

       @arg zone: Zone with/-out Band (C{scalar} or C{str}).
       @kwarg band: Optional I{longitudinal/polar} Band letter (C{str}).
       @kwarg hemipole: Optional hemisphere/pole letter (C{str}).
       @kwarg Error: Optional error to raise, overriding the default
                     C{ValueError}.

       @return: 3-Tuple (C{zone, Band, hemisphere/pole}) as (C{int, str,
                'N'|'S'}) where C{zone} is C{0} for UPS or C{1..60} for
                UTM and C{Band} is C{'A'..'Z'} I{NOT} checked for valid
                UTM/UPS bands.

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
            elif zone.upper() in _UPS_BANDS:  # single letter
                B = zone
                z = _UPS_ZONE

        if _UTMUPS_ZONE_MIN <= z <= _UTMUPS_ZONE_MAX:
            hp = hemipole[:1].upper()
            if hp in _NS_ or not hp:
                z = Zone(z)
                B = Band(B.upper())
                if B.isalpha():
                    return z, B, (hp or _hemi(B, _N_))
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
    x = wrap360(lon + _180_0)  # use wrap360 to get ...
    z = int(x) // 6 + 1  # ... longitudinal UTM zone [1, 60] and ...
    lon = x - _180_0  # ... lon [-180, 180) i.e. -180 <= lon < 180
    return Zone(z), wrap90(lat), lon


class UtmUpsBase(_NamedBase):
    '''(INTERNAL) Base class for L{Utm} and L{Ups} coordinates.
    '''
    _band        =  NN     # latitude band letter ('A..Z')
    _Bands       =  NN     # valid Band letters, see L{Utm} and L{Ups}
    _convergence =  None   # meridian conversion (C{degrees})
    _datum       = _WGS84  # L{Datum}
    _easting     = _0_0    # Easting, see B{C{falsed}} (C{meter})
    _Error       =  None   # I{Must be overloaded}, see function C{notOverloaded}
    _falsed      =  True   # falsed easting and northing (C{bool})
    _hemisphere  =  NN     # hemisphere ('N' or 'S'), different from UPS pole
    _latlon      =  None   # cached toLatLon (C{LatLon} or C{._toLLEB})
    _northing    = _0_0    # Northing, see B{C{falsed}} (C{meter})
    _scale       =  None   # grid or point scale factor (C{scalar}) or C{None}
#   _scale0      = _K0     # central scale factor (C{scalar})
    _ups         =  None   # cached toUps (L{Ups})
    _utm         =  None   # cached toUtm (L{Utm})

    def __init__(self, easting, northing, band=NN, datum=None, falsed=True,
                                          convergence=None, scale=None):
        '''(INTERNAL) New L{UtmUpsBase}.
        '''
        E = self._Error
        if not E:
            notOverloaded(self, callername='_Error')

        self._easting  = Easting(easting,   Error=E)
        self._northing = Northing(northing, Error=E)

        if band:
            self._band1(band)

        if datum not in (None, self._datum):
            self._datum = _ellipsoidal_datum(datum)  # raiser=True, name=band

        if not falsed:
            self._falsed = False

        if convergence is not self._convergence:
            self._convergence = Scalar(convergence=convergence, Error=E)
        if scale is not self._scale:
            self._scale = Scalar(scale=scale, Error=E)

    def __repr__(self):
        return self.toRepr(B=True)

    def __str__(self):
        return self.toStr()

    def _band1(self, band):
        '''(INTERNAL) Re/set the latitudinal or polar band.
        '''
        if band:
            _xinstanceof(str, band=band)
#           if not self._Bands:
#               notOverloaded(self, callername='_Bands')
            if band not in self._Bands:
                t = _or(*sorted(set(map(repr, self._Bands))))
                raise self._Error(band=band, txt=_not_(t))
            self._band = band
        elif self._band:  # reset
            self._band = NN

    @property_RO
    def convergence(self):
        '''Get the meridian convergence (C{degrees}) or C{None}
           if not available.
        '''
        return self._convergence

    @property_doc_(''' the (ellipsoidal) datum of this coordinate.''')
    def datum(self):
        '''Get the datum (L{Datum}).
        '''
        return self._datum

    @datum.setter  # PYCHOK setter!
    def datum(self, datum):
        '''Set the (ellipsoidal) datum L{Datum}, L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}).
        '''
        d = _ellipsoidal_datum(datum)
        if d != self.datum:
            self._update(True)
            self._datum = d

    @Property_RO
    def easting(self):
        '''Get the easting (C{meter}).
        '''
        return self._easting

    @Property_RO
    def eastingnorthing(self):
        '''Get easting and northing (L{EasNor2Tuple}C{(easting, northing)}).
        '''
        return EasNor2Tuple(self.easting, self.northing)

    def eastingnorthing2(self, falsed=True):
        '''Return easting and northing, falsed or unfalsed.

           @kwarg falsed: If C{True} return easting and northing falsed
                          (C{bool}), otherwise unfalsed.

           @return: An L{EasNor2Tuple}C{(easting, northing)} in C{meter}.
        '''
        e, n = self.falsed2
        if self.falsed and not falsed:
            e, n = neg_(e, n)
        elif falsed and not self.falsed:
            pass
        else:
            e = n = _0_0
        return EasNor2Tuple(Easting( e + self.easting,  Error=self._Error),
                            Northing(n + self.northing, Error=self._Error))

    @Property_RO
    def _epsg(self):
        '''(INTERNAL) Cache for method L{toEpsg}.
        '''
        return _MODS.epsg.Epsg(self)

    @Property_RO
    def falsed(self):
        '''Get easting and northing falsed (C{bool}).
        '''
        return self._falsed

    @Property_RO
    def falsed2(self):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}, see function C{notOverloaded}.
        '''
        notOverloaded(self)

    @property_RO
    def hemisphere(self):
        '''Get the hemisphere (C{str}, 'N'|'S').
        '''
        if not self._hemisphere:
            self._toLLEB()
        return self._hemisphere

    def _latlon5(self, LatLon, **LatLon_kwds):
        '''(INTERNAL) Get cached C{._toLLEB} as B{C{LatLon}} instance.
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

    def _latlon5args(self, ll, _toBand, unfalse, *other):
        '''(INTERNAL) See C{._toLLEB} methods, functions C{ups.toUps8} and C{utm._toXtm8}
        '''
        ll._toLLEB_args = (unfalse,) + other
        if unfalse:
            if not self._band:
                self._band = _toBand(ll.lat, ll.lon)
            if not self._hemisphere:
                self._hemisphere = _hemi(ll.lat)
        self._latlon = ll

    @Property_RO
    def northing(self):
        '''Get the northing (C{meter}).
        '''
        return self._northing

    @Property_RO
    def scale(self):
        '''Get the grid scale (C{float}) or C{None}.
        '''
        return self._scale

    @Property_RO
    def scale0(self):
        '''Get the central scale factor (C{float}).
        '''
        return self._scale0

    @deprecated_method
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
        return self._epsg

    def _toLLEB(self, **kwds):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}.
        '''
        notOverloaded(self, **kwds)

    def _toRepr(self, fmt, B, cs, prec, sep):  # PYCHOK expected
        '''(INTERNAL) Return a representation for this ETM/UTM/UPS coordinate.
        '''
        t =  self.toStr(prec=prec, sep=None, B=B, cs=cs)  # hemipole
        T = 'ZHENCS'[:len(t)]
        return _xzipairs(T, t, sep=sep, fmt=fmt)

    def _toStr(self, hemipole, B, cs, prec, sep):
        '''(INTERNAL) Return a string for this ETM/UTM/UPS coordinate.
        '''
        z = NN(Fmt.zone(self.zone), (self.band if B else NN))  # PYCHOK band
        t = (z, hemipole) + _fstrENH2(self, prec, None)[0]
        if cs:
            prec = cs if isint(cs) else 8  # for backward compatibility
            t += (_n_a_ if self.convergence is None else
                    degDMS(self.convergence, prec=prec, pos=_PLUS_),
                  _n_a_ if self.scale is None else
                      fstr(self.scale, prec=prec))
        return t if sep is None else sep.join(t)


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
    def _UTMUPS5(strUTMUPS, UPS, band, sep):
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

    return _parseX(_UTMUPS5, strUTMUPS, UPS, band, sep,
                             strUTMUPS=strUTMUPS, Error=Error)


__all__ += _ALL_DOCS(UtmUpsBase)

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
