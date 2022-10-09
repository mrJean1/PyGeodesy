
# -*- coding: utf-8 -*-

u'''(INTERNAL) Private class C{UtmUpsBase}, functions and constants
for L{epsg}, L{etm}, L{mgrs}, L{ups} and L{utm}.
'''

from pygeodesy.basics import isint, isscalar, isstr, map1, neg_, \
                            _xinstanceof, _xsubclassof
from pygeodesy.constants import _float, _0_0, _0_5, _N_90_0, _180_0
from pygeodesy.datums import _ellipsoidal_datum, _WGS84
from pygeodesy.dms import degDMS, parseDMS2
from pygeodesy.ellipsoidalBase import LatLonEllipsoidalBase as _LLEB
from pygeodesy.errors import _or, ParseError, _parseX, _ValueError, \
                             _xkwds, _xkwds_get, _xkwds_not
from pygeodesy.interns import NN, _A_, _B_, _COMMA_, _Error_, \
                             _gamma_, _n_a_, _not_, _N_, _NS_, _PLUS_, \
                             _scale_, _SPACE_, _Y_, _Z_, _UNDER
from pygeodesy.lazily import _ALL_DOCS, _ALL_MODS as _MODS
from pygeodesy.named import _NamedBase, nameof, notOverloaded, _xnamed
from pygeodesy.namedTuples import EasNor2Tuple, LatLonDatum5Tuple
from pygeodesy.props import deprecated_method, property_doc_, _update_all, \
                            deprecated_property_RO, Property_RO, property_RO
from pygeodesy.streprs import Fmt, fstr, _fstrENH2, _xattrs, _xzipairs
from pygeodesy.units import Band, Easting, Northing, Scalar, Zone
from pygeodesy.utily import wrap360

__all__ = ()
__version__ = '22.10.05'

_UPS_BANDS = _A_, _B_, _Y_, _Z_  # UPS polar bands SE, SW, NE, NW
# _UTM_BANDS = _MODS.utm._Bands

_UTM_LAT_MAX  = _float( 84)          # PYCHOK for export (C{degrees})
_UTM_LAT_MIN  = _float(-80)          # PYCHOK for export (C{degrees})

_UPS_LAT_MAX  = _UTM_LAT_MAX - _0_5  # PYCHOK includes 30' UTM overlap
_UPS_LAT_MIN  = _UTM_LAT_MIN + _0_5  # PYCHOK includes 30' UTM overlap

_UPS_LATS = {_A_: _N_90_0, _Y_: _UTM_LAT_MAX,  # UPS band bottom latitudes,
             _B_: _N_90_0, _Z_: _UTM_LAT_MAX}  # PYCHOK see .Mgrs.bandLatitude

_UTM_ZONE_MAX     =  60  # PYCHOK for export
_UTM_ZONE_MIN     =   1  # PYCHOK for export
_UTM_ZONE_OFF_MAX =  60  # PYCHOK max Central meridian offset (C{degrees})

_UPS_ZONE     = _UTM_ZONE_MIN - 1     # PYCHOK for export
_UPS_ZONE_STR =  Fmt.zone(_UPS_ZONE)  # PYCHOK for export

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
                B =  zone
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

        raise ValueError  # _invalid_
    except (AttributeError, IndexError, TypeError, ValueError) as x:
        raise Error(zone=zone, band=B, hemipole=hemipole, cause=x)


def _to3zll(lat, lon):  # imported by .ups, .utm
    '''Wrap lat- and longitude and determine UTM zone.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).

       @return: 3-Tuple (C{zone, lat, lon}) as (C{int}, C{degrees90},
                C{degrees180}) where C{zone} is C{1..60} for UTM.
    '''
    x = wrap360(lon + _180_0)  # use wrap360 to get ...
    z = int(x) // 6 + 1  # ... longitudinal UTM zone [1, 60] and ...
    return Zone(z), lat, (x - _180_0)  # ... -180 <= lon < 180


class UtmUpsBase(_NamedBase):
    '''(INTERNAL) Base class for L{Utm} and L{Ups} coordinates.
    '''
    _band       =  NN     # latitude band letter ('A..Z')
    _Bands      =  NN     # valid Band letters, see L{Utm} and L{Ups}
    _datum      = _WGS84  # L{Datum}
    _easting    = _0_0    # Easting, see B{C{falsed}} (C{meter})
    _Error      =  None   # I{Must be overloaded}, see function C{notOverloaded}
    _falsed     =  True   # falsed easting and northing (C{bool})
    _gamma      =  None   # meridian conversion (C{degrees})
    _hemisphere =  NN     # hemisphere ('N' or 'S'), different from UPS pole
    _latlon     =  None   # cached toLatLon (C{LatLon} or C{._toLLEB})
    _northing   = _0_0    # Northing, see B{C{falsed}} (C{meter})
    _scale      =  None   # grid or point scale factor (C{scalar}) or C{None}
#   _scale0     = _K0     # central scale factor (C{scalar})
    _ups        =  None   # cached toUps (L{Ups})
    _utm        =  None   # cached toUtm (L{Utm})

    def __init__(self, easting, northing, band=NN, datum=None, falsed=True,
                                          gamma=None, scale=None, **convergence):
        '''(INTERNAL) New L{UtmUpsBase}.
        '''
        E = self._Error
        if not E:
            notOverloaded(self, callername=_UNDER(_Error_))

        self._easting  = Easting(easting,   Error=E)
        self._northing = Northing(northing, Error=E)

        if band:
            self._band1(band)

        if datum not in (None, self._datum):
            self._datum = _ellipsoidal_datum(datum)  # raiser=_datum_, name=band

        if not falsed:
            self._falsed = False

        if convergence:  # for backward compatibility
            gamma = _xkwds_get(convergence, convergence=gamma)
        if gamma is not self._gamma:
            self._gamma = Scalar(gamma=gamma, Error=E)
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
#               notOverloaded(self, callername=_UNDER('Bands'))
            if band not in self._Bands:
                t = _or(*sorted(set(map(repr, self._Bands))))
                raise self._Error(band=band, txt=_not_(t))
            self._band = band
        elif self._band:  # reset
            self._band = NN

    @deprecated_property_RO
    def convergence(self):
        '''DEPRECATED, use property C{gamma}.'''
        return self.gamma

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
        if self._datum != d:
            _update_all(self)
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

    @Property_RO
    def gamma(self):
        '''Get the meridian convergence (C{degrees}) or C{None}
           if not available.
        '''
        return self._gamma

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
                                  ll.gamma, ll.scale)
        else:
            _xsubclassof(_LLEB, LatLon=LatLon)
            kwds = _xkwds(LatLon_kwds, datum=ll.datum)
            r = _xattrs(LatLon(ll.lat, ll.lon, **kwds),
                               ll, _UNDER(_gamma_), _UNDER(_scale_))
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
    def _lowerleft(self):  # by .ellipsoidalBase._lowerleft
        '''Get this UTM or UPS C{un}-centered (L{Utm} or L{Ups}) to its C{lowerleft}.
        '''
        return _lowerleft(self, 0)

    @Property_RO
    def _mgrs(self):
        '''(INTERNAL) Cache for method L{toMgrs}.
        '''
        return _toMgrs(self)

    @Property_RO
    def _mgrs_lowerleft(self):
        '''(INTERNAL) Cache for method L{toMgrs}, I{un}-centered.
        '''
        utmups = self._lowerleft
        return self._mgrs if utmups is self else _toMgrs(utmups)

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

    def toMgrs(self, center=False):
        '''Convert this UTM/UPS coordinate to an MGRS grid reference.

           @kwarg center: If C{True}, I{un}-center this UTM or UPS to
                          its C{lowerleft} (C{bool}) or by C{B{center}
                          meter} (C{scalar}).

           @return: The MGRS grid reference (L{Mgrs}).

           @see: Function L{pygeodesy.toMgrs} in module L{mgrs} for more details.

           @note: If not specified, the I{latitudinal} C{band} is computed from
                  the (geodetic) latitude and the C{datum}.
        '''
        return self._mgrs if center in (False, 0, _0_0) else (
               self._mgrs_lowerleft if center in (True,) else
              _toMgrs(_lowerleft(self, center)))  # PYCHOK indent

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
            t += (_n_a_ if self.gamma is None else
                    degDMS(self.gamma, prec=prec, pos=_PLUS_),
                  _n_a_ if self.scale is None else
                      fstr(self.scale, prec=prec))
        return t if sep is None else sep.join(t)


def _lowerleft(utmups, center):  # by .ellipsoidalBase._lowerleft
    '''(INTERNAL) I{Un}-center a B{C{utmups}} to its C{lowerleft} by
       C{B{center} meter} or by a I{guess} if B{C{center}} is C{0}.
    '''
    if center:
        e = n = -center
    else:
        c = 5  # center
        for _ in range(3):
            c *= 10  # 50, 500, 5000
            t  = c * 2
            e = int(utmups.easting  % t)
            n = int(utmups.northing % t)
            if (e == c and n in (c, c - 1)) or \
               (n == c and e in (c, c - 1)):
                break
        else:
            return utmups  # unchanged

    r = _xkwds_not(None, datum=utmups.datum,
                         gamma=utmups.gamma,
                         scale=utmups.scale)
    return utmups.classof(utmups.zone, utmups.hemisphere,
                          utmups.easting - e, utmups.northing - n,
                          band=utmups.band, falsed=utmups.falsed, **r)


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
            raise ValueError(_not_(_UPS_ZONE_STR))

        u = u.replace(sep, _SPACE_).strip().split()
        if len(u) < 4:
            raise ValueError(_not_(sep))

        z, h = u[:2]
        if h[:1].upper() not in _NS_:
            raise ValueError(_SPACE_(h, _not_(_NS_)))

        if z.isdigit():
            z, B = int(z), band
        else:
            for i in range(len(z)):
                if not z[i].isdigit():
                    # int('') raises ValueError
                    z, B = int(z[:i]), z[i:]
                    break
            else:
                raise ValueError(z)

        e, n = map(float, u[2:4])
        return z, h.upper(), e, n, B.upper()

    return _parseX(_UTMUPS5, strUTMUPS, UPS, band, sep,
                             strUTMUPS=strUTMUPS, Error=Error)


def _toMgrs(utmups):
    '''(INTERNAL) Convert a L{Utm} or L{Ups} to an L{Mgrs} instance.
    '''
    return _MODS.mgrs.toMgrs(utmups, datum=utmups.datum, name=utmups.name)


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
