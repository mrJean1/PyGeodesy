
# -*- coding: utf-8 -*-

u'''Military Grid Reference System (MGRS/NATO) references.

Military Grid Reference System (MGRS/NATO) classes L{Mgrs} and L{MGRSError}
and functions L{parseMGRS} and L{toMgrs}.

Pure Python implementation of MGRS / UTM / UPS conversion functions using
an ellipsoidal earth model, in part transcoded from JavaScript originals by
I{(C) Chris Veness 2014-2016} published under the same MIT Licence**, see
U{MGRS<https://www.Movable-Type.co.UK/scripts/latlong-utm-mgrs.html>} and
U{Module mgrs<https://www.Movable-Type.co.UK/scripts/geodesy/docs/module-mgrs.html>}.

The MGRS/NATO grid references provide geocoordinate references covering the entire
globe, based on UTM and UPS projections.

MGRS references comprise a grid zone designation (GZD), a 100 km grid (square)
tile identification, and an easting and northing (in C{meter}).  The GZD consists
of a longitudinal zone (or column) number and latitudinal band (row) letter.  Each
zone (column) is 6 degrees wide and each band (row) is 8 degrees high, except the
top one 'X' is 12 degrees tall.

See also the U{United States National Grid<https://www.FGDC.gov/standards/projects/
FGDC-standards-projects/usng/fgdc_std_011_2001_usng.pdf>} and U{Military Grid
Reference System<https://WikiPedia.org/wiki/Military_grid_reference_system>}.
'''

from pygeodesy.basics import halfs2, _xinstanceof
from pygeodesy.datums import _ellipsoidal_datum, _WGS84
from pygeodesy.errors import _AssertionError, MGRSError, _parseX, \
                             _ValueError, _xkwds
from pygeodesy.interns import NN, _A_, _AtoZnoIO_, _band_, _B_, _COMMASPACE_, \
                             _datum_, _easting_, _northing_, _not_, _SPACE_, \
                             _splituple, _W_, _Y_, _Z_, _zone_, _0_, \
                             _0_5, _N_90_0
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import _NamedBase, _NamedTuple, _Pass, _xnamed
from pygeodesy.namedTuples import EasNor2Tuple, UtmUps5Tuple
from pygeodesy.props import deprecated_property_RO, Property_RO
from pygeodesy.streprs import enstr2, Fmt, _xzipairs
from pygeodesy.units import Easting, Northing, Str, _100km
from pygeodesy.units import _1um, _2000km  # PYCHOK used!
from pygeodesy.ups import _hemi, toUps8, Ups, _UPS_ZONE
from pygeodesy.utm import toUtm8, _to3zBlat, Utm, _UTM_LAT_MAX, \
                         _UTM_ZONE_MAX, _UTM_ZONE_MIN

from math import log10

__all__ = _ALL_LAZY.mgrs
__version__ = '22.07.23'

# <https://GitHub.com/hrbrmstr/mgrs/blob/master/src/mgrs.c>
_AtoPx_ = _AtoZnoIO_.tillP
_B2UPS  = {_A_: _N_90_0, _Y_: _UTM_LAT_MAX,  # UPS band bottom latitudes,
           _B_: _N_90_0, _Z_: _UTM_LAT_MAX}  # PYCHOK see .utm._to3zBlat
_FeUPS  = {_A_: 8, _B_: 20, _Y_:  8, _Z_: 20}  # falsed offsets (C{_100kms})
_FnUPS  = {_A_: 8, _B_:  8, _Y_: 13, _Z_: 13}  # falsed offsets (C{_100kms})
_JtoZx_ = 'JKLPQRSTUXYZ'  # _AtoZnoDEIMNOVW.fromJ
# 100 km grid tile UTM column (E) letters, repeating every third zone
_LeUTM  = _AtoZnoIO_.tillH, _AtoZnoIO_.fromJ.tillR, _AtoZnoIO_.fromS  # grid E colums
# 100 km grid tile UPS column (E) letters for each polar zone
_LeUPS  = {_A_: _JtoZx_, _B_: 'ABCFGHJKLPQR', _Y_: _JtoZx_, _Z_: 'ABCFGHJ'}
# 100 km grid tile UTM and UPS row (N) letters, repeating every other zone
_LnUTM  = _AtoZnoIO_.tillV, _AtoZnoIO_.fromF.tillV + _AtoZnoIO_.tillE  # grid N rows
_LnUPS  = {_A_: _AtoZnoIO_, _B_: _AtoZnoIO_, _Y_: _AtoPx_, _Z_: _AtoPx_}
_polar_ = _SPACE_('polar', _zone_)


class _RE(object):
    '''(INTERNAL) Lazily compiled regex-es.
    '''
    @Property_RO
    def MGRS_UPS(self):  # split an MGRS polar string "BEN1235..." into 3 parts
        import re  # PYCHOK warning locale.Error
        return re.compile('([ABYZ]{1})([A-Z]{2})([0-9]+)', re.IGNORECASE)

    @Property_RO
    def MGRS_UTM(self):  # split an MGRS string "12BEN1235..." into 3 parts
        import re  # PYCHOK warning locale.Error
        return re.compile('([0-9]{1,2}[C-X]{1})([A-Z]{2})([0-9]+)', re.IGNORECASE)

    @Property_RO
    def ZBEN_UPS(self):  # split an MGRS polar string "BEN" into 2 parts
        import re  # PYCHOK warning locale.Error
        return re.compile('([ABYZ]{1})([A-Z]{2})', re.IGNORECASE)

    @Property_RO
    def ZBEN_UTM(self):  # split an MGRS string "12BEN" into 2 parts
        import re  # PYCHOK warning locale.Error
        return re.compile('([0-9]{1,2}[C-X]{1})([A-Z]{2})', re.IGNORECASE)

_RE = _RE()  # PYCHOK singleton


class Mgrs(_NamedBase):
    '''Military Grid Reference System (MGRS/NATO) references,
       with method to convert to UTM coordinates.
    '''
    _band       =  NN     # latitudinal (C..X) or polar (ABYZ) band
    _bandLat    =  None   # band latitude (C{degrees90} or C{None})
    _datum      = _WGS84  # Datum (L{Datum})
    _easting    =  0      # Easting (C{meter}), within 100 km grid tile
    _EN         =  NN     # EN digraph (C{str}), 100 km grid tile
    _northing   =  0      # Northing (C{meter}), within 100 km grid tile
    _resolution =  0      # MGRS cell size (C{meter})
    _zone       =  0      # longitudinal or polar zone (C{int}), 0..60

    def __init__(self, zone, en100k, easting, northing, band=NN,
                             datum=_WGS84, resolution=0, name=NN):
        '''New L{Mgrs} Military grid reference.

           @arg zone: 6° longitudinal zone (C{int}), 1..60 covering 180°W..180°E.
           @arg en100k: Two-letter EN digraph (C{str}), 100 km grid tile using
                        I{only} the I{AA} or I{MGRS-New} (row) U{lettering scheme
                        <http://Wikipedia.org/wiki/Military_Grid_Reference_System>}.
           @arg easting: Easting (C{meter}), within 100 km grid tile.
           @arg northing: Northing (C{meter}), within 100 km grid tile.
           @kwarg band: Optional, 8° I{latitudinal} band (C{str}), 'C'|..|'W'|'X'
                        covering 80°S..84°N (without 'I' and 'O') or I{polar}
                        region 'A'|'B' at the south or 'Y'|'Z' at the north pole.
           @kwarg datum: Optional this reference's datum (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg resolution: Optional resolution, cell size (C{meter}) or C{0}.
           @kwarg name: Optional name (C{str}).

           @raise MGRSError: Invalid MGRS grid reference, B{C{zone}}, B{C{en100k}},
                             B{C{easting}}, B{C{northing}} or B{C{band}}.

           @raise TypeError: Invalid B{C{datum}}.

           @example:

            >>> from pygeodesy import Mgrs
            >>> m = Mgrs('31U', 'DQ', 48251, 11932)  # 31U DQ 48251 11932
        '''
        if name:
            self.name = name

        self._zone, self._band, self._bandLat = _to3zBlat(zone, band, Error=MGRSError)
        try:
            en = str(en100k)
            if len(en) != 2 or not en.isalpha():
                raise IndexError  # caught below
            self._EN = en.upper()
            self._EN2m()  # check E and N
        except (IndexError, TypeError, ValueError):
            raise MGRSError(en100k=en100k, zone=zone)

        self._easting  = Easting(easting,   Error=MGRSError)
        self._northing = Northing(northing, Error=MGRSError)
        if datum not in (None, self._datum):
            self._datum = _ellipsoidal_datum(datum, name=name)

        if resolution:
            self.resolution = resolution

    @Property_RO
    def band(self):
        '''Get the I{latitudinal} C{'C'|..|'W'|'X'} or I{polar}
           C{'A'|'B'|'Y'|'Z'}) band letter (C{str}).
        '''
        return self._band

    @Property_RO
    def bandLatitude(self):
        '''Get the band latitude (C{degrees90}).
        '''
        a = self._bandLat
        if a is None:
            if self.isUPS:
                self._bandLat = a = _B2UPS[self.band]
            else:  # must have been set
                raise _AssertionError(band=self.band, Lat=a)
        return a

    @Property_RO
    def datum(self):
        '''Get the datum (L{Datum}).
        '''
        return self._datum

    @Property_RO
    def EN(self):
        '''Get the 2-character grid EN digraph (C{str}).
        '''
        return self._EN

    digraph = EN

    @deprecated_property_RO
    def en100k(self):  # PYCHOK no cover
        '''DEPRECATED, use property C{EN}.'''
        return self.EN

    def _EN2m(self):
        '''(INTERNAL) Convert grid letters to east-/northing meter.
        '''
        EN = self.EN
        if self.isUTM:
            i = self.zone - 1
            # get easting from E column (note, +1 because
            # eastings start at 166e3 due to 500 km false origin)
            e = _LeUTM[i % 3].index(EN[0]) + 1
            # similarly, get northing from N row
            n = _LnUTM[i % 2].index(EN[1])
        elif self.isUPS:
            B =  self.band
            e = _LeUPS[B].index(EN[0]) + _FeUPS[B]
            n = _LnUPS[B].index(EN[1]) + _FnUPS[B]
        else:
            raise _AssertionError(zone=self.zone)
        return float(e * _100km), float(n * _100km)  # meter

    @Property_RO
    def easting(self):
        '''Gets the easting (C{meter} within grid tile).
        '''
        return self._easting

    @Property_RO
    def eastingnorthing(self):
        '''Get easting and northing (L{EasNor2Tuple}C{(easting, northing)}).
        '''
        return EasNor2Tuple(self.easting, self.northing)

    @Property_RO
    def isUPS(self):
        '''Is this MGRS in a (polar) UPS zone (C{bool}).
        '''
        return self._zone == _UPS_ZONE

    @Property_RO
    def isUTM(self):
        '''Is this MGRS in a (non-polar) UTM zone (C{bool}).
        '''
        return _UTM_ZONE_MIN <= self._zone <= _UTM_ZONE_MAX

    @Property_RO
    def northing(self):
        '''Get the northing (C{meter} within grid tile).
        '''
        return self._northing

    @Property_RO
    def northingBottom(self):
        '''Get northing of the band bottom (C{meter}), extended
           to include entirety of bottom-most 100 km tile.
        '''
        a = self.bandLatitude
        u = toUtm8(a, 0, datum=self.datum, Utm=None) if self.isUTM else \
            toUps8(a, 0, datum=self.datum, Ups=None)
        return int(u.northing / _100km) * _100km

    def parse(self, strMGRS, name=NN):
        '''Parse a string to a similar L{Mgrs} instance.

           @arg strMGRS: The MGRS reference (C{str}),
                         see function L{parseMGRS}.
           @kwarg name: Optional instance name (C{str}),
                        overriding this name.

           @return: The similar instance (L{Mgrs}).

           @raise MGRSError: Invalid B{C{strMGRS}}.
        '''
        return parseMGRS(strMGRS, datum=self.datum, Mgrs=self.classof,
                                  name=name or self.name)

    @property
    def resolution(self):
        '''Get the resolution (C{meter}).
        '''
        return self._resolution

    @resolution.setter  # PYCHOK setter!
    def resolution(self, resolution):
        '''Set the resolution of this L{Mgrs} instance, as cell size (C{meter}, power of 10).

           @raise MGRSError: Invalid B{C{resolution}}.
        '''
        try:
            if resolution and resolution > 0:
                r = int(log10(resolution))
                if r > 5:
                    raise ValueError
                r = 10**r
            else:
                r = 0
        except (ValueError, TypeError):
            raise MGRSError(resolution=resolution)
        self._resolution = r

    @Property_RO
    def tile(self):
        '''Get the MGRS grid tile size (C{meter}).
        '''
        assert _MODS.utmups._MGRS_TILE is _100km
        return _100km

    def toLatLon(self, LatLon=None, center=True, **toLatLon_kwds):
        '''Convert this MGRS grid reference to a UTM coordinate.

           @kwarg LatLon: Optional, ellipsoidal class to return the
                          geodetic point (C{LatLon}) or C{None}.
           @kwarg center: Optionally, return the grid's center or
                          lower left corner (C{bool}).
           @kwarg toLatLon_kwds: Optional, additional L{Utm.toLatLon}
                                 and B{C{LatLon}} keyword arguments.

           @return: A B{C{LatLon}} instance or if C{B{LatLon} is None}
                    a L{LatLonDatum5Tuple}C{(lat, lon, datum,
                    convergence, scale)}.

           @raise TypeError: If B{C{LatLon}} is not ellipsoidal.

           @raise UTMError: Invalid meridional radius or H-value.

           @see: Methods L{Mgrs.toUtm} and L{Utm.toLatLon}.
        '''
        u = self.toUtmUps(center=center)
        return u.toLatLon(LatLon=LatLon, **toLatLon_kwds)

    def toRepr(self, prec=10, fmt=Fmt.SQUARE, sep=_COMMASPACE_):  # PYCHOK expected
        '''Return a string representation of this MGRS grid reference.

           @kwarg prec: Number of digits (C{int}), 4:Km, 10:m.
           @kwarg fmt: Enclosing backets format (C{str}).
           @kwarg sep: Separator between name:values (C{str}).
           @return: This Mgrs as "[Z:00B, G:EN, E:meter, N:meter]" (C{str}).
        '''
        t = self.toStr(prec=prec, sep=None)
        return _xzipairs('ZGEN', t, sep=sep, fmt=fmt)

    def toStr(self, prec=10, sep=_SPACE_):  # PYCHOK expected
        '''Return a string representation of this MGRS grid reference.

           Note that MGRS grid references are truncated, not rounded
           (unlike UTM coordinates).

           @kwarg prec: Number of digits (C{int}), 4:Km, 10:m.
           @kwarg sep: Optional separator to join (C{str}) or C{None}
                       to return an unjoined C{tuple} of C{str}s.

           @return: This Mgrs as "00B EN easting northing" (C{str}).

           @raise ValueError: Invalid B{C{prec}}.

           @example:

            >>> m = Mgrs(31, 'DQ', 48251, 11932, band='U')
            >>> m.toStr()  # '31U DQ 48251 11932'
        '''
        zB = self.zoneB
        t  = enstr2(self._easting, self._northing, prec, zB, self.EN)
        return t if sep is None else sep.join(t)

    def toUps(self, Ups=Ups, center=False):
        '''Convert this MGRS grid reference to a UPS coordinate.

           @kwarg Ups: Optional class to return the UPS coordinate
                       (L{Ups}) or C{None}.
           @kwarg center: Optionally, center easting and northing
                          by the resolution (C{bool}).

           @return: A B{C{Ups}} instance or if C{B{Ups} is None}
                    a L{UtmUps5Tuple}C{(zone, hemipole, easting,
                    northing, band)}.

           @raise MGRSError: This MGRS is a I{non-polar} UTM reference.
        '''
        if self.isUTM:
            raise MGRSError(zoneB=self.zoneB, txt=_not_(_polar_))
        return self._toUtmUps(Ups, center)

    def toUtm(self, Utm=Utm, center=False):
        '''Convert this MGRS grid reference to a UTM coordinate.

           @kwarg Utm: Optional class to return the UTM coordinate
                       (L{Utm}) or C{None}.
           @kwarg center: Optionally, center easting and northing
                          by the resolution (C{bool}).

           @return: A B{C{Utm}} instance or if C{B{Utm} is None}
                    a L{UtmUps5Tuple}C{(zone, hemipole, easting,
                    northing, band)}.

           @raise MGRSError: This MGRS is a I{polar} UPS reference.
        '''
        if self.isUPS:
            raise MGRSError(zoneB=self.zoneB, txt=_polar_)
        return self._toUtmUps(Utm, center)

    def toUtmUps(self, Utm=Utm, Ups=Ups, center=False):
        '''Convert this MGRS grid reference to a UTM or UPS coordinate.

           @kwarg Utm: Optional class to return the UTM coordinate
                       (L{Utm}) or C{None}.
           @kwarg Ups: Optional class to return the UPS coordinate
                       (L{Utm}) or C{None}.
           @kwarg center: Optionally, center easting and northing
                          by the resolution (C{bool}).

           @return: A B{C{Utm}} or B{C{Ups}} instance or if C{B{Utm}
                    or B{Ups} is None} a L{UtmUps5Tuple}C{(zone,
                    hemipole, easting, northing, band)}.
        '''
        return self._toUtmUps((Utm if self.isUTM else
                              (Ups if self.isUPS else None)), center)

    def _toUtmUps(self, U, center):
        '''(INTERNAL) Helper for C{.toUps} and C{.toUtm}.
        '''
        e, n = self._EN2m()
        # 100 km grid tile row letters repeat every 2,000 km north;
        # add enough 2,000 km blocks to get into required band
        e += self.easting
        n += self.northing
        if self.isUTM:
            b = (self.northingBottom - n) / _2000km
            if b > 0:
                b  = int(b + 1)
                b  = min(b, (3 if self.band == _W_ else 4))
                n += b * _2000km
        if center:
            c = self.resolution * _0_5
            if c:
                e += c
                n += c
        z =  self.zone
        h = _hemi(self.bandLatitude)  # if self.band < _N_
        B =  self.band
        m =  self.name
        return UtmUps5Tuple(z, h, e, n, B, name=m, Error=MGRSError) if U is None \
                     else U(z, h, e, n, B, name=m, datum=self.datum)

    @Property_RO
    def zone(self):
        '''Get the I{longitudinal} zone (C{int}, 1..60 or 0 for polar).
        '''
        return self._zone

    @Property_RO
    def zoneB(self):
        '''Get the I{longitudinal} zone digits and I{latitudinal} band letter (C{str}).
        '''
        return NN(Fmt.zone(self.zone), self.band)


class Mgrs4Tuple(_NamedTuple):
    '''4-Tuple C{(zone, digraph, easting, northing)}, C{zone} and
       C{digraph} as C{str}, C{easting} and C{northing} in C{meter}.
    '''
    _Names_ = (_zone_, 'digraph', _easting_, _northing_)
    _Units_ = ( Str,    Str,       Easting,   Northing)

    def __new__(cls, z, di, e, n, Error=MGRSError, name=NN):
        if Error is not None:
            e = Easting( e, Error=Error)
            n = Northing(n, Error=Error)
        return _NamedTuple.__new__(cls, z, di, e, n, name=name)

    def to6Tuple(self, band, datum):
        '''Extend this L{Mgrs4Tuple} to a L{Mgrs6Tuple}.

           @arg band: The band to add (C{str} or C{None}).
           @arg datum: The datum to add (L{Datum} or C{None}).

           @return: An L{Mgrs6Tuple}C{(zone, digraph, easting,
                    northing, band, datum)}.
        '''
        return self._xtend(Mgrs6Tuple, band, datum)


class Mgrs6Tuple(_NamedTuple):  # XXX only used in the line above
    '''6-Tuple C{(zone, digraph, easting, northing, band, datum)},
       C{zone}, C{digraph} and C{band} as C{str}, C{easting} and
       C{northing} in C{meter} and C{datum} a L{Datum}.
    '''
    _Names_ = Mgrs4Tuple._Names_ + (_band_, _datum_)
    _Units_ = Mgrs4Tuple._Units_ + ( Str,   _Pass)


def parseMGRS(strMGRS, datum=_WGS84, Mgrs=Mgrs, name=NN):
    '''Parse a string representing a MGRS grid reference,
       consisting of C{"zoneBand, grid, easting, northing"}.

       @arg strMGRS: MGRS grid reference (C{str}).
       @kwarg datum: Optional datum to use (L{Datum}).
       @kwarg Mgrs: Optional class to return the MGRS grid
                    reference (L{Mgrs}) or C{None}.
       @kwarg name: Optional B{C{Mgrs}} name (C{str}).

       @return: The MGRS grid reference as B{C{Mgrs}} or if
                C{B{Mgrs} is None} as an L{Mgrs4Tuple}C{(zone,
                digraph, easting, northing)}.

       @raise MGRSError: Invalid B{C{strMGRS}}.

       @example:

        >>> m = parseMGRS('31U DQ 48251 11932')
        >>> str(m)  # '31U DQ 48251 11932'
        >>> m = parseMGRS('31UDQ4825111932')
        >>> repr(m)  # [Z:31U, G:DQ, E:48251, N:11932]
        >>> m = mgrs.parseMGRS('42SXD0970538646')
        >>> str(m)  # '42S XD 09705 38646'
        >>> m = mgrs.parseMGRS('42SXD9738')  # Km
        >>> str(m)  # '42S XD 97000 38000'
        >>> m = mgrs.parseMGRS('YUB17770380')  # polar
        >>> str(m)  # '00Y UB 17770 03800'
    '''
    def _mg(s, re_UTM, re_UPS):  # return re.match groups
        s = s.lstrip(_0_)
        m = re_UTM.match(s)
        if m:
            return m.groups()
        m = re_UPS.match(s)
        if m:
            m = m.groups()
            t = '00' + m[0]
            return (t,) + m[1:]
        raise ValueError

    def _s2m(g):  # e or n string to float meter
        # convert to meter if less than 5 digits
        m = g + '00000'
        return float(m[:5])

    def _MGRS(strMGRS, datum, Mgrs, name):
        m = _splituple(strMGRS)
        if len(m) == 1:  # [01]BEN1234512345'
            m = _mg(m[0], _RE.MGRS_UTM, _RE.MGRS_UPS)
            m = m[:2] + halfs2(m[2])
        elif len(m) == 2:  # [01]BEN 1234512345'
            m = _mg(m[0], _RE.ZBEN_UTM, _RE.ZBEN_UPS) + halfs2(m[1])
        elif len(m) == 3:  # [01]BEN 12345 12345'
            m = _mg(m[0], _RE.ZBEN_UTM, _RE.ZBEN_UPS) + m[1:]
        if len(m) != 4:  # [01]B EN 12345 12345
            raise ValueError
        e, n = map(_s2m, m[2:])

        zB, EN = m[0], m[1].upper()
        if Mgrs is None:
            r = Mgrs4Tuple(zB, EN, e, n, name=name)
        else:
            p = max(map(len, m[2:]))  # precision 1..5
            m = 10**(5 - p)  # resolution in meter
            r = Mgrs(zB, EN, e, n, datum=datum, resolution=m, name=name)
        return r

    return _parseX(_MGRS, strMGRS, datum, Mgrs, name,
                          strMGRS=strMGRS, Error=MGRSError)


def toMgrs(utmups, Mgrs=Mgrs, name=NN, **Mgrs_kwds):
    '''Convert a UTM or UPS coordinate to an MGRS grid reference.

       @arg utmups: A UTM or UPS coordinate (L{Utm}, L{Etm} or L{Ups]).
       @kwarg Mgrs: Optional class to return the MGRS grid reference
                    (L{Mgrs}) or C{None}.
       @kwarg name: Optional B{C{Mgrs}} name (C{str}).
       @kwarg Mgrs_kwds: Optional, additional B{C{Mgrs}} keyword
                         arguments, ignored if C{B{Mgrs} is None}.

       @return: The MGRS grid reference as B{C{Mgrs}} or if
                C{B{Mgrs} is None} as an L{Mgrs6Tuple}C{(zone,
                digraph, easting, northing, band, datum)}.

       @raise MGRSError: Invalid B{C{utmups}}.

       @raise TypeError: If B{C{utmups}} is not L{Utm} nor L{Etm}
                         nor L{Ups}.

       @example:

        >>> u = Utm(31, 'N', 448251, 5411932)
        >>> m = u.toMgrs()  # 31U DQ 48251 11932
    '''
#   _MODS.utmups.utmupsValidate(utmups, MGRS=True, Error-MGRSError)
    _xinstanceof(Utm, Ups, utmups=utmups)  # Utm, Etm, Ups
    try:
        e, n =  utmups.eastingnorthing2(falsed=True)
        E, e = _um100km2(e)
        N, n = _um100km2(n)
        B, z =  utmups.band, utmups.zone
        if _UTM_ZONE_MIN <= z <= _UTM_ZONE_MAX:
            i = z - 1
            # columns in zone 1 are A-H, zone 2 J-R, zone 3 S-Z,
            # then repeating every 3rd zone (note E-1 because
            # eastings start at 166e3 due to 500km false origin)
            EN  = _LeUTM[i % 3][E - 1]
            # rows in even zones are A-V, in odd zones are F-E
            EN += _LnUTM[i % 2][N % len(_LnUTM[0])]
        elif z == _UPS_ZONE:
            EN  = _LeUPS[B][E - _FeUPS[B]]
            EN += _LnUPS[B][N - _FnUPS[B]]
        else:
            raise _ValueError(zone=z)
    except (IndexError, TypeError, ValueError) as x:
        raise MGRSError(B=B, E=E, N=N, utmups=utmups, txt=str(x))

    if Mgrs is None:
        r = Mgrs4Tuple(Fmt.zone(z), EN, e, n).to6Tuple(B, utmups.datum)
    else:
        kwds = _xkwds(Mgrs_kwds, band=B, datum=utmups.datum)
        r = Mgrs(z, EN, e, n, **kwds)
    return _xnamed(r, name or utmups.name)


def _um100km2(m):
    '''(INTERNAL) An MGRS east-/northing truncated to micrometer (um)
       precision and to grid tile C{M} and C{m}eter within the tile.
    '''
    m = int(m / _1um) * _1um  # micrometer
    M, m = divmod(m, _100km)
    return int(M), m


if __name__ == '__main__':

    from pygeodesy.ellipsoidalVincenty import LatLon
    from pygeodesy.lazily import printf

    e = n = 0
    for lat in range(-89, 90, 1):
        for lon in range(-173, 180, 1):
            m = LatLon(lat, lon).toMgrs()
            t = m.toLatLon(LatLon=LatLon)
            d = max(abs(t.lat - lat), abs(t.lon - lon))
            if d > 1e-9:
                e += 1
                printf('%s: %s %s vs %s %.6e', e, m, t.latlon, (float(lat), float(lon)), d)
            n += 1
    p = e * 100.0 / n
    printf('%s/%s failures (%.2f%%)', (e if e else 'no'), n, p)

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
