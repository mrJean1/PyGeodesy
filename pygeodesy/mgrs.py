
# -*- coding: utf-8 -*-

u'''Military Grid Reference System (MGRS/NATO) references.

Classes L{Mgrs}, L{Mgrs4Tuple} and L{Mgrs6Tuple} and functions L{parseMGRS}
and L{toMgrs}.

Pure Python implementation of MGRS, UTM and UPS conversions covering the entire
I{ellipsoidal} earth, transcoded from I{Chris Veness}' JavaScript originals U{MGRS
<https://www.Movable-Type.co.UK/scripts/latlong-utm-mgrs.html>} and U{Module mgrs
<https://www.Movable-Type.co.UK/scripts/geodesy/docs/module-mgrs.html>} and from
I{Charles Karney}'s C++ class U{MGRS<https://GeographicLib.SourceForge.io/C++/doc/
classGeographicLib_1_1MGRS.html>}.

MGRS references comprise a grid zone designation (GZD), a 100 Km grid (square)
tile identification and an easting and northing (in C{meter}).  The GZD consists
of a longitudinal zone (or column) I{number} and latitudinal band (row) I{letter}
in the UTM region between 80°S and 84°N.  Each zone (column) is 6° wide and each
band (row) is 8° high, except top band 'X' is 12° tall.  In UPS polar regions
below 80°S and above 84°N the GZD contains only a single I{letter}, C{'A'} or
C{'B'} near the south and C{'Y'} or C{'Z'} around the north pole (for west
respectively east longitudes).

See also the U{United States National Grid<https://www.FGDC.gov/standards/projects/
FGDC-standards-projects/usng/fgdc_std_011_2001_usng.pdf>} and U{Military Grid
Reference System<https://WikiPedia.org/wiki/Military_grid_reference_system>}.

See module L{pygeodesy.ups} for env variable C{PYGEODESY_UPS_POLES} determining
the UPS encoding I{at} the south and north pole.

Set env variable C{PYGEODESY_GEOCONVERT} to the (fully qualified) path of the
C{GeoConvert} executable to run this module as I{python[3] -m pygeodesy.mgrs}
and compare the MGRS results with those from I{Karney}'s utility U{GeoConvert
<https://GeographicLib.sourceforge.io/C++/doc/GeoConvert.1.html>}.
'''

from pygeodesy.basics import halfs2, _splituple, _xinstanceof
# from pygeodesy.constants import _0_5  # from .units
from pygeodesy.datums import _ellipsoidal_datum, _WGS84
from pygeodesy.errors import _AssertionError, MGRSError, _parseX, \
                             _ValueError, _xkwds
from pygeodesy.interns import NN, _0_, _A_, _AtoZnoIO_, _band_, _B_, \
                             _COMMASPACE_, _datum_, _easting_, _invalid_, \
                             _northing_, _not_, _SPACE_, _W_, _Y_, _Z_, _zone_
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS, _PYGEODESY_GEOCONVERT_
from pygeodesy.named import _NamedBase, _NamedTuple, _Pass, _xnamed
from pygeodesy.namedTuples import EasNor2Tuple, UtmUps5Tuple
from pygeodesy.props import deprecated_property_RO, property_RO, Property_RO
from pygeodesy.streprs import enstr2, _enstr2m3, Fmt, _resolution10, _xzipairs
from pygeodesy.units import Easting, Northing, Str, _100km,  _0_5
from pygeodesy.units import _1um, _2000km  # PYCHOK used!
from pygeodesy.ups import _hemi, toUps8, Ups, _UPS_ZONE
from pygeodesy.utm import toUtm8, _to3zBlat, Utm, _UTM_ZONE_MAX, _UTM_ZONE_MIN
# from pygeodesy.utmupsBase import _UTM_ZONE_MAX, _UTM_ZONE_MIN  # from .utm

__all__ = _ALL_LAZY.mgrs
__version__ = '23.09.14'

_AN_    = 'AN'  # default south pole grid tile and band B
_AtoPx_ = _AtoZnoIO_.tillP
# <https://GitHub.com/hrbrmstr/mgrs/blob/master/src/mgrs.c>
_FeUPS  = {_A_: 8, _B_: 20, _Y_:  8, _Z_: 20}  # falsed offsets (C{_100kms})
_FnUPS  = {_A_: 8, _B_:  8, _Y_: 13, _Z_: 13}  # falsed offsets (C{_100kms})
_JtoZx_ = 'JKLPQRSTUXYZZ'  # _AtoZnoDEIMNOVW.fromJ, duplicate Z
# 100 Km grid tile UTM column (E) letters, repeating every third zone
_LeUTM  = _AtoZnoIO_.tillH, _AtoZnoIO_.fromJ.tillR, _AtoZnoIO_.fromS  # grid E colums
# 100 Km grid tile UPS column (E) letters for each polar zone
_LeUPS  = {_A_: _JtoZx_, _B_: 'ABCFGHJKLPQR', _Y_: _JtoZx_, _Z_: 'ABCFGHJ'}
# 100 Km grid tile UTM and UPS row (N) letters, repeating every other zone
_LnUTM  = _AtoZnoIO_.tillV, _AtoZnoIO_.fromF.tillV + _AtoZnoIO_.tillE  # grid N rows
_LnUPS  = {_A_: _AtoZnoIO_, _B_: _AtoZnoIO_, _Y_: _AtoPx_, _Z_: _AtoPx_}
_polar_ = _SPACE_('polar', _zone_)


class Mgrs(_NamedBase):
    '''Military Grid Reference System (MGRS/NATO) references,
       with method to convert to UTM coordinates.
    '''
    _band       =  NN     # latitudinal (C..X) or polar (ABYZ) band
    _bandLat    =  None   # band latitude (C{degrees90} or C{None})
    _datum      = _WGS84  # Datum (L{Datum})
    _easting    =  0      # Easting (C{meter}), within 100 Km grid tile
    _EN         =  NN     # EN digraph (C{str}), 100 Km grid tile
    _northing   =  0      # Northing (C{meter}), within 100 Km grid tile
    _resolution =  0      # from L{parseMGRS}, centering (C{meter})
    _zone       =  0      # longitudinal or polar zone (C{int}), 0..60

    def __init__(self, zone=0, EN=NN, easting=0, northing=0, band=NN,
                               datum=_WGS84, resolution=0, name=NN):
        '''New L{Mgrs} Military grid reference.

           @arg zone: The 6° I{longitudinal} zone (C{int}), 1..60 covering
                      180°W..180°E or C{0} for I{polar} regions or (C{str})
                      with the zone number and I{latitudinal} band letter.
           @arg EN: Two-letter EN digraph (C{str}), grid tile I{using only}
                    the I{AA} aka I{MGRS-New} (row) U{lettering scheme
                    <http://Wikipedia.org/wiki/Military_Grid_Reference_System>}.
           @kwarg easting: Easting (C{meter}), within 100 Km grid tile.
           @kwarg northing: Northing (C{meter}), within 100 Km grid tile.
           @kwarg band: Optional, I{latitudinal} band or I{polar} region letter
                        (C{str}), 'C'|..|'X' covering 80°S..84°N (no 'I'|'O'),
                        'A'|'B' at the south or 'Y'|'Z' at the north pole.
           @kwarg datum: This reference's datum (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
           @kwarg resolution: Optional resolution (C{meter}), C{0} for default.
           @kwarg name: Optional name (C{str}).

           @raise MGRSError: Invalid B{C{zone}}, B{C{EN}}, B{C{easting}},
                             B{C{northing}}, B{C{band}} or B{C{resolution}}.

           @raise TypeError: Invalid B{C{datum}}.

           @example:

            >>> from pygeodesy import Mgrs
            >>> m = Mgrs('31U', 'DQ', 48251, 11932)  # 31U DQ 48251 11932
            >>> m = Mgrs()  # defaults to south pole
            >>> m.toLatLon()  # ... lat=-90.0, lon=0.0, datum=...
        '''
        if name:
            self.name = name

        if not (zone or EN or band):
            EN, band = _AN_, _B_  # default, south pole
        try:
            self._zone, self._band, self._bandLat = _to3zBlat(zone, band, Error=MGRSError)
            en = str(EN)
            if len(en) != 2 or not en.isalpha():
                raise ValueError  # caught below
            self._EN = en.upper()
            _ = self._EN2m  # check E and N
        except (IndexError, KeyError, TypeError, ValueError):
            raise MGRSError(band=band, EN=EN, zone=zone)

        self._easting  = Easting(easting,   Error=MGRSError)
        self._northing = Northing(northing, Error=MGRSError)
        if datum not in (None, Mgrs._datum):
            self._datum = _ellipsoidal_datum(datum, name=name)  # XXX raiser=_datum_

        if resolution:
            self.resolution = resolution

    def __str__(self):
        return self.toStr(sep=_SPACE_)  # for backward compatibility

    @property_RO
    def band(self):
        '''Get the I{latitudinal} band C{'C'|..|'X'} (no C{'I'|'O'})
           or I{polar} region C{'A'|'B'|'Y'|'Z'}) letter (C{str}).
        '''
        return self._band

    @Property_RO
    def bandLatitude(self):
        '''Get the band latitude (C{degrees90}).
        '''
        return self._bandLat

    @Property_RO
    def datum(self):
        '''Get the datum (L{Datum}).
        '''
        return self._datum

    @deprecated_property_RO
    def digraph(self):
        '''DEPRECATED, use property C{EN}.'''
        return self.EN

    @property_RO
    def EN(self):
        '''Get the 2-letter grid tile (C{str}).
        '''
        return self._EN

    @deprecated_property_RO
    def en100k(self):
        '''DEPRECATED, use property C{EN}.'''
        return self.EN

    @Property_RO
    def _EN2m(self):
        '''(INTERNAL) Get the grid 2-tuple (easting, northing) in C{meter}.

           @note: Raises AssertionError, IndexError or KeyError: Invalid
                  C{zone} number, C{EN} letter or I{polar} region letter.
        '''
        EN = self.EN
        if self.isUTM:
            i = self.zone - 1
            # get easting from the E column (note, +1 because
            # easting starts at 166e3 due to 500 Km falsing)
            e = _LeUTM[i % 3].index(EN[0]) + 1
            # similarly, get northing from the N row
            n = _LnUTM[i % 2].index(EN[1])
        elif self.isUPS:
            B =  self.band
            e = _LeUPS[B].index(EN[0]) + _FeUPS[B]
            n = _LnUPS[B].index(EN[1]) + _FnUPS[B]
        else:
            raise _AssertionError(zone=self.zone)
        return float(e * _100km), float(n * _100km)  # meter

    @property_RO
    def easting(self):
        '''Get the easting (C{meter} within grid tile).
        '''
        return self._easting

    @Property_RO
    def eastingnorthing(self):
        '''Get easting and northing (L{EasNor2Tuple}C{(easting, northing)})
           I{within} the MGRS grid tile, both in C{meter}.
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

    @property_RO
    def northing(self):
        '''Get the northing (C{meter} within grid tile).
        '''
        return self._northing

    @Property_RO
    def northingBottom(self):
        '''Get the northing of the band bottom (C{meter}).
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
        '''Get the MGRS resolution (C{meter}, power of 10)
           or C{0} if undefined.
        '''
        return self._resolution

    @resolution.setter  # PYCHOK setter!
    def resolution(self, resolution):
        '''Set the MGRS resolution (C{meter}, power of 10)
           or C{0} to undefine and disable UPS/UTM centering.

           @raise MGRSError: Invalid B{C{resolution}}, over
                             C{1.e+5} or under C{1.e-6}.
        '''
        if resolution:  # and resolution > 0
            r = _resolution10(resolution, Error=MGRSError)
        else:
            r = 0
        if self._resolution != r:
            self._resolution = r

    @Property_RO
    def tilesize(self):
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
                    a L{LatLonDatum5Tuple}C{(lat, lon, datum, gamma,
                    scale)}.

           @raise TypeError: If B{C{LatLon}} is not ellipsoidal.

           @raise UTMError: Invalid meridional radius or H-value.

           @see: Methods L{Mgrs.toUtm} and L{Utm.toLatLon}.
        '''
        u = self.toUtmUps(center=center)
        return u.toLatLon(LatLon=LatLon, **toLatLon_kwds)

    def toRepr(self, fmt=Fmt.SQUARE, sep=_COMMASPACE_, **prec):  # PYCHOK expected
        '''Return a string representation of this MGRS grid reference.

           @kwarg fmt: Enclosing backets format (C{str}).
           @kwarg sep: Separator between name:values (C{str}).
           @kwarg prec: Precision (C{int}), see method L{Mgrs.toStr}.

           @return: This Mgrs as "[Z:[dd]B, G:EN, E:easting, N:northing]"
                    (C{str}), with C{B{sep} ", "}.

           @note: MGRS grid references are truncated, not rounded (unlike
                  UTM/UPS coordinates).

           @raise ValueError: Invalid B{C{prec}}.
        '''
        t = self.toStr(sep=None, **prec)
        return _xzipairs('ZGEN', t, sep=sep, fmt=fmt)

    def toStr(self, prec=0, sep=NN):  # PYCHOK expected
        '''Return this MGRS grid reference as a string.

           @kwarg prec: Precision, the number of I{decimal} digits (C{int}) or if
                        negative, the number of I{units to drop}, like MGRS U{PRECISION
                        <https://GeographicLib.SourceForge.io/C++/doc/GeoConvert.1.html#PRECISION>}.
           @kwarg sep: Optional separator to join (C{str}) or C{None} to return an unjoined
                       3-C{tuple} of C{str}s.

           @return: This Mgrs as 4-tuple C{("dd]B", "EN", "easting", "northing")} if C{B{sep}=NN}
                    or "[dd]B EN easting northing" (C{str}) with C{B{sep} " "}.

           @note: Both C{easting} and C{northing} strings are C{NN} or missing if C{B{prec} <= -5}.

           @note: MGRS grid references are truncated, not rounded (unlike UTM/UPS).

           @raise ValueError: Invalid B{C{prec}}.

           @example:

            >>> from pygeodesy import Mgrs, NN, parseMGRS
            >>> m = Mgrs(31, 'DQ', 48251, 11932, band='U')
            >>> m.toStr()  # '31U DQ 48251 11932'
            >>> m = parseMGRS('BAN1234567890')
            >>> str(m)  # 'B AN 12345 67890'
            >>> m.toStr()  # 'BAN1234567890'
            >>> m.toStr(prec=-2)  # 'BAN123678'
        '''
        zB = self.zoneB
        t  = enstr2(self._easting, self._northing, prec, zB, self.EN)
        return t if sep is None else sep.join(t).rstrip()

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
        e, n = self._EN2m
        e += self.easting
        n += self.northing
        if self.isUTM:
            # 100 Km row letters repeat every 2,000 Km north;
            # add 2,000 Km blocks to get into required band
            b = (self.northingBottom - n) / _2000km
            if b > 0:
                b  = int(b) + 1
                b  = min(b, (3 if self.band == _W_ else 4))
                n += b * _2000km
        if center:
            c = self.resolution
            if c:
                c *= _0_5
                e +=  c
                n +=  c
        z =  self.zone
        h = _hemi(self.bandLatitude)  # _S_ if self.band < _N_ else _N_
        B =  self.band
        m =  self.name
        return UtmUps5Tuple(z, h, e, n, B, name=m, Error=MGRSError) if U is None \
                     else U(z, h, e, n, B, name=m, datum=self.datum)

    @property_RO
    def zone(self):
        '''Get the I{longitudinal} zone (C{int}), 1..60 or 0 for I{polar}.
        '''
        return self._zone

    @Property_RO
    def zoneB(self):
        '''Get the I{polar} region letter or the I{longitudinal} zone digits
           plus I{latitudinal} band letter (C{str}).
        '''
        return self.band if self.isUPS else NN(Fmt.zone(self.zone), self.band)


class Mgrs4Tuple(_NamedTuple):
    '''4-Tuple C{(zone, EN, easting, northing)}, C{zone} and grid
       tile C{EN} as C{str}, C{easting} and C{northing} in C{meter}.

       @note: The C{zone} consists of either the I{longitudinal} zone
              number plus the I{latitudinal} band letter or only the
              I{polar} region letter.
    '''
    _Names_ = (_zone_, 'EN', _easting_, _northing_)
    _Units_ = ( Str,    Str,  Easting,   Northing)

    @deprecated_property_RO
    def digraph(self):
        '''DEPRECATED, use attribute C{EN}.'''
        return self.EN  # PYCHOK or [1]

    def toMgrs(self, **Mgrs_and_kwds):
        '''Return this L{Mgrs4Tuple} as an L{Mgrs} instance.
        '''
        return self.to6Tuple(NN, _WGS84).toMgrs(**Mgrs_and_kwds)

    def to6Tuple(self, band=NN, datum=_WGS84):
        '''Extend this L{Mgrs4Tuple} to a L{Mgrs6Tuple}.

           @kwarg band: The band (C{str}).
           @kwarg datum: The datum (L{Datum}).

           @return: An L{Mgrs6Tuple}C{(zone, EN, easting,
                    northing, band, datum)}.
        '''
        z = self.zone  # PYCHOK or [0]
        B = z[-1:]
        if B.isalpha():
            z = z[:-1] or Fmt.zone(0)
            t = Mgrs6Tuple(z, self.EN, self.easting, self.northing,  # PYCHOK attrs
                              band or B, datum, name=self.name)
        else:
            t = self._xtend(Mgrs6Tuple, band, datum)
        return t


class Mgrs6Tuple(_NamedTuple):  # XXX only used above
    '''6-Tuple C{(zone, EN, easting, northing, band, datum)}, with
       C{zone}, grid tile C{EN} and C{band} as C{str}, C{easting}
       and C{northing} in C{meter} and C{datum} a L{Datum}.

       @note: The C{zone} is the I{longitudinal} zone C{"01".."60"}
              or C{"00"} for I{polar} regions and C{band} is the
              I{latitudinal} band or I{polar} region letter.
    '''
    _Names_ = Mgrs4Tuple._Names_ + (_band_, _datum_)
    _Units_ = Mgrs4Tuple._Units_ + ( Str,   _Pass)

    @deprecated_property_RO
    def digraph(self):
        '''DEPRECATED, use attribute C{EN}.'''
        return self.EN  # PYCHOK or [1]

    def toMgrs(self, Mgrs=Mgrs, **Mgrs_kwds):
        '''Return this L{Mgrs6Tuple} as an L{Mgrs} instance.
        '''
        kwds = dict(self.items())
        if self.name:
            kwds.update(name=self.name)
        if Mgrs_kwds:
            kwds.update(Mgrs_kwds)
        return Mgrs(**kwds)


class _RE(object):
    '''(INTERNAL) Lazily compiled C{re}gex-es to parse MGRS strings.
    '''
    _EN = '([A-Z]{2})'            # 2-letter grid tile designation
    _en = '([0-9]+)'              # easting_northing digits, 2-10+
    _pB = '([ABYZ]{1})'           # polar region letter, pseudo-zone 0
    _zB = '([0-9]{1,2}[C-X]{1})'  # zone number and band letter, no I|O

    @Property_RO
    def pB_EN(self):  # split polar "BEN" into 2 parts
        import re  # PYCHOK warning locale.Error
        return re.compile(_RE._pB + _RE._EN, re.IGNORECASE)

    @Property_RO
    def pB_EN_en(self):  # split polar "BEN1235..." into 3 parts
        import re  # PYCHOK warning locale.Error
        return re.compile(_RE._pB + _RE._EN + _RE._en, re.IGNORECASE)

    @Property_RO
    def zB_EN(self):  # split "1[2]BEN" into 2 parts
        import re  # PYCHOK warning locale.Error
        return re.compile(_RE._zB + _RE._EN, re.IGNORECASE)

    @Property_RO
    def zB_EN_en(self):  # split "1[2]BEN1235..." into 3 parts
        import re  # PYCHOK warning locale.Error
        return re.compile(_RE._zB + _RE._EN + _RE._en, re.IGNORECASE)

_RE = _RE()  # PYCHOK singleton


def parseMGRS(strMGRS, datum=_WGS84, Mgrs=Mgrs, name=NN):
    '''Parse a string representing a MGRS grid reference,
       consisting of C{"[zone]Band, EN, easting, northing"}.

       @arg strMGRS: MGRS grid reference (C{str}).
       @kwarg datum: Optional datum to use (L{Datum}).
       @kwarg Mgrs: Optional class to return the MGRS grid
                    reference (L{Mgrs}) or C{None}.
       @kwarg name: Optional B{C{Mgrs}} name (C{str}).

       @return: The MGRS grid reference as B{C{Mgrs}} or if
                C{B{Mgrs} is None} as an L{Mgrs4Tuple}C{(zone,
                EN, easting, northing)}.

       @raise MGRSError: Invalid B{C{strMGRS}}.

       @example:

        >>> m = parseMGRS('31U DQ 48251 11932')
        >>> str(m)  # '31U DQ 48251 11932'
        >>> m = parseMGRS('31UDQ4825111932')
        >>> repr(m)  # [Z:31U, G:DQ, E:48251, N:11932]
        >>> m = parseMGRS('42SXD0970538646')
        >>> str(m)  # '42S XD 09705 38646'
        >>> m = parseMGRS('42SXD9738')  # Km
        >>> str(m)  # '42S XD 97000 38000'
        >>> m = parseMGRS('YUB17770380')  # polar
        >>> str(m)  # 'Y UB 17770 03800'
    '''
    def _mg(s, re_UTM, re_UPS):  # return re.match groups
        m = re_UTM.match(s)
        if m:
            return m.groups()
        m = re_UPS.match(s.lstrip(_0_))
        if m:
            return m.groups()
#           m = m.groups()
#           t = '00' + m[0]
#           return (t,) + m[1:]
        raise ValueError(_SPACE_(repr(s), _invalid_))

    def _MGRS(strMGRS, datum, Mgrs, name):
        m = _splituple(strMGRS.strip())
        if len(m) == 1:  # [01]BEN1234512345'
            m = _mg(m[0], _RE.zB_EN_en, _RE.pB_EN_en)
            m = m[:2] + halfs2(m[2])
        elif len(m) == 2:  # [01]BEN 1234512345'
            m = _mg(m[0], _RE.zB_EN, _RE.pB_EN) + halfs2(m[1])
        elif len(m) == 3:  # [01]BEN 12345 12345'
            m = _mg(m[0], _RE.zB_EN, _RE.pB_EN) + m[1:]
        if len(m) != 4:  # [01]B EN 12345 12345
            raise ValueError

        zB, EN = m[0].upper(), m[1].upper()
        if zB[-1:] in 'IO':
            raise ValueError(_SPACE_(repr(m[0]), _invalid_))
        e, n, m = _enstr2m3(*m[2:])

        if Mgrs is None:
            r = Mgrs4Tuple(zB, EN, e, n, name=name)
            _ = r.toMgrs(resolution=m)  # validate
        else:
            r = Mgrs(zB, EN, e, n, datum=datum, resolution=m, name=name)
        return r

    return _parseX(_MGRS, strMGRS, datum, Mgrs, name,
                          strMGRS=strMGRS, Error=MGRSError)


def toMgrs(utmups, Mgrs=Mgrs, name=NN, **Mgrs_kwds):
    '''Convert a UTM or UPS coordinate to an MGRS grid reference.

       @arg utmups: A UTM or UPS coordinate (L{Utm}, L{Etm} or L{Ups}).
       @kwarg Mgrs: Optional class to return the MGRS grid reference
                    (L{Mgrs}) or C{None}.
       @kwarg name: Optional B{C{Mgrs}} name (C{str}).
       @kwarg Mgrs_kwds: Optional, additional B{C{Mgrs}} keyword
                         arguments, ignored if C{B{Mgrs} is None}.

       @return: The MGRS grid reference as B{C{Mgrs}} or if
                C{B{Mgrs} is None} as an L{Mgrs6Tuple}C{(zone,
                EN, easting, northing, band, datum)}.

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
        raise MGRSError(B=B, E=E, N=N, utmups=utmups, cause=x)

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

    from pygeodesy.ellipsoidalVincenty import fabs, LatLon
    from pygeodesy.lazily import _getenv, printf

#   from math import fabs  # from .ellipsoidalVincenty
    from os import access as _access, linesep as _NL, X_OK as _X_OK

    # <https://GeographicLib.sourceforge.io/C++/doc/GeoConvert.1.html>
    _GeoConvert = _getenv(_PYGEODESY_GEOCONVERT_, '/opt/local/bin/GeoConvert')
    if _access(_GeoConvert, _X_OK):
        GC_m = _GeoConvert, '-m'  # -m converts latlon to MGRS
        printf(' using: %s ...', _SPACE_.join(GC_m))
        from pygeodesy.solveBase import _popen2
    else:
        GC_m = _popen2 = None

    e = n = 0
    try:
        for lat in range(-90, 91, 1):
            printf('%6s: lat %s ...', n, lat, end=NN, flush=True)
            nl = _NL
            for lon in range(-180, 181, 1):
                m = LatLon(lat, lon).toMgrs()
                if _popen2:
                    t = '%s %s' % (lat, lon)
                    g = _popen2(GC_m, stdin=t)[1]
                    t =  m.toStr()  # sep=NN
                    if t != g:
                        e += 1
                        printf('%s%6s: %s: %r vs %r (lon %s)', nl, -e, m, t, g, lon)
                        nl = NN
                t = m.toLatLon(LatLon=LatLon)
                d = max(fabs(t.lat - lat), fabs(t.lon - lon))
                if d > 1e-9 and -90 < lat < 90 and -180 < lon < 180:
                    e += 1
                    printf('%s%6s: %s: %s vs %s %.6e', nl, -e, m, t.latlon, (float(lat), float(lon)), d)
                    nl = NN
                n += 1
            if nl:
                printf(' OK')
    except KeyboardInterrupt:
        printf(nl)

    p = e * 100.0 / n
    printf('%6s: %s errors (%.2f%%)', n, (e if e else 'no'), p)

# **) MIT License
#
# Copyright (C) 2016-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
