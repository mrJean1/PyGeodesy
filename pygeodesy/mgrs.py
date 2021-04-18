
# -*- coding: utf-8 -*-

u'''Military Grid Reference System (MGRS/NATO) classes L{Mgrs} and
L{MGRSError} and functions L{parseMGRS} and L{toMgrs}.

Pure Python implementation of MGRS / UTM conversion functions using
an ellipsoidal earth model, transcribed from JavaScript originals by
I{(C) Chris Veness 2014-2016} published under the same MIT Licence**, see
U{MGRS<https://www.Movable-Type.co.UK/scripts/latlong-utm-mgrs.html>} and
U{Module mgrs<https://www.Movable-Type.co.UK/scripts/geodesy/docs/module-mgrs.html>}.

The MGRS/NATO grid references provides geocoordinate references
covering the entire globe, based on UTM projections.

MGRS references comprise a grid zone designation, a 100 km square
identification, and an easting and northing (in metres).

Depending on requirements, some parts of the reference may be omitted
(implied), and easting/northing may be given to varying resolution.

See also U{United States National Grid
<https://www.FGDC.gov/standards/projects/FGDC-standards-projects/usng/fgdc_std_011_2001_usng.pdf>}
and U{Military Grid Reference System<https://WikiPedia.org/wiki/Military_grid_reference_system>}.

@newfield example: Example, Examples
'''

from pygeodesy.basics import halfs2, _xinstanceof
from pygeodesy.datums import _ellipsoidal_datum, _WGS84
from pygeodesy.errors import _parseX, _ValueError, _xkwds
from pygeodesy.interns import NN, _AtoZnoIO_, _band_, \
                             _COMMASPACE_, _datum_, _easting_, \
                             _northing_, _SPACE_, _zone_
from pygeodesy.interns import _0_0  # PYCHOK used!
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import _NamedBase, _NamedTuple, _Pass, _xnamed
from pygeodesy.namedTuples import UtmUps5Tuple
from pygeodesy.props import Property_RO
from pygeodesy.streprs import enstr2, Fmt, _xzipairs
from pygeodesy.units import Easting, Northing, Str, _100km
from pygeodesy.units import _2000km  # PYCHOK used!
from pygeodesy.utm import toUtm8, _to3zBlat, Utm
from pygeodesy.utmupsBase import _hemi


__all__ = _ALL_LAZY.mgrs
__version__ = '21.04.15'

# 100 km grid square column (‘e’) letters repeat every third zone
_Le100k = _AtoZnoIO_.tillH, _AtoZnoIO_.fromJ.tillR, _AtoZnoIO_.fromS  # grid E colums
# 100 km grid square row (‘n’) letters repeat every other zone
_Ln100k = _AtoZnoIO_.tillV, _AtoZnoIO_.fromF.tillV + _AtoZnoIO_.tillE  # grid N rows


class _RE(object):
    '''(INTERNAL) Lazily compiled regex-es.
    '''
    @Property_RO
    def MGRS(self):  # split an MGRS string "12ABC1235..." into 3 parts
        import re  # PYCHOK warning locale.Error
        return re.compile('([0-9]{1,2}[C-X]{1})([A-Z]{2})([0-9]+)', re.IGNORECASE)

    @Property_RO
    def ZBG(self):  # split an MGRS string "12ABC1235..." into 3 parts
        import re  # PYCHOK warning locale.Error
        return re.compile('([0-9]{1,2}[C-X]{1})([A-Z]{2})', re.IGNORECASE)

_RE = _RE()  # PYCHOK singleton


class MGRSError(_ValueError):
    '''Military Grid Reference System (MGRS) parse or other L{Mgrs} issue.
    '''
    pass


class Mgrs(_NamedBase):
    '''Military Grid Reference System (MGRS/NATO) references,
       with method to convert to UTM coordinates.
    '''
    _band     =  NN     # latitudinal band (C..X)
    _bandLat  =  None   # band latitude (C{degrees90} or C{None})
    _datum    = _WGS84  # Datum (L{Datum})
    _easting  =  0      # Easting (C{meter}), within 100 km grid square
    _en100k   =  NN     # grid EN digraph (C{str}), 100 km grid square
    _northing =  0      # Northing (C{meter}), within 100 km grid square
    _zone     =  0      # longitudinal zone (C{int}), 1..60

    def __init__(self, zone, en100k, easting, northing,
                             band=NN, datum=_WGS84, name=NN):
        '''New L{Mgrs} Military grid reference.

           @arg zone: 6° longitudinal zone (C{int}), 1..60 covering 180°W..180°E.
           @arg en100k: Two-letter EN digraph (C{str}), 100 km grid square.
           @arg easting: Easting (C{meter}), within 100 km grid square.
           @arg northing: Northing (C{meter}), within 100 km grid square.
           @kwarg band: Optional 8° latitudinal band (C{str}), C..X covering 80°S..84°N.
           @kwarg datum: Optional this reference's datum (L{Datum}, L{Ellipsoid},
                         L{Ellipsoid2} or L{a_f2Tuple}).
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

        self._zone, self._band, self._bandLat = _to3zBlat(zone, band, MGRSError)

        try:
            en = str(en100k).upper()
            if len(en) != 2:
                raise IndexError  # caught below
            self._en100k = en
            self._en100k2m()
        except IndexError:
            raise MGRSError(en100k=en100k)

        self._easting  = Easting(easting,   Error=MGRSError)
        self._northing = Northing(northing, Error=MGRSError)
        if datum not in (None, self._datum):
            self._datum = _ellipsoidal_datum(datum, name=name)

    def _en100k2m(self):
        # check and convert grid letters to meter
        z = self._zone - 1
        # get easting specified by e100k (note, +1 because
        # eastings start at 166e3 due to 500 km false origin)
        e = float(_Le100k[z % 3].index(self._en100k[0]) + 1) * _100km  # meter
        # similarly, get northing specified by n100k
        n = float(_Ln100k[z % 2].index(self._en100k[1])) * _100km  # meter
        return e, n

    @Property_RO
    def band(self):
        '''Get the latitudinal band (C{str, 'A'|'B'..'Y'|'Z'}).
        '''
        return self._band

    @Property_RO
    def bandLatitude(self):
        '''Get the band latitude (C{degrees90} or C{None}).
        '''
        return self._bandLat

    @Property_RO
    def datum(self):
        '''Get the datum (L{Datum}).
        '''
        return self._datum

    @Property_RO
    def en100k(self):
        '''Get the 2-character grid EN digraph (C{str}).
        '''
        return self._en100k

    digraph = en100k

    @Property_RO
    def easting(self):
        '''Gets the easting (C{meter}).
        '''
        return self._easting

    @Property_RO
    def northing(self):
        '''Get the northing (C{meter}).
        '''
        return self._northing

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

    def toRepr(self, prec=10, fmt=Fmt.SQUARE, sep=_COMMASPACE_):  # PYCHOK expected
        '''Return a string representation of this MGRS grid reference.

           @kwarg prec: Optional number of digits (C{int}), 4:km, 10:m.
           @kwarg fmt: Optional enclosing backets format (C{str}).
           @kwarg sep: Optional separator between name:values (C{str}).

           @return: This Mgrs as "[Z:00B, G:EN, E:meter, N:meter]" (C{str}).
        '''
        t = self.toStr(prec=prec, sep=None)
        return _xzipairs('ZGEN', t, sep=sep, fmt=fmt)

    def toStr(self, prec=10, sep=_SPACE_):  # PYCHOK expected
        '''Return a string representation of this MGRS grid reference.

           Note that MGRS grid references are truncated, not rounded
           (unlike UTM coordinates).

           @kwarg prec: Optional number of digits (C{int}), 4:km, 10:m.
           @kwarg sep: Optional separator to join (C{str}) or C{None}
                       to return an unjoined C{tuple} of C{str}s.

           @return: This Mgrs as "00B EN easting northing" (C{str}).

           @raise ValueError: Invalid B{C{prec}}.

           @example:

            >>> m = Mgrs(31, 'DQ', 48251, 11932, band='U')
            >>> m.toStr()  # '31U DQ 48251 11932'
        '''
        t = NN(Fmt.zone(self._zone), self._band)
        t = enstr2(self._easting, self._northing, prec, t, self._en100k)
        return t if sep is None else sep.join(t)

    def toUtm(self, Utm=Utm):
        '''Convert this MGRS grid reference to a UTM coordinate.

           @kwarg Utm: Optional class to return the UTM coordinate
                       (L{Utm}) or C{None}.

           @return: The UTM coordinate (L{Utm}) or a
                    L{UtmUps5Tuple}C{(zone, hemipole, easting,
                    northing, band)} if B{C{Utm}} is C{None}.

           @example:

            >>> m = Mgrs('31U', 'DQ', 448251, 11932)
            >>> u = m.toUtm()  # 31 N 448251 5411932
            >>> u = m.toUtm(None)  # 31 N 448251 5411932 U
        '''
        r, u = self._utmups5utm2
        if Utm is None:
            u = r
        elif Utm is not u.__class__:
            u = Utm(r.zone, r.hemipole, r.easting, r.northing,
                    band=r.band, datum=self.datum, name=r.name)
        return u

    @Property_RO
    def _utmups5utm2(self):
        '''(INTERNAL) Cache for L{toUtm}.
        '''
        # get northing of the band bottom, extended to
        # include entirety of bottom-most 100 km square
        n  = toUtm8(self.bandLatitude, _0_0, datum=self.datum).northing
        nb = int(n / _100km) * _100km

        e, n = self._en100k2m()
        # 100 km grid square row letters repeat every 2,000 km north;
        # add enough 2,000 km blocks to get into required band
        e += self.easting
        n += self.northing
        while n < nb:
            n += _2000km

        z =  self.zone
        h = _hemi(self.bandLatitude)  # if self.band < _N_
        B =  self.band
        return (UtmUps5Tuple(z, h, e, n, B, Error=MGRSError, name=self.name),
                Utm(         z, h, e, n, band=B, datum=self.datum, name=self.name))

    @Property_RO
    def zone(self):
        '''Get the longitudal zone (C{int}, 1..60).
        '''
        return self._zone


class Mgrs4Tuple(_NamedTuple):
    '''4-Tuple C{(zone, digraph, easting, northing)}, C{zone} and
       C{digraph} as C{str}, C{easting} and C{northing} in C{meter}.
    '''
    _Names_ = (_zone_, 'digraph', _easting_, _northing_)
    _Units_ = ( Str,    Str,       Easting,   Northing)

    def __new__(cls, z, di, e, n, Error=MGRSError):
        if Error is not None:
            e = Easting( e, Error=Error)
            n = Northing(n, Error=Error)
        return _NamedTuple.__new__(cls, z, di, e, n)

    def to6Tuple(self, band, datum):
        '''Extend this L{Mgrs4Tuple} to a L{Mgrs6Tuple}.

           @arg band: The band to add (C{str} or C{None}).
           @arg datum: The datum to add (L{Datum} or C{None}).

           @return: A L{Mgrs6Tuple}C{(zone, digraph, easting,
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

       @return: The MGRS grid reference (B{C{Mgrs}}) or an
                L{Mgrs4Tuple}C{(zone, digraph, easting, northing)}
                if C{B{Mgrs}=None}.

       @raise MGRSError: Invalid B{C{strMGRS}}.

       @example:

        >>> m = parseMGRS('31U DQ 48251 11932')
        >>> str(m)  # 31U DQ 48251 11932
        >>> m = parseMGRS('31UDQ4825111932')
        >>> repr(m)  # [Z:31U, G:DQ, E:48251, N:11932]
        >>> m = mgrs.parseMGRS('42SXD0970538646')
        >>> str(m)  # 42S XD 09705 38646
        >>> m = mgrs.parseMGRS('42SXD9738')  # Km
        >>> str(m)  # 42S XD 97000 38000
    '''
    def _mg(RE, s):  # return re.match groups
        m = RE.match(s)
        if not m:
            raise ValueError
        return m.groups()

    def _s2m(g):  # e or n string to float meter
        # convert to meter if less than 5 digits
        m = g + '00000'
        return float(m[:5])

    def _MGRS_(strMGRS, datum, Mgrs, name):
        m = tuple(strMGRS.replace(',', ' ').strip().split())
        if len(m) == 1:  # 01ABC1234512345'
            m = _mg(_RE.MGRS, m[0])
            m = m[:2] + halfs2(m[2])
        elif len(m) == 2:  # 01ABC 1234512345'
            m = _mg(_RE.ZBG, m[0]) + halfs2(m[1])
        elif len(m) == 3:  # 01ABC 12345 12345'
            m = _mg(_RE.ZBG, m[0]) + m[1:]
        if len(m) != 4:  # 01A BC 1234 12345
            raise ValueError
        e, n = map(_s2m, m[2:])

        z, EN = m[0], m[1].upper()
        r = Mgrs4Tuple(z, EN, e, n, name=name) if Mgrs is None else \
          _xnamed(Mgrs(z, EN, e, n, datum=datum), name)
        return r

    return _parseX(_MGRS_, strMGRS, datum, Mgrs, name,
                           strMGRS=strMGRS, Error=MGRSError)


def toMgrs(utm, Mgrs=Mgrs, name=NN, **Mgrs_kwds):
    '''Convert a UTM coordinate to an MGRS grid reference.

       @arg utm: A UTM coordinate (L{Utm} or L{Etm}).
       @kwarg Mgrs: Optional class to return the MGRS grid
                    reference (L{Mgrs}) or C{None}.
       @kwarg name: Optional B{C{Mgrs}} name (C{str}).
       @kwarg Mgrs_kwds: Optional, additional B{C{Mgrs}} keyword
                         arguments, ignored if C{B{Mgrs}=None}.

       @return: The MGRS grid reference (B{C{Mgrs}}) or an
                L{Mgrs6Tuple}C{(zone, digraph, easting, northing,
                band, datum)} if C{B{Mgrs}=None}.

       @raise TypeError: If B{C{utm}} is not L{Utm} nor L{Etm}.

       @raise MGRSError: Invalid B{C{utm}}.

       @example:

        >>> u = Utm(31, 'N', 448251, 5411932)
        >>> m = u.toMgrs()  # 31U DQ 48251 11932
    '''
    _xinstanceof(Utm, utm=utm)  # Utm, Etm

    e, n = utm.eastingnorthing2(falsed=True)
    # truncate east-/northing to within 100 km grid square
    # XXX add rounding to nm precision?
    E, e = divmod(e, _100km)
    N, n = divmod(n, _100km)

    # columns in zone 1 are A-H, zone 2 J-R, zone 3 S-Z, then
    # repeating every 3rd zone (note -1 because eastings start
    # at 166e3 due to 500km false origin)
    z = utm.zone - 1
    en = (_Le100k[z % 3][int(E) - 1] +
          # rows in even zones are A-V, in odd zones are F-E
          _Ln100k[z % 2][int(N) % len(_Ln100k[0])])

    if Mgrs is None:
        r = Mgrs4Tuple(utm.zone, en, e, n).to6Tuple(utm.band, utm.datum)
    else:
        kwds = _xkwds(Mgrs_kwds, band=utm.band, datum=utm.datum)
        r = Mgrs(utm.zone, en, e, n, **kwds)
    return _xnamed(r, name or utm.name)

# **) MIT License
#
# Copyright (C) 2016-2021 -- mrJean1 at Gmail -- All Rights Reserved.
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
