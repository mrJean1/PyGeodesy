
# -*- coding: utf-8 -*-

u'''Military Grid Reference System (MGRS/NATO) class L{Mgrs} and
functions L{parseMGRS} and L{toMgrs}.

Pure Python implementation of MGRS / UTM conversion functions using
an ellipsoidal earth model, transcribed from JavaScript originals by
I{(C) Chris Veness 2014-2016} published under the same MIT Licence**, see
U{MGRS<http://www.Movable-Type.co.UK/scripts/latlong-utm-mgrs.html>} and
U{Module mgrs<http://www.Movable-Type.co.UK/scripts/geodesy/docs/module-mgrs.html>}.

The MGRS/NATO grid references provides geocoordinate references
covering the entire globe, based on UTM projections.

MGRS references comprise a grid zone designation, a 100 km square
identification, and an easting and northing (in metres).

Depending on requirements, some parts of the reference may be omitted
(implied), and easting/northing may be given to varying resolution.

See also U{United States National Grid
<http://www.FGDC.gov/standards/projects/FGDC-standards-projects/usng/fgdc_std_011_2001_usng.pdf>}
and U{Military Grid Reference System<http://WikiPedia.org/wiki/Military_grid_reference_system>}.

@newfield example: Example, Examples
'''

from bases import Based, _xattrs, _xnamed
from datum import Datums
from utily import enStr2, halfs, property_RO
from utm import toUtm, Utm, _toZBlat3

from math import log10
import re  # PYCHOK warning locale.Error

# all public contants, classes and functions
__all__ = ('Mgrs',  # classes
           'parseMGRS', 'toMgrs')  # functions
__version__ = '18.10.10'

_100km  =  100e3  #: (INTERNAL) 100 km in meter.
_2000km = 2000e3  #: (INTERNAL) 2,000 km in meter.

# 100 km grid square column (‘e’) letters repeat every third zone
_Le100k = 'ABCDEFGH', 'JKLMNPQR', 'STUVWXYZ'  #: (INTERNAL) Grid E colums.
# 100 km grid square row (‘n’) letters repeat every other zone
_Ln100k = 'ABCDEFGHJKLMNPQRSTUV', 'FGHJKLMNPQRSTUVABCDE'  #: (INTERNAL) Grid N rows.

# split an MGRS string "12ABC1235..." into 3 parts
_MGRSre = re.compile('([0-9]{1,2}[C-X]{1})([A-Z]{2})([0-9]+)', re.IGNORECASE)  #: (INTERNAL) Regex.
_GZDre  = re.compile('([0-9]{1,2}[C-X]{1})', re.IGNORECASE)  #: (INTERNAL) Regex.


class Mgrs(Based):
    '''Military Grid Reference System (MGRS/NATO) references,
       with method to convert to UTM coordinates.
    '''
    _band     = ''  #: (INTERNAL) Latitudinal band (C..X).
    _bandLat  = None  #: (INTERNAL) Band latitude (C{degrees90} or C{None}).
    _datum    = Datums.WGS84  #: (INTERNAL) Datum (L{Datum}).
    _easting  = 0   #: (INTERNAL) Easting (C{meter}), within 100 km grid square.
    _en100k   = ''  #: (INTERNAL) Grid EN digraph (C{str}), 100 km grid square.
    _northing = 0   #: (INTERNAL) Northing (C{meter}), within 100 km grid square.
    _zone     = 0   #: (INTERNAL) Longitudinal zone (C{int}), 1..60.

    def __init__(self, zone, en100k, easting, northing,
                             band='', datum=Datums.WGS84, name=''):
        '''New MGRS grid reference.

           @param zone: 6° longitudinal zone (C{int}), 1..60 covering 180°W..180°E.
           @param en100k: Two-letter EN digraph (C{str}), 100 km grid square.
           @param easting: Easting (C{meter}), within 100 km grid square.
           @param northing: Northing (C{meter}), within 100 km grid square.
           @keyword band: Optional 8° latitudinal band (C{str}), C..X covering 80°S..84°N.
           @keyword datum: Optional this reference's datum (L{Datum}).
           @keyword name: Optional name (C{str}).

           @raise ValueError: Invalid MGRS grid reference, I{zone},
                              I{en100k} or I{band}.

           @example:

           >>> from mgrs import Mgrs
           >>> m = Mgrs('31U', 'DQ', 48251, 11932)  # 31U DQ 48251 11932
        '''
        if name:
            self.name = name

        self._zone, self._band, self._bandLat = _toZBlat3(zone, band, True)

        try:
            en = str(en100k).upper()
            if len(en) != 2:
                raise IndexError  # caught below
            self._en100k = en
            self._en100k2m()
        except IndexError:
            raise ValueError('%s invalid: %r' ('en100k', en100k))

        self._easting, self._northing = float(easting), float(northing)

        if self._datum != datum:
            self._datum = datum

    def _en100k2m(self):
        # check and convert grid letters to meter
        z = self._zone - 1
        # get easting specified by e100k (note, +1 because
        # eastings start at 166e3 due to 500 km false origin)
        e = float(_Le100k[z % 3].index(self._en100k[0]) + 1) * _100km  # metres
        # similarly, get northing specified by n100k
        n = float(_Ln100k[z % 2].index(self._en100k[1])) * _100km  # metres
        return e, n

    def _xcopy(self, *attrs):
        '''(INTERNAL) Make copy with add'l, subclass attributes.
        '''
        return _xattrs(self.classof(self.zone, self.en100k,
                                    self.easting, self.northing,
                                    band=self.band, datum=self.datum),
                       self, *attrs)

    @property_RO
    def band(self):
        '''Get the latitudinal band (A..Z C{str}).
        '''
        return self._band

    @property_RO
    def bandLatitude(self):
        '''Get the band latitude (C{degrees90} or C{None}).
        '''
        return self._bandLat

    def copy(self):
        '''Copy this MGRS reference.

           @return: The copy (L{Mgrs} or subclass thereof).
        '''
        return self._xcopy()

    @property_RO
    def datum(self):
        '''Get the datum (L{Datum}).
        '''
        return self._datum

    @property_RO
    def en100k(self):
        '''Get the 2-character grid (EN digraph C{str}).
        '''
        return self._en100k

    @property_RO
    def easting(self):
        '''Gets the easting (C{meter}).
        '''
        return self._easting

    @property_RO
    def northing(self):
        '''Get the northing (C{meter}).
        '''
        return self._northing

    def parse(self, strMGRS):
        '''Parse a string to a MGRS grid reference.

           @param strMGRS: MGRS grid reference (C{str}).

           @return: MGRS reference (L{Mgrs}).

           @raise ValueError: Invalid I{strMGRS}.

           @see: Function L{parseMGRS} in this module L{mgrs}.
        '''
        return parseMGRS(strMGRS, datum=self.datum)

    def toStr(self, prec=10, sep=' '):  # PYCHOK expected
        '''Return a string representation of this MGRS grid reference.

           Note that MGRS grid references are truncated, not rounded
           (unlike UTM coordinates).

           @keyword prec: Optional number of digits (C{int}), 4:km, 10:m.
           @keyword sep: Optional separator to join (C{str}).

           @return: This Mgrs as "00B EN easting northing" (C{str}).

           @raise ValueError: Invalid I{prec}.

           @example:

           >>> m = Mgrs(31, 'DQ', 48251, 11932, band='U')
           >>> m.toStr()  # '31U DQ 48251 11932'
        '''
        t = enStr2(self._easting, self._northing, prec,
                   '%02d%s' % (self._zone, self._band), self._en100k)
        return sep.join(t)

    def toStr2(self, prec=10, fmt='[%s]', sep=', '):  # PYCHOK expected
        '''Return a string representation of this MGRS grid reference.

           @keyword prec: Optional number of digits (C{int}), 4:km, 10:m.
           @keyword fmt: Optional enclosing backets format (C{str}).
           @keyword sep: Optional separator between name:values (C{str}).

           @return: This Mgrs as "[Z:00B, G:EN, E:meter, N:meter]" (C{str}).
        '''
        t = self.toStr(prec=prec, sep=' ').split()
        return fmt % (sep.join('%s:%s' % t for t in zip('ZGEN', t)),)

    def toUtm(self, Utm=Utm):
        '''Convert this MGRS grid reference to a UTM coordinate.

           @keyword Utm: Optional (sub-)class to use for the UTM
                         coordinate (L{Utm}) or C{None}.

           @return: The UTM coordinate (L{Utm}) or 4-tuple (zone,
                    hemisphere, easting, northing) if I{Utm} is C{None}.

           @example:

           >>> m = Mgrs('31U', 'DQ', 448251, 11932)
           >>> u = m.toUtm()  # 31 N 448251 5411932
        '''
        # get northing of the band bottom, extended to
        # include entirety of bottom-most 100 km square
        n = toUtm(self._bandLat, 0, datum=self._datum).northing
        nb = int(n / _100km) * _100km

        e, n = self._en100k2m()
        # 100 km grid square row letters repeat every 2,000 km north;
        # add enough 2,000 km blocks to get into required band
        e += self._easting
        n += self._northing
        while n < nb:
            n += _2000km

        h = 'S' if self._bandLat < 0 else 'N'  # if self._band < 'N'

        return (self.zone, h, e, n) if Utm is None else _xnamed(Utm(
                self.zone, h, e, n, band=self.band, datum=self.datum), self.name)

    @property_RO
    def zone(self):
        '''Get the longitudal zone (1..60 C{int}).
        '''
        return self._zone


def parseMGRS(strMGRS, datum=Datums.WGS84, Mgrs=Mgrs, name=''):
    '''Parse a string representing a MGRS grid reference,
       consisting of zoneBand, grid, easting and northing.

       @param strMGRS: MGRS grid reference (C{str}).
       @keyword datum: Optional datum to use (L{Datum}).
       @keyword Mgrs: Optional (sub-)class to use for the MGRS
                      grid reference (L{Mgrs}) or C{None}.
       @keyword name: Optional I{Mgrs} name (C{str}).

       @return: The MGRS grid reference (L{Mgrs}) or 4-tuple (zone,
                ENdigraph, easting, northing) if I{Mgrs} is C{None}.

       @raise ValueError: Invalid I{strMGRS}.

       @example:

       >>> m = parseMGRS('31U DQ 48251 11932')
       >>> str(m)  # 31U DQ 48251 11932
       >>> m = parseMGRS('31UDQ4825111932')
       >>> repr(m)  # [Z:31U, G:DQ, E:48251, N:11932]
    '''
    def _mg(cre, s):  # return re.match groups
        m = cre.match(s)
        if not m:
            raise ValueError
        return m.groups()

    def _s2m(g):  # e or n string to meter
        f = float(g)
        if f > 0:
            x = int(log10(f))
            if 0 <= x < 4:  # at least 5 digits
                f *= (10000, 1000, 100, 10)[x]
        return f

    m = tuple(strMGRS.strip().replace(',', ' ').split())
    try:
        if len(m) == 1:  # 01ABC1234512345'
            m = _mg(_MGRSre, m[0])
            m = m[:2] + halfs(m[2])
        elif len(m) == 2:  # 01ABC 1234512345'
            m = _mg(_GZDre, m[0]) + halfs(m[1])
        elif len(m) == 3:  # 01ABC 12345 12345'
            m = _mg(_GZDre, m[0]) + m[1:]
        if len(m) != 4:  # 01A BC 1234 12345
            raise ValueError
        e, n = map(_s2m, m[2:])
    except ValueError:
        raise ValueError('%s invalid: %r' % ('strMGRS', strMGRS))

    z, EN = m[0], m[1].upper()
    return (z, EN, e, n) if Mgrs is None else _xnamed(Mgrs(
            z, EN, e, n, datum=datum), name)


def toMgrs(utm, Mgrs=Mgrs, name=''):
    '''Convert a UTM coordinate to an MGRS grid reference.

       @param utm: A UTM coordinate (L{Utm}).
       @keyword Mgrs: Optional (sub-)class to use for the MGRS
                      grid reference (L{Mgrs}).

       @return: The MGRS grid reference (L{Mgrs}).
       @keyword name: Optional I{Mgrs} name (C{str}).

       @raise TypeError: If I{utm} is not L{Utm}.

       @raise ValueError: Invalid I{utm}.

       @example:

       >>> u = Utm(31, 'N', 448251, 5411932)
       >>> m = u.toMgrs()  # 31U DQ 48251 11932
    '''
    if not isinstance(utm, Utm):
        raise TypeError('%s not Utm: %s' % ('utm', type(utm)))

    # truncate east-/northing to within 100 km grid square
    # XXX add rounding to nm precision?
    E, e = divmod(utm.easting, _100km)
    N, n = divmod(utm.northing, _100km)

    # columns in zone 1 are A-H, zone 2 J-R, zone 3 S-Z, then
    # repeating every 3rd zone (note -1 because eastings start
    # at 166e3 due to 500km false origin)
    z = utm.zone - 1
    en = (_Le100k[z % 3][int(E) - 1] +
          # rows in even zones are A-V, in odd zones are F-E
          _Ln100k[z % 2][int(N) % len(_Ln100k[0])])

    return _xnamed(Mgrs(utm.zone, en, e, n, band=utm.band, datum=utm.datum),
                   name or utm.name)

# **) MIT License
#
# Copyright (C) 2016-2018 -- mrJean1 at Gmail dot com
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
