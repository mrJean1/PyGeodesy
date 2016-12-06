
# -*- coding: utf-8 -*-

# Python implementation of MGRS / UTM conversion functions using
# an ellipsoidal earth model.  Transcribed from JavaScript originals
# by (C) Chris Veness 2014-2016 published under the same MIT Licence**,
# see <http://www.movable-type.co.uk/scripts/latlong-utm-mgrs.html>
# and <http://www.movable-type.co.uk/scripts/geodesy/docs/module-mgrs.html>

from bases import _Base
from datum import Datums
from utils import halfs
from utm   import toUtm, Utm, _toZBL

import re  # PYCHOK warning locale.Error

# Military Grid Reference System (MGRS/NATO) grid references provides
# geocoordinate references covering the entire globe, based on UTM
# projections.

# MGRS references comprise a grid zone designation, a 100 km square
# identification, and an easting and northing (in metres).

# Depending on requirements, some parts of the reference may be omitted
# (implied), and easting/northing may be given to varying resolution.

# Qv <http://www.fgdc.gov/standards/projects/FGDC-standards-projects/usng/fgdc_std_011_2001_usng.pdf>
# and <https://en.wikipedia.org/wiki/Military_grid_reference_system>

# all public contants, classes and functions
__all__ = ('Mgrs',  # classes
           'parseMGRS', 'toMgrs')  # functions
__version__ = '16.11.11'

_100km  =  100e3  # 100 km in meter
_2000km = 2000e3  # 2,000 km in meter

# 100 km grid square column (‘e’) letters repeat every third zone
_Le100k = 'ABCDEFGH', 'JKLMNPQR', 'STUVWXYZ'
# 100 km grid square row (‘n’) letters repeat every other zone
_Ln100k = 'ABCDEFGHJKLMNPQRSTUV', 'FGHJKLMNPQRSTUVABCDE'

# split an MGRS string "12ABC1235..." into 3 parts
_MGRSre = re.compile('(\d{1,2}[C-X]{1})([A-Z]{2})(\d+)', re.IGNORECASE)
_GZDre  = re.compile('(\d{1,2}[C-X]{1})', re.IGNORECASE)


class Mgrs(_Base):
    '''Military Grid Reference System (MGRS/NATO) references,
       with method to convert to UTM coordinates.
    '''
    _band     = ''
    _bandLat  = None
    _datum    = Datums.WGS84
    _en100k   = ''
    _easting  = 0
    _northing = 0
    _zone     = 0

    def __init__(self, zone, en100k, easting, northing,
                             band='', datum=Datums.WGS84):
        '''Create an Mgrs grid reference object.

           @param  {number} zone - 6° longitudinal zone (1..60 covering 180°W..180°E).
           @param  {string} band - 8° latitudinal band (C..X covering 80°S..84°N).
           @param  {string} en100k - Two-letter (EN digraph) 100 km grid square.
           @param  {meter} easting - Easting in metres within 100 km grid square.
           @param  {meter} northing - Northing in metres within 100 km grid square.
           @param  {Datum} [datum=Datums.WGS84] - This reference's datum.

           @throws {Error} Invalid MGRS grid reference.

           @example
           from mgrs import Mgrs
           m = Mgrs('31U', 'DQ', 48251, 11932)  # 31U DQ 48251 11932
        '''
        self._zone, self._band, self._bandLat = _toZBL(zone, band, True)

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

    @property
    def band(self):
        '''Return latitudinal band (A..Z).'''
        return self._band

    @property
    def bandLatitude(self):
        '''Return band latitude in degrees90 or None.'''
        return self._bandLat

    @property
    def datum(self):
        '''Return the datum.'''
        return self._datum

    @property
    def en100k(self):
        '''Return 2-character grid.'''
        return self._en100k

    @property
    def easting(self):
        '''Return easting in meter.'''
        return self._easting

    @property
    def northing(self):
        '''Return northing in meter.'''
        return self._northing

    def parse(self, strMGRS):
        '''Parse a string to a MGRS grid reference.

           For details, see function parseMGRS in
           this module mrgs.
        '''
        return parseMGRS(strMGRS, datum=self.datum)

    def toStr(self, prec=10, sep=' '):  # PYCHOK expected
        '''Returns a string representation of this MGRS grid reference.

           Note that MGRS grid references are truncated, not rounded
           (unlike UTM coordinates).

           @param {number} [prec=10] - Number of digits (4:km, 10:m).
           @param {string} [sep=' '] - Separator to join.

           @returns {string} Mgrs as string "00B EN meter meter".

           @example
           m = Mgrs(31, 'DQ', 48251, 11932, band='U')
           m.toStr()  # '31U DQ 48251 11932'
        '''
        w = prec // 2
        if 1 > w or w > 5:
            raise ValueError('%s invalid: %r' % ('prec', prec))
        p = (0, 1e-4, 1e-3, 1e-2, 1e-1, 1)[w]  # 10 ** (5 - w)

        t = ['%02d%s' % (self._zone, self._band), self._en100k,
             '%0*d' % (w, int(self._easting * p)),
             '%0*d' % (w, int(self._northing * p))]
        return sep.join(t)

    def toStr2(self, prec=10, fmt='[%s]', sep=', '):  # PYCHOK expected
        '''Returns a string representation of this MGRS grid reference.

           @param {number} [prec=10] - Number of digits (4:km, 10:m).
           @param {string} [fmt='[%s]'] - Enclosing backets format.
           @param {string} [sep=', '] - Separator between name:values.

           @returns {string} This Mgrs as "[Z:00B, G:EN, E:meter, N:meter]".
        '''
        t = self.toStr(prec=prec, sep=' ').split()
        return fmt % (sep.join('%s:%s' % t for t in zip('ZGEN', t)),)

    def toUtm(self):
        '''Convert this MGRS grid reference to a UTM coordinate.

           @returns {Utm} The UTM coordinate equivalent to
                          the MGRS grid reference.

           @example
           m = Mgrs('31U', 'DQ', 448251, 11932)
           u = m.toUtm()  # 31 N 448251 5411932
        '''
        # get northing of the band bottom, extended to
        # include entirety of bottom-most 100 km square
        n = toUtm(self._bandLat, 0, datum=self._datum).northing
        nb = int(n / _100km) * _100km

        e, n = self._en100k2m()
        # 100 km grid square row letters repeat every 2,000 km north;
        # add enough 2,000 km blocks to get into required band
        n += self._northing
        while n < nb:
            n += _2000km

        h = 'S' if self._bandLat < 0 else 'N'  # if self._band < 'N'

        return Utm(self._zone, h, e + self._easting, n, band=self._band, datum=self._datum)

    @property
    def zone(self):
        '''Return longitudal zone (1..60).'''
        return self._zone


def parseMGRS(strMGRS, datum=Datums.WGS84):
    '''Parse a string representing a MGRS grid reference,
       consisting of zoneBand, grid, easting and northing.

       @param {string} strMGRS - MGRS grid reference.
       @param {Datum} [datum=Datums.WGS84] - The datum to use.

       @returns {Mgrs} The MGRS grid reference.

       @throws {ValueError} Invalid strMGRS.

       @example
       m = parseMGRS('31U DQ 48251 11932')
       str(m)  # 31U DQ 48251 11932
       m = parseMGRS('31UDQ4825111932')
       repr(m)  # [Z:31U, G:DQ, E:48251, N:11932]
    '''
    def _mg(cre, s):  # return re.match groups
        m = cre.match(s)
        if not m:
            raise ValueError
        return m.groups()

    def _s2m(g):  # e or n string to meter
        f = float(g)
        n = int(f)
        if n:
            n = len(str(n))
        if n < 5:
            f *= pow(10, 5 - n)
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

    return Mgrs(m[0], m[1].upper(), e, n, datum=datum)


def toMgrs(utm):
    '''Convert a UTM coordinate to an MGRS grid reference.

       @param {Utm} utm - A UTM coordinate.

       @returns {Mgrs} The MGRS grid reference.

       @throws {TypeError} Invalid UTM coordinate.

       @example
       u = Utm(31, 'N', 448251, 5411932)
       m = u.toMgrs()  # 31U DQ 48251 11932
    '''
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

    return Mgrs(utm.zone, en, e, n, band=utm.band, datum=utm.datum)

# **) MIT License
#
# Copyright (c) 2016-2017 -- mrJean1@Gmail.com
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
