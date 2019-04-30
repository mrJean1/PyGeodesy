
# -*- coding: utf-8 -*-

u'''Classes L{Epsg} and L{EPSGError} and functions to L{encode} and L{decode2}
U{EPSG (European Petroleum Survey Group)<http://www.EPSG-Registry.org>} codes
from and to U{UTM
<http://WikiPedia.org/wiki/Universal_Transverse_Mercator_coordinate_system>} and
U{UPS<http://WikiPedia.org/wiki/Universal_polar_stereographic_coordinate_system>}
zones.

A pure Python implementation transcribed from C++ class U{UTMUPS
<http://GeographicLib.SourceForge.io/html/classGeographicLib_1_1UTMUPS.html>}
by I{Charles Karney}, including coverage of UPS as zone C{0}.
'''

from bases import _Named
from ellipsoidalBase import _to3zBhp, _UPS_ZONE, \
                            _UTM_ZONE_MIN, _UTM_ZONE_MAX, \
                            _UTMUPS_ZONE_INVALID
from lazily import _ALL_LAZY
from utily import property_RO, _Strs
from ups import Ups
from utm import Utm

# all public contants, classes and functions
__all__ = _ALL_LAZY.epsg + ('decode2', 'encode')
__version__ = '19.04.26'

# _EPSG_INVALID = _UTMUPS_ZONE_INVALID
_EPSG_N_01 = 32601  # EPSG code for UTM zone 01 N
_EPSG_N_60 = 32660  # EPSG code for UTM zone 60 N
_EPSG_N    = 32661  # EPSG code for UPS pole N

_EPSG_S_01 = 32701  # EPSG code for UTM zone 01 S
_EPSG_S_60 = 32760  # EPSG code for UTM zone 60 S
_EPSG_S    = 32761  # EPSG code for UPS pole S


class EPSGError(ValueError):
    '''EPSG encode or decode error.
    '''
    pass


class Epsg(int, _Named):
    '''U{EPSG<http://www.EPSG-Registry.org>} class, sub-class of C{int}.
    '''
    _band       = ''
    _epsg       = None
    _hemisphere = ''
    _utmups     = None
    _zone       = _UTMUPS_ZONE_INVALID

    def __new__(cls, eisu):
        '''New L{Epsg} from a UTM/USP coordinate or other EPSG.

           @param eisu: Other (L{Epsg}, C{int}, C{str}, L{Utm} or L{Ups}).

           @return: New L{Epsg}.

           @raise TypeError: Invalid I{eisu}.

           @raise EPSGError: Invalid I{eisu}.
        '''
        if isinstance(eisu, Epsg):
            self = int.__new__(cls, int(eisu))
            self._band       = eisu.band
            self._epsg       = eisu.epsg
            self._hemisphere = eisu.hemisphere
            self._utmups     = eisu.utmups
            self._zone       = eisu.zone
            if eisu.name:
                self.name = eisu.name

        elif isinstance(eisu, int):
            self = int.__new__(cls, eisu)
            self._epsg = eisu
            self._zone, self._hemisphere = decode2(eisu)

        elif isinstance(eisu, _Strs):
            self = encode(eisu)

        elif isinstance(eisu, (Utm, Ups)):
            u = eisu
            self = encode(u.zone, hemipole=u.hemisphere, band=u.band)  # PYCHOK **kwds
            self._utmups = u
            if u.name:
                self._name = u.name
        else:
            raise TypeError('%s invalid: %r' % ('eisu', eisu))

        return self

    def __repr__(self):
        return '%s(%s)' % (self.name or self.classname, self)

    def __str__(self):
        return int.__str__(self)

    @property_RO
    def band(self):
        '''Get the UTM/UPS Band (C{'A'|'B'|'C'|'D'..'W'|'X'|'Y'|'Z'} or C{""}).
        '''
        return self._band

    @property_RO
    def hemisphere(self):
        '''Get the UTM/UPS hemisphere (C{'N'|'S'}).
        '''
        return self._hemisphere

    @property_RO
    def utmups(self):
        '''Get the UTM/TMP original (L{Utm}, L{Ups}).
        '''
        return self._utmups

    def utmupsStr(self, B=False):
        '''Get the UTM/TMP zone, band and hemisphere (C{str}).
        '''
        z = '%02d' % (self.zone,)
        b = self.band if B else ''
        h = s = self.hemisphere
        if h:
            s = ' '
        return ''.join(map(str, (z, b, s, h)))

    @property_RO
    def zone(self):
        '''Get the UTM/UPS zone (C{int}, C{1..60} for UTM, C{0} for UPS).
        '''
        return self._zone


def decode2(epsg):
    '''Determine the UTM/USP zone and hemisphere from a given
       U{EPSG<http://www.EPSG-Registry.org>}.

       @param epsg: The EPSG (L{Epsg}, C{str} or C{scalar}).

       @return: 2-Tuple (C{zone, 'N'|'S'}) as (C{int}, C{str})
                where C{zone} is C{1..60} for UTM or C{0} for UPS.

       @raise EPSGError: Invalid I{epsg}.

       @note: Coverage of UPS follows Karney's function U{UTMUPS::DecodeEPSG
              <http://GeographicLib.SourceForge.io/html/classGeographicLib_1_1UTMUPS.html>}.
    '''
    if isinstance(epsg, Epsg):
        return epsg.zone, epsg.hemisphere

    try:
        e = int(epsg)
        if _EPSG_N_01 <= e <= _EPSG_N_60:
            return (e - _EPSG_N_01 + _UTM_ZONE_MIN), 'N'

        elif _EPSG_S_01 <= e <= _EPSG_S_60:
            return (e - _EPSG_S_01 + _UTM_ZONE_MIN), 'S'

        elif e == _EPSG_N:
            return _UPS_ZONE, 'N'

        elif e == _EPSG_S:
            return _UPS_ZONE, 'S'

    except (TypeError, ValueError):
        pass
    raise EPSGError('%s invalid: %r' % ('epsg', epsg))


def encode(zone, hemipole='', band=''):
    '''Determine the U{EPSG<http://www.EPSG-Registry.org>} code for
       a given UTM/UPS zone number, hemisphere/pole and/or Band.

       @param zone: The (longitudinal) UTM zone (C{int}, 1..60) or UPS
                    zone (C{int}, 0) or UTM zone with/-out (latitudinal)
                    Band letter (C{str}, '01C'..'60X') or UPS zone
                    with/-out (polar) Band letter (C{str}, '00A', '00B',
                    '00Y' or '00Z').
       @keyword hemipole: UTM/UPS hemisphere or UPS projection top/center
                          pole (C{str}, C{'N[orth]'} or C{'S[outh]'}).
       @keyword band: Optional (latitudinal) UTM or (polar) UPS Band
                      letter (C{str}).

       @return: C{EPSG} code (L{Epsg}).

       @raise EPSGError: Invalid I{zone}, I{hemipole} or I{band}.

       @note: Coverage of UPS follows Karney's function U{UTMUPS::EncodeEPSG
              <http://GeographicLib.SourceForge.io/html/classGeographicLib_1_1UTMUPS.html>}.
    '''
    try:
        z, B, hp = _to3zBhp(zone, band, hemipole=hemipole)  # in .ellipsoidalBase
        if hp not in ('N', 'S'):
            raise ValueError
    except ValueError:
        raise EPSGError('%s, %s or %s invalid: %r' %
                        ('zone', 'hemipole', 'band', (zone, hp or hemipole, B or band)))

    if _UTM_ZONE_MIN <= z <= _UTM_ZONE_MAX:
        e = z - _UTM_ZONE_MIN + (_EPSG_N_01 if hp == 'N' else _EPSG_S_01)
    elif z == _UPS_ZONE:
        e = _EPSG_N if hp == 'N' else _EPSG_S
    else:
        raise EPSGError('%s invalid: %r' % ('zone', zone))

    e = Epsg(e)
    e._band = B
    # e._hemisphere = hp
    return e

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
