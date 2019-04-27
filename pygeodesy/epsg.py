
# -*- coding: utf-8 -*-

u'''Functions L{encodeEPSG} and L{decodeEPSG2} provide en- and decoding
of I{EPSG (European Petroleum Survery Group)} codes from and to U{UTM
<http://WikiPedia.org/wiki/Universal_Transverse_Mercator_coordinate_system>} and
U{UPS<http://WikiPedia.org/wiki/Universal_polar_stereographic_coordinate_system>}
zones.

A pure Python implementation transcribed from C++ class U{UTMUPS
<http://GeographicLib.SourceForge.io/html/classGeographicLib_1_1UTMUPS.html>}
by I{Charles Karney}.
'''

from ellipsoidalBase import _to3zBhp, _UPS_ZONE, \
                            _UTM_ZONE_MIN, _UTM_ZONE_MAX
from lazily import _ALL_LAZY

# all public contants, classes and functions
__all__ = _ALL_LAZY.epsg
__version__ = '19.04.21'

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


def decodeEPSG2(epsg):
    '''Determine the UTM/USP zone number and hemisphere from a
       given I{EPSG (European Petroleum Survery Group)} code.

       @param epsg: The EPSG code (C{str} or C{scalar}).

       @return: 2-Tuple (C{zone, 'N'|'S'}) as (C{int}, C{str})
                where C{zone} is C{1..60} for UTM or C{0} for UPS.

       @raise UTMUPSError: Invalid I{epsg}.
    '''
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


def encodeEPSG(zone, hemipole='', band=''):
    '''Determine the I{EPSG (European Petroleum Survery Group)} code
       for a given UTM/UPS zone number and hemipole or Band.

       @param zone: The (longitudinal) UTM zone (C{int}, 1..60) or
                    zone with/-out (latitudinal) Band letter (C{str},
                    '01C'..'60X') or C{0} for UPS.
       @keyword hemipole: UTM hemisphere or UPS projection top/center
                          pole (C{str}, C{'N[orth]'} or C{'S[outh]'}).
       @keyword band: Optional (latitudinal) UTM Band letter (C{str}).

       @return: C{EPSG} code (C{int}).

       @raise EPSGError: Invalid I{zone}, I{band} or I{pole}.
    '''
    z, B, hp = _to3zBhp(zone, band, hemipole=hemipole)  # in .ellipsoidalBase

    if hp in ('N', 'S'):
        if _UTM_ZONE_MIN <= z <= _UTM_ZONE_MAX:
            return z - _UTM_ZONE_MIN + (_EPSG_N_01 if hp == 'N' else _EPSG_S_01)

        elif z == _UPS_ZONE:
            return _EPSG_N if hp == 'N' else _EPSG_S

        raise EPSGError('%s invalid: %r' % ('zone', zone))
    else:
        raise EPSGError('%s or %s invalid: %r' %
                        ('hemipole', 'band', (hemipole, B or band)))

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
