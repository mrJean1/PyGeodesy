
# -*- coding: utf-8 -*-

u'''Classes L{Epsg} and L{EPSGError} and functions to L{encode} and L{decode2}
I{European Petroleum Survey Group} (U{EPSG<https://www.EPSG-Registry.org>}) codes
from and to U{UTM
<https://WikiPedia.org/wiki/Universal_Transverse_Mercator_coordinate_system>} and
U{UPS<https://WikiPedia.org/wiki/Universal_polar_stereographic_coordinate_system>}
zones.

A pure Python implementation transcribed from I{Charles Karney}'s C++ class
U{UTMUPS<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1UTMUPS.html>},
including coverage of UPS as zone C{0}.
'''

from pygeodesy.basics import isint, isstr, _xinstanceof
from pygeodesy.errors import _ValueError
from pygeodesy.interns import NN, _N_, _NS_, _S_, _SPACE_
from pygeodesy.lazily import _ALL_LAZY, _ALL_OTHER
from pygeodesy.namedTuples import UtmUps2Tuple
from pygeodesy.props import Property_RO
from pygeodesy.streprs import Fmt
from pygeodesy.units import Int
from pygeodesy.ups import Ups
from pygeodesy.utm import Utm
from pygeodesy.utmupsBase import _to3zBhp, _UPS_ZONE, _UTM_ZONE_MIN, \
                                 _UTM_ZONE_MAX, _UTMUPS_ZONE_INVALID

__all__ = _ALL_LAZY.epsg
__version__ = '21.01.20'

# _EPSG_INVALID = _UTMUPS_ZONE_INVALID
_EPSG_N_01 = 32601  # EPSG code for UTM zone 01 N
_EPSG_N_60 = 32660  # EPSG code for UTM zone 60 N
_EPSG_N    = 32661  # EPSG code for UPS pole N

_EPSG_S_01 = 32701  # EPSG code for UTM zone 01 S
_EPSG_S_60 = 32760  # EPSG code for UTM zone 60 S
_EPSG_S    = 32761  # EPSG code for UPS pole S


class Epsg(Int):
    '''U{EPSG<https://www.EPSG-Registry.org>} class, a named C{int}.
    '''
    _band       =  NN
    _epsg       =  None
    _hemisphere =  NN
    _utmups     =  None
    _zone       = _UTMUPS_ZONE_INVALID

    def __new__(cls, eisu, name=NN):
        '''New L{Epsg} (I{European Petroleum Survey Group}) code from a
           UTM/USP coordinate or other EPSG code.

           @arg eisu: Other code (L{Epsg}, C{int}, C{str}, L{Utm} or L{Ups}).

           @return: New L{Epsg}.

           @raise TypeError: Invalid B{C{eisu}}.

           @raise EPSGError: Invalid B{C{eisu}}.
        '''
        if isinstance(eisu, Epsg):
            self = int.__new__(cls, int(eisu))
            self._band       = eisu.band
            self._epsg       = self  # XXX eisu
            self._hemisphere = eisu.hemisphere
            self._utmups     = eisu.utmups
            self._zone       = eisu.zone
            if eisu.name:
                self.name = eisu.name

        elif isint(eisu):
            self = int.__new__(cls, eisu)
            self._epsg = eisu
            self._zone, self._hemisphere = decode2(eisu)  # PYCHOK UtmUps2Tuple

        elif isstr(eisu):
            self = encode(eisu)

        else:
            u = eisu
            _xinstanceof(Utm, Ups, eisu=u)
            self = encode(u.zone, hemipole=u.hemisphere, band=u.band)  # PYCHOK **kwds
            self._utmups = u
            if u.name:
                self.name = u.name

        if name:
            self.name = name
        return self

    def __repr__(self):
        return Fmt.PAREN(self.named, int.__repr__(self))

    def __str__(self):
        return int.__str__(self)

    @Property_RO
    def band(self):
        '''Get the (latitudinal) UTM/UPS Band (C{'A'|'B'|'C'|'D'..'W'|'X'|'Y'|'Z'} or C{""}).
        '''
        return self._band

    @Property_RO
    def hemisphere(self):
        '''Get the UTM/UPS hemisphere/-pole (C{'N'|'S'}).
        '''
        return self._hemisphere

    @Property_RO
    def utmups(self):
        '''Get the UTM/UPS original (L{Utm}, L{Ups}).
        '''
        return self._utmups

    def utmupsStr(self, B=False):
        '''Get the UTM/UPS zone, band and hemisphere/-pole (C{str}).
        '''
        b = self.band if B else NN
        h = s = self.hemisphere
        if h:
            s = _SPACE_
        return NN(Fmt.zone(self.zone), b, s, h)

    @Property_RO
    def zone(self):
        '''Get the (longitudinal) UTM/UPS zone (C{int}, C{1..60} for UTM, C{0} for UPS).
        '''
        return self._zone


class EPSGError(_ValueError):
    '''EPSG encode, decode or other L{Epsg} issue.
    '''
    pass


def decode2(epsg):
    '''Determine the UTM/USP zone and hemisphere from a given
       U{EPSG<https://www.EPSG-Registry.org>}.

       @arg epsg: The EPSG (L{Epsg}, C{str} or C{scalar}).

       @return: A L{UtmUps2Tuple}C{(zone, hemipole)}.

       @raise EPSGError: Invalid B{C{epsg}}.

       @note: Coverage of UPS as zone C{0} follows I{Karney}'s function U{UTMUPS::DecodeEPSG
              <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1UTMUPS.html>}.
    '''
    if isinstance(epsg, Epsg):
        z, h = epsg.zone, epsg.hemisphere

    else:
        try:
            e = int(epsg)  # int(long) OK
            if _EPSG_N_01 <= e <= _EPSG_N_60:
                z, h = int(e - _EPSG_N_01 + _UTM_ZONE_MIN), _N_

            elif _EPSG_S_01 <= e <= _EPSG_S_60:
                z, h = int(e - _EPSG_S_01 + _UTM_ZONE_MIN), _S_

            elif e == _EPSG_N:
                z, h = _UPS_ZONE, _N_

            elif e == _EPSG_S:
                z, h = _UPS_ZONE, _S_

            else:
                raise ValueError
        except (TypeError, ValueError) as x:
            raise EPSGError(epsg=epsg, txt=str(x))

    return UtmUps2Tuple(z, h)


def encode(zone, hemipole=NN, band=NN):
    '''Determine the U{EPSG<https://www.EPSG-Registry.org>} code for
       a given UTM/UPS zone number, hemisphere/pole and/or Band.

       @arg zone: The (longitudinal) UTM zone (C{int}, 1..60) or UPS
                  zone (C{int}, 0) or UTM zone with/-out (latitudinal)
                  Band letter (C{str}, '01C'..'60X') or UPS zone
                  with/-out (polar) Band letter (C{str}, '00A', '00B',
                  '00Y' or '00Z').
       @kwarg hemipole: UTM/UPS hemisphere or UPS projection top/center
                        pole (C{str}, C{'N[orth]'} or C{'S[outh]'}).
       @kwarg band: Optional (latitudinal) UTM or (polar) UPS Band
                    letter (C{str}).

       @return: C{EPSG} code (L{Epsg}).

       @raise EPSGError: Invalid B{C{zone}}, B{C{hemipole}} or B{C{band}}.

       @note: Coverage of UPS as zone C{0} follows I{Karney}'s function U{UTMUPS::EncodeEPSG
              <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1UTMUPS.html>}.
    '''
    try:
        z, B, hp = _to3zBhp(zone, band, hemipole=hemipole)  # in .ellipsoidalBase
        if hp not in _NS_:
            raise ValueError
    except (TypeError, ValueError) as x:
        raise EPSGError(zone=zone, hemipole=hemipole, band=band, txt=str(x))

    if _UTM_ZONE_MIN <= z <= _UTM_ZONE_MAX:
        e = z - _UTM_ZONE_MIN + (_EPSG_N_01 if hp == _N_ else _EPSG_S_01)
    elif z == _UPS_ZONE:
        e = _EPSG_N if hp == _N_ else _EPSG_S
    else:
        raise EPSGError(zone=zone)

    e = Epsg(e)
    e._band = B
    # e._hemisphere = hp
    return e


__all__ += _ALL_OTHER(decode2, encode)

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
