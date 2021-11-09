
# -*- coding: utf-8 -*-

u'''World Geographic Reference System (WGRS) en-/decoding.

Classes L{Georef} and L{WGRSError} and several functions to encode,
decode and inspect WGRS references.

Transcoded from I{Charles Karney}'s C++ class U{Georef
<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1Georef.html>},
but with modified C{precision} and extended with C{height} and C{radius}.  See
also U{World Geographic Reference System
<https://WikiPedia.org/wiki/World_Geographic_Reference_System>}.
'''

from pygeodesy.basics import isstr
from pygeodesy.dms import parse3llh  # parseDMS2
from pygeodesy.errors import _ValueError, _xkwds
from pygeodesy.interns import EPS1_2, NN, _AtoZnoIO_, _float, \
                             _height_, _radius_, _SPACE_, \
                             _0to9_, _0_5, _0_001, _1_0, _2_0, \
                             _60_0, _90_0, _1000_0
from pygeodesy.lazily import _ALL_LAZY, _ALL_OTHER
from pygeodesy.named import nameof
from pygeodesy.namedTuples import LatLon2Tuple, LatLonPrec3Tuple
from pygeodesy.props import Property_RO
from pygeodesy.streprs import Fmt, _0wd
from pygeodesy.units import Height, Int, Lat, Lon, Precision_, \
                            Radius, Scalar_, Str, _xStrError
from pygeodesy.utily import ft2m, m2ft, m2NM

from math import floor

__all__ = _ALL_LAZY.wgrs
__version__ = '21.07.31'

_Base    =  10
_BaseLen =  4
_DegChar = _AtoZnoIO_.tillQ
_Digits  = _0to9_
_Height_ =  Height.__name__
_INV_    = 'INV'
_LatOrig = -90
_LatTile = _AtoZnoIO_.tillM
_LonOrig = -180
_LonTile = _AtoZnoIO_
_M_      =  60000000000  # == 60_000_000_000 == 60 * pow(10, 9)
_MaxPrec =  11
_Radius_ =  Radius.__name__
_Tile    =  15  # tile size in degrees

_MaxLen  = _BaseLen + 2 * _MaxPrec
_MinLen  = _BaseLen - 2

_LatOrig_M_ = _LatOrig * _M_
_LonOrig_M_ = _LonOrig * _M_

_float_Tile   = _float(_Tile)
_LatOrig_Tile = _float(_LatOrig) / _Tile
_LonOrig_Tile = _float(_LonOrig) / _Tile


def _2divmod3(x, Orig_M_):
    '''(INTERNAL) Convert B{C{x}} to 3_tuple C{(tile, modulo, fraction)}/
    '''
    i = int(floor(x * _M_))
    i, x = divmod(i - Orig_M_, _M_)
    xt, xd = divmod(i, _Tile)
    return xt, xd, x


def _2fllh(lat, lon, height=None):
    '''(INTERNAL) Convert lat, lon, height.
    '''
    # lat, lon = parseDMS2(lat, lon)
    return (Lat(lat, Error=WGRSError),
            Lon(lon, Error=WGRSError), height)


def _2geostr2(georef):
    '''(INTERNAL) Check a georef string.
    '''
    try:
        n, geostr = len(georef), georef.upper()
        p, o = divmod(n, 2)
        if o or n < _MinLen or n > _MaxLen \
             or geostr[:3] == _INV_ \
             or not geostr.isalnum():
            raise ValueError
        return geostr, _2Precision(p - 1)

    except (AttributeError, TypeError, ValueError) as x:
        raise WGRSError(Georef.__name__, georef, txt=str(x))


def _2Precision(precision):
    '''(INTERNAL) Return a L{Precision_} instance.
    '''
    return Precision_(precision, Error=WGRSError, low=0, high=_MaxPrec)


class WGRSError(_ValueError):
    '''World Geographic Reference System (WGRS) encode, decode or other L{Georef} issue.
    '''
    pass


class Georef(Str):
    '''Georef class, a named C{str}.
    '''
    # no str.__init__ in Python 3
    def __new__(cls, cll, precision=3, name=NN):
        '''New L{Georef} from an other L{Georef} instance or georef
           C{str} or from a C{LatLon} instance or lat-/longitude C{str}.

           @arg cll: Cell or location (L{Georef} or C{str}, C{LatLon}
                     or C{str}).
           @kwarg precision: Optional, the desired georef resolution
                             and length (C{int} 0..11), see function
                             L{wgrs.encode} for more details.
           @kwarg name: Optional name (C{str}).

           @return: New L{Georef}.

           @raise RangeError: Invalid B{C{cll}} lat- or longitude.

           @raise TypeError: Invalid B{C{cll}}.

           @raise WGRSError: INValid or non-alphanumeric B{C{cll}}.
        '''
        ll = p = None

        if isinstance(cll, Georef):
            g, p = _2geostr2(str(cll))

        elif isstr(cll):
            if ',' in cll:
                lat, lon, h = _2fllh(*parse3llh(cll, height=None))
                g  = encode(lat, lon, precision=precision, height=h)  # PYCHOK false
                ll = lat, lon  # original lat, lon
            else:
                g = cll.upper()

        else:  # assume LatLon
            try:
                lat, lon, h = _2fllh(cll.lat, cll.lon)
                h = getattr(cll, _height_, h)
            except AttributeError:
                raise _xStrError(Georef, cll=cll)  # Error=WGRSError
            g  = encode(lat, lon, precision=precision, height=h)  # PYCHOK false
            ll = lat, lon  # original lat, lon

        self = Str.__new__(cls, g, name=name or nameof(cll))
        self._latlon    = ll
        self._precision = p
        return self

    @Property_RO
    def decoded3(self):
        '''Get this georef's attributes (L{LatLonPrec3Tuple}).
        '''
        lat, lon = self.latlon
        return LatLonPrec3Tuple(lat, lon, self.precision, name=self.name)

    @Property_RO
    def decoded5(self):
        '''Get this georef's attributes (L{LatLonPrec5Tuple}) with
           height and radius set to C{None} if missing.
        '''
        return self.decoded3.to5Tuple(self.height, self.radius)

    @Property_RO
    def _decoded5(self):
        '''(INTERNAL) Initial L{LatLonPrec5Tuple}.
        '''
        return decode5(self)

    @Property_RO
    def height(self):
        '''Get this georef's height in C{meter} or C{None} if missing.
        '''
        return self._decoded5.height

    @Property_RO
    def latlon(self):
        '''Get this georef's (center) lat- and longitude (L{LatLon2Tuple}).
        '''
        lat, lon = self._latlon or self._decoded5[:2]
        return LatLon2Tuple(lat, lon, name=self.name)

    @Property_RO
    def latlonheight(self):
        '''Get this georef's (center) lat-, longitude and height (L{LatLon3Tuple}),
           with height set to C{0} if missing.
        '''
        return self.latlon.to3Tuple(self.height or 0)

    @Property_RO
    def precision(self):
        '''Get this georef's precision (C{int}).
        '''
        p = self._precision
        return self._decoded5.precision if p is None else p

    @Property_RO
    def radius(self):
        '''Get this georef's radius in C{meter} or C{None} if missing.
        '''
        return self._decoded5.radius

    def toLatLon(self, LatLon=None, height=None, **LatLon_kwds):
        '''Return (the center of) this georef cell as an instance
           of the supplied C{LatLon} class.

           @kwarg LatLon: Class to use (C{LatLon}) or C{None}.
           @kwarg height: Optional height ({meter}).
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                               arguments, ignored if C{B{LatLon} is None}.

           @return: This georef location (B{C{LatLon}}) or if B{C{LatLon}}
                    is C{None}, a L{LatLon3Tuple}C{(lat, lon, height)}.

           @raise TypeError: Invalid B{C{LatLon}} or B{C{LatLon_kwds}}.
        '''
        if LatLon is None:
            r = self.latlonheight if height is None else \
                self.latlon.to3Tuple(height)
        else:
            h = (self.height or 0) if height is None else height
            r =  LatLon(*self.latlon, height=h, **_xkwds(LatLon_kwds, name=self.name))
        return r


def decode3(georef, center=True):
    '''Decode a C{georef} to lat-, longitude and precision.

       @arg georef: To be decoded (L{Georef} or C{str}).
       @kwarg center: If C{True} the center, otherwise the south-west,
                      lower-left corner (C{bool}).

       @return: A L{LatLonPrec3Tuple}C{(lat, lon, precision)}.

       @raise WGRSError: Invalid B{C{georef}}, INValid, non-alphanumeric
                         or odd length B{C{georef}}.
    '''
    def _digit(ll, g, i, m):
        d = _Digits.find(g[i])
        if d < 0 or d >= m:
            raise _Error(i)
        return ll * m + d

    def _Error(i):
        return WGRSError(Fmt.SQUARE(georef=i), georef)

    def _index(chars, g, i):
        k = chars.find(g[i])
        if k < 0:
            raise _Error(i)
        return k

    g, precision = _2geostr2(georef)
    lon = _index(_LonTile, g, 0) + _LonOrig_Tile
    lat = _index(_LatTile, g, 1) + _LatOrig_Tile

    u = _1_0
    if precision > 0:
        lon = lon * _Tile + _index(_DegChar, g, 2)
        lat = lat * _Tile + _index(_DegChar, g, 3)
        m, p = 6, precision - 1
        for i in range(_BaseLen, _BaseLen + p):
            lon = _digit(lon, g, i,     m)
            lat = _digit(lat, g, i + p, m)
            u *= m
            m  = _Base
        u *= _Tile

    if center:
        lon = lon * _2_0 + _1_0
        lat = lat * _2_0 + _1_0
        u *= _2_0
    u = _Tile / u
    return LatLonPrec3Tuple(Lat(lat * u, Error=WGRSError),
                            Lon(lon * u, Error=WGRSError),
                            precision, name=nameof(georef))


def decode5(georef, center=True):
    '''Decode a C{georef} to lat-, longitude, precision, height and radius.

       @arg georef: To be decoded (L{Georef} or C{str}).
       @kwarg center: If C{True} the center, otherwise the south-west,
                      lower-left corner (C{bool}).

       @return: A L{LatLonPrec5Tuple}C{(lat, lon,
                precision, height, radius)} where C{height} and/or
                C{radius} are C{None} if missing.

       @raise WGRSError: Invalid B{C{georef}}, INValid, non-alphanumeric
                         or odd length B{C{georef}}.
    '''
    def _h2m(kft, name):
        return Height(ft2m(kft * _1000_0), name=name, Error=WGRSError)

    def _r2m(NM, name):
        return Radius(NM / m2NM(1), name=name, Error=WGRSError)

    def _split2(g, name, _2m):
        i = max(g.find(name[0]), g.rfind(name[0]))
        if i > _BaseLen:
            return g[:i], _2m(int(g[i+1:]), _SPACE_(georef, name))
        else:
            return g, None

    g = Str(georef, Error=WGRSError)

    g, h = _split2(g, _Height_, _h2m)  # H is last
    g, r = _split2(g, _Radius_, _r2m)  # R before H

    return decode3(g, center=center).to5Tuple(h, r)


def encode(lat, lon, precision=3, height=None, radius=None):  # MCCABE 14
    '''Encode a lat-/longitude as a C{georef} of the given precision.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).
       @kwarg precision: Optional, the desired C{georef} resolution and length
                         (C{int} 0..11).
       @kwarg height: Optional, height in C{meter}, see U{Designation of area
                      <https://WikiPedia.org/wiki/World_Geographic_Reference_System>}.
       @kwarg radius: Optional, radius in C{meter}, see U{Designation of area
                      <https://WikiPedia.org/wiki/World_Geographic_Reference_System>}.

       @return: The C{georef} (C{str}).

       @raise RangeError: Invalid B{C{lat}} or B{C{lon}}.

       @raise WGRSError: Invalid B{C{precision}}, B{C{height}} or B{C{radius}}.

       @note: The B{C{precision}} value differs from U{Georef<https://
              GeographicLib.SourceForge.io/html/classGeographicLib_1_1Georef.html>}.
              The C{georef} length is M{2 * (precision + 1)} and the
              C{georef} resolution is I{15°} for B{C{precision}} 0, I{1°}
              for 1, I{1′} for 2, I{0.1′} for 3, I{0.01′} for 4, ...
              M{10**(2 - precision)}.
    '''
    def _option(name, m, m2_, K):
        f = Scalar_(m, name=name, Error=WGRSError)
        return NN(name[0].upper(), int(m2_(f * K) + _0_5))

    p = _2Precision(precision)

    lat, lon, _ = _2fllh(lat, lon)
    if lat == _90_0:
        lat *= EPS1_2

    xt, xd, x = _2divmod3(lon, _LonOrig_M_)
    yt, yd, y = _2divmod3(lat, _LatOrig_M_)

    g = _LonTile[xt], _LatTile[yt]
    if p > 0:
        g += _DegChar[xd], _DegChar[yd]
        p -= 1
        if p > 0:
            d =  pow(_Base, _MaxPrec - p)
            x = _0wd(p, x // d)
            y = _0wd(p, y // d)
            g += x, y

    if radius is not None:  # R before H
        g += _option(_radius_, radius, m2NM, _1_0),
    if height is not None:  # H is last
        g += _option(_height_, height, m2ft, _0_001),

    return NN.join(g)  # XXX Georef(''.join(g))


def precision(res):
    '''Determine the L{Georef} precision to meet a required (geographic)
       resolution.

       @arg res: The required resolution (C{degrees}).

       @return: The L{Georef} precision (C{int} 0..11).

       @raise ValueError: Invalid B{C{res}}.

       @see: Function L{wgrs.encode} for more C{precision} details.
    '''
    r = Scalar_(res=res)
    for p in range(_MaxPrec):
        if resolution(p) <= r:
            return p
    return _MaxPrec


def resolution(prec):
    '''Determine the (geographic) resolution of a given L{Georef} precision.

       @arg prec: The given precision (C{int}).

       @return: The (geographic) resolution (C{degrees}).

       @raise ValueError: Invalid B{C{prec}}.

       @see: Function L{wgrs.encode} for more C{precision} details.
    '''
    p = Int(prec=prec, Error=WGRSError)
    if p > 1:
        r = _1_0 / (_60_0 * pow(_Base, min(p, _MaxPrec) - 1))
    elif p < 1:
        r = _float_Tile
    else:
        r = _1_0
    return r


__all__ += _ALL_OTHER(decode3, decode5,  # functions
                      encode, precision, resolution)

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
