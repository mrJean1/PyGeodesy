
# -*- coding: utf-8 -*-

u'''Classes L{Georef} and L{WGRSError} and several functions to encode,
decode and inspect I{World Geographic Reference System (WGRS)} references.

Transcribed from I{Charles Karney}'s C++ class U{Georef
<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1Georef.html>},
but with modified C{precision} and extended with C{height} and C{radius}.  See
also U{World Geographic Reference System
<https://WikiPedia.org/wiki/World_Geographic_Reference_System>}.
'''

from pygeodesy.basics import EPS1_2, isstr, property_RO
from pygeodesy.dms import parse3llh  # parseDMS2
from pygeodesy.errors import _ValueError
from pygeodesy.interns import _height_, _item_sq, _Missing, \
                               NN, _prec_, _radius_, _res_
from pygeodesy.lazily import _ALL_LAZY, _ALL_OTHER
from pygeodesy.named import LatLon2Tuple, LatLonPrec3Tuple, \
                            nameof, _xnamed
from pygeodesy.units import Height, Int, Lat, Lon, Precision_, \
                            Radius, Scalar_, Str, _xStrError
from pygeodesy.utily import ft2m, m2ft, m2NM

from math import floor

__all__ = _ALL_LAZY.wgrs
__version__ = '20.07.08'

_Base    = 10
_BaseLen =  4
_DegChar = 'ABCDEFGHJKLMNPQ'
_Digits  = '0123456789'
_LatOrig = -90
_LatTile = 'ABCDEFGHJKLM'
_LonOrig = -180
_LonTile = 'ABCDEFGHJKLMNPQRSTUVWXYZ'
_M_      = 60000000000  # == 60_000_000_000 == 60 * pow(10, 9)
_MaxPrec = 11
_Tile    = 15  # tile size in degrees

_MaxLen  = _BaseLen + 2 * _MaxPrec
_MinLen  = _BaseLen - 2

_LatOrig_M_ = _LatOrig * _M_
_LonOrig_M_ = _LonOrig * _M_

_LatOrig_Tile = float(_LatOrig) / _Tile
_LonOrig_Tile = float(_LonOrig) / _Tile


def _2divmod3(x, Orig_M_):
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
             or geostr[:3] == 'INV' \
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
    '''Georef class, a named c{str}.
    '''
    _height    = _Missing  # meter
    _latlon    =  None     # cached latlon property
    _precision =  None
    _radius    = _Missing  # meter

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
        h = None

        if isinstance(cll, Georef):
            g, p = _2geostr2(str(cll))
            self = Str.__new__(cls, g)
            self._latlon = LatLon2Tuple(*cll._latlon)
            self._precision = p  # cll._precision
            if cll._name:
                self._name = cll._name

        elif isstr(cll):
            if ',' in cll:
                lat, lon, h = _2fllh(*parse3llh(cll, height=None))
                g = encode(lat, lon, precision=precision, height=h)  # PYCHOK false
                self = Str.__new__(cls, g)
                self._latlon = LatLon2Tuple(lat, lon)
            else:
                self = Str.__new__(cls, cll.upper())
                self._decode()

        else:  # assume LatLon
            try:
                lat, lon, h = _2fllh(cll.lat, cll.lon)
                h = getattr(cll, _height_, h)
            except AttributeError:
                raise _xStrError(Georef, cll=cll)  # Error=WGRSError
            g = encode(lat, lon, precision=precision, height=h)  # PYCHOK false
            self = Str.__new__(cls, g)
            self._latlon = LatLon2Tuple(lat, lon)

        if h not in (None, _Missing):
            self._height = Height(h)
        if self._precision is None:
            self._precision = _2Precision(precision)

        if name:
            self.name = name
        return self

    def _decode(self):
        # cache all decoded attrs
        lat, lon, p, h, r = decode5(self)  # PYCHOK LatLonPrec5Tuple
        if self._latlon is None:
            self._latlon = LatLon2Tuple(lat, lon)
        if self._precision is None:
            self._precision = p
        if self._height is _Missing:
            self._height = h
        if self._radius is _Missing:
            self._radius = r

    @property_RO
    def height(self):
        '''Get this georef's height in C{meter} or C{None} if missing.
        '''
        if self._height is _Missing:
            self._decode()
        return self._height

    @property_RO
    def latlon(self):
        '''Get this georef's (center) lat- and longitude (L{LatLon2Tuple}).
        '''
        if self._latlon is None:
            self._decode()
        return self._latlon

    @property_RO
    def precision(self):
        '''Get this georef's precision (C{int}).
        '''
        if self._precision is None:
            self._decode()
        return self._precision

    @property_RO
    def radius(self):
        '''Get this georef's radius in C{meter} or C{None} if missing.
        '''
        if self._radius is _Missing:
            self._decode()
        return self._radius

    def toLatLon(self, LatLon=None, height=None, **LatLon_kwds):
        '''Return (the center of) this georef cell as an instance
           of the supplied C{LatLon} class.

           @kwarg LatLon: Class to use (C{LatLon}) or C{None}.
           @kwarg height: Optional height ({meter}).
           @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                               arguments, ignored if B{C{LatLon=None}}.

           @return: This georef location (B{C{LatLon}}) or if B{C{LatLon}}
                    is C{None}, a L{LatLon3Tuple}C{(lat, lon, height)}.

           @raise TypeError: Invalid B{C{LatLon}} or B{C{LatLon_kwds}}.
        '''
        h = (self.height or 0) if height is None else height
        r = self.latlon.to3Tuple(h) if LatLon is None else \
                 LatLon(*self.latlon, height=h, **LatLon_kwds)
        return self._xnamed(r)


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
        return WGRSError(_item_sq(georef=i), georef)

    def _index(chars, g, i):
        k = chars.find(g[i])
        if k < 0:
            raise _Error(i)
        return k

    g, precision = _2geostr2(georef)
    lon = _index(_LonTile, g, 0) + _LonOrig_Tile
    lat = _index(_LatTile, g, 1) + _LatOrig_Tile

    u = 1.0
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
        lon = lon * 2 + 1
        lat = lat * 2 + 1
        u *= 2
    u = _Tile / u
    r = LatLonPrec3Tuple(Lat(lat * u, Error=WGRSError),
                         Lon(lon * u, Error=WGRSError), precision)
    return _xnamed(r, nameof(georef))


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
        return Height(ft2m(kft * 1000.0), name=name, Error=WGRSError)

    def _r2m(NM, name):
        return Radius(NM / m2NM(1), name=name, Error=WGRSError)

    def _split2(g, name, _2m):
        i = max(g.find(name[0]), g.rfind(name[0]))
        if i > _BaseLen:
            return g[:i], _2m(int(g[i+1:]), '%s %s' % (georef, name))
        else:
            return g, None

    g = Str(georef, Error=WGRSError)

    g, h = _split2(g, 'Height', _h2m)  # H is last
    g, r = _split2(g, 'Radius', _r2m)  # R before H

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
        return '%s%d' % (name[0].upper(), int(m2_(f * K) + 0.5))

    def _pstr(p, x):
        return '%0*d' % (p, x)

    p = _2Precision(precision)

    lat, lon, _ = _2fllh(lat, lon)
    if lat == 90:
        lat *= EPS1_2

    xt, xd, x = _2divmod3(lon, _LonOrig_M_)
    yt, yd, y = _2divmod3(lat, _LatOrig_M_)

    g = _LonTile[xt], _LatTile[yt]
    if p > 0:
        g += _DegChar[xd], _DegChar[yd]
        p -= 1
        if p > 0:
            d = pow(_Base, _MaxPrec - p)
            x = _pstr(p, x // d)
            y = _pstr(p, y // d)
            g += x, y

    if radius is not None:  # R before H
        g += _option(_radius_, radius, m2NM, 1.0),
    if height is not None:  # H is last
        g += _option(_height_, height, m2ft, 1e-3),

    return NN.join(g)  # XXX Georef(''.join(g))


def precision(res):
    '''Determine the L{Georef} precision to meet a required (geographic)
       resolution.

       @arg res: The required resolution (C{degrees}).

       @return: The L{Georef} precision (C{int} 0..11).

       @raise ValueError: Invalid B{C{res}}.

       @see: Function L{wgrs.encode} for more C{precision} details.
    '''
    r = Scalar_(res, name=_res_)
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
    p = Int(prec, name=_prec_, Error=WGRSError)
    if p > 1:
        r = 1.0 / (60.0 * pow(_Base, min(p, _MaxPrec) - 1))
    elif p < 1:
        r = float(_Tile)
    else:
        r = 1.0
    return r


__all__ += _ALL_OTHER(decode3, decode5,  # functions
                      encode, precision, resolution)

# **) MIT License
#
# Copyright (C) 2016-2020 -- mrJean1 at Gmail -- All Rights Reserved.
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
