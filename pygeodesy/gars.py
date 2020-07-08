
# -*- coding: utf-8 -*-

u'''Classes L{Garef} and L{GARSError} and several functions to encode,
decode and inspect I{Global Area Reference System} (GARS) references.

Transcribed from C++ class U{GARS
<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1GARS.html>}
by I{Charles Karney}.  See also U{Global Area Reference System
<https://WikiPedia.org/wiki/Global_Area_Reference_System>} and U{NGA (GARS)
<https://Earth-Info.NGA.mil/GandG/coordsys/grids/gars.html>}.

@newfield example: Example, Examples
'''

from pygeodesy.basics import EPS1_2, isstr, property_RO
from pygeodesy.dms import parse3llh  # parseDMS2
from pygeodesy.errors import _ValueError
from pygeodesy.interns import _item_sq, NN, _prec_, _res_
from pygeodesy.lazily import _ALL_LAZY, _ALL_OTHER
from pygeodesy.named import LatLon2Tuple, LatLonPrec3Tuple, nameof, \
                           _xnamed
from pygeodesy.units import Int, Lat, Lon, Precision_, Scalar_, \
                            Str, _xStrError

from math import floor

__all__ = _ALL_LAZY.gars
__version__ = '20.07.08'

_Digits  = '0123456789'
_LatLen  =    2
_LatOrig =  -90
_Letters = 'ABCDEFGHJKLMNPQRSTUVWXYZ'
_LonLen  =    3
_LonOrig = -180
_MaxPrec =    2

_MinLen = _LonLen + _LatLen
_MaxLen = _MinLen + _MaxPrec

_M1 = 2
_M2 = 2
_M3 = 3
_M_ = _M1 * _M2 * _M3

_LatOrig_M_ = _LatOrig * _M_
_LonOrig_M_ = _LonOrig * _M_

_LatOrig_M1   = _LatOrig * _M1
_LonOrig_M1_1 = _LonOrig * _M1 - 1

_Resolutions = tuple(1.0 / _ for _ in (_M1, _M1 * _M2, _M_))


def _2divmod2(ll, Orig_M_):
    x = int(floor(ll * _M_)) - Orig_M_
    i = (x * _M1) // _M_
    x -= i * _M_ // _M1
    return i, x


def _2fll(lat, lon, *unused):
    '''(INTERNAL) Convert lat, lon.
    '''
    # lat, lon = parseDMS2(lat, lon)
    return (Lat(lat, Error=GARSError),
            Lon(lon, Error=GARSError))


# def _2Garef(garef):
#     '''(INTERNAL) Check or create a L{Garef} instance.
#     '''
#     if not isinstance(garef, Garef):
#         try:
#             garef = Garef(garef)
#         except (TypeError, ValueError):
#             raise _xStrError(Garef, Str, garef=garef)
#     return garef


def _2garstr2(garef):
    '''(INTERNAL) Check a garef string.
    '''
    try:
        n, garstr = len(garef), garef.upper()
        if n < _MinLen or n > _MaxLen \
                       or garstr[:3] == 'INV' \
                       or not garstr.isalnum():
            raise ValueError
        return garstr, _2Precision(n - _MinLen)

    except (AttributeError, TypeError, ValueError) as x:
        raise GARSError(Garef.__name__, garef, txt=str(x))


def _2Precision(precision):
    '''(INTERNAL) Return a L{Precision_} instance.
    '''
    return Precision_(precision, Error=GARSError, low=0, high=_MaxPrec)


class GARSError(_ValueError):
    '''Global Area Reference System (GARS) encode, decode or other L{Garef} issue.
    '''
    pass


class Garef(Str):
    '''Garef class, a named C{str}.
    '''
    _latlon    = None  # cached latlon property
    _precision = None

    # no str.__init__ in Python 3
    def __new__(cls, cll, precision=1, name=NN):
        '''New L{Garef} from an other L{Garef} instance or garef
           C{str} or from a C{LatLon} instance or lat-/longitude C{str}.

           @arg cll: Cell or location (L{Garef} or C{str}, C{LatLon}
                     or C{str}).
           @kwarg precision: Optional, the desired garef resolution
                             and length (C{int} 0..2), see function
                             L{gars.encode} for more details.
           @kwarg name: Optional name (C{str}).

           @return: New L{Garef}.

           @raise RangeError: Invalid B{C{cll}} lat- or longitude.

           @raise TypeError: Invalid B{C{cll}}.

           @raise GARSError: INValid or non-alphanumeric B{C{cll}}.
        '''
        if isinstance(cll, Garef):
            g, p = _2garstr2(str(cll))
            self = Str.__new__(cls, g)
            self._latlon = LatLon2Tuple(*cll._latlon)
            self._name = cll._name
            self._precision = p  # cll._precision

        elif isstr(cll):
            if ',' in cll:
                lat, lon = _2fll(*parse3llh(cll))
                cll = encode(lat, lon, precision=precision)  # PYCHOK false
                self = Str.__new__(cls, cll)
                self._latlon = LatLon2Tuple(lat, lon)
            else:
                self = Str.__new__(cls, cll.upper())
                self._decode()

        else:  # assume LatLon
            try:
                lat, lon = _2fll(cll.lat, cll.lon)
            except AttributeError:
                raise _xStrError(Garef, cll=cll)  # Error=GARSError
            cll = encode(lat, lon, precision=precision)  # PYCHOK false
            self = Str.__new__(cls, cll)
            self._latlon = LatLon2Tuple(lat, lon)

        if self._precision is None:
            self._precision = _2Precision(precision)

        if name:
            self.name = name
        return self

    def _decode(self):
        # cache all decoded attrs
        lat, lon, p = decode3(self)
        if self._latlon is None:
            self._latlon = LatLon2Tuple(lat, lon)
        if self._precision is None:
            self._precision = p

    @property_RO
    def latlon(self):
        '''Get this garef's (center) lat- and longitude (L{LatLon2Tuple}).
        '''
        if self._latlon is None:
            self._decode()
        return self._latlon

    @property_RO
    def precision(self):
        '''Get this garef's precision (C{int}).
        '''
        if self._precision is None:
            self._decode()
        return self._precision

    def toLatLon(self, LatLon, **kwds):
        '''Return (the center of) this garef cell as an instance
           of the supplied C{LatLon} class.

           @arg LatLon: Class to use (C{LatLon}).
           @kwarg kwds: Optional keyword arguments for B{C{LatLon}}.

           @return: This garef location (B{C{LatLon}}).

           @raise GARSError: Invalid B{C{LatLon}}.
        '''
        if LatLon is None:
            raise GARSError(LatLon=LatLon)

        return self._xnamed(LatLon(*self.latlon, **kwds))


def decode3(garef, center=True):
    '''Decode a C{garef} to lat-, longitude and precision.

       @arg garef: To be decoded (L{Garef} or C{str}).
       @kwarg center: If C{True} the center, otherwise the south-west,
                      lower-left corner (C{bool}).

       @return: A L{LatLonPrec3Tuple}C{(lat, lon, precision)}.

       @raise GARSError: Invalid B{C{garef}}, INValid, non-alphanumeric
                         or bad length B{C{garef}}.
    '''
    def _Error(i):
        return GARSError(garef=_item_sq(repr(garef), i))

    def _ll(chars, g, i, j, lo, hi):
        ll, b = 0, len(chars)
        for i in range(i, j):
            d = chars.find(g[i])
            if d < 0:
                raise _Error(i)
            ll = ll * b + d
        if ll < lo or ll > hi:
            raise _Error(j)
        return ll

    def _ll2(lon, lat, g, i, m):
        d = _Digits.find(g[i])
        if d < 1 or d > m * m:
            raise _Error(i)
        d, r = divmod(d - 1, m)
        lon = lon * m + r
        lat = lat * m + (m - 1 - d)
        return lon, lat

    g, precision = _2garstr2(garef)

    lon = _ll(_Digits,  g,       0, _LonLen, 1, 720) + _LonOrig_M1_1
    lat = _ll(_Letters, g, _LonLen, _MinLen, 0, 359) + _LatOrig_M1
    if precision > 0:
        lon, lat = _ll2(lon, lat, g, _MinLen, _M2)
        if precision > 1:
            lon, lat = _ll2(lon, lat, g, _MinLen + 1, _M3)

    if center:  # ll = (ll * 2 + 1) / 2
        lon += 0.5
        lat += 0.5

    r = _Resolutions[precision]  # == 1.0 / unit
    r = LatLonPrec3Tuple(Lat(lat * r, Error=GARSError),
                         Lon(lon * r, Error=GARSError), precision)
    return _xnamed(r, nameof(garef))


def encode(lat, lon, precision=1):  # MCCABE 14
    '''Encode a lat-/longitude as a C{garef} of the given precision.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).
       @kwarg precision: Optional, the desired C{garef} resolution
                         and length (C{int} 0..2).

       @return: The C{garef} (C{str}).

       @raise RangeError: Invalid B{C{lat}} or B{C{lon}}.

       @raise GARSError: Invalid B{C{precision}}.

       @note: The C{garef} length is M{precision + 5} and the C{garef}
              resolution is B{30′} for B{C{precision}} 0, B{15′} for 1
              and B{5′} for 2, respectively.
    '''
    def _digit(x, y, m):
        return _Digits[m * (m - y - 1) + x + 1],

    def _str(chars, x, n):
        s, b = [], len(chars)
        for i in range(n):
            x, i = divmod(x, b)
            s.append(chars[i])
        return tuple(reversed(s))

    p = _2Precision(precision)

    lat, lon = _2fll(lat, lon)
    if lat == 90:
        lat *= EPS1_2

    ix, x = _2divmod2(lon, _LonOrig_M_)
    iy, y = _2divmod2(lat, _LatOrig_M_)

    g = _str(_Digits, ix + 1, _LonLen) + _str(_Letters, iy, _LatLen)
    if p > 0:
        ix, x = divmod(x, _M3)
        iy, y = divmod(y, _M3)
        g += _digit(ix, iy, _M2)
        if p > 1:
            g += _digit(x, y, _M3)

    return NN.join(g)


def precision(res):
    '''Determine the L{Garef} precision to meet a required (geographic)
       resolution.

       @arg res: The required resolution (C{degrees}).

       @return: The L{Garef} precision (C{int} 0..2).

       @raise ValueError: Invalid B{C{res}}.

       @see: Function L{gars.encode} for more C{precision} details.
    '''
    r = Scalar_(res, name=_res_)
    for p in range(_MaxPrec):
        if resolution(p) <= r:
            return p
    return _MaxPrec


def resolution(prec):
    '''Determine the (geographic) resolution of a given L{Garef} precision.

       @arg prec: The given precision (C{int}).

       @return: The (geographic) resolution (C{degrees}).

       @raise ValueError: Invalid B{C{prec}}.

       @see: Function L{gars.encode} for more C{precision} details.
    '''
    p = Int(prec, name=_prec_, Error=GARSError)
    return _Resolutions[max(0, min(p, _MaxPrec))]


__all__ += _ALL_OTHER(decode3,  # functions
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
