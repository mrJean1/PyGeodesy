
# -*- coding: utf-8 -*-

u'''Base L{geodesicx} classes, functions and constants.

Copyright (C) Charles Karney (2012-2021) <Charles@Karney.com>
and licensed under the MIT/X11 License.  For more information,
see U{GeographicLib<https://GeographicLib.SourceForge.io>}.
'''

# from pygeodesy.basics import isodd  # from .karney
# from pygeodesy.errors import _or  # from .karney
from pygeodesy.interns import MIN as _MIN, _not_, \
                             _0_0, _1_0, _N_1_0, _2_0
from pygeodesy.karney import GeodesicError, isodd, _or, \
                            _NamedBase, Property, _sum2_
from pygeodesy.lazily import _ALL_DOCS
# from pygeodesy.named import _NamedBase  # from .karney
# from pygeodesy.props import Property  # from .karney

from math import sqrt, ldexp as _ldexp

__all__ = ()
__version__ = '21.11.24'

# valid C{nC4}s and C{C4Order}s, see _xnC4 below
_nC4s = {24: 2900, 27: 4032, 30: 5425}
# underflow guard, we require _TINY * EPS > 0, _TINY + EPS == EPS
_TINY = sqrt(_MIN)  # PYCHOK exported


class Caps(object):  # PYCHOK
    '''(INTERNAL) Overriden by C{Caps} below.
    '''
    EMPTY          =  0        # formerly ka NONE
    LATITUDE       =  1 << 7   # compute latitude C{lat2} (0x80)
    LONGITUDE      =  1 << 8   # compute longitude C{lon2}
    AZIMUTH        =  1 << 9   # azimuths C{azi1} and C{azi2}
    DISTANCE       =  1 << 10  # compute distance C{s12}
    DISTANCE_IN    =  1 << 11  # allow distance C{s12} in Direct
    REDUCEDLENGTH  =  1 << 12  # compute reduced length C{m12}
    GEODESICSCALE  =  1 << 13  # compute geodesic scales C{M12} and C{M21}
    AREA           =  1 << 14  # compute area C{S12} (0x4000)

    STANDARD       =  AZIMUTH | DISTANCE | DISTANCE_IN | LATITUDE | LONGITUDE
    ALL            =  0x7F80   # without LONG_UNROLL, REVERSE2 and _DEBUG_*

    LONG_UNROLL    =  1 << 15  # unroll C{lon2} in GeodesicExact.Direct
    REVERSE2       =  1 << 16  # reverse C{azi2}

    _ANGLE_ONLY    =  1 << 17  # angular distance C{a12} only
    _SALPs_CALPs   =  1 << 18  # (INTERNAL) GeodesicExact._GenInverse

    _DEBUG_AREA    =  1 << 19  # (INTERNAL) include Line details
    _DEBUG_DIRECT  =  1 << 20  # (INTERNAL) include Direct details
    _DEBUG_INVERSE =  1 << 21  # (INTERNAL) include Inverse details
    _DEBUG_LINE    =  1 << 22  # (INTERNAL) include Line details
    _DEBUG_ALL     = _DEBUG_AREA | _DEBUG_DIRECT | _DEBUG_INVERSE | \
                     _DEBUG_LINE | _ANGLE_ONLY | _SALPs_CALPs

    _OUT_ALL       =  ALL
    _OUTMASK       =  ALL | LONG_UNROLL | REVERSE2  | _DEBUG_ALL

    _AZIMUTH_DISTANCE                     = AZIMUTH | DISTANCE
    _AZIMUTH_LATITUDE_LONGITUDE           = AZIMUTH | LATITUDE | LONGITUDE
    _LINE                                 = AZIMUTH | LATITUDE | LONG_UNROLL
    _REDUCEDLENGTH_GEODESICSCALE          = REDUCEDLENGTH | GEODESICSCALE
    _DISTANCE_REDUCEDLENGTH_GEODESICSCALE = REDUCEDLENGTH | GEODESICSCALE | DISTANCE

Caps = Caps()  # PYCHOK singleton
'''I{Enum}-style masks to be bit-C{or}'ed to specify geodesic
capabilities (C{caps}) and expected results (C{outmask}).

C{AREA} - compute area C{S12},

C{AZIMUTH} - include azimuths C{azi1} and C{azi2},

C{DISTANCE} - compute distance C{s12},

C{DISTANCE_IN} - allow distance C{s12} in C{.Direct},

C{EMPTY} - nothing, formerly aka C{NONE},

C{GEODESICSCALE} - compute geodesic scales C{M12} and C{M21},

C{LATITUDE} - compute latitude C{lat2},

C{LONGITUDE} - compute longitude C{lon2},

C{LONG_UNROLL} - unroll C{lon2} in C{.Direct},

C{REDUCEDLENGTH} - compute reduced length C{m12},

C{REVERSE2} - reverse C{azi2},

and C{ALL} - all of the above.

C{STANDARD} = C{AZIMUTH | DISTANCE | DISTANCE_IN | LATITUDE | LONGITUDE}'''


class _GeodesicBase(_NamedBase):  # in .geodsolve
    '''(INTERNAL) Base class for C{[_]Geodesic*Exact}.
    '''
    ALL           = Caps.ALL
    AREA          = Caps.AREA
    AZIMUTH       = Caps.AZIMUTH
    DISTANCE      = Caps.DISTANCE
    DISTANCE_IN   = Caps.DISTANCE_IN
    EMPTY         = Caps.EMPTY  # aka NONE
    GEODESICSCALE = Caps.GEODESICSCALE
    LATITUDE      = Caps.LATITUDE
    LONGITUDE     = Caps.LONGITUDE
    LONG_UNROLL   = Caps.LONG_UNROLL
    REDUCEDLENGTH = Caps.REDUCEDLENGTH
    STANDARD      = Caps.STANDARD

    _debug        = 0  # or Caps._DEBUG_...

#   def toRepr(self, prec=6, sep=_COMMASPACE_, **unused):  # PYCHOK signature
#       '''Return this C{GeodesicExact*} items string.
#
#          @kwarg prec: The C{float} precision, number of decimal digits (0..9).
#                       Trailing zero decimals are stripped for B{C{prec}} values
#                       of 1 and above, but kept for negative B{C{prec}} values.
#          @kwarg sep: Optional separator to join (C{str}).
#
#          @return: C{GeodesicExact*} (C{str}).
#       '''
#       return Fmt.PAREN(self.named, self.toStr(prec=prec, sep=sep))

    @Property
    def debug(self):
        '''Get the C{debug} option (C{bool}).
        '''
        return bool(self._debug)

    @debug.setter  # PYCHOK setter!
    def debug(self, debug):
        '''Set the C{debug} option.

           @arg debug: Include more details in results (C{bool}).
        '''
        self._debug = Caps._DEBUG_ALL if debug else 0


def _all_caps(_caps, caps):  # PYCHOK shared
    '''(INTERNAL) Check all available capabilities: C{True}
       if I{all} B{C{caps}} are available in B{C{_caps}},
       C{False} otherwise (C{bool}).
    '''
    caps &= Caps._OUT_ALL
    return (_caps & caps) == caps


def _coSeries(c4s, sx, cx):  # PYCHOK shared .geodesicx.gx and -.gxline
    '''(INTERNAL) I{Karney}'s cosine series expansion using U{Clenshaw
       summation<https://WikiPedia.org/wiki/Clenshaw_algorithm>}.
    '''
    ar = _2_0 * (cx - sx) * (cx + sx)  # 2 * cos(2 * x)
    y0 = t0 = y1 = t1 = _0_0
    i = len(c4s)  # c4s = list(c4s)
    if isodd(i):
        i -= 1
        y0 = c4s[i]  # c4s.pop()
    for i in range(i - 1, 0, -2):  # reversed
        # y1 = ar * y0 - y1 + c4s.pop()
        # y0 = ar * y1 - y0 + c4s.pop()
        y1, t1 = _sum2_(ar * y0, ar * t0, -y1, -t1, c4s[i])
        y0, t0 = _sum2_(ar * y1, ar * t1, -y0, -t0, c4s[i - 1])
    s, _ = _sum2_(_1_0, cx * y0, cx * t0, -cx * y1, -cx * t1, _N_1_0)
    return s  # cx * (y0 - y1)


_f = float  # in _f2 and .geodesicx._C4_24, _27 and _30


def _f2(hi, lo):  # in .geodesicx._C4_24, _27 and _30
    '''(INTERNAL) For C{_coeffs}.
    '''
    return _ldexp(_f(hi), 52) + _f(lo)


def _polynomial(x, c4s, i, j):  # PYCHOK shared
    '''(INTERNAL) Like C{GeographicLib.Math.hpp.polyval} but with a
       different signature and cascaded summation as C{karney._fsum2}.

       @return: M{sum(c4s[k] * x**(j - k - 1) for k in range(i, j)}
    '''
    s, t = c4s[i], _0_0
    for i in range(i + 1, j):
        s, t = _sum2_(s * x, t * x, c4s[i])
    return s  # + t


def _sincos12(sin1, cos1, sin2, cos2, noneg=False):  # PYCHOK shared
    '''(INTERNAL) Compute the C{sin12} and C{cos12} of
       M{ang12 = atan2(sin2, cos2) - atan2(sin1, cos1)}.

       Use C{-sin1} to get C{sin12} and C{cos12} of the sum
       M{ang12 = atan2(sin2, cos2) + atan2(sin1, cos1)}.

       @kwarg noneg: Limit C{sin12} to non-negative (C{bool}).

       @return: 2-Tuple C{(sin12, cos12)}.
    '''
    s = sin2 * cos1 - sin1 * cos2
    c = cos2 * cos1 + sin1 * sin2
    if noneg and s < 0:
        s = _0_0  # max(s, _0_0)
    return s, c


def _xnC4(**name_nC4):
    '''(INTERNAL) Validate C{C4Order}.
    '''
    n, nC4 = name_nC4.popitem()
    if nC4 not in _nC4s or not isinstance(nC4, int):
        raise GeodesicError(n, nC4, txt=_not_(_or(*map(str, _nC4s))))
    return _nC4s[nC4]


__all__ += _ALL_DOCS(Caps)

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
