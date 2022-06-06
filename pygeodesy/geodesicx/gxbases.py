
# -*- coding: utf-8 -*-

u'''Base L{geodesicx} classes, functions and constants.

Copyright (C) U{Charles Karney<mailto:Charles@Karney.com>} (2008-2022)
and licensed under the MIT/X11 License.  For more information, see the
U{GeographicLib<https://GeographicLib.SourceForge.io>} documentation.
'''

# from pygeodesy.basics import isodd  # from .karney
# from pygeodesy.errors import _or  # from .karney
from pygeodesy.interns import MIN as _MIN, _not_, _0_0, _2_0
from pygeodesy.karney import _CapsBase, GeodesicError, \
                              isodd, _hypot, _or, _sum2_
# from pygeodesy.props import Property  # from .karney

from math import sqrt, ldexp as _ldexp

__all__ = ()
__version__ = '22.06.01'

# valid C{nC4}s and C{C4order}s, see _xnC4 below
_nC4s = {24: 2900, 27: 4032, 30: 5425}
# underflow guard, we require _TINY * EPS > 0, _TINY + EPS == EPS
_TINY = sqrt(_MIN)  # PYCHOK exported
# assert (_TINY * EPS) > 0 and (_TINY + EPS) == EPS


class _GeodesicBase(_CapsBase):  # in .geodsolve
    '''(INTERNAL) Base class for C{[_]Geodesic*Exact}.
    '''
#   def toRepr(self, prec=6, sep=_COMMASPACE_, **unused):  # PYCHOK signature
#       '''Return this C{GeodesicExact*} items string.
#
#          @kwarg prec: The C{float} precision, number of decimal digits (0..9).
#                       Trailing zero decimals are stripped for B{C{prec}} values
#                       of 1 and above, but kept for negative B{C{prec}} values.
#          @kwarg sep: Separator to join (C{str}).
#
#          @return: C{GeodesicExact*} (C{str}).
#       '''
#       return Fmt.PAREN(self.named, self.toStr(prec=prec, sep=sep))
    pass


def _cosSeries(c4s, sx, cx):  # PYCHOK shared .geodesicx.gx and -.gxline
    '''(INTERNAL) I{Karney}'s cosine series expansion using U{Clenshaw
       summation<https://WikiPedia.org/wiki/Clenshaw_algorithm>}.
    '''
    ar = _2_0 * (cx - sx) * (cx + sx)  # 2 * cos(2 * x)
    y0 = t0 = y1 = t1 = _0_0
    i  = len(c4s)  # c4s = list(c4s)
    if isodd(i):
        i -= 1
        y0 = c4s[i]  # c4s.pop()
    _s2_ = _sum2_
    for i in range(i - 1, 0, -2):  # reversed
        # y1 = ar * y0 - y1 + c4s.pop()
        # y0 = ar * y1 - y0 + c4s.pop()
        y1, t1 = _s2_(ar * y0, ar * t0, -y1, -t1, c4s[i])
        y0, t0 = _s2_(ar * y1, ar * t1, -y0, -t0, c4s[i - 1])
    s, _ = _s2_(cx *  y0, cx * t0, -cx * y1, -cx * t1)
    return s  # cx * (y0  -     y1)


_f = float  # in _f2 and .geodesicx._C4_24, _27 and _30


def _f2(hi, lo):  # in .geodesicx._C4_24, _27 and _30
    '''(INTERNAL) For C{_coeffs}.
    '''
    return _ldexp(_f(hi), 52) + _f(lo)


def _sincos12(sin1, cos1, sin2, cos2, sineg0=False):  # PYCHOK shared
    '''(INTERNAL) Compute the sine and cosine of angle
       M{ang12 = atan2(sin2, cos2) - atan2(sin1, cos1)}.

       Use C{-sin1} to get C{sin12} and C{cos12} of the sum
       M{ang12 = atan2(sin2, cos2) + atan2(sin1, cos1)}.

       @kwarg sineg0: If C{True}, make negative C{sin12} zero (C{bool}).

       @return: 2-Tuple C{(sin12, cos12)}.
    '''
    s = sin2 * cos1 - sin1 * cos2
    c = cos2 * cos1 + sin1 * sin2
    if sineg0 and s < 0:
        s = _0_0  # max(s, _0_0) or NEG0?
    return s, c


def _sin1cos2(sin1, cos1, sin2, cos2):  # PYCHOK shared
    '''(INTERNAL) Compute the C{sin1 * cos2} sine and its cosine.

       @return: 2-Tuple C{(sin1 * cos2, hypot(sin1 * sin2, cos1)}.
    '''
    s =        sin1 * cos2
    c = _hypot(sin1 * sin2, cos1)
    return s, c  # _norm2(s, c)


def _xnC4(**name_nC4):
    '''(INTERNAL) Validate C{C4order}.
    '''
    n, nC4 = name_nC4.popitem()
    if nC4 not in _nC4s or not isinstance(nC4, int):
        raise GeodesicError(n, nC4, txt=_not_(_or(*map(str, _nC4s))))
    return _nC4s[nC4]


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
