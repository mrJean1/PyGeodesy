
# -*- coding: utf-8 -*-

u'''(INTERNAL) Private L{geodesicx} base class, functions and constants.

Copyright (C) U{Charles Karney<mailto:Karney@Alum.MIT.edu>} (2008-2024)
and licensed under the MIT/X11 License.  For more information, see the
U{GeographicLib<https://GeographicLib.SourceForge.io>} documentation.
'''

from pygeodesy.basics import isodd,  _MODS
from pygeodesy.constants import _EPSmin as _TINY, _0_0, isfinite
from pygeodesy.errors import _or, _xkwds_item2
from pygeodesy.fmath import hypot as _hypot
# from pygeodesy.interns import _numpy_  # _MODS
from pygeodesy.karney import _CapsBase, GeodesicError, _2cos2x, \
                             _norm2, _sincos2d, _sum3
# from pygeodesy.lazily import _ALL_MODS as _MODS  # from .basics

from math import fabs, ldexp as _ldexp

__all__ = ()
__version__ = '25.06.01'

# valid C{nC4}s and C{C4order}s, see _xnC4 below
_nC4s = {24: 2900, 27: 4032, 30: 5425}
# assert (_TINY * EPS) > 0 and (_TINY + EPS) == EPS  # underflow guard


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


class _Gfloats(dict):
    '''(INTERNAL) Numpy or "Unique" floats.
    '''
    n   = 0  # total number of floats
    nC4 = 0

    def __init__(self, nC4):  # PYCHOK signature
        self.nC4 = nC4

    def __call__(self, fs):
        '''Return a C{numpy.array} or C{tuple} of C{float}s.
        '''
        np = _MODS.imported(_MODS.interns._numpy_)
        if np:  # use numpy, already imported
            cs = np.array(fs, dtype=float)
        else:
            _f = self.setdefault  # avoid duplicates
            cs = tuple(_f(f, f) for f in map(float, fs))  # PYCHOK as attr
        self.n += len(fs)
        return cs


def _cosSeries(c4s, sx, cx):  # PYCHOK shared .geodesicx.gx and -.gxline
    '''(INTERNAL) I{Karney}'s cosine series expansion using U{Clenshaw
       summation<https://WikiPedia.org/wiki/Clenshaw_algorithm>}.
    '''
    ar  = _2cos2x(cx, sx)
    y0  =  t0 = y1 = t1 = _0_0
    c4  =  list(c4s)
    _c4 =  c4.pop
    if isodd(len(c4)):
        y0 = _c4()
    while c4:
        # y1 = ar * y0 - y1 + c4.pop()
        # y0 = ar * y1 - y0 + c4.pop()
        y1, t1, _ = _sum3(-y1, -t1, ar * y0, ar * t0, _c4())
        y0, t0, _ = _sum3(-y0, -t0, ar * y1, ar * t1, _c4())
    # s  = (y0 - y1) * cx
    s, t, _ = _sum3(cx * y0, _0_0, cx * t0, -cx * y1, -cx * t1)
    return s + t

#   Y0, Y1 = Fsum(), Fsum()
#   ar  = _2cos2x(cx, sx)
#   c4  =  list(c4s)
#   _c4 =  c4.pop
#   if isodd(len(c4)):
#       Y0 += _c4()
#   while c4:
#       # y1 = ar * y0 - y1 + c4.pop()
#       # y0 = ar * y1 - y0 + c4.pop()
#       Y1 = Y0 * ar - Y1 + _c4()
#       Y0 = Y1 * ar - Y0 + _c4()
#   # s  = (y0 - y1) * cx
#   return float((Y0 - Y1) * cx)


_f = float  # in _f2 and .geodesicx._C4_24, _27 and _30


def _f2(hi, lo):  # in .geodesicx._C4_24, _27 and _30
    '''(INTERNAL) For C{_coeffs}.
    '''
    return _ldexp(_f(hi), 52) + _f(lo)


def _sincos12(sin1, cos1, sin2, cos2, sineg0=False):
    '''(INTERNAL) Compute the sine and cosine of angle
       M{ang12 = atan2(sin2, cos2) - atan2(sin1, cos1)}.

       Negate C{sin1} to get C{sin12} and C{cos12} of the sum
       M{ang12 = atan2(sin2, cos2) + atan2(sin1, cos1)}.

       @kwarg sineg0: If C{True}, make negative C{sin12} zero (C{bool}).

       @return: 2-Tuple C{(sin12, cos12)}.
    '''
    s = sin2 * cos1 - sin1 * cos2
    c = cos2 * cos1 + sin1 * sin2
    if sineg0 and s < 0:
        s = _0_0  # max(s, _0_0) or NEG0?
    return s, c


def _sin1cos2(sin1, cos1, sin2, cos2):
    '''(INTERNAL) Compute the C{sin1 * cos2} sine and its cosine.

       @return: 2-Tuple C{(sin1 * cos2, hypot(sin1 * sin2, cos1)}.
    '''
    s =        sin1 * cos2
    c = _hypot(sin1 * sin2, cos1)
    return s, c  # _norm2(s, c)


def _sinf1cos2d(lat, f1):
    '''(INTERNAL) See C{GeodesicExact} and C{_GeodesicLineExact}.
    '''
    sbet, cbet = _sincos2d(lat)
    # ensure cbet1 = +epsilon at poles; doing the fix on beta means
    # that sig12 will be <= 2*tiny for two points at the same pole
    sbet, cbet = _norm2(sbet * f1, cbet)
    return sbet, (cbet if fabs(cbet) > _TINY else _TINY)


def _toNAN(outmask, *args):
    '''(INTERNAL) Is any C{arg} not finite?
    '''
    if (outmask & _CapsBase.NONFINITONAN):  # Caps.NONFINITONAN
        for arg in args:
            if not isfinite(arg):
                return True
    return False


def _xnC4(**name_nC4):
    '''(INTERNAL) Validate C{C4order}.
    '''
    n, nC4 = _xkwds_item2(name_nC4)
    if nC4 not in _nC4s or not isinstance(nC4, int):
        t = map(str, _nC4s)
        raise GeodesicError(n, nC4, txt_not_=_or(*t))
    return _nC4s[nC4]


# **) MIT License
#
# Copyright (C) 2016-2025 -- mrJean1 at Gmail -- All Rights Reserved.
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
