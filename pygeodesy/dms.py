
# -*- coding: utf-8 -*-

u'''Parsers and formatters of angles in degrees, minutes and seconds.

Functions to parse and format bearing, compass, lat- and longitudes in various
forms of degrees, minutes and seconds with and without symbols and suffix,
including parsing of C{decimal} and C{sexagecimal} degrees.

After I{(C) Chris Veness 2011-2015} published under the same MIT Licence**, see
U{Latitude/Longitude<https://www.Movable-Type.co.UK/scripts/latlong.html>} and
U{Vector-based geodesy<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}.
'''

from pygeodesy.basics import copysign0, isodd, issequence, isstr, map2, neg
from pygeodesy.errors import ParseError, _parseX,  RangeError, \
                            _rangerrors, _ValueError
from pygeodesy.interns import _COMMA_, _NE_, _NSEW_, _NW_, _SE_  # PYCHOK used!
from pygeodesy.interns import NN, _d_, _deg_, _DEGREE_, _degrees_, _DOT_, \
                             _e_, _E_, _f_, _g_, _MINUS_, _PLUSMINUS_, _EW_, \
                             _MINUTE_, _N_, _NS_, _PERCENTDOTSTAR_, _PLUS_, \
                             _radians_, _S_, _SECOND_, _SPACE_, _SW_, _W_, \
                             _0_, _0_5, _60_0, _360_0, _3600_0
from pygeodesy.interns import _QUOTE1_, _QUOTE2_  # PYCHOK used!
from pygeodesy.lazily import _ALL_LAZY, _sys_version_info2  # PYCHOK used!
from pygeodesy.streprs import Fmt, fstr, fstrzs, _0wpF

from math import modf, radians
try:
    from string import letters as _LETTERS
except ImportError:  # Python 3+
    from string import ascii_letters as _LETTERS

__all__ = _ALL_LAZY.dms
__version__ = '21.11.01'

(F_D,    # unsigned format "deg°" plus suffix N, S, E or W
 F_DM,   # unsigned format "deg°min′" plus suffix
 F_DMS,  # unsigned format "deg°min′sec″" plus suffix
 F_DEG,  # unsigned format "[D]DD" plus suffix without symbol
 F_MIN,  # unsigned format "[D]DDMM" plus suffix without symbols
 F_SEC,  # unsigned format "[D]DDMMSS" plus suffix without symbols
 F_D60,  # unsigned format "[D]DD.MMSS" plus suffix without symbols
 F__E,   # unsigned format "%E" plus suffix without symbol
 F__F,   # unsigned format "%F" plus suffix without symbol
 F__G,   # unsigned format "%G" plus suffix without symbol
 F_RAD,  # convert degrees to radians and format unsigned "RR" plus suffix
 ) = _F_s = (_d_, 'dm', 'dms', _deg_, 'min', 'sec', 'd60', _e_, _f_, _g_, 'rad')

(F_D_,    # signed format "-/deg°" without suffix N, S, E or W
 F_DM_ ,  # signed format "-/deg°min′" without suffix
 F_DMS_,  # signed format "-/deg°min′sec″" without suffix
 F_DEG_,  # signed format "-/[D]DD" without suffix and symbol
 F_MIN_,  # signed format "-/[D]DDMM" without suffix and symbols
 F_SEC_,  # signed format "-/[D]DDMMSS" without suffix and symbols
 F_D60_,  # signed format "-/[D]DD.MMSS" without suffix and symbols
 F__E_ ,  # signed format "-/%E" without suffix and symbol
 F__F_,   # signed format "-/%F" without suffix and symbol
 F__G_,   # signed format "-/%G" without suffix and symbol
 F_RAD_,  # convert degrees to radians and format as signed "-/RR" without suffix
 ) = _F_s_ = (NN(_MINUS_, _) for _ in _F_s)  # PYCHOK unused

(F_D__,    # signed format "-/+deg°" without suffix N, S, E or W
 F_DM__,   # signed format "-/+deg°min′" without suffix
 F_DMS__,  # signed format "-/+deg°min′sec″" without suffix
 F_DEG__,  # signed format "-/+[D]DD" without suffix and symbol
 F_MIN__,  # signed format "-/+[D]DDMM" without suffix and symbols
 F_SEC__,  # signed format "-/+[D]DDMMSS" without suffix and symbols
 F_D60__,  # signed format "-/+[D]DD.MMSS" without suffix and symbols
 F__E__,   # signed format "-/+%E" without suffix and symbol
 F__F__,   # signed format "-/+%F" without suffix and symbol
 F__G__,   # signed format "-/+%G" without suffix and symbol
 F_RAD__,  # convert degrees to radians and format signed "-/+RR" without suffix
 ) = _F_s__ = (NN(_PLUS_, _) for _ in _F_s)  # PYCHOK unused

S_DEG = _DEGREE_  # degrees "°" symbol
S_MIN = _MINUTE_  # minutes "′" symbol
S_SEC = _SECOND_  # seconds "″" symbol
S_RAD =  NN   # radians symbol ""
S_SEP =  NN   # separator between deg, min and sec ""
S_NUL =  NN   # empty string

_F_case = {F_D:   F_D,   F_DEG: F_D,   _degrees_: F_D,
           F_DM:  F_DM,  F_MIN: F_DM,  'deg+min': F_DM,
           F_D60: F_D60, F_RAD: F_RAD, _radians_: F_RAD,
           F__E:  F__E,  F__F:  F__F,   F__G:     F__G}

_F_prec = {F_D:   6, F_DM:  4, F_DMS: 2,  # default precs
           F_DEG: 6, F_MIN: 4, F_SEC: 2, F_D60: 0,
           F__E:  8, F__F:  8, F__G:  8, F_RAD: 5}

_F_symb = {F_DEG, F_MIN, F_SEC, F_D60, F__E, F__F, F__G}  # set({}) pychok -Tb

# note ord(_DEGREE_) == ord('°') == 176, ord('˚') == 730
_S_norm = {_DEGREE_: S_DEG, '˚': S_DEG, '^': S_DEG,  # normalized DMS symbols
           _MINUTE_: S_MIN, '’': S_MIN, _QUOTE1_: S_MIN,
           _SECOND_: S_SEC, '”': S_SEC, _QUOTE2_: S_SEC}
assert S_DEG in _S_norm and S_MIN in _S_norm and S_SEC in _S_norm

_beyond_     = 'beyond'
_SEXAGECIMUL =  1.e4  # decimal C{D.MMSS} into sexagecimal C{DMMSS}


def _toDMS(deg, form, prec, sep, ddd, suff):  # MCCABE 15 by .units.py
    '''(INTERNAL) Convert degrees to C{str}, with/-out sign and/or suffix.
    '''
    def _DMS3(d, ddd, p, w):
        d, s = divmod(round(d * _3600_0, p), _3600_0)
        m, s = divmod(s, _60_0)
        return (_0wpF(ddd, 0, d),
                _0wpF(  2, 0, m),
                _0wpF(w+2, p, s))

    try:
        deg = float(deg)
    except (TypeError, ValueError) as x:
        raise _ValueError(deg=deg, txt=str(x))

    form = form.lower()
    sign = form[:1]
    if sign in _PLUSMINUS_:
        form = form[1:]
    else:
        sign = S_NUL

    if prec is None:
        z = p = _F_prec.get(form, 6)
    else:
        z = int(prec)
        p = abs(z)
    w = p + (1 if p else 0)
    d = abs(deg)

    if form in _F_symb:
        s_deg = s_min = s_sec = S_NUL  # no symbols
    else:
        s_deg, s_min, s_sec = S_DEG, S_MIN, S_SEC

    F = _F_case.get(form, F_DMS)
    if F is F_DMS:  # 'deg+min+sec'
        d, m, s = _DMS3(d, ddd, p, w)
        t = NN(d, s_deg, sep,
               m, s_min, sep, s)
        s = s_sec

    elif F is F_DM:  # 'deg+min'
        d, m = divmod(round(d * _60_0, p), _60_0)
        t = NN(_0wpF(ddd, 0, d), s_deg, sep,
               _0wpF(w+2, p, m))
        s = s_min

    elif F is F_D:  # 'deg'
        t = _0wpF(ddd+w, p, d)
        s = s_deg

    elif F is F_RAD:
        t = NN(_PERCENTDOTSTAR_, 'F') % (p, radians(d))
        s = S_RAD

    elif F is F_D60:  # 'deg.MMSSss'
        d, m, s = _DMS3(d, ddd, p, w)
        if z > 1:
            s = fstrzs(s)
            z = 0  # skip fstrzs below
        t = NN(d, _DOT_, m, *s.split(_DOT_))
        s = S_NUL

    else:  # F in (F__E, F__F, F__G)
        t = NN(_PERCENTDOTSTAR_, F) % (p, d)
        s = S_NUL

    if z > 1:
        t = fstrzs(t, ap1z=F is F__G)

    if sign:
        if deg < 0:
            t = _MINUS_ + t
        elif deg > 0 and sign == _PLUS_:
            t = _PLUS_ + t
    elif suff:  # and deg:  # zero suffix?
        s += sep + suff
    return t + s


def bearingDMS(bearing, form=F_D, prec=None, sep=S_SEP):
    '''Convert bearing to a string.

       @arg bearing: Bearing from North (compass C{degrees360}).
       @kwarg form: Optional B{C{bearing}} format (C{str} or L{F_D},
                    L{F_DM}, L{F_DMS}, L{F_DEG}, L{F_MIN}, L{F_SEC},
                    L{F_D60}, L{F__E}, L{F__F}, L{F__G}, L{F_RAD},
                    L{F_D_}, L{F_DM_}, L{F_DMS_}, L{F_DEG_}, L{F_MIN_},
                    L{F_SEC_}, L{F_D60_}, L{F__E_}, L{F__F_}, L{F__G_},
                    L{F_RAD_}, L{F_D__}, L{F_DM__}, L{F_DMS__}, L{F_DEG__},
                    L{F_MIN__}, L{F_SEC__}, L{F_D60__}, L{F__E__},
                    L{F__F__}, L{F__G__} or L{F_RAD__}).
       @kwarg prec: Optional number of decimal digits (0..9 or C{None}
                    for default).  Trailing zero decimals are stripped
                    for B{C{prec}} values of 1 and above, but kept for
                    negative B{C{prec}}.
       @kwarg sep: Optional separator between degrees, minutes, seconds
                   and suffix (C{str}).

       @return: Compass degrees per the specified B{C{form}} (C{str}).

       @JSname: I{toBrng}.
    '''
    return _toDMS(bearing % _360_0, form, prec, sep, 1, NN)


def _clipped_(angle, limit, units):
    '''(INTERNAL) Helper for C{clipDegrees} and C{clipRadians}.
    '''
    c = min(limit, max(-limit, angle))
    if c != angle and _rangerrors:
        t = _SPACE_(fstr(angle, prec=6, ints=True), _beyond_,
                    copysign0(limit, angle), units)
        raise RangeError(t, txt=None)
    return c


def clipDegrees(deg, limit):
    '''Clip a lat- or longitude to the given range.

       @arg deg: Unclipped lat- or longitude (C{degrees}).
       @arg limit: Valid C{-/+B{limit}} range (C{degrees}).

       @return: Clipped value (C{degrees}).

       @raise RangeError: If B{C{deg}} outside the valid C{-/+B{limit}}
                          range and L{pygeodesy.rangerrors} set to C{True}.
    '''
    return _clipped_(deg, limit, _degrees_) if limit and limit > 0 else deg


def clipRadians(rad, limit):
    '''Clip a lat- or longitude to the given range.

       @arg rad: Unclipped lat- or longitude (C{radians}).
       @arg limit: Valid C{-/+B{limit}} range (C{radians}).

       @return: Clipped value (C{radians}).

       @raise RangeError: If B{C{rad}} outside the valid C{-/+B{limit}}
                          range and  L{pygeodesy.rangerrors} set to C{True}.
    '''
    return _clipped_(rad, limit, _radians_) if limit and limit > 0 else rad


def compassDMS(bearing, form=F_D, prec=None, sep=S_SEP):
    '''Convert bearing to a string suffixed with compass point.

       @arg bearing: Bearing from North (compass C{degrees360}).
       @kwarg form: Optional B{C{bearing}} format (C{str} or L{F_D},
                    L{F_DM}, L{F_DMS}, L{F_DEG}, L{F_MIN}, L{F_SEC},
                    L{F_D60}, L{F__E}, L{F__F}, L{F__G}, L{F_RAD},
                    L{F_D_}, L{F_DM_}, L{F_DMS_}, L{F_DEG_}, L{F_MIN_},
                    L{F_SEC_}, L{F_D60_}, L{F__E_}, L{F__F_}, L{F__G_},
                    L{F_RAD_}, L{F_D__}, L{F_DM__}, L{F_DMS__}, L{F_DEG__},
                    L{F_MIN__}, L{F_SEC__}, L{F_D60__}, L{F__E__},
                    L{F__F__}, L{F__G__} or L{F_RAD__}).
       @kwarg prec: Optional number of decimal digits (0..9 or C{None}
                    for default).  Trailing zero decimals are stripped
                    for B{C{prec}} values of 1 and above, but kept for
                    negative B{C{prec}}.
       @kwarg sep: Optional separator between degrees, minutes, seconds
                   and suffix (C{str}).

       @return: Compass degrees and point in the specified form (C{str}).
    '''
    return _toDMS(bearing % _360_0, form, prec, sep, 1, compassPoint(bearing))


def compassPoint(bearing, prec=3):
    '''Convert bearing to a compass point.

       @arg bearing: Bearing from North (compass C{degrees360}).
       @kwarg prec: Optional precision (1 for cardinal or basic winds,
                    2 for intercardinal or ordinal or principal winds,
                    3 for secondary-intercardinal or half-winds or
                    4 for quarter-winds).

       @return: Compass point (1-, 2-, 3- or 4-letter C{str}).

       @raise ValueError: Invalid B{C{prec}}.

       @see: U{Dms.compassPoint
             <https://GitHub.com/ChrisVeness/geodesy/blob/master/dms.js>}
             and U{Compass rose<https://WikiPedia.org/wiki/Compass_rose>}.

       @example:

        >>> p = compassPoint(24, 1)  # 'N'
        >>> p = compassPoint(24, 2)  # 'NE'
        >>> p = compassPoint(24, 3)  # 'NNE'
        >>> p = compassPoint(24)     # 'NNE'
        >>> p = compassPoint(11, 4)  # 'NbE'
        >>> p = compassPoint(30, 4)  # 'NEbN'

        >>> p = compassPoint(11.249)  # 'N'
        >>> p = compassPoint(11.25)   # 'NNE'
        >>> p = compassPoint(-11.25)  # 'N'
        >>> p = compassPoint(348.749) # 'NNW'
    '''
    try:  # m = 2 << prec; x = 32 // m
        m, x = _MOD_X[prec]
    except KeyError:
        raise _ValueError(prec=prec)
    # not round(), i.e. half-even rounding in Python 3+,
    # but round-away-from-zero as int(b + 0.5) iff b is
    # non-negative, otherwise int(b + copysign0(_0_5, b))
    q = int((bearing % _360_0) * m / _360_0 + _0_5) % m
    return _WINDS[q * x]


_MOD_X = {1: (4, 8), 2: (8, 4), 3: (16, 2), 4: (32, 1)}  # [prec]
_WINDS = (_N_, 'NbE', 'NNE', 'NEbN', _NE_, 'NEbE', 'ENE', 'EbN',
          _E_, 'EbS', 'ESE', 'SEbE', _SE_, 'SEbS', 'SSE', 'SbE',
          _S_, 'SbW', 'SSW', 'SWbS', _SW_, 'SWbW', 'WSW', 'WbS',
          _W_, 'WbN', 'WNW', 'NWbW', _NW_, 'NWbN', 'NNW', 'NbW')  # cardinals


def degDMS(deg, prec=6, s_D=S_DEG, s_M=S_MIN, s_S=S_SEC, neg=_MINUS_, pos=NN):
    '''Convert degrees to a string in degrees, minutes I{or} seconds.

       @arg deg: Value in degrees (C{scalar}).
       @kwarg prec: Optional number of decimal digits (0..9 or C{None}
                    for default).  Trailing zero decimals are stripped
                    for B{C{prec}} values of 1 and above, but kept for
                    negative B{C{prec}}.
       @kwarg s_D: Symbol for degrees (C{str}).
       @kwarg s_M: Symbol for minutes (C{str}) or C{""}.
       @kwarg s_S: Symbol for seconds (C{str}) or C{""}.
       @kwarg neg: Optional sign for negative (C{'-'}).
       @kwarg pos: Optional sign for positive (C{''}).

       @return: I{Either} degrees, minutes I{or} seconds (C{str}).
    '''
    try:
        deg = float(deg)
    except (TypeError, ValueError) as x:
        raise _ValueError(deg=deg, txt=str(x))

    d, s = abs(deg), s_D
    if d < 1:
        if s_M:
            d *= _60_0
            if d < 1 and s_S:
                d *= _60_0
                s = s_S
            else:
                s = s_M
        elif s_S:
            d *= _3600_0
            s = s_S

    z = int(prec)
    t = Fmt.F(d, prec=abs(z))
    if z > 1:
        t = fstrzs(t)
    n = neg if deg < 0 else pos
    return NN(n, t, s)


def latDMS(deg, form=F_DMS, prec=2, sep=S_SEP):
    '''Convert latitude to a string, optionally suffixed with N or S.

       @arg deg: Latitude to be formatted (C{degrees}).
       @kwarg form: Optional B{C{deg}} format (C{str} or L{F_D},
                    L{F_DM}, L{F_DMS}, L{F_DEG}, L{F_MIN}, L{F_SEC},
                    L{F_D60}, L{F__E}, L{F__F}, L{F__G}, L{F_RAD},
                    L{F_D_}, L{F_DM_}, L{F_DMS_}, L{F_DEG_}, L{F_MIN_},
                    L{F_SEC_}, L{F_D60_}, L{F__E_}, L{F__F_}, L{F__G_},
                    L{F_RAD_}, L{F_D__}, L{F_DM__}, L{F_DMS__}, L{F_DEG__},
                    L{F_MIN__}, L{F_SEC__}, L{F_D60__}, L{F__E__},
                    L{F__F__}, L{F__G__} or L{F_RAD__}).
       @kwarg prec: Optional number of decimal digits (0..9 or C{None}
                    for default).  Trailing zero decimals are stripped
                    for B{C{prec}} values of 1 and above, but kept for
                    negative B{C{prec}}.
       @kwarg sep: Optional separator between degrees, minutes, seconds
                   and suffix (C{str}).

       @return: Degrees in the specified form (C{str}).

       @JSname: I{toLat}.
    '''
    return _toDMS(deg, form, prec, sep, 2, _S_ if deg < 0 else _N_)


def latlonDMS(lls, form=F_DMS, prec=None, sep=None):
    '''Convert one or more C{LatLon} instances to strings.

       @arg lls: Single or a list, sequence, tuple, etc. (C{LatLon}s).
       @kwarg form: Optional B{C{deg}} format (C{str} or L{F_D},
                    L{F_DM}, L{F_DMS}, L{F_DEG}, L{F_MIN}, L{F_SEC},
                    L{F_D60}, L{F__E}, L{F__F}, L{F__G}, L{F_RAD},
                    L{F_D_}, L{F_DM_}, L{F_DMS_}, L{F_DEG_}, L{F_MIN_},
                    L{F_SEC_}, L{F_D60_}, L{F__E_}, L{F__F_}, L{F__G_},
                    L{F_RAD_}, L{F_D__}, L{F_DM__}, L{F_DMS__}, L{F_DEG__},
                    L{F_MIN__}, L{F_SEC__}, L{F_D60__}, L{F__E__},
                    L{F__F__}, L{F__G__} or L{F_RAD__}).
       @kwarg prec: Optional number of decimal digits (0..9 or C{None}
                    for default).  Trailing zero decimals are stripped
                    for B{C{prec}} values of 1 and above, but kept for
                    negative B{C{prec}}.
       @kwarg sep: Separator joining B{C{lls}} (C{str} or C{None}).

       @return: A C{str} or C{tuple} of B{C{sep}} is C{None} or C{NN}.
    '''
    if issequence(lls):
        t = tuple(ll.toStr(form=form, prec=prec) for ll in lls)
        if sep:
            t = sep.join(t)
    else:
        t = lls.toStr(form=form, prec=prec)
    return t


def lonDMS(deg, form=F_DMS, prec=2, sep=S_SEP):
    '''Convert longitude to a string, optionally suffixed with E or W.

       @arg deg: Longitude to be formatted (C{degrees}).
       @kwarg form: Optional B{C{deg}} format (C{str} or L{F_D},
                    L{F_DM}, L{F_DMS}, L{F_DEG}, L{F_MIN}, L{F_SEC},
                    L{F_D60}, L{F__E}, L{F__F}, L{F__G}, L{F_RAD},
                    L{F_D_}, L{F_DM_}, L{F_DMS_}, L{F_DEG_}, L{F_MIN_},
                    L{F_SEC_}, L{F_D60_}, L{F__E_}, L{F__F_}, L{F__G_},
                    L{F_RAD_}, L{F_D__}, L{F_DM__}, L{F_DMS__}, L{F_DEG__},
                    L{F_MIN__}, L{F_SEC__}, L{F_D60__}, L{F__E__},
                    L{F__F__}, L{F__G__} or L{F_RAD__}).
       @kwarg prec: Optional number of decimal digits (0..9 or C{None}
                    for default).  Trailing zero decimals are stripped
                    for B{C{prec}} values of 1 and above, but kept for
                    negative B{C{prec}}.
       @kwarg sep: Optional separator between degrees, minutes, seconds
                   and suffix (C{str}).

       @return: Degrees in the specified form (C{str}).

       @JSname: I{toLon}.
    '''
    return _toDMS(deg, form, prec, sep, 3, _W_ if deg < 0 else _E_)


if _sys_version_info2 < (3,):  # Python 2-

    def normDMS(strDMS, norm=None):
        '''Normalize all degree, minute and second DMS symbols in a string
           to the default DMS symbols %r, %r and %r.

           @arg strDMS: DMS (C{str}).
           @kwarg norm: Optional replacement symbol (C{str}) or C{None} for
                        the default DMS symbols.  Use C{B{norm}=""} or
                        C{B{norm}=pygeodesy.NN} to remove all DMS symbols.

           @return: Normalized DMS (C{str}).
        '''
        if norm is None:
            for s, S in _S_norm.items():
                strDMS = strDMS.replace(s, S)
        else:
            for s in _S_norm.keys():
                strDMS = strDMS.replace(s, norm)
            strDMS = strDMS.rstrip(norm)
        return strDMS

else:  # Python 3+

    def normDMS(strDMS, norm=None):
        '''Normalize all degree, minute and second DMS symbols in a string
           to the default DMS symbols %r, %r and %r.

           @arg strDMS: DMS (C{str}).
           @kwarg norm: Optional replacement symbol (C{str}) or C{None} for
                        the default DMS symbols.  Use C{B{norm}=""} or
                        C{B{norm}=pygeodesy.NN} to remove all DMS symbols.

           @return: Normalized DMS (C{str}).
        '''
        S_norm = _S_norm

        if norm is None:
            S_norm_get = S_norm.get
            t = NN.join(S_norm_get(s, s) for s in strDMS)

        elif norm:  # replace all DMS symbols
            t = NN.join((norm if s in S_norm else s) for s in strDMS)
            t = t.rstrip(norm)

        else:  # remove all DMS symbols
            t = NN.join(s for s in strDMS if s not in S_norm)

        return t


if __debug__:  # no __doc__ at -O and -OO
    normDMS.__doc__ %= (S_DEG, S_MIN, S_SEC)


def parseDDDMMSS(strDDDMMSS, suffix=_NSEW_, sep=S_SEP, clip=0, sexagecimal=False):
    '''Parse a lat- or longitude represention form [D]DDMMSS in degrees.

       @arg strDDDMMSS: Degrees in any of several forms (C{str}) and
                        types (C{float}, C{int}, other).
       @kwarg suffix: Optional, valid compass points (C{str}, C{tuple}).
       @kwarg sep: Optional separator between "[D]DD", "MM" and "SS" (%r).
       @kwarg clip: Optionally, limit value to range C{-/+B{clip}} (C{degrees}).
       @kwarg sexagecimal: If C{True}, convert C{float(D.MMSS)} to C{base-60}
                           "MM" and "SS" digits.

       @return: Degrees (C{float}).

       @raise ParseError: Invalid B{C{strDDDMMSS}} or B{C{clip}} or the
                          B{C{strDDDMMSS}} form is incompatible with the
                          suffixed or B{C{suffix}} compass point.

       @raise RangeError: Value of B{C{strDDDMMSS}} outside the valid
                          C{-/+B{clip}} range and L{pygeodesy.rangerrors}
                          set to C{True}.

       @note: Type C{str} values "[D]DD", "[D]DDMM" and "[D]DDMMSS" for
              B{C{strDDDMMSS}} are parsed properly only if I{either}
              unsigned and suffixed with a valid, compatible, C{cardinal}
              L{compassPoint} I{or} if unsigned or signed, unsuffixed and
              with keyword argument B{C{suffix}} set to B{%r}, B{%r} or a
              compatible L{compassPoint}.

       @note: Unlike function L{parseDMS}, type C{float}, C{int} and
              other non-C{str} B{C{strDDDMMSS}} values are interpreted
              form [D]DDMMSS.  For example, C{int(1230)} is returned as
              12.5 and I{not 1230.0} degrees.  However, C{int(345)} is
              considered form "DDD" 345 I{and not "DDMM" 0345}, unless
              B{C{suffix}} specifies compass point B{%r}.  Also,
              C{float(15.0523)} is returned as 15.0523 decimal degrees
              and I{not 15°5′23″ sexagecimal}.  To consider the latter,
              use C{float(15.0523)} and C{B{sexagecimal}=True}.

       @see: Functions L{parseDMS}, L{parseDMS2} and L{parse3llh}.
    '''
    def _DDDMMSS_(strDDDMMSS, suffix, sep, clip):
        S = suffix.upper()
        if isstr(strDDDMMSS):
            t = strDDDMMSS.strip()
            if sep:
                t = t.replace(sep, NN).strip()

            s = t[:1]   # sign or digit
            P = t[-1:]  # compass point, digit or dot

            t = t.lstrip(_PLUSMINUS_).rstrip(S).strip()
            f = t.split(_DOT_)
            d = len(f[0])
            f = NN.join(f)
            if 1 < d < 8 and f.isdigit() and (
                                 (P in S and s.isdigit()) or
                            (P.isdigit() and s in '-0123456789+'  # PYCHOK indent
                                         and S in ((_NS_, _EW_) + _WINDS))):
                # check [D]DDMMSS form and compass point
                X = _EW_ if isodd(d) else _NS_
                if not (P in X or (S in X and (P.isdigit() or P == _DOT_))):
                    t = 'DDDMMSS'[0 if isodd(d) else 1:d | 1], X[:1], X[1:]
                    raise ParseError('form %s applies %s-%s' % t)
                f = 0  # fraction
            else:  # try other forms
                return _DMS2deg(strDDDMMSS, S, sep, clip)

        else:  # float or int to [D]DDMMSS[.fff]
            f = float(strDDDMMSS)
            if sexagecimal:
                f *= _SEXAGECIMUL
            s = _MINUS_ if f < 0 else NN
            P = _0_  # anything except _SW_
            f, i = modf(abs(f))
            t = Fmt.f(i, prec=0)  # str(i) == 'i.0'
            d = len(t)
            # bump number of digits to match
            # the given, valid compass point
            if S in (_NS_ if isodd(d) else _EW_):
                t = _0_ + t
                d += 1
            #   P = S
            # elif d > 1:
            #   P = (_EW_ if isodd(d) else _NS_)[0]

        if d < 4:  # [D]DD[.ddd]
            if f:
                t = float(t) + f
            t = t, 0, 0
        else:
            f += float(t[d-2:])
            if d < 6:  # [D]DDMM[.mmm]
                t = t[:d-2], f, 0
            else:  # [D]DDMMSS[.sss]
                t = t[:d-4], t[d-4:d-2], f
        d = _dms2deg(s, P, *map2(float, t))

        return clipDegrees(d, float(clip)) if clip else d

    return _parseX(_DDDMMSS_, strDDDMMSS, suffix, sep, clip,
                              strDDDMMSS=strDDDMMSS, suffix=suffix)


if __debug__:  # no __doc__ at -O and -OO
    parseDDDMMSS.__doc__ %= (S_SEP, _NS_, _EW_, _NS_)


def _dms2deg(s, P, deg, min, sec):
    '''(INTERNAL) Helper for C{parseDDDMMSS} and C{_DMS2deg}.
    '''
    deg += (min + (sec / _60_0)) / _60_0
    if s == _MINUS_ or P in _SW_:
        deg = neg(deg)
    return deg


def _DMS2deg(strDMS, suffix, sep, clip):
    '''(INTERNAL) Helper for C{parseDDDMMSS} and C{parseDMS}.
    '''
    try:
        d = float(strDMS)

    except (TypeError, ValueError):
        strDMS = strDMS.strip()

        t = strDMS.lstrip(_PLUSMINUS_).rstrip(suffix.upper()).strip()
        if sep:  # remove all DMS symbols
            t = normDMS(t.replace(sep, _SPACE_), norm=NN)
        else:  # replace all DMS symbols
            t = normDMS(t, norm=_SPACE_)
        t =  map2(float, t.strip().split()) + (0, 0)
        d = _dms2deg(strDMS[:1], strDMS[-1:], *t[:3])

    return clipDegrees(d, float(clip)) if clip else d


def parseDMS(strDMS, suffix=_NSEW_, sep=S_SEP, clip=0):  # MCCABE 14
    '''Parse a lat- or longitude representation in C{degrees}.

       This is very flexible on formats, allowing signed decimal
       degrees, degrees and minutes or degrees minutes and seconds
       optionally suffixed by a cardinal compass point.

       A variety of symbols, separators and suffixes are accepted,
       for example "3°37′09″W".  Minutes and seconds may be omitted.

       @arg strDMS: Degrees in any of several forms (C{str}) and
                    types (C{float}, C{int}, other).
       @kwarg suffix: Optional, valid compass points (C{str}, C{tuple}).
       @kwarg sep: Optional separator between deg°, min′ and sec″ (C{''}).
       @kwarg clip: Optionally, limit value to range C{-/+B{clip}} (C{degrees}).

       @return: Degrees (C{float}).

       @raise ParseError: Invalid B{C{strDMS}} or B{C{clip}}.

       @raise RangeError: Value of B{C{strDMS}} outside the valid C{-/+B{clip}}
                          range and L{pygeodesy.rangerrors} set to C{True}.

       @note: Unlike function L{parseDDDMMSS}, type C{float}, C{int}
              and other non-C{str} B{C{strDMS}} values are considered
              as decimal degrees.  For example, C{int(1230)} is returned
              as 1230.0 I{and not as 12.5} degrees and C{float(345)} as
              345.0 I{and not as 3.75} degrees!

       @see: Functions L{parseDDDMMSS}, L{parseDMS2} and L{parse3llh}.
    '''
    return _parseX(_DMS2deg, strDMS, suffix, sep, clip, strDMS=strDMS, suffix=suffix)


def parseDMS2(strLat, strLon, sep=S_SEP, clipLat=90, clipLon=180):
    '''Parse a lat- and a longitude representions C{"lat, lon"} in C{degrees}.

       @arg strLat: Latitude in any of several forms (C{str} or C{degrees}).
       @arg strLon: Longitude in any of several forms (C{str} or C{degrees}).
       @kwarg sep: Optional separator between deg°, min′ and sec″ (C{''}).
       @kwarg clipLat: Limit latitude to range C{-/+B{clipLat}} (C{degrees}).
       @kwarg clipLon: Limit longitude to range C{-/+B{clipLon}} (C{degrees}).

       @return: A L{LatLon2Tuple}C{(lat, lon)} in C{degrees}.

       @raise ParseError: Invalid B{C{strLat}} or B{C{strLon}}.

       @raise RangeError: Value of B{C{strLat}} or B{C{strLon}} outside the
                          valid C{-/+B{clipLat}} or C{-/+B{clipLon}} range
                          and L{pygeodesy.rangerrors} set to C{True}.

       @note: See the B{Notes} at function L{parseDMS}.

       @see: Functions L{parseDDDMMSS}, L{parseDMS} and L{parse3llh}.
    '''
    from pygeodesy.namedTuples import LatLon2Tuple  # avoid circluar import

    return LatLon2Tuple(parseDMS(strLat, suffix=_NS_, sep=sep, clip=clipLat),
                        parseDMS(strLon, suffix=_EW_, sep=sep, clip=clipLon))


def parse3llh(strllh, height=0, sep=_COMMA_, clipLat=90, clipLon=180):
    '''Parse a string C{"lat lon [h]"} representing lat-, longitude in
       C{degrees} and optional height in C{meter}.

       The lat- and longitude value must be separated by a separator
       character.  If height is present it must follow, separated by
       another separator.

       The lat- and longitude values may be swapped, provided at least
       one ends with the proper compass point.

       @arg strllh: Latitude, longitude[, height] (C{str}, ...).
       @kwarg height: Optional, default height (C{meter}) or C{None}.
       @kwarg sep: Optional separator between C{"lat lon [h]"} (C{str}).
       @kwarg clipLat: Limit latitude to range C{-/+B{clipLat}} (C{degrees}).
       @kwarg clipLon: Limit longitude to range C{-/+B{clipLon}} (C{degrees}).

       @return: A L{LatLon3Tuple}C{(lat, lon, height)} in C{degrees},
                C{degrees} and C{float}.

       @raise RangeError: Lat- or longitude value of B{C{strllh}} outside
                          the valid C{-/+B{clipLat}} or C{-/+B{clipLon}}
                          range and L{pygeodesy.rangerrors} set to C{True}.

       @raise ValueError: Invalid B{C{strllh}} or B{C{height}}.

       @note: See the B{Notes} at function L{parseDMS}.

       @see: Functions L{parseDDDMMSS}, L{parseDMS} and L{parseDMS2}.

       @example:

        >>> parse3llh('000°00′05.31″W, 51° 28′ 40.12″ N')
        (51.4778°N, 000.0015°W, 0)
    '''
    from pygeodesy.namedTuples import LatLon3Tuple  # avoid circluar import

    def _3llh_(strllh, height, sep):
        ll = strllh.strip().split(sep)
        if len(ll) > 2:  # XXX interpret height unit
            h = float(ll.pop(2).rstrip(_LETTERS).strip())
        else:
            h = height  # None from wgrs.Georef.__new__
        if len(ll) != 2:
            raise ValueError

        a, b = [_.strip() for _ in ll]  # PYCHOK false
        if a[-1:] in _EW_ or b[-1:] in _NS_:
            a, b = b, a
        return LatLon3Tuple(parseDMS(a, suffix=_NS_, clip=clipLat),
                            parseDMS(b, suffix=_EW_, clip=clipLon), h)

    return _parseX(_3llh_, strllh, height, sep, strllh=strllh)


def parseRad(strRad, suffix=_NSEW_, clip=0):
    '''Parse a string representing angle in C{radians}.

       @arg strRad: Degrees in any of several forms (C{str} or C{radians}).
       @kwarg suffix: Optional, valid compass points (C{str}, C{tuple}).
       @kwarg clip: Optionally, limit value to range C{-/+B{clip}} (C{radians}).

       @return: Radians (C{float}).

       @raise ParseError: Invalid B{C{strRad}} or B{C{clip}}.

       @raise RangeError: Value of B{C{strRad}} outside the valid C{-/+B{clip}}
                          range and L{pygeodesy.rangerrors} set to C{True}.
    '''
    def _Rad_(strRad, suffix, clip):
        try:
            r = float(strRad)

        except (TypeError, ValueError):
            strRad = strRad.strip()

            r = float(strRad.lstrip(_PLUSMINUS_).rstrip(suffix.upper()).strip())
            if strRad[:1] == _MINUS_ or strRad[-1:] in _SW_:
                r = neg(r)

        return clipRadians(r, float(clip)) if clip else r

    return _parseX(_Rad_, strRad, suffix, clip, strRad=strRad, suffix=suffix)


def precision(form, prec=None):
    '''Set the default precison for a given F_ form.

       @arg form: L{F_D}, L{F_DM}, L{F_DMS}, L{F_DEG}, L{F_MIN},
                  L{F_SEC}, L{F_D60}, L{F__E}, L{F__F}, L{F__G}
                  or L{F_RAD} (C{str}).
       @kwarg prec: Optional number of decimal digits (0..9 or C{None}
                    for default).  Trailing zero decimals are stripped
                    for B{C{prec}} values of 1 and above, but kept for
                    negative B{C{prec}}.

       @return: Previous precision for the B{C{form}} (C{int}).

       @raise ValueError: Invalid B{C{form}} or B{C{prec}} or B{C{prec}}
                          outside the valid range C{-/+9}.
    '''
    try:
        p = _F_prec[form]
    except KeyError:
        raise _ValueError(form=form)

    if prec is not None:
        from pygeodesy.units import Precision_
        _F_prec[form] = Precision_(prec=prec, low=-9, high=9)

    return p


def toDMS(deg, form=F_DMS, prec=2, sep=S_SEP, ddd=2, neg=_MINUS_, pos=NN):
    '''Convert I{signed} C{degrees} to string, without suffix.

       @arg deg: Degrees to be formatted (C{degrees}).
       @kwarg form: Optional B{C{deg}} format (C{str} or L{F_D},
                    L{F_DM}, L{F_DMS}, L{F_DEG}, L{F_MIN}, L{F_SEC},
                    L{F_D60}, L{F__E}, L{F__F}, L{F__G}, L{F_RAD},
                    L{F_D_}, L{F_DM_}, L{F_DMS_}, L{F_DEG_}, L{F_MIN_},
                    L{F_SEC_}, L{F_D60_}, L{F__E_}, L{F__F_}, L{F__G_},
                    L{F_RAD_}, L{F_D__}, L{F_DM__}, L{F_DMS__}, L{F_DEG__},
                    L{F_MIN__}, L{F_SEC__}, L{F_D60__}, L{F__E__},
                    L{F__F__}, L{F__G__} or L{F_RAD__}).
       @kwarg prec: Optional number of decimal digits (0..9 or C{None}
                    for default).  Trailing zero decimals are stripped
                    for B{C{prec}} values of 1 and above, but kept for
                    negative B{C{prec}}.
       @kwarg sep: Optional separator between degrees, minutes and
                   seconds (C{str}).
       @kwarg ddd: Optional number of digits for deg° (2 or 3).
       @kwarg neg: Optional sign for negative degrees (C{'-'}).
       @kwarg pos: Optional sign for positive degrees (C{''}).

       @return: Degrees in the specified form (C{str}).
    '''
    t = _toDMS(deg, form, prec, sep, ddd, NN)
    if deg and form[:1] not in _PLUSMINUS_:
        t = NN((neg if deg < 0 else (pos if deg > 0 else NN)), t)
    return t

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
