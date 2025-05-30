
# -*- coding: utf-8 -*-

u'''Parsers and formatters of angles in degrees, minutes and seconds or radians.

Functions to parse and format bearing, compass, lat- and longitudes in various forms of
degrees, minutes and seconds with or without degrees, minute and second symbols plus a
compass point suffix, including parsing of C{decimal} and C{sexagecimal} degrees.

Set env variable C{PYGEODESY_FMT_FORM} to any C{F_...} form to override default C{F_DMS}
formatting of lat- and longitudes or to an empty string to restore the default.

After I{(C) Chris Veness 2011-2024} published under the same MIT Licence**, see
U{Latitude/Longitude<https://www.Movable-Type.co.UK/scripts/latlong.html>} and
U{Vector-based geodesy<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}.

@var F_D:   Format degrees as unsigned "deg°" with symbol, plus compass point suffix C{N, S, E} or C{W} (C{str}).
@var F_DM:  Format degrees as unsigned "deg°min′" with symbols, plus suffix (C{str}).
@var F_DMS: Format degrees as unsigned "deg°min′sec″" with symbols, plus suffix (C{str}).
@var F_DEG: Format degrees as unsigned "[D]DD" I{without} symbol, plus suffix (C{str}).
@var F_MIN: Format degrees as unsigned "[D]DDMM" I{without} symbols, plus suffix (C{str}).
@var F_SEC: Format degrees as unsigned "[D]DDMMSS" I{without} symbols, plus suffix (C{str}).
@var F_D60: Format degrees as unsigned "[D]DD.MMSS" C{sexagecimal} I{without} symbols, plus suffix (C{str}).
@var F__E:  Format degrees as unsigned "%E" I{without} symbols, plus suffix (C{str}).
@var F__F:  Format degrees as unsigned "%F" I{without} symbols, plus suffix (C{str}).
@var F__G:  Format degrees as unsigned "%G" I{without} symbols, plus suffix (C{str}).
@var F_RAD: Convert degrees to radians and format as unsigned "RR" with symbol, plus suffix (C{str}).

@var F_D_:   Format degrees as signed "-/deg°" with symbol, I{without} suffix (C{str}).
@var F_DM_:  Format degrees as signed "-/deg°min′" with symbols, I{without} suffix (C{str}).
@var F_DMS_: Format degrees as signed "-/deg°min′sec″" with symbols, I{without} suffix (C{str}).
@var F_DEG_: Format degrees as signed "-/[D]DD" I{without} symbol, I{without} suffix (C{str}).
@var F_MIN_: Format degrees as signed "-/[D]DDMM" I{without} symbols, I{without} suffix (C{str}).
@var F_SEC_: Format degrees as signed "-/[D]DDMMSS" I{without} symbols, I{without} suffix (C{str}).
@var F_D60_: Format degrees as signed "-/[D]DD.MMSS" C{sexagecimal} I{without} symbols, I{without} suffix (C{str}).
@var F__E_:  Format degrees as signed "-/%E" I{without} symbols, I{without} suffix (C{str}).
@var F__F_:  Format degrees as signed "-/%F" I{without} symbols, I{without} suffix (C{str}).
@var F__G_:  Format degrees as signed "-/%G" I{without} symbols, I{without} suffix (C{str}).
@var F_RAD_: Convert degrees to radians and format as signed "-/RR" I{without} symbol, I{without} suffix (C{str}).

@var F_D__:   Format degrees as signed "-/+deg°" with symbol, I{without} suffix (C{str}).
@var F_DM__:  Format degrees as signed "-/+deg°min′" with symbols, I{without} suffix (C{str}).
@var F_DMS__: Format degrees as signed "-/+deg°min′sec″" with symbols, I{without} suffix (C{str}).
@var F_DEG__: Format degrees as signed "-/+[D]DD" I{without} symbol, I{without} suffix (C{str}).
@var F_MIN__: Format degrees as signed "-/+[D]DDMM" I{without} symbols, without suffix (C{str}).
@var F_SEC__: Format degrees as signed "-/+[D]DDMMSS" I{without} symbols, I{without} suffix (C{str}).
@var F_D60__: Format degrees as signed "-/+[D]DD.MMSS" C{sexagecimal} I{without} symbols, I{without} suffix (C{str}).
@var F__E__:  Format degrees as signed "-/+%E" I{without} symbols, I{without} suffix (C{str}).
@var F__F__:  Format degrees as signed "-/+%F" I{without} symbols, I{without} suffix (C{str}).
@var F__G__:  Format degrees as signed "-/+%G" I{without} symbols, I{without} suffix (C{str}).
@var F_RAD__: Convert degrees to radians and format as signed "-/+RR" I{without} symbol, I{without} suffix (C{str}).

@var S_DEG: Degrees symbol, default C{"°"}
@var S_MIN: Minutes symbol, default C{"′"} aka I{PRIME}
@var S_SEC: Seconds symbol, default C{"″"} aka I{DOUBLE_PRIME}
@var S_RAD: Radians symbol, default C{""} aka L{pygeodesy.NN}
@var S_DMS: If C{True} include, otherwise cancel all DMS symbols, default C{True}.
@var S_SEP: Separator between C{deg°|min′|sec″|suffix}, default C{""} aka L{pygeodesy.NN}

@note: In Python 2-, L{S_DEG}, L{S_MIN}, L{S_SEC}, L{S_RAD} and L{S_SEP} may be multi-byte,
       non-ascii characters and if so, I{not} C{unicode}.
'''

from pygeodesy.basics import copysign0, isLatLon, isodd, issequence, isstr, \
                             neg as _neg,  typename  # in .ups
from pygeodesy.constants import _umod_360, _0_0, _0_5, _60_0, _360_0, _3600_0
from pygeodesy.errors import ParseError, RangeError, _TypeError, _ValueError, \
                            _parseX, rangerrors, _xError, _xkwds,  _envPYGEODESY
# from pygeodesy.internals import _envPYGEODESY, typename  # from .errors
from pygeodesy.interns import NN, _COMMA_, _d_, _DASH_, _deg_, _degrees_, _DOT_, \
                             _0_, _e_, _E_, _EW_, _f_, _F_, _g_, _MINUS_, _N_, \
                             _NE_, _NS_, _NSEW_, _NW_, _PERCENTDOTSTAR_, _PLUS_, \
                             _PLUSMINUS_, _QUOTE1_, _QUOTE2_, _radians_, _S_, \
                             _SE_, _SPACE_, _SW_, _W_
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS
# from pygeodesy.namedTuples import LatLon2Tuple  # _MODS
# from pygeodesy.props import _throwarning  # _MODS
from pygeodesy.streprs import Fmt, fstr, fstrzs, _0wpF
# from pygeodesy.units import Precision_  # _MODS
# from pygeodesy.utily import _Wrap  # _MODS

from math import fabs, modf, radians
try:
    from string import letters as _LETTERS
except ImportError:  # Python 3+
    from string import ascii_letters as _LETTERS

__all__ = _ALL_LAZY.dms
__version__ = '25.04.14'

_beyond_      = 'beyond'
_deg_min_     = 'deg+min'
_sexagecimal_ = 'sexagecimal'
_SEXAGECIMUL  =  1.e4  # sexagecimal C{D.MMSSss} into decimal C{DMMSS.ss}

F_D,   F_DM,   F_DMS,   F_DEG,   F_MIN,   F_SEC,   F_D60,   F__E,   F__F,   F__G,   F_RAD   = _F_s = (
 _d_,   'dm',   'dms',   _deg_,   'min',   'sec',   'd60',    _e_,    _f_,    _g_,   'rad')
F_D_,  F_DM_,  F_DMS_,  F_DEG_,  F_MIN_,  F_SEC_,  F_D60_,  F__E_,  F__F_,  F__G_,  F_RAD_  = (NN(
 _MINUS_, _) for _ in _F_s)
F_D__, F_DM__, F_DMS__, F_DEG__, F_MIN__, F_SEC__, F_D60__, F__E__, F__F__, F__G__, F_RAD__ = (NN(
 _PLUS_, _)  for _ in _F_s)
del _F_s
_F_DMS  = _envPYGEODESY('FMT_FORM') or F_DMS

_F_case = {F_D:   F_D,   F_DEG: F_D,   _degrees_: F_D,  # unsigned _F_s
           F_DM:  F_DM,  F_MIN: F_DM,  _deg_min_: F_DM,
           F_D60: F_D60, F_RAD: F_RAD, _radians_: F_RAD,
           F__E:  F__E,  F__F:  F__F,   F__G:     F__G}  # default F_DMS
_F_prec = {F_D:   6, F_DM:  4, F_DMS: 2,  # default precs
           F_DEG: 6, F_MIN: 4, F_SEC: 2, F_D60: 0,
           F__E:  8, F__F:  8, F__G:  8, F_RAD: 5}
_F_symb = set((F_D, F_DM, F_DMS, _deg_min_))  # == {} pychok -Tb

S_DEG = _DEGREES_ = '°'  # ord() = 176
S_MIN = _MINUTES_ = '′'  # PRIME
S_SEC = _SECONDS_ = '″'  # DOUBLE_PRIME
S_RAD = _RADIANS_ = NN   # PYCHOK radians symbol ""
S_DMS =  True  # include DMS symbols
S_SEP =  NN    # separator between deg|min|sec|suffix ""
S_NUL =  NN    # empty string, kept INTERNAL

# note: ord(_DEGREES_) == ord('°') == 176, ord('˚') == 730
_S_norm = {S_DEG: _DEGREES_, '˚': _DEGREES_, '^': _DEGREES_,  # _d_: _DEGREES_,
           S_MIN: _MINUTES_, '’': _MINUTES_, _QUOTE1_: _MINUTES_,  # _r_: _RADIANS_
           S_SEC: _SECONDS_, '”': _SECONDS_, _QUOTE2_: _SECONDS_}

_WINDS = (_N_, 'NbE', 'NNE', 'NEbN', _NE_, 'NEbE', 'ENE', 'EbN',
          _E_, 'EbS', 'ESE', 'SEbE', _SE_, 'SEbS', 'SSE', 'SbE',
          _S_, 'SbW', 'SSW', 'SWbS', _SW_, 'SWbW', 'WSW', 'WbS',
          _W_, 'WbN', 'WNW', 'NWbW', _NW_, 'NWbN', 'NNW', 'NbW')


def _D603(sep, s_D=_DOT_, s_M=None, s_S=S_NUL, s_DMS=S_DMS, **unused):
    '''(INTERNAL) Get the overridden or default pseudo-C{DMS} symbols.
    '''
    if s_DMS:
        M = sep if s_M is None else s_M
        return s_D, (M or S_NUL), s_S
    else:  # no overriden symbols
        return _DOT_, sep, S_NUL


def _DMS3(form, s_D=S_DEG, s_M=S_MIN, s_S=S_SEC, s_DMS=S_DMS, **unused):
    '''(INTERNAL) Get the overridden or default C{DMS} symbols.
    '''
    return (s_D, s_M, s_S) if s_DMS and form in _F_symb else (S_NUL, S_NUL, S_NUL)


def _dms3(d, ddd, p, w):
    '''(INTERNAL) Format C{d} as (deg, min, sec) C{str}s with leading zeros.
    '''
    d =    round(d * _3600_0, p)
    d, s = divmod(d, _3600_0)
    m, s = divmod(s, _60_0)
    return (_0wpF(ddd, 0, d),
            _0wpF(  2, 0, m),
            _0wpF(w+2, p, s))


def _DR2(s_D=S_NUL, s_R=S_RAD, **unused):
    '''(INTERNAL) Get the overridden or default C{D} and C{RAD} symbols.
    '''
    return s_D, s_R


def _fstrzs(t, **unused):
    '''(INTERNAL) Pass-thru version of C{.streprs.fstrzs}.
    '''
    return t


def _split3(strDMS, suffix=_NSEW_):
    '''(INTERNAL) Return sign, stripped B{C{strDMS}} and compass point.
    '''
    t = strDMS.strip()
    s = t[:1]   # sign or digit
    P = t[-1:]  # compass point or digit or dot
    t = t.lstrip(_PLUSMINUS_).rstrip(suffix)
    return s, t.strip(), P


def _toDMS(deg, form, prec, sep, ddd, P, s_D_M_S):  # MCCABE 15 in .units
    '''(INTERNAL) Convert C{deg} to C{str}, with/-out sign, DMS symbols and/or suffix.
    '''
    try:
        F   = None
        deg = float(deg)

        f, s = form, form[:1]
        if s in _PLUSMINUS_:  # signed
            n = _MINUS_ if deg < 0 else (
                _PLUS_  if deg > 0 and s == _PLUS_ else NN)
            P = NN  # no suffix
            f = f.lstrip(_PLUSMINUS_)
        else:  # suffixed
            n = NN
            if sep and P:  # no sep if no suffix
                P = NN(sep, P)
        try:
            F = _F_case[f]  # .strip()
        except KeyError:
            f =  f.lower()  # .strip()
            F = _F_case.get(f, F_DMS)

        if prec is None:
            z = p = _F_prec.get(F, 6)
        else:
            z = int(prec)
            p = abs(z)
        w = p + (1 if p else 0)
        z = fstrzs if z > 1 else _fstrzs
        d = fabs(deg)

        if F is F_DMS:  # 'deg+min+sec', default
            D, M, S = _DMS3(f, **s_D_M_S)
            d, m, s = _dms3(d, ddd, p, w)
            t = NN(n, d,  D, sep,
                      m,  M, sep,
                    z(s), S, P)

        elif F is F_DM:  # 'deg+min'
            D, M, _ = _DMS3(f, **s_D_M_S)
            d, m = divmod(round(d * _60_0, p), _60_0)
            t = NN(n, _0wpF(ddd, 0, d),  D, sep,
                    z(_0wpF(w+2, p, m)), M, P)

        elif F is F_D:  # 'deg'
            D, _, _ = _DMS3(f, **s_D_M_S)
            t = NN(n, z(_0wpF(w+ddd, p, d)), D, P)

        elif F is F_D60:  # 'deg.MM|SSss|'
            D, M, S = _D603(sep, **s_D_M_S)
            d, m, s = _dms3(d, ddd, p, w)
            t = z(s).split(_DOT_) + [S, P]
            t = NN(n, d, D, m, M, *t)

        elif F is F_RAD:
            _, R = _DR2(**s_D_M_S)
            r = NN(_PERCENTDOTSTAR_, _F_) % (p, radians(d))
            t = NN(n, z(r), R, P)

        else:  # F in (F__E, F__F, F__G)
            D, _ = _DR2(**s_D_M_S)
            d = NN(_PERCENTDOTSTAR_, F) % (p, d)  # XXX f?
            t = NN(n, z(d, ap1z=F is F__G), D, P)

    except Exception as x:
        d = {} if F is None else dict(F_=F, P=P)
        raise _xError(x, deg=deg, form=form, prec=prec, **d)

    return t  # NOT unicode in Python 2-


def bearingDMS(bearing, form=F_D, prec=None, sep=S_SEP, **s_D_M_S):
    '''Convert bearing to a string (without compass point suffix).

       @arg bearing: Bearing from North (compass C{degrees360}).
       @kwarg form: Format specifier for B{C{deg}} (C{str} or C{F_...}).
       @kwarg prec: Number of decimal digits (0..9 or C{None} for default).
       @kwarg sep: Separator between degrees, minutes, seconds, suffix (C{str}).
       @kwarg s_D_M_S: Optional keyword arguments to override any or cancel all
                       DMS symbol suffixes, see function L{pygeodesy.toDMS}.

       @return: Compass degrees per the specified B{C{form}} (C{str}).

       @see: Function L{pygeodesy.toDMS} and its B{Notes} for further details.
    '''
    return _toDMS(_umod_360(bearing), form, prec, sep, 1, NN, s_D_M_S)


def _clip(angle, limit, units):
    '''(INTERNAL) Helper for C{clipDegrees} and C{clipRadians}.
    '''
    c = min(limit, max(-limit, angle))
    if c != angle and rangerrors():
        t = _SPACE_(fstr(angle, prec=6, ints=True), _beyond_,
                    copysign0(limit, angle), units)
        raise RangeError(t, txt=None)
    return c


def clipDegrees(deg, limit):
    '''Clip a lat- or longitude to the given range.

       @arg deg: Unclipped lat- or longitude (C{scalar degrees}).
       @arg limit: Valid C{-/+B{limit}} range (C{degrees}).

       @return: Clipped value (C{degrees}).

       @raise RangeError: If B{C{deg}} outside the valid C{-/+B{limit}} range
                          and L{rangerrors<pygeodesy.rangerrors>} is C{True}.
    '''
    return _clip(deg, limit, _degrees_) if limit and limit > 0 else deg


def clipRadians(rad, limit):
    '''Clip a lat- or longitude to the given range.

       @arg rad: Unclipped lat- or longitude (C{radians}).
       @arg limit: Valid C{-/+B{limit}} range (C{radians}).

       @return: Clipped value (C{radians}).

       @raise RangeError: If B{C{rad}} outside the valid C{-/+B{limit}} range
                          and L{rangerrors<pygeodesy.rangerrors>} is C{True}.
    '''
    return _clip(rad, limit, _radians_) if limit and limit > 0 else rad


def compassDMS(bearing, form=F_D, prec=None, sep=S_SEP, **s_D_M_S):
    '''Convert bearing to a string suffixed with compass point.

       @arg bearing: Bearing from North (compass C{degrees360}).
       @kwarg form: Format specifier for B{C{deg}} (C{str} or C{F_...}).
       @kwarg prec: Number of decimal digits (0..9 or C{None} for default).
       @kwarg sep: Separator between degrees, minutes, seconds, suffix (C{str}).
       @kwarg s_D_M_S: Optional keyword arguments to override any or cancel all
                       DMS symbol suffixes, see function L{pygeodesy.toDMS}.

       @return: Compass degrees and point in the specified form (C{str}).

       @see: Function L{pygeodesy.toDMS} and its B{Notes} for further details.
    '''
    b = _umod_360(bearing)
    return _toDMS(b, form, prec, sep, 1, compassPoint(b), s_D_M_S)


def compassPoint(bearing, prec=3):
    '''Convert a C{bearing} from North to a compass point.

       @arg bearing: Bearing (compass C{degrees360}).
       @kwarg prec: Precision, number of compass point characters:
                    1 for cardinal or basic winds,
                    2 for intercardinal or ordinal or principal winds,
                    3 for secondary-intercardinal or half-winds or
                    4 for quarter-winds).

       @return: Compass point (1-, 2-, 3- or 4-letter C{str}).

       @raise ValueError: Invalid B{C{bearing}} or B{C{prec}}.

       @see: U{Dms.compassPoint
             <https://GitHub.com/ChrisVeness/geodesy/blob/master/dms.js>}
             and U{Compass rose<https://WikiPedia.org/wiki/Compass_rose>}.
    '''
    try:  # like .streprs.enstr2
        b = _umod_360(bearing)
        p = _MODS.units.Precision_(prec, low=1, high=4) \
             if prec != 3 else int(prec)
        m =  2 << p
        w = 32 // m  # if _isin(m, 4, 8, 16, 32)
        # not round(b), half-even rounding in Python 3+, but
        # round-away-from-zero as int(b + copysign0(_0_5, b))
        w *= int(b * m / _360_0 + _0_5) % m
        return _WINDS[w]
    except Exception as x:
        raise _xError(x, bearing=bearing, prec=prec)


def degDMS(deg, prec=6, s_D=S_DEG, s_M=S_MIN, s_S=S_SEC, neg=_MINUS_, pos=NN):
    '''Convert degrees to a string in degrees, minutes I{or} seconds.

       @arg deg: Value in degrees (C{scalar degrees}).
       @kwarg prec: Number of decimal digits (0..9 or C{None} for default).
                    Trailing zero decimals are stripped for C{B{prec}=1}
                    and above, but kept for negative B{C{prec}}.
       @kwarg s_D: D symbol for degrees (C{str}).
       @kwarg s_M: M symbol for minutes (C{str}) or C{""} aka L{pygeodesy.NN}.
       @kwarg s_S: S symbol for seconds (C{str}) or C{""} aka L{pygeodesy.NN}.
       @kwarg neg: Optional sign for negative (C{'-'}).
       @kwarg pos: Optional sign for positive (C{""}) aka L{pygeodesy.NN}.

       @return: I{Either} degrees, minutes I{or} seconds (C{str}).

       @see: Function L{pygeodesy.toDMS}.
    '''
    try:
        deg = float(deg)
    except Exception as x:
        raise _xError(x, deg=deg, prec=prec)

    d, s = fabs(deg), s_D
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
    return NN(n, t, s)  # NOT unicode in Python 2-


def latDMS(deg, form=_F_DMS, prec=None, sep=S_SEP, **s_D_M_S):
    '''Convert latitude to a string, optionally suffixed with N or S.

       @arg deg: Latitude to be formatted (C{scalar degrees}).
       @kwarg form: Format specifier for B{C{deg}} (C{str} or C{F_...}).
       @kwarg prec: Number of decimal digits (0..9 or C{None} for default).
       @kwarg sep: Separator between degrees, minutes, seconds, suffix (C{str}).
       @kwarg s_D_M_S: Optional keyword arguments to override any or cancel all
                       DMS symbol suffixes, see function L{pygeodesy.toDMS}.

       @return: Degrees in the specified form (C{str}).

       @see: Function L{pygeodesy.toDMS} and its B{Notes} for further details.
    '''
    P = _S_ if deg < 0 else _N_
    return _toDMS(deg, form, prec, sep, 2, P, s_D_M_S)


def latlonDMS(lls, **m_form_prec_sep_s_D_M_S):
    '''Convert one or more C{LatLon} instances to strings.

       @arg lls: Single (C{LatLon}) or an iterable (C{LatLon}s).
       @kwarg m_form_prec_sep_s_D_M_S: Optional keyword arguments C{B{m}eter},
                C{B{form}at}, C{B{prec}ision}, B{C{s_D}}, B{C{s_M}}, B{C{s_S}},
                B{C{s_DMS}} and I{DEPRECATED}  C{B{sep}=None}, see function
                L{pygeodesy.toDMS} and method C{LatLon.toStr} for more details.

       @return: A C{tuple} of C{str}s if B{C{lls}} is a list, sequence, tuple, etc.
                of C{LatLon}s or a single C{str} if B{C{lls}} is a single C{LatLon}.

       @see: Functions L{pygeodesy.latlonDMS_}, L{pygeodesy.latDMS}, L{pygeodesy.lonDMS}
             and L{pygeodesy.toDMS} and method C{LatLon.toStr}.

       @note: Keyword argument C{B{sep}=None} to join the resturned C{tuple} has
              been I{DEPRECATED}, use C{B{sep}.join(B{latlonDMS_}(...))} instead.
    '''
    sep, kwds = _latlonDMS_sep2(latlonDMS, **m_form_prec_sep_s_D_M_S)
    if isLatLon(lls):
        t = lls.toStr(**kwds)
    elif issequence(lls):
        t = tuple(ll.toStr(**kwds) for ll in lls)
        if sep:  # XXX DEPRECATED, to be removed
            t = sep.join(t)
    else:
        raise _TypeError(lls=lls, **m_form_prec_sep_s_D_M_S)
    return t


def latlonDMS_(*lls, **m_form_prec_sep_s_D_M_S):
    '''Convert one or more C{LatLon} instances to strings.

       @arg lls: The instances (C{LatLon}s), all positional arguments.
       @kwarg m_form_prec_sep_s_D_M_S: Optional keyword arguments C{B{m}eter},
                C{B{form}at}, C{B{prec}ision}, B{C{s_D}}, B{C{s_M}}, B{C{s_S}},
                B{C{s_DMS}} and I{DEPRECATED}  C{B{sep}=None}, see function
                L{pygeodesy.toDMS} and method C{LatLon.toStr} for more details.

       @return: A C{tuple} of C{str}s if 2 or more C{LatLon} instances
                or a single C{str} if only a single C{LatLon} instance
                is given in B{C{lls}}.

       @see: Functions L{pygeodesy.latlonDMS}, L{pygeodesy.latDMS} and
             L{pygeodesy.lonDMS} and L{pygeodesy.toDMS} and method
             C{LatLon.toStr}.

       @note: Keyword argument C{B{sep}=None} to join the resturned C{tuple} has
              been I{DEPRECATED}, use C{B{sep}.join(B{latlonDMS_}(...))} instead.
    '''
    sep, kwds = _latlonDMS_sep2(latlonDMS, **m_form_prec_sep_s_D_M_S)
    if not lls:
        raise _ValueError(lls=lls, **m_form_prec_sep_s_D_M_S)
    elif len(lls) < 2:
        lls, sep = lls[0], None
    t = latlonDMS(lls, **kwds)
    return sep.join(t) if sep else t


def _latlonDMS_sep2(where, sep=None, **kwds):
    '''DEPRECATED, use %r.join(%s(...)) instead.'''
    if sep:  # PYCHOK no cover
        i = _MODS.inters
        k =  Fmt.EQUAL(sep=repr(sep))
        k = _SPACE_(i._keyword_, i._arg_, k, i._of_)
        n =  typename(where)
        t = _latlonDMS_sep2.__doc__ % (sep, n)
        _MODS.props._throwarning(k, n, t)
    return sep, kwds


def lonDMS(deg, form=_F_DMS, prec=None, sep=S_SEP, **s_D_M_S):
    '''Convert longitude to a string, optionally suffixed with E or W.

       @arg deg: Longitude to be formatted (C{scalar degrees}).
       @kwarg form: Format specifier for B{C{deg}} (C{str} or C{F_...}).
       @kwarg prec: Number of decimal digits (0..9 or C{None} for default).
       @kwarg sep: Separator between degrees, minutes, seconds, suffix (C{str}).
       @kwarg s_D_M_S: Optional keyword arguments to override any or cancel all
                       DMS symbol suffixes.

       @return: Degrees in the specified form (C{str}).

       @see: Function L{pygeodesy.toDMS} and its B{Notes} for further details.
    '''
    P = _W_ if deg < 0 else _E_
    return _toDMS(deg, form, prec, sep, 3, P, s_D_M_S)


def normDMS(strDMS, norm=None, **s_D_M_S):
    '''Normalize all degrees, minutes and seconds (DMS) I{symbol} suffixes in
       a string to the default symbols L{S_DEG}, L{S_MIN}, L{S_SEC}.

       @arg strDMS: Original DMS string (C{str}).
       @kwarg norm: Optional replacement symbol (C{str}) or C{None} for the default
                    DMS symbol).  Use C{B{norm}=""} to remove all DMS symbols.
       @kwarg s_D_M_S: Optional, alternate DMS symbol suffixes C{B{s_D}=}L{S_DEG},
                       C{B{s_M}=}L{S_MIN}, C{B{s_S}=}L{S_SEC} and C{B{s_R}=}L{S_RAD}
                       for radians, each to be replaced by B{C{norm}}.

       @return: Normalized DMS (C{str}).
    '''
    def _s2S2(s_D=S_DEG, s_M=S_MIN, s_S=S_SEC, s_R=S_RAD):
        d = {s_D: S_DEG, s_M: S_MIN, s_S: S_SEC, s_R: S_RAD}
        for s, S in _xkwds(d, **_S_norm).items():
            if s:
                yield s, S

    # XXX strDMS isn't unicode in Python 2- and looping
    # thru strDMS will yield each byte, hence the loop
    # thru _s2S2 and replacing the DMS symbols in strDMS

    if norm is None:  # back to default DMS
        for s, S in _s2S2(**s_D_M_S):
            if s != S:
                strDMS = strDMS.replace(s, S)

    else:  # replace or remove all DMS
        n = norm or NN
        for s, _ in _s2S2(**s_D_M_S):
            if s != n:
                strDMS = strDMS.replace(s, n)
        if n:
            strDMS = strDMS.rstrip(n)  # XXX not .strip?

    return strDMS  # NOT unicode in Python 2-


def parseDDDMMSS(strDDDMMSS, suffix=_NSEW_, sep=S_SEP, clip=0, sexagecimal=False):
    '''Parse a lat- or longitude represention forms as [D]DDMMSS in degrees.

       @arg strDDDMMSS: Degrees in any of several forms (C{str}) and types (C{float},
                        C{int}, other).
       @kwarg suffix: Optional compass points (C{str}), valid in B{C{strDDDMMSS}}.
       @kwarg sep: Optional separator between "[D]DD", "MM", "SS", B{C{suffix}}
                   (L{S_SEP}) in B{C{strDDDMMSS}}.
       @kwarg clip: Optionally, limit value to range C{-/+B{clip}} (C{degrees}).
       @kwarg sexagecimal: If C{True}, convert C{"D.MMSS"} or C{float(D.MMSS)} to
                           C{base-60} "MM" and "SS" digits.  See C{form}s L{F_D60},
                           L{F_D60_} and L{F_D60__}.

       @return: Degrees (C{float}).

       @raise ParseError: Invalid B{C{strDDDMMSS}} or B{C{clip}} or the form of
                          B{C{strDDDMMSS}} is incompatible with or invalid for the
                          given B{C{suffix}} compass point(s).

       @raise RangeError: Value of B{C{strDDDMMSS}} outside the valid C{-/+B{clip}}
                          range and L{rangerrors<pygeodesy.rangerrors>} is C{True}.

       @note: Type C{str} values "[D]DD", "[D]DDMM", "[D]DDMMSS" and "[D]DD.MMSS"
              for B{C{strDDDMMSS}} are parsed properly only if I{either} unsigned
              and suffixed with a valid, compatible, C{cardinal} L{compassPoint}
              I{or} signed I{or} unsigned, unsuffixed and with keyword argument
              B{C{suffix}="NS"}, B{C{suffix}="EW"} or a compatible L{compassPoint}.

       @note: Unlike function L{parseDMS}, type C{float}, C{int} and other non-C{str}
              B{C{strDDDMMSS}} values are interpreted as C{form} [D]DDMMSS or
              [D]DD.MMSS.  For example, C{int(1230)} is returned as 12.5 and I{not
              1230.0} degrees.  However, C{int(345)} is considered C{form} "DDD"
              345 I{and not "DDMM" 0345}, unless B{C{suffix}} specifies the compass
              point.  Also, C{float(15.0523)} is returned as 15.0523 decimal
              degrees and I{not 15°5′23″ sexagecimal}.  To consider the latter, use
              C{float(15.0523)} or C{"15.0523"} and specify the keyword argument
              C{B{sexagecimal}=True}.

       @see: Functions L{pygeodesy.parseDMS}, L{pygeodesy.parseDMS2} and
             L{pygeodesy.parse3llh}.
    '''
    def _DDDMMSS(strDDDMMSS, suffix, sep, clip, sexagecimal):
        S = suffix.upper() or _N_  # bool('' in _NE_) is True
        if isstr(strDDDMMSS):
            t =  strDDDMMSS.replace(sep, NN) if sep else strDDDMMSS
            n, s, t, f, P, X = _DDDMMSS6(t, S)
            if X:
                pass
            elif not sexagecimal:  # try other forms
                return _DMS2deg(strDDDMMSS, S, sep, clip, {})

            if sexagecimal:  # move decimal dot from ...
                n += 4  # ... [D]DD.MMSSs to [D]DDMMSS.s
                if n < 6:
                    t = _SPACE_(_sexagecimal_, 'digits')
                    raise ParseError(t, n)
                z = n - len(f)  # zeros to append
                t = (f + (_0_ * z)) if z > 0 else _DOT_(f[:n], f[n:])
            f = _0_0  # fraction
        else:  # float or int to [D]DDMMSS[.fff]
            f = float(strDDDMMSS)
            n, s, t, f, P = _DDDMMSS5(f, S, sexagecimal)

        if n < 4:  # [D]DD[.ddd]
            if t:
                f += float(t)
            t = f,
        else:
            n -= 2
            f += float(t[n:])
            if n < 4:  # [D]DDMM[.mmm]
                t = float(t[:n]), f
            else:  # [D]DDMMSS[.sss]
                d = n - 2
                t = float(t[:d]), float(t[d:n]), f
        d = _dms2deg(s, P, *t)
        return clipDegrees(d, float(clip)) if clip else d

    return _parseX(_DDDMMSS, strDDDMMSS, suffix, sep, clip, sexagecimal,
                             strDDDMMSS=strDDDMMSS, suffix=suffix, sexagecimal=sexagecimal)


def _DDDMMSS5(f, S, sexagecimal):
    '''(INTERNAL) Partial C{parseDDDMMSS} of C{float}.
    '''
    if sexagecimal:
        f *= _SEXAGECIMUL
        m = 6
    else:
        m = 0
    s = P = _PLUS_  # anything except NN, _S_, _SW_, _W_
    if f < 0:
        f = -f
        s = _MINUS_
    f, i = modf(f)   # returns ...
    t = str(int(i))  # ... float(i)
    n = len(t)  # number of digits to ...
    if n < m:  # ... required min or ...
        t = (_0_ * (m - n)) + t
    # ... match the given compass point
    elif S in (_NS_ if isodd(n) else _EW_):
        t = _0_ + t
    #   P = S
    # elif n > 1:
    #   P = (_EW_ if isodd(n) else _NS_)[0]
    return len(t), s, t, f, P  # float(f) fraction


def _DDDMMSS6(t, S):
    '''(INTERNAL) Partial C{parseDDDMMSS} of C{str}.
    '''
    s, t, P = _split3(t, S)
    f = t.split(_DOT_)
    if len(f) > 2:
        raise ParseError('dots', len(f) - 1)
    n = len(f[0])
    f = NN.join(f)
    if 1 < n < 8 and f.isdigit() and (  # dddN/S/E/W or ddd or +/-ddd
                         (P in S and  s.isdigit()) or
                    (P.isdigit() and  S in _WINDS  # PYCHOK indent
                                 and (s in _PLUSMINUS_ or s.isdigit()))):
        # check [D]DDMMSS form and compass point
        X = _EW_ if isodd(n) else _NS_
        if not (P in X or (S in X and (P.isdigit() or P == _DOT_))):
            t =  typename(parseDDDMMSS)[(5 if isodd(n) else 6):]
            t = _SPACE_('form', t, 'applies', _DASH_.join(X))
            raise ParseError(t)
    else:
        X = NN
    return n, s, t, f, P, X  # str(f) fraction


def _DMS2deg(strDMS, suffix, sep, clip, s_D_M_S):
    '''(INTERNAL) Helper for C{parseDDDMMSS} and C{parseDMS}.
    '''
    try:
        d = float(strDMS)
    except (TypeError, ValueError):
        s, t, P = _split3(strDMS, suffix.upper())
        if sep:  # remove all DMS symbols
            t = t.replace(sep, _SPACE_)
            t = normDMS(t, norm=NN, **s_D_M_S)
        else:  # replace all DMS symbols
            t = normDMS(t, norm=_SPACE_, **s_D_M_S)
        t =  t.strip().split()
        d = _dms2deg(s, P, *map(float, t))
    return clipDegrees(d, float(clip)) if clip else d


def _dms2deg(n, P, d, m=_0_0, s=_0_0):
    '''(INTERNAL) Helper for C{parseDDDMMSS}, C{parseRad} and C{_DMS2deg}.
    '''
    m += s / _60_0
    d += m / _60_0
    if n == _MINUS_ or (P and P in _SW_):
        d = _neg(d)
    return d


def parseDMS(strDMS, suffix=_NSEW_, sep=S_SEP, clip=0, **s_D_M_S):  # MCCABE 14
    '''Parse a lat- or longitude representation in C{degrees}.

       This is very flexible on formats, allowing signed decimal degrees, degrees
       and minutes or degrees minutes and seconds optionally suffixed by a cardinal
       compass point.

       A variety of symbols, separators and suffixes are accepted, for example
       "3°37′09″W".  Minutes and seconds may be omitted.

       @arg strDMS: Degrees in any of several forms (C{str}) and types (C{float},
                    C{int}, other).
       @kwarg suffix: Optional, valid compass points (C{str}, C{tuple}).
       @kwarg sep: Optional separator between C{deg°}, C{min′}, C{sec″},
                   B{C{suffix}} (C{''}).
       @kwarg clip: Optionally, limit value to range C{-/+B{clip}} (C{degrees}).
       @kwarg s_D_M_S: Optional, alternate DMS symbol suffixes for degrees
                       C{B{s_D}=}L{S_DEG}, minutes C{B{s_M}=}L{S_MIN} and
                       seconds C{B{s_S}=}L{S_SEC}, see function L{pygeodesy.toDMS}.

       @return: Degrees (C{float}).

       @raise ParseError: Invalid B{C{strDMS}} or B{C{clip}}.

       @raise RangeError: Value of B{C{strDMS}} outside the valid C{-/+B{clip}} range
                          and L{rangerrors<pygeodesy.rangerrors>} is C{True}.

       @note: Unlike function L{parseDDDMMSS}, type C{float}, C{int} and other non-C{str}
              B{C{strDMS}} values are considered decimal (and not sexagecimal) degrees.
              For example, C{int(1230)} is returned as 1230.0 I{and not as 12.5} degrees
              and C{float(345)} as 345.0 I{and not as 3.75} degrees!

       @see: Functions L{pygeodesy.parseDDDMMSS}, L{pygeodesy.parseDMS2},
             L{pygeodesy.parse3llh} and L{pygeodesy.toDMS}.
    '''
    return _parseX(_DMS2deg, strDMS, suffix, sep, clip, s_D_M_S, strDMS=strDMS, suffix=suffix)


def parseDMS2(strLat, strLon, sep=S_SEP, clipLat=90, clipLon=180, wrap=False, **s_D_M_S):
    '''Parse a lat- and a longitude representions C{"lat, lon"} in C{degrees}.

       @arg strLat: Latitude in any of several forms (C{str} or C{degrees}).
       @arg strLon: Longitude in any of several forms (C{str} or C{degrees}).
       @kwarg sep: Optional separator between deg°, min′, sec″, suffix (C{''}).
       @kwarg clipLat: Limit latitude to range C{-/+B{clipLat}} (C{degrees}).
       @kwarg clipLon: Limit longitude to range C{-/+B{clipLon}} (C{degrees}).
       @kwarg wrap: If C{True}, wrap or I{normalize} the lat- and longitude,
                    overriding B{C{clipLat}} and B{C{clipLon}} (C{bool}).
       @kwarg s_D_M_S: Optional, alternate DMS symbol suffixes for degrees
                       C{B{s_D}=}L{S_DEG}, minutes C{B{s_M}=}L{S_MIN} and
                       seconds C{B{s_S}=}L{S_SEC}, see function L{pygeodesy.toDMS}.

       @return: A L{LatLon2Tuple}C{(lat, lon)} in C{degrees}.

       @raise ParseError: Invalid B{C{strLat}} or B{C{strLon}}.

       @raise RangeError: Value of B{C{strLat}} or B{C{strLon}} outside the
                          valid C{-/+B{clipLat}} or C{-/+B{clipLon}} range
                          and L{rangerrors<pygeodesy.rangerrors>} is C{True}.

       @note: See the B{Notes} at function L{parseDMS}.

       @see: Functions L{pygeodesy.parseDDDMMSS}, L{pygeodesy.parseDMS},
             L{pygeodesy.parse3llh} and L{pygeodesy.toDMS}.
    '''
    return _2Tuple(strLat, strLon, clipLat, clipLon, wrap, sep=sep, **s_D_M_S)


def _2Tuple(strLat, strLon, clipLat, clipLon, wrap, **kwds):
    '''(INTERNAL) Helper for C{parseDMS2} and C{parse3llh}.
    '''
    if wrap:
        _W = _MODS.utily._Wrap
        ll = _W.latlon(parseDMS(strLat, suffix=_NS_, **kwds),
                       parseDMS(strLon, suffix=_EW_, **kwds))
    else:
        # if wrap is None:
        #     clipLat = clipLon = 0
        ll = (parseDMS(strLat, suffix=_NS_, clip=clipLat, **kwds),
              parseDMS(strLon, suffix=_EW_, clip=clipLon, **kwds))
    return _MODS.namedTuples.LatLon2Tuple(*ll)


def parse3llh(strllh, height=0, sep=_COMMA_, clipLat=90, clipLon=180, wrap=False, **s_D_M_S):
    '''Parse a string C{"lat, lon [, h]"} representing lat-, longitude in
       C{degrees} and optional height in C{meter}.

       The lat- and longitude value must be separated by a separator
       character.  If height is present it must follow, separated by
       another separator.

       The lat- and longitude values may be swapped, provided at least
       one ends with the proper compass point.

       @arg strllh: Latitude, longitude[, height] (C{str}, ...).
       @kwarg height: Optional, default height (C{meter}) or C{None}.
       @kwarg sep: Optional separator between C{"lat lon [h] suffix"} (C{str}).
       @kwarg clipLat: Limit latitude to range C{-/+B{clipLat}} (C{degrees}).
       @kwarg clipLon: Limit longitude to range C{-/+B{clipLon}} (C{degrees}).
       @kwarg wrap: If C{True}, wrap or I{normalize} the lat- and longitude,
                    overriding B{C{clipLat}} and B{C{clipLon}} (C{bool}).
       @kwarg s_D_M_S: Optional, alternate DMS symbol suffixes for degrees
                       C{B{s_D}=}L{S_DEG}, minutes C{B{s_M}=}L{S_MIN} and
                       seconds C{B{s_S}=}L{S_SEC}, see function L{pygeodesy.toDMS}.

       @return: A L{LatLon3Tuple}C{(lat, lon, height)} in C{degrees},
                C{degrees} and C{float}.

       @raise RangeError: Lat- or longitude value of B{C{strllh}} outside the
                          valid C{-/+B{clipLat}} or C{-/+B{clipLon}} range
                          and L{rangerrors<pygeodesy.rangerrors>} is C{True}.

       @raise ValueError: Invalid B{C{strllh}} or B{C{height}}.

       @note: See the B{Notes} at function L{parseDMS}.

       @see: Functions L{pygeodesy.parseDDDMMSS}, L{pygeodesy.parseDMS},
             L{pygeodesy.parseDMS2} and L{pygeodesy.toDMS}.
    '''
    def _3llh(strllh, height, sep, wrap):
        ll = strllh.strip().split(sep)
        if len(ll) > 2:  # XXX interpret height unit
            h = float(ll.pop(2).rstrip(_LETTERS + _SPACE_))
        else:
            h = height  # None from wgrs.Georef.__new__
        if len(ll) != 2:
            raise ValueError

        a, b = [_.strip() for _ in ll]  # PYCHOK false
        if a[-1:] in _EW_ or b[-1:] in _NS_:
            a, b = b, a
        t = _2Tuple(a, b, clipLat, clipLon, wrap, **s_D_M_S)
        return t.to3Tuple(h)

    return _parseX(_3llh, strllh, height, sep, wrap, strllh=strllh)


def parseRad(strRad, suffix=_NSEW_, clip=0):
    '''Parse a string representing angle in C{radians}.

       @arg strRad: Degrees in any of several forms (C{str} or C{radians}).
       @kwarg suffix: Optional, valid compass points (C{str}, C{tuple}).
       @kwarg clip: Optionally, limit value to range C{-/+B{clip}} (C{radians}).

       @return: Radians (C{float}).

       @raise ParseError: Invalid B{C{strRad}} or B{C{clip}}.

       @raise RangeError: Value of B{C{strRad}} outside the valid C{-/+B{clip}}
                          range and L{rangerrors<pygeodesy.rangerrors>} is C{True}.
    '''
    def _Rad(strRad, suffix, clip):
        try:
            r = float(strRad)
        except (TypeError, ValueError):
            s, t, P = _split3(strRad, suffix.upper())
            r = _dms2deg(s, P, float(t))
        return clipRadians(r, float(clip)) if clip else r

    return _parseX(_Rad, strRad, suffix, clip, strRad=strRad, suffix=suffix)


def precision(form, prec=None):
    '''Set the default precison for a given F_ form.

       @arg form: L{F_D}, L{F_DM}, L{F_DMS}, L{F_DEG}, L{F_MIN}, L{F_SEC},
                  L{F_D60}, L{F__E}, L{F__F}, L{F__G} or L{F_RAD} (C{str}).
       @kwarg prec: Number of decimal digits (0..9 or C{None} for default).
                    Trailing zero decimals are stripped for C{B{prec}=1}
                    and above, but kept for negative B{C{prec}}.

       @return: Previous precision for the B{C{form}} (C{int}).

       @raise ValueError: Invalid B{C{form}} or B{C{prec}} or B{C{prec}}
                          outside range C{-9..+9}.
    '''
    try:
        p = _F_prec[form]
    except KeyError:
        raise _ValueError(form=form)
    if prec is not None:  # set as default
        P_ = _MODS.units.Precision_
        _F_prec[form] = P_(prec, low=-9, high=9)
    return p


def toDMS(deg, form=_F_DMS, prec=2, sep=S_SEP, ddd=2, neg=_MINUS_, pos=_PLUS_, **s_D_M_S):
    '''Convert I{signed} C{degrees} to string, without suffix.

       @arg deg: Degrees to be formatted (C{scalar degrees}).
       @kwarg form: Format specifier for B{C{deg}} (C{str} or L{F_D}, L{F_DM}, L{F_DMS},
                    L{F_DEG}, L{F_MIN}, L{F_SEC}, L{F_D60}, L{F__E}, L{F__F}, L{F__G},
                    L{F_RAD}, L{F_D_}, L{F_DM_}, L{F_DMS_}, L{F_DEG_}, L{F_MIN_}, L{F_SEC_},
                    L{F_D60_}, L{F__E_}, L{F__F_}, L{F__G_}, L{F_RAD_}, L{F_D__}, L{F_DM__},
                    L{F_DMS__}, L{F_DEG__}, L{F_MIN__}, L{F_SEC__}, L{F_D60__}, L{F__E__},
                    L{F__F__}, L{F__G__} or L{F_RAD__}).
       @kwarg prec: Number of decimal digits (0..9 or C{None} for default).  Trailing zero
                    decimals are stripped for C{B{prec}=1} and above, but kept for negative
                    B{C{prec}}.
       @kwarg sep: Separator between degrees, minutes, seconds, suffix (C{str}).
       @kwarg ddd: Number of (integer) digits for B{C{deg}°} (2 or 3).
       @kwarg neg: Prefix for negative B{C{deg}} (C{'-'}).
       @kwarg pos: Prefix for positive B{C{deg}} and signed B{C{form}} (C{'+'}).
       @kwarg s_D_M_S: Optional keyword arguments C{B{s_D}=}L{S_DEG}, C{B{s_M}=}L{S_MIN}
                       C{B{s_S}=}L{S_SEC}, C{B{s_R}=}L{S_RAD} and C{B{s_DMS}=True} to
                       override any or cancel all DMS symbol suffixes.  See the B{Notes}
                       below.

       @return: Degrees in the specified form (C{str}).

       @note: The degrees, minutes and seconds (DMS) symbol can be overridden in this and
              other C{*DMS} functions by using optional keyword argments C{B{s_D}="d"},
              C{B{s_M}="'"} respectively C{B{s_S}='"'}.  Using B{C{s_DMS}=None} cancels
              all C{DMS} symbols to C{""} aka L{pygeodesy.NN}.

       @note: Sexagecimal format B{C{F_D60}} supports overridable pseudo-DMS symbols
              positioned at C{"[D]DD<B{s_D}>MM<B{s_M}>SS<B{s_S}>"} with defaults
              C{B{s_D}="."}, C{B{s_M}=B{sep}} and C{B{s_S}=}L{pygeodesy.NN}.

       @note: Formats B{C{F__E}}, B{C{F__F}} and B{C{F__G}} is extended with a C{D}-only
              symbol suffix if defined with keyword argument C{B{s_D}=}L{pygeodesy.NN}.
              Likewise for B{C{F_RAD}} forms with keyword argument C{B{s_R}=}L{S_RAD}.

       @see: Function L{pygeodesy.degDMS}
    '''
    s =  form[:1]
    f =  form[1:] if s in _PLUSMINUS_ else form
    t = _toDMS(deg, f, prec, sep, ddd, NN, s_D_M_S)  # unsigned and -suffixed
    if deg < 0 and neg:
        t = neg + t
    elif deg > 0 and s == _PLUS_ and pos:
        t = pos + t
    return t

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
