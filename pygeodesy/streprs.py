
# -*- coding: utf-8 -*-

u'''Floating point and other formatting utilities.

'''

from pygeodesy.basics import isint, isscalar
from pygeodesy.errors import _AttributeError, _IsnotError, _TypeError, \
                             _ValueError, _xkwds_pop
from pygeodesy.interns import NN, MISSING, _BAR_, _COMMASPACE_, _DOT_, _E_, \
                             _EQUAL_, _H_, _N_, _name_, _not_, _PERCENT_, \
                             _OKd_, _scalar_, _SPACE_, _STAR_, _UNDER_, _0_, \
                             _0_0, _0_001, _0_01, _0_1, _1_0
from pygeodesy.interns import _convergence_, _distant_, _e_, _EPS0_, \
                              _EQUALSPACED_, _exceeds_, _f_, _F_, _g_  # PYCHOK used!
from pygeodesy.lazily import _ALL_LAZY

__all__ = _ALL_LAZY.streprs
__version__ = '21.11.15'

_E_4_E0 = (1e-4, _0_001, _0_01, _0_1, _1_0)


def _pct(fmt):
    '''(INTERNAL) Prefix C{%} if needed.
    '''
    return fmt if _PERCENT_ in fmt else NN(_PERCENT_, fmt)


def _0wd(*w_i):  # in .osgr, .wgrs
    '''(INTERNAL) Int formatter'.
    '''
    return '%0*d' % w_i


def _0wpF(*w_p_f):  # in .dms, .osgr
    '''(INTERNAL) Float deg, min, sec formatter'.
    '''
    return '%0*.*f' % w_p_f  # XXX was F


class _Fmt(str):  # in .streprs
    '''(INTERNAL) Callable formatting.
    '''
    name = NN

    def __call__(self, *name_value_, **name_value):
        '''Format a C{name=value} pair or C{name, value} pair
           or just a single C{value}.
        '''
        for n, v in name_value.items():
            break
        else:
            if len(name_value_) > 1:
                n, v = name_value_[:2]
            elif name_value_:
                n, v = NN, name_value_[0]
            else:
                n, v = NN, MISSING
        t = str.__mod__(self, v)
        return NN(n, t) if n else t

#   def __mod__(self, arg, **unused):
#       '''Regular C{%} operator.
#       '''
#       return str.__mod__(self, arg)


class Fstr(str):
    '''(INTERNAL) C{float} format.
    '''
    name = NN

    def __call__(self, flt, prec=None, ints=False):
        '''Format the B{C{flt}} like function L{fstr}.
        '''
        # see also function C{fstr} if isscalar case below
        t = str.__mod__(_pct(self), flt) if prec is None else next(
           _streprs(prec, (flt,), self, ints, True, None))
        return t

    def __mod__(self, arg, **unused):
        '''Regular C{%} operator.

           @arg arg: A C{scalar} value to be formatted (either
                     the C{scalar}, or a 1-tuple C{(scalar,)},
                     or 2-tuple C{(prec, scalar)}.

           @raise TypeError: Non-scalar B{C{arg}} value.

           @raise ValueError: Invalid B{C{arg}}.
        '''
        def _error(arg):
            n = _DOT_(Fstr.__name__, self.name or self)
            return _SPACE_(n, _PERCENT_, repr(arg))

        prec = 6  # default std %f and %F
        if isinstance(arg, (tuple, list)):
            n = len(arg)
            if n == 1:
                arg = arg[0]
            elif n == 2:
                prec, arg = arg
            else:
                raise _ValueError(_error(arg))

        if not isscalar(arg):
            raise _TypeError(_error(arg))
        return self(arg, prec=prec)


class Fmt(object):
    '''Formatting options.
    '''
    ANGLE       = _Fmt('<%s>')
    COLON       = _Fmt(':%s')
#   COLONSPACE  = _Fmt(': %s')  # == _COLONSPACE_(n, v)
#   COMMASPACE  = _Fmt(', %s')  # == _COMMASPACE_(n, v)
    convergence = _Fmt(_convergence_('(%g)'))
    CURLY       = _Fmt('{%s}')  # BRACES
    distant     = _Fmt(_distant_('(%.3g)'))
    DOT         = _Fmt('.%s')   # == NN(_DOT_, n)
    e           = Fstr(_e_)
    E           = Fstr(_E_)
    EPS0        = _Fmt(_SPACE_(_EPS0_, '(%g)'))
    EQUAL       = _Fmt(_EQUAL_(NN, '%s'))
    EQUALSPACED = _Fmt(_EQUALSPACED_(NN, '%s'))
    exceeds_eps = _Fmt(_exceeds_('eps', '(%g)'))
    f           = Fstr(_f_)
    F           = Fstr(_F_)
    g           = Fstr(_g_)
    G           = Fstr('G')
    h           = Fstr('%+.*f')  # height, .streprs.hstr
    limit       = _Fmt(' %s limit')  # .units
    LOPEN       = _Fmt('(%s]')  # left-open range (L, R]
    PAREN       = _Fmt('(%s)')
    PARENSPACED = _Fmt(' (%s)')
    QUOTE2      = _Fmt('"%s"')
    ROPEN       = _Fmt('[%s)')  # right-open range [L, R)
#   SPACE       = _Fmt(' %s')   # == _SPACE_(n, v)
    SQUARE      = _Fmt('[%s]')  # BRACKETS
    TAG         =  ANGLE
    TAGEND      = _Fmt('</%s>')
    zone        = _Fmt('%02d')  # .epsg, .mgrs, .utmupsBase

    def __init__(self):
        for n, a in self.__class__.__dict__.items():
            if isinstance(a, (Fstr, _Fmt)):
                setattr(a, _name_, n)

Fmt          = Fmt()  # PYCHOK singleton
Fmt.__name__ = Fmt.__class__.__name__

_DOTSTAR_ = Fmt.DOT(_STAR_)
# formats %G and %g drop all trailing zeros and the
# decimal point, making the float appear as an int
_Gg     = (Fmt.G, Fmt.g)
_FfEeGg = (Fmt.F, Fmt.f, Fmt.E, Fmt.e) + _Gg  # float formats
_Fspec_ = NN('[%[<flags>][<width>]', _DOTSTAR_, ']', _BAR_.join(_FfEeGg))  # in testStreprs


def _streprs(prec, objs, fmt, ints, force, strepr):
    '''(INTERNAL) Helper for C{fstr}, C{pairs}, C{reprs} and C{strs}
    '''
    # <https://docs.Python.org/3/library/stdtypes.html#printf-style-string-formatting>
    if fmt in _FfEeGg:
        fGg = fmt in _Gg
        fmt = NN(_PERCENT_, _DOT_, abs(prec), fmt)

    elif fmt.startswith(_PERCENT_):
        fGg = False
        try:  # to make sure fmt is valid
            f = fmt.replace(_DOTSTAR_, Fmt.DOT(abs(prec)))
            _ = f % (_0_0,)
        except (TypeError, ValueError):
            raise _ValueError(fmt=fmt, txt=_not_(repr(_DOTSTAR_)))
        fmt = f

    else:
        raise _ValueError(fmt=fmt, txt=_not_(repr(_Fspec_)))

    for o in objs:
        if force or isinstance(o, float):
            t = fmt % (float(o),)
            if ints and (isint(o, both=True) or  # for ...
                         # corner case testLcc lon0=-96.0
                         t.rstrip(_0_).endswith(_DOT_)):
                t = t.split(_DOT_)[0]
            elif prec > 1:
                t = fstrzs(t, ap1z=fGg)
        elif strepr:
            t = strepr(o)
        else:
            raise _IsnotError(_scalar_, floats=o)
        yield t


def anstr(name, OKd=_OKd_, sub=_UNDER_):
    '''Make a valid name of alphanumeric and OKd characters.

       @arg name: The original name (C{str}).
       @kwarg OKd: Other acceptable characters (C{str}).
       @kwarg sub: Substitute for invalid charactes (C{str}).

       @return: The modified name (C{str}).

       @note: Leading and trailing whitespace characters are removed,
              intermediate whitespace characters are coalesced and
              substituted.
    '''
    s = n = str(name).strip()
    for c in n:
        if not (c.isalnum() or c in OKd or c in sub):
            s = s.replace(c, _SPACE_)
    return sub.join(s.strip().split())


def attrs(inst, *names, **pairs_kwds):  # prec=6, fmt=Fmt.F, ints=False, Nones=True, sep=_EQUAL_
    '''Get instance attributes as I{name=value} strings, with C{float}s
       formatted by function L{fstr}.

       @arg inst: The instance (any C{type}).
       @arg names: The attribute names (C{str}s).
       @kwarg pairs_kwds: Keyword argument for function L{pairs}, except
                          C{B{Nones}=True} to in- or exclude missing or
                          or C{None}-valued attributes.

       @return: A C{tuple(sep.join(t) for t in zip(names, reprs(values, ...)))}
                of C{str}s.
    '''
    Nones = _xkwds_pop(pairs_kwds, Nones=True)

    def items():
        for n in names:
            v = getattr(inst, n, None)
            if Nones or v is not None:
                yield n, v

    return pairs(items(), **pairs_kwds)


def enstr2(easting, northing, prec, *extras):
    '''Return easting, northing string representations.

       @arg easting: Easting from false easting (C{meter}).
       @arg northing: Northing from from false northing (C{meter}).
       @arg prec: Precision in number of digits (C{int}, [1..5]).
       @arg extras: Optional leading items (C{str}s).

       @return: B{C{extras}} + 2-Tuple C{(eastingStr, northingStr)}.

       @raise ValueError: Invalid B{C{prec}}.
    '''
    w = int(prec) // 2
    if not 0 < w < 6:
        raise _ValueError(prec=prec)
    p = _E_4_E0[w - 1]  # 10**(5 - w)
    return extras + (_0wd(w, int(easting  * p)),
                     _0wd(w, int(northing * p)))


def fstr(floats, prec=6, fmt=Fmt.F, ints=False, sep=_COMMASPACE_, strepr=None):
    '''Convert one or more floats to string, optionally stripped of trailing zero decimals.

       @arg floats: Single or a list, sequence, tuple, etc. (C{scalar}s).
       @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                    Trailing zero decimals are stripped if B{C{prec}} is
                    positive, but kept for negative B{C{prec}} values.  In
                    addition, trailing decimal zeros are stripped for U{alternate,
                    form '#'<https://docs.Python.org/3/library/stdtypes.html
                    #printf-style-string-formatting>}.
       @kwarg fmt: Optional, C{float} format (C{str}).
       @kwarg ints: Optionally, remove the decimal dot for C{int} values (C{bool}).
       @kwarg sep: Separator joining the B{C{floats}} (C{str}).
       @kwarg strepr: Optional callable to format non-C{floats} (typically
                      C{repr}, C{str}) or C{None} to raise a TypeError.

       @return: The C{sep.join(strs(floats, ...)} joined (C{str}) or single
                C{strs((floats,), ...)} (C{str}) if B{C{floats}} is C{scalar}.
    '''
    if isscalar(floats):  # see Fstr.__call__ above
        return next(_streprs(prec, (floats,), fmt, ints, True, strepr))
    else:
        return sep.join(_streprs(prec, floats, fmt, ints, True, strepr))


def _fstrENH2(inst, prec, m):  # in .css, .lcc, .utmupsBase
    # (INTERNAL) For C{Css.} and C{Lcc.} C{toRepr} and C{toStr} and C{UtmUpsBase._toStr}.
    t = inst.easting, inst.northing
    t = tuple(_streprs(prec, t, Fmt.F, False, True, None))
    T = _E_, _N_
    if m is not None and abs(inst.height):  # abs(self.height) > EPS
        t +=  hstr(inst.height, prec=-2, m=m),
        T += _H_,
    return t, T


def _fstrLL0(inst, prec, toRepr):  # in .azimuthal, .css
    # (INTERNAL) For C{_AzimuthalBase.} and C{CassiniSoldner.} C{toStr} and C{toRepr}.
    t = tuple(_streprs(prec, inst.latlon0, Fmt.F, False, True, None))
    if toRepr:
        n = inst.name
        if n:
            t += Fmt.EQUAL(_name_, repr(n)),
        t = Fmt.PAREN(inst.classname, _COMMASPACE_.join(t))
    return t


def fstrzs(efstr, ap1z=False):
    '''Strip trailing zero decimals from a C{float} string.

       @arg efstr: Float with or without exponent (C{str}).
       @kwarg ap1z: Append the decimal point and one zero decimal
                    if the B{C{efstr}} is all digits (C{bool}).

       @return: Float (C{str}).
    '''
    s = efstr.find(_DOT_)
    if s >= 0:
        e = efstr.rfind(Fmt.e)
        if e < 0:
            e = efstr.rfind(Fmt.E)
            if e < 0:
                e = len(efstr)
        s += 2  # keep 1st _DOT_ + _0_
        if s < e and efstr[e-1] == _0_:
            efstr = NN(efstr[:s], efstr[s:e].rstrip(_0_), efstr[e:])

    elif ap1z:
        # %.G and %.g formats may drop the decimal
        # point and all trailing zeros, ...
        if efstr.isdigit():
            efstr += _DOT_ + _0_  # ... append or ...
        else:  # ... insert one dot and zero
            e = efstr.rfind(Fmt.e)
            if e < 0:
                e = efstr.rfind(Fmt.E)
            if e > 0:
                efstr = NN(efstr[:e], _DOT_, _0_, efstr[e:])

    return efstr


def hstr(height, prec=2, fmt=Fmt.h, ints=False, m=NN):
    '''Return a string for the height value.

       @arg height: Height value (C{float}).
       @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                    Trailing zero decimals are stripped if B{C{prec}} is
                    positive, but kept for negative B{C{prec}} values.
       @kwarg fmt: Optional, C{float} format (C{str}).
       @kwarg ints: Optionally, remove the decimal dot for C{int} values (C{bool}).
       @kwarg m: Optional unit of the height (C{str}).
    '''
    h = next(_streprs(prec, (height,), fmt, ints, True, None))
    return NN(h, str(m)) if m else h


def instr(inst, *args, **kwds):
    '''Return the string representation of an instantiation.

       @arg inst: The instance (any C{type}).
       @arg args: Optional positional arguments.
       @kwarg kwds: Optional keyword arguments.

       @return: Representation (C{str}).
    '''
    from pygeodesy.named import classname
    return unstr(classname(inst), *args, **kwds)


def pairs(items, prec=6, fmt=Fmt.F, ints=False, sep=_EQUAL_):
    '''Convert items to I{name=value} strings, with C{float}s handled like L{fstr}.

       @arg items: Name-value pairs (C{dict} or 2-{tuple}s of any C{type}s).
       @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                    Trailing zero decimals are stripped if B{C{prec}} is
                    positive, but kept for negative B{C{prec}} values.
       @kwarg fmt: Optional, C{float} format (C{str}).
       @kwarg ints: Optionally, remove the decimal dot for C{int} values (C{bool}).
       @kwarg sep: Separator joining I{names} and I{values} (C{str}).

       @return: A C{tuple(sep.join(t) for t in zip(names, reprs(values,...)))}
                of C{str}s.
    '''
    try:
        if isinstance(items, dict):
            items = sorted(items.items())
        elif not isinstance(items, (list, tuple)):
            items = tuple(items)
        # can't unzip empty items tuple, list, etc.
        n, v = zip(*items) if items else ((), ())
    except (TypeError, ValueError):
        raise _IsnotError(dict.__name__, '2-tuples', items=items)
    v = _streprs(prec, v, fmt, ints, False, repr)
    return tuple(sep.join(t) for t in zip(map(str, n), v))


def reprs(objs, prec=6, fmt=Fmt.F, ints=False):
    '''Convert objects to C{repr} strings, with C{float}s handled like L{fstr}.

       @arg objs: List, sequence, tuple, etc. (any C{type}s).
       @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                    Trailing zero decimals are stripped if B{C{prec}} is
                    positive, but kept for negative B{C{prec}} values.
       @kwarg fmt: Optional, C{float} format (C{str}).
       @kwarg ints: Optionally, remove the decimal dot for C{int} values (C{bool}).

       @return: A C{tuple(map(fstr|repr, objs))} of C{str}s.
    '''
    return tuple(_streprs(prec, objs, fmt, ints, False, repr)) if objs else ()


def strs(objs, prec=6, fmt=Fmt.F, ints=False):
    '''Convert objects to C{str} strings, with C{float}s handled like L{fstr}.

       @arg objs: List, sequence, tuple, etc. (any C{type}s).
       @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                    Trailing zero decimals are stripped if B{C{prec}} is
                    positive, but kept for negative B{C{prec}} values.
       @kwarg fmt: Optional, C{float} format (C{str}).
       @kwarg ints: Optionally, remove the decimal dot for C{int} values (C{bool}).

       @return: A C{tuple(map(fstr|str, objs))} of C{str}s.
    '''
    return tuple(_streprs(prec, objs, fmt, ints, False, str)) if objs else ()


def unstr(named, *args, **kwds):
    '''Return the string representation of an invokation.

       @arg named: Function, method or class name (C{str}).
       @arg args: Optional positional arguments.
       @kwarg kwds: Optional keyword arguments.

       @return: Representation (C{str}).
    '''
    t = reprs(args, fmt=Fmt.g) + pairs(sorted(kwds.items()))
    return Fmt.PAREN(named, _COMMASPACE_.join(t))


def _boolkwds(inst, **name_value_pairs):  # in .frechet, .hausdorff, .heights
    '''(INTERNAL) Set applicable C{bool} properties/attributes.
    '''
    for n, v in name_value_pairs.items():
        b = getattr(inst, n, None)
        if b is None:  # invalid bool attr
            t = _SPACE_(_EQUAL_(n, repr(v)), 'for', inst.__class__.__name__)  # XXX .classname
            raise _AttributeError(t, txt=_not_('applicable'))
        if v in (False, True) and v != b:
            setattr(inst, NN(_UNDER_, n), v)


def _xattrs(insto, other, *attrs):
    '''(INTERNAL) Copy attribute values from B{C{other}} to B{C{insto}}.

       @arg insto: Object to copy attribute values to (any C{type}).
       @arg other: Object to copy attribute values from (any C{type}).
       @arg attrs: One or more attribute names (C{str}s).

       @return: Object B{C{insto}}, updated.

       @raise AttributeError: An B{C{attrs}} doesn't exist
                              or is not settable.
    '''
    def _getattr(o, a):
        if hasattr(o, a):
            return getattr(o, a)
        try:
            n = o._DOT_(a)
        except AttributeError:
            n = Fmt.DOT(a)
        raise _AttributeError(o, name=n)

    for a in attrs:
        s = _getattr(other, a)
        g = _getattr(insto, a)
        if (g is None and s is not None) or g != s:
            setattr(insto, a, s)  # not settable?
    return insto


def _xzipairs(lefts, rights, sep=_COMMASPACE_, fmt=NN, pair_fmt=Fmt.COLON):
    '''(INTERNAL) Zip C{lefts} and C{rights} into a C{str}.
    '''
    t = sep.join(pair_fmt(*t) for t in zip(lefts, rights))
    return (fmt % (t,)) if fmt else t

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
