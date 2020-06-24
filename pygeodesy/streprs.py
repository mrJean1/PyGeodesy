
# -*- coding: utf-8 -*-

u'''Floating point and other formatting utilities.

'''

from pygeodesy.basics import isint, isscalar
from pygeodesy.errors import _IsnotError, _ValueError, _xkwds_pop
from pygeodesy.lazily import _ALL_LAZY

# all public contants, classes and functions
__all__ = _ALL_LAZY.streprs
__version__ = '20.06.24'

# formats %G and %.g drop all trailing zeros and the
# decimal point making the float appear as an int
_Gg     = ('G', 'g')
_EeFfGg = ('F', 'f', 'E', 'e') + _Gg  # float formats
_Fmt    =  'F'


def _streprs(prec, objs, fmt, ints, force, strepr):
    '''(INTERNAL) Helper for C{fstr}, C{pairs}, C{reprs} and C{strs}
    '''
    if fmt in _EeFfGg:
        fGg = fmt in _Gg
        fmt = '%.' + str(abs(prec)) + fmt
    elif fmt.startswith('%'):
        fGg = False
        fmt = fmt.replace('*', str(abs(prec)))
    else:
        t = '[%s]%s' % ('%.*', '|'.join(_EeFfGg))
        raise _ValueError(fmt=fmt, txt='not %r' % (t,))

    for o in objs:
        if force or isinstance(o, float):
            t = fmt % (float(o),)
            if ints and (isint(o, both=True) or  # for ...
                         # corner case testLcc lon0=-96.0
                         t.rstrip('0').endswith('.')):
                t = t.split('.')[0]
            elif prec > 1:
                t = fstrzs(t, ap1z=fGg)
        elif strepr:
            t = strepr(o)
        else:
            raise _IsnotError('scalar', floats=o)
        yield t


def anstr(name, OKd='._-', sub='_'):
    '''Make a valid name of alphanumeric and OKd characters.

       @arg name: The original name (C{str}).
       @kwarg OKd: Other acceptable characters (C{str}).
       @kwarg sub: Substitute for invalid charactes (C{str}).

       @return: The modified name (C{str}).

       @note: Leading and trailing whitespace characters are removed
              and intermediate whitespace characters are coalesced
              and substituted.
    '''
    s = n = str(name).strip()
    for c in n:
        if not (c.isalnum() or c in OKd or c in sub):
            s = s.replace(c, ' ')
    return sub.join(s.strip().split())


def attrs(inst, *names, **kwds):  # prec=6, fmt='F', ints=False, Nones=True, sep='='
    '''Get instance attributes as I{name=value} strings, with C{float}s handled like L{fstr}.

       @arg inst: The instance (any C{type}).
       @arg names: The attribute names (C{str}s).
       @kwarg kwds: Keyword argument for function L{pairs}, except
                    B{C{Nones=True}} to in-/exclude missing or
                    attributes with a C{None} I{value}.

       @return: A C{tuple(sep.join(t) for t in zip(names, reprs(values, ...)))}
                of C{str}s.
    '''
    Nones = _xkwds_pop(kwds, Nones=True)

    def items():
        for n in names:
            v = getattr(inst, n, None)
            if Nones or v is not None:
                yield n, v

    return pairs(items(), **kwds)


def enstr2(easting, northing, prec, *extras):
    '''Return easting, northing string representations.

       @arg easting: Easting from false easting (C{meter}).
       @arg northing: Northing from from false northing (C{meter}).
       @arg prec: Precision in number of digits (C{int}).
       @arg extras: Optional leading items (C{str}s).

       @return: B{C{extras}} + 2-Tuple C{(eastingStr, northingStr)}.

       @raise ValueError: Invalid B{C{prec}}.
    '''
    w = prec // 2
    try:
        p10 = (1e-4, 1e-3, 1e-2, 1e-1, 1)[w - 1]  # 10**(5 - w)
    except IndexError:
        raise _ValueError(prec=prec)
    return extras + ('%0*d' % (w, int(easting  * p10)),
                     '%0*d' % (w, int(northing * p10)))


def fstr(floats, prec=6, fmt=_Fmt, ints=False, sep=', ', strepr=None):
    '''Convert one or more floats to string, optionally stripped of trailing zero decimals.

       @arg floats: Single or a list, sequence, tuple, etc. (C{scalar}s).
       @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                    Trailing zero decimals are stripped if B{C{prec}} is
                    positive, but kept for negative B{C{prec}} values.
       @kwarg fmt: Optional, float format (C{str}).
       @kwarg ints: Optionally, remove the decimal dot (C{bool}).
       @kwarg sep: Separator joining the B{C{floats}} (C{str}).
       @kwarg strepr: Optional callable to format non-C{floats} (typically
                      C{repr}, C{str}) or C{None} to raise a TypeError.

       @return: The C{sep.join(strs(floats, ...)} joined (C{str}) or single
                C{strs((floats,), ...)} (C{str}) if B{C{floats}} is C{scalar}.
    '''
    if isscalar(floats):
        floats = (floats,)
    return sep.join(_streprs(prec, floats, fmt, ints, True, strepr))


def fstrzs(efstr, ap1z=False):
    '''Strip trailing zero decimals from a C{float} string.

       @arg efstr: Float with or without exponent (C{str}).
       @kwarg ap1z: Append the decimal point and one zero decimal
                    if the B{C{efstr}} is all digits (C{bool}).

       @return: Float (C{str}).
    '''
    s = efstr.find('.')
    if s >= 0:
        e = efstr.rfind('e')
        if e < 0:
            e = efstr.rfind('E')
            if e < 0:
                e = len(efstr)
        s += 2  # keep 1st '.0'
        if s < e and efstr[e-1] == '0':
            efstr = efstr[:s] + efstr[s:e].rstrip('0') + efstr[e:]

    elif ap1z:
        # %.G and %.g formats drop the decimal
        # point and all trailing zeros, ...
        if efstr.isdigit():
            efstr += '.0'  # ... append or ...
        else:  # ... insert one dot and zero
            e = efstr.rfind('e')
            if e < 0:
                e = efstr.rfind('E')
            if e > 0:
                efstr = efstr[:e] + '.0' + efstr[e:]

    return efstr


def instr(inst, *args, **kwds):
    '''Return the string representation of an instantiation.

       @arg inst: The instance (any C{type}).
       @arg args: Optional positional arguments.
       @kwarg kwds: Optional keyword arguments.

       @return: Representation (C{str}).
    '''
    from pygeodesy.named import classname
    return unstr(classname(inst), *args, **kwds)


def pairs(items, prec=6, fmt=_Fmt, ints=False, sep='='):
    '''Convert items to I{name=value} strings, with C{float}s handled like L{fstr}.

       @arg items: Name-value pairs (C{dict} or 2-{tuple}s of any C{type}s).
       @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                    Trailing zero decimals are stripped if B{C{prec}} is
                    positive, but kept for negative B{C{prec}} values.
       @kwarg fmt: Optional, float format (C{str}).
       @kwarg ints: Optionally, remove the decimal dot (C{bool}).
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
        raise _IsnotError('dict', '2-tuples', items=items)
    v = _streprs(prec, v, fmt, ints, False, repr)
    return tuple(sep.join(t) for t in zip(map(str, n), v))


def reprs(objs, prec=6, fmt=_Fmt, ints=False):
    '''Convert objects to C{repr} strings, with C{float}s handled like L{fstr}.

       @arg objs: List, sequence, tuple, etc. (any C{type}s).
       @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                    Trailing zero decimals are stripped if B{C{prec}} is
                    positive, but kept for negative B{C{prec}} values.
       @kwarg fmt: Optional, float format (C{str}).
       @kwarg ints: Optionally, remove the decimal dot (C{bool}).

       @return: A C{tuple(map(fstr|repr, objs))} of C{str}s.
    '''
    return tuple(_streprs(prec, objs, fmt, ints, False, repr)) if objs else ()


def strs(objs, prec=6, fmt=_Fmt, ints=False):
    '''Convert objects to C{str} strings, with C{float}s handled like L{fstr}.

       @arg objs: List, sequence, tuple, etc. (any C{type}s).
       @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                    Trailing zero decimals are stripped if B{C{prec}} is
                    positive, but kept for negative B{C{prec}} values.
       @kwarg fmt: Optional, float format (C{str}).
       @kwarg ints: Optionally, remove the decimal dot (C{bool}).

       @return: A C{tuple(map(fstr|str, objs))} of C{str}s.
    '''
    return tuple(_streprs(prec, objs, fmt, ints, False, str)) if objs else ()


def unstr(name, *args, **kwds):
    '''Return the string representation of an invokation.

       @arg name: Function, method or class name (C{str}).
       @arg args: Optional positional arguments.
       @kwarg kwds: Optional keyword arguments.

       @return: Representation (C{str}).
    '''
    t = reprs(args, fmt='g') + pairs(sorted(kwds.items()))
    return '%s(%s)' % (name, ', '.join(t))

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
