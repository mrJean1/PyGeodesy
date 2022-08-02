
# -*- coding: utf-8 -*-

u'''Some, basic definitions, functions and dependencies.

Use env variable C{PYGEODESY_XPACKAGES} to avoid import of dependencies
C{geographiclib}, C{numpy} and/or C{scipy}.  Set C{PYGEODESY_XPACKAGES}
to a comma-separated list of package names to be excluded from import.
'''
# make sure int/int division yields float quotient
from __future__ import division
division = 1 / 2  # .albers, .azimuthal, etc., .utily
if not division:
    raise ImportError('%s 1/2 == %s' % ('division', division))
del division

from pygeodesy.errors import _AttributeError, _ImportError, _TypeError, \
                             _TypesError, _ValueError, _xError, _xkwds_get
from pygeodesy.interns import EPS0, INF, INT0, MISSING, NAN, NEG0, NINF, NN, \
                             _by_, _DOT_, _enquote, _EQUAL_, _in_, _INF_, \
                             _invalid_, _N_A_, _name_, _NAN_, _NINF_, _SPACE_, \
                             _splituple, _UNDER_, _version_, _0_0, _0_5, \
                             _1_0, _360_0
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS, _FOR_DOCS, \
                             _getenv, _sys_version_info2

from copy import copy as _copy, deepcopy as _deepcopy
from math import copysign as _copysign, isinf, isnan

__all__ = _ALL_LAZY.basics
__version__ = '22.08.02'

_below_               = 'below'
_cannot_              = 'cannot'
_ELLIPSIS4_           = '....'
_INF_NAN_NINF         = {INF: _INF_, NAN: _NAN_, NINF: _NINF_}
_odd_                 = 'odd'
_required_            = 'required'
_PYGEODESY_XPACKAGES_ = 'PYGEODESY_XPACKAGES'
_XPACKAGES            = _splituple(_getenv(_PYGEODESY_XPACKAGES_, NN))

try:  # Luciano Ramalho, "Fluent Python", page 395, O'Reilly, 2016
    from numbers import Integral as _Ints, Real as _Scalars
except ImportError:
    try:
        _Ints = int, long  # int objects (C{tuple})
    except NameError:  # Python 3+
        _Ints = int,  # int objects (C{tuple})
    _Scalars = _Ints + (float,)

try:
    try:  # use C{from collections.abc import ...} in Python 3.9+
        from collections.abc import Sequence as _Sequence  # in .points
    except ImportError:  # no .abc in Python 3.8- and 2.7-
        from collections import Sequence as _Sequence  # in .points
    if isinstance([], _Sequence) and isinstance((), _Sequence):
        #                        and isinstance(range(1), _Sequence):
        _Seqs = _Sequence
    else:
        raise ImportError  # _AssertionError
except ImportError:
    _Sequence = tuple  # immutable for .points._Basequence
    _Seqs     = list, _Sequence  # , range for function len2 below

try:
    _Bytes = unicode, bytearray  # PYCHOK expected
    _Strs  = basestring, str  # XXX , bytes

    def _NOP(x):
        '''NOP, pass thru.'''
        return x

    str2ub = ub2str = _NOP  # avoids UnicodeDecodeError

    def _Xstr(exc):  # PYCHOK no cover
        '''I{Invoke only with caught ImportError} B{C{exc}}.

           C{... "cannot import name _distributor_init" ...}

           only for C{numpy}, C{scipy} import errors occurring
           on arm64 Apple Silicon running macOS' Python 2.7.16?
        '''
        t = str(exc)
        if '_distributor_init' in t:
            from sys import exc_info
            from traceback import extract_tb
            tb = exc_info()[2]  # 3-tuple (type, value, traceback)
            t4 = extract_tb(tb, 1)[0]  # 4-tuple (file, line, name, 'import ...')
            t = _SPACE_(_cannot_, t4[3] or _N_A_)
            del tb, t4
        return t

except NameError:  # Python 3+
    from pygeodesy.interns import _utf_8_

    _Bytes = bytes, bytearray
    _Strs  = str,  # tuple
    _Xstr  = str

    def str2ub(sb):
        '''Convert C{str} to C{unicode bytes}.
        '''
        if isinstance(sb, _Strs):
            sb = sb.encode(_utf_8_)
        return sb

    def ub2str(ub):
        '''Convert C{unicode bytes} to C{str}.
        '''
        if isinstance(ub, _Bytes):
            ub = str(ub.decode(_utf_8_))
        return ub


def clips(sb, limit=50, white=NN):
    '''Clip a string to the given length limit.

       @arg sb: String (C{str} or C{bytes}).
       @kwarg limit: Length limit (C{int}).
       @kwarg white: Optionally, replace all whitespace (C{str}).

       @return: The clipped or unclipped B{C{sb}}.
    '''
    T = type(sb)
    if len(sb) > limit > 8:
        h  = limit // 2
        sb = T(_ELLIPSIS4_).join((sb[:h], sb[-h:]))
    if white:  # replace whitespace
        sb = T(white).join(sb.split())
    return sb


def copysign0(x, y):
    '''Like C{math.copysign(x, y)} except C{zero}, I{unsigned}.

       @return: C{math.copysign(B{x}, B{y})} if B{C{x}} else
                C{type(B{x})(0)}.
    '''
    return _copysign(x, (y if y else 0)) if x else copytype(0, x)


def copytype(x, y):
    '''Return the value of B{x} as C{type} of C{y}.

       @return: C{type(B{y})(B{x})}.
    '''
    return type(y)(x if x else _0_0)


def halfs2(str2):
    '''Split a string in 2 halfs.

       @arg str2: String to split (C{str}).

       @return: 2-Tuple (_1st, _2nd) half (C{str}).

       @raise ValueError: Zero or odd C{len}(B{C{str2}}).
    '''
    h, r = divmod(len(str2), 2)
    if r or not h:
        raise _ValueError(str2=str2, txt=_odd_)
    return str2[:h], str2[h:]


def isbool(obj):
    '''Check whether an object is C{bool}ean.

       @arg obj: The object (any C{type}).

       @return: C{True} if B{C{obj}} is C{bool}ean,
                C{False} otherwise.
    '''
    return isinstance(obj, bool)  # and (obj is True
#                                     or obj is False)


if _FOR_DOCS:  # XXX avoid epidoc Python 2.7 error
    from inspect import isclass as _isclass

    def isclass(obj):
        '''Return C{True} if B{C{obj}} is a C{class}.

           @see: Python's C{inspect.isclass}.
        '''
        return _isclass(obj)
else:
    from inspect import isclass  # PYCHOK re-import


try:
    from math import isclose as _isclose
except ImportError:  # Python 3.4-

    def _isclose(a, b, rel_tol=1e-9, abs_tol=0):
        '''Mimick Python 3.5+ C{math.isclose}.
        '''
        t, d = abs_tol, abs(a - b)
        if d > t:
            r = max(abs(a), abs(b)) * rel_tol
            t = max(r, t)
        return d <= t


def isclose(a, b, rel_tol=1e-12, abs_tol=EPS0):
    '''Like C{math.isclose}, but with defaults such
       that C{isclose(0, EPS0)} is C{True} by default.
    '''
    return _isclose(a, b, rel_tol=rel_tol, abs_tol=abs_tol)


def iscomplex(obj):
    '''Check whether an object is C{complex}.

       @arg obj: The object (any C{type}).

       @return: C{True} if B{C{obj}} is C{complex},
                C{False} otherwise.
    '''
    # hasattr('conjugate'), hasattr('real') and hasattr('imag')
    return isinstance(obj, complex)  # numbers.Complex?


try:
    from math import isfinite as _isfinite  # in .ellipsoids, .fsums, .karney
except ImportError:  # Python 3.1-

    def _isfinite(x):
        '''Mimick Python 3.2+ C{math.isfinite}.
        '''
        return not (isinf(x) or isnan(x))


def isfinite(obj):
    '''Check a finite C{scalar} or C{complex} value.

       @arg obj: Value (C{scalar} or C{complex}).

       @return: C{False} if B{C{obj}} is C{INF}, C{NINF}
                or C{NAN}, C{True} otherwise.

       @raise TypeError: Non-scalar and non-complex B{C{obj}}.
    '''
    try:
        return (obj not in _INF_NAN_NINF) and _isfinite(obj)
    except Exception as x:
        if iscomplex(obj):  # _isfinite(complex) thows TypeError
            return isfinite(obj.real) and isfinite(obj.imag)
        P = _MODS.streprs.Fmt.PAREN
        raise _xError(x, P(isfinite.__name__, obj))


try:
    isidentifier = str.isidentifier  # Python 3, must be str
except AttributeError:  # Python 2-

    def isidentifier(obj):
        '''Return C{True} if B{C{obj}} is a valid Python identifier.
        '''
        return bool(obj and obj.replace(_UNDER_, NN).isalnum()
                        and not obj[:1].isdigit())

# from math import isinf


def isint(obj, both=False):
    '''Check for C{int} type or an integer C{float} value.

       @arg obj: The object (any C{type}).
       @kwarg both: If C{true}, check C{float} and L{Fsum}
                    type and value (C{bool}).

       @return: C{True} if B{C{obj}} is C{int} or an integer
                C{float} or L{Fsum}, C{False} otherwise.

       @note: C{isint(True)} or C{isint(False)} returns
              C{False} (and no longer C{True}).
    '''
    if isinstance(obj, _Ints) and not isbool(obj):
        return True
    elif both:  # and isinstance(obj, (float, Fsum)):  # NOT _Scalars!
        try:
            return obj.is_integer()
        except AttributeError:
            pass  # XXX float(int(obj)) == obj?
    return False


def isint0(obj, both=False):
    '''Check for L{INT0} or C{int(0)} value.

       @arg obj: The object (any C{type}).
       @kwarg both: If C{true}, also check C{float(0)} (C{bool}).

       @return: C{True} if B{C{obj}} is L{INT0}, C{int(0)} or
                C{float(0)}, C{False} otherwise.
    '''
    return (obj is INT0 or obj is int(0) or bool(both and
       (not obj) and isint(obj, both=True))) and not isbool(obj)


try:
    from keyword import iskeyword  # Python 2.7+
except ImportError:

    def iskeyword(unused):
        '''Not Implemented.  Return C{False}, always.
        '''
        return False

# from math import isnan


def isnear0(x, eps0=EPS0):
    '''Is B{C{x}} near zero?

       @arg x: Value (C{scalar}).
       @kwarg eps0: Near-zero (C{EPS0}).

       @return: C{True} if C{abs(B{x}) < B{eps0}},
                C{False} otherwise.

       @see: Function L{isnon0}.
    '''
    return eps0 > x > -eps0


def isnear1(x, eps0=EPS0):
    '''Is B{C{x}} near one?

       @arg x: Value (C{scalar}).
       @kwarg eps0: Near-zero (C{EPS0}).

       @return: C{isnear0(B{x} - 1)}.

       @see: Function L{isnear0}.
    '''
    return isnear0(x - _1_0, eps0=eps0)


def isneg0(x):
    '''Check for L{NEG0}, negative C{0.0}.

       @arg x: Value (C{scalar}).

       @return: C{True} if B{C{x}} is C{NEG0} or C{-0.0},
                C{False} otherwise.
    '''
    return x in (_0_0, NEG0) and _copysign(1, x) < 0
#                            and str(x).startswith(_MINUS_)


def isninf(x):
    '''Check for L{NINF}, negative C{INF}.

       @arg x: Value (C{scalar}).

       @return: C{True} if B{C{x}} is C{NINF} or C{-inf},
                C{False} otherwise.
    '''
    return x is NINF or ((not isfinite(x)) and x < 0)


def isnon0(x, eps0=EPS0):
    '''Is B{C{x}} non-zero?

       @arg x: Value (C{scalar}).
       @kwarg eps0: Near-zero (C{EPS0}).

       @return: C{True} if C{abs(B{x}) > B{eps0}},
                C{False} otherwise.

       @see: Function L{isnear0}.
    '''
    return eps0 <= x or x <= -eps0  # not isnear0


def isodd(x):
    '''Is B{C{x}} odd?

       @arg x: Value (C{scalar}).

       @return: C{True} if B{C{x}} is odd,
                C{False} otherwise.
    '''
    return bool(int(x) & 1)


def isscalar(obj):
    '''Check for scalar types.

       @arg obj: The object (any C{type}).

       @return: C{True} if B{C{obj}} is C{scalar}, C{False} otherwise.
    '''
    return isinstance(obj, _Scalars) and not isbool(obj)


def issequence(obj, *excls):
    '''Check for sequence types.

       @arg obj: The object (any C{type}).
       @arg excls: Classes to exclude (C{type}), all positional.

       @note: Excluding C{tuple} implies excluding C{namedtuple}.

       @return: C{True} if B{C{obj}} is a sequence, C{False} otherwise.
    '''
    return isinstance(obj, _Seqs) and not (excls and isinstance(obj, excls))


def isstr(obj):
    '''Check for string types.

       @arg obj: The object (any C{type}).

       @return: C{True} if B{C{obj}} is C{str}, C{False} otherwise.
    '''
    return isinstance(obj, _Strs)


def issubclassof(Sub, *Supers):
    '''Check whether a class is a sub-class of some other class(es).

       @arg Sub: The sub-class (C{class}).
       @arg Supers: One or more C(super) classes (C{class}).

       @return: C{True} if B{C{Sub}} is a sub-class of any B{C{Supers}},
                C{False} if not (C{bool}) or C{None} if B{C{Sub}} is not
                a class or if no B{C{Supers}} are given or none of those
                are a class.
    '''
    if isclass(Sub):
        t = tuple(S for S in Supers if isclass(S))
        if t:
            return bool(issubclass(Sub, t))
    return None


def istuplist(obj, minum=0):
    '''Check for tuple or list types and minumal length.

       @arg obj: The object (any C{type}).
       @kwarg minum: Minimal C{len} required C({int}).

       @return: C{True} if B{C{obj}} is C{tuple} or C{list} with
                C{len} at least B{C{minum}}, C{False} otherwise.
    '''
    return type(obj) in (tuple, list) and len(obj) >= (minum or 0)


def len2(items):
    '''Make built-in function L{len} work for generators, iterators,
       etc. since those can only be started exactly once.

       @arg items: Generator, iterator, list, range, tuple, etc.

       @return: 2-Tuple C{(n, items)} of the number of items (C{int})
                and the items (C{list} or C{tuple}).
    '''
    if not isinstance(items, _Seqs):  # NOT hasattr(items, '__len__'):
        items = list(items)
    return len(items), items


def map1(fun1, *xs):  # XXX map_
    '''Apply each argument to a single-argument function and
       return a C{tuple} of results.

       @arg fun1: 1-Arg function to apply (C{callable}).
       @arg xs: Arguments to apply (C{any positional}).

       @return: Function results (C{tuple}).
    '''
    return tuple(map(fun1, xs))  # note xs, not *xs


def map2(func, *xs):
    '''Apply arguments to a function and return a C{tuple} of results.

       Unlike Python 2's built-in L{map}, Python 3+ L{map} returns a
       L{map} object, an iterator-like object which generates the
       results only once.  Converting the L{map} object to a tuple
       maintains the Python 2 behavior.

       @arg func: Function to apply (C{callable}).
       @arg xs: Arguments to apply (C{list, tuple, ...}).

       @return: Function results (C{tuple}).
    '''
    return tuple(map(func, *xs))  # note *xs, not xs


def neg(x):
    '''Negate C{x} unless C{zero} or C{NEG0}.

       @return: C{-B{x}} if B{C{x}} else C{0.0}.
    '''
    return -x if x else _0_0


def neg_(*xs):
    '''Negate all of C{xs} with L{neg}.

       @return: A C{tuple(neg(x) for x in B{xs})}.
    '''
    return tuple(map(neg, xs))  # like map1


try:
    from math import remainder
except ImportError:  # Python 3.6-
    from math import fmod as _fmod

    def remainder(x, y):
        '''Mimick Python 3.7+ C{math.remainder}.
        '''
        if isnan(y):
            x =  NAN
        elif x and not isnan(x):
            y =  abs(y)
            x = _fmod(x, y)
            h = _0_5 * y
            if x < -h:
                x += y
            elif x >= h:
                x -= y
        return x  # keep signed 0.0


def signBit(x):
    '''Return C{signbit(B{x})}, like C++.

       @return: C{True} if C{B{x} < 0} or C{NEG0} (C{bool}).
    '''
    return x < 0 or isneg0(x)


def _signOf(x, off):
    '''(INTERNAL) Return the sign of B{C{x}} versus B{C{off}}.
    '''
    return 1 if x > off else (-1 if x < off else 0)


def signOf(x):
    '''Return sign of C{x} as C{int}.

       @return: -1, 0 or +1 (C{int}).
    '''
    try:
        s = x.signOf()  # Fsum instance?
    except AttributeError:
        s = _signOf(x, 0)
    return s


def splice(iterable, n=2, **fill):
    '''Split an iterable into C{n} slices.

       @arg iterable: Items to be spliced (C{list}, C{tuple}, ...).
       @kwarg n: Number of slices to generate (C{int}).
       @kwarg fill: Optional fill value for missing items.

       @return: A generator for each of B{C{n}} slices,
                M{iterable[i::n] for i=0..n}.

       @raise TypeError: Invalid B{C{n}}.

       @note: Each generated slice is a C{tuple} or a C{list},
              the latter only if the B{C{iterable}} is a C{list}.

       @example:

        >>> from pygeodesy import splice

        >>> a, b = splice(range(10))
        >>> a, b
        ((0, 2, 4, 6, 8), (1, 3, 5, 7, 9))

        >>> a, b, c = splice(range(10), n=3)
        >>> a, b, c
        ((0, 3, 6, 9), (1, 4, 7), (2, 5, 8))

        >>> a, b, c = splice(range(10), n=3, fill=-1)
        >>> a, b, c
        ((0, 3, 6, 9), (1, 4, 7, -1), (2, 5, 8, -1))

        >>> tuple(splice(list(range(9)), n=5))
        ([0, 5], [1, 6], [2, 7], [3, 8], [4])

        >>> splice(range(9), n=1)
        <generator object splice at 0x0...>
    '''
    if not isint(n):
        raise _TypeError(n=n)

    t = iterable
    if not isinstance(t, (list, tuple)):
        t = tuple(t)  # force tuple, also for PyPy3

    if n > 1:
        if fill:
            fill = _xkwds_get(fill, fill=MISSING)
            if fill is not MISSING:
                m = len(t) % n
                if m > 0:  # fill with same type
                    t += type(t)((fill,)) * (n - m)
        for i in range(n):
            # XXX t[i::n] chokes PyChecker
            yield t[slice(i, None, n)]
    else:
        yield t


def _umod_360(deg):
    '''(INTERNAL) Non-negative C{deg} modulo 360, basic C{.utily.wrap360}.
    '''
    return (deg % _360_0) or _0_0


def unsigned0(x):
    '''Return C{0.0} unsigned.

       @return: C{B{x}} if B{C{x}} else C{0.0}.
    '''
    return x if x else _0_0


def _xcopy(inst, deep=False):
    '''(INTERNAL) Copy an object, shallow or deep.

       @arg inst: The object to copy (any C{type}).
       @kwarg deep: If C{True} make a deep, otherwise
                    a shallow copy (C{bool}).

       @return: The copy of B{C{inst}}.
    '''
    return _deepcopy(inst) if deep else _copy(inst)


def _xdup(inst, **items):
    '''(INTERNAL) Duplicate an object, replacing some attributes.

       @arg inst: The object to copy (any C{type}).
       @kwarg items: Attributes to be changed (C{any}).

       @return: Shallow duplicate of B{C{inst}} with modified
                attributes, if any B{C{items}}.

       @raise AttributeError: Some B{C{items}} invalid.
    '''
    d = _xcopy(inst, deep=False)
    for n, v in items.items():
        if not hasattr(d, n):
            from pygeodesy.named import classname
            t = _SPACE_(_DOT_(classname(inst), n), _invalid_)
            raise _AttributeError(txt=t, this=inst, **items)
        setattr(d, n, v)
    return d


def _xgeographiclib(where, *required):
    '''(INTERNAL) Import C{geographiclib} and check required version
    '''
    try:
        _xpackage(_xgeographiclib)
        import geographiclib
    except ImportError as x:
        raise _xImportError(x, where)
    return _xversion(geographiclib, where, *required)


def _xImportError(x, where, **name):
    '''(INTERNAL) Embellish an C{ImportError}.
    '''
    t = _SPACE_(_required_, _by_, _xwhere(where, **name))
    return _ImportError(_Xstr(x), txt=t)


def _xinstanceof(*Types, **name_value_pairs):
    '''(INTERNAL) Check C{Types} of all C{name=value} pairs.

       @arg Types: One or more classes or types (C{class}),
                   all positional.
       @kwarg name_value_pairs: One or more C{B{name}=value} pairs
                                with the C{value} to be checked.

       @raise TypeError: One of the B{C{name_value_pairs}} is not
                         an instance of any of the B{C{Types}}.
    '''
    for n, v in name_value_pairs.items():
        if not isinstance(v, Types):
            raise _TypesError(n, v, *Types)


def _xnumpy(where, *required):
    '''(INTERNAL) Import C{numpy} and check required version
    '''
    try:
        _xpackage(_xnumpy)
        import numpy
    except ImportError as x:
        raise _xImportError(x, where)
    return _xversion(numpy, where, *required)


def _xpackage(_xpkg):
    '''(INTERNAL) Check dependency to be excluded.
    '''
    n = _xpkg.__name__[2:]
    if n in _XPACKAGES:
        x = _SPACE_(n, _in_, _PYGEODESY_XPACKAGES_)
        e = _enquote(_getenv(_PYGEODESY_XPACKAGES_, NN))
        raise ImportError(_EQUAL_(x, e))


def _xor(x, *xs):
    '''(INTERNAL) Exclusive-or C{x} and C{xs}.
    '''
    for x_ in xs:
        x ^= x_
    return x


def _xscipy(where, *required):
    '''(INTERNAL) Import C{scipy} and check required version
    '''
    try:
        _xpackage(_xscipy)
        import scipy
    except ImportError as x:
        raise _xImportError(x, where)
    return _xversion(scipy, where, *required)


def _xsubclassof(*Classes, **name_value_pairs):
    '''(INTERNAL) Check (super) class of all C{name=value} pairs.

       @arg Classes: One or more classes or types (C{class}),
                     all positional.
       @kwarg name_value_pairs: One or more C{B{name}=value} pairs
                                with the C{value} to be checked.

       @raise TypeError: One of the B{C{name_value_pairs}} is not
                         a (sub-)class of any of the B{C{Classes}}.
    '''
    for n, v in name_value_pairs.items():
        if not issubclassof(v, *Classes):
            raise _TypesError(n, v, *Classes)


def _xversion(package, where, *required, **name):  # in .karney
    '''(INTERNAL) Check the C{package} version vs B{C{required}}.
    '''
    n = len(required)
    if n:
        t = _xversion_info(package)
        if t[:n] < required:
            t = _SPACE_(package.__name__, _version_, _DOT_(*t),
                       _below_, _DOT_(*required),
                       _required_, _by_, _xwhere(where, **name))
            raise ImportError(t)
    return package


def _xversion_info(package):  # in .karney
    '''(INTERNAL) Get the C{package.__version_info__} as a 2- or
       3-tuple C{(major, minor, revision)} if C{int}s.
    '''
    try:
        t = package.__version_info__
    except AttributeError:
        t = package.__version__.strip()
        t = t.replace(_DOT_, _SPACE_).split()[:3]
    return map2(int, t)


def _xwhere(where, **name):
    '''(INTERNAL) Get the fully qualified name.
    '''
    m = _MODS.named.modulename(where, prefixed=True)
    n =  name.get(_name_, NN)
    if n:
        m = _DOT_(m, n)
    return m


if _sys_version_info2 < (3, 10):
    _zip = zip  # PYCHOK exported
else:  # Python 3.10+

    def _zip(*args):
        return zip(*args, strict=True)

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
