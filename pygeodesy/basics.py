
# -*- coding: utf-8 -*-

u'''Some, basic definitions and functions.
'''
# make sure int/int division yields float quotient
from __future__ import division
division = 1 / 2  # .albers, .azimuthal, etc., .utily
if not division:
    raise ImportError('%s 1/2 == %s' % ('division', division))
del division

from pygeodesy.errors import _AttributeError, _ImportError, _TypeError, \
                             _TypesError, _ValueError, _xkwds_get
from pygeodesy.interns import EPS0, MISSING, NEG0, NN, _by_, _DOT_, \
                             _invalid_, _N_A_, _name_, _not_finite_, \
                             _not_scalar_, _SPACE_, _UNDER_, _utf_8_, \
                             _version_, _0_0, _1_0
from pygeodesy.lazily import _ALL_LAZY, _FOR_DOCS

from copy import copy as _copy, deepcopy as _deepcopy
from math import copysign as _copysign, isinf, isnan
try:
    from math import isfinite as _isfinite  # in .ellipsoids, .fsums
except ImportError:  # Python 2-

    def _isfinite(x):  # mimick math.isfinite
        return not (isinf(x) or isnan(x))

__all__ = _ALL_LAZY.basics
__version__ = '22.03.28'

_below_     = 'below'
_ELLIPSIS4_ = '....'
_odd_       = 'odd'
_required_  = 'required'

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
            t = _SPACE_('cannot', t4[3] or _N_A_)
            del tb, t4
        return t

except NameError:  # Python 3+
    _Bytes = bytes, bytearray
    _Strs  = str,  # tuple
    _Xstr  = str


def clips(bstr, limit=50, white=NN):
    '''Clip a string to the given length limit.

       @arg bstr: String (C{bytes} or C{str}).
       @kwarg limit: Length limit (C{int}).
       @kwarg white: Optionally, replace all whitespace (C{str}).

       @return: The clipped or unclipped B{C{bstr}}.
    '''
    T = type(bstr)
    if len(bstr) > limit > 8:
        h = limit // 2
        bstr = T(_ELLIPSIS4_).join((bstr[:h], bstr[-h:]))
    if white:  # replace whitespace
        bstr = T(white).join(bstr.split())
    return bstr


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


def iscomplex(obj):
    '''Check whether an object is C{complex}.

       @arg obj: The object (any C{type}).

       @return: C{True} if B{C{obj}} is C{complex},
                C{False} otherwise.
    '''
    return isinstance(obj, complex)  # numbers.Complex?


def isfinite(obj):
    '''Check for C{Inf} and C{NaN} values.

       @arg obj: Value (C{scalar}).

       @return: C{False} if B{C{obj}} is C{INF} or C{NAN},
                C{True} otherwise.

       @raise TypeError: Non-scalar B{C{obj}}.
    '''
    try:
        return _isfinite(obj)
    except TypeError as x:
        raise _TypeError(_not_scalar_, obj, txt=str(x))
    except Exception as x:  # ValueError
        raise _ValueError(_not_finite_, obj, txt=str(x))


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
    return isinstance(obj, _Scalars)


def issequence(obj, *excluded):
    '''Check for sequence types.

       @arg obj: The object (any C{type}).
       @arg excluded: Optional exclusions (C{type}).

       @note: Excluding C{tuple} implies excluding C{namedtuple}.

       @return: C{True} if B{C{obj}} is a sequence, C{False} otherwise.
    '''
    return False if (excluded and isinstance(obj,  excluded)) else \
                                  isinstance(obj, _Seqs)


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


def ub2str(ub):
    '''Convert C{unicode} or C{bytes} to C{str}.
    '''
    if isinstance(ub, _Bytes):
        ub = str(ub.decode(_utf_8_))
    return ub


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
        import numpy
    except ImportError as x:
        raise _xImportError(x, where)
    return _xversion(numpy, where, *required)


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
        t = map2(int, package.__version__.split(_DOT_))
        if t[:n] < required:
            t = _SPACE_(package.__name__, _version_, _DOT_(*t),
                       _below_, _DOT_(*required),
                       _required_, _by_, _xwhere(where, **name))
            raise ImportError(t)
    return package


def _xwhere(where, **name):
    '''(INTERNAL) Get the fully qualified name.
    '''
    from pygeodesy.named import modulename
    m = modulename(where, prefixed=True)
    n = name.get(_name_, NN)
    if n:
        m = _DOT_(m, n)
    return m

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
