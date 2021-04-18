
# -*- coding: utf-8 -*-

u'''Basic definitions and functions.
'''
# make sure int/int division yields float quotient
from __future__ import division
division = 1 / 2  # .albers, .datums, .ellipsoidalVincenty, .ellipsoids,
if not division:  # .elliptic, .etm, .fmath, .formy, .lcc, .osgr, .utily
    raise ImportError('%s 1/2 == %s' % ('division', division))
del division

from pygeodesy.errors import _IsnotError, _TypesError, _ValueError
from pygeodesy.interns import NEG0, NN, _by_, _DOT_, _name_, \
                             _SPACE_, _UNDER_, _utf_8_, _version_, _0_0
from pygeodesy.lazily import _ALL_LAZY, _FOR_DOCS

from copy import copy as _copy, deepcopy as _deepcopy
from inspect import isclass as _isclass
from math import copysign as _copysign, isinf, isnan

__all__ = _ALL_LAZY.basics
__version__ = '21.03.28'

try:  # Luciano Ramalho, "Fluent Python", page 395, O'Reilly, 2016
    from numbers import Integral as _Ints  # int objects
except ImportError:  # PYCHOK no cover
    try:
        _Ints = int, long  # int objects (C{tuple})
    except NameError:  # Python 3+
        _Ints = int,  # int objects (C{tuple})

try:  # similarly ...
    from numbers import Real as _Scalars  # scalar objects
except ImportError:  # PYCHOK no cover
    try:
        _Scalars = int, long, float  # scalar objects (C{tuple})
    except NameError:
        _Scalars = int, float  # scalar objects (C{tuple})

try:
    try:  # use C{from collections.abc import ...} in Python 3.9+
        from collections.abc import Sequence as _Sequence  # imported by .points
    except ImportError:  # no .abc in Python 2.7-
        from collections import Sequence as _Sequence  # imported by .points
    if isinstance([], _Sequence) and isinstance((), _Sequence):
        #                        and isinstance(range(1), _Sequence):
        _Seqs = _Sequence
    else:
        raise ImportError  # _AssertionError
except ImportError:  # PYCHOK no cover
    _Sequence = tuple  # immutable for .points._Basequence
    _Seqs     = list, _Sequence  # , range for function len2 below

try:
    _Bytes = unicode, bytearray  # PYCHOK expected
    _Strs  = basestring, str
except NameError:  # Python 3+
    _Bytes = bytes, bytearray
    _Strs  = str,


def clips(bstr, limit=50, white=NN):
    '''Clip a string to the given length limit.

       @arg bstr: String (C{bytes} or C{str}).
       @kwarg limit: Length limit (C{int}).
       @kwarg white: Optionally, replace all whitespace (C{str}).

       @return: The clipped or unclipped B{C{bstr}}.
    '''
    if len(bstr) > limit > 8:
        h = limit // 2
        bstr = NN(bstr[:h], type(bstr)('....'), bstr[-h:])
    if white:  # replace whitespace
        bstr = type(bstr)(white).join(bstr.split())
    return bstr


def copysign(x, y):
    '''Like C{math.copysign(x, y)} except C{zero}, I{unsigned} .

       @return: C{math.copysign(B{x}, B{y})} if B{C{x}} else C{0.0}.
    '''
    return _copysign(x, y) if x else _0_0


def halfs2(str2):
    '''Split a string in 2 halfs.

       @arg str2: String to split (C{str}).

       @return: 2-Tuple (_1st, _2nd) half (C{str}).

       @raise ValueError: Zero or odd C{len}(B{C{str2}}).
    '''
    h, r = divmod(len(str2), 2)
    if r or not h:
        raise _ValueError(str2=str2, txt='odd')
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
    def isclass(obj):
        '''Return C{True} if B{C{obj}} is a C{class}.

           @see: Python's C{inspect.isclass}.
        '''
        return _isclass(obj)
else:
    isclass = _isclass


try:
    from math import isfinite  # new in Python 3+
except ImportError:

    def isfinite(obj):
        '''Check for C{Inf} and C{NaN} values.

           @arg obj: Value (C{scalar}).

           @return: C{False} if B{C{obj}} is C{INF} or C{NAN},
                    C{True} otherwise.

           @raise TypeError: Non-scalar B{C{obj}}.
        '''
        if not isscalar(obj):
            raise _IsnotError(isscalar.__name__, obj=obj)
        return not (isinf(obj) or isnan(obj))


try:
    isidentifier = str.isidentifier  # Python 3, must be str
except AttributeError:  # Python 2-

    def isidentifier(obj):
        '''Return C{True} if B{C{obj}} is a valid Python identifier.
        '''
        return True if (obj and obj.replace(_UNDER_, NN).isalnum()
                            and not obj[:1].isdigit()) else False


def isint(obj, both=False):
    '''Check for C{int} type or an integer C{float} value.

       @arg obj: The object (any C{type}).
       @kwarg both: Optionally, check C{float} type and value (C{bool}).

       @return: C{True} if B{C{obj}} is C{int} or an integer
                C{float}, C{False} otherwise.

       @note: C{isint(True)} or C{isint(False)} returns C{False} (and
              no longer C{True}).
    '''
    if both and isinstance(obj, float):  # NOT _Scalars!
        try:
            return obj.is_integer()
        except AttributeError:  # PYCHOK no cover
            return False  # XXX float(int(obj)) == obj?
    return isinstance(obj, _Ints) and not isbool(obj)


try:
    from keyword import iskeyword  # Python 2.7+
except ImportError:

    def iskeyword(unused):
        '''Not Implemented.  Return C{False}, always.
        '''
        return False


def isneg0(obj):
    '''Check for L{NEG0}, negative C{0.0}.

       @arg obj: Value (C{scalar}).

       @return: C{True} if B{C{obj}} is C{NEG0} or C{-0.0},
                C{False} otherwise.
    '''
    return obj in (_0_0, NEG0) and _copysign(1, obj) < 0
#                              and str(obj).startswith('-')


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
    '''Check whether a class is a sub-class of some class(es).

       @arg Sub: The sub-class (C{class}).
       @arg Supers: One or more C(super) classes (C{class}).

       @return: C{True} if B{C{Sub}} is a sub-class of any
                B{C{Supers}}, C{False} otherwise (C{bool}).
    '''
    if isclass(Sub):
        for S in Supers:  # any()
            if isclass(S) and issubclass(Sub, S):
                return True
    return False


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
       maintains Python 2 behavior.

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
    return tuple(map(neg, xs))  # see map1


def ub2str(ub):
    '''Convert C{unicode} or C{bytes} to C{str}.
    '''
    if isinstance(ub, _Bytes):
        ub = str(ub.decode(_utf_8_))
    return ub


def _xcopy(inst, deep=False):
    '''(INTERNAL) Copy an object, shallow or deep.

       @arg inst: The object to copy (any C{type}).
       @kwarg deep: If C{True} make a deep, otherwise
                    a shallow copy (C{bool}).

       @return: The copy of B{C{inst}}.
    '''
    return _deepcopy(inst) if deep else _copy(inst)


def _xinstanceof(*Types, **name_value_pairs):
    '''(INTERNAL) Check C{Types} of all C{name=value} pairs.

       @arg Types: One or more classes or types (C{class}).
       @kwarg name_value_pairs: One or more C{B{name}=value} pairs
                                with the C{value} to be checked.

       @raise TypeError: At least one of the B{C{name_value_pairs}}
                         is not any of the B{C{Types}}.
    '''
    for n, v in name_value_pairs.items():
        if not isinstance(v, Types):
            raise _TypesError(n, v, *Types)


def _xnumpy(where, *required):
    '''(INTERNAL) Import C{numpy} and check required version
    '''
    import numpy
    return _xversion(numpy, where, *required)


def _xscipy(where, *required):
    '''(INTERNAL) Import C{scipy} and check required version
    '''
    import scipy
    return _xversion(scipy, where, *required)


def _xsubclassof(Class, **name_value_pairs):
    '''(INTERNAL) Check super C{Class} of all C{name=value} pairs.

       @arg Class: A class or type (C{class}).
       @kwarg name_value_pairs: One or more C{B{name}=value} pairs
                                with the C{value} to be checked.

       @raise TypeError: At least one of the B{C{name_value_pairs}}
                         is not a sub-class of B{C{Class}}.
    '''
    for n, v in name_value_pairs.items():
        if not issubclassof(v, Class):
            raise _TypesError(n, v, Class)


def _xversion(package, where, *required, **name):  # in .karney
    '''(INTERNAL) Check the C{package} version vs B{C{required}}.
    '''
    t = map2(int, package.__version__.split(_DOT_)[:2])
    if t < required:
        from pygeodesy.named import modulename
        m = modulename(where, prefixed=True)
        n = name.get(_name_, NN)
        if n:
            m = _DOT_(m, n)
        t = _SPACE_(package.__name__, _version_, _DOT_.join_(*t),
                   'below', _DOT_.join_(*required),
                   'required', _by_, m)
        raise ImportError(t)
    return package

# **) MIT License
#
# Copyright (C) 2016-2021 -- mrJean1 at Gmail -- All Rights Reserved.
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
