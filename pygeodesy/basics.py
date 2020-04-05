
# -*- coding: utf-8 -*-

u'''Basic constants, definitions and functions.

'''
from pygeodesy.lazily import _ALL_LAZY

from copy import copy as _copy, deepcopy as _deepcopy
from inspect import isclass
from math import copysign, isinf, isnan
from sys import float_info as _float_info

# all public contants, classes and functions
__all__ = _ALL_LAZY.basics
__version__ = '20.04.02'

try:  # Luciano Ramalho, "Fluent Python", page 395, O'Reilly, 2016
    from numbers import Integral as _Ints  #: (INTERNAL) Int objects
except ImportError:  # PYCHOK no cover
    try:  # _Ints imported by .utily
        _Ints = int, long  #: (INTERNAL) Int objects (C{tuple})
    except NameError:  # Python 3+
        _Ints = int,  #: (INTERNAL) Int objects (C{tuple})

try:  # similarly ...
    from numbers import Real as _Scalars  #: (INTERNAL) Scalar objects
except ImportError:  # PYCHOK no cover
    try:
        _Scalars = int, long, float  #: (INTERNAL) Scalar objects (C{tuple})
    except NameError:
        _Scalars = int, float  #: (INTERNAL) Scalar objects (C{tuple})

try:
    try:  # use C{from collections.abc import ...} in Python 3.9+
        from collections.abc import Sequence as _Sequence  # imported by .points
    except ImportError:  # no .abc in Python 2.7-
        from collections import Sequence as _Sequence  # imported by .points
    if isinstance([], _Sequence) and isinstance((), _Sequence):
        #                        and isinstance(range(1), _Sequence):
        _Seqs = _Sequence
    else:
        raise ImportError  # AssertionError
except ImportError:  # PYCHOK no cover
    _Sequence = tuple  # immutable for .points._Basequence
    _Seqs     = list, _Sequence  # , range for function len2 below

try:
    _Strs = basestring, str
except NameError:  # Python 3+
    _Strs = str,

try:
    EPS    = _float_info.epsilon   #: System's epsilon (C{float})
    MANTIS = _float_info.mant_dig  #: System's mantissa bits (C{int})
    MAX    = _float_info.max       #: System's float max (C{float})
    MIN    = _float_info.min       #: System's float min (C{float})
except AttributeError:  # PYCHOK no cover
    EPS    = 2.220446049250313e-16  #: Epsilon (C{float}) 2**-52?
    MANTIS = 53  #: Mantissa bits ≈53 (C{int})
    MAX    = pow(2.0,  1023) * (2 - EPS)  #: Float max (C{float}) ≈10**308, 2**1024?
    MIN    = pow(2.0, -1021)  # Float min (C{float}) ≈10**-308, 2**-1021?
EPS_2  = EPS / 2.0    #: M{EPS / 2}   ≈1.110223024625e-16 (C{float})
EPS1   = 1.0 - EPS    #: M{1 - EPS}   ≈0.9999999999999998 (C{float})
EPS1_2 = 1.0 - EPS_2  #: M{1 - EPS_2} ≈0.9999999999999999 (C{float})
# _1EPS  = 1.0 + EPS  #: M{1 + EPS}   ≈1.0000000000000002 (C{float})

INF  = float('inf')  #: Infinity (C{float}), see function C{isinf}, C{isfinite}
NAN  = float('nan')  #: Not-A-Number (C{float}), see function C{isnan}
NEG0 = -0.0          #: Negative 0.0 (C{float}), see function C{isneg0}

OK = 'OK'  # OK for test like I{if ... is OK: ...}

# R_M moved here to avoid circular imports
R_M = 6371008.771415  #: Mean, spherical earth radius (C{meter}).

_limiterrors = True      # imported by .formy
_MISSING     = object()  # singleton, imported by .wgrs


class LenError(ValueError):
    '''Error raised for mis-matching C{len} values.
    '''
    def __init__(self, where, **lens):  # Error=ValueError
        '''New L{LenError}.

           @arg where: Object with C{.__name__} attribute (C{class}, C{method}, or C{function}).
           @kwarg lens: Two or more C{name=len(name)} pairs (C{keyword arguments}).
        '''
        ns, vs = zip(*sorted(lens.items()))
        ns = ', '.join(ns)
        vs = ' vs '.join(map(str, vs))
        t  = where.__name__, ns, 'len', vs
        ValueError.__init__(self, '%s(%s) %s: %s' % t)


class LimitError(ValueError):
    '''Error raised for lat- or longitudinal deltas exceeding
       the B{C{limit}} in functions L{equirectangular} and
       L{equirectangular_} and C{nearestOn*} and C{simplify*}
       functions or methods.
    '''
    pass


class _Adict(dict):
    '''(INTERNAL) Basic C{dict} with key I{and} attribute
       access to the items.
    '''
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            return dict.__getattr__(self, name)


def _bkwds(inst, kwds, Error):
    '''(INTERNAL) Set applicable C{bool} attributes.
    '''
    for k, v in kwds.items():
        b = getattr(inst, k, None)
        if b is None:
            t = k, v, inst.__class__.__name__
            raise Error('not applicable: %s=%r for %s' % t)
        if v in (False, True) and v != b:
            setattr(inst, '_' + k, v)


def clips(bstr, limit=50, white=''):
    '''Clip a string to the given length limit.

       @arg bstr: String (C{bytes} or C{str}).
       @kwarg limit: Length limit (C{int}).
       @kwarg white: Whitespace replacement (C{str}).

       @return: Un-/clipped B{C{bstr}}.
    '''
    if len(bstr) > limit > 8:
        h = limit // 2
        bstr = bstr[:h] + type(bstr)('....') + bstr[-h:]
    if white:  # replace whitespace
        bstr = type(bstr)(white).join(bstr.split())
    return bstr


def halfs2(str2):
    '''Split a string in 2 halfs.

       @arg str2: String to split (C{str}).

       @return: 2-Tuple (1st, 2nd) half (C{str}).

       @raise ValueError: Zero or odd C{len}(B{str2}).
    '''
    h, r = divmod(len(str2), 2)
    if r or not h:
        raise ValueError('%s invalid: %r' % ('str2', str2))
    return str2[:h], str2[h:]


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
            raise _isnotError(isscalar.__name__, obj=obj)
        return not (isinf(obj) or isnan(obj))


def isint(obj, both=False):
    '''Check for integer type or an integer C{float}.

       @arg obj: The object (any C{type}).
       @kwarg both: Optionally, check both type and value (C{bool}).

       @return: C{True} if B{C{obj}} is C{int}, C{False} otherwise.
    '''
    if both and isinstance(obj, float):  # NOT _Scalars!
        try:
            return obj.is_integer()
        except AttributeError:
            return False  # XXX float(int(obj)) == obj?
    return isinstance(obj, _Ints)


def isneg0(obj):
    '''Check for L{NEG0}, negative 0.0.

       @arg obj: Value (C{scalar}).

       @return: C{True} if B{C{obj}} is C{NEG0} or -0.0,
                C{False} otherwise.
    '''
    return obj in (0.0, NEG0) and copysign(1, obj) < 0
#                             and str(obj).rstrip('0') == '-0.'


def _isnotError(*names, **pair):  # Error=TypeError, name=value
    '''(INTERNAL) Format a C{TypeError} for a C{name=value} pair.
    '''
    Error = pair.pop('Error', TypeError)
    for n, v in pair.items():
        break
    else:
        n, v = 'pair', 'N/A'
    t = ' or ' .join(names)
    return Error('%s not %s: %r' % (n, t, v))


def isscalar(obj):
    '''Check for scalar types.

       @arg obj: The object (any C{type}).

       @return: C{True} if B{C{obj}} is C{scalar}, C{False} otherwise.
    '''
    return isinstance(obj, _Scalars)


def issequence(obj, *excluded):
    '''Check for sequence types.

       @arg obj: The object (any C{type}).
       @arg excluded: Optional, exclusions (C{type}).

       @note: Excluding C{tuple} implies excluding C{namedtuple}.

       @return: C{True} if B{C{obj}} is a sequence, C{False} otherwise.
    '''
    if excluded:
        return isinstance(obj, _Seqs) and not \
               isinstance(obj, excluded)
    else:
        return isinstance(obj, _Seqs)


def isstr(obj):
    '''Check for string types.

       @arg obj: The object (any C{type}).

       @return: C{True} if B{C{obj}} is C{str}, C{False} otherwise.
    '''
    return isinstance(obj, _Strs)


def issubclassof(sub, sup):
    '''Check whether a class is a subclass of a super class.

       @arg sub: The subclass (C{class}).
       @arg sup: The super class (C{class}).

       @return: C{True} if B{C{sub}} is a subclass of B{C{sup}}.
    '''
    return isclass(sub) and isclass(sup) and issubclass(sub, sup)


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


def limiterrors(raiser=None):
    '''Get/set the raising of limit errors.

       @kwarg raiser: Choose C{True} to throw or C{False} to
                      ignore L{LimitError} exceptions.  Use
                      C{None} to leave the setting unchanged.

       @return: Previous setting (C{bool}).
    '''
    global _limiterrors
    t = _limiterrors
    if raiser in (True, False):
        _limiterrors = raiser
    return t


def map1(func, *xs):  # XXX map_
    '''Apply each argument to a single-argument function and
       return a C{tuple} of results.

       @arg func: Function to apply (C{callable}).
       @arg xs: Arguments to apply (C{any positional}).

       @return: Function results (C{tuple}).
    '''
    return tuple(map(func, xs))


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
    return tuple(map(func, *xs))


def property_doc_(doc):
    '''Decorator for a property with documentation.

       @arg doc: The property documentation (C{str}).

       @example:

       >>> @property_doc_("documentation text.")
       >>> def name(self):
       >>>     ...
       >>>
       >>> @name.setter
       >>> def name(self, value):
       >>>     ...
    '''
    # See Luciano Ramalho, "Fluent Python", page 212ff, O'Reilly, 2016,
    # "Parameterized Decorators", especially Example 7-23.  Also, see
    # <https://Python-3-Patterns-Idioms-Test.ReadTheDocs.io/en/latest/PythonDecorators.html>

    def _property(method):
        '''(INTERNAL) Return C{method} as documented C{property.getter}.
        '''
        t = 'get and set' if doc.startswith(' ') else ''
        return property(method, None, None, 'Property to ' + t + doc)

    return _property


class property_RO(property):
    # No __doc__ on purpose

    def __init__(self, method):  # PYCHOK signature
        '''New immutable, read-only L{property_RO}.

           @arg method: The callable to be decorated as C{property.getter}.

           @note: Like standard Python C{property} without a C{property.setter},
                  but with a more descriptive error message when set.
        '''
        # U{Descriptor HowTo Guide<https://docs.Python.org/3/howto/descriptor.html>}
        def immutable(inst, value):
            '''Throws an C{AttributeError}, always.
            '''
            t = immutable.__name__, inst, method.__name__, value
            raise AttributeError('%s property: %r.%s = %r' % t)

        property.__init__(self, method, immutable, None,  method.__doc__ or 'N/A')


# def property_RO(method):  # OBSOLETE
#     '''An immutable property (C{Read Only}).
#
#        @arg method: The callable to be decorated as C{property.getter}.
#
#        @note: Like standard Python C{property} without a C{property.setter},
#               but with a more descriptive error message when set.
#     '''
#     def Read_Only(inst, value):
#         '''Throws an C{AttributeError}, always.
#         '''
#         t = Read_Only.__name__, inst, method.__name__, value
#         raise AttributeError('%s property: %r.%s = %r' % t)
#
#     return property(method, Read_Only, None, method.__doc__ or 'N/A')


def scalar(value, low=EPS, high=1.0, name='scalar', Error=ValueError):
    '''Validate a scalar.

       @arg value: The value (C{scalar}).
       @kwarg low: Optional lower bound (C{scalar}).
       @kwarg high: Optional upper bound (C{scalar}).
       @kwarg name: Optional name of value (C{str}).
       @kwarg Error: Exception to raise (C{ValueError}).

       @return: New value (C{type} of B{C{low}}).

       @raise TypeError: Non-scalar B{C{value}}.

       @raise Error: Out-of-bounds B{C{value}}.
    '''
    if not isscalar(value):
        raise _isnotError(scalar.__name__, **{name: value})
    try:
        if low is None:
            v = float(value)
        else:
            v = type(low)(value)
            if low > v or v > high:
                raise ValueError
    except (TypeError, ValueError):
        raise _isnotError('valid', Error=Error, **{name: value})
    return v


def splice(iterable, n=2, fill=_MISSING):
    '''Split an iterable into C{n} slices.

       @arg iterable: Items to be spliced (C{list}, C{tuple}, ...).
       @kwarg n: Number of slices to generate (C{int}).
       @kwarg fill: Fill value for missing items.

       @return: Generator of B{C{n}} slices M{iterable[i::n] for i=0..n}.

       @note: Each generated slice is a C{tuple} or a C{list},
              the latter only if the B{C{iterable}} is a C{list}.

       @raise ValueError: Non-C{int} or non-positive B{C{n}}.

       @example:

       >>> from pygeodesy import splice

       >>> a, b = splice(range(10))
       >>> a, b
       ((0, 2, 4, 6, 8), (1, 3, 5, 7, 9))

       >>> a, b, c = splice(range(10), n=3)
       >>> a, b, c
       ((0, 3, 6, 9), (1, 4, 7], [2, 5, 8))

       >>> a, b, c = splice(range(10), n=3, fill=-1)
       >>> a, b, c
       ((0, 3, 6, 9), (1, 4, 7, -1), (2, 5, 8, -1))

       >>> list(splice(range(12), n=5))
       [(0, 5, 10), (1, 6, 11), (2, 7), (3, 8), (4, 9)]

       >>> splice(range(9), n=1)
       <generator object splice at 0x0...>
    '''
    if not (isinstance(n, _Ints) and n > 0):
        raise ValueError('%s %s=%s' % ('splice', 'n', n))

    t = iterable
    if not isinstance(t, (list, tuple)):
        t = tuple(t)  # force tuple, also for PyPy3
    if n > 1:
        if fill is not _MISSING:
            m = len(t) % n
            if m > 0:  # fill with same type
                t += type(t)((fill,)) * (n - m)
        for i in range(n):
            yield t[i::n]  # [i:None:n] pychok -Tb ...
    else:
        yield t


def _TypeError(*Types, **pairs):
    '''(INTERNAL) Check C{Types} of all C{name=value} pairs.
    '''
    for n, v in pairs.items():
        if not isinstance(v, Types):
            t = ' or '.join(t.__name__ for t in Types)
            # first letter of Type name I{pronounced} as vowel
            a = 'an' if t[:1].lower() in 'aeinoux' else 'a'
            raise TypeError('%s not %s %s: %r' % (n, a, t, v))


def _xcopy(inst, deep=False):
    '''(INTERNAL) Copy an instance, shallow or deep.

       @arg inst: The instance to copy (C{_Named}).
       @kwarg deep: If C{True} make a deep, otherwise
                    shallow copy (C{bool}).

       @return: The copy (C{This class} or subclass thereof).
    '''
    return _deepcopy(inst) if deep else _copy(inst)


def _xkwds(kwds, **dflts):
    '''(INTERNAL) Override C{dflts} with C{kwds}.
    '''
    d = dflts
    if kwds:
        d = d.copy()
        d.update(kwds)
    return d

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
