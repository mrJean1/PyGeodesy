
# -*- coding: utf-8 -*-

u'''Basic constants, definitions and functions.
'''
from pygeodesy.errors import _AttributeError, _IsnotError, \
                             _TypesError, _ValueError
from pygeodesy.interns import _COMMA_SPACE_, NN, _UNDERSCORE_
from pygeodesy.lazily import _ALL_LAZY, _FOR_DOCS

from copy import copy as _copy, deepcopy as _deepcopy
from inspect import isclass
from math import copysign, isinf, isnan, pi as PI
from sys import float_info as _float_info

__all__ = _ALL_LAZY.basics
__version__ = '20.07.28'

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
        raise ImportError  # _AssertionError
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
EPS2   = EPS * 2.0    #: M{EPS * 2}   ≈4.440892098501e-16 (C{float})

INF  = float('inf')  #: Infinity (C{float}), see function C{isinf}, C{isfinite}
NAN  = float('nan')  #: Not-A-Number (C{float}), see function C{isnan}
NEG0 = -0.0          #: Negative 0.0 (C{float}), see function C{isneg0}

PI2  = PI * 2.0  #: Two PI, M{PI * 2} aka Tau (C{float})  # PYCHOK expected
PI_2 = PI / 2.0  #: Half PI, M{PI / 2} (C{float})
PI_4 = PI / 4.0  #: Quarter PI, M{PI / 4} (C{float})

R_M  = 6371008.771415  #: Mean, spherical earth radius (C{meter}).


def _bkwds(inst, Error=AttributeError, **name_value_pairs):  # in .frechet, .hausdorff, .heights
    '''(INTERNAL) Set applicable C{bool} properties/attributes.
    '''
    for n, v in name_value_pairs.items():
        b = getattr(inst, n, None)
        if b is None:  # invalid bool attr
            t = n, v, inst.__class__.__name__  # XXX .classname
            raise Error('not applicable: %s=%r for %s' % t)
        if v in (False, True) and v != b:
            setattr(inst, _UNDERSCORE_ + n, v)


def clips(bstr, limit=50, white=NN):
    '''Clip a string to the given length limit.

       @arg bstr: String (C{bytes} or C{str}).
       @kwarg limit: Length limit (C{int}).
       @kwarg white: Optionally, replace all whitespace (C{str}).

       @return: The clipped or unclipped B{C{bstr}}.
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

       @return: 2-Tuple (_1st, _2nd) half (C{str}).

       @raise ValueError: Zero or odd C{len}(B{C{str2}}).
    '''
    h, r = divmod(len(str2), 2)
    if r or not h:
        raise _ValueError(str2=str2, txt='odd')
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
            raise _IsnotError(isscalar.__name__, obj=obj)
        return not (isinf(obj) or isnan(obj))


def isint(obj, both=False):
    '''Check for C{int} type or an integer C{float} value.

       @arg obj: The object (any C{type}).
       @kwarg both: Optionally, check C{float} type and value (C{bool}).

       @return: C{True} if B{C{obj}} is C{int} or an integer
                C{float}, C{False} otherwise.
    '''
    if both and isinstance(obj, float):  # NOT _Scalars!
        try:
            return obj.is_integer()
        except AttributeError:  # PYCHOK no cover
            return False  # XXX float(int(obj)) == obj?
    return isinstance(obj, _Ints)


def isneg0(obj):
    '''Check for L{NEG0}, negative 0.0.

       @arg obj: Value (C{scalar}).

       @return: C{True} if B{C{obj}} is C{NEG0} or -0.0,
                C{False} otherwise.
    '''
    return obj in (0.0, NEG0) and copysign(1, obj) < 0
#                             and str(obj).rstrip(_0_) == '-0.'


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


def issubclassof(Sub, Super):
    '''Check whether a class is a sub-class of a super class.

       @arg Sub: The sub-class (C{class}).
       @arg Super: The super class (C{class}).

       @return: C{True} if B{C{Sub}} is a sub-class of B{C{Super}},
                C{False} otherwise (C{bool}).
    '''
    return isclass(Sub) and isclass(Super) and issubclass(Sub, Super)


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


def map1(func, *xs):  # XXX map_
    '''Apply each argument to a single-argument function and
       return a C{tuple} of results.

       @arg func: Function to apply (C{callable}).
       @arg xs: Arguments to apply (C{any positional}).

       @return: Function results (C{tuple}).
    '''
    return tuple(map(func, xs))  # note, NO *xs


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
    return tuple(map(func, *xs))  # note, *xs


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
        t = 'get and set' if doc.startswith(' ') else NN
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
        if _FOR_DOCS and method.__doc__:
            self.__doc__ = method.__doc__
        self.name = method.__name__  # == self.fget.__name__

        # U{Descriptor HowTo Guide<https://docs.Python.org/3/howto/descriptor.html>}
        def immutable(inst, value):
            '''Throws an C{AttributeError}, always.
            '''
            t = immutable.__name__, inst, method.__name__, value
            raise _AttributeError('%s property: %r.%s = %r' % t)

        property.__init__(self, method, immutable, None, method.__doc__ or 'N/A')


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
        raise _AttributeError('.%s' % (a,), o)

    for a in attrs:
        s = _getattr(other, a)
        g = _getattr(insto, a)
        if (g is None and s is not None) or g != s:
            setattr(insto, a, s)  # not settable?
    return insto


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
       @kwarg name_value_pairs: One or more B{C{name=value}} pairs
                                with the C{value} to be checked.

       @raise TypeError: At least one of the B{C{name_value_pairs}}
                         is not any of the B{C{Types}}.
    '''
    for n, v in name_value_pairs.items():
        if not isinstance(v, Types):
            raise _TypesError(n, v, *Types)


def _xkwds(kwds, **dflts):
    '''(INTERNAL) Override C{dflts} with C{kwds}.
    '''
    d = dflts
    if kwds:
        d = _copy(d)
        d.update(kwds)
    return d


def _xsubclassof(Class, **name_value_pairs):
    '''(INTERNAL) Check super C{Class} of all C{name=value} pairs.

       @arg Class: A class or type (C{class}).
       @kwarg name_value_pairs: One or more B{C{name=value}} pairs
                                with the C{value} to be checked.

       @raise TypeError: At least one of the B{C{name_value_pairs}}
                         is not a sub-class of B{C{Class}}.
    '''
    for n, v in name_value_pairs.items():
        if not issubclassof(v, Class):
            raise _TypesError(n, v, Class)


def _xzipairs(lefts, rights, sep=_COMMA_SPACE_, fmt=NN, pair='%s:%s'):
    '''(INTERNAL) Zip C{lefts} and C{rights} into a C{str}.
    '''
    t = sep.join(pair % t for t in zip(lefts, rights))
    if fmt:
        t = fmt % (t,)
    return t

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
