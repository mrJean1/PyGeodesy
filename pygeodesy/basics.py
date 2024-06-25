
# -*- coding: utf-8 -*-

u'''Some, basic definitions, functions and dependencies.

Use env variable C{PYGEODESY_XPACKAGES} to avoid import of dependencies
C{geographiclib}, C{numpy} and/or C{scipy}.  Set C{PYGEODESY_XPACKAGES}
to a comma-separated list of package names to be excluded from import.
'''
# make sure int/int division yields float quotient
from __future__ import division
division = 1 / 2  # .albers, .azimuthal, .constants, etc., .utily
if not division:
    raise ImportError('%s 1/2 == %s' % ('division', division))
del division

# from pygeodesy.cartesianBase import CartesianBase  # _MODS
# from pygeodesy.constants import isneg0, NEG0  # _MODS
from pygeodesy.errors import _AttributeError, _ImportError, _NotImplementedError, \
                             _TypeError, _TypesError, _ValueError, _xAssertionError, \
                             _xkwds_get1
from pygeodesy.internals import _0_0, _enquote, _passarg, _version_info
from pygeodesy.interns import MISSING, NN, _1_, _by_, _COMMA_, _DOT_, _DEPRECATED_, \
                             _ELLIPSIS4_, _EQUAL_, _in_, _invalid_, _N_A_, _not_, \
                             _not_scalar_, _odd_, _SPACE_, _UNDER_, _version_
# from pygeodesy.latlonBase import LatLonBase  # _MODS
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS, _FOR_DOCS, _getenv, \
                             LazyImportError, _sys_version_info2
# from pygeodesy.named import classname, modulename, _name__  # _MODS
# from pygeodesy.nvectorBase import NvectorBase  # _MODS
# from pygeodesy.props import _update_all  # _MODS
# from pygeodesy.streprs import Fmt  # _MODS
# from pygeodesy.unitsBase import _NamedUnit, Str  # _MODS

from copy import copy as _copy, deepcopy as _deepcopy
from math import copysign as _copysign
import inspect as _inspect

__all__ = _ALL_LAZY.basics
__version__ = '24.06.15'

_below_               = 'below'
_list_tuple_types     = (list, tuple)
_PYGEODESY_XPACKAGES_ = 'PYGEODESY_XPACKAGES'
_required_            = 'required'

try:  # Luciano Ramalho, "Fluent Python", O'Reilly, 2016 p. 395, 2022 p. 577+
    from numbers import Integral as _Ints, Real as _Scalars  # .units
except ImportError:
    try:
        _Ints = int, long  # int objects (C{tuple})
    except NameError:  # Python 3+
        _Ints = int,  # int objects (C{tuple})
    _Scalars = (float,) + _Ints

try:
    try:  # use C{from collections.abc import ...} in Python 3.9+
        from collections.abc import Sequence as _Sequence  # in .points
    except ImportError:  # no .abc in Python 3.8- and 2.7-
        from collections import Sequence as _Sequence  # in .points
    if isinstance([], _Sequence) and isinstance((), _Sequence):
        #                        and isinstance(range(1), _Sequence):
        _Seqs = _Sequence
    else:
        raise ImportError()  # _AssertionError
except ImportError:
    _Sequence = tuple  # immutable for .points._Basequence
    _Seqs     = list, _Sequence  # range for function len2 below

try:
    _Bytes = unicode, bytearray  # PYCHOK expected
    _Strs  = basestring, str  # XXX , bytes
    str2ub = ub2str = _passarg  # avoids UnicodeDecodeError

    def _Xstr(exc):  # PYCHOK no cover
        '''I{Invoke only with caught ImportError} B{C{exc}}.

           C{... "can't import name _distributor_init" ...}

           only for C{numpy}, C{scipy} import errors occurring
           on arm64 Apple Silicon running macOS' Python 2.7.16?
        '''
        t = str(exc)
        if '_distributor_init' in t:
            from sys import exc_info
            from traceback import extract_tb
            tb = exc_info()[2]  # 3-tuple (type, value, traceback)
            t4 = extract_tb(tb, 1)[0]  # 4-tuple (file, line, name, 'import ...')
            t = _SPACE_("can't", t4[3] or _N_A_)
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


def _args_kwds_count2(func, exelf=True):
    '''(INTERNAL) Get a C{func}'s args and kwds count as 2-tuple
       C{(nargs, nkwds)}, including arg C{self} for methods.

       @kwarg exelf: If C{True}, exclude C{self} in the C{args}
                     of a method (C{bool}).
    '''
    try:
        a = k = 0
        for _, p in _inspect.signature(func).parameters.items():
            if p.kind is p.POSITIONAL_OR_KEYWORD:
                if p.default is p.empty:
                    a += 1
                else:
                    k += 1
    except AttributeError:  # .signature new Python 3+
        s = _inspect.getargspec(func)
        k = len(s.defaults or ())
        a = len(s.args) - k
    if exelf and a > 0 and _inspect.ismethod(func):
        a -= 1
    return a, k


def _args_kwds_names(func, splast=False):
    '''(INTERNAL) Get a C{func}'s args and kwds names, including
       C{self} for methods.

       @kwarg splast: If C{True}, split the last keyword argument
                      at UNDERscores (C{bool}).

       @note: Python 2 may I{not} include the C{*args} nor the
              C{**kwds} names.
    '''
    try:
        args_kwds = _inspect.signature(func).parameters.keys()
    except AttributeError:  # .signature new Python 3+
        args_kwds = _inspect.getargspec(func).args
    if splast and args_kwds:
        args_kwds = list(args_kwds)
        t = args_kwds[-1:]
        if t:
            s = t[0].strip(_UNDER_).split(_UNDER_)
            if len(s) > 1 or s != t:
                args_kwds += s
    return tuple(args_kwds)


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


def _enumereverse(iterable):
    '''(INTERNAL) Reversed C{enumberate}.
    '''
    for j in _reverange(len(iterable)):
        yield j, iterable[j]


def halfs2(str2):
    '''Split a string in 2 halfs.

       @arg str2: String to split (C{str}).

       @return: 2-Tuple C{(_1st, _2nd)} half (C{str}).

       @raise ValueError: Zero or odd C{len(B{str2})}.
    '''
    h, r = divmod(len(str2), 2)
    if r or not h:
        raise _ValueError(str2=str2, txt=_odd_)
    return str2[:h], str2[h:]


def int1s(x):
    '''Count the number of 1-bits in an C{int}, I{unsigned}.

       @note: C{int1s(-B{x}) == int1s(abs(B{x}))}.
    '''
    try:
        return x.bit_count()  # Python 3.10+
    except AttributeError:
        # bin(-x) = '-' + bin(abs(x))
        return bin(x).count(_1_)


def isbool(obj):
    '''Is B{C{obj}}ect a C{bool}ean?

       @arg obj: The object (any C{type}).

       @return: C{True} if C{bool}ean, C{False} otherwise.
    '''
    return isinstance(obj, bool)  # and (obj is False
#                                     or obj is True)

assert not (isbool(1) or isbool(0) or isbool(None))  # PYCHOK 2


def isCartesian(obj, ellipsoidal=None):
    '''Is B{C{obj}}ect some C{Cartesian}?

       @arg obj: The object (any C{type}).
       @kwarg ellipsoidal: If C{None}, return the type of any C{Cartesian},
                           if C{True}, only an ellipsoidal C{Cartesian type}
                           or if C{False}, only a spherical C{Cartesian type}.

       @return: C{type(B{obj}} if a C{Cartesian} of the required type, C{False}
                if a C{Cartesian} of an other type or {None} otherwise.
    '''
    if ellipsoidal is not None:
        try:
            return obj.ellipsoidalCartesian if ellipsoidal else obj.sphericalCartesian
        except AttributeError:
            return None
    return isinstanceof(obj, _MODS.cartesianBase.CartesianBase)


if _FOR_DOCS:  # XXX avoid epydoc Python 2.7 error

    def isclass(obj):
        '''Is B{C{obj}}ect a C{Class} or C{type}?
        '''
        return _inspect.isclass(obj)
else:
    isclass = _inspect.isclass


def iscomplex(obj, both=False):
    '''Is B{C{obj}}ect a C{complex} or complex literal C{str}?

       @arg obj: The object (any C{type}).
       @kwarg both: If C{True}, check complex C{str} (C{bool}).

       @return: C{True} if C{complex}, C{False} otherwise.
    '''
    try:  # hasattr('conjugate', 'real' and 'imag')
        return isinstance(obj, complex) or bool(both and isstr(obj) and
               isinstance(complex(obj), complex))  # numbers.Complex?
    except (TypeError, ValueError):
        return False


def isDEPRECATED(obj):
    '''Is B{C{obj}}ect a C{DEPRECATED} class, method or function?

       @return: C{True} if C{DEPRECATED}, {False} if not or
                C{None} if undetermined.
    '''
    try:  # XXX inspect.getdoc(obj) or obj.__doc__
        doc = obj.__doc__.lstrip()
        return bool(doc and doc.startswith(_DEPRECATED_))
    except AttributeError:
        return None


def isfloat(obj, both=False):
    '''Is B{C{obj}}ect a C{float} or float literal C{str}?

       @arg obj: The object (any C{type}).
       @kwarg both: If C{True}, check float C{str} (C{bool}).

       @return: C{True} if C{float}, C{False} otherwise.
    '''
    try:
        return isinstance(obj, float) or bool(both and
               isstr(obj) and isinstance(float(obj), float))
    except (TypeError, ValueError):
        return False


try:
    isidentifier = str.isidentifier  # Python 3, must be str
except AttributeError:  # Python 2-

    def isidentifier(obj):
        '''Is B{C{obj}}ect a Python identifier?
        '''
        return bool(obj and isstr(obj)
                        and obj.replace(_UNDER_, NN).isalnum()
                        and not obj[:1].isdigit())


def isinstanceof(obj, *Classes):
    '''Is B{C{obj}}ect an instance of one of the C{Classes}?

       @arg obj: The object (any C{type}).
       @arg Classes: One or more classes (C{Class}).

       @return: C{type(B{obj}} if one of the B{C{Classes}},
                C{None} otherwise.
    '''
    return type(obj) if isinstance(obj, Classes) else None


def isint(obj, both=False):
    '''Is B{C{obj}}ect an C{int} or integer C{float} value?

       @arg obj: The object (any C{type}).
       @kwarg both: If C{True}, check C{float} and L{Fsum}
                    type and value (C{bool}).

       @return: C{True} if C{int} or I{integer} C{float}
                or L{Fsum}, C{False} otherwise.

       @note: Both C{isint(True)} and C{isint(False)} return
              C{False} (and no longer C{True}).
    '''
    if isinstance(obj, _Ints):
        return not isbool(obj)
    elif both:  # and isinstance(obj, (float, Fsum))
        try:  # NOT , _Scalars) to include Fsum!
            return obj.is_integer()
        except AttributeError:
            pass  # XXX float(int(obj)) == obj?
    return False


def isiterable(obj):
    '''Is B{C{obj}}ect C{iterable}?

       @arg obj: The object (any C{type}).

       @return: C{True} if C{iterable}, C{False} otherwise.
    '''
    # <https://PyPI.org/project/isiterable/>
    return hasattr(obj, '__iter__')  # map, range, set


def isiterablen(obj):
    '''Is B{C{obj}}ect C{iterable} and has C{len}gth?

       @arg obj: The object (any C{type}).

       @return: C{True} if C{iterable} with C{len}gth, C{False} otherwise.
    '''
    return hasattr(obj, '__len__') and hasattr(obj, '__getitem__')


try:
    from keyword import iskeyword  # Python 2.7+
except ImportError:

    def iskeyword(unused):
        '''Not Implemented, C{False} always.
        '''
        return False


def isLatLon(obj, ellipsoidal=None):
    '''Is B{C{obj}}ect some C{LatLon}?

       @arg obj: The object (any C{type}).
       @kwarg ellipsoidal: If C{None}, return the type of any C{LatLon},
                           if C{True}, only an ellipsoidal C{LatLon type}
                           or if C{False}, only a spherical C{LatLon type}.

       @return: C{type(B{obj}} if a C{LatLon} of the required type, C{False}
                if a C{LatLon} of an other type or {None} otherwise.
    '''
    if ellipsoidal is not None:
        try:
            return obj.ellipsoidalLatLon if ellipsoidal else obj.sphericalLatLon
        except AttributeError:
            return None
    return isinstanceof(obj, _MODS.latlonBase.LatLonBase)


def islistuple(obj, minum=0):
    '''Is B{C{obj}}ect a C{list} or C{tuple} with non-zero length?

       @arg obj: The object (any C{type}).
       @kwarg minum: Minimal C{len} required C({int}).

       @return: C{True} if a C{list} or C{tuple} with C{len} at
                least B{C{minum}}, C{False} otherwise.
    '''
    return isinstance(obj, _list_tuple_types) and len(obj) >= minum


def isNvector(obj, ellipsoidal=None):
    '''Is B{C{obj}}ect some C{Nvector}?

       @arg obj: The object (any C{type}).
       @kwarg ellipsoidal: If C{None}, return the type of any C{Nvector},
                           if C{True}, only an ellipsoidal C{Nvector type}
                           or if C{False}, only a spherical C{Nvector type}.

       @return: C{type(B{obj}} if an C{Nvector} of the required type, C{False}
                if an C{Nvector} of an other type or {None} otherwise.
    '''
    if ellipsoidal is not None:
        try:
            return obj.ellipsoidalNvector if ellipsoidal else obj.sphericalNvector
        except AttributeError:
            return None
    return isinstanceof(obj, _MODS.nvectorBase.NvectorBase)


def isodd(x):
    '''Is B{C{x}} odd?

       @arg x: Value (C{scalar}).

       @return: C{True} if odd, C{False} otherwise.
    '''
    return bool(int(x) & 1)  # == bool(int(x) % 2)


def isscalar(obj, both=False):
    '''Is B{C{obj}}ect an C{int} or integer C{float} value?

       @arg obj: The object (any C{type}).
       @kwarg both: If C{True}, check L{Fsum<Fsum.residual>}.

       @return: C{True} if C{int}, C{float} or L{Fsum} with
                zero residual, C{False} otherwise.
    '''
    if isinstance(obj, _Scalars):
        return not isbool(obj)
    elif both:  # and isinstance(obj, Fsum)
        try:
            return bool(obj.residual == 0)
        except (AttributeError, TypeError):
            pass  # XXX float(int(obj)) == obj?
    return False


def issequence(obj, *excls):
    '''Is B{C{obj}}ect some sequence type?

       @arg obj: The object (any C{type}).
       @arg excls: Classes to exclude (C{type}), all positional.

       @note: Excluding C{tuple} implies excluding C{namedtuple}.

       @return: C{True} if a sequence, C{False} otherwise.
    '''
    return isinstance(obj, _Seqs) and not (excls and isinstance(obj, excls))


def isstr(obj):
    '''Is B{C{obj}}ect some string type?

       @arg obj: The object (any C{type}).

       @return: C{True} if a C{str}, C{bytes}, ...,
                C{False} otherwise.
    '''
    return isinstance(obj, _Strs)


def issubclassof(Sub, *Supers):
    '''Is B{C{Sub}} a class and sub-class of some other class(es)?

       @arg Sub: The sub-class (C{Class}).
       @arg Supers: One or more C(super) classes (C{Class}).

       @return: C{True} if a sub-class of any B{C{Supers}}, C{False}
                if not (C{bool}) or C{None} if not a class or if no
                B{C{Supers}} are given or none of those are a class.
    '''
    if isclass(Sub):
        t = tuple(S for S in Supers if isclass(S))
        if t:
            return bool(issubclass(Sub, t))
    return None


def itemsorted(adict, *items_args, **asorted_reverse):
    '''Return the items of C{B{adict}} sorted I{alphabetically,
       case-insensitively} and in I{ascending} order.

       @arg items_args: Optional positional argument(s) for method
                        C{B{adict}.items(B*{items_args})}.
       @kwarg asorted_reverse: Use C{B{asorted}=False} for I{alphabetical,
                      case-sensitive} sorting and C{B{reverse}=True} for
                      sorting in C{descending} order.
    '''
    def _ins(item):  # functools.cmp_to_key
        k, v = item
        return k.lower()

    def _reverse_key(asorted=True, reverse=False):
        return dict(reverse=reverse, key=_ins if asorted else None)

    items = adict.items(*items_args) if items_args else adict.items()
    return sorted(items, **_reverse_key(**asorted_reverse))


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
    '''Call a single-argument function to each B{C{xs}}
       and return a C{tuple} of results.

       @arg fun1: 1-Arg function (C{callable}).
       @arg xs: Arguments (C{any positional}).

       @return: Function results (C{tuple}).
    '''
    return tuple(map(fun1, xs))


def map2(fun, *xs):
    '''Like Python's B{C{map}} but returning a C{tuple} of results.

       Unlike Python 2's built-in L{map}, Python 3+ L{map} returns a
       L{map} object, an iterator-like object which generates the
       results only once.  Converting the L{map} object to a tuple
       maintains the Python 2 behavior.

       @arg fun: Function (C{callable}).
       @arg xs: Arguments (C{all positional}).

       @return: Function results (C{tuple}).
    '''
    return tuple(map(fun, *xs))


def neg(x, neg0=None):
    '''Negate C{x} and optionally, negate C{0.0} and C{-0.0}.

       @kwarg neg0: Defines the return value for zero C{B{x}}: if C{None}
                    return C{0.0}, if C{True} return C{NEG0 if B{x}=0.0}
                    and C{0.0 if B{x}=NEG0} or if C{False} return C{B{x}}
                    I{as-is} (C{bool} or C{None}).

       @return: C{-B{x} if B{x} else 0.0, NEG0 or B{x}}.
    '''
    return (-x) if x else (
           _0_0 if neg0 is None else (
             x  if not neg0 else (
           _0_0 if signBit(x) else _MODS.constants.
           NEG0)))  # PYCHOK indent


def neg_(*xs):
    '''Negate all C{xs} with L{neg}.

       @return: A C{map(neg, B{xs})}.
    '''
    return map(neg, xs)


def _neg0(x):
    '''(INTERNAL) Return C{NEG0 if x < 0 else _0_0},
       unlike C{_copysign_0_0} which returns C{_N_0_0}.
    '''
    return _MODS.constants.NEG0 if x < 0 else _0_0


def _req_d_by(where, **name):
    '''(INTERNAL) Get the fully qualified name.
    '''
    m = _MODS.named
    n =  m._name__(**name)
    m =  m.modulename(where, prefixed=True)
    if n:
        m = _DOT_(m, n)
    return _SPACE_(_required_, _by_, m)


def _reverange(n, stop=-1, step=-1):
    '''(INTERNAL) Reversed range yielding C{n-1, n-1-step, ..., stop+1}.
    '''
    return range(n - 1, stop, step)


def signBit(x):
    '''Return C{signbit(B{x})}, like C++.

       @return: C{True} if C{B{x} < 0} or C{NEG0} (C{bool}).
    '''
    return x < 0 or _MODS.constants.isneg0(x)


def _signOf(x, ref):  # in .fsums
    '''(INTERNAL) Return the sign of B{C{x}} versus B{C{ref}}.
    '''
    return (-1) if x < ref else (+1 if x > ref else 0)


def signOf(x):
    '''Return sign of C{x} as C{int}.

       @return: -1, 0 or +1 (C{int}).
    '''
    try:
        s = x.signOf()  # Fsum instance?
    except AttributeError:
        s =  _signOf(x, 0)
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

    t = _xiterablen(iterable)
    if not isinstance(t, _list_tuple_types):
        t = tuple(t)

    if n > 1:
        if fill:
            fill = _xkwds_get1(fill, fill=MISSING)
            if fill is not MISSING:
                m = len(t) % n
                if m > 0:  # same type fill
                    t = t + type(t)((fill,) * (n - m))
        for i in range(n):
            # XXX t[i::n] chokes PyChecker
            yield t[slice(i, None, n)]
    else:
        yield t  # 1 slice, all


def _splituple(strs, *sep_splits):  # in .mgrs, .osgr, .webmercator
    '''(INTERNAL) Split a C{comma}- or C{whitespace}-separated
       string into a C{tuple} of stripped strings.
    '''
    t = (strs.split(*sep_splits) if sep_splits else
         strs.replace(_COMMA_, _SPACE_).split()) if strs else ()
    return tuple(s.strip() for s in t if s)


def unsigned0(x):
    '''Unsign if C{0.0}.

       @return: C{B{x}} if B{C{x}} else C{0.0}.
    '''
    return x if x else _0_0


def _xcopy(obj, deep=False):
    '''(INTERNAL) Copy an object, shallow or deep.

       @arg obj: The object to copy (any C{type}).
       @kwarg deep: If C{True} make a deep, otherwise
                    a shallow copy (C{bool}).

       @return: The copy of B{C{obj}}.
    '''
    return _deepcopy(obj) if deep else _copy(obj)


def _xdup(obj, deep=False, **items):
    '''(INTERNAL) Duplicate an object, replacing some attributes.

       @arg obj: The object to copy (any C{type}).
       @kwarg deep: If C{True} copy deep, otherwise shallow.
       @kwarg items: Attributes to be changed (C{any}).

       @return: A duplicate of B{C{obj}} with modified
                attributes, if any B{C{items}}.

       @raise AttributeError: Some B{C{items}} invalid.
    '''
    d = _xcopy(obj, deep=deep)
    for n, v in items.items():
        if getattr(d, n, v) != v:
            setattr(d, n, v)
        elif not hasattr(d, n):
            t = _MODS.named.classname(obj)
            t = _SPACE_(_DOT_(t, n), _invalid_)
            raise _AttributeError(txt=t, obj=obj, **items)
#   if items:
#       _MODS.props._update_all(d)
    return d


def _xgeographiclib(where, *required):
    '''(INTERNAL) Import C{geographiclib} and check required version.
    '''
    try:
        _xpackage(_xgeographiclib)
        import geographiclib
    except ImportError as x:
        raise _xImportError(x, where, Error=LazyImportError)
    return _xversion(geographiclib, where, *required)


def _xImportError(exc, where, Error=_ImportError, **name):
    '''(INTERNAL) Embellish an C{Lazy/ImportError}.
    '''
    t = _req_d_by(where, **name)
    return Error(_Xstr(exc), txt=t, cause=exc)


def _xinstanceof(*Types, **names_values):
    '''(INTERNAL) Check C{Types} of all C{name=value} pairs.

       @arg Types: One or more classes or types (C{class}), all
                   positional.
       @kwarg names_values: One or more C{B{name}=value} pairs
                            with the C{value} to be checked.

       @raise TypeError: One B{C{names_values}} pair is not an
                         instance of any of the B{C{Types}}.
    '''
    if not (Types and names_values):
        raise _xAssertionError(_xinstanceof, *Types, **names_values)

    for n, v in names_values.items():
        if not isinstance(v, Types):
            raise _TypesError(n, v, *Types)


def _xiterable(obj):
    '''(INTERNAL) Return C{obj} if iterable, otherwise raise C{TypeError}.
    '''
    return obj if isiterable(obj) else _xiterror(obj, _xiterable)  # PYCHOK None


def _xiterablen(obj):
    '''(INTERNAL) Return C{obj} if iterable with C{__len__}, otherwise raise C{TypeError}.
    '''
    return obj if isiterablen(obj) else _xiterror(obj, _xiterablen)  # PYCHOK None


def _xiterror(obj, _xwhich):
    '''(INTERNAL) Helper for C{_xinterable} and C{_xiterablen}.
    '''
    t = _not_(_xwhich.__name__[2:])  # _dunder_nameof
    raise _TypeError(repr(obj), txt=t)


def _xnumpy(where, *required):
    '''(INTERNAL) Import C{numpy} and check required version.
    '''
    try:
        _xpackage(_xnumpy)
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


def _xpackage(_xpkg):
    '''(INTERNAL) Check dependency to be excluded.
    '''
    n = _xpkg.__name__[2:]  # _dunder_nameof
    if n in _XPACKAGES:
        x = _SPACE_(n, _in_, _PYGEODESY_XPACKAGES_)
        e = _enquote(_getenv(_PYGEODESY_XPACKAGES_, NN))
        raise ImportError(_EQUAL_(x, e))


def _xscalar(**names_values):
    '''(INTERNAL) Check all C{name=value} pairs to be C{scalar}.
    '''
    for n, v in names_values.items():
        if not isscalar(v):
            raise _TypeError(n, v, txt=_not_scalar_)


def _xscipy(where, *required):
    '''(INTERNAL) Import C{scipy} and check required version.
    '''
    try:
        _xpackage(_xscipy)
        import scipy
    except ImportError as x:
        raise _xImportError(x, where)
    return _xversion(scipy, where, *required)


def _xsubclassof(*Classes, **names_values):
    '''(INTERNAL) Check (super) class of all C{name=value} pairs.

       @arg Classes: One or more classes or types (C{class}), all
                     positional.
       @kwarg names_values: One or more C{B{name}=value} pairs
                            with the C{value} to be checked.

       @raise TypeError: One B{C{names_values}} pair is not a
                         (sub-)class of any of the B{C{Classes}}.
    '''
    if not (Classes and names_values):
        raise _xAssertionError(_xsubclassof, *Classes, **names_values)

    for n, v in names_values.items():
        if not issubclassof(v, *Classes):
            raise _TypesError(n, v, *Classes)


def _xversion(package, where, *required, **name):
    '''(INTERNAL) Check the C{package} version vs B{C{required}}.
    '''
    if required:
        t = _version_info(package)
        if t[:len(required)] < required:
            t = _SPACE_(package.__name__,  # _dunder_nameof
                       _version_, _DOT_(*t),
                       _below_, _DOT_(*required),
                       _req_d_by(where, **name))
            raise ImportError(t)
    return package


def _xzip(*args, **strict):  # PYCHOK no cover
    '''(INTERNAL) Standard C{zip(..., strict=True)}.
    '''
    s = _xkwds_get1(strict, strict=True)
    if s:
        if _zip is zip:  # < (3, 10)
            t = _MODS.streprs.unstr(_xzip, *args, strict=s)
            raise _NotImplementedError(t, txt=None)
        return _zip(*args)
    return zip(*args)


if _sys_version_info2 < (3, 10):  # see .errors
    _zip = zip  # PYCHOK exported
else:  # Python 3.10+

    def _zip(*args):
        return zip(*args, strict=True)

_XPACKAGES = _splituple(_getenv(_PYGEODESY_XPACKAGES_, NN).lower())

# **) MIT License
#
# Copyright (C) 2016-2024 -- mrJean1 at Gmail -- All Rights Reserved.
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
