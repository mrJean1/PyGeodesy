
# -*- coding: utf-8 -*-

u'''Errors, exceptions and exception chaining.

Error, exception classes and functions to format PyGeodesy errors,
including the setting of I{exception chaining} in Python 3+.

By default, I{exception chaining} is turned I{off}.  To enable
I{exception chaining}, use command line option C{python -X dev}
I{OR} set env var C{PYTHONDEVMODE} to C{1} or any non-empyty
string I{OR} set env var C{PYGEODESY_EXCEPTION_CHAINING=std}
or any other non-empty string.
'''
from pygeodesy.interns import MISSING, NN, _a_,_an_, _and_, \
                             _COLON_, _COMMA_, _COMMASPACE_, \
                             _datum_, _ellipsoidal_, _EQUAL_, \
                             _incompatible_, _invalid_, _len_, \
                             _name_, _no_, _not_, _or_, _SPACE_, \
                             _specified_, _UNDER_, _value_
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS, \
                             _getenv, _PYTHON_X_DEV

__all__ = _ALL_LAZY.errors  # _ALL_DOCS('_InvalidError', '_IsnotError')
__version__ = '22.04.09'

_default_    = 'default'
_kwargs_     = 'kwargs'
_limiterrors =  True  # imported by .formy
_multiple_   = 'multiple'
_name_value_ =  repr('name=value')
_rangerrors  =  True  # imported by .dms
__vs__       = ' vs '  # _vsSPACED_
_with_       = 'with'

try:
    _exception_chaining = None  # not available
    _ = Exception().__cause__   # Python 3+ exception chaining

    if _PYTHON_X_DEV or _getenv('PYGEODESY_EXCEPTION_CHAINING', NN):  # == _std_
        _exception_chaining = True  # turned on, std
        raise AttributeError  # allow exception chaining

    _exception_chaining = False  # turned off

    def _error_chain(inst, other=None):
        '''(INTERNAL) Set or avoid Python 3+ exception chaining.

           Setting C{inst.__cause__ = None} is equivalent to syntax
           C{raise Error(...) from None} to avoid exception chaining.

           @arg inst: An error instance (C{Exception}).
           @kwarg other: The previous error instance (C{Exception}) or
                         C{None} to avoid exception chaining.

           @see: Alex Martelli, et.al., "Python in a Nutshell", 3rd Ed., page 163,
                 O'Reilly, 2017, U{PEP-3134<https://www.Python.org/dev/peps/pep-3134>},
                 U{here <https://StackOverflow.com/questions/17091520/how-can-i-more-
                 easily-suppress-previous-exceptions-when-i-raise-my-own-exception>}
                 and U{here<https://StackOverflow.com/questions/1350671/
                 inner-exception-with-traceback-in-python>}.
        '''
        inst.__cause__ = other  # None, no exception chaining
        return inst

except AttributeError:  # Python 2+

    def _error_chain(inst, **unused):  # PYCHOK expected
        return inst  # no-op


class _AssertionError(AssertionError):
    '''(INTERNAL) Format an C{AssertionError} without exception chaining.
    '''
    def __init__(self, *name_value, **txt_name_values):  # txt=_invalid_
        _error_init(AssertionError, self, name_value, **txt_name_values)


class _AttributeError(AttributeError):
    '''(INTERNAL) Format an C{AttributeError} without exception chaining.
    '''
    def __init__(self, *name_value, **txt_name_values):  # txt=_invalid_
        _error_init(AttributeError, self, name_value, **txt_name_values)


class _ImportError(ImportError):
    '''(INTERNAL) Format an C{ImportError} without exception chaining.
    '''
    def __init__(self, *name_value, **txt_name_values):  # txt=_invalid_
        _error_init(ImportError, self, name_value, **txt_name_values)


class _IndexError(IndexError):
    '''(INTERNAL) Format an C{IndexError} without exception chaining.
    '''
    def __init__(self, *name_value, **txt_name_values):  # txt=_invalid_
        _error_init(IndexError, self, name_value, **txt_name_values)


class _NameError(NameError):
    '''(INTERNAL) Format a C{NameError} without exception chaining.
    '''
    def __init__(self, *name_value, **txt_name_values):  # txt=_invalid_
        _error_init(NameError, self, name_value, **txt_name_values)


class _NotImplementedError(NotImplementedError):
    '''(INTERNAL) Format a C{NotImplementedError} without exception chaining.
    '''
    def __init__(self, *name_value, **txt_name_values):  # txt=_invalid_
        _error_init(NotImplementedError, self, name_value, **txt_name_values)


class _OverflowError(OverflowError):
    '''(INTERNAL) Format an C{OverflowError} without exception chaining.
    '''
    def __init__(self, *name_value, **txt_name_values):  # txt=_invalid_
        _error_init(OverflowError, self, name_value, **txt_name_values)


class _TypeError(TypeError):
    '''(INTERNAL) Format a C{TypeError} without exception chaining.
    '''
    def __init__(self, *name_value, **txt_name_values):
        _error_init(TypeError, self, name_value, fmt_name_value='type(%s) (%r)',
                                               **txt_name_values)


class _TypesError(_TypeError):
    '''(INTERNAL) Format a C{TypeError} without exception chaining.
    '''
    def __init__(self, name, value, *Types):
        t = _not_(_an(_or(*(t.__name__ for t in Types))))
        _TypeError.__init__(self, name, value, txt=t)


class _ValueError(ValueError):
    '''(INTERNAL) Format a C{ValueError} without exception chaining.
    '''
    def __init__(self, *name_value, **txt_name_values):  # name, value, txt=_invalid_
        _error_init(ValueError, self, name_value, **txt_name_values)


class _ZeroDivisionError(ZeroDivisionError):
    '''(INTERNAL) Format a C{ZeroDivisionError} without exception chaining.
    '''
    def __init__(self, *name_value, **txt_name_values):  # name, value, txt=_invalid_
        _error_init(ZeroDivisionError, self, name_value, **txt_name_values)


class CrossError(_ValueError):
    '''Error raised for zero or near-zero vectorial cross products,
       occurring for coincident or colinear points, paths or bearings.
    '''
    pass


class IntersectionError(_ValueError):  # in .ellipsoidalBaseDI, .formy, ...
    '''Error raised for path or circle intersection issues.
    '''
    def __init__(self, *args, **kwds):  # txt=_invalid_
        '''New L{IntersectionError}.
        '''
        if args:
            _ValueError.__init__(self, _SPACE_(*args), **kwds)
        else:
            _ValueError.__init__(self, **kwds)


class LenError(_ValueError):  # in .ecef, .fmath, .heights, .iters, .named
    '''Error raised for mis-matching C{len} values.
    '''
    def __init__(self, where, **lens_txt):  # txt=None
        '''New L{LenError}.

           @arg where: Object with C{.__name__} attribute
                       (C{class}, C{method}, or C{function}).
           @kwarg lens_txt: Two or more C{name=len(name)} pairs
                            (C{keyword arguments}).
        '''
        PAREN = _MODS.streprs.Fmt.PAREN
        x = _xkwds_pop(lens_txt, txt=_invalid_)
        ns, vs = zip(*sorted(lens_txt.items()))
        ns = _COMMASPACE_.join(ns)
        vs = __vs__.join(map(str, vs))
        t  = _SPACE_(PAREN(where.__name__, ns), _len_, vs)
        _ValueError.__init__(self, t, txt=x)


class LimitError(_ValueError):
    '''Error raised for lat- or longitudinal deltas exceeding
       the B{C{limit}} in functions L{pygeodesy.equirectangular} and
       L{pygeodesy.equirectangular_} and several C{nearestOn*} and
       C{simplify*} functions or methods.
    '''
    pass


class NumPyError(_ValueError):
    '''Error raised for C{NumPy} errors.
    '''
    pass


class ParseError(_ValueError):  # in .dms, .elevations, .utmupsBase
    '''Error parsing degrees, radians or several other formats.
    '''
    pass


class PointsError(_ValueError):  # in .clipy, .frechet, ...
    '''Error for an insufficient number of points.
    '''
    pass


class RangeError(_ValueError):  # in .dms, .ellipsoidalBase, .geoids, .units, .ups, .utm, .utmups
    '''Error raised for lat- or longitude values outside the B{C{clip}}, B{C{clipLat}},
       B{C{clipLon}} or B{C{limit}} range in function L{pygeodesy.clipDegrees},
       L{pygeodesy.clipRadians}, L{pygeodesy.parse3llh}, L{pygeodesy.parseDMS},
       L{pygeodesy.parseDMS2} or L{pygeodesy.parseRad}.

       @see: Function L{pygeodesy.rangerrors}.
    '''
    pass


class TriangleError(_ValueError):  # in .resections, .vector2d
    '''Error raised for triangle, inter- or resection issues.
    '''
    pass


class SciPyError(PointsError):
    '''Error raised for C{SciPy} errors.
    '''
    pass


class SciPyWarning(PointsError):
    '''Error thrown for C{SciPy} warnings.

       To raise C{SciPy} warnings as L{SciPyWarning} exceptions, Python
       C{warnings} must be filtered as U{warnings.filterwarnings('error')
       <https://docs.Python.org/3/library/warnings.html#the-warnings-filter>}
       I{prior to} C{import scipy} OR by setting env var U{PYTHONWARNINGS
       <https://docs.Python.org/3/using/cmdline.html#envvar-PYTHONWARNINGS>}
       OR by invoking C{python} with command line option U{-W<https://docs.
       Python.org/3/using/cmdline.html#cmdoption-w>} set to C{-W error}.
    '''
    pass


class TRFError(_ValueError):  # in .ellipsoidalBase, .trf, .units
    '''Terrestrial Reference Frame (TRF), L{Epoch}, L{RefFrame}
       or L{RefFrame} conversion issue.
    '''
    pass


class UnitError(_ValueError):  # in .named, .units
    '''Default exception for L{units} issues.
    '''
    pass


class VectorError(_ValueError):  # in .nvectorBase, .vector3d, .vector3dBase
    '''L{Vector3d}, C{Cartesian*} or C{*Nvector} issues.
    '''
    pass


def _an(noun):
    '''(INTERNAL) Prepend an article to a noun based
       on the pronounciation of the first letter.
    '''
    a = _an_ if noun[:1].lower() in 'aeinoux' else _a_
    return _SPACE_(a, noun)


def _and(*words):
    '''(INTERNAL) Join C{words} with C{", "} and C{" and "}.
    '''
    return _and_or(_and_, *words)


def _and_or(last, *words):
    '''(INTERNAL) Join C{words} with C{", "} and C{B{last}}.
    '''
    t, w = NN, list(words)
    if w:
        t = w.pop()
        if w:
            w = _COMMASPACE_.join(w)
            t = _SPACE_(w, last, t)
    return t


def crosserrors(raiser=None):
    '''Report or ignore vectorial cross product errors.

       @kwarg raiser: Use C{True} to throw or C{False} to ignore
                      L{CrossError} exceptions.  Use C{None} to
                      leave the setting unchanged.

       @return: Previous setting (C{bool}).

       @see: Property C{Vector3d[Base].crosserrors}.
    '''
    B = _MODS.vector3dBase.Vector3dBase
    t =  B._crosserrors  # XXX class attr!
    if raiser in (True, False):
        B._crosserrors = raiser
    return t


def _error_init(Error, inst, name_value, fmt_name_value='%s (%r)',
                                         txt=_invalid_, **name_values):  # by .lazily
    '''(INTERNAL) Format an error text and initialize an C{Error} instance.

       @arg Error: The error super-class (C{Exception}).
       @arg inst: Sub-class instance to be initialized (C{_Exception}).
       @arg name_value: Either just a value or several name, value, ...
                        positional arguments (C{str}, any C{type}), in
                        particular for name conflicts with keyword
                        arguments of C{error_init} or which can't be
                        used as C{name=value} keyword arguments.
       @kwarg fmt_name_value: Format for (name, value) (C{str}).
       @kwarg txt: Optional explanation of the error (C{str}).
       @kwarg name_values: One or more C{B{name}=value} pairs overriding
                           any B{C{name_value}} positional arguments.
    '''
    if name_values:
        t = _or(*sorted(fmt_name_value % t for t in name_values.items()))
    elif len(name_value) > 1:
        t = _or(*sorted(fmt_name_value % t for t in zip(name_value[0::2],
                                                        name_value[1::2])))
    elif name_value:
        t = str(name_value[0])
    else:
        t = _SPACE_(_name_value_, str(MISSING))

    if txt is not None:
        c = _COMMA_ if _COLON_ in t else _COLON_
        x =  str(txt) or _invalid_
        t =  NN(t, c, _SPACE_, x)
#   else:
#       x = NN  # XXX or t?
    Error.__init__(inst, t)
#   inst.__x_txt__ = x  # hold explanation
    _error_chain(inst)  # no Python 3+ exception chaining
    _error_under(inst)


def _error_under(inst):
    '''(INTERNAL) Remove leading underscore from instance' class name.
    '''
    n = inst.__class__.__name__
    if n.startswith(_UNDER_):
        inst.__class__.__name__ = n.lstrip(_UNDER_)
    return inst


def exception_chaining(error=None):
    '''Get the previous exception's or exception chaining setting.

       @kwarg error: An error instance (C{Exception}) or C{None}.

       @return: If C{B{error} is None}, return C{True} if exception
                chaining is enabled for PyGeodesy errors, C{False}
                if turned off and C{None} if not available.  If
                B{C{error}} is not C{None}, return the previous,
                chained error or C{None} otherwise.

       @note: Set env var C{PYGEODESY_EXCEPTION_CHAINING} to any
              non-empty value prior to C{import pygeodesy} to
              enable exception chaining for C{pygeodesy} errors.
    '''
    return _exception_chaining if error is None else \
            getattr(error, '__cause__', None)


def _incompatible(this):
    '''(INTERNAL) Format an incompatible text.
    '''
    return _SPACE_(_incompatible_, _with_, this)


def _InvalidError(Error=_ValueError, **txt_name_values):  # txt=_invalid_, name=value [, ...]
    '''(INTERNAL) Create an C{Error} instance.

       @kwarg Error: The error class or sub-class (C{Exception}).
       @kwarg txt_name_values: One or more C{B{name}=value} pairs
                               and optionally, a C{B{txt}=...}
                               keyword argument to override the
                               default C{B{txt}='invalid'}

       @return: An B{C{Error}} instance.
    '''
    try:
        e = Error(**txt_name_values)
    except TypeError:  # std *Error takes no keyword arguments
        e = _ValueError(**txt_name_values)
        e = Error(str(e))
        _error_chain(e)
        _error_under(e)
    return e


def _IsnotError(*nouns, **name_value_Error):  # name=value [, Error=TypeError]
    '''Create a C{TypeError} for an invalid C{name=value} type.

       @arg nouns: One or more expected class or type names,
                   usually nouns (C{str}).
       @kwarg name_value_Error: One C{B{name}=value} pair and
                                optionally, an C{B{Error}=...}
                                keyword argument to override
                                the default C{B{Error}=TypeError}.

       @return: A C{TypeError} or an B{C{Error}} instance.
    '''
    Error = _xkwds_pop(name_value_Error, Error=TypeError)
    n, v  = _xkwds_popitem(name_value_Error) if name_value_Error else (
                          _name_value_, MISSING)  # XXX else tuple(...)
    v = _MODS.streprs.Fmt.PAREN(repr(v))
    t = _or(*nouns) or _specified_
    if len(nouns) > 1:
        t = _an(t)
    e = Error(_SPACE_(n, v, _not_, t))
    _error_chain(e)
    _error_under(e)
    return e


def limiterrors(raiser=None):
    '''Get/set the throwing of L{LimitError}s.

       @kwarg raiser: Choose C{True} to raise or C{False} to
                      ignore L{LimitError} exceptions.  Use
                      C{None} to leave the setting unchanged.

       @return: Previous setting (C{bool}).
    '''
    global _limiterrors
    t = _limiterrors
    if raiser in (True, False):
        _limiterrors = raiser
    return t


def _or(*words):
    '''(INTERNAL) Join C{words} with C{", "} and C{" or "}.
    '''
    return _and_or(_or_, *words)


def _parseX(parser, *args, **name_values_Error):  # name=value[, ..., Error=ParseError]
    '''(INTERNAL) Invoke a parser and handle exceptions.

       @arg parser: The parser (C{callable}).
       @arg args: Any parser positional arguments (any C{type}s).
       @kwarg name_values_Error: One or more C{B{name}=value} pairs
                                 and optionally, an C{B{Error}=...}
                                 keyword argument to override the
                                 default C{B{Error}=ParseError}.

       @return: Parser result.

       @raise ParseError: Or the specified C{B{Error}=...}.

       @raise RangeError: If that error occurred.
    '''
    try:
        return parser(*args)

    except RangeError as x:
        E, t = type(x), str(x)
    except (AttributeError, IndexError, TypeError, ValueError) as x:
        E, t = ParseError, str(x)
    raise _InvalidError(**_xkwds(name_values_Error, Error=E, txt=t))


def rangerrors(raiser=None):
    '''Get/set the throwing of L{RangeError}s.

       @kwarg raiser: Choose C{True} to raise or C{False} to ignore
                      L{RangeError} exceptions.  Use C{None} to leave
                      the setting unchanged.

       @return: Previous setting (C{bool}).
    '''
    global _rangerrors
    t = _rangerrors
    if raiser in (True, False):
        _rangerrors = raiser
    return t


def _SciPyIssue(x, *extras):  # PYCHOK no cover
    if isinstance(x, (RuntimeWarning, UserWarning)):
        X = SciPyWarning
    else:
        X = SciPyError  # PYCHOK not really
    t = _SPACE_(str(x).strip(), *extras)
    return _error_chain(X(t), other=x)


def _xdatum(datum1, datum2, Error=None):
    '''(INTERNAL) Check for datum, ellipsoid or rhumb mis-match.
    '''
    if Error:
        E1, E2 = datum1.ellipsoid, datum2.ellipsoid
        if E1 != E2:
            raise Error(E2.named2, txt=_incompatible(E1.named2))
    elif datum1 != datum2:
        t = _SPACE_(_datum_, repr(datum1.name), _not_, repr(datum2.name))
        raise _AssertionError(t)


def _xellipsoidal(**name_value):
    '''(INTERNAL) Check an I{ellipsoidal} item.

       @return: The B{C{value}} if ellipsoidal.

       @raise TypeError: Not ellipsoidal B{C{value}}.
    '''
    try:
        for n, v in name_value.items():
            if v.isEllipsoidal:
                return v
            break
        else:
            n = v = MISSING
    except AttributeError:
        pass
    raise _TypeError(n, v, txt=_not_(_ellipsoidal_))


# map certain C{Exception} classes to the C{_Error}
_X2Error = {AssertionError:      _AssertionError,
            AttributeError:      _AttributeError,
            ImportError:         _ImportError,
            IndexError:          _IndexError,
            NameError:           _NameError,
            NotImplementedError: _NotImplementedError,
            OverflowError:       _OverflowError,
            TypeError:           _TypeError,
            ValueError:          _ValueError,
            ZeroDivisionError:   _ZeroDivisionError}
_xErrors = _TypeError, _TypeError, _ValueError


def _xError(x, *args, **kwds):
    '''(INTERNAL) Embellish an exception.

       @arg x: The exception instance (usually, C{_Error}).
       @arg args: Embelishments (C{any}).
       @arg kwds: Embelishments (C{any}).
    '''
    X = x.__class__
    t = str(x)
    try:  # C{_Error} style
        return X(txt=t, *args, **kwds)
    except TypeError:  # no keyword arguments
        pass
    # not an C{_Error}, format as C{_Error}
    t = str(_ValueError(txt=t, *args, **kwds))
    return x if _exception_chaining else X(t)


def _xError2(x):  # in .fsums
    '''(INTERNAL) Map an exception to 2-tuple (C{_Error} class, error C{txt}).

       @arg x: The exception instance (usually, C{Exception}).
    '''
    X =  x.__class__
    E = _X2Error.get(X, X)
    if E is X and not isinstance(x, _xErrors):
        E = _NotImplementedError
        t =  repr(x)
    else:
        t =  str(x)
    return E, t


try:
    _ = {}.__or__  # {} | {}  # Python 3.9+

    def _xkwds(kwds, **dflts):
        '''(INTERNAL) Override C{dflts} with specified C{kwds}.
        '''
        return (dflts | kwds) if kwds else dflts

except AttributeError:

    from copy import copy as _copy

    def _xkwds(kwds, **dflts):  # PYCHOK expected
        '''(INTERNAL) Override C{dflts} with specified C{kwds}.
        '''
        d = dflts
        if kwds:
            d = _copy(d)
            d.update(kwds)
        return d


def _xkwds_Error(where, kwds, name_txt, txt=_default_):
    # Helper for _xkwds_get, _xkwds_pop and _xkwds_popitem below
    pairs = _MODS.streprs.pairs
    f = _COMMASPACE_.join(pairs(kwds) + pairs(name_txt))
    f = _MODS.streprs.Fmt.PAREN(where.__name__, f)
    t = _multiple_ if name_txt else _no_
    t = _SPACE_(t, _EQUAL_(_name_, txt), _kwargs_)
    return _AssertionError(f, txt=t)


def _xkwds_get(kwds, **name_default):
    '''(INTERNAL) Get a C{kwds} value by C{name}, or the C{default}.
    '''
    if len(name_default) == 1:
        for n, d in name_default.items():
            return kwds.get(n, d)

    raise _xkwds_Error(_xkwds_get, kwds, name_default)


def _xkwds_get_(kwds, **names_defaults):
    '''(INTERNAL) Yield each C{kwds} value by C{name} or its C{default}
       in I{alphabetically sorted} order.
    '''
    for n, d in sorted(names_defaults.items()):
        yield kwds.get(n, d)


def _xkwds_not(*args, **kwds):
    '''(INTERNAL) Return C{kwds} with a value not in C{args}.
    '''
    return dict((n, v) for n, v in kwds.items() if v not in args)


def _xkwds_pop(kwds, **name_default):
    '''(INTERNAL) Pop a C{kwds} value by C{name}, or the C{default}.
    '''
    if len(name_default) == 1:
        return kwds.pop(*name_default.popitem())

    raise _xkwds_Error(_xkwds_pop, kwds, name_default)


def _xkwds_popitem(name_value):
    '''(INTERNAL) Return exactly one C{(name, value)} item.
    '''
    if len(name_value) == 1:  # XXX TypeError
        return name_value.popitem()  # XXX AttributeError

    raise _xkwds_Error(_xkwds_popitem, (), name_value, txt=_value_)

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
