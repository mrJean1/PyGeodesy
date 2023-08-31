
# -*- coding: utf-8 -*-

u'''Errors, exceptions, exception formatting and exception chaining.

Error, exception classes and functions to format PyGeodesy errors,
including the setting of I{exception chaining} in Python 3+.

By default, I{exception chaining} is turned I{off}.  To enable I{exception
chaining}, use command line option C{python -X dev} I{OR} set env variable
C{PYTHONDEVMODE=1} or to any non-empyty string I{OR} set env variable
C{PYGEODESY_EXCEPTION_CHAINING=std} or to any non-empty string.
'''

from pygeodesy.interns import MISSING, NN, _a_, _an_, _and_, _clip_, \
                             _COLON_, _COLONSPACE_, _COMMASPACE_, _datum_, \
                             _ellipsoidal_, _EQUAL_, _incompatible_, _invalid_, \
                             _len_, _name_, _no_, _not_, _or_, _SPACE_, \
                             _specified_, _UNDER_, _value_, _vs_, _with_
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS, _getenv, \
                             _pairs, _PYTHON_X_DEV

__all__ = _ALL_LAZY.errors  # _ALL_DOCS('_InvalidError', '_IsnotError')  _under
__version__ = '23.08.30'

_box_        = 'box'
_default_    = 'default'
_kwargs_     = 'kwargs'  # XXX _kwds_?
_limiterrors =  True  # imported by .formy
_multiple_   = 'multiple'
_name_value_ =  repr('name=value')
_rangerrors  =  True  # imported by .dms
_region_     = 'region'
_vs__        = _SPACE_(NN, _vs_, NN)

try:
    _exception_chaining = None  # not available
    _ = Exception().__cause__   # Python 3+ exception chaining

    if _PYTHON_X_DEV or _getenv('PYGEODESY_EXCEPTION_CHAINING', NN):  # == _std_
        _exception_chaining = True  # turned on, std
        raise AttributeError  # allow exception chaining

    _exception_chaining = False  # turned off

    def _error_cause(inst, cause=None):
        '''(INTERNAL) Set or avoid Python 3+ exception chaining.

           Setting C{inst.__cause__ = None} is equivalent to syntax
           C{raise Error(...) from None} to avoid exception chaining.

           @arg inst: An error instance (I{caught} C{Exception}).
           @kwarg cause: A previous error instance (I{caught} C{Exception})
                         or C{None} to avoid exception chaining.

           @see: Alex Martelli, et.al., "Python in a Nutshell", 3rd Ed., page 163,
                 O'Reilly, 2017, U{PEP-3134<https://www.Python.org/dev/peps/pep-3134>},
                 U{here<https://StackOverflow.com/questions/17091520/how-can-i-more-
                 easily-suppress-previous-exceptions-when-i-raise-my-own-exception>}
                 and U{here<https://StackOverflow.com/questions/1350671/
                 inner-exception-with-traceback-in-python>}.
        '''
        inst.__cause__ = cause  # None, no exception chaining
        return inst

except AttributeError:  # Python 2+

    def _error_cause(inst, **unused):  # PYCHOK expected
        return inst  # no-op


class _AssertionError(AssertionError):
    '''(INTERNAL) Format an C{AssertionError} with/-out exception chaining.
    '''
    def __init__(self, *args, **kwds):
        _error_init(AssertionError, self, args, **kwds)


class _AttributeError(AttributeError):
    '''(INTERNAL) Format an C{AttributeError} with/-out exception chaining.
    '''
    def __init__(self, *args, **kwds):
        _error_init(AttributeError, self, args, **kwds)


class _ImportError(ImportError):
    '''(INTERNAL) Format an C{ImportError} with/-out exception chaining.
    '''
    def __init__(self, *args, **kwds):
        _error_init(ImportError, self, args, **kwds)


class _IndexError(IndexError):
    '''(INTERNAL) Format an C{IndexError} with/-out exception chaining.
    '''
    def __init__(self, *args, **kwds):
        _error_init(IndexError, self, args, **kwds)


class _KeyError(KeyError):
    '''(INTERNAL) Format a C{KeyError} with/-out exception chaining.
    '''
    def __init__(self, *args, **kwds):  # txt=_invalid_
        _error_init(KeyError, self, args, **kwds)


class _NameError(NameError):
    '''(INTERNAL) Format a C{NameError} with/-out exception chaining.
    '''
    def __init__(self, *args, **kwds):
        _error_init(NameError, self, args, **kwds)


class _NotImplementedError(NotImplementedError):
    '''(INTERNAL) Format a C{NotImplementedError} with/-out exception chaining.
    '''
    def __init__(self, *args, **kwds):
        _error_init(NotImplementedError, self, args, **kwds)


class _OverflowError(OverflowError):
    '''(INTERNAL) Format an C{OverflowError} with/-out exception chaining.
    '''
    def __init__(self, *args, **kwds):  # txt=_invalid_
        _error_init(OverflowError, self, args, **kwds)


class _TypeError(TypeError):
    '''(INTERNAL) Format a C{TypeError} with/-out exception chaining.
    '''
    def __init__(self, *args, **kwds):
        _error_init(TypeError, self, args, fmt_name_value='type(%s) (%r)', **kwds)


class _TypesError(_TypeError):
    '''(INTERNAL) Format a C{TypeError} with/-out exception chaining.
    '''
    def __init__(self, name, value, *Types, **kwds):
        t = _not_(_an(_or(*(t.__name__ for t in Types))))
        _TypeError.__init__(self, name, value, txt=t, **kwds)


class _ValueError(ValueError):
    '''(INTERNAL) Format a C{ValueError} with/-out exception chaining.
    '''
    def __init__(self, *args, **kwds):  # ..., cause=None, txt=_invalid_, ...
        _error_init(ValueError, self, args, **kwds)


class _ZeroDivisionError(ZeroDivisionError):
    '''(INTERNAL) Format a C{ZeroDivisionError} with/-out exception chaining.
    '''
    def __init__(self, *args, **kwds):
        _error_init(ZeroDivisionError, self, args, **kwds)


class AuxError(_ValueError):
    '''Error raised for an L{rhumbaux} C{Aux}, C{AuxDLat} or C{AuxLat} issue.
    '''
    pass


class ClipError(_ValueError):
    '''Clip box or clip region issue.
    '''
    def __init__(self, *name_n_corners, **txt_cause):
        '''New L{ClipError}.

           @arg name_n_corners: Either just a name (C{str}) or
                                name, number, corners (C{str},
                                C{int}, C{tuple}).
           @kwarg txt_cause: Optional C{B{txt}=str} explanation
                             of the error and C{B{cause}=None}
                             for exception chaining.
        '''
        if len(name_n_corners) == 3:
            t, n, v = name_n_corners
            n = _SPACE_(t, _clip_, (_box_ if n == 2 else _region_))
            name_n_corners = n, v
        _ValueError.__init__(self, *name_n_corners, **txt_cause)


class CrossError(_ValueError):
    '''Error raised for zero or near-zero vectorial cross products,
       occurring for coincident or colinear points, lines or bearings.
    '''
    pass


class GeodesicError(_ValueError):
    '''Error raised for L{pygeodesy.geodesicx} lack of convergence
       or other L{pygeodesy.geodesicx} or L{pygeodesy.karney} issues.
    '''
    pass


class IntersectionError(_ValueError):  # in .ellipsoidalBaseDI, .formy, ...
    '''Error raised for line or circle intersection issues.
    '''
    def __init__(self, *args, **kwds):
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
        def _ns_vs_txt_x(cause=None, txt=_invalid_, **kwds):
            ns, vs = zip(*itemsorted(kwds))  # unzip
            return ns, vs, txt, cause

        ns, vs, txt, x = _ns_vs_txt_x(**lens_txt)
        ns = _COMMASPACE_.join(ns)
        t  = _MODS.streprs.Fmt.PAREN(where.__name__, ns)
        vs = _vs__.join(map(str, vs))
        t  = _SPACE_(t, _len_, vs)
        _ValueError.__init__(self, t, txt=txt, cause=x)


class LimitError(_ValueError):
    '''Error raised for lat- or longitudinal values or deltas exceeding
       the given B{C{limit}} in functions L{pygeodesy.equirectangular},
       L{pygeodesy.equirectangular_}, C{nearestOn*} and C{simplify*}
       or methods with C{limit} or C{options} keyword arguments.

       @see: Subclass L{UnitError}.
    '''
    pass


class MGRSError(_ValueError):
    '''Military Grid Reference System (MGRS) parse or other L{Mgrs} issue.
    '''
    pass


class NumPyError(_ValueError):
    '''Error raised for C{NumPy} issues.
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


class RangeError(_ValueError):
    '''Error raised for lat- or longitude values outside the B{C{clip}},
       B{C{clipLat}}, B{C{clipLon}} in functions L{pygeodesy.parse3llh},
       L{pygeodesy.parseDMS}, L{pygeodesy.parseDMS2} and L{pygeodesy.parseRad}
       or the given B{C{limit}} in functions L{pygeodesy.clipDegrees} and
       L{pygeodesy.clipRadians}.

       @see: Function L{pygeodesy.rangerrors}.
    '''
    pass


class RhumbError(_ValueError):
    '''Error raised for L{pygeodesy.rhumbaux}, L{pygeodesy.rhumbsolve}
       or L{pygeodesy.rhumbx} issues.
    '''
    pass


class TriangleError(_ValueError):  # in .resections, .vector2d
    '''Error raised for triangle, inter- or resection issues.
    '''
    pass


class SciPyError(PointsError):
    '''Error raised for C{SciPy} issues.
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


class UnitError(LimitError):  # in .named, .units
    '''Default exception for L{units} issues for a value exceeding the
       C{low} or C{high} limit.
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


def _error_init(Error, inst, args, fmt_name_value='%s (%r)', txt=NN,
                                   cause=None, **kwds):  # by .lazily
    '''(INTERNAL) Format an error text and initialize an C{Error} instance.

       @arg Error: The error super-class (C{Exception}).
       @arg inst: Sub-class instance to be __init__-ed (C{_Exception}).
       @arg args: Either just a value or several name, value, ...
                  positional arguments (C{str}, any C{type}), in
                  particular for name conflicts with keyword
                  arguments of C{error_init} or which can't be
                  given as C{name=value} keyword arguments.
       @kwarg fmt_name_value: Format for (name, value) (C{str}).
       @kwarg txt: Optional explanation of the error (C{str}).
       @kwarg cause: Optional, caught error (L{Exception}), for
                     exception chaining (supported in Python 3+).
       @kwarg kwds: Additional C{B{name}=value} pairs, if any.
    '''
    def _fmtuple(pairs):
        return tuple(fmt_name_value % t for t in pairs)

    t, n = (), len(args)
    if n > 2:
        s = _MODS.basics.isodd(n)
        t = _fmtuple(zip(args[0::2], args[1::2]))
        if s:  # XXX _xzip(..., strict=s)
            t += args[-1:]
    elif n == 2:
        t = (fmt_name_value % args),
    elif n:  # == 1
        t =  str(args[0]),

    if kwds:
        t += _fmtuple(itemsorted(kwds))
    t = _or(*t) if t else _SPACE_(_name_value_, MISSING)

    if txt is not None:
        x =  str(txt) or (str(cause) if cause else _invalid_)
        C = _COMMASPACE_ if _COLON_ in t else _COLONSPACE_
        t =  C(t, x)
#   else:  # LenError, _xzip, .dms, .heights, .vector2d
#       x =  NN  # XXX or t?
    Error.__init__(inst, t)
#   inst.__x_txt__ = x  # hold explanation
    _error_cause(inst, cause=cause if _exception_chaining else None)
    _error_under(inst)


def _error_under(inst):
    '''(INTERNAL) Remove leading underscore from instance' class name.
    '''
    n = inst.__class__.__name__
    if n.startswith(_UNDER_):
        inst.__class__.__name__ = n.lstrip(_UNDER_)
    return inst


def exception_chaining(error=None):
    '''Get an error's I{cause} or the exception chaining setting.

       @kwarg error: An error instance (C{Exception}) or C{None}.

       @return: If C{B{error} is None}, return C{True} if exception
                chaining is enabled for PyGeodesy errors, C{False}
                if turned off and C{None} if not available.  If
                B{C{error}} is not C{None}, return it's error
                I{cause} or C{None}.

       @note: To enable exception chaining for C{pygeodesy} errors,
              set env var C{PYGEODESY_EXCEPTION_CHAINING} to any
              non-empty value prior to C{import pygeodesy}.
    '''
    return _exception_chaining if error is None else \
            getattr(error, '__cause__', None)


def _incompatible(this):
    '''(INTERNAL) Format an C{"incompatible with ..."} text.
    '''
    return _SPACE_(_incompatible_, _with_, this)


def _InvalidError(Error=_ValueError, **txt_name_values_cause):  # txt=_invalid_, name=value [, ...]
    '''(INTERNAL) Create an C{Error} instance.

       @kwarg Error: The error class or sub-class (C{Exception}).
       @kwarg txt_name_values: One or more C{B{name}=value} pairs
                  and optionally, keyword argument C{B{txt}=str}
                  to override the default C{B{txt}='invalid'} and
                  C{B{cause}=None} for exception chaining.

       @return: An B{C{Error}} instance.
    '''
    return _XError(Error, **txt_name_values_cause)


def isError(obj):
    '''Check a (caught) exception.

       @arg obj: The exception C({Exception}).

       @return: C{True} if B{C{obj}} is a C{pygeodesy} error,
                C{False} if B{C{obj}} is a standard Python error
                of C{None} if neither.
    '''
    return True  if isinstance(obj, _XErrors)   else (
           False if isinstance(obj,  Exception) else None)


def _IsnotError(*nouns, **name_value_Error_cause):  # name=value [, Error=TypeError, cause=None]
    '''Create a C{TypeError} for an invalid C{name=value} type.

       @arg nouns: One or more expected class or type names, usually nouns (C{str}).
       @kwarg name_value_Error_cause: One C{B{name}=value} pair and optionally,
                   keyword argument C{B{Error}=TypeError} to override the default
                   and C{B{cause}=None} for exception chaining.

       @return: A C{TypeError} or an B{C{Error}} instance.
    '''
    x, Error = _xkwds_pop_(name_value_Error_cause, cause=None, Error=TypeError)
    n, v =  _xkwds_popitem(name_value_Error_cause) if name_value_Error_cause \
                    else (_name_value_, MISSING)  # XXX else tuple(...)
    n = _MODS.streprs.Fmt.PARENSPACED(n, repr(v))
    t = _not_(_an(_or(*nouns)) if nouns else _specified_)
    return _XError(Error, n, txt=t, cause=x)


def itemsorted(adict, *args, **asorted_reverse):
    '''Return the items of C{B{adict}} sorted I{alphabetically, case-insensitively}
       and in I{ascending} order.

       @arg args: Optional argument(s) for method C{B{adict}.items(B*{args})}.
       @kwarg asorted_reverse: Use keyword argument C{B{asorted}=False} for
                      I{case-sensitive} sorting and C{B{reverse}=True} for
                      results in C{descending} order.
    '''
    def _un(item):
        return item[0].lower()

    # see .rhumb.Rhumb and ._RhumbLine
    a, r = _xkwds_get_(asorted_reverse, asorted=True, reverse=False) \
                    if asorted_reverse  else   (True,         False)
    items = adict.items(*args) if args else adict.items()
    return sorted(items, reverse=r, key=_un if a else None)


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
       @arg args: Any B{C{parser}} arguments (any C{type}s).
       @kwarg name_values_Error: Any C{B{name}=value} pairs and
                   optionally, C{B{Error}=ParseError} keyword
                   argument to override the default.

       @return: Parser result.

       @raise ParseError: Or the specified C{B{Error}}.
    '''
    try:
        return parser(*args)
    except Exception as x:
        E = _xkwds_pop(name_values_Error, Error=type(x) if isError(x) else
                                                ParseError)
        raise _XError(E, **_xkwds(name_values_Error, cause=x))


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
        Error = SciPyWarning
    else:
        Error = SciPyError  # PYCHOK not really
    t = _SPACE_(str(x).strip(), *extras)
    return Error(t, cause=x)


def _xattr(obj, **name_default):  # see .strerprs._xattrs
    '''(INTERNAL) Get an C{obj}'s attribute by C{name}.
    '''
    if len(name_default) == 1:
        for n, d in name_default.items():
            return getattr(obj, n, d)
    raise _xkwds_Error(_xattr, {}, name_default)


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


def _XError(Error, *args, **kwds):
    '''(INTERNAL) Format an C{Error} or C{_Error}.
    '''
    try:  # C{_Error} style
        return Error(*args, **kwds)
    except TypeError:  # no keyword arguments
        pass
    e = _ValueError(*args, **kwds)
    E =  Error(str(e))
    if _exception_chaining:
        _error_cause(E, cause=e.__cause__)  # PYCHOK OK
    return E


_XErrors = _TypeError, _ValueError
# map certain C{Exception} classes to the C{_Error}
_X2Error = {AssertionError:      _AssertionError,
            AttributeError:      _AttributeError,
            ImportError:         _ImportError,
            IndexError:          _IndexError,
            KeyError:            _KeyError,
            NameError:           _NameError,
            NotImplementedError: _NotImplementedError,
            OverflowError:       _OverflowError,
            TypeError:           _TypeError,
            ValueError:          _ValueError,
            ZeroDivisionError:   _ZeroDivisionError}


def _xError(x, *args, **kwds):
    '''(INTERNAL) Embellish a (caught) exception.

       @arg x: The exception (usually, C{_Error}).
       @arg args: Embelishments (C{any}).
       @kwarg kwds: Embelishments (C{any}).
    '''
    return _XError(type(x), *args, **_xkwds(kwds, cause=x))


def _xError2(x):  # in .fsums
    '''(INTERNAL) Map an exception to 2-tuple (C{_Error} class, error C{txt}).

       @arg x: The exception instance (usually, C{Exception}).
    '''
    X =  type(x)
    E = _X2Error.get(X, X)
    if E is X and not isError(x):
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


def _xkwds_bool(inst, **kwds):  # in .frechet, .hausdorff, .heights
    '''(INTERNAL) Set applicable C{bool} properties/attributes.
    '''
    for n, v in kwds.items():
        b = getattr(inst, n, None)
        if b is None:  # invalid bool attr
            t = _SPACE_(_EQUAL_(n, repr(v)), 'for', inst.__class__.__name__)  # XXX .classname
            raise _AttributeError(t, txt=_not_('applicable'))
        if v in (False, True) and v != b:
            setattr(inst, NN(_UNDER_, n), v)


def _xkwds_Error(where, kwds, name_txt, txt=_default_):
    # Helper for _xkwds_get, _xkwds_pop and _xkwds_popitem below
    f = _COMMASPACE_.join(_pairs(kwds) + _pairs(name_txt))
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
    '''(INTERNAL) Yield each C{kwds} value or its C{default}
       in I{case-insensitive, alphabetical} order.
    '''
    for n, d in itemsorted(names_defaults):
        yield kwds.get(n, d)


def _xkwds_not(*args, **kwds):
    '''(INTERNAL) Return C{kwds} with a value not in C{args}.
    '''
    return dict((n, v) for n, v in kwds.items() if v not in args)


def _xkwds_pop(kwds, **name_default):
    '''(INTERNAL) Pop a C{kwds} value by C{name}, or the C{default}.
    '''
    if len(name_default) == 1:
        for n, d in name_default.items():
            return kwds.pop(n, d)
    raise _xkwds_Error(_xkwds_pop, kwds, name_default)


def _xkwds_pop_(kwds, **names_defaults):
    '''(INTERNAL) Pop and yield each C{kwds} value or its C{default}
       in I{case-insensitive, alphabetical} order.
    '''
    for n, d in itemsorted(names_defaults):
        yield kwds.pop(n, d)


def _xkwds_popitem(name_value):
    '''(INTERNAL) Return exactly one C{(name, value)} item.
    '''
    if len(name_value) == 1:  # XXX TypeError
        return name_value.popitem()  # XXX AttributeError
    raise _xkwds_Error(_xkwds_popitem, (), name_value, txt=_value_)


def _Xorder(_Coeffs, Error, **Xorder):  # in .ktm, .rhumbBase
    '''(INTERNAL) Validate C{RAorder} or C{TMorder}.
    '''
    X, m = Xorder.popitem()
    if m in _Coeffs and _MODS.basics.isint(m):
        return m
    t = sorted(map(str, _Coeffs.keys()))
    raise Error(X, m, txt=_not_(_or(*t)))


def _xzip(*args, **strict):  # PYCHOK no cover
    '''(INTERNAL) Standard C{zip(..., strict=True)}.
    '''
    s = _xkwds_get(strict, strict=True)
    if s:
        _zip = _MODS.basics._zip
        if _zip is zip:  # < (3, 10)
            t = _MODS.streprs.unstr(_xzip.__name__, *args, strict=s)
            raise _NotImplementedError(t, txt=None)
        return _zip(*args)
    return zip(*args)

# **) MIT License
#
# Copyright (C) 2016-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
