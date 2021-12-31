
# -*- coding: utf-8 -*-

u'''Im-/mutable, caching or memoizing properties and deprecation decorators.

To enable C{DeprecationWarning}s from C{PyGeodesy}, set environment
variable C{PYGEODESY_WARNINGS} to a non-empty string and run C{python}
with command line option C{-X dev} or one or the C{-W} choices, see
function L{DeprecationWarnings} below.
'''
from pygeodesy.errors import _AssertionError, _AttributeError, \
                             _xkwds, _xkwds_get
from pygeodesy.interns import NN, _DOT_, _EQUALSPACED_, _immutable_, \
                             _invalid_, MISSING, _N_A_, _NLNL_, \
                             _SPACE_, _UNDER_
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS, _getenv, \
                             _FOR_DOCS, _PYTHON_X_DEV, _sys

from functools import wraps as _wraps

__all__ = _ALL_LAZY.props
__version__ =  '21.12.28'

_DEPRECATED_ = 'DEPRECATED'
_dont_use_   = _DEPRECATED_ + ", don't use."
_has_been_   = 'has been'
_Warnings    =  0
_W_DEV       = (bool(_sys.warnoptions) or _PYTHON_X_DEV) \
                and _getenv('PYGEODESY_WARNINGS', NN)


def _hasProperty(inst, name, *Classes):
    '''(INTERNAL) Check whether C{inst} has a C{P/property/_RO} by this C{name}.
    '''
    ps = Classes if Classes else _PropertyBase
    for c in inst.__class__.__mro__[:-1]:
        p = c.__dict__.get(name, None)
        if isinstance(p, ps) and p.name == name:
            return True
    return False


def _update_all(inst, *attrs, **Base):
    '''(INTERNAL) Zap all I{cached} L{property_RO}s, L{Property_RO}s
       and the named C{attrs} from an instance' C{__dict__}.

       @return: The number of updates (C{int}), if any.
    '''
    B = _xkwds_get(Base, Base=_PropertyBase)
    d = inst.__dict__
    u = len(d)

    for c in inst.__class__.__mro__[:-1]:
        for n, p in c.__dict__.items():
            if isinstance(p, B) and p.name == n:
                p._update(inst, c)

    p = d.pop
    for a in attrs:  # PYCHOK no cover
        if hasattr(inst, a):
            p(a, None)
        else:
            n = _MODS.named.classname(inst, prefixed=True)
            a = _DOT_(n, _SPACE_(a, _invalid_))
            raise _AssertionError(a, txt=repr(inst))

    return u - len(d)  # of updates


class _PropertyBase(property):
    '''(INTERNAL) Base class for C{P/property/_RO}.
    '''
    def __init__(self, method, fget, fset, doc=NN):

        if not callable(method):
            self.getter(method)  # PYCHOK no cover

        self.method = method
        self.name   = method.__name__
        d = doc or method.__doc__
        if _FOR_DOCS and d:
            self.__doc__ = d   # PYCHOK no cover

        property.__init__(self, fget, fset, self._fdel, d or _N_A_)

    def _fdel(self, inst):
        '''Zap the I{cached/memoized} C{property} value.
        '''
        self._update(inst, None)   # PYCHOK no cover

    def _fget(self, inst):
        '''Get and I{cache/memoize} the C{property} value.
        '''
        try:  # to get the value cached in instance' __dict__
            return inst.__dict__[self.name]
        except KeyError:
            # cache the value in the instance' __dict__
            inst.__dict__[self.name] = val = self.method(inst)
            return val

    def _fset_error(self, inst, val):
        '''Throws an C{AttributeError}, always.
        '''
        n = _MODS.named.classname(inst)
        n = _DOT_(n, self.name)
        n = _EQUALSPACED_(n, repr(val))
        raise self._Error(_immutable_, n, None)

    def _update(self, inst, *unused):
        '''(INTERNAL) Zap the I{cached/memoized} C{inst.__dict__[name]} item.
        '''
        inst.__dict__.pop(self.name, None)  # name, NOT _name

    def deleter(self, fdel):
        '''Throws an C{AttributeError}, always.
        '''
        raise self._Error(_invalid_, self.deleter, fdel)

    def getter(self, fget):
        '''Throws an C{AttributeError}, always.
        '''
        raise self._Error(_invalid_, self.getter, fget)

    def setter(self, fset):
        '''Throws an C{AttributeError}, always.
        '''
        raise self._Error(_immutable_, self.setter, fset)

    def _Error(self, kind, nameter, farg):
        '''(INTERNAL) Return an C{AttributeError} instance.
        '''
        if farg:
            n = _DOT_(self.name, nameter.__name__)
            n = _SPACE_(n, farg.__name__)
        else:
            n = nameter
        e = _SPACE_(kind, _MODS.named.classname(self))
        return _AttributeError(e, txt=n)


class Property_RO(_PropertyBase):
    # No __doc__ on purpose
    def __init__(self, method, doc=NN):  # PYCHOK expected
        '''New I{immutable}, I{caching}, I{memoizing} C{property} I{Factory}
           to be used as C{decorator}.

           @arg method: The callable being decorated as this C{property}'s C{getter},
                        to be invoked only once.
           @kwarg doc: Optional property documentation (C{str}).

           @note: Like standard Python C{property} without a C{setter}, but with
                  a more descriptive error message when set.

           @see: Python 3's U{functools.cached_property<https://docs.Python.org/3/
                 library/functools.html#functools.cached_property>} and U{-.cache
                 <https://Docs.Python.org/3/library/functools.html#functools.cache>}
                 to I{cache} or I{memoize} the property value.

           @see: Luciano Ramalho, "Fluent Python", page 636, O'Reilly, 2016,
                 "Coding a Property Factory", especially Example 19-24 and U{class
                 Property<https://docs.Python.org/3/howto/descriptor.html>}.
        '''
        _fget = method if _FOR_DOCS else self._fget  # XXX force method.__doc__ to epydoc
        _PropertyBase.__init__(self, method, _fget, self._fset_error)

    def __get__(self, inst, *unused):  # PYCHOK 2 vs 3 args
        if inst is None:
            return self
        try:  # to get the cached value immediately
            return inst.__dict__[self.name]
        except (AttributeError, KeyError):
            return self._fget(inst)


class Property(Property_RO):
    # No __doc__ on purpose
    __init__ = Property_RO.__init__
    '''New I{mutable}, I{caching}, I{memoizing} C{property} I{Factory}
       to be used as C{decorator}.

       @see: L{Property_RO} for more details.

       @note: Unless and until the C{setter} is defined, this L{Property} behaves
              like an I{immutable}, I{caching}, I{memoizing} L{Property_RO}.
    '''

    def setter(self, method):
        '''Make this C{Property} I{mutable}.

           @arg method: The callable being decorated as this C{Property}'s C{setter}.

           @note: Setting a new property value always clears the previously I{cached}
                  or I{memoized} value I{after} invoking the B{C{method}}.
        '''
        if not callable(method):
            _PropertyBase.setter(self, method)  # PYCHOK no cover

        if _FOR_DOCS:  # XXX force method.__doc__ into epydoc
            _PropertyBase.__init__(self, self.method, self.method, method)
        else:

            def _fset(inst, val):
                '''Set and I{cache}, I{memoize} the C{property} value.
                '''
                method(inst, val)
                self._update(inst)  # un-cache this item

            # class Property <https://docs.Python.org/3/howto/descriptor.html>
            _PropertyBase.__init__(self, self.method, self._fget, _fset)
        return self


class property_RO(_PropertyBase):
    # No __doc__ on purpose
    _uname = NN

    def __init__(self, method, doc=NN):  # PYCHOK expected
        '''New I{immutable}, standard C{property} to be used as C{decorator}.

           @arg method: The callable being decorated as C{property}'s C{getter}.
           @kwarg doc: Optional property documentation (C{str}).

           @note: Like standard Python C{property} without a setter, but with
                  a more descriptive error message when set.

           @see: L{Property_RO}.
        '''
        _PropertyBase.__init__(self, method, method, self._fset_error, doc=doc)
        self._uname = NN(_UNDER_, self.name)  # actual attr UNDER<name>

    def _update(self, inst, *Clas):  # PYCHOK signature
        '''(INTERNAL) Zap the I{cached} C{B{inst}.__dict__[_name]} item.
        '''
        uname = self._uname
        if uname in inst.__dict__:
            if Clas:  # overrides inst.__class__
                d = Clas[0].__dict__.get(uname, MISSING)
            else:
                d = getattr(inst.__class__, uname, MISSING)
#               if d is MISSING:  # XXX superfluous
#                   for c in inst.__class__.__mro__[:-1]:
#                       if uname in c.__dict__:
#                           d = c.__dict__[uname]
#                           break
            if d is None:  # remove inst value
                inst.__dict__.pop(uname)


def property_doc_(doc):
    '''Decorator for a standard C{property} with basic documentation.

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

    def _documented_property(method):
        '''(INTERNAL) Return the documented C{property}.
        '''
        t = 'get and set' if doc.startswith(_SPACE_) else NN
        return property(method, None, None, NN('Property to ', t, doc))

    return _documented_property


def _deprecated(call, kind, qual_d):
    '''(INTERNAL) Decorator for DEPRECATED functions, methods, etc.

       @see: Brett Slatkin, "Effective Python", page 105, 2nd ed,
             Addison-Wesley, 2019.
    '''
    doc = _docof(call)

    @_wraps(call)  # PYCHOK self?
    def _deprecated_call(*args, **kwds):
        if qual_d:  # function
            q = qual_d
        elif args:  # method
            q = _qualified(args[0], call.__name__)
        else:  # PYCHOK no cover
            q = call.__name__
        _throwarning(kind, q, doc)
        return call(*args, **kwds)

    return _deprecated_call


def deprecated_class(cls_or_class):
    '''Use inside __new__ or __init__ of a DEPRECATED class.

       @arg cls_or_class: The class (C{cls} or C{Class}).

       @note: NOT a decorator!
    '''
    if _W_DEV:
        q = _DOT_(cls_or_class.__module__, cls_or_class.__name__)
        _throwarning('class', q, cls_or_class.__doc__)


def deprecated_function(call):
    '''Decorator for a DEPRECATED function.

       @arg call: The deprecated function (C{callable}).

       @return: The B{C{call}} DEPRECATED.
    '''
    return _deprecated(call, 'function', _DOT_(
                       call.__module__, call.__name__)) if \
           _W_DEV else call


def deprecated_method(call):
    '''Decorator for a DEPRECATED method.

       @arg call: The deprecated method (C{callable}).

       @return: The B{C{call}} DEPRECATED.
    '''
    return _deprecated(call, 'method', NN) if _W_DEV else call


def _deprecated_module(name):  # PYCHOK no cover
    '''(INTERNAL) Callable within a DEPRECATED module.
    '''
    if _W_DEV:
        _throwarning('module', name, _dont_use_)


def deprecated_Property_RO(method):
    '''Decorator for a DEPRECATED L{Property_RO}.

       @arg method: The C{Property_RO.fget} method (C{callable}).

       @return: The B{C{method}} DEPRECATED.
    '''
    return _deprecated_RO(method, Property_RO)


def deprecated_property_RO(method):
    '''Decorator for a DEPRECATED L{property_RO}.

       @arg method: The C{property_RO.fget} method (C{callable}).

       @return: The B{C{method}} DEPRECATED.
    '''
    return _deprecated_RO(method, property_RO)


def _deprecated_RO(method, _RO):
    '''(INTERNAL) Create a DEPRECATED C{property_RO} or C{Property_RO}.
    '''
    doc = _docof(method)

    if _W_DEV:

        class _Deprecated_RO(_PropertyBase):
            __doc__ = doc

            def __init__(self, method):
                _PropertyBase.__init__(self, method, self._fget, self._fset_error, doc=doc)

            def _fget(self, inst):  # PYCHOK no cover
                q = _qualified(inst, self.name)
                _throwarning(_RO.__name__, q, doc)
                return self.method(inst)

        return _Deprecated_RO(method)
    else:  # PYCHOK no cover
        return _RO(method, doc=doc)


def DeprecationWarnings():
    '''Have any C{DeprecationWarning}s been reported or raised?

       @return: The number of C{DeprecationWarning}s (C{int}) so
                far or C{None} if not enabled.

       @note: To get C{DeprecationWarning}s if any, run C{python}
              with environment variable C{PYGEODESY_WARNINGS} set
              to a non-empty string I{AND} use C{python[3]} command
              line option C{-X dev}, C{-W always} or C{-W error}, etc.
    '''
    return _Warnings if _W_DEV else None


def _docof(obj):
    '''(INTERNAL) Get uniform DEPRECATED __doc__ string.
    '''
    try:
        d = obj.__doc__.strip()
        i = d.find(_DEPRECATED_)
    except AttributeError:
        i = -1
    return _DOT_(_DEPRECATED_, NN) if i < 0 else d[i:]


def _qualified(inst, name):
    '''(INTERNAL) Fully qualify a name.
    '''
    # _DOT_(inst.classname, name), not _DOT_(inst.named4, name)
    c =  inst.__class__
    q = _DOT_(c.__module__, c.__name__, name)
    return q


def _throwarning(kind, name, doc, **stacklevel):  # stacklevel=3
    '''(INTERNAL) Report or raise a C{DeprecationWarning}.
    '''
    from warnings import warn

    line =  doc.split(_NLNL_, 1)[0].split()
    name = _MODS.streprs.Fmt.CURLY(L=name)
    text = _SPACE_(kind, name, _has_been_, *line)
    kwds = _xkwds(stacklevel, stacklevel=3)
    # XXX invoke warn or raise DeprecationWarning(text)
    warn(text, category=DeprecationWarning, **kwds)

    global _Warnings
    _Warnings += 1

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
