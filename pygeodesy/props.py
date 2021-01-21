
# -*- coding: utf-8 -*-

u'''Caching, immutable and mutable properties and decorators.
'''
from pygeodesy.errors import _AssertionError, _AttributeError
from pygeodesy.interns import NN, _DOT_, _EQUALSPACED_, \
                             _immutable_, _invalid_, _N_A_, \
                             _SPACE_, _UNDER_
from pygeodesy.lazily import _ALL_LAZY, _FOR_DOCS

__all__ = _ALL_LAZY.props
__version__ = '21.01.21'


def _hasProperty(inst, name, *Classes):
    '''(INTERNAL) Check whether C{inst} has a C{P/property/_RO} by this C{name}.
    '''
    ps = Classes if Classes else _PropertyBase
    for c in inst.__class__.__mro__[:-1]:
        p = c.__dict__.get(name, None)
        if isinstance(p, ps) and p.name == name:
            return True
    return False


def _update_all(inst, *attrs):
    '''(INTERNAL) Zap all I{cached} L{property_RO}s, L{Property_RO}s
       and the named C{attrs} from an instance' C{__dict__}.

       @return: The number of updates (C{int}), if any.
    '''
    d = inst.__dict__
    u = len(d)

    for c in inst.__class__.__mro__[:-1]:
        for n, p in c.__dict__.items():
            if isinstance(p, _PropertyBase) and p.name == n:
                p._update(inst, c)

    for a in attrs:
        if hasattr(inst, a):
            d.pop(a, None)
        else:  # PYCHOK no cover
            from pygeodesy.named import classname
            n =  classname(inst, prefixed=True)
            a = _DOT_(n, _SPACE_(a, _invalid_))
            raise _AssertionError(a, txt=repr(inst))

    return u - len(d)


class _PropertyBase(property):
    '''(INTERNAL) Base class for C{P/property/_RO}.
    '''
    def __init__(self, method, fget, fset):

        if not callable(method):
            self.getter(method)  # PYCHOK no cover

        self.method = method
        self.name   = method.__name__
        if _FOR_DOCS and method.__doc__:
            self.__doc__ = method.__doc__   # PYCHOK no cover

        property.__init__(self, fget, fset, self._fdel, method.__doc__ or _N_A_)

    def _fdel(self, inst):  # deleter
        self._update(inst, None)   # PYCHOK no cover

    def _fget(self, inst):
        try:  # to get the value cached in instance' __dict__
            return inst.__dict__[self.name]
        except KeyError:
            # cache the value in the instance' __dict__
            inst.__dict__[self.name] = val = self.method(inst)
            return val

    def _fset_error(self, inst, val):
        '''Throws an C{AttributeError}, always.
        '''
        from pygeodesy.named import classname
        n = _DOT_(classname(inst), self.name)
        self._immutable_error(_EQUALSPACED_(n, repr(val)))

    def _immutable_error(self, name):
        from pygeodesy.named import classname
        e = _SPACE_(_immutable_, classname(self))
        raise _AttributeError(e, txt=name)

    def _invalid_error(self, name):
        from pygeodesy.named import classname
        e = _SPACE_(_invalid_, classname(self))
        raise _AttributeError(e, txt=name)

    def _update(self, inst, *unused):
        '''(INTERNAL) Zap the I{cached} C{inst.__dict__[name]} item.
        '''
        inst.__dict__.pop(self.name, None)

    def deleter(self, fdel):
        '''Throws an C{AttributeError}, always.
        '''
        n = _DOT_(self.name, self.deleter.__name__)
        self._invalid_error(_SPACE_(n, fdel.__name__))

    def getter(self, fget):
        '''Throws an C{AttributeError}, always.
        '''
        n = _DOT_(self.name, self.getter.__name__)
        self._invalid_error(_SPACE_(n, fget.__name__))

    def setter(self, fset):
        '''Throws an C{AttributeError}, always.
        '''
        n = _DOT_(self.name, self.setter.__name__)
        self._immutable_error(_SPACE_(n, fset.__name__))


class Property_RO(_PropertyBase):
    # No __doc__ on purpose
    def __init__(self, method):  # PYCHOK expected
        '''New I{immutable}, I{caching} C{property} I{Factory} to be used as C{decorator}.

           @arg method: The callable being decorated as this C{property}'s getter,
                        to be invoked only once.

           @note: Like standard Python C{property} without a setter, but with
                  a more descriptive error message when set.

           @see: Python 3's U{functools.cached_property<https://docs.Python.org/3/
                 library/functools.html#functools.cached_property>} and U{-.cache
                 <https://Docs.Python.org/3/library/functools.html#functools.cache>}
                 to I{memoize} the property value.

           @see: Luciano Ramalho, "Fluent Python", page 636, O'Reilly, 2016,
                 "Coding a Property Factory", especially Example 19-24 and U{class
                 Property<https://docs.Python.org/3/howto/descriptor.html>}.
        '''
        _PropertyBase.__init__(self, method, self._fget, self._fset_error)

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
    '''New I{mutable}, I{caching} C{property} I{Factory} to be used as C{decorator}.

       @see: L{Property_RO} for some more details.

       @note: Unless and until the C{setter} is defined, this L{Property} behaves
              like an I{immutable}, I{caching} L{Property_RO}.
    '''

    def setter(self, method):
        '''Make this C{Property} I{mutable}.

           @arg method: The callable being decorated as this C{Property}'s setter.
        '''
        if not callable(method):
            _PropertyBase.setter(self, method)  # PYCHOK no cover

        def _fset(inst, val):
            self._update(inst)  # uncache this item
            method(inst, val)

        # class Property <https://docs.Python.org/3/howto/descriptor.html>
        _PropertyBase.__init__(self, self.method, self._fget, _fset)
        return self


class property_RO(_PropertyBase):
    # No __doc__ on purpose
    def __init__(self, method):  # PYCHOK expected
        '''New I{immutable}, standard C{property} to be used as C{decorator}.

           @arg method: The callable being decorated as C{property}'s getter.

           @note: Like standard Python C{property} without a setter, but with
                  a more descriptive error message when set.

           @see: L{Property_RO}.
        '''
        _PropertyBase.__init__(self, method, method, self._fset_error)

    def _update(self, inst, clas=None):  # PYCHOK signature
        '''(INTERNAL) Zap the I{cached} C{inst.__dict__[_name]} item.
        '''
        c = clas or inst.__class__
        if c:
            _n = NN(_UNDER_, self.name)
            if c.__dict__.get(_n, False) is None:  # in (None, ()):
                inst.__dict__.pop(_n, None)


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
