
# -*- coding: utf-8 -*-

u'''Basic C{Float}, C{Int} and C{Str}ing units classes.
'''

from pygeodesy.errors import UnitError, _XError, _xkwds_popitem
from pygeodesy.interns import NN, _degrees_, _degrees2_, _invalid_, \
                             _meter_, _radians_, _radians2_, \
                             _radius_, _UNDER_,  _std_  # PYCHOK used!
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import modulename, _Named, property_doc_
# from pygeodesy.props import property_doc_  # from .named
from pygeodesy.streprs import Fmt, fstr

__all__ = _ALL_LAZY.unitsBase
__version__ = '23.05.10'


class _NamedUnit(_Named):
    '''(INTERNAL) Base class for C{units}.
    '''
    _std_repr = True  # set below
    _units    = None

    @property_doc_(' standard C{repr} or named C{toRepr} representation.')
    def std_repr(self):
        '''Get the representation (C{bool}, C{True} means standard).
        '''
        return self._std_repr

    @std_repr.setter  # PYCHOK setter!
    def std_repr(self, std):
        '''Set the representation (C{True} or C{"std"} for standard).
        '''
        self._std_repr = std in (True, _std_)

    def _toRepr(self, value):
        '''(INTERNAL) Representation "<name> (<value>)" or "<classname>(<value>)".
        '''
        return Fmt.PARENSPACED(self.name, value) if self.name else \
               Fmt.PAREN( self.classname, value)

    @property_doc_(' units name.')
    def units(self):
        '''Get the units name (C{str}).
        '''
        if self._units is None:
            self._units = self.classname
        return self._units

    @units.setter  # PYCHOK setter!
    def units(self, units):
        '''Set the units name for this instance (C{str} or C{None} for default).
        '''
        self._units = None if units is None else str(units)


class Float(float, _NamedUnit):
    '''Named C{float}.
    '''
    # _std_repr = True  # set below

    def __new__(cls, arg=None, name=NN, Error=UnitError, **name_arg):
        '''New C{Ffloat} instance.

           @kwarg arg: The value (any C{type} convertable to C{float}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the default
                         L{UnitError}.
           @kwarg name_arg: Optional C{name=arg} keyword argument, inlieu
                            of B{C{name}} and B{C{arg}}.

           @returns: A C{Float} instance.

           @raise Error: Invalid B{C{arg}}.
        '''
        if name_arg:
            name, arg = _xkwds_popitem(name_arg)
        try:
            self = float.__new__(cls, arg)
            if name:
                _NamedUnit.name.fset(self, name)  # see _Named.name
        except Exception as x:  # XXX not ... as x:
            raise _Error(cls, arg, name, Error, x=x)
        return self

    def __repr__(self):  # to avoid MRO(float)
        '''Return a representation of this C{Float}.

           @see: Method C{Float.toRepr} and property C{Float.std_repr}.

           @note: Use C{env} variable C{PYGEODESY_FLOAT_STD_REPR=std} prior
                  to C{import pygeodesy} to get the standard C{repr} or set
                  property C{std_repr=False} to always get the named C{toRepr}
                  representation.
        '''
        return self.toRepr(std=self._std_repr)

    def __str__(self):  # to avoid MRO(float)
        '''Return this C{Float} as standard C{str}.
        '''
        # XXX must use super(Float, self)... since super()...
        # only works for Python 3+ and float.__str__(self)
        # invokes .__repr__(self); calling self.toRepr(std=True)
        # super(Float, self).__repr__() mimicks this behavior
        # XXX the default number of decimals is 10-12 when using
        # float.__str__(self) with both python 3.8+ and 2.7-, but
        # float.__repr__(self) shows DIG decimals in python2.7!
        # return super(Float, self).__repr__()  # see .test.testCss.py
        return float.__str__(self)  # always _std_str_

    def toRepr(self, std=False, **prec_fmt_ints):  # PYCHOK prec=8, ...
        '''Return a representation of this C{Float}.

           @kwarg std: If C{True} return the standard C{repr},
                       otherwise the named representation (C{bool}).

           @see: Methods L{Float.__repr__}, L{Float.toStr} and function
                 L{pygeodesy.fstr} for more documentation.
        '''
        # XXX must use super(Float, self)... since
        # super()... only works for Python 3+
        # return super(Float, self).__repr__() if std else \
        return float.__repr__(self) if std else \
               self._toRepr(self.toStr(**prec_fmt_ints))

    def toStr(self, prec=12, fmt=Fmt.g, ints=False):  # PYCHOK prec=8, ...
        '''Format this C{Float} as C{str}.

           @see: Method L{Float.__repr__} and function L{pygeodesy.fstr}
                 for more documentation.
        '''
        return fstr(self, prec=prec, fmt=fmt, ints=ints)


class Int(int, _NamedUnit):
    '''Named C{int}.
    '''
    # _std_repr = True  # set below

    def __new__(cls, arg=None, name=NN, Error=UnitError, **name_arg):
        '''New C{Int} instance.

           @kwarg arg: The value (any C{type} convertable to C{float}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the
                         default L{UnitError}.
           @kwarg name_arg: Optional C{name=arg} keyword argument,
                            inlieu of B{C{name}} and B{C{arg}}.

           @returns: An C{Int} instance.

           @raise Error: Invalid B{C{arg}}.
        '''
        if name_arg:
            name, arg = _xkwds_popitem(name_arg)
        try:
            self = int.__new__(cls, arg)
            if name:
                _NamedUnit.name.fset(self, name)  # see _Named.name
        except Exception as x:  # XXX not ... as x:
            raise _Error(cls, arg, name, Error, x=x)
        return self

    def __repr__(self):  # to avoid MRO(int)
        '''Return a representation of this named C{int}.

           @see: Method C{Int.toRepr} and property C{Int.std_repr}.

           @note: Use C{env} variable C{PYGEODESY_INT_STD_REPR=std}
                  prior to C{import pygeodesy} to get the standard
                  C{repr} or set property C{std_repr=False} to always
                  get the named C{toRepr} representation.
        '''
        return self.toRepr(std=self._std_repr)

    def __str__(self):  # to avoid MRO(int)
        '''Return this C{Int} as standard C{str}.
        '''
        return self.toStr()

    def toRepr(self, std=False, **unused):  # PYCHOK **unused
        '''Return a representation of this C{Int}.

           @kwarg std: If C{True} return the standard C{repr},
                       otherwise the named representation (C{bool}).

           @see: Method L{Int.__repr__} for more documentation.
        '''
        r = int.__repr__(self)  # self.toStr()
        return r if std else self._toRepr(r)

    def toStr(self, **unused):  # PYCHOK **unused
        '''Return this C{Int} as standard C{str}.

           @see: Method L{Int.__repr__} for more documentation.
        '''
        # XXX must use '%d' % (self,) since
        # int.__str__(self) fails with 3.8+
        return '%d' % (self,)


class Radius(Float):
    '''Named C{float} representing a radius, conventionally in C{meter}.
    '''
    def __new__(cls, arg=None, name=_radius_, **Error_name_arg):
        '''New L{Radius} instance, see L{Float}.
        '''
        return Float.__new__(cls, arg=arg, name=name, **Error_name_arg)


class Str(str, _NamedUnit):
    '''Named, callable C{str}.
    '''
    # _std_repr = True  # set below

    def __new__(cls, arg=None, name=NN, Error=UnitError, **name_arg):
        '''New  C{Str} instance.

           @kwarg cls: This class (C{Str} or sub-class).
           @kwarg arg: The value (any C{type} convertable to C{str}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the
                         default (C{ValueError}).
           @kwarg name_arg: Optional C{name=arg} keyword argument,
                            inlieu of B{C{name}} and B{C{arg}}.

           @returns: A L{Str} instance.

           @raise Error: Invalid B{C{arg}}.

           @see: Callable, not-nameable class L{pygeodesy.Str_}.
        '''
        if name_arg:
            name, arg = _xkwds_popitem(name_arg)
        try:
            self = str.__new__(cls, arg)
            if name:
                _NamedUnit.name.fset(self, name)  # see _Named.name
        except Exception as x:  # XXX not ... as x:
            raise _Error(cls, arg, name, Error, x=x)
        return self

    def __repr__(self):
        '''Return a representation of this C{Str}.

           @see: Method C{Str.toRepr} and property C{Str.std_repr}.

           @note: Use C{env} variable C{PYGEODESY_STR_STD_REPR=std}
                  prior to C{import pygeodesy} to get the standard
                  C{repr} or set property C{std_repr=False} to always
                  get the named C{toRepr} representation.
        '''
        return self.toRepr(std=self._std_repr)  # see .test/testGars.py

    def __str__(self):
        '''Return this C{Str} as standard C{str}.
        '''
        return self.toStr()

    def join_(self, *args, **name_Error):
        '''Join all positional B{C{args}} like C{self.join(B{args})}.

           @return: All B{C{args}} joined by this instance (L{Str_}).

           @note: An other L{Str} instance is returned to make the
                  result re-callable.
        '''
        return Str(str.join(self, map(str, args)), **name_Error)  # re-callable

    __call__ = join_

    def toRepr(self, std=False, **unused):  # PYCHOK **unused
        '''Return a representation of this C{Str}.

           @kwarg std: If C{True} return the standard C{repr},
                       otherwise the named representation (C{bool}).

           @see: Method L{Str.__repr__} for more documentation.
        '''
        # must use super(Str, self).. since
        # super()... only works for Python 3+ and
        # str.__repr__(self) fails with Python 3.8+
        r = super(Str, self).__repr__()
        return r if std else self._toRepr(r)

    def toStr(self, **unused):  # PYCHOK **unused
        '''Return this C{Str} as standard C{str}.
        '''
        # must use super(Str, self)... since
        # super()... only works for Python 3+ and
        # str.__str__(self) fails with Python 3.8+
        return super(Str, self).__str__()


_Str_degrees  = Str(_degrees_)   # PYCHOK in .frechet, .hausdorff
_Str_degrees2 = Str(_degrees2_)  # PYCHOK in .frechet, .hausdorff
_Str_meter    = Str(_meter_)     # PYCHOK in .frechet, .hausdorff
_Str_NN       = Str(NN)          # PYCHOK in .frechet, .hausdorff
_Str_radians  = Str(_radians_)   # PYCHOK in .frechet, .hausdorff
_Str_radians2 = Str(_radians2_)  # PYCHOK in .frechet, .hausdorff


def _Error(clas, arg, name, Error, txt=_invalid_, x=None):
    '''(INTERNAL) Return an error with explanation.

       @arg clas: The C{units} class or sub-class.
       @arg arg: The original C{unit} value.
       @arg name: The instance name (C{str}).
       @arg Error: The Error class to use (C{Excetion}).
       @kwarg txt: An explanation of the error )C{str}).
       @kwarg x: Caught exception, used for exception
                 chaining (iff enabled in Python3+).

       @returns: An B{C{Error}} instance.
    '''
    if x is not None:  # caught exception, cause
        if Error is UnitError:  # and isError(x)
            Error = type(x)  # i.e. if not overridden
        if txt is _invalid_:
            txt = str(x)  # i.e. if not overridden
    n = name if name else modulename(clas).lstrip(_UNDER_)
    return _XError(Error, n, arg, txt=txt, cause=x)


__all__ += _ALL_DOCS(_NamedUnit)

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
