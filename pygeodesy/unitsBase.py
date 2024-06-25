
# -*- coding: utf-8 -*-

u'''Basic C{Float}, C{Int} and C{Str}ing units classes.
'''

from pygeodesy.basics import isstr, issubclassof, _xsubclassof
from pygeodesy.errors import _IsnotError, _UnexpectedError, UnitError, _XError
from pygeodesy.interns import NN, _degrees_, _degrees2_, _invalid_, _meter_, \
                             _radians_, _radians2_, _radius_, _UNDER_, _units_, \
                             _std_  # PYCHOK used!
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import modulename, _Named,  property_doc_
# from pygeodesy.props import property_doc_  # from .named
from pygeodesy.streprs import Fmt, fstr

__all__ = _ALL_LAZY.unitsBase
__version__ = '24.06.15'


class _NamedUnit(_Named):
    '''(INTERNAL) Base class for C{units}.
    '''
    _std_repr = True  # set below
    _units    = None

    def __new__(cls, typ, arg, name, Error=UnitError, **name_arg):
        '''(INTERNAL) Return a named C{typ.__new__(cls, arg)} instance.
        '''
        if name_arg:
            name, arg = _NamedUnit._arg_name_arg2(arg, name, **name_arg)
        try:  # assert typ in cls.__mro__
            self = typ.__new__(cls, arg)
            if name:
                self.name = name
        except Exception as x:
            raise _NamedUnit._Error(cls, arg, name, Error, cause=x)
        return self

    @staticmethod
    def _arg_name_arg2(arg, name=NN, name__=None, **name_arg):  # in .units
        '''(INTERNAL) Get the 2-tuple C{(name, arg)}.
        '''
        if name_arg:
            if len(name_arg) > 1:
                raise _UnexpectedError(**name_arg)
            for name, arg in name_arg.items():  # next(iter(.items()))
                break
        elif name:
            pass
        elif name__ is not None:
            name = name__.__name__
        return name, arg

    @staticmethod  # PYCHOK unused suffix
    def _Error(cls, arg, name, Error=UnitError, suffix=NN,  # unused
                                txt=_invalid_, cause=None, **name_arg):
        '''(INTERNAL) Return a C{_NamedUnit} error with explanation.

           @returns: An B{C{Error}} instance.
        '''
        kwds, x = {}, cause
        if x is not None:  # caught exception
            if Error is UnitError:  # and isError(x)
                Error = type(x)  # i.e. not overridden
            if txt is _invalid_:
                txt = str(x)  # i.e. not overridden
            kwds.update(cause=x)
        if name_arg:
            try:
                name, arg = _NamedUnit._arg_name_arg2(arg, name, **name_arg)
            except Exception:  # ignore, same error?
                kwds.update(name_arg)
        n = name if name else modulename(cls).lstrip(_UNDER_)
        return _XError(Error, n, arg, txt=txt, **kwds)

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

    def __new__(cls, arg=None, name=NN, **Error_name_arg):
        '''New, named C{Ffloat}.

           @kwarg arg: The value (any C{type} acceptable to C{float}).
           @kwarg name: Optional name (C{str}).
           @kwarg Error_name_arg: Optional C{B{Error}=UnitError} to raise
                        and optional C{name=arg} keyword argument, inlieu
                        of separate B{C{arg}} and B{C{name}} ones.

           @returns: A named C{Float}.

           @raise Error: Invalid B{C{arg}}.
        '''
        return _NamedUnit.__new__(cls, float, arg, name, **Error_name_arg)

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
        # must use super(Float, self)... since super()... only works
        # for Python 3+ and float.__str__(self) invokes .__repr__(self);
        # calling self.toRepr(std=True) super(Float, self).__repr__()
        # mimicks this bhavior

        # XXX the default number of decimals is 10-12 when using
        # float.__str__(self) with both python 3.8+ and 2.7-, but
        # float.__repr__(self) shows DIG decimals in python2.7!
        # return super(Float, self).__repr__()  # see .testCss.py
        return float.__str__(self)  # always _std_str_

    def toRepr(self, std=False, **prec_fmt_ints):  # PYCHOK prec=8, ...
        '''Return a representation of this C{Float}.

           @kwarg std: If C{True} return the standard C{repr},
                       otherwise the named representation (C{bool}).

           @see: Function L{fstr<pygeodesy.streprs.fstr>} and methods
                 L{Float.__repr__}, L{Float.toStr} for further details.
        '''
        # must use super(Float, self)... since super()... only works for
        # Python 3+; return super(Float, self).__repr__() if std else \
        return float.__repr__(self) if std else \
               self._toRepr(self.toStr(**prec_fmt_ints))

    def toStr(self, prec=12, fmt=Fmt.g, ints=False):  # PYCHOK prec=8, ...
        '''Format this C{Float} as C{str}.

           @see: Function L{fstr<pygeodesy.streprs.fstr>} and method
                 L{Float.__repr__} and for further information.
        '''
        return fstr(self, prec=prec, fmt=fmt, ints=ints)


class Int(int, _NamedUnit):
    '''Named C{int}.
    '''
    # _std_repr = True  # set below

    def __new__(cls, arg=None, name=NN, **Error_name_arg):
        '''New, named C{Int}.

           @kwarg arg: The value (any C{type} acceptable to C{int}).
           @kwarg name: Optional name (C{str}).
           @kwarg Error_name_arg: Optional C{B{Error}=UnitError} to raise
                        and optional C{name=arg} keyword argument, inlieu
                        of separate B{C{arg}} and B{C{name}} ones.

           @returns: A named C{Int}.

           @raise Error: Invalid B{C{arg}}.
        '''
        return _NamedUnit.__new__(cls, int, arg, name, **Error_name_arg)

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

    def __new__(cls, arg=None, name=NN, **Error_name_arg):
        '''New, named and callable C{Str}.

           @kwarg cls: This class (C{Str} or sub-class).
           @kwarg arg: The value (any C{type} acceptable to C{str}).
           @kwarg name: Optional name (C{str}).
           @kwarg Error_name_arg: Optional C{B{Error}=UnitError} to raise
                        and optional C{name=arg} keyword argument, inlieu
                        of separate B{C{arg}} and B{C{name}} ones.

           @returns: A named L{Str}.

           @raise Error: Invalid B{C{arg}}.

           @see: Callable, not nameable class L{Str_<pygeodesy.interns.Str_>}.
        '''
        return _NamedUnit.__new__(cls, str, arg, name, **Error_name_arg)

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
        # must use super(Str, self)... since super()... only works
        # for Python 3+ and str.__repr__(self) fails in Python 3.8+
        r = super(Str, self).__repr__()
        return r if std else self._toRepr(r)

    def toStr(self, **unused):  # PYCHOK **unused
        '''Return this C{Str} as standard C{str}.
        '''
        # must use super(Str, self)... since super()... only works
        # for Python 3+ and str.__repr__(self) fails in Python 3.8+
        return super(Str, self).__str__()


_Str_degrees  = Str(_degrees_)   # PYCHOK in .frechet, .hausdorff
_Str_degrees2 = Str(_degrees2_)  # PYCHOK in .frechet, .hausdorff
_Str_meter    = Str(_meter_)     # PYCHOK in .frechet, .hausdorff
_Str_NN       = Str(NN)          # PYCHOK in .frechet, .hausdorff
_Str_radians  = Str(_radians_)   # PYCHOK in .frechet, .hausdorff
_Str_radians2 = Str(_radians2_)  # PYCHOK in .frechet, .hausdorff


def _xUnit(units, Base):  # PYCHOK in .frechet,  .hausdorff
    '''(INTERNAL) Get C{Unit} from C{units} or C{name}, ortherwise C{Base}.
    '''
    _xsubclassof(_NamedUnit, Base=Base)
    U = getattr(_MODS, units.capitalize(), Base) if isstr(units) else units
    return U if issubclassof(U, Base) else Base


def _xUnits(units, Base=_NamedUnit):  # in .frechet, .hausdorff
    '''(INTERNAL) Set property C{units} as C{Unit} or C{Str}.
    '''
    _xsubclassof(_NamedUnit, Base=Base)
    if issubclassof(units, Base):
        U = units
    elif isstr(units):
        U = Str(units, name=_units_)  # XXX Str to _Pass and for backward compatibility
    else:
        raise _IsnotError(Base, Str, str, units=units)
    return U


__all__ += _ALL_DOCS(_NamedUnit)

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
