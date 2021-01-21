
# -*- coding: utf-8 -*-

u'''Sub-classes C{Float}, C{Int} and C{Str} from basic
C{float}, C{int} respectively C{str} to named units as
L{Degrees}, L{Feet}, L{Meter}, L{Radians}, etc.

@newfield example: Example, Examples
'''

from pygeodesy.basics import isstr, issubclassof
from pygeodesy.dms import F__F, F__F_, parseDMS, parseRad, \
                          S_NUL, S_SEP, _toDMS
from pygeodesy.errors import _IsnotError, RangeError, TRFError, \
                              UnitError, _xkwds_popitem
from pygeodesy.interns import EPS, EPS1, NN, PI, PI_2, _band_, \
                             _bearing_, _degrees_, _degrees2_, \
                             _distance_, _E_, _easting_, _epoch_, \
                             _EW_, _feet_, _height_, _invalid_, _N_, \
                             _lam_, _lat_, _LatLon_, _lon_, _meter_, \
                             _meter2_, _northing_, _NS_, _NSEW_, \
                             _number_, _PERCENT_, _phi_, _precision_, \
                             _radians_, _radians2_, _radius_, _S_, \
                             _scalar_, _SPACE_, _UNDER_, _units_, \
                             _W_, _zone_, _0_0, _0_001
from pygeodesy.interns import _std_  # PYCHOK used!
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import modulename, _Named
from pygeodesy.props import Property_RO, property_doc_
from pygeodesy.streprs import Fmt, fstr

from math import radians

__all__ = _ALL_LAZY.units
__version__ = '21.01.19'


class _NamedUnit(_Named):
    '''(INTERNAL) Base class for C{units}.
    '''
    _std_repr = True  # set below
    _units    = None

    def _toRepr(self, value):
        '''(INTERNAL) Representation "<name> (<value>)" or "<classname>(<value>)".
        '''
        t = NN(self.name, _SPACE_) if self.name else self.classname
        return Fmt.PAREN(t, value)

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
        '''New named C{float} instance.

           @kwarg arg: The value (any C{type} convertable to C{float}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the
                         default L{UnitError}.
           @kwarg name_arg: Optional C{name=arg} keyword argument,
                            inlieu of B{C{name}} and B{C{arg}}.

           @returns: A C{Float} instance.

           @raise Error: Invalid B{C{arg}}.
        '''
        if name_arg:
            name, arg = _xkwds_popitem(name_arg)
        try:
            self = float.__new__(cls, arg)
            if name:
                _NamedUnit.name.fset(self, name)  # see _Named.name
            return self

        except (TypeError, ValueError) as x:  # XXX not ... as x:
            raise _Error(cls, arg, name=name, Error=Error, txt=str(x))

    def __repr__(self):  # to avoid MRO(float)
        '''Return a representation of this C{Float}.

           @see: Method C{Float.toRepr} and property C{Float.std_repr}.

           @note: Use C{env} variable C{PYGEODESY_FLOAT_STD_REPR=std} prior
                  to C{import pygeodesy} to get the standard C{repr} or use
                  C{-=named} for the named C{toRepr} representation.
        '''
        return self.toRepr(std=self._std_repr)

    def __str__(self):  # to avoid MRO(float)
        '''Return this C{Float} as standard C{str}.
        '''
        # XXX must use super(Float, self)... since super()...
        # only works for Python 3+ and float.__str__(self)
        # invokes .__repr__(self); calling self.toRepr(std=True)
        # super(Float, self).__repr__() mimicks this behavior
        return super(Float, self).__repr__()  # see .test.testCss.py

    def toRepr(self, prec=12, fmt=Fmt.g, ints=False, std=False):  # PYCHOK prec=8, ...
        '''Return a representation of this C{Float}.

           @kwarg std: Use the standard C{repr} or the named
                       representation (C{bool}).

           @see: Function L{fstr} for more documentation.
        '''
        # XXX must use super(Float, self)... since
        # super()... only works for Python 3+
        return super(Float, self).__repr__() if std else \
               self._toRepr(fstr(self, prec=prec, fmt=fmt, ints=ints))

    def toStr(self, prec=12, fmt=Fmt.g, ints=False):  # PYCHOK prec=8, ...
        '''Format this C{Float} as C{str}.

           @see: Function L{fstr} for more documentation.
        '''
        return fstr(self, prec=prec, fmt=fmt, ints=ints)


class Float_(Float):
    '''Named C{float} with optional C{low} and C{high} limit.
    '''
    def __new__(cls, arg=None, name=NN, Error=UnitError, low=EPS, high=None, **name_arg):
        '''New named C{float} instance with limits.

           @arg cls: This class (C{Float_} or sub-class).
           @kwarg arg: The value (any C{type} convertable to C{float}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the default
                         L{UnitError}.
           @kwarg low: Optional lower B{C{arg}} limit (C{float} or C{None}).
           @kwarg high: Optional upper B{C{arg}} limit (C{float} or C{None}).
           @kwarg name_arg: Optional C{name=arg} keyword argument, inlieu of
                            B{C{name}} and B{C{arg}}.

           @returns: A C{Float_} instance.

           @raise Error: Invalid B{C{arg}} or B{C{arg}} below B{C{low}} or
                         above B{C{high}}.
        '''
        if name_arg:
            name, arg = _xkwds_popitem(name_arg)
        self = Float.__new__(cls, arg=arg, name=name, Error=Error)
        if (low is not None) and self < low:
            txt = Fmt.limit(below=Fmt.g(low, prec=6, ints=isinstance(self, Epoch)))
        elif (high is not None) and self > high:
            txt = Fmt.limit(above=Fmt.g(high, prec=6, ints=isinstance(self, Epoch)))
        else:
            return self
        raise _Error(cls, arg, name=name, Error=Error, txt=txt)


class Int(int, _NamedUnit):
    '''Named C{int}.
    '''
    # _std_repr = True  # set below

    def __new__(cls, arg=None, name=NN, Error=UnitError, **name_arg):
        '''New named C{int} instance.

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
            return self

        except (TypeError, ValueError) as x:  # XXX not ... as x:
            raise _Error(cls, arg, name=name, Error=Error, txt=str(x))

    def __repr__(self):  # to avoid MRO(int)
        '''Return a representation of this named C{int}.

           @see: Method C{Int.toRepr} and property C{Int.std_repr}.

           @note: Use C{env} variable C{PYGEODESY_INT_STD_REPR=std}
                  prior to C{import pygeodesy} to get the standard
                  C{repr} or use C{-=named} for the named C{toRepr}
                  representation.
        '''
        return self.toRepr(std=self._std_repr)

    def __str__(self):  # to avoid MRO(int)
        '''Return this C{Int} as standard C{str}.
        '''
        return self.toStr()

    def toRepr(self, std=False, **unused):  # PYCHOK **unused
        '''Return the representation of thisb is True C{Int}.

           @kwarg std: Use the standard C{repr} or the named
                       representation (C{bool}).
        '''
        r = int.__repr__(self)  # self.toStr()
        return r if std else self._toRepr(r)

    def toStr(self, **unused):  # PYCHOK **unused
        '''Return this C{Int} as standard C{str}.
        '''
        # XXX must use '%d' % (self,) since
        # int.__str__(self) fails with 3.8+
        return '%d' % (self,)


class Int_(Int):
    '''Named C{int} with optional limits C{low} and C{high}.
    '''
    def __new__(cls, arg=None, name=NN, Error=UnitError, low=0, high=None, **name_arg):
        '''New named C{int} instance with limits.

           @kwarg cls: This class (C{Int_} or sub-class).
           @arg arg: The value (any C{type} convertable to C{int}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the default
                         C{UnitError}.
           @kwarg low: Optional lower B{C{arg}} limit (C{float} or C{None}).
           @kwarg high: Optional upper B{C{arg}} limit (C{float} or C{None}).
           @kwarg name_arg: Optional C{name=arg} keyword argument,
                            inlieu of B{C{name}} and B{C{arg}}.

           @returns: An L{Int_} instance.

           @raise Error: Invalid B{C{arg}} or B{C{arg}} below B{C{low}} or
                         above B{C{high}}.
        '''
        if name_arg:
            name, arg = _xkwds_popitem(name_arg)
        self = Int.__new__(cls, arg=arg, name=name, Error=Error)
        if (low is not None) and self < low:
            txt = Fmt.limit(below=low)
        elif (high is not None) and self > high:
            txt = Fmt.limit(above=high)
        else:
            return self
        raise _Error(cls, arg, name=name, Error=Error, txt=txt)


class Bool(Int, _NamedUnit):
    '''Named C{bool}, a sub-class of C{int} like Python's C{bool}.
    '''
    # _std_repr = True  # set below
    _bool_True_or_False = None

    def __new__(cls, arg=None, name=NN, Error=UnitError, **name_arg):
        '''New named C{bool} instance.

           @kwarg cls: This class (C{Bool} or sub-class).
           @kwarg arg: The value (any C{type} convertable to C{bool}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the
                         default L{UnitError}.
           @kwarg name_arg: Optional C{name=arg} keyword argument,
                            inlieu of B{C{name}} and B{C{arg}}.

           @returns: A L{Bool}, a C{bool}-like instance.

           @raise Error: Invalid B{C{arg}}.
        '''
        if name_arg:
            name, arg = _xkwds_popitem(name_arg)
        try:
            b = bool(arg)
        except (TypeError, ValueError) as x:  # XXX not ... as x:
            raise _Error(cls, arg, name=name, Error=Error, txt=str(x))

        self = Int.__new__(cls, b, name=name, Error=Error)
        self._bool_True_or_False = b
        return self

    # <https://StackOverflow.com/questions/9787890/assign-class-boolean-value-in-python>
    def __bool__(self):  # PYCHOK Python 3+
        return self._bool_True_or_False
    __nonzero__ = __bool__  # PYCHOK Python 2-

    def toRepr(self, std=False, **unused):  # PYCHOK **unused
        '''Return the representation of this C{Bool}.

           @kwarg std: Use the standard C{repr} or the named
                       representation (C{bool}).

           @note: Use C{env} variable C{PYGEODESY_BOOL_STD_REPR=std}
                  prior to C{import pygeodesy} to get the standard
                  C{repr} or use C{-=named} for the named C{toRepr}
                  representation.
        '''
        r = repr(self._bool_True_or_False)  # self.toStr()
        return r if std else self._toRepr(r)

    def toStr(self, **unused):  # PYCHOK **unused
        '''Return this C{Bool} as standard C{str}.
        '''
        return str(self._bool_True_or_False)


class Str(str, _NamedUnit):
    '''Named C{str}.
    '''
    # _std_repr = True  # set below

    def __new__(cls, arg=None, name=NN, Error=UnitError, **name_arg):
        '''New named C{str} instance.

           @kwarg cls: This class (C{Str} or sub-class).
           @kwarg arg: The value (any C{type} convertable to C{str}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding
                         the default (C{ValueError}).
           @kwarg name_arg: Optional C{name=arg} keyword argument,
                            inlieu of B{C{name}} and B{C{arg}}.

           @returns: A L{Str} instance.

           @raise Error: Invalid B{C{arg}}.
        '''
        if name_arg:
            name, arg = _xkwds_popitem(name_arg)
        try:
            self = str.__new__(cls, arg)
            if name:
                _NamedUnit.name.fset(self, name)  # see _Named.name
            return self

        except (TypeError, ValueError) as x:  # XXX not ... as x:
            raise _Error(cls, arg, name=name, Error=Error, txt=str(x))

    def __repr__(self):
        '''Return a representation of this C{Str}.

           @see: Method C{Str.toRepr} and property C{Str.std_repr}.

           @note: Use C{env} variable C{PYGEODESY_STR_STD_REPR=std}
                  prior to C{import pygeodesy} to get the standard
                  C{repr} or use C{-=named} for the named C{toRepr}
                  representation.
        '''
        return self.toRepr(std=self._std_repr)  # see .test/testGars.py

    def __str__(self):
        '''Return this C{Str} as standard C{str}.
        '''
        return self.toStr()

    def toRepr(self, std=False, **unused):  # PYCHOK **unused
        '''Return the named representation of this C{Str}.

           @kwarg std: Use the standard C{repr} or the named
                       representation (C{bool}).
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
_Str_meter    = Str(_meter_)     # PYCHOK in .frechet, .hausdorff
_Str_NN       = Str(NN)          # PYCHOK in .frechet, .hausdorff
_Str_radians  = Str(_radians_)   # PYCHOK in .frechet, .hausdorff
_Str_radians2 = Str(_radians2_)  # PYCHOK in .frechet, .hausdorff


class Band(Str):
    '''Named C{str} representing a UTM/UPS band letter, unchecked.
    '''
    def __new__(cls, arg=None, name=_band_, **Error_name_arg):
        '''See L{Str}.
        '''
        return Str.__new__(cls, arg=arg, name=name, **Error_name_arg)


class Degrees(Float):
    '''Named C{float} representing a coordinate in C{degrees}, optionally clipped.
    '''
    _ddd_ = 1  # default for .dms._toDMS
    _suf_ = S_NUL, S_NUL, S_NUL
    _sep_ = S_SEP

    def __new__(cls, arg=None, name=_degrees_, Error=UnitError, suffix=_NSEW_, clip=0, **name_arg):
        '''New named C{Degrees} instance.

           @arg cls: This class (C{Degrees} or sub-class).
           @kwarg arg: The value (any C{type} convertable to C{float} or
                       parsable by L{parseDMS}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the default
                         L{UnitError}.
           @kwarg suffix: Optional, valid compass direction suffixes (C{NSEW}).
           @kwarg clip: Optional B{C{arg}} range B{C{-clip..+clip}}
                        (C{degrees} or C{0} or C{None} for unclipped).
           @kwarg name_arg: Optional C{name=arg} keyword argument,
                            inlieu of B{C{name}} and B{C{arg}}.

           @returns: A C{Degrees} instance.

           @raise Error: Invalid B{C{arg}} or B{C{abs(deg)}} outside the
                         B{C{clip}} range and L{rangerrors} set to C{True}.
        '''
        if name_arg:
            name, arg = _xkwds_popitem(name_arg)
        try:
            return Float.__new__(cls, parseDMS(arg, suffix=suffix, clip=clip),
                                      name=name, Error=Error)
        except RangeError as x:
            t, E = str(x), type(x)
        except (TypeError, ValueError) as x:
            t, E = str(x), Error
        raise _Error(cls, arg, name=name, Error=E, txt=t)

    def toRepr(self, prec=None, fmt=F__F_, ints=False, std=False):  # PYCHOK prec=8, ...
        return Float.toRepr(self, std=std) if std else \
               self._toRepr(self.toStr(prec=prec, fmt=fmt, ints=ints))

    def toStr(self, prec=None, fmt=F__F_, ints=False):  # PYCHOK prec=8, ...
        if fmt.startswith(_PERCENT_):  # use regular formatting
            p = 8 if prec is None else prec
            return fstr(self, prec=p, fmt=fmt, ints=ints, sep=self._sep_)
        else:
            return _toDMS(self, fmt, prec, self._sep_, self._ddd_,
                                           self._suf_[0 if self > 0 else
                                                     (1 if self < 0 else 2)])


class Degrees_(Degrees):
    '''Named C{Degrees} representing a coordinate in C{degrees} with optional limits C{low} and C{high}.
    '''
    def __new__(cls, arg=None, name=_degrees_, Error=UnitError, suffix=_NSEW_, low=None, high=None, **name_arg):
        '''New named C{Degrees} instance.

           @arg cls: This class (C{Degrees_} or sub-class).
           @kwarg arg: The value (any C{type} convertable to C{float} or
                       parsable by L{parseDMS}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the default
                         L{UnitError}.
           @kwarg suffix: Optional, valid compass direction suffixes (C{NSEW}).
           @kwarg low: Optional lower B{C{arg}} limit (C{float} or C{None}).
           @kwarg high: Optional upper B{C{arg}} limit (C{float} or C{None}).
           @kwarg name_arg: Optional C{name=arg} keyword argument,
                            inlieu of B{C{name}} and B{C{arg}}.

           @returns: A C{Degrees} instance.

           @raise Error: Invalid B{C{arg}} or B{C{abs(deg)}} below B{C{low}}
                         or above B{C{high}}.
        '''
        if name_arg:
            name, arg = _xkwds_popitem(name_arg)
        self = Degrees.__new__(cls, arg=arg, name=name, Error=Error, suffix=suffix, clip=0)
        if (low is not None) and self < low:
            txt = Fmt.limit(below=low)
        elif (high is not None) and self > high:
            txt = Fmt.limit(above=high)
        else:
            return self
        raise _Error(cls, arg, name=name, Error=Error, txt=txt)


class Degrees2(Float):
    '''Named C{float} representing a distance in C{degrees squared}.
    '''
    def __new__(cls, arg=None, name=_degrees2_, **Error_name_arg):
        '''See L{Float}.
        '''
        return Float.__new__(cls, arg=arg, name=name, **Error_name_arg)


class Radians(Float):
    '''Named C{float} representing a coordinate in C{radians}, optionally clipped.
    '''
    def __new__(cls, arg=None, name=_radians_, Error=UnitError, suffix=_NSEW_, clip=0, **name_arg):
        '''New named C{Radians} instance.

           @arg cls: This class (C{Radians} or sub-class).
           @kwarg arg: The value (any C{type} convertable to C{float} or
                       parsable by L{parseRad}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the default
                         L{UnitError}.
           @kwarg suffix: Optional, valid compass direction suffixes (C{NSEW}).
           @kwarg clip: Optional B{C{arg}} range B{C{-clip..+clip}}
                        (C{radians} or C{0} or C{None} for unclipped).
           @kwarg name_arg: Optional C{name=arg} keyword argument,
                            inlieu of B{C{name}} and B{C{arg}}.

           @returns: A C{Radians} instance.

           @raise Error: Invalid B{C{arg}} or B{C{abs(deg)}} outside the
                         B{C{clip}} range and L{rangerrors} set to C{True}.
        '''
        if name_arg:
            name, arg = _xkwds_popitem(name_arg)
        try:
            return Float.__new__(cls, parseRad(arg, suffix=suffix, clip=clip),
                                      name=name, Error=Error)
        except RangeError as x:
            t, E = str(x), type(x)
        except (TypeError, ValueError) as x:
            t, E = str(x), Error
        raise _Error(cls, arg, name=name, Error=E, txt=t)

    def toRepr(self, prec=8, fmt=F__F, ints=False, std=False):  # PYCHOK prec=8, ...
        return Float.toRepr(self, prec=prec, fmt=fmt, ints=ints, std=std)

    def toStr(self, prec=8, fmt=F__F, ints=False):  # PYCHOK prec=8, ...
        return fstr(self, prec=prec, fmt=fmt, ints=ints)


class Radians2(Float_):
    '''Named C{float} representing a distance in C{radians squared}.
    '''
    def __new__(cls, arg=None, name=_radians2_, Error=UnitError, **name_arg):
        '''See L{Float_}.
        '''
        return Float_.__new__(cls, arg=arg, name=name, Error=Error, low=_0_0, **name_arg)


class Bearing(Degrees):
    '''Named C{float} representing a bearing in compass C{degrees} from (true) North.
    '''
    _ddd_ =  1
    _suf_ = _N_ * 3  # always suffix N

    def __new__(cls, arg=None, name=_bearing_, Error=UnitError, clip=0, **name_arg):
        '''See L{Degrees}.
        '''
        if name_arg:
            name, arg = _xkwds_popitem(name_arg)
        d = Degrees.__new__(cls, arg=arg, name=name, Error=Error, suffix=_N_, clip=clip)
        b = d % 360
        return d if b == d else Degrees.__new__(cls, arg=b, name=name, Error=Error)


class Bearing_(Radians):
    '''Named C{float} representing a bearing in C{radians} from compass C{degrees} from (true) North.
    '''
    def __new__(cls, arg=None, name=_bearing_, clip=0, **Error_name_arg):
        '''See L{Bearing} and L{Radians}.
        '''
        d = Bearing.__new__(cls, arg=arg, name=name, clip=clip, **Error_name_arg)
        return Radians.__new__(cls, radians(d), name=name)


class Distance(Float):
    '''Named C{float} representing a distance, conventionally in C{meter}.
    '''
    def __new__(cls, arg=None, name=_distance_, **Error_name_arg):
        '''See L{Distance}.
        '''
        return Float.__new__(cls, arg=arg, name=name, **Error_name_arg)


class Distance_(Float_):
    '''Named C{float} with optional C{low} and C{high} limits representing a distance, conventionally in C{meter}.
    '''
    def __new__(cls, arg=None, name=_distance_, low=EPS, **high_Error_name_arg):
        '''See L{Distance_}.
        '''
        return Float_.__new__(cls, arg=arg, name=name, low=low, **high_Error_name_arg)


class Easting(Float):
    '''Named C{float} representing an easting, conventionally in C{meter}.
    '''
    def __new__(cls, arg=None, name=_easting_, Error=UnitError, falsed=False, osgr=False, **name_arg):
        '''New named C{Easting} instance.

           @arg cls: This class (C{Easting} or sub-class).
           @kwarg arg: The value (any C{type} convertable to C{float}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the default
                         L{UnitError}.
           @kwarg falsed: The B{C{arg}} value includes false origin (C{bool}).
           @kwarg osgr: Check B{C{arg}} as an OSGR easting (C{bool}).
           @kwarg name_arg: Optional C{name=arg} keyword argument,
                            inlieu of B{C{name}} and B{C{arg}}.

           @returns: An C{Easting} instance.

           @raise Error: Invalid B{C{arg}} or negative, falsed B{C{arg}}.
        '''
        if name_arg:
            name, arg = _xkwds_popitem(name_arg)
        self = Float.__new__(cls, arg=arg, name=name, Error=Error)
        if osgr and (self < 0 or self > 700e3):  # like Veness
            raise _Error(cls, arg, name=name, Error=Error)
        elif falsed and self < 0:
            raise _Error(cls, arg, name=name, Error=Error, txt='negative, falsed')
        return self


class Epoch(Float_):  # by .ellipsoidalBase
    '''Named C{epoch} with optional C{low} and C{high} limits representing a fractional calendar year.
    '''

    _std_repr = False

    def __new__(cls, arg=None, name=_epoch_, Error=TRFError, low=1900, high=9999, **name_arg):
        '''See L{Float_}.
        '''
        if name_arg:
            name, arg = _xkwds_popitem(name_arg)
        return arg if isinstance(arg, Epoch) else \
               Float_.__new__(cls, arg=arg, name=name, Error=Error, low=low, high=high)

    def toRepr(self, std=False, **unused):  # PYCHOK prec=3, fmt=Fmt.F, ints=True
        '''Return a representation of this C{Epoch}.

           @kwarg std: Use the standard C{repr} or the named
                       representation (C{bool}).

           @see: Function L{fstr} for more documentation.
        '''
        return Float_.toRepr(self, prec=-3, fmt=Fmt.F, ints=True, std=std)

    def toStr(self, **unused):  # PYCHOK prec=3, fmt=Fmt.F, ints=True
        '''Format this C{Epoch} as C{str}.

           @see: Function L{fstr} for more documentation.
        '''
        return Float_.toStr(self, prec=-3, fmt=Fmt.F, ints=True)

    __str__ = toStr  # PYCHOK default '%.3F', without trailing zeros and decimal point


class Feet(Float):
    '''Named C{float} representing a distance or length in C{feet}.
    '''
    def __new__(cls, arg=None, name=_feet_, **Error_name_arg):
        '''See L{Float}.
        '''
        return Float.__new__(cls, arg=arg, name=name, **Error_name_arg)


class FIx(Float_):
    '''A named I{Fractional Index}, an C{int} or C{float} index into
       a C{list} or C{tuple} of C{points}, typically.  A C{float}
       I{Fractional Index} C{fi} represents a location on the edge
       between C{points[int(fi)]} and C{points[(int(fi) + 1) %
       len(points)]}.
    '''
    _fin = 0

    def __new__(cls, fi, fin=None, **name_Error):
        '''New I{Fractional Index} in a C{list} or C{tuple} of points.

           @arg fi: The fractional index (C{float} or C{int}).
           @kwarg fin: Optional C{len}, the number of C{points}, the index
                       C{[n]} wrapped to C{[0]} (C{int} or C{None}).
           @kwarg name_Error: Optional keyword argument C{name=<name>}
                              and/or C{Error=<Exception>}.

           @return: The B{C{fi}} (named L{FIx}).

           @note: The returned B{C{fi}} may exceed the B{C{flen}} of
                  the original C{points} in certain open/closed cases.

           @see: Method L{fractional} or function L{pygeodesy.fractional}.
        '''
        n = Int_(fin=fin, low=0) if fin else None
        f = Float_.__new__(cls, fi, low=_0_0, high=n, **name_Error)
        i = int(f)
        r = f - float(i)
        if r < EPS:  # see .points._fractional
            f = Float_.__new__(cls, i)
        elif r > EPS1:
            f = Float_.__new__(cls, i + 1, high=n, **name_Error)
        if n:
            f._fin = n
        return f

    @Property_RO
    def fin(self):
        '''Get the given C{len}, the index C{[n]} wrapped to C{[0]} (C{int}).
        '''
        return self._fin

    def fractional(self, points, **LatLon_LatLon_kwds):
        '''Return the point at this I{Fractional Index}.

           @arg points: The points (C{LatLon}[], L{Numpy2LatLon}[],
                        L{Tuple2LatLon}[] or C{other}[]).
           @kwarg **LatLon_LatLon_kwds: Optional class to return the
                                        I{intermediate}, I{fractional}
                                        point (C{LatLon}) or C{None}
                                        and optional B{C{LatLon}} keyword
                                        arguments thereof.

           @return: See function L{pygeodesy.fractional}.

           @raise IndexError: Fractional index invalid for B{C{points}}
                              or B{C{points}} not subscriptable or not
                              closed.
        '''
        from pygeodesy.points import fractional
        # fi = 0 if self == self.fin else self
        return fractional(points, self, **LatLon_LatLon_kwds)


class Height(Float):  # here to avoid circular import
    '''Named C{float} representing a height, conventionally in C{meter}.
    '''
    def __new__(cls, arg=None, name=_height_, **Error_name_arg):
        '''See L{Float}.
        '''
        return Float.__new__(cls, arg=arg, name=name, **Error_name_arg)


class Lam(Radians):
    '''Named C{float} representing a longitude in C{radians}.
    '''
    def __new__(cls, arg=None, name=_lam_, clip=PI, **Error_name_arg):
        '''See L{Radians}.
        '''
        return Radians.__new__(cls, arg=arg, name=name, suffix=_EW_, clip=clip, **Error_name_arg)


class Lam_(Lam):
    '''Named C{float} representing a longitude in C{radians} converted from C{degrees}.
    '''
    def __new__(cls, arg=None, name=_lon_, clip=180, **Error_name_arg):
        '''See L{Degrees} and L{Radians}.
        '''
        d = Lam.__new__(cls, arg=arg, name=name, clip=clip, **Error_name_arg)
        return Radians.__new__(cls, radians(d), name=name)


class Lat(Degrees):
    '''Named C{float} representing a latitude in C{degrees}.
    '''
    _ddd_ =  2
    _suf_ = _N_, _S_, S_NUL  # no zero suffix

    def __new__(cls, arg=None, name=_lat_, clip=90, **Error_name_arg):
        '''See L{Degrees}.
        '''
        return Degrees.__new__(cls, arg=arg, name=name, suffix=_NS_, clip=clip, **Error_name_arg)


class Lat_(Degrees_):
    '''Named C{float} representing a latitude in C{degrees} within limits C{low} and C{high}.
    '''
    _ddd_ =  2
    _suf_ = _N_, _S_, S_NUL  # no zero suffix

    def __new__(cls, arg=None, name=_lat_, low=-90, high=90, **Error_name_arg):
        '''See L{Degrees_}.
        '''
        return Degrees_.__new__(cls, arg=arg, name=name, suffix=_NS_, low=low, high=high, **Error_name_arg)


class Lon(Degrees):
    '''Named C{float} representing a longitude in C{degrees}.
    '''
    _ddd_ =  3
    _suf_ = _E_, _W_, S_NUL  # no zero suffix

    def __new__(cls, arg=None, name=_lon_, clip=180, **Error_name_arg):
        '''See L{Degrees}.
        '''
        return Degrees.__new__(cls, arg=arg, name=name, suffix=_EW_, clip=clip, **Error_name_arg)


class Lon_(Degrees_):
    '''Named C{float} representing a longitude in C{degrees} within limits C{low} and C{high}.
    '''
    _ddd_ =  3
    _suf_ = _E_, _W_, S_NUL  # no zero suffix

    def __new__(cls, arg=None, name=_lon_, low=-180, high=180, **Error_name_arg):
        '''See L{Degrees_}.
        '''
        return Degrees_.__new__(cls, arg=arg, name=name, suffix=_EW_, low=low, high=high, **Error_name_arg)


class Meter(Float):
    '''Named C{float} representing a distance or length in C{meter}.
    '''
    def __new__(cls, arg=None, name=_meter_, **Error_name_arg):
        '''See L{Float}.
        '''
        return Float.__new__(cls, arg=arg, name=name, **Error_name_arg)


_1mm    = Meter(   _1mm=_0_001)  # PYCHOK 1 millimeter in .ellipsoidal...
_10um   = Meter(  _10um= 1e-5)   # PYCHOK 0.01 millimeter in .osgr
_100km  = Meter( _100km= 1e+5)   # PYCHOK 100 kilometer in .formy, .mgrs, .osgr
_2000km = Meter(_2000km= 2e+6)   # PYCHOK 2,000 kilometer in .mgrs


class Meter2(Float_):
    '''Named C{float} representing an area in C{meter squared}.
    '''
    def __new__(cls, arg=None, name=_meter2_, Error=UnitError, **name_arg):
        '''See L{Float_}.
        '''
        return Float_.__new__(cls, arg=arg, name=name, Error=Error, low=_0_0, **name_arg)


class Meter3(Float_):
    '''Named C{float} representing a volume in C{meter cubed}.
    '''
    def __new__(cls, arg=None, name='meter3', Error=UnitError, **name_arg):
        '''See L{Float_}.
        '''
        return Float_.__new__(cls, arg=arg, name=name, Error=Error, low=_0_0, **name_arg)


class Northing(Float):
    '''Named C{float} representing a northing, conventionally in C{meter}.
    '''
    def __new__(cls, arg=None, name=_northing_, Error=UnitError, falsed=False, osgr=False, **name_arg):
        '''New named C{Northing} instance.

           @arg cls: This class (C{Northing} or sub-class).
           @kwarg arg: The value (any C{type} convertable to C{float}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the default
                         L{UnitError}.
           @kwarg falsed: The B{C{arg}} value includes false origin (C{bool}).
           @kwarg osgr: Check B{C{arg}} as an OSGR northing (C{bool}).
           @kwarg name_arg: Optional C{name=arg} keyword argument,
                            inlieu of B{C{name}} and B{C{arg}}.

           @returns: A C{Northing} instance.

           @raise Error: Invalid B{C{arg}} or negative, falsed B{C{arg}}.
        '''
        if name_arg:
            name, arg = _xkwds_popitem(name_arg)
        self = Float.__new__(cls, arg=arg, name=name, Error=Error)
        if osgr and (self < 0 or self > 1300e3):  # like Veness
            raise _Error(cls, arg, name=name, Error=Error)
        elif falsed and self < 0:
            raise _Error(cls, arg, name=name, Error=Error, txt='negative, falsed')
        return self


class Number_(Int_):
    '''Named C{int} representing a non-negative number.
    '''
    def __new__(cls, arg=None, name=_number_, low=0, **high_Error_name_arg):
        '''See L{Int_}.
        '''
        return Int_.__new__(cls, arg=arg, name=name, low=low, **high_Error_name_arg)


class Phi(Radians):
    '''Named C{float} representing a latitude in C{radians}.
    '''
    def __new__(cls, arg=None, name=_phi_, clip=PI_2, **Error_name_arg):
        '''See L{Radians}.
        '''
        return Radians.__new__(cls, arg=arg, name=name, suffix=_NS_, clip=clip, **Error_name_arg)


class Phi_(Phi):
    '''Named C{float} representing a latitude in C{radians} converted from C{degrees}.
    '''
    def __new__(cls, arg=None, name=_lat_, Error=UnitError, clip=90, **name_arg):
        '''See L{Degrees} and L{Radians}.
        '''
        if name_arg:
            name, arg = _xkwds_popitem(name_arg)
        d = Phi.__new__(cls, arg=arg, name=name, Error=Error, clip=clip)
        return Radians.__new__(cls, arg=radians(d), name=name, Error=Error)


class Precision_(Int_):
    '''Named C{int} with optional C{low} and C{high} limits representing a precision.
    '''
    def __new__(cls, arg=None, name=_precision_, low=0, **high_Error_name_arg):
        '''See L{Int_}.
        '''
        return Int_.__new__(cls, arg=arg, name=name, low=low, **high_Error_name_arg)


class Radius(Float):
    '''Named C{float} representing a radius, conventionally in C{meter}.
    '''
    def __new__(cls, arg=None, name=_radius_, **Error_name_arg):
        '''See L{Float}.
        '''
        return Float.__new__(cls, arg=arg, name=name, **Error_name_arg)


class Radius_(Float_):
    '''Named C{float} with optional C{low} and C{high} limits representing a radius, conventionally in C{meter}.
    '''
    def __new__(cls, arg=None, name=_radius_, low=EPS, **high_Error_name_arg):
        '''See L{Float}.
        '''
        return Float_.__new__(cls, arg=arg, name=name, low=low, **high_Error_name_arg)


class Scalar(Float):
    '''Named C{float} representing a factor, fraction, scale, etc.
    '''
    def __new__(cls, arg=None, name=_scalar_, **Error_name_arg):
        '''See L{Float}.
        '''
        return Float.__new__(cls, arg=arg, name=name, **Error_name_arg)


class Scalar_(Float_):
    '''Named C{float} with optional C{low} and C{high} limits representing a factor, fraction, scale, etc.
    '''
    def __new__(cls, arg=None, name=_scalar_, low=_0_0, **high_Error_name_arg):
        '''See L{Float_}.
        '''
        return Float_.__new__(cls, arg=arg, name=name, low=low, **high_Error_name_arg)


class Zone(Int):
    '''Named C{int} representing a UTM/UPS zone number.
    '''
    def __new__(cls, arg=None, name=_zone_, **Error_name_arg):
        '''See L{Int}
        '''
        # usually low=_UTMUPS_ZONE_MIN, high=_UTMUPS_ZONE_MAX
        return Int_.__new__(cls, arg=arg, name=name, **Error_name_arg)


def _Error(clas, arg, name=NN, Error=UnitError, txt=_invalid_):
    '''(INTERNAL) Return an error with explanation.

       @arg clas: The C{units} class or sub-class.
       @arg arg: The original C{unit} value.
       @kwarg name: The instance name (C{str}).
       @kwarg Error: Optional error, overriding the default
                     L{UnitError}.
       @kwarg txt: Optional explanation of the error (C{str}).

       @returns: An B{C{Error}} instance.
    '''
    n = name if name else modulename(clas).lstrip(_UNDER_)
    return Error(n, arg, txt=txt)


def _xStrError(*Refs, **name_value_Error):
    '''(INTERNAL) Create a C{TypeError} for C{Garef}, C{Geohash}, C{Wgrs}.
    '''
    r = tuple(r.__name__ for r in Refs) + (Str.__name__, _LatLon_, 'LatLon*Tuple')
    return _IsnotError(*r, **name_value_Error)


def _xUnit(units, Base):  # in .frechet,  .hausdorff
    '''(INTERNAL) Get C{Unit} from C{Unit} or C{name}, ortherwise C{Base}.
    '''
    if not issubclassof(Base, _NamedUnit):
        raise _IsnotError(_NamedUnit.__name__, Base=Base)
    U = globals().get(units.capitalize(), Base) if isstr(units) else (
                      units if issubclassof(units, Base) else Base)
    return U if issubclassof(U, Base) else Base


def _xUnits(units, Base=_NamedUnit):  # in .frechet, .hausdorff
    '''(INTERNAL) Set property C{units} as C{Unit} or C{Str}.
    '''
    if not issubclassof(Base, _NamedUnit):
        raise _IsnotError(_NamedUnit.__name__, Base=Base)
    elif issubclassof(units, Base):
        return units
    elif isstr(units):
        return Str(units, name=_units_)  # XXX Str to _Pass and for backward compatibility
    else:
        raise _IsnotError(Base.__name__, Str.__name__, str.__name__, units=units)


def _std_repr(*classes):
    '''(INTERNAL) Use standard C{repr} or named C{toRepr}.
    '''
    from os import environ as _environ
    for C in classes:
        if hasattr(C, _std_repr.__name__):  # PYCHOK del _std_repr
            env = 'PYGEODESY_%s_STD_REPR' % (C.__name__.upper(),)
            if _environ.get(env, _std_).lower() != _std_:
                C._std_repr = False

_std_repr(Bool, Float, Int, Meter, Str)  # PYCHOK expected
del _std_repr

__all__ += _ALL_DOCS(_NamedUnit)

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
