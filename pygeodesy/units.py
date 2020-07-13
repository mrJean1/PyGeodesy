
# -*- coding: utf-8 -*-

u'''Sub-classes C{Float}, C{Int} and C{Str} from basic
C{float}, C{int} respectively C{str} to named units as
L{Degrees}, L{Feet}, L{Meter}, L{Radians}, etc.

@newfield example: Example, Examples
'''

from pygeodesy.basics import EPS, PI, PI_2, property_doc_
from pygeodesy.dms import F__F, F__F_, parseDMS, parseRad, \
                          S_NUL, S_SEP, _toDMS
from pygeodesy.errors import _IsnotError, RangeError, _ValueError
from pygeodesy.interns import _band_, _bearing_, _degrees_, _E_, _easting_, \
                              _EW_, _feet_, _height_, _invalid_, _lam_, \
                              _lat_, _LatLon_, _lon_, _meter_, _N_, NN, _northing_, \
                              _NS_, _NSEW_, _number_, _PERCENT_, _phi_, \
                              _precision_, _radians_, _radius_, _S_, _scalar_, \
                              _SPACE_, _std_, _UNDERSCORE_, _W_, _zone_  # PYCHOK used!
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import modulename, _Named
from pygeodesy.streprs import fstr, _g

from math import radians

__all__ = _ALL_LAZY.units
__version__ = '20.07.12'


class UnitError(_ValueError):
    '''Default exception for L{units} issues.
    '''
    pass


class _NamedUnit(_Named):
    '''(INTERNAL) Base class for C{units}.
    '''
    _std_repr = True  # set below
    _units    = None

    def _toRepr(self, value):
        '''(INTERNAL) Representation "<name> (<value>)" or "<classname>(<value>)".
        '''
        t = (self.name, _SPACE_) if self.name else (self.classname,)
        return ''.join(t + ('(', str(value), ')'))

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
            self._units = self.classname.lower()
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

    def __new__(cls, arg, name=NN, Error=UnitError):
        '''New named C{float} instance.

           @arg arg: The value (any C{type} convertable to C{float}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the
                         default L{UnitError}.

           @returns: A C{Float} instance.

           @raise Error: Invalid B{C{arg}}.
        '''
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

           @note: Use env variable {PYGEODESY_FLOAT_STD_REPR=std}
                  to get the standard C{repr} or use C{-=named}
                  for the named C{toRepr}.
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

    def toRepr(self, prec=12, fmt=_g, ints=False, std=False):  # PYCHOK prec=8, ...
        '''Return a representation of this named C{float}.

           @kwarg std: Use the standard C{repr} or the named
                       representation (C{bool}).

           @see: Function L{fstr} for more documentation.
        '''
        # XXX must use super(Float, self)... since
        # super()... only works for Python 3+
        return super(Float, self).__repr__() if std else \
               self._toRepr(fstr(self, prec=prec, fmt=fmt, ints=ints))

    def toStr(self, prec=12, fmt=_g, ints=False):  # PYCHOK prec=8, ...
        '''Format this C{Float} as C{str}.

           @see: Function L{fstr} for more documentation.
        '''
        return fstr(self, prec=prec, fmt=fmt, ints=ints)


class Float_(Float):
    '''Named C{float} with optional C{low} and C{high} limit.
    '''
    def __new__(cls, arg, name=NN, Error=UnitError, low=EPS, high=None):
        '''New named C{float} instance with limits.

           @arg cls: This class (C{Float_} or sub-class).
           @arg arg: The value (any C{type} convertable to C{float}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the default
                         L{UnitError}.
           @kwarg low: Optional lower B{C{arg}} limit (C{float} or C{None}).
           @kwarg high: Optional upper B{C{arg}} limit (C{float} or C{None}).

           @returns: A C{Float_} instance.

           @raise Error: Invalid B{C{arg}} or B{C{arg}} below B{C{low}} or
                         above B{C{high}}.
        '''
        self = Float.__new__(cls, arg, name=name, Error=Error)
        if (low is not None) and self < low:
            raise _Error(cls, arg, name=name, Error=Error, txt='below %.6G limit' % (low,))
        if (high is not None) and self > high:
            raise _Error(cls, arg, name=name, Error=Error, txt='above %.6G limit' % (high,))
        return self


class Int(int, _NamedUnit):
    '''Named C{int}.
    '''
    # _std_repr = True  # set below

    def __new__(cls, arg, name=NN, Error=UnitError):
        '''New named {int} instance.

           @arg arg: The value (any C{type} convertable to C{int}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the
                         default L{UnitError}.

           @returns: An C{Int} instance.

           @raise Error: Invalid B{C{arg}}.
        '''
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

           @note: Use env variable {PYGEODESY_INT_STD_REPR=std}
                  to get the standard C{repr} or use C{...=named}
                  for the named C{toRepr}.
        '''
        return self.toRepr(std=self._std_repr)

    def __str__(self):  # to avoid MRO(int)
        '''Return this C{Int} as standard C{str}.
        '''
        return self.toStr()

    def toRepr(self, std=False, **unused):  # PYCHOK **unused
        '''Return the representation of this named C{int}.

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
    def __new__(cls, arg, name=NN, Error=UnitError, low=0, high=None):
        '''New named C{int} instance with limits.

           @arg cls: This class (C{Int_} or sub-class).
           @arg arg: The value (any C{type} convertable to C{int}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the default
                         C{UnitError}.
           @kwarg low: Optional lower B{C{arg}} limit (C{float} or C{None}).
           @kwarg high: Optional upper B{C{arg}} limit (C{float} or C{None}).

           @returns: An L{Int_} instance.

           @raise Error: Invalid B{C{arg}} or B{C{arg}} below B{C{low}} or
                         above B{C{high}}.
        '''
        self = Int.__new__(cls, arg, name=name, Error=Error)
        if (low is not None) and self < low:
            raise _Error(cls, arg, name=name, Error=Error, txt='below %s limit' % (low,))
        if (high is not None) and self > high:
            raise _Error(cls, arg, name=name, Error=Error, txt='above %s limit' % (high,))
        return self


class Str(str, _NamedUnit):
    '''Named C{str}.
    '''
    # _std_repr = True  # set below

    def __new__(cls, arg, name=NN, Error=UnitError):
        '''New named C{str} instance.

           @arg arg: The value (any C{type} convertable to C{str}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding
                         the default (C{ValueError}).

           @returns: A L{Str} instance.

           @raise Error: Invalid B{C{arg}}.
        '''
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

           @note: Use env variable {PYGEODESY_STR_STD_REPR=std}
                  to get the standard C{repr} or use C{...=named}
                  for the named C{toRepr}.
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


class Band(Str):
    '''Named C{str} representing a UTM/UPS band latter, unchecked.
    '''
    def __new__(cls, arg, name=_band_, Error=UnitError):
        '''See L{Str}.
        '''
        return Str.__new__(cls, arg, name=name, Error=Error)


class Degrees(Float):
    '''Named C{float} representing a coordinate in C{degrees}, optionally clipped.
    '''
    _ddd_ = 1  # default for .dms._toDMS
    _suf_ = S_NUL, S_NUL, S_NUL
    _sep_ = S_SEP

    def __new__(cls, arg, name=_degrees_, Error=UnitError, suffix=_NSEW_, clip=0):
        '''New named C{Degrees} instance.

           @arg cls: This class (C{Degrees} or sub-class).
           @arg arg: The value (any C{type} convertable to C{float} or
                     parsable by L{parseDMS}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the default
                         L{UnitError}.
           @kwarg suffix: Optional, valid compass direction suffixes (C{NSEW}).
           @kwarg clip: Optional B{C{arg}} range B{C{-clip..+clip}}
                        (C{degrees} or C{0} or C{None} for unclipped).

           @returns: A C{Degrees} instance.

           @raise Error: Invalid B{C{arg}} or B{C{abs(deg)}} outside the
                         B{C{clip}} range and L{rangerrors} set to C{True}.
        '''
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


class Bearing(Degrees):
    '''Named C{float} representing a bearing in compass C{degrees}.
    '''
    _ddd_ =  1
    _suf_ = _N_ * 3  # always suffix N

    def __new__(cls, arg, name=_bearing_, Error=UnitError, clip=0):
        '''See L{Degrees}.
        '''
        d = Degrees.__new__(cls, arg, name=name, Error=Error, suffix=_N_, clip=clip)
        b = d % 360
        return d if b == d else Degrees.__new__(cls, b, name=name)


class Radians(Float):
    '''Named C{float} representing a coordinate in C{radians}, optionally clipped.
    '''
    def __new__(cls, arg, name=_radians_, Error=UnitError, suffix=_NSEW_, clip=0):
        '''New named C{Radians} instance.

           @arg cls: This class (C{Radians} or sub-class).
           @arg arg: The value (any C{type} convertable to C{float} or
                     parsable by L{parseRad}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the default
                         L{UnitError}.
           @kwarg suffix: Optional, valid compass direction suffixes (C{NSEW}).
           @kwarg clip: Optional B{C{arg}} range B{C{-clip..+clip}}
                        (C{radians} or C{0} or C{None} for unclipped).

           @returns: A C{Radians} instance.

           @raise Error: Invalid B{C{arg}} or B{C{abs(deg)}} outside the
                         B{C{clip}} range and L{rangerrors} set to C{True}.
        '''
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


class Bearing_(Radians):
    '''Named C{float} representing a bearing in C{radians} from compass C{degrees}.
    '''
    def __new__(cls, arg, name=_bearing_, Error=UnitError, clip=0):
        '''See L{Bearing} and L{Radians}.
        '''
        d = Bearing.__new__(cls, arg, name=name, Error=Error, clip=clip)
        return Radians.__new__(cls, radians(d), name=name)


class Distance(Float):
    '''Named C{float} representing a distance, conventionally in C{meter}.
    '''
    def __new__(cls, arg, name='distance', Error=UnitError):
        '''See L{Float}.
        '''
        return Float.__new__(cls, arg, name=name, Error=Error)


class Easting(Float):
    '''Named C{float} representing an easting.
    '''
    def __new__(cls, arg, name=_easting_, Error=UnitError, falsed=False, osgr=False):
        '''New named C{Easting} instance.

           @arg cls: This class (C{Easting} or sub-class).
           @arg arg: The value (any C{type} convertable to C{float}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the default
                         L{UnitError}.
           @kwarg falsed: The B{C{arg}} value includes false origin (C{bool}).
           @kwarg osgr: Check B{C{arg}} as an OSGR easting (C{bool}).

           @returns: An C{Easting} instance.

           @raise Error: Invalid B{C{arg}} or negative, falsed B{C{arg}}.
        '''
        self = Float.__new__(cls, arg, name=name, Error=Error)
        if osgr and (self < 0 or self > 700e3):  # like Veness
            raise _Error(cls, arg, name=name, Error=Error)
        elif falsed and self < 0:
            raise _Error(cls, arg, name=name, Error=Error, txt='negative, falsed')
        return self


class Feet(Float):
    '''Named C{float} representing a distance or length in C{feet}.
    '''
    def __new__(cls, arg, name=_feet_, Error=UnitError):
        '''See L{Float}.
        '''
        return Float.__new__(cls, arg, name=name, Error=Error)


class Height(Float):  # here to avoid circular import
    '''Named C{float} representing a height, conventionally in C{meter}.
    '''
    def __new__(cls, arg, name=_height_, Error=UnitError):
        '''See L{Float}.
        '''
        return Float.__new__(cls, arg, name=name, Error=Error)


class Lam(Radians):
    '''Named C{float} representing a longitude in C{radians}.
    '''
    def __new__(cls, arg, name=_lam_, Error=UnitError, clip=PI):
        '''See L{Radians}.
        '''
        return Radians.__new__(cls, arg, name=name, Error=Error, suffix=_EW_, clip=clip)


class Lam_(Lam):
    '''Named C{float} representing a longitude in C{radians} converted from C{degrees}.
    '''
    def __new__(cls, arg, name=_lon_, Error=UnitError, clip=180):
        '''See L{Degrees} and L{Radians}.
        '''
        d = Lam.__new__(cls, arg, name=name, Error=Error, clip=clip)
        return Radians.__new__(cls, radians(d), name=name)


class Lat(Degrees):
    '''Named C{float} representing a latitude in C{degrees}.
    '''
    _ddd_ =  2
    _suf_ = _N_, _S_, S_NUL  # no zero suffix

    def __new__(cls, arg, name=_lat_, Error=UnitError, clip=90):
        '''See L{Degrees}.
        '''
        return Degrees.__new__(cls, arg, name=name, Error=Error, suffix=_NS_, clip=clip)


class Lon(Float):
    '''Named C{float} representing a longitude in C{degrees}.
    '''
    _ddd_ =  3
    _suf_ = _E_, _W_, S_NUL  # no zero suffix

    def __new__(cls, arg, name=_lon_, Error=UnitError, clip=180):
        '''See L{Degrees}.
        '''
        return Degrees.__new__(cls, arg, name=name, Error=Error, suffix=_EW_, clip=clip)


class Meter(Float):
    '''Named C{float} representing a distance or length in C{meter}.
    '''
    def __new__(cls, arg, name=_meter_, Error=UnitError):
        '''See L{Float}.
        '''
        return Float.__new__(cls, arg, name=name, Error=Error)


class Northing(Float):
    '''Named C{float} representing a northing.
    '''
    def __new__(cls, arg, name=_northing_, Error=UnitError, falsed=False, osgr=False):
        '''New named C{Northing} instance.

           @arg cls: This class (C{Northing} or sub-class).
           @arg arg: The value (any C{type} convertable to C{float}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the default
                         L{UnitError}.
           @kwarg falsed: The B{C{arg}} value includes false origin (C{bool}).
           @kwarg osgr: Check B{C{arg}} as an OSGR northing (C{bool}).

           @returns: A C{Northing} instance.

           @raise Error: Invalid B{C{arg}} or negative, falsed B{C{arg}}.
        '''
        self = Float.__new__(cls, arg, name=name, Error=Error)
        if osgr and (self < 0 or self > 1300e3):  # like Veness
            raise _Error(cls, arg, name=name, Error=Error)
        elif falsed and self < 0:
            raise _Error(cls, arg, name=name, Error=Error, txt='negative, falsed')
        return self


class Number_(Int_):
    '''Named C{int} with optional limits C{low} and C{high} representing a non-negtive number.
    '''
    def __new__(cls, arg, name=_number_, Error=UnitError, low=0, high=None):
        '''See L{Int_}.
        '''
        return Int_.__new__(cls, arg, name=name, Error=Error, low=low, high=high)


class Phi(Radians):
    '''Named C{float} representing a latitude in C{radians}.
    '''
    def __new__(cls, arg, name=_phi_, Error=UnitError, clip=PI_2):
        '''See L{Radians}.
        '''
        return Radians.__new__(cls, arg, name=name, Error=Error, suffix=_NS_, clip=clip)


class Phi_(Phi):
    '''Named C{float} representing a latitude in C{radians} converted from C{degrees}.
    '''
    def __new__(cls, arg, name=_lat_, Error=UnitError, clip=90):
        '''See L{Degrees} and L{Radians}.
        '''
        d = Phi.__new__(cls, arg, name=name, Error=Error, clip=clip)
        return Radians.__new__(cls, radians(d), name=name)


class Precision_(Int_):
    '''Named C{int} with optional C{low} and C{high} limits representing a precision.
    '''
    def __new__(cls, arg, name=_precision_, Error=UnitError, low=0, high=None):
        '''See L{Int_}.
        '''
        return Int_.__new__(cls, arg, name=name, Error=Error, low=low, high=high)


class Radius(Float):
    '''Named C{float} representing a radius, conventionally in C{meter}.
    '''
    def __new__(cls, arg, name=_radius_, Error=UnitError):
        '''See L{Float}.
        '''
        return Float.__new__(cls, arg, name=name, Error=Error)


class Radius_(Float_):
    '''Named C{float} with optional C{low} and C{high} limits representing a radius, conventionally in C{meter}.
    '''
    def __new__(cls, arg, name=_radius_, Error=UnitError, low=EPS, high=None):
        '''See L{Float}.
        '''
        return Float_.__new__(cls, arg, name=name, Error=Error, low=low, high=high)


class Scalar(Float):
    '''Named C{float} representing a factor, fraction, scale, etc.
    '''
    def __new__(cls, arg, name=_scalar_, Error=UnitError):
        '''See L{Float}.
        '''
        return Float.__new__(cls, arg, name=name, Error=Error)


class Scalar_(Float_):
    '''Named C{float} with optional C{low} and C{high} limits representing a factor, fraction, scale, etc.
    '''
    def __new__(cls, arg, name=_scalar_, Error=UnitError, low=0.0, high=None):
        '''See L{Float_}.
        '''
        return Float_.__new__(cls, arg, name=name, Error=Error, low=low, high=high)


class Zone(Int):
    '''Named C{int} representing a UTM/UPS zone number.
    '''
    def __new__(cls, arg, name=_zone_, Error=UnitError):
        '''See L{Int}
        '''
        # usually low=_UTMUPS_ZONE_MIN, high=_UTMUPS_ZONE_MAX
        return Int_.__new__(cls, arg, name=name, Error=Error)


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
    n = name if name else modulename(clas).lstrip(_UNDERSCORE_)
    return Error(n, arg, txt=txt)


def _xStrError(*Refs, **name_value_Error):
    '''(INTERNAL) Create a C{TypeError} for C{Garef}, C{Geohash}, C{Wgrs}.
    '''
    r = tuple(r.__name__ for r in Refs) + (Str.__name__, _LatLon_, 'LatLon*Tuple')
    return _IsnotError(*r, **name_value_Error)


def _std_repr(*classes):
    '''(INTERNAL) Use standard C{repr} or named C{toRepr}.
    '''
    from os import environ as _environ
    for C in classes:
        if hasattr(C, _std_repr.__name__):  # PYCHOK del _std_repr
            env = 'PYGEODESY_%s_STD_REPR' % (C.__name__.upper(),)
            if _environ.get(env, _std_).lower() != _std_:
                C._std_repr = False

_std_repr(Float, Int, Str)  # PYCHOK expected
del _std_repr

__all__ += _ALL_DOCS(_NamedUnit)

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
