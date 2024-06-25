
# -*- coding: utf-8 -*-

u'''Various named units, all sub-classes of C{Float}, C{Int} or C{Str} from
basic C{float}, C{int} respectively C{str} to named units as L{Degrees},
L{Feet}, L{Meter}, L{Radians}, etc.
'''

from pygeodesy.basics import isscalar, issubclassof, signOf
from pygeodesy.constants import EPS, EPS1, PI, PI2, PI_2, _umod_360, _0_0, \
                               _0_001,  _0_5, INT0  # PYCHOK for .mgrs, .namedTuples
from pygeodesy.dms import F__F, F__F_, S_NUL, S_SEP, parseDMS, parseRad, _toDMS
from pygeodesy.errors import _AssertionError, TRFError, UnitError, _xattr, _xcallable
from pygeodesy.interns import NN, _band_, _bearing_, _COMMASPACE_, _degrees_, \
                             _degrees2_, _distance_, _E_, _easting_, _epoch_, _EW_, \
                             _feet_, _height_, _lam_, _lat_, _LatLon_, _lon_, \
                             _meter_, _meter2_, _N_, _negative_, _northing_, _NS_, \
                             _NSEW_, _number_, _PERCENT_, _phi_, _precision_, \
                             _radians_, _radians2_, _radius_, _S_, _scalar_, \
                             _units_, _W_, _zone_,  _std_  # PYCHOK used!
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS, _getenv
# from pygeodesy.named import _name__  # _MODS
from pygeodesy.props import Property_RO
# from pygeodesy.streprs import Fmt, fstr  # from .unitsBase
from pygeodesy.unitsBase import Float, Int, _NamedUnit, Radius, Str,  Fmt, fstr

from math import degrees, radians

__all__ = _ALL_LAZY.units
__version__ = '24.06.15'


class Float_(Float):
    '''Named C{float} with optional C{low} and C{high} limit.
    '''
    def __new__(cls, arg=None, name=NN, low=EPS, high=None, **Error_name_arg):
        '''New, named C{Float_}, see L{Float}.

           @arg cls: This class (C{Float_} or sub-class).
           @kwarg arg: The value (any C{type} convertable to C{float}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg low: Optional lower B{C{arg}} limit (C{float} or C{None}).
           @kwarg high: Optional upper B{C{arg}} limit (C{float} or C{None}).

           @returns: A named C{Float_}.

           @raise Error: Invalid B{C{arg}} or B{C{arg}} below B{C{low}} or above B{C{high}}.
        '''
        self = Float.__new__(cls, arg=arg, name=name, **Error_name_arg)
        t = _xlimits(self, low, high, g=True)
        if t:
            raise _NamedUnit._Error(cls, arg, name, txt=t, **Error_name_arg)
        return self


class Int_(Int):
    '''Named C{int} with optional limits C{low} and C{high}.
    '''
    def __new__(cls, arg=None, name=NN, low=0, high=None, **Error_name_arg):
        '''New, named C{int} instance with limits, see L{Int}.

           @kwarg cls: This class (C{Int_} or sub-class).
           @arg arg: The value (any C{type} convertable to C{int}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg low: Optional lower B{C{arg}} limit (C{int} or C{None}).
           @kwarg high: Optional upper B{C{arg}} limit (C{int} or C{None}).

           @returns: A named L{Int_}.

           @raise Error: Invalid B{C{arg}} or B{C{arg}} below B{C{low}} or above B{C{high}}.
        '''
        self = Int.__new__(cls, arg=arg, name=name, **Error_name_arg)
        t = _xlimits(self, low, high)
        if t:
            raise _NamedUnit._Error(cls, arg, name, txt=t, **Error_name_arg)
        return self


class Bool(Int, _NamedUnit):
    '''Named C{bool}, a sub-class of C{int} like Python's C{bool}.
    '''
    # _std_repr = True  # set below
    _bool_True_or_False = None

    def __new__(cls, arg=None, name=NN, Error=UnitError, **name_arg):
        '''New, named C{Bool}.

           @kwarg cls: This class (C{Bool} or sub-class).
           @kwarg arg: The value (any C{type} convertable to C{bool}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg Error: Optional error to raise, overriding the default L{UnitError}.
           @kwarg name_arg: Optional C{name=arg} keyword argument, inlieu of separate
                            B{C{arg}} and B{C{name}} ones.

           @returns: A named L{Bool}, C{bool}-like.

           @raise Error: Invalid B{C{arg}}.
        '''
        if name_arg:
            name, arg = _NamedUnit._arg_name_arg2(arg, **name_arg)
        try:
            b = bool(arg)
        except Exception as x:
            raise _NamedUnit._Error(cls, arg, name, Error, cause=x)

        self = Int.__new__(cls, b, name=name, Error=Error)
        self._bool_True_or_False = b
        return self

    # <https://StackOverflow.com/questions/9787890/assign-class-boolean-value-in-python>
    def __bool__(self):  # PYCHOK Python 3+
        return self._bool_True_or_False

    __nonzero__ = __bool__  # PYCHOK Python 2-

    def toRepr(self, std=False, **unused):  # PYCHOK **unused
        '''Return a representation of this C{Bool}.

           @kwarg std: Use the standard C{repr} or the named representation (C{bool}).

           @note: Use C{env} variable C{PYGEODESY_BOOL_STD_REPR=std} prior to C{import
                  pygeodesy} to get the standard C{repr} or set property C{std_repr=False}
                  to always get the named C{toRepr} representation.
        '''
        r = repr(self._bool_True_or_False)
        return r if std else self._toRepr(r)

    def toStr(self, **unused):  # PYCHOK **unused
        '''Return this C{Bool} as standard C{str}.
        '''
        return str(self._bool_True_or_False)


class Band(Str):
    '''Named C{str} representing a UTM/UPS band letter, unchecked.
    '''
    def __new__(cls, arg=None, name=_band_, **Error_name_arg):
        '''New, named L{Band}, see L{Str}.
        '''
        return Str.__new__(cls, arg=arg, name=name, **Error_name_arg)


class Degrees(Float):
    '''Named C{float} representing a coordinate in C{degrees}, optionally clipped.
    '''
    _ddd_ =  1  # default for .dms._toDMS
    _sep_ =  S_SEP
    _suf_ = (S_NUL,) * 3

    def __new__(cls, arg=None, name=_degrees_, suffix=_NSEW_, clip=0, wrap=None, Error=UnitError, **name_arg):
        '''New C{Degrees} instance, see L{Float}.

           @arg cls: This class (C{Degrees} or sub-class).
           @kwarg arg: The value (any C{type} convertable to C{float} or parsable by
                       function L{parseDMS<pygeodesy.dms.parseDMS>}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg suffix: Optional, valid compass direction suffixes (C{NSEW}).
           @kwarg clip: Optional B{C{arg}} range B{C{-clip..+clip}} (C{degrees} or C{0}
                        or C{None} for unclipped).
           @kwarg wrap: Optionally adjust the B{C{arg}} value (L{wrap90<pygeodesy.wrap90>},
                        L{wrap180<pygeodesy.wrap180>} or L{wrap360<pygeodesy.wrap360>}).
           @kwarg Error: Optional error to raise, overriding the default L{UnitError}.
           @kwarg name_arg: Optional C{name=arg} keyword argument, inlieu of separate
                            B{C{arg}} and B{C{name}} ones.

           @returns: A C{Degrees} instance.

           @raise Error: Invalid B{C{arg}} or B{C{abs(arg)}} outside the B{C{clip}}
                         range and L{rangerrors<pygeodesy.rangerrors>} is C{True}.
        '''
        if name_arg:
            name, arg = _NamedUnit._arg_name_arg2(arg, name, **name_arg)
        try:
            arg = parseDMS(arg, suffix=suffix, clip=clip)
            if wrap:
                _xcallable(wrap=wrap)
                arg = wrap(arg)
            self = Float.__new__(cls, arg=arg, name=name, Error=Error)
        except Exception as x:
            raise _NamedUnit._Error(cls, arg, name, Error, cause=x)
        return self

    def toDegrees(self):
        '''Convert C{Degrees} to C{Degrees}.
        '''
        return self

    def toRadians(self):
        '''Convert C{Degrees} to C{Radians}.
        '''
        return Radians(radians(self), name=self.name)

    def toRepr(self, std=False, **prec_fmt_ints):  # PYCHOK prec=8, ...
        '''Return a representation of this C{Degrees}.

           @kwarg std: If C{True} return the standard C{repr}, otherwise
                       the named representation (C{bool}).

           @see: Methods L{Degrees.toStr}, L{Float.toRepr} and function
                 L{pygeodesy.toDMS} for futher C{prec_fmt_ints} details.
        '''
        return Float.toRepr(self, std=std, **prec_fmt_ints)

    def toStr(self, prec=None, fmt=F__F_, ints=False, **s_D_M_S):  # PYCHOK prec=8, ...
        '''Return this C{Degrees} as standard C{str}.

           @see: Function L{pygeodesy.toDMS} for futher details.
        '''
        if fmt.startswith(_PERCENT_):  # use regular formatting
            p = 8 if prec is None else prec
            return fstr(self, prec=p, fmt=fmt, ints=ints, sep=self._sep_)
        else:
            s = self._suf_[signOf(self) + 1]
            return _toDMS(self, fmt, prec, self._sep_, self._ddd_, s, s_D_M_S)


class Degrees_(Degrees):
    '''Named C{Degrees} representing a coordinate in C{degrees} with optional limits C{low} and C{high}.
    '''
    def __new__(cls, arg=None, name=_degrees_, low=None, high=None, **suffix_Error_name_arg):
        '''New, named C{Degrees_}, see  L{Degrees} and L{Float}.

           @arg cls: This class (C{Degrees_} or sub-class).
           @kwarg arg: The value (any C{type} convertable to C{float} or parsable by
                       function L{parseDMS<pygeodesy.dms.parseDMS>}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg low: Optional lower B{C{arg}} limit (C{float} or C{None}).
           @kwarg high: Optional upper B{C{arg}} limit (C{float} or C{None}).

           @returns: A named C{Degrees}.

           @raise Error: Invalid B{C{arg}} or B{C{arg}} below B{C{low}} or above B{C{high}}.
        '''
        self = Degrees.__new__(cls, arg=arg, name=name, clip=0, **suffix_Error_name_arg)
        t = _xlimits(self, low, high)
        if t:
            raise _NamedUnit._Error(cls, arg, name, txt=t, **suffix_Error_name_arg)
        return self


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
    def __new__(cls, arg=None, name=_radians_, suffix=_NSEW_, clip=0, Error=UnitError, **name_arg):
        '''New, named C{Radians}, see L{Float}.

           @arg cls: This class (C{Radians} or sub-class).
           @kwarg arg: The value (any C{type} convertable to C{float} or parsable by
                       L{pygeodesy.parseRad}).
           @kwarg name: Optional instance name (C{str}).
           @kwarg suffix: Optional, valid compass direction suffixes (C{NSEW}).
           @kwarg clip: Optional B{C{arg}} range B{C{-clip..+clip}} (C{radians} or C{0}
                        or C{None} for unclipped).
           @kwarg Error: Optional error to raise, overriding the default L{UnitError}.
           @kwarg name_arg: Optional C{name=arg} keyword argument, inlieu of separate
                            B{C{arg}} and B{C{name}} ones.

           @returns: A named C{Radians}.

           @raise Error: Invalid B{C{arg}} or B{C{abs(arg)}} outside the B{C{clip}}
                         range and L{rangerrors<pygeodesy.rangerrors>} is C{True}.
        '''
        if name_arg:
            name, arg = _NamedUnit._arg_name_arg2(arg, name, **name_arg)
        try:
            arg = parseRad(arg, suffix=suffix, clip=clip)
            return Float.__new__(cls, arg, name=name, Error=Error)
        except Exception as x:
            raise _NamedUnit._Error(cls, arg, name, Error, cause=x)

    def toDegrees(self):
        '''Convert C{Radians} to C{Degrees}.
        '''
        return Degrees(degrees(self), name=self.name)

    def toRadians(self):
        '''Convert C{Radians} to C{Radians}.
        '''
        return self

    def toRepr(self, std=False, **prec_fmt_ints):  # PYCHOK prec=8, ...
        '''Return a representation of this C{Radians}.

           @kwarg std: If C{True} return the standard C{repr}, otherwise
                       the named representation (C{bool}).

           @see: Methods L{Radians.toStr}, L{Float.toRepr} and function
                 L{pygeodesy.toDMS} for more documentation.
        '''
        return Float.toRepr(self, std=std, **prec_fmt_ints)

    def toStr(self, prec=8, fmt=F__F, ints=False):  # PYCHOK prec=8, ...
        '''Return this C{Radians} as standard C{str}.

           @see: Function L{pygeodesy.fstr} for keyword argument details.
        '''
        return fstr(self, prec=prec, fmt=fmt, ints=ints)


class Radians_(Radians):
    '''Named C{float} representing a coordinate in C{radians} with optional limits C{low} and C{high}.
    '''
    def __new__(cls, arg=None, name=_radians_, low=_0_0, high=PI2, **suffix_Error_name_arg):
        '''New, named C{Radians_}, see L{Radians}.

           @arg cls: This class (C{Radians_} or sub-class).
           @kwarg arg: The value (any C{type} convertable to C{float} or parsable by
                       function L{parseRad<pygeodesy.dms.parseRad>}).
           @kwarg name: Optional name (C{str}).
           @kwarg low: Optional lower B{C{arg}} limit (C{float} or C{None}).
           @kwarg high: Optional upper B{C{arg}} limit (C{float} or C{None}).

           @returns: A named C{Radians_}.

           @raise Error: Invalid B{C{arg}} or B{C{arg}} below B{C{low}} or above B{C{high}}.
        '''
        self = Radians.__new__(cls, arg=arg, name=name, **suffix_Error_name_arg)
        t = _xlimits(self, low, high)
        if t:
            raise _NamedUnit._Error(cls, arg, name, txt=t, **suffix_Error_name_arg)
        return self


class Radians2(Float_):
    '''Named C{float} representing a distance in C{radians squared}.
    '''
    def __new__(cls, arg=None, name=_radians2_, **Error_name_arg):
        '''New, named L{Radians2}, see L{Float_}.
        '''
        return Float_.__new__(cls, arg=arg, name=name, low=_0_0, **Error_name_arg)


class Bearing(Degrees):
    '''Named C{float} representing a bearing in compass C{degrees} from (true) North.
    '''
    _ddd_ =  1
    _suf_ = _N_ * 3  # always suffix N

    def __new__(cls, arg=None, name=_bearing_, clip=0, **Error_name_arg):
        '''New, named L{Bearing}, see L{Degrees}.
        '''
        d =  Degrees.__new__(cls, arg=arg, name=name, suffix=_N_, clip=clip, **Error_name_arg)
        b = _umod_360(d)  # 0 <= b < 360
        return d if b == d else Degrees.__new__(cls, arg=b, name=d.name)


class Bearing_(Radians):
    '''Named C{float} representing a bearing in C{radians} from compass C{degrees} from (true) North.
    '''
    def __new__(cls, arg=None, **name_clip_Error_name_arg):
        '''New, named L{Bearing_}, see L{Bearing} and L{Radians}.
        '''
        d = Bearing.__new__(cls, arg=arg, **name_clip_Error_name_arg)
        return Radians.__new__(cls, radians(d), name=d.name)


class Distance(Float):
    '''Named C{float} representing a distance, conventionally in C{meter}.
    '''
    def __new__(cls, arg=None, name=_distance_, **Error_name_arg):
        '''New, named L{Distance}, see L{Float}.
        '''
        return Float.__new__(cls, arg=arg, name=name, **Error_name_arg)


class Distance_(Float_):
    '''Named C{float} with optional C{low} and C{high} limits representing a distance, conventionally in C{meter}.
    '''
    def __new__(cls, arg=None, name=_distance_, **low_high_Error_name_arg):
        '''New L{Distance_} instance, see L{Float}.
        '''
        return Float_.__new__(cls, arg=arg, name=name, **low_high_Error_name_arg)


class _EasNorBase(Float):
    '''(INTERNAL) L{Easting} and L{Northing} base class.
    '''
    def __new__(cls, arg, name, falsed, high, **Error_name_arg):
        self = Float.__new__(cls, arg=arg, name=name, **Error_name_arg)
        low  = self < 0
        if (high is not None) and (low or self > high):  # like Veness
            t = _negative_ if low else Fmt.limit(above=high)
        elif low and falsed:
            t = _COMMASPACE_(_negative_, 'falsed')
        else:
            return self
        raise _NamedUnit._Error(cls, arg, name, txt=t, **Error_name_arg)


class Easting(_EasNorBase):
    '''Named C{float} representing an easting, conventionally in C{meter}.
    '''
    def __new__(cls, arg=None, name=_easting_, falsed=False, high=None, **Error_name_arg):
        '''New, named C{Easting} or C{Easting of Point}, see L{Float}.

           @arg cls: This class (C{Easting} or sub-class).
           @kwarg arg: The value (any C{type} convertable to C{float}).
           @kwarg name: Optional name (C{str}).
           @kwarg falsed: If C{True}, the B{C{arg}} value is falsed (C{bool}).
           @kwarg high: Optional upper B{C{arg}} limit (C{scalar} or C{None}).

           @returns: A named C{Easting}.

           @raise Error: Invalid B{C{arg}}, above B{C{high}} or negative, falsed B{C{arg}}.
        '''
        return _EasNorBase.__new__(cls, arg, name, falsed, high, **Error_name_arg)


class Epoch(Float_):  # in .ellipsoidalBase, .trf
    '''Named C{epoch} with optional C{low} and C{high} limits representing a fractional
       calendar year.
    '''
    _std_repr = False

    def __new__(cls, arg=None, name=_epoch_, low=1900, high=9999, Error=TRFError, **name_arg):
        '''New, named L{Epoch}, see L{Float_}.
        '''
        if name_arg:
            name, arg = _NamedUnit._arg_name_arg2(arg, name, **name_arg)
        return arg if isinstance(arg, Epoch) else Float_.__new__(cls,
               arg=arg, name=name, Error=Error, low=low, high=high)

    def toRepr(self, prec=3, std=False, **unused):  # PYCHOK fmt=Fmt.F, ints=True
        '''Return a representation of this C{Epoch}.

           @kwarg std: Use the standard C{repr} or the named
                       representation (C{bool}).

           @see: Method L{Float.toRepr} for more documentation.
        '''
        return Float_.toRepr(self, prec=prec, std=std)  # fmt=Fmt.F, ints=True

    def toStr(self, prec=3, **unused):  # PYCHOK fmt=Fmt.F, nts=True
        '''Format this C{Epoch} as C{str}.

           @see: Function L{pygeodesy.fstr} for more documentation.
        '''
        return fstr(self, prec=prec, fmt=Fmt.F, ints=True)

    __str__ = toStr  # PYCHOK default '%.3F', with trailing zeros and decimal point


class Feet(Float):
    '''Named C{float} representing a distance or length in C{feet}.
    '''
    def __new__(cls, arg=None, name=_feet_, **Error_name_arg):
        '''New, named L{Feet}, see L{Float}.
        '''
        return Float.__new__(cls, arg=arg, name=name, **Error_name_arg)


class FIx(Float_):
    '''A named I{Fractional Index}, an C{int} or C{float} index into a C{list}
       or C{tuple} of C{points}, typically.  A C{float} I{Fractional Index}
       C{fi} represents a location on the edge between C{points[int(fi)]} and
       C{points[(int(fi) + 1) % len(points)]}.
    '''
    _fin = 0

    def __new__(cls, fi, fin=None, Error=UnitError, **name):
        '''New, named I{Fractional Index} in a C{list} or C{tuple} of points.

           @arg fi: The fractional index (C{float} or C{int}).
           @kwarg fin: Optional C{len}, the number of C{points}, the index
                       C{[n]} wrapped to C{[0]} (C{int} or C{None}).
           @kwarg Error: Optional error to raise.
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: A named B{C{fi}} (L{FIx}).

           @note: The returned B{C{fi}} may exceed the B{C{len}}, number of
                  original C{points} in certain open/closed cases.

           @see: Method L{fractional} or function L{pygeodesy.fractional}.
        '''
        _ = _MODS.named._name__(name) if name else NN  # check error
        n =  Int_(fin=fin, low=0) if fin else None
        f =  Float_.__new__(cls, fi, low=_0_0, high=n, Error=Error, **name)
        i =  int(f)
        r =  f - float(i)
        if r < EPS:  # see .points._fractional
            f = Float_.__new__(cls, i, low=_0_0, Error=Error, **name)
        elif r > EPS1:
            f = Float_.__new__(cls, i + 1, high=n, Error=Error, **name)
        if n:  # non-zero and non-None
            f._fin = n
        return f

    @Property_RO
    def fin(self):
        '''Get the given C{len}, the index C{[n]} wrapped to C{[0]} (C{int}).
        '''
        return self._fin

    def fractional(self, points, wrap=None, **LatLon_or_Vector_and_kwds):
        '''Return the point at this I{Fractional Index}.

           @arg points: The points (C{LatLon}[], L{Numpy2LatLon}[], L{Tuple2LatLon}[] or
                        C{other}[]).
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the B{C{points}}
                        (C{bool}) or C{None} for backward compatible L{LatLon2Tuple} or
                        B{C{LatLon}} with I{averaged} lat- and longitudes.
           @kwarg LatLon_or_Vector_and_kwds: Optional C{B{LatLon}=None} I{or} C{B{Vector}=None}
                         to return the I{intermediate}, I{fractional} point and optional,
                         additional B{C{LatLon}} I{or} B{C{Vector}} keyword arguments, see
                         function L{fractional<pygeodesy.points.fractional>}.

           @return: See function L{fractional<pygeodesy.points.fractional>}.

           @raise IndexError: In fractional index invalid or B{C{points}} not
                              subscriptable or not closed.

           @raise TypeError: Invalid B{C{LatLon}}, B{C{Vector}} or B{C{kwds}} argument.

           @see: Function L{pygeodesy.fractional}.
        '''
        # fi = 0 if self == self.fin else self
        return _MODS.points.fractional(points, self, wrap=wrap, **LatLon_or_Vector_and_kwds)


def _fi_j2(f, n):  # PYCHOK in .ellipsoidalBaseDI, .vector3d
    # Get 2-tuple (C{fi}, C{j})
    i = int(f)  # like L{FIx}
    if not 0 <= i < n:
        raise _AssertionError(i=i, n=n, f=f, r=f - float(i))
    return FIx(fi=f, fin=n), (i + 1) % n


class Height(Float):  # here to avoid circular import
    '''Named C{float} representing a height, conventionally in C{meter}.
    '''
    def __new__(cls, arg=None, name=_height_, **Error_name_arg):
        '''New, named L{Height}, see L{Float}.
        '''
        return Float.__new__(cls, arg=arg, name=name, **Error_name_arg)


class Height_(Float_):  # here to avoid circular import
    '''Named C{float} with optional C{low} and C{high} limits representing a height, conventionally in C{meter}.
    '''
    def __new__(cls, arg=None, name=_height_, **low_high_Error_name_arg):
        '''New, named L{Height}, see L{Float}.
        '''
        return Float_.__new__(cls, arg=arg, name=name, **low_high_Error_name_arg)


class HeightX(Height):
    '''Like L{Height}, used to distinguish the interpolated height
       from an original L{Height} at a clip intersection.
    '''
    pass


def _heigHt(inst, height):
    '''(INTERNAL) Override the C{inst}ance' height.
    '''
    return inst.height if height is None else Height(height)


class Lam(Radians):
    '''Named C{float} representing a longitude in C{radians}.
    '''
    def __new__(cls, arg=None, name=_lam_, clip=PI, **Error_name_arg):
        '''New, named L{Lam}, see L{Radians}.
        '''
        return Radians.__new__(cls, arg=arg, name=name, suffix=_EW_, clip=clip, **Error_name_arg)


class Lamd(Lam):
    '''Named C{float} representing a longitude in C{radians} converted from C{degrees}.
    '''
    def __new__(cls, arg=None, name=_lon_, clip=180, **Error_name_arg):
        '''New, named L{Lamd}, see L{Lam} and L{Radians}.
        '''
        d = Degrees(arg=arg, name=name, clip=clip, **Error_name_arg)
        return Lam.__new__(cls, radians(d), clip=radians(clip), name=d.name)


class Lat(Degrees):
    '''Named C{float} representing a latitude in C{degrees}.
    '''
    _ddd_ =  2
    _suf_ = _S_, S_NUL, _N_  # no zero suffix

    def __new__(cls, arg=None, name=_lat_, clip=90, **Error_name_arg):
        '''New, named L{Lat}, see L{Degrees}.
        '''
        return Degrees.__new__(cls, arg=arg, name=name, suffix=_NS_, clip=clip, **Error_name_arg)


class Lat_(Degrees_):
    '''Named C{float} representing a latitude in C{degrees} within limits C{low} and C{high}.
    '''
    _ddd_ =  2
    _suf_ = _S_, S_NUL, _N_  # no zero suffix

    def __new__(cls, arg=None, name=_lat_, low=-90, high=90, **Error_name_arg):
        '''See L{Degrees_}.
        '''
        return Degrees_.__new__(cls, arg=arg, name=name, suffix=_NS_, low=low, high=high, **Error_name_arg)


class Lon(Degrees):
    '''Named C{float} representing a longitude in C{degrees}.
    '''
    _ddd_ =  3
    _suf_ = _W_, S_NUL, _E_  # no zero suffix

    def __new__(cls, arg=None, name=_lon_, clip=180, **Error_name_arg):
        '''New, named L{Lon}, see L{Degrees}.
        '''
        return Degrees.__new__(cls, arg=arg, name=name, suffix=_EW_, clip=clip, **Error_name_arg)


class Lon_(Degrees_):
    '''Named C{float} representing a longitude in C{degrees} within limits C{low} and C{high}.
    '''
    _ddd_ =  3
    _suf_ = _W_, S_NUL, _E_  # no zero suffix

    def __new__(cls, arg=None, name=_lon_, low=-180, high=180, **Error_name_arg):
        '''New, named L{Lon_}, see L{Lon} and L{Degrees_}.
        '''
        return Degrees_.__new__(cls, arg=arg, name=name, suffix=_EW_, low=low, high=high, **Error_name_arg)


class Meter(Float):
    '''Named C{float} representing a distance or length in C{meter}.
    '''
    def __new__(cls, arg=None, name=_meter_, **Error_name_arg):
        '''New, named L{Meter}, see L{Float}.
        '''
        return Float.__new__(cls, arg=arg, name=name, **Error_name_arg)

    def __repr__(self):
        '''Return a representation of this C{Meter}.

           @see: Method C{Str.toRepr} and property C{Str.std_repr}.

           @note: Use C{env} variable C{PYGEODESY_METER_STD_REPR=std} prior to C{import
                  pygeodesy} to get the standard C{repr} or set property C{std_repr=False}
                  to always get the named C{toRepr} representation.
        '''
        return self.toRepr(std=self._std_repr)


# _1Å   = Meter(     _Å= 1e-10)  # PyCHOK 1 Ångstrōm
_1um    = Meter(   _1um= 1.e-6)  # PYCHOK 1 micrometer in .mgrs
_10um   = Meter(  _10um= 1.e-5)  # PYCHOK 10 micrometer in .osgr
_1mm    = Meter(   _1mm=_0_001)  # PYCHOK 1 millimeter in .ellipsoidal...
_100km  = Meter( _100km= 1.e+5)  # PYCHOK 100 kilometer in .formy, .mgrs, .osgr, .sphericalBase
_2000km = Meter(_2000km= 2.e+6)  # PYCHOK 2,000 kilometer in .mgrs


class Meter_(Float_):
    '''Named C{float} representing a distance or length in C{meter}.
    '''
    def __new__(cls, arg=None, name=_meter_, low=_0_0, **high_Error_name_arg):
        '''New, named L{Meter_}, see L{Meter} and L{Float_}.
        '''
        return Float_.__new__(cls, arg=arg, name=name, low=low, **high_Error_name_arg)


class Meter2(Float_):
    '''Named C{float} representing an area in C{meter squared}.
    '''
    def __new__(cls, arg=None, name=_meter2_, **Error_name_arg):
        '''New, named L{Meter2}, see L{Float_}.
        '''
        return Float_.__new__(cls, arg=arg, name=name, low=_0_0, **Error_name_arg)


class Meter3(Float_):
    '''Named C{float} representing a volume in C{meter cubed}.
    '''
    def __new__(cls, arg=None, name='meter3', **Error_name_arg):
        '''New, named L{Meter3}, see L{Float_}.
        '''
        return Float_.__new__(cls, arg=arg, name=name, low=_0_0, **Error_name_arg)


class Northing(_EasNorBase):
    '''Named C{float} representing a northing, conventionally in C{meter}.
    '''
    def __new__(cls, arg=None, name=_northing_, falsed=False, high=None, **Error_name_arg):
        '''New, named C{Northing} or C{Northing of point}, see L{Float}.

           @arg cls: This class (C{Northing} or sub-class).
           @kwarg arg: The value (any C{type} convertable to C{float}).
           @kwarg name: Optional name (C{str}).
           @kwarg falsed: If C{True}, the B{C{arg}} value is falsed (C{bool}).
           @kwarg high: Optional upper B{C{arg}} limit (C{scalar} or C{None}).

           @returns: A named C{Northing}.

           @raise Error: Invalid B{C{arg}}, above B{C{high}} or negative, falsed B{C{arg}}.
        '''
        return _EasNorBase.__new__(cls, arg, name, falsed, high, **Error_name_arg)


class Number_(Int_):
    '''Named C{int} representing a non-negative number.
    '''
    def __new__(cls, arg=None, name=_number_, **low_high_Error_name_arg):
        '''New, named L{Number_}, see L{Int_}.
        '''
        return Int_.__new__(cls, arg=arg, name=name, **low_high_Error_name_arg)


class Phi(Radians):
    '''Named C{float} representing a latitude in C{radians}.
    '''
    def __new__(cls, arg=None, name=_phi_, clip=PI_2, **Error_name_arg):
        '''New, named L{Phi}, see L{Radians}.
        '''
        return Radians.__new__(cls, arg=arg, name=name, suffix=_NS_, clip=clip, **Error_name_arg)


class Phid(Phi):
    '''Named C{float} representing a latitude in C{radians} converted from C{degrees}.
    '''
    def __new__(cls, arg=None, name=_lat_, clip=90, **Error_name_arg):
        '''New, named L{Phid}, see L{Phi} and L{Radians}.
        '''
        d = Degrees(arg=arg, name=name, clip=clip, **Error_name_arg)
        return Phi.__new__(cls, arg=radians(d), clip=radians(clip), name=d.name)


class Precision_(Int_):
    '''Named C{int} with optional C{low} and C{high} limits representing a precision.
    '''
    def __new__(cls, arg=None, name=_precision_, **low_high_Error_name_arg):
        '''New, named L{Precision_}, see L{Int_}.
        '''
        return Int_.__new__(cls, arg=arg, name=name, **low_high_Error_name_arg)


class Radius_(Float_):
    '''Named C{float} with optional C{low} and C{high} limits representing a radius, conventionally in C{meter}.
    '''
    def __new__(cls, arg=None, name=_radius_, **low_high_Error_name_arg):
        '''New, named L{Radius_}, see L{Radius} and L{Float}.
        '''
        return Float_.__new__(cls, arg=arg, name=name, **low_high_Error_name_arg)


class Scalar(Float):
    '''Named C{float} representing a factor, fraction, scale, etc.
    '''
    def __new__(cls, arg=None, name=_scalar_, **Error_name_arg):
        '''New, named L{Scalar}, see L{Float}.
        '''
        return Float.__new__(cls, arg=arg, name=name, **Error_name_arg)


class Scalar_(Float_):
    '''Named C{float} with optional C{low} and C{high} limits representing a factor, fraction, scale, etc.
    '''
    def __new__(cls, arg=None, name=_scalar_, low=_0_0, **high_Error_name_arg):
        '''New, named L{Scalar_}, see L{Scalar} and L{Float_}.
        '''
        return Float_.__new__(cls, arg=arg, name=name, low=low, **high_Error_name_arg)


class Zone(Int):
    '''Named C{int} representing a UTM/UPS zone number.
    '''
    def __new__(cls, arg=None, name=_zone_, **Error_name_arg):
        '''New, named L{Zone}, see L{Int}
        '''
        # usually low=_UTMUPS_ZONE_MIN, high=_UTMUPS_ZONE_MAX
        return Int_.__new__(cls, arg=arg, name=name, **Error_name_arg)


_Scalars =  Float, Float_, Scalar, Scalar_
_Degrees = (Bearing, Bearing_, Degrees, Degrees_) + _Scalars
_Meters  = (Distance, Distance_, Meter, Meter_) + _Scalars
_Radians = (Radians, Radians_) + _Scalars  # PYCHOK unused
_Radii   = _Meters + (Radius, Radius_)


def _isDegrees(obj):
    # Check for valid degrees types.
    return isinstance(obj, _Degrees) or _isScalar(obj)


def _isHeight(obj):
    # Check for valid heigth types.
    return isinstance(obj, _Meters) or _isScalar(obj)


def _isMeter(obj):
    # Check for valid meter types.
    return isinstance(obj, _Meters) or _isScalar(obj)


def _isRadius(obj):
    # Check for valid earth radius types.
    return isinstance(obj, _Radii) or _isScalar(obj)


def _isScalar(obj):
    # Check for pure scalar types.
    return isscalar(obj) and not isinstance(obj, _NamedUnit)


def _toUnit(Unit, arg, name=NN, **Error):
    '''(INTERNAL) Wrap C{arg} in a C{name}d C{Unit}.
    '''
    if not (issubclassof(Unit, _NamedUnit) and isinstance(arg, Unit) and
                                              _xattr(arg, name=NN) == name):
        arg = Unit(arg, name=name, **Error)
    return arg


def _xlimits(arg, low, high, g=False):
    '''(INTERNAL) Check C{low <= unit <= high}.
    '''
    if (low is not None) and arg < low:
        if g:
            low = Fmt.g(low, prec=6, ints=isinstance(arg, Epoch))
        t = Fmt.limit(below=low)
    elif (high is not None) and arg > high:
        if g:
            high = Fmt.g(high, prec=6, ints=isinstance(arg, Epoch))
        t = Fmt.limit(above=high)
    else:
        t = NN
    return t


def _std_repr(*Classes):
    '''(INTERNAL) Use standard C{repr} or named C{toRepr}.
    '''
    for C in Classes:
        if hasattr(C, _std_repr.__name__):  # PYCHOK del _std_repr
            env = 'PYGEODESY_%s_STD_REPR' % (C.__name__.upper(),)
            if _getenv(env, _std_).lower() != _std_:
                C._std_repr = False

_std_repr(Bearing, Bool, Degrees, Float, Int, Meter, Radians, Str)  # PYCHOK expected
del _std_repr

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
