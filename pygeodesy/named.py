
# -*- coding: utf-8 -*-

u'''(INTERNAL) Classes C{_Named}, C{_NamedDict} and C{_NamedTuple},
all with nameable instances and several subclasses thereof.

In addition, the items in a C{_NamedDict} are accessable as attributes
and the items in a C{_NamedTuple} can be named to be accessable as
attributes, similar to standard Python C{namedtuple}s.

Results previously returned as tuples by C{pygeodesy} functions and
class methods are now instances of some C{...Tuple} class, all
sub-classes of C{_NamedTuple} defined here.

@newfield example: Example, Examples
'''

# update imported names under if __name__ == '__main__':
from pygeodesy.basics import isstr, issubclassof, property_doc_, \
                             property_RO, _xcopy, _xinstanceof, _xkwds
from pygeodesy.errors import _AssertionError, _AttributeError, _incompatible, \
                             _IndexError, _IsnotError, LenError, _NameError, \
                             _NotImplementedError, \
                             _TypeError, _TypesError, _ValueError  # PYCHOK indent
from pygeodesy.interns import _angle_, _AT_, _COLON_, _COLON_SPACE_, _COMMA_SPACE_, \
                              _CURLY_, _datum_, _distance_, _doesn_t_exist_, _DOT_, \
                              _dot_, _DUNDER_, _dunder_name, _easting_, _EQUAL_, _h_, \
                              _height_, _item_ps, _item_sq, _lam_, _lat_, _lon_, \
                              _name_, NN, _northing_, _number_, _other_, _PARENTH_, \
                              _phi_, _points_, _precision_, _radius_, _UNDERSCORE_, \
                              _valid_, _x_, _y_, _z_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.streprs import attrs, _Fmt, pairs, reprs, unstr

__all__ = _ALL_LAZY.named
__version__ = '20.08.04'

# __DUNDER gets mangled in class
_final_     = 'final'
_immutable_ = 'immutable'
_initial_   = 'initial'
_name       = '_name'
_Names_     = '_Names_'


def _xnamed(inst, name):
    '''(INTERNAL) Set the instance' C{.name = }B{C{name}}.

       @arg inst: The instance (C{_Named}).
       @arg name: The name (C{str}).

       @return: The B{C{inst}}, named if not named before.
    '''
    if name and isinstance(inst, _Named) and not inst.name:
        inst.name = name
    return inst


class _Named(object):
    '''(INTERNAL) Root class for named objects.
    '''
    _name        = NN     #: (INTERNAL) name (C{str})
    _classnaming = False  #: (INTERNAL) prefixed (C{bool})

    def __repr__(self):
        '''Default C{repr(self)}.
        '''
        return '<%s at %#x>' % (self, id(self))

    def __str__(self):
        '''Default C{str(self)}.
        '''
        return self.named2

    def attrs(self, *names, **kwds):
        '''Join attributes as C{name=value} pairs.

           @arg names: The attribute names (C{str}s).
           @kwarg kwds: Keyword argument for function L{attrs}.

           @return: All C{name=value} pairs joined (C{str}).
        '''
        return _COMMA_SPACE_.join(attrs(self, *names, **kwds))

    @property_RO
    def classname(self):
        '''Get this object's C{[module.]class} name (C{str}), see
           property C{.classnaming} and function C{classnaming}.
        '''
        return classname(self, prefixed=self._classnaming)

    @property_doc_(''' the class naming (C{bool}).''')
    def classnaming(self):
        '''Get the class naming (C{bool}), see function C{classnaming}.
        '''
        return self._classnaming

    @classnaming.setter  # PYCHOK setter!
    def classnaming(self, prefixed):
        '''Set the class naming for C{[module.].class} names.

           @arg prefixed: Include the module name (C{bool}).
        '''
        self._classnaming = bool(prefixed)

    def classof(self, *args, **kwds):
        '''Create another instance of this very class.

           @arg args: Optional, positional arguments.
           @kwarg kwds: Optional, keyword arguments.

           @return: New instance (B{self.__class__}).
        '''
        return _xnamed(self.__class__(*args, **kwds), self.name)

    def copy(self, deep=False):
        '''Make a shallow or deep copy of this instance.

           @kwarg deep: If C{True} make a deep, otherwise
                          a shallow copy (C{bool}).

           @return: The copy (C{This class} or subclass thereof).
        '''
        return _xcopy(self, deep=deep)

    def _dot_(self, name):
        '''(INTERNAL) Period-join C{self.name} and C{name}.
        '''
        return _dot_(self.name, name)

    @property_doc_(''' the name (C{str}).''')
    def name(self):
        '''Get the name (C{str}).
        '''
        return self._name

    @name.setter  # PYCHOK setter!
    def name(self, name):
        '''Set the name.

           @arg name: New name (C{str}).
        '''
        self._name = str(name)
        # to set the name from a sub-class, use
        # self.name = name or
        # _Named.name.fset(self, name), but not
        # _Named(self).name = name

    @property_RO
    def named(self):
        '''Get the name I{or} class name or C{""} (C{str}).
        '''
        return self.name or self.classname

    @property_RO
    def named2(self):
        '''Get the class name I{and/or} the name or C{""} (C{str}).
        '''
        n, c = self.name, self.classname
        return ('%s %r' % (c, n)) if c and n else (c or n)

    def toRepr(self, **unused):  # PYCHOK expected
        '''Default C{repr(self)}.
        '''
        return repr(self)

    toStr2 = toRepr  # PYCHOK for backward compatibility
    '''DEPRECATED, used method C{toRepr}.'''

    def toStr(self, **unused):  # PYCHOK expected
        '''Default C{str(self)}.
        '''
        return str(self)

    def _xnamed(self, inst, name=NN):
        '''(INTERNAL) Set the instance' C{.name = self.name}.

           @arg inst: The instance (C{_Named}).
           @kwarg name: Optional name, overriding C{self.name} (C{str}).

           @return: The B{C{inst}}, named if not named before.
        '''
        n = name or self.name
        return _xnamed(inst, n) if n else inst

    def _xrenamed(self, inst):
        '''(INTERNAL) Rename the instance' C{.name = self.name}.

           @arg inst: The instance (C{_Named}).

           @return: The B{C{inst}}, named if not named before.
        '''
        if not isinstance(inst, _Named):
            raise _IsnotError(_valid_, inst=inst)

        if inst.name != self.name:
            inst.name = self.name
        return inst


class _NamedBase(_Named):
    '''(INTERNAL) Base class with name.
    '''

    def __repr__(self):
        return self.toRepr()

    def __str__(self):
        return self.toStr()

    def _update(self, updated, *attrs):
        '''(INTERNAL) Zap cached instance attributes.
        '''
        if updated and attrs:
            for a in attrs:  # zap attrs to None
                if getattr(self, a, None) is not None:
                    setattr(self, a, None)
                elif not hasattr(self, a):
                    raise _AssertionError('.%s invalid: %r' % (a, self))

#   def notImplemented(self, attr):
#       '''Raise error for a missing method, function or attribute.
#
#          @arg attr: Attribute name (C{str}).
#
#          @raise NotImplementedError: No such attribute.
#       '''
#       c = self.__class__.__name__
#       return NotImplementedError(_dot_(c, attr))

    def others(self, other, name=_other_, up=1):  # see .points.LatLon_.others
        '''Check this and an other instance for type compatibility.

           @arg other: The other instance (any C{type}).
           @kwarg name: Optional, name for other (C{str}).
           @kwarg up: Number of call stack frames up (C{int}).

           @return: The B{C{other}} if compatible.

           @raise TypeError: Mismatch of the B{C{other}} and this
                             C{class} or C{type}.
        '''
        if not (isinstance(self, other.__class__) or
                isinstance(other, self.__class__)):
            n = _callname(name, classname(self, prefixed=True), self.name, up=up + 1)
            raise _TypeError(name, other, txt=_incompatible(n))
        return other

    def toRepr(self, **kwds):  # PYCHOK expected
        '''(INTERNAL) I{Could be overloaded}.

           @kwarg kwds: Optional, keyword arguments.

           @return: C{toStr}() with keyword arguments (as C{str}).
        '''
        t = self.toStr(**kwds).lstrip('([{').rstrip('}])')
        return _item_ps(self.classname, t)  # XXX (self.named, t)

#   def toRepr(self, **kwds)
#       if kwds:
#           s = NN.join(reprs((self,), **kwds))
#       else:  # super().__repr__ only for Python 3+
#           s = super(self.__class__, self).__repr__()
#       return _item_ps(self.named, s)  # clips(s)

    def toStr(self, **kwds):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}.

           @raise _AssertionError: Always, see function L{notOverloaded}.
        '''
        notOverloaded(self, self.toStr, **kwds)

#   def toStr(self, **kwds):
#       if kwds:
#           s = NN.join(strs((self,), **kwds))
#       else:  # super().__str__ only for Python 3+
#           s = super(self.__class__, self).__str__()
#       return s


class _NamedDict(dict, _Named):
    '''(INTERNAL) Named C{dict} with key I{and} attribute
       access to the items.
    '''

    def __init__(self, *args, **kwds):
        if args:  # args override kwds
            if len(args) != 1:
                t = unstr(self.classname, *args, **kwds)
                raise _ValueError(args=len(args), txt=t)
            kwds = _xkwds(dict(args[0]), **kwds)
        if _name_ in kwds:
            _Named.name.fset(self, kwds.pop(_name_))  # see _Named.name
        dict.__init__(self, kwds)

    def __delattr__(self, name):
        if name in dict.keys(self):
            dict.pop(name)
        elif name in (_name_, _name):
            dict.__setattr__(self, name, NN)
        else:
            dict.__delattr__(self, name)

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            if name == _name_:
                return _Named.name.fget(self)
            return dict.__getattr__(self, name)

    def __getitem__(self, key):
        if key == _name_:
            raise KeyError(_item_sq(self.classname, key))
        return dict.__getitem__(self, key)

    def __repr__(self):
        return self.toRepr()

    def __setattr__(self, name, value):
        if name in dict.keys(self):
            dict.__setitem__(self, name, value)  # self[name] = value
        else:
            dict.__setattr__(self, name, value)

    def __setitem__(self, key, value):
        if key == _name_:
            raise KeyError('%s = %r' % (_item_sq(self.classname, key), value))
        dict.__setitem__(self, key, value)

    def __str__(self):
        return self.toStr()

    def toRepr(self, prec=6, fmt=_Fmt):  # PYCHOK _Named
        '''Like C{repr(dict)} but with C{name} and  C{floats} formatting by C{fstr}.
        '''
        t = pairs(self.items(), prec=prec, fmt=fmt, sep=_EQUAL_)
        return _item_ps(self.name, _COMMA_SPACE_.join(sorted(t)))

    toStr2 = toRepr  # PYCHOK for backward compatibility
    '''DEPRECATED, use method C{toRepr}.'''

    def toStr(self, prec=6, fmt=_Fmt):  # PYCHOK _Named
        '''Like C{str(dict)} but with C{floats} formatting by C{fstr}.
        '''
        t = pairs(self.items(), prec=prec, fmt=fmt, sep=_COLON_SPACE_)
        return _CURLY_ % (_COMMA_SPACE_.join(sorted(t)),)


class _NamedEnum(_NamedDict):
    '''(INTERNAL) Enum-like C{_NamedDict} with attribute access
       restricted to valid keys.
    '''
    _item_Classes = ()

    def __init__(self, name, *Classes):
        '''New C{_NamedEnum}.

           @arg name: Name (C{str}).
           @arg Classes: One or more classes or types acceptable
                         as enum values (C{class}- or C{type}s).
        '''
        if Classes:
            self._item_Classes = Classes
        if name:
            _Named.name.fset(self, name)  # see _Named.name

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            if name == _name_:
                return _NamedDict.name.fget(self)
        raise _AttributeError(item=self._dot_(name), txt=_doesn_t_exist_)

    def __repr__(self):
        return self.toRepr()

    def __str__(self):
        return self.toStr()

    def _assert(self, **kwds):
        '''(INTERNAL) Check names against given, registered names.
        '''
        for a, v in kwds.items():
            assert self[a] is v and getattr(self, a) \
                                and self.find(v) == a

    def find(self, item):
        '''Find a registered item.

           @arg item: The item to look for (any C{type}).

           @return: If found the B{C{item}}'s name (C{str}), C{None} otherwise.
        '''
        for k, v in self.items():
            if v is item:
                return k
        return None

    def register(self, item):
        '''Registed a new item.

           @arg item: The item (any C{type}).

           @return: The item name (C{str}).

           @raise NameError: An B{C{item}} already registered with
                             that name or the B{C{item}} has no, an
                             empty or an invalid name.

           @raise TypeError: The B{C{item}} type invalid.
        '''
        try:
            n = item.name
            if not (n and n.replace(_UNDERSCORE_, NN).isalnum() and isstr(n)):
                raise ValueError
        except (AttributeError, ValueError, TypeError) as x:
            raise _NameError(_dot_('item', _name_), item, txt=str(x))
        if n in self:
            raise _NameError(self._dot_(n), item, txt='exists')
        if not (self._item_Classes and isinstance(item, self._item_Classes)):
            raise _TypesError(self._dot_(n), item, *self._item_Classes)
        self[n] = item

    def unregister(self, name_or_item):
        '''Remove a registered item.

           @arg name_or_item: Name (C{str}) of or the item (any C{type}).

           @return: The unregistered item.

           @raise NameError: No item with that B{C{name}}.

           @raise ValueError: No such item.
        '''
        name = self.find(name_or_item)
        if name is None:
            if not isstr(name_or_item):
                raise _ValueError(name_or_item=name_or_item)
            name = name_or_item
        try:
            item = dict.pop(self, name)
        except KeyError:
            raise _NameError(item=self._dot_(name), txt=_doesn_t_exist_)
        item._enum = None
        return item

    def toRepr(self, prec=6, fmt=_Fmt, sep=',\n'):  # PYCHOK _NamedDict
        '''Like C{repr(dict)} but with C{name} and C{floats} formatting by C{fstr}.
        '''
        t = sorted((self._dot_(n), v) for n, v in self.items())
        return sep.join(pairs(t, prec=prec, fmt=fmt, sep=_COLON_SPACE_))

    toStr2 = toRepr  # PYCHOK for backward compatibility
    '''DEPRECATED, use method C{toRepr}.'''

    def toStr(self, *unused):  # PYCHOK _NamedDict
        '''Like C{str(dict)} but with C{floats} formatting by C{fstr}.
        '''
        return self._dot_(', .'.join(sorted(self.keys())))


class _NamedEnumItem(_NamedBase):
    '''(INTERNAL) Base class for items in a C{_NamedEnum} registery.
    '''
    _enum = None

    def __ne__(self, other):
        '''Compare this and an other item.

           @return: C{True} if different, C{False} otherwise.
        '''
        return not self.__eq__(other)

    def _instr(self, prec, *attrs, **kwds):
        '''(INTERNAL) Format, used by C{Conic}, C{Ellipsoid}, C{Transform}.
        '''
        t = ('%s=%r' % (_name_, self.name),)
        if attrs:
            t += pairs(((a, getattr(self, a)) for a in attrs),
                       prec=prec, ints=True)
        if kwds:
            t += pairs(kwds, prec=prec)
        return _COMMA_SPACE_.join(t)

    @property_doc_(''' the I{registered} name (C{str}).''')
    def name(self):
        '''Get the I{registered} name (C{str}).
        '''
        return self._name

    @name.setter  # PYCHOK setter!
    def name(self, name):
        '''Set the name, unless already registered.

           @arg name: New name (C{str}).
        '''
        if self._enum:
            raise _NameError(str(name), self, txt='registered')  # XXX _TypeError
        self._name = str(name)

    def _register(self, enum, name):
        '''(INTERNAL) Add this item as B{C{enum.name}}.

           @note: Don't register if name is empty or doesn't
                  start with a letter.
        '''
        if name:
            self.name = name
            if name[:1].isalpha():  # '_...' not registered
                enum.register(self)
                self._enum = enum

    def unregister(self):
        '''Remove this instance from its C{_NamedEnum} registry.

           @raise _AssertionError: Mismatch of this and registered item.

           @raise NameError: This item is unregistered.
        '''
        enum = self._enum
        if enum and self.name and self.name in enum:
            item = enum.unregister(self.name)
            if item is not self:
                raise _AssertionError('%r vs %r' % (item, self))


class _NamedTuple(tuple, _Named):
    '''(INTERNAL) Named C{tuple} with index I{and} attribute
       access to the items.

       @note: This class is similar to Python's C{namedtuple},
              but statically defined, lighter and limited.
    '''
    _iteration = None  #: (INTERNAL) Iteration number (C{int} or C{None}).
    _Names_    = ()    # must be tuple of 2 or more attr names

    def __new__(cls, *args):
        '''New L{_NamedTuple} initialized with B{C{positional}} arguments.
        '''
        self = tuple.__new__(cls, args)
        ns = self._Names_
        if not (isinstance(ns, tuple) and len(ns) > 1):  # XXX > 0
            raise _TypeError(_dot_(self.classname, _Names_), ns)
        if len(ns) != len(args) or not ns:
            raise LenError(cls, args=len(args), ns=len(ns))
        if _name_ in ns:
            t = unstr(_dot_(self.classname, _Names_), *ns)
            raise _NameError(_name_, _name_, txt=t)
        return self

    def __delattr__(self, name):
        if name in self._Names_:
            raise _TypeError('del', _dot_(self.classname, name), txt=_immutable_)
        elif name in (_name_, _name):
            _Named.__setattr__(self, name, '')  # XXX _Named.name.fset(self, '')
        else:
            tuple.__delattr__(self, name)

    def __getattr__(self, name):
        try:
            return tuple.__getitem__(self, self._Names_.index(name))
        except IndexError:
            raise _IndexError(_dot_(self.classname, '<name>'), name)
        except ValueError:
            return tuple.__getattribute__(self, name)

    def __getitem__(self, item):  # index, slice, etc.
        return tuple.__getitem__(self, item)

    def __repr__(self):
        return self.toRepr()

    def __setattr__(self, name, value):
        if name in self._Names_:
            raise _TypeError(_dot_(self.classname, name), value, txt=_immutable_)
        elif name in (_name_, _name):
            _Named.__setattr__(self, name, value)  # XXX _Named.name.fset(self, value)
        else:
            tuple.__setattr__(self, name, value)

    def __str__(self):
        return self.toStr()

    def items(self):
        '''Get the items as C{name, value} pairs (C{2-tuple}s).
        '''
        for i, n in enumerate(self._Names_):
            yield n, tuple.__getitem__(self, i)

    iteritems = items

    @property_RO
    def iteration(self):
        '''Get the iteration number (C{int} or C{None} if not available/applicable).
        '''
        return self._iteration

    def _xtend(self, NamedTuple, *items):
        '''(INTERNAL) Extend this C{_Tuple} with C{items} to an other C{NamedTuple}.
        '''
        if not (issubclassof(NamedTuple, _NamedTuple) and
               (len(self._Names_) + len(items)) == len(NamedTuple._Names_)
                and self._Names_ == NamedTuple._Names_[:len(self)]):
            raise TypeError('%s%r vs %s%r' % (self.classname, self._Names_,
                            NamedTuple.__name__, NamedTuple._Names_))
        return self._xnamed(NamedTuple(*(self + items)))

    def toRepr(self, prec=6, sep=_COMMA_SPACE_, **unused):  # PYCHOK signature
        '''Return the -Tuple items as C{name=value} string(s).

           @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                        Trailing zero decimals are stripped for B{C{prec}} values
                        of 1 and above, but kept for negative B{C{prec}} values.
           @kwarg sep: Optional separator to join (C{str}).

           @return: Tuple items (C{str}).
        '''
        return _item_ps(self.named, sep.join(pairs(self.items(), prec=prec)))

    toStr2 = toRepr  # PYCHOK for backward compatibility
    '''DEPRECATED, use method C{toRepr}.'''

    def toStr(self, prec=6, sep=_COMMA_SPACE_, **unused):  # PYCHOK signature
        '''Return the -Tuple items as string(s).

           @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                        Trailing zero decimals are stripped for B{C{prec}} values
                        of 1 and above, but kept for negative B{C{prec}} values.
           @kwarg sep: Optional separator to join (C{str}).

           @return: Tuple items (C{str}).
        '''
        return _PARENTH_ % (sep.join(reprs(self, prec=prec)),)


class Bearing2Tuple(_NamedTuple):
    '''2-Tuple C{(initial, final)} bearings, both in compass C{degrees360}.
    '''
    _Names_ = (_initial_, _final_)


class Bounds2Tuple(_NamedTuple):  # .geohash.py, .latlonBase.py, .points.py
    '''2-Tuple C{(latlonSW, latlonNE)} with the bounds' lower-left and
       upper-right corner as C{LatLon} instance.
    '''
    _Names_ = ('latlonSW', 'latlonNE')


class Bounds4Tuple(_NamedTuple):  # .geohash.py, .points.py
    '''4-Tuple C{(latS, lonW, latN, lonE)} with the bounds' lower-left
       C{(LatS, LowW)} and upper-right C{(latN, lonE)} corner lat- and
       longitudes.
    '''
    _Names_ = ('latS', 'lonW', 'latN', 'lonE')


class Destination2Tuple(_NamedTuple):  # .ellipsoidalKarney.py, -Vincenty.py
    '''2-Tuple C{(destination, final)}, C{destination} in C{LatLon}
       and C{final} bearing in compass C{degrees360}.
    '''
    _Names_ = ('destination', _final_)


class Destination3Tuple(_NamedTuple):  # .karney.py
    '''3-Tuple C{(lat, lon, final)}, destination C{lat}, C{lon} in
       and C{final} bearing in C{degrees}.
    '''
    _Names_ = (_lat_, _lon_, _final_)


class Distance2Tuple(_NamedTuple):  # .datum.py, .ellipsoidalBase.py
    '''2-Tuple C{(distance, initial)}, C{distance} in C{meter} and
       C{initial} bearing in compass C{degrees360}.
    '''
    _Names_ = (_distance_, _initial_)


class Distance3Tuple(_NamedTuple):  # .ellipsoidalKarney.py, -Vincenty.py
    '''3-Tuple C{(distance, initial, final)}, C{distance} in C{meter}
       and C{initial} and C{final} bearing, both in compass C{degrees360}.
    '''
    _Names_ = (_distance_, _initial_, _final_)


class Distance4Tuple(_NamedTuple):  # .formy.py, .points.py
    '''4-Tuple C{(distance2, delta_lat, delta_lon, unroll_lon2)} with
       the distance in C{degrees squared}, the latitudinal C{delta_lat}
       = B{C{lat2}}-B{C{lat1}}, the wrapped, unrolled and adjusted
       longitudinal C{delta_lon} = B{C{lon2}}-B{C{lon1}} and the
       C{unroll_lon2} unrollment for B{C{lon2}}.

       @note: Use Function L{degrees2m} to convert C{degrees squared}
              to C{meter} as M{degrees2m(sqrt(distance2), ...)} or
              M{degrees2m(hypot(delta_lat, delta_lon), ...)}.
    '''
    _Names_ = ('distance2', 'delta_lat', 'delta_lon', 'unroll_lon2')


class EasNor2Tuple(_NamedTuple):  # .css.py, .osgr.py, .ups.py, .utm.py, .utmupsBase.py
    '''2-Tuple C{(easting, northing)}, both in C{meter}.
    '''
    _Names_ = (_easting_, _northing_)


class EasNor3Tuple(_NamedTuple):  # .css.py, .lcc.py
    '''3-Tuple C{(easting, northing, height)}, all in C{meter}.
    '''
    _Names_ = (_easting_, _northing_, _height_)


class LatLon2Tuple(_NamedTuple):
    '''2-Tuple C{(lat, lon)} in C{degrees90} and C{degrees180}.
    '''
    _Names_ = (_lat_, _lon_)

    def to3Tuple(self, height):
        '''Extend this L{LatLon2Tuple} to a L{LatLon3Tuple}.

           @arg height: The height to add (C{scalar}).

           @return: A L{LatLon3Tuple}C{(lat, lon, height)}.

           @raise ValueError: Invalid B{C{height}}.
        '''
        from pygeodesy.units import Height
        return self._xtend(LatLon3Tuple, Height(height))

    def to4Tuple(self, height, datum):
        '''Extend this L{LatLon2Tuple} to a L{LatLon4Tuple}.

           @arg height: The height to add (C{scalar}).
           @arg datum: The datum to add (C{Datum}).

           @return: A L{LatLon4Tuple}C{(lat, lon, height, datum)}.

           @raise TypeError: If B{C{datum}} not a C{Datum}.

           @raise ValueError: Invalid B{C{height}}.
        '''
        return self.to3Tuple(height).to4Tuple(datum)


class LatLon3Tuple(_NamedTuple):
    '''3-Tuple C{(lat, lon, height)} in C{degrees90}, C{degrees180}
       and C{meter}.
    '''
    _Names_ = (_lat_, _lon_, _height_)

    def to4Tuple(self, datum):
        '''Extend this L{LatLon3Tuple} to a L{LatLon4Tuple}.

           @arg datum: The datum to add (C{Datum}).

           @return: A L{LatLon4Tuple}C{(lat, lon, height, datum)}.

           @raise TypeError: If B{C{datum}} not a C{Datum}.
        '''
        from pygeodesy.datum import Datum
        _xinstanceof(Datum, datum=datum)
        return self._xtend(LatLon4Tuple, datum)


class LatLon4Tuple(_NamedTuple):  # .cartesianBase.py, .css.py, .ecef.py, .lcc.py
    '''4-Tuple C{(lat, lon, height, datum)} in C{degrees90},
       C{degrees180}, C{meter} and L{Datum}.
    '''
    _Names_ = (_lat_, _lon_, _height_, _datum_)


class LatLonDatum3Tuple(_NamedTuple):  # .lcc.py, .osgr.py
    '''3-Tuple C{(lat, lon, datum)} in C{degrees90}, C{degrees180}
       and L{Datum}.
    '''
    _Names_ = (_lat_, _lon_, _datum_)


class LatLonPrec3Tuple(_NamedTuple):  # .gars.py, .wgrs.py
    '''3-Tuple C{(lat, lon, precision)} in C{degrees}, C{degrees}
       and C{int}.
    '''
    _Names_ = (_lat_, _lon_, _precision_)

    def to5Tuple(self, height, radius):
        '''Extend this L{LatLonPrec3Tuple} to a L{LatLonPrec5Tuple}.

           @arg height: The height to add (C{float} or C{None}).
           @arg radius: The radius to add (C{float} or C{None}).

           @return: A L{LatLonPrec5Tuple}C{(lat, lon, precision,
                    height, radius)}.
        '''
        return self._xtend(LatLonPrec5Tuple, height, radius)


class LatLonPrec5Tuple(_NamedTuple):  # .wgrs.py
    '''5-Tuple C{(lat, lon, precision, height, radius)} in C{degrees},
       C{degrees}, C{int} and C{height} or C{radius} in C{meter} (or
       C{None} if missing).
    '''
    _Names_ = (_lat_, _lon_, _precision_, _height_, _radius_)


class NearestOn3Tuple(_NamedTuple):  # .points.py, .sphericalTrigonometry.py
    '''3-Tuple C{(closest, distance, angle)} of the C{closest}
       point on the polygon, either a C{LatLon} instance or a
       L{LatLon3Tuple}C{(lat, lon, height)} and the C{distance}
       and C{angle} to the C{closest} point are in C{meter}
       respectively compass C{degrees360}.
    '''
    _Names_ = ('closest', _distance_, _angle_)


class PhiLam2Tuple(_NamedTuple):  # .frechet.py, .hausdorff.py, .latlonBase.py, .points.py, .vector3d.py
    '''2-Tuple C{(phi, lam)} with latitude C{phi} in C{radians[PI_2]}
       and longitude C{lam} in C{radians[PI]}.

       @note: Using C{phi/lambda} for lat-/longitude in C{radians}
              follows Chris Veness' U{convention
              <https://www.Movable-Type.co.UK/scripts/latlong.html>}.
    '''
    _Names_ = (_phi_, _lam_)

    def to3Tuple(self, height):
        '''Extend this L{PhiLam2Tuple} to a L{PhiLam3Tuple}.

           @arg height: The height to add (C{scalar}).

           @return: A L{PhiLam3Tuple}C{(phi, lam, height)}.

           @raise ValueError: Invalid B{C{height}}.
        '''
        from pygeodesy.units import Height
        return self._xtend(PhiLam3Tuple, Height(height))

    def to4Tuple(self, height, datum):
        '''Extend this L{PhiLam2Tuple} to a L{PhiLam4Tuple}.

           @arg height: The height to add (C{scalar}).
           @arg datum: The datum to add (C{Datum}).

           @return: A L{PhiLam4Tuple}C{(phi, lam, height, datum)}.

           @raise TypeError: If B{C{datum}} not a C{Datum}.

           @raise ValueError: Invalid B{C{height}}.
        '''
        return self.to3Tuple(height).to4Tuple(datum)


class PhiLam3Tuple(_NamedTuple):  # .nvector.py, extends -2Tuple
    '''3-Tuple C{(phi, lam, height)} with latitude C{phi} in
       C{radians[PI_2]}, longitude C{lam} in C{radians[PI]} and
       C{height} in C{meter}.

       @note: Using C{phi/lambda} for lat-/longitude in C{radians}
              follows Chris Veness' U{convention
              <https://www.Movable-Type.co.UK/scripts/latlong.html>}.
    '''
    _Names_ = (_phi_, _lam_, _height_)

    def to4Tuple(self, datum):
        '''Extend this L{PhiLam3Tuple} to a L{PhiLam4Tuple}.

           @arg datum: The datum to add (C{Datum}).

           @return: A L{PhiLam4Tuple}C{(phi, lam, height, datum)}.

           @raise TypeError: If B{C{datum}} not a C{Datum}.
        '''
        from pygeodesy.datum import Datum
        _xinstanceof(Datum, datum=datum)
        return self._xtend(PhiLam4Tuple, datum)


class PhiLam4Tuple(_NamedTuple):  # extends -3Tuple
    '''4-Tuple C{(phi, lam, height, datum)} with latitude C{phi} in
       C{radians[PI_2]}, longitude C{lam} in C{radians[PI]}, C{height}
       in C{meter} and L{Datum}.

       @note: Using C{phi/lambda} for lat-/longitude in C{radians}
              follows Chris Veness' U{convention
              <https://www.Movable-Type.co.UK/scripts/latlong.html>}.
    '''
    _Names_ = (_phi_, _lam_, _height_, _datum_)


class Points2Tuple(_NamedTuple):  # .formy.py, .latlonBase.py
    '''2-Tuple C{(number, points)} with the C{number} of points
       and -possible reduced- C{list} or C{tuple} of C{points}.
    '''
    _Names_ = (_number_, _points_)


class Vector3Tuple(_NamedTuple):
    '''3-Tuple C{(x, y, z)} of (geocentric) components, all in
       C{meter} or C{units}.
    '''
    _Names_ = (_x_, _y_, _z_)

    def to4Tuple(self, h):
        '''Extend this L{Vector3Tuple} to a L{Vector4Tuple}.

           @arg h: The height to add (C{scalar}).

           @return: A L{Vector4Tuple}C{(x, y, z, h)}.

           @raise ValueError: Invalid B{C{h}}.
        '''
        from pygeodesy.units import Height
        return self._xtend(Vector4Tuple, Height(h, name=_h_))


class Vector4Tuple(_NamedTuple):  # .nvector.py
    '''4-Tuple C{(x, y, z, h)} of (geocentric) components, all
       in C{meter} or C{units}.
    '''
    _Names_ = (_x_, _y_, _z_, _h_)


def callername(up=1, dflt=NN, source=False):
    '''Get the name of the calling callable.

       @kwarg up: Number of call stack frames up (C{int}).
       @kwarg dflt: Default return value (C{any}).
       @kwarg source: Include source file name and line
                      number (C{bool}).

       @return: Name of the non-internal callable (C{str})
                or B{C{dflt}} if none found.
    '''
    try:
        from sys import _getframe
        for u in range(up, up + 32):
            f = _getframe(u)
            n = f.f_code.co_name
            if n and (n.startswith(_DUNDER_) or
                  not n.startswith(_UNDERSCORE_)):
                if source:
                    from os.path import basename
                    n = NN.join((n, _AT_, basename(f.f_code.co_filename),
                                    _COLON_, str(f.f_lineno)))
                return n
    except (AttributeError, ImportError):
        pass
    return dflt


def _callname(name, class_name, self_name, up=1):  # imported by .points
    '''(INTERNAL) Assemble the name for an invokation.
    '''
    n, c = class_name, callername(up=up + 1)
    if c:
        n = _dot_(n, _item_ps(c, name))
    if self_name:
        n = '%s %r' % (n, self_name)
    return n


def classname(inst, prefixed=None):
    '''Return the instance' class name optionally prefixed with the
       module name.

       @arg inst: The object (any C{type}).
       @kwarg prefixed: Include the module name (C{bool}), see
                        function C{classnaming}.

       @return: The B{C{inst}}'s C{[module.]class} name (C{str}).
    '''
    return modulename(inst.__class__, prefixed or
             (getattr(inst, 'classnaming', classnaming()) if
                                      prefixed is None else False))


def classnaming(prefixed=None):
    '''Get/set the default class naming for C{[module.]class} names.

       @kwarg prefixed: Include the module name (C{bool}).

       @return: Previous class naming setting (C{bool}).
    '''
    t = _Named._classnaming
    if prefixed in (True, False):
        _Named._classnaming = prefixed
    return t


def modulename(clas, prefixed=None):
    '''Return the class name optionally prefixed with the
       module name.

       @arg clas: The class (any C{class}).
       @kwarg prefixed: Include the module name (C{bool}), see
                        function C{classnaming}.

       @return: The B{C{class}}'s C{[module.]class} name (C{str}).
    '''
    try:
        n = clas.__name__
    except AttributeError:
        n = '--'
    if prefixed or (classnaming() if prefixed is None else False):
        try:
            m = clas.__module__.rsplit(_DOT_, 1)
            n = _dot_(*(m[1:] + [n]))
        except AttributeError:
            pass
    return n


def nameof(inst):
    '''Get the name of an instance.

       @arg inst: The object (any C{type}).

       @return: The instance' name (C{str}) or C{""}.
    '''
    return getattr(inst, _name_, NN)


def _notError(inst, name, args, kwds):  # PYCHOK no cover
    '''(INTERNAL) Format an error message.
    '''
    n = _dot_(classname(inst, prefixed=True), _dunder_name(name, name))
    m = _COMMA_SPACE_.join(modulename(c, prefixed=True) for c in inst.__class__.__mro__[1:-1])
    t = '%s, MRO(%s)' % (unstr(n, *args, **kwds), m)
    return t


def notImplemented(inst, name, *args, **kwds):  # PYCHOK no cover
    '''Raise a C{NotImplementedError} for a missing method or property.

       @arg inst: Instance (C{any}).
       @arg name: Method, property or name (C{str} or C{callable}).
       @arg args: Method or property positional arguments (any C{type}s).
       @arg kwds: Method or property keyword arguments (any C{type}s).
    '''
    t = _notError(inst, name, args, kwds)
    raise _NotImplementedError(t, txt=notImplemented.__name__)


def notOverloaded(inst, name, *args, **kwds):  # PYCHOK no cover
    '''Raise an C{AssertionError} for a method or property not overloaded.

       @arg inst: Instance (C{any}).
       @arg name: Method, property or name (C{str} or C{callable}).
       @arg args: Method or property positional arguments (any C{type}s).
       @arg kwds: Method or property keyword arguments (any C{type}s).
    '''
    t = _notError(inst, name, args, kwds)
    raise _AssertionError(t, txt=notOverloaded.__name__)


__all__ += _ALL_DOCS(_Named,
                     _NamedBase,  # _NamedDict,
                     _NamedEnum, _NamedEnumItem,  # _NamedTuple,
                      Bearing2Tuple, Bounds2Tuple, Bounds4Tuple,
                      Destination2Tuple, Destination3Tuple,
                      Distance2Tuple, Distance3Tuple, Distance4Tuple,
                      EasNor2Tuple, EasNor3Tuple,
                      LatLon2Tuple, LatLon3Tuple, LatLon4Tuple,
                      LatLonDatum3Tuple, LatLonPrec3Tuple, LatLonPrec5Tuple,
                      NearestOn3Tuple,
                      PhiLam2Tuple, PhiLam3Tuple, PhiLam4Tuple, Points2Tuple,
                      Vector3Tuple, Vector4Tuple)

if __name__ == '__main__':

    from sys import argv, exit  # PYCHOK shadows exit

    from pygeodesy.lazily import _FOR_DOCS

    if not _FOR_DOCS:
        exit('%s\n' % (' '.join('usage: env PYGEODESY_FOR_DOCS=1 python -m'.split() + argv),))

    ls = locals()
    for n in __all__:
        if n not in ls:
            raise NameError('%s %r not in %s' % ('__all__', n, _dot_('named', 'locals')))
    for n, o in ls.items():
        if n not in __all__ and not n.startswith(_UNDERSCORE_) \
                            and getattr(o, '__module__', '') == __name__:
            raise NameError('%s %r not in %s' % ('locals', n, _dot_('named', '__all__')))

    print('%s: %s vs %s OK' % (argv[0], '__all__', 'locals'))

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

# % env PYGEODESY_FOR_DOCS=1 python -m pygeodesy.named
# all 71 locals OK
