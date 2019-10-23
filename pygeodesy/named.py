
# -*- coding: utf-8 -*-

u'''(INTERNAL) Classes C{_Named}, C{_NamedDict}, C{_NamedInt},
C{_NamedStr} and C{_NamedTuple} with nameable instances and several
subclasses thereof.

In addition, the items in a C{_NamedDict} are accessable as attributes
and the items in a C{_NamedTuple} can be named to be accessable as
attributes, similar to standard Python C{namedtuple}s.

Results previously returned as tuples by C{pygeodesy} functions and
class methods are now instances of some C{...Tuple} class, all
sub-classes of C{_NamedTuple} defined here.

@newfield example: Example, Examples
'''

from pygeodesy.fmath import fStr, _IsNotError
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.utily import issubclassof, property_RO, _Strs, unStr

# XXX 'FsumDelta2Tuple' is in _ALL_LAZY.named
__all__ = _ALL_LAZY.named + _ALL_DOCS(  # '_Named', '_NamedBase',
         # '_NamedDict', '_NamedEnum', '_NamedEnumItem',
         # '_NamedInt', '_NamedStr', '_NamedTuple',
         'Bearing2Tuple', 'Bounds2Tuple', 'Bounds4Tuple',
         'ClipCS3Tuple', 'ClipSH3Tuple', 'Curvature2Tuple',
         'Destination2Tuple',
         'Distance2Tuple', 'Distance3Tuple', 'Distance4Tuple',
         'EasNor2Tuple', 'EasNor3Tuple',
         'EasNorAziRk4Tuple', 'EasNorExact4Tuple', 'EasNorRadius3Tuple',
         'Elevation2Tuple', 'GeoidHeight2Tuple', 'GeoidHeight5Tuple',
         'LatLon2Tuple', 'LatLon3Tuple', 'LatLon4Tuple',
         'LatLonAziRk4Tuple', 'LatLonDatum3Tuple', 'LatLonDatum5Tuple',
         'LatLonExact4Tuple', 'LatLonPrec3Tuple', 'LatLonPrec5Tuple',
         'Mgrs4Tuple', 'Mgrs6Tuple',
         'NearestOn3Tuple', 'NearestOn5Tuple',
         'Ned3Tuple', 'Neighbors8Dict',
         'PhiLam2Tuple', 'PhiLam3Tuple',
         'Point3Tuple', 'Points2Tuple', 'Shape2Tuple',
         'UtmUps2Tuple', 'UtmUps4Tuple', 'UtmUps5Tuple', 'UtmUps8Tuple',
         'UtmUpsLatLon5Tuple',
         'Vector3Tuple', 'Vector4Tuple')
__version__ = '19.10.19'

_NAME_ = 'name'  # __NAME gets mangled in class


def _xattrs(inst, other, *attrs):
    '''(INTERNAL) Copy attribute values from B{C{other}} to B{C{inst}}.

       @param inst: Instance to copy attribute values to.
       @param other: Instance to copy attribute values from.
       @param attrs: Attribute names (C{str})s.

       @return: The B{C{inst}}, updated.

       @raise AttributeError: An B{C{attr}} doesn't exist or is not settable.
    '''
    def _getattr(o, a):
        if hasattr(o, a):
            return getattr(o, a)
        raise AttributeError('.%s invalid: %r' % (a, o))

    for a in attrs:
        s = _getattr(other, a)
        if _getattr(inst, a) != s:
            setattr(inst, a, s)  # not settable?
    return inst


def _xnamed(inst, name):
    '''(INTERNAL) Set the instance' C{.name = }B{C{name}}.

       @param inst: The instance (C{_Named}).
       @param name: The name (C{str}).

       @return: The B{C{inst}}, named if not named before.
    '''
    if name and isinstance(inst, _Named) and not inst.name:
        inst.name = name
    return inst


class _Named(object):
    '''(INTERNAL) Root class for named objects.
    '''
    _name        = ''     #: (INTERNAL) name (C{str})
    _classnaming = False  #: (INTERNAL) prefixed (C{bool})

    def __repr__(self):
        return '<%s at %#x>' % (self.named2, id(self))

    __str__ = __repr__

    def _xcopy(self, *attrs):
        '''(INTERNAL) I{Must be overloaded}.
        '''
        self._notOverloaded(self._xcopy.__name__, *attrs)

    @property_RO
    def classname(self):
        '''Get this object's C{[module.]class} name (C{str}), see
           property C{classnaming} and function C{classnaming}.
        '''
        return classname(self, prefixed=self._classnaming)

    @property
    def classnaming(self):
        '''Get the class naming (C{bool}).
        '''
        return self._classnaming

    @classnaming.setter  # PYCHOK setter!
    def classnaming(self, prefixed):
        '''Set the class naming for C{[module.].class} names.

           @param prefixed: Include the module name (C{bool}).
        '''
        self._classnaming = bool(prefixed)

    def classof(self, *args, **kwds):
        '''Create another instance of this very class.

           @param args: Optional, positional arguments.
           @keyword kwds: Optional, keyword arguments.

           @return: New instance (B{self.__class__}).
        '''
        return _xnamed(self.__class__(*args, **kwds), self.name)

    def copy(self):
        '''Make a copy of this instance.

           @return: The copy (C{This class} or subclass thereof).
        '''
        return self._xcopy()

    __copy__     = copy
#   __deepcopy__ = copy

    @property
    def name(self):
        '''Get the name (C{str}).
        '''
        return self._name

    @name.setter  # PYCHOK setter!
    def name(self, name):
        '''Set the name.

           @param name: New name (C{str}).
        '''
        self._name = str(name)
        # to set the name ffrom a sub-class, use
        # self.name = name or
        # _Named.name.fset(self, name), but
        # not _Named(self).name = name

    @property_RO
    def named(self):
        '''Get the name I{or} class name or C{""} (C{str}).
        '''
        return self.name or self.classname

    @property_RO
    def named2(self):
        '''Get the class name I{and/or} thename or C{""} (C{str}).
        '''
        n, c = self.name, self.classname
        return ('%s %r' % (c, n)) if c and n else (c or n)

    def _notOverloaded(self, name, *args, **kwds):
        '''Raise an error for a method or property not overloaded.
        '''
        n = '%s %s.%s' % (self._notOverloaded.__name__, self.classname, name)
        raise AssertionError(unStr(n, *args, **kwds))

    def toStr(self, **unused):
        '''Default C{str(self)}.
        '''
        return str(self)

    def toStr2(self, **unused):
        '''Default C{repr(self)}.
        '''
        return repr(self)

    def _xrenamed(self, inst):
        '''(INTERNAL) Rename the instance' C{.name = self.name}.

           @param inst: The instance (C{_Named}).

           @return: The B{C{inst}}, named if not named before.
        '''
        if not isinstance(inst, _Named):
            raise _IsNotError('valid', inst=inst)

        if inst.name != self.name:
            inst.name = self.name
        return inst

    def _xnamed(self, inst, name=''):
        '''(INTERNAL) Set the instance' C{.name = self.name}.

           @param inst: The instance (C{_Named}).
           @keyword name: Optional name, overriding C{self.name} (C{str}).

           @return: The B{C{inst}}, named if not named before.
        '''
        n = name or self.name
        return _xnamed(inst, n) if n else inst


class _NamedBase(_Named):
    '''(INTERNAL) Base class with name.
    '''

    def __repr__(self):
        return self.toStr2()

    def __str__(self):
        return self.toStr()

    def _update(self, unused):
        '''(INTERNAL) To be overloaded.
        '''
        pass

#   def notImplemented(self, attr):
#       '''Raise error for a missing method, function or attribute.
#
#          @param attr: Attribute name (C{str}).
#
#          @raise NotImplementedError: No such attribute.
#       '''
#       c = self.__class__.__name__
#       return NotImplementedError('%s.%s' % (c, attr))

    def others(self, other, name='other'):
        '''Check this and an other instance for type compatiblility.

           @param other: The other instance (any C{type}).
           @keyword name: Optional, name for other (C{str}).

           @return: C{None}.

           @raise TypeError: Mismatch of this and B{C{other}} C{type}.
        '''
        if not (isinstance(self, other.__class__) or
                isinstance(other, self.__class__)):
            raise TypeError('type(%s) mismatch: %s vs %s' % (name,
                             classname(other), self.classname))

    def toStr(self, **kwds):
        '''(INTERNAL) Must be overloaded.

           @param kwds: Optional, keyword arguments.
        '''
        self._notOverloaded(self.toStr.__name__, **kwds)

#   def toStr(self, *args, **kwds):
#       if kwds or args:
#           t = list('%s=%r' % t for t in sorted(kwds.items()))
#           if args:
#               t = map(repr, args) + t
#           s = ', '.join(t)
#       else:
#           s = super(self.__class__, self).__str__()
#       return s

    def toStr2(self, **kwds):
        '''(INTERNAL) To be overloaded.

           @keyword kwds: Optional, keyword arguments.

           @return: L{toStr}() plus keyword arguments (as C{str}).
        '''
        t = self.toStr(**kwds).lstrip('([{').rstrip('}])')
        return '%s(%s)' % (self.classname, t)

#   def toStr2(self, *args, **kwds)
#       if args or kwds:
#           s = self.toStr(*args, **kwds)
#       else:
#           s = super(self.__class__, self).__repr__()
#       return '%s(%s)' % (self.named, s)  # clip(s)


class _NamedDict(dict, _Named):
    '''(INTERNAL) Named C{dict} with key I{and} attribute
       access to the items.
    '''

    def __init__(self, *args, **kwds):
        if kwds:
            if _NAME_ in kwds:
                _Named.name.fset(self, kwds.pop(_NAME_))  # see _Named.name
            dict.__init__(self, kwds)
        if args:
            raise ValueError('%s(%s) invalid: %r' %
                             (self.classname, 'args', args,))

    def __delattr__(self, name):
        if name in dict.keys(self):
            dict.pop(name)
        elif name == _NAME_:
            dict.__setattr__(self, name, '')
        else:
            dict.__delattr__(self, name)

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            if name == _NAME_:
                return _Named.name.fget(self)
            return dict.__getattr__(self, name)

    def __getitem__(self, key):
        if key == _NAME_:
            raise KeyError('%s[%r]' % (self.classname, key))
        return dict.__getitem__(self, key)

    def __repr__(self):
        t = ', '.join('%s=%r' % t for t in self.items())
        return '%s(%s)' % (self.name, t)

    def __setattr__(self, name, value):
        if name in dict.keys(self):
            dict.__setitem__(self, name, value)  # self[name] = value
        else:
            dict.__setattr__(self, name, value)

    def __setitem__(self, key, value):
        if key == _NAME_:
            raise KeyError('%s[%r] = %r' % (self.classname, key, value))
        dict.__setitem__(self, key, value)

    def __str__(self):
        return dict.__repr__(self)  # dict.__str__(self) fails


class Neighbors8Dict(_NamedDict):
    '''8-Dict C{(N, NE, E, SE, S, SW, W, NW)} of L{Geohash}es,
       providing key I{and} attribute access to the items.
    '''
    _Keys_ = ('N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW')


class _NamedEnum(_NamedDict):
    '''(INTERNAL) Enum-like C{_NamedDict} with attribute access
       restricted to valid keys.
    '''
    _item_classes = ()

    def __init__(self, name, *classes):
        '''New C{_NamedEnum}.

           @param name: Name (C{str}).
        '''
        if classes:
            self._item_classes = classes
        if name:
            _Named.name.fset(self, name)  # see _Named.name

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            if name == _NAME_:
                return _NamedDict.name.fget(self)
        raise AttributeError("%s.%s doesn't exist" % (self.named, name))

    def __repr__(self):
        return '\n'.join('%s.%s: %r,' % (self.name, n, v) for n, v in sorted(self.items()))

    def __str__(self):
        return (self.name + '.') + ', .'.join(sorted(self.keys()))

    def _assert(self, **kwds):
        '''(INTERNAL) Check names against given, registered names.
        '''
        for a, v in kwds.items():
            assert self[a] is v and getattr(self, a) \
                                and self.find(v) == a

    def find(self, item):
        '''Find a registered item.

           @param item: The item to look for (any C{type}).

           @return: If found the B{C{item}}'s name (C{str}), C{None} otherwise.
        '''
        for k, v in self.items():
            if v is item:
                return k
        return None

    def register(self, item):
        '''Registed a new item.

           @param item: The item (any C{type}).

           @return: The item name (C{str}).

           @raise NameError: An B{C{item}} already registered with
                             that name or the B{C{item}} has no, an
                             empty or an invalid name.

           @raise TypeError: The B{C{item}} type invalid.
        '''
        if not (self._item_classes and isinstance(item, self._item_classes)):
            raise TypeError('%s.%s: %r' % ('item', 'type', item))

        try:
            n = item.name
            if not (n and n.replace('_', '').isalnum() and isinstance(n, _Strs)):
                raise ValueError
        except (AttributeError, ValueError):
            raise NameError('%s.%s %s: %r' % ('item', 'name', 'invalid', item))
        if n in self:
            raise NameError('%s.%s %s: %r' % (self.named, n, 'exists', item))
        self[n] = item

    def unregister(self, name_or_item):
        '''Remove a registered item.

           @param name_or_item: Name (C{str}) of or the item (any C{type}).

           @return: The unregistered item.

           @raise NameError: No item with that B{C{name}}.

           @raise ValueError: No such item.
        '''
        name = self.find(name_or_item)
        if name is None:
            if not isinstance(name_or_item, _Strs):
                raise ValueError('no such %r' % (name_or_item,))
            name = name_or_item
        try:
            item = dict.pop(self, name)
        except KeyError:
            raise NameError('no %s.%r' % (self.named, name))
        item._enum = None
        return item


class _NamedEnumItem(_NamedBase):
    '''(INTERNAL) Base class for items in a C{_NamedEnum} registery.
    '''
    _enum = None

    def __ne__(self, other):
        '''Compare this and an other item.

           @return: C{True} if different, C{False} otherwise.
        '''
        return not self.__eq__(other)

    def _fStr(self, prec, *attrs, **others):
        '''(INTERNAL) Format.
        '''
        t = fStr([getattr(self, a) for a in attrs], prec=prec, sep=' ', ints=True)
        t = ['%s=%s' % (a, v) for a, v in zip(attrs, t.split())]
        if others:
            t += ['%s=%s' % (a, v) for a, v in sorted(others.items())]
        return ', '.join(['name=%r' % (self.name,)] + t)

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

    @property
    def name(self):
        '''Get the I{registered, immutable} name (C{str}).
        '''
        return self._name

    @name.setter  # PYCHOK setter!
    def name(self, name):
        '''Set the name, unless already registered.

           @param name: New name (C{str}).
        '''
        if self._enum:
            raise ValueError('%s, %s: %r' % ('registered', 'immutable', self))
        self._name = str(name)

    def unregister(self):
        '''Remove this instance from its C{_NamedEnum} registry.

           @raise AssertionError: Mismatch of this and registered item.

           @raise NameError: This item is unregistered.
        '''
        enum = self._enum
        if enum and self.name and self.name in enum:
            item = enum.unregister(self.name)
            if item is not self:
                raise AssertionError('%r vs %r' % (item, self))


class _NamedInt(int, _Named):
    '''(INTERNAL) Named C{int}.
    '''
    def __repr__(self):
        return "%s(%s)" % (self.named, int.__repr__(self))

    def __str__(self):
        return int.__str__(self)


class _NamedStr(str, _Named):
    '''(INTERNAL) Named C{str}.
    '''
    def __repr__(self):
        return "%s(%s)" % (self.named, str.__repr__(self))

    def __str__(self):
        return str.__str__(self)


class _NamedTuple(tuple, _Named):
    '''(INTERNAL) Named C{tuple} with index I{and} attribute
       access to the items.

       @note: This class is similar to Python's C{namedtuple},
              but limited and lighter in setup and weight.
    '''
    _Names_ = ()  # must be tuple of 2 or more attr names

    def __new__(cls, *args):
        '''New L{_NamedTuple} initialized with B{C{positional}} arguments.
        '''
        self = tuple.__new__(cls, args)
        ns = self._Names_
        if not (isinstance(ns, tuple) and len(ns) > 1):
            raise TypeError('%s.%s inot valid: %r' %
                            (self.classname, '_Names_', ns))
        if len(ns) != len(args) or not ns:
            raise ValueError('%s(%s) not valid: %r[%s] vs %s' %
                             (self.classname, 'args',
                              args, len(args), len(ns)))
        if _NAME_ in ns:
            raise NameError('%s.%s not valid: %r' %
                             (self.classname, '_Names_', _NAME_))
        return self

    def __delattr__(self, name):
        if name in self._Names_:
            raise TypeError('%s not mutable: %s .%s' %
                            (self.classname, 'del', name))
        elif name == _NAME_:
            tuple.__setattr__(self, name, '')
        else:
            tuple.__delattr__(self, name)

    def __getattr__(self, name):
        try:
            return tuple.__getitem__(self, self._Names_.index(name))
        except IndexError:
            raise IndexError('%s.%s not valid: %r' %
                             (self.classname, '<name>', name))
        except ValueError:
            return tuple.__getattribute__(self, name)

    def __getitem__(self, item):  # index, slice, etc.
        return tuple.__getitem__(self, item)

    def __repr__(self):
        t = ', '.join('%s=%r' % t for t in self.items())
        return '%s(%s)' % (self.name, t)

    def __setattr__(self, name, value):
        if name in self._Names_:
            raise TypeError('%s not mutable: %s=%r' %
                            (self.classname, name, value))
        else:
            tuple.__setattr__(self, name, value)

    def __str__(self):
        return tuple.__repr__(self)  # tuple.__str__(self) fails

    def items(self):
        '''Get the items as C{name, value} pairs (C{2-tuple}s).
        '''
        for i, n in enumerate(self._Names_):
            yield n, tuple.__getitem__(self, i)

    iteritems = items

    def _xtend(self, namedTuple, *items):
        '''(INTERNAL) Extend this C{_Tuple} with C{items} to an other C{namedTuple}.
        '''
        if not (issubclassof(namedTuple, _NamedTuple) and
               (len(self._Names_) + len(items)) == len(namedTuple._Names_)
                and self._Names_ == namedTuple._Names_[:len(self)]):
            raise TypeError('%s%r vs %s%r' % (self.classname, self._Names_,
                            namedTuple.__name__, namedTuple._Names_))
        return self._xnamed(namedTuple(*(self + items)))


class Bearing2Tuple(_NamedTuple):
    '''2-Tuple C{(initial, final)} bearings, both in compass C{degrees360}.
    '''
    _Names_ = ('initial', 'final')


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


class ClipCS3Tuple(_NamedTuple):  # .clipy.py
    '''3-Tuple C{(start, end, index)} for each edge of a I{clipped}
       path with the C{start} and C{end} points (C{LatLon}) of the
       portion of the edge inside or on the clip box and the C{index}
       (C{int}) of the edge in the original path.
    '''
    _Names_ = ('start', 'end', 'index')


class ClipSH3Tuple(_NamedTuple):  # .clipy.py
    '''3-Tuple C{(start, end, original)} for each edge of a I{clipped}
       polygon, the C{start} and C{end} points (C{LatLon}) of the
       portion of the edge inside or on the clip region and the
       C{original} indicates whether the edge is part of the original
       polygon or part of the clip region (C{bool}).
    '''
    _Names_ = ('start', 'end', 'original')


class Curvature2Tuple(_NamedTuple):  # .datum.py
    '''2-Tuple C{(meridional, prime_vertical)} of radii of curvature,
       both in C{meter}.
    '''
    _Names_ = ('meridional', 'prime_vertical')


class Destination2Tuple(_NamedTuple):  # .ellipsoidalKarney.py, -Vincenty.py
    '''2-Tuple C{(destination, final)}, C{destination} in C{LatLon}
       and C{final} bearing in compass C{degrees360}.
    '''
    _Names_ = ('destination', 'final')


class Distance2Tuple(_NamedTuple):  # .datum.py, .ellipsoidalBase.py
    '''2-Tuple C{(distance, initial)}, C{distance} in C{meter} and
       C{initial} bearing in compass C{degrees360}.
    '''
    _Names_ = ('distance', 'initial')


class Distance3Tuple(_NamedTuple):  # .ellipsoidalKarney.py, -Vincenty.py
    '''3-Tuple C{(distance, initial, final)}, C{distance} in C{meter}
       and C{initial} and C{final} bearing, both in compass C{degrees360}.
    '''
    _Names_ = ('distance', 'initial', 'final')


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


class EasNor2Tuple(_NamedTuple):  # .css.py, .osgr.py
    '''2-Tuple C{(easting, northing)}, both in C{meter}.
    '''
    _Names_ = ('easting', 'northing')


class EasNor3Tuple(_NamedTuple):  # .css.py, .lcc.py
    '''3-Tuple C{(easting, northing, height)}, all in C{meter}.
    '''
    _Names_ = ('easting', 'northing', 'height')


class EasNorAziRk4Tuple(_NamedTuple):
    '''4-Tuple C{(easting, northing, azimuth, reciprocal)} for the
       Cassini-Soldner location with C{easting} and C{northing} in
       C{meters}, the C{azimuth} of easting direction azimuth and
       C{reciprocal} the reciprocal of azimuthal northing scale
       reciprocal, both in C{degrees}.
    '''
    _Names_ = ('easting', 'northing', 'azimuth', 'reciprocal')


class EasNorExact4Tuple(_NamedTuple):  # .etm.py
    '''4-Tuple C{(easting, northing, convergence, scale)} in
       C{meter}, C{meter}, {degrees} and C{scalar}.
    '''
    _Names_ = ('easting', 'northing', 'convergence', 'scale')


class EasNorRadius3Tuple(_NamedTuple):  # .webmercator.py
    '''3-Tuple C{(easting, northing, radius)}, all in C{meter}.
    '''
    _Names_ = ('easting', 'northing', 'radius')


class Elevation2Tuple(_NamedTuple):  # .elevations.py
    '''2-Tuple C{(elevation, data_source)} in C{meter} and C{str}.
    '''
    _Names_ = ('elevation', 'data_source')


class GeoidHeight2Tuple(_NamedTuple):  # .elevations.py
    '''2-Tuple C{(height, model_name)}, geoid C{height} in C{meter}
       and C{model_name} as C{str}.
    '''
    _Names_ = ('height', 'model_name')


class GeoidHeight5Tuple(_NamedTuple):  # .geoids.py
    '''5-Tuple C{(lat, lon, egm84, egm96, egm2008)} for U{GeoidHeights.dat
       <https://SourceForge.net/projects/geographiclib/files/testdata/>}
       tests with the heights for 3 different EGM grids with C{-90.0 <=
       lat <= 90.0} and C{-180.0 <= lon <= 180.0} degrees (and C{lon}
       converted from the original C{0.0 <= EasterLon <= 360.0}).
    '''
    _Names_ = ('lat', 'lon', 'egm84', 'egm96', 'egm2008')


class LatLon2Tuple(_NamedTuple):
    '''2-Tuple C{(lat, lon)} in C{degrees[90]} and C{degrees[180]}.
    '''
    _Names_ = ('lat', 'lon')

    def _3Tuple(self, height):
        '''(INTERNAL) Extend to a L{LatLon3Tuple}.
        '''
        return self._xtend(LatLon3Tuple, height)

    def _4Tuple(self, height, datum):
        '''(INTERNAL) Extend to a L{LatLon4Tuple}.
        '''
        return self._xtend(LatLon4Tuple, height, datum)


class LatLon3Tuple(_NamedTuple):
    '''3-Tuple C{(lat, lon, height)} in C{degrees[90]}, C{degrees[180]}
       and C{meter}.
    '''
    _Names_ = ('lat', 'lon', 'height')

    def _4Tuple(self, datum):
        '''(INTERNAL) Extend to a L{LatLon4Tuple}.
        '''
        return self._xtend(LatLon4Tuple, datum)


class LatLon4Tuple(_NamedTuple):  # .css.py
    '''4-Tuple C{(lat, lon, height, datum)} in C{degrees90},
       C{degrees180}, C{meter} and L{Datum}.
    '''
    _Names_ = ('lat', 'lon', 'height', 'datum')


class LatLonAziRk4Tuple(_NamedTuple):  # .css.py
    '''4-Tuple C{(lat, lon, azimuth, reciprocal)}, all in C{degrees}
       where C{azimuth} is the azimuth of easting direction and
       C{reciprocal} the reciprocal of azimuthal northing scale.
    '''
    _Names_ = ('lat', 'lon', 'azimuth', 'reciprocal')


class LatLonDatum3Tuple(_NamedTuple):  # .lcc.py, .osgr.py
    '''3-Tuple C{(lat, lon, datum)} in C{degrees[90]}, C{degrees[180]}
       and L{Datum}.
    '''
    _Names_ = ('lat', 'lon', 'datum')


class LatLonDatum5Tuple(_NamedTuple):  # .ups.py, .utm.py
    '''5-Tuple C{(lat, lon, datum, convergence, scale)} in
       C{degrees[90]}, C{degrees[180]}, L{Datum}, {degrees}
       and C{float}.
    '''
    _Names_ = ('lat', 'lon', 'datum', 'convergence', 'scale')


class LatLonExact4Tuple(_NamedTuple):  # .etm.py
    '''4-Tuple C{(lat, lon, convergence, scale)} in C{degrees},
       C{degrees180}, C{degrees} and C{sclar}.
    '''
    _Names_ = ('lat', 'lon', 'convergence', 'scale')


class LatLonPrec3Tuple(_NamedTuple):  # .gars.py, .wgrs.py
    '''3-Tuple C{(lat, lon, precision)} in C{degrees}, C{degrees}
       and C{int}.
    '''
    _Names_ = ('lat', 'lon', 'precision')


class LatLonPrec5Tuple(_NamedTuple):  # .wgrs.py
    '''5-Tuple C{(lat, lon, precision, height, radius)} in C{degrees},
       C{degrees}, C{int} and C{height} or C{radius} in C{meter} (or
       C{None} if missing).
    '''
    _Names_ = ('lat', 'lon', 'precision', 'height', 'radius')


class Mgrs4Tuple(_NamedTuple):  # see .mgrs.py
    '''4-Tuple C{(zone, digraph, easting, northing)}, C{zone} and
       C{digraph} as C{str}, C{easting} and C{northing} in C{meter}.
    '''
    _Names_ = ('zone', 'digraph', 'easting', 'northing')


class Mgrs6Tuple(_NamedTuple):
    '''6-Tuple C{(zone, digraph, easting, northing, band, datum)},
       C{zone}, C{digraph} and C{band} as C{str}, C{easting} and
       C{northing} in C{meter} and C{datum} a L{Datum}.
    '''
    _Names_ = ('zone', 'digraph', 'easting', 'northing', 'band', 'datum')


class NearestOn3Tuple(_NamedTuple):  # .points.py, .sphericalTrigonometry.py
    '''3-Tuple C{(closest, distance, angle)} of the C{closest}
       point on the polygon, either a C{LatLon} instance or a
       L{LatLon3Tuple}C{(lat, lon, height)} and the C{distance}
       and C{angle} to the C{closest} point are in C{meter}
       respectively compass C{degrees360}.
    '''
    _Names_ = ('closest', 'distance', 'angle')


class NearestOn5Tuple(_NamedTuple):  # .points.py
    '''5-Tuple C{(lat, lon, distance, angle, height)} all in C{degrees},
       except C{height}.  The C{distance} is the L{equirectangular_}
       distance between the closest and the reference B{C{point}} in
       C{degrees}.  The C{angle} from the reference B{C{point}} to
       the closest point is in compass C{degrees360}, see function
       L{compassAngle}.  The C{height} is the (interpolated) height
       at the closest point in C{meter} or C{0}.
    '''
    _Names_ = ('lat', 'lon', 'distance', 'angle', 'height')


class Ned3Tuple(_NamedTuple):  # .ellipsoidalNvector.py
    '''3-Tuple C{(north, east, down)}, all in C{degrees}.
    '''
    _Names_ = ('north', 'east', 'down')


class PhiLam2Tuple(_NamedTuple):  # .frechet.py, .hausdorff.py, .latlonBase.py, .points.py, .vector3d.py
    '''2-Tuple C{(phi, lam)} with latitude C{phi} in C{radians[PI_2]}
       and longitude C{lam} in C{radians[PI]}.

       @note: Using C{phi/lambda} for lat-/longitude in radians
              follows Chris Veness' U{convention
              <https://www.Movable-Type.co.UK/scripts/latlong.html>}.
    '''
    _Names_ = ('phi', 'lam')

    def _3Tuple(self, height):
        '''(INTERNAL) Extend to a L{PhiLam3Tuple}.
        '''
        return self._xtend(PhiLam3Tuple, height)


class PhiLam3Tuple(_NamedTuple):  # .nvector.py
    '''3-Tuple C{(phi, lam, height)} with latitude C{phi} in
       C{radians[PI_2]}, longitude C{lam} in C{radians[PI]} and
       C{height} in C{meter}.

       @note: Using C{phi/lamda} for lat-/longitude in radians
              follows Chris Veness' U{convention
              <https://www.Movable-Type.co.UK/scripts/latlong.html>}.
    '''
    _Names_ = ('phi', 'lam', 'height')


class Point3Tuple(_NamedTuple):  # .points.py
    '''3-Tuple C{(x, y, ll)} in C{meter}, C{meter} and C{LatLon}.
    '''
    _Names_ = ('x', 'y', 'll')


class Points2Tuple(_NamedTuple):  # .formy.py
    '''2-Tuple C{(number, points)} with the C{number} of points
       and -possible reduced- C{list} or C{tuple} of C{points}.
    '''
    _Names_ = ('number', 'points')


class Shape2Tuple(_NamedTuple):  # .points.py
    '''2-Tuple C{(nrows, ncols)}, the number of rows and columns,
       both C{int}.
    '''
    _Names_ = ('nrows', 'ncols')


class UtmUps2Tuple(_NamedTuple):  # .epsg.py
    '''2-Tuple C{(zone, hemipole)} as C{int} and C{str}, where
       C{zone} is C{1..60} for UTM or C{0} for UPS and C{hemipole}
       C{'N'|'S'} is the UTM hemisphere or the UPS pole.
    '''
    _Names_ = ('zone', 'hemipole')


class UtmUps4Tuple(_NamedTuple):  # .mgrs.py
    '''4-Tuple C{(zone, hemipole, easting, northing)} as C{str},
       C{str}, C{meter} and C{meter}.
    '''
    _Names_ = ('zone', 'hemipole', 'easting', 'northing')


class UtmUps5Tuple(_NamedTuple):  # .ups.py, .utm.py, .utmups.py
    '''5-Tuple C{(zone, hemipole, easting, northing, band)} as C{int},
       C{str}, C{meter}, C{meter} and C{band} letter, where C{zone}
       is C{1..60} for UTM or C{0} for UPS, C{hemipole} C{'N'|'S'} is
       the UTM hemisphere or the UPS pole and {band} is C{""} or the
       (longitudinal) UTM band C{'C'|'D'..'W'|'X'} or the (polar) UPS
       band C{'A'|'B'|'Y'|'Z'}.
    '''
    _Names_ = ('zone', 'hemipole', 'easting', 'northing', 'band')


class UtmUps8Tuple(_NamedTuple):  # .ups.py, .utm.py, .utmups.py
    '''8-Tuple C{(zone, hemipole, easting, northing, band, datum,
       convergence, scale)} as C{int}, C{str}, C{meter}, C{meter},
       C{band} letter, C{Datum}, C{degrees} and C{float}, where
       C{zone} is C{1..60} for UTM or C{0} for UPS, C{hemipole}
       C{'N'|'S'} is the UTM hemisphere or the UPS pole and {band}
       is C{""} or the (longitudinal) UTM band C{'C'|'D'..'W'|'X'}
       or the (polar) UPS band C{'A'|'B'|'Y'|'Z'}.
    '''
    _Names_ = ('zone', 'hemipole', 'easting', 'northing',
               'band', 'datum', 'convergence', 'scale')


class UtmUpsLatLon5Tuple(_NamedTuple):  # .ups.py, .utm.py, .utmups.py
    '''5-Tuple C{(zone, band, hemipole, lat, lon)} as C{int},
       C{str}, C{str}, C{degrees90} and C{degrees180}, where
       C{zone} is C{1..60} for UTM or C{0} for UPS, {band} is
       C{""} or the (longitudinal) UTM band C{'C'|'D'..'W'|'X'}
       or (polar) UPS band C{'A'|'B'|'Y'|'Z'} and C{hemipole}
       C{'N'|'S'} is the UTM hemisphere or the UPS pole.
    '''
    _Names_ = ('zone', 'band', 'hemipole', 'lat', 'lon')


class Vector3Tuple(_NamedTuple):
    '''3-Tuple C{(x, y, z)} of (geocentric) components, all in
       C{meter} or C{units}.
    '''
    _Names_ = ('x', 'y', 'z')

    def _4Tuple(self, h):
        '''(INTERNAL) Extend to a L{Vector4Tuple}.
        '''
        return self._xtend(Vector4Tuple, h)


class Vector4Tuple(_NamedTuple):  # .nvector.py
    '''4-Tuple C{(x, y, z, h)} of (geocentric) components, all in
       C{meter} or C{units}.
    '''
    _Names_ = ('x', 'y', 'z', 'h')


def classname(inst, prefixed=None):
    '''Return the instance' class name optionally prefixed with the
       module name.

       @param inst: The object (any C{type}).
       @keyword prefixed: Include the module name (C{bool}), see
                          function L{classnaming}.

       @return: The B{C{inst}}'s C{[module.]class} name (C{str}).
    '''
    try:
        n = inst.__class__.__name__
    except AttributeError:
        n = 'Nn'
    if prefixed or (getattr(inst, 'classnaming', _Named._classnaming)
                    if prefixed is None else False):
        try:
            m = inst.__module__
            n = '.'.join(m.split('.')[-1:] + [n])
        except AttributeError:
            pass
    return n


def classnaming(prefixed=None):
    '''Get/set the default naming for C{[module.]class} names.

       @keyword prefixed: Include the module name (C{bool}).

       @return: Previous class naming setting (C{bool}).
    '''
    t = _Named._classnaming
    if prefixed in (True, False):
        _Named._classnaming = prefixed
    return t


def inStr(inst, *args, **kwds):
    '''Return the string representation of an instance.

       @param inst: The instance (any C{type}).
       @param args: Optional positional arguments.
       @keyword kwds: Optional keyword arguments.

       @return: The B{C{inst}}'s representation (C{str}).
    '''
    return unStr(classname(inst), *args, **kwds)


def nameof(inst):
    '''Get the name of an instance.

       @param inst: The object (any C{type}).

       @return: The instance' name (C{str}) or C{""}.
    '''
    return getattr(inst, 'name', '')


if __name__ == '__main__':

    from pygeodesy.lazily import _FOR_DOCS
    if not _FOR_DOCS:
        raise NameError('usage: %s' % ('env PYGEODESY_FOR_DOCS=1 python ...',))

    ls = set(locals().keys()) - \
         set(('fStr', 'ls', 'n', 'property_RO', 'issubclassof', 'unStr'))
    for n in __all__:
        if n not in ls:
            raise NameError('%s %r not in %s' % ('__all__', n, 'locals'))
    for n in ls:
        if n not in __all__ and not n.startswith('_'):
            raise NameError('%s %r not in %s' % ('locals', n, '__all__'))

    print('all %s %s OK' % (len(ls), 'locals'))

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
