
# -*- coding: utf-8 -*-

u'''(INTERNAL) Nameable class instances.

Classes C{_Named}, C{_NamedDict}, C{_NamedEnum}, C{_NamedEnumItem} and
C{_NamedTuple} and several subclasses thereof, all with nameable instances.

The items in a C{_NamedDict} are accessable as attributes and the items
in a C{_NamedTuple} are named to be accessable as attributes, similar to
standard Python C{namedtuple}s.

@see: Module L{pygeodesy.namedTuples} for (most of) the C{Named-Tuples}.
'''

from pygeodesy.basics import isclass, isidentifier, iskeyword, isstr, \
                             issubclassof, len2, _xcopy, _xdup, _zip
from pygeodesy.errors import _AssertionError, _AttributeError, _incompatible, \
                             _IndexError, _IsnotError, itemsorted, LenError, \
                             _NameError, _NotImplementedError, _TypeError, \
                             _TypesError, _ValueError, UnitError, \
                             _xkwds, _xkwds_popitem
from pygeodesy.interns import NN, _at_, _AT_, _COLON_, _COLONSPACE_, _COMMA_, \
                             _COMMASPACE_, _doesn_t_exist_, _DOT_, _DUNDER_, \
                             _EQUAL_, _EQUALSPACED_, _exists_, _immutable_, _name_, \
                             _NL_, _NN_, _not_, _other_, _s_, _SPACE_, _std_, \
                             _UNDER_, _valid_, _vs_,  _dunder_nameof, _UNDER
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _caller3, _getenv
from pygeodesy.props import _allPropertiesOf_n, deprecated_method, Property_RO, \
                            _hasProperty, property_doc_, property_RO, \
                            _update_all, _update_attrs
from pygeodesy.streprs import attrs, Fmt, lrstrip, pairs, reprs, unstr

__all__ = _ALL_LAZY.named
__version__ = '22.10.30'

_COMMANL_           = _COMMA_ + _NL_
_COMMASPACEDOT_     = _COMMASPACE_ + _DOT_
_del_               = 'del'
_item_              = 'item'
_MRO_               = 'MRO'
# __DUNDER gets mangled in class
_name               = _UNDER(_name_)
_Names_             = '_Names_'
_registered_        = 'registered'  # PYCHOK used!
_std_NotImplemented = _getenv('PYGEODESY_NOTIMPLEMENTED', NN).lower() == _std_
_Units_             = '_Units_'


def _xjoined_(prefix, name):
    '''(INTERNAL) Join C{pref} and non-empty C{name}.
    '''
    return _SPACE_(prefix, repr(name)) if name and prefix else (prefix or name)


def _xnamed(inst, name, force=False):
    '''(INTERNAL) Set the instance' C{.name = B{name}}.

       @arg inst: The instance (C{_Named}).
       @arg name: The name (C{str}).
       @kwarg force: Force name change (C{bool}).

       @return: The B{C{inst}}, named if B{C{force}}d or
                not named before.
    '''
    if name and isinstance(inst, _Named):
        if not inst.name:
            inst.name = name
        elif force:
            inst.rename(name)
    return inst


def _xother3(inst, other, name=_other_, up=1, **name_other):
    '''(INTERNAL) Get C{name} and C{up} for a named C{other}.
    '''
    if name_other:  # and not other and len(name_other) == 1
        name, other = _xkwds_popitem(name_other)
    elif other and len(other) == 1:
        other = other[0]
    else:
        raise _AssertionError(name, other, txt=classname(inst, prefixed=True))
    return other, name, up


def _xotherError(inst, other, name=_other_, up=1):
    '''(INTERNAL) Return a C{_TypeError} for an incompatible, named C{other}.
    '''
    n = _callname(name, classname(inst, prefixed=True), inst.name, up=up + 1)
    return _TypeError(name, other, txt=_incompatible(n))


def _xvalid(name, _OK=False):
    '''(INTERNAL) Check valid attribute name C{name}.
    '''
    return True if (name and isstr(name)
                         and name != _name_
                         and (_OK or not name.startswith(_UNDER_))
                         and (not iskeyword(name))
                         and isidentifier(name)) else False


class _Dict(dict):
    '''(INTERNAL) An C{dict} with both key I{and}
       attribute access to the C{dict} items.
    '''
    def __getattr__(self, name):
        '''Get an attribute or item by B{C{name}}.
        '''
        try:
            return self[name]
        except KeyError:
            pass
        name = _DOT_(self.__class__.__name__, name)
        raise _AttributeError(item=name, txt=_doesn_t_exist_)

    def __repr__(self):
        '''Default C{repr(self)}.
        '''
        return self.toRepr()

    def __str__(self):
        '''Default C{str(self)}.
        '''
        return self.toStr()

    def set_(self, **items):  # PYCHOK signature
        '''Add one or several new items or replace existing ones.

           @kwarg items: One or more C{name=value} pairs.
        '''
        dict.update(self, items)

    def toRepr(self, **prec_fmt):  # PYCHOK signature
        '''Like C{repr(dict)} but with C{name} prefix and with
           C{floats} formatted by function L{pygeodesy.fstr}.
        '''
        n = self.get(_name_, classname(self))
        return Fmt.PAREN(n, self._toT(_EQUAL_, **prec_fmt))

    def toStr(self, **prec_fmt):  # PYCHOK signature
        '''Like C{str(dict)} but with C{floats} formatted by
           function L{pygeodesy.fstr}.
        '''
        return Fmt.CURLY(self._toT(_COLONSPACE_, **prec_fmt))

    def _toT(self, sep, **kwds):
        '''(INTERNAL) Helper for C{.toRepr} and C{.toStr}, also
           in C{_NamedDict} below.
        '''
        kwds = _xkwds(kwds, prec=6, fmt=Fmt.F, sep=sep)
        return _COMMASPACE_.join(pairs(itemsorted(self), **kwds))


class _Named(object):
    '''(INTERNAL) Root class for named objects.
    '''
    _iteration   = None   # iteration number (C{int}) or C{None}
    _name        = NN     # name (C{str})
    _classnaming = False  # prefixed (C{bool})
#   _updates     = 0      # OBSOLETE Property/property updates

    def __imatmul__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)  # PYCHOK Python 3.5+

    def __matmul__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)  # PYCHOK Python 3.5+

    def __repr__(self):
        '''Default C{repr(self)}.
        '''
        return Fmt.ANGLE(_SPACE_(self, _at_, hex(id(self))))

    def __rmatmul__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)  # PYCHOK Python 3.5+

    def __str__(self):
        '''Default C{str(self)}.
        '''
        return self.named2

    def attrs(self, *names, **sep_COMMASPACE__Nones_True__pairs_kwds):
        '''Join named attributes as I{name=value} strings, with C{float}s formatted by
           function L{pygeodesy.fstr}.

           @arg names: The attribute names, all positional (C{str}).
           @kwarg sep_COMMASPACE__Nones_True__pairs_kwds: Keyword argument for function
                  L{pygeodesy.pairs}, except C{B{sep}=", "} and C{B{Nones}=True} to
                  in- or exclude missing or C{None}-valued attributes.

           @return: All C{name=value} pairs, joined by B{C{sep}} (C{str}).

           @see: Functions L{pygeodesy.attrs}, L{pygeodesy.fstr} and L{pygeodesy.pairs}.
        '''
        def _sep_kwds(sep=_COMMASPACE_, **kwds):
            return sep, kwds

        sep, kwds = _sep_kwds(**sep_COMMASPACE__Nones_True__pairs_kwds)
        return sep.join(attrs(self, *names, **kwds))

    @Property_RO
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
        '''Set the class naming for C{[module.].class} names (C{bool}),
           to C{True} to include the module name.
        '''
        b = bool(prefixed)
        if self._classnaming != b:
            self._classnaming = b
            _update_attrs(self, *_Named_Property_ROs)

    def classof(self, *args, **kwds):
        '''Create another instance of this very class.

           @arg args: Optional, positional arguments.
           @kwarg kwds: Optional, keyword arguments.

           @return: New instance (B{self.__class__}).
        '''
        return _xnamed(self.__class__(*args, **kwds), self.name)

    def copy(self, deep=False, name=NN):
        '''Make a shallow or deep copy of this instance.

           @kwarg deep: If C{True} make a deep, otherwise
                        a shallow copy (C{bool}).
           @kwarg name: Optional, non-empty name (C{str}).

           @return: The copy (C{This class} or sub-class thereof).
        '''
        c = _xcopy(self, deep=deep)
        if name:
            c.rename(name)
        return c

    def _DOT_(self, *names):
        '''(INTERNAL) Period-join C{self.name} and C{names}.
        '''
        return _DOT_(self.name, *names)

    def dup(self, name=NN, **items):
        '''Duplicate this instance, replacing some items.

           @kwarg name: Optional, non-empty name (C{str}).
           @kwarg items: Attributes to be changed (C{any}).

           @return: The duplicate (C{This class} or sub-class thereof).

           @raise AttributeError: Some B{C{items}} invalid.
        '''
        d = _xdup(self, **items)
        if name:
            d.rename(name=name)
        return d

    def _instr(self, name, prec, *attrs, **props_kwds):
        '''(INTERNAL) Format, used by C{Conic}, C{Ellipsoid}, C{Transform}, C{Triaxial}.
        '''
        def _props_kwds(props=(), **kwds):
            return props, kwds

        t = Fmt.EQUAL(_name_, repr(name or self.name)),
        if attrs:
            t += pairs(((a, getattr(self, a)) for a in attrs),
                       prec=prec, ints=True)
        props, kwds =_props_kwds(**props_kwds)
        if props:
            t += pairs(((p.name, getattr(self, p.name)) for p in props),
                       prec=prec, ints=True)
        if kwds:
            t += pairs(kwds, prec=prec)
        return _COMMASPACE_.join(t[1:] if name is None else t)

    @property_RO
    def iteration(self):  # see .karney.GDict
        '''Get the most recent iteration number (C{int}) or C{None}
           if not available or not applicable.

           @note: The interation number may be an aggregate number over
                  several, nested functions.
        '''
        return self._iteration

    def methodname(self, which):
        '''Get a method C{[module.]class.method} name of this object (C{str}).

           @arg which: The method (C{callable}).
        '''
        return _DOT_(self.classname, which.__name__ if callable(which) else _NN_)

    @property_doc_(''' the name (C{str}).''')
    def name(self):
        '''Get the name (C{str}).
        '''
        return self._name

    @name.setter  # PYCHOK setter!
    def name(self, name):
        '''Set the name (C{str}).

           @raise NameError: Can't rename, use method L{rename}.
        '''
        m, n = self._name, str(name)
        if not m:
            self._name = n
        elif n != m:
            n =  repr(n)
            c =  self.classname
            t = _DOT_(c, Fmt.PAREN(self.rename.__name__, n))
            m =  Fmt.PAREN(_SPACE_('was', repr(m)))
            n = _DOT_(c, _EQUALSPACED_(_name_, n))
            n = _SPACE_(n, m)
            raise _NameError(_SPACE_('use', t), txt=_not_(n))
        # to set the name from a sub-class, use
        #  self.name = name or
        # _Named.name.fset(self, name), but NOT
        # _Named(self).name = name

    @Property_RO
    def named(self):
        '''Get the name I{or} class name or C{""} (C{str}).
        '''
        return self.name or self.classname

    @Property_RO
    def named2(self):
        '''Get the C{class} name I{and/or} the name or C{""} (C{str}).
        '''
        return _xjoined_(self.classname, self.name)

    @Property_RO
    def named3(self):
        '''Get the I{prefixed} C{class} name I{and/or} the name or C{""} (C{str}).
        '''
        return _xjoined_(classname(self, prefixed=True), self.name)

    @Property_RO
    def named4(self):
        '''Get the C{package.module.class} name I{and/or} the name or C{""} (C{str}).
        '''
        return _xjoined_(_DOT_(self.__module__, self.__class__.__name__),  self.name)

    def rename(self, name):
        '''Change the name.

           @arg name: The new name (C{str}).

           @return: The old name (C{str}).
        '''
        m, n = self._name, str(name)
        if n != m:
            self._name = n
            _update_attrs(self, *_Named_Property_ROs)
        return m

    def toRepr(self, **unused):  # PYCHOK no cover
        '''Default C{repr(self)}.
        '''
        return repr(self)

    def toStr(self, **unused):  # PYCHOK no cover
        '''Default C{str(self)}.
        '''
        return str(self)

    @deprecated_method
    def toStr2(self, **kwds):  # PYCHOK no cover
        '''DEPRECATED, use method C{toRepr}.'''
        return self.toRepr(**kwds)

    def _unstr(self, which, *args, **kwds):
        '''(INTERNAL) Return the string representation of a method
           invokation of this instance: C{str(self).method(...)}

           @see: Function L{pygeodesy.unstr}.
        '''
        n = _DOT_(self, which.__name__ if callable(which) else _NN_)
        return unstr(n, *args, **kwds)

    def _xnamed(self, inst, name=NN, force=False):
        '''(INTERNAL) Set the instance' C{.name = self.name}.

           @arg inst: The instance (C{_Named}).
           @kwarg name: Optional name, overriding C{self.name} (C{str}).
           @kwarg force: Force name change (C{bool}).

           @return: The B{C{inst}}, named if not named before.
        '''
        return _xnamed(inst, name or self.name, force=force)

    def _xrenamed(self, inst):
        '''(INTERNAL) Rename the instance' C{.name = self.name}.

           @arg inst: The instance (C{_Named}).

           @return: The B{C{inst}}, named if not named before.

           @raise TypeError: Not C{isinstance(B{inst}, _Named)}.
        '''
        if not isinstance(inst, _Named):
            raise _IsnotError(_valid_, inst=inst)

        inst.rename(self.name)
        return inst

_Named_Property_ROs = _allPropertiesOf_n(5, _Named, Property_RO)  # PYCHOK once


class _NamedBase(_Named):
    '''(INTERNAL) Base class with name.
    '''
    def __repr__(self):
        '''Default C{repr(self)}.
        '''
        return self.toRepr()

    def __str__(self):
        '''Default C{str(self)}.
        '''
        return self.toStr()

#   def notImplemented(self, attr):
#       '''Raise error for a missing method, function or attribute.
#
#          @arg attr: Attribute name (C{str}).
#
#          @raise NotImplementedError: No such attribute.
#       '''
#       c = self.__class__.__name__
#       return NotImplementedError(_DOT_(c, attr))

    def others(self, *other, **name_other_up):  # see .points.LatLon_.others
        '''Refined class comparison, invoked as C{.others(other=other)},
           C{.others(name=other)} or C{.others(other, name='other')}.

           @arg other: The other instance (any C{type}).
           @kwarg name_other_up: Overriding C{name=other} and C{up=1}
                                 keyword arguments.

           @return: The B{C{other}} iff compatible with this instance's
                    C{class} or C{type}.

           @raise TypeError: Mismatch of the B{C{other}} and this
                             instance's C{class} or C{type}.
        '''
        if other:  # most common, just one arg B{C{other}}
            other0 = other[0]
            if isinstance(other0, self.__class__) or \
               isinstance(self, other0.__class__):
                return other0

        other, name, up = _xother3(self, other, **name_other_up)
        if isinstance(self, other.__class__) or \
           isinstance(other, self.__class__):
            return _xnamed(other, name)

        raise _xotherError(self, other, name=name, up=up + 1)

    def toRepr(self, **kwds):  # PYCHOK expected
        '''(INTERNAL) I{Could be overloaded}.

           @kwarg kwds: Optional, C{toStr} keyword arguments.

           @return: C{toStr}() with keyword arguments (as C{str}).
        '''
        t = lrstrip(self.toStr(**kwds))
#       if self.name:
#           t =  NN(Fmt.EQUAL(name=repr(self.name)), sep, t)
        return Fmt.PAREN(self.classname, t)  # XXX (self.named, t)

#   def toRepr(self, **kwds)
#       if kwds:
#           s = NN.join(reprs((self,), **kwds))
#       else:  # super().__repr__ only for Python 3+
#           s = super(self.__class__, self).__repr__()
#       return Fmt.PAREN(self.named, s)  # clips(s)

    def toStr(self, **kwds):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}, see function C{notOverloaded}.
        '''
        notOverloaded(self, **kwds)

#   def toStr(self, **kwds):
#       if kwds:
#           s = NN.join(strs((self,), **kwds))
#       else:  # super().__str__ only for Python 3+
#           s = super(self.__class__, self).__str__()
#       return s

    def _update(self, updated, *attrs, **setters):
        '''(INTERNAL) Zap cached instance attributes and overwrite C{__dict__} or L{Property_RO} values.
        '''
        u = _update_all(self, *attrs) if updated else 0
        if setters:
            d = self.__dict__
            # double-check that setters are Property_RO's
            for n, v in setters.items():
                if n in d or _hasProperty(self, n, Property_RO):
                    d[n] = v
                else:
                    raise _AssertionError(n, v, txt=repr(self))
            u += len(setters)
        return u


class _NamedDict(_Dict, _Named):
    '''(INTERNAL) Named C{dict} with key I{and} attribute
       access to the items.
    '''
    def __init__(self, *args, **kwds):
        if args:  # args override kwds
            if len(args) != 1:
                t = unstr(self.classname, *args, **kwds)  # PYCHOK no cover
                raise _ValueError(args=len(args), txt=t)
            kwds = _xkwds(dict(args[0]), **kwds)
        if _name_ in kwds:
            _Named.name.fset(self, kwds.pop(_name_))  # see _Named.name
        _Dict.__init__(self, kwds)

    def __delattr__(self, name):
        '''Delete an attribute or item by B{C{name}}.
        '''
        if name in _Dict.keys(self):
            _Dict.pop(name)
        elif name in (_name_, _name):
            # _Dict.__setattr__(self, name, NN)
            _Named.rename(self, NN)
        else:
            _Dict.__delattr__(self, name)

    def __getattr__(self, name):
        '''Get an attribute or item by B{C{name}}.
        '''
        try:
            return self[name]
        except KeyError:
            if name == _name_:
                return _Named.name.fget(self)
        raise _AttributeError(item=self._DOT_(name), txt=_doesn_t_exist_)

    def __getitem__(self, key):
        '''Get the value of an item by B{C{key}}.
        '''
        if key == _name_:
            raise KeyError(Fmt.SQUARE(self.classname, key))
        return _Dict.__getitem__(self, key)

    def __setattr__(self, name, value):
        '''Set attribute or item B{C{name}} to B{C{value}}.
        '''
        if name in _Dict.keys(self):
            _Dict.__setitem__(self, name, value)  # self[name] = value
        else:
            _Dict.__setattr__(self, name, value)

    def __setitem__(self, key, value):
        '''Set item B{C{key}} to B{C{value}}.
        '''
        if key == _name_:
            raise KeyError(_EQUAL_(Fmt.SQUARE(self.classname, key), repr(value)))
        _Dict.__setitem__(self, key, value)

    def toRepr(self, **prec_fmt):  # PYCHOK signature
        '''Like C{repr(dict)} but with C{name} prefix and with
           C{floats} formatted by function L{pygeodesy.fstr}.
        '''
        return Fmt.PAREN(self.name, self._toT(_EQUAL_, **prec_fmt))

    def toStr(self, **prec_fmt):  # PYCHOK signature
        '''Like C{str(dict)} but with C{floats} formatted by
           function L{pygeodesy.fstr}.
        '''
        return Fmt.CURLY(self._toT(_COLONSPACE_, **prec_fmt))


class _NamedEnum(_NamedDict):
    '''(INTERNAL) Enum-like C{_NamedDict} with attribute access
       restricted to valid keys.
    '''
    _item_Classes = ()

    def __init__(self, Class, *Classes, **name):
        '''New C{_NamedEnum}.

           @arg Class: Initial class or type acceptable as items
                       values (C{type}).
           @arg Classes: Additional, acceptable classes or C{type}s.
        '''
        self._item_Classes = (Class,) + Classes
        n = name.get(_name_, NN) or NN(Class.__name__, _s_)
        if n and _xvalid(n, _OK=True):
            _Named.name.fset(self, n)  # see _Named.name

    def __getattr__(self, name):
        '''Get the value of an attribute or item by B{C{name}}.
        '''
        try:
            return self[name]
        except KeyError:
            if name == _name_:
                return _NamedDict.name.fget(self)
        raise _AttributeError(item=self._DOT_(name), txt=_doesn_t_exist_)

    def __repr__(self):
        '''Default C{repr(self)}.
        '''
        return self.toRepr()

    def __str__(self):
        '''Default C{str(self)}.
        '''
        return self.toStr()

    def _assert(self, **kwds):
        '''(INTERNAL) Check attribute name against given, registered name.
        '''
        for n, v in kwds.items():
            if isinstance(v, _LazyNamedEnumItem):  # property
                assert n is v.name
                # assert not hasattr(self.__class__, n)
                setattr(self.__class__, n, v)
            elif isinstance(v, self._item_Classes):  # PYCHOK no cover
                assert self[n] is v and getattr(self, n) \
                                    and self.find(v) == n
            else:
                raise _TypeError(v, name=n)

    def find(self, item, dflt=None):
        '''Find a registered item.

           @arg item: The item to look for (any C{type}).
           @kwarg dflt: Value to return if not found (any C{type}).

           @return: The B{C{item}}'s name if found (C{str}), or C{{dflt}} if
                    there is no such I{registered} B{C{item}}.
        '''
        for k, v in _Dict.items(self):  # XXX not self.items()
            if v is item:
                return k
        return dflt

    def get(self, name, dflt=None):
        '''Get the value of a I{registered} item.

           @arg name: The name of the item (C{str}).
           @kwarg dflt: Value to return (any C{type}).

           @return: The item with B{C{name}} if found, or B{C{dflt}} if
                    there is no item I{registered} with that B{C{name}}.
        '''
        # getattr needed to instantiate L{_LazyNamedEnumItem}
        return getattr(self, name, dflt)

    def items(self, all=False, asorted=False):
        '''Yield all or only the I{registered} items.

           @kwarg all: Use C{True} to yield {all} items or C{False}
                       for only the currently I{registered} ones.
           @kwarg asorted: If C{True}, yield the items sorted in
                           I{alphabetical, case-insensitive} order.
        '''
        if all:  # instantiate any remaining L{_LazyNamedEnumItem} ...
            # ... and remove the L{_LazyNamedEnumItem} from the class
            for n in tuple(n for n, p in self.__class__.__dict__.items()
                                      if isinstance(p, _LazyNamedEnumItem)):
                _ = getattr(self, n)

        return itemsorted(_Dict, self) if asorted else _Dict.items(self)

    def keys(self, **all_asorted):
        '''Yield the keys (C{str}) of all or only the I{registered} items,
           optionally sorted I{alphabetically} and I{case-insensitively}.

           @kwarg all_asorted: See method C{items}..
        '''
        for k, _ in self.items(**all_asorted):
            yield k

    def popitem(self):
        '''Remove I{an, any} curretly I{registed} item.

           @return: The removed item.
        '''
        return self._zapitem(*_Dict.pop(self))

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
            if not (n and isstr(n) and isidentifier(n)):
                raise ValueError
        except (AttributeError, ValueError, TypeError) as x:
            raise _NameError(_DOT_(_item_, _name_), item, cause=x)
        if n in self:
            raise _NameError(self._DOT_(n), item, txt=_exists_)
        if not (self._item_Classes and isinstance(item, self._item_Classes)):
            raise _TypesError(self._DOT_(n), item, *self._item_Classes)
        self[n] = item

    def unregister(self, name_or_item):
        '''Remove a I{registered} item.

           @arg name_or_item: Name (C{str}) or the item (any C{type}).

           @return: The unregistered item.

           @raise NameError: No item with that B{C{name}}.

           @raise ValueError: No such item.
        '''
        if isstr(name_or_item):
            name = name_or_item
        else:
            name = self.find(name_or_item)
        try:
            item = _Dict.pop(self, name)
        except KeyError:
            raise _NameError(item=self._DOT_(name), txt=_doesn_t_exist_)
        return self._zapitem(name, item)

    pop = unregister

    def toRepr(self, prec=6, fmt=Fmt.F, sep=_COMMANL_, **all_asorted):  # PYCHOK _NamedDict
        '''Like C{repr(dict)} but C{name}s optionally sorted and
           C{floats} formatted by function L{pygeodesy.fstr}.
        '''
        t = ((self._DOT_(n), v) for n, v in self.items(**all_asorted))
        return sep.join(pairs(t, prec=prec, fmt=fmt, sep=_COLONSPACE_))

    def toStr(self, *unused, **all_asorted):  # PYCHOK _NamedDict
        '''Return a string with all C{name}s, optionally sorted.
        '''
        return self._DOT_(_COMMASPACEDOT_.join(self.keys(**all_asorted)))

    def values(self, **all_asorted):
        '''Yield the value (C{type}) of all or only the I{registered} items,
           optionally sorted I{alphabetically} and I{case-insensitively}.

           @kwarg all_asorted: See method C{items}.
        '''
        for _, v in self.items(**all_asorted):
            yield v

    def _zapitem(self, name, item):
        # remove _LazyNamedEnumItem property value if still present
        if self.__dict__.get(name, None) is item:
            self.__dict__.pop(name)  # [name] = None
        item._enum = None
        return item


class _LazyNamedEnumItem(property_RO):  # XXX or descriptor?
    '''(INTERNAL) Lazily instantiated L{_NamedEnumItem}.
    '''
    pass


def _lazyNamedEnumItem(name, *args, **kwds):
    '''(INTERNAL) L{_LazyNamedEnumItem} property-like factory.

       @see: Luciano Ramalho, "Fluent Python", page 636, O'Reilly, 2016,
             "Coding a Property Factory", especially Example 19-24.
    '''
    def _fget(inst):
        # assert isinstance(inst, _NamedEnum)
        try:  # get the item from the instance' __dict__
            item = inst.__dict__[name]
        except KeyError:
            # instantiate an _NamedEnumItem, it self-registers
            item = inst._Lazy(*args, **_xkwds(kwds, name=name))
            # assert inst[name] is item  # MUST be registered
            # store the item in the instance' __dict__
            inst.__dict__[name] = item
            # remove the property from the registry class, such that
            # (a) the property no longer overrides the instance' item
            # in inst.__dict__ and (b) _NamedEnum.items(all=True) only
            # sees un-instantiated ones to be instantiated
            p = getattr(inst.__class__, name, None)
            if isinstance(p, _LazyNamedEnumItem):
                delattr(inst.__class__, name)
        # assert isinstance(item, _NamedEnumItem)
        return item

    p = _LazyNamedEnumItem(_fget)
    p.name = name
    return p


class _NamedEnumItem(_NamedBase):
    '''(INTERNAL) Base class for items in a C{_NamedEnum} registery.
    '''
    _enum = None

#   def __ne__(self, other):  # XXX fails for Lcc.conic = conic!
#       '''Compare this and an other item.
#
#          @return: C{True} if different, C{False} otherwise.
#       '''
#       return not self.__eq__(other)

    @property_doc_(''' the I{registered} name (C{str}).''')
    def name(self):
        '''Get the I{registered} name (C{str}).
        '''
        return self._name

    @name.setter  # PYCHOK setter!
    def name(self, name):
        '''Set the name, unless already registered (C{str}).
        '''
        if self._enum:
            raise _NameError(str(name), self, txt=_registered_)  # XXX _TypeError
        self._name = str(name)

    def _register(self, enum, name):
        '''(INTERNAL) Add this item as B{C{enum.name}}.

           @note: Don't register if name is empty or doesn't
                  start with a letter.
        '''
        if name and _xvalid(name, _OK=True):
            self.name = name
            if name[:1].isalpha():  # '_...' not registered
                enum.register(self)
                self._enum = enum

    def unregister(self):
        '''Remove this instance from its C{_NamedEnum} registry.

           @raise AssertionError: Mismatch of this and registered item.

           @raise NameError: This item is unregistered.
        '''
        enum = self._enum
        if enum and self.name and self.name in enum:
            item = enum.unregister(self.name)
            if item is not self:
                t = _SPACE_(repr(item), _vs_, repr(self))  # PYCHOK no cover
                raise _AssertionError(t)


class _NamedTuple(tuple, _Named):
    '''(INTERNAL) Base for named C{tuple}s with both index I{and}
       attribute name access to the items.

       @note: This class is similar to Python's C{namedtuple},
              but statically defined, lighter and limited.
    '''
    _Names_ = ()  # item names, non-identifier, no leading underscore
    '''Tuple specifying the C{name} of each C{Named-Tuple} item.

       @note: Specify at least 2 item names.
    '''
    _Units_    = ()    # .units classes
    '''Tuple defining the C{units} of the value of each C{Named-Tuple} item.

       @note: The C{len(_Units_)} must match C{len(_Names_)}.
    '''
    _validated = False  # set to True I{per sub-class!}

    def __new__(cls, arg, *args, **iteration_name):
        '''New L{_NamedTuple} initialized with B{C{positional}} arguments.

           @arg arg: Tuple items (C{tuple}, C{list}, ...) or first tuple
                     item of several more in other positional arguments.
           @arg args: Tuple items (C{any}), all positional arguments.
           @kwarg iteration_name: Only keyword arguments C{B{iteration}=None}
                                  and C{B{name}=NN} are used, any other are
                                  I{silently} ignored.

           @raise LenError: Unequal number of positional arguments and
                            number of item C{_Names_} or C{_Units_}.

           @raise TypeError: The C{_Names_} or C{_Units_} attribute is
                             not a C{tuple} of at least 2 items.

           @raise ValueError: Item name is not a C{str} or valid C{identifier}
                              or starts with C{underscore}.
        '''
        n, args = len2(((arg,) + args) if args else arg)
        self = tuple.__new__(cls, args)
        if not self._validated:
            self._validate()

        N = len(self._Names_)
        if n != N:
            raise LenError(self.__class__, args=n, _Names_=N)

        if iteration_name:
            self._kwdself(**iteration_name)
        return self

    def __delattr__(self, name):
        '''Delete an attribute by B{C{name}}.

           @note: Items can not be deleted.
        '''
        if name in self._Names_:
            raise _TypeError(_del_, _DOT_(self.classname, name), txt=_immutable_)
        elif name in (_name_, _name):
            _Named.__setattr__(self, name, NN)  # XXX _Named.name.fset(self, NN)
        else:
            tuple.__delattr__(self, name)

    def __getattr__(self, name):
        '''Get the value of an attribute or item by B{C{name}}.
        '''
        try:
            return tuple.__getitem__(self, self._Names_.index(name))
        except IndexError:
            raise _IndexError(_DOT_(self.classname, Fmt.ANGLE(_name_)), name)
        except ValueError:  # e.g. _iteration
            return tuple.__getattribute__(self, name)

#   def __getitem__(self, index):  # index, slice, etc.
#       '''Get the item(s) at an B{C{index}} or slice.
#       '''
#       return tuple.__getitem__(self, index)

    def __repr__(self):
        '''Default C{repr(self)}.
        '''
        return self.toRepr()

    def __setattr__(self, name, value):
        '''Set attribute or item B{C{name}} to B{C{value}}.
        '''
        if name in self._Names_:
            raise _TypeError(_DOT_(self.classname, name), value, txt=_immutable_)
        elif name in (_name_, _name):
            _Named.__setattr__(self, name, value)  # XXX _Named.name.fset(self, value)
        else:  # e.g. _iteration
            tuple.__setattr__(self, name, value)

    def __str__(self):
        '''Default C{repr(self)}.
        '''
        return self.toStr()

    def dup(self, name=NN, **items):
        '''Duplicate this tuple replacing one or more items.

           @kwarg name: Optional new name (C{str}).
           @kwarg items: Items to be replaced (C{name=value} pairs), if any.

           @return: A copy of this tuple with B{C{items}}.

           @raise NameError: Some B{C{items}} invalid.
        '''
        tl = list(self)
        if items:
            _ix = self._Names_.index
            try:
                for n, v in items.items():
                    tl[_ix(n)] = v
            except ValueError:  # bad item name
                raise _NameError(_DOT_(self.classname, n), v, this=self)
        return self.classof(*tl, name=name or self.name)

    def items(self):
        '''Yield the items, each as a C{(name, value)} pair (C{2-tuple}).

           @see: Method C{.units}.
        '''
        for n, v in _zip(self._Names_, self):  # strict=True
            yield n, v

    iteritems = items

    def _kwdself(self, iteration=None, name=NN, **unused):
        '''(INTERNAL) Set C{__new__} keyword arguments.
        '''
        if iteration is not None:
            self._iteration = iteration
        if name:
            self.name = name

    def toRepr(self, prec=6, sep=_COMMASPACE_, fmt=Fmt.F, **unused):  # PYCHOK signature
        '''Return this C{Named-Tuple} items as C{name=value} string(s).

           @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                        Trailing zero decimals are stripped for B{C{prec}} values
                        of 1 and above, but kept for negative B{C{prec}} values.
           @kwarg sep: Separator to join (C{str}).
           @kwarg fmt: Optional, C{float} format (C{str}).

           @return: Tuple items (C{str}).
        '''
        t = pairs(self.items(), prec=prec, fmt=fmt)
#       if self.name:
#           t = (Fmt.EQUAL(name=repr(self.name)),) + t
        return Fmt.PAREN(self.named, sep.join(t))  # XXX (self.classname, sep.join(t))

    def toStr(self, prec=6, sep=_COMMASPACE_, fmt=Fmt.F, **unused):  # PYCHOK signature
        '''Return this C{Named-Tuple} items as string(s).

           @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                        Trailing zero decimals are stripped for B{C{prec}} values
                        of 1 and above, but kept for negative B{C{prec}} values.
           @kwarg sep: Separator to join (C{str}).
           @kwarg fmt: Optional C{float} format (C{str}).

           @return: Tuple items (C{str}).
        '''
        return Fmt.PAREN(sep.join(reprs(self, prec=prec, fmt=fmt)))

    def toUnits(self, Error=UnitError):  # overloaded in .frechet, .hausdorff
        '''Return a copy of this C{Named-Tuple} with each item value wrapped
           as an instance of its L{units} class.

           @kwarg Error: Error to raise for L{units} issues (C{UnitError}).

           @return: A duplicate of this C{Named-Tuple} (C{C{Named-Tuple}}).

           @raise Error: Invalid C{Named-Tuple} item or L{units} class.
        '''
        t = (v for _, v in self.units(Error=Error))
        return self.classof(*tuple(t))

    def units(self, Error=UnitError):
        '''Yield the items, each as a C{(name, value}) pair (C{2-tuple}) with
           the value wrapped as an instance of its L{units} class.

           @kwarg Error: Error to raise for L{units} issues (C{UnitError}).

           @raise Error: Invalid C{Named-Tuple} item or L{units} class.

           @see: Method C{.items}.
        '''
        for n, v, U in _zip(self._Names_, self, self._Units_):  # strict=True
            if not (v is None or U is None
                              or (isclass(U) and
                                  isinstance(v, U) and
                                  hasattr(v, _name_) and
                                  v.name == n)):  # PYCHOK indent
                v = U(v, name=n, Error=Error)
            yield n, v

    iterunits = units

    def _validate(self, _OK=False):  # see .EcefMatrix
        '''(INTERNAL) One-time check of C{_Names_} and C{_Units_}
           for each C{_NamedUnit} I{sub-class separately}.
        '''
        ns = self._Names_
        if not (isinstance(ns, tuple) and len(ns) > 1):  # XXX > 0
            raise _TypeError(_DOT_(self.classname, _Names_), ns)
        for i, n in enumerate(ns):
            if not _xvalid(n, _OK=_OK):
                t = Fmt.SQUARE(_Names_=i)  # PYCHOK no cover
                raise _ValueError(_DOT_(self.classname, t), n)

        us = self._Units_
        if not isinstance(us, tuple):
            raise _TypeError(_DOT_(self.classname, _Units_), us)
        if len(us) != len(ns):
            raise LenError(self.__class__, _Units_=len(us), _Names_=len(ns))
        for i, u in enumerate(us):
            if not (u is None or callable(u)):
                t = Fmt.SQUARE(_Units_=i)  # PYCHOK no cover
                raise _TypeError(_DOT_(self.classname, t), u)

        self.__class__._validated = True

    def _xtend(self, xTuple, *items):
        '''(INTERNAL) Extend this C{Named-Tuple} with C{items} to an other B{C{xTuple}}.
        '''
        if (issubclassof(xTuple, _NamedTuple) and
            (len(self._Names_) + len(items)) == len(xTuple._Names_) and
                 self._Names_ == xTuple._Names_[:len(self)]):
            return self._xnamed(xTuple(self + items))  # *(self + items)
        c = NN(self.classname,  repr(self._Names_))  # PYCHOK no cover
        x = NN(xTuple.__name__, repr(xTuple._Names_))  # PYCHOK no cover
        raise TypeError(_SPACE_(c, _vs_, x))


def callername(up=1, dflt=NN, source=False, underOK=False):
    '''Get the name of the invoking callable.

       @kwarg up: Number of call stack frames up (C{int}).
       @kwarg dflt: Default return value (C{any}).
       @kwarg source: Include source file name and line
                      number (C{bool}).
       @kwarg underOK: Private, internal callables are OK (C{bool}).

       @return: The callable name (C{str}) or B{C{dflt}} if none found.
    '''
    try:  # see .lazily._caller3
        for u in range(up, up + 32):
            n, f, s = _caller3(u)
            if n and (underOK or n.startswith(_DUNDER_) or
                             not n.startswith(_UNDER_)):
                if source:
                    n = NN(n, _AT_, f, _COLON_, str(s))
                return n
    except (AttributeError, ValueError):
        pass
    return dflt


def _callname(name, class_name, self_name, up=1):
    '''(INTERNAL) Assemble the name for an invokation.
    '''
    n, c = class_name, callername(up=up + 1)
    if c:
        n = _DOT_(n, Fmt.PAREN(c, name))
    if self_name:
        n = _SPACE_(n, repr(self_name))
    return n


def classname(inst, prefixed=None):
    '''Return the instance' class name optionally prefixed with the
       module name.

       @arg inst: The object (any C{type}).
       @kwarg prefixed: Include the module name (C{bool}), see
                        function C{classnaming}.

       @return: The B{C{inst}}'s C{[module.]class} name (C{str}).
    '''
    if prefixed is None:
        prefixed = getattr(inst, classnaming.__name__, prefixed)
    return modulename(inst.__class__, prefixed=prefixed)


def classnaming(prefixed=None):
    '''Get/set the default class naming for C{[module.]class} names.

       @kwarg prefixed: Include the module name (C{bool}).

       @return: Previous class naming setting (C{bool}).
    '''
    t = _Named._classnaming
    if prefixed in (True, False):
        _Named._classnaming = prefixed
    return t


def modulename(clas, prefixed=None):  # in .basics._xversion
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
        n = '__name__'  # _DUNDER_(NN, _name_, NN)
    if prefixed or (classnaming() if prefixed is None else False):
        try:
            m =  clas.__module__.rsplit(_DOT_, 1)
            n = _DOT_.join(m[1:] + [n])
        except AttributeError:
            pass
    return n


def nameof(inst):
    '''Get the name of an instance.

       @arg inst: The object (any C{type}).

       @return: The instance' name (C{str}) or C{""}.
    '''
    n = getattr(inst, _name_, NN)
    if not n:  # and isinstance(inst, property):
        try:
            n = inst.fget.__name__
        except AttributeError:
            n = NN
    return n


def _notError(inst, name, args, kwds):  # PYCHOK no cover
    '''(INTERNAL) Format an error message.
    '''
    n = _DOT_(classname(inst, prefixed=True), _dunder_nameof(name, name))
    m = _COMMASPACE_.join(modulename(c, prefixed=True) for c in inst.__class__.__mro__[1:-1])
    return _COMMASPACE_(unstr(n, *args, **kwds), Fmt.PAREN(_MRO_, m))


def _NotImplemented(inst, *other, **kwds):
    '''(INTERNAL) Raise a C{__special__} error or return C{NotImplemented},
       but only if env variable C{PYGEODESY_NOTIMPLEMENTED=std}.
    '''
    if _std_NotImplemented:
        return NotImplemented
    n = _DOT_(classname(inst), callername(up=2, underOK=True))  # source=True
    raise _NotImplementedError(unstr(n, *other, **kwds), txt=repr(inst))


def notImplemented(inst, *args, **kwds):  # PYCHOK no cover
    '''Raise a C{NotImplementedError} for a missing method or property.

       @arg inst: Instance (C{any}).
       @arg args: Method or property positional arguments (any C{type}s).
       @arg kwds: Method or property keyword arguments (any C{type}s).
    '''
    n =  kwds.pop(callername.__name__, NN) or callername(up=2)
    t = _notError(inst, n, args, kwds)
    raise _NotImplementedError(t, txt=notImplemented.__name__.replace('I', ' i'))


def notOverloaded(inst, *args, **kwds):  # PYCHOK no cover
    '''Raise an C{AssertionError} for a method or property not overloaded.

       @arg inst: Instance (C{any}).
       @arg args: Method or property positional arguments (any C{type}s).
       @arg kwds: Method or property keyword arguments (any C{type}s).
    '''
    n =  kwds.pop(callername.__name__, NN) or callername(up=2)
    t = _notError(inst, n, args, kwds)
    raise _AssertionError(t, txt=notOverloaded.__name__.replace('O', ' o'))


def _Pass(arg, **unused):  # PYCHOK no cover
    '''(INTERNAL) I{Pass-thru} class for C{_NamedTuple._Units_}.
    '''
    return arg


__all__ += _ALL_DOCS(_Named,
                     _NamedBase,  # _NamedDict,
                     _NamedEnum, _NamedEnumItem,
                     _NamedTuple)

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

# % env PYGEODESY_FOR_DOCS=1 python -m pygeodesy.named
# all 71 locals OK
