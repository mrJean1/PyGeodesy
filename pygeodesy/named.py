
# -*- coding: utf-8 -*-

u'''(INTERNAL) Nameable class instances.

Classes C{_Named}, C{_NamedDict}, C{_NamedEnum}, C{_NamedEnumItem} and
C{_NamedTuple} and several subclasses thereof, all with nameable instances.

The items in a C{_NamedDict} are accessable as attributes and the items
in a C{_NamedTuple} are named to be accessable as attributes, similar to
standard Python C{namedtuple}s.

@see: Module L{pygeodesy.namedTuples} for (most of) the C{Named-Tuples}.
'''

from pygeodesy.basics import isidentifier, iskeyword, isstr, itemsorted, len2, \
                            _xcopy, _xdup, _xinstanceof, _xsubclassof, _zip
from pygeodesy.errors import _AssertionError, _AttributeError, _incompatible, \
                             _IndexError, _KeyError, LenError, _NameError, \
                             _NotImplementedError, _TypeError, _TypesError, \
                             _UnexpectedError, UnitError, _ValueError, \
                             _xattr, _xkwds, _xkwds_item2, _xkwds_pop2
from pygeodesy.internals import _caller3, _dunder_nameof, _isPyPy, _sizeof, _under
from pygeodesy.interns import MISSING, NN, _AT_, _COLON_, _COLONSPACE_, _COMMA_, \
                             _COMMASPACE_, _doesn_t_exist_, _DOT_, _DUNDER_, \
                             _dunder_name_, _EQUAL_, _exists_, _immutable_, _name_, \
                             _NL_, _NN_, _no_, _other_, _s_, _SPACE_, _std_, \
                             _UNDER_, _vs_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS, _getenv
from pygeodesy.props import _allPropertiesOf_n, deprecated_method, _hasProperty, \
                            _update_all, property_doc_, Property_RO, property_RO, \
                            _update_attrs
from pygeodesy.streprs import attrs, Fmt, lrstrip, pairs, reprs, unstr
# from pygeodesy.units import _toUnit  # _MODS

__all__ = _ALL_LAZY.named
__version__ = '24.06.24'

_COMMANL_           = _COMMA_ + _NL_
_COMMASPACEDOT_     = _COMMASPACE_ + _DOT_
_del_               = 'del'
_item_              = 'item'
_MRO_               = 'MRO'
# __DUNDER gets mangled in class
_name               = _under(_name_)
_Names_             = '_Names_'
_registered_        = 'registered'  # PYCHOK used!
_std_NotImplemented = _getenv('PYGEODESY_NOTIMPLEMENTED', NN).lower() == _std_
_such_              = 'such'
_Units_             = '_Units_'
_UP                 =  2


class ADict(dict):
    '''A C{dict} with both key I{and} attribute access to
       the C{dict} items.
    '''
    _iteration = None  # Iteration number (C{int}) or C{None}

    def __getattr__(self, name):
        '''Get the value of an item by B{C{name}}.
        '''
        try:
            return self[name]
        except KeyError:
            if name == _name_:
                return NN
        raise self._AttributeError(name)

    def __repr__(self):
        '''Default C{repr(self)}.
        '''
        return self.toRepr()

    def __setattr__(self, name, value):
        '''Set the value of a I{known} item by B{C{name}}.
        '''
        try:
            if self[name] != value:
                self[name] = value
        except KeyError:
            dict.__setattr__(self, name, value)

    def __str__(self):
        '''Default C{str(self)}.
        '''
        return self.toStr()

    def _AttributeError(self, name):
        '''(INTERNAL) Create an C{AttributeError}.
        '''
        if _DOT_ not in name:  # NOT classname(self)!
            name = _DOT_(self.__class__.__name__, name)
        return _AttributeError(item=name, txt=_doesn_t_exist_)

    @property_RO
    def iteration(self):  # see .named._NamedBase
        '''Get the iteration number (C{int}) or
           C{None} if not available/applicable.
        '''
        return self._iteration

    def set_(self, iteration=None, **items):  # PYCHOK signature
        '''Add one or several new items or replace existing ones.

           @kwarg iteration: Optional C{iteration} (C{int}).
           @kwarg items: One or more C{name=value} pairs.
        '''
        if iteration is not None:
            self._iteration = iteration
        if items:
            dict.update(self, items)
        return self  # in RhumbLineBase.Intersecant2, _PseudoRhumbLine.Position

    def toRepr(self, **prec_fmt):
        '''Like C{repr(dict)} but with C{name} prefix and with
           C{floats} formatted by function L{pygeodesy.fstr}.
        '''
        n = _xattr(self, name=NN) or self.__class__.__name__
        return Fmt.PAREN(n, self._toT(_EQUAL_, **prec_fmt))

    def toStr(self, **prec_fmt):
        '''Like C{str(dict)} but with C{floats} formatted by
           function L{pygeodesy.fstr}.
        '''
        return Fmt.CURLY(self._toT(_COLONSPACE_, **prec_fmt))

    def _toT(self, sep, **kwds):
        '''(INTERNAL) Helper for C{.toRepr} and C{.toStr}.
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
        return Fmt.repr_at(self)

    def __rmatmul__(self, other):  # PYCHOK no cover
        '''Not implemented.'''
        return _NotImplemented(self, other)  # PYCHOK Python 3.5+

    def __str__(self):
        '''Default C{str(self)}.
        '''
        return self.named2

    def attrs(self, *names, **sep_Nones_pairs_kwds):
        '''Join named attributes as I{name=value} strings, with C{float}s formatted by
           function L{pygeodesy.fstr}.

           @arg names: The attribute names, all positional (C{str}).
           @kwarg sep_Nones_pairs_kwds: Keyword arguments for function L{pygeodesy.pairs},
                      except C{B{sep}=", "} and C{B{Nones}=True} to in-/exclude missing
                      or C{None}-valued attributes.

           @return: All C{name=value} pairs, joined by B{C{sep}} (C{str}).

           @see: Functions L{pygeodesy.attrs}, L{pygeodesy.pairs} and L{pygeodesy.fstr}
        '''
        sep, kwds = _xkwds_pop2(sep_Nones_pairs_kwds, sep=_COMMASPACE_)
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
        '''Set the class naming for C{[module.].class} names (C{bool})
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

    def copy(self, deep=False, **name):
        '''Make a shallow or deep copy of this instance.

           @kwarg deep: If C{True} make a deep, otherwise
                        a shallow copy (C{bool}).
           @kwarg name: Optional, non-empty C{B{name}=NN} (C{str}).

           @return: The copy (C{This class}).
        '''
        c = _xcopy(self, deep=deep)
        if name:
            _ = c.rename(name)
        return c

    def _DOT_(self, *names):
        '''(INTERNAL) Period-join C{self.name} and C{names}.
        '''
        return _DOT_(self.name, *names)

    def dup(self, deep=False, **items):
        '''Duplicate this instance, replacing some attributes.

           @kwarg deep: If C{True} duplicate deep, otherwise shallow.
           @kwarg items: Attributes to be changed (C{any}), including
                         optional C{B{name}} (C{str}).

           @return: The duplicate (C{This class}).

           @raise AttributeError: Some B{C{items}} invalid.
        '''
        n = self.name
        m, items = _xkwds_pop2(items, name=n)
        d = _xdup(self, deep=deep, **items)
        if m != n:
            d.rename(m)  # zap _Named_Property_ROs
#       if items:
#           _update_all(d)
        return d

    def _instr(self, *attrs, **fmt_prec_props_sep_name__kwds):
        '''(INTERNAL) Format, used by C{Conic}, C{Ellipsoid}, C{Geodesic...}, C{Transform}, C{Triaxial}.
        '''
        def _fmt_prec_props_kwds(fmt=Fmt.F, prec=6, props=(), sep=_COMMASPACE_, **kwds):
            return fmt, prec, props, sep, kwds

        name, kwds = _name2__(**fmt_prec_props_sep_name__kwds)
        fmt, prec, props, sep, kwds = _fmt_prec_props_kwds(**kwds)

        t = () if name is None else (Fmt.EQUAL(name=repr(name or self.name)),)
        if attrs:
            t += pairs(((a, getattr(self, a)) for a in attrs),
                       prec=prec, ints=True, fmt=fmt)
        if props:
            t += pairs(((p.name, getattr(self, p.name)) for p in props),
                       prec=prec, ints=True)
        if kwds:
            t += pairs(kwds, prec=prec)
        return sep.join(t) if sep else t

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
        '''Set the C{B{name}=NN} (C{str}).

           @raise NameError: Can't rename, use method L{rename}.
        '''
        m, n = self._name, _name__(name)
        if not m:
            self._name = n
        elif n != m:
            n =  repr(n)
            c =  self.classname
            t = _DOT_(c, Fmt.PAREN(self.rename.__name__, n))
            n = _DOT_(c, Fmt.EQUALSPACED(name=n))
            m =  Fmt.PAREN(_SPACE_('was', repr(m)))
            n = _SPACE_(n, m)
            raise _NameError(n, txt=_SPACE_('use', t))
        # to set the name from a sub-class, use
        #   self.name = name or
        #  _Named.name.fset(self, name), but NOT
        #  _Named(self).name = name

    def _name__(self, name):  # usually **name
        '''(INTERNAL) Get the C{name} or this C{name}.
        '''
        return _name__(name, _or_nameof=self)  # nameof(self)

    def _name1__(self, kwds):
        '''(INTERNAL) Resolve and set the C{B{name}=NN}.
        '''
        return _name1__(kwds, _or_nameof=self.name) if self.name else kwds

    @Property_RO
    def named(self):
        '''Get the name I{or} class name or C{""} (C{str}).
        '''
        return self.name or self.classname

#   @Property_RO
#   def named_(self):
#       '''Get the C{class} name I{and/or} the str(name) or C{""} (C{str}).
#       '''
#       return _xjoined_(self.classname, self.name, enquote=False)

    @Property_RO
    def named2(self):
        '''Get the C{class} name I{and/or} the repr(name) or C{""} (C{str}).
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

    def _notImplemented(self, *args, **kwds):
        '''(INTERNAL) See function L{notImplemented}.
        '''
        notImplemented(self, *args, **_xkwds(kwds, up=_UP + 1))

    def _notOverloaded(self, *args, **kwds):
        '''(INTERNAL) See function L{notOverloaded}.
        '''
        notOverloaded(self, *args, **_xkwds(kwds, up=_UP + 1))

    def rename(self, name):
        '''Change the name.

           @arg name: The new name (C{str}).

           @return: The previous name (C{str}).
        '''
        m, n = self._name, _name__(name)
        if n != m:
            self._name = n
            _update_attrs(self, *_Named_Property_ROs)
        return m

    def renamed(self, name):
        '''Change the name.

           @arg name: The new name (C{str}).

           @return: This instance (C{str}).
        '''
        _ = self.rename(name)
        return self

    @property_RO
    def sizeof(self):
        '''Get the current size in C{bytes} of this instance (C{int}).
        '''
        return _sizeof(self)

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
        '''DEPRECATED on 23.10.07, use method C{toRepr}.'''
        return self.toRepr(**kwds)

#   def _unstr(self, which, *args, **kwds):
#       '''(INTERNAL) Return the string representation of a method
#          invokation of this instance: C{str(self).method(...)}
#
#          @see: Function L{pygeodesy.unstr}.
#       '''
#       return _DOT_(self, unstr(which, *args, **kwds))

    def _xnamed(self, inst, name=NN, **force):
        '''(INTERNAL) Set the instance' C{.name = self.name}.

           @arg inst: The instance (C{_Named}).
           @kwarg name: The name (C{str}).
           @kwarg force: If C{True}, force rename (C{bool}).

           @return: The B{C{inst}}, renamed if B{C{force}}d
                    or if not named before.
        '''
        return _xnamed(inst, name, **force)

    def _xrenamed(self, inst):
        '''(INTERNAL) Rename the instance' C{.name = self.name}.

           @arg inst: The instance (C{_Named}).

           @return: The B{C{inst}}, named if not named before.

           @raise TypeError: Not C{isinstance(B{inst}, _Named)}.
        '''
        _xinstanceof(_Named, inst=inst)  # assert
        return inst.renamed(self.name)

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

    def others(self, *other, **name_other_up):
        '''Refined class comparison, invoked as C{.others(other)},
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
        '''I{Must be overloaded}.'''
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


class _NamedDict(ADict, _Named):
    '''(INTERNAL) Named C{dict} with key I{and} attribute access
       to the items.
    '''
    def __init__(self, *args, **kwds):
        if args:  # is args[0] a dict?
            if len(args) != 1:  # or not isinstance(args[0], dict)
                kwds = _name1__(kwds)
                t = unstr(self.classname, *args, **kwds)  # PYCHOK no cover
                raise _ValueError(args=len(args), txt=t)
            kwds = _xkwds(dict(args[0]), **kwds)  # args[0] overrides kwds
        n, kwds = _name2__(**kwds)
        if n:
            _Named.name.fset(self, n)  # see _Named.name
        ADict.__init__(self, kwds)

    def __delattr__(self, name):
        '''Delete an attribute or item by B{C{name}}.
        '''
        if name in self:  # in ADict.keys(self):
            ADict.pop(self, name)
        elif name in (_name_, _name):
            # ADict.__setattr__(self, name, NN)
            _Named.rename(self, NN)
        else:
            ADict.__delattr__(self, name)

    def __getattr__(self, name):
        '''Get the value of an item by B{C{name}}.
        '''
        try:
            return self[name]
        except KeyError:
            if name == _name_:
                return _Named.name.fget(self)
        raise ADict._AttributeError(self, self._DOT_(name))

    def __getitem__(self, key):
        '''Get the value of an item by B{C{key}}.
        '''
        if key == _name_:
            raise self._KeyError(key)
        return ADict.__getitem__(self, key)

    def _KeyError(self, key, *value):  # PYCHOK no cover
        '''(INTERNAL) Create a C{KeyError}.
        '''
        n = self.name or self.__class__.__name__
        t = Fmt.SQUARE(n, key)
        if value:
            t = Fmt.EQUALSPACED(t, *value)
        return _KeyError(t)

    def __setattr__(self, name, value):
        '''Set attribute or item B{C{name}} to B{C{value}}.
        '''
        if name in self:  # in ADict.keys(self)
            ADict.__setitem__(self, name, value)  # self[name] = value
        else:
            ADict.__setattr__(self, name, value)

    def __setitem__(self, key, value):
        '''Set item B{C{key}} to B{C{value}}.
        '''
        if key == _name_:
            raise self._KeyError(key, repr(value))
        ADict.__setitem__(self, key, value)


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
        n = _name__(**name) or NN(Class.__name__, _s_)  # _dunder_nameof
        if n and _xvalid(n, underOK=True):
            _Named.name.fset(self, n)  # see _Named.name

    def __getattr__(self, name):
        '''Get the value of an attribute or item by B{C{name}}.
        '''
        return _NamedDict.__getattr__(self, name)

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
        pypy = _isPyPy()
        _isa =  isinstance
        for n, v in kwds.items():
            if _isa(v, _LazyNamedEnumItem):  # property
                assert (n == v.name) if pypy else (n is v.name)
                # assert not hasattr(self.__class__, n)
                setattr(self.__class__, n, v)
            elif _isa(v, self._item_Classes):  # PYCHOK no cover
                assert self[n] is v and getattr(self, n) \
                                    and self.find(v) == n
            else:
                raise _TypeError(v, name=n)

    def find(self, item, dflt=None, all=False):
        '''Find a registered item.

           @arg item: The item to look for (any C{type}).
           @kwarg dflt: Value to return if not found (any C{type}).
           @kwarg all: Use C{True} to search I{all} items or C{False} only
                       the currently I{registered} ones (C{bool}).

           @return: The B{C{item}}'s name if found (C{str}), or C{{dflt}}
                    if there is no such B{C{item}}.
        '''
        for k, v in self.items(all=all):  # or ADict.items(self)
            if v is item:
                return k
        return dflt

    def get(self, name, dflt=None):
        '''Get the value of a I{registered} item.

           @arg name: The name of the item (C{str}).
           @kwarg dflt: Value to return (any C{type}).

           @return: The item with B{C{name}} if found, or B{C{dflt}} if
                    there is no I{registered} item with that B{C{name}}.
        '''
        # getattr needed to instantiate L{_LazyNamedEnumItem}
        return getattr(self, name, dflt)

    def items(self, all=False, asorted=False):
        '''Yield all or only the I{registered} items.

           @kwarg all: Use C{True} to yield I{all} items or C{False} for
                       only the currently I{registered} ones (C{bool}).
           @kwarg asorted: If C{True}, yield the items in I{alphabetical,
                           case-insensitive} order (C{bool}).
        '''
        if all:  # instantiate any remaining L{_LazyNamedEnumItem}
            _isa = isinstance
            for n, p in tuple(self.__class__.__dict__.items()):
                if _isa(p, _LazyNamedEnumItem):
                    _ = getattr(self, n)
        return itemsorted(self) if asorted else ADict.items(self)

    def keys(self, **all_asorted):
        '''Yield the name (C{str}) of I{all} or only the currently I{registered}
           items, optionally sorted I{alphabetically, case-insensitively}.

           @kwarg all_asorted: See method C{items}.
        '''
        for k, _ in self.items(**all_asorted):
            yield k

    def popitem(self):
        '''Remove I{an, any} currently I{registed} item.

           @return: The removed item.
        '''
        return self._zapitem(*ADict.popitem(self))

    def register(self, item):
        '''Registed one new item or I{all} or I{any} unregistered ones.

           @arg item: The item (any C{type}) or B{I{all}} or B{C{any}}.

           @return: The item name (C{str}) or C("all") or C{"any"}.

           @raise NameError: An B{C{item}} with that name is already
                             registered the B{C{item}} has no or an
                             invalid name.

           @raise TypeError: The B{C{item}} type invalid.
        '''
        if item is all or item is any:
            _ = self.items(all=True)
            n = item.__name__
        else:
            try:
                n = item.name
                if not (n and isstr(n) and isidentifier(n)):
                    raise ValueError()
            except (AttributeError, ValueError, TypeError) as x:
                n = _DOT_(_item_, _name_)
                raise _NameError(n, item, cause=x)
            if n in self:
                t = _SPACE_(_item_, self._DOT_(n), _exists_)
                raise _NameError(t, txt=repr(item))
            if not isinstance(item, self._item_Classes):  # _xinstanceof
                n = self._DOT_(n)
                raise _TypesError(n, item, *self._item_Classes)
            self[n] = item
        return n

    def unregister(self, name_or_item):
        '''Remove a I{registered} item.

           @arg name_or_item: Name (C{str}) or the item (any C{type}).

           @return: The unregistered item.

           @raise AttributeError: No such B{C{item}}.

           @raise NameError: No item with that B{C{name}}.
        '''
        if isstr(name_or_item):
            name = name_or_item
        else:
            name = self.find(name_or_item, dflt=MISSING)  # all=True?
            if name is MISSING:
                t = _SPACE_(_no_, _such_, self.classname, _item_)
                raise _AttributeError(t, txt=repr(name_or_item))
        try:
            item = ADict.pop(self, name)
        except KeyError:
            raise _NameError(item=self._DOT_(name), txt=_doesn_t_exist_)
        return self._zapitem(name, item)

    pop = unregister

    def toRepr(self, prec=6, fmt=Fmt.F, sep=_COMMANL_, **all_asorted):  # PYCHOK _NamedDict, ADict
        '''Like C{repr(dict)} but C{name}s optionally sorted and
           C{floats} formatted by function L{pygeodesy.fstr}.
        '''
        t = ((self._DOT_(n), v) for n, v in self.items(**all_asorted))
        return sep.join(pairs(t, prec=prec, fmt=fmt, sep=_COLONSPACE_))

    def toStr(self, *unused, **all_asorted):  # PYCHOK _NamedDict, ADict
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

       @see: Luciano Ramalho, "Fluent Python", O'Reilly, Example
             19-24, 2016 p. 636 or Example 22-28, 2022 p. 869+
    '''
    def _fget(inst):
        # assert isinstance(inst, _NamedEnum)
        try:  # get the item from the instance' __dict__
            # item = inst.__dict__[name]  # ... or ADict
            item = inst[name]
        except KeyError:
            # instantiate an _NamedEnumItem, it self-registers
            item = inst._Lazy(*args, **_xkwds(kwds, name=name))
            # assert inst[name] is item  # MUST be registered
            # store the item in the instance' __dict__ ...
            # inst.__dict__[name] = item  # ... or update the
            inst.update({name: item})  # ... ADict for Triaxials
            # remove the property from the registry class, such that
            # (a) the property no longer overrides the instance' item
            # in inst.__dict__ and (b) _NamedEnum.items(all=True) only
            # sees any un-instantiated ones yet to be instantiated
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
        name = _name__(name) or _NN_
        if self._enum:
            raise _NameError(name, self, txt=_registered_)  # _TypeError
        if name:
            self._name = name

    def _register(self, enum, name):
        '''(INTERNAL) Add this item as B{C{enum.name}}.

           @note: Don't register if name is empty or doesn't
                  start with a letter.
        '''
        name = _name__(name)
        if name and _xvalid(name, underOK=True):
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
            if item is not self:  # PYCHOK no cover
                t = _SPACE_(repr(item), _vs_, repr(self))
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
    _Units_ = ()    # .units classes
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
            i, name = _xkwds_pop2(iteration_name, iteration=None)
            if i is not None:
                self._iteration = i
            if name:
                self.name = name
        return self

    def __delattr__(self, name):
        '''Delete an attribute by B{C{name}}.

           @note: Items can not be deleted.
        '''
        if name in self._Names_:
            t = _SPACE_(_del_, self._DOT_(name))
            raise _TypeError(t, txt=_immutable_)
        elif name in (_name_, _name):
            _Named.__setattr__(self, name, NN)  # XXX _Named.name.fset(self, NN)
        else:
            tuple.__delattr__(self, name)

    def __getattr__(self, name):
        '''Get the value of an attribute or item by B{C{name}}.
        '''
        try:
            return tuple.__getitem__(self, self._Names_.index(name))
        except IndexError as x:
            raise _IndexError(self._DOT_(name), cause=x)
        except ValueError:  # e.g. _iteration
            return tuple.__getattr__(self, name)  # __getattribute__

#   def __getitem__(self, index):  # index, slice, etc.
#       '''Get the item(s) at an B{C{index}} or slice.
#       '''
#       return tuple.__getitem__(self, index)

    def __hash__(self):
        return tuple.__hash__(self)

    def __repr__(self):
        '''Default C{repr(self)}.
        '''
        return self.toRepr()

    def __setattr__(self, name, value):
        '''Set attribute or item B{C{name}} to B{C{value}}.
        '''
        if name in self._Names_:
            t = Fmt.EQUALSPACED(self._DOT_(name), repr(value))
            raise _TypeError(t, txt=_immutable_)
        elif name in (_name_, _name):
            _Named.__setattr__(self, name, value)  # XXX _Named.name.fset(self, value)
        else:  # e.g. _iteration
            tuple.__setattr__(self, name, value)

    def __str__(self):
        '''Default C{repr(self)}.
        '''
        return self.toStr()

    def _DOT_(self, *names):
        '''(INTERNAL) Period-join C{self.classname} and C{names}.
        '''
        return _DOT_(self.classname, *names)

    def dup(self, name=NN, **items):
        '''Duplicate this tuple replacing one or more items.

           @kwarg name: Optional new name (C{str}).
           @kwarg items: Items to be replaced (C{name=value} pairs), if any.

           @return: A copy of this tuple with B{C{items}}.

           @raise NameError: Some B{C{items}} invalid.
        '''
        t = list(self)
        U = self._Units_
        if items:
            _ix =  self._Names_.index
            _2U = _MODS.units._toUnit
            try:
                for n, v in items.items():
                    i = _ix(n)
                    t[i] = _2U(U[i], v, name=n)
            except ValueError:  # bad item name
                raise _NameError(self._DOT_(n), v, this=self)
        return self.classof(*t).reUnit(*U, name=name)

    def items(self):
        '''Yield the items, each as a C{(name, value)} pair (C{2-tuple}).

           @see: Method C{.units}.
        '''
        for n, v in _zip(self._Names_, self):  # strict=True
            yield n, v

    iteritems = items

    def reUnit(self, *Units, **name):
        '''Replace some of this C{Named-Tuple}'s C{Units}.

           @arg Units: One or more C{Unit} classes, all positional.
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: This instance with updated C{Units}.

           @note: This C{Named-Tuple}'s values are I{not updated}.
        '''
        U = self._Units_
        n = min(len(U), len(Units))
        if n:
            R = Units + U[n:]
            if R != U:
                self._Units_ = R
        return self.renamed(name) if name else self

    def toRepr(self, prec=6, sep=_COMMASPACE_, fmt=Fmt.F, **unused):  # PYCHOK signature
        '''Return this C{Named-Tuple} items as C{name=value} string(s).

           @kwarg prec: The C{float} precision, number of decimal digits (0..9).
                        Trailing zero decimals are stripped for B{C{prec}} values
                        of 1 and above, but kept for negative B{C{prec}} values.
           @kwarg sep: Separator to join (C{str}).
           @kwarg fmt: Optional C{float} format (C{letter}).

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
           @kwarg fmt: Optional C{float} format (C{letter}).

           @return: Tuple items (C{str}).
        '''
        return Fmt.PAREN(sep.join(reprs(self, prec=prec, fmt=fmt)))

    def toUnits(self, Error=UnitError, **name):  # overloaded in .frechet, .hausdorff
        '''Return a copy of this C{Named-Tuple} with each item value wrapped
           as an instance of its L{units} class.

           @kwarg Error: Error to raise for L{units} issues (C{UnitError}).
           @kwarg name: Optional C{B{name}=NN} (C{str}).

           @return: A duplicate of this C{Named-Tuple} (C{C{Named-Tuple}}).

           @raise Error: Invalid C{Named-Tuple} item or L{units} class.
        '''
        t = tuple(v for _, v in self.units(Error=Error))
        return self.classof(*t).reUnit(*self._Units_, **name)

    def units(self, **Error):
        '''Yield the items, each as a C{2-tuple (name, value}) with the
           value wrapped as an instance of its L{units} class.

           @kwarg Error: Optional C{B{Error}=UnitError} to raise.

           @raise Error: Invalid C{Named-Tuple} item or L{units} class.

           @see: Method C{.items}.
        '''
        _2U = _MODS.units._toUnit
        for n, v, U in _zip(self._Names_, self, self._Units_):  # strict=True
            yield n, _2U(U, v, name=n, **Error)

    iterunits = units

    def _validate(self, underOK=False):  # see .EcefMatrix
        '''(INTERNAL) One-time check of C{_Names_} and C{_Units_}
           for each C{_NamedUnit} I{sub-class separately}.
        '''
        ns = self._Names_
        if not (isinstance(ns, tuple) and len(ns) > 1):  # XXX > 0
            raise _TypeError(self._DOT_(_Names_), ns)
        for i, n in enumerate(ns):
            if not _xvalid(n, underOK=underOK):
                t = Fmt.SQUARE(_Names_=i)  # PYCHOK no cover
                raise _ValueError(self._DOT_(t), n)

        us = self._Units_
        if not isinstance(us, tuple):
            raise _TypeError(self._DOT_(_Units_), us)
        if len(us) != len(ns):
            raise LenError(self.__class__, _Units_=len(us), _Names_=len(ns))
        for i, u in enumerate(us):
            if not (u is None or callable(u)):
                t = Fmt.SQUARE(_Units_=i)  # PYCHOK no cover
                raise _TypeError(self._DOT_(t), u)

        self.__class__._validated = True

    def _xtend(self, xTuple, *items, **name):
        '''(INTERNAL) Extend this C{Named-Tuple} with C{items} to an other B{C{xTuple}}.
        '''
        _xsubclassof(_NamedTuple, xTuple=xTuple)
        if len(xTuple._Names_)       != (len(self._Names_) + len(items)) or \
               xTuple._Names_[:len(self)] != self._Names_:  # PYCHOK no cover
            c = NN(self.classname,  repr(self._Names_))
            x = NN(xTuple.__name__, repr(xTuple._Names_))
            raise TypeError(_SPACE_(c, _vs_, x))
        t = self + items
        return xTuple(t, name=self._name__(name))  # .reUnit(*self._Units_)


def callername(up=1, dflt=NN, source=False, underOK=False):
    '''Get the name of the invoking callable.

       @kwarg up: Number of call stack frames up (C{int}).
       @kwarg dflt: Default return value (C{any}).
       @kwarg source: Include source file name and line number (C{bool}).
       @kwarg underOK: If C{True}, private, internal callables are OK,
                       otherwise private callables are skipped (C{bool}).

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


def _callername2(args, callername=NN, source=False, underOK=False, up=_UP, **kwds):
    '''(INTERNAL) Extract C{callername}, C{source}, C{underOK} and C{up} from C{kwds}.
    '''
    n = callername or _MODS.named.callername(up=up + 1, source=source,
                                        underOK=underOK or bool(args or kwds))
    return n, kwds


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
        n =  clas.__name__
    except AttributeError:
        n = _dunder_name_
    if prefixed or (classnaming() if prefixed is None else False):
        try:
            m =  clas.__module__.rsplit(_DOT_, 1)
            n = _DOT_.join(m[1:] + [n])
        except AttributeError:
            pass
    return n


# def _name__(name=NN, name__=None, _or_nameof=None, **kwds):
#     '''(INTERNAL) Get single keyword argument C{B{name}=NN|None}.
#     '''
#     if kwds:  # "unexpected keyword arguments ..."
#         m = _MODS.errors
#         raise m._UnexpectedError(**kwds)
#     if name:  # is given
#         n = _name__(**name) if isinstance(name, dict) else str(name)
#     elif name__ is not None:
#         n = getattr(name__, _dunder_name_, NN)  # _xattr(name__, __name__=NN)
#     else:
#         n = name  # NN or None or {} or any False type
#     if _or_nameof is not None and not n:
#         n = getattr(_or_nameof, _name_, NN) # _xattr(_or_nameof, name=NN)
#     return n  # str or None or {}


def _name__(name=NN, **kwds):
    '''(INTERNAL) Get single keyword argument C{B{name}=NN|None}.
    '''
    if name or kwds:
        name, kwds = _name2__(name, **kwds)
        if kwds:  # "unexpected keyword arguments ..."
            raise _UnexpectedError(**kwds)
    return name if name or name is None else NN


def _name1__(kwds_name, **name__or_nameof):
    '''(INTERNAL) Resolve and set the C{B{name}=NN}.
    '''
    if kwds_name or name__or_nameof:
        n, kwds_name = _name2__(kwds_name, **name__or_nameof)
        kwds_name.update(name=n)
    return kwds_name


def _name2__(name=NN, name__=None, _or_nameof=None, **kwds):
    '''(INTERNAL) Get the C{B{name}=NN|None} and other C{kwds}.
    '''
    if name:  # is given
        if isinstance(name, dict):
            kwds.update(name)  # kwds = _xkwds(kwds, **name)?
            n, kwds = _name2__(**kwds)
        else:
            n = str(name)
    elif name__ is not None:
        n = _dunder_nameof(name__, NN)
    else:
        n = name if name is None else NN
    if _or_nameof is not None and not n:
        n = _xattr(_or_nameof, name=NN)  # nameof
    return n, kwds  # (str or None or {}), dict


def nameof(inst):
    '''Get the name of an instance.

       @arg inst: The object (any C{type}).

       @return: The instance' name (C{str}) or C{""}.
    '''
    n = _xattr(inst, name=NN)
    if not n:  # and isinstance(inst, property):
        try:
            n = inst.fget.__name__
        except AttributeError:
            n = NN
    return n


def _notDecap(where):
    '''De-Capitalize C{where.__name__}.
    '''
    n = where.__name__
    c = n[3].lower()  # len(_not_)
    return NN(n[:3], _SPACE_, c, n[4:])


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
    n, kwds = _callername2(other, **kwds)  # source=True
    t = unstr(_DOT_(classname(inst), n), *other, **kwds)
    raise _NotImplementedError(t, txt=repr(inst))


def notImplemented(inst, *args, **kwds):  # PYCHOK no cover
    '''Raise a C{NotImplementedError} for a missing instance method or
       property or for a missing caller feature.

       @arg inst: Caller instance (C{any}) or C{None} for function.
       @arg args: Method or property positional arguments (any C{type}s).
       @arg kwds: Method or property keyword arguments (any C{type}s),
                  except C{B{callername}=NN}, C{B{underOK}=False} and
                  C{B{up}=2}.
    '''
    n, kwds = _callername2(args, **kwds)
    t = _notError(inst, n, args, kwds) if inst else unstr(n, *args, **kwds)
    raise _NotImplementedError(t, txt=_notDecap(notImplemented))


def notOverloaded(inst, *args, **kwds):  # PYCHOK no cover
    '''Raise an C{AssertionError} for a method or property not overloaded.

       @arg inst: Instance (C{any}).
       @arg args: Method or property positional arguments (any C{type}s).
       @arg kwds: Method or property keyword arguments (any C{type}s),
                  except C{B{callername}=NN}, C{B{underOK}=False} and
                  C{B{up}=2}.
    '''
    n, kwds = _callername2(args, **kwds)
    t = _notError(inst, n, args, kwds)
    raise _AssertionError(t, txt=_notDecap(notOverloaded))


def _Pass(arg, **unused):  # PYCHOK no cover
    '''(INTERNAL) I{Pass-thru} class for C{_NamedTuple._Units_}.
    '''
    return arg


def _xjoined_(prefix, name=NN, enquote=True, **name__or_nameof):
    '''(INTERNAL) Join C{prefix} and non-empty C{name}.
    '''
    if name__or_nameof:
        name = _name__(name, **name__or_nameof)
    if name and prefix:
        if enquote:
            name = repr(name)
        t = _SPACE_(prefix, name)
    else:
        t =  prefix or name
    return t


def _xnamed(inst, name=NN, force=False, **name__or_nameof):
    '''(INTERNAL) Set the instance' C{.name = B{name}}.

       @arg inst: The instance (C{_Named}).
       @kwarg name: The name (C{str}).
       @kwarg force: If C{True}, force rename (C{bool}).

       @return: The B{C{inst}}, renamed if B{C{force}}d
                or if not named before.
    '''
    if name__or_nameof:
        name = _name__(name, **name__or_nameof)
    if name and isinstance(inst, _Named):
        if not inst.name:
            inst.name = name
        elif force:
            inst.rename(name)
    return inst


def _xother3(inst, other, name=_other_, up=1, **name_other):
    '''(INTERNAL) Get C{name} and C{up} for a named C{other}.
    '''
    if name_other:  # and other is None
        name, other = _xkwds_item2(name_other)
    elif other and len(other) == 1:
        name, other = _name__(name), other[0]
    else:
        raise _AssertionError(name, other, txt=classname(inst, prefixed=True))
    return other, name, up


def _xotherError(inst, other, name=_other_, up=1):
    '''(INTERNAL) Return a C{_TypeError} for an incompatible, named C{other}.
    '''
    n = _callname(name, classname(inst, prefixed=True), inst.name, up=up + 1)
    return _TypeError(name, other, txt=_incompatible(n))


def _xvalid(name, underOK=False):
    '''(INTERNAL) Check valid attribute name C{name}.
    '''
    return bool(name and isstr(name)
                     and name != _name_
                     and (underOK or not name.startswith(_UNDER_))
                     and (not iskeyword(name))
                     and isidentifier(name))


__all__ += _ALL_DOCS(_Named,
                     _NamedBase,  # _NamedDict,
                     _NamedEnum, _NamedEnumItem,
                     _NamedTuple)

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
