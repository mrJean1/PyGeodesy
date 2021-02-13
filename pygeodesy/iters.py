
# -*- coding: utf-8 -*-

u'''Iterator classes L{LatLon2PsxyIter} and L{PointsIter} to iterate
over iterables, lists, sets, tuples, etc. with optional loop-back to
the initial items, skipping of duplicate items and copying of the
iterated items.

@newfield example: Example, Examples
'''

from pygeodesy.basics import issubclassof, len2, map2
from pygeodesy.errors import _IndexError, LenError, PointsError
from pygeodesy.interns import NN, _few_, _points_, _too_, _0_, _1_0
from pygeodesy.lazily import _ALL_LAZY, _ALL_DOCS
from pygeodesy.named import _Named
from pygeodesy.namedTuples import Point3Tuple, Points2Tuple
from pygeodesy.props import property_RO
from pygeodesy.streprs import Fmt
from pygeodesy.units import Int, Radius
from pygeodesy.utily import degrees2m, isNumpy2, isTuple2, wrap90, wrap180

__all__ = _ALL_LAZY.iters
__version__ = '21.02.12'

_items_ = 'items'


class _BaseIter(_Named):
    '''(INTERNAL) Iterator over items with loop-back and de-duplication.

       @see: Luciano Ramalho, "Fluent Python", page 418+, O'Reilly, 2016.
    '''
    _closed =  True
    _copies = ()
    _dedup  =  False
    _Error  =  LenError
    _items  =  None
    _loop   = ()
    _name   = _items_
    _prev   =  object()

    def __init__(self, items, loop=0, dedup=False, Error=None, name=NN):
        '''New iterator over an iterable of B{C{items}}.

           @arg items: Iterable (any C{type}).
           @kwarg loop: Number of loop-back items, also initial enumerate
                        and iterate index (non-negative C{int}).
           @kwarg dedup: Skip duplicate items (C{bool}).
           @kwarg Error: Error to raise (L{LenError}).
           @kwarg name: Optional name (C{str}).

           @raise Error: Insufficient number of B{C{items}}.
        '''
        if dedup:
            self._dedup = True
        if issubclassof(Error, Exception):
            self._Error = Error
        if name:
            self.rename(name)

        if isinstance(items, (list, tuple)):  # range in Python 2
            self._items = items
        self._iter = iter(items)
        self._indx = -1
        if Int(loop) > 0:
            try:
                self._loop = tuple(self.next for _ in range(loop))
                if self.loop != loop:
                    raise RuntimeError  # force Error
            except (RuntimeError, StopIteration):
                raise self._Error(self.name, self.loop, txt=_too_(_few_))

    @property_RO
    def copies(self):
        '''Get the saved copies, if any (C{tuple} or C{list}) and only I{once}.
        '''
        cs = self._copies
        if cs:
            self._copies = ()
        return cs

    @property_RO
    def dedup(self):
        '''Get the de-duplication setting (C{bool}).
        '''
        return self._dedup

    def enumerate(self, closed=False, copies=False, dedup=False):
        '''Yield all items, each as a 2-tuple C{(index, item)}.

           @kwarg closed: Look back to the first B{C{point(s)}}.
           @kwarg copies: Make a copy of all B{C{items}} (C{bool}).
           @kwarg dedup: Set de-duplication in loop-back (C{bool}).
        '''
        for item in self.iterate(closed=closed, copies=copies, dedup=dedup):
            yield self._indx, item

    def __getitem__(self, index):
        '''Get the item(s) at the given B{C{index}} or C{slice}.

           @raise IndexError: Invalid B{C{index}}, beyond B{C{loop}}.
        '''
        t = self._items or self._copies or self._loop
        try:  # Luciano Ramalho, "Fluent Python", page 293+, O'Reilly, 2016.
            if isinstance(index, slice):
                return t[index.start:index.stop:index.step]
            else:
                return t[index]
        except IndexError as x:
            t = Fmt.SQUARE(self.name, index)
            raise _IndexError(str(x), txt=t)

    def __iter__(self):  # PYCHOK no cover
        '''Make this iterator C{iterable}.
        '''
        # Luciano Ramalho, "Fluent Python", page 421, O'Reilly, 2016.
        return self.iterate()  # XXX or self?

    def iterate(self, closed=False, copies=False, dedup=False):
        '''Yield all items, each as C{item}.

           @kwarg closed: Look back to the first B{C{point(s)}}.
           @kwarg copies: Make a copy of all B{C{items}} (C{bool}).
           @kwarg dedup: Set de-duplication in loop-back (C{bool}).

           @raise Error: Using C{B{closed}=True} without B{C{loop}}-back.
        '''
        if closed and not self.loop:
            raise self._Error(closed=closed, loop=self.loop)

        if copies:
            if self._items:
                self._copies = self._items
                self._items  = copy_ = None
            else:
                self._copies = list(self._loop)
                copy_ = self._copies.append
        else:  # del B{C{items}} reference
            self._items = copy_ = None

        self._closed = closed
        if self._iter:
            try:
                next_ = self.next_
                if copy_:
                    while True:
                        item = next_(dedup=dedup)
                        copy_(item)
                        yield item
                else:
                    while True:
                        yield next_(dedup=dedup)
            except StopIteration:
                self._iter = ()  # del self._iter, prevent re-iterate

    @property_RO
    def loop(self):
        '''Get the B{C{loop}} setting (C{int}), always C{0} in loop-back.
        '''
        return len(self._loop)

    @property_RO
    def next(self):
        '''Get the next item.
        '''
        return self._next_dedup() if self._dedup else self._next(False)

#   __next__  # NO __next__ AND __iter__ ... see Ramalho, page 426

    def next_(self, dedup=False):
        '''Return the next item.

           @kwarg dedup: Set de-duplication for loop-back (C{bool}).
        '''
        return self._next_dedup() if self._dedup else self._next(dedup)

    def _next(self, dedup):
        '''Return the next item, regardless.

           @arg dedup: Set de-duplication for loop-back (C{bool}).
        '''
        try:
            self._indx += 1
            self._prev  = item = next(self._iter)
            return item
        except StopIteration:
            pass
        if self._closed and self._loop:  # loop back
            self._iter  = iter(self._loop)
            self._loop  = ()
            self._indx  = 0
            self._dedup = bool(dedup or self._dedup)
        return next(self._iter)

    def _next_dedup(self):
        '''Return the next item, different from the previous one.
        '''
        prev = self._prev
        item = self._next(True)
        while item == prev:
            item = self._next(True)
        return item


class PointsIter(_BaseIter):
    '''Iterator for C{points} with optional loop-back and copies.
    '''
    _base  = None
    _Error = PointsError

    def __init__(self, points, loop=0, base=None):
        '''New L{PointsIter} iterator.

           @arg points: C{Iterable} or C{list}, C{sequence}, C{set}, C{tuple},
                        etc. (C{point}s).
           @kwarg loop: Number of loop-back points, also initial C{enumerate}
                        and C{iterate} index (non-negative C{int}).
           @kwarg base: Optional B{C{points}} instance for type checking (C{any}).

           @raise PointsError: Insufficient number of B{C{points}}.

           @raise TypeError: Some B{C{points}} are not B{C{base}}.
       '''
        _BaseIter.__init__(self, points, loop=loop, name=_points_)

        if base and not (isNumpy2(points) or isTuple2(points)):
            self._base = base

    def enumerate(self, closed=False, copies=False):  # PYCHOK signature
        '''Iterate and yield each point as a 2-tuple C{(index, point)}.

           @kwarg closed: Look back to the first B{C{point(s)}}, de-dup'ed (C{bool}).
           @kwarg copies: Save a copy of all B{C{points}} (C{bool}).

           @raise PointsError: Insufficient number of B{C{points}} or using
                               C{B{closed}=True} without B{C{loop}}-back.

           @raise TypeError: Some B{C{points}} are not B{C{base}}-compatible.
        '''
        for p in self.iterate(closed=closed, copies=copies):
            yield self._indx, p

    def iterate(self, closed=False, copies=False):  # PYCHOK signature
        '''Iterate through all B{C{points}} starting at index C{loop}.

           @kwarg closed: Look back to the first B{C{point(s)}}, de-dup'ed (C{bool}).
           @kwarg copies: Save a copy of all B{C{points}} (C{bool}).

           @raise PointsError: Insufficient number of B{C{points}} or using
                               C{B{closed}=True} without B{C{loop}}-back.

           @raise TypeError: Some B{C{points}} are not B{C{base}}-compatible.
        '''
        if self._base:
            base_ = self._base.others
            fmt_  = Fmt.SQUARE(points=0).replace
        else:
            base_ = fmt_ = None

        n = self.loop if self._iter else 0
        for p in _BaseIter.iterate(self, closed=closed, copies=copies, dedup=closed):
            if base_:
                base_(p, name=fmt_(_0_, str(self._indx)), up=2)
            yield p
            n += 1
        if n < (4 if closed else 2):
            raise self._Error(self.name, n, txt=_too_(_few_))


class LatLon2PsxyIter(PointsIter):
    '''Iterate and convert for C{points} with optional loop-back and copies.
    '''
    _deg2m  = None
    _radius = None  # keep degrees
    _wrap   = True

    def __init__(self, points, loop=0, base=None, wrap=True, radius=None):
        '''New L{LatLon2PsxyIter} iterator.

           @note: The C{LatLon} latitude is considered the I{pseudo-y} and
                  longitude the I{pseudo-x} coordinate, like L{LatLon2psxy}.

           @arg points: C{Iterable} or C{list}, C{sequence}, C{set}, C{tuple},
                        etc. (C{LatLon}[]).
           @kwarg loop: Number of loop-back points, also initial C{enumerate}
                        and C{iterate} index (non-negative C{int}).
           @kwarg base: Optional B{C{points}} instance for type checking (C{any}).
           @kwarg wrap: Wrap lat- and longitudes (C{bool}).
           @kwarg radius: Mean earth radius (C{meter}) for conversion from
                          C{degrees} to C{meter} (or C{radians} if C{B{radius}=1}).

           @raise PointsError: Insufficient number of B{C{points}}.

           @raise TypeError: Some B{C{points}} are not B{C{base}}-compatible.
        '''
        PointsIter.__init__(self, points, loop=loop, base=base)
        if not wrap:
            self._wrap = False
        if radius:
            self._radius = r = Radius(radius)
            self._deg2m  = degrees2m(_1_0, r)

    def __getitem__(self, index):
        '''Get the point(s) at the given B{C{index}} or C{slice}.

           @raise IndexError: Invalid B{C{index}}, beyond B{C{loop}}.
        '''
        ll = PointsIter.__getitem__(self, index)
        if isinstance(index, slice):
            return map2(self._point3Tuple, ll)
        else:
            return self._point3Tuple(ll)

    def enumerate(self, closed=False, copies=False):  # PYCHOK signature
        '''Iterate and yield each point as a 2-tuple C{(index, L{Point3Tuple})}.

           @kwarg closed: Look back to the first B{C{point(s)}}, de-dup'ed (C{bool}).
           @kwarg copies: Save a copy of all B{C{points}} (C{bool}).

           @raise PointsError: Insufficient number of B{C{points}} or using
                               C{B{closed}=True} without B{C{loop}}-back.

           @raise TypeError: Some B{C{points}} are not B{C{base}}-compatible.
        '''
        return PointsIter.enumerate(self, closed=closed, copies=copies)

    def iterate(self, closed=False, copies=False):  # PYCHOK signature
        '''Iterate the B{C{points}} starting at index B{C{loop}} and
           yield each as a L{Point3Tuple}C{(x, y, ll)}.

           @kwarg closed: Loop back to the first B{C{point(s)}}, de-dup'ed (C{bool}).
           @kwarg copies: Save a copy of all B{C{points}} (C{bool}).

           @raise PointsError: Insufficient number of B{C{points}} or using
                               C{B{closed}=True} without B{C{loop}}-back.

           @raise TypeError: Some B{C{points}} are not B{C{base}}-compatible.
        '''
        if self._deg2m not in (None, _1_0):
            _point3Tuple = self._point3Tuple
        elif self._wrap:
            def _point3Tuple(ll):
                return Point3Tuple(wrap180(ll.lon), wrap90(ll.lat), ll)
        else:
            def _point3Tuple(ll):  # PYCHOK redef
                return Point3Tuple(ll.lon, ll.lat, ll)

        for ll in PointsIter.iterate(self, closed=closed, copies=copies):
            yield _point3Tuple(ll)

    def _point3Tuple(self, ll):
        '''(INTERNAL) Create a L{Point3Tuple} for point B{C{ll}}.
        '''
        x, y = ll.lon, ll.lat  # note, x, y = lon, lat
        if self._wrap:
            x, y = wrap180(x), wrap90(y)
        d = self._deg2m
        if d:  # convert degrees
            x *= d
            y *= d
        return Point3Tuple(x, y, ll)


def _imdex2(closed, n):  # PYCHOK by .clipy
    '''(INTERNAL) Return first and second index of C{range(B{n})}.
    '''
    return (n-1, 0) if closed else (0, 1)


def points2(points, closed=True, base=None, Error=PointsError):
    '''Check a path or polygon represented by points.

       @arg points: The path or polygon points (C{LatLon}[])
       @kwarg closed: Optionally, consider the polygon closed,
                      ignoring any duplicate or closing final
                      B{C{points}} (C{bool}).
       @kwarg base: Optionally, check all B{C{points}} against
                    this base class, if C{None} don't check.
       @kwarg Error: Exception to raise (C{ValueError}).

       @return: A L{Points2Tuple}C{(number, points)} with the number
                of points and the points C{list} or C{tuple}.

       @raise PointsError: Insufficient number of B{C{points}}.

       @raise TypeError: Some B{C{points}} are not B{C{base}}
                         compatible.
    '''
    n, points = len2(points)

    if closed:
        # remove duplicate or closing final points
        while n > 1 and points[n-1] in (points[0], points[n-2]):
            n -= 1
        # XXX following line is unneeded if points
        # are always indexed as ... i in range(n)
        points = points[:n]  # XXX numpy.array slice is a view!

    if n < (3 if closed else 1):
        raise Error(points=n, txt=_too_(_few_))

    if base and not (isNumpy2(points) or isTuple2(points)):
        for i in range(n):
            base.others(points[i], name=Fmt.SQUARE(points=i))

    return Points2Tuple(n, points)


__all__ += _ALL_DOCS(_BaseIter)

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
