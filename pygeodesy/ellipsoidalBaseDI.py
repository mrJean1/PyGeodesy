
# -*- coding: utf-8 -*-

u'''(INTERNAL) Ellipsoidal direct/inverse geodesy base class
C{LatLonEllipsoidalBaseDI} and functions.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division

from pygeodesy.basics import issubclassof
from pygeodesy.ellipsoidalBase import LatLonEllipsoidalBase, Property_RO
from pygeodesy.errors import _AssertionError, IntersectionError, _IsnotError, \
                             _ValueError, _xellipsoidal, _xError, _xkwds_not
from pygeodesy.fmath import favg, fmean_, fsum_
from pygeodesy.formy import _radical2
from pygeodesy.interns import _ellipsoidal_  # PYCHOK used!
from pygeodesy.interns import EPS, PI, _datum_, _epoch_, _exceed_PI_radians_, \
                             _height_, _no_, _near_concentric_, _reframe_, \
                             _too_, _0_0
from pygeodesy.lazily import _ALL_DOCS
from pygeodesy.namedTuples import Bearing2Tuple, Destination2Tuple, \
                           Intersection3Tuple, _LL4Tuple
# from pygeodesy.props import Property_RO  # from .ellipsoidalBase
from pygeodesy.streprs import Fmt
from pygeodesy.units import Height, Radius_, Scalar, _1mm as _TOL_M
from pygeodesy.utily import m2degrees, unroll180, wrap90, wrap180, wrap360

__all__ = ()
__version__ = '21.06.28'

_TRIPS = 17  # _intersect3, _intersects2, _nearestOn interations, 6 is sufficient


class LatLonEllipsoidalBaseDI(LatLonEllipsoidalBase):
    '''(INTERNAL) Base class for C{ellipsoidal*.LatLon} classes
       with I{overloaded} C{Direct} and C{Inverse} methods.
    '''

    def bearingTo2(self, other, wrap=False):
        '''Compute the initial and final bearing (forward and reverse
           azimuth) from this to an other point, using this C{Inverse}
           method.  See methods L{initialBearingTo} and  L{finalBearingTo}
           for more details.

           @arg other: The other point (C{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: A L{Bearing2Tuple}C{(initial, final)}.

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.
        '''
        r = self._Inverse(other, wrap)
        return Bearing2Tuple(r.initial, r.final, name=self.name)

    def destination(self, distance, bearing, height=None):
        '''Compute the destination point after having travelled for
           the given distance from this point along a geodesic given
           by an initial bearing, using this C{Direct} method.  See
           method L{destination2} for more details.

           @arg distance: Distance (C{meter}).
           @arg bearing: Initial bearing in (compass C{degrees360}).
           @kwarg height: Optional height, overriding the default
                          height (C{meter}, same units as C{distance}).

           @return: The destination point (C{LatLon}).
        '''
        return self._Direct(distance, bearing, self.classof, height).destination

    def destination2(self, distance, bearing, height=None):
        '''Compute the destination point and the final bearing (reverse
           azimuth) after having travelled for the given distance from
           this point along a geodesic given by an initial bearing,
           using this C{Direct} method.

           The distance must be in the same units as this point's datum
           axes, conventionally C{meter}.  The distance is measured on
           the surface of the ellipsoid, ignoring this point's height.

           The initial and final bearing (forward and reverse azimuth)
           are in compass C{degrees360}.

           The destination point's height and datum are set to this
           point's height and datum, unless the former is overridden.

           @arg distance: Distance (C{meter}).
           @arg bearing: Initial bearing (compass C{degrees360}).
           @kwarg height: Optional height, overriding the default
                          height (C{meter}, same units as C{distance}).

           @return: A L{Destination2Tuple}C{(destination, final)}.
        '''
        r = self._Direct(distance, bearing, self.classof, height)
        return self._xnamed(r)

    def _Direct(self, distance, bearing, LL, height):  # overloaded by I{Vincenty}
        '''(INTERNAL) I{Karney}'s C{Direct} method.

           @return: A L{Destination2Tuple}C{(destination, final)} or
                    a L{Destination3Tuple}C{(lat, lon, final)} if
                    B{C{LL}} is C{None}.
        '''
        g = self.geodesic
        r = g.Direct3(self.lat, self.lon, bearing, distance)
        if LL:
            r = self._Direct2Tuple(LL, height, r)
        return r

    def _Direct2Tuple(self, LL, height, r):
        '''(INTERNAL) Helper for C{._Direct} result L{Destination2Tuple}.
        '''
        h = self.height if height is None else height
        d = LL(wrap90(r.lat), wrap180(r.lon), height=h, datum=self.datum, name=self.name,
                 **_xkwds_not(None, epoch=self.epoch, reframe=self.reframe))
        return Destination2Tuple(d, wrap360(r.final))

    def distanceTo(self, other, wrap=False, **unused):  # ignore radius=R_M
        '''Compute the distance between this and an other point
           along a geodesic, using this C{Inverse} method. See method
           L{distanceTo3} for more details.

           @arg other: The other point (C{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Distance (C{meter}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.
        '''
        return self._Inverse(other, wrap, azis=False).distance

    def distanceTo3(self, other, wrap=False):
        '''Compute the distance, the initial and final bearing along
           a geodesic between this and an other point, using this
           C{Inverse} method.

           The distance is in the same units as this point's datum axes,
           conventionally meter.  The distance is measured on the surface
           of the ellipsoid, ignoring this point's height.

           The initial and final bearing (forward and reverse azimuth)
           are in compass C{degrees360} from North.

           @arg other: Destination point (C{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: A L{Distance3Tuple}C{(distance, initial, final)}.

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.
        '''
        return self._xnamed(self._Inverse(other, wrap))

    def finalBearingOn(self, distance, bearing):
        '''Compute the final bearing (reverse azimuth) after having
           travelled for the given distance along a geodesic given by
           an initial bearing from this point, using this C{Direct}
           method.  See method L{destination2} for more details.

           @arg distance: Distance (C{meter}).
           @arg bearing: Initial bearing (compass C{degrees360}).

           @return: Final bearing (compass C{degrees360}).
        '''
        return self._Direct(distance, bearing, None, None).final

    def finalBearingTo(self, other, wrap=False):
        '''Compute the final bearing (reverse azimuth) after having
           travelled along a geodesic from this point to an other
           point, using this C{Inverse} method.  See method
           L{distanceTo3} for more details.

           @arg other: The other point (C{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Final bearing (compass C{degrees360}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.
        '''
        return self._Inverse(other, wrap).final

    @Property_RO
    def geodesic(self):  # overloaded by I{Karney}'s, N/A for I{Vincenty}
        '''N/A, invalid (C{None} I{always}).
        '''
        return None  # PYCHOK no cover

    def initialBearingTo(self, other, wrap=False):
        '''Compute the initial bearing (forward azimuth) to travel
           along a geodesic from this point to an other point,
           using this C{Inverse} method.  See method L{distanceTo3}
           for more details.

           @arg other: The other point (C{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Initial bearing (compass C{degrees360}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.
        '''
        return self._Inverse(other, wrap).initial

    def intermediateTo(self, other, fraction, height=None, wrap=False):
        '''Return the point at given fraction along the geodesic between
           this and an other point, using this C{Direct} and C{Inverse}
           methods.

           @arg other: The other point (C{LatLon}).
           @arg fraction: Fraction between both points ranging from
                          0, meaning this to 1, the other point (C{float}).
           @kwarg height: Optional height, overriding the fractional
                          height (C{meter}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Intermediate point (C{LatLon}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @raise UnitError: Invalid B{C{fraction}} or B{C{height}}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.

           @see: Methods L{distanceTo3} and L{destination}.
        '''
        t = self.distanceTo3(other, wrap=wrap)
        f = Scalar(fraction=fraction)
        h = self._havg(other, f=f) if height is None else Height(height)
        return self.destination(t.distance * f, t.initial, height=h)

    def _Inverse(self, other, wrap, **unused):  # azis=False, overloaded by I{Vincenty}
        '''(INTERNAL) I{Karney}'s C{Inverse} method.

           @return: A L{Distance3Tuple}C{(distance, initial, final)}.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's
                              L{Datum} ellipsoids are not compatible.
        '''
        _ = self.ellipsoids(other)
        g = self.geodesic
        _, lon = unroll180(self.lon, other.lon, wrap=wrap)
        return g.Inverse3(self.lat, self.lon, other.lat, lon)


def _Equidistant00(equidistant, p1):
    '''(INTERNAL) Get an C{Equidistant*(0, 0, ...)} instance.
    '''
    from pygeodesy.azimuthal import _Equidistants

    if equidistant is None or not callable(equidistant):
        equidistant = p1.Equidistant
    elif not issubclassof(equidistant, *_Equidistants):  # PYCHOK no cover
        t = tuple(_.__name__ for _ in _Equidistants)
        raise _IsnotError(*t, equidistant=equidistant)
    return equidistant(0, 0, p1.datum)


def _intersect3(s1, end1, s2, end2, height=None, wrap=True,  # MCCABE?
                equidistant=None, tol=_TOL_M, LatLon=None, **LatLon_kwds):
    '''(INTERNAL) Intersect two (ellipsoidal) path, see L{_intersection}
       above, separated to allow callers to embellish any exceptions.
    '''
    from pygeodesy.sphericalTrigonometry import _intersect as _si, LatLon as _LLS
    from pygeodesy.vector3d import _intersect3d3 as _vi3

    # E = s1.ellipsoids(s2)
    # assert E == s1.ellispoids(e1) == s1.ellipsoids(e2)

    e1 = s1.others(end1=end1)
    e2 = s1.others(end2=end2)

    if wrap:  # unroll180 == .karney._unroll2
        e1 = _unrollon(s1, e1)
        e2 = _unrollon(s2, e2)

    # get the azimuthal equidistant projection
    A = _Equidistant00(equidistant, s1)
    e = max(tol, EPS)

    # gu-/estimate initial intersection, spherically ...
    t = _si(_LLS(s1.lat, s1.lon, height=s1.height),
            _LLS(e1.lat, e1.lon, height=e1.height),
            _LLS(s2.lat, s2.lon, height=s2.height),
            _LLS(e2.lat, e2.lon, height=e2.height),
            height=height, wrap=False, LatLon=_LLS)  # unrolled already
    h, n = t.height, t.name

    # ... and iterate as Karney describes, @see:
    # LatLonEllipsoidalBase.LatLon.intersections2
    c = None  # force first d == c to False
    for i in range(1, _TRIPS):
        A.reset(t.lat, t.lon)  # gu-/estimate as origin
        # convert start and end points to projection
        # space and compute an intersection there
        v, o1, o2 = _vi3(A.forward(s1.lat, s1.lon),
                         A.forward(e1.lat, e1.lon),
                         A.forward(s2.lat, s2.lon),
                         A.forward(e2.lat, e2.lon),
                         eps=e, useZ=False)
        # convert intersection back to geodetic
        t, d = A._reverse2(v.x, v.y)
        # break if below tolerance or if unchanged
        if d < e or d == c:
            t._iteration = i
            break
        c = d
    else:
        raise IntersectionError(_no_(Fmt.convergence(e)))

    r = _LL4Tuple(t.lat, t.lon, h, t.datum, LatLon, LatLon_kwds, inst=s1, name=n)
    r._iteration = t._iteration  # _NamedTuple._iteration
    return Intersection3Tuple(r, o1, o2)


def _intersection3(start1, end1, start2, end2, height=None, wrap=True,
                   equidistant=None, tol=_TOL_M, LatLon=None, **LatLon_kwds):
    '''(INTERNAL) Iteratively compute the intersection point of two paths,
       each defined by an (ellipsoidal) start and end point.
    '''
    s1 = _xellipsoidal(start1=start1)
    s2 = s1.others(start2=start2)

    try:
        return _intersect3(s1, end1, s2, end2, height=height, wrap=wrap,
                                          equidistant=equidistant, tol=tol,
                                               LatLon=LatLon, **LatLon_kwds)
    except (TypeError, ValueError) as x:
        raise _xError(x, start1=start1, end1=end1, start2=start2, end2=end2)


def _intersections2(center1, radius1, center2, radius2, height=None, wrap=True,
                    equidistant=None, tol=_TOL_M, LatLon=None, **LatLon_kwds):
    '''(INTERNAL) Iteratively compute the intersection points of two circles,
       each defined by an (ellipsoidal) center point and a radius.
    '''
    c1 = _xellipsoidal(center1=center1)
    c2 = c1.others(center2=center2)

    try:
        return _intersects2(c1, radius1, c2, radius2, height=height, wrap=wrap,
                                                 equidistant=equidistant, tol=tol,
                                                      LatLon=LatLon, **LatLon_kwds)
    except (TypeError, ValueError) as x:
        raise _xError(x, center1=center1, radius1=radius1,
                         center2=center2, radius2=radius2)


def _intersects2(c1, radius1, c2, radius2, height=None, wrap=True,  # MCCABE 15
                 equidistant=None, tol=_TOL_M, LatLon=None, **LatLon_kwds):
    '''(INTERNAL) Intersect two (ellipsoidal) circles, see L{_intersections2}
       above, separated to allow callers to embellish any exceptions.
    '''
    from pygeodesy.sphericalTrigonometry import _intersects2 as _si2, LatLon as _LLS
    from pygeodesy.vector3d import _intersects2 as _vi2

    def _latlon4(t, h, n, c):
        r = _LL4Tuple(t.lat, t.lon, h, t.datum, LatLon, LatLon_kwds, inst=c, name=n)
        r._iteration = t.iteration  # ._iteration for tests
        return r

    r1 = Radius_(radius1=radius1)
    r2 = Radius_(radius2=radius2)

    E = c1.ellipsoids(c2)
    # get the azimuthal equidistant projection
    A = _Equidistant00(equidistant, c1)

    if r1 < r2:
        c1, c2 = c2, c1
        r1, r2 = r2, r1

    if r1 > (min(E.b, E.a) * PI):
        raise _ValueError(_exceed_PI_radians_)

    if wrap:  # unroll180 == .karney._unroll2
        c2 = _unrollon(c1, c2)

    # distance between centers and radii are
    # measured along the ellipsoid's surface
    m = c1.distanceTo(c2, wrap=False)  # meter
    if m < max(r1 - r2, EPS):
        raise IntersectionError(_near_concentric_)
    if fsum_(r1, r2, -m) < 0:
        raise IntersectionError(_too_(Fmt.distant(m)))

    f = _radical2(m, r1, r2).ratio  # "radical fraction"
    r = E.rocMean(favg(c1.lat, c2.lat, f=f))
    e = max(m2degrees(tol, radius=r), EPS)

    # gu-/estimate initial intersections, spherically ...
    t1, t2 = _si2(_LLS(c1.lat, c1.lon, height=c1.height), r1,
                  _LLS(c2.lat, c2.lon, height=c2.height), r2,
                   radius=r, height=height, wrap=False, too_d=m)  # unrolled already
    h, n = t1.height, t1.name

    # ... and iterate as Karney describes, @see:
    # LatLonEllipsoidalBase.LatLon.intersections2
    ts, ta = [], None
    for t in ((t1,) if t1 is t2 else (t1, t2)):
        c = None  # force first d == c to False
        for i in range(1, _TRIPS):
            A.reset(t.lat, t.lon)  # gu-/estimate as origin
            # convert centers to projection space
            t1 = A.forward(c1.lat, c1.lon)
            t2 = A.forward(c2.lat, c2.lon)
            # compute intersections in projection space
            v1, v2 = _vi2(t1, r1,  # XXX * t1.scale?,
                          t2, r2,  # XXX * t2.scale?,
                          sphere=False, too_d=m)
            # convert intersections back to geodetic
            t1, d1 = A._reverse2(v1.x, v1.y)
            if v1 is v2:  # abutting
                t, d = t1, d1  # PYCHOK no cover
            else:
                t2, d2 = A._reverse2(v2.x, v2.y)
                # consider only the closer intersection
                t, d = (t1, d1) if d1 < d2 else (t2, d2)
            # break if below tolerance or if unchanged
            if d < e or d == c:
                t._iteration = i  # _NamedTuple._iteration
                ts.append(t)
                if v1 is v2:  # abutting
                    ta = t  # PYCHOK no coves
                break
            c = d
        else:
            raise IntersectionError(_no_(Fmt.convergence(tol)))

    if ta:  # abutting circles
        pass  # PYCHOK no cover
    elif len(ts) == 2:
        return (_latlon4(ts[0], h, n, c1),
                _latlon4(ts[1], h, n, c2))
    elif len(ts) == 1:  # PYCHOK no cover
        ta = ts[0]  # assume abutting
    else:  # PYCHOK no cover
        raise _AssertionError(ts=ts)
    r = _latlon4(ta, h, n, c1)
    return r, r


def _nearestOn(point, point1, point2, within=True, height=None, wrap=True,
               equidistant=None, tol=_TOL_M, LatLon=None, **LatLon_kwds):
    '''(INTERNAL) Get closest point, imported by .ellipsoidalExact,
       -GeodSolve, -Karney and -Vincenty to embellish exceptions.
    '''
    try:
        p = _xellipsoidal(point=point)
        return _nearestOne(p, point1, point2, within=within, height=height, wrap=wrap,
                           equidistant=equidistant, tol=tol, LatLon=LatLon, **LatLon_kwds)
    except (TypeError, ValueError) as x:
        raise _xError(x, point=point, point1=point1, point2=point2)


def _nearestOne(p, point1, point2, within=True, height=None, wrap=True,
                equidistant=None, tol=_TOL_M, LatLon=None, **LatLon_kwds):
    '''(INTERNAL) Get closest point, like L{_intersects2} above,
       separated to allow callers to embellish any exceptions.
    '''
    from pygeodesy.sphericalNvector import LatLon as _LLS
    from pygeodesy.vector3d import _nearestOn as _vnOn, Vector3d

    def _v(t, h):
        return Vector3d(t.x, t.y, h)

    p1 = p.others(point1=point1)
    p2 = p.others(point2=point2)

    _ = p.ellipsoids(p1)
    E = p.ellipsoids(p2)

    if wrap:
        p1 = _unrollon(p,  p1)
        p2 = _unrollon(p,  p2)
        p2 = _unrollon(p1, p2)

    r = E.rocMean(fmean_(p.lat, p1.lat, p2.lat))
    e = max(m2degrees(tol, radius=r), EPS)

    # get the azimuthal equidistant projection
    A = _Equidistant00(equidistant, p)

    # gu-/estimate initial nearestOn, spherically ... wrap=False, only!
    t = _LLS(p.lat,  p.lon,  height=p.height).nearestOn(
        _LLS(p1.lat, p1.lon, height=p1.height),
        _LLS(p2.lat, p2.lon, height=p2.height), within=within, height=height)
    n = t.name
    if height is False:  # PYCHOK no cover
        h  = t.height  # use heights as Z
        h1 = p1.height
        h2 = p2.height
    else:
        h = h1 = h2 = _0_0

    # ... and iterate as Karney describes, @see:
    # LatLonEllipsoidalBase.LatLon.intersections2
    c = None  # force first d == c to False
    # closest to origin, .z to interpolate height
    p = Vector3d(0, 0, h)
    for i in range(1, _TRIPS):
        A.reset(t.lat, t.lon)  # gu-/estimate as origin
        # convert points to projection space
        # and compute the nearest one there
        v = _vnOn(p, _v(A.forward(p1.lat, p1.lon), h1),
                     _v(A.forward(p2.lat, p2.lon), h2),
                      within=within)
        # convert nearest one back to geodetic
        t, d = A._reverse2(v.x, v.y)
        # break if below tolerance or if unchanged
        if d < e or d == c:
            t._iteration = i  # _NamedTuple._iteration
            if height is False:
                h = v.z  # nearest interpolated
            break
        c = d
    else:
        raise _ValueError(_no_(Fmt.convergence(tol)))

    r = _LL4Tuple(t.lat, t.lon, h, t.datum, LatLon, LatLon_kwds, inst=p, name=n)
    r._iteration = t.iteration  # ._iteration for tests
    return r


def _unrollon(p1, p2):  # unroll180 == .karney._unroll2
    # wrap, unroll and replace longitude if different
    _, lon = unroll180(p1.lon, p2.lon, wrap=True)
    if abs(lon - p2.lon) > EPS:
        p2 = p2.classof(p2.lat, lon, **_xkwds_not(None,
                                   height=getattr(p2, _height_,  None),
                                    datum=getattr(p2, _datum_,   None),
                                    epoch=getattr(p2, _epoch_,   None),
                                  reframe=getattr(p2, _reframe_, None)))
    return p2


__all__ += _ALL_DOCS(LatLonEllipsoidalBaseDI)

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
