
# -*- coding: utf-8 -*-

u'''(INTERNAL) Private, ellipsoidal Direct/Inverse geodesy base
class C{LatLonEllipsoidalBaseDI} and functions.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import isscalar, issubclassof
from pygeodesy.constants import EPS, MAX, PI, PI2, PI_4, isnear0, isnear1, \
                               _0_0, _0_5, _1_5, _3_0
from pygeodesy.ellipsoidalBase import LatLonEllipsoidalBase, Property_RO, \
                                      property_RO, _TOL_M
from pygeodesy.errors import _AssertionError, IntersectionError, _IsnotError, \
                             _or, _ValueError, _xellipsoidal, _xError, _xkwds_not
from pygeodesy.fmath import favg, fmean_
from pygeodesy.fsums import Fmt, fsumf_
from pygeodesy.formy import opposing, _radical2
from pygeodesy.interns import _antipodal_, _concentric_, _exceed_PI_radians_, \
                              _low_, _near_, _SPACE_, _too_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.namedTuples import Bearing2Tuple, Destination2Tuple, \
                                  Intersection3Tuple, NearestOn2Tuple, \
                                  NearestOn8Tuple, _LL4Tuple
# from pygeodesy.props import Property_RO, property_RO  # from .ellipsoidalBase
# from pygeodesy.streprs import Fmt  # from .fsums
from pygeodesy.units import _fi_j2, Radius_, Scalar
from pygeodesy.utily import m2km, unroll180, _unrollon, _unrollon3, \
                           _Wrap, wrap360

from math import degrees, radians

__all__ = _ALL_LAZY.ellipsoidalBaseDI
__version__ = '23.05.15'

_polar__  = 'polar?'
_B2END    = _1_5  # _intersect3 bearing to pseudo-end point factor
_TRIPS    =  33   # _intersect3, _intersects2, _nearestOn interations, 6..9 sufficient?


class LatLonEllipsoidalBaseDI(LatLonEllipsoidalBase):
    '''(INTERNAL) Base class for C{ellipsoidal*.LatLon} classes
       with I{overloaded} C{Direct} and C{Inverse} methods.
    '''

    def bearingTo2(self, other, wrap=False):
        '''Compute the initial and final bearing (forward and reverse
           azimuth) from this to an other point, using this C{Inverse}
           method.  See methods L{initialBearingTo} and L{finalBearingTo}
           for more details.

           @arg other: The other point (C{LatLon}).
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{other}} point (C{bool}).

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
        d = LL(*_Wrap.latlon(r.lat, r.lon), height=h, **_xkwds_not(None,
                                            datum=self.datum, name=self.name,
                                            epoch=self.epoch, reframe=self.reframe))
        return Destination2Tuple(d, wrap360(r.final))

    def distanceTo(self, other, wrap=False, **unused):  # ignore radius=R_M
        '''Compute the distance between this and an other point along
           a geodesic, using this C{Inverse} method. See method
           L{distanceTo3} for more details.

           @arg other: The other point (C{LatLon}).
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{other}} point (C{bool}).

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
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{other}} point (C{bool}).

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
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                        the B{C{other}} point (C{bool}).

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
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                        the B{C{other}} point (C{bool}).

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
           @arg fraction: Fraction between both points (C{scalar}, 0.0
                          at this and 1.0 at the other point.
           @kwarg height: Optional height, overriding the fractional
                          height (C{meter}).
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{other}} point (C{bool}).

           @return: Intermediate point (C{LatLon}).

           @raise TypeError: The B{C{other}} point is not C{LatLon}.

           @raise UnitError: Invalid B{C{fraction}} or B{C{height}}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.

           @see: Methods L{distanceTo3}, L{destination}, C{midpointTo} and
                 C{rhumbMidpointTo}.
        '''
        f = Scalar(fraction=fraction)
        if isnear0(f):
            r = self
        elif isnear1(f) and not wrap:
            r = self.others(other)
        else:  # negative fraction OK
            t = self.distanceTo3(other, wrap=wrap)
            h = self._havg(other, f=f, h=height)
            r = self.destination(t.distance * f, t.initial, height=h)
        return r

    def _Inverse(self, other, wrap, **unused):  # azis=False, overloaded by I{Vincenty}
        '''(INTERNAL) I{Karney}'s C{Inverse} method.

           @return: A L{Distance3Tuple}C{(distance, initial, final)}.
        '''
        _ = self.ellipsoids(other)
        g = self.geodesic
        _, lon = unroll180(self.lon, other.lon, wrap=wrap)
        return g.Inverse3(self.lat, self.lon, other.lat, lon)

    def nearestOn8(self, points, closed=False, height=None, wrap=False,
                                               equidistant=None, tol=_TOL_M):
        '''I{Iteratively} locate the point on a path or polygon closest
           to this point.

           @arg points: The path or polygon points (C{LatLon}[]).
           @kwarg closed: Optionally, close the polygon (C{bool}).
           @kwarg height: Optional height, overriding the height of this and
                          all other points (C{meter}, conventionally).  If
                          B{C{height}} is C{None}, the height of each point
                          is taken into account for distances.
           @kwarg wrap: If C{True}, wrap or I{normalize} and unroll the
                        B{C{points}} (C{bool}).

           @return: A L{NearestOn8Tuple}C{(closest, distance, fi, j, start,
                    end, initial, final)} with C{distance} in C{meter},
                    conventionally and with the C{closest}, the C{start}
                    the C{end} point each an instance of this C{LatLon}.

           @raise PointsError: Insufficient number of B{C{points}}.

           @raise TypeError: Some B{C{points}} or B{C{equidistant}} invalid.

           @raise ValueError: Some B{C{points}}' datum or ellipsoid incompatible
                              or no convergence for the given B{C{tol}}.

           @see: Function L{pygeodesy.nearestOn6} and method C{nearestOn6}.
        '''
        _d3 =  self.distanceTo3  # Distance3Tuple
        _n3 = _nearestOn3
        try:
            Ps =  self.PointsIter(points, loop=1, wrap=wrap)
            p1 =  c = s = e = Ps[0]
            _  =  self.ellipsoids(p1)
            c3 = _d3(c, wrap=wrap)  # XXX wrap=False?

        except (TypeError, ValueError) as x:
            raise _xError(x, Fmt.SQUARE(points=0), p1, this=self, tol=tol,
                             closed=closed, height=height, wrap=wrap)

        # get the azimuthal equidistant projection, once
        A = _Equidistant00(equidistant, c)
        b = _Box(c, c3.distance)
        m = f = i = 0  # p1..p2 == points[i]..[j]

        kwds = dict(within=True, height=height, tol=tol,
                    LatLon=self.classof,  # this LatLon
                    datum=self.datum, epoch=self.epoch, reframe=self.reframe)
        try:
            for j, p2 in Ps.enumerate(closed=closed):
                if wrap and j != 0:  # not Ps.looped
                    p2 = _unrollon(p1, p2)
                # skip edge if no overlap with box around closest
                if j < 4 or b.overlaps(p1.lat, p1.lon, p2.lat, p2.lon):
                    p, t, _ = _n3(self, p1, p2, A, **kwds)
                    d3 = _d3(p, wrap=False)  # already unrolled
                    if d3.distance < c3.distance:
                        c3, c, s, e, f = d3, p, p1, p2, (i + t)
                        b = _Box(c, c3.distance)
                        m =  max(m, c.iteration)
                p1, i = p2, j

        except (TypeError, ValueError) as x:
            raise _xError(x, Fmt.SQUARE(points=i), p1,
                             Fmt.SQUARE(points=j), p2, this=self, tol=tol,
                             closed=closed, height=height, wrap=wrap)

        f, j = _fi_j2(f, len(Ps))  # like .vector3d.nearestOn6

        n = self.nearestOn8.__name__
        c.rename(n)
        if s is not c:
            s = s.copy(name=n)
        if e is not c:
            e = e.copy(name=n)
        return NearestOn8Tuple(c, c3.distance, f, j, s, e, c3.initial, c3.final,
                               iteration=m)  # ._iteration for tests


class _Box(object):
    '''Bounding box around a C{LatLon} point.

       @see: Function C{_box4} in .clipy.py.
    '''
    def __init__(self, center, distance):
        '''New L{_Box} around point.

           @arg center: The center point (C{LatLon}).
           @arg distance: Radius, half-size of the box
                          (C{meter}, conventionally)
        '''
        E = center.ellipsoid()
        d = degrees(distance / max(E.a, E.b)) + _0_5  # some margin
        self._N = center.lat + d
        self._S = center.lat - d
        self._E = center.lon + d
        self._W = center.lon - d

    def overlaps(self, lat1, lon1, lat2, lon2):
        '''Check whether this box overlaps a line between 2 points.

           @arg lat1: Latitude of first point (C{degrees}).
           @arg lon1: Longitude of first point (C{degrees}).
           @arg lat2: Latitude of second point (C{degrees}).
           @arg lon2: Longitude of second point (C{degrees}).

           @return: C{False} if there is certainly no overlap,
                    C{True} otherwise (C{bool}).
        '''
        if lat1 > lat2:
            lat1, lat2 = lat2, lat1
        if lat1 > self._N or lat2 < self._S:
            return False
        if lon1 > lon2:
            lon1, lon2 = lon2, lon1
        if lon1 > self._E or lon2 < self._W:
            return False
        return True


class _Tol(object):
    '''Handle a tolerance in C{meter} as C{degrees} and C{meter}.
    '''
    _deg  = 0    # tol in degrees
    _lat  = 0
    _m    = 0    # tol in meter
    _min  = MAX  # degrees
    _prev = None
    _r    = 0

    def __init__(self, tol_m, E, lat, *lats):
        '''New L{_Tol}.

           @arg tol_m: Tolerance (C{meter}, only).
           @arg E: Earth ellipsoid (L{Ellipsoid}).
           @arg lat: Latitude (C{degrees}).
           @arg lats: Additional latitudes (C{degrees}).
        '''
        self._lat = fmean_(lat, *lats) if lats else lat
        self._r   = max(EPS, E.rocMean(self._lat))
        self._m   = max(EPS, tol_m)
        self._deg = max(EPS, degrees(self._m / self._r))  # avoid m2degrees!

    @property_RO
    def degrees(self):
        '''Get this tolerance in C{degrees}.
        '''
        return self._deg

    def degrees2m(self, deg):
        '''Convert B{C{deg}} to meter at the same C{lat} and earth radius.
        '''
        return self.radius * radians(deg) / PI2  # avoid degrees2m!

    def degError(self, Error=_ValueError):
        '''Compose an error with the C{deg}rees minimum.
        '''
        return self.mError(self.degrees2m(self._min), Error=Error)

    def done(self, deg):
        '''Check C{deg} vs. tolerance and previous value.
        '''
        if deg < self._deg or deg == self._prev:
            return True
        self._min  = min(self._min, deg)
        self._prev = deg
        return False

    @property_RO
    def lat(self):
        '''Get the mean latitude  in C{degrees}.
        '''
        return self._lat

    def mError(self, m, Error=_ValueError):
        '''Compose an error with B{C{m}}eter minimum.
        '''
        t = _SPACE_(Fmt.tolerance(self.meter), _too_(_low_))
        if m2km(m) > self.meter:
            t = _or(t, _antipodal_, _near_(_polar__))
        return Error(Fmt.no_convergence(m), txt=t)

    @property_RO
    def meter(self):
        '''Get this tolerance in C{meter}.
        '''
        return self._m

    @property_RO
    def radius(self):
        '''Get the earth radius in C{meter}.
        '''
        return self._r

    def reset(self):
        '''Reset tolerances.
        '''
        self._min  = MAX
        self._prev = None


def _Equidistant00(equidistant, p1):
    '''(INTERNAL) Get an C{Equidistant*(0, 0, ...)} instance.
    '''
    if equidistant is None or not callable(equidistant):
        equidistant = p1.Equidistant
    elif not issubclassof(equidistant, *_MODS.azimuthal._Equidistants):  # PYCHOK no cover
        t = tuple(_.__name__ for _  in  _MODS.azimuthal._Equidistants)
        raise _IsnotError(*t, equidistant=equidistant)
    return equidistant(0, 0, p1.datum)


def _intersect3(s1, end1, s2, end2, height=None, wrap=False,  # MCCABE 16  was=True
                equidistant=None, tol=_TOL_M, LatLon=None, **LatLon_kwds):
    '''(INTERNAL) Intersect two (ellipsoidal) lines, see ellipsoidal method
       L{intersection3}, separated to allow callers to embellish any exceptions.
    '''
    _LLS = _MODS.sphericalTrigonometry.LatLon
    _si  = _MODS.sphericalTrigonometry._intersect
    _vi3 = _MODS.vector3d._intersect3d3

    def _b_d(s, e, w, t, h=_0_0):
        # compute opposing and distance
        t = s.classof(t.lat, t.lon, height=h, name=t.name)
        t = s.distanceTo3(t, wrap=w)  # Distance3Tuple
        b = opposing(e, t.initial)  # "before" start
        return b, t.distance

    def _b_e(s, e, w, t):
        # compute an end point along the initial bearing
        # about 1.5 times the distance to the gu-/estimate, at
        # least 1/8 and at most 3/8 of the earth perimeter like
        # radians in .sphericalTrigonometry._int3d2 and bearing
        # comparison in .sphericalTrigonometry._intb
        b, d = _b_d(s, e, w, t, h=t.height)
        m = s.ellipsoid().R2x * PI_4  # authalic exact
        d = min(max(d * _B2END, m), m * _3_0)
        e = s.destination(d, e)
        return b, (_unrollon(s, e) if w else e)

    def _e_ll(s, e, w, **end):
        # return 2-tuple (end, False if bearing else True)
        ll = not isscalar(e)
        if ll:
            e = s.others(**end)
            if w:  # unroll180 == .karney._unroll2
                e = _unrollon(s, e)
        return e, ll

    def _llS(ll):  # return an C{_LLS} instance
        return _LLS(ll.lat, ll.lon, height=ll.height)

    def _o(o, b, n, s, t, e):
        # determine C{o}utside before, on or after start point
        if not o:  # intersection may be on start
            if _MODS.latlonBase._isequalTo(s, t, eps=e.degrees):
                return o
        return -n if b else n

    E = s1.ellipsoids(s2)

    e1, ll1 = _e_ll(s1, end1, wrap, end1=end1)
    e2, ll2 = _e_ll(s2, end2, wrap, end2=end2)

    e = _Tol(tol, E, s1.lat, (e1.lat if ll1 else s1.lat),
                     s2.lat, (e2.lat if ll2 else s2.lat))

    # get the azimuthal equidistant projection
    A = _Equidistant00(equidistant, s1)

    # gu-/estimate initial intersection, spherically ...
    t = _si(_llS(s1), (_llS(e1) if ll1 else e1),
            _llS(s2), (_llS(e2) if ll2 else e2),
            height=height, wrap=False, LatLon=_LLS)  # unrolled already
    h, n = t.height, t.name

    if not ll1:
        b1, e1 = _b_e(s1, e1, wrap, t)
    if not ll2:
        b2, e2 = _b_e(s2, e2, wrap, t)

    # ... and iterate as Karney describes, for references
    # @see: Function L{ellipsoidalKarney.intersection3}.
    for i in range(1, _TRIPS):
        A.reset(t.lat, t.lon)  # gu-/estimate as origin
        # convert start and end points to projection
        # space and compute an intersection there
        v, o1, o2 = _vi3(*A._forwards(s1, e1, s2, e2),
                          eps=e.meter, useZ=False)
        # convert intersection back to geodetic
        t, d = A._reverse2(v)
        if e.done(d):  # below tol or unchanged?
            t._iteration = i
            break
    else:
        raise e.degError(Error=IntersectionError)

    # like .sphericalTrigonometry._intersect, if this intersection
    # is "before" the first point, use the antipodal intersection
    if not (ll1 or ll2):  # end1 and end2 are an initial bearing
        b1, _ = _b_d(s1, end1, wrap, t)
        if b1:
            t  = t.antipodal()
            b1 = not b1
        b2, _ = _b_d(s2, end2, wrap, t)

    r = _LL4Tuple(t.lat, t.lon, h, t.datum, LatLon, LatLon_kwds, inst=s1,
                                            iteration=t._iteration, name=n)
    return Intersection3Tuple(r, (o1 if ll1 else _o(o1, b1, 1, s1, t, e)),
                                 (o2 if ll2 else _o(o2, b2, 2, s2, t, e)))


def _intersection3(start1, end1, start2, end2, height=None, wrap=False,  # was=True
                   equidistant=None, tol=_TOL_M, LatLon=None, **LatLon_kwds):
    '''(INTERNAL) Iteratively compute the intersection point of two lines,
       each defined by two (ellipsoidal) points or an (ellipsoidal) start
       point and an initial bearing from North.
    '''
    s1 = _xellipsoidal(start1=start1)
    s2 =  s1.others(start2=start2)
    try:
        return _intersect3(s1, end1, s2, end2, height=height, wrap=wrap,
                                          equidistant=equidistant, tol=tol,
                                               LatLon=LatLon, **LatLon_kwds)
    except (TypeError, ValueError) as x:
        raise _xError(x, start1=start1, end1=end1, start2=start2, end2=end2)


def _intersections2(center1, radius1, center2, radius2, height=None, wrap=False,  # was=True
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


def _intersects2(c1, radius1, c2, radius2, height=None, wrap=False,  # MCCABE 16 was=True
                 equidistant=None, tol=_TOL_M, LatLon=None, **LatLon_kwds):
    '''(INTERNAL) Intersect two (ellipsoidal) circles, see L{_intersections2}
       above, separated to allow callers to embellish any exceptions.
    '''
    _LLS = _MODS.sphericalTrigonometry.LatLon
    _si2 = _MODS.sphericalTrigonometry._intersects2
    _vi2 = _MODS.vector3d._intersects2

    def _latlon4(t, h, n, c):
        return _LL4Tuple(t.lat, t.lon, h, t.datum, LatLon, LatLon_kwds, inst=c,
                                                   iteration=t.iteration, name=n)

    def _llS(ll):  # return an C{_LLS} instance
        return _LLS(ll.lat, ll.lon, height=ll.height)

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
        raise IntersectionError(_near_(_concentric_))  # XXX ConcentricError?
    if fsumf_(r1, r2, -m) < 0:
        raise IntersectionError(_too_(Fmt.distant(m)))

    f = _radical2(m, r1, r2).ratio  # "radical fraction"
    e = _Tol(tol, E, favg(c1.lat, c2.lat, f=f))

    # gu-/estimate initial intersections, spherically ...
    t1, t2 = _si2(_llS(c1), r1, _llS(c2), r2, radius=e.radius,
                   height=height, wrap=False, too_d=m)  # unrolled already
    h, n = t1.height, t1.name

    # ... and iterate as Karney describes, for references
    # @see: Function L{ellipsoidalKarney.intersections2}.
    ts, ta = [], None
    for t in ((t1,) if t1 is t2 else (t1, t2)):
        for i in range(1, _TRIPS):
            A.reset(t.lat, t.lon)  # gu-/estimate as origin
            # convert centers to projection space and
            # compute the intersections there
            t1, t2 = A._forwards(c1, c2)
            v1, v2 = _vi2(t1, r1,  # XXX * t1.scale?,
                          t2, r2,  # XXX * t2.scale?,
                          sphere=False, too_d=m)
            # convert intersections back to geodetic
            if v1 is v2:  # abutting
                t,  d  = A._reverse2(v1)
            else:  # consider the closer intersection
                t1, d1 = A._reverse2(v1)
                t2, d2 = A._reverse2(v2)
                t, d = (t1, d1) if d1 < d2 else (t2, d2)
            if e.done(d):  # below tol or unchanged?
                t._iteration = i  # _NamedTuple._iteration
                ts.append(t)
                if v1 is v2:  # abutting
                    ta = t
                break
        else:
            raise e.degError(Error=IntersectionError)
        e.reset()

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


def _nearestOn2(p, point1, point2, within=True, height=None, wrap=False,  # was=True
                   equidistant=None, tol=_TOL_M, **LatLon_and_kwds):
    '''(INTERNAL) Closest point and fraction, like L{_intersects2} above,
       separated to allow callers to embellish any exceptions.
    '''
    p1 = p.others(point1=point1)
    p2 = p.others(point2=point2)

    _ = p.ellipsoids(p1)
#   E = p.ellipsoids(p2)  # done in _nearestOn3

    # get the azimuthal equidistant projection
    A = _Equidistant00(equidistant, p)

    p1, p2, _ = _unrollon3(p, p1, p2, wrap)  # XXX don't unroll?
    r,  f,  _ = _nearestOn3(p, p1, p2, A, within=within, height=height,
                                          tol=tol, **LatLon_and_kwds)
    return NearestOn2Tuple(r, f)


def _nearestOn3(p, p1, p2, A, within=True, height=None, tol=_TOL_M,
                              LatLon=None, **LatLon_kwds):
    # Only in function C{_nearestOn2} and method C{nearestOn8} above
    _LLS = _MODS.sphericalNvector.LatLon
    _V3d = _MODS.vector3d.Vector3d
    _vn2 = _MODS.vector3d._nearestOn2

    def _llS(ll):  # return an C{_LLS} instance
        return _LLS(ll.lat, ll.lon, height=ll.height)

    E =  p.ellipsoids(p2)
    e = _Tol(tol, E, p.lat, p1.lat, p2.lat)

    # gu-/estimate initial nearestOn, spherically ... wrap=False, only!
    # using sphericalNvector.LatLon.nearestOn for within=False support
    t = _llS(p).nearestOn(_llS(p1), _llS(p2), within=within,
                                              height=height)
    n, h = t.name, t.height
    if height is None:
        h1 = p1.height  # use heights as pseudo-Z in projection space
        h2 = p2.height  # to be included in the closest function
        h0 = favg(h1, h2)
    else:  # ignore heights in distances, Z=0
        h0 = h1 = h2 = _0_0

    # ... and iterate to find the closest (to the origin with .z
    # to interpolate height) as Karney describes, for references
    # @see: Function L{ellipsoidalKarney.nearestOn}.
    vp, f = _V3d(_0_0, _0_0, h0), None
    for i in range(1, _TRIPS):
        A.reset(t.lat, t.lon)  # gu-/estimate as origin
        # convert points to projection space and compute
        # the nearest one (and its height) there
        s, t = A._forwards(p1, p2)
        v, f = _vn2(vp, _V3d(s.x, s.y, h1),
                        _V3d(t.x, t.y, h2), within=within)
        # convert nearest one back to geodetic
        t, d = A._reverse2(v)
        if e.done(d):  # below tol or unchanged?
            t._iteration = i  # _NamedTuple._iteration
            break
    else:
        raise e.degError()

    if height is None:
        h = v.z  # nearest
    elif isscalar(height):
        h = height
    r = _LL4Tuple(t.lat, t.lon, h, t.datum, LatLon, LatLon_kwds, inst=p,
                                            iteration=t.iteration, name=n)
    return r, f, e  # fraction or None


__all__ += _ALL_DOCS(LatLonEllipsoidalBaseDI)

# **) MIT License
#
# Copyright (C) 2016-2023 -- mrJean1 at Gmail -- All Rights Reserved.
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
