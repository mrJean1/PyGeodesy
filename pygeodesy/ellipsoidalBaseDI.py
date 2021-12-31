
# -*- coding: utf-8 -*-

u'''(INTERNAL) Ellipsoidal direct/inverse geodesy base class
C{LatLonEllipsoidalBaseDI} and functions.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # PYCHOK semicolon

from pygeodesy.basics import isnear0, isnear1, isscalar, issubclassof
from pygeodesy.ellipsoidalBase import LatLonEllipsoidalBase, Property_RO, \
                                      property_RO, _TOL_M
from pygeodesy.errors import _AssertionError, IntersectionError, _IsnotError, \
                             _ValueError, _xellipsoidal, _xError, _xkwds_not
from pygeodesy.fmath import favg, fmean_, Fmt, fsum_
from pygeodesy.formy import opposing, _radical2
from pygeodesy.interns import EPS, PI, PI2, PI_4, _antipodal_, _COMMASPACE_, \
                             _concentric_, _exceed_PI_radians_, _near_, _no_, \
                             _SPACE_, _too_, _0_0, _0_5, _1_5, _3_0
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.namedTuples import Bearing2Tuple, Destination2Tuple, \
                                  Intersection3Tuple, NearestOn2Tuple, \
                                  NearestOn8Tuple, _LL4Tuple
# from pygeodesy.props import Property_RO, property_RO  # from .ellipsoidalBase
# from pygeodesy.streprs import Fmt  # from .fmath
from pygeodesy.units import _fi_j2, Height, Radius_, Scalar
from pygeodesy.utily import m2km, unroll180, _unrollon, wrap90, wrap180, wrap360

from math import degrees, radians

__all__ = _ALL_LAZY.ellipsoidalBaseDI
__version__ = '21.12.28'

_polar__    = 'polar?'
_tolerance_ = 'tolerance'
_too_low_   = _too_('low')
_B2END      = _1_5  # _intersect3 bearing to pseudo-end point factor
_TRIPS      = 33    # _intersect3, _intersects2, _nearestOn interations, 6..9 sufficient?


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
           @arg fraction: Fraction between both points (C{scalar},
                          0.0 at this and 1.0 at the other point.
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
        f = Scalar(fraction=fraction)
        if isnear0(f):
            r = self
        elif isnear1(f) and not wrap:
            r = self.others(other)
        else:  # negative fraction OK
            t = self.distanceTo3(other, wrap=wrap)
            h = self._havg(other, f=f) if height is None else Height(height)
            r = self.destination(t.distance * f, t.initial, height=h)
        return r

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

    def nearestOn8(self, points, closed=False, height=None, wrap=False,
                                               equidistant=None, tol=_TOL_M):
        '''Iteratively locate the point on a path or polygon closest
           to this point.

           @arg points: The path or polygon points (C{LatLon}[]).
           @kwarg closed: Optionally, close the polygon (C{bool}).
           @kwarg height: Optional height, overriding the height of this and
                          all other points (C{meter}, conventionally).  If
                          B{C{height}} is C{None}, the height of each point
                          is taken into account for distances.

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
        D3 = self.distanceTo3  # Distance3Tuple

        try:
            Ps = self.PointsIter(points, loop=1)
            p1 = c = s = e = Ps[0]
            _  = self.ellipsoids(p1)
            c3 = D3(c, wrap=wrap)  # XXX wrap=False?

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
                if wrap and j != 0:
                    p2 = _unrollon(p1, p2)
                # skip edge if no overlap with box around closest
                if j < 4 or b.overlaps(p1.lat, p1.lon, p2.lat, p2.lon):
                    p, t, _ = _nearestOn3_(self, p1, p2, A, **kwds)
                    d3 = D3(p, wrap=False)  # already unrolled
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
        r = NearestOn8Tuple(c, c3.distance, f, j, s, e, c3.initial, c3.final)
        r._iteration = m  # ._iteration for tests
        return r


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
        non_ = ((lat1 > self._N or lat2 < self._S) if lat1 < lat2 else
                (lat2 > self._N or lat1 < self._S)) or \
               ((lon1 > self._E or lon2 < self._W) if lon1 < lon2 else
                (lon2 > self._E or lon1 < self._W))
        return not non_


class _Tol(object):
    '''Handle a tolerance in C{meter} as C{degrees} and C{meter}.
    '''
    _deg = 0
    _lat = 0
    _m   = 0
    _r   = 0

    def __init__(self, tol_m, E, lat, *lats):
        '''New L{_Tol}.

           @arg tol_m: Tolerance (C{meter}, only).
           @arg E: Earth ellispoid (L{Ellipsoid}).
           @arg lat: Latitude (C{degrees}).
           @arg lats: Additional latitudes (C{degrees}).
        '''
        if lats:
            lat = fmean_(lat, *lats)
        self._lat = lat
        self._m   = max(EPS, tol_m)
        self._r   = max(EPS, E.rocMean(lat))
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

    def degError(self, deg, Error=_ValueError):
        '''Compose an error with C{deg}rees minimum.
        '''
        return self.mError(self.degrees2m(deg), Error=Error)

    @property_RO
    def lat(self):
        '''Get the mean latitude  in C{degrees}.
        '''
        return self._lat

    def mError(self, m, Error=_ValueError):
        '''Compose an error with B{C{m}}eter minimum.
        '''
        t = _SPACE_(_tolerance_, Fmt.PAREN(Fmt.g(self.meter)), _too_low_)
        if m2km(m) > self.meter:
            t = _COMMASPACE_(t, _antipodal_, _near_(_polar__))
        return Error(_no_(Fmt.convergence(m)), txt=t)

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


def _Equidistant00(equidistant, p1):
    '''(INTERNAL) Get an C{Equidistant*(0, 0, ...)} instance.
    '''
    if equidistant is None or not callable(equidistant):
        equidistant = p1.Equidistant
    elif not issubclassof(equidistant, *_MODS.azimuthal._Equidistants):  # PYCHOK no cover
        t = tuple(_.__name__ for _  in  _MODS.azimuthal._Equidistants)
        raise _IsnotError(*t, equidistant=equidistant)
    return equidistant(0, 0, p1.datum)


def _intersect3(s1, end1, s2, end2, height=None, wrap=True,  # MCCABE 16
                equidistant=None, tol=_TOL_M, LatLon=None, **LatLon_kwds):
    '''(INTERNAL) Intersect two (ellipsoidal) path, see ellipsoidal method
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
    t = _si(_LLS(s1.lat, s1.lon, height=s1.height),
           (_LLS(e1.lat, e1.lon, height=e1.height) if ll1 else e1),
            _LLS(s2.lat, s2.lon, height=s2.height),
           (_LLS(e2.lat, e2.lon, height=e2.height) if ll2 else e2),
            height=height, wrap=False, LatLon=_LLS)  # unrolled already
    h, n = t.height, t.name

    if not ll1:
        b1, e1 = _b_e(s1, e1, wrap, t)
    if not ll2:
        b2, e2 = _b_e(s2, e2, wrap, t)

    # ... and iterate as Karney describes, @see:
    # LatLonEllipsoidalBase.LatLon.intersections2
    c = m = None  # force first d == c to False
    for i in range(1, _TRIPS):
        A.reset(t.lat, t.lon)  # gu-/estimate as origin
        # convert start and end points to projection
        # space and compute an intersection there
        v, o1, o2 = _vi3(A.forward(s1.lat, s1.lon),
                         A.forward(e1.lat, e1.lon),
                         A.forward(s2.lat, s2.lon),
                         A.forward(e2.lat, e2.lon),
                         eps=e.meter, useZ=False)
        # convert intersection back to geodetic
        t, d = A._reverse2(v.x, v.y)
        # break if below tolerance or if unchanged
        if d < e.degrees or d == c:
            t._iteration = i
            break
        if m is None or m > d:
            m = d  # min tol in degrees
        c = d
    else:
        raise e.degError(m, Error=IntersectionError)

    # like .sphericalTrigonometry._intersect, if this intersection
    # is "before" the first point, use the antipodal intersection
    if not (ll1 or ll2):  # end1 and end2 are an initial bearing
        b1, _ = _b_d(s1, end1, wrap, t)
        if b1:
            t  = t.antipodal()
            b1 = not b1
        b2, _ = _b_d(s2, end2, wrap, t)

    r = _LL4Tuple(t.lat, t.lon, h, t.datum, LatLon, LatLon_kwds, inst=s1, name=n)
    r._iteration = t._iteration  # _NamedTuple._iteration
    return Intersection3Tuple(r, (o1 if ll1 else _o(o1, b1, 1, s1, t, e)),
                                 (o2 if ll2 else _o(o2, b2, 2, s2, t, e)))


def _intersection3(start1, end1, start2, end2, height=None, wrap=True,
                   equidistant=None, tol=_TOL_M, LatLon=None, **LatLon_kwds):
    '''(INTERNAL) Iteratively compute the intersection point of two paths,
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


def _intersects2(c1, radius1, c2, radius2, height=None, wrap=True,  # MCCABE 16
                 equidistant=None, tol=_TOL_M, LatLon=None, **LatLon_kwds):
    '''(INTERNAL) Intersect two (ellipsoidal) circles, see L{_intersections2}
       above, separated to allow callers to embellish any exceptions.
    '''
    _LLS = _MODS.sphericalTrigonometry.LatLon
    _si2 = _MODS.sphericalTrigonometry._intersects2
    _vi2 = _MODS.vector3d._intersects2

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
        raise IntersectionError(_near_(_concentric_))
    if fsum_(r1, r2, -m) < 0:
        raise IntersectionError(_too_(Fmt.distant(m)))

    f = _radical2(m, r1, r2).ratio  # "radical fraction"
    e = _Tol(tol, E, favg(c1.lat, c2.lat, f=f))

    # gu-/estimate initial intersections, spherically ...
    t1, t2 = _si2(_LLS(c1.lat, c1.lon, height=c1.height), r1,
                  _LLS(c2.lat, c2.lon, height=c2.height), r2,
                   radius=e.radius, height=height, wrap=False, too_d=m)  # unrolled already
    h, n = t1.height, t1.name

    # ... and iterate as Karney describes, @see:
    # LatLonEllipsoidalBase.LatLon.intersections2
    ts, ta = [], None
    for t in ((t1,) if t1 is t2 else (t1, t2)):
        c = m = None  # force first d == c to False
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
            if d < e.degrees or d == c:
                t._iteration = i  # _NamedTuple._iteration
                ts.append(t)
                if v1 is v2:  # abutting
                    ta = t  # PYCHOK no coves
                break
            c = d
            if m is None or m > d:
                m = d  # min tol in degrees
        else:
            raise e.degError(m, Error=IntersectionError)

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


def _nearestOn2(p, point1, point2, within=True, height=None, wrap=True,
                   equidistant=None, tol=_TOL_M, **LatLon_and_kwds):
    '''(INTERNAL) Closest point and fraction, like L{_intersects2} above,
       separated to allow callers to embellish any exceptions.
    '''
    p1 = p.others(point1=point1)
    p2 = p.others(point2=point2)

    _ = p.ellipsoids(p1)
#   E = p.ellipsoids(p2)  # done in __nearestOn2__

    # get the azimuthal equidistant projection
    A = _Equidistant00(equidistant, p)

    if wrap:
        p1 = _unrollon(p,  p1)  # XXX do not unroll?
        p2 = _unrollon(p,  p2)  # XXX do not unroll?
        p2 = _unrollon(p1, p2)

    r, f, e = _nearestOn3_(p, p1, p2, A, within=within, height=height,
                                         tol=tol, **LatLon_and_kwds)
    return NearestOn2Tuple(r, f)


def _nearestOn3_(p, p1, p2, A, within=True, height=None, tol=_TOL_M,
                               LatLon=None, **LatLon_kwds):
    # Only function C{_nearestOn2} and method C{nearestOn8} above
    _LLS   = _MODS.sphericalNvector.LatLon
    _vnOn2 = _MODS.vector3d._nearestOn2
    _V3d   = _MODS.vector3d.Vector3d

    def _v3d(t, h):
        return _V3d(t.x, t.y, h)

    E =  p.ellipsoids(p2)
    e = _Tol(tol, E, p.lat, p1.lat, p2.lat)

    # gu-/estimate initial nearestOn, spherically ... wrap=False, only!
    # using sphericalNvector.LatLon.nearestOn for within=False support
    t = _LLS(p.lat,  p.lon,  height=p.height).nearestOn(
        _LLS(p1.lat, p1.lon, height=p1.height),
        _LLS(p2.lat, p2.lon, height=p2.height), within=within, height=height)
    n, h = t.name, t.height
    if height is None:
        h1 = p1.height  # use heights as pseudo-Z in projection space
        h2 = p2.height  # to be included in the closest function
        h0 = favg(h1, h2)
    else:  # ignore heights in distances, Z=0
        h0 = h1 = h2 = _0_0

    # ... and iterate as Karney describes, @see:
    # LatLonEllipsoidalBase.LatLon.intersections2
    c = m = f = None  # force first d == c to False
    # closest to origin, .z to interpolate height
    vp = _V3d(_0_0, _0_0, h0)
    for i in range(1, _TRIPS):
        A.reset(t.lat, t.lon)  # gu-/estimate as origin
        # convert points to projection space
        # and compute the nearest one there
        v, f = _vnOn2(vp, _v3d(A.forward(p1.lat, p1.lon), h1),
                          _v3d(A.forward(p2.lat, p2.lon), h2),
                           within=within)
        # convert nearest one back to geodetic
        t, d = A._reverse2(v.x, v.y)
        # break if below tolerance or if unchanged
        if d < e.degrees or d == c:
            t._iteration = i  # _NamedTuple._iteration
            break
        c = d
        if m is None or m > d:
            m = d  # min tol in degrees
    else:
        raise e.degError(m)

    if height is None:
        h = v.z  # nearest
    elif isscalar(height):
        h = height
    r = _LL4Tuple(t.lat, t.lon, h, t.datum, LatLon, LatLon_kwds, inst=p, name=n)
    r._iteration = t.iteration  # for tests
    return r, f, e  # fraction or None


__all__ += _ALL_DOCS(LatLonEllipsoidalBaseDI)

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
