
# -*- coding: utf-8 -*-

u'''Ellipsoidal I{Thaddeus Vincenty}'s geodetic (lat-/longitude) L{LatLon},
geocentric (ECEF) L{Cartesian} and L{VincentyError} classes and functions
L{areaOf}, L{intersections2}, L{nearestOn} and L{perimeterOf}.

Pure Python implementation of geodesy tools for ellipsoidal earth models,
transcribed from JavaScript originals by I{(C) Chris Veness 2005-2016}
and published under the same MIT Licence**, see U{Vincenty geodesics
<https://www.Movable-Type.co.UK/scripts/LatLongVincenty.html>}.  More at
U{geographiclib<https://PyPI.org/project/geographiclib>} and
U{GeoPy<https://PyPI.org/project/geopy>}.

Calculate geodesic distance between two points using the U{Vincenty
<https://WikiPedia.org/wiki/Vincenty's_formulae>} formulae and one of
several ellipsoidal earth models.  The default model is WGS-84, the
most accurate and widely used globally-applicable model for the earth
ellipsoid.

Other ellipsoids offering a better fit to the local geoid include Airy
(1830) in the UK, Clarke (1880) in Africa, International 1924 in much of
Europe, and GRS-67 in South America.  North America (NAD83) and Australia
(GDA) use GRS-80, which is equivalent to the WGS-84 model.

Great-circle distance uses a spherical model of the earth with the mean
earth radius defined by the International Union of Geodesy and Geophysics
(IUGG) as M{(2 * a + b) / 3 = 6371008.7714150598} meter or approx. 6371009
meter (for WGS-84, resulting in an error of up to about 0.5%).

Here's an example usage of C{ellipsoidalVincenty}:

    >>> from pygeodesy.ellipsoidalVincenty import LatLon
    >>> Newport_RI = LatLon(41.49008, -71.312796)
    >>> Cleveland_OH = LatLon(41.499498, -81.695391)
    >>> Newport_RI.distanceTo(Cleveland_OH)
    866,455.4329158525  # meter

You can change the ellipsoid model used by the Vincenty formulae
as follows:

    >>> from pygeodesy import Datums
    >>> from pygeodesy.ellipsoidalVincenty import LatLon
    >>> p = LatLon(0, 0, datum=Datums.OSGB36)

or by converting to anothor datum:

    >>> p = p.toDatum(Datums.OSGB36)

@newfield example: Example, Examples
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division

from pygeodesy.datums import Datums
from pygeodesy.ellipsoidalBase import _intermediateTo, _intersections2, \
                                       CartesianEllipsoidalBase, \
                                       LatLonEllipsoidalBase, _nearestOn
from pygeodesy.errors import _ValueError, _xellipsoidal, _xkwds
from pygeodesy.fmath import fpolynomial, hypot, hypot1
from pygeodesy.interns import EPS0, EPS, NN, _ambiguous_, _convergence_, _no_, \
                             _to_, _0_0, _1_0, _2_0, _3_0, _4_0, _6_0, _16_0
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_OTHER
from pygeodesy.namedTuples import Bearing2Tuple, Destination2Tuple, \
                                  Distance3Tuple
from pygeodesy.points import ispolar  # PYCHOK exported
from pygeodesy.props import property_doc_, property_RO
from pygeodesy.units import Number_, Scalar_, _1mm as _TOL_M
from pygeodesy.utily import atan2b, degrees90, degrees180, \
                            sincos2, unroll180

from math import atan2, cos, radians, tan

__all__ = _ALL_LAZY.ellipsoidalVincenty
__version__ = '21.01.14'

_antipodal_ = 'antipodal '  # trailing _SPACE_
_limit_     = 'limit'  # PYCHOK used!


class VincentyError(_ValueError):
    '''Error raised from I{Vincenty}'s C{direct} and C{inverse} methods
       for coincident points or lack of convergence.
    '''
    pass


class Cartesian(CartesianEllipsoidalBase):
    '''Extended to convert geocentric, L{Cartesian} points to
       Vincenty-based, ellipsoidal, geodetic L{LatLon}.
    '''
    _Ecef = None  # preferred C{EcefVeness} class

    @property_RO
    def Ecef(self):
        '''Get the ECEF I{class} (L{EcefVeness}).
        '''
        if Cartesian._Ecef is None:
            from pygeodesy.ecef import EcefVeness
            Cartesian._Ecef = EcefVeness
        return Cartesian._Ecef

    def toLatLon(self, **LatLon_datum_kwds):  # PYCHOK LatLon=LatLon, datum=None
        '''Convert this cartesian point to a C{Vincenty}-based
           geodetic point.

           @kwarg LatLon_datum_kwds: Optional L{LatLon}, B{C{datum}} and
                                     other keyword arguments, ignored if
                                     C{B{LatLon}=None}.  Use
                                     B{C{LatLon=...}} to override this
                                     L{LatLon} class or specify
                                     C{B{LatLon}=None}.

           @return: The geodetic point (L{LatLon}) or if B{C{LatLon}} is
                    C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height,
                    C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{LatLon}}, B{C{datum}} or other
                             B{C{LatLon_datum_kwds}}.
        '''
        kwds = _xkwds(LatLon_datum_kwds, LatLon=LatLon, datum=self.datum)
        return CartesianEllipsoidalBase.toLatLon(self, **kwds)


class LatLon(LatLonEllipsoidalBase):
    '''Using the formulae devised by U{I{Thaddeus Vincenty (1975)}
       <https://WikiPedia.org/wiki/Vincenty's_formulae>} for an (oblate)
       ellipsoidal model of the earth to compute the geodesic distance
       and bearings between two given points or the destination point
       given an start point and (initial) bearing.

       Set the earth model to be used with the keyword argument
       datum.  The default is Datums.WGS84, which is the most globally
       accurate.  For other models, see the Datums in module datum.

       Note: This implementation of I{Vincenty} methods may not converge
       for some valid points, raising a L{VincentyError}.  In that case,
       a result may be obtained by increasing the tolerance C{epsilon}
       and/or iteration C{limit}, see properties L{LatLon.epsilon} and
       L{LatLon.iterations}.
    '''
    _Ecef       = None     # preferred C{EcefVeness} class
    _epsilon    = 1.0e-12  # about 0.006 mm
    _iteration  = 0        # iteration number
    _iterations = 100      # vs Veness' 500

    def bearingTo(self, other, wrap=False):  # PYCHOK no cover
        '''DEPRECATED, use method C{initialBearingTo}.
        '''
        return self.initialBearingTo(other, wrap=wrap)

    def bearingTo2(self, other, wrap=False):
        '''Compute the initial and final bearing (forward and reverse
           azimuth) from this to an other point, using I{Vincenty}'s
           C{inverse} method.  See methods L{initialBearingTo} and
           L{finalBearingTo} for more details.

           @arg other: The other point (L{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: A L{Bearing2Tuple}C{(initial, final)}.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.

           @raise VincentyError: Vincenty fails to converge for the current
                                 L{LatLon.epsilon} and L{LatLon.iterations}
                                 limit and/or if this and the B{C{other}}
                                 point are coincident or near-antipodal.
        '''
        r = self._inverse(other, True, wrap)
        return self._xnamed(Bearing2Tuple(r.initial, r.final))

    def destination(self, distance, bearing, height=None):
        '''Compute the destination point after having travelled
           for the given distance from this point along a geodesic
           given by an initial bearing, using I{Vincenty}'s C{direct}
           method.  See method L{destination2} for more details.

           @arg distance: Distance (C{meter}).
           @arg bearing: Initial bearing (compass C{degrees360}).
           @kwarg height: Optional height, overriding the default
                          height (C{meter}, same units as C{distance}).

           @return: The destination point (L{LatLon}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise VincentyError: Vincenty fails to converge for the current
                                 L{LatLon.epsilon} and L{LatLon.iterations}
                                 limit.

           @example:

           >>> p = LatLon(-37.95103, 144.42487)
           >>> d = p.destination(54972.271, 306.86816)  # 37.6528°S, 143.9265°E
        '''
        r = self._direct(distance, bearing, True, height=height)
        return self._xnamed(r.destination)

    def destination2(self, distance, bearing, height=None):
        '''Compute the destination point and the final bearing (reverse
           azimuth) after having travelled for the given distance from
           this point along a geodesic given by an initial bearing,
           using I{Vincenty}'s C{direct} method.

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
                          height (C{meter}, same units as B{C{distance}}).

           @return: A L{Destination2Tuple}C{(destination, final)}.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise VincentyError: Vincenty fails to converge for the current
                                 L{LatLon.epsilon} and L{LatLon.iterations}
                                 limit.

           @example:

           >>> p = LatLon(-37.95103, 144.42487)
           >>> b = 306.86816
           >>> d, f = p.destination2(54972.271, b)
           >>> d
           LatLon(37°39′10.14″S, 143°55′35.39″E)  # 37.652818°S, 143.926498°E
           >>> f
           307.1736313846706
        '''
        r = self._direct(distance, bearing, True, height=height)
        return self._xnamed(r)

    def distanceTo(self, other, wrap=False, **unused):  # ignore radius=R_M
        '''Compute the distance between this and an other point
           along a geodesic, using I{Vincenty}'s C{inverse} method.
           See method L{distanceTo3} for more details.

           @arg other: The other point (L{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Distance (C{meter}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.

           @raise VincentyError: Vincenty fails to converge for the current
                                 L{LatLon.epsilon} and L{LatLon.iterations}
                                 limit and/or if this and the B{C{other}}
                                 point are coincident or near-antipodal.

           @example:

           >>> p = LatLon(50.06632, -5.71475)
           >>> q = LatLon(58.64402, -3.07009)
           >>> d = p.distanceTo(q)  # 969,954.166 m
        '''
        return self._inverse(other, False, wrap).distance

    def distanceTo3(self, other, wrap=False):
        '''Compute the distance, the initial and final bearing along a
           geodesic between this and an other point, using Vincenty's
           C{inverse} method.

           The distance is in the same units as this point's datum axes,
           conventially meter.  The distance is measured on the surface
           of the ellipsoid, ignoring this point's height.

           The initial and final bearing (forward and reverse azimuth)
           are in compass C{degrees360} from North.

           @arg other: Destination point (L{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: A L{Distance3Tuple}C{(distance, initial, final)}.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.

           @raise VincentyError: Vincenty fails to converge for the current
                                 L{LatLon.epsilon} and L{LatLon.iterations}
                                 limit and/or if this and the B{C{other}}
                                 point are coincident or near-antipodal.
        '''
        return self._xnamed(self._inverse(other, True, wrap))

    @property_RO
    def Ecef(self):
        '''Get the ECEF I{class} (L{EcefVeness}).
        '''
        if LatLon._Ecef is None:
            from pygeodesy.ecef import EcefVeness
            LatLon._Ecef = EcefVeness
        return LatLon._Ecef

    @property_doc_(''' the convergence epsilon (C{scalar}).''')
    def epsilon(self):
        '''Get the convergence epsilon (C{scalar}).
        '''
        return self._epsilon

    @epsilon.setter  # PYCHOK setter!
    def epsilon(self, eps):
        '''Set the convergence epsilon.

           @arg eps: New epsilon (C{scalar}).

           @raise TypeError: Non-scalar B{C{eps}}.

           @raise ValueError: Out of bounds B{C{eps}}.
        '''
        self._epsilon = Scalar_(epsilon=eps)

    def finalBearingOn(self, distance, bearing):
        '''Compute the final bearing (reverse azimuth) after having
           travelled for the given distance along a geodesic given
           by an initial bearing from this point, using Vincenty's
           C{direct} method.  See method L{destination2} for more details.

           @arg distance: Distance (C{meter}).
           @arg bearing: Initial bearing (compass C{degrees360}).

           @return: Final bearing (compass C{degrees360}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.

           @raise VincentyError: Vincenty fails to converge for the current
                                 L{LatLon.epsilon} and L{LatLon.iterations}
                                 limit.

           @example:

           >>> p = LatLon(-37.95103, 144.42487)
           >>> b = 306.86816
           >>> f = p.finalBearingOn(54972.271, b)  # 307.1736
        '''
        return self._direct(distance, bearing, False).final

    def finalBearingTo(self, other, wrap=False):
        '''Compute the final bearing (reverse azimuth) after having
           travelled along a geodesic from this point to an other
           point, using I{Vincenty}'s C{inverse} method.  See method
           L{distanceTo3} for more details.

           @arg other: The other point (L{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Final bearing (compass C{degrees360}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.

           @raise VincentyError: Vincenty fails to converge for the current
                                 L{LatLon.epsilon} and L{LatLon.iterations}
                                 limit and/or if this and the B{C{other}}
                                 point are coincident or near-antipodal.

           @example:

           >>> p = new LatLon(50.06632, -5.71475)
           >>> q = new LatLon(58.64402, -3.07009)
           >>> f = p.finalBearingTo(q)  # 11.2972°

           >>> p = LatLon(52.205, 0.119)
           >>> q = LatLon(48.857, 2.351)
           >>> f = p.finalBearingTo(q)  # 157.9
        '''
        return self._inverse(other, True, wrap).final

    def initialBearingTo(self, other, wrap=False):
        '''Compute the initial bearing (forward azimuth) to travel
           along a geodesic from this point to an other point,
           using I{Vincenty}'s C{inverse} method.  See method
           L{distanceTo3} for more details.

           @arg other: The other point (L{LatLon}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Initial bearing (compass C{degrees360}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.

           @raise VincentyError: Vincenty fails to converge for the current
                                 L{LatLon.epsilon} and L{LatLon.iterations}
                                 limit and/or if this and the B{C{other}}
                                 point are coincident or near-antipodal.

           @example:

           >>> p = LatLon(50.06632, -5.71475)
           >>> q = LatLon(58.64402, -3.07009)
           >>> b = p.initialBearingTo(q)  # 9.141877°

           >>> p = LatLon(52.205, 0.119)
           >>> q = LatLon(48.857, 2.351)
           >>> b = p.initialBearingTo(q)  # 156.11064°

           @JSname: I{bearingTo}.
        '''
        return self._inverse(other, True, wrap).initial

    def intermediateTo(self, other, fraction, height=None, wrap=False):
        '''Return the point at given fraction along the geodesic between
           this and an other point, using I{Vincenty}'s C{direct} and
           C{inverse} methods.

           @arg other: The other point (L{LatLon}).
           @arg fraction: Fraction between both points ranging from
                          0, meaning this to 1, the other point (C{float}).
           @kwarg height: Optional height, overriding the fractional
                          height (C{meter}).
           @kwarg wrap: Wrap and unroll longitudes (C{bool}).

           @return: Intermediate point (L{LatLon}).

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise UnitError: Invalid B{C{fraction}} or B{C{height}}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.

           @raise VincentyError: Vincenty fails to converge for the current
                                 L{LatLon.epsilon} and L{LatLon.iterations}
                                 limit and/or if this and the B{C{other}}
                                 point are coincident or near-antipodal.

           @see: Methods L{distanceTo3} and L{destination}.

           @example:

           >>> p = ellipsoidalVincenty.LatLon(52.205, 0.119)
           >>> q = ellipsoidalVincenty.LatLon(48.857, 2.351)
           >>> i = p.intermediateTo(q, 0.25)  # 51.372275°N, 000.707253°E
        '''
        return _intermediateTo(self, other, fraction, height, wrap)

    @property_doc_(''' the iteration limit (C{int}).''')
    def iterations(self):
        '''Get the iteration limit (C{int}).
        '''
        return self._iterations

    @iterations.setter  # PYCHOK setter!
    def iterations(self, limit):
        '''Set the iteration limit.

           @arg limit: New iteration limit (C{int}).

           @raise TypeError: Non-scalar B{C{limit}}.

           @raise ValueError: Out-of-bounds B{C{limit}}.
        '''
        self._iterations = Number_(limit, name=_limit_, low=4, high=500)

    def toCartesian(self, **Cartesian_datum_kwds):  # PYCHOK Cartesian=Cartesian, datum=None
        '''Convert this point to C{Vincenty}-based cartesian (ECEF)
           coordinates.

           @kwarg Cartesian_datum_kwds: Optional L{Cartesian}, B{C{datum}}
                                        and other keyword arguments, ignored
                                        if B{C{Cartesian=None}}.  Use
                                        B{C{Cartesian=...}} to override this
                                        L{Cartesian} class or specify
                                        B{C{Cartesian=None}}.

           @return: The cartesian point (L{Cartesian}) or if B{C{Cartesian}}
                    is C{None}, an L{Ecef9Tuple}C{(x, y, z, lat, lon, height,
                    C, M, datum)} with C{C} and C{M} if available.

           @raise TypeError: Invalid B{C{Cartesian}}, B{C{datum}} or other
                             B{C{Cartesian_datum_kwds}}.
        '''
        kwds = _xkwds(Cartesian_datum_kwds, Cartesian=Cartesian,
                                                datum=self.datum)
        return LatLonEllipsoidalBase.toCartesian(self, **kwds)

    def _direct(self, distance, bearing, llr, height=None):
        '''(INTERNAL) Direct Vincenty method.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.

           @raise VincentyError: Vincenty fails to converge for the current
                                 L{LatLon.epsilon} and L{LatLon.iterations}
                                 limit.
        '''
        E = self.ellipsoid()

        c1, s1, t1 = _r3(self.lat, E.f)

        i = radians(bearing)  # initial bearing (forward azimuth)
        si, ci = sincos2(i)
        s12 = atan2(t1, ci) * _2_0

        sa = c1 * si
        c2a = _1_0 - sa**2
        if c2a < EPS:
            c2a = _0_0
            A, B = _1_0, _0_0
        else:  # e22 == (a / b)**2 - 1
            A, B = _p2(c2a * E.e22)

        s = d = distance / (E.b * A)
        for self._iteration in range(1, self._iterations + 1):
            ss, cs = sincos2(s)
            c2sm = cos(s12 + s)
            s_, s = s, d + _ds(B, cs, ss, c2sm)
            if abs(s - s_) < self._epsilon:
                break
        else:
            raise VincentyError(_no_(_convergence_), txt=repr(self))  # self.toRepr()

        t = s1 * ss - c1 * cs * ci
        # final bearing (reverse azimuth +/- 180)
        r = atan2b(sa, -t)

        if llr:
            # destination latitude in [-90, 90)
            a = degrees90(atan2(s1 * cs + c1 * ss * ci, E.b_a * hypot(sa, t)))
            # destination longitude in [-180, 180)
            t = radians(self.lon) - _dl(E.f, c2a, sa, s, cs, ss, c2sm)
            b = degrees180(atan2(ss * si, c1 * cs - s1 * ss * ci) + t)
            h = self.height if height is None else height
            d = self.classof(a, b, height=h, datum=self.datum)
        else:
            d = None
        return Destination2Tuple(d, r)

    def _inverse(self, other, azis, wrap):
        '''(INTERNAL) Inverse Vincenty method.

           @raise TypeError: The B{C{other}} point is not L{LatLon}.

           @raise ValueError: If this and the B{C{other}} point's L{Datum}
                              ellipsoids are not compatible.

           @raise VincentyError: Vincenty fails to converge for the current
                                 L{LatLon.epsilon} and L{LatLon.iterations}
                                 limit and/or if this and the B{C{other}}
                                 point are coincident or near-antipodal.
        '''
        E = self.ellipsoids(other)

        c1, s1, _ = _r3(self.lat, E.f)
        c2, s2, _ = _r3(other.lat, E.f)

        c1c2, s1c2 = c1 * c2, s1 * c2
        c1s2, s1s2 = c1 * s2, s1 * s2

        dl, _ = unroll180(self.lon, other.lon, wrap=wrap)
        ll = dl = radians(dl)
        for self._iteration in range(1, self._iterations + 1):
            ll_ = ll
            sll, cll = sincos2(ll)

            ss = hypot(c2 * sll, c1s2 - s1c2 * cll)
            if ss < EPS:  # coincident or antipodal, ...
                if self.isantipodeTo(other, eps=self._epsilon):
                    t = '%r %s%s %r' % (self, _antipodal_, _to_, other)
                    raise VincentyError(_ambiguous_, txt=t)
                # return zeros like Karney, but unlike Veness
                return Distance3Tuple(_0_0, 0, 0)

            cs = s1s2 + c1c2 * cll
            s = atan2(ss, cs)

            sa = c1c2 * sll / ss
            c2a = _1_0 - sa**2
            if abs(c2a) < EPS0:
                c2a = _0_0  # equatorial line
                ll = dl + E.f * sa * s
            else:
                c2sm = cs - 2 * s1s2 / c2a
                ll = dl + _dl(E.f, c2a, sa, s, cs, ss, c2sm)

            if abs(ll - ll_) < self._epsilon:
                break
#           # omitted and applied only after failure to converge below, see footnote
#           # under Inverse at <https://WikiPedia.org/wiki/Vincenty's_formulae>
#           # <https://GitHub.com/ChrisVeness/geodesy/blob/master/latlon-vincenty.js>
#           elif abs(ll) > PI and self.isantipodeTo(other, eps=self._epsilon):
#              raise VincentyError('%s, %r %sto %r' % ('ambiguous', self,
#                                  _antipodal_, other))
        else:
            t = _antipodal_ if self.isantipodeTo(other, eps=self._epsilon) else NN
            t = '%r %s%s %r' % (self, t, _to_, other)
            raise VincentyError(_no_(_convergence_), txt=t)

        if c2a:  # e22 == (a / b)**2 - 1
            A, B = _p2(c2a * E.e22)
            s = A * (s - _ds(B, cs, ss, c2sm))

        b = E.b
#       if self.height or other.height:
#           b += self._havg(other)
        d = b * s

        if azis:  # forward and reverse azimuth
            sll, cll = sincos2(ll)
            f = atan2b(c2 * sll,  c1s2 - s1c2 * cll)
            r = atan2b(c1 * sll, -s1c2 + c1s2 * cll)
        else:
            f = r = _0_0
        return Distance3Tuple(d, f, r)


def _c2sm2(c2sm):
    '''(INTERNAL) 2 * c2sm**2 - 1.
    '''
    return _2_0 * c2sm**2 - _1_0


def _dl(f, c2a, sa, s, cs, ss, c2sm):
    '''(INTERNAL) Dl.
    '''
    C = f / _16_0 * c2a * (_4_0 + f * (_4_0 - _3_0 * c2a))
    return (_1_0 - C) * f * sa * (s + C * ss * (c2sm +
                                      C * cs * _c2sm2(c2sm)))


def _ds(B, cs, ss, c2sm):
    '''(INTERNAL) Ds.
    '''
    if B:
        c2sm2 = _c2sm2(c2sm)
        ss2 = (_4_0 * ss**2 - _3_0) * (_2_0 * c2sm2 - _1_0)
        B *= ss * (c2sm + B / _4_0 * (c2sm2 * cs -
                          B / _6_0 *  c2sm  * ss2))
    return B


def _p2(u2):  # e'2 WGS84 = 0.00673949674227643
    '''(INTERNAL) Compute A, B polynomials.
    '''
    A = fpolynomial(u2, 16384, 4096, -768, 320, -175) / 16384
    B = fpolynomial(u2,     0,  256, -128,  74,  -47) / 1024
    return A, B


def _r3(a, f):
    '''(INTERNAL) Reduced cos, sin, tan.
    '''
    t = (_1_0 - f) * tan(radians(a))
    c =  _1_0 / hypot1(t)
    s =  t * c
    return c, s, t


def areaOf(points, datum=Datums.WGS84, wrap=True):  # PYCHOK no cover
    '''DEPRECATED, use function C{ellipsoidalKarney.areaOf}.
    '''
    from pygeodesy.ellipsoidalKarney import areaOf
    return areaOf(points, datum=datum, wrap=wrap)


def _Equidistant(equidistant):
    # (INTERNAL) Get the C{Equidistant} class.
    if equidistant is None or not callable(equidistant):
        from pygeodesy.azimuthal import Equidistant as equidistant
    return equidistant


def intersections2(center1, radius1, center2, radius2, height=None, wrap=True,
                   equidistant=None, tol=_TOL_M, LatLon=LatLon, **LatLon_kwds):
    '''Iteratively compute the intersection points of two circles each defined
       by an (ellipsoidal) center point and a radius.

       @arg center1: Center of the first circle (L{LatLon}).
       @arg radius1: Radius of the first circle (C{meter}, conventionally).
       @arg center2: Center of the second circle (L{LatLon}).
       @arg radius2: Radius of the second circle (C{meter}, same units as
                     B{C{radius2}}).
       @kwarg height: Optional height for the intersection points,
                      overriding the "radical height" at the "radical
                      line" between both centers (C{meter}) or C{None}.
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).
       @kwarg equidistant: An azimuthal equidistant projection (class
                           L{EquidistantKarney} or function L{equidistant}),
                           or C{None} for L{Equidistant}.
       @kwarg tol: Convergence tolerance (C{meter}, same units as B{C{radius1}}
                   and B{C{radius2}}).
       @kwarg LatLon: Optional class to return the intersection points
                      (L{LatLon}) or C{None}.
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                           arguments, ignored if C{B{LatLon}=None}.

       @return: 2-Tuple of the intersection points, each a B{C{LatLon}}
                instance or L{LatLon4Tuple}C{(lat, lon, height, datum)}
                if B{C{LatLon}} is C{None}.  For abutting circles, both
                intersection points are the same instance.

       @raise IntersectionError: Concentric, antipodal, invalid or
                                 non-intersecting circles or no
                                 convergence for the B{C{tol}}.

       @raise ImportError: Package U{geographiclib
                           <https://PyPI.org/project/geographiclib>}
                           not installed or not found.

       @raise TypeError: Invalid or non-ellipsoidal B{C{center1}} or B{C{center2}}
                         or invalid B{C{equidistant}}.

       @raise UnitError: Invalid B{C{radius1}}, B{C{radius2}} or B{C{height}}.

       @see: U{The B{ellipsoidal} case<https://GIS.StackExchange.com/questions/48937/
             calculating-intersection-of-two-circles>}, U{Karney's paper
             <https://ArXiv.org/pdf/1102.1215.pdf>}, pp 20-21, section B{I{14. MARITIME BOUNDARIES}},
             U{circle-circle<https://MathWorld.Wolfram.com/Circle-CircleIntersection.html>} and
             U{sphere-sphere<https://MathWorld.Wolfram.com/Sphere-SphereIntersection.html>}
             intersections.
    '''
    E = _Equidistant(equidistant)
    return _intersections2(center1, radius1, center2, radius2, height=height, wrap=wrap,
                                    equidistant=E, tol=tol, LatLon=LatLon, **LatLon_kwds)


def nearestOn(point, point1, point2, within=True, height=None, wrap=False,
              equidistant=None, tol=_TOL_M, LatLon=LatLon, **LatLon_kwds):
    '''Iteratively locate the closest point on the arc between two other points.

       @arg point: Reference point (C{LatLon}).
       @arg point1: Start point of the arc (C{LatLon}).
       @arg point2: End point of the arc (C{LatLon}).
       @kwarg within: If C{True} return the closest point I{between}
                      B{C{point1}} and B{C{point2}}, otherwise the
                      closest point elsewhere on the arc (C{bool}).
       @kwarg height: Optional height for the closest point (C{meter})
                      or C{None}.
       @kwarg wrap: Wrap and unroll longitudes (C{bool}).
       @kwarg equidistant: An azimuthal equidistant projection (class
                           L{EquidistantKarney} or function L{equidistant}),
                           or C{None} for L{Equidistant}.
       @kwarg tol: Convergence tolerance (C{meter}).
       @kwarg LatLon: Optional class to return the closest point
                      (L{LatLon}) or C{None}.
       @kwarg LatLon_kwds: Optional, additional B{C{LatLon}} keyword
                           arguments, ignored if C{B{LatLon}=None}.

       @return: Closest point (B{C{LatLon}}).

       @raise ImportError: Package U{geographiclib
                           <https://PyPI.org/project/geographiclib>}
                           not installed or not found.

       @raise TypeError: Invalid or non-ellipsoidal B{C{point}}, B{C{point1}}
                         or B{C{point2}} or invalid B{C{equidistant}}.
    '''
    p = _xellipsoidal(point=point)
    p1 = p.others(point1=point1)
    p2 = p.others(point2=point2)
    E = _Equidistant(equidistant)
    return _nearestOn(p, p1, p2, within=within, height=height, wrap=wrap,
                      equidistant=E, tol=tol, LatLon=LatLon, **LatLon_kwds)


def perimeterOf(points, closed=False, datum=Datums.WGS84, wrap=True):  # PYCHOK no cover
    '''DEPRECATED, use function C{ellipsoidalKarney.perimeterOf}.

       @raise ImportError: Package U{geographiclib
                           <https://PyPI.org/project/geographiclib>}
                           not installed or not found.
    '''
    from pygeodesy.ellipsoidalKarney import perimeterOf
    return perimeterOf(points, closed=closed, datum=datum, wrap=wrap)


__all__ += _ALL_OTHER(Cartesian, LatLon,
                      intersections2, ispolar,  # from .points
                      nearestOn) + _ALL_DOCS(perimeterOf)

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
