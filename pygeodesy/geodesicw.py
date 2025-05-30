
# -*- coding: utf-8 -*-

u'''Wrapper around Python classes C{geodesic.Geodesic} and C{geodesicline.GeodesicLine} from
I{Karney}'s Python package U{geographiclib<https://PyPI.org/project/geographiclib>}, provided
that package is installed.

The I{wrapped} class methods return a L{GDict} instance offering access to the C{dict} items
either by C{key} or by C{attribute} name.

With env variable C{PYGEODESY_GEOGRAPHICLIB} left undefined or set to C{"2"}, this module and modules
L{pygeodesy.geodesici}, L{pygeodesy.geodesicx} and L{pygeodesy.karney} will use U{GeographicLib 2.0
<https://GeographicLib.SourceForge.io/C++/doc/>} transcoding, otherwise C{1.52} or older.
'''

from pygeodesy.basics import _copysign, _xinstanceof
from pygeodesy.constants import EPS, NAN, _EPSqrt as _TOL, _0_5
from pygeodesy.datums import _earth_datum, _WGS84,  _EWGS84
# from pygeodesy.dms import F_D  # from .latlonBase
# from pygeodesy.ellipsoids import _EWGS84  # from .datums
from pygeodesy.errors import _AssertionError, GeodesicError, \
                              IntersectionError
from pygeodesy.fsums import Fsum,  Fmt, unstr
from pygeodesy.internals import typename, _under
from pygeodesy.interns import NN, _DOT_, _SPACE_, _to_, _too_
from pygeodesy.karney import _atan2d, Caps, Direct9Tuple, GDict, \
                              Inverse10Tuple, _kWrapped
from pygeodesy.latlonBase import LatLonBase as _LLB,  F_D, Radius_
from pygeodesy.lazily import _ALL_LAZY, _ALL_MODS as _MODS
from pygeodesy.named import callername, classname, _name1__, _name2__
from pygeodesy.namedTuples import Destination3Tuple, Distance3Tuple
from pygeodesy.props import Property, Property_RO, property_RO, \
                            property_ROver
# from pygeodesy.streprs import Fmt, unstr  # from .fsums
# from pygeodesy.units import Radius_  # from .latlonBase
from pygeodesy.utily import _unrollon, _Wrap, wrap360,  fabs  # PYCHOK used!

from contextlib import contextmanager
# from math import fabs  # from .utily

__all__ = _ALL_LAZY.geodesicw
__version__ = '25.05.28'

_plumb_ = 'plumb'
_TRIPS  =  65


class _gWrapped(_kWrapped):
    '''(INTERNAL) Wrapper for some of I{Karney}'s U{geographiclib
        <https://PyPI.org/project/geographiclib>} classes.
    '''

    @property_ROver  # MCCABE 24
    def Geodesic(self):
        '''Get the I{wrapped} C{geodesic.Geodesic} class from I{Karney}'s Python
           U{geographiclib<https://GitHub.com/geographiclib/geographiclib-python>},
           provided the latter is installed.
        '''
        _Geodesic = self.geographiclib.Geodesic
        if not (Caps.LATITUDE      == _Geodesic.LATITUDE and
                Caps.LONGITUDE     == _Geodesic.LONGITUDE and
                Caps.AZIMUTH       == _Geodesic.AZIMUTH and
                Caps.DISTANCE      == _Geodesic.DISTANCE and
                Caps.DISTANCE_IN   == _Geodesic.DISTANCE_IN and
                Caps.REDUCEDLENGTH == _Geodesic.REDUCEDLENGTH and
                Caps.GEODESICSCALE == _Geodesic.GEODESICSCALE and
                Caps.AREA          == _Geodesic.AREA and
                Caps.STANDARD      == _Geodesic.STANDARD and
                Caps.ALL           == _Geodesic.ALL):
            raise _AssertionError(Caps=bin(Caps.ALL),
                                  Geodesic=bin(_Geodesic.ALL))

        class Geodesic(_Geodesic):
            '''I{Wrapper} for I{Karney}'s Python U{geodesic.Geodesic
               <https://PyPI.org/project/geographiclib>} class.
            '''
            _datum        = _WGS84
            _debug        =  0  # like .geodesicx.bases._GeodesicBase
            LINE_OFF      =  0  # in .azimuthal._GnomonicBase and .css.CassiniSoldner
            _name         =  NN
            STANDARD_LINE =  Caps.STANDARD_LINE

            def __init__(self, a_ellipsoid=_EWGS84, f=None, **name):  # PYCHOK signature
                '''New I{wrapped} C{geodesic.Geodesic} instance.

                   @arg a_ellipsoid: The equatorial radius I{a} (C{meter}, conventionally),
                                     an ellipsoid (L{Ellipsoid}) or a datum (L{Datum}).
                   @arg f: The ellipsoid's flattening (C{scalar}), required if B{C{a_ellipsoid})
                           is C{meter}, ignored otherwise.
                   @kwarg name: Optional C{B{name}=NN} (C{str}).
                '''
                _earth_datum(self, a_ellipsoid, f=f, **name)  # raiser=NN
                E = self.ellipsoid
                with _wargs(self, *E.a_f, **name) as args:
                    _Geodesic.__init__(self, *args)
                if name:
                    self._name, _ = _name2__(name, _or_nameof=E)

            def Area(self, polyline=False, **name):  # like GeodesicExact.Area
                '''Return a C{PolygonArea} instance with method C{Compute} extended.
                '''
                _AreaBase = _wrapped._PolygonArea  # in .karney._kwrapped

                class _PolygonArea(_AreaBase):
                    # def __init__(self, *earth_polyline):
                    #     _PolygonArea.__init__(self, *earth_polyline)

                    def Compute(self, reverse=False, sign=True, polar=False):
                        '''Use C{B{polar}=True} to adjust the area, see function
                           L{areaOf<pygeodesy.geodesicx.gxarea.areaOf>}.
                        '''
                        n, p, a = _AreaBase.Compute(self, reverse=reverse, sign=sign)
                        if polar:  # see .geodesicx.gxarea.GeodesicAreaExact._reduced
                            a += _copysign(self.area0 * _0_5 * n, a)
                        return n, p, a

                A = _PolygonArea(self, polyline)
                A.name = _name2__(name, _or_nameof=self)
                return A

            def ArcDirect(self, lat1, lon1, azi1, a12, outmask=Caps.STANDARD):  # PYCHOK no cover
                '''Return the C{_Geodesic.ArcDirect} result as L{GDict}.
                '''
                with _wargs(self, lat1, lon1, azi1, a12, outmask) as args:
                    d = _Geodesic.ArcDirect(self, *args)
                return GDict(d)

            def ArcDirectLine(self, lat1, lon1, azi1, a12, caps=Caps.STANDARD_LINE, **name):  # PYCHOK no cover
                '''Return the C{_Geodesic.ArcDirectLine} as I{wrapped} C{GeodesicLine}.
                '''
                return self._GenDirectLine(lat1, lon1, azi1, True, a12, caps, **name)

            @property_RO
            def datum(self):
                '''Get this geodesic's datum (C{Datum}).
                '''
                return self._datum

            @Property
            def debug(self):
                '''Get the C{debug} option (C{bool}).
                '''
                return bool(self._debug)

            @debug.setter  # PYCHOK setter!
            def debug(self, debug):
                '''Set the C{debug} option (C{bool}) to include more
                   details in L{GDict} results.
                '''
                self._debug = Caps._DEBUG_ALL if debug else 0

            def Direct(self, lat1, lon1, azi1, s12=0, outmask=Caps.STANDARD):
                '''Return the C{_Geodesic.Direct} result as L{GDict}.
                '''
                with _wargs(self, lat1, lon1, azi1, s12, outmask) as args:
                    d = _Geodesic.Direct(self, *args)
                return GDict(d)

            def Direct3(self, lat1, lon1, azi1, s12):  # PYCHOK outmask
                '''Return the destination lat, lon and reverse azimuth
                   in C{degrees} as L{Destination3Tuple}.
                '''
                d = self.Direct(lat1, lon1, azi1, s12, outmask=Caps._DIRECT3)
                return Destination3Tuple(d.lat2, d.lon2, d.azi2)

            def _DirectLine(self, ll1, azi12, s12=0, **caps_name):
                '''(INTERNAL) Short-cut version.
                '''
                return self.DirectLine(ll1.lat, ll1.lon, azi12, s12, **caps_name)

            def DirectLine(self, lat1, lon1, azi1, s12, caps=Caps.STANDARD_LINE, **name):
                '''Return the C{_Geodesic.DirectLine} as I{wrapped} C{GeodesicLine}.
                '''
                return self._GenDirectLine(lat1, lon1, azi1, False, s12, caps, **name)

            @Property_RO
            def ellipsoid(self):
                '''Get this geodesic's ellipsoid (C{Ellipsoid}).
                '''
                return self.datum.ellipsoid

            @property_RO
            def f1(self):  # in .css.CassiniSoldner.reset
                '''Get the geodesic's ellipsoid's I{1 - flattening} (C{float}).
                '''
                return getattr(self, _under(Geodesic.f1.name), self.ellipsoid.f1)

            def _GDictDirect(self, lat, lon, azi, arcmode, s12_a12, outmask=Caps.STANDARD):
                '''(INTERNAL) Get C{_Geodesic._GenDirect} result as C{GDict}.
                '''
                with _wargs(self, lat, lon, azi, arcmode, s12_a12, outmask) as args:
                    t = _Geodesic._GenDirect(self, *args)
                return Direct9Tuple(t).toGDict()  # *t

            def _GDictInverse(self, lat1, lon1, lat2, lon2, outmask=Caps.STANDARD):
                '''(INTERNAL) Get C{_Geodesic._GenInverse} result as L{Inverse10Tuple}.
                '''
                with _wargs(self, lat1, lon1, lat2, lon2, outmask) as args:
                    t = _Geodesic._GenInverse(self, *args)
                return Inverse10Tuple(t).toGDict(lon1=lon1, lon2=lon2)  # *t

            def _GenDirectLine(self, lat1, lon1, azi1, arcmode, s12_a12, *caps, **name):
                '''(INTERNAL) Invoked by C{_Geodesic.DirectLine} and C{-.ArcDirectLine},
                   returning the result as a I{wrapped} C{GeodesicLine}.
                '''
                with _wargs(self, lat1, lon1, azi1, arcmode, s12_a12, *caps, **name) as args:
                    t = _Geodesic._GenDirectLine(self, *args)
                return self._Line13(t, **name)

            def _Inverse(self, ll1, ll2, wrap, **outmask):
                '''(INTERNAL) Short-cut version, see .ellipsoidalBaseDI.intersecant2.
                '''
                if wrap:
                    ll2 = _unrollon(ll1, _Wrap.point(ll2))
                return self.Inverse(ll1.lat, ll1.lon, ll2.lat, ll2.lon, **outmask)

            def Inverse(self, lat1, lon1, lat2, lon2, outmask=Caps.STANDARD):
                '''Return the C{_Geodesic.Inverse} result as L{GDict}.
                '''
                with _wargs(self, lat1, lon1, lat2, lon2, outmask) as args:
                    d = _Geodesic.Inverse(self, *args)
                return GDict(d)

            def Inverse1(self, lat1, lon1, lat2, lon2, wrap=False):
                '''Return the non-negative, I{angular} distance in C{degrees}.

                   @kwarg wrap: If C{True}, wrap or I{normalize} and unroll
                                B{C{lat2}} and BC{lon2}} (C{bool}).
                '''
                # see .FrechetKarney.distance, .HausdorffKarney._distance
                # and .HeightIDWkarney._distances
                if wrap:
                    _, lat2, lon2 = _Wrap.latlon3(lat1, lat2, lon2, True)  # _Geodesic.LONG_UNROLL
                r = self.Inverse(lat1, lon1, lat2, lon2)
                # XXX _Geodesic.DISTANCE needed for 'a12'?
                return fabs(r.a12)

            def Inverse3(self, lat1, lon1, lat2, lon2):  # PYCHOK outmask
                '''Return the distance in C{meter} and the forward and reverse
                   azimuths in C{degrees} as L{Distance3Tuple}.
                '''
                r = self.Inverse(lat1, lon1, lat2, lon2, outmask=Caps._INVERSE3)
                return Distance3Tuple(r.s12, wrap360(r.azi1), wrap360(r.azi2))

            def _InverseLine(self, ll1, ll2, wrap, **caps_name):
                '''(INTERNAL) Short-cut version.
                '''
                if wrap:
                    ll2 = _unrollon(ll1, _Wrap.point(ll2))
                return self.InverseLine(ll1.lat, ll1.lon, ll2.lat, ll2.lon, **caps_name)

            def InverseLine(self, lat1, lon1, lat2, lon2, caps=Caps.STANDARD_LINE, **name):
                '''Return the C{_Geodesic.InverseLine} as I{wrapped} C{GeodesicLine}.
                '''
                with _wargs(self, lat1, lon1, lat2, lon2, caps, **name) as args:
                    t = _Geodesic.InverseLine(self, *args)
                return self._Line13(t, **name)

            def Line(self, lat1, lon1, azi1, caps=Caps.STANDARD_LINE, **name):
                '''Set up a I{wrapped} C{GeodesicLine} to compute several points
                   along a single, I{wrapped} (this) geodesic.
                '''
                return _wrapped.GeodesicLine(self, lat1, lon1, azi1, caps=caps, **name)

            def _Line13(self, t, **name):
                '''(INTERNAL) Wrap C{_GeodesicLine}, add distance and arc length
                   to reference point 3.
                '''
                gl = _wrapped.GeodesicLine(self, t.lat1, t.lon1, t.azi1, caps=t.caps,
                                                 salp1=t.salp1, calp1=t.calp1, **name)
                gl.a13, gl.s13 = t.a13, t.s13
                return gl

            @property_RO
            def name(self):
                '''Get the name (C{str}).
                '''
                return self._name

            Polygon = Area
            WGS84   = None  # _EWGS84.geodesicw recusion

        # Geodesic.ArcDirect.__doc__   = _Geodesic.ArcDirect.__doc__
        # Geodesic.Direct.__doc__      = _Geodesic.Direct.__doc__
        # Geodesic.Inverse.__doc__     = _Geodesic.Inverse.__doc__
        # Geodesic.InverseLine.__doc__ = _Geodesic.InverseLinr.__doc__
        # Geodesic.Line.__doc__        = _Geodesic.Line.__doc__
        return Geodesic  # overwrite property_ROver

    @property_ROver  # MCCABE 16
    def GeodesicLine(self):
        '''Get the I{wrapped} C{geodesicline.GeodesicLine} class from I{Karney}'s
           Python U{geographiclib<https://GitHub.com/geographiclib/geographiclib-python>},
           provided the latter is installed.
        '''
        _GeodesicLine = self.geographiclib.GeodesicLine

        class GeodesicLine(_GeodesicLine):
            '''I{Wrapper} for I{Karney}'s Python U{geodesicline.GeodesicLine
               <https://PyPI.org/project/geographiclib>} class.
            '''
            _geodesic = None
            _name     = NN

            def __init__(self, geodesic, lat1, lon1, azi1, **caps_name_):  # salp1=NAN, calp1=NAN
                '''New I{wrapped} C{geodesicline.GeodesicLine} instance.

                   @arg geodesic: A I{wrapped} C{Geodesic} instance.
                   @arg lat1: Latitude of the first points (C{degrees}).
                   @arg lon1: Longitude of the first points (C{degrees}).
                   @arg azi1: Azimuth at the first points (compass C{degrees360}).
                   @kwarg caps_name_: Optional keyword arguments C{B{caps}=Caps.STANDARD},
                               a bit-or'ed combination of L{Caps<pygeodesy.karney.Caps>}
                               values specifying the capabilities the C{GeodesicLine} instance
                               should possess, an optional C{B{name}=NN} plus C{salp1=NAN} and
                               C{calp1=NAN} for I{INTERNAL} use.
                '''
                _xinstanceof(_wrapped.Geodesic, geodesic=geodesic)
                with _wargs(self, geodesic, lat1, lon1, azi1, **caps_name_) as args:
                    name, caps_ = _name2__(caps_name_, _or_nameof=geodesic)
                    _GeodesicLine.__init__(self, *args, **caps_)  # XXX avoid updates?
                    if name:
                        self._name = name
                self._geodesic = geodesic

            @Property_RO
            def a1(self):
                '''Get the I{equatorial arc} (C{degrees}), the arc length between
                   the northward equatorial crossing and point C{(lat1, lon1)}.

                   @see: U{EquatorialArc<https://GeographicLib.SourceForge.io/
                         C++/doc/classGeographicLib_1_1GeodesicLine.html>}
                '''
                try:
                    return _atan2d(self._ssig1, self._csig1)
                except AttributeError:
                    return NAN  # see .geodesicx.gxline._GeodesicLineExact

            equatorarc = a1

            def Arc(self):
                '''Return the angular distance to point 3 (C{degrees} or C{NAN}).
                '''
                return self.a13

            def ArcPosition(self, a12, outmask=Caps.STANDARD):
                '''Return the position at C{B{a12} degrees} on this line.

                   @arg a12: Angular distance from this line's first point
                             (C{degrees}).

                   @see: Method L{Position} for further details.
                '''
                with _wargs(self, a12, outmask) as args:
                    d = _GeodesicLine.ArcPosition(self, *args)
                return GDict(d)

            @Property_RO
            def azi0(self):  # see .css.CassiniSoldner.forward4
                '''Get the I{equatorial azimuth} (C{degrees}), the azimuth of the
                   geodesic line as it crosses the equator in a northward direction.

                   @see: U{EquatorialAzimuth<https://GeographicLib.SourceForge.io/
                         C++/doc/classGeographicLib_1_1GeodesicLine.html>}
                '''
                try:
                    return _atan2d(self._salp0, self._calp0)
                except AttributeError:
                    return NAN  # see .geodesicx.gxline._GeodesicLineExact

            equatorazimuth = azi0

            def Distance(self):
                '''Return the distance to reference point 3 (C{meter} or C{NAN}).
                '''
                return self.s13

            @property_RO
            def geodesic(self):
                '''Get the I{wrapped} geodesic (L{Geodesic}).
                '''
                return self._geodesic

            def Intersecant2(self, lat0, lon0, radius, tol=_TOL):
                '''Compute the intersection(s) of this geodesic line and a circle.

                   @arg lat0: Latitude of the circle center (C{degrees}).
                   @arg lon0: Longitude of the circle center (C{degrees}).
                   @arg radius: Radius of the circle (C{meter}, conventionally).
                   @kwarg tol: Convergence tolerance (C{scalar}).

                   @return: 2-Tuple C{(P, Q)} with both intersections points (representing
                            a geodesic chord), each a L{GDict} from method L{Position} and
                            extended to 14 items C{lat1, lon1, azi1, lat2, lon2, azi2, a12,
                            s12, lat0, lon0, azi0, a02, s02, at} with the circle center
                            C{lat0}, C{lon0}, azimuth C{azi0} at the intersection, distance
                            C{a02} in C{degrees} and C{s02} in C{meter} along the geodesic
                            from the circle center to the intersection C{lat2, lon2} and
                            the angle C{at} between the geodesic and this line at the
                            intersection.  The I{geodesic} azimuth at the intersection is
                            C{(at + azi2)}.  If this line is tangential to the circle, both
                            intersections are the same L{GDict} instance.

                   @raise IntersectionError: The circle and this geodesic line do not
                                             intersect.

                   @raise UnitError: Invalid B{C{radius}}.
                '''
                return _Intersecant2(self, lat0, lon0, radius, tol=tol)

            def PlumbTo(self, lat0, lon0, est=None, tol=_TOL):
                '''Compute the I{perpendicular} intersection of this geodesic line
                   with a geodesic from the given point.

                   @arg lat0: Latitude of the point (C{degrees}).
                   @arg lon0: Longitude of the point (C{degrees}).
                   @kwarg est: Optional, initial estimate for the distance C{s12} of
                               the intersection I{along} this geodesic line (C{meter}).
                   @kwarg tol: Convergence tolerance (C(meter)).

                   @return: The intersection point on this geodesic line, a L{GDict}
                            from method L{Position} extended to 14 items C{lat1, lon1,
                            azi1, lat2, lon2, azi2, a12, s12, lat0, lon0, azi0, a02,
                            s02, at} with C{a02} and C{s02} the distance in C{degrees}
                            and C{meter} from the given point C{lat0, lon0} to the
                            intersection C{lat2, lon2}, azimuth C{azi0} at the given
                            point and the (perpendicular) angle C{at} between the
                            geodesic and this line at the intersection point.  The
                            geodesic azimuth at the intersection is C{(at + azi2)}.
                            See method L{Position} for further details.

                   @see: Methods C{Intersecant2}, C{Intersection} and C{Position}.
                '''
                return _PlumbTo(self, lat0, lon0, est=est, tol=tol)

            def Position(self, s12, outmask=Caps.STANDARD):
                '''Return the position at distance C{B{s12} meter} on this line.

                   @arg s12: Distance from this line's first point (C{meter}).
                   @kwarg outmask: Bit-or'ed combination of L{Caps<pygeodesy.karney.Caps>}
                                   values specifying the quantities to be returned.

                   @return: A L{GDict} with up to 12 items C{lat1, lon1, azi1, lat2,
                            lon2, azi2, m12, a12, s12, M12, M21, S12} with C{lat1},
                            C{lon1}, C{azi1} and arc length C{a12} always included,
                            except when C{a12=NAN}.
                '''
                with _wargs(self, s12, outmask) as args:
                    d = _GeodesicLine.Position(self, *args)
                return GDict(d)

        # GeodesicLine.ArcPosition.__doc__ = _GeodesicLine.ArcPosition.__doc__
        # GeodesicLine.Position.__doc__    = _GeodesicLine.Position.__doc__
        return GeodesicLine  # overwrite property_ROver

    @property_ROver
    def Geodesic_WGS84(self):
        '''Get the I{wrapped} C{Geodesic(WGS84)} singleton, provided the
           U{geographiclib<https://PyPI.org/project/geographiclib>} package
           is installed, otherwise an C{ImportError}.
        '''
        return _EWGS84.geodesicw  # overwrite property_ROver

_wrapped = _gWrapped()  # PYCHOK singleton, .ellipsoids, .test/base.py


class Geodesic(_gWrapped):  # overwritten by 1st instance
    '''I{Wrapper} around I{Karney}'s class U{geographiclib.geodesic.Geodesic
       <https://GeographicLib.SourceForge.io/Python/doc/code.html>}.
    '''
    def __new__(unused, a_ellipsoid=_EWGS84, f=None, **name):
        '''Return a I{wrapped} C{geodesic.Geodesic} instance from I{Karney}'s
           Python U{geographiclib<https://PyPI.org/project/geographiclib>},
           provide the latter is installed, otherwise an C{ImportError}.

           @arg a_ellipsoid: An ellipsoid (L{Ellipsoid}) or datum (L{Datum})
                  or the equatorial radius I{a} of the ellipsoid (C{meter}).
           @arg f: The flattening of the ellipsoid (C{scalar}), required if
                   B{C{a_ellipsoid}}) is C{meter}, ignored otherwise.
           @kwarg name: Optional C{B{name}=NN} (C{str}).
        '''
        g = _wrapped.Geodesic(a_ellipsoid, f=f, **name)
        _MODS.geodesicw.Geodesic = type(g)  # overwrite class
        return g


class GeodesicLine(_gWrapped):  # overwritten by 1st instance
    '''I{Wrapper} around I{Karney}'s class U{geographiclib.geodesicline.GeodesicLine
       <https://GeographicLib.SourceForge.io/Python/doc/code.html>}.
    '''
    def __new__(unused, geodesic, lat1, lon1, azi1, caps=Caps.STANDARD_LINE, **name):
        '''Return a I{wrapped} C{geodesicline.GeodesicLine} instance from I{Karney}'s
           Python U{geographiclib<https://PyPI.org/project/geographiclib>}, provided
           the latter is installed, otherwise an C{ImportError}.

           @arg geodesic: A I{wrapped} L{Geodesic} instance.
           @arg lat1: Latitude of the first points (C{degrees}).
           @arg lon1: Longitude of the first points (C{degrees}).
           @arg azi1: Azimuth at the first points (compass C{degrees360}).
           @kwarg caps: Optional, bit-or'ed combination of L{Caps<pygeodesy.karney.Caps>}
                        values specifying the capabilities the C{GeodesicLine} instance
                        should possess, i.e., which quantities can be returned by methods
                        C{GeodesicLine.Position} and C{GeodesicLine.ArcPosition}.
           @kwarg name: Optional C{B{name}=NN} (C{str}).
        '''
        gl = _wrapped.GeodesicLine(geodesic, lat1, lon1, azi1, caps=caps, **name)
        _MODS.geodesicw.GeodesicLine = type(gl)  # overwrite class
        return gl


def Geodesic_WGS84():
    '''Get the I{wrapped} L{Geodesic}C{(WGS84)} singleton, provided
       U{geographiclib<https://PyPI.org/project/geographiclib>} is
       installed, otherwise an C{ImportError}.
    '''
    return _wrapped.Geodesic_WGS84


class _wargs(object):  # see also .formy._idllmn6, .latlonBase._toCartesian3, .vector2d._numpy
    '''(INTERNAL) C{geographiclib} arguments and exception handler.
    '''
    @contextmanager  # <https://www.Python.org/dev/peps/pep-0343/> Examples
    def __call__(self, inst, *args, **kwds):
        '''(INTERNAL) Yield C{tuple(B{args})} with any errors raised
           as L{GeodesicError} embellished with all B{C{kwds}}.
        '''
        try:
            yield args
        except Exception as x:
            u = _DOT_(classname(inst), callername(up=2, underOK=True))
            raise GeodesicError(unstr(u, *args, **_name1__(kwds)), cause=x)

_wargs = _wargs()  # PYCHOK singleton


def _Intersecant2(gl, lat0, lon0, radius, tol=_TOL, form=F_D):  # MCCABE in LatLonEllipsoidalBaseDI.intersecant2, .geodesicx.gxline.Intersecant2
    # (INTERNAL) Return the intersections of a circle at C{lat0, lon0}
    # and a geodesic line as a 2-Tuple C{(P, Q)}, each a C{GDict}.
    r  = Radius_(radius)
    n  = typename(_Intersecant2)[1:]
    _P = gl.Position
    _I = gl.geodesic.Inverse

    def _R3(s):
        # radius, intersection, etc. at distance C{s}
        P = _P(s)
        d = _I(lat0, lon0, P.lat2, P.lon2)
        return fabs(d.s12), P, d

    def _bisect2(s, c, Rc, r, tol, _R3):
        _s = Fsum(c).fsumf_
        for i in range(_TRIPS):
            b = _s(s)
            Rb, P, d = _R3(b)
            if Rb > r:
                break
        else:  # b >>> s and c >>> s
            raise ValueError(Fmt.no_convergence(b, s))
        # Rb > r > Rc
        for i in range(_TRIPS):  # 47-48
            s = (b + c) * _0_5
            R, P, d = _R3(s)
            if Rb > R > r:
                b, Rb = s, R
            elif Rc < R < r:
                c, Rc = s, R
#           else:
#               break
            t = fabs(b - c)
            if t < tol:  # or fabs(R - r) < tol:
                break
        else:  # t = min(t, fabs(R - r))
            raise ValueError(Fmt.no_convergence(t, tol))
        i += C.iteration  # combine iterations
        P.set_(lat0=lat0, lon0=lon0, azi0=d.azi1, iteration=i,
               a02=d.a12, s02=d.s12, at=d.azi2 - P.azi2, name=n)
        return P, s

    # get the perpendicular intersection of 2 geodesics,
    # one the plumb, pseudo-rhumb line to the other
    C = _PlumbTo(gl, lat0, lon0, tol=tol)
    try:
        a = fabs(C.s02)  # distance between centers
        if a < r:
            c =  C.s12  # distance along pseudo-rhumb line
            h = _copysign(r, c)  # past half chord length
            P, p = _bisect2( h, c, a, r, tol, _R3)
            Q, q = _bisect2(-h, c, a, r, tol, _R3)
            if fabs(p - q) < max(EPS, tol):
                Q = P
        elif a > r:
            raise ValueError(_too_(Fmt.distant(a)))
        else:  # tangential
            P = Q = C
    except Exception as x:
        t = _LLB(C.lat2, C.lon2).toStr(form=form)
        t = _SPACE_(x, _plumb_, _to_, Fmt.PAREN(t))
        raise IntersectionError(t, txt=None, cause=x)

    return P, Q


def _PlumbTo(gl, lat0, lon0, est=None, tol=_TOL):
    # (INTERNAL) Return the I{perpendicular} intersection of
    # a geodesic line C{gl} and geodesic from C{(lat0, lon0)}.
    pl = _MODS.rhumb.bases._PseudoRhumbLine(gl)
    return pl.PlumbTo(lat0, lon0, exact=gl.geodesic,
                                  est=est, tol=tol)

# **) MIT License
#
# Copyright (C) 2016-2025 -- mrJean1 at Gmail -- All Rights Reserved.
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
