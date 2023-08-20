
# -*- coding: utf-8 -*-

u'''Wrapper around Python classes C{geodesic.Geodesic} and C{geodesicline.GeodesicLine} from
I{Karney}'s Python package U{geographiclib<https://PyPI.org/project/geographiclib>}, provided
that package is installed.

The I{wrapped} class methods return a L{GDict} instance offering access to the C{dict} items
either by C{key} or by C{attribute} name.

With env variable C{PYGEODESY_GEOGRAPHICLIB} left undefined or set to C{"2"}, this module,
L{pygeodesy.geodesicx} and L{pygeodesy.karney} will use U{GeographicLib 2.0
<https://GeographicLib.SourceForge.io/C++/doc/>} transcoding, otherwise C{1.52} or older.
'''

# from pygeodesy.basics import _xinstanceof  # from .karney
# from pygeodesy.constants import NAN  # from .karney
# from pygeodesy.datums import _a_ellipsoid  # from .karney
# from pygeodesy.ellipsoids import _EWGS84  # from .karney
# from pygeodesy.errors import _xkwds  # from .karney
from pygeodesy.interns import NN, _DOT_, _under
from pygeodesy.karney import _atan2d, Caps, Direct9Tuple, GDict, \
                             _EWGS84, GeodesicError, Inverse10Tuple, \
                             _kWrapped,  _a_ellipsoid, _ALL_LAZY, NAN, \
                             _xinstanceof, _xkwds  # PYCHOK used!
# from pygeodesy.lazily import _ALL_LAZY  # from .karney
from pygeodesy.named import callername, classname,  unstr
from pygeodesy.namedTuples import Destination3Tuple, Distance3Tuple
from pygeodesy.props import Property, Property_RO
# from pygeodesy.streps import unstr  # from .named
from pygeodesy.utily import _Wrap, wrap360,  fabs  # PYCHOK used!

from contextlib import contextmanager
# from math import fabs  # from .utily

__all__ = _ALL_LAZY.geodesicw
__version__ = '23.08.20'


class _gWrapped(_kWrapped):
    ''''(INTERNAL) Wrapper for some of I{Karney}'s U{geographiclib
        <https://PyPI.org/project/geographiclib>} classes.
    '''

    @Property_RO  # MCCABE 24
    def Geodesic(self):
        '''Get the I{wrapped} C{geodesic.Geodesic} class from I{Karney}'s Python
           U{geographiclib<https://GitHub.com/geographiclib/geographiclib-python>},
           provided the latter is installed.
        '''
        _Geodesic = self.geographiclib.Geodesic
        # assert Caps._STD == _Geodesic.STANDARD

        class Geodesic(_Geodesic):
            '''I{Wrapper} for I{Karney}'s Python U{geodesic.Geodesic
               <https://PyPI.org/project/geographiclib>} class.
            '''
            _debug   =  0  # like .geodesicx.bases._GeodesicBase
            _E       = _EWGS84
            LINE_OFF =  0  # in .azimuthal._GnomonicBase and .css.CassiniSoldner

            def __init__(self, a_ellipsoid=_EWGS84, f=None, name=NN):  # PYCHOK signature
                '''New I{wrapped} C{geodesic.Geodesic} instance.

                   @arg a_ellipsoid: An ellipsoid (L{Ellipsoid}) or datum (L{Datum})
                          or the equatorial radius I{a} of the ellipsoid (C{meter}).
                   @arg f: The flattening of the ellipsoid (C{scalar}), ignored if
                           B{C{a_ellipsoid}) is not specified as C{meter}.
                   @kwarg name: Optional ellipsoid name (C{str}), ignored like B{C{f}}.
                '''
                if a_ellipsoid not in (Geodesic._E, None):  # spherical OK
                    self._E = _a_ellipsoid(a_ellipsoid, f, name=name)  # raiser=NN
                with _wargs(self, *self.ellipsoid.a_f, name=name) as args:
                    _Geodesic.__init__(self, *args)

            def ArcDirect(self, lat1, lon1, azi1, a12, outmask=Caps._STD):
                '''Return the C{_Geodesic.ArcDirect} result as L{GDict}.
                '''
                with _wargs(self, lat1, lon1, azi1, a12, outmask) as args:
                    d = _Geodesic.ArcDirect(self, *args)
                return GDict(d)

            def ArcDirectLine(self, lat1, lon1, azi1, a12, caps=Caps._STD_LINE):
                '''Return the C{_Geodesic.ArcDirectLine} as I{wrapped} C{GeodesicLine}.
                '''
                return self._GenDirectLine(lat1, lon1, azi1, True, a12, caps)

            Area = _Geodesic.Polygon  # like GeodesicExact.Area

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

            def Direct(self, lat1, lon1, azi1, s12, outmask=Caps._STD):
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

            def DirectLine(self, lat1, lon1, azi1, s12, caps=Caps._STD_LINE):
                '''Return the C{_Geodesic.DirectLine} as I{wrapped} C{GeodesicLine}.
                '''
                return self._GenDirectLine(lat1, lon1, azi1, False, s12, caps)

            @Property_RO
            def ellipsoid(self):
                '''Get this geodesic's ellipsoid (C{Ellipsoid[2]}).
                '''
                return self._E

            @Property_RO
            def f1(self):  # in .css.CassiniSoldner.reset
                '''Get the geodesic's ellipsoid's I{1 - flattening} (C{float}).
                '''
                return getattr(self, _under(Geodesic.f1.name), self.ellipsoid.f1)

            def _GDictDirect(self, lat, lon, azi, arcmode, s12_a12, outmask=Caps._STD):
                '''(INTERNAL) Get C{_Geodesic._GenDirect} result as C{GDict}.
                '''
                with _wargs(self, lat, lon, azi, arcmode, s12_a12, outmask) as args:
                    t = _Geodesic._GenDirect(self, *args)
                return Direct9Tuple(t).toGDict()  # *t

            def _GDictInverse(self, lat1, lon1, lat2, lon2, outmask=Caps._STD):
                '''(INTERNAL) Get C{_Geodesic._GenInverse} result as L{Inverse10Tuple}.
                '''
                with _wargs(self, lat1, lon1, lat2, lon2, outmask) as args:
                    t = _Geodesic._GenInverse(self, *args)
                return Inverse10Tuple(t).toGDict(lon1=lon1, lon2=lon2)  # *t

            def _GenDirectLine(self, lat1, lon1, azi1, arcmode, s12_a12, *caps):
                '''(INTERNAL) Invoked by C{_Geodesic.DirectLine} and C{-.ArcDirectLine},
                   returning the result as a I{wrapped} C{GeodesicLine}.
                '''
                with _wargs(self, lat1, lon1, azi1, arcmode, s12_a12, *caps) as args:
                    t = _Geodesic._GenDirectLine(self, *args)
                return self._Line13(t)

            def Inverse(self, lat1, lon1, lat2, lon2, outmask=Caps._STD):
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

            def InverseLine(self, lat1, lon1, lat2, lon2, caps=Caps._STD_LINE):
                '''Return the C{_Geodesic.InverseLine} as I{wrapped} C{GeodesicLine}.
                '''
                with _wargs(self, lat1, lon1, lat2, lon2, caps) as args:
                    t = _Geodesic.InverseLine(self, *args)
                return self._Line13(t)

            def Line(self, lat1, lon1, azi1, caps=Caps._STD_LINE):
                '''Set up a I{wrapped} C{GeodesicLine} to compute several points
                   along a single, I{wrapped} (this) geodesic.
                '''
                return _wrapped.GeodesicLine(self, lat1, lon1, azi1, caps=caps)

            def _Line13(self, t):
                '''(INTERNAL) Wrap C{_GeodesicLine}, add distance and arc length
                   to reference point 3.
                '''
                rl = _wrapped.GeodesicLine(self, t.lat1, t.lon1, t.azi1, caps=t.caps,
                                           salp1=t.salp1, calp1=t.calp1)
                rl.a13, rl.s13 = t.a13, t.s13
                return rl

#           Polygon = _Geodesic.Polygon

        # Geodesic.ArcDirect.__doc__   = _Geodesic.ArcDirect.__doc__
        # Geodesic.Direct.__doc__      = _Geodesic.Direct.__doc__
        # Geodesic.Inverse.__doc__     = _Geodesic.Inverse.__doc__
        # Geodesic.InverseLine.__doc__ = _Geodesic.InverseLinr.__doc__
        # Geodesic.Line.__doc__        = _Geodesic.Line.__doc__
        return Geodesic

    @Property_RO  # MCCABE 16
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
            def __init__(self, geodesic, lat1, lon1, azi1, **caps):  # caps, salp1=NAN, calp1=NAN
                '''New I{wrapped} C{geodesicline.GeodesicLine} instance.

                   @arg geodesic: A I{wrapped} C{Geodesic} instance.
                   @arg lat1: Latitude of the first points (C{degrees}).
                   @arg lon1: Longitude of the first points (C{degrees}).
                   @arg azi1: Azimuth at the first points (compass C{degrees360}).
                   @kwarg caps: Optional, bit-or'ed combination of L{Caps} values
                                specifying the capabilities the C{GeodesicLine}
                                instance should possess, i.e., which quantities can
                                be returned by calls to C{GeodesicLine.Position}
                                and C{GeodesicLine.ArcPosition}.
                '''
                _xinstanceof(_wrapped.Geodesic, geodesic=geodesic)
                with _wargs(self, geodesic, lat1, lon1, azi1, **caps) as args:
                    _GeodesicLine.__init__(self, *args, **caps)

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
                '''Return the arc length to reference point 3 (C{degrees} or C{NAN}).
                '''
                return self.a13

            def ArcPosition(self, a12, outmask=Caps._STD):
                '''Return the position at arc length C{B{a12} degrees} on this line.
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

            def Position(self, s12, outmask=Caps._STD):
                '''Return the position at distance C{B{s12} meter} on this line.
                '''
                with _wargs(self, s12, outmask) as args:
                    d = _GeodesicLine.Position(self, *args)
                return GDict(d)

        # GeodesicLine.ArcPosition.__doc__ = _GeodesicLine.ArcPosition.__doc__
        # GeodesicLine.Position.__doc__    = _GeodesicLine.Position.__doc__
        return GeodesicLine

    @Property_RO
    def Geodesic_WGS84(self):
        '''Get the I{wrapped} C{Geodesic(WGS84)} singleton, provided the
           U{geographiclib<https://PyPI.org/project/geographiclib>} package
           is installed, otherwise an C{ImportError}.
        '''
        return _EWGS84.geodesic

_wrapped = _gWrapped()  # PYCHOK singleton, .ellipsoids, .test/base.py


def Geodesic(a_ellipsoid, f=None, name=NN):
    '''Return a I{wrapped} C{geodesic.Geodesic} instance from I{Karney}'s
       Python U{geographiclib<https://PyPI.org/project/geographiclib>},
       provide the latter is installed, otherwise an C{ImportError}.

       @arg a_ellipsoid: An ellipsoid (L{Ellipsoid}) or datum (L{Datum})
              or the equatorial radius I{a} of the ellipsoid (C{meter}).
       @arg f: The flattening of the ellipsoid (C{scalar}), ignored if
               B{C{a_ellipsoid}}) is not specified as C{meter}.
       @kwarg name: Optional ellipsoid name (C{str}), ignored like B{C{f}}.
    '''
    return _wrapped.Geodesic(a_ellipsoid, f=f, name=name)


def GeodesicLine(geodesic, lat1, lon1, azi1, caps=Caps._STD_LINE):
    '''Return a I{wrapped} C{geodesicline.GeodesicLine} instance from I{Karney}'s
       Python U{geographiclib<https://PyPI.org/project/geographiclib>}, provided
       the latter is installed, otherwise an C{ImportError}.

       @arg geodesic: A I{wrapped} L{Geodesic} instance.
       @arg lat1: Latitude of the first points (C{degrees}).
       @arg lon1: Longitude of the first points (C{degrees}).
       @arg azi1: Azimuth at the first points (compass C{degrees360}).
       @kwarg caps: Optional, bit-or'ed combination of L{Caps} values
                    specifying the capabilities the C{GeodesicLine}
                    instance should possess, i.e., which quantities can
                    be returned by calls to C{GeodesicLine.Position}
                    and C{GeodesicLine.ArcPosition}.
    '''
    return _wrapped.GeodesicLine(geodesic, lat1, lon1, azi1, caps=caps)


def Geodesic_WGS84():
    '''Get the I{wrapped} L{Geodesic}C{(WGS84)} singleton, provided
       U{geographiclib<https://PyPI.org/project/geographiclib>} is
       installed, otherwise an C{ImportError}.
    '''
    return _wrapped.Geodesic_WGS84


class _wargs(object):  # see also .vector2d._numpy
    '''(INTERNAL) C{geographiclib} caller, catching exceptions.
    '''
    @contextmanager  # <https://www.python.org/dev/peps/pep-0343/> Examples
    def __call__(self, inst, *args, **kwds):
        '''(INTERNAL) Yield C{tuple(B{args})} with any errors raised as L{NumPyError}.
        '''
        try:
            yield args
        except (AttributeError, TypeError, ValueError) as x:
            n = _DOT_(classname(inst), callername(up=3, underOK=True))
            raise GeodesicError(unstr(n, *args, **kwds), cause=x)

_wargs = _wargs()  # PYCHOK singleton


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
