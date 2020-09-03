
# -*- coding: utf-8 -*-

u'''Classes L{Datum} and L{Transform} and registries thereof L{Datums}
and L{Transforms}, respectively.

Pure Python implementation of geodesy tools for ellipsoidal earth models,
including datums and ellipsoid parameters for different geographic coordinate
systems and methods for converting between them and to cartesian coordinates.
Transcribed from JavaScript originals by I{(C) Chris Veness 2005-2016} and
published under the same MIT Licence**, see U{latlon-ellipsoidal.js
<https://www.Movable-Type.co.UK/scripts/geodesy/docs/latlon-ellipsoidal.js.html>}.

Historical geodetic datums: a latitude/longitude point defines a geographic
location on or above/below the earth’s surface, measured in degrees from
the equator and the International Reference Meridian and meters above the
ellipsoid, and based on a given datum.  The datum is based on a reference
ellipsoid and tied to geodetic survey reference points.

Modern geodesy is generally based on the WGS84 datum (as used for instance
by GPS systems), but previously various reference ellipsoids and datum
references were used.

The UK Ordnance Survey National Grid References are still based on the otherwise
historical OSGB36 datum, q.v. U{Ordnance Survey 'A guide to coordinate systems
in Great Britain', Section 6<https://www.OrdnanceSurvey.co.UK/docs/support/
guide-coordinate-systems-great-britain.pdf>}.

@newfield example: Example, Examples

@var Datums.BD72: Datum(name='BD72', ellipsoid=Ellipsoids.Intl1924, transform=Transforms.BD72)
@var Datums.DHDN: Datum(name='DHDN', ellipsoid=Ellipsoids.Bessel1841, transform=Transforms.DHDN)
@var Datums.ED50: Datum(name='ED50', ellipsoid=Ellipsoids.Intl1924, transform=Transforms.ED50)
@var Datums.GDA2020: Datum(name='GDA2020', ellipsoid=Ellipsoids.GRS80, transform=Transforms.WGS84)
@var Datums.GRS80: Datum(name='GRS80', ellipsoid=Ellipsoids.GRS80, transform=Transforms.WGS84)
@var Datums.Irl1975: Datum(name='Irl1975', ellipsoid=Ellipsoids.AiryModified, transform=Transforms.Irl1975)
@var Datums.Krassovski1940: Datum(name='Krassovski1940', ellipsoid=Ellipsoids.Krassovski1940, transform=Transforms.Krassovski1940)
@var Datums.Krassowsky1940: Datum(name='Krassowsky1940', ellipsoid=Ellipsoids.Krassowsky1940, transform=Transforms.Krassowsky1940)
@var Datums.MGI: Datum(name='MGI', ellipsoid=Ellipsoids.Bessel1841, transform=Transforms.MGI)
@var Datums.NAD27: Datum(name='NAD27', ellipsoid=Ellipsoids.Clarke1866, transform=Transforms.NAD27)
@var Datums.NAD83: Datum(name='NAD83', ellipsoid=Ellipsoids.GRS80, transform=Transforms.NAD83)
@var Datums.NTF: Datum(name='NTF', ellipsoid=Ellipsoids.Clarke1880IGN, transform=Transforms.NTF)
@var Datums.OSGB36: Datum(name='OSGB36', ellipsoid=Ellipsoids.Airy1830, transform=Transforms.OSGB36)
@var Datums.Potsdam: Datum(name='Potsdam', ellipsoid=Ellipsoids.Bessel1841, transform=Transforms.Bessel1841)
@var Datums.Sphere: Datum(name='Sphere', ellipsoid=Ellipsoids.Sphere, transform=Transforms.WGS84)
@var Datums.TokyoJapan: Datum(name='TokyoJapan', ellipsoid=Ellipsoids.Bessel1841, transform=Transforms.TokyoJapan)
@var Datums.WGS72: Datum(name='WGS72', ellipsoid=Ellipsoids.WGS72, transform=Transforms.WGS72)
@var Datums.WGS84: Datum(name='WGS84', ellipsoid=Ellipsoids.WGS84, transform=Transforms.WGS84)

@var Transforms.BD72: Transform(name='BD72', tx=106.86863, ty=-52.29778, tz=103.72389, rx=-0, ry=-0, rz=-0.00001, s=1.2727, s1=1, sx=-0.33657, sy=-0.45696, sz=-1.84218)
@var Transforms.Bessel1841: Transform(name='Bessel1841', tx=-582, ty=-105, tz=-414, rx=-0.00001, ry=-0, rz=0.00001, s=-8.3, s1=0.99999, sx=-1.04, sy=-0.35, sz=3.08)
@var Transforms.Clarke1866: Transform(name='Clarke1866', tx=8, ty=-160, tz=-176, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
@var Transforms.DHDN: Transform(name='DHDN', tx=-591.28, ty=-81.35, tz=-396.39, rx=0.00001, ry=-0, rz=-0.00001, s=-9.82, s1=0.99999, sx=1.477, sy=-0.0736, sz=-1.458)
@var Transforms.ED50: Transform(name='ED50', tx=89.5, ty=93.8, tz=123.1, rx=0, ry=0, rz=0, s=-1.2, s1=1, sx=0, sy=0, sz=0.156)
@var Transforms.Identity: Transform(name='Identity', tx=0, ty=0, tz=0, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
@var Transforms.Irl1965: Transform(name='Irl1965', tx=-482.53, ty=130.596, tz=-564.557, rx=0.00001, ry=0, rz=0, s=-8.15, s1=0.99999, sx=1.042, sy=0.214, sz=0.631)
@var Transforms.Irl1975: Transform(name='Irl1975', tx=-482.53, ty=130.596, tz=-564.557, rx=-0.00001, ry=-0, rz=-0, s=-1.1, s1=1, sx=-1.042, sy=-0.214, sz=-0.631)
@var Transforms.Krassovski1940: Transform(name='Krassovski1940', tx=-24, ty=123, tz=94, rx=-0, ry=0, rz=0, s=-2.423, s1=1, sx=-0.02, sy=0.26, sz=0.13)
@var Transforms.Krassowsky1940: Transform(name='Krassowsky1940', tx=-24, ty=123, tz=94, rx=-0, ry=0, rz=0, s=-2.423, s1=1, sx=-0.02, sy=0.26, sz=0.13)
@var Transforms.MGI: Transform(name='MGI', tx=-577.326, ty=-90.129, tz=-463.92, rx=0.00002, ry=0.00001, rz=0.00003, s=-2.423, s1=1, sx=5.137, sy=1.474, sz=5.297)
@var Transforms.NAD27: Transform(name='NAD27', tx=8, ty=-160, tz=-176, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
@var Transforms.NAD83: Transform(name='NAD83', tx=1.004, ty=-1.91, tz=-0.515, rx=0, ry=0, rz=0, s=-0.0015, s1=1, sx=0.0267, sy=0.00034, sz=0.011)
@var Transforms.NTF: Transform(name='NTF', tx=-168, ty=-60, tz=320, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
@var Transforms.OSGB36: Transform(name='OSGB36', tx=-446.448, ty=125.157, tz=-542.06, rx=-0, ry=-0, rz=-0, s=20.4894, s1=1.00002, sx=-0.1502, sy=-0.247, sz=-0.8421)
@var Transforms.TokyoJapan: Transform(name='TokyoJapan', tx=148, ty=-507, tz=-685, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
@var Transforms.WGS72: Transform(name='WGS72', tx=0, ty=0, tz=-4.5, rx=0, ry=0, rz=0, s=-0.22, s1=1, sx=0, sy=0, sz=0.554)
@var Transforms.WGS84: Transform(name='WGS84', tx=0, ty=0, tz=0, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
'''

# make sure int/int division yields float quotient
from __future__ import division
division = 1 / 2  # double check int division, see .ellipsoids, .fmath, .utily
if not division:
    raise ImportError('%s 1/2 == %d' % ('division', division))
del division

from pygeodesy.basics import isscalar, property_RO, _xinstanceof
from pygeodesy.ellipsoids import a_f2Tuple, _4Ecef, Ellipsoid, \
                                 Ellipsoid2, Ellipsoids
from pygeodesy.errors import _IsnotError
from pygeodesy.fmath import _2_3rd, cbrt, cbrt2, fdot, fpowers, Fsum, \
                             fsum_, hypot1, hypot2, sqrt3  # PYCHOK _2_3rd
from pygeodesy.interns import _COMMA_SPACE_, _ellipsoid_, _flt, _name_, NN, \
                              _spherical_, _transform_, _UNDERSCORE_
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import _NamedEnum, _NamedEnumItem, Vector3Tuple
from pygeodesy.units import Radius_, Scalar

from math import radians

__all__ = _ALL_LAZY.datums
__version__ = '20.09.01'


def _r_s2(s):
    '''(INTERNAL) rotation in C{radians} and C{degree seconds}.
    '''
    return _flt(radians(s / 3600.0)), _flt(s)


class Transform(_NamedEnumItem):
    '''Helmert transformation.
    '''
    tx = 0  #: X translation (C{meter}).
    ty = 0  #: Y translation (C{meter}).
    tz = 0  #: Z translation (C{meter}).

    rx = 0  #: X rotation (C{radians}).
    ry = 0  #: Y rotation (C{radians}).
    rz = 0  #: Z rotation (C{radians}).

    s  = 0  #: Scale ppm (C{float}).
    s1 = 1  #: Scale + 1 (C{float}).

    sx = 0  #: X rotation (degree seconds).
    sy = 0  #: Y rotation (degree seconds).
    sz = 0  #: Z rotation (degree seconds).

    def __init__(self, name=NN, tx=0, ty=0, tz=0,
                                sx=0, sy=0, sz=0, s=0):
        '''New L{Transform}.

           @kwarg name: Optional, unique name (C{str}).
           @kwarg tx: Optional X translation (C{meter}).
           @kwarg ty: Optional Y translation (C{meter}).
           @kwarg tz: Optional Z translation (C{meter}).
           @kwarg s: Optional scale ppm (C{float}).
           @kwarg sx: Optional X rotation (C{degree seconds}).
           @kwarg sy: Optional Y rotation (C{degree seconds}).
           @kwarg sz: Optional Z rotation (C{degree seconds}).

           @raise NameError: Transform with that B{C{name}} already exists.
        '''
        if tx:
            self.tx = _flt(tx)
        if ty:
            self.ty = _flt(ty)
        if tz:
            self.tz = _flt(tz)
        if sx:  # secs to rads
            self.rx, self.sx = _r_s2(sx)
        if sy:
            self.ry, self.sy = _r_s2(sy)
        if sz:
            self.rz, self.sz = _r_s2(sz)
        if s:
            self.s  = _flt(s)
            self.s1 = _flt(s * 1e-6 + 1)  # normalize ppm to (s + 1)

        self._register(Transforms, name)

    def __eq__(self, other):
        '''Compare this and an other transform.

           @arg other: The other transform (L{Transform}).

           @return: C{True} if equal, C{False} otherwise.
        '''
        return self is other or (isinstance(other, Transform) and
                                 self.tx == other.tx and
                                 self.ty == other.ty and
                                 self.tz == other.tz and
                                 self.rx == other.rx and
                                 self.ry == other.ry and
                                 self.rz == other.rz and
                                 self.s  == other.s)

    def inverse(self, name=NN):
        '''Return the inverse of this transform.

           @kwarg name: Optional, unique name (C{str}).

           @return: Inverse (Transform).

           @raise NameError: Transform with that B{C{name}} already exists.
        '''
        return Transform(name=name or (self.name + 'Inverse'),
                         tx=-self.tx, ty=-self.ty, tz=-self.tz,
                         sx=-self.sx, sy=-self.sy, sz=-self.sz, s=-self.s)

    def toStr(self, prec=5):  # PYCHOK expected
        '''Return this transform as a string.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).

           @return: Transform attributes (C{str}).
        '''
        return self._instr(prec, 'tx', 'ty', 'tz',
                                 'rx', 'ry', 'rz', 's', 's1',
                                 'sx', 'sy', 'sz')

    def transform(self, x, y, z, inverse=False):
        '''Transform a (geocentric) Cartesian point, forward or inverse.

           @arg x: X coordinate (C{meter}).
           @arg y: Y coordinate (C{meter}).
           @arg z: Z coordinate (C{meter}).
           @kwarg inverse: Optional direction, forward or inverse (C{bool}).

           @return: A L{Vector3Tuple}C{(x, y, z)}, transformed.
        '''
        if inverse:
            _xyz = -1, -x, -y, -z
            _s1 = self.s1 - 2  # == -(1 - s * 1e-6)) == -(1 - (s1 - 1))
        else:
            _xyz =  1,  x,  y,  z
            _s1  = self.s1
        # x', y', z' = (.tx + x * .s1 - y * .rz + z * .ry,
        #               .ty + x * .rz + y * .s1 - z * .rx,
        #               .tz - x * .ry + y * .rx + z * .s1)
        r = Vector3Tuple(fdot(_xyz, self.tx,      _s1, -self.rz,  self.ry),
                         fdot(_xyz, self.ty,  self.rz,      _s1, -self.rx),
                         fdot(_xyz, self.tz, -self.ry,  self.rx,      _s1))
        return self._xnamed(r)


Transforms = _NamedEnum('Transforms', Transform)  #: Registered transforms.
# <https://WikiPedia.org/wiki/Helmert_transformation> from WGS84
Transforms._assert(
    BD72           = Transform('BD72', tx=106.868628, ty=-52.297783, tz=103.723893,
                     # <https://www.NGI.BE/FR/FR4-4.shtm> ETRS89 == WG84
                     # <https://GeoRepository.com/transformation_15929/BD72-to-WGS-84-3.html>
                                       sx=-0.33657,   sy= -0.456955, sz= -1.84218,
                                        s= 1.2727),
    Bessel1841     = Transform('Bessel1841', tx=-582.0,  ty=-105.0, tz=-414.0,
                                             sx=  -1.04, sy= -0.35, sz=   3.08,
                                              s=  -8.3),
    Clarke1866     = Transform('Clarke1866', tx=8, ty=-160, tz=-176),
    DHDN           = Transform('DHDN', tx=-591.28,  ty=-81.35,   tz=-396.39,
                                       sx=   1.477, sy= -0.0736, sz=  -1.458,
                                        s=  -9.82),  # Germany
    ED50           = Transform('ED50', tx=89.5, ty=93.8, tz=123.1,
                     # <https://GeoNet.ESRI.com/thread/36583> sz=-0.156
                     # <https://GitHub.com/ChrisVeness/geodesy/blob/master/latlon-ellipsoidal.js>
                     # <https://www.Gov.UK/guidance/oil-and-gas-petroleum-operations-notices#pon-4>
                                                         sz=  0.156, s=-1.2),
    Identity       = Transform('Identity'),
    Irl1965        = Transform('Irl1965', tx=-482.530, ty=130.596, tz=-564.557,
                                          sx=   1.042, sy=  0.214, sz=   0.631,
                                           s=  -8.15),
    Irl1975        = Transform('Irl1975', tx=-482.530, ty=130.596, tz=-564.557,
                     # XXX rotation signs may be opposite, to be checked
                                          sx=  -1.042, sy= -0.214, sz=  -0.631,
                                           s=  -1.1),
    Krassovski1940 = Transform('Krassovski1940', tx=-24.0,  ty=123.0,  tz=94.0,
                                                 sx= -0.02, sy=  0.26, sz= 0.13,
                                                  s= -2.423),  # spelling
    Krassowsky1940 = Transform('Krassowsky1940', tx=-24.0,  ty=123.0,  tz=94.0,
                                                 sx= -0.02, sy=  0.26, sz= 0.13,
                                                  s= -2.423),  # spelling
    MGI            = Transform('MGI', tx=-577.326, ty=-90.129, tz=-463.920,
                                      sx=   5.137, sy=  1.474, sz=   5.297,
                                       s=  -2.423),  # Austria
    NAD27          = Transform('NAD27', tx=8, ty=-160, tz=-176),
    NAD83          = Transform('NAD83', tx= 1.004,  ty=-1.910,   tz=-0.515,
                                        sx= 0.0267, sy= 0.00034, sz= 0.011,
                                         s=-0.0015),
    NTF            = Transform('NTF', tx=-168, ty= -60, tz=320),  # XXX verify
    OSGB36         = Transform('OSGB36', tx=-446.448,  ty=125.157,  tz=-542.060,
                                         sx=  -0.1502, sy= -0.2470, sz=  -0.8421,
                                          s=  20.4894),
    TokyoJapan     = Transform('TokyoJapan', tx=148, ty=-507, tz=-685),
    WGS72          = Transform('WGS72', tz=-4.5, sz=0.554, s=-0.22),
    WGS84          = Transform('WGS84'),  # unity
)


class Datum(_NamedEnumItem):
    '''Ellipsoid and transform parameters for an earth model.
    '''
    _ellipsoid = Ellipsoids.WGS84  #: (INTERNAL) Default ellipsoid (L{Ellipsoid}, L{Ellipsoid2}).
    _exactTM   = None              #: (INTERNAL) L{ExactTransverseMercator} projection.
    _transform = Transforms.WGS84  #: (INTERNAL) Default transform (L{Transform}).

    def __init__(self, ellipsoid, transform=None, name=NN):
        '''New L{Datum}.

           @arg ellipsoid: The ellipsoid (L{Ellipsoid} or L{Ellipsoid2}).
           @kwarg transform: Optional transform (L{Transform}).
           @kwarg name: Optional, unique name (C{str}).

           @raise NameError: Datum with that B{C{name}} already exists.

           @raise TypeError: If B{C{ellipsoid}} is not an L{Ellipsoid}
                             nor L{Ellipsoid2} or B{C{transform}} is
                             not a L{Transform}.
        '''
        self._ellipsoid = ellipsoid or Datum._ellipsoid
        _xinstanceof(Ellipsoid, ellipsoid=self.ellipsoid)

        self._transform = transform or Datum._transform
        _xinstanceof(Transform, transform=self.transform)

        self._register(Datums, name or self.transform.name or self.ellipsoid.name)

    def __eq__(self, other):
        '''Compare this and an other datum.

           @arg other: The other datum (L{Datum}).

           @return: C{True} if equal, C{False} otherwise.
        '''
        return self is other or (isinstance(other, Datum) and
                                 self.ellipsoid == other.ellipsoid and
                                 self.transform == other.transform)

    def ecef(self, Ecef=None):
        '''Return U{ECEF<https://WikiPedia.org/wiki/ECEF>} converter.

           @kwarg Ecef: ECEF class to use (L{EcefKarney}, L{EcefVeness}
                        or L{EcefYou}).

           @return: An ECEF converter for this C{datum} (L{EcefKarney},
                    L{EcefVeness} or L{EcefYou}).

           @raise TypeError: Invalid B{C{Ecef}}.
        '''
        return _4Ecef(self, Ecef)

    @property_RO
    def ellipsoid(self):
        '''Get this datum's ellipsoid (L{Ellipsoid} or L{Ellipsoid2}).
        '''
        return self._ellipsoid

    @property_RO
    def exactTM(self):
        '''Get the C{ExactTM} projection (L{ExactTransverseMercator}).
        '''
        if self._exactTM is None:
            from pygeodesy.etm import ExactTransverseMercator
            self._exactTM = ExactTransverseMercator(datum=self)
        return self._exactTM

    @property_RO
    def isEllipsoidal(self):
        '''Check whether this datum is ellipsoidal (C{bool}).
        '''
        return self._ellipsoid.isEllipsoidal

    @property_RO
    def isSpherical(self):
        '''Check whether this datum is spherical (C{bool}).
        '''
        return self._ellipsoid.isSpherical

    def toStr(self, **unused):  # PYCHOK expected
        '''Return this datum as a string.

           @return: Datum attributes (C{str}).
        '''
        t = ['%s=%r' % (_name_, self.named)]
        for a in (_ellipsoid_, _transform_):
            v = getattr(self, a)
            t.append('%s=%ss.%s' % (a, v.classname, v.name))
        return _COMMA_SPACE_.join(t)

    @property_RO
    def transform(self):
        '''Get this datum's transform (L{Transform}).
        '''
        return self._transform


def _En2(arg, name):
    '''(INTERNAL) Helper for C{_ellipsoid} amd C{_ellipsoidal_datum}.
    '''
    if isinstance(arg, (Ellipsoid, Ellipsoid2)):
        E = arg
        n = _UNDERSCORE_ + (name or E.name)
    elif isinstance(arg, Datum):
        E = arg.ellipsoid
        n = _UNDERSCORE_ + (name or arg.name)
    elif isinstance(arg, a_f2Tuple):
        n = _UNDERSCORE_ + (name or arg.name)
        E = Ellipsoid(arg.a, arg.b, name=n)
    elif isinstance(arg, (tuple, list)) and len(arg) == 2:
        n = _UNDERSCORE_ + (name or getattr(arg, _name_, NN))
        a_f = a_f2Tuple(*arg)
        E = Ellipsoid(a_f.a, a_f.b, name=n)  # PYCHOK .a
    else:
        E, n = None, NN
    return E, n


def _ellipsoid(ellipsoid, name=NN):  # in .trf
    '''(INTERNAL) Create an L{Ellipsoid} or L{Ellipsoid2} from L{datum} or C{a_f2Tuple}.
    '''
    E, _ = _En2(ellipsoid, name)
    if not E:
        _xinstanceof(Ellipsoid, Ellipsoid2, a_f2Tuple, Datum, ellipsoid=ellipsoid)
    return E


def _ellipsoidal_datum(a_f, name=NN):
    '''(INTERNAL) Create a L{Datum} from an L{Ellipsoid} or L{Ellipsoid2} or C{a_f2Tuple}.
    '''
    if isinstance(a_f, Datum):
        return a_f
    E, n = _En2(a_f, name)
    if not E:
        _xinstanceof(Datum, Ellipsoid, Ellipsoid2, a_f2Tuple, datum=a_f)
    return Datum(E, transform=Transforms.Identity, name=n)


def _spherical_datum(radius, name=NN, raiser=False):
    '''(INTERNAL) Create a L{Datum} from an L{Ellipsoid}, L{Ellipsoid2} or scalar earth C{radius}.
    '''
    try:
        d = _ellipsoidal_datum(radius, name=name)
    except TypeError:
        d = None
    if d is None:
        if not isscalar(radius):
            _xinstanceof(Datum, Ellipsoid, Ellipsoid2, a_f2Tuple, Scalar, datum=radius)
        n = _UNDERSCORE_ + name
        r = Radius_(radius, Error=TypeError)
        E = Ellipsoid(r, r, name=n)
        d = Datum(E, transform=Transforms.Identity, name=n)
    elif raiser and not d.isSpherical:  # raiser if no spherical
        raise _IsnotError(_spherical_, datum=radius)
    return d


Datums = _NamedEnum('Datums', Datum)      #: Registered datums.
# Datums with associated ellipsoid and Helmert transform parameters
# to convert from WGS84 into the given datum.  More are available at
# <https://Earth-Info.NGA.mil/GandG/coordsys/datums/NATO_DT.pdf> and
# <XXX://www.FieldenMaps.info/cconv/web/cconv_params.js>.
Datums._assert(
    # Belgian Datum 1972, based on Hayford ellipsoid.
    # <https://NL.WikiPedia.org/wiki/Belgian_Datum_1972>
    # <https://SpatialReference.org/ref/sr-org/7718/html/>
    BD72           = Datum(Ellipsoids.Intl1924, Transforms.BD72),

    # Germany <https://WikiPedia.org/wiki/Bessel-Ellipsoid>
    #         <https://WikiPedia.org/wiki/Helmert_transformation>
    DHDN           = Datum(Ellipsoids.Bessel1841, Transforms.DHDN),

    # <https://www.Gov.UK/guidance/oil-and-gas-petroleum-operations-notices#pon-4>
    ED50           = Datum(Ellipsoids.Intl1924, Transforms.ED50),

    # Australia <https://ICSM.Gov.AU/datum/gda2020-and-gda94-technical-manuals>
#   ADG66          = Datum(Ellipsoids.ANS, Transforms.WGS84, name='ADG66'),  # XXX Transform?
#   ADG84          = Datum(Ellipsoids.ANS, Transforms.WGS84, name='ADG84'),  # XXX Transform?
#   GDA94          = Datum(Ellipsoids.GRS80, Transforms.WGS84, name='GDA94'),
    GDA2020        = Datum(Ellipsoids.GRS80, Transforms.WGS84, name='GDA2020'),  # XXX Transform?

    # <https://WikiPedia.org/wiki/GRS_80>
    GRS80          = Datum(Ellipsoids.GRS80, Transforms.WGS84, name='GRS80'),

    # <https://OSI.IE/OSI/media/OSI/Content/Publications/transformations_booklet.pdf>
    Irl1975        = Datum(Ellipsoids.AiryModified, Transforms.Irl1975),

    # Germany <https://WikiPedia.org/wiki/Helmert_transformation>
    Krassovski1940 = Datum(Ellipsoids.Krassovski1940, Transforms.Krassovski1940),  # XXX spelling?
    Krassowsky1940 = Datum(Ellipsoids.Krassowsky1940, Transforms.Krassowsky1940),  # XXX spelling?

    # Austria <https://DE.WikiPedia.org/wiki/Datum_Austria>
    MGI            = Datum(Ellipsoids.Bessel1841, Transforms.MGI),

    # <https://WikiPedia.org/wiki/Helmert_transformation>
    NAD27          = Datum(Ellipsoids.Clarke1866, Transforms.NAD27),

    # NAD83 (2009) == WGS84 - <https://www.UVM.edu/giv/resources/WGS84_NAD83.pdf>
    # (If you *really* must convert WGS84<->NAD83, you need more than this!)
    NAD83          = Datum(Ellipsoids.GRS80, Transforms.NAD83),

    #  Nouvelle Triangulation Francaise (Paris)  XXX verify
    NTF            = Datum(Ellipsoids.Clarke1880IGN, Transforms.NTF),

    # <https://www.OrdnanceSurvey.co.UK/docs/support/guide-coordinate-systems-great-britain.pdf>
    OSGB36         = Datum(Ellipsoids.Airy1830, Transforms.OSGB36),

    # Germany <https://WikiPedia.org/wiki/Helmert_transformation>
    Potsdam        = Datum(Ellipsoids.Bessel1841, Transforms.Bessel1841, name='Potsdam'),

    # XXX psuedo-ellipsoids for spherical LatLon
    Sphere         = Datum(Ellipsoids.Sphere, Transforms.WGS84, name='Sphere'),

    # <https://www.GeoCachingToolbox.com?page=datumEllipsoidDetails>
    TokyoJapan     = Datum(Ellipsoids.Bessel1841, Transforms.TokyoJapan),

    # <https://www.ICAO.int/safety/pbn/documentation/eurocontrol/eurocontrol%20wgs%2084%20implementation%20manual.pdf>
    WGS72          = Datum(Ellipsoids.WGS72, Transforms.WGS72),

    WGS84          = Datum(Ellipsoids.WGS84, Transforms.WGS84),
)

if __name__ == '__main__':

    # __doc__ of this file
    for e in (Datums, Transforms):
        t = [NN] + repr(e).split('\n')
        print('\n@var '.join(i.strip(',') for i in t))

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

# % python -m pygeodesy.datums
