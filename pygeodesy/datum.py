
# -*- coding: utf-8 -*-

u'''Classes L{Datum}, L{Ellipsoid} and L{Transform} and registries thereof.

Pure Python implementation of geodesy tools for ellipsoidal earth models,
including datums and ellipsoid parameters for different geographic coordinate
systems and methods for converting between them and to cartesian coordinates.
Transcribed from JavaScript originals by I{(C) Chris Veness 2005-2016} and
published under the same MIT Licence**, see U{http://www.movable-type.co.uk/
scripts/geodesy/docs/latlon-ellipsoidal.js.html}.

Historical geodetic datums: a latitude/longitude point defines a geographic
location on or above/below the earth’s surface, measured in degrees from
the equator and the International Reference Meridian and meters above the
ellipsoid, and based on a given datum.  The datum is based on a reference
ellipsoid and tied to geodetic survey reference points.

Modern geodesy is generally based on the WGS84 datum (as used for instance
by GPS systems), but previously various reference ellipsoids and datum
references were used.

The UK Ordnance Survey National Grid References are still based on the
otherwise historical OSGB36 datum, q.v. Ordnance Survey 'A guide to
coordinate systems in Great Britain' Section 6 U{http://www.ordnancesurvey
.co.uk/docs/support/guide-coordinate-systems-great-britain.pdf} and also
U{http://www.ordnancesurvey.co.uk/blog/2014/12/2}.

@newfield example: Example, Examples
'''

# force int division to yield float quotient
from __future__ import division as _
if not 1/2:  # PYCHOK 1/2 == 0
    raise ImportError('1/2 == %d' % (1/2,))

from bases import Base, Named
from utils import R_M, cbrt, cbrt2, fdot, fStr, \
                  m2km, m2NM, m2SM, radians

from math import atanh, sqrt

R_M  = R_M        #: Mean, spherical earth radius (meter).
R_KM = m2km(R_M)  #: Mean, spherical earth radius (kilo meter).
R_NM = m2NM(R_M)  #: Mean, spherical earth radius (nautical miles).
R_SM = m2SM(R_M)  #: Mean, spherical earth radius (statute miles).
# <http://www.edwilliams.org/avform.htm#XTE>
# R_AV = 6366707.0194937  #: Aviation earth radius (meter)
# R_?? = 6372797.560856   #: XXX some other earth radius?

# all public contants, classes and functions
__all__ = ('R_KM', 'R_M', 'R_NM', 'R_SM',  # constants
           'Datum',  'Ellipsoid',  'Transform',  # classes
           'Datums', 'Ellipsoids', 'Transforms')  # enum-like
__version__ = '17.06.25'


class _Enum(dict, Named):
    '''(INTERNAL) Enum-like dict sub-class.
    '''
    def __init__(self, name):
        '''New Enum.

           @param name: Name (string).
        '''
        self._name = name

    def __getattr__(self, attr):
        try:
            return self[attr]
        except KeyError:
            raise AttributeError("%s.%s doesn't exist" % (self._name, attr))

    def __repr__(self):
        return '\n'.join('%s.%s: %r' % (self._name, n, v) for n, v in sorted(self.items()))

    def __str__(self):
        return self._name + ', '.join(sorted('.' + n for n in self.keys()))

    def _assert(self, **kwds):
        '''(INTERNAL) Check names against given names.
        '''
        for a, v in kwds.items():
            assert getattr(self, a) is v

    def unregister(self, name):
        '''Remove a registered instance.

           @param name: Name of the instance (string).

           @return: The unregistered instance.

           @raise NameError: No instance with that name.
        '''
        try:
            return dict.pop(self, name)
        except KeyError:
            raise NameError('no %s %r' % (self._name, name))


Datums     = _Enum('Datums')      #: Registered datums (L{_Enum}).
Ellipsoids = _Enum('Ellipsoids')  #: Registered ellipsoids (L{_Enum}).
Transforms = _Enum('Transforms')  #: Registered transforms (L{_Enum}).


class _Based(Base, Named):
    '''(INTERNAL) Base class.
    '''
    def __ne__(self, other):
        '''Compare this and an other ellipsoid.

           @return: True if different (bool).
        '''
        return not self.__eq__(other)

    def _fStr(self, prec, *attrs, **others):
        '''(INTERNAL) Format.
        '''
        t = fStr([getattr(self, a) for a in attrs], prec=prec, sep=' ', ints=True)
        t = ['%s=%s' % (a, v) for a, v in zip(attrs, t.split())]
        if others:
            t += ['%s=%s' % (a, v) for a, v in sorted(others.items())]
        return ', '.join(['name=%r' % (self.name,)] + t)

    def _register(self, enum, name):
        '''(INTERNAL) Add this as enum.name.
        '''
        if name:
            self._name = name
            if name[:1].isalpha():
                if name in enum:
                    raise NameError('%s.%s exists' % (enum.name, name))
                enum[name] = self


class Ellipsoid(_Based):
    '''Ellipsoid with semi-major, semi-minor axis, inverse flattening
       and a number of other pre-computed, frequently used values.
    '''
    a    = 0  #: Semi-major, equatorial axis (meter).
    b    = 0  #: Semi-minor, polar axis (meter): a * (f - 1) / f.
    # pre-computed, frequently used values
    a2   = 0  #: (1 / a**2) (float).
    ab   = 1  #: (a / b) = 1 / (1 - f) (float).
    e    = 0  #: 1st Eccentricity: sqrt(1 - (b / a)**2)) (float).
    e2   = 0  #: 1st Eccentricity squared: f * (2 - f) = (a**2 - b**2) / a**2 (float).
    e4   = 0  #: e2**2 (float).
    e12  = 1  #: (1 - e2) (float).
    e22  = 0  #: 2nd Eccentricity squared: e2 / (1 - e2) = ab**2 - 1 (float).
    f    = 0  #: Flattening: (a - b) / a (float).
    f_   = 0  #: Inverse flattening: a / (a - b) = 1 /f (float).
    n    = 0  #: 3rd Flattening: f / (2 - f) = (a - b) / (a + b) (float).
    # radii from <http://wikipedia.org/wiki/Earth_radius>
    R    = 0  #: Mean radius: (2 * a + b) / 3 per IUGG definition (meter).
    Rm   = 0  #: Mean radius: sqrt(a * b) (meter).
    R2   = 0  #: Authalic radius: sqrt((a**2 + b**2 * atanh(e) / e) / 2) (meter).
    R3   = 0  #: Volumetric radius: cbrt(a * a * b) (meter).
    Rr   = 0  #: Rectifying radius: ((a**3/2 + b**3/2) / 2)**2/3 (meter).

    _A      = None  #: (INTERNAL) meridian radius
    _Alpha6 = None  #: (INTERNAL) 6th-order Krüger Alpha series
    _Beta6  = None  #: (INTERNAL) 6th-order Krüger Beta series
    _Mabcd  = None  #: (INTERNAL) OSGB meridional coefficients

    def __init__(self, a, b, f_, name=''):
        '''New ellipsoid.

           @param a: Semi-major, equatorial axis (meter).
           @param b: Semi-minor, polar axis (meter).
           @param f_: Inverse flattening: a / (a - b) (float >>> 1).
           @keyword name: Optional, unique name (string).

           @raise NameError: If ellipsoid name already exists.
        '''
        self.a = a = float(a)  # major half-axis in meter
        if not b:  # get b from a and f_
            self.b = b = a * (f_ - 1) / float(f_)
        else:  # get f_ from a and b if not spherical
            self.b = b = float(b)  # minor half-axis in meter
            if not f_ and a > b:
                f_ = a / (a - b)
        if f_ > 0 and a > b:
            self.f_ = f_ = float(f_)  # inverse flattening
            self.f  = f  = 1 / f_  # flattening
            self.n  = n  = f / (2 - f)  # 3rd flattening for utm
            self.e2 = e2 = f * (2 - f)  # 1st eccentricity squared
            self.e4 = e2 * e2  # for Nvector.Cartesian.toNvector
            self.e  = sqrt(e2)  # eccentricity for utm
            self.e12 = 1 - e2  # for Nvector.Cartesian.toNvector and utm
            self.e22 = e2 / (1 - e2)  # 2nd eccentricity squared
            self.ab = a / b  # for Nvector.toCartesian
            self.R = (2 * a + b) / 3  # per IUGG definition for WGS84
            self.Rm = sqrt(a * b)  # mean radius
            self.R2 = sqrt((a * a + b * b * atanh(self.e) / self.e) * 0.5)  # authalic radius
            self.R3 = cbrt(a * a * b)  # volumetric radius
            self.Rr = cbrt2((pow(a, 1.5) + pow(b, 1.5)) * 0.5)  # rectifying radius
        else:
            self.R = self.Rm = self.R2 = self.R3 = self.Rr = self.b = b = a
            f_ = f = n = 0
        self.a2 = 1 / (a * a)  # for Nvector.Cartesian.toNvector

        d = a - b
        self._ab_90 = d / 90  # for radiusAt below

        # some sanity checks to catch mistakes
        if d < 0 or min(a, b) < 1:
            raise AssertionError('%s: %s=%0.9f vs %s=%0.9f' % (name,
                                 'a', a, 'b', b))
        t = d / a
        if abs(f - t) > 1e-8:
            raise AssertionError('%s: %s=%.9e vs %s=%.9e' % (name,
                                 'f', f, '(a-b)/a', t))
        t = d / (a + b)
        if abs(n - t) > 1e-8:
            raise AssertionError('%s: %s=%.9e vs %s=%.9e' % (name,
                                 'n', n, '(a-b)/(a+b)', t))
        t = self.ab ** 2 - 1
        if abs(self.e22 - t) > 1e-8:
            raise AssertionError('%s: %s=%.9e vs %s=%.9e' % (name,
                                 'e22', self.e22, 'ab**2-1', t))

        self._register(Ellipsoids, name)

    def __eq__(self, other):
        '''Compare this and an other ellipsoid.

           @param other: The other ellipsoid (L{Ellipsoid}).

           @return: True if equal (bool).
        '''
        return self is other or (isinstance(other, Ellipsoid) and
                                 self.a == other.a and
                                 self.b == other.b)

    @property
    def A(self):
        '''Get the meridional radius (meter).
        '''
        if self._A is None:
            n = self.n
            n2_4 = n * n / 4
            n4_64 = n2_4 * n2_4 / 4
            # A = a / (1 + n) * (1 + n**2 / 4 + n**4 / 64 + n**6 / 256 +
            #                   n**8 * 25 / 16384 + n**10 * 49 / 65536)
            self._A = self.a / (1 + n) * (1 + n2_4) * (1 + n4_64)
        return self._A

    @property
    def Alpha6(self):
        '''Get the 6th-order Krüger Alpha series (7-tuple, 1-origin).
        '''
        if self._Alpha6 is None:
            self._Alpha6 = self._K6(
                # XXX i/i quotients require  from __future__ import division
                # n    n**2   n**3      n**4         n**5            n**6
                (1/2, -2/3,   5/16,    41/180,    -127/288,       7891/37800),
                     (13/48, -3/5,    557/1440,    281/630,   -1983433/1935360),  # PYCHOK expected
                            (61/240, -103/140,   15061/26880,   167603/181440),  # PYCHOK expected
                                   (49561/161280, -179/168,    6601661/7257600),  # PYCHOK expected
                                                (34729/80640, -3418889/1995840),  # PYCHOK expected
                                                            (212378941/319334400,))  # PYCHOK expected
        return self._Alpha6

    @property
    def Beta6(self):
        '''Get the 6th-order Krüger Beta series (7-tuple, 1-origin).
        '''
        if self._Beta6 is None:
            self._Beta6 = self._K6(
                # XXX i/i quotients require  from __future__ import division
                # n    n**2  n**3     n**4        n**5            n**6
                (1/2, -2/3, 37/96,   -1/360,    -81/512,      96199/604800),
                      (1/48, 1/15, -437/1440,    46/105,   -1118711/3870720),  # PYCHOK expected
                           (17/480, -37/840,   -209/4480,      5569/90720),  # PYCHOK expected
                                  (4397/161280, -11/504,    -830251/7257600),  # PYCHOK expected
                                              (4583/161280, -108847/3991680),  # PYCHOK expected
                                                          (20648693/638668800,))  # PYCHOK expected
        return self._Beta6

    def e2s2(self, s):
        '''Compute the norm sqrt(1 - e2 * s**2).

           @param s: S value (scalar).

           @return: Norm (float).
        '''
        return sqrt(1 - self.e2 * s * s)

    @property
    def isEllipsoidal(self):
        '''Check whether this model is ellipsoidal (bool).
        '''
        return self.a > self.R > self.b

    @property
    def isSpherical(self):
        '''Check whether this model is spherical (bool).
        '''
        return self.a == self.R == self.b

    def _K6(self, *fs6):
        '''(INTERNAL) Compute the 6th-order Krüger Alpha or Beta
           series per Karney 2011, 'Transverse Mercator with an
           accuracy of a few nanometers', page 7, equations 35
           and 36, see <http://arxiv.org/pdf/1002.1417v3.pdf>.

           @param fs6: 6-Tuple of coefficent tuples.

           @return: 6th-Order Krüger (7-tuple, 1-origin).
        '''
        ns = [self.n]  # 3rd flattening: n, n**2, ... n**6
        for i in range(len(fs6) - 1):
            ns.append(ns[0] * ns[i])

        k6 = [0]  # 1-origin
        for fs in fs6:
            i = len(ns) - len(fs)
            k6.append(fdot(fs, *ns[i:]))

        return tuple(k6)

    @property
    def Mabcd(self):
        '''Get the OSGR meridional coefficients, Airy130 only (4-tuple).
        '''
        if self._Mabcd is None:
            n = self.n
            n2 = n * n
            n3 = n * n2
            # XXX i/i quotients require  from __future__ import division
            self._Mabcd = (fdot((1, n, n2, n3), 1, 1,  5/4, 5/4),
                           fdot(   (n, n2, n3), 3, 3, 21/8),
                           fdot(      (n2, n3), 15/8, 15/8),
                                   35/24 * n3)
        return self._Mabcd

    def radiusAt(self, lat):
        '''Approximate the ellipsoid radius at the given
           latitude in degrees by trivial interpolation.

           @param lat: Latitude (degrees90).

           @return: Radius at that latitude (meter).
        '''
        # r = major - (major - minor) * |lat| / 90
        return self.a - self._ab_90 * min(abs(lat), 90)

    def toStr(self, prec=8):  # PYCHOK expected
        '''Return this ellipsoid as a string.

           @keyword prec: Number of decimals, unstripped (int).

           @return: Ellipsoid attributes (string).
        '''
        return self._fStr(prec, 'a', 'b', 'f_', 'f', 'e2', 'e22',
                                'R', 'Rm', 'R2', 'R3', 'Rr')


# <http://www.gnu.org/software/gama/manual/html_node/Supported-ellipsoids.html>
Ellipsoids._assert(  # <http://wikipedia.org/wiki/Earth_ellipsoid>
    Airy1830       = Ellipsoid(6377563.396, 6356256.909,       299.3249646,   'Airy1830'),
    AiryModified   = Ellipsoid(6377340.189, 6356034.448,       299.3249646,   'AiryModified'),
    Australia1966  = Ellipsoid(6378160.0,   6356774.719,       298.25,        'Australia1966'),
    Bessel1841     = Ellipsoid(6377397.155, 6356078.963,       299.152815351, 'Bessel1841'),  # XXX 299.1528128
    Clarke1866     = Ellipsoid(6378206.4,   6356583.8,         294.978698214, 'Clarke1866'),
    Clarke1880IGN  = Ellipsoid(6378249.2,   6356515.0,         293.466021294, 'Clarke1880IGN'),  # XXX confirm
    CPM1799        = Ellipsoid(6375738.7,   6356671.92557493,  334.39,        'CPM1799'),  # Comm. des Poids et Mesures
    Delambre1810   = Ellipsoid(6376428.0,   6355957.92616372,  311.5,         'Delambre1810'),  # Belgium
    Engelis1985    = Ellipsoid(6378136.05,  6356751.32272154,  298.2566,      'Engelis1985'),
    Everest1969    = Ellipsoid(6377295.664, 6356094.667915,    300.8017,      'Everest1969'),
    Fisher1968     = Ellipsoid(6378150.0,   6356768.33724438,  298.3,         'Fisher1968'),
    GRS67          = Ellipsoid(6378160.0,   6356774.516,       298.247167427, 'GRS67'),  # Lucerne
    GRS80          = Ellipsoid(6378137.0,   6356752.314140347, 298.257222101, 'GRS80'),  # ITRS, ETRS89
    Helmert1906    = Ellipsoid(6378200.0,   6356818.16962789,  298.3,         'Helmert1906'),
    IERS1989       = Ellipsoid(6378136.0,   6356751.302,       298.257,       'IERS1989'),
    IERS2003       = Ellipsoid(6378136.6,   6356751.85797165,  298.25642,     'IERS2003'),
    Intl1924       = Ellipsoid(6378388.0,   6356911.946,       297.0,         'Intl1924'),  # aka Hayford
    Intl1967       = Ellipsoid(6378157.5,   6356772.2,         298.24961539,  'Intl1967'),  # New Int'l
    Krassovsky1940 = Ellipsoid(6378245.0,   6356863.019,       298.3,         'Krassovsky1940'),
    Maupertuis1738 = Ellipsoid(6397300.0,   6363806.28272251,  191.0,         'Maupertuis1738'),  # France
    NWL1965        = Ellipsoid(6378145.0,   6356759.76948868,  298.25,        'NWL1965'),  # Naval Weapons Lab.
    Plessis1817    = Ellipsoid(6397523.0,   6355863.0,         153.56512242,  'Plessis1817'),  # France
    SGS85          = Ellipsoid(6378136.0,   6356751.30156878,  298.257,       'SGS85'),  # Soviet Geodetic System
    WGS60          = Ellipsoid(6378165.0,   6356783.28695944,  298.3,         'WGS60'),
    WGS66          = Ellipsoid(6378145.0,   6356759.76948868,  298.25,        'WGS66'),
    WGS72          = Ellipsoid(6378135.0,   6356750.52,        298.26,        'WGS72'),
    WGS84          = Ellipsoid(6378137.0,   6356752.31425,     298.257223563, 'WGS84'),  # GPS
    Sphere         = Ellipsoid(R_M,         R_M,                 0.0,         'Sphere'),  # pseudo
)


class Transform(_Based):
    '''Helmert transformation.
    '''
    tx = 0  #: X translation (meter).
    ty = 0  #: Y translation (meter).
    tz = 0  #: Z translation (meter).

    rx = 0  #: X rotation (radians).
    ry = 0  #: Y rotation (radians).
    rz = 0  #: Z rotation (radians).

    s  = 0  #: Scale ppm (float).
    s1 = 1  #: Scale + 1 (float).

    sx = 0  #: X rotation (degree seconds).
    sy = 0  #: Y rotation (degree seconds).
    sz = 0  #: Z rotation (degree seconds).

    def __init__(self, name='', tx=0, ty=0, tz=0,
                                sx=0, sy=0, sz=0, s=0):
        '''New transform.

           @keyword name: Optional, unique name (string).
           @keyword tx: X translation (meter).
           @keyword ty: Y translation (meter).
           @keyword tz: Z translation (meter).
           @keyword s: Scale ppm (float).
           @keyword sx: X rotation (degree seconds).
           @keyword sy: Y rotation (degree seconds).
           @keyword sz: Z rotation (degree seconds).

           @raise NameError: If transform name already exists.
        '''
        if tx:
            self.tx = float(tx)
        if ty:
            self.ty = float(ty)
        if tz:
            self.tz = float(tz)
        if sx:  # secs to rads
            self.rx = radians(sx / 3600.0)
            self.sx = sx
        if sy:
            self.ry = radians(sy / 3600.0)
            self.sy = sy
        if sz:
            self.rz = radians(sz / 3600.0)
            self.sz = sz
        if s:
            self.s = float(s)
            self.s1 = s * 1.e-6 + 1  # normalize ppm to (s + 1)

        self._register(Transforms, name)

    def __eq__(self, other):
        '''Compare this and an other transform.

           @param other: The other transform (L{Transform}).

           @return: True if equal (bool).
        '''
        return self is other or (isinstance(other, Transform) and
                                 self.tx == other.tx and
                                 self.ty == other.ty and
                                 self.tz == other.tz and
                                 self.rx == other.rx and
                                 self.ry == other.ry and
                                 self.rz == other.rz and
                                 self.s  == other.s)

    def inverse(self, name=''):
        '''Return the inverse of this transform.

           @keyword name: Optional, unique name (string).

           @return: Inverse (Transform).

           @raise NameError: If transform name already exists.
        '''
        return Transform(name=name or 'Inverse' + self.name,
                         tx=-self.tx, ty=-self.ty, tz=-self.tz,
                         sx=-self.sx, sy=-self.sy, sz=-self.sz, s=-self.s)

    def toStr(self, prec=4):  # PYCHOK expected
        '''Return this transform as a string.

           @keyword prec: Number of decimals, unstripped (int).

           @return: Transform attributes (string).
        '''
        return self._fStr(prec, 'tx', 'ty', 'tz',
                                'rx', 'ry', 'rz', 's', 's1',
                                'sx', 'sy', 'sz')

    def transform(self, x, y, z, inverse=False):
        '''Transform a (geocentric) Cartesian point, forward or inverse.

           @param x: X coordinate (meter).
           @param y: Y coordinate (meter).
           @param z: Z coordinate (meter).
           @keyword inverse: Direction, forward or inverse (bool).

           @return: 3-Tuple (x, y, z) transformed.
        '''
        if inverse:
            xyz = -1, -x, -y, -z
            _s1 = self.s1 - 2  # negative inverse: -(1 - s * 1.e-6)
        else:
            xyz =  1,  x,  y,  z
            _s1 = self.s1
        # x', y', z' = (.tx + x * .s1 - y * .rz + z * .ry,
        #               .ty + x * .rz + y * .s1 - z * .rx,
        #               .tz - x * .ry + y * .rx + z * .s1)
        return (fdot(xyz, self.tx,      _s1, -self.rz,  self.ry),
                fdot(xyz, self.ty,  self.rz,      _s1, -self.rx),
                fdot(xyz, self.tz, -self.ry,  self.rx,      _s1))


# <http://wikipedia.org/wiki/Helmert_transformation> from WGS84
Transforms._assert(
    BD72           = Transform('BD72', tx=106.868628, ty=-52.297783, tz=103.723893,
                     # <http://www.ngi.be/FR/FR4-4.shtm> ETRS89 == WG84
                     # <http://georepository.com/transformation_15929/BD72-to-WGS-84-3.html>
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
                     # <http://geonet.esri.com/thread/36583> sz=-0.156
                     # <http://github.com/chrisveness/geodesy/blob/master/latlon-ellipsoidal.js>
                     # <http://www.gov.uk/guidance/oil-and-gas-petroleum-operations-notices#pon-4>
                                                         sz=  0.156, s=-1.2),
    Irl1965        = Transform('Irl1965', tx=-482.530, ty=130.596, tz=-564.557,
                                          sx=   1.042, sy=  0.214, sz=   0.631,
                                           s=  -8.15),
    Irl1975        = Transform('Irl1975', tx=-482.530, ty=130.596, tz=-564.557,
                     # XXX rotation signs may be opposite, to be checked
                                          sx=  -1.042, sy= -0.214, sz=  -0.631,
                                           s=  -1.1),
    Krassovsky1940 = Transform('Krassovsky1940', tx=-24.0,  ty=123.0,  tz=94.0,
                                                 sx= -0.02, sy=  0.26, sz= 0.13,
                                                  s= -2.423),
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


class Datum(_Based):
    '''Ellipsoid and transform parameters for an earth model.
    '''
    _ellipsoid = Ellipsoids.WGS84  #: (INTERNAL) Default ellipsoid (L{Ellipsoid}).
    _transform = Transforms.WGS84  #: (INTERNAL) Default transform (L{Transform}).

    def __init__(self, ellipsoid, transform=None, name=''):
        '''New datum.

           @param ellipsoid: The ellipsoid (L{Ellipsoid}).
           @keyword transform: The transform (L{Transform}).
           @keyword name: Optional, unique name (string).

           @raise NameError: If datum name already exists.

           @raise TypeError: If ellipsoid is not an L{Ellipsoid}
                             or transform is not a L{Transform}.
        '''
        self._ellipsoid = ellipsoid or Datum._ellipsoid
        if not isinstance(self.ellipsoid, Ellipsoid):
            raise TypeError('%s not an %s: %r' % ('ellipsoid', Ellipsoid.__name__, self.ellipsoid))

        self._transform = transform or Datum._transform
        if not isinstance(self.transform, Transform):
            raise TypeError('%s not a %s: %r' % ('transform', Transform.__name__, self.transform))

        self._register(Datums, name or self.transform.name or self.ellipsoid.name)

    def __eq__(self, other):
        '''Compare this and an other datum.

           @param other: The other datum (L{Datum}).

           @return: True if equal (bool)
        '''
        return self is other or (isinstance(other, Datum) and
                                 self.ellipsoid == other.ellipsoid and
                                 self.transform == other.transform)

    @property
    def ellipsoid(self):
        '''Get this datum's ellipsoid (L{Ellipsoid}).
        '''
        return self._ellipsoid

    @property
    def isEllipsoidal(self):
        '''Check whether this datum is ellipsoidal (bool).
        '''
        return self._ellipsoid.isEllipsoidal

    @property
    def isSpherical(self):
        '''Check whether this datum is spherical (bool).
        '''
        return self._ellipsoid.isSpherical

    def toStr(self, **unused):  # PYCHOK expected
        '''Return this datum as a string.

           @return: Datum attributes (string).
        '''
        t = []
        for a in ('ellipsoid', 'transform'):
            v = getattr(self, a)
            t.append('%s=%ss.%s' % (a, v.__class__.__name__, v.name))
        return ', '.join(['name=%r' % (self.name,)] + t)

    @property
    def transform(self):
        '''Get this datum's transform (L{Transform}).
        '''
        return self._transform


# Datums with associated ellipsoid and Helmert transform parameters
# to convert from WGS84 into the given datum.  More are available at
# <http://earth-info.nga.mil/GandG/coordsys/datums/NATO_DT.pdf> and
# <http://www.fieldenmaps.info/cconv/web/cconv_params.js>.
Datums._assert(
    # Belgian Datum 1972, based on Hayford ellipsoid.
    # <http://nl.wikipedia.org/wiki/Belgian_Datum_1972>
    # <http://spatialreference.org/ref/sr-org/belge-1972-belgian-
    #         lambert-72-corrected-transformation-parameters/>
    BD72           = Datum(Ellipsoids.Intl1924, Transforms.BD72),
    # Germany <http://wikipedia.org/wiki/Bessel-Ellipsoid>
    #         <http://wikipedia.org/wiki/Helmert_transformation>
    DHDN           = Datum(Ellipsoids.Bessel1841, Transforms.DHDN),

    # <http://www.gov.uk/guidance/oil-and-gas-petroleum-operations-notices#pon-4>
    ED50           = Datum(Ellipsoids.Intl1924, Transforms.ED50),

    # <http://wikipedia.org/wiki/GRS_80>
    GRS80          = Datum(Ellipsoids.GRS80, Transforms.WGS84, name='GRS80'),

    # <http://osi.ie/OSI/media/OSI/Content/Publications/transformations_booklet.pdf>
    Irl1975        = Datum(Ellipsoids.AiryModified, Transforms.Irl1975),

    # Germany <http://wikipedia.org/wiki/Helmert_transformation>
    Krassovsky1940 = Datum(Ellipsoids.Krassovsky1940, Transforms.Krassovsky1940),

    # Austria <http://de.wikipedia.org/wiki/Datum_Austria>
    MGI            = Datum(Ellipsoids.Bessel1841, Transforms.MGI),

    # <http://wikipedia.org/wiki/Helmert_transformation>
    NAD27          = Datum(Ellipsoids.Clarke1866, Transforms.NAD27),

    # NAD83 (2009) == WGS84 - <http://www.uvm.edu/giv/resources/WGS84_NAD83.pdf>
    # (If you *really* must convert WGS84<->NAD83, you need more than this!)
    NAD83          = Datum(Ellipsoids.GRS80, Transforms.NAD83),

    #  Nouvelle Triangulation Francaise (Paris)  XXX verify
    NTF            = Datum(Ellipsoids.Clarke1880IGN, Transforms.NTF),

    # <http://www.ordnancesurvey.co.uk/docs/support/guide-coordinate-systems-great-britain.pdf>
    OSGB36         = Datum(Ellipsoids.Airy1830, Transforms.OSGB36),

    # Germany <http://wikipedia.org/wiki/Helmert_transformation>
    Potsdam        = Datum(Ellipsoids.Bessel1841, Transforms.Bessel1841, name='Potsdam'),

    # XXX psuedo-ellipsoids for spherical LatLon
    Sphere         = Datum(Ellipsoids.Sphere, Transforms.WGS84, name='Sphere'),

    # <http://www.geocachingtoolbox.com?page=datumEllipsoidDetails>
    TokyoJapan     = Datum(Ellipsoids.Bessel1841, Transforms.TokyoJapan),

    # <http://www.icao.int/safety/pbn/documentation/eurocontrol/eurocontrol%20wgs%2084%20implementation%20manual.pdf>
    WGS72          = Datum(Ellipsoids.WGS72, Transforms.WGS72),

    WGS84          = Datum(Ellipsoids.WGS84, Transforms.WGS84),
)


if __name__ == '__main__':

    # print all
    for e in (Ellipsoids, Transforms, Datums):
        print('\n%r' % (e,))

    E = Datums.WGS84.ellipsoid
    e = (E.a - E.b) / (E.a + E.b) - E.n
    t = (E.toStr(prec=10),
        'A=%r, e=%.10f, f=1/%.10f, n=%.10f(%.10e)' % (E.A, E.e, 1/E.f, E.n, e),
        'Alpha6=%r' % (E.Alpha6,),
        'Beta6=%r' % (E.Beta6,))
    print('\n%s: %s' % (E.name, ',\n       '.join(t)))

# **) MIT License
#
# Copyright (C) 2016-2017 -- mrJean1 at Gmail dot com
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

# Typical result (on macOS 10.12.3, 10,12,4 and 10.12.5 Sierra)

# Ellipsoids.Airy1830: Ellipsoid(name='Airy1830', a=6377563.396, b=6356256.909, f_=299.3249646, f=0.00334085, e2=0.00667054, e22=0.00671533, R=6370461.23366667, Rm=6366901.23988196, R2=6370459.65458944, R3=6370453.30986645, Rr=6366914.60880589)
# Ellipsoids.AiryModified: Ellipsoid(name='AiryModified', a=6377340.189, b=6356034.448, f_=299.3249646, f=0.00334085, e2=0.00667054, e22=0.00671533, R=6370238.27533333, Rm=6366678.40619415, R2=6370236.69636116, R3=6370230.35181066, Rr=6366691.7746498)
# Ellipsoids.Australia1966: Ellipsoid(name='Australia1966', a=6378160.0, b=6356774.719, f_=298.25, f=0.00335289, e2=0.00669454, e22=0.00673966, R=6371031.573, Rm=6367458.38162583, R2=6371029.98238815, R3=6371023.59117818, Rr=6367471.84843391)
# Ellipsoids.Bessel1841: Ellipsoid(name='Bessel1841', a=6377397.155, b=6356078.963, f_=299.15281535, f=0.00334277, e2=0.00667437, e22=0.00671922, R=6370291.091, Rm=6366729.13634557, R2=6370289.51018729, R3=6370283.15827603, Rr=6366742.52032409)
# Ellipsoids.CPM1799: Ellipsoid(name='CPM1799', a=6375738.7, b=6356671.92557493, f_=334.39, f=0.00299052, e2=0.0059721, e22=0.00600798, R=6369383.10852498, Rm=6366198.17466371, R2=6369381.8434158, R3=6369376.76247021, Rr=6366208.88184734)
# Ellipsoids.Clarke1866: Ellipsoid(name='Clarke1866', a=6378206.4, b=6356583.8, f_=294.97869821, f=0.00339008, e2=0.00676866, e22=0.00681478, R=6370998.86666667, Rm=6367385.92165547, R2=6370997.240633, R3=6370990.70659881, Rr=6367399.68916895)
# Ellipsoids.Clarke1880IGN: Ellipsoid(name='Clarke1880IGN', a=6378249.2, b=6356515.0, f_=293.46602129, f=0.00340755, e2=0.00680349, e22=0.00685009, R=6371004.46666667, Rm=6367372.82664821, R2=6371002.82383111, R3=6370996.22212394, Rr=6367386.73667251)
# Ellipsoids.Delambre1810: Ellipsoid(name='Delambre1810', a=6376428.0, b=6355957.92616372, f_=311.5, f=0.00321027, e2=0.00641024, e22=0.0064516, R=6369604.64205457, Rm=6366184.7355549, R2=6369603.18419749, R3=6369597.32739068, Rr=6366197.07684267)
# Ellipsoids.Engelis1985: Ellipsoid(name='Engelis1985', a=6378136.05, b=6356751.32272154, f_=298.2566, f=0.00335282, e2=0.00669439, e22=0.00673951, R=6371007.80757385, Rm=6367434.70891814, R2=6371006.21707085, R3=6370999.82613572, Rr=6367448.17507892)
# Ellipsoids.Everest1969: Ellipsoid(name='Everest1969', a=6377295.664, b=6356094.667915, f_=300.8017, f=0.00332445, e2=0.00663785, e22=0.0066822, R=6370228.665305, Rm=6366686.3410779, R2=6370227.10178534, R3=6370220.81951617, Rr=6366699.57839424)
# Ellipsoids.Fisher1968: Ellipsoid(name='Fisher1968', a=6378150.0, b=6356768.33724438, f_=298.3, f=0.00335233, e2=0.00669342, e22=0.00673853, R=6371022.77908146, Rm=6367450.19377421, R2=6371021.18903735, R3=6371014.79995034, Rr=6367463.65604301)
# Ellipsoids.GRS67: Ellipsoid(name='GRS67', a=6378160.0, b=6356774.516, f_=298.24716743, f=0.00335292, e2=0.00669461, e22=0.00673973, R=6371031.50533333, Rm=6367458.27995524, R2=6371029.91470873, R3=6371023.52335984, Rr=6367471.74701921)
# Ellipsoids.GRS80: Ellipsoid(name='GRS80', a=6378137.0, b=6356752.31414035, f_=298.2572221, f=0.00335281, e2=0.00669438, e22=0.0067395, R=6371008.77138012, Rm=6367435.67966369, R2=6371007.18088351, R3=6371000.78997413, Rr=6367449.14577025)
# Ellipsoids.Helmert1906: Ellipsoid(name='Helmert1906', a=6378200.0, b=6356818.16962789, f_=298.3, f=0.00335233, e2=0.00669342, e22=0.00673853, R=6371072.7232093, Rm=6367500.10989561, R2=6371071.13315272, R3=6371064.74401563, Rr=6367513.57226994)
# Ellipsoids.IERS1989: Ellipsoid(name='IERS1989', a=6378136.0, b=6356751.302, f_=298.257, f=0.00335281, e2=0.00669438, e22=0.0067395, R=6371007.76733333, Rm=6367434.6735819, R2=6371006.17690648, R3=6370999.78591702, Rr=6367448.13970588)
# Ellipsoids.IERS2003: Ellipsoid(name='IERS2003', a=6378136.6, b=6356751.85797165, f_=298.25642, f=0.00335282, e2=0.0066944, e22=0.00673951, R=6371008.35265722, Rm=6367435.25153158, R2=6371006.76215217, R3=6371000.37120876, Rr=6367448.71770978)
# Ellipsoids.Intl1924: Ellipsoid(name='Intl1924', a=6378388.0, b=6356911.946, f_=297.0, f=0.003367, e2=0.00672267, e22=0.00676817, R=6371229.31533333, Rm=6367640.91900784, R2=6371227.71127046, R3=6371221.26583212, Rr=6367654.49999285)
# Ellipsoids.Intl1967: Ellipsoid(name='Intl1967', a=6378157.5, b=6356772.2, f_=298.24961539, f=0.0033529, e2=0.00669455, e22=0.00673967, R=6371029.06666667, Rm=6367455.87210634, R2=6371027.4760839, R3=6371021.08482752, Rr=6367469.33894366)
# Ellipsoids.Krassovsky1940: Ellipsoid(name='Krassovsky1940', a=6378245.0, b=6356863.019, f_=298.3, f=0.00335233, e2=0.00669342, e22=0.00673853, R=6371117.673, Rm=6367545.03451854, R2=6371116.08297003, R3=6371109.69375021, Rr=6367558.49698756)
# Ellipsoids.Maupertuis1738: Ellipsoid(name='Maupertuis1738', a=6397300.0, b=6363806.28272251, f_=191.0, f=0.0052356, e2=0.01044379, e22=0.01055402, R=6386135.42757417, Rm=6380531.16381863, R2=6386131.54144846, R3=6386115.88628229, Rr=6380564.13011364)
# Ellipsoids.NWL1965: Ellipsoid(name='NWL1965', a=6378145.0, b=6356759.76948868, f_=298.25, f=0.00335289, e2=0.00669454, e22=0.00673966, R=6371016.58982956, Rm=6367443.40689145, R2=6371014.999254, R3=6371008.60802666, Rr=6367456.87366762)
# Ellipsoids.Plessis1817: Ellipsoid(name='Plessis1817', a=6397523.0, b=6355863.0, f_=153.56512242, f=0.0065119, e2=0.01298139, e22=0.01315212, R=6383636.33333333, Rm=6376658.97844232, R2=6383630.32549925, R3=6383606.08096947, Rr=6376710.01073346)
# Ellipsoids.SGS85: Ellipsoid(name='SGS85', a=6378136.0, b=6356751.30156878, f_=298.257, f=0.00335281, e2=0.00669438, e22=0.0067395, R=6371007.76718959, Rm=6367434.67336593, R2=6371006.17669087, R3=6370999.78577296, Rr=6367448.13949045)
# Ellipsoids.Sphere: Ellipsoid(name='Sphere', a=6371008.771415, b=6371008.771415, f_=0, f=0, e2=0, e22=0, R=6371008.771415, Rm=6371008.771415, R2=6371008.771415, R3=6371008.771415, Rr=6371008.771415)
# Ellipsoids.WGS60: Ellipsoid(name='WGS60', a=6378165.0, b=6356783.28695944, f_=298.3, f=0.00335233, e2=0.00669342, e22=0.00673853, R=6371037.76231981, Rm=6367465.16861063, R2=6371036.17227197, R3=6371029.78316993, Rr=6367478.6309111)
# Ellipsoids.WGS66: Ellipsoid(name='WGS66', a=6378145.0, b=6356759.76948868, f_=298.25, f=0.00335289, e2=0.00669454, e22=0.00673966, R=6371016.58982956, Rm=6367443.40689145, R2=6371014.999254, R3=6371008.60802666, Rr=6367456.87366762)
# Ellipsoids.WGS72: Ellipsoid(name='WGS72', a=6378135.0, b=6356750.52, f_=298.26, f=0.00335278, e2=0.00669432, e22=0.00673943, R=6371006.84, Rm=6367433.78276368, R2=6371005.24953082, R3=6370998.85874532, Rr=6367447.24861499)
# Ellipsoids.WGS84: Ellipsoid(name='WGS84', a=6378137.0, b=6356752.31425, f_=298.25722356, f=0.00335281, e2=0.00669438, e22=0.0067395, R=6371008.77141667, Rm=6367435.67971861, R2=6371007.18092088, R3=6371000.79001076, Rr=6367449.14582503)

# Transforms.BD72: Transform(name='BD72', tx=106.8686, ty=-52.2978, tz=103.7239, rx=-0.0, ry=-0.0, rz=-0.0, s=1.2727, s1=1.0, sx=-0.3366, sy=-0.457, sz=-1.8422)
# Transforms.Bessel1841: Transform(name='Bessel1841', tx=-582.0, ty=-105.0, tz=-414.0, rx=-0.0, ry=-0.0, rz=0.0, s=-8.3, s1=1.0, sx=-1.04, sy=-0.35, sz=3.08)
# Transforms.Clarke1866: Transform(name='Clarke1866', tx=8.0, ty=-160.0, tz=-176.0, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
# Transforms.DHDN: Transform(name='DHDN', tx=-591.28, ty=-81.35, tz=-396.39, rx=0.0, ry=-0.0, rz=-0.0, s=-9.82, s1=1.0, sx=1.477, sy=-0.0736, sz=-1.458)
# Transforms.ED50: Transform(name='ED50', tx=89.5, ty=93.8, tz=123.1, rx=0, ry=0, rz=0.0, s=-1.2, s1=1.0, sx=0, sy=0, sz=0.156)
# Transforms.Irl1965: Transform(name='Irl1965', tx=-482.53, ty=130.596, tz=-564.557, rx=0.0, ry=0.0, rz=0.0, s=-8.15, s1=1.0, sx=1.042, sy=0.214, sz=0.631)
# Transforms.Irl1975: Transform(name='Irl1975', tx=-482.53, ty=130.596, tz=-564.557, rx=-0.0, ry=-0.0, rz=-0.0, s=-1.1, s1=1.0, sx=-1.042, sy=-0.214, sz=-0.631)
# Transforms.Krassovsky1940: Transform(name='Krassovsky1940', tx=-24.0, ty=123.0, tz=94.0, rx=-0.0, ry=0.0, rz=0.0, s=-2.423, s1=1.0, sx=-0.02, sy=0.26, sz=0.13)
# Transforms.MGI: Transform(name='MGI', tx=-577.326, ty=-90.129, tz=-463.92, rx=0.0, ry=0.0, rz=0.0, s=-2.423, s1=1.0, sx=5.137, sy=1.474, sz=5.297)
# Transforms.NAD27: Transform(name='NAD27', tx=8.0, ty=-160.0, tz=-176.0, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
# Transforms.NAD83: Transform(name='NAD83', tx=1.004, ty=-1.91, tz=-0.515, rx=0.0, ry=0.0, rz=0.0, s=-0.0015, s1=1.0, sx=0.0267, sy=0.0003, sz=0.011)
# Transforms.NTF: Transform(name='NTF', tx=-168.0, ty=-60.0, tz=320.0, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
# Transforms.OSGB36: Transform(name='OSGB36', tx=-446.448, ty=125.157, tz=-542.06, rx=-0.0, ry=-0.0, rz=-0.0, s=20.4894, s1=1.0, sx=-0.1502, sy=-0.247, sz=-0.8421)
# Transforms.TokyoJapan: Transform(name='TokyoJapan', tx=148.0, ty=-507.0, tz=-685.0, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)
# Transforms.WGS72: Transform(name='WGS72', tx=0, ty=0, tz=-4.5, rx=0, ry=0, rz=0.0, s=-0.22, s1=1.0, sx=0, sy=0, sz=0.554)
# Transforms.WGS84: Transform(name='WGS84', tx=0, ty=0, tz=0, rx=0, ry=0, rz=0, s=0, s1=1, sx=0, sy=0, sz=0)

# Datums.BD72: Datum(name='BD72', ellipsoid=Ellipsoids.Intl1924, transform=Transforms.BD72)
# Datums.DHDN: Datum(name='DHDN', ellipsoid=Ellipsoids.Bessel1841, transform=Transforms.DHDN)
# Datums.ED50: Datum(name='ED50', ellipsoid=Ellipsoids.Intl1924, transform=Transforms.ED50)
# Datums.GRS80: Datum(name='GRS80', ellipsoid=Ellipsoids.GRS80, transform=Transforms.WGS84)
# Datums.Irl1975: Datum(name='Irl1975', ellipsoid=Ellipsoids.AiryModified, transform=Transforms.Irl1975)
# Datums.Krassovsky1940: Datum(name='Krassovsky1940', ellipsoid=Ellipsoids.Krassovsky1940, transform=Transforms.Krassovsky1940)
# Datums.MGI: Datum(name='MGI', ellipsoid=Ellipsoids.Bessel1841, transform=Transforms.MGI)
# Datums.NAD27: Datum(name='NAD27', ellipsoid=Ellipsoids.Clarke1866, transform=Transforms.NAD27)
# Datums.NAD83: Datum(name='NAD83', ellipsoid=Ellipsoids.GRS80, transform=Transforms.NAD83)
# Datums.NTF: Datum(name='NTF', ellipsoid=Ellipsoids.Clarke1880IGN, transform=Transforms.NTF)
# Datums.OSGB36: Datum(name='OSGB36', ellipsoid=Ellipsoids.Airy1830, transform=Transforms.OSGB36)
# Datums.Potsdam: Datum(name='Potsdam', ellipsoid=Ellipsoids.Bessel1841, transform=Transforms.Bessel1841)
# Datums.Sphere: Datum(name='Sphere', ellipsoid=Ellipsoids.Sphere, transform=Transforms.WGS84)
# Datums.TokyoJapan: Datum(name='TokyoJapan', ellipsoid=Ellipsoids.Bessel1841, transform=Transforms.TokyoJapan)
# Datums.WGS72: Datum(name='WGS72', ellipsoid=Ellipsoids.WGS72, transform=Transforms.WGS72)
# Datums.WGS84: Datum(name='WGS84', ellipsoid=Ellipsoids.WGS84, transform=Transforms.WGS84)

# WGS84: name='WGS84', a=6378137.0, b=6356752.3142499998, f_=298.257223563, f=0.0033528107, e2=0.00669438, e22=0.0067394967, R=6371008.7714166669, Rm=6367435.6797186071, R2=6371007.180920884, R3=6371000.7900107643, Rr=6367449.1458250266,
#        A=6367449.145823415, e=0.0818191908, f=1/298.2572235630, n=0.0016792204(-3.7914875232e-13),
#        Alpha6=(0, 0.0008377318206244698, 7.608527773572307e-07, 1.1976455033294527e-09, 2.4291706072013587e-12, 5.711757677865804e-15, 1.4911177312583895e-17),
#        Beta6=(0, 0.0008377321640579486, 5.905870152220203e-08, 1.6734826652839968e-10, 2.1647980400627059e-13, 3.7879780461686053e-16, 7.2487488906941545e-19)
