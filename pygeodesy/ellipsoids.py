
# -*- coding: utf-8 -*-

u'''Classes L{a_f2Tuple}, L{Ellipsoid}, L{Ellipsoid2}, an L{Ellipsoids} registry
and a dozen functions to convert I{equatorial} radius, I{polar} radius,
I{eccentricities}, I{flattenings} and I{inverse flattening}.

See module L{datums} for more information and other details.

@newfield example: Example, Examples

@var Ellipsoids.Airy1830: Ellipsoid(name='Airy1830', a=6377563.396, b=6356256.90923729, f_=299.3249646, f=0.00334085, f2=0.00335205, n=0.00167322, e=0.08167337, e2=0.00667054, e22=0.00671533, e32=0.00334643, L=10001126.0807165, R1=6370461.23374576, R2=6370459.65470808, R3=6370453.30994572)
@var Ellipsoids.AiryModified: Ellipsoid(name='AiryModified', a=6377340.189, b=6356034.44793853, f_=299.3249646, f=0.00334085, f2=0.00335205, n=0.00167322, e=0.08167337, e2=0.00667054, e22=0.00671533, e32=0.00334643, L=10000776.05340819, R1=6370238.27531284, R2=6370236.69633043, R3=6370230.35179012)
@var Ellipsoids.Australia1966: Ellipsoid(name='Australia1966', a=6378160, b=6356774.71919531, f_=298.25, f=0.00335289, f2=0.00336417, n=0.00167926, e=0.08182018, e2=0.00669454, e22=0.00673966, e32=0.00335851, L=10002001.39064442, R1=6371031.5730651, R2=6371029.9824858, R3=6371023.59124343)
@var Ellipsoids.Bessel1841: Ellipsoid(name='Bessel1841', a=6377397.155, b=6356078.962818, f_=299.1528128, f=0.00334277, f2=0.00335398, n=0.00167418, e=0.08169683, e2=0.00667437, e22=0.00671922, e32=0.00334836, L=10000855.76443237, R1=6370291.09093933, R2=6370289.51012659, R3=6370283.15821522)
@var Ellipsoids.CPM1799: Ellipsoid(name='CPM1799', a=6375738.7, b=6356671.92557493, f_=334.39, f=0.00299052, f2=0.00299949, n=0.0014975, e=0.07727934, e2=0.0059721, e22=0.00600798, e32=0.00299499, L=10000017.52721564, R1=6369383.10852498, R2=6369381.8434158, R3=6369376.76247021)
@var Ellipsoids.Clarke1866: Ellipsoid(name='Clarke1866', a=6378206.4, b=6356583.8, f_=294.97869821, f=0.00339008, f2=0.00340161, n=0.00169792, e=0.08227185, e2=0.00676866, e22=0.00681478, e32=0.00339582, L=10001888.04298286, R1=6370998.86666667, R2=6370997.240633, R3=6370990.70659881)
@var Ellipsoids.Clarke1880: Ellipsoid(name='Clarke1880', a=6378249.145, b=6356514.86954978, f_=293.465, f=0.00340756, f2=0.00341921, n=0.00170669, e=0.0824834, e2=0.00680351, e22=0.00685012, e32=0.00341337, L=10001867.55164747, R1=6371004.38651659, R2=6371002.74366963, R3=6370996.1419165)
@var Ellipsoids.Clarke1880IGN: Ellipsoid(name='Clarke1880IGN', a=6378249.2, b=6356515, f_=293.46602129, f=0.00340755, f2=0.0034192, n=0.00170668, e=0.08248326, e2=0.00680349, e22=0.00685009, e32=0.00341336, L=10001867.69724906, R1=6371004.46666667, R2=6371002.82383112, R3=6370996.22212394)
@var Ellipsoids.Clarke1880Mod: Ellipsoid(name='Clarke1880Mod', a=6378249.145, b=6356514.96582849, f_=293.4663, f=0.00340755, f2=0.0034192, n=0.00170668, e=0.08248322, e2=0.00680348, e22=0.00685009, e32=0.00341335, L=10001867.62720001, R1=6371004.4186095, R2=6371002.77577708, R3=6370996.17408252)
@var Ellipsoids.Delambre1810: Ellipsoid(name='Delambre1810', a=6376428, b=6355957.92616372, f_=311.5, f=0.00321027, f2=0.00322061, n=0.00160772, e=0.08006397, e2=0.00641024, e22=0.0064516, e32=0.00321543, L=9999998.98395793, R1=6369604.64205457, R2=6369603.18419749, R3=6369597.32739068)
@var Ellipsoids.Engelis1985: Ellipsoid(name='Engelis1985', a=6378136.05, b=6356751.32272154, f_=298.2566, f=0.00335282, f2=0.0033641, n=0.00167922, e=0.08181928, e2=0.00669439, e22=0.00673951, e32=0.00335844, L=10001964.20447208, R1=6371007.80757385, R2=6371006.21707085, R3=6370999.82613572)
@var Ellipsoids.Everest1969: Ellipsoid(name='Everest1969', a=6377295.664, b=6356094.667915, f_=300.8017, f=0.00332445, f2=0.00333554, n=0.00166499, e=0.08147298, e2=0.00663785, e22=0.0066822, e32=0.00332998, L=10000788.3115495, R1=6370228.665305, R2=6370227.10178537, R3=6370220.81951617)
@var Ellipsoids.Fisher1968: Ellipsoid(name='Fisher1968', a=6378150, b=6356768.33724438, f_=298.3, f=0.00335233, f2=0.00336361, n=0.00167898, e=0.08181333, e2=0.00669342, e22=0.00673853, e32=0.00335795, L=10001988.52191361, R1=6371022.77908146, R2=6371021.18903735, R3=6371014.79995034)
@var Ellipsoids.GEM10C: Ellipsoid(name='GEM10C', a=6378137, b=6356752.31424783, f_=298.2572236, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181919, e2=0.00669438, e22=0.0067395, e32=0.00335843, L=10001965.7293148, R1=6371008.77141594, R2=6371007.18091936, R3=6371000.79001004)
@var Ellipsoids.GRS67: Ellipsoid(name='GRS67', a=6378160, b=6356774.51609071, f_=298.24716743, f=0.00335292, f2=0.0033642, n=0.00167928, e=0.08182057, e2=0.00669461, e22=0.00673973, e32=0.00335854, L=10002001.2312605, R1=6371031.50536357, R2=6371029.91475409, R3=6371023.52339014)
@var Ellipsoids.GRS80: Ellipsoid(name='GRS80', a=6378137, b=6356752.31414035, f_=298.2572221, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181919, e2=0.00669438, e22=0.0067395, e32=0.00335843, L=10001965.72923046, R1=6371008.77138012, R2=6371007.18088351, R3=6371000.78997413)
@var Ellipsoids.Helmert1906: Ellipsoid(name='Helmert1906', a=6378200, b=6356818.16962789, f_=298.3, f=0.00335233, f2=0.00336361, n=0.00167898, e=0.08181333, e2=0.00669342, e22=0.00673853, e32=0.00335795, L=10002066.93013953, R1=6371072.7232093, R2=6371071.13315272, R3=6371064.74401563)
@var Ellipsoids.IERS1989: Ellipsoid(name='IERS1989', a=6378136, b=6356751.30156878, f_=298.257, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181922, e2=0.00669438, e22=0.0067395, e32=0.00335843, L=10001964.14856985, R1=6371007.76718959, R2=6371006.17669088, R3=6370999.78577296)
@var Ellipsoids.IERS1992TOPEX: Ellipsoid(name='IERS1992TOPEX', a=6378136.3, b=6356751.61659215, f_=298.25722356, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181919, e2=0.00669438, e22=0.0067395, e32=0.00335843, L=10001964.63159783, R1=6371008.07219738, R2=6371006.48170097, R3=6371000.09079235)
@var Ellipsoids.IERS2003: Ellipsoid(name='IERS2003', a=6378136.6, b=6356751.85797165, f_=298.25642, f=0.00335282, f2=0.0033641, n=0.00167922, e=0.0818193, e2=0.0066944, e22=0.00673951, e32=0.00335844, L=10001965.05683465, R1=6371008.35265722, R2=6371006.76215217, R3=6371000.37120876)
@var Ellipsoids.Intl1924: Ellipsoid(name='Intl1924', a=6378388, b=6356911.94612795, f_=297, f=0.003367, f2=0.00337838, n=0.00168634, e=0.08199189, e2=0.00672267, e22=0.00676817, e32=0.00337267, L=10002288.29898944, R1=6371229.31537598, R2=6371227.71133444, R3=6371221.26587487)
@var Ellipsoids.Intl1967: Ellipsoid(name='Intl1967', a=6378157.5, b=6356772.2, f_=298.24961539, f=0.0033529, f2=0.00336418, n=0.00167926, e=0.08182023, e2=0.00669455, e22=0.00673967, e32=0.00335852, L=10001997.44859308, R1=6371029.06666667, R2=6371027.47608389, R3=6371021.08482752)
@var Ellipsoids.Krassovski1940: Ellipsoid(name='Krassovski1940', a=6378245, b=6356863.01877305, f_=298.3, f=0.00335233, f2=0.00336361, n=0.00167898, e=0.08181333, e2=0.00669342, e22=0.00673853, e32=0.00335795, L=10002137.49754285, R1=6371117.67292435, R2=6371116.08285656, R3=6371109.69367439)
@var Ellipsoids.Krassowsky1940: Ellipsoid(name='Krassowsky1940', a=6378245, b=6356863.01877305, f_=298.3, f=0.00335233, f2=0.00336361, n=0.00167898, e=0.08181333, e2=0.00669342, e22=0.00673853, e32=0.00335795, L=10002137.49754285, R1=6371117.67292435, R2=6371116.08285656, R3=6371109.69367439)
@var Ellipsoids.Maupertuis1738: Ellipsoid(name='Maupertuis1738', a=6397300, b=6363806.28272251, f_=191, f=0.0052356, f2=0.00526316, n=0.00262467, e=0.10219488, e2=0.01044379, e22=0.01055402, e32=0.00524931, L=10022566.69846922, R1=6386135.42757417, R2=6386131.54144847, R3=6386115.88628229)
@var Ellipsoids.Mercury1960: Ellipsoid(name='Mercury1960', a=6378166, b=6356784.28360711, f_=298.3, f=0.00335233, f2=0.00336361, n=0.00167898, e=0.08181333, e2=0.00669342, e22=0.00673853, e32=0.00335795, L=10002013.61254591, R1=6371038.76120237, R2=6371037.17115427, R3=6371030.78205124)
@var Ellipsoids.Mercury1968Mod: Ellipsoid(name='Mercury1968Mod', a=6378150, b=6356768.33724438, f_=298.3, f=0.00335233, f2=0.00336361, n=0.00167898, e=0.08181333, e2=0.00669342, e22=0.00673853, e32=0.00335795, L=10001988.52191361, R1=6371022.77908146, R2=6371021.18903735, R3=6371014.79995034)
@var Ellipsoids.NWL1965: Ellipsoid(name='NWL1965', a=6378145, b=6356759.76948868, f_=298.25, f=0.00335289, f2=0.00336417, n=0.00167926, e=0.08182018, e2=0.00669454, e22=0.00673966, e32=0.00335851, L=10001977.86818326, R1=6371016.58982956, R2=6371014.999254, R3=6371008.60802666)
@var Ellipsoids.OSU86F: Ellipsoid(name='OSU86F', a=6378136.2, b=6356751.51693008, f_=298.2572236, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181919, e2=0.00669438, e22=0.0067395, e32=0.00335843, L=10001964.47478349, R1=6371007.97231003, R2=6371006.38181364, R3=6370999.99090512)
@var Ellipsoids.OSU91A: Ellipsoid(name='OSU91A', a=6378136.3, b=6356751.6165948, f_=298.2572236, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181919, e2=0.00669438, e22=0.0067395, e32=0.00335843, L=10001964.63159991, R1=6371008.07219827, R2=6371006.48170186, R3=6371000.09079324)
@var Ellipsoids.Plessis1817: Ellipsoid(name='Plessis1817', a=6376523, b=6355862.93325557, f_=308.64, f=0.00324002, f2=0.00325055, n=0.00162264, e=0.08043347, e2=0.00646954, e22=0.00651167, e32=0.00324527, L=9999999.1100364, R1=6369636.31108519, R2=6369634.82608583, R3=6369628.85999667)
@var Ellipsoids.SGS85: Ellipsoid(name='SGS85', a=6378136, b=6356751.30156878, f_=298.257, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181922, e2=0.00669438, e22=0.0067395, e32=0.00335843, L=10001964.14856985, R1=6371007.76718959, R2=6371006.17669087, R3=6370999.78577296)
@var Ellipsoids.SoAmerican1969: Ellipsoid(name='SoAmerican1969', a=6378160, b=6356774.71919531, f_=298.25, f=0.00335289, f2=0.00336417, n=0.00167926, e=0.08182018, e2=0.00669454, e22=0.00673966, e32=0.00335851, L=10002001.39064442, R1=6371031.5730651, R2=6371029.98248581, R3=6371023.59124343)
@var Ellipsoids.Sphere: Ellipsoid(name='Sphere', a=6371008.771415, b=6371008.771415, f_=0, f=0, f2=0, n=0, e=0, e2=0, e22=0, e32=0, L=10007557.17611675, R1=6371008.771415, R2=6371008.771415, R3=6371008.771415)
@var Ellipsoids.SphereAuthalic: Ellipsoid(name='SphereAuthalic', a=6371000, b=6371000, f_=0, f=0, f2=0, n=0, e=0, e2=0, e22=0, e32=0, L=10007543.39801029, R1=6371000, R2=6371000, R3=6371000)
@var Ellipsoids.SpherePopular: Ellipsoid(name='SpherePopular', a=6378137, b=6378137, f_=0, f=0, f2=0, n=0, e=0, e2=0, e22=0, e32=0, L=10018754.17139462, R1=6378137, R2=6378137, R3=6378137)
@var Ellipsoids.Struve1860: Ellipsoid(name='Struve1860', a=6378298.3, b=6356657.14266956, f_=294.73, f=0.00339294, f2=0.00340449, n=0.00169935, e=0.0823065, e2=0.00677436, e22=0.00682056, e32=0.00339869, L=10002017.83655713, R1=6371084.58088985, R2=6371082.95208988, R3=6371076.40691418)
@var Ellipsoids.WGS60: Ellipsoid(name='WGS60', a=6378165, b=6356783.28695944, f_=298.3, f=0.00335233, f2=0.00336361, n=0.00167898, e=0.08181333, e2=0.00669342, e22=0.00673853, e32=0.00335795, L=10002012.04438139, R1=6371037.76231981, R2=6371036.17227197, R3=6371029.78316993)
@var Ellipsoids.WGS66: Ellipsoid(name='WGS66', a=6378145, b=6356759.76948868, f_=298.25, f=0.00335289, f2=0.00336417, n=0.00167926, e=0.08182018, e2=0.00669454, e22=0.00673966, e32=0.00335851, L=10001977.86818326, R1=6371016.58982956, R2=6371014.999254, R3=6371008.60802666)
@var Ellipsoids.WGS72: Ellipsoid(name='WGS72', a=6378135, b=6356750.52001609, f_=298.26, f=0.00335278, f2=0.00336406, n=0.0016792, e=0.08181881, e2=0.00669432, e22=0.00673943, e32=0.0033584, L=10001962.74919858, R1=6371006.84000536, R2=6371005.24953886, R3=6370998.85875069)
@var Ellipsoids.WGS84: Ellipsoid(name='WGS84', a=6378137, b=6356752.31424518, f_=298.25722356, f=0.00335281, f2=0.00336409, n=0.00167922, e=0.08181919, e2=0.00669438, e22=0.0067395, e32=0.00335843, L=10001965.72931272, R1=6371008.77141506, R2=6371007.18091847, R3=6371000.79000915)
'''

# make sure int/int division yields float quotient
from __future__ import division
division = 1 / 2  # double check int division, see .datums, .fmath, .utily
if not division:
    raise ImportError('%s 1/2 == %d' % ('division', division))
del division

from pygeodesy.basics import EPS, EPS1, INF, PI2, PI_2, R_M, \
                             property_doc_, property_RO, _xinstanceof
from pygeodesy.errors import _AssertionError, _ValueError
from pygeodesy.fmath import _2_3rd, cbrt, cbrt2, fdot, fpowers, Fsum, \
                             fsum_, hypot1, hypot2, sqrt3  # PYCHOK _2_3rd
from pygeodesy.interns import _flt, _lat_, _lat0_, _lat1_, _lon0_, _lon1_, \
                              _n_a_, _negative_, NN, _x_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY
from pygeodesy.named import Distance2Tuple, _NamedEnum, _NamedEnumItem, _NamedTuple
from pygeodesy.streprs import _e, instr, _Fmt, fstr
from pygeodesy.units import Bearing_, Distance, Float, Float_, Lam_, Lat, \
                            Phi, Phi_, Radius, Radius_, Scalar
from pygeodesy.utily import atan2d, degrees90, degrees360, degrees2m, \
                            m2degrees, m2km, m2NM, m2SM

from math import asinh, atan, atan2, atanh, copysign, cos, \
                 degrees, exp, hypot, sin, sinh, sqrt, tan

R_M  = Radius(R_M,       name='R_M')   #: Mean (spherical) earth radius (C{meter}).
R_MA = Radius(6378137.0, name='R_MA')  #: Equatorial (major) earth radius (C{meter}) WGS84, EPSG:3785.
R_MB = Radius(6356752.0, name='R_MB')  #: Polar (minor) earth radius (C{meter}) WGS84, EPSG:3785.
R_KM = Radius(m2km(R_M), name='R_KM')  #: Mean (spherical) earth radius (C{KM}, kilo meter).
R_NM = Radius(m2NM(R_M), name='R_NM')  #: Mean (spherical) earth radius (C{NM}, nautical miles).
R_SM = Radius(m2SM(R_M), name='R_SM')  #: Mean (spherical) earth radius (C{SM}, statute miles).
# See <https://www.EdWilliams.org/avform.htm>,
# <https://www.DTIC.mil/dtic/tr/fulltext/u2/a216843.pdf> and
# <https://GitHub.com/NASA/MultiDop/blob/master/src/share/man/man3/geog_lib.3>
# based on International Standard Nautical Mile of 1,852 meter (1' latitude)
R_FM = Radius(6371000.0,       name='R_FM')  #: Former FAI Sphere earth radius (C{meter}).
R_VM = Radius(6366707.0194937, name='R_VM')  #: Aviation/Navigation earth radius (C{meter}).
# R_ = Radius(6372797.560856,  name='R_')    #: XXX some other earth radius???

__all__ = _ALL_LAZY.ellipsoids
__version__ = '20.09.02'

_PI4   = PI2 * 2          # PYCHOK used!
_PI4_3 = PI2 * _2_3rd     # PYCHOK used!
_TOL   = sqrt(EPS * 0.1)  # for Ellipsoid.estauf, imported by .ups, testEllipdoidal.py


def _atand(y_x):
    '''Return C{atan} angle in C{degrees}.
    '''
    return atan2d(y_x, 1)


def _aux(lat, inverse, auxLat, clip=90):
    '''Return named auxiliary latitude degrees.
    '''
    return Lat(lat, clip=clip, name=_lat_ if inverse else auxLat)


def _c2_s2(phi):
    '''Return 2-tuple C{(cos(phi)**2, sin(phi)**2)}.
    '''
    c2 = cos(phi)**2
    return c2, (1 - c2)


def _4Ecef(this, Ecef):
    '''Return an ECEF converter.
    '''
    from pygeodesy.ecef import EcefKarney, EcefVeness, EcefYou

    if Ecef is None:
        Ecef = EcefKarney
    else:
        _xinstanceof(EcefKarney, EcefVeness, EcefYou, Ecef=Ecef)
    return Ecef(this, name=this.name)  # datum or ellipsoid


class a_f2Tuple(_NamedTuple):
    '''2-Tuple C{(a, f)} specifying an ellipsoid by I{equatorial}
       radius (C{meter}) and scalar I{flattening}.

       @note: M{abs(f) < EPS} is forced to M{f = 0}, spherical.

       @see: Class L{Ellipsoid2}.
    '''
    _Names_ = ('a', 'f')

    def __new__(cls, a, f):
        a = Radius_(_flt(a), name=Ellipsoid.a.name)
        f = Float_(_flt(f), name=Ellipsoid.f.name, low=None, high=EPS1)
        if abs(f) < EPS:  # force spherical
            f = Float(_flt(0), name=Ellipsoid.f.name)
        return _NamedTuple.__new__(cls, a, f)

    @property_RO
    def b(self):
        '''Get the I{polar} radius (C{meter}), M{a * (1 - f)}.
        '''
        return a_f2b(self.a, self.f)  # PYCHOK .a and .f

    @property_RO
    def f_(self):
        '''Get the I{inverse} flattening (C{float}), M{1 / f} == M{a / (a - b)}.
        '''
        return f2f_(self.f)  # PYCHOK .f


class Curvature2Tuple(_NamedTuple):
    '''2-Tuple C{(meridional, prime_vertical)} of radii of curvature,
       both in C{meter}.
    '''
    _Names_ = ('meridional', 'prime_vertical')


class Ellipsoid(_NamedEnumItem):
    '''Ellipsoid with I{equatorial} and I{polar} radius, I{flattening},
       I{inverse flattening} and other, often used, cached attributes,
       supporting I{spherical} and I{oblate} and I{prolate} ellipsoidal models.
    '''
    _a  = 0  #: Equatorial (major) radius, semi-axis (C{meter}).
    _b  = 0  #: Polar (minor) radius, semi-axis (C{meter}): a * (f - 1) / f
    _f  = 0  #: (1st) flattening: (a - b) / a
    _f_ = 0  #: (INTERNAL) Inverse flattening: 1 / f = a / (a - b)

    _f2 = None  #: 2nd flattening: (a - b) / b  # rarely used
    _n  = None  #: 3rd flattening: f / (2 - f) = (a - b) / (a + b)  # for .A and .utm

    _a2  = None  #: a**2
    _a2_ = None  #: (INTERNAL) (1 / a**2)  for .ellipsiodalNvector.Cartesian.toNvector
    _a_b = None  #: (a / b) = 1 / (1 - f)  for .ellipsoidalNvector.Nvector.toCartesian
    _b2  = None  #: b**2
    _b_a = None  #: (b / a) = 1 - f  for .formy, .R2, .Rgeocentric below

    # curvatures <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>
    _a2_b = None  #: (INTERNAL) Meridional radius of curvature at poles: a**2 / b (C{meter})
    _b2_a = None  #: (INTERNAL) Meridional radius of curvature at equator: b**2 / a (C{meter})

    # eccentricities
    _e    = None  #: (1st) eccentricity: sqrt(1 - (b / a)**2))  # for utm
    _e2   = None  #: (1st) eccentricity squared: f * (2 - f) = 1 - (b / a)**2
    _e22  = None  #: 2nd eccentricity squared: e2 / (1 - e2) = (a / b)**2 - 1
    _e32  = None  #: 3rd eccentricity squared: e2 / (2 - e2) = (a**2 - b**2 / (a**2 + b**2)
    _e12  = None  #: 1 - e2 = (1 - f)**2 for .ellipsoidalNvector.Cartesian.toNvector, .ecef .utm
    _e4   = None  #: e2=**4 = e2**2 for .ellipsoidalNvector.Cartesian.toNvector, .ecef
    _es_c = None  #: M{(1 - f) * exp(es_atanh(1))}

    # fixed earth radii from <https://WikiPedia.org/wiki/Earth_radius>
    _A  = None  #: (INTERNAL) UTM meridional radius
    _L  = None  #: (INTERNAL) Quarter meridian: b * Elliptic(-e2 / (1 - e2)).E (C{meter})
    _R1 = None  #: (INTERNAL) Mean earth radius: (2 * a + b) / 3 per IUGG definition (C{meter})
    _R2 = None  #: (INTERNAL) Authalic radius: sqrt((a**2 + b**2 * atanh(e) / e) / 2) (C{meter})
#   _c  = None  #: (INTERNAL) Authalic radius: equ (60) in Karney's "Algorithms for Geodesics"
    _R3 = None  #: (INTERNAL) Volumetric radius: (a * a * b)**(1/3) (C{meter})
    _Rb = None  #: (INTERNAL) Biaxial mean earth radius: sqrt((a**2 * b**2) / 2) (C{meter})
    _Rg = None  #: (INTERNAL) Geometric mean earth radius: sqrt(a * b) (C{meter})
    _Rr = None  #: (INTERNAL) Rectifying radius: ((a**(3/2) + b**(3/2)) / 2)**(2/3) (C{meter})

    _ab_90  = None  #: (a - b) / 90  # for .Rlat below
    _area   = None  #: (INTERNAL) Surface area: 4 * PI * R2**2
    _b2_a2  = None  #: b**2 / a**2 = (b / a)**2 = _b_a**2 = (1 - f)**2
    _volume = None  #: (INTERNAL) Volume: 4 / 3 * PI * a**2 * b

    _AlphaKs = None  #: (INTERNAL) Up to 8th-order Krüger Alpha series
    _BetaKs  = None  #: (INTERNAL) Up to 8th-order Krüger Beta series
    _KsOrder = 8     #: (INTERNAL) Krüger series order (4, 6 or 8)
    _Mabcd   = None  #: (INTERNAL) OSGR meridional coefficients

    _albersC  = None  #: (INTERNAL) Cached L{AlbersEqualAreaCylindrical} instance
    _elliptic = None  #: (INTERNAL) Cached elliptic function L{Elliptic} instance
    _geodesic = None  #: (INTERNAL) Cached C{karney._wrapped_.Geodesic} instance
    _Math     = None  #: (INTERNAL) Cached C{geographiclib.geomath.Math} module

    def __init__(self, a, b=None, f_=None, name=NN):
        '''New L{Ellipsoid} from I{equatorial} and I{polar} radius or
           I{equatorial} radius and I{inverse flattening}.

           @arg a: Equatorial (major) radius, semi-axis (C{meter}).
           @arg b: Optional, polar (minor) radius, semi-axis (C{meter}).
           @arg f_: Inverse flattening: M{a / (a - b)} (C{float} >>> 1.0).
           @kwarg name: Optional, unique name (C{str}).

           @raise NameError: Ellipsoid with that B{C{name}} already exists.

           @raise ValueError: Invalid B{C{a}}, B{C{b}} or B{C{f_}}.

           @note: M{abs(1 / f_) < EPS} is considered M{1 / f_ = 0}, spherical.
        '''
        try:
            a = Radius_(_flt(a), name=Ellipsoid.a.name)  # low=EPS
            if b:
                b = Radius_(_flt(b), name=Ellipsoid.b.name)  # low=EPS
                f = a_b2f(a, b)
                if f_ is None:
                    f_ = f2f_(f)
            elif f_:
                b = a_f_2b(a, f_)  # a * (f_ - 1) / f_
                b = Radius_(_flt(b), name=Ellipsoid.b.name)  # low=EPS
                f = a_b2f(a, b)
            else:  # only a, spherical
                f = f_ = 0
                b = a  # superfluous

            if abs(f) < EPS or a == b or not f_:  # spherical
                b = a
                f = f_ = 0
            elif f > EPS1:  # sanity check
                raise ValueError

        except (TypeError, ValueError) as x:
            raise _ValueError(instr(self, a=a, b=b, f_=f_), txt=str(x))

        self._a  = a
        self._b  = b
        self._f  = Float(_flt(f),  name=Ellipsoid.f.name)
        self._f_ = Float(_flt(f_), name=Ellipsoid.f_.name)

        self._register(Ellipsoids, name)

        if f and f_:  # see .test/testEllipsoidal.py
            self._assert(1 / f,  f_=f_, eps=_TOL)
            self._assert(1 / f_, f=f,   eps=_TOL)
        self._assert(self.b2_a2, e12=self.e12, eps=EPS)

    def __eq__(self, other):
        '''Compare this and an other ellipsoid.

           @arg other: The other ellipsoid (L{Ellipsoid} or L{Ellipsoid2}).

           @return: C{True} if equal, C{False} otherwise.
        '''
        return self is other or (isinstance(other, Ellipsoid) and
                                 self.a == other.a and
                                (self.b == other.b or
                                 self.f == other.f))

    def _Kseries(self, *AB8Ks):
        '''(INTERNAL) Compute the 4-, 6- or 8-th order I{Krüger} Alpha
           or Beta series coefficients per I{Karney} 2011, 'Transverse
           Mercator with an accuracy of a few nanometers', U{page 7,
           equations 35 and 36<https://Arxiv.org/pdf/1002.1417v3.pdf>}.

           @arg AB8Ks: 8-Tuple of 8-th order I{Krüger} Alpha or Beta series
                       coefficient tuples.

           @return: I{Krüger} series coefficients (C{KsOrder}-tuple).

           @see: I{Karney}s 30-th order U{TMseries30
                 <https://GeographicLib.SourceForge.io/html/tmseries30.html>}.
        '''
        k = self.KsOrder
        ns = fpowers(self.n, k)
        return tuple(fdot(AB8Ks[i][:k-i], *ns[i:]) for i in range(k))

    @property_RO
    def a(self):
        '''Get the equatorial (major) radius, semi-axis (C{meter}).
        '''
        return self._a

    equatoradius = a  # = Requatorial

    @property_RO
    def a2(self):
        '''Get the equatorial (major) radius I{squared} (C{float}), M{a**2}.
        '''
        if self._a2 is None:
            self._a2 = Float(self.a**2, name=Ellipsoid.a2.name)
        return self._a2

    @property_RO
    def a2_(self):
        '''Get the inverse of the equatorial (major) radius I{squared} (C{float}), M{1 / a**2}.
        '''
        if self._a2_ is None:
            self._a2_ = Float(1 / self.a2, name=Ellipsoid.a2_.name)
        return self._a2_  # (1 / a**2)

    @property_RO
    def a_b(self):
        '''Get ratio M{a / b} (C{float}).
        '''
        if self._a_b is None:
            self._a_b = Float(self.a / self.b if self.f else 1, name=Ellipsoid.a_b.name)
        return self._a_b

    @property_RO
    def a2_b(self):
        '''Get the polar meridional radius of curvature (C{meter}), M{a**2 / b}, see C{rocPolar}.

           @see: U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}
                 and U{Moritz, H. (1980), Geodetic Reference System 1980
                 <https://WikiPedia.org/wiki/Earth_radius#cite_note-Moritz-2>}.

           @note: Symbol C{c} is used by IUGG and IERS for the U{polar radius of
                  curvature<https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        if self._a2_b is None:
            self._a2_b = Radius(self.a2 / self.b if self.f else self.a, name=Ellipsoid.a2_b.name)
        return self._a2_b

    @property_RO
    def a_f(self):
        '''Get equatorial radius and flattening (L{a_f2Tuple}).
        '''
        return self._xnamed(a_f2Tuple(self.a, self.f))

    @property_RO
    def A(self):
        '''Get the UTM meridional radius (C{meter}).
        '''
        if self._A is None:
            # <https://GeographicLib.SourceForge.io/html/transversemercator.html>
            A, n = self.a, self.n
            if n and self.f:
                A = A / (1 + n) * (fsum_(65536, 16384 * n**2,
                                                 1024 * n**4,
                                                  256 * n**6,
                                                  100 * n**8,
                                                   49 * n**10) / 65536)
            # <https://www.MyGeodesy.id.AU/documents/Karney-Krueger%20equations.pdf>
            # A = self.a / (1 + n) * (fhorner(n**2, 16384, 4096, 256, 64, 25) / 16384)
            self._A = Radius(A, name=Ellipsoid.A.name)
        return self._A

    @property_RO
    def area(self):
        '''Get the ellipsoid's surface area (C{meter**2}), M{4 * PI * R2**2} or M{4 * PI * a**2}.
        '''
        if self._area is None:
            r = self.R2 if self.f else self.a  # authalic radius
            self._area = Float(_PI4 * r**2, name=Ellipsoid.area.name)
        return self._area

    @property_RO
    def AlphaKs(self):
        '''Get the I{Krüger} U{Alpha series coefficients<https://GeographicLib.SourceForge.io/html/tmseries30.html>} (C{KsOrder}-tuple).
        '''
        if self._AlphaKs is None:
            self._AlphaKs = self._Kseries(  # XXX int/int quotients may require  from __future__ import division
                # n    n**2   n**3      n**4         n**5            n**6                 n**7                     n**8
                (1/2, -2/3,   5/16,    41/180,    -127/288,       7891/37800,         72161/387072,        -18975107/50803200),
                     (13/48, -3/5,    557/1440,    281/630,   -1983433/1935360,       13769/28800,         148003883/174182400),      # PYCHOK unaligned
                            (61/240, -103/140,   15061/26880,   167603/181440,    -67102379/29030400,       79682431/79833600),       # PYCHOK unaligned
                                   (49561/161280, -179/168,    6601661/7257600,       97445/49896,      -40176129013/7664025600),     # PYCHOK unaligned
                                                (34729/80640, -3418889/1995840,    14644087/9123840,      2605413599/622702080),      # PYCHOK unaligned
                                                            (212378941/319334400, -30705481/10378368,   175214326799/58118860800),    # PYCHOK unaligned
                                                                                (1522256789/1383782400, -16759934899/3113510400),     # PYCHOK unaligned
                                                                                                      (1424729850961/743921418240,))  # PYCHOK unaligned
        return self._AlphaKs

    def _assert(self, val, eps=_TOL, f0=0, **name_value):
        '''(INTERNAL) Assert a C{name=value} vs C{val}.
        '''
        for n, v in name_value.items():
            if abs(v - val) > eps:
                raise _AssertionError('%s=%.12g vs %.12g' % (self._dot_(n), v, val))
            return Float(v if self.f else f0, name=n)
        raise _AssertionError(self._dot_(self._assert.__name__), name_value)

    def auxAuthalic(self, lat, inverse=False):
        '''Compute the I{authalic} auxiliary latitude or inverse thereof.

           @arg lat: The geodetic (or I{authalic}) latitude (C{degrees90}).
           @kwarg inverse: If C{True}, B{C{lat}} is the I{authalic} and
                           return the geodetic latitude (C{bool}).

           @return: The I{authalic} (or geodetic) latitude in C{degrees90}.

           @see: U{Inverse-/AuthalicLatitude<https://geographiclib.sourceforge.io/
                 html/classGeographicLib_1_1Ellipsoid.html>}, U{Authalic latitude
                 <https://WikiPedia.org/wiki/Latitude#Authalic_latitude>}, and
                 U{Snyder<https://pubs.er.USGS.gov/djvu/PP/PP_1395.pdf>}, p 16.

        '''
        if self.isEllipsoidal:
            if self._albersC is None:
                from pygeodesy.albers import AlbersEqualAreaCylindrical as _AC
                self._albersC = _AC(datum=self, name=self.name)

            f = self._albersC._tanf if inverse else self._albersC._txif  # PYCHOK attr
            lat = _atand(f(tan(Phi_(lat))))  # PYCHOK attr
        return _aux(lat, inverse, Ellipsoid.auxAuthalic.__name__)

    def auxConformal(self, lat, inverse=False):
        '''Compute the I{conformal} auxiliary latitude or inverse thereof.

           @arg lat: The geodetic (or I{conformal}) latitude (C{degrees90}).
           @kwarg inverse: If C{True}, B{C{lat}} is the I{conformal} and
                           return the geodetic latitude (C{bool}).

           @return: The I{conformal} (or geodetic) latitude in C{degrees90}.

           @see: U{Inverse-/ConformalLatitude<https://geographiclib.sourceforge.io/
                 html/classGeographicLib_1_1Ellipsoid.html>}, U{Conformal latitude
                 <https://WikiPedia.org/wiki/Latitude#Conformal_latitude>}, and
                 U{Snyder<https://pubs.er.USGS.gov/djvu/PP/PP_1395.pdf>}, pp 15-16.
        '''
        if self.isEllipsoidal:
            f = self.es_tauf if inverse else self.es_taupf  # PYCHOK attr
            lat = _atand(f(tan(Phi_(lat))))  # PYCHOK attr
        return _aux(lat, inverse, Ellipsoid.auxConformal.__name__)

    def auxGeocentric(self, lat, inverse=False):
        '''Compute the I{geocentric} auxiliary latitude or inverse thereof.

           @arg lat: The geodetic (or I{geocentric}) latitude (C{degrees90}).
           @kwarg inverse: If C{True}, B{C{lat}} is the geocentric and
                           return the I{geocentric} latitude (C{bool}).

           @return: The I{geocentric} (or geodetic) latitude in C{degrees90}.

           @see: U{Inverse-/GeocentricLatitude<https://geographiclib.sourceforge.io/
                 html/classGeographicLib_1_1Ellipsoid.html>}, U{Geocentric latitude
                 <https://WikiPedia.org/wiki/Latitude#Geocentric_latitude>}, and
                 U{Snyder<<https://pubs.er.USGS.gov/djvu/PP/PP_1395.pdf>}, pp 17-18.
        '''
        if self.isEllipsoidal:
            f = (1 / self.b2_a2) if inverse else self.b2_a2
            lat = _atand(f * tan(Phi_(lat)))
        return _aux(lat, inverse, Ellipsoid.auxGeocentric.__name__)

    def auxIsometric(self, lat, inverse=False):
        '''Compute the I{isometric} auxiliary latitude or inverse thereof.

           @arg lat: The geodetic (or I{isometric}) latitude (C{degrees}).
           @kwarg inverse: If C{True}, B{C{lat}} is the I{isometric} and
                           return the geodetic latitude (C{bool}).

           @return: The I{isometric} (or geodetic) latitude in C{degrees}.

           @note: The I{isometric} latitude for geodetic C{+/-90} is far
                  outside the C{[-90..+90]} range but the inverse
                  thereof is the original geodetic latitude.

           @see: U{Inverse-/IsometricLatitude<https://geographiclib.sourceforge.io/
                 html/classGeographicLib_1_1Ellipsoid.html>}, U{Isometric latitude
                 <https://WikiPedia.org/wiki/Latitude#Isometric_latitude>}, and
                 U{Snyder<https://pubs.er.USGS.gov/djvu/PP/PP_1395.pdf>}, pp 15-16.
        '''
        if self.isEllipsoidal:
            r = Phi_(lat, clip=0)
            lat = degrees( atan(self.es_tauf(sinh(r))) if inverse else
                          asinh(self.es_taupf(tan(r))))
        # clip=0, since auxIsometric(+/-90) is far outside [-90..+90]
        return _aux(lat, inverse, Ellipsoid.auxIsometric.__name__, clip=0)

    def auxParametric(self, lat, inverse=False):
        '''Compute the I{parametric} auxiliary latitude or inverse thereof.

           @arg lat: The geodetic (or I{parametric}) latitude (C{degrees90}).
           @kwarg inverse: If C{True}, B{C{lat}} is the I{parametric} and
                           return the geodetic latitude (C{bool}).

           @return: The I{parametric} (or geodetic) latitude in C{degrees90}.

           @see: U{Inverse-/ParametricLatitude<https://geographiclib.sourceforge.io/
                 html/classGeographicLib_1_1Ellipsoid.html>}, U{Parametric latitude
                 <https://WikiPedia.org/wiki/Latitude#Parametric_(or_reduced)_latitude>},
                 and U{Snyder<https://pubs.er.USGS.gov/djvu/PP/PP_1395.pdf>}, p 18.
        '''
        if self.isEllipsoidal:
            f = (1 / self.b_a) if inverse else self.b_a
            lat = _atand(f * tan(Phi_(lat)))
        return _aux(lat, inverse, Ellipsoid.auxParametric.__name__)

    auxReduced = auxParametric  # synonyms

    def auxRectifying(self, lat, inverse=False):
        '''Compute the I{rectifying} auxiliary latitude or inverse thereof.

           @arg lat: The geodetic (or I{rectifying}) latitude (C{degrees90}).
           @kwarg inverse: If C{True}, B{C{lat}} is the I{rectifying} and
                           return the geodetic latitude (C{bool}).

           @return: The I{rectifying} (or geodetic) latitude in C{degrees90}.

           @see: U{Inverse-/RectifyingLatitude<https://geographiclib.sourceforge.io/
                 html/classGeographicLib_1_1Ellipsoid.html>}, U{Rectifying latitude
                 <https://WikiPedia.org/wiki/Latitude#Rectifying_latitude>}, and
                 U{Snyder<https://pubs.er.USGS.gov/djvu/PP/PP_1395.pdf>}, pp 16-17.
        '''
        lat = Lat(lat)
        if 0 < abs(lat) < 90 and self.isEllipsoidal:
            if inverse:
                e = self._elliptic_e22
                lat = degrees90(e.fEinv(e.cE * lat / 90))
                lat = self.auxParametric(lat, inverse=True)
            else:
                lat = 90 * self.Llat(lat) / self.L
        return _aux(lat, inverse, Ellipsoid.auxRectifying.__name__)

    @property_RO
    def b(self):
        '''Get the polar (minor) radius, semi-axis (C{meter}).
        '''
        return self._b

    polaradius = b  # = Rpolar

    @property_RO
    def b_a(self):
        '''Get ratio M{b / a} (C{float}), M{1 - f}.
        '''
        if self._b_a is None:
            self._b_a = self._assert(self.b / self.a, b_a=1 - self.f, f0=1)
        return self._b_a

    @property_RO
    def b2(self):
        '''Get the polar (minor) radius I{squared} (C{float}), M{b**2}.
        '''
        if self._b2 is None:
            self._b2 = Float(self.b**2, name=Ellipsoid.b2.name)
        return self._b2

    @property_RO
    def b2_a(self):
        '''Get the equatorial meridional radius of curvature (C{meter}), M{b**2 / a}.

           @see: U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        if self._b2_a is None:
            self._b2_a = Radius(self.b2 / self.a if self.f else self.b, name=Ellipsoid.b2_a.name)
        return self._b2_a

    @property_RO
    def b2_a2(self):
        '''Get ratio M{(b / a)**2} (C{float}), M{(1 - f)**2}, see C{e12}.
        '''
        if self._b2_a2 is None:
            self._b2_a2 = Float(self.b_a**2 if self.f else 1, name=Ellipsoid.b2_a2.name)
        return self._b2_a2

    @property_RO
    def BetaKs(self):
        '''Get the I{Krüger} U{Beta series coefficients<https://GeographicLib.SourceForge.io/html/tmseries30.html>} (C{KsOrder}-tuple).
        '''
        if self._BetaKs is None:
            self._BetaKs = self._Kseries(  # XXX int/int quotients may require  from __future__ import division
                # n    n**2  n**3     n**4        n**5            n**6                 n**7                   n**8
                (1/2, -2/3, 37/96,   -1/360,    -81/512,      96199/604800,     -5406467/38707200,      7944359/67737600),
                      (1/48, 1/15, -437/1440,    46/105,   -1118711/3870720,       51841/1209600,      24749483/348364800),       # PYCHOK unaligned
                           (17/480, -37/840,   -209/4480,      5569/90720,       9261899/58060800,     -6457463/17740800),        # PYCHOK unaligned
                                  (4397/161280, -11/504,    -830251/7257600,      466511/2494800,     324154477/7664025600),      # PYCHOK unaligned
                                              (4583/161280, -108847/3991680,    -8005831/63866880,     22894433/124540416),       # PYCHOK unaligned
                                                          (20648693/638668800, -16363163/518918400, -2204645983/12915302400),     # PYCHOK unaligned
                                                                              (219941297/5535129600, -497323811/12454041600),     # PYCHOK unaligned
                                                                                                  (191773887257/3719607091200,))  # PYCHOK unaligned
        return self._BetaKs

    @property_RO
    def c(self):
        '''Get the I{authalic} earth radius (C{meter}), see C{R2}.

           @note: Symbol C{c} in U{equation 60
                  <https://Link.Springer.com/article/10.1007%2Fs00190-012-0578-z>}.

           @note: Symbol C{c} is used by IUGG and IERS for the U{polar radius of
                  curvature<https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>},
                  see C{rocPolar}.
        '''
        return self.R2  # if self._R2 is None else self._R2

    def degrees2m(self, deg, lat=0):
        '''Convert an angle to the distance along the equator or along
           a parallel at the given latitude.

           @arg deg: Angle (C{degrees}).
           @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

           @return: Angle (C{degrees}).

           @raise RangeError: Latitude B{C{lat}} outside valid range
                              and L{rangerrors} set to C{True}.

           @raise ValueError: Invalid B{C{deg}} or B{C{lat}}.
        '''
        return degrees2m(deg, radius=self.a, lat=lat)

    def distance2(self, lat0, lon0, lat1, lon1):
        '''Approximate the distance and (initial) bearing between two
           points based on the U{local, flat earth approximation
           <https://www.EdWilliams.org/avform.htm#flat>} aka U{Hubeny
           <https://www.OVG.AT/de/vgi/files/pdf/3781/>} formula.

           Suitable only for distances of several hundred Km or Miles
           and only between points not near-polar.

           @arg lat0: From latitude (C{degrees}).
           @arg lon0: From longitude (C{degrees}).
           @arg lat1: To latitude (C{degrees}).
           @arg lon1: To longitude (C{degrees}).

           @return: A L{Distance2Tuple}C{(distance, initial)} with C{distance}
                    in same units as this ellipsoid's axes.

           @note: The meridional and prime_vertical radii of curvature are
                  taken and scaled I{at the initial latitude}, see C{roc2}.

           @see: Function L{flatLocal}/L{hubeny}.
        '''
        phi0 = Phi_(lat0, name=_lat0_)
        m, n = self.roc2_(phi0, scaled=True)
        m *= Phi_(lat1, name=_lat1_) - phi0
        n *= Lam_(lon1, name=_lon1_) - Lam_(lon0, name=_lon0_)
        return Distance2Tuple(hypot(m, n), degrees360(atan2(n, m)))

    @property_RO
    def e(self):
        '''Get the (1st) eccentricity (C{float}), M{sqrt(1 - (b / a)**2))}.
        '''
        if self._e is None:
            self._e = Float(sqrt(abs(self.e2)) if self.f else 0, name=Ellipsoid.e.name)
        return self._e

    eccentricity = e  # eccentricity

    @property_RO
    def e2(self):
        '''Get the (1st) eccentricity I{squared} (C{float}), M{f * (2 - f) == 1 - (b / a)**2}.
        '''
        if self._e2 is None:
            self._e2 = self._assert(1 - (self.b / self.a)**2, e2=f2e2(self.f))
        return self._e2

    eccentricity2 = e2  # eccentricity squared

    @property_RO
    def e22(self):
        '''Get the 2nd eccentricity I{squared} (C{float}), M{e2 / (1 - e2) == (a / b)**2 - 1}.
        '''
        if self._e22 is None:
            self._e22 = self._assert(self.a_b**2 - 1, e22=f2e22(self.f))
        return self._e22

    eccentricity2nd2 = e2  # second eccentricity squared

    @property_RO
    def e32(self):
        '''Get the 3rd eccentricity I{squared} (C{float}), M{e2 / (2 - e2) == (a**2 - b**2) / (a**2 + b**2)}.
        '''
        if self._e32 is None:
            self._e32 = self._assert((self.a2 - self.b2) / (self.a2 + self.b2), e32=f2e32(self.f))
        return self._e32

    eccentricity3rd2 = e32  # third eccentricity squared

    @property_RO
    def e12(self):
        '''Get M{1 - e**2} (C{float}), M{(1 - f)**2}, M{b**2 / a**2}, see C{b2_a2}.
        '''
        if self._e12 is None:
            self._e12 = self._assert((1 - self.f)**2, e12=1 - self.e2, f0=1)
        return self._e12  # 1 - e2

    @property_RO
    def e4(self):
        '''Get the (1st) eccentricity to 4th power (C{float}), M{e**4 == e2**2}.
        '''
        if self._e4 is None:
            self._e4 = Float(self.e2**2 if self.f else 0, name=Ellipsoid.e4.name)
        return self._e4

    def ecef(self, Ecef=None):
        '''Return U{ECEF<https://WikiPedia.org/wiki/ECEF>} converter.

           @kwarg Ecef: ECEF class to use (L{EcefKarney}, L{EcefVeness}
                        or L{EcefYou}).

           @return: An ECEF converter for this C{ellipsoid} (L{EcefKarney},
                    L{EcefVeness} or L{EcefYou}).

           @raise TypeError: Invalid B{C{Ecef}}.
        '''
        return _4Ecef(self, Ecef)

    @property_RO
    def _elliptic_e22(self):
        '''(INTERNAL) Elliptic function for C{auxRectifying}, C{L}, C{Llat}.
        '''
        if self._elliptic is None:
            from pygeodesy.elliptic import Elliptic
            self._elliptic = Elliptic(-self.e22)
        return self._elliptic

    def e2s(self, s):
        '''Compute norm M{sqrt(1 - e2 * s**2)}.

           @arg s: S value (C{scalar}).

           @return: Norm (C{float}).

           @raise ValueError: Invalid B{C{s}}.
        '''
        return sqrt(self.e2s2(s))

    def e2s2(self, s):
        '''Compute M{1 - e2 * s**2}.

           @arg s: S value (C{scalar}).

           @return: Result (C{float}).

           @raise ValueError: Invalid B{C{s}}.
        '''
        try:
            r = 1 - self.e2 * Scalar(s, name='s')**2
            if r < 0:
                raise ValueError(_negative_)
        except (TypeError, ValueError) as x:
            raise _ValueError(self._dot_(Ellipsoid.e2s2.__name__), s, txt=str(x))
        return r

    @property_RO
    def es(self):
        '''Get the (1st) eccentricity I{signed} (C{float}).
        '''
        # note, self.e is always non-negative
        return Float(copysign(self.e, self.f), name=Ellipsoid.es.name)  # see .ups

    def es_atanh(self, x):
        '''Compute M{es * atanh(es * x)} where I{es} is the I{signed}
           (1st) eccentricity.

           @raise ValueError: Invalid B{C{x}}.

           @see: Function U{Math::eatanhe<https://GeographicLib.SourceForge.io/
                 html/classGeographicLib_1_1Math.html>}.
        '''
        # note, self.e is always non-negative
        if self.f > 0:  # oblate
            r = self.e * atanh(self.e * Scalar(x, name=_x_))
        elif self.f < 0:  # prolate
            r = -self.e * atan(self.e * Scalar(x, name=_x_))
        else:
            r = 0
        return r

    @property_RO
    def es_c(self):
        '''Get M{(1 - f) * exp(es_atanh(1))} (C{float}), M{b_a * exp(es_atanh(1))}.
        '''
        if self._es_c is None:
            self._es_c = Float(self.b_a * exp(self.es_atanh(1.0)) if self.f else 1, name=Ellipsoid.es_c.name)
        return self._es_c

    def es_tauf(self, taup):
        '''Compute I{Karney}'s U{equations (19), (20) and (21)
           <https://ArXiv.org/abs/1002.1417>}.

           @see: U{Math::tauf<https://GeographicLib.SourceForge.io/
                 html/classGeographicLib_1_1Math.html>}.
        '''
        T_  = Scalar(taup, name='taup')
        tol = max(abs(T_), 1) * _TOL
        e12 = self.e12
        t = T_ / e12
        T = Fsum(t)
        for _ in range(9):
            a = self.es_taupf(t)
            d = (T_ - a) * (1 + e12 * t**2) / (e12 * hypot1(t) * hypot1(a))
            t, d = T.fsum2_(d)
            if abs(d) < tol:
                break
        return t

    def es_taupf(self, tau):
        '''Compute I{Karney}'s U{equations (7), (8) and (9)
           <https://ArXiv.org/abs/1002.1417>}.

           @see: U{Math::taupf<https://GeographicLib.SourceForge.io/
                 html/classGeographicLib_1_1Math.html>}.
        '''
        T = Scalar(tau, name='tau')
        t = hypot1(T)
        s = sinh(self.es_atanh(T / t))
        return hypot1(s) * T - s * t

    @property_RO
    def f(self):
        '''Get the flattening (C{float}), M{(a - b) / a}, C{0} for spherical.
        '''
        return self._f

    flattening = f

    @property_RO
    def f_(self):
        '''Get the I{inverse} flattening (C{float}), M{1 / f} == M{a / (a - b)}, C{0} for spherical.
        '''
        return self._f_

    @property_RO
    def f2(self):
        '''Get the 2nd flattening (C{float}), M{(a - b) / b == f / (1 - f)}, C{0} for spherical.
        '''
        if self._f2 is None:
            self._f2 = self._assert(self.a_b - 1, f2=f2f2(self.f))
        return self._f2

    flattening2nd = f2

    @property_RO
    def geodesic(self):
        '''Get this ellipsoid's I{wrapped Karney} U{Geodesic
           <https://GeographicLib.SourceForge.io/html/python/code.html>},
           provided the U{geographiclib<https://PyPI.org/project/geographiclib>}
           package is installed.
        '''
        if self._geodesic is None:
            # if not self.isEllipsoidal:
            #     raise _IsnotError('ellipsoidal', ellipsoid=self)
            from pygeodesy.karney import _wrapped
            self._geodesic = _wrapped.Geodesic(self.a, self.f)
        return self._geodesic

    @property_RO
    def _geodesic_Math2(self):
        '''(INTERNAL) Get this ellipsoid's I{wrapped Karney} C{Geodesic}
           and I{Karney}'s C{Math} class.
        '''
        if Ellipsoid._Math is None:
            from pygeodesy.karney import _wrapped
            Ellipsoid._Math = _wrapped.Math
        return self.geodesic, Ellipsoid._Math

    def _hubeny2_(self, phi2, phi1, lam21):
        '''(INTERNAL) like function L{flatLocal_}/L{hubeny_} but
           returning the I{angular} distance in C{radians squared}.
        '''
        m, n = self.roc2_((phi2 + phi1) * 0.5, scaled=True)
        return hypot2(m * (phi2 - phi1), n * lam21) * self.a2_

    @property_RO
    def isEllipsoidal(self):
        '''Is this model I{ellipsoidal} (C{bool})?
        '''
        return self.f != 0

    @property_RO
    def isOblate(self):
        '''Is this ellipsoid I{oblate} (C{bool})?  I{Prolate} or
           spherical otherwise.
        '''
        return self.f > 0

    @property_RO
    def isProlate(self):
        '''Is this ellipsoid I{prolate} (C{bool})?  I{Oblate} or
           spherical otherwise.
        '''
        return self.f < 0

    @property_RO
    def isSpherical(self):
        '''Is this model I{spherical} (C{bool})?
        '''
        return self.f == 0

    @property_doc_(''' the I{Krüger} series' order (C{int}).''')
    def KsOrder(self):
        '''Get the Krüger series order (C{int} 4, 6 or 8).
        '''
        return self._KsOrder

    @KsOrder.setter  # PYCHOK setter!
    def KsOrder(self, order):
        '''Set the I{Krüger} series' order.

           @arg order: New I{Krüger} series' order (C{int} 4, 6 or 8).

           @raise ValueError: Invalid B{C{order}}.
        '''
        if order not in (4, 6, 8):
            raise _ValueError(order=order)
        if order != self._KsOrder:
            if self._AlphaKs:
                self._AlphaKs = None
            if self._BetaKs:
                self._BetaKs = None
            self._KsOrder = order

    @property_RO
    def L(self):
        '''Get the I{quarter meridian} C{L} (aka polar distance) the distance
           along a meridian between the equator and a pole (C{meter}),
           M{b * Elliptic(-e2 / (1 - e2)).E} or M{a * PI / 2}.
        '''
        if self._L is None:
            if self.f:  # complete integral 2nd ...
                # kind: Elliptic(-e2 / (1 - e2)).E
                L = self.b * self._elliptic_e22.cE
            else:  # spherical
                L = self.a * PI_2
            self._L = Distance(L, name=Ellipsoid.L.name)
        return self._L

    quarteradius = L  # -meridian
    '''DEPRECATED, use C{L}.'''

    def Llat(self, lat):
        '''Return the I{meridional length}, the distance along a meridian
           between the equator and the given latitude (C{Meter}).

           @arg lat: Geodetic latitude (C{degrees90}).

           @return: The meridional length at B{C{lat}}, negative on southern
                    hemisphere (C{meter}).
        '''
        L = self._elliptic_e22.fEd(self.auxParametric(lat)) if self.f else Phi_(lat)
        return Distance(self.b * L, name=Ellipsoid.Llat.__name__)

    Lmeridian = Llat  # meridional distance

    def m2degrees(self, meter, lat=0):
        '''Convert distance to angle along equator or along
           a parallel at an other latitude.

           @arg meter: Distance (C{meter}).
           @kwarg lat: Parallel latitude (C{degrees90}, C{str}).

           @return: Angle (C{degrees}).

           @raise RangeError: Latitude B{C{lat}} outside valid range
                              and L{rangerrors} set to C{True}.

           @raise ValueError: Invalid B{C{meter}} or B{C{lat}}.
       '''
        return m2degrees(meter, radius=self.a, lat=lat)

    @property_RO
    def Mabcd(self):
        '''Get the OSGR meridional coefficients (C{4-Tuple}), C{Airy130} only.
        '''
        if self._Mabcd is None:
            n, n2, n3 = fpowers(self.n, 3)  # PYCHOK false!
            self._Mabcd = (fdot((1, n, n2, n3), 4, 4, 5, 5) / 4,
                           fdot(   (n, n2, n3), 24, 24, 21) / 8,
                           fdot(      (n2, n3), 15, 15) / 8,
                                      35 * n3 / 24)
        return self._Mabcd

    majoradius = a
    '''DEPRECATED, use C{a} or C{Requatorial}.'''
    minoradius = b
    '''DEPRECATED, use C{b} or C{Rpolar}.'''

    @property_RO
    def n(self):
        '''Get the 3rd flattening (C{float}), M{f / (2 - f) == (a - b) / (a + b)}.
        '''
        if self._n is None:
            self._n = self._assert((self.a - self.b) / (self.a + self.b), n=f2n(self.f))
        return self._n

    flattening3rd = n

    @property_RO
    def R1(self):
        '''Get the I{mean} earth radius per IUGG (C{meter}), M{(2 * a + b) / 3}.

           @see: U{Earth radius<https://WikiPedia.org/wiki/Earth_radius>}
                 and method C{Rgeometric}.
        '''
        if self._R1 is None:
            self._R1 = Radius((self.a * 2 + self.b) / 3 if self.f else self.a, name=Ellipsoid.R1.name)
        return self._R1

    Rmean = R1

    @property_RO
    def R2(self):
        '''Get the I{authalic} earth radius (C{meter}), M{sqrt((a**2 + b**2 * atanh(e) / e) / 2)}.

           @see: U{Earth radius<https://WikiPedia.org/wiki/Earth_radius>}, C{area} and
                 U{c<https://Link.Springer.com/article/10.1007%2Fs00190-012-0578-z>}.
        '''
        if self._R2 is None:
            if self.f:
                r = self.b2_a2
                if self.f > 0:  # self.e2 > 0
                    r *= atanh(self.e) / self.e
                elif self.f < 0:  # self.e2 < 0
                    r *= atan(self.e) / self.e
                r = self.a * sqrt((1 + r) * 0.5)
            else:
                r = self.a
            self._R2 = Radius(r, name=Ellipsoid.R2.name)
        return self._R2

    Rauthalic = R2

    @property_RO
    def R3(self):
        '''Get the I{volumetric} earth radius (C{meter}), M{(a * a * b)**(1/3)}.

           @see: U{Earth radius<https://WikiPedia.org/wiki/Earth_radius>}.
        '''
        if self._R3 is None:
            self._R3 = Radius(cbrt(self.a2 * self.b) if self.f else self.a, name=Ellipsoid.R3.name)
        return self._R3

    Rvolumetric = R3

    @property_RO
    def Rbiaxial(self):
        '''Get the I{biaxial}, I{quadratic} mean earth radius (C{meter}), M{sqrt((a**2 + b**2) / 2)}.
        '''
        if self._Rb is None:
            self._Rb = Radius(sqrt((self.a2 + self.b2) / 2) if self.f else self.a, name=Ellipsoid.Rbiaxial.name)
        return self._Rb

    Requatorial = a  # for consistent naming

    def Rgeocentric(self, lat):
        '''Compute the I{geocentric} earth radius at the given latitude.

           @arg lat: Latitude (C{degrees90}).

           @return: Geocentric earth radius (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.

           @see: U{Geocentric Radius
                 <https://WikiPedia.org/wiki/Earth_radius#Geocentric_radius>}
        '''
        c2, s2 = _c2_s2(Phi_(lat))
        b2_a2_s2 = self.b2_a2 * s2
        # R_ == sqrt((a2**2 * c2 + b2**2 * s2) / (a2 * c2 + b2 * s2))
        #    == sqrt(a2**2 * (c2 + (b2 / a2)**2 * s2) / (a2 * (c2 + b2 / a2 * s2)))
        #    == sqrt(a2 * (c2 + (b2 / a2)**2 * s2) / (c2 + (b2 / a2) * s2))
        #    == a * sqrt((c2 + b2_a2 * b2_a2 * s2) / (c2 + b2_a2 * s2))
        #    == a * sqrt((c2 + b2_a2 * b2_a2_s2) / (c2 + b2_a2_s2))
        return Radius(self.a * sqrt((c2 + self.b2_a2 * b2_a2_s2) / (c2 + b2_a2_s2)) if self.f else self.a,
                      name=Ellipsoid.Rgeocentric.__name__)

    @property_RO
    def Rgeometric(self):
        '''Get the I{geometric} mean earth radius (C{meter}), M{sqrt(a * b)}.

           @see: Method C{R1}.
        '''
        if self._Rg is None:
            self._Rg = Radius(sqrt(self.a * self.b) if self.f else self.a, name=Ellipsoid.Rgeometric.name)
        return self._Rg

    Rs = Rgeometric  # for backward compatibility
    '''DEPRECATED, use C{Rgeometric}.'''
    Rpolar = b  # for consistent naming

    def Rlat(self, lat):
        '''Approximate the earth radius at the given latitude.

           @arg lat: Latitude (C{degrees90}).

           @return: Approximate earth radius (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.
        '''
        if self._ab_90 is None:
            self._ab_90 = (self.a - self.b) / 90.0
        # r = a - (a - b) * |lat| / 90
        r = self.a
        if lat:
            r -= self._ab_90 * min(abs(Lat(lat, clip=0)), 90)
        return Radius(r, name=Ellipsoid.Rlat.__name__)

    Rquadratic = Rbiaxial  # synonyms

    @property_RO
    def Rrectifying(self):
        '''Get the I{rectifying} earth radius (C{meter}), M{((a**(3/2) + b**(3/2)) / 2)**(2/3)}.

           @see: U{Earth radius<https://WikiPedia.org/wiki/Earth_radius>}.
        '''
        if self._Rr is None:
            self._Rr = Radius(cbrt2((sqrt3(self.a) + sqrt3(self.b)) * 0.5) if self.f else self.a,
                              name=Ellipsoid.Rrectifying.name)
        return self._Rr

    Rr = Rrectifying
    '''DEPRECATED, use C{Rrectifying}.'''

    def roc2(self, lat, scaled=False):
        '''Compute the I{meridional} and I{prime-vertical}, I{normal}
           radii of curvature at the given latitude.

           @arg lat: Latitude (C{degrees90}).
           @kwarg scaled: Scale prime_vertical by B{C{cos(radians(lat))}} (C{bool}).

           @return: An L{Curvature2Tuple}C{(meridional, prime_vertical)} with
                    the radii of curvature.

           @raise ValueError: Invalid B{C{lat}}.

           @see: Method C{roc2_}, U{Local, flat earth approximation
                 <https://www.EdWilliams.org/avform.htm#flat>} and
                 U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        return self.roc2_(Phi_(lat), scaled=scaled)

    def roc2_(self, phi, scaled=False):
        '''Compute the I{meridional} and I{prime-vertical}, I{normal}
           radii of curvature at the given latitude.

           @arg phi: Latitude (C{radians}).
           @kwarg scaled: Scale prime_vertical by B{C{cos(phi)}} (C{bool}).

           @return: An L{Curvature2Tuple}C{(meridional, prime_vertical)} with
                    the radii of curvature.

           @raise ValueError: Invalid B{C{phi}}.

           @see: Method C{roc2}, U{Local, flat earth approximation
                 <https://www.EdWilliams.org/avform.htm#flat>} and meridional
                 and prime vertical U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        a = abs(Phi(phi))
        r = self.e2s2(sin(a) if a < PI_2 else 1)
        if r < EPS:
            m = n = 0  # PYCHOK attr
        elif r < EPS1:
            n = self.a / sqrt(r)
            m = n * self.e12 / r  # PYCHOK attr
        else:
            n = self.a
            m = n * self.e12
        if scaled:
            n *= cos(a) if a < PI_2 else 0
        return Curvature2Tuple(Radius(m, name=Ellipsoid.rocMeridional.__name__),
                               Radius(n, name=Ellipsoid.rocPrimeVertical.__name__))

    def rocBearing(self, lat, bearing):
        '''Compute the I{directional} radius of curvature at a
           given latitude and compass direction.

           @arg lat: Latitude (C{degrees90}).
           @arg bearing: Direction (compass C{degrees360}).

           @return: Directional radius of curvature (C{meter}).

           @raise RangeError: Latitude B{C{lat}} outside valid range
                              and L{rangerrors} set to C{True}.

           @raise ValueError: Invalid B{C{lat}} or B{C{bearing}}.

           @see: U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}
        '''
        c2, s2 = _c2_s2(Bearing_(bearing))
        m, n = self.roc2_(Phi_(lat))
        if n < m:  # == n / (c2 * n / m + s2)
            c2 *= n / m
        elif m < n:  # == m / (c2 + s2 * m / n)
            s2 *= m / n
            n = m
        return Radius(n / (c2 + s2), name=Ellipsoid.rocBearing.__name__)  # == 1 / (c2 / m + s2 / n)

    def rocGauss(self, lat):
        '''Compute the I{Gaussian} radius of curvature at the given latitude.

           @arg lat: Latitude (C{degrees90}).

           @return: Gaussian radius of curvature (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.

           @see: U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}
        '''
        # using ...
        #    m, n = self.roc2_(Phi_(lat))
        #    return sqrt(m * n)
        # ... requires 1 or 2 sqrt
        c2, s2 = _c2_s2(Phi_(lat))
        return Radius(self.b / (c2 + self.b2_a2 * s2), name=Ellipsoid.rocGauss.__name__)

    def rocMean(self, lat):
        '''Compute the I{mean} radius of curvature at the given latitude.

           @arg lat: Latitude (C{degrees90}).

           @return: Mean radius of curvature (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.

           @see: U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}
        '''
        m, n = self.roc2_(Phi_(lat))
        return Radius(2 * m * n / (m + n), name=Ellipsoid.rocMean.__name__)  # == 2 / (1 / m + 1 / n)

    def rocMeridional(self, lat):
        '''Compute the I{meridional} radius of curvature at the given latitude.

           @arg lat: Latitude (C{degrees90}).

           @return: Meridional radius of curvature (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.

           @see: U{Local, flat earth approximation
                 <https://www.EdWilliams.org/avform.htm#flat>} and
                 U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        return self.roc2_(Phi_(lat)).meridional

    @property_RO
    def rocPolar(self):
        '''Get the I{polar} radius of curvature (C{meter}), see C{a2_b}.
        '''
        return self.a2_b if self._a2_b is None else self._a2_b

    def rocPrimeVertical(self, lat):
        '''Compute the I{prime-vertical}, I{normal} radius of curvature at
           the given latitude, aka the transverse radius of curvature.

           @arg lat: Latitude (C{degrees90}).

           @return: Prime-vertical radius of curvature (C{meter}).

           @raise ValueError: Invalid B{C{lat}}.

           @see: U{Local, flat earth approximation
                 <https://www.EdWilliams.org/avform.htm#flat>} and
                 U{Radii of Curvature
                 <https://WikiPedia.org/wiki/Earth_radius#Radii_of_curvature>}.
        '''
        return self.roc2_(Phi_(lat)).prime_vertical

    rocTransverse = rocPrimeVertical  # synonyms

    def toStr(self, prec=8):  # PYCHOK expected
        '''Return this ellipsoid as a text string.

           @kwarg prec: Optional number of decimals, unstripped (C{int}).

           @return: Ellipsoid attributes (C{str}).
        '''
        return self._instr(prec, 'a', 'b', 'f_', 'f', 'f2', 'n',
                                 'e', 'e2', 'e22', 'e32', 'L',
                                 'R1', 'R2', 'R3')

    @property_RO
    def volume(self):
        '''Get the ellipsoid's I{volume} (C{meter**3}), M{4 / 3 * PI * a**2 * b}.
        '''
        if self._volume is None:
            self._volume = Float(_PI4_3 * self.a2 * self.b, name=Ellipsoid.volume.name)
        return self._volume


class Ellipsoid2(Ellipsoid):
    '''Like L{Ellipsoid}, but specified by I{equatorial} radius and I{flattening}.
    '''
    def __init__(self, a, f, name=NN):
        '''New L{Ellipsoid2}.

           @arg a: Equatorial (major) radius, semi-axis (C{meter}).
           @arg f: Flattening: (C{float} < 1.0).
           @kwarg name: Optional, unique name (C{str}).

           @raise NameError: Ellipsoid with that B{C{name}} already exists.

           @raise ValueError: Invalid B{C{a}} or B{C{f}}.

           @note: M{abs(f) < EPS} is forced to M{f = 0}, spherical.
        '''
        try:
            a = Float_(a, low=EPS, high=None)  # like Radius_
            f = Float_(f, low=None, high=EPS1)
            Ellipsoid.__init__(self, a, a_f2b(a, f), name=name)
        except (TypeError, ValueError) as x:
            raise _ValueError(instr(self, a=a, f=f), txt=str(x))


def a_b2f(a, b):
    '''Return C{f}, the I{flattening} for a given I{equatorial} and I{polar} radius.

       @arg a: Equatorial (major) radius (C{scalar} > 0).
       @arg b: Polar (minor) radius (C{scalar} > 0).

       @return: The flattening (C{float} or C{0}).
    '''
    f = 0 if a < EPS or b < EPS else (a - b) / a
    return Float(0 if abs(f) < EPS else f, name=Ellipsoid.f.name)


def a_b2f_(a, b):
    '''Return C{f_}, the I{inverse flattening} for a given I{equatorial} and I{polar} radius.

       @arg a: Equatorial (major) radius (C{scalar} > 0).
       @arg b: Polar (minor) radius (C{scalar} > 0).

       @return: The inverse flattening (C{float} or C{0}).
    '''
    f_ = 0 if b < EPS or a < EPS else (a - b) / b
    return Float(0 if abs(f_) < EPS else f_, name=Ellipsoid.f_.name)


def a_f2b(a, f):
    '''Return C{b}, the I{polar} radius for a given I{equatorial} radius and I{flattening}.

       @arg a: Equatorial (major) radius (C{scalar} > 0).
       @arg f: Flattening (C{scalar} < 1).

       @return: The polar (minor) radius (C{float}).
    '''
    return Float(a if abs(f) < EPS else a * (1 - f), name=Ellipsoid.b.name)


def a_f_2b(a, f_):
    '''Return C{b}, the I{polar} radius for a given I{equatorial} radius and I{inverse flattening}.

       @arg a: Equatorial (major) radius (C{scalar} > 0).
       @arg f_: Inverse flattening (C{scalar} >>> 1).

       @return: The polar (minor) radius (C{float}).
    '''
    return Float(a if abs(f_) < EPS else a * (f_ - 1) / f_, name=Ellipsoid.b.name)


def b_f2a(b, f):
    '''Return C{a}, the I{equatorial} radius for a given I{polar} radius and I{flattening}.

       @arg b: Polar (minor) radius (C{scalar} > 0).
       @arg f: Flattening (C{scalar} < 1).

       @return: The equatorial (major) radius (C{float}).
    '''
    f = 1 - f
    return Float(b if abs(f) < EPS else b / f, name=Ellipsoid.a.name)


def b_f_2a(b, f_):
    '''Return C{a}, the I{equatorial} radius for a given I{polar} radius and I{inverse flattening}.

       @arg b: Polar (minor) radius (C{scalar} > 0).
       @arg f_: Inverse flattening (C{scalar} >>> 1).

       @return: The equatorial (major) radius (C{float}).
    '''
    return Float(b if abs(f_) < EPS else b * (f_ + 1), name=Ellipsoid.a.name)


def f2e2(f):
    '''Return C{e2}, the I{eccentricity squared} for a given I{flattening}.

       @arg f: Flattening (C{scalar} < 1).

       @return: The (1st) eccentricity I{squared} (C{float} < 1).

       @note: The result C{e2} will be negative for I{prolate} ellipsoids.

       @see: U{Eccentricity conversions<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1Ellipsoid.html>} and U{Flattening
             <https://WikiPedia.org/wiki/Flattening>}.
    '''
    return Float(f * (2 - f), name=Ellipsoid.e2.name)


def f2e22(f):
    '''Return C{e22}, the I{2nd eccentricity squared} for a given I{flattening}.

       @arg f: Flattening (C{scalar} < 1).

       @return: The 2nd eccentricity I{squared} (C{float} > -1 or C{INF}).

       @see: U{Eccentricity conversions<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1Ellipsoid.html>}.
    '''
    # e2 / (1 - e2) == f * (2 - f) / (1 - f)**2
    t = (1 - f)**2
    return Float(INF if t < EPS else f2e2(f) / t, name=Ellipsoid.e22.name)  # PYCHOK type


def f2e32(f):
    '''Return C{e32}, the I{3rd eccentricity squared} for a given I{flattening}.

       @arg f: Flattening (C{scalar} < 1).

       @return: The 3rd eccentricity I{squared} (C{float}).

       @see: U{Eccentricity conversions<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1Ellipsoid.html>}.
    '''
    # e2 / (2 - e2) == f * (2 - f) / (1 + (1 - f)**2)
    e2 = f2e2(f)
    return Float(e2 / (2 - e2), name=Ellipsoid.e32.name)


def f_2f(f_):
    '''Return C{f}, the I{flattening} for a given I{inverse flattening}.

       @arg f_: Inverse flattening (C{scalar} >>> 1).

       @return: The flattening (C{float} or C{0}).
    '''
    return Float(0 if abs(f_) < EPS else 1.0 / f_, name=Ellipsoid.f.name)  # PYCHOK type


def f2f_(f):
    '''Return C{f_}, the I{inverse flattening} for a given I{flattening}.

       @arg f: Flattening (C{scalar} < 1).

       @return: The inverse flattening (C{float} or C{0}).
    '''
    return Float(0 if abs(f) < EPS else 1.0 / f, name=Ellipsoid.f_.name)  # PYCHOK type


def f2f2(f):
    '''Return C{f2}, the I{2nd flattening} for a given I{flattening}.

       @arg f: Flattening (C{scalar} < 1).

       @return: The 2nd flattening (C{float} or C{INF}).

       @see: U{Eccentricity conversions<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1Ellipsoid.html>} and U{Flattening
             <https://WikiPedia.org/wiki/Flattening>}.
    '''
    t = 1 - f
    return Float(INF if abs(t) < EPS else f / t, name=Ellipsoid.f2.name)  # PYCHOK type


def f2n(f):
    '''Return C{n}, the I{3rd flattening} for a given I{flattening}.

       @arg f: Flattening (C{scalar} < 1).

       @return: The 3rd flattening (-1 < C{float} < 1).

       @see: U{Eccentricity conversions<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1Ellipsoid.html>} and U{Flattening
             <https://WikiPedia.org/wiki/Flattening>}.
    '''
    return Float(f / (2 - f) if f else 0, name=Ellipsoid.n.name)


def n2e2(n):
    '''Return C{e2}, the I{eccentricity squared} for a given I{3rd flattening}.

       @arg n: The 3rd flattening (-1 < C{scalar} < 1).

       @return: The (1st) eccentricity I{squared} (C{float} or -INF).

       @note: The result C{e2} will be negative for I{prolate} ellipsoids.

       @see: U{Flattening<https://WikiPedia.org/wiki/Flattening>}.
    '''
    t = (1 + n)**2
    return Float(-INF if t < EPS else 4 * n / t, name=Ellipsoid.e2.name)


def n2f(n):
    '''Return C{f}, the I{flattening} for a given I{3rd flattening}.

       @arg n: The 3rd flattening (-1 <= C{scalar} <= 1).

       @return: The flattening (C{float} or -INF).

       @see: U{Eccentricity conversions<https://GeographicLib.SourceForge.io/
             html/classGeographicLib_1_1Ellipsoid.html>} and U{Flattening
             <https://WikiPedia.org/wiki/Flattening>}.
    '''
    t = 1 + n
    return Float(-INF if t < EPS else (0 if abs(n) < EPS else 2 * n / t), name=Ellipsoid.f.name)


Ellipsoids = _NamedEnum('Ellipsoids', Ellipsoid)  #: Registered ellipsoids.
# <https://www.GNU.org/software/gama/manual/html_node/Supported-ellipsoids.html>
# <https://w3.Energistics.org/archive/Epicentre/Epicentre_v3.0/DataModel/
#         LogicalDictionary/StandardValues/ellipsoid.html>
# <https://kb.OSU.edu/dspace/handle/1811/77986>
Ellipsoids._assert(  # <https://WikiPedia.org/wiki/Earth_ellipsoid>
    Airy1830       = Ellipsoid(6377563.396, None,              299.3249646,   'Airy1830'),   # b=6356256.909
    AiryModified   = Ellipsoid(6377340.189, None,              299.3249646,   'AiryModified'),  # b=6356034.448
#   ANS            = Ellipsoid(6378160.0,   None,              298.25,        'ANS'),  # b=6356774.719
    Australia1966  = Ellipsoid(6378160.0,   None,              298.25,        'Australia1966'),  # b=6356774.719
#   Bessel1841     = Ellipsoid(6377397.155, 6356078.963,       299.152815351, 'Bessel1841'),
    Bessel1841     = Ellipsoid(6377397.155, 6356078.962818,    299.1528128,   'Bessel1841'),
    Clarke1866     = Ellipsoid(6378206.4,   6356583.8,         294.978698214, 'Clarke1866'),
    Clarke1880     = Ellipsoid(6378249.145, 6356514.86954978,  293.465,       'Clarke1880'),
    Clarke1880IGN  = Ellipsoid(6378249.2,   6356515.0,         293.466021294, 'Clarke1880IGN'),
    Clarke1880Mod  = Ellipsoid(6378249.145, 6356514.96582849,  293.4663,      'Clarke1880Mod'),
    CPM1799        = Ellipsoid(6375738.7,   6356671.92557493,  334.39,        'CPM1799'),  # Comm. des Poids et Mesures
    Delambre1810   = Ellipsoid(6376428.0,   6355957.92616372,  311.5,         'Delambre1810'),  # Belgium
    Engelis1985    = Ellipsoid(6378136.05,  6356751.32272154,  298.2566,      'Engelis1985'),
    Everest1969    = Ellipsoid(6377295.664, 6356094.667915,    300.8017,      'Everest1969'),
    Fisher1968     = Ellipsoid(6378150.0,   6356768.33724438,  298.3,         'Fisher1968'),
    GEM10C         = Ellipsoid(6378137.0,   6356752.31424783,  298.2572236,   'GEM10C'),
    GRS67          = Ellipsoid(6378160.0,   None,              298.247167427, 'GRS67'),  # Lucerne b=6356774.516
    GRS80          = Ellipsoid(6378137.0,   6356752.314140347, 298.257222101, 'GRS80'),  # ITRS, ETRS89
    Helmert1906    = Ellipsoid(6378200.0,   6356818.16962789,  298.3,         'Helmert1906'),
    IERS1989       = Ellipsoid(6378136.0,   None,              298.257,       'IERS1989'),  # b=6356751.302
    IERS1992TOPEX  = Ellipsoid(6378136.3,   6356751.61659215,  298.257223563, 'IERS1992TOPEX'),  # IERS/TOPEX/Poseidon/McCarthy
    IERS2003       = Ellipsoid(6378136.6,   6356751.85797165,  298.25642,     'IERS2003'),
    Intl1924       = Ellipsoid(6378388.0,   None,              297.0,         'Intl1924'),  # aka Hayford b=6356911.946
    Intl1967       = Ellipsoid(6378157.5,   6356772.2,         298.24961539,  'Intl1967'),  # New Int'l
    Krassovski1940 = Ellipsoid(6378245.0,   6356863.01877305,  298.3,         'Krassovski1940'),  # spelling
    Krassowsky1940 = Ellipsoid(6378245.0,   6356863.01877305,  298.3,         'Krassowsky1940'),  # spelling
    Maupertuis1738 = Ellipsoid(6397300.0,   6363806.28272251,  191.0,         'Maupertuis1738'),  # France
    Mercury1960    = Ellipsoid(6378166.0,   6356784.28360711,  298.3,         'Mercury1960'),
    Mercury1968Mod = Ellipsoid(6378150.0,   6356768.33724438,  298.3,         'Mercury1968Mod'),
    NWL1965        = Ellipsoid(6378145.0,   6356759.76948868,  298.25,        'NWL1965'),  # Naval Weapons Lab.
    OSU86F         = Ellipsoid(6378136.2,   6356751.51693008,  298.2572236,   'OSU86F'),
    OSU91A         = Ellipsoid(6378136.3,   6356751.6165948,   298.2572236,   'OSU91A'),
#   Plessis1817    = Ellipsoid(6397523.0,   6355863.0,         153.56512242,  'Plessis1817'),  # XXX incorrect?
    Plessis1817    = Ellipsoid(6376523.0,   6355862.93325557,  308.64,        'Plessis1817'),  # XXX IGN France 1972
    SGS85          = Ellipsoid(6378136.0,   6356751.30156878,  298.257,       'SGS85'),  # Soviet Geodetic System
    SoAmerican1969 = Ellipsoid(6378160.0,   6356774.71919531,  298.25,        'SoAmerican1969'),  # South American
    Struve1860     = Ellipsoid(6378298.3,   6356657.14266956,  294.73,        'Struve1860'),
    WGS60          = Ellipsoid(6378165.0,   6356783.28695944,  298.3,         'WGS60'),
    WGS66          = Ellipsoid(6378145.0,   6356759.76948868,  298.25,        'WGS66'),
    WGS72          = Ellipsoid(6378135.0,   None,              298.26,        'WGS72'),  # b=6356750.52
    WGS84          = Ellipsoid(6378137.0,   None,              298.257223563, 'WGS84'),  # GPS b=6356752.31425
#   Prolate        = Ellipsoid(6356752.3,   6378137.0,         None,          'Prolate'),
    Sphere         = Ellipsoid(R_M,         R_M,                 0.0,         'Sphere'),  # pseudo
    SphereAuthalic = Ellipsoid(R_FM,        R_FM,                0.0,         'SphereAuthalic'),  # pseudo
    SpherePopular  = Ellipsoid(R_MA,        R_MA,                0.0,         'SpherePopular'),  # EPSG:3857 Spheroid
)

__all__ += _ALL_DOCS(Curvature2Tuple)

if __name__ == '__main__':

    for E in (Ellipsoids.WGS84, Ellipsoids.GRS80,  # NAD83,
              Ellipsoids.Sphere, Ellipsoids.SpherePopular,
              Ellipsoid(6356752.3, 6378137.0, name='_Prolate')):
        if E.f:
            f_ = 'f_=1/%.10F' % (1 / E.f,)
        else:
            f_ = 'f_=' + _n_a_
        e = (E.a - E.b) / (E.a + E.b) - E.n
        t = (E.toStr(prec=10),
            'A=%r, e=%s, %s, n=%s(%s)' % (E.A, fstr(E.e, prec=13, fmt=_e),
                                          f_,  fstr(E.n, prec=13, fmt=_Fmt),
                                               fstr(e,   prec=3,  fmt=_e),),
            '%s=(%s)'   % (Ellipsoid.AlphaKs.name, fstr(E.AlphaKs, prec=20),),
            '%s= (%s)'  % (Ellipsoid.BetaKs.name,  fstr(E.BetaKs,  prec=20),),
            '%s= %s'    % ('KsOrder',                   E.KsOrder),
            '%s=  (%s)' % (Ellipsoid.Mabcd.name,   fstr(E.Mabcd,   prec=20),))
        print('\nEllipsoid.%s: %s' % (E.name, ',\n    '.join(t)))

    # __doc__ of this file
    for e in (Ellipsoids,):
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

# % python -m pygeodesy.ellipsoids

# Ellipsoid.WGS84: name='WGS84', a=6378137, b=6356752.3142451793, f_=298.257223563, f=0.0033528107, f2=0.0033640898, n=0.0016792204, e=0.0818191908, e2=0.00669438, e22=0.0067394967, e32=0.0033584313, L=10001965.7293127235, R1=6371008.7714150595, R2=6371007.1809184738, R3=6371000.790009154,
#     A=6367449.145823414, e=8.1819190842622e-02, f_=1/298.2572235630, n=0.0016792203864(-2.168e-19),
#     AlphaKs=(0.00083773182062447786, 0.00000076085277735726, 0.00000000119764550324, 0.00000000000242917068, 0.00000000000000571182, 0.0000000000000000148, 0.00000000000000000004, 0.0),
#     BetaKs= (0.00083773216405795667, 0.0000000590587015222, 0.00000000016734826653, 0.00000000000021647981, 0.00000000000000037879, 0.00000000000000000072, 0.0, 0.0),
#     KsOrder= 8,
#     Mabcd=  (1.00168275103155868244, 0.00504613293193333871, 0.00000529596776243457, 0.00000000690525779769)

# Ellipsoid.GRS80: name='GRS80', a=6378137, b=6356752.3141403468, f_=298.257222101, f=0.0033528107, f2=0.0033640898, n=0.0016792204, e=0.081819191, e2=0.00669438, e22=0.0067394968, e32=0.0033584313, L=10001965.7292304561, R1=6371008.7713801153, R2=6371007.1808835138, R3=6371000.7899741307,
#     A=6367449.145771043, e=8.1819191042833e-02, f_=1/298.2572221009, n=0.0016792203946(0.0e+00),
#     AlphaKs=(0.00083773182472890429, 0.00000076085278481561, 0.00000000119764552086, 0.00000000000242917073, 0.00000000000000571182, 0.0000000000000000148, 0.00000000000000000004, 0.0),
#     BetaKs= (0.0008377321681623882, 0.00000005905870210374, 0.000000000167348269, 0.00000000000021647982, 0.00000000000000037879, 0.00000000000000000072, 0.0, 0.0),
#     KsOrder= 8,
#     Mabcd=  (1.00168275103983916985, 0.0050461329567537995, 0.00000529596781448937, 0.00000000690525789941)

# Ellipsoid.Sphere: name='Sphere', a=6371008.7714149999, b=6371008.7714149999, f_=0, f=0, f2=0, n=0, e=0, e2=0, e22=0, e32=0, L=10007557.1761167478, R1=6371008.7714149999, R2=6371008.7714149999, R3=6371008.7714149999,
#     A=6371008.771415, e=0.0e+00, f_=n/a, n=0.0(0.0e+00),
#     AlphaKs=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
#     BetaKs= (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
#     KsOrder= 8,
#     Mabcd=  (1.0, 0.0, 0.0, 0.0)

# Ellipsoid.SpherePopular: name='SpherePopular', a=6378137, b=6378137, f_=0, f=0, f2=0, n=0, e=0, e2=0, e22=0, e32=0, L=10018754.171394622, R1=6378137, R2=6378137, R3=6378137,
#     A=6378137.0, e=0.0e+00, f_=n/a, n=0.0(0.0e+00),
#     AlphaKs=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
#     BetaKs= (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
#     KsOrder= 8,
#     Mabcd=  (1.0, 0.0, 0.0, 0.0)

# Ellipsoid._Prolate: name='_Prolate', a=6356752.2999999998, b=6378137, f_=-297.2570248823, f=-0.0033640921, f2=-0.0033528129, n=-0.0016792215, e=0.0820944654, e2=-0.0067395013, e22=-0.0066943844, e32=-0.0033584335, L=10001965.718133986, R1=6363880.5333333341, R2=6363878.9318593405, R3=6363872.5548945786,
#     A=6367449.138706811, e=8.2094465430912e-02, f_=1/-297.2570248823, n=-0.0016792215069(4.337e-19),
#     AlphaKs=(-0.00084149208789184528, 0.00000076653583100163, -0.00000000120934745864, 0.00000000000245762907, -0.00000000000000578865, 0.00000000000000001502, -0.00000000000000000004, 0.0),
#     BetaKs= (-0.00084149243499239313, 0.00000005842742972611, -0.00000000016804906023, 0.00000000000021706319, -0.00000000000000038002, 0.00000000000000000073, -0.0, 0.0),
#     KsOrder= 8,
#     Mabcd=  (0.99832429730542537882, -0.00502921759545196367, 0.00000527821842317521, -0.00000000690527162046)
