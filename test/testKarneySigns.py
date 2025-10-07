
# -*- coding: utf-8 -*-

# Sign and some geodesic tests from Karney's Python U{test/test_sign.py
# <https://GitHub.com/geographiclib/geographiclib-python/tree/main/geographiclib/>}.

__all__ = ('Tests',)
__version__ = '25.10.05'

from bases import TestsBase

from pygeodesy import EPS, INF, NAN, signBit  # atan2d, sincos2d
from pygeodesy.geodesicx import GeodesicExact as Geodesic
from pygeodesy.karney import _around, _atan2d, _diff182, \
                             _norm180, _sincos2d, _sum2

g_WGS84 = Geodesic()


def _sign(x):
    return '-' if signBit(x) else '+'


class Tests(TestsBase):
    """Sign test suite"""

    def _T2(self, which, *tuples):
        i = 0
        self.test(which.__name__, '...', '...', nl=1)
        for t in tuples:
            i += 1
            yield ('test_' + str(i)), t

    def test_eqv(self, n, y, x, **kwds):
        if not self.test(n, y, x, **kwds):  # pass
            if signBit(y) != signBit(x):
                self.test(n, _sign(y), _sign(x))

    def test_AngDiff(self):
        """Test special cases of AngDiff"""
        E128 = EPS * 128
        for n, (a, b, x) in self._T2(self.test_AngDiff,
               (  +0.0,   +0.0, +0.0),
               (  +0.0,   -0.0, -0.0),
               (  -0.0,   +0.0, +0.0),
               (  -0.0,   -0.0, +0.0),
               (  +5.0, +365.0, +0.0),
               (+365.0,   +5.0, -0.0),
               (+  5.0, +185.0, +180.0),
               (+185.0,   +5.0, -180.0),
               ( EPS  , +180.0, +180.0),
               (-EPS  , +180.0, -180.0),
               ( EPS  , -180.0, +180.0),
               (-EPS  , -180.0, -180.0),
               (E128+138, -164, 58-E128)):
            d, _ = _diff182(a, b)
            self.test_eqv(n, d, x)

    def test_AngNormalize(self):
        """Test special cases of AngNormalize"""
        for n, (d, x) in self._T2(self.test_AngNormalize,
               (-900.0, -180.0),
               (-720.0,   -0.0),
               (-540.0, -180.0),
               (-360.0,   -0.0),
               (-180.0, -180.0),
               (  -0.0,   -0.0),
               (  +0.0,   +0.0),
               ( 180.0, +180.0),
               ( 360.0,   +0.0),
               ( 540.0, +180.0),
               ( 720.0,   +0.0),
               ( 900.0, +180.0)):
            self.test_eqv(n, _norm180(d), x)

    def test_AngRound(self):
        """Test special cases for AngRound"""
        for n, (a, x) in self._T2(self.test_AngRound,
               (      EPS/32,       EPS/32),
               (     -EPS/32,      -EPS/32),
               (     -EPS/64, -0.0        ),  # PYCHOK E202
               (        -0.0, -0.0        ),  # PYCHOK E202
               (         0.0, +0.0        ),  # PYCHOK E202
               (      EPS/64, +0.0        ),  # PYCHOK E202
               ((1-2*EPS)/64, (1-2*EPS)/64),
               ((1-EPS  )/64,  1.0     /64),  # PYCHOK E202
               ((1-EPS/2)/64,  1.0     /64),
               ((1-EPS/4)/64,  1.0     /64),
               ( 1.0     /64,  1.0     /64),
               ((1+EPS/2)/64,  1.0     /64),
               ((1+EPS  )/64,  1.0     /64),  # PYCHOK E202
               ((1+2*EPS)/64, (1+2*EPS)/64),
               ((1-EPS  )/32, (1-EPS  )/32),  # PYCHOK E202
               ((1-EPS/2)/32,  1.0     /32),
               ((1-EPS/4)/32,  1.0     /32),
               ( 1.0     /32,  1.0     /32),
               ((1+EPS/2)/32,  1.0     /32),
               ((1+EPS  )/32, (1+EPS  )/32),  # PYCHOK E202
               ((1-EPS  )/16, (1-EPS  )/16),  # PYCHOK E202
               ((1-EPS/2)/16, (1-EPS/2)/16),
               ((1-EPS/4)/16,  1.0     /16),
               ( 1.0     /16,  1.0     /16),
               ((1+EPS/4)/16,  1.0     /16),
               ((1+EPS/2)/16,  1.0     /16),
               ((1+EPS  )/16, (1+EPS  )/16),  # PYCHOK E202
               ((1-EPS  )/ 8, (1-EPS  )/ 8),  # PYCHOK E202
               ((1-EPS/2)/ 8, (1-EPS/2)/ 8),
               ((1-EPS/4)/ 8,  1.0     / 8),
               ((1+EPS/2)/ 8,  1.0     / 8),
               ((1+EPS  )/ 8, (1+EPS  )/ 8),  # PYCHOK E202
               ( 1-EPS      ,  1-EPS      ),  # PYCHOK E202
               ( 1-EPS/2    ,  1-EPS/2    ),  # PYCHOK E202
               ( 1-EPS/4    ,  1.0        ),  # PYCHOK E202
               ( 1.0        ,  1.0        ),  # PYCHOK E202
               ( 1+EPS/4    ,  1.0        ),  # PYCHOK E202
               ( 1+EPS/2    ,  1.0        ),  # PYCHOK E202
               ( 1+EPS      ,  1+  EPS    ),  # PYCHOK E202
               ( 90.0-64*EPS,  90-64*EPS  ),  # PYCHOK E202
               ( 90.0-32*EPS,  90.0       ),  # PYCHOK E202
               ( 90.0       ,  90.0       )):  # PYCHOK E202
            self.test_eqv(n, _around(a), x)

    def test_antipodal(self):
        """How does the exact antipodal equatorial path go N/S + E/W"""
        for n, (lat1, lat2, lon2, azi1, azi2) in self._T2(self.test_antipodal,
               (+0.0, +0.0, +180,   +0.0, +180.0),
               (-0.0, -0.0, +180, +180.0,   +0.0),
               (+0.0, +0.0, -180,   -0.0, -180.0),
               (-0.0, -0.0, -180, -180.0,   -0.0)):
            r = g_WGS84.Inverse(lat1, 0.0, lat2, lon2)
            self.test_eqv(n, r.azi1, azi1)
            self.test_eqv(n, r.azi2, azi2)

    def test_antipodal_prolate(self):
        """Antipodal points on the equator with prolate ellipsoid"""
        g = Geodesic(6.4e6, -1 / 300.0)
        for n, (lon2, azi) in self._T2(self.test_antipodal_prolate,
               (+180, +90.0),
               (-180, -90.0)):
            r = g.Inverse(0.0, 0.0, 0.0, lon2)
            self.test_eqv(n, r.azi1, azi)
            self.test_eqv(n, r.azi2, azi)

    def test_azimuth_0_180(self):
        """Azimuths = +/-0 and +/-180 for the direct problem"""
        for n, (azi1, lon2, azi2) in self._T2(self.test_azimuth_0_180,
               (+0.0, +180.0, +180.0),
               (-0.0, -180.0, -180.0),
               (+180, +180.0,   +0.0),
               (-180, -180.0,   -0.0)):
            r = g_WGS84.Direct(0.0, 0.0, azi1, 15e6,
                               Geodesic.STANDARD | Geodesic.LONG_UNROLL)
            self.test_eqv(n, r.lon2, lon2)
            self.test_eqv(n, r.azi2, azi2)

    def test_equatorial_coincident(self):
        """Azimuth with coincident point on equator"""
        for n, (lat1, lat2, azi) in self._T2(self.test_equatorial_coincident,
               (+0.0, -0.0, 180.0),
               (-0.0, +0.0,   0.0)):
            r = g_WGS84.Inverse(lat1, 0.0, lat2, 0.0)
            self.test_eqv(n, r.azi1, azi)
            self.test_eqv(n, r.azi2, azi)

    def test_equatorial_NS(self):
        """Does the nearly antipodal equatorial solution go north or south?"""
        for n, (lat1, lat2, azi1, azi2) in self._T2(self.test_equatorial_NS,
               (+0.0, +0.0,  56., 124.),
               (-0.0, -0.0, 124.,  56.)):
            r = g_WGS84.Inverse(lat1, 0.0, lat2, 179.5)
            self.test_eqv(n, r.azi1, azi1, known=True, prec=2)
            self.test_eqv(n, r.azi2, azi2, known=True, prec=2)

    def test_atan2d(self):
        """Test special cases for atan2d"""
        eps, k = 7e-16, True
        for n, (a, b, x) in self._T2(self.test_atan2d,
               ( eps, -1.0,  180 - _atan2d(eps, 1.0)),  # 179.9..4 vs 7 Intel
               (+0.0, -0.0, +180.0),
               (-0.0, -0.0, -180.0),
               (+0.0, +0.0,   +0.0),
               (-0.0, +0.0,   -0.0),
               (+0.0, -1.0, +180.0),
               (-0.0, -1.0, -180.0),
               (+0.0, +1.0,   +0.0),
               (-0.0, +1.0,   -0.0),
               (-1.0, +0.0,  -90.0),
               (-1.0, -0.0,  -90.0),
               (+1.0, +0.0,  +90.0),
               (+1.0, -0.0,  +90.0),
               (+1.0, -INF, +180.0),
               (-1.0, -INF, -180.0),
               (+1.0,  INF,   +0.0),
               (-1.0,  INF,   -0.0),
               ( INF, +1.0,  +90.0),
               ( INF, -1.0,  +90.0),
               (-INF, +1.0,  -90.0),
               (-INF, -1.0,  -90.0),
               ( INF, -INF, +135.0),
               (-INF, -INF, -135.0),
               ( INF,  INF,  +45.0),
               (-INF,  INF,  -45.0),
               ( NAN, +1.0,    NAN),
               (+1.0,  NAN,    NAN)):
            self.test_eqv(n, _atan2d(a, b), x, known=k)
            k = False

    def test_sincosd(self):
        """Test special cases for sincosd"""
        for n, (d, x, y) in self._T2(self.test_sincosd,
               (-810.0, -1.0, +0.0),
               (-720.0, -0.0, +1.0),
               (-630.0, +1.0, +0.0),
               (-540.0, -0.0, -1.0),
               (-450.0, -1.0, +0.0),
               (-360.0, -0.0, +1.0),
               (-270.0, +1.0, +0.0),
               (-180.0, -0.0, -1.0),
               (- 90.0, -1.0, +0.0),
               (-  0.0, -0.0, +1.0),
               (+  0.0, +0.0, +1.0),
               (+ 90.0, +1.0, +0.0),
               (+180.0, +0.0, -1.0),
               (+270.0, -1.0, +0.0),
               (+360.0, +0.0, +1.0),
               (+450.0, +1.0, +0.0),
               (+540.0, +0.0, -1.0),
               (+630.0, -1.0, +0.0),
               (+720.0, +0.0, +1.0),
               (+810.0, +1.0, +0.0),
               (-  INF,  NAN,  NAN),
               (   INF,  NAN,  NAN),
               (   NAN,  NAN,  NAN)):
            s, c = _sincos2d(d)
            self.test_eqv(n, s, x, prec=16)  # error=(s - x)
            self.test_eqv(n, c, y)

        s1, c1 = _sincos2d(         9.0)
        s2, c2 = _sincos2d(        81.0)
        s3, c3 = _sincos2d(-123456789.0)
        self.test_eqv(n, s1,  c2, prec=16)
        self.test_eqv(n, c1,  s2)
        self.test_eqv(n, s1,  s3, prec=16)
        self.test_eqv(n, c1, -c3)

    def test_sum2(self):
        """Test special cases of sum"""
        for n, (a, b, x) in self._T2(self.test_sum2,
               (+9.0, -9.0, +0.0),
               (-9.0, +9.0, +0.0),
               (-0.0, +0.0, +0.0),
               (+0.0, -0.0, +0.0),
               (-0.0, -0.0, -0.0),
               (+0.0, +0.0, +0.0)):
            s, _ = _sum2(a, b)
            self.test_eqv(n, s, x)


if __name__ == '__main__':

    t = Tests(__file__, __version__)

    t.test_AngDiff()
    t.test_AngNormalize()
    t.test_AngRound()

    t.test_antipodal()
    t.test_antipodal_prolate()
    t.test_azimuth_0_180()

    t.test_equatorial_coincident()
    t.test_equatorial_NS()

    t.test_atan2d()
    t.test_sincosd()
    t.test_sum2()

    t.results()
    t.exit()
