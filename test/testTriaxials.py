
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '22.10.22'

from base import startswith, TestsBase

from pygeodesy import Ellipsoids, F_DMS, JacobiConformal, map2, PI_2, PI_4, Triaxial, triaxials, \
                      Vector3d as V3, sincos2d_
from math import radians

# <https://www.researchgate.net/profile/Maxim-Nyrtsov/publication/
# 299693481_Jacobi_Conformal_Projection_of_the_Triaxial_Ellipsoid_New_Projection_for_Mapping_of_Small_Celestial_Bodies/links/57ee6b2308ae8da3ce499cfc/
#           Jacobi-Conformal-Projection-of-the-Triaxial-Ellipsoid-New-Projection-for-Mapping-of-Small-Celestial-Bodies.pdf>


class Tests(TestsBase):

    def testJacobiConformal(self, module):
        self.subtitle(module, JacobiConformal.__name__)

        n = JacobiConformal.__name__
        # <https://GeographicLib.sourceforge.io/1.52/jacobi.html>
        j = JacobiConformal(6378137+35, 6378137-35, 6356752, name='Test')
        self.test(n, repr(j), "JacobiConformal(name='Test', a=6378172, b=6378102, c=6356752, e2ab=", known=startswith)

        n = j.xR.__name__
        x = j.xR(PI_2)
        self.test(n, x, '1.572092804', known=round(x, 7) == 1.5720928)

        n = j.yR.__name__
        y = j.yR(PI_2)
        self.test(n, y, '4.246581015',  known=round(y, 7) == 4.2465810)

        n = j.xyR2.__name__ + '.toDegrees'
        p = j.xyR2(PI_2, PI_2)
        t = p.toDegrees()
        self.test(n, t, '(90.074283, 243.31117)', known=map2(int, t) == (90, 243))
        t = p.toDegrees(form=F_DMS)
        self.test(n, t, "('90°04′27.42″', '243°18′40.21″')", known=True)
        n = j.area.name
        a = j.area
        self.test(n, int(a), '510065604942135', known=int(a * 1e-19) == 510)
        n = j.area_p.__name__
        p = j.area_p()
        self.test(n, int(p), '510065609807745', known=int(p * 1e-19) == 510)
        e = abs((a - p) / a)
        self.test('error', e, '9.54e-09', fmt='%.2e')
        n = j.volume.name
        p = j.volume
        self.test(n, p, '1.083207e+21', fmt='%.6e')

        n = JacobiConformal.__name__
        j = JacobiConformal(267.5, 147, 104.5, name='Itokawa25134')
        self.test(n, repr(j), "JacobiConformal(name='Itokawa25134', a=267.5, b=147, c=104.5, e2ab=", known=startswith, nl=1)

        n = j.xyR2.__name__
        p = j.xyR2(PI_4, 0)
        self.test(n, p, '(0.0, 0.61539)')

        n = p.toDegrees.__name__
        t = p.toDegrees()
        self.test(n, t, '(0.0, 35.259243)', known=True)
        t = p.toDegrees(form=F_DMS)
        self.test(n, t, "('00°00′00.0″', '35°15′33.27″')", known=True)

        n = JacobiConformal.xyQ2.name
        q = j.xyQ2
        self.test(n, q, '(3.13215, 1.42547)')

        n = q.toDegrees.__name__
        t = q.toDegrees()
        self.test(n, t, '(179.4589659, 81.673412)', known=True)
        t = q.toDegrees(form=F_DMS)
        self.test(n, t, "('179°27′32.28″', '81°40′24.28″')", known=True)

    def testTriaxial(self, module):
        self.subtitle(module, Triaxial.__name__)

        # <http://OJS.BBWPublisher.com/index.php/JWA/article/view/74>
        n = Triaxial.__name__
        t = Triaxial(6378388, 6378318, 6356911.9461, name='Test')  # a - b = 70
        self.test(n, repr(t), "Triaxial(name='Test', a=6378388, b=6378318, c=6356911.9461, e2ab=", known=startswith)

        n = t.forwardBetaOmega.__name__
        p = t.forwardBetaOmega(radians(30), radians(40), 1200)
        self.test(n, p, '(4234607.381429, 3551286.590486, 3176009.080037)', nl=1)
        v = V3(p)
        p = t.forwardBetaOmega(radians(30), radians(40))
        self.test(n, p, '(4233813.533025, 3550620.827453, 3175409.655093)')
        d = v.minus_(*p).length
        self.test('length', d, '1196.97367069', known=abs(d - 1200) < 5)

        n = t.forwardBetaOmega_.__name__
        p = t.forwardBetaOmega_(*sincos2d_(30, 40))  # h=0
        self.test(n, p, '(4233813.533025, 3550620.827453, 3175409.655093)')

        n = t.reverseLatLon.__name__
        q = t.reverseLatLon(*p)
        self.test(n, q, '(30.051881, 39.984967, 0.0)', nl=1)
        n = t.forwardLatLon.__name__
        q = t.forwardLatLon(*q)
        self.test(n, q, p)

        n = t.reverseBetaOmega.__name__
        p = t.reverseBetaOmega(4235882.4602, 3554249.4108, 3171030.2321)  # Ex-2 p 82 T 2

        self.test(n, p, '(0.520687, 0.698121, 12892.55755)', nl=1)  # (30, 40, 12000)
        p = t.reverseBetaOmega(4233721.2616, 3554717.2818, 3173743.2226)  # Ex-2 p 82 T 2
        self.test(n, p, '(0.521162, 0.698437, 12918.032538)')  # (30°01′38.2729″, 40°01′05.2057″, 1384.1361)
        n = p.toDegrees.__name__
        d = p.toDegrees()
        self.test(n, d, '(29.860398, 40.017494, 12918.032538)', known=True)
        d = p.toDegrees(form=F_DMS)
        self.test(n, d, "('29°51′37.43″', '40°01′02.98″', 12918.032538)", known=True)  # (30°01′38.2729″, 40°01′05.2057″)

        n = t.hartzell.__name__
        pov, los = V3(1e7, 1e7, 1e7), V3(-0.7274, -0.3637, -0.5819)
        self.test(n, t.hartzell(pov, los).toStr(prec=6), '(1125622.950901, 5562811.47545, 2900742.363389, 12200109.384877)', nl=1)
        self.test(n, t.hartzell(pov).toStr(prec=6),      '(3678403.581531, 3678403.581531, 3678403.581531, 10949326.181734)')
        E = Ellipsoids.WGS84  # earth as testFormy
        e = Triaxial(E.a, E.a, E.b, name=E.name)
        self.test(n, e.hartzell(pov, los).toStr(prec=6), '(1125440.234789, 5562720.117395, 2900596.195524, 12200360.575077)')
        self.test(n, e.hartzell(pov).toStr(prec=6),      '(3678289.79469, 3678289.79469, 3678289.79469, 10949523.266323)')

        t = Triaxial(3, 2, 1)  # Eberly
        n = t.height4.__name__
        p = t.height4(2, 4, 3)
        self.test(n, p.toStr(prec=6), '(1.206423, 1.61288, 0.433517, 3.593736)', nl=1)
        self.test(n, p.iteration, p.iteration)
        p = t.height4(-2, -4, -3)
        self.test(n, p.toStr(prec=6), '(-1.206423, -1.61288, -0.433517, 3.593736)')
        p = t.height4(0, 4, 3)
        self.test(n, p.toStr(prec=6), '(0.0, 1.746769, 0.487031, 3.375213)')
        p = t.height4(2, 0, 3)
        self.test(n, p.toStr(prec=6), '(1.563196, 0.0, 0.853517, 2.190477)')
        p = t.height4(2, 4, 0)
        self.test(n, p.toStr(prec=6), '(1.297504, 1.803267, 0.0, 2.306326)')


if __name__ == '__main__':

    t = Tests(__file__, __version__)
    t.testJacobiConformal(triaxials)
    t.testTriaxial(triaxials)
    t.results()
    t.exit()
