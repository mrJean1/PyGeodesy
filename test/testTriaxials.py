
# -*- coding: utf-8 -*-

# Test L{triaxials} module.
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # noqa: E702 ;

__all__ = ('Tests',)
__version__ = '25.12.12'

from bases import Geod3Solve, numpy, random, startswith, TestsBase

from pygeodesy import EPS4, F_DEG_, F_DMS, PI_2, PI_4, \
                      ConformalSphere, Conformal, Degrees, Ellipsoids, \
                      fstr, LLK, Los, map1, map2, signBit, sincos2d_, \
                      Triaxial, Triaxial_, Triaxials, Triaxial3s, \
                      triaxum5, Vector3d
from math import radians


def _r(s):
    r = random()
    if r < 0.05:
        r = 0
    elif r > 0.5:
        r = -r
    return r * s

# <https://www.ResearchGate.net/profile/Maxim-Nyrtsov/publication/299693481/links/57ee6b2308ae8da3ce499cfc/
#  Jacobi-Conformal-Projection-of-the-Triaxial-Ellipsoid-New-Projection-for-Mapping-of-Small-Celestial-Bodies.pdf>


class Tests(TestsBase):

    def testHartzell(self, module, LatLon):
        self.subtitle(module, 'Hartzell')

        T = Triaxial(6378388, 6378318, 6356911.9461, name='Test')  # a - b = 70
        n = T.hartzell4.__name__
        U = Triaxial_(T.c, T.a, T.b, name='Un')

        p = Vector3d(10.1e6, 10.2e6, 10.3e6)  # 10+ km
        v = Vector3d(-0.7274, -0.3637, -0.5819)
        t = T.hartzell4(p, v)
        self.test(n, t.toStr(prec=6), '(884268.349816, 5592134.174908, 2927668.068131, 12669388.912805)')
        self.test(n, T.sideOf(t), 0)
        t = T.hartzell4(p)
        self.test(n, t.toStr(prec=6), '(3642143.609933, 3678204.437754, 3714265.265575, 11296443.179278)')
        self.test(n, T.sideOf(t), 0)

        t = U.hartzell4(p, v)
        self.test(n, t.toStr(prec=6), '(888679.181482, 5594339.590741, 2931196.612187, 12663325.092381)', nl=1)
        self.test(n, U.sideOf(t), 0)
        t = U.hartzell4(p)
        self.test(n, t.toStr(prec=6), '(3642304.092727, 3678366.509487, 3714428.926247, 11296162.453809)')
        self.test(n, U.sideOf(t, eps=EPS4), 0)

        T = Ellipsoids.WGS84.toTriaxial()
        t = T.hartzell4(p, v)
        self.test(n, t.toStr(prec=6), '(884080.396945, 5592040.198472, 2927517.711001, 12669647.302276)', nl=1)
        self.test(n, T.sideOf(t), 0)
        t = T.hartzell4(p)
        self.test(n, t.toStr(prec=6), '(3642031.283571, 3678090.99925, 3714150.714929, 11296639.666827)')
        self.test(n, T.sideOf(t), 0)

        p = LatLon(30, 60, height=100e3)
        n = p.hartzell.__name__
        t = p.hartzell(Los(45, -45))
        self.test(n, t.toStr(prec=3), '30°38′27.119″N, 060°44′36.777″E, +142549.69m', nl=1)
        self.test(n, t.toStr(form=F_DEG_), '30.640866, 060.743549, +142549.69m')
        c = p.toCartesian()
        self.test(n, c.toStr(prec=3), '[2807429.59, 4862610.688, 3220373.735]')
        t = c.hartzell(Los(45, -45))
        self.test(n, t.toStr(prec=3), '[2684238.298, 4791786.806, 3231700.636]')
        t = (t - c).length
        self.test(n, t, '142549.6943849337', known=int(t)==142549)

        t = p.hartzell(los=True)
        self.test(n, t.toStr(prec=3), '30°00′00.0″N, 060°00′00.0″E, +100000.00m', nl=1)
        self.test(n, t.toStr(form=F_DEG_), '30.0, 060.0, +100000.00m')
        c = p.toCartesian()
        self.test(n, c.toStr(prec=3), '[2807429.59, 4862610.688, 3220373.735]')
        t = c.hartzell(los=True)
        self.test(n, t.toStr(prec=3), '[2764128.32, 4787610.688, 3170373.735]')
        t = (t - c).length
        self.test(n, t, '100000.0', known=abs(t - 100000.0) < 1e-6)

    def testConformal(self, module):
        self.subtitle(module, Conformal.__name__)

        n = Conformal.__name__
        # <https://GeographicLib.SourceForge.io/1.52/jacobi.html>
        J = Conformal(6378137+35, 6378137-35, 6356752, name='Test')
        self.test(n, repr(J), "%s(name='Test', a=6378172, b=6378102, c=6356752, e2ab=" % (n,), known=startswith)

        n = J.xR.__name__
        x = J.xR(PI_2)
        self.test(n, x, '1.572092804', known=round(x, 7) == 1.5720928)

        n = J.yR.__name__
        y = J.yR(PI_2)
        self.test(n, y, '4.246581015',  known=round(y, 7) == 4.2465810)

        n = J.xyR2.__name__ + '.toDegrees'
        p = J.xyR2(PI_2, PI_2)
        t = p.toDegrees()
        self.test(n, t, '(90.074283, 243.31117)', known=map2(int, t) == (90, 243))
        t = p.toDegrees(form=F_DMS)
        self.test(n, t, "('90°04′27.42″N', '243°18′40.21″E')", known=True)
        n = J.area.name
        a = J.area
        self.test(n, int(a), '510065604942135', known=int(a * 1e-19) == 510)
        n = J.area_p.__name__
        p = J.area_p()
        self.test(n, int(p), '510065609807745', known=int(p * 1e-19) == 510)
        e = abs((a - p) / a)
        self.test('error', e, '9.54e-09', fmt='%.2e')
        n = J.volume.name
        p = J.volume
        self.test(n, p, '1.083207e+21', fmt='%.6e')

        n = Conformal.__name__
        J = Conformal(267.5, 147, 104.5, name='Itokawa25134')
        self.test(n, repr(J), "%s(name='Itokawa25134', a=267.5, b=147, c=104.5, e2ab=" % (n,), known=startswith, nl=1)

        n = J.xyR2.__name__
        p = J.xyR2(PI_4, 0)
        self.test(n, p, '(0.0, 0.61539)')

        n = p.toDegrees.__name__
        t = p.toDegrees()
        self.test(n, t, '(0.0, 35.259243)', known=True)
        t = p.toDegrees(form=F_DMS)
        self.test(n, t, "('00°00′00.0″N', '035°15′33.27″E')", known=True)

        n = Conformal.xyQR2.name
        q = J.xyQR2
        self.test(n, q, '(3.13215, 1.42547)')

        n = q.toDegrees.__name__
        t = q.toDegrees()
        self.test(n, t, '(179.4589659, 81.673412)', known=True)
        t = q.toDegrees(form=F_DMS)
        self.test(n, t, "('179°27′32.28″N', '081°40′24.28″E')", known=True)

    def testConformalSphere(self, module):
        self.subtitle(module, ConformalSphere.__name__)

        n = ConformalSphere.__name__
        # <https://GeographicLib.SourceForge.io/1.52/jacobi.html>
        J = ConformalSphere(6378137+35, ab=1, bc=2, name='Test')
        self.test(n, repr(J), "%s(name='Test', a=6378172, ab=1, bc=2, e2ab=0, " % (n,), known=startswith)

        n = J.xR.__name__
        x = J.xR(PI_2)
        self.test(n, x, '1.73391688526', known=round(x, 7) == 1.7339169)

        n = J.yR.__name__
        y = J.yR(PI_2)
        self.test(n, y, '2.02895910275',  known=round(y, 7) == 2.0289591)

        n = J.xyR2.__name__ + '.toDegrees'
        p = J.xyR2(PI_2, PI_2)
        t = p.toDegrees()
        self.test(n, t, '(99.34612, 116.250793)', known=map2(int, t) == (99, 116))
        t = p.toDegrees(form=F_DMS)
        self.test(n, t, "('99°20′46.03″N', '116°15′02.86″E')", known=True)
        n = J.area.name
        a = J.area
        self.test(n, int(a), '511213503913540', known=int(a * 1e-19) == 511)
        n = J.area_p.__name__
        p = J.area_p()
        self.test(n, int(p), '511213503913540', known=int(p * 1e-19) == 511)
        e = abs((a - p) / a)
        self.test('error', e, '0.00e+00', fmt='%.2e')
        n = J.volume.name
        p = J.volume
        self.test(n, p, '1.086869e+21', fmt='%.6e')

        n = ConformalSphere.__name__
        J = ConformalSphere(267.5, 147, 104.5, name='Itokawa25134')
        self.test(n, repr(J), "%s(name='Itokawa25134', a=267.5, ab=147, bc=104.5, e2ab=0, " % (n,), known=startswith, nl=1)

        n = J.xyR2.__name__
        p = J.xyR2(PI_4, 0)
        self.test(n, p, '(0.0, 0.818354)')

        n = p.toDegrees.__name__
        t = p.toDegrees()
        self.test(n, t, '(0.0, 46.888217)', known=True)
        t = p.toDegrees(form=F_DMS)
        self.test(n, t, "('00°00′00.0″N', '046°53′17.58″E')", known=True)

        n = ConformalSphere.xyQR2.name
        q = J.xyQR2
        self.test(n, q, '(1.933157, 1.788429)')

        n = q.toDegrees.__name__
        t = q.toDegrees()
        self.test(n, t, '(179.4589659, 81.673412)', known=True)
        t = q.toDegrees(form=F_DMS)
        self.test(n, t, "('110°45′42.27″N', '102°28′10.04″E')", known=True)

    def testConformal3(self, module):
        n = module.Conformal3.__name__
        self.subtitle(module, n)

        # <https://GeographicLib.SourceForge.io/C++/doc/Cart3Convert.1.html>
        T = TX = module.Conformal3(Triaxials.WGS84_3)
        self.test(n, T, "name='WGS84_3', a=6378171.36, b=6378101.61, c=6356751.84, ", known=True)

        n = T.forwardBetOmg.__name__
        t = T.forwardBetOmg(33.3, 44.4, M=True, unit=Degrees)
        self.test(n, t.toStr(prec=5), "(-5077726.43188, 3922574.86203, 0, 1.19703, 'CONFORMAL')")  # Conformal3Proj
        n = T.reverseBetOmg.__name__
        t = T.reverseBetOmg(t.x, t.y, M=True)
        self.test(n, t.toDegrees(), "(Degrees(33.3), Degrees(44.4), None, 1.197032", known=startswith)  # Conformal3Proj

#       n = T.forwardSphere3.__name__
#       t, d, _ = T.forwardSphere3(33.3, 44.4, M=True, unit=Degrees)
#       self.test(n, t, "(-5077726.431881, 3922574.862033, 0, 1.1970322930911, 'CONFORMAL')", known=TBD)  # Conformal3Proj
#       n = T.reverseSphere.__name__
#       t = T.reverseSphere(t, dir3d=d, M=True)
#       self.test(n, t.toDegrees(), "(Degrees(33.3), Degrees(44.4), None, 1.197032", known=startswith)  # Conformal3Proj

        n = module.Conformal3.__name__
        T = module.Conformal3(Triaxials.WGS84_3r)  # rounded
        self.test(n, T, "name='WGS84_3r', a=6378172, b=6378102, c=6356752, ", known=startswith, nl=1)
        n = T.forwardBetOmg.__name__
        t = T.forwardBetOmg(Degrees(33.3), Degrees(44.4), M=True)
        self.test(n, t.toStr(prec=3), '(-5077732.396, 3922571.859, ', known=startswith)  # Conformal3Proj
        n = T.reverseBetOmg.__name__
        t = T.reverseBetOmg(t.x, t.y, M=True)
        self.test(n, t.toDegrees(), "(Degrees(33.3), Degrees(44.4), None, 1.197034, 'CONFORMAL')")

        n = T.forwardOther.__name__
        t = T.forwardOther(TX, 33.3, 44.4, M=True, unit=Degrees)
        self.test(n, t.toDegrees(0, prec=6), "(33.299887, 44.399927, 0.000263, 1.000046, 'CONFORMAL')")  # Conformal3Proj -tx
        n = T.reverseOther.__name__
        t = T.reverseOther(TX, t.bet, t.omg, M=True)
        self.test(n, t.toDegrees(0), "(Degrees(33.3), Degrees(44.4), Degrees(0.00026262), 0.999954, 'CONFORMAL')")  # Conformal3Proj -tx -r

        n = module.Conformal3.__name__
        c = 6378137 * (1 - 1 / (298257223563 / 1000000000))
        T = module.Conformal3(6378172, 6378102, c, name='WGS84+/-35')
        self.test(n, T, "name='WGS84+/-35', a=6378172, b=6378102, c=6356752.314245179, ", known=startswith, nl=1)
        n = T.forwardBetOmg.__name__
        t = T.forwardBetOmg(Degrees(33.3), Degrees(44.4), M=True)
        self.test(n, t.toStr(prec=3), '(-5077732.419, 3922572.019, 0', known=startswith)  # Conformal3Proj -t 6378172 6378102 6356752.314245179
        n = T.reverseBetOmg.__name__
        t = T.reverseBetOmg(t.x, t.y, M=True)
        self.test(n, t.toDegrees(), "(Degrees(33.3), Degrees(44.4), None, 1.197", known=startswith)

    def testGeod3Solve(self, Geodesic3Solve):
        g3S = Geodesic3Solve()
        n = g3S.__class__.__name__
        self.test(n, g3S, 'Geod3Solve=', known=startswith, nl=1)
        n = g3S.Direct.__name__
        t = g3S.Direct(40.57193, -54.38111, 3.20824, 15347602)
        self.test(n, t, '{a12: 137.869863, alp1: 3.20824, alp2: ', known=startswith)
        t = t.toGeod3Solve8Tuple()
        self.test(n, t, '(40.57193, ', known=startswith)

        n = g3S.Inverse.__name__
        t = g3S.Inverse(40.57193, -54.38111, 1.355292, 123.419711)
        self.test(n, t, '{a12: 137.869863, alp1: 3.20824, alp2: ', known=startswith)
        t = t.toGeod3Solve8Tuple()
        self.test(n, t, '(40.57193, ', known=startswith)

        gl3S = g3S.Line(40.57193, -54.38111, 3.20824)
        n = gl3S.__class__.__name__
        self.test(n, gl3S.toStr(), "alp1=", known=startswith)
        n = gl3S.Position.__name__
        t = gl3S.Position(15347602)
        self.test(n, t, '{a12: 137.869863, alp1: 3.20824, alp2: ', known=startswith)
        t = t.toGeod3Solve8Tuple()
        self.test(n, t, '(40.57193, ', known=startswith)

        gl3S = g3S.InverseLine(40.57193, -54.38111, 1.355292, 123.419711)
        n = gl3S.__class__.__name__
        self.test(n, gl3S.toStr(), "alp1=", known=startswith)
        n = gl3S.Position.__name__
        t = gl3S.Position(15347602)
        self.test(n, t, '{a12: 137.869863, alp1: 3.20824, alp2: ', known=startswith)
        t = t.toGeod3Solve8Tuple()
        self.test(n, t, '(40.57193, ', known=startswith)

    def testTriaxial3(self, module):
        n = module.Triaxial3.__name__
        self.subtitle(module, n)

#       E = Ellipsoids.WGS84  # earth as testFormy

        # <https://GeographicLib.SourceForge.io/C++/doc/Cart3Convert.1.html>
        T = module.Triaxial3(3, 2, 1)
        self.test(n, repr(T), "Triaxial3(name='', a=3, b=2, c=1, k2=0.375, kp2=0.625, ", known=startswith)

        n = T.reverseLatLon.__name__
        t = T.reverseLatLon(1, 2, 3)
        self.test(n, t.toDegrees(0), "(Degrees(58.69140449), Degrees(75.11263103), None, 2.586065, 'ELLIPSOIDAL')")
        n = T.forwardLatLon.__name__
        t = T.forwardLatLon(t.bet, t.omg, h=t.h)
        self.test(n, t, "(1.0, 2.0, 3.0, 2.586065, 'ELLIPSOIDAL')")

        # <https://GeographicLib.SourceForge.io/C++/doc/Cart3Convert.1.html>
#       T = module.Triaxial3(3, 2, 1)
#       self.test(n, repr(T), "Triaxial3(name='', a=3, b=2, c=1, k2=0.375, kp2=0.625, ", known=startswith, nl=1)
        for llk, x in ((LLK.ELLIPSOIDAL,   "(58.691404, 75.112631, None, 2.586065, "),  # Cart3Convert
                       (LLK.GEOCENTRIC,    "(29.12663, 56.916602, None, 2.391078, "),  # Cart3Convert
                       (LLK.GEOCENTRIC_X,  "(28.478775, 123.624552, None, 2.391078, "),  # Cart3Convert
                       (LLK.GEODETIC,      "(68.626017, 73.851827, None, 2.391078, "),  # Cart3Convert
# not WGS_3.           (LLK.GEODETIC_LON0, "(68.626017, 58.921827, None, 2.391078, "),  # Cart3Convert
                       (LLK.GEODETIC_X,    "(5.817652, 159.397221, None, 2.391078, "),  # Cart3Convert
                       (LLK.PARAMETRIC,    "(50.658091, 66.523762, None, 2.391078, "),  # Cart3Convert
                       (LLK.PARAMETRIC_X,  "(14.628136, 143.06191, None, 2.391078, ")):  # Cart3Convert
            n = str(llk)
            t = T.reverseLatLon(1, 2, 3, llk=llk)
            self.test(n, t.toDegrees(0, fmt=F_DEG_), x, known=startswith, nl=1)
            t = T.forwardLatLon(t.bet, t.omg, height=t.h) if llk in LLK._NOIDAL else \
                T.forwardLatLon(t.phi, t.lam, height=t.h, llk=llk)
            self.test(n, t, "(1.0, 2.0, 3.0, 2.", known=startswith)

        n = module.Triaxial3.__name__
        # <http://OJS.BBWPublisher.com/index.php/JWA/article/view/74>  # Bektas
        T = module.Triaxial3(6378388, 6378318, 6356911.9461, name='Bektas')  # a - b = 70
        self.test(n, repr(T), "Triaxial3(name='Bektas', a=6378388, b=6378318, c=6356911.9461, k2=", known=startswith, nl=1)

        n = T.forwardBetaOmega.__name__
        q = T.forwardBetaOmega(30, 40, 1200, unit=Degrees)
        self.test(n, q, "(4234607.381429, 3551286.590486, 3176009.080037, 1200.0, 'ELLIPSOIDAL')")
        v = Vector3d(q.xyz3)
        p = T.forwardBetaOmega(radians(30), radians(40))
        self.test(n, p, "(4233813.533025, 3550620.827453, 3175409.655093, 0, 'ELLIPSOIDAL')")  # C++
        d = v.minus_(*p[:3]).length
        self.test('length', d, '1196.973671', known=abs(d - 1200) < 5, prec=6)
        n = T.reverseBetaOmega.__name__
        t = T.reverseBetaOmega(4234607.381429, 3551286.590486, 3176009.080037)
        self.test(n, t.toDegrees(0), "(Degrees(30.0), Degrees(40.0), None, 1200.0, 'ELLIPSOIDAL')")
        t = T.reverseBetaOmega(4233813.533025, 3550620.827453, 3175409.655093)
        self.test(n, t.toDegrees(0), "(Degrees(30.0), Degrees(40.0), None, 0.0, 'ELLIPSOIDAL')")  # C++

        n = T.forwardCartesian.__name__
        t = T.forwardCartesian(q, normal=True)
        self.test(n, t, "(4233813.533151, 3550620.827558, 3175409.654809, 1196.973671, 'ELLIPSOIDAL')", known=abs(t.h - 1200) < 5, nl=1)
        f = T.forwardCartesian(q, normal=False)
        self.test(n, f, "(4239665.951888, 3553574.566129, 3164352.410834, 12911.309173, 'ELLIPSOIDAL')", known=abs(t.h - 1200) < 100)

        n = T.reverseCartesian.__name__
        t = T.reverseCartesian(t, height=t.h)  # normal=True
        self.test(n, t, "(4234607.381429, 3551286.590486, 3176009.080037, ", known=startswith, nl=1)  # == q above
        f = T.reverseCartesian(f, height=f.h, normal=False)  # off by 0.27%
        d = v.minus_(f[:3]).length * 100 / v.length
        self.test(n, f, q, known=abs(d) < 0.5)

        n = T.forwardBetaOmega_.__name__
        p = T.forwardBetaOmega_(*sincos2d_(30, 40))  # h=0
        self.test(n, p, "(4233813.533025, 3550620.827453, 3175409.655093, 0, 'ELLIPSOIDAL')", nl=1)
        n = T.reverseLatLon.__name__
        q = T.reverseLatLon(p)
        self.test(n, q.toDegrees(0), "(Degrees(30.0), Degrees(40.0), None, ", known=startswith)
        n = T.forwardLatLon.__name__
        q = T.forwardLatLon(q.bet, q.omg)
        self.test(n, q, p)

        n = T.reverseBetaOmega.__name__
        p = T.reverseBetaOmega(4235882.4602, 3554249.4108, 3171030.2321)  # Ex-2 p 82 T 2
        self.test(n, p.toDegrees(0), "(Degrees(29.94812666), Degrees(40.01497072), None, 1203.037176, 'ELLIPSOIDAL')", nl=1)  # (30, 40, 12000)
        t = p.toDegrees(0, fmt=F_DMS)  # sep=_COMMASPACE_
        self.test(n, t, "(29°56′53.26″, 40°00′53.89″, None, 1203.0", known=startswith)
        p = T.reverseBetaOmega(4233721.2616, 3554717.2818, 3173743.2226)  # Ex-2 p 82 T 2
        self.test(n, p.toDegrees(0), "(Degrees(29.97539672), Degrees(40.03311872), None, 1387.637345, 'ELLIPSOIDAL')")
#       self.test(n, p.toRadians(0), "(Radians(0.52316937),  Radians(0.69870973),  None, 1387.637345, 'ELLIPSOIDAL')")
        t = p.toDegrees(0, fmt=F_DMS)  # sep=_COMMASPACE_
        self.test(n, t, "(29°58′31.43″, 40°01′59.23″, None, 1387.6", known=startswith)  # (30°01′38.2729″, 40°01′05.2057″)

        n = T.height4.__name__
        t = T.height4(3909863.9271, 3909778.123, 3170932.5016)  # Bektas
        self.test(n, t, '(3909251.554667, 3909165.750567, 3170432.501602, 999.999996)', nl=1)
        ct = T.forwardCartesian(t, normal=True)
        self.test(n, ct, '(3909251.554667, 3909165.750567, 3170432.501602, 0, None)')
        ct = T.forwardCartesian(t, normal=False)
        self.test(n, ct, '(3909251.554667, 3909165.750567, 3170432.501602, 0, None)')  # XXX

        n = 'JFK-SIN'
        # <https://GeographicLib.SourceForge.io/C++/doc/Geod3Solve.1.html>
        # echo 40:38:23N 073:46:44W-19.43W X 01:21:33N 103:59:22E-19.43W | \  # note 19.43W
        # tr X '\n' | Cart3Convert -G | Cart3Convert -E -r | tr '\n' ' ' | Geod3Solve -i -: -p 0
        # 003:12:29.7 177:28:59.5 15347602 == 3.20824 177.48319 15347602
        T = Triaxial3s.WGS84_35  # module.Triaxial3(Triaxials.WGS84_35)  # a - b = 70
        self.test(n, repr(T), "Triaxial3(name='WGS84_35', a=6378172, b=6378102, c=6356752.314245179, k2=0.9967", known=startswith, nl=1)
        d = Degrees('73 46 44W') - Degrees('19.43W')  # XXX T.lon0?
        self.test(n, d, '-54.34', prec=5, known=startswith)
        f = T.forwardLatLon(Degrees('40 38 23N'), Degrees(d), llk=LLK.GEODETIC)
        self.test(n, f, "(2824949.36608, -3938333.736799, 4132149.896611, 0, 'GEODETIC')")  # 2824949.425 -3938333.819 4132149.574
        r = T.reverseLatLon(f)
        self.test(n, r.toDegrees(0, fmt=F_DMS), "(40°38′23.0″, 54°20′56.0″, None, 0", known=startswith)  # 0.0, 'ELLIPSOIDAL')")
        t = T.reverseLatLon(f, llk=LLK.ELLIPSOIDAL)
        self.test(n, t.toDegrees(0), "(Degrees(40.57193395), Degrees(-54.38110954), None, 0", known=startswith)  # 0.0, 'ELLIPSOIDAL')")  # 40.57193215 -54.38110906

        d = Degrees('103 59 22E') - Degrees('19.43W')  # XXX T.lon0?
        self.test(n, d, '123.4', prec=2, known=startswith, nl=1)
        f = T.forwardLatLon(Degrees('1 21 33N'), Degrees(d), llk=LLK.GEODETIC)
        self.test(n, f, "(-3511912.82574, 5322047.492059, 150275.382099, 0, 'GEODETIC')")  # -3511912.826 5322047.492 150275.367
        r = T.reverseLatLon(f)  # llk=LLK.GEODETIC)
        self.test(n, r.toDegrees(0, fmt=F_DMS), "(1°21′33.0″, 123°25′10.0″, None, ", known=startswith)
        t = T.reverseLatLon(f, llk=LLK.ELLIPSOIDAL)
        self.test(n, t.toDegrees(0, prec=6), "(1.355287, 123.419709, None, 0", known=startswith)  # 1.35528738 123.41970939

        # echo 40:38:23N 073:46:44W-14.93W X 01:21:33N 103:59:22E-14.93W | \  # note 14.93W
        # tr X '\n' | Cart3Convert -G | Cart3Convert -E -r | tr '\n' ' ' | Geod3Solve -i -: -p 0
        # 003:12:53.4 177:29:00.3 15347592 == 3.21482 177.48341 15347592
        f = T.forwardLatLon(Degrees('40 38 23N'), Degrees('73 46 44W'), llk=LLK.GEODETIC_LON0)
        self.test(n, f, "(2507237.249613, -4147833.171672, 4132151.785141, 0, 'GEODETIC_LON0')")  # 2507237.302 -4147833.258 4132151.463
        r = T.reverseLatLon(f, llk=LLK.ELLIPSOIDAL)
        self.test(n, r.toDegrees(0), "(Degrees(40.56616585), Degrees(-58.87899203), ", known=startswith)  # 40.56616414 -58.87899158

        f = T.forwardLatLon(Degrees('1 21 33N'), Degrees('103 59 22E'), llk=LLK.GEODETIC_LON0)
        self.test(n, f, "(-3083516.921703, 5581181.10656, 150275.496647, 0, 'GEODETIC_LON0')")  # -3083516.922 5581181.107 150275.482
        r = T.reverseLatLon(f, llk=LLK.ELLIPSOIDAL)
        self.test(n, r.toDegrees(0), "(Degrees(1.35513418), Degrees(118.9196884), ", known=startswith)  # 1.35513411 118.91968840

        for n, llk in sorted(LLK.items()):
            if llk is not LLK.CONFORMAL:
                ct, d3 = T.random2(llk, True)
                self.test(n, ct, ct, nl=1)
                self.test(n, d3, d3)
                if llk is LLK.ELLIPSOIDAL:
                    r    = T.reverseBetOmgAlp(ct, dir3d=d3)
                    f, d = T.forwardBetOmgAlp2(r.bet, r.omg, r.alp)
                else:
                    r    = T.reversePhiLamZet(ct, dir3d=d3, llk=llk)
                    f, d = T.forwardPhiLamZet2(r.phi, r.lam, r.zet, llk=llk)
                self.test(n, r, r)
                self.test(n, f, ct)
                self.test(n, d, d3, known=True)
                _ = T.normed2(ct, d3)  # PYCHOK coverage

        self.test(T.Lon0.name, T.Lon0, -14.93, nl=1)
        n = T.forwardLatLon.__name__
        t = T.forwardLatLon(0, T.Lon0)  # bi- to triaxial lon
        self.test(n, t, "(6162853.284268, -1643246.23441, 0.0, 0, 'ELLIPSOIDAL')")
        n = T.reverseLatLon.__name__
        t = T.reverseLatLon(t)
        self.test(n, t.toDegrees(0), "(Degrees(0.0), Degrees(-14.93), None, ", known=startswith)
        t = t.omg - T.Lon0
        self.test(n, t.degrees0, 0, known=abs(t.degrees0) < 1e-13)  # tri- to biaxial lon

    def testTriaxial5(self, module):
        n = Triaxial.__name__
        self.subtitle(module, n)

        E = Ellipsoids.WGS84  # earth as testFormy

        # <http://OJS.BBWPublisher.com/index.php/JWA/article/view/74>  # Bektas
        T = Triaxial(6378388, 6378318, 6356911.9461, name='Test')  # a - b = 70
        self.test(n, repr(T), "Triaxial(name='Test', a=6378388, b=6378318, c=6356911.9461, e2ab=", known=startswith)
        U = Triaxial_(T.c, T.a, T.b, name='Un')
        self.test(n, repr(U), "Triaxial_(name='Un', a=6356911.9461, b=6378388, c=6378318, e2ab=", known=startswith)

#       e = EPS4
#       for a in range(91):
#           for b in range(91):
#               sa, ca, sb, cb = sincos2d_(a, b)
#               sU = U.sideOf(*U._radialTo3(sa, ca, sb, cb))
#               sT = T.sideOf(*T._radialTo3(sa, ca, sb, cb))
#               s = abs(sU - sT)
#               if s > e:
#                   self.test('sU-sT%r' % ((a, b),), s, s)
#                   e = s

        n = T.forwardBetaOmega.__name__
        q = T.forwardBetaOmega(radians(30), radians(40), 1200)
        self.test(n, q, '(4234607.381429, 3551286.590486, 3176009.080037)', nl=1)
        v = Vector3d(q.xyz3)
        p = T.forwardBetaOmega(radians(30), radians(40))
        self.test(n, p, '(4233813.533025, 3550620.827453, 3175409.655093)')
        d = v.minus_(*p[:3]).length
        self.test('length', d, '1196.973671', known=abs(d - 1200) < 5, prec=6)

        n = T.forwardCartesian.__name__
        t = T.forwardCartesian(q, normal=True)
        self.test(n, t, '(4233813.533151, 3550620.827558, 3175409.654809, 1196.973671)', known=abs(t.h - 1200) < 5, nl=1)
        f = T.forwardCartesian(q, normal=False)
        self.test(n, f, '(4239665.951888, 3553574.566129, 3164352.410834, 12911.309173)', known=abs(t.h - 1200) < 100)

        n = T.reverseCartesian.__name__
        t = T.reverseCartesian(*t)  # normal=True
        self.test(n, t, q, nl=1)  # q above
        f = T.reverseCartesian(*f, normal=False)  # off by 0.27%
        d = v.minus_(f[:3]).length * 100 / v.length
        self.test(n, f, q, known=abs(d) < 0.5)

        n = T.forwardBetaOmega_.__name__
        p = T.forwardBetaOmega_(*sincos2d_(30, 40))  # h=0
        self.test(n, p, '(4233813.533025, 3550620.827453, 3175409.655093)', nl=1)

        n = T.reverseLatLon.__name__
        q = T.reverseLatLon(p)
        self.test(n, q, '(30.051881, 39.984967, 0.0)', nl=1)
        n = T.forwardLatLon.__name__
        q = T.forwardLatLon(*q)
        self.test(n, q, p)

        n = T.reverseBetaOmega.__name__
        p = T.reverseBetaOmega(4235882.4602, 3554249.4108, 3171030.2321)  # Ex-2 p 82 T 2

        self.test(n, p, '(0.520687, 0.698121, 12892.55755)', nl=1)  # (30, 40, 12000)
        p = T.reverseBetaOmega(4233721.2616, 3554717.2818, 3173743.2226)  # Ex-2 p 82 T 2
        self.test(n, p, '(0.521162, 0.698437, 12918.032538)')  # (30°01′38.2729″, 40°01′05.2057″, 1384.1361)
        n = p.toDegrees.__name__
        d = p.toDegrees()
        self.test(n, d, '(29.860398, 40.017494, 12918.032538)', known=True)
        d = p.toDegrees(form=F_DMS)
        self.test(n, d, "('29°51′37.43″', '40°01′02.98″', 12918.032538)", known=True)  # (30°01′38.2729″, 40°01′05.2057″)

        n = T.height4.__name__
        t = T.height4(3909863.9271, 3909778.123, 3170932.5016)  # Bektas
        self.test(n, t, '(3909251.554667, 3909165.750567, 3170432.501602, 999.999996)', nl=1)

        T = Triaxial(3, 2, 1)  # Eberly
        p = T.height4(2, 4, 3)
        self.test(n, p.toStr(prec=6), '(1.206423, 1.61288, 0.433517, 3.593736)', nl=1)
        self.test(n, p.iteration, p.iteration)
        self.test(n, T.sideOf(p), 0)
        p = T.height4(-2, -4, -3)
        self.test(n, p.toStr(prec=6), '(-1.206423, -1.61288, -0.433517, 3.593736)')
        p = T.height4(0, 4, 3)
        self.test(n, p.toStr(prec=6), '(0.0, 1.746769, 0.487031, 3.375213)')
        p = T.height4(2, 0, 3)
        self.test(n, p.toStr(prec=6), '(1.563196, 0.0, 0.853517, 2.190477)')
        p = T.height4(2, 4, 0)
        self.test(n, p.toStr(prec=6), '(1.297504, 1.803267, 0.0, 2.306326)', nt=1)

        b = signBit.__name__
        for x in (-2, 0, 2):
            for y in (-4, 0, 4):
                for z in (-3, 0, 3):
                    p = T.height4(x, y, z)
                    s = '%s %s' % (str(p), p.iteration)
                    xyz = str((x, y, z))
                    self.test(n + xyz, s, s)  # decoration
                    self.test(b + xyz, map1(signBit, p.x, p.y, p.z),
                                       map1(signBit,   x,   y,   z), nt=1)

        a, b, c = 6, 5, 4
        t  = T.height4(a, b, c)
        t += t.iteration,
        x, y, z, d, i = t
        f  = module._plumbTo5
        n  = f.__name__
        s  = f(a, b, c, T)  # ordered
        self.test(f.__name__, fstr(t, prec=3, ints=True), fstr(s, prec=3, ints=True))
        u  = f(a, c, b, Triaxial_(T.a, T.c, T.b))
        s  = x, z, y, d, i
        self.test(n, fstr(u, prec=3, ints=True), fstr(s, prec=3, ints=True))
        u  = f(b, a, c, Triaxial_(T.b, T.a, T.c))
        s  = y, x, z, d, i
        self.test(n, fstr(u, prec=3, ints=True), fstr(s, prec=3, ints=True))
        u  = f(b, c, a, Triaxial_(T.b, T.c, T.a))
        s  = y, z, x, d, i
        self.test(n, fstr(u, prec=3, ints=True), fstr(s, prec=3, ints=True))
        u  = f(c, a, b, Triaxial_(T.c, T.a, T.b))
        s  = z, x, y, d, i
        self.test(n, fstr(u, prec=3, ints=True), fstr(s, prec=3, ints=True))
        u  = f(c, b, a, Triaxial_(T.c, T.b, T.a))
        s  = z, y, x, d, i
        self.test(n, fstr(u, prec=3, ints=True), fstr(s, prec=3, ints=True), nt=1)

#       n = U.height4.__name__
        x, y, z, _s = U.a * 2, U.b * 2, U.c * 2, signBit
        for _ in range(256):
            xyz = _r(x), _r(y), _r(z)
            s = n + str(xyz)
            p = U.height4(*xyz)
            t = '%s %s' % (str(p), p.iteration)
            self.test(s, t, t)
            sp   = map2(_s, p[:3])
            sxyz = map2(_s, xyz)
            if sp != sxyz:
                self.test(s, sp, sxyz)

#       for xyz in ((0, 0, 0), (0, 0, 1), (0, 1, 0),
#                   (0, 1, 1), (1, 0, 0), (1, 0, 1),
#                   (1, 1, 0), (1, 1, 1)):
#           p = T.height4(*xyz)
#           s = '%s %s' % (str(p), p.iteration)
#           self.test(n + str(xyz), s, s)

        U = Triaxial_(3, 2, 1)  # coverage
        p = U.height4(2, 4, 3, normal=False)
        self.test(n, p.toStr(prec=6), '(0.545455, 1.090909, 0.818182, 3.916483)', nl=1)  # ???
        self.test(n, p.iteration, p.iteration)
        self.test(n, U.sideOf(p), 0)

        U = Triaxial_(2, 3, 1)  # coverage
        p = U.height4(4, 2, 3, normal=False)
        self.test(n, p.toStr(prec=6), '(1.090909, 0.545455, 0.818182, 3.916483)', nl=1)  # ???
        self.test(n, p.iteration, p.iteration)
        self.test(n, U.sideOf(p), 0)

        U = Triaxial_(2, 2, 2)  # coverage
        p = U.height4(2, 3, 4, normal=False)
        self.test(n, p.toStr(prec=6), '(0.742781, 1.114172, 1.485563, 3.385165)', nl=1)  # ???
        self.test(n, p.iteration, p.iteration)
        self.test(n, U.sideOf(p), 0)

        n = E.toTriaxial.__name__
        T = E.toTriaxial()
        t = str(T)
        self.test(n, t, t, nl=1)
        n = T.toEllipsoid.__name__
        e = T.toEllipsoid(name='_')
        t = str(e)
        self.test(n, t, t)
        T = Triaxial(3, 2, 2)
        e = T.toEllipsoid()
        t = str(e)
        self.test(n, t, t)

        self.test('Triaxials', len(Triaxials.items(all=True)), 14, nl=1)
        for t in Triaxials.values(all=True):
            self.test(t.name, t, getattr(Triaxials, t.name))

    def testTriaxum5(self):
        # <https://Jekel.me/2020/Least-Squares-Ellipsoid-Fit/>
        ps, v = [], Vector3d(1.2, 0.2, 0.9)
        for d in range(0, 180, 9):
            su, cu, sv, cv = sincos2d_(d, d * 0.5)
            ps.append(v.times_(cu * sv, su * sv, cv))
        n = triaxum5.__name__
        t = triaxum5(ps, useZ=True)
        self.test(n, t, '(1.2, 0.2, 0.9, 3, ', known=startswith, nl=1)
        t = triaxum5(ps, useZ=False)
        self.test(n, t, '(1.244625, 0.145582, 0.0, 2, ', known=startswith)


if __name__ == '__main__':

    from pygeodesy.ellipsoidalVincenty import LatLon
    import pygeodesy.triaxials.triaxial5 as triaxial5
    import pygeodesy.triaxials.conformal3 as conformal3
    import pygeodesy.triaxials.triaxial3 as triaxial3

    t = Tests(__file__, __version__)
    t.testHartzell(triaxial5, LatLon)
    t.testConformal(triaxial5)
    t.testConformalSphere(triaxial5)
    t.testTriaxial5(triaxial5)
    t.testTriaxial3(triaxial3)
    t.testConformal3(conformal3)
    if numpy:
        t.testTriaxum5()
    else:
        t.skip(triaxum5.__name__, 2)
    if Geod3Solve:
        from pygeodesy.geod3solve import Geodesic3Solve
        t.testGeod3Solve(Geodesic3Solve)
    t.results()
    t.exit()
