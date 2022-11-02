
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '22.11.02'

from base import random, startswith, TestsBase

from pygeodesy import EPS4, PI_2, PI_4, Ellipsoids, F_DMS, fstr, JacobiConformal, map1, map2, \
                      signBit, sincos2d_, Triaxial, Triaxial_, triaxials, Vector3d
from math import radians


def _r(s):
    r = random()
    if r < 0.05:
        r = 0
    elif r > 0.5:
        r = -r
    return r * s

# <https://www.researchgate.net/profile/Maxim-Nyrtsov/publication/
# 299693481_Jacobi_Conformal_Projection_of_the_Triaxial_Ellipsoid_New_Projection_for_Mapping_of_Small_Celestial_Bodies/links/57ee6b2308ae8da3ce499cfc/
#           Jacobi-Conformal-Projection-of-the-Triaxial-Ellipsoid-New-Projection-for-Mapping-of-Small-Celestial-Bodies.pdf>


class Tests(TestsBase):

    def testJacobiConformal(self, module):
        self.subtitle(module, JacobiConformal.__name__)

        n = JacobiConformal.__name__
        # <https://GeographicLib.sourceforge.io/1.52/jacobi.html>
        J = JacobiConformal(6378137+35, 6378137-35, 6356752, name='Test')
        self.test(n, repr(J), "JacobiConformal(name='Test', a=6378172, b=6378102, c=6356752, e2ab=", known=startswith)

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
        self.test(n, t, "('90°04′27.42″', '243°18′40.21″')", known=True)
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

        n = JacobiConformal.__name__
        J = JacobiConformal(267.5, 147, 104.5, name='Itokawa25134')
        self.test(n, repr(J), "JacobiConformal(name='Itokawa25134', a=267.5, b=147, c=104.5, e2ab=", known=startswith, nl=1)

        n = J.xyR2.__name__
        p = J.xyR2(PI_4, 0)
        self.test(n, p, '(0.0, 0.61539)')

        n = p.toDegrees.__name__
        t = p.toDegrees()
        self.test(n, t, '(0.0, 35.259243)', known=True)
        t = p.toDegrees(form=F_DMS)
        self.test(n, t, "('00°00′00.0″', '35°15′33.27″')", known=True)

        n = JacobiConformal.xyQ2.name
        q = J.xyQ2
        self.test(n, q, '(3.13215, 1.42547)')

        n = q.toDegrees.__name__
        t = q.toDegrees()
        self.test(n, t, '(179.4589659, 81.673412)', known=True)
        t = q.toDegrees(form=F_DMS)
        self.test(n, t, "('179°27′32.28″', '81°40′24.28″')", known=True)

    def testTriaxial(self, module):
        self.subtitle(module, Triaxial.__name__)

        E = Ellipsoids.WGS84  # earth as testFormy

        # <http://OJS.BBWPublisher.com/index.php/JWA/article/view/74>
        n = Triaxial.__name__
        T = Triaxial(6378388, 6378318, 6356911.9461, name='Test')  # a - b = 70
        self.test(n, repr(T), "Triaxial(name='Test', a=6378388, b=6378318, c=6356911.9461, e2ab=", known=startswith)
        n = Triaxial_.__name__
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
        v = Vector3d(q)
        p = T.forwardBetaOmega(radians(30), radians(40))
        self.test(n, p, '(4233813.533025, 3550620.827453, 3175409.655093)')
        d = v.minus_(*p).length
        self.test('length', d, '1196.973671', known=abs(d - 1200) < 5, prec=6)

        n = T.forwardCartesian.__name__
        t = T.forwardCartesian(*q, normal=True)
        self.test(n, t, '(4233813.533151, 3550620.827558, 3175409.654809, 1196.973671)', known=abs(t.h - 1200) < 5, nl=1)
        f = T.forwardCartesian(*q, normal=False)
        self.test(n, f, '(4239665.951888, 3553574.566129, 3164352.410834, 12911.309173)', known=abs(t.h - 1200) < 100)

        n = T.reverseCartesian.__name__
        t = T.reverseCartesian(*t)  # normal=True
        self.test(n, t, q, nl=1)  # q above
        f = T.reverseCartesian(*f, normal=False)  # off by 0.27%
        d = v.minus_(*f).length * 100 / v.length
        self.test(n, f, q, known=abs(d) < 0.5)

        n = T.forwardBetaOmega_.__name__
        p = T.forwardBetaOmega_(*sincos2d_(30, 40))  # h=0
        self.test(n, p, '(4233813.533025, 3550620.827453, 3175409.655093)', nl=1)

        n = T.reverseLatLon.__name__
        q = T.reverseLatLon(*p)
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

        n = T.hartzell4.__name__
        p = Vector3d(10.1e6, 10.2e6, 10.3e6)  # 10+ km
        v = Vector3d(-0.7274, -0.3637, -0.5819)
        t = T.hartzell4(p, v)
        self.test(n, t.toStr(prec=6), '(884268.349816, 5592134.174908, 2927668.068131, 12669388.912805)', nl=1)
        self.test(n, T.sideOf(t), 0)
        t = T.hartzell4(p)
        self.test(n, t.toStr(prec=6), '(3642143.609933, 3678204.437754, 3714265.265575, 11296443.179278)')
        self.test(n, T.sideOf(t), 0)

        e = Triaxial(E.a, E.a, E.b, name=E.name)
        self.test(n, e.hartzell4(p, v).toStr(prec=6), '(884080.396945, 5592040.198472, 2927517.711001, 12669647.302276)')
        self.test(n, e.hartzell4(p).toStr(prec=6),    '(3642031.283571, 3678090.99925, 3714150.714929, 11296639.666827)')

        t = U.hartzell4(p, v)
        self.test(n, t.toStr(prec=6), '(5585791.305438, 2917519.825775, 871582.610876, 12686828.717519)', nl=1)
        self.test(n, U.sideOf(t), 0)
        t = U.hartzell4(p)
        self.test(n, t.toStr(prec=6), '(3678286.263029, 3714347.893058, 3642224.632999, 11296301.449205)')
        self.test(n, U.sideOf(t, eps=EPS4), 0)

        T = Triaxial(3, 2, 1)  # Eberly
        n = T.height4.__name__
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

        t  = T.height4(6, 5, 4)
        t += t.iteration,
        f  = module._normalTo5
        s  = f(6, 5, 4, T)  # ordered
        self.test(f.__name__, fstr(t, prec=3, ints=True), fstr(s, prec=3, ints=True))
        f  = module._normalTo5u
        u  = f(5, 6, 4, 2, 3, 1)  # swap x and y
        s  = t[1], t[0], t[2], t[3], t[4]
        self.test(f.__name__, fstr(u, prec=3, ints=True), fstr(s, prec=3, ints=True))
        u  = f(4, 5, 6, 1, 2, 3)  # swap x and z
        s  = t[2], t[1], t[0], t[3], t[4]
        self.test(f.__name__, fstr(u, prec=3, ints=True), fstr(s, prec=3, ints=True))
        u  = f(5, 4, 6, 2, 1, 3)  # swap x and z, then x and y
        s  = t[1], t[2], t[0], t[3], t[4]
        self.test(f.__name__, fstr(u, prec=3, ints=True), fstr(s, prec=3, ints=True), nt=1)

#       n = U.height4.__name__
        a, b, c, _s = U.a * 2, U.b * 2, U.c * 2, signBit
        for _ in range(256):
            xyz = _r(a), _r(b), _r(c)
            p = U.height4(*xyz)
            s = '%s %s' % (str(p), p.iteration)
            self.test(n + str(xyz), s, s)
            sp = _s(p.x), _s(p.y), _s(p.z)
            sxyz = map2(_s, xyz)
            if sp != sxyz:
                self.test(n + str(xyz), sp, sxyz)

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


if __name__ == '__main__':

    t = Tests(__file__, __version__)
    t.testJacobiConformal(triaxials)
    t.testTriaxial(triaxials)
    t.results()
    t.exit()
