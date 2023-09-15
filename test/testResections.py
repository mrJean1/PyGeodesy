
# -*- coding: utf-8 -*-

# Test L{resections} module.

__all__ = ('Tests',)
__version__ = '23.09.11'

from bases import endswith, isWindows, TestsBase

from pygeodesy import EPS0, PI, PI_4, cassini, collins5, fstr, pierlot, pierlotx, \
                      ResectionError, snellius3, tienstra7, triAngle, triAngle5, \
                      triArea, triSide, triSide2, triSide4, Vector3d, wildberger3, \
                      collins, tienstra, triAngle4  # DEPRECATED
from math import degrees


class Tests(TestsBase):

    def testResections(self, C_, V3d):

        # Tienstra example <http://MesaMike.org/geocache/GC1B0Q9/resection-methods.pdf>
        A, B, C = V3d(1000, 5300), V3d(2200, 6300), V3d(3100, 5000)  # DDD.MMSS 115.5220 109.3045

        p = cassini(A, C, B, 109.5125, 115.0889)  # note B center, alpha and beta definition
        self.test(cassini.__name__, p.toStr(prec=4), '(2128.3903, 5578.1443, 0)')
        p = C_(A).cassini(C, B, 109.5125, 115.0889)
        self.test(cassini.__name__, p.toRepr(prec=4), 'Cartesian_(2128.3903, 5578.1443, 0)')

        t = collins5(A, C, B, 109.5125, 115.0889)  # note B center, alpha and beta definition
        self.test(collins5.__name__, t.pointP.toStr(prec=4), '(2128.3903, 5578.1443, 0)', nl=1)
        self.test(collins5.__name__, t.pointH.toStr(prec=4), '(1830.5948, 2576.2429, 0)')
        r = fstr(t[2:], prec=4)
        self.test(collins5.__name__, r, '1581.1388, 1562.0499, 2121.3203')
        t = C_(A).collins5(C, B, 109.5125, 115.0889)
        self.test(collins5.__name__, t.pointP.toRepr(prec=4), 'Cartesian_(2128.3903, 5578.1443, 0)')
        self.test(collins5.__name__, t.pointH.toRepr(prec=4), 'Cartesian_(1830.5948, 2576.2429, 0)')
        self.test(collins5.__name__, fstr(t[2:], prec=4), r)
        self.test(collins5.__name__, fstr(t[2:], prec=4), r)
        self.test(collins.__name__, C_(A).collins(C, B, 109.5125, 115.0889), t, nl=1)  # DEPRECATED

        p = pierlot(C, B, A, 115.0889, 109.5125)  # note CCW order, alpha12 and alpha23 definition
        self.test(pierlot.__name__, p.toStr(prec=4), '(2128.3903, 5578.1443, 0)', nl=1)
        p = C_(C).pierlot(B, A, 115.0889, 109.5125)
        self.test(pierlot.__name__, p.toRepr(prec=4), 'Cartesian_(2128.3903, 5578.1443, 0)')
        p = C_(C).pierlot(B, A, 115.0889, 109.5125, useZ=True)  # _zidw coverage
        self.test(pierlot.__name__, p.toRepr(prec=4), 'Cartesian_(2128.3903, 5578.1443, 0.0)')

        p = pierlotx(C, B, A, -115.0889, 0, 109.5125)  # note CCW order, alpha1, alpha2, alpha3 definition
        self.test(pierlotx.__name__, p.toStr(prec=4), '(2128.3903, 5578.1443, 0)', nl=1)
        p = C_(C).pierlotx(B, A, -115.0889, 0, 109.5125)
        self.test(pierlotx.__name__, p.toRepr(prec=4), 'Cartesian_(2128.3903, 5578.1443, 0)')
        p = C_(C).pierlotx(B, A, -115.0889, 0, 109.5125, useZ=True)  # _zidw coverage
        self.test(pierlotx.__name__, p.toRepr(prec=4), 'Cartesian_(2128.3903, 5578.1443, 0.0)')

        t = tienstra7(A, B, C, 115.0889, None, 109.5125)  # note alpha, beta and gamma definition
        self.test(tienstra7.__name__, t.pointP.toStr(prec=4), '(2128.3903, 5578.1443, 0)', nl=1)
        r = fstr(t[1:], prec=4)
        self.test(tienstra7.__name__, r, '47.9357, 84.8896, 47.1747, 1581.1388, 2121.3203, 1562.0499')
        t = C_(A).tienstra7(B, C, 115.0889, None, 109.5125)
        self.test(tienstra7.__name__, t.pointP.toRepr(prec=4), 'Cartesian_(2128.3903, 5578.1443, 0)')
        self.test(tienstra7.__name__, fstr(t[1:], prec=4), r)
        self.test(tienstra.__name__, C_(A).tienstra(B, C, 115.0889, None, 109.5125), t, nl=1)  # DEPRECATED

        p = cassini(A, C, B, 109.3, 115.1)  # note B center, alpha and beta definition
        self.test(cassini.__name__, p.toStr(prec=4), '(2129.3018, 5575.8016, 0)', nl=1)

        t = collins5(A, C, B, 109.3, 115.1)  # note B center, alpha and beta definition
        self.test(collins5.__name__, t.pointP.toStr(prec=4), '(2129.3018, 5575.8016, 0)', nl=1)
        self.test(collins5.__name__, t.pointH.toStr(prec=4), '(1835.1911, 2563.0708, 0)')
        self.test(collins5.__name__, fstr(t[2:], prec=4), '1581.1388, 1562.0499, 2121.3203')
        self.test(collins.__name__, collins(A, C, B, 109.3, 115.1), t, nl=1)  # DEPRECATED

        p = pierlot(C, B, A, 115.1, 109.3)  # note CCW order, alpha12 and alpha23 definition
        self.test(pierlot.__name__, p.toStr(prec=4), '(2129.3018, 5575.8016, 0)', nl=1)
        try:
            p = pierlot(C, B, A, 115.1, 109.3, eps=0).toStr(prec=4)
        except ResectionError as x:
            p = str(x)
        self.test(pierlot.__name__, p, 'eps not positive', known=endswith)

        p = pierlotx(C, B, A, -115.1, 0, 109.3)  # note CCW order, alpha1, alpha2, alpha3 definition
        self.test(pierlotx.__name__, p.toStr(prec=4), '(2129.3018, 5575.8016, 0)', nl=1)
        p = pierlotx(A, B, C, -115.1, 0, 109.3)
        self.test(pierlotx.__name__, p.toStr(prec=4), '(2128.2026, 4708.1218, 0)')
        p = pierlotx(B, C, A, 0, 0, 109.3)
        self.test(pierlotx.__name__, p.toStr(prec=4), '(1969.0673, 6633.5695, 0)')
        p = pierlotx(C, A, B, 0, 0, 109.3)
        self.test(pierlotx.__name__, p.toStr(prec=4), '(2438.0239, 5094.568, 0)')

        t = tienstra7(A, B, C, 115.1, 135.6, 109.3)  # note alpha, beta and gamma definition
        self.test(tienstra7.__name__, t.pointP.toStr(prec=4), '(2129.3018, 5575.8016, 0)',nl=1)
        self.test(tienstra7.__name__, fstr(t[1:], prec=4), '47.9357, 84.8896, 47.1747, 1581.1388, 2121.3203, 1562.0499')
        self.test(tienstra.__name__, tienstra(A, B, C, 115.1, beta=135.6, gamma=109.3), t, nl=1)  # DEPRECATED

    def testSnellius(self):
        n = snellius3.__name__

        c = triSide(10, 30, PI_4)
        rA, rB, rC, _ = triAngle4(10, 30, c)
        t = triSide4(rA, rB, c)
        self.test(triSide4.__name__, t, '(10.0, 30.0, 0.785398, 8.840862)', nl=1)

        t = snellius3(10, 30, 45, degrees(rA), degrees(rA / 2))
        self.test(n, t, '(17.54582, 38.564239, 46.317675)', nl=1)

        rA, rB, rC, _ = triAngle4(320, 435, 598)
        t = snellius3(320, 435, degrees(rC), 30, 15)
        self.test(n, t, '(844.880591, 571.107418, 835.462796)')

        rA, rB, rC, _ = triAngle4(100, 100, 100)
        t = snellius3(100, 100, degrees(rC), 30, 20)
        self.test(n, t, '(128.557522, 100.0, 187.938524)')

        rA, rB, rC, _ = triAngle4(435, 320, 600)
        t = snellius3(435, 320, degrees(rC), 15, 30)
        self.test(n, t, '(567.480866, 847.344375, 832.446688)')

        rA, rB, rC, _ = triAngle4(1716, 924, 1056)
        t = snellius3(1716, 924, degrees(rC), 0.0, 14.5)  # alpha 0.0 OK
        self.test(n, t, '(4064.197388, 3652.539386, 4988.197388)', known=isWindows)  # (0.0, 3652....)

    def testTri(self):  # for coverage
        t = triAngle(1, 2, 3)
        self.test(triAngle.__name__, t, PI, prec=9, nl=1)
        t = triAngle4(1, 2, EPS0 / 2)  # DEPRECATED
        self.test(triAngle4.__name__, t, '(1.570796, 1.570796, 0.0, 0.0)', known=True)
        t = triAngle5(1, 2, EPS0 / 2)
        self.test(triAngle5.__name__, t, '(1.570796, 1.570796, 0.0, 0.0, 0.0)', known=True)
        t = triAngle5(4, 13, 15)
        self.test(triAngle5.__name__, t, '(0.24871, 0.927295, 1.965587, 1.5, 24.0)', known=True)
        t = triArea(4, 13, 15)
        self.test(triArea.__name__, t, '24.0')
        t = triSide2(0, 2, 0)
        self.test(triSide2.__name__, t, '(2.0, 0.0)')
        t = triSide2(0, 2, PI)
        self.test(triSide2.__name__, t, '(2.0, 3.141593)')

    def testWildberger(self):
        n = wildberger3.__name__

        t = wildberger3(10, 30, 23.994498, 17.139272, 8.569636)
        self.test(n, t, '(17.54582, 38.56424, 46.317675)', nl=1)

        t = wildberger3(320, 435, 598, 30, 15)
        self.test(n, t, '(844.880591, 571.107418, 835.462796)')

        t = wildberger3(100, 100, 100, 30, 20)
        self.test(n, t, '(128.557522, 100.0, 187.938524)')

        t = wildberger3(435, 320, 600, 15, 30)
        self.test(n, t, '(567.480866, 847.344375, 832.446688)')

        t = wildberger3(1716, 924, 1056, 0.0000001, 14.5)
        t = str(t)
        self.test(n, t, '(4064.197343, 3652.539342, 4988.197355)', known=t in
                       ('(4064.197353, 3652.539342, 4988.197355)',
                        '(4064.197363, 3652.539342, 4988.197355)'))


if __name__ == '__main__':

    from pygeodesy.cartesianBase import CartesianBase

    class Cartesian_(CartesianBase):
        pass

    t = Tests(__file__, __version__)
    t.testResections(Cartesian_, Vector3d)
    t.testSnellius()
    t.testWildberger()
    t.testTri()
    t.results()
    t.exit()
