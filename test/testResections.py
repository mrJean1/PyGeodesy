
# -*- coding: utf-8 -*-

# Test module attributes.

__all__ = ('Tests',)
__version__ = '21.10.29'

from base import isWindows, TestsBase

from pygeodesy import PI_4, cassini, collins, fstr, pierlot, snellius3, tienstra, \
                      triAngle4, triSide, triSide4, Vector3d, wildberger3
from math import degrees


class Tests(TestsBase):

    def testResections(self, C_, V3d):

        # Tienstra example <http://MesaMike.org/geocache/GC1B0Q9/resection-methods.pdf>
        A, B, C = V3d(1000, 5300), V3d(2200, 6300), V3d(3100, 5000)  # DDD.MMSS 115.5220 109.3045

        p = cassini(A, C, B, 109.5125, 115.0889)  # note B center, alpha and beta definition
        self.test(cassini.__name__, p.toStr(prec=4), '(2128.3903, 5578.1443, 0)')
        p = C_(A).cassini(C, B, 109.5125, 115.0889)
        self.test(cassini.__name__, p.toRepr(prec=4), 'Cartesian_(2128.3903, 5578.1443, 0.0)')

        t = collins(A, C, B, 109.5125, 115.0889)  # note B center, alpha and beta definition
        self.test(collins.__name__, t.pointP.toStr(prec=4), '(2128.3903, 5578.1443, 0)')
        self.test(collins.__name__, t.pointH.toStr(prec=4), '(1830.5948, 2576.2429, 0)')
        r = fstr(t[2:], prec=4)
        self.test(collins.__name__, r, '1581.1388, 1562.0499, 2121.3203')
        t = C_(A).collins(C, B, 109.5125, 115.0889)
        self.test(collins.__name__, t.pointP.toRepr(prec=4), 'Cartesian_(2128.3903, 5578.1443, 0.0)')
        self.test(collins.__name__, t.pointH.toRepr(prec=4), 'Cartesian_(1830.5948, 2576.2429, 0.0)')
        self.test(collins.__name__, fstr(t[2:], prec=4), r)

        p = pierlot(C, B, A, 115.0889, 109.5125)  # note CCW order, alpha and beta definition
        self.test(pierlot.__name__, p.toStr(prec=4), '(2128.3903, 5578.1443, 0)')
        p = C_(C).pierlot(B, A, 115.0889, 109.5125)
        self.test(pierlot.__name__, p.toRepr(prec=4), 'Cartesian_(2128.3903, 5578.1443, 0.0)')

        t = tienstra(A, B, C, 115.0889, None, 109.5125)  # note alpha, beta and gamma definition
        self.test(tienstra.__name__, t.pointP.toStr(prec=4), '(2128.3903, 5578.1443, 0)')
        r = fstr(t[1:], prec=4)
        self.test(tienstra.__name__, r, '47.9357, 84.8896, 47.1747, 1581.1388, 2121.3203, 1562.0499')
        t = C_(A).tienstra(B, C, 115.0889, None, 109.5125)
        self.test(tienstra.__name__, t.pointP.toRepr(prec=4), 'Cartesian_(2128.3903, 5578.1443, 0.0)')
        self.test(tienstra.__name__, fstr(t[1:], prec=4), r)

        p = cassini(A, C, B, 109.3, 115.1)  # note B center, alpha and beta definition
        self.test(cassini.__name__, p.toStr(prec=4), '(2129.3018, 5575.8016, 0)')

        t = collins(A, C, B, 109.3, 115.1)  # note B center, alpha and beta definition
        self.test(collins.__name__, t.pointP.toStr(prec=4), '(2129.3018, 5575.8016, 0)')
        self.test(collins.__name__, t.pointH.toStr(prec=4), '(1835.1911, 2563.0708, 0)')
        self.test(collins.__name__, fstr(t[2:], prec=4), '1581.1388, 1562.0499, 2121.3203')

        p = pierlot(C, B, A, 115.1, 109.3)  # note CCW order, alpha and beta definition
        self.test(pierlot.__name__, p.toStr(prec=4), '(2129.3018, 5575.8016, 0)')

        t = tienstra(A, B, C, 115.1, 135.6, 109.3)  # note alpha, beta and gamma definition
        self.test(tienstra.__name__, t.pointP.toStr(prec=4), '(2129.3018, 5575.8016, 0)')
        self.test(tienstra.__name__, fstr(t[1:], prec=4), '47.9357, 84.8896, 47.1747, 1581.1388, 2121.3203, 1562.0499')

    def testSnellius(self):
        n = snellius3.__name__

        c = triSide(10, 30, PI_4)
        rA, rB, rC, _ = triAngle4(10, 30, c)
        t = triSide4(rA, rB, c)
        self.test(triSide4.__name__, t, '(10.0, 30.0, 0.785398, 8.840862)', nl=1)

        t = snellius3(10, 30, 45, degrees(rA), degrees(rA / 2))
        self.test(n, t, '(17.54582, 38.564239, 46.317675)')

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
        self.test(n, t, '(4064.197313, 3652.539342, 4988.197355)')


if __name__ == '__main__':

    from pygeodesy.cartesianBase import CartesianBase

    class Cartesian_(CartesianBase):
        pass

    t = Tests(__file__, __version__)
    t.testResections(Cartesian_, Vector3d)
    t.testSnellius()
    t.testWildberger()
    t.results()
    t.exit()
