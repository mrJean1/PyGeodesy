
# -*- coding: utf-8 -*-

# Test module attributes.

__all__ = ('Tests',)
__version__ = '21.10.23'

from base import TestsBase

from pygeodesy import cassini, collins, fstr, tienstra, \
                      Vector3d as V3d


class Tests(TestsBase):

    def testResections(self):

        # Tienstra example <http://MesaMike.org/geocache/GC1B0Q9/resection-methods.pdf>
        A, B, C = V3d(1000, 5300), V3d(2200, 6300), V3d(3100, 5000)  # DDD.MMSS 115.5220 109.3045

        p = cassini(A, C, B, 109.5125, 115.0889)  # note B center, alpha and beta definition
        self.test(cassini.__name__, p.toStr(prec=4), '(2128.3903, 5578.1443, 0)')
        t = collins(A, C, B, 109.5125, 115.0889)  # note B center, alpha and beta definition
        self.test(collins.__name__, t.pointP.toStr(prec=4), '(2128.3903, 5578.1443, 0)')
        self.test(collins.__name__, t.pointH.toStr(prec=4), '(1830.5948, 2576.2429, 0)')
        self.test(collins.__name__, fstr(t[2:], prec=4), '1581.1388, 1562.0499, 2121.3203')
        t = tienstra(A, B, C, 115.0889, None, 109.5125)  # note alpha, beta and gamma definition
        self.test(tienstra.__name__, t.pointP.toStr(prec=4), '(2128.3903, 5578.1443, 0)')
        self.test(tienstra.__name__, fstr(t[1:], prec=4), '47.9357, 84.8896, 47.1747, 1581.1388, 2121.3203, 1562.0499')

        p = cassini(A, C, B, 109.3, 115.1)  # note B center, alpha and beta definition
        self.test(cassini.__name__, p.toStr(prec=4), '(2129.3018, 5575.8016, 0)')
        t = collins(A, C, B, 109.3, 115.1)  # note B center, alpha and beta definition
        self.test(collins.__name__, t.pointP.toStr(prec=4), '(2129.3018, 5575.8016, 0)')
        self.test(collins.__name__, t.pointH.toStr(prec=4), '(1835.1911, 2563.0708, 0)')
        self.test(collins.__name__, fstr(t[2:], prec=4), '1581.1388, 1562.0499, 2121.3203')
        t = tienstra(A, B, C, 115.1, 135.6, 109.3)  # note alpha, beta and gamma definition
        self.test(tienstra.__name__, t.pointP.toStr(prec=4), '(2129.3018, 5575.8016, 0)')
        self.test(tienstra.__name__, fstr(t[1:], prec=4), '47.9357, 84.8896, 47.1747, 1581.1388, 2121.3203, 1562.0499')


if __name__ == '__main__':

    t = Tests(__file__, __version__)
    t.testResections()
    t.results()
    t.exit()
