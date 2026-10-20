
# -*- coding: utf-8 -*-

# Some L{lom} C{LabordeOblique} Mercator tests.

__all__ = ('Tests',)
__version__ = '26.10.20'

from bases import TestsBase, startswith

from pygeodesy import Ellipsoids as Es, grades2degrees, \
                      lom, LabordeOblique

from math import fabs


def _err(f, x, y):
    return max(fabs(f.x - x), fabs(f.y - y))


def _str(t):
    return str(t).split('LabordeOblique')[0] + '...'


class Tests(TestsBase):

    def testLaborde(self):

        p = LabordeOblique(grades2degrees(21), grades2degrees(51.5969213),
                           azimuth=grades2degrees(-21), earth=Es.Intl1924,
                           kC=0.9995, x0=4e5, y0=8e5)  # IOGP page 6
        self.test('laborde', p, '18.9, 46.4372, 341.1, ', known=startswith)
        self.test('latlonC', p.latlonC, '(18.9, 46.437229)')
        self.test('philamC', p.philamC, '(0.329867, 0.810483)')
        self.test('azimuthC', p.azimuthC, '-18.90', prec=2)
        self.test('kC', p.kC, '0.9995')
        self.test('latlonParis', p.latlonParis, '(48.8051, 2.337229)')
        self.test('philamParis', p.philamParis, '(0.85181, 0.040792)')
        self.test('xy0', p.xy0, '(400000.0, 800000.0)')

        f = p.forward('16 11 23.280S', '44 27 27.260E')
        self.test('forward', _str(f), '(188367.887184, 1098721.824251, 0.998343, -0.046949, -16.1898, 44.457572, 0.9995, ...', known=startswith, nl=1)
        r = p.reverse(f.x, f.y)
        self.test('reverse', _str(r), '(188367.887184, 1098721.824251, 0.959983, -0.033272, -16.1898, 44.457572, 0.9995, ...', known=startswith)

        r = p.reverse(188333.848, 1098841.091)
        self.test('reverse', _str(r), '(188333.848, 1098841.091, 0.959988, -0.033277, -16.188719, 44.457265, 0.9995, ...', known=startswith, nl=1)

        f = p.forward('18 54S', '47 30E')  # IOGP page 66, Comparing the Hotine Oblique Mercator (variant B) ...
        self.test('forward', _str(f), '(511921.054, 799665.521, ...', error=_err(f, 511921.054, 799665.521), known=True, nl=1)  # IOGP N/X? 799665.521, E/Y? 511921.054, ... XXX
        f = p.forward('16 12S', '44 24E')
        self.test('forward', _str(f), '(182184.982, 1097651.4447, ...', error=_err(f, 182184.982, 1097651.4447), known=True)  # IOGP N/X? 1097651.447, E/Y? 182184.982, ... XXX
        f = p.forward('25 40S', '45 18E')
        self.test('forward', _str(f), '(285294.334, 50636.222, ...', error=_err(f, 285294.334, 50636.222), known=True)  # IOGP N/X? 50636.222, E/Y? 285294.334, ... XXX
        f = p.forward('12 0S', '49 12E')
        self.test('forward', _str(f), '(701354.05, 1561109.146, ...', error=_err(f, 701354.05, 1561109.146), known=True)  # IOGP N/X? 1561109.146, E/Y? 701354.056, ... XXX


if __name__ == '__main__':

    t = Tests(__file__, __version__, lom)
    t.testLaborde()
    t.results()
    t.exit()
