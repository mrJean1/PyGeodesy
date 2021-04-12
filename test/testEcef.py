
# -*- coding: utf-8 -*-

u'''Test Ecef conversions.
'''

__all__ = ('Tests',)
__version__ = '21.04.06'

from base import TestsBase

from pygeodesy import Datums, EcefError, EcefFarrell21, EcefFarrell22, EcefKarney, \
                      EcefMatrix, EcefSudano, EcefVeness, EcefYou, Ellipsoids, \
                      fstr, latDMS, lonDMS, parse3llh, sphericalNvector  # deprecated.nvector

from math import radians


def _known(t, lat, h, d_h=10.0):
    # if lat is off, so is height
    return abs(t.lat - lat) < 0.1 and abs(t.height - h) < d_h


class Tests(TestsBase):

    def testEcef(self, Ecef):

        self.test(Ecef.__name__, '...', '...', nl=1)

        Karney    = Ecef is EcefKarney
        Sudano    = Ecef is EcefSudano
        Farrell_2 = Ecef is EcefFarrell22

        g = Ecef(Datums.WGS84, name='Test')
        self.test('name', g.name, 'Test')

        t = g.toRepr()
        self.test('toStr', t, g.classname, known=True)

        t = Ecef(g.a, g.f, name=g.name)  # coverage
        self.test('a, f', t, t.classname, known=True)

        self.testCopy(g)

        # <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1Geocentric.html>
        t = g.forward(27.99, 86.93, 8820)  # Mt Everest
        self.test('forward', fstr(t[3:6], prec=2), '27.99, 86.93, 8820.0')
        self.test('forward', fstr(t[0:3], prec=1), '302271.4, 5635928.4, 2979666.1')
        self.test('name', t.name, 'Test')

        t = g.reverse(302271.4, 5635928.4, 2979666.1)
        self.test('reverse', fstr(t[0:3], prec=1), '302271.4, 5635928.4, 2979666.1')
        self.test('reverse', fstr(t[3:6], prec=2), '27.99, 86.93, 8820.01', known=Sudano and _known(t, 27.99, 8820))
        self.test('case', t.C, 2 if Karney else (6 if Sudano else 1))
        self.test('iteration', t.iteration, t.iteration)
        self.test('name', t.name, 'Test')

        # <https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1Geocentric.html>
        t = g.reverse(302e3, 5636e3, 2980e3)
        self.test('reverse', fstr(t[0:3], prec=1), '302000.0, 5636000.0, 2980000.0')
        self.test('reverse', fstr(t[3:6], prec=2), '27.99, 86.93, 9027.03', known=Sudano and _known(t, 27.99, 9027))
        self.test('case', t.C, 2 if Karney else (6 if Sudano else 1))
        self.test('iteration', t.iteration, t.iteration)

        t = g.forward(27.99, 86.93, 8820.0)
        self.test('forward', fstr(t[3:6], prec=2), '27.99, 86.93, 8820.0')
        self.test('forward', fstr(t[0:3], prec=2), '302271.43, 5635928.37, 2979666.13')

        # <https://GeographicLib.SourceForge.io/html/CartConvert.1.html>
        t = g.forward(33.3, 44.4, 6000)
        self.test('forward', fstr(t[3:6], prec=2), '33.3, 44.4, 6000.0')
        self.test('forward', fstr(t[0:3], prec=2), '3816209.6, 3737108.55, 3485109.57')

        t = g.reverse(3816209.6, 3737108.55, 3485109.57)
        self.test('reverse', fstr(t[0:3], prec=2), '3816209.6, 3737108.55, 3485109.57')
        self.test('reverse', fstr(t[3:6], prec=3), '33.3, 44.4, 5999.996', known=Sudano and _known(t, 33.3, 5999))
        self.test('case', t.C, 2 if Karney else (6 if Sudano else 1))
        self.test('iteration', t.iteration, t.iteration)

        # <https://GeographicLib.SourceForge.io/html/CartConvert.1.html>
        t = g.reverse(30000, 30000, 0)
        self.test('reverse', fstr(t[0:3], prec=1), '30000.0, 30000.0, 0.0')
        self.test('reverse', fstr(t[3:6], prec=3), '6.483, 45.0, -6335709.726', known=not Karney)
        self.test('case', t.C, 3 if Karney else (1 if Sudano else 1))
        self.test('iteration', t.iteration, t.iteration)

        t = g.forward(6.483, 45.0, -6335709.726)
        self.test('forward', fstr(t[3:6], prec=3), '6.483, 45.0, -6335709.726')
        self.test('forward', fstr(t[0:3], prec=1), '30000.0, 30000.0, -0.0', known=True)

        # <https://Search.ProQuest.com/docview/847292978> pp 113-114
        t = g.reverse(-2578.0e3, -504.9e3, 5792.9e3)
        self.test('Vermeille', t.lonVermeille,      '-168.919', fmt='%.3f')
        t = g.reverse(-2577.1e3, -498.1e3, 5793.9e3)
        self.test('Vermeille', t.lonVermeille + 360, '190.939', fmt='%.3f')

        # Rey-Jer You <https://www.ResearchGate.net/publication/240359424>
        for i, (x, y, z, h) in enumerate(((-2259148.993, 3912960.837, 4488055.516, 1e3),
                                          (-2259502.546, 3913573.210, 4488762.622, 2e3),
                                          (-2259856.100, 3914185.582, 4489469.729, 3e3),
                                          (-2260209.653, 3914797.955, 4490176.836, 4e3),
                                          (-2262330.973, 3918472.189, 4494419.477, 1e4),
                                          (-2265866.507, 3924595.914, 4501490.544, 2e4),
                                          (-2294150.778, 3973585.709, 4558059.087, 1e5),
                                          (-2541638.152, 4402246.414, 5053033.834, 8e5),
                                          (-2612348.830, 4524720.901, 5194455.190, 1e6))):
            i = '-%d' % (i + 1,)
            r = '45.0, 120.0, %.1f' % (h,)  # Zero and First order columns
            t = g.reverse(x, y, z)
            k = Sudano and _known(t, 45, h) or Farrell_2 and _known(t, 45, h, d_h=0.01)
            self.test('reverse' + i, fstr(t[3:6], prec=3), r, known=k)
            f = g.forward(t.lat, t.lon, t.height)
            self.test('forward' + i, fstr(f[0:3], prec=1), fstr((x, y, z), prec=1), known=k)
            self.test('xyzh' + i, fstr(f.xyzh, prec=1), fstr((x, y, z, h), prec=1), known=k)
            f = f.phi, f.lam
            t = radians(t.lat), radians(t.lon)
            self.test('philam' + i, fstr(f, prec=4), fstr(t, prec=4))

        # <https://www.ResearchGate.net/publication/3709199_An_exact_conversion_from_an_Earth-centered_coordinate_system_to_latitude_longitude_and_altitude>
        t = g.reverse(4588301.55696757, 0, 4558059.086984613)
        self.test('sudano', fstr(t[3:6], prec=3), '45.0, 0.0, 100000.0', known=Sudano)

        # <https://www.OrdnanceSurvey.co.UK/documents/resources/guide-coordinate-systems-great-britain.pdf> pp 47
        g = Ecef(Ellipsoids.GRS80, name='OS-UK')
        self.test('name', g.name, 'OS-UK')

        t = g.forward(parse3llh('''53°36′43.1653"N, 001°39′51.9920"W, 299.800'''))
        self.test('forward', fstr(t[3:6], prec=8), '53.61199036, -1.66444222, 299.8')
        self.test('forward', fstr(t[0:3], prec=2), '3790644.9, -110149.21, 5111482.97')

        t = g.reverse(3790644.9, -110149.21, 5111482.97)
        self.test('reverse', fstr(t[0:3], prec=2), '3790644.9, -110149.21, 5111482.97')
        self.test('reverse', fstr(t[3:5], prec=8), '53.61199036, -1.66444223', known=Sudano)
        self.test('reverse.lat', latDMS(t.lat, prec=4), '53°36′43.1653″N', known=Sudano)
        self.test('reverse.lon', lonDMS(t.lon, prec=4), '001°39′51.992″W')
        self.test('reverse.height', fstr(t.height, prec=-3), '299.800', known=Sudano)
        self.test('case', t.C, 2 if Karney else (7 if Sudano else 1))
        self.test('iteration', t.iteration, t.iteration)

        t = g.forward(parse3llh('34 0 0.00174N, 117 20 0.84965W, 251.702'))  # near LA, CA
        self.test('forward', fstr(t[0:3], prec=6), '-2430601.827685, -4702442.703125, 3546587.358103')
        t = g.reverse(t)
        self.test('reverse', fstr(t[3:6], prec=6), '34.0, -117.333569, 251.702', known=not Karney)

        # coverage
        self.test(EcefError.__name__, g.reverse(0, 0, 0), '(0.0, 0.0, ...)', known=True)
        try:
            self.test(EcefError.__name__, g.forward(None, None, None), EcefError.__name__)
        except Exception as x:
            self.test(EcefError.__name__, str(x).split(':')[0], 'lat (None), lon (None) ...', known=True)
        try:
            self.test(Ecef.__name__, Ecef(None), EcefError.__name__)
        except Exception as x:
            self.test(Ecef.__name__, str(x), Ecef.__name__, known=True)

    def testEcefMatrix(self):

        self.test(EcefMatrix.__name__, '...', '...', nl=1)

        # index order in .multiply
        t = tuple(r * 3 + c for r in range(3) for c in range(3))
        self.test('index', t, '(0, 1, 2, 3, 4, 5, 6, 7, 8)')

        M = EcefMatrix(*t)
        self.test('matrix', fstr(M, prec=0), '0, 1, 2, 3, 4, 5, 6, 7, 8')
        t = M.multiply(M)
        self.test('multiply', fstr(t, prec=0), '45, 54, 63, 54, 66, 78, 63, 78, 93')

        self.testCopy(M)

        I = [0] * 9  # PYCHOK I
        I[0] = I[4] = I[8] = 1
        I = EcefMatrix(*I)  # PYCHOK I
        self.test('matrix', fstr(I, prec=0), '1, 0, 0, 0, 1, 0, 0, 0, 1')
        t = I.multiply(I)
        self.test('multiply', fstr(t, prec=0), '1, 0, 0, 0, 1, 0, 0, 0, 1')

        self.testCopy(I)

    def testLatLonEcef(self, module, Ecef):

        C = module.Cartesian
        self.test(module.__name__, C.__name__, C.__name__, nl=1)
        self.test('Ecef', C(0, 0, 0).Ecef, Ecef)

        LL = module.LatLon
        self.test(module.__name__, LL.__name__, LL.__name__)
        ll = LL(48.833, 2.333, name='Paris')
        self.test('Ecef', ll.Ecef, Ecef)

        t = ll.toEcef()
        self.test('forward', fstr(t[3:6], prec=3), '48.833, 2.333, 0.0')
        self.test('forward', fstr(t[0:3], prec=2), '4202946.8, 171232.47, 4778354.17' if ll.isEllipsoidal
                                              else '4190278.55, 170716.35, 4796058.21')
        self.test('name', t.name, 'Paris')

        e = ll.datum.ecef().reverse(t)
        self.test('reverse', fstr(e[3:6], prec=3), '48.833, 2.333, 0.0')
        self.test('name', e.name, 'Paris')

        ll = e.toLatLon(LL)
        self.test('toLatLon', repr(ll), 'LatLon(48°49′58.8″N, 002°19′58.8″E, +0.00m)' if ll.isEllipsoidal
                                   else 'LatLon(48°49′58.8″N, 002°19′58.8″E)')
        self.test('name', ll.name, 'Paris')
        self.test('Ecef', ll.Ecef, Ecef)

        t = e.toLatLon(LatLon=None)
        self.test('to4Tuple', t.classname, 'LatLon4Tuple')
        self.test('to4Tuple', repr(t), 'Paris(lat=48.833, lon=2.333, height=0.0, datum=%r)' % (t.datum,))

        t = e.toLatLon(LatLon=None, datum=None)
        self.test('to3Tuple', t.classname, 'LatLon3Tuple')
        self.test('to3Tuple', repr(t), 'Paris(lat=48.833, lon=2.333, height=0.0)')

        # XXX import NvectorBase from ellipsoidal- or sphericalNvector since nvectorBase is private
        v = e.toVector(getattr(module, 'Nvector', sphericalNvector.NvectorBase))  # DEPRECATED nvector.Nvector
        self.test('toVector', str(v), '(4202946.79528, 171232.46613, 4778354.17)' if ll.isEllipsoidal
                                 else '(4190278.55277, 170716.34863, 4796058.20898)')
        self.test('name', v.name, 'Paris')

        c = e.toCartesian(C)
        self.test('forward', c.toStr(prec=2), '[4202946.8, 171232.47, 4778354.17]' if ll.isEllipsoidal
                                         else '[4190278.55, 170716.35, 4796058.21]')
        self.test('Ecef', c.Ecef, Ecef)


if __name__ == '__main__':

    from pygeodesy import ellipsoidalKarney, ellipsoidalNvector, ellipsoidalVincenty, \
                          sphericalTrigonometry

    t = Tests(__file__, __version__)
    t.testEcef(EcefKarney)
    t.testEcef(EcefFarrell21)
    t.testEcef(EcefFarrell22)
    t.testEcef(EcefVeness)
    t.testEcef(EcefSudano)
    t.testEcef(EcefYou)
    t.testEcefMatrix()
    t.testLatLonEcef(ellipsoidalKarney, EcefKarney)
    t.testLatLonEcef(ellipsoidalNvector, EcefVeness)
    t.testLatLonEcef(ellipsoidalVincenty, EcefVeness)
    t.testLatLonEcef(sphericalNvector, EcefKarney)
    t.testLatLonEcef(sphericalTrigonometry, EcefKarney)
    t.results()
    t.exit()
