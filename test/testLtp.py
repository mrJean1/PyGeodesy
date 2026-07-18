
# -*- coding: utf-8 -*-

# Test L{ltp} I{local tangent plane} classes.

__all__ = ('Tests',)
__version__ = '26.07.17'

from bases import startswith, TestsBase

from pygeodesy import Aer, Attitude, EcefFarrell21, EcefFarrell22, EcefKarney, \
                      EcefVeness, EcefSudano, Ecef9Tuple, EcefUPC, EcefYou, \
                      Enu, EPS0, Frustum, fstr, LatLon_, LocalCartesian, \
                      Local9Tuple, Ltp, LqRD, Ned, tyr3d, XyzLocal
from pygeodesy.deprecated import EcefCartesian  # DEPRECATED, use L{LocalCartesian}

from math import fabs


def _absdiff(a, *b):
    for a, b in zip(a, b):
        yield abs(a - b)


class Tests(TestsBase):

    def testLtp(self, Ltp, **kwds):

        self.test(Ltp.__name__, kwds, kwds, nl=1)

        # <https://GeographicLib.SourceForge.io/C++/doc/CartConvert.1.html>
        c = Ltp(33, 44, 20, name='Test', **kwds)
        self.test('name', c.name, 'Test')
        t = c.toRepr()
        self.test('toStr', t, c.classname, known=True)
        Sudano = c.ecef.__class__ is EcefSudano

        self.testCopy(c)
        t = Ltp(c)  # like c.copy() or c.dup()
        self.test('New', t, c)

        t = c.forward(33.3, 44.4, 6000)
        self.test('forward', fstr(t[0:3], prec=2), '37288.97, 33374.29, 5783.65')  # 5783.64
        self.test('name', c.name, 'Test')

        t = c.reverse(37288.97, 33374.29, 5783.65)
        self.test('reverse', fstr(t[3:6], prec=2), '33.3, 44.4, 6000.0', known=Sudano)
        self.test('name', c.name, 'Test')

        # <https://SourceForge.net/p/geographiclib/code/ci/release/tree/examples/example-LocalCartesian.cpp>
        c.reset(48 + 50 / 60.0, 2 + 20 / 60.0, name='Paris')
        self.test('name', c.name, 'Paris')
        self.test(c.name, fstr((c.lat0, c.lon0, c.height0), prec=3), '48.833, 2.333, 0.0')

        t = c.forward(LatLon_(50.9, 1.8, name='Calais'))  # Local9Tuple
        self.test('forward', fstr(t[0:3], prec=2), '-37518.64, 229949.65, -4260.43')
        self.test('name', t.name, 'Calais')

        t = c.reverse(-37518.64, 229949.65, -4260.43)  # Local9Tuple
        self.test('reverse', fstr(t[3:6], prec=2), '50.9, 1.8, -0.0', known=Sudano)
        self.test('name', t.name, 'Paris')

        t = c.reverse(-38e3, 230e3, -4e3)
        self.test('reverse', fstr(t[0:3], prec=1), '-38000.0, 230000.0, -4000.0')
        self.test('reverse', fstr(t[3:6], prec=2), '50.9, 1.79, 264.92', known=Sudano)

        t = c.forward(50.9, 1.79, 264.92)  # Local9Tuple
        self.test('forward', fstr(t[0:3], prec=1), '-38223.7, 229964.2, -4000.0')

        # <https://www.MathWorks.com/help/map/ref/enu2geodetic.html>
        Z = Ltp(46.017, 7.750, 1673, name='Zermatt', **kwds)
        self.test('Zermatt', Z.toStr(), c.classname, known=True, nl=1)
        Sudano = Z.ecef.__class__ is EcefSudano
        M = XyzLocal(-7134.8, -4556.3, 2852.4)  # Matterhorn XYZ
        t = Z.reverse(M).toLatLon(datum=None)  # Matterhorn Xyz to LatLon
        self.test('Matterhorn', t.toStr(prec=3), '(45.976, 7.658, 4531.01)', known=Sudano)
        self.test('xyz', Z.forward(t).xyz.toStr(prec=1), '(-7134.8, -4556.3, 2852.4)', known=Sudano)

        t = Z._local2ecef(M, nine=False)  # coverage
        self.test('_local2ecef', fstr(t, prec=3), '4403757.602, 592124.536, 4566652.082', known=Sudano)
        t = Z._local2ecef(M, nine=True)  # coverage
        self.test('_local2ecef', t.toStr(prec=3), Ecef9Tuple.__name__, known=True)
        self.test('_local2ecef', t.__class__.__name__, Ecef9Tuple.__name__)
        t = Z._ecef2local(t, None, {})  # coverage
        self.test('_ecef2local', t.toStr(prec=3), Local9Tuple.__name__, known=True)
        self.test('_ecef2local', t.__class__.__name__, Local9Tuple.__name__)

        t = M.toXyz()  # Matterhorn XYZ/ENU
        self.test('Xyz', t.toStr(prec=3), '(-7134.8, -4556.3, 2852.4, None)', known=Sudano)
        t = Aer(238.08, 18.744, 8876.8).toXyz()  # Matterhorn Aer
        self.test('Aer', t.toStr(prec=3), '(-7134.912, -4444.548, 2852.474, None)', known=Sudano)

        t = M.toEnu(Enu).toXyz()  # Matterhorn XYZ/ENU
        self.test('Enu', t.toStr(prec=3), '(-7134.8, -4556.3, 2852.4, None)', known=Sudano)
        t = Ned(-4556.3, -7134.8, -2852.4)  # Matterhorn Aer
        self.test('Ned', t.toStr(prec=3), '[-4556.3, -7134.8, -2852.4]', known=Sudano)
        t = t.toXyz(Enu)
        self.test('Enu', t.toStr(prec=3), '[-7134.8, -4556.3, 2852.4]', known=Sudano)
        t = t.toXyz(Ned)
        self.test('Ned', t.toStr(prec=3), '[-4556.3, -7134.8, -2852.4]', known=Sudano)

        # courtesy ThomasGreenhill <https://github.com/mrJean1/PyGeodesy/issues/77>
        c = LocalCartesian(90.0, 57.3)
        t = c.reverse(0, 0, 0).toLatLon()
        self.test('lon00', t, '(90.0, 57.3, 0.0, ', known=startswith, nl=1)
        t = c.reverse(0, 0, 0, lon00=3.75).toLatLon()
        self.test('lon00', t, '(90.0, 3.75, 0.0, ', known=startswith)

        f = Frustum(90, 90)  # ltp=None
        self.test(Frustum.__name__, f, '90.0, 90.0', nl=1)
        self.test('hfov', f.hfov, 90.0)
        self.test('vfov', f.vfov, 90.0)

        p = f.footprint5(1000, -90)
        t = p.toStr()
        self.test('footprint', t, t)
        t = p.xyzLocal5().toStr()
        self.test('footprint.xyzLocal5', t, t)
        t = p.toLatLon5(ltp=Ltp(), LatLon=LatLon_).toStr()
        self.test('footprint.toLatLon5', t, t)

        t = f.footprint5(1000, -90, 0, 44.99).toStr()
        self.test('footprint', t, t, nl=1)

        f = Frustum(45, 45)  # ltp=None
        t = f.footprint5(1000, -90, 0, 22.5).toStr()
        self.test('footprint', t, t)
        f = Frustum(45, 45)  # ltp=None
        t = f.footprint5(1000, -179, 0, 22.5).toStr()
        self.test('footprint', t, t)

        a = Attitude(0, tilt=350, roll=340, yaw=-30, name='test')
        self.test(Attitude.__name__, a, '(0.0, -10.0, 330.0, -20.0)', nl=1)
        self.test('alt',  a.alt,  '0.0')
        self.test('tilt', a.tilt, '-10.0')
        self.test('roll', a.roll, '-20.0')
        self.test('yaw',  a.yaw,  '330.0')
        self.test('matrix', a.matrix, '((0.8137976813493737, -0.4409696105298823, -0.3785223063697926),'
                                      ' (0.4698463103929541, 0.8825641192593856, -0.01802831123629725),'
                                      ' (0.3420201433256688, -0.16317591116653488, 0.9254165783983233))', known=True)
        self.test('rotate', a.rotate(1, 1, 1), '(-0.005694, 1.334382, 1.104261)')
        d = tyr3d(tilt=0, yaw=0, roll=0)
        self.test(tyr3d.__name__, d, '(0.0, 0.0, 0.0)')
        d = tyr3d(tilt=90, yaw=0, roll=0)
        self.test(tyr3d.__name__, d, '(0.0, -2.0, 0.0)')
        d = tyr3d(tilt=0, yaw=90, roll=0)
        self.test(tyr3d.__name__, d, '(0.0, -2.0, 0.0)')
        d = tyr3d(tilt=0, yaw=0, roll=90)
        self.test(tyr3d.__name__, d, '(0.0, 0.0, -2.0)')

    def testLqRD(self, LqRD, eps, **kwds):

        def _x(t, x):
            a, b = map(float, x[1:-1].replace(',', '').split()[:2])
            x += ', diff (%.3f, %.3f)' % (t.x - a, t.y - b)
            return x

        self.test(LqRD.__name__, kwds, kwds, nl=1)

        c = LqRD(name='Test', **kwds)
        self.test('name', c.name, 'Test')
        t = c.toRepr()
        self.test('toStr', t, c.classname, known=True)

        t = c.region4()
        self.test('region4', t.toRepr(), 'RD region (latS=50.0, lonW=2.0, latN=56.0, lonE=8.0)')
        t = c.region4(asRD=True)
        self.test('region4RD', t.toRepr(), 'RD region (minRDx=', known=startswith)  # coverage

        self.testCopy(c)
        t = LqRD(c, name=c.name)  # like c.copy() or c.dup()
        self.test('New', t, c, nt=1)

        t = c.forward(c.lat0, c.lon0, height=c.height0_ETRS)  # 52.15616056, 5.38763889
        self.test('forward', t.xyz.toStr(prec=3), _x(t, '(155000.0, 463000.0, 43.0)'), known=True)
        t = c.reverse(*t.xyz)
        self.test('reverse', t.latlonheight.toStr(prec=6), '(52.156161, 5.387639, 43.0)')

        r = c.region4()
        t = c.forward(r.latC, r.lonC)
        self.test('forward', t.xyz.toStr(prec=6), _x(t, '(196139.431611, 557179.026084, -826.266581)'), known=True)  # -27.8, -107.9
        t = c.reverse(*t.xyz)
        self.test('reverse', t.latlonheight.toStr(prec=3), '(53.0, 6.0, 0.0)', known=fabs(t.height) < eps)

        t = c.forward(r.latN, r.lonE)
        self.test('forward', t.xyz.toStr(prec=6), _x(t, '(318159.695836, 894090.746829, -16626.936575)'), known=True)  # -223.7, -634.2
        t = c.reverse(*t.xyz)
        self.test('reverse', t.latlonheight.toStr(prec=6), '(56.0, 8.0, 0.0)', known=fabs(t.height) < eps)

        t = c.forward(r.latN, r.lonW)
        self.test('forward', t.xyz.toStr(prec=2), _x(t, '(-56496.660094, 896141.991394, -18180.08)'), known=True)  # -225, -686.3
        t = c.reverse(*t.xyz)
        self.test('reverse', t.latlonheight.toStr(prec=6), '(56.0, 2.0, 0.0)', known=fabs(t.height) < eps)

        t = c.forward(r.latS, r.lonE)
        self.test('forward', t.xyz.toStr(prec=6), _x(t, '(342349.924377, 226555.86418, -7131.753557)'), known=True)  # -119.6, -0.4
        t = c.reverse(*t.xyz)
        self.test('reverse', t.latlonheight.toStr(prec=6), '(50.0, 8.0, 0.0)', known=fabs(t.height) < eps)

        t = c.forward(r.latS, r.lonW)
        self.test('forward', t.xyz.toStr(prec=2), _x(t, '(-87853.983835, 228817.835659, -8916.47)'), known=True)  # -116.1, 34.9
        t = c.reverse(*t.xyz)
        self.test('reverse', t.latlonheight.toStr(prec=6), '(50.0, 2.0, 0.0)', known=fabs(t.height) < eps)

        t = c.forward(51.728601274, 4.712120126, 301.7981)
        self.test('forward', t.xyz.toStr(prec=2), _x(t, '(108360.88, 415757.28, 258.01)'), known=True)
        t = c.reverse(*t.xyz)
        self.test('reverse', t.latlonheight.toStr(prec=2), '(51.73, 4.71, 301.8)')


if __name__ == '__main__':

    from pygeodesy import Ellipsoids

    E_WGS84 = Ellipsoids.WGS84

    t = Tests(__file__, __version__)
    t.testLtp(EcefCartesian)  # DEPRECATED, use L{LocalCartesian}
    t.testLtp(LocalCartesian)
    t.testLtp(Ltp)
    t.testLtp(Ltp, ecef=EcefKarney(E_WGS84))

    t.testLtp(Ltp, ecef=EcefFarrell21(E_WGS84))
    t.testLtp(Ltp, ecef=EcefFarrell22(E_WGS84))
    t.testLtp(Ltp, ecef=EcefVeness(E_WGS84))
    t.testLtp(Ltp, ecef=EcefSudano(E_WGS84))
    t.testLtp(Ltp, ecef=EcefYou(E_WGS84))

    t.testLqRD(LqRD, EPS0)
    t.testLqRD(LqRD, 1e-8, ecef=EcefUPC(Ellipsoids.Bessel1841))
    t.testLqRD(LqRD, 1e-9, ecef=EcefUPC(Ellipsoids.GRS80))

    t.results()
    t.exit()
