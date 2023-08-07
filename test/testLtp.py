
# -*- coding: utf-8 -*-

# Test L{ltp} I{local tangent plane} classes.

__all__ = ('Tests',)
__version__ = '23.08.07'

from bases import startswith, TestsBase

from pygeodesy import Aer, Attitude, ChLV, ChLVa, ChLVe, \
                      EcefFarrell21, EcefFarrell22, EcefKarney, \
                      EcefVeness, EcefSudano, Ecef9Tuple, EcefYou, \
                      Enu, Frustum, fstr, LatLon_, LocalCartesian, \
                      Local9Tuple, Ltp, Ned, tyr3d, XyzLocal, latDMS, lonDMS, \
                      EcefCartesian  # DEPRECATED, use L{LocalCartesian}


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

    def testChLV(self, ChLV_):

        def _trim(t, neg0=False):
            t = str(t)
            if neg0:
                t = t.replace(' -0.0,', ' 0.0,')
            i = t.find('EcefMatrix')
            return t[:i] + '...'

        self.test(ChLV_.__name__, '...', '...', nl=1)

        c = ChLV_(name='Test')
        self.test('name', c.name, 'Test')

        t = c.forward(46.95108, 7.438637)
        self.test('ChLV_', t.ltp.__class__.__name__, ChLV_.__name__)

        self.test('forward1', _trim(t), '(-72.039994, -147.361444, -49.552111, 46.95108, 7.438637, 0.0, ' if t.isChLV
                                   else '(0.329415, -0.292702, -49.554242, 46.95108, 7.438637, 0.0, ' if t.isChLVa
                                   else '(-72.031251, -147.344948, -49.554242, 46.95108, 7.438637, 0.0, ', known=startswith, nl=1)
        self.test('Y, X, h_', (t.Y, t.X, t.h_), (t.x, t.y, t.z))
        self.test('EN2_LV95', t.EN2_LV95, '(2599927.960006, 1199852.638556)' if t.isChLV
                                     else '(2600000.329415, 1199999.707298)' if t.isChLVa
                                     else '(2599927.968749, 1199852.655052)')
        self.test('yx2_LV03', t.yx2_LV03, '(599927.960006, 199852.638556)' if t.isChLV
                                     else '(600000.329415, 199999.707298)' if t.isChLVa
                                     else '(599927.968749, 199852.655052)')
        r = c.reverse(t)  # .Y, t.X, t.h_)
        self.test('reverse1', _trim(r), '(-72.039994, -147.361444, -49.552111, 46.95108, 7.438637, 0.0, ' if r.isChLV
                                   else '(0.329415, -0.292702, -49.554242, 46.951078, 7.438642, -0.004239, ' if r.isChLVa
                                   else '(-72.031251, -147.344948, -49.554242, 46.95108, 7.438637, 0.000002, ', known=startswith)

        r = c.reverse(700000, 100000, 600) if ChLV_ is ChLVa else c.reverse(2700000, 1200000, 600)
        self.test('reverse2', _trim(r), '(100000.0, 0.0, 600.0, 46.944873, 8.752874, 1431.948128, ' if r.isChLV
                                   else '(100000.0, -100000.0, 600.0, 46.044127, 8.730499, 650.554, ' if r.isChLVa
                                   else '(100000.0, 0.0, 600.0, 46.944869, 8.753274, 648.29, ', known=startswith, nl=1)
        t = c.forward(r.lat, r.lon, r.height)
        self.test('forward2', _trim(t, neg0=True), '(100000.0, 0.0, 600.0, 46.944873, 8.752874, 1431.948128, ' if t.isChLV
                                              else '(99999.933937, -100000.44412, 600.003469, 46.044127, 8.730499, 650.554, ' if t.isChLVa
                                              else '(100000.000001, 0.0, 600.012265, 46.944869, 8.753274, 648.29, ', known=startswith)
        self.test('Y, X, h_', (t.Y, t.X, t.h_), (t.x, t.y, t.z))
        self.test('EN2_LV95', t.EN2_LV95, '(2700000.0, 1200000.0)' if t.isChLV
                                     else '(2699999.933937, 1099999.55588)' if t.isChLVa
                                     else '(2700000.000001, 1200000.0)')
        self.test('yx2_LV03', t.yx2_LV03, '(700000.0, 200000.0)' if t.isChLV
                                     else '(699999.933937, 99999.55588)' if t.isChLVa
                                     else '(700000.000001, 200000.0)')

        t = c.forward('46 2 38.87', '8 43 49.79', 650.60)
        self.test('forward3', _trim(t), '(99920.639806, -100148.24791, -967.661696, 46.044131, 8.730497, 650.6, ' if t.isChLV
                                   else '(99999.763621, -100000.026905, 600.049476, 46.044131, 8.730497, 650.6, ' if t.isChLVa
                                   else '(99914.740455, -100135.079382, 600.049476, 46.044131, 8.730497, 650.6, ', known=startswith, nl=1)
        self.test('Y, X, h_', (t.Y, t.X, t.h_), (t.x, t.y, t.z))
        self.test('EN2_LV95', t.EN2_LV95, '(2699920.639806, 1099851.75209)' if t.isChLV
                                     else '(2699999.763621, 1099999.973095)' if t.isChLVa
                                     else '(2699914.740455, 1099864.920618)')
        self.test('yx2_LV03', t.yx2_LV03, '(699920.639806, 99851.75209)' if t.isChLV
                                     else '(699999.763621, 99999.973095)' if t.isChLVa
                                     else '(699914.740455, 99864.920618)')
        r = c.reverse(t)  # .Y, t.X, t.h_)
        self.test('reverse3', _trim(r), '(99920.639806, -100148.24791, -967.661696, 46.044131, 8.730497, 650.6, ' if r.isChLV
                                   else '(99999.763621, -100000.026905, 600.049476, 46.044127, 8.730496, 650.603479, ' if r.isChLVa
                                   else '(99914.740455, -100135.079382, 600.049476, 46.044131, 8.730497, 650.607608, ', known=startswith)

        t = c.forward('''47° 03' 28.95659233"''', '''8° 29' 11.11127154"''')  # Rigi
        self.test('forward4', _trim(t), '(79527.502386, 12274.804229, -556.312155, 47.058043, 8.48642, 0.0, ' if t.isChLV
                                   else '(79602.736359, 12421.914221, -48.257243, 47.058043, 8.48642, 0.0, ' if t.isChLVa
                                   else '(79520.050001, 12273.439993, -48.257243, 47.058043, 8.48642, 0.0, ', known=startswith, nl=1)
        self.test('Y, X, h_', (t.Y, t.X, t.h_), (t.x, t.y, t.z))
        self.test('EN2_LV95', t.EN2_LV95, '(2679527.502386, 1212274.804229)' if t.isChLV
                                     else '(2679602.736359, 1212421.914221)' if t.isChLVa
                                     else '(2679520.050001, 1212273.439993)')
        self.test('yx2_LV03', t.yx2_LV03, '(679527.502386, 212274.804229)' if t.isChLV
                                     else '(679602.736359, 212421.914221)' if t.isChLVa
                                     else '(679520.050001, 212273.439993)')
        r = c.reverse(t)  # .Y, t.X, t.h_)
        self.test('reverse4', _trim(r), '(79527.502386, 12274.804229, -556.312155, 47.058043, 8.48642, ' if t.isChLV
                                   else '(79602.736359, 12421.914221, -48.257243, 47.058038, 8.486421, 0.00853, ' if r.isChLVa
                                   else '(79520.050001, 12273.439993, -48.257243, 47.058043, 8.48642, 0.012933, ', known=startswith, nt=1)

        _ = '''swisstopo GeoSuite, 10/2/2022 6:56 PM, REFRAME
        <https://www.SwissTopo.admin.ch/en/maps-data-online/calculation-services/reframe.html>
        Change planimetry:      yes
                Input:          Plane coordinates LV03 (CH1903)
                                Format: Military (Berne=600/200 km)
                Output:         Global coordinates ETRF93/CHTRF95 (ETRS89/CHTRS95/≈WGS84)
                                Format: Geographic Lon/Lat/Alt in [°]  (LTOP ED)
        Change altimetry:       yes
                Input:          National leveling network LN02 (leveled heights)
                Output:         Ellipsoidal heights (GRS80)
        Parameters and datasets used:
                1: Htrans       HTRANS (swisstopo 2005)
                2: CHGeo        CHGeo2004 (swisstopo 2004)
                3: Fineltra     CHENyx06 (swisstopo 03.2007)
                4: GPSRef       CHTRF95/CH1903+ (swisstopo 1997)
        '''
        i = 1  # XXX (lon, lat) all wrong the ChLVe, but height is correct
        for e, n, h_, t in ((690112.2987, 167571.01,   2316.17,  ( 8.61588452888461, 46.6532932166132, 2367.70367439557)),
                            (811208.67,   184302.34,  1495.279,  (10.2046874099038,  46.7763864160066, 1546.62256775331)),
                            (584293.29,   249277.33,   890.93,   ( 7.23058814336896, 47.3941353658685,  940.206020880491)),
                            (783515.164,  187452.333, 1595.183,  ( 9.84351348200287, 46.8129146868336, 1645.59174753912)),
                            (533184.98,   152611.74,   410.638,  ( 6.56789363190021, 46.5214647527993,  460.491580192)),
                            (700000,      100000,      600,      ( 8.73049738375132, 46.0441209510033,  650.012275829911))):
            r = c.reverse(e, n, h_, name='REFRAME' + str(i))  # LV03 falsed!
            d = max(_absdiff(t, r.lon, r.lat))
            self.test(r.name, (r.lon, r.lat, r.height), t, known=d < 0.005)  # 5 mm
            d = ' %.3g m' % (d,)
            self.test(r.name, d, d)
            i += 1

        # <https://eMuseum.GGGS.CH/literatur-lv/liste-Dateien/1967_Bolliger_a.pdf> p 151
        if ChLV_ is ChLVe:
            r, g = c.reverse(*ChLV.false2(130689.06, -83702.04, LV95=True), name='Bollinger', gamma=True)
        else:
            r, g = c.reverse(*ChLV.false2(130689.06, -83702.04, LV95=True), name='Bollinger'), None
        self.test(r.name + '.lat', latDMS(r.lat, prec=5),  "46°11′12.03969″N", known=abs(r.lat - 46.18667769) < 0.007, nl=1)
        self.test(r.name + '.lon', lonDMS(r.lon, prec=5), "009°07′57.29664″E", known=abs(r.lon -  9.13258239) < 0.007)
        if g is not None:
            self.test(r.name + '.gamma', g, g, nt=1)

        # <https://www.SwissTopo.admin.ch/en/maps-data-online/calculation-services/navref.html>
        for e, n, h_, t, m in ((623706.01,   96020.08,  1687.567, (46.015314107,  7.74720273,   1742.449),   'Zermatt'),
                               (602030.68,  191775.03,   897.915, (46.878427446,  7.466218797,   953.0973),  'Zimmerwald'),
                               (617306.3,   268507.30,   456.064, (47.568308647,  7.6695485844,  897.4060),  'Chrischona'),
                               (776668.105, 265372.681, 1042.624, (47.5163201698, 9.7843787372, 3868.8599),  'Pfaender'),
                               (497313.292, 145625.438, 1207.434, (46.4555740686, 6.1033162697, 2313.94497), 'La Givrine'),
                               (722758.81,   87649.67,  1636.6,   (45.9311184445, 9.0216020471, 3855.4332),  'Monte Generoso')):
            r = c.reverse(e, n, h_, name=m)
            d = max(_absdiff(t, r.lat, r.lon))
            self.test(r.name, (r.lat, r.lon, r.height), t, known=d < 0.005)  # 5 mm
            d = ' %.3g m' % (d,)
            self.test(r.name, d, d)


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

    t.testChLV(ChLV)
    t.testChLV(ChLVa)
    t.testChLV(ChLVe)

    t.results()
    t.exit()
