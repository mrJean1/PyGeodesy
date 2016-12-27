
# -*- coding: utf-8 -*-

# Test OSGR functions and methods.

__version__ = '16.12.07'

if __name__ == '__main__':

    from tests import Tests as _Tests

    from geodesy import Datums, ellipsoidalNvector as eN, F_DMS, osgr

    class Tests(_Tests):

        def testOSgr(self):

            # check convertDatum and back
            p = eN.LatLon(51.4778, -0.0016, datum=Datums.WGS84)
            self.test('WGS84', p, '51.4778°N, 000.0016°W')
            r = p.convertDatum(Datums.OSGB36)
            r.height = 0
            self.test('OSGB36', r, '51.477284°N, 000.00002°E')
            r = r.convertDatum(Datums.WGS84)
            r.height = 0
            self.test('WGS84', r, '51.4778°N, 000.0016°W')

            g = osgr.Osgr(651409.903, 313177.270)
            self.test('OSgr1', str(g), 'TG 51409 13177')
            self.test('OSgr1', repr(g), '[G:TG, E:51409, N:13177]')

            p = g.toLatLon(eN.LatLon)
            p.height = 0
            self.test('toLatLon1', p.toStr(F_DMS), '52°39′28.72″N, 001°42′57.74″E', known=True)
            self.test('toLatLon1', p, '52.657977°N, 001.716038°E', known=True)
            r = p.toOsgr()
            self.test('toOsgr1', r.toStr(0), '651409.903, 313177.270', known=True)

            p = g.toLatLon(eN.LatLon, datum=Datums.OSGB36)
            self.test('toLatLon2', p.toStr(F_DMS), '52°39′27.25″N, 001°43′04.47″E', known=True)
            self.test('toLatLon2', p, '52.657568°N, 001.717908°E', known=True)
            r = osgr.toOsgr(p)
            self.test('toOsgr2', r.toStr(0), '651409,313177', known=True)

            p = eN.LatLon(52.65798, 1.71605)
            r = osgr.toOsgr(p)  # TG 51409 13177
            self.test('toOsgr3', str(r), 'TG 51409 13177')

            r = osgr.toOsgr(52.65757, 1.71791, datum=Datums.OSGB36)
            self.test('toOsgr4', str(r), 'TG 51409 13177')

            g = osgr.parseOSGR('TG 48251 11932')
            self.test('OSGR1', str(g), 'TG 48251 11932')
            self.test('OSGR1', repr(g), '[G:TG, E:48251, N:11932]')

            g = osgr.parseOSGR('TG51409 13177')
            self.test('OSGR2', str(g), 'TG 51409 13177')
            self.test('OSGR2', repr(g), '[G:TG, E:51409, N:13177]')

            g = osgr.parseOSGR('TG5140913177')
            self.test('OSGR3', str(g), 'TG 51409 13177')
            self.test('OSGR3', repr(g), '[G:TG, E:51409, N:13177]')

            g = osgr.parseOSGR('651409,313177')
            self.test('OSGR4', str(g), 'TG 51409 13177')
            self.test('OSGR4', repr(g), '[G:TG, E:51409, N:13177]')

            self.test('OSGR5', g.toStr(prec=0), '651409,313177')
            self.test('OSGR5', g.toStr2(prec=-3), '[OSGR:651409.000,313177.000]')

            r = osgr.parseOSGR(g.toStr(prec=-3))
            self.test('OSGR6', r.toStr(prec=0), '651409,313177')

    t = Tests(__file__, __version__, osgr)
    t.testOSgr()
    t.results()
    t.exit()

    # Typical test results (on MacOS X):

    # testing geodesy.osgr version 16.11.11
    # test 1 WGS84: 51.4778°N, 000.0016°W
    # test 2 OSGB36: 51.477284°N, 000.00002°E
    # test 3 WGS84: 51.4778°N, 000.0016°W
    # test 4 OSgr1: TG 51409 13177
    # test 5 OSgr1: [G:TG, E:51409, N:13177]
    # test 6 toLatLon1: 52°39′28.72″N, 001°43′00.63″E  FAILED, KNOWN, expected 52°39′28.72″N, 001°42′57.74″E
    # test 7 toLatLon1: 52.657979°N, 001.716843°E  FAILED, KNOWN, expected 52.657977°N, 001.716038°E
    # test 8 toOsgr1: 651463,313180  FAILED, KNOWN, expected 651409.903, 313177.270
    # test 9 toLatLon2: 52°39′27.25″N, 001°43′07.37″E  FAILED, KNOWN, expected 52°39′27.25″N, 001°43′04.47″E
    # test 10 toLatLon2: 52.65757°N, 001.718713°E  FAILED, KNOWN, expected 52.657568°N, 001.717908°E
    # test 11 toOsgr2: 651463,313180  FAILED, KNOWN, expected 651409,313177
    # test 12 toOsgr3: TG 51409 13177
    # test 13 toOsgr4: TG 51409 13177
    # test 14 OSGR1: TG 48251 11932
    # test 15 OSGR1: [G:TG, E:48251, N:11932]
    # test 16 OSGR2: TG 51409 13177
    # test 17 OSGR2: [G:TG, E:51409, N:13177]
    # test 18 OSGR3: TG 51409 13177
    # test 19 OSGR3: [G:TG, E:51409, N:13177]
    # test 20 OSGR4: TG 51409 13177
    # test 21 OSGR4: [G:TG, E:51409, N:13177]
    # test 22 OSGR5: 651409,313177
    # test 23 OSGR5: [OSGR:651409.000,313177.000]
    # test 24 OSGR6: 651409,313177
    # 6 geodesy.osgr tests (25.0%) FAILED, incl. 6 KNOWN (Python 2.7.10 64bit)

    # testing geodesy.osgr version 16.11.11
    # test 1 WGS84: 51.4778°N, 000.0016°W
    # test 2 OSGB36: 51.477284°N, 000.00002°E
    # test 3 WGS84: 51.4778°N, 000.0016°W
    # test 4 OSgr1: TG 51409 13177
    # test 5 OSgr1: [G:TG, E:51409, N:13177]
    # test 6 toLatLon1: 52°39′28.72″N, 001°43′00.63″E  FAILED, KNOWN, expected 52°39′28.72″N, 001°42′57.74″E
    # test 7 toLatLon1: 52.657979°N, 001.716843°E  FAILED, KNOWN, expected 52.657977°N, 001.716038°E
    # test 8 toOsgr1: 651463,313180  FAILED, KNOWN, expected 651409.903, 313177.270
    # test 9 toLatLon2: 52°39′27.25″N, 001°43′07.37″E  FAILED, KNOWN, expected 52°39′27.25″N, 001°43′04.47″E
    # test 10 toLatLon2: 52.65757°N, 001.718713°E  FAILED, KNOWN, expected 52.657568°N, 001.717908°E
    # test 11 toOsgr2: 651463,313180  FAILED, KNOWN, expected 651409,313177
    # test 12 toOsgr3: TG 51409 13177
    # test 13 toOsgr4: TG 51409 13177
    # test 14 OSGR1: TG 48251 11932
    # test 15 OSGR1: [G:TG, E:48251, N:11932]
    # test 16 OSGR2: TG 51409 13177
    # test 17 OSGR2: [G:TG, E:51409, N:13177]
    # test 18 OSGR3: TG 51409 13177
    # test 19 OSGR3: [G:TG, E:51409, N:13177]
    # test 20 OSGR4: TG 51409 13177
    # test 21 OSGR4: [G:TG, E:51409, N:13177]
    # test 22 OSGR5: 651409,313177
    # test 23 OSGR5: [OSGR:651409.000,313177.000]
    # test 24 OSGR6: 651409,313177
    # 6 geodesy.osgr tests (25.0%) FAILED, incl. 6 KNOWN (Python 3.6.0 64bit)
