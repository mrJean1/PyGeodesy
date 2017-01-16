
# -*- coding: utf-8 -*-

# Test MGRS functions and methods.

__version__ = '17.01.12'

if __name__ == '__main__':

    from tests import Tests as _Tests

    from geodesy import ellipsoidalVincenty, mgrs

    LatLon = ellipsoidalVincenty.LatLon

    class Tests(_Tests):

        def testMgrs(self):

            m = mgrs.Mgrs('31U', 'DQ', 48251, 11932)
            self.test('Mgrs1', str(m), '31U DQ 48251 11932')
            self.test('Mgrs1', repr(m), '[Z:31U, G:DQ, E:48251, N:11932]')

            m = mgrs.parseMGRS('31U DQ 48251, 11932')
            self.test('Mgrs2', str(m), '31U DQ 48251 11932')
            self.test('Mgrs2', repr(m), '[Z:31U, G:DQ, E:48251, N:11932]')

            m = mgrs.parseMGRS('31UDQ4825111932')
            self.test('Mgrs3', str(m), '31U DQ 48251 11932')
            self.test('Mgrs3', repr(m), '[Z:31U, G:DQ, E:48251, N:11932]')

            u = m.toUtm()
            self.test('toUtm1', str(u), '31 N 448251 5411932')
            self.test('toUtm1', repr(u), '[Z:31, H:N, E:448251, N:5411932]')

            m = u.toMgrs()
            self.test('toMgrs', str(m), '31U DQ 48251 11932')

            for lat, lon, x in ((60.0,  1.0, '31V CG 88455 53097'),  # southern Norway
                                (60.0,  3.0, '32V JM 65640 66593'),
                                (60.0,  9.0, '32V NM 00000 51411'),
                                (76.0,  1.0, '31X DE 45999 36099'),  # Svalbard
                                (76.0, 13.0, '33X VE 45999 36099'),
                                (76.0, 25.0, '35X ME 45999 36099'),
                                (76.0, 37.0, '37X DE 45999 36099')):
                p = LatLon(lat, lon)
                m = p.toUtm().toMgrs()
                self.test('toUtm(%s).toMgrs' % (p,), m, x)

    t = Tests(__file__, __version__, mgrs)
    t.testMgrs()
    t.results()
    t.exit()

    # Typical test results (on MacOS X):

    # testing geodesy.mgrs version 17.01.15
    # test 1 Mgrs1: 31U DQ 48251 11932
    # test 2 Mgrs1: [Z:31U, G:DQ, E:48251, N:11932]
    # test 3 Mgrs2: 31U DQ 48251 11932
    # test 4 Mgrs2: [Z:31U, G:DQ, E:48251, N:11932]
    # test 5 Mgrs3: 31U DQ 48251 11932
    # test 6 Mgrs3: [Z:31U, G:DQ, E:48251, N:11932]
    # test 7 toUtm1: 31 N 448251 5411932
    # test 8 toUtm1: [Z:31, H:N, E:448251, N:5411932]
    # test 9 toMgrs: 31U DQ 48251 11932
    # test 10 toUtm(60.0°N, 001.0°E).toMgrs: 31V CG 88455 53097
    # test 11 toUtm(60.0°N, 003.0°E).toMgrs: 32V JM 65640 66593
    # test 12 toUtm(60.0°N, 009.0°E).toMgrs: 32V NM 00000 51411
    # test 13 toUtm(76.0°N, 001.0°E).toMgrs: 31X DE 45999 36099
    # test 14 toUtm(76.0°N, 013.0°E).toMgrs: 33X VE 45999 36099
    # test 15 toUtm(76.0°N, 025.0°E).toMgrs: 35X ME 45999 36099
    # test 16 toUtm(76.0°N, 037.0°E).toMgrs: 37X DE 45999 36099
    # all geodesy.mgrs tests passed (Python 2.7.13 64bit)

    # testing geodesy.mgrs version 17.01.15
    # test 1 Mgrs1: 31U DQ 48251 11932
    # test 2 Mgrs1: [Z:31U, G:DQ, E:48251, N:11932]
    # test 3 Mgrs2: 31U DQ 48251 11932
    # test 4 Mgrs2: [Z:31U, G:DQ, E:48251, N:11932]
    # test 5 Mgrs3: 31U DQ 48251 11932
    # test 6 Mgrs3: [Z:31U, G:DQ, E:48251, N:11932]
    # test 7 toUtm1: 31 N 448251 5411932
    # test 8 toUtm1: [Z:31, H:N, E:448251, N:5411932]
    # test 9 toMgrs: 31U DQ 48251 11932
    # test 10 toUtm(60.0°N, 001.0°E).toMgrs: 31V CG 88455 53097
    # test 11 toUtm(60.0°N, 003.0°E).toMgrs: 32V JM 65640 66593
    # test 12 toUtm(60.0°N, 009.0°E).toMgrs: 32V NM 00000 51411
    # test 13 toUtm(76.0°N, 001.0°E).toMgrs: 31X DE 45999 36099
    # test 14 toUtm(76.0°N, 013.0°E).toMgrs: 33X VE 45999 36099
    # test 15 toUtm(76.0°N, 025.0°E).toMgrs: 35X ME 45999 36099
    # test 16 toUtm(76.0°N, 037.0°E).toMgrs: 37X DE 45999 36099
    # all geodesy.mgrs tests passed (Python 3.6.0 64bit)
