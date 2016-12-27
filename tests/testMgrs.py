
# -*- coding: utf-8 -*-

# Test MGRS functions and methods.

__version__ = '16.12.07'

if __name__ == '__main__':

    from tests import Tests as _Tests

    from geodesy import mgrs

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

    t = Tests(__file__, __version__, mgrs)
    t.testMgrs()
    t.results()
    t.exit()

    # Typical test results (on MacOS X):

    # testing geodesy.mgrs version 16.11.11
    # test 1 Mgrs1: 31U DQ 48251 11932
    # test 2 Mgrs1: [Z:31U, G:DQ, E:48251, N:11932]
    # test 3 Mgrs2: 31U DQ 48251 11932
    # test 4 Mgrs2: [Z:31U, G:DQ, E:48251, N:11932]
    # test 5 Mgrs3: 31U DQ 48251 11932
    # test 6 Mgrs3: [Z:31U, G:DQ, E:48251, N:11932]
    # test 7 toUtm1: 31 N 448251 5411932
    # test 8 toUtm1: [Z:31, H:N, E:448251, N:5411932]
    # test 9 toMgrs: 31U DQ 48251 11932
    # all geodesy.mgrs tests passed (Python 2.7.10 64bit)

    # testing geodesy.mgrs version 16.11.11
    # test 1 Mgrs1: 31U DQ 48251 11932
    # test 2 Mgrs1: [Z:31U, G:DQ, E:48251, N:11932]
    # test 3 Mgrs2: 31U DQ 48251 11932
    # test 4 Mgrs2: [Z:31U, G:DQ, E:48251, N:11932]
    # test 5 Mgrs3: 31U DQ 48251 11932
    # test 6 Mgrs3: [Z:31U, G:DQ, E:48251, N:11932]
    # test 7 toUtm1: 31 N 448251 5411932
    # test 8 toUtm1: [Z:31, H:N, E:448251, N:5411932]
    # test 9 toMgrs: 31U DQ 48251 11932
    # all geodesy.mgrs tests passed (Python 3.6.0 64bit)
