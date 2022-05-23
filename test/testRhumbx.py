
# -*- coding: utf-8 -*-

u'''Some basic C{rhumbx} vs C++ C{RhumbSolve} tests.
'''
__all__ = ('Tests',)
__version__ = '22.05.22'

from base import startswith, TestsBase

from pygeodesy import NN, Caps, classname, Ellipsoid, GDict, \
                      latDMS, lonDMS, parseDMS, parseDMS2, RhumbLine
from pygeodesy.interns import DIG, _DOT_

_C = ':'
_G = '%%.%sg' % (DIG,)


def _key2(t2):
    n, v = t2
    return n.lower(), v


class Tests(TestsBase):

    def testDiffs(self, name, r, rx, nl, e=1e-13):
        for n, v in sorted(r.items(), key=_key2):
            x = rx.get(n, None)
            if x is not None:
                k = (abs(v - x) / (x or 1)) < e  # rel error
                self.test(_DOT_(name, n), v, x, fmt=_G, known=k, nl=nl)
                nl = 0

    def testDirect(self, E, debug=False):
        self.subtitle(rhumbx, 'DirectX vs ...')
        R = E.rhumbx

        R.exact = False
        r = R.Direct(40.6, -73.8, -92.38889, 12782581.068)
        self.testDiffs(R.Direct.__name__, r, GDict(lat1=40.6, lat2=35.79999,
                                                   lon1=-73.8, lon2=140.23651,
                                                   azi12=-92.38889, s12=12782581.068), 0, e=1e-5)
        R.debug = debug
        R.exact = True
        # GeographicLib.RhumbSolve example results
        rX = GDict(a=6378137, b=6356752.31424518,  # WGS84
                   f=3.35281066474748e-03, f1=9.966471893352525e-01,
                   e=8.181919084262149e-02, e2=6.694379990141316e-03,
                   L=10001965.7293127, k2=-0.00673949674227643,
                   lat1=40.6, lat2=71.688899882813, azi12=51.0,
                   lon1=-73.8, lon2=0.255519824423374, m12=0, M12=1, M21=1,
                   mu1=40.457426625097, mu2=71.6026636715225, mu12=31.1452370464256,
                   psi1=44.2483764794879, salp=0.777145961456971, calp=0.629320391049837,
                   s12=5500000, S12=44054189889953.1)

        r = R.Direct(40.6, -73.8, 51, 5.5e6, R.ALL)  # from JFK about NE
        self.testDiffs('GDict', r, rX, 1, e=1e-11)  # Windows lon2=0.2555..23445
#       self.test('iteration', r.iteration, r.iteration)

        Rl = R.Line(40.6, -73.8, 51, 5.5e6)  # coverage
        n = classname(Rl)
        t = Rl.toStr()
        self.test(n, t, t, nl=1)
        self.test(classname(R), Rl.rhumb.toRepr(), R.toRepr())
#       r = Rl.Position(Rl.s13, Rl.ALL)  # coverage
#       t = r.toRepr()
#       self.test(n, t, t)

        r = R.Direct(40.6, -73.8, -92.38889, 12782581.068)  # coverage
        self.testDiffs(R.Direct.__name__, r, GDict(lat1=40.6, lat2=35.8,
                                                   lon1=-73.8, lon2=140.3,
                                                   azi12=-92.38889, s12=12782581.068), 1, e=1e-5)

        r = R.Direct7(40.6, -73.8, -92.38889, 12782581.068)  # coverage
        t = str(r)
        self.test(R.Direct7.__name__, t, t)
        t = str(r.toDirect9Tuple())  # coverage
        self.test(r.toDirect9Tuple.__name__, t, t)

        t = R.Line(40.6, -73.8, 51, R.STANDARD)  # coverage
        t = str(r)
        self.test(R.Line.__name__, t, t)  # == DirectLine

        r = RhumbLine(R, 40.6, -73.8, 51, name='Test')  # coverage
        self.test(RhumbLine.__name__, r, R.Line(40.6, -73.8, 51), nl=1)
        r = R.DirectLine(35.8, 140.3, -51)  # coverage
        self.test(R.DirectLine.__name__, r, R.Line(35.8, 140.3, -51),)

    def testInverse(self, E, debug=False):
        self.subtitle(rhumbx, 'InverseX vs ...')
        R = E.rhumbx

        R.exact = False
        r = R.Inverse(40.6, -73.8, 35.8, 140.3, R.ALL)  # JFK to Tokyo Narita
        self.testDiffs(R.Inverse.__name__, r, GDict(lat1=40.6, lat2=35.8,
                                                    lon1=-73.8, lon2=140.3,
                                                    azi12=-92.38889, s12=12776983.7488964,
                                                    sS12=-63760642939073), 0, e=1e-5)

        # GeographicLib.RhumbSolve example results
        rX = GDict(a=6378137, b=6356752.31424518,  # WGS84
                   f=3.35281066474748e-03, f1=9.966471893352525e-01,
                   e=8.181919084262149e-02, e2=6.694379990141316e-03,
                   L=10001965.7293127, k2=-0.00673949674227643,
                   lat1=40.6, lat2=51.6,  azi12=77.7683897102557,
                   lon1=-73.8, lon2=-0.5, lon12=73.3, m12=0, M12=1, M21=1,
                   psi1=44.2483764794879, psi2=60.1387327632216, psi12=15.8903562837337,
                   s12=5771083.38332803, S12=37362984299663.1)

        R.debug = debug
        R.exact = True
        r = R.Inverse(40.6, -73.8, 51.6, -0.5, R.ALL)  # JFK to LHR
        self.testDiffs('GDict', r, rX, 1)

        # <https://GeographicLib.SourceForge.io/C++/doc/RhumbSolve.1.html>
        r = R.Inverse(*(parseDMS2('40:38:23N', '073:46:44W', sep=':') +  # JFK
                        parseDMS2('01:21:33N', '103:59:22E', sep=':')))  # Chanhi Airport, Singapore
        self.testDiffs(R.Inverse.__name__, r, GDict(lat1=40.639722, lat2=1.359167,
                                                    lon=-73.778889, lon2=103.989444,
                                                    azi12=parseDMS('103:34:58.2', sep=':'),
                                                    s12=18523563), 1, e=1e-6)

        # <https://GeographicLib.SourceForge.io/C++/doc/RhumbSolve.1.html>
        L = R.Line(*(parseDMS2('40:38:23N', '073:46:44W', sep=':') +
                     (parseDMS('103:34:58.2', sep=':'),)))
        for i, x in enumerate(('40:38:23.0N 073:46:44.0W 0',
                               '36:24:30.3N 051:28:26.4W 9817078307821',
                               '32:10:26.8N 030:20:57.3W 18224745682005',
                               '27:56:13.2N 010:10:54.2W 25358020327741',
                               '23:41:50.1N 009:12:45.5E 31321269267102',
                               '19:27:18.7N 027:59:22.1E 36195163180159',
                               '15:12:40.2N 046:17:01.1E 40041499143669',
                               '10:57:55.9N 064:12:52.8E 42906570007050',
                               '06:43:07.3N 081:53:28.8E 44823504180200',
                               '02:28:16.2N 099:24:54.5E 45813843358737',
                               '01:46:36.0S 116:52:59.7E 45888525219677')):
            lat2, lon2, S12 = x.split()
            n = '%d,000 Km ' % (i,)
            r = L.Position(i * 2e6, outmask=Caps.LATITUDE_LONGITUDE_AREA)
            self.test(n + 'lat2', latDMS(r.lat2, prec=1, s_D=_C, s_M=_C, s_S=NN), lat2, nl=0 if i else 1)
            self.test(n + 'lon2', lonDMS(r.lon2, prec=1, s_D=_C, s_M=_C, s_S=NN), lon2)
            S = int(S12)
            e = (abs(r.S12 - S) / S) if S else 0
            self.test(n + 'S12 ', int(r.S12), S12, known=e < 2e-3)

        P = Ellipsoid(E.b, E.a, name='_Prolate').rhumbx  # '_...' for iOS
        t = str(P.Inverse(40.6, -73.8, 51.6, -0.5))  # coverage
        self.test(P.Inverse.__name__, t, t, nl=1)
        r = P.Inverse7(40.6, -73.8, 51.6, -0.5)  # coverage
        t = str(r)
        self.test(P.Inverse7.__name__, t, t)
        t = str(r.toInverse10Tuple())
        self.test(r.toInverse10Tuple.__name__, t, t)
        rl = R.InverseLine(51.6, -0.5, 40.6, -73.8)
        t = str(rl.azi12)
        self.test(R.InverseLine.__name__, t, t)

    def testRhumb_Line(self, E):
        R = E.rhumbx
        self.test('E.rhumbx.exact', R.exact, True, nl=1)
        R.exact = x = 0
        self.test('E.rhumbx.exact', R.exact, bool(x))
        self.test('E.rhumbx', repr(R), '''Rhumb(RAorder=6, TMorder=6, ellipsoid=Ellipsoid(name='WGS84',''', known=startswith)
        Rl = R.Line(1, 2, 3)
        R.exact = x = 1
        self.test('E.rhumbx.exact', R.exact, bool(x), nl=1)
        self.test('E.rhumbx.Line.exact', Rl.exact, bool(x))
        self.test('E.rhumbx.Line', repr(Rl), '''RhumbLine(azi12=3.0, exact=True, lat1=1.0, lon1=2.0, rhumb=Rhumb(RAorder=6, TMorder=6, ellipsoid=Ellipsoid(name='WGS84',''', known=startswith)
        t = R.orders(4, 8)
        self.test(R.orders.__name__, str(t), '(6, 6)')
        t = R.orders(6, 6)
        self.test(R.orders.__name__, str(t), '(4, 8)')


if __name__ == '__main__':

    from pygeodesy import Ellipsoids, rhumbx
    from sys import argv

    _debug = '-d' in argv or '--debug' in argv

    E = Ellipsoids.WGS84
    t = Tests(__file__, __version__, rhumbx)

    t.testDirect(E, debug=_debug)
    t.testInverse(E, debug=_debug)
    t.testRhumb_Line(E)

    t.results()
    t.exit()
