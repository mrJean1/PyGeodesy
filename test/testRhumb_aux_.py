
# -*- coding: utf-8 -*-

# Some basic C{rhumb.aux_} vs C++ C{RhumbSolve} tests.

__all__ = ('Tests',)
__version__ = '24.08.30'

from bases import coverage, _fLate, RhumbSolve, startswith, TestsBase

from pygeodesy import NN, Caps, classname, DIG, Ellipsoid, Ellipsoids, GDict, \
                      isfinite, itemsorted, latDMS, lonDMS, parseDMS, parseDMS2, \
                      RhumbAux, RhumbLineAux, RhumbLineSolve, R_M, \
                      Fwelford, fremainder, pairs
from pygeodesy.interns import _COLON_, _COMMASPACE_, _DOT_
from pygeodesy.rhumb import aux_ as rhumbaux  # NOT deprecated.rhumbaux

_C = _COLON_
_G = '%%.%sg' % (DIG,)


class Tests(TestsBase):

    def testDiffs(self, name, r, rx, nl, e=1e-13, known=False):
        for n, v in itemsorted(r):
            x = rx.get(n, None)
            if x is not None:
                r = abs(v - x) / abs(x or 1)
                k = known or int(v) == int(x) or r < e  # rel error
                self.test(_DOT_(name, n), v, x, fmt=_G, error=r, known=k, nl=nl)
                nl = 0

    def testDirect(self, E, debug=False):
        self.subtitle(rhumbaux, 'DirectX vs ...')
        R = E.rhumbaux

        R.exact = False
        r = R.Direct(40.6, -73.8, -92.38889, 12782581.068)
        self.testDiffs(R.Direct.__name__, r, GDict(lat1=40.6, lat2=35.79999595,
                                                   lon1=-73.8, lon2=140.30000410,
                                                   azi12=-92.38889, s12=12782581.068), 0, e=90)  # XXX neg!
        R.debug = debug
        R.exact = True
        # GeographicLib.RhumbSolve example, _DEBUG_ALL results
        rX = GDict(a=6378137, b=6356752.31424518,  # WGS84
                   f=3.35281066474748e-03, f1=9.966471893352525e-01,
                   e=8.181919084262149e-02, e2=6.694379990141316e-03,
                   L=10001965.7293127, k2=-0.00673949674227643,
                   lat1=40.6, lat2=71.688899882813, azi12=51.0,
                   lon1=-73.8, lon2=0.255519824423359,  # m12=0, M12=1, M21=1,
                   mu1=40.457426625097, mu2=71.6026636715225, mu12=31.1452370464256,
                   psi1=44.2483764794879, salp=0.777145961456971, calp=0.629320391049837,
                   s12=5500000, S12=44095641862956.1)

        r = R.Direct(40.6, -73.8, 51, 5.5e6, outmask=R.ALL)  # from JFK about NE
        self.testDiffs('GDict', r, rX, 1, e=1e-11)  # Windows lon2=0.2555..23445
#       self.test('iteration', r.iteration, r.iteration)

        rl = R.Line(40.6, -73.8, 51)  # coverage
        n = classname(rl)
        t = rl.toStr()
        self.test(n, t, t, nl=1)
        self.test(classname(R), rl.rhumb.toRepr(), R.toRepr())
        self.test(rl.__class__.isLoxodrome.name, rl.isLoxodrome, True)

        a = 51.0
        r = rl.ArcPosition(a, rl.ALL)  # coverage
        self.testDiffs(rl.ArcPosition.__name__, r, GDict(lat1=40.6, lat2=72.635128,
                                                   lon1=-73.8, lon2=4.068528, S12=46665957571716.4,
                                                   azi12=51.0, a12=a, s12=5667780.579944), 1, e=1e-5)
        s = r.s12
        r = rl.Position(s, rl.ALL)  # coverage
        self.testDiffs(rl.Position.__name__, r, GDict(lat1=40.6, lat2=72.635128,
                                                lon1=-73.8, lon2=4.068528, S12=46665957571716.4,
                                                azi12=51.0, a12=a, s12=s), 1, e=1e-5)

        s = 12782581.068
        r = R.Direct(40.6, -73.8, -92.38889, s)  # coverage
        self.testDiffs(R.Direct.__name__, r, GDict(lat1=40.6, lat2=35.8,
                                                   lon1=-73.8, lon2=140.3,
                                                   azi12=-92.38889, a12=115.02062, s12=s), 1, e=1e-5)
        a = r.a12
        r = R.ArcDirect(40.6, -73.8, -92.38889, a)  # coverage
        self.testDiffs(R.ArcDirect.__name__, r, GDict(lat1=40.6, lat2=35.8,
                                                      lon1=-73.8, lon2=140.3,
                                                      azi12=-92.38889, a12=a, s12=s), 1, e=1e-5)

        r = R.Direct8(40.6, -73.8, -92.38889, 12782581.068)  # coverage
        t = str(r)
        self.test(R.Direct8.__name__, t, t, nl=1)
        t = str(r.toDirect9Tuple())  # coverage
        self.test(r.toDirect9Tuple.__name__, t, t)

        t = R.Line(40.6, -73.8, 51, caps=R.STANDARD)  # coverage
        t = str(r)
        self.test(R.Line.__name__, t, t)  # == DirectLine

        r = RhumbLineAux(R, 40.6, -73.8, 51, name='Test')  # coverage
        self.test(RhumbLineAux.__name__, r, R.Line(40.6, -73.8, 51), nl=1)
        r = R.DirectLine(35.8, 140.3, -51)  # coverage
        self.test(R.DirectLine.__name__, r, R.Line(35.8, 140.3, -51),)

        if RhumbSolve:
            S = E.rhumbsolve
            S.reverse2 = True
            t = S.Direct3(40.6, -73.8, 51, 5.5e6).toStr()
            self.test('Direct3', t, '(71.6889, 0.25552, 231.0)')
            S.reverse2 = False
            t = S.Direct3(40.6, -73.8, 51, 5.5e6).toStr()
            self.test('Direct3', t, '(71.6889, 0.25552, 51.0)')

            s = S.Direct(40.6, -73.8, 51, 5.5e6, outmask=S.ALL)
            self.testDiffs('RhumbSolve', rX, s, 1, e=9)  # XXX FIX
            self.test('iteration', s.iteration, r.iteration)
            # extreme ob- and prolate
            for f in range(-7, 10):  # -9, -8 throw an Ellipsoid.e21 AssertionError
                try:
                    f *= 0.1
                    r = R.classof(E.a, f).Direct(40.6, -73.8, 51, 5.5e6, outmask=R.ALL)
                    s = S.classof(E.a, f).Direct(40.6, -73.8, 51, 5.5e6, outmask=S.ALL)
                    self.testDiffs(_fLate(f), r, s, 1, e=9)  # XXX FIX
                except AssertionError:  # eps for f < -0.7
                    pass

    def testInverse(self, E, debug=False):
        self.subtitle(rhumbaux, 'InverseX vs ...')
        R = E.rhumbaux

        R.exact = False
        r = R.Inverse(40.6, -73.8, 35.8, 140.3, outmask=R.ALL)  # JFK to Tokyo Narita
        self.testDiffs(R.Inverse.__name__, r, GDict(lat1=40.6, lat2=35.8,
                                                    lon1=-73.8, lon2=140.3,
                                                    azi12=-92.38889, s12=1282.19384,  # 12782581.0676842,
                                                    S12=21207525604650.8), 0, e=1e-2)  # XXX neg!

        # GeographicLib.RhumbSolve example, _DEBUG_ALL results
        rX = GDict(a=6378137, b=6356752.31424518,  # WGS84
                   f=3.35281066474748e-03, f1=9.966471893352525e-01,
                   e=8.181919084262149e-02, e2=6.694379990141316e-03,
                   L=10001965.7293127, k2=-0.00673949674227643,
                   lat1=40.6, lat2=51.6,  azi12=77.7683897102557,
                   lon1=-73.8, lon2=-0.5, lon12=73.3,  # m12=0, M12=1, M21=1,
                   psi1=44.2483764794879, psi2=60.1387327632216, psi12=15.8903562837337,
                   s12=5771083.38332803, S12=37395209100030.4)

        R.debug = debug
        R.exact = True
        r = R.Inverse(40.6, -73.8, 51.6, -0.5, outmask=R.ALL)  # JFK to LHR
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

        P = Ellipsoid(E.b, E.a, name='_Prolate').rhumbaux  # '_...' for iOS
        t = str(P.Inverse(40.6, -73.8, 51.6, -0.5))  # coverage
        self.test(P.Inverse.__name__, t, t, nl=1)
        r = P.Inverse8(40.6, -73.8, 51.6, -0.5)  # coverage
        t = str(r)
        self.test(P.Inverse8.__name__, t, t)
        t = str(r.toInverse10Tuple())
        self.test(r.toInverse10Tuple.__name__, t, t)
        rl = R.InverseLine(51.6, -0.5, 40.6, -73.8)
        t = str(rl.azi12)
        self.test(R.InverseLine.__name__, t, t)

        if RhumbSolve:
            S = E.rhumbsolve
            t = S.Inverse1(40.6, -73.8, 51.6, -0.5)
            self.test('Inverse1', t, '51.9295425', prec=7, nl=1)
            t = S.Inverse3(40.6, -73.8, 51.6, -0.5).toStr()
            self.test('Inverse3', t, '(5771083.383328, 77.76839, 77.76839)')

            s = S.Inverse(40.6, -73.8, 51.6, -0.5, outmask=S.ALL)
            self.testDiffs('RhumbSolve', rX, s, 1, e=1e-5, known=True)  # XXX FIX
            self.test('iteration', s.iteration, s.iteration)

            # extreme ob- and prolate
            for f in range(-7, 10):  # -9, -8 throw an Ellipsoid.e21 AssertionError
                try:
                    f *= 0.1
                    r = R.classof(E.a, f).Inverse(40.6, -73.8, 51.6, -0.5, outmask=R.ALL)
                    s = S.classof(E.a, f).Inverse(40.6, -73.8, 51.6, -0.5, outmask=S.ALL)
                    self.testDiffs(_fLate(f), s, r, 1, e=1e-5, known=True)  # XXX FIX
                except AssertionError:  # eps for f < -0.7
                    pass

    def testRhumbLine(self, E):
        R = E.rhumbaux
        self.test('R.exact', R.exact, True, nl=1)
        R.exact = x = 0
        self.test('R.exact', R.exact, bool(x))
        self.test('R', repr(R), '''RhumbAux(RAorder=None, TMorder=6, ellipsoid=Ellipsoid(name='WGS84',''', known=startswith)
        rl = R.Line(1, 2, 3)
        R.exact = x = 1
        self.test('R.exact', R.exact, bool(x), nl=1)
        self.test('R.Line.exact', rl.exact, bool(x))
        self.test('R.Line', repr(rl), '''RhumbLineAux(TMorder=6, azi12=3.0, exact=True, lat1=1.0, lon1=2.0, rhumb=RhumbAux(RAorder=None, TMorder=6, ellipsoid=Ellipsoid(name='WGS84',''', known=startswith)
#       t = R.orders(4, 8)
#       self.test(R.orders.__name__, str(t), '(6, 6)')
#       t = R.orders(6, 6)
#       self.test(R.orders.__name__, str(t), '(4, 8)', nt=1)

        R.TMorder = 7
        for exact in (False, True, None):
            nonexact = exact is None
            # <https://SourceForge.net/p/geographiclib/discussion/1026620/thread/2ddc295e/>
            r = R.Line(30, 0, 45, caps=R.LINE_OFF)
            for est in (None,) if nonexact else (1e6, None):
                p = r.PlumbTo(60, 0, exact=exact, est=est)
                t = p.toRepr()
                self.test('PlumbTo(exact=%s, est=%s)' % (exact, est), t, t, nl=1)
                if nonexact:
                    self.test('a02',   p.a02,      17.798332, prec=6)
                    self.test('s02',   p.s02, 1977981.142985, prec=6)
                    self.test('s12',   p.s12, 2169465.957531, prec=6)
                    self.test('azi02', p.azi02,   135.000,    prec=3)
                else:
                    self.test('a02',   p.a02,      17.967658, prec=6)
                    self.test('s02',   p.s02, 1997960.116871, prec=6)
                    self.test('s12',   p.s12, 3083112.636236, prec=6)
                    self.test('azi0',  p.azi0,    113.736,    prec=3)
                    azi2 = p.at + p.azi12
                    self.test('azi2',    azi2,    135.000,    prec=3)
                self.test('iteration', p.iteration, p.iteration)

            t = str(r.Intersecant2(60, 0, radius=r.degrees2m(30)))
            self.test('Intersecant2', t, t, nl=1)

            # stats
            m = n = j = 0
            s = Fwelford()
            r = R.Line(20, 0, 0, caps=R.LINE_OFF)
            for d in range(0, 361, 60 if coverage else 6):
                r.azi12 = d
                p = r.PlumbTo(0, 40, exact=exact)
                i = p.iteration
#               t = p.toRepr()
#               self.test('at %d PlumbTo' % (d,), t, t)
                t = r.distance2(p.lat2, p.lon2)
                z = t.initial
#               t = t.toRepr()
#               self.test('at %d distance2' % (d,), t, t)
#               self.test('at %d iteration' % (d,), i, i)
                j = max(i, j)
                if isfinite(z):  # XXX
                    # azi difference with d
                    if z < 0:
                        z += 360
                    z = fremainder(z - d, 360)
                    if z > 90:
                        z -= 180
                    elif z < -90:
                        z += 180
                    m  = max(m, z)
                    n  = min(n, z)
                    s += z
            d =  dict(exact=exact, min=n, mean=s.fmean(), stdev=s.fstdev(), max=m, iteration=j)
            t = _COMMASPACE_(*pairs(d, prec=6))
            z =  p.azi02 if nonexact else p.azi0
            self.test('azi0*=%.3f' % (z,), t, t)
        t = r.xTM.toRepr()  # coverage
        self.test('xTM', t, t, nl=1)

        # <https://www.MathWorks.com/help/map/ref/rhxrh.html>
        R = RhumbAux(R_M, 0)  # sphere
#       R.TMorder = 5
        r = R.Line(10, -56,  35)
        s = R.Line( 0, -10, 310)
        p = r.Intersection(s)
        t = p.toRepr()
        self.test('Intersection', t, '(26.9774, -43.4088)', nl=1, known=True)
        t = r.PlumbTo(p.lat2, p.lon2).toRepr()
        self.test('PlumbTo', t, t)
        t = s.PlumbTo(p.lat2, p.lon2).toRepr()
        self.test('PlumbTo', t, t)
        t = r.xTM.toRepr()  # coverage
        self.test('xTM', t, t)

        # <https://lost-contact.MIT.edu/afs/inf.ed.ac.uk/group/teaching/matlab-help/R2018a/help/map/calculate-intersection-of-rhumb-line-tracks.html>
        R = RhumbAux(R_M, 0)  # sphere
#       R.TMorder = 7
        r = R.Line(37, -76,  90)
        s = R.Line(15, -17, 315)
        p = r.Intersection(s)
        t = p.toRepr()
        self.test('Intersection', t, '(37.0, -41.7028)', nl=1, known=True)
        t = r.PlumbTo(p.lat2, p.lon2).toRepr()
        self.test('PlumbTo', t, t)
        t = s.PlumbTo(p.lat2, p.lon2).toRepr()
        self.test('PlumbTo', t, t)
        t = r.xTM.toRepr()  # coverage
        self.test('xTM', t, t)

        if RhumbSolve:  # coverage
            for S in (E.rhumbsolve.Line(1, 2, 3), RhumbLineSolve(E.rhumbsolve, 1, 2, 3)):
                t = S.toStr()
                self.test('toStr', t, t, nl=1)
                self.test('lat1', S.lat1, S.lat1)
                self.test('lon1', S.lon1, S.lon1)
                self.test('a', S.a, S.a)
                self.test('f', S.f, S.f)
                t = S.Position(1e6).toStr()
                self.test('Position', t, t)
                S.prec = t = 9
                self.test('prec', S.prec, t)
                S.reverse2 = t = True
                self.test('reverse2', S.reverse2, t)
                S.unroll = t = True
                self.test('unroll', S.unroll, t)
                S.verbose = t = True
                self.test('verbose', S.verbose, t)


if __name__ == '__main__':

    from sys import argv

    _debug = '-d' in argv or '--debug' in argv

    E = Ellipsoids.WGS84
    t = Tests(__file__, __version__, rhumbaux)

    try:
        t.testDirect(E, debug=_debug)
    except ImportError as x:
        t.skip(str(x), n=20)
    try:
        t.testInverse(E, debug=_debug)
    except ImportError as x:
        t.skip(str(x), n=60)
    try:
        t.testRhumbLine(E)
    except ImportError as x:
        t.skip(str(x), n=4)

    t.results()
    t.exit()
