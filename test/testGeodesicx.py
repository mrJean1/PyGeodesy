
# -*- coding: utf-8 -*-

# Some basic C{geodesicx} vs C++ C{GeographicLib}, C{GeodSolve}
# and Python C{geographiclib} tests.

__all__ = ('Tests',)
__version__ = '25.01.05'

from bases import _fLate, GeodSolve, geographiclib, isPython2, TestsBase

from pygeodesy import classname, DIG, Ellipsoid, GDict, GeodesicLineExact, \
                      itemsorted, map2, NN
from pygeodesy.interns import _DOT_

_G = '%%.%sg' % (DIG,)


def _t3(n_p_a):
    return map2(int, map(round, n_p_a))


class Tests(TestsBase):

    def testDiffs(self, name, r, rx, nl, e=1e-13):
        for n, v in itemsorted(r):
            x = rx.get(n, None)
            if x is not None:
                k = (abs(v - x) / (x or 1)) < e  # rel error
                self.test(_DOT_(name, n), v, x, fmt=_G, known=k, nl=nl)
                nl = 0

    def testDirect(self, geodesicx, E, debug=False):
        self.subtitle(geodesicx, 'DirectX vs ...')
        # GeographicLib.GeodesicExact example, results from GeoSolve Direct
        rCpp = GDict(lat1=40.6, lon1=-73.799999999999997, azi1=51.0,
                     lat2=51.884564505606761, lon2=-1.141172861200829, azi2=107.189397162605886,
                     s12=5500000.0, a12=49.475527463251467, m12=4844148.703101486,
                     M12=0.65091056699808603, M21=0.65122865892196558, S12=39735075134877.094)

        gX = E.geodesicx
        gX.debug = debug
        rX = gX.Direct(40.6, -73.8, 51, 5.5e6, gX.ALL)
        self.testDiffs('C++X', rX, rCpp, 0)
        self.test('iteration', rX.iteration, rX.iteration)
        gs = gX,

        if geographiclib:
            gP = E.geodesic
            rP = gP.Direct(40.6, -73.8, 51, 5.5e6, gP.ALL)
            self.testDiffs('Python', rX, rP, 1)
            self.test('iteration', rP.iteration, rP.iteration)
            gs += gP,

        if GeodSolve:
            gS = E.geodsolve
            rS = gS.Direct(40.6, -73.8, 51, 5.5e6, gS.ALL)
            self.testDiffs('GeodSolve', rX, rS, 1)
            self.test('iteration', rS.iteration, rS.iteration)
            gs += gS,

            # extreme ob- and prolate
            for f in range(-7, 10):  # -9, -8 throw an Ellipsoid.e21 AssertionError
                try:
                    f *= 0.1
                    rX = gX.classof(E.a, f).Direct(40.6, -73.8, 51, 5.5e6, gX.ALL)
                    rS = gS.classof(E.a, f).Direct(40.6, -73.8, 51, 5.5e6, gS.ALL)
                    self.testDiffs(_fLate(f), rX, rS, 1, e=1e-4)  # XXX S12 in GeodSolve 2.2
                except AssertionError:  # eps for f < -0.7
                    pass

        glX = gX.DirectLine(40.6, -73.8, 51, 5.5e6)  # coverage
        n = classname(glX)
        S = glX.toStr()
        self.test(n, S, S, nl=1)
        self.test(classname(gX), glX.geodesic.toRepr(), gX.toRepr())
        t = glX._GenPosition(False, glX.s13, glX.ALL)  # coverage
        S = t.toRepr()
        self.test(n, S, S)

        for g in gs:  # coverage
            g.debug = False  # coverage
            n = classname(g)
            r = g._GDictDirect(40.6, -73.8, 51, False, 5.5e6)
            S = r.toStr()
            self.test(n, S, S, nl=1)
            t = r.toDirect9Tuple()
            S = t.toStr()
            self.test(n, S, S)
            self.test(n, t.toGDict(), r, known=True)

        t = gX.ArcDirect(40.6, -73.8, 51, 49.8)  # coverage
        self.testDiffs(gX.ArcDirect.__name__, t, GDict(lat1=40.6, lon1=-73.8,
                                                       lat2=51.7876867, lon2=-0.641731,
                                                       azi1=51.0, azi2=107.5820825,
                                                       a12=49.8, s12=5536073.734393), 1, e=1e-7)
        self.test('iteration', t.iteration, t.iteration)

        t = gX.ArcDirectLine(40.6, -73.8, 51, 49.8, gX.STANDARD)  # coverage
        s = str(t)
        self.test(gX.ArcDirectLine.__name__, s, s, nl=1)
        self.test('iteration', t.iteration, t.iteration)

        t = GeodesicLineExact(gX, 40.6, -73.8, 51, caps=gX.ALL)  # coverage
        self.test(GeodesicLineExact.__name__, t, gX.Line(40.6, -73.8, 51), nl=1)
        self.test('iteration', t.iteration, t.iteration)

    def testInverse(self, geodesicx, E, debug=False):
        self.subtitle(geodesicx, 'InverseX vs ...')
        # GeographicLib.GeodesicExact example, results from instrumented GeoSolve Inverse
        rCpp = GDict(a=6378137, f=3.35281066474748e-03,
                     e=8.181919084262149e-02, e2=6.694379990141316e-03,
                     ep2=6.739496742276434e-03, f1=9.966471893352525e-01,
                     c2x=4.058973249931476e+13, c2=4.058973249931476e+13,
                     sig12=8.716402960622249e-01,
                     slam12=9.578224948453149e-01, clam12=2.873605198497121e-01,
                     sbet1=-7.82676535388952e-01, cbet1=6.224286633434762e-01,
                     sbet2=-6.49513671781464e-01, cbet2=7.603499129801756e-01,
                     ssig1=-9.716339593508425e-01, csig1=-2.364898497530187e-01, dn1=1.002062122905103,
                     ssig2=-8.063222953935967e-01, csig2=5.914764204524144e-01, dn2=1.001420580015174,
                     salp1=7.793257472478793e-01, calp1=6.266190067948263e-01,
                     salp2=9.520131366060522e-01, calp2=-3.060571641532119e-01,
                     somg12=9.583181997903982e-01, comg12=2.857030415492466e-01,
                     A4=1.299900770381802e+11, B41=-1.574868160309433e-01, B42=3.939014756725822e-01,
                     k2=4.373072977143198e-03, alp12=-9.882559303867356e-01,
                     v=3.919390853535099e-17, dv=1.605102632009351,
                     iter=3, eFk2=-4.373072977143198e-03, eFa2=-6.739496742276434e-03,
                     eFcD=7.841136961631642e-01, eFcE=1.572512222993233e+00, eFcH=7.836521215162365e-01,
                     lat1=40.6, lon1=-73.799999999999997, azi1=51.198882845579824,
                     lat2=51.6, lon2=-0.5, azi2=107.821776735514248,
                     s12=5551759.4003186841, a12=49.941310217899037,
                     m12=4877684.6027061976, M12=0.64472969205948238,
                     M21=0.64504567852134398, S12=40041368848742.531)

        gX = E.geodesicx
        gX.debug = debug
        rX = gX.Inverse(40.6, -73.8, 51.6, -0.5, outmask=gX.ALL)
        self.testDiffs('C++X', rX, rCpp, 0)
        self.test('iteration', rX.iteration, rX.iteration)
        gs = gX,

        if geographiclib:
            gP = E.geodesic
            rP = gP.Inverse(40.6, -73.8, 51.6, -0.5, outmask=gP.ALL)
            self.testDiffs('Python', rX, rP, 1)
            self.test('iteration', rP.iteration, rP.iteration)
            gs += gP,

        if GeodSolve:
            gS = E.geodsolve
            rS = gS.Inverse(40.6, -73.8, 51.6, -0.5, outmask=gS.ALL)
            self.testDiffs('GeodSolve', rX, rS, 1)
            self.test('iteration', rS.iteration, rS.iteration)
            gs += gS,
            # extreme ob- and prolate
            for f in range(-7, 10):  # -9, -8 throw an Ellipsoid.e21 AssertionError
                try:
                    f *= 0.1
                    rX = gX.classof(E.a, f).Inverse(40.6, -73.8, 51.6, -0.5, outmask=gX.ALL)
                    rS = gS.classof(E.a, f).Inverse(40.6, -73.8, 51.6, -0.5, outmask=gS.ALL)
                    self.testDiffs(_fLate(f), rX, rS, 1, e=1e-4)  # XXX S12 in GeodSolve 2.2
                except AssertionError:  # eps for f < -0.7
                    pass

        glX = gX.InverseLine(40.6, -73.8, 51.6, -0.5)  # coverage
        n = classname(glX)
        R = glX.toRepr()
        self.test(n, R, R, nl=1)
        self.test(classname(gX), glX.geodesic.toStr(), gX.toStr())
        t = glX._GenPosition(True, glX.a13, glX.ALL)  # coverage
        R = t.toRepr()
        self.test(n, R, R)

        for g in gs:  # coverage
            g.debug = False  # coverage
            n = classname(g)
            r = g._GDictInverse(40.6, -73.8, 51.6, -0.5)
            R = r.toRepr()
            self.test(n, R, R, nl=1)
            t = r.toInverse10Tuple()
            R = t.toRepr()
            self.test(n, R, R)
            self.test(n, t.toGDict(), r, known=True)

        gX = Ellipsoid(E.b, E.a, name='_Prolate').geodesicx  # '_...' for iOS
        t = str(gX.Inverse(40.6, -73.8, 51.6, -0.5))  # coverage
        self.test(gX.Inverse.__name__, t, t, nl=1)
        t = str(gX.Inverse1(40.6, -73.8, 51.6, -0.5))  # coverage
        self.test(gX.Inverse1.__name__, t, t)

    def testPlumbTo(self, geodesicx, E, debug=False):
        self.subtitle(geodesicx, 'PlumbTo ...')
        # see Baselga Moreno, S. & Martinez-Llario, J.C, "Intersection and point-to-line
        # solutions for geodesics on the ellipsoid",  "3. Minimum point-to-line distance ..."
        # <https://riunet.UPV.ES/bitstream/handle/10251/122902/Revised_Manuscript.pdf>
        gX = E.geodesicx  # == geodesicx.GeodesicExact(E)
        gX.debug = debug
        g = gX.InverseLine(52, 5, 51.4, 6).PlumbTo(52, 5.5)
        # {a02: 0.213783, a12: 0.222926, at: -270.0, azi0: -136.00267, azi1: 133.603738, azi2: 133.808744,
        # lat0: 52, lat1: 52.0, lat2: 51.846089, lon0: 5.5, lon1: 5.0, lon2: 5.260428, s02: 23767.724184, s12: 24784.288415}
        self.test('lat2',      g.lat2,   51.846089, prec=6)  # == 51°50′45.9212″N  51°50'45.9212"
        self.test('lon2',      g.lon2,    5.260428, prec=6)  # ==  5°15′37.5426″E   5°15'37.5426"
        self.test('s12',       g.s12, 24784.288415, prec=6)  # vs       Table 4  sAX 24,784.2886
        self.test('at',        g.at,   -270.0, prec=6)
        self.test('iteration', g.iteration, 6, known=True, nt=1)
        g = gX.InverseLine(42, 29, 39, -77).PlumbTo(64, -22)
        # {a02: 9.086172, a12: 35.341373, at: 270.0, azi0: 179.771122, azi1: -50.693753, azi2: -90.174696,
        # lat0: 64, lat1: 42.0, lat2: 54.928531, lon0: -22, lon1: 29.0, lon2: -21.937291, s02: 1010585.998837, s12: 3928788.572003}
        self.test('lat2',      g.lat2,     54.928531, prec=6)  # == 54°55′42.7116″N  54°55'42.7134"
        self.test('lon2',      g.lon2,    -21.937291, prec=6)  # == 21°56′14.2476″W -21°56'14.2477"
        self.test('s12',       g.s12, 3928788.572003, prec=6)  # vs    Table 5  sAX 3,928,857.7554
        self.test('at',        g.at,      270.0, prec=6)
        self.test('iteration', g.iteration, 27, known=True, nt=1)
        g = gX.InverseLine(42, 29, -35, -70).PlumbTo(64, -22)  # 12,200 Km > 10K Km, too long!
        # {a02: 35.339953, a12: 9.725372, at: 270.864557, azi0: 118.886127, azi1: -112.646672, azi2: -119.919402,
        # lat0: 64, lat1: 42.0, lat2: 37.673899, lon0: -22, lon1: 29.0, lon2: 17.677028, s02: 3928936.714834, s12: 1080484.867328}
        self.test('lat2',      g.lat2,     37.976217, prec=6, known=int(g.lat2)==37)  # == 37°58′34.3812″N  37°58'41.2236"
        self.test('lon2',      g.lon2,     18.344820, prec=6, known=int(g.lon2)==18)  # == 18°20′41.3520″E  18°20'56.6279"
        self.test('s12',       g.s12, 1012790.599291, prec=6,  # varies 5+ meter        vs    Table 6  sAX 1,012,443.9063
                                                              known=1012789<g.s12<1012798)
        self.test('at',        g.at,      270.005437, prec=6, known=int(g.at)  ==270)
        self.test('iteration', g.iteration, 128, known=True)

    def testPolygon(self, module, g, nC4=NN, K=False):
        self.subtitle(module, 'Polygon' + str(nC4))

        if nC4:
            g.C4order = nC4
        p = g.Polygon()  # GeodesicAreaExact(g)

        p.AddPoint( 0,    0)
        self.test('Compute', _t3(p.Compute()), '(1, 0, 0)')  # coverage
        p.AddEdge( 90, 1000)
        p.AddEdge(  0, 1000)
        p.AddEdge(-90, 1000)
        self.test('AddEdges',  _t3(p.Compute()),           '(4, 4000, 1000000)', known=K)
        self.test('TestEdge',  _t3(p.TestEdge(180, 1000)), '(5, 4000, 1000000)', known=K)

        n = p.Clear()  # .Reset()
        self.test('Clear', n, 0, known=n is None)
        self.test('TestPoint', _t3(p.TestPoint(52, 0)), '(1, 0, 0)')  # coverage

        p.AddPoint( 52,   0)  # London
        p.AddPoint( 41, -74)  # New York
        p.AddPoint(-23, -43)  # Rio de Janeiro
        p.AddPoint(-26,  28)  # Johannesburg
        self.test('AddPoints',  _t3(p.Compute()),        '(4, 29506941, 65690027591346)', known=isPython2)
        self.test('TestPoint',  _t3(p.TestPoint(52, 0)), '(5, 29506941, 65690027591346)', known=isPython2)

        if not K:  # coverage, but not for pygeodesy.karney.PolygonArea
            t = p.toStr()
            self.test(p.toStr.__name__, t, t)


if __name__ == '__main__':

    from pygeodesy import Ellipsoids, geodesicw, geodesicx, geodsolve
    from sys import argv

    _debug = '-d' in argv or '--debug' in argv

    t = Tests(__file__, __version__, geodesicx)

    E = Ellipsoids.WGS84

    t.testDirect(geodesicx, E, debug=_debug)

    t.testInverse(geodesicx, E, debug=_debug)

    t.testPolygon(geodesicx, E.geodesicx, nC4=24)
    t.testPolygon(geodesicx, E.geodesicx, nC4=27)
    t.testPolygon(geodesicx, E.geodesicx, nC4=30)

    if GeodSolve:
        t.testPolygon(geodsolve, E.geodsolve)

    if geographiclib:
        t.testPolygon(geodesicw, E.geodesic, K=True)  # XXX geographiclib 1.49 issue?

    t.testPlumbTo(geodesicx, E, debug=_debug)

    t.results()
    t.exit()
