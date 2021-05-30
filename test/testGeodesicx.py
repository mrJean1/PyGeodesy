
# -*- coding: utf-8 -*-

u'''Some basic C{geodsicx} vs C++ C{GeographicLib}, C{GeodSolve}
    and Python C{geographiclib} tests.
'''
__all__ = ('Tests',)
__version__ = '21.05.25'

from base import geographiclib, GeodSolve, TestsBase

from pygeodesy import classname, map2, NN
from pygeodesy.interns import DIG, _DOT_

_G = '%%.%sg' % (DIG,)


def _fLate(f):
    t = 'proLate' if f < 0 else \
        ('obLate' if f > 0 else 'sphere')
    return 'f(%.1f)%s' % (f, t)


def _key2(t2):
    n, v = t2
    return n.lower(), v


def _t3(n_p_a):
    return map2(int, map(round, n_p_a))


class Tests(TestsBase):

    def testDiffs(self, name, r, rx, nl):
        for n, v in sorted(r.items(), key=_key2):
            x = rx.get(n, None)
            if x is not None:
                k = (abs(v - x) / (x or 1)) < 1e-13  # rel error
                self.test(_DOT_(name, n), v, x, fmt=_G, known=k, nl=nl)
            nl = 0

    def testDirect(self, E):
        self.subtitle(geodesicx, 'DirectX vs ...')
        # GeographicLib.GeodesicExact example, results from GeoSolve Direct
        rCpp = dict(lat1=40.6, lon1=-73.799999999999997, azi1=51.0,
                    lat2=51.884564505606761, lon2=-1.141172861200829, azi2=107.189397162605886,
                    s12=5500000.0, a12=49.475527463251467, m12=4844148.703101486,
                    M12=0.65091056699808603, M21=0.65122865892196558, S12=39735075134877.094)

        gX = E.geodesicx
        rX = gX.Direct(40.6, -73.8, 51, 5.5e6, gX.ALL)
        self.testDiffs('C++', rX, rCpp, 0)
        gs = gX,

        if geographiclib:
            gP = E.geodesic
            rP = gP.Direct(40.6, -73.8, 51, 5.5e6, gP.ALL)
            self.testDiffs('Python', rX, rP, 1)
            gs += gP,

        if GeodSolve:
            gS = E.geodsolve
            rS = gS.Direct(40.6, -73.8, 51, 5.5e6, gS.ALL)
            self.testDiffs('GeodSolve', rX, rS, 1)
            gs += gS,

            # extreme ob- and prolate
            for f in range(-7, 10):  # -9, -8 throw an Ellipsoid.e12 AssertionError
                try:
                    f *= 0.1
                    rX = gX.classof(E.a, f).Direct(40.6, -73.8, 51, 5.5e6, gX.ALL)
                    rS = gS.classof(E.a, f).Direct(40.6, -73.8, 51, 5.5e6, gS.ALL)
                    self.testDiffs(_fLate(f), rX, rS, 1)
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

    def testInverse(self, E):
        self.subtitle(geodesicx, 'InverseX vs ...')
        # GeographicLib.GeodesicExact example, results from instrumented GeoSolve Inverse
        rCpp = dict(a=6378137, f=3.35281066474748e-03,
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
        rX = gX.Inverse(40.6, -73.8, 51.6, -0.5, gX.ALL)
        self.testDiffs('C++', rX, rCpp, 0)
        gs = gX,

        if geographiclib:
            gP = E.geodesic
            rP = gP.Inverse(40.6, -73.8, 51.6, -0.5, gP.ALL)
            self.testDiffs('Python', rX, rP, 1)
            gs += gP,

        if GeodSolve:
            gS = E.geodsolve
            rS = gS.Inverse(40.6, -73.8, 51.6, -0.5, gS.ALL)
            self.testDiffs('GeodSolve', rX, rS, 1)
            gs += gS,
            # extreme ob- and prolate
            for f in range(-7, 10):  # -9, -8 throw an Ellipsoid.e12 AssertionError
                try:
                    f *= 0.1
                    rX = gX.classof(E.a, f).Inverse(40.6, -73.8, 51.6, -0.5, gX.ALL)
                    rS = gS.classof(E.a, f).Inverse(40.6, -73.8, 51.6, -0.5, gS.ALL)
                    self.testDiffs(_fLate(f), rX, rS, 1)
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

    def testPolygon(self, module, g, nC4=NN, K=False):
        self.subtitle(module, 'Polygon' + str(nC4))

        if nC4:
            g.C4Order = nC4
        p = g.Polygon()  # GeodesicAreaExact(g)

        p.AddPoint( 0,    0)
        p.AddEdge( 90, 1000)
        p.AddEdge(  0, 1000)
        p.AddEdge(-90, 1000)
        self.test('AddEdges',  _t3(p.Compute()),           (4, 4000, 1000000), known=K)
        self.test('TestEdge',  _t3(p.TestEdge(180, 1000)), (5, 4000, 1000000), known=K)

        n = p.Clear()  # .Reset()
        self.test('Clear', n, 0, known=n is None)

        p.AddPoint( 52,   0)  # London
        p.AddPoint( 41, -74)  # New York
        p.AddPoint(-23, -43)  # Rio de Janeiro
        p.AddPoint(-26,  28)  # Johannesburg
        self.test('AddPoints',  _t3(p.Compute()),        (4, 29506941, 65690027591346))
        self.test('TestPoint',  _t3(p.TestPoint(52, 0)), (5, 29506941, 65690027591346))


if __name__ == '__main__':

    from pygeodesy import Ellipsoids, geodesicx, geodsolve, karney

    t = Tests(__file__, __version__, geodesicx)

    E = Ellipsoids.WGS84

    t.testDirect(E)

    t.testInverse(E)

    t.testPolygon(geodesicx, E.geodesicx, nC4=24)
    t.testPolygon(geodesicx, E.geodesicx, nC4=27)
    t.testPolygon(geodesicx, E.geodesicx, nC4=30)

    if GeodSolve:
        t.testPolygon(geodsolve, E.geodsolve)

    if geographiclib:
        t.testPolygon(karney, E.geodesic, K=True)  # XXX geographiclib 1.49 issue?

    t.results()
    t.exit()
