
# -*- coding: utf-8 -*-

# Test module L{ellipses}.

__all__ = ('Tests',)
__version__ = '26.03.25'

from bases import startswith, TestsBase

from pygeodesy import Ellipse, Ellipsoids, EPS, elliperim, polar2d


class Tests(TestsBase):

    def testEllipse(self):

        def _a(name, n=14):
            return Ellipse.__name__ + '.' + (name + ' ' * 12)[:n]

        def _n(name):
            return _a('perimeter' + name)

        a = 6378172.0
        for b, x in ((6378102.0, '40075016.6858801'),
                     (a * 0.9,   '38097844.6222377'),
                     (a / 2,     '30897294.5'),
                     (a / 4,     '273573'),
                     (a / 8,     '26106'),
                     (EPS,       '25512')):
            m = 1.0 - (b / a)**2  # e2
            t = a, b, (b / a), m
            self.test('a, b, b/a, m', t, '(6378172.0, ', known=startswith, nl=1)
            E = Ellipse(a, b)
            for p, n in ((E.perimeter2k,  '2k'),
                         (E.perimeter2k_, '2k_'),  # scipy?
                         (E.perimeterAGM, 'AGM'),
                         (E.perimeterHGK, 'HGK'),
                         (E.perimeterCR,  'CR'),
                         (E.perimeterGK,  'GK'),
                         (E.perimeter2RC, '2RC'),
                         (E.perimeter2R,  '2R')):
                if p is not None:
                    self.test(_n(n), p, x, known=startswith, prec=9)

            p = E.perimeter4Arc3
            t = str(p).lstrip('(').rstrip(')')
            self.test(_n('4Arc3'), t, x[:1], known=startswith)

            n = elliperim.__name__ + ' DEPRECATED  '  # for backward compatibility
            x = x.split('.')[0]  # str(int(float(x)))
            self.test(n, elliperim(a, b), x, known=startswith, nl=1)
            self.test(n, elliperim(a, a), 40075236.597, prec=3)
            self.test(n, elliperim(a, 0), a * 4, prec=1)
            self.test(n, elliperim(0, b), b * 4, prec=1)
            self.test(n, elliperim(0, 0), '0.0')

            self.test(_a('apses2'), E.apses2, E.apses2, nl=1)
            m =  Ellipse(b, a).arc
            n = _a('arc')
            x =  E.arc(45)
            self.test(n, x, x, prec=6)
            self.test(n, m(90, 45), x, prec=6)
            x = E.arc(270)
            self.test(n, x, x, prec=6)
            self.test(n, m(360, 90), x, prec=6)

            self.test(_a('area'), E.area, E.area)
            self.test(_a('e '),   E.e, E.e)
            self.test(_a('c '),   E.c, E.c)
            self.test(_a('p '),   E.p, E.p)
            x = E.point(45)
            self.test(_a('point 45'), x, x)
            x = E.polar2d(45)
            self.test(_a('polar2d 45'), x, x)
            self.test(_a('R1'), E.R1, E.R1)
            self.test(_a('R2'), E.R2, E.R2)
            x = E.Roc(60)
            self.test(_a('Roc 60'), x, x)
            x = E.Roc_(-1, 1)
            self.test(_a('Roc_(-1, 1)'), x, x)
            self.test(_a('Rrectifying'), E.Rrectifying, E.Rrectifying)
            x = E.sideOf(0, 1)
            self.test(_a('sideOf'), x, x)
            self.test(_a('sideOf(c, p)'), E.sideOf(E.c, E.p), 0.0)
            x = E.slope(60)
            self.test(_a('slope 60'), x, x)

            x = E.hartzell4(a, b)
            self.test(_a('hartzell4'), x, x)
            x = E.height4(a, b)
            self.test(_a('height4'), x, x)
            x = E.normal3d(45)
            self.test(_a('normal3d'), x, x)
            x = E.normal4(60)
            self.test(_a('normal4'), x, x, nt=1)

            x = E.toRepr(terse=0)
            self.test(_a('toRepr(terse=0)', n=-12), x, x)
            x = E.toEllipsoid().toStr()
            self.test(_a('toEllipsoid', n=-12), x, x)
            x = E.toTriaxial_().toStr()
            self.test(_a('toTriaxial_', n=-12), x, x)

            if not E.isFlat:  # coverage
                cw = tuple(E.points(20, ccw=False, ended=True))
                n  = sum(int(bool(E.sideOf(*t))) for t in cw)
                self.total += len(cw)
                self.test(_a('sideOfs'), n, 0, nl=1)
                self.test(_a('sideOf point'),  E.sideOf(*E.point(45)), '0.0')
                ccw = tuple(E.points(20, ccw=True, ended=True))
                self.test(_a('cw vs ccw'), len(cw), len(ccw))
                t = str(cw == tuple(reversed(ccw)))
                self.test(_a("cw vs ccw'd"), t, True)
                n = 0
                for d in range(0, 360, 15):
                    x, y, _, _ = E.normal4(d, height=10000)
                    _, p = polar2d(x, y, *E.point(d))
                    s = p - E.slope(d)
                    if abs(s - 90) > 1e-11:
                        self.test(_a('slopes at %d' % (d,)), s, 90.0)
                        n += 1
                    else:
                        self.total += 1
                self.test(_a('slopes'), n, 0)

        self.test(_a('isCircular'), E.isCircular, E.isCircular, nl=1)
        self.test(_a('isFlat'),     E.isFlat, E.isFlat)
        self.test(_a('isOblate'),   E.isOblate, E.isOblate)
        self.test(_a('isProlate'),  E.isProlate, E.isProlate, nt=1)

        for n, W in Ellipsoids.items(all=True, asorted=True):
            self.test(n + '.polarimeter/2k', W.polarimeter, W.toEllipse().perimeter2k, prec=6)


if __name__ == '__main__':

    t = Tests(__file__, __version__)
    t.testEllipse()
    t.results()
    t.exit()
