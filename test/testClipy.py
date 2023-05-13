
# -*- coding: utf-8 -*-

# Test L{clipy} module.

__all__ = ('Tests',)
__version__ = '23.03.27'

from bases import TestsBase

from pygeodesy import F_D, F__F_, boundsOf, clipCS4, \
               ClipError, clipFHP4, ClipFHP4Tuple, \
               clipGH4, clipLB6, clipSH, clipSH3


def _lles3(*lles):
    ll1, _ = lles[-1]
    for ll2, e in lles:
        yield ll1, ll2, e
        ll1 = ll2


def _llsTrue3(lls):
    ll1 = lls[-1]
    for ll2 in lls:
        yield ll1, ll2, True
        ll1 = ll2


def _4Ts(*cs):
    for i, c in enumerate(cs):
        for xy in c.split(','):  # _COMMA_
            x, y = map(int, xy.strip().split())
            yield ClipFHP4Tuple(y, x, 0, i)


class Tests(TestsBase):

    def testClipy(self, module):  # MCCABE 15

        self.subtitle(module)
        LatLon = module.LatLon

        # <https://www.CS.Helsinki.Fi/group/goa/viewing/leikkaus/example1.html>
        ll, ur = LatLon(0, 0), LatLon(10, 10)
        ps = LatLon(3, -5), LatLon(5, 5), LatLon(5, 5), LatLon(9, 15)
        t = tuple(clipCS4(ps, ll, ur, closed=True, inull=True))
        p1, p2, i, j = t[0]  # closing edge
        self.test('clipCS4.p1', p1, '07.5°N, 010.0°E')
        self.test('clipCS4.p2', p2, '04.5°N, 000.0°E')
        self.test('clipCS4.i', i, 3)
        self.test('clipCS4.j', j, 0)
        p1, p2, i, j = t[2]  # only null edge
        self.test('clipCS4.p1', p1, '05.0°N, 005.0°E')
        self.test('clipCS4.p2', p2, '05.0°N, 005.0°E')
        self.test('clipCS4.i', i, 1)
        self.test('clipCS4.j', j, 2)

        t = tuple(clipLB6(ps, ll, ur, closed=True, inull=True))
        p1, p2, i, fi, fj, j = t[0]  # closing edge
        self.test('clipLB6.p1', p1, '07.5°N, 010.0°E')
        self.test('clipLB6.p2', p2, '04.5°N, 000.0°E')
        self.test('clipLB6.i',   i, 3)
        self.test('clipLB6.fi', fi, 3.25, prec=2)
        self.test('clipLB6.fi', fi.fractional(ps, LatLon=LatLon), p1)
        self.test('clipLB6.fj', fj, 3.75, prec=2)
        self.test('clipLB6.fj', fj.fractional(ps, LatLon=LatLon), p2)
        self.test('clipLB6.j',   j, 0)
        self.test('clipLB6.fin', fi.fin, 4)
        p1, p2, i, fi, fj, j = t[2]  # only null edge
        self.test('clipLB6.p1', p1, '05.0°N, 005.0°E')
        self.test('clipLB6.p2', p2, '05.0°N, 005.0°E')
        self.test('clipLB6.i',   i, 1)
        self.test('clipLB6.fi', fi, 1.00, prec=2)
        self.test('clipLB6.fi', fi.fractional(ps, LatLon=LatLon), p1)
        self.test('clipLB6.fj', fj, 2.00, prec=2)
        self.test('clipLB6.fj', fj.fractional(ps, LatLon=LatLon), p2)
        self.test('clipLB6.j',   j, 2)
        self.test('clipLB6.fin', fi.fin, 4)

        ll, ur = LatLon(60, 70), LatLon(70, 130)
        ps = LatLon(20, 30), LatLon(80, 170)
        for p1, p2, i, j in clipCS4(ps, ll, ur):
            self.test('clipCS4.p1', p1, '60.0°N, 123.333333°E')
            self.test('clipCS4.p2', p2, '62.857143°N, 130.0°E')
            self.test('clipCS4.i', i, 0)
            self.test('clipCS4.j', j, 1)

        for p1, p2, i, fi, fj, j in clipLB6(ps, ll, ur):
            self.test('clipLB6.p1', p1, '60.0°N, 123.333333°E')
            self.test('clipLB6.p2', p2, '62.857143°N, 130.0°E')
            self.test('clipLB6.i',   i, 0)
            self.test('clipLB6.fi', fi, 0.666667, prec=6)
            self.test('clipLB6.fi', fi.fractional(ps, LatLon=LatLon), p1)
            self.test('clipLB6.fj', fj, 0.714286, prec=6)
            self.test('clipLB6.fj', fj.fractional(ps, LatLon=LatLon), p2)
            self.test('clipLB6.j',   j, 1)
            self.test('clipLB6.fin', fi.fin, 0)

        ll, ur = LatLon(15, 15), LatLon(20, 20)
        ps = LatLon(15, 10), LatLon(25, 20), LatLon(20, 30)
        for p1, p2, i, j in clipCS4(ps, ll, ur, closed=True, inull=False):
            self.test('clipCS4.p1', p1, '17.5°N, 020.0°E')
            self.test('clipCS4.p2', p2, '16.25°N, 015.0°E')
            self.test('clipCS4.i', i, 2)  # closing edge
            self.test('clipCS4.j', j, 0)

        for p1, p2, i, fi, fj, j in clipLB6(ps, ll, ur, closed=True, inull=False):
            self.test('clipLB6.p1', p1, '17.5°N, 020.0°E')
            self.test('clipLB6.p2', p2, '16.25°N, 015.0°E')
            self.test('clipLB6.i',   i, 2)  # closing edge
            self.test('clipLB6.fi', fi, 2.500, prec=3)  # closing edge
            self.test('clipLB6.fi', fi.fractional(ps, LatLon=LatLon), p1)
            self.test('clipLB6.fj', fj, 2.750, prec=3)  # closing edge
            self.test('clipLB6.fj', fj.fractional(ps, LatLon=LatLon), p2)
            self.test('clipLB6.j',   j, 0)  # closing edge
            self.test('clipLB6.fin', fi.fin, 3)

        # ps = LatLon(15, 10), LatLon(25, 20), LatLon(20, 30)
        sh = tuple(clipSH(ps, (ll, ur)))
        self.testClipSH_('clipSH1.', sh, [LatLon(20, 20), LatLon(17.5, 20), LatLon(16.25, 15), LatLon(20, 15)])

        # ps = LatLon(15, 10), LatLon(25, 20), LatLon(20, 30)
        cs = LatLon(30, 10), LatLon(30, 30), LatLon(10, 20)
        for r in ('', 'reversed.'):
            sh = tuple(clipSH(ps, cs, closed=False))
            self.testClipSH_('clipSH2.' + r, sh, [LatLon(18.571, 24.286), LatLon(16.667, 16.667), LatLon(20, 15), LatLon(25, 20), LatLon(22, 26)])
            cs = tuple(reversed(cs))

            # ps = LatLon(15, 10), LatLon(25, 20), LatLon(20, 30)
            # cs = LatLon(30, 10), LatLon(30, 30), LatLon(10, 20)
            sh = clipSH3(ps, cs, closed=True)
            for r3, x3 in zip(sh, _lles3((LatLon(16.667, 16.667), True),
                                         (LatLon(20, 15), False),
                                         (LatLon(25, 20), True),
                                         (LatLon(22, 26), True),
                                         (LatLon(18.571, 24.286), False))):
                self.testClipSH_('clipSH3.' + r, r3[:2], x3[:2])
                self.test('clipSH3.edge.' + r, r3[2], x3[2])

        # ps = LatLon(15, 10), LatLon(25, 20), LatLon(20, 30)
        cs = LatLon(30, -10), LatLon(30, -30), LatLon(10, -20)  # all outside
        for r in ('', '.reversed'):
            sh = tuple(clipSH(ps, cs, closed=True))
            self.test('clipSH.allout' + r, sh, ())
            sh = tuple(clipSH3(ps, cs, closed=True))
            self.test('clipSH3.allout' + r, sh, ())
            cs = tuple(reversed(cs))

        # ps = LatLon(15, 10), LatLon(25, 20), LatLon(20, 30)
        cs = LatLon(15, 10), LatLon(30, 10), LatLon(30, 30), LatLon(15, 30)  # all inside
        for r in ('', 'reversed.'):
            sh = tuple(clipSH(ps, cs, closed=True))
            self.testClipSH_('clipSH.allin.' + r, sh, (ps[-1],) + ps)  # close ps
            sh = tuple(clipSH3(ps, cs, closed=False))
            for r3, x3 in zip(sh, _llsTrue3(ps)):
                self.testClipSH_('clipSH3.allin.' + r, r3[:2], x3[:2])
                self.test('clipSH3.edge.' + r, r3[2], x3[2])
            cs = tuple(reversed(cs))

        # ps = LatLon(15, 10), LatLon(25, 20), LatLon(20, 30)
        cs = LatLon(10, 10), LatLon(20, 20), LatLon(10, 20), LatLon(20, 10)  # warped
        for r in ('', 'reversed.'):
            try:
                # use list to force exception, see
                # <https://RickardLindberg.me/writing/bitten-by-python-generators>
                sh = list(clipSH(ps, cs))
                t = ClipError.__name__
            except ClipError as x:
                t = sh = str(x)  # .split(':')[0]
            self.test('clipSH.warped' + r, sh, t)
            try:
                # use list to force exception, see
                # <https://RickardLindberg.me/writing/bitten-by-python-generators>
                sh = list(clipSH3(ps, cs))
                t = ClipError.__name__
            except ClipError as x:
                t = sh = str(x)  # .split(':')[0]
            self.test('clipSH3.warped' + r, sh, t)
            cs = tuple(reversed(cs))

        cs = boundsOf(cs)
        ps = boundsOf(ps)
        self.test(boundsOf.__name__, cs, '(10.0, 10.0, 20.0, 20.0)')
        self.test(boundsOf.__name__, ps, '(15.0, 10.0, 25.0, 30.0)')
        self.test('enclosures', cs.enclosures(ps), '(5.0, 0.0, -5.0, -10.0)')
        self.test('overlap', cs.overlap(ps), '(15.0, 10.0, 20.0, 20.0)')

        i = 3  # <https://GitHub.com/mdabdk/sutherland-hodgman>
        for ps, cs, x in (([(-1,  1), ( 1, 1), (1, -1), (-1, -1)], [(0, 0), (0, 2), (2, 2), (2, 0)],
                                                                   '0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0'),
                          ([(-1, -1), (-1, 1), (1,  1), ( 1, -1)], [(2, 0), (0, 0), (0, 2), (2, 2)],
                                                                   '0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0'),
                          ([(0, 0), (2, 1), (2, 0)],               [(1, 0.5), (3, 1.5), (3, 0.5)],
                                                                   '1.0, 0.5, 2.0, 1.0, 2.0, 0.5'),
                          ([(0, 3), (0.5, 0.5), (3, 0), (0.5, -0.5), (0, -3), (-0.5, -0.5), (-3, 0), (-0.5, 0.5)],
                                                                   [(-2, -2), (-2, 2), (2, 2), (2, -2)],
                                                                   '-0.2, 2.0, 0.2, 2.0, 0.5, 0.5, 2.0, 0.2, 2.0, -0.2, 0.5, -0.5, 0.2, -2.0, -0.2, -2.0, -0.5, -0.5, -2.0, -0.2, -2.0, 0.2, -0.5, 0.5'),
                          ([(0, 3), (0.5, 0.5), (3, 0), (0.5, -0.5), (0, -3), (-0.5, -0.5), (-3, 0), (-0.5, 0.5)],
                                                                   [(0, 2), (2, -2), (-2, -2)],
                                                                   '-0.33, 1.33, 0.0, 2.0, 0.33, 1.33, 0.5, 0.5, 0.78, 0.44, 1.18, -0.36, 0.5, -0.5, 0.2, -2.0, -0.2, -2.0, -0.5, -0.5, -1.18, -0.36, -0.78, 0.44, -0.5, 0.5'),
                          ([(86, -174), (86, -173), (75, -82), (73, -65), (77, -30)],
                                                                   [(28, -95), (28, -94), (29, -94), (29, -95)],  # CCW
                                                                   None),  # courtesy of U{christophelebrun<https://GitHub.com/mrJean1/PyGeodesy/issues/61>}
                          ([(86, -174), (86, -173), (75, -82), (73, -65), (77, -30)],
                                                                   [(29, -95), (29, -94), (28, -94), (28, -95)],  # CW
                                                                   None),  # courtesy of U{christophelebrun<https://GitHub.com/mrJean1/PyGeodesy/issues/61>}
                          ([(28, -95), (28, -94), (29, -94)],      [(30, -94), (29, -95), (30, -95)], None),
                          ([(30, -94), (29, -95), (30, -95)],      [(28, -95), (28, -94), (29, -94)], None)):
            sh = clipSH((LatLon(*ll) for ll in ps),
                        (LatLon(*ll) for ll in cs))
            sh = ', '.join(ll.toStr(form=F__F_, prec=2) for ll in sh) or None
            i += 1
            self.test('clipSH' + str(i), sh, x)

        p = LatLon(0, 0, height=1.),  LatLon(7, 5, height=2.), LatLon(0, 10, height=3.)  # (0, 0)
        q = LatLon(10, 0, height=1.), LatLon(3, 5, height=2.), LatLon(10, 10, height=3.)  # (5, 0)

        r = clipGH4(p, q, raiser=True)
        self.test(clipGH4.__name__, tuple(r), '(ClipGH4Tuple(lat=5.0, lon=3.571429, height=1.714286, clipid=0), ClipGH4Tuple(lat=7.0, lon=5.0, height=2.0, clipid=0), ClipGH4Tuple(lat=5.0, lon=6.428571, height=2.285714, clipid=0), ClipGH4Tuple(lat=3.0, lon=5.0, height=2.0, clipid=0))', nl=1)

        r = clipFHP4(p, q)
        self.test(clipFHP4.__name__, tuple(r), '(ClipFHP4Tuple(lat=7.0, lon=5.0, height=2.0, clipid=0), ClipFHP4Tuple(lat=5.0, lon=6.428571, height=2.285714, clipid=0), ClipFHP4Tuple(lat=3.0, lon=5.0, height=2.0, clipid=0), ClipFHP4Tuple(lat=5.0, lon=3.571429, height=1.714286, clipid=0))', nl=1)

        # see <https://www.sciencedirect.com/science/article/pii/S259014861930007X?via%3Dihub>, Figs <https://ars.els-cdn.com/content/image/1-s2.0-S259014861930007X-mmc1.zip>
        p = LatLon(0, 0), LatLon(4, 2), LatLon(8, 2), LatLon(6, 7), LatLon(8, 10), LatLon(4, 12), LatLon(1, 12), LatLon(-2, 15), LatLon(0, 9)  # Fig 8
        q = LatLon(6, 12), LatLon(6, 3), LatLon(2, 1), LatLon(-1,4), LatLon(3, 8), LatLon(-1, 12)
        r = clipFHP4(p, q)
        self.test('Fig 8', tuple(r), '(ClipFHP4Tuple(lat=4.0, lon=12.0, height=0.0, clipid=0), ClipFHP4Tuple(lat=-1.0, lon=12.0, height=0.0, clipid=0), ClipFHP4Tuple(lat=3.0, lon=8.0, height=0.0, clipid=0), ClipFHP4Tuple(lat=0.0, lon=5.0, height=0.0, clipid=0), ClipFHP4Tuple(lat=0.0, lon=3.0, height=0.0, clipid=0), ClipFHP4Tuple(lat=2.0, lon=1.0, height=0.0, clipid=0), ClipFHP4Tuple(lat=6.0, lon=3.0, height=0.0, clipid=0), ClipFHP4Tuple(lat=6.0, lon=11.0, height=0.0, clipid=0))', nl=1)

        p = LatLon(0, 0), LatLon(0, 6), LatLon(3, 3), LatLon(6, 6), LatLon(9, 3), LatLon(11, 1), LatLon(13, 3), LatLon(13, 0)  # Fig 14
        q = LatLon(0, 0), LatLon(3, 3), LatLon(6, 0), LatLon(9, 3), LatLon(11, 5), LatLon(13, 3), LatLon(13, 6), LatLon(0, 6)
        r = clipFHP4(p, q)
        self.test('Fig 14', tuple(r), '(ClipFHP4Tuple(lat=3.0, lon=3.0, height=0.0, clipid=0), ClipFHP4Tuple(lat=0.0, lon=0.0, height=0.0, clipid=0), ClipFHP4Tuple(lat=0.0, lon=6.0, height=0.0, clipid=0), ClipFHP4Tuple(lat=6.0, lon=6.0, height=0.0, clipid=1), ClipFHP4Tuple(lat=9.0, lon=3.0, height=0.0, clipid=1), ClipFHP4Tuple(lat=6.0, lon=0.0, height=0.0, clipid=1), ClipFHP4Tuple(lat=3.0, lon=3.0, height=0.0, clipid=1))')

        p = LatLon(0, 0), LatLon(0, 4), LatLon(2, 4), LatLon(1, 3), LatLon(1, 1), LatLon(2, 0)  # Fig 15
        q = LatLon(0, 4), LatLon(1, 3), LatLon(1, 1), LatLon(0, 0), LatLon(2, 0), LatLon(2, 4)
        r = clipFHP4(p, q)
        self.test('Fig 15', tuple(r), '(ClipFHP4Tuple(lat=1.0, lon=3.0, height=0.0, clipid=0), ClipFHP4Tuple(lat=0.0, lon=4.0, height=0.0, clipid=0), ClipFHP4Tuple(lat=2.0, lon=4.0, height=0.0, clipid=0), ClipFHP4Tuple(lat=2.0, lon=0.0, height=0.0, clipid=1), ClipFHP4Tuple(lat=0.0, lon=0.0, height=0.0, clipid=1), ClipFHP4Tuple(lat=1.0, lon=1.0, height=0.0, clipid=1))')

        p = LatLon(0, 0), LatLon(0, 4), LatLon(4, 4), LatLon(4, 0)  # Fig 16
        q = LatLon(2, 2), LatLon(0, 2), LatLon(0, 6), LatLon(-2, 6), LatLon(-2, -2), LatLon(0, -2)
        r = clipFHP4(p, q)
        self.test('Fig 16', tuple(r), '(ClipFHP4Tuple(lat=0.0, lon=0.0, height=0.0, clipid=0), ClipFHP4Tuple(lat=1.0, lon=0.0, height=0.0, clipid=0), ClipFHP4Tuple(lat=2.0, lon=2.0, height=0.0, clipid=0), ClipFHP4Tuple(lat=0.0, lon=2.0, height=0.0, clipid=0))')

        p = _4Ts('-10 -10, 10 -10, 10 10, -10 10',  # Fig 18
                 '-8 -8, 8 -8, 8 8, -8 8',
                 '-6 6, -6 -6, 6 -6, 6 6',
                 '4 4, 4 -4, -4 -4, -4 4',
                 '20 -8, 20 8, 12 8, 12 -8',
                 '14 -6, 14 6, 18 6, 18 -6')
        q = _4Ts('-2 -2, -2 -12, -12 -12, -12 -2',
                 '0 0, -14 0, -14 -14, 0 -14',
                 '2 2, -16 2, -16 -16, 2 -16',
                 '-4 4, -18 4, -18 18, -4 18',
                 '-6 6, -16 6, -16 16, -6 16',
                 '-8 8, -14 8, -14 14, -8 14',
                 '-10 10, -12 10, -12 12, -10 12',
                 '6 -8, 6 8, 20 8, 20 -8',
                 '8 -6, 8 6, 18 6, 18 -6',
                 '16 -4, 16 4, 10 4, 10 -4')
        r = clipFHP4(p, q)
        t = tuple((int(x), int(y), clipid) for y, x, _, clipid in r)
        self.test('Fig 18', t, ' '.join('''((-10, -10, 0), (-10, -2, 0), (-8, -2, 0), (-8, -8, 0), (-2, -8, 0), (-2, -10, 0),
                                            (2, -10, 1), (2, -8, 1), (0, -8, 1), (0, -10, 1),
                                            (10, -6, 2), (8, -6, 2), (8, -8, 2), (10, -8, 2),
                                            (10, 8, 3), (8, 8, 3), (8, 6, 3), (10, 6, 3),
                                            (-6, 10, 4), (-6, 8, 4), (-4, 8, 4), (-4, 10, 4),
                                            (-10, 10, 5), (-10, 8, 5), (-8, 8, 5), (-8, 10, 5),
                                            (-10, 4, 6), (-8, 4, 6), (-8, 6, 6), (-10, 6, 6),
                                            (-10, 0, 7), (-8, 0, 7), (-8, 2, 7), (-10, 2, 7),
                                            (-6, 6, 8), (-4, 6, 8), (-4, 4, 8), (-6, 4, 8),
                                            (-6, 0, 9), (-4, 0, 9), (-4, 2, 9), (-6, 2, 9),
                                            (-6, -6, 10), (-2, -6, 10), (-2, -4, 10), (-4, -4, 10), (-4, -2, 10), (-6, -2, 10),
                                            (2, -6, 11), (2, -4, 11), (0, -4, 11), (0, -6, 11),
                                            (12, 8, 12), (20, 8, 12), (20, -8, 12), (12, -8, 12), (12, -6, 12), (18, -6, 12), (18, 6, 12), (12, 6, 12),
                                            (12, -4, 13), (14, -4, 13), (14, 4, 13), (12, 4, 13))'''.split()))

    def testClipSH_(self, text, sh, lls):
        self.test(text + 'len', len(sh), len(lls))
        for i, (r, x) in enumerate(zip(sh, lls)):
            self.test(text + str(i), r.toStr(F_D, prec=3), x.toStr(F_D))
            self.test(text + 'LL', isinstance(r, type(x)), True)


if __name__ == '__main__':

    from pygeodesy import ellipsoidalNvector, ellipsoidalVincenty, \
                          sphericalNvector, sphericalTrigonometry

    t = Tests(__file__, __version__)
    t.testClipy(ellipsoidalNvector)
    t.testClipy(ellipsoidalVincenty)
    t.testClipy(sphericalNvector)
    t.testClipy(sphericalTrigonometry)
    t.results()
    t.exit()
