
# -*- coding: utf-8 -*-

# Test gars module.

__all__ = ('Tests',)
__version__ = '20.01.22'

from base import TestsBase

from pygeodesy import degDMS, fstr, gars, Garef, S_MIN


class Tests(TestsBase):

    def testCodec3(self, g, x, prec=4):
        self.test('codec3', Garef(g), g)
        t = gars.decode3(g)
        self.test('decode3', fstr(t, prec=prec), x)
        self.test('encode', gars.encode(*t), g)

    def testGars(self, LL):

        # Karney's geographiclib/1.49/examples/example-GARS.cpp
        # <https://SourceForge.net/p/geographiclib/code/ci/release/tree/examples/example-GARS.cpp>
        g = Garef('57.64911, 10.40744', precision=2)
        self.test('Garef', g, '381NH45')
        self.test('Garef', g.toStr(), '381NH45')
        self.test('Garef', g.toRepr(), "Garef('381NH45')")
        self.test('Garef', g.toRepr(std=True), "'381NH45'")
        self.test_('Garef', repr(g), "Garef('381NH45')", "'381NH45'")  # PYGEODESY_NAMEDSTR_REPR
        self.test('Garef.precision', g.precision, 2)
        self.testCopy(g)

        self.test('Garef.latlon', fstr(g.latlon, prec=5), '57.64911, 10.40744')
        t = g.toLatLon(LL)
        self.test('Garef.toLatLon', repr(t), 'LatLon(57°38′56.8″N, 010°24′26.78″E)')
        self.testCodec3(g, '57.625, 10.375, 2.0', prec=4)
        t = Garef(t, precision=2, name='self')
        self.test('Garef(LatLon)', t, g)
        self.testCopy(t)

        for t in range(-1, 4):
            r = gars.resolution(t)
            p = gars.precision(r)
            self.test('precision', t, p, known=t < 0 or t > 2)
            b = degDMS(r, prec=0, s_D='', s_S='')  # only S_MIN
            x = ('30' + S_MIN) if p < 1 else (
                ('15' + S_MIN) if p < 2 else ('5' + S_MIN))
            self.test('resolution', b, x)  # also to test degDMS


if __name__ == '__main__':

    from pygeodesy import ellipsoidalVincenty

    t = Tests(__file__, __version__, gars)
    t.testGars(ellipsoidalVincenty.LatLon)
    t.results()
    t.exit()
