
# -*- coding: utf-8 -*-

# Test geohash module.

__all__ = ('Tests',)
__version__ = '20.04.24'

from base import TestsBase

from pygeodesy import classname, fstr, geohash, Geohash


class Tests(TestsBase):

    def testGeohash(self, LL):
        cn = classname(LL(0, 0))

        g = Geohash('geek')
        self.test('Geohash', str(g), 'geek')
        self.test('Geohash', g.toStr(), 'geek')
        self.test('Geohash', repr(g), "'geek'")  # std_repr
        self.test('Geohash', g.toRepr(), "Geohash('geek')")
        self.test('Geohash', Geohash(g), 'geek')
        self.test('bounds', g.bounds(LL), '(%s(65°23′26.25″N, 017°55′46.88″W), %s(65°33′59.06″N, 017°34′41.25″W))' % (cn, cn))
        self.test('toLatLon', g.toLatLon(LL), '65.478516°N, 017.753906°W')
        self.test('latlon', fstr(g.latlon, prec=7), '65.4785156, -17.7539062')
        self.test('ab', fstr(g.ab, prec=7), '1.1428157, -0.3098641')
        self.testCopy(g)

        g = Geohash(LL(65.390625, -17.929689), precision=9)
        self.test('Geohash', g, 'geehpbpbp')
        self.test('toLatLon', g.toLatLon(LL), '65.390625°N, 017.929689°W')
        self.test('latlon', fstr(g.latlon, prec=8), '65.390625, -17.929689')
        self.test('ab', fstr(g.ab, prec=7), '1.1412817, -0.3129321')
        self.test('decode', geohash.decode(g), "('65.390646', '-17.929709')")
        self.test('decode_error', fstr(geohash.decode_error(g), fmt='%*e'), '2.145767e-05, 2.145767e-05')
        self.test('distance1', g.distance1('geehpb'), '2758.887', fmt='%.3f')
        self.test('distance2', g.distance2('geehpb'),  '682.760', fmt='%.3f')
        self.test('distance3', g.distance3('geehpb'),  '397.404', fmt='%.3f')
        self.test('sizes', fstr(g.sizes, prec=1), '4.8, 4.8')
        self.testCopy(g)

        for d in (g.neighbors, geohash.neighbors(g)):
            self.test('N',  d.N,  g.N)
            self.test('NE', d.NE, g.NE)
            self.test('E',  d.E,  g.E)
            self.test('SE', d.SE, g.SE)
            self.test('S',  d.S,  g.S)
            self.test('SW', d.SW, g.SW)
            self.test('W',  d.W,  g.W)
            self.test('NW', d.NW, g.NW)

            self.test('N',  d['N'],  g.N)
            self.test('NE', d['NE'], g.NE)
            self.test('E',  d['E'],  g.E)
            self.test('SE', d['SE'], g.SE)
            self.test('S',  d['S'],  g.S)
            self.test('SW', d['SW'], g.SW)
            self.test('W',  d['W'],  g.W)
            self.test('NW', d['NW'], g.NW)

        b = geohash.bounds('u120fxw')
        self.test('bounds', fstr(b, prec=8), '52.20428467, 0.11810303, 52.20565796, 0.11947632')
        d = geohash.decode('u120fxw')
        self.test('decode', fstr(d, prec=4), '52.205, 0.1188')

        for g in ('u120fxw', 'geek', 'fur', 'geehpbpbp', 'u4pruydqqvj8', 'bgr96qxvpd46', '0123456789', 'zzzzzz'):
            self.test('encode-decode', geohash.encode(*geohash.decode(g)), g)

        for p in range(8, 13):
            g = Geohash(LL(57.64911, 10.40744), precision=p)  # Jutland, Denamrk
            self.test('Geohash', g, 'u4pruydqqvj8'[:p], )
            self.test('N.E.S.W', g.N.E.S.W == g, 'True')
            self.test('E.S.W.N', g.E.S.W.N == g, 'True')
            self.test('S.W.N.E', g.S.W.N.E == g, 'True')
            self.test('W.N.E.S', g.W.N.E.S == g, 'True')
            self.test('N.E.S.S.W.W.N.N.E.S', g.N.E.S.S.W.W.N.N.E.S == g, True)  # MCCABE Law of Demeter

        self.test('encode', geohash.encode(52.205, 0.1188), 'u120fxw')
        self.test('decode', geohash.decode('u120fxw'), "('52.205', '0.1188')")
        self.test('decode_error', fstr(geohash.decode_error('u120fxw'), fmt='%*e'), '6.866455e-04, 6.866455e-04')
        self.test('distance1', geohash.distance1('u120fxw', 'u120fxws0'), '486.710', fmt='%.3f')
        self.test('distance2', geohash.distance2('u120fxw', 'u120fxws0'),   '3.374', fmt='%.3f')
        self.test('distance3', geohash.distance3('u120fxw', 'u120fxws0'),   '2.798', fmt='%.3f')
        self.test('sizes', fstr(geohash.sizes('u120fxw'), prec=1), '153.0, 153.0')

        g = Geohash('52.5009, 13.354')
        self.test('Geohash', g, 'u336xv')
        e = geohash.encode(52.5009, 13.354)
        self.test('encode', e, 'u336xv')
        self.test('equal', g == e, True)
        self.test('sizes', fstr(geohash.sizes(g), prec=1), '610.0, 1220.0')

        self.test('encode', geohash.encode(69.6, -45.7), 'fur')
        self.test('decode', geohash.decode('fur'), "('69.6', '-45.7')")
        self.test('decode', geohash.decode('fu'), "('70.3', '-51')")
        self.test('decode', geohash.decode('f'), "('68', '-68')")
        self.test('decode_error', geohash.decode_error('fur'), '(0.703125, 0.703125)')
        self.test('decode_error', geohash.decode_error('fu'), '(2.8125, 5.625)')
        self.test('decode_error', geohash.decode_error('f'), '(22.5, 22.5)')

        for t in range(0, 14):
            r = geohash.resolution2(t, t)
            p = geohash.precision(*r)
            self.test('precision', t, p, known=t < 1 or t > 12)
            b = fstr(r, prec=t + 1)
            self.test('resolution', b, b)  # just to show


if __name__ == '__main__':

    from pygeodesy import ellipsoidalVincenty

    t = Tests(__file__, __version__, geohash)
    t.testGeohash(ellipsoidalVincenty.LatLon)
    t.results()
    t.exit()
