
# -*- coding: utf-8 -*-

# Test geohash module.

__all__ = ('Tests',)
__version__ = '17.06.23'

from base import TestsBase

from pygeodesy import geohash, Geohash, ellipsoidalVincenty


class Tests(TestsBase):

    def testGeohash(self):
        # geohash module tests
        LL = ellipsoidalVincenty.LatLon

        g = Geohash('geek')
        self.test('Geohash', repr(g), "Geohash('geek')")
        self.test('Geohash', g, 'geek')
        self.test('Geohash', Geohash(g), 'geek')
        self.test('bounds', g.bounds(LL), '(LatLon(65°23′26.25″N, 017°55′46.88″W), LatLon(65°33′59.06″N, 017°34′41.25″W))')
        self.test('toLatLon', g.toLatLon(LL), '65.478516°N, 017.753906°W')
        self.test('latlon', g.latlon, '(65.478515625, -17.75390625)')

        g = Geohash(LL(65.390625, -17.929689), precision=9)
        self.test('Geohash', g, 'geehpbpbp')
        self.test('latlon', g.latlon, '(65.390625, -17.929689)')
        self.test('toLatLon', g.toLatLon(LL), '65.390625°N, 017.929689°W')
        self.test('decode', geohash.decode(g), "('65.390646', '-17.929709')")
        self.test('decode_error', geohash.decode_error(g), '(2.1457672119140625e-05, 2.1457672119140625e-05)')
        self.test('distance1', round(g.distance1('geehpb'), 3), '2758.887')
        self.test('distance2', round(g.distance2('geehpb'), 3), '676.254')
        self.test('distance3', round(g.distance3('geehpb'), 3), '397.404')

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
        self.test('decode_error', geohash.decode_error('u120fxw'), '(0.0006866455078125, 0.0006866455078125)')
        self.test('distance1', round(geohash.distance1('u120fxw', 'u120fxws0'), 3), '486.71')
        self.test('distance2', round(geohash.distance2('u120fxw', 'u120fxws0'), 3), '3.374')
        self.test('distance3', round(geohash.distance3('u120fxw', 'u120fxws0'), 3), '2.798')

        self.test('encode', geohash.encode(69.6, -45.7), 'fur')
        self.test('decode', geohash.decode('fur'), "('69.6', '-45.7')")
        self.test('decode', geohash.decode('fu'), "('70.3', '-51')")
        self.test('decode', geohash.decode('f'), "('68', '-68')")
        self.test('decode_error', geohash.decode_error('fur'), '(0.703125, 0.703125)')
        self.test('decode_error', geohash.decode_error('fu'), '(2.8125, 5.625)')
        self.test('decode_error', geohash.decode_error('f'), '(22.5, 22.5)')


if __name__ == '__main__':

    t = Tests(__file__, __version__, geohash)
    t.testGeohash()
    t.results(nl=0)
    t.exit()
