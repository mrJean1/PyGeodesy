
# -*- coding: utf-8 -*-

# Test geohash module.

__all__ = ('Tests',)
__version__ = '17.04.22'

from tests import Tests as _Tests

from pygeodesy import geohash, Geohash, ellipsoidalVincenty


class Tests(_Tests):

    def testGeohash(self):
        # geohash module tests
        _LL = ellipsoidalVincenty.LatLon

        g = Geohash('geek')
        self.test('Geohash', repr(g), "Geohash('geek')")
        self.test('Geohash', str(g), 'geek')
        self.test('Geohash', Geohash(g), 'geek')
        self.test('bounds', g.bounds(_LL), '(LatLon(65°23′26.25″N, 017°55′46.88″W), LatLon(65°33′59.06″N, 017°34′41.25″W))')
        self.test('toLatLon', str(g.toLatLon(_LL)), '65.478516°N, 017.753906°W')
        self.test('latlon', str(g.latlon), '(65.478515625, -17.75390625)')

        for p in range(9, 13):
            g = Geohash(_LL(57.64911, 10.40744), precision=p)  # Jutland, Denamrk
            self.test('Geohash', g, 'u4pruydqqvj8'[:p], )
            self.test('N.E.S.W', g.N.E.S.W == g, 'True')
            self.test('E.S.W.N', g.E.S.W.N == g, 'True')
            self.test('S.W.N.E', g.S.W.N.E == g, 'True')
            self.test('W.N.E.S', g.W.N.E.S == g, 'True')
            self.test('N.E.S.S.W.W.N.N.E.S', g.N.E.S.S.W.W.N.N.E.S == g, 'True')  # PYCHOK expected

        self.test('encode', geohash.encode(52.205, 0.1188), 'u120fxw')
        self.test('decode', geohash.decode('u120fxw'), "('52.205', '0.1188')")
        self.test('decode_error', geohash.decode_error('u120fxw'), '(0.0006866455078125, 0.0006866455078125)')

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
    t.results()
    t.exit()
