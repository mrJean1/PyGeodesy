
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '19.05.09'

from base import TestsBase

from pygeodesy import decodeEPSG2, encodeEPSG, epsg, toUtmUps8


class Tests(TestsBase):

    def testEpsg(self, decode, encode):

        for p in ('N', 'S'):
            z = epsg._UPS_ZONE
            e = encode(z, hemipole=p)
            d = decode(e)
            self.test(str(z) + ' ' + p, str(d), (z, p))
            for z in range(epsg._UTM_ZONE_MIN, epsg._UTM_ZONE_MAX + 1):
                e = encode(z, hemipole=p)
                d = decode(e)
                self.test(str(z) + ' ' + p, str(d), (z, p))

    def testEpsgTMcoord(self, n, lat, lon):
        u = toUtmUps8(lat, lon)
        e = u.toEpsg()
        s = ' '.join(u.toStr(B=True).split()[:2])
        self.test('TMcoord ' + str(n), e.utmupsStr(B=True), s)


if __name__ == '__main__':

    from testTMcoords import _TMcoords

    t = Tests(__file__, __version__)
    t.testEpsg(epsg.decode2, epsg.encode)
    t.testEpsg(decodeEPSG2, encodeEPSG)  # DEPRECATED

    for n, coord in enumerate(_TMcoords.readlines()):
        t.testEpsgTMcoord(n + 1, *map(float, coord.split()[:2]))

    t.results()
    t.exit()
