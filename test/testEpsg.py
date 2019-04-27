
# -*- coding: utf-8 -*-

# Test base classes.

__all__ = ('Tests',)
__version__ = '19.04.21'

from base import TestsBase

from pygeodesy import decodeEPSG2, encodeEPSG, epsg


class Tests(TestsBase):

    def testEpsg(self):

        for p in ('N', 'S'):
            z = epsg._UPS_ZONE
            e = encodeEPSG(z, hemipole=p)
            d = decodeEPSG2(e)
            self.test(str(z) + ' ' + p, str(d), (z, p))

            for z in range(epsg._UTM_ZONE_MIN, epsg._UTM_ZONE_MAX + 1):
                e = encodeEPSG(z, hemipole=p)
                d = decodeEPSG2(e)
                self.test(str(z) + ' ' + p, str(d), (z, p))


if __name__ == '__main__':

    t = Tests(__file__, __version__)
    t.testEpsg()
    t.results()
    t.exit()
