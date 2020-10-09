
# -*- coding: utf-8 -*-

# Test LatLon base classes.

__all__ = ('Tests',)
__version__ = '20.10.08'

from base import TestsBase

from pygeodesy import F_D, F_DMS, precision


class Tests(TestsBase):

    def testBases(self, module, LatLon):

        self.subtitle(module, LatLon.__name__)

        p = LatLon(50.06632, -5.71475)
        self.test('lat, lon', p, '50.06632°N, 005.71475°W')
        q = LatLon('50°03′59″N', """005°42'53"W""")
        self.test('lat, lon', q, '50.066389°N, 005.714722°W')

        p = LatLon(52.205, 0.119)
        q = p.copy()  # LatLon(52.205, 0.119)
        self.test('isequalTo',  q.isequalTo(p), True)
        self.test('isequalTo3', q.isequalTo3(p), True)

        self.test('latlon',       q.latlon, p.latlon)
        self.test('latlonheight', q.latlonheight, p.latlonheight)
        self.test('phimlam',       q.philam, p.philam)
        self.test('phimlamheight', q.philamheight, p.philamheight)

        t = precision(F_DMS, 0)
        p = LatLon(51.4778, -0.0016)
        self.test('toStr', p.toStr(), '''51°28'40"N, 000°00'06"W''')
        self.test('toStr', p.toStr(F_D), '51.4778°N, 000.0016°W')

        self.test('precision', precision(F_DMS), 0)
        p = LatLon(51.4778, -0.0016, 42)
        self.test('toStr', p.toStr(), '''51°28'40"N, 000°00'06"W, +42.00m''')
        precision(F_DMS, t)  # restore

        q.latlon = p.latlon
        self.test('isequalTo',  q.isequalTo(p), True)
        self.test('isequalTo3', q.isequalTo3(p), False)

        self.test('latlon', q.latlon, p.latlon)
        self.test('phimlam', q.philam, p.philam)

        q.latlon = p.latlon + (p.height,)
        self.test('isequalTo', q.isequalTo(p), True)
        self.test('isequalTo3', q.isequalTo3(p), True)

        self.test('latlon',       q.latlon, p.latlon)
        self.test('latlonheight', q.latlonheight, p.latlonheight)
        self.test('phimlam',       q.philam, p.philam)
        self.test('phimlamheight', q.philamheight, p.philamheight)

        self.test('latlon',       repr(q.latlon), repr(p.latlon))
        self.test('latlonheight', repr(q.latlonheight), repr(p.latlonheight))
        self.test('phimlam',       repr(q.philam), repr(p.philam))
        self.test('phimlamheight', repr(q.philamheight), repr(p.philamheight))


if __name__ == '__main__':

    from pygeodesy import ellipsoidalKarney, ellipsoidalNvector, \
                          ellipsoidalVincenty, sphericalNvector, \
                          sphericalTrigonometry, LazyImportError

    t = Tests(__file__, __version__)

    t.testBases(ellipsoidalKarney, ellipsoidalKarney.LatLon)
    t.testBases(ellipsoidalNvector, ellipsoidalNvector.LatLon)
    t.testBases(ellipsoidalVincenty, ellipsoidalVincenty.LatLon)

    t.testBases(sphericalNvector, sphericalNvector.LatLon)
    t.testBases(sphericalTrigonometry, sphericalTrigonometry.LatLon)

    try:  # (INTERNAL) modules not explicitly exported
        from pygeodesy import ellipsoidalBase, latlonBase, \
                              nvectorBase, sphericalBase

        t.testBases(latlonBase, latlonBase.LatLonBase)
        t.testBases(nvectorBase, nvectorBase.LatLonNvectorBase)
        t.testBases(ellipsoidalBase, ellipsoidalBase.LatLonEllipsoidalBase)
        t.testBases(sphericalBase, sphericalBase.LatLonSphericalBase)
    except LazyImportError:
        pass

    t.results()
    t.exit()
