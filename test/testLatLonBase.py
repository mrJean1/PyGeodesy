
# -*- coding: utf-8 -*-

# Test LatLon base classes.

__all__ = ('Tests',)
__version__ = '21.07.01'

from base import GeodSolve, TestsBase

from pygeodesy import clips, F_D, F_DMS, precision


class Tests(TestsBase):

    def testBases(self, module, LatLon, Base=False, Sph=False):

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

        if not Base:  # coverage
            t = clips(p.rhumbLine(30).toStr(), limit=150)
            self.test('rhumbLine', t, t, nl=1)
            q = LatLon(50.964, 1.853)
            t = clips(p.rhumbLine(q).toStr(), limit=150)
            self.test('rhumbLine', t, t)

            # same rhumb test as testSpherical, invoking rhumbx.Rhumb[Line]
            # with exact=False for all tests and the differences between
            # the ellipsoidal and spherical cases are due to the ellipsoid.
            p = LatLon(51.127, 1.338)
            z = p.rhumbAzimuthTo(q, exact=Sph)
            self.test('rhumbAzimuthTo', z, 107.563 if Sph else 116.661, fmt='%.3f')

            d = p.rhumbDestination(40300, 116.7, exact=Sph)
            self.test('rhumbDestination', d, '50.964155°N, 001.853°E' if Sph
                                        else '50.964234°N, 001.851383°E')
            self.test('rhumbDestination', isinstance(d, LatLon), True)

            d = p.rhumbDistanceTo(q, exact=Sph)  # force rhumbx.Rhumb[Line]
            self.test('rhumbDistanceTo', d, 42186.1 if Sph else 40413.1, fmt='%.1f')

            m = p.rhumbMidpointTo(q, exact=Sph)  # ditto
            self.test('rhumbMidpointo-0.5', m, '51.069759°N, 001.625988°E' if Sph
                                          else '51.045501°N, 001.595726°E')
            self.test('rhumbMidpointo', isinstance(m, LatLon), True)
            m = p.rhumbMidpointTo(q, exact=Sph, fraction=0)  # ditto
            self.test('rhumbMidpointo-0.0', m, '51.127°N, 001.338°E')
            m = p.rhumbMidpointTo(q, exact=Sph, fraction=0.25)  # ditto
            self.test('rhumbMidpointo-0.25', m, '51.09838°N, 001.482038°E' if Sph
                                           else '51.08625°N, 001.46692°E')
            m = p.rhumbMidpointTo(q, exact=Sph, fraction=0.75)  # ditto
            self.test('rhumbMidpointo-0.75', m, '51.041139°N, 001.769848°E' if Sph
                                           else '51.00475°N, 001.724419°E')
            m = p.rhumbMidpointTo(q, exact=Sph, fraction=1)  # ditto
            self.test('rhumbMidpointo-1.0', m, '51.012519°N, 001.913619°E' if Sph
                                          else '50.964°N, 001.853°E')
            m = p.rhumbMidpointTo(q, exact=Sph, fraction=2)  # ditto
            self.test('rhumbMidpointo-2.0', m, '50.898038°N, 002.48782°E' if Sph
                                          else '50.800995°N, 002.366201°E')


if __name__ == '__main__':

    from pygeodesy import ellipsoidalExact, ellipsoidalKarney, \
                          ellipsoidalNvector, ellipsoidalVincenty, \
                          sphericalNvector, sphericalTrigonometry, \
                          LazyImportError

    t = Tests(__file__, __version__)

    t.testBases(sphericalNvector, sphericalNvector.LatLon, Sph=True)
    t.testBases(sphericalTrigonometry, sphericalTrigonometry.LatLon, Sph=True)

    t.testBases(ellipsoidalNvector, ellipsoidalNvector.LatLon)
    t.testBases(ellipsoidalVincenty, ellipsoidalVincenty.LatLon)
    t.testBases(ellipsoidalKarney, ellipsoidalKarney.LatLon)
    t.testBases(ellipsoidalExact, ellipsoidalExact.LatLon)

    if GeodSolve:
        from pygeodesy import ellipsoidalGeodSolve
        t.testBases(ellipsoidalGeodSolve, ellipsoidalGeodSolve.LatLon)

    try:  # (INTERNAL) modules not explicitly exported
        from pygeodesy import ellipsoidalBase, ellipsoidalBaseDI, \
                              latlonBase, nvectorBase, sphericalBase

        t.testBases(ellipsoidalBase, ellipsoidalBase.LatLonEllipsoidalBase)
        t.testBases(ellipsoidalBaseDI, ellipsoidalBaseDI.LatLonEllipsoidalBaseDI)
        t.testBases(latlonBase, latlonBase.LatLonBase, Base=True)
        t.testBases(nvectorBase, nvectorBase.LatLonNvectorBase, Base=True)
        t.testBases(sphericalBase, sphericalBase.LatLonSphericalBase, Sph=True)
    except LazyImportError:
        pass

    t.results()
    t.exit()
