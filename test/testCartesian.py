
# -*- coding: utf-8 -*-

# Test module attributes.

__all__ = ('Tests',)
__version__ = '20.03.07'

from base import geographiclib, TestsBase

from pygeodesy import Datums, degrees  # PYCHOK expected
from pygeodesy.cartesianBase import CartesianBase
from pygeodesy.ecef import Ecef9Tuple
from pygeodesy.named import LatLon2Tuple, LatLon3Tuple, LatLon4Tuple, \
                            PhiLam2Tuple, PhiLam3Tuple, PhiLam4Tuple, \
                                          Vector3Tuple, Vector4Tuple  # PYCHOK hanging


class Tests(TestsBase):

    def testCartesian(self, module, Sph=False, Nv=True):  # MCCABE 45

        self.subtitle(module, 'Cartesian')

        Cartesian = module.Cartesian
        LatLon    = module.LatLon
        Nvector   = module.Nvector if Nv else Vector4Tuple

        datum  = Datums.Sphere if Sph else Datums.WGS84
        datum2 = None          if Sph else Datums.WGS72
        # <https://www.Movable-Type.co.UK/scripts/geodesy/docs/
        #        latlon-nvector-ellipsoidal.js.html#line309>
        c = Cartesian(3980581, 97, 4966825, datum=datum)
        self.test('Cartesian0', c.toStr(prec=0), '[3980581, 97, 4966825]')
        self.test('Cartesian4', c.toStr(prec=4), '[3980581.0, 97.0, 4966825.0]')

        self.test('isEllipsoidal', c.isEllipsoidal, not Sph)
        self.test('isSpherical',   c.isSpherical,       Sph)
        self.testCopy(c)

        n = c.toNvector()  # (x=0.622818, y=0.00002, z=0.782367, h=0.242887)
        t = n.classname  # Nvector.__name__
        if Nv:
            self.test(t, repr(n), 'Nvector(0.62538, 0.00002, 0.78032, -5918.38)' if Sph
                             else 'Nvector(0.62282, 0.00002, 0.78237, +0.24)')
            self.test(t+'3', n.toStr(prec=3), '(0.625, 0.0, 0.78, -5918.38)' if Sph
                                         else '(0.623, 0.0, 0.782, +0.24)')
            self.test(t+'6', n.toStr(prec=6), '(0.625377, 0.000015, 0.780323, -5918.38)' if Sph
                                         else '(0.622818, 0.000015, 0.782367, +0.24)')  # PYCHOK attribute
        else:
            self.test(t, repr(n), '(x=0.6253769790183048, y=1.5239375097448227e-05, z=0.7803227754472505, h=-5918.3802583276365)' if Sph
                             else '(x=0.6228177647454303, y=1.517701139112776e-05, z=0.782366941841975, h=0.24288680875513333)', known=True)

        for ll in ((50.0379, 8.5622),  # FRA
                   (51.47,   0.4543),  # LHR
                   # <https://www.EdWilliams.org/avform.htm#XTE>
                   (degrees(0.709186), -degrees(1.287762)),  # JFK
                   (33.+57./60, -(118.+24./60)),  # LAX
                   # <https://GeographicLib.SourceForge.io/html/python/examples.html>
                   (-41.32, 174.81),  # WNZ, Wellington, NZ
                   (40.96,    5.50),  # SAL, Salamanca, Spain
                   (40.1,   116.6),   # BJS, Beijing Airport
                   (37.6,  -122.4)):  # SFO
            if datum2:
                t = c.convertDatum(datum2).convertDatum(datum)
                self.test('convertDatum', t, c)  # PYCHOK attribute
            p = LatLon(*ll)
            q = p.toCartesian().toLatLon()
            t = str(q)
            self.test('LatLon', t, p, known=t.endswith('m'))  # PYCHOK attribute

        # c = Cartesian(3980581, 97, 4966825, datum=datum)
        t = c.copy()
        self.test('copy', t.isequalTo(c), True)
        self.test('__eq__', t == c, True)
        self.test('__ne__', t != c, False)

        if hasattr(Cartesian, 'convertRefFrame'):
            pass  # PYCHOK attribute

        for B in (False, True):  # check return types
            t = c.__class__
            self.test('Cartesian', t, t)
            # self.testReturnType(c.Ecef,             Ecef,       c.Ecef.__name__)
            self.testReturnType(c.latlon,            LatLon2Tuple, 'latlon')
            self.testReturnType(c.latlonheight,      LatLon3Tuple, 'latlonheight')
            self.testReturnType(c.latlonheightdatum, LatLon4Tuple, 'latlonheightdatum')
            self.testReturnType(c.isequalTo(c),      bool,         'isequalTo')
            self.testReturnType(c.philam,            PhiLam2Tuple, 'philam')
            self.testReturnType(c.philamheight,      PhiLam3Tuple, 'philamheight')
            self.testReturnType(c.philamheightdatum, PhiLam4Tuple, 'philamheightdatum')
            self.testReturnType(c.to3llh(),          LatLon4Tuple, 'to3llh')
            self.testReturnType(c.toEcef(),          Ecef9Tuple,   'toEcef')
            self.testReturnType(c.toLatLon(),        Ecef9Tuple if B else LatLon, 'toLatLon')
            self.testReturnType(c.toNvector(),       Vector4Tuple if B else Nvector, 'toNvector')
            self.testReturnType(c.xyz,               Vector3Tuple, 'xyz')
            c = CartesianBase(c)  # PYCHOK attribute

    def testReturnType(self, inst, clas, name):
        self.test(name, type(inst), clas)  # type(inst).__name__ == clas.__name__


if __name__ == '__main__':

    from pygeodesy import ellipsoidalNvector, ellipsoidalVincenty, \
                          sphericalNvector, sphericalTrigonometry

    t = Tests(__file__, __version__)

    if geographiclib:
        from pygeodesy import ellipsoidalKarney
        t.testCartesian(ellipsoidalKarney, Sph=False, Nv=False)

    t.testCartesian(ellipsoidalNvector, Sph=False, Nv=True)
    t.testCartesian(ellipsoidalVincenty, Sph=False, Nv=False)

    t.testCartesian(sphericalNvector, Sph=True, Nv=True)
    t.testCartesian(sphericalTrigonometry, Sph=True, Nv=False)

    t.results()
    t.exit()
