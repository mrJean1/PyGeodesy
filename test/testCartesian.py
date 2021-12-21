
# -*- coding: utf-8 -*-

# Test module attributes.

__all__ = ('Tests',)
__version__ = '21.12.20'

from base import GeodSolve, geographiclib, isPython35, TestsBase

from pygeodesy import R_M, classname, Datums, degrees, fstr, Height, \
                modulename, Transforms  # PYCHOK expected
from pygeodesy.cartesianBase import CartesianBase
from pygeodesy.ecef import Ecef9Tuple
from pygeodesy.named import _std_NotImplemented
from pygeodesy.namedTuples import LatLon2Tuple, LatLon3Tuple, LatLon4Tuple, \
                                  PhiLam2Tuple, PhiLam3Tuple, PhiLam4Tuple, \
                                                Vector3Tuple, Vector4Tuple  # PYCHOK hanging


class Tests(TestsBase):

    def testCartesian(self, module, Sph=False, Nv=False, X=False):  # MCCABE 45

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

        if datum2:
            d = c.convertDatum(datum2)
            t = d.convertDatum(datum)
            self.test('convertDatum', t, c)  # PYCHOK attribute
            if isPython35:
                # using eval avoids SyntaxError with Python 3.4- ...
                t = eval('c @ datum2')
                self.test('__matmul__', t, d)
                if _std_NotImplemented:  # PYGEODESY_NOTIMPLEMENTED=std
                    t = eval('datum2 @ c')
                    self.test('__rmatmul__', t, d)
                # ... however t = eval("d @= ...") throws a SyntaxError
                # d @= datum2
                # self.test('__imatmul__', d, c)
                t = eval('d @ Transforms.Identity')
                self.test('__matmul__', t, d)

        self.test('height',  c.height, '-5918.380258' if Sph else '0.242887', prec=6)
        self.test('height4', c.height4().toStr(prec=1), '(3984282.2, 97.1, 4971443.2, -5918.4)' if Sph
                                                   else '(3980580.8, 97.0, 4966824.8, 0.2)')
        self.test('height4', c.height4(Cartesian=Cartesian, height=0).toStr(prec=1), '[3984282.2, 97.1, 4971443.2]' if Sph
                                                                                else '[3980580.8, 97.0, 4966824.8]')

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
            n = fstr(n, fmt='g', prec=12)
            self.test(t, n, '0.625376979018, 1.52393750974e-05, 0.780322775447, -5918.38025833' if Sph
                       else '0.622817764745, 1.51770113911e-05, 0.782366941842, 0.242886808456')

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
            self.testReturnType(c.height4(),         Vector4Tuple, 'height4')
            self.testReturnType(c.isequalTo(c),      bool,         'isequalTo')
            self.testReturnType(c.philam,            PhiLam2Tuple, 'philam')
            self.testReturnType(c.philamheight,      PhiLam3Tuple, 'philamheight')
            self.testReturnType(c.philamheightdatum, PhiLam4Tuple, 'philamheightdatum')
            self.testReturnType(c.latlonheight,      LatLon3Tuple, 'latlonheight')
            self.testReturnType(c.toEcef(),          Ecef9Tuple,   'toEcef')
            self.testReturnType(c.toLatLon(),        Ecef9Tuple if B else LatLon, 'toLatLon')
            self.testReturnType(c.toNvector(),       Vector4Tuple if B else Nvector, 'toNvector')
            self.testReturnType(c.xyz,               Vector3Tuple, 'xyz')
            c = CartesianBase(c)  # PYCHOK attribute

        if hasattr(Cartesian, 'intersections2'):
            # <https://GIS.StackExchange.com/questions/48937/calculating-intersection-of-two-circles>
            c = Cartesian(-0.00323306, -0.7915, 0.61116)
            n = classname(c, prefixed=True) + '.intersections2'
            self.test(n, c.toLatLon(height=0), '37.673442°N, 090.234036°W' if Sph
                                          else '89.998941°N, 090.234036°W')  # XXX?
            d = Cartesian(-0.0134464, -0.807775, 0.589337)
            self.test(n, d.toLatLon(height=0), '36.109987°N, 090.95367°W' if Sph
                                          else '89.99892°N, 090.95367°W')  # XXX?
            if Sph:
                x, y = c.intersections2(0.0312705, d, 0.0421788, radius=None)  # radii in radians
                self.test(n, x.toStr(prec=6), '[-0.032779, -0.784769, 0.61892]')  # -0.0327606, -0.784759, 0.618935
                self.test(n, x.toLatLon(height=0), '38.237342°N, 092.391779°W')  # 38.23838°N, 092.390487°W
                if y is not x:
                    self.test(n, y.toStr(prec=6), '[0.025768, -0.798347, 0.601646]')  # 0.0257661, -0.798332, 0.601666
                    self.test(n, y.toLatLon(height=0), '36.987868°N, 088.151309°W')  # 36.98931°N, 088.151425°W
                try:
                    from pygeodesy import trilaterate3d2  # with earth ... equivalent to Cartesian.intersections2?
                    n = modulename(trilaterate3d2, prefixed=True)
                    i, j = trilaterate3d2(c, 0.0312705, d, 0.0421788, Cartesian(0, 0, 0), 1)  # radians
                    self.test(n, i.toStr(prec=6), '[-0.032761, -0.784757, 0.618937]', known=x.minus(i).length < 5e-5)
                    self.test(n, j.toStr(prec=6),  '[0.025768, -0.798331, 0.601668]', known=y.minus(j).length < 5e-5)
                except ImportError as x:
                    self.skip(str(x), n=2)
            else:
                x, y = c.intersections2(0.0312705, d, 0.0421788, sphere=True)
                self.test(n, x.toStr(prec=6), '[-0.0035, -0.791926, 0.610589]')
                self.test(n, x.toLatLon(height=0), '89.998941°N, 090.253237°W')
                self.test(n, y.toStr(prec=6), '0.0312613')  # radius

        try:
            from pygeodesy.vector3d import intersections2
            n = modulename(intersections2, prefixed=True)
            u = Vector3Tuple(-0.00323306, -0.7915, 0.61116)
            v = Vector3Tuple(-0.0134464, -0.807775, 0.589337)
            c, r = intersections2(u, 0.0312705, v, 0.0421788, sphere=True)
            self.test(n, c.toStr(prec=6), '(-0.0035, -0.791926, 0.610589)')
            self.test(n, r.toStr(prec=6), '0.0312613', known=True)  # XXX G and g formats may add 1 decimal
            v1, v2 = intersections2(u, 0.0312705, v, 0.0421788, sphere=False)
            self.test(n, v1.toStr(prec=6), '(-0.021973, -0.766467, 0.0)')
            if v2 is not v1:
                self.test(n, v2.toStr(prec=6), '(0.027459, -0.797488, 0.0)')
        except ImportError as x:
            self.skip(str(x), n=4)

    def testReturnType(self, inst, clas, name):
        self.test(name, type(inst), clas)  # type(inst).__name__ == clas.__name__


if __name__ == '__main__':

    from pygeodesy import ellipsoidalExact, ellipsoidalNvector, ellipsoidalVincenty, \
                          sphericalNvector, sphericalTrigonometry

    t = Tests(__file__, __version__)

    t.testCartesian(sphericalNvector, Sph=True, Nv=True)
    t.testCartesian(sphericalTrigonometry, Sph=True)

    t.testCartesian(ellipsoidalNvector, Nv=True)
    t.testCartesian(ellipsoidalVincenty)

    if geographiclib:
        from pygeodesy import ellipsoidalKarney
        t.testCartesian(ellipsoidalKarney)

    if GeodSolve:
        from pygeodesy import ellipsoidalGeodSolve
        t.testCartesian(ellipsoidalGeodSolve)

    t.testCartesian(ellipsoidalExact, X=True)

    t.results()
    t.exit()
