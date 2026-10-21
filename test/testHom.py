
# -*- coding: utf-8 -*-

# Some L{hom} C{HotineOblique} Mercator tests.

__all__ = ('Tests',)
__version__ = '26.10.21'

from bases import isPython3, TestsBase, startswith

from pygeodesy import Ellipsoids as Es, fstr, hom, grades2degrees, \
                      HotineOblique1, HotineOblique2, HotineOblique3, \
                      HotineObliqueA, HotineObliqueB, latDMS, lonDMS


from math import fabs


def _err(f, x, y):
    return max(fabs(f.x - x), fabs(f.y - y))


def _str(t):
    return str(t).split('HotineOblique')[0] + '...'


class Tests(TestsBase):

    def testHotine123(self):

        h = HotineOblique1('47 30N, 122 18W', '25 42N, 80 12W', k0=0.9996)  # Snyder page 272
        self.test('1hotine', h, '(47.5, -122.3), (25.7, -80.2), k0=0.9996, radius=', known=startswith)
        self.test('1latlon1', h.latlon1, '(47.5, -122.3)')
        self.test('1latlon2', h.latlon2, '(25.7, -80.2)')
        f = h.forward(38.8092128, -76.8707953)
        self.test('1forward', f, '(4655575.277446, -1248887.591514, 4655575.277446, -1248887.591514, 38.809213, -76.870795, 1.018883, ', known=startswith)
        r = h.reverse(f.x, f.y)
        self.test('1reverse', r, '(4655575.277446, -1248887.591514, 4655575.277446, -1248887.591514, 38.809213, -76.870795, 0.9996, ', known=startswith)
        self.test('1roundtrip', fstr(r[:4], prec=6), fstr(f[:4], prec=6))

        h = HotineOblique1((45, 0), (0, -90), radius=1, k0=1)  # page 272
        self.test('1hotine', h, '(45.0, 0.0), (0.0, -90.0), radius=1.0', nl=1)
        f = h.forward(-30, 120)
        self.test('1forward', f.toStr(prec=9), '(-2.420133502, -0.047402646, -2.420133502, -0.047402646, -30.0, 120.0, 1.001123716, ', known=startswith)
        r = h.reverse(f.x, f.y)
        self.test('1reverse', r.toStr(prec=9), '(-2.420133502, -0.047402646, -2.420133502, -0.047402646, -30.0, 120.0, 1.0, ', known=startswith)
        self.test('1roundtrip', fstr(r[:4], prec=6), fstr(f[:4], prec=6))
        self.test('1poleP', h.latlonP, '(45.0, 180.0)')
        self.test('1pole0', h.latlon0, '(0.0, -90.0)')

        h = HotineOblique2(20, -68.6557771, azimuth=48.8062990, radius=1, k0=1)  # page 273
        self.test('2hotine', h, '20.0, -68.6558, 48.8063, radius=1.0', nl=1)
        self.test('1latlonC', h.latlonC, '(20.0, -68.655777)')
        self.test('2latlonP', h.latlonP, '(45.0, 180.0)')

        h = HotineOblique1((45, 0), (0, -90), radius=1)  # page 273
        self.test('1hotine', h, '(45.0, 0.0), (0.0, -90.0), radius=1.0', nl=1)
        r = h.reverse(-2.4201335, -0.0474026)
        self.test('1reverse', r.toStr(prec=3), '(-2.42, -0.047, -2.42, -0.047, -30.0, 120.0, 1.0, ', known=startswith)
        f = h.forward(r.lat, r.lon)
        self.test('1forward', f.toStr(prec=3), '(-2.42, -0.047, -2.42, -0.047, -30.0, 120.0, 1.001, ', known=startswith)

        h = HotineOblique3(20.0, -68.6557, radius=1)
        self.test('3hotine', h, '20.0, -68.6557, radius=1.0')

    def testHotineAB(self):

        h = HotineObliqueA('47 30N, 122 18W', '25 42N, 80 12W', earth=Es.Clarke1866, k0=0.9996, lat0=40, x0=4e6, y0=5e5)  # Snyder page 274, "Alternate A"
        self.test('Ahotine1', h, '(47.5, -122.3), (25.7, -80.2), ', known=startswith, nl=1)
        self.test('AmapAB', h.ABmap.toRepr(), 'HotineObliqueA(A=6379333.23054, B=1.001173, D=1.304333, D2=1.701284, E=1.002186, G=1.088777, uC=0.0, zC=-0.993906, a0=0.698132, b0=-0.984976, g0=-0.697885, t0=0.468343, x0=4000000.0, y0=500000.0)')
        self.test('Alatlon0', h.latlon0, '(40.0, -56.434963)')
        self.test('Alatlon1', h.latlon1, '(47.5, -122.3)')
        self.test('Alatlon2', h.latlon2, '(25.7, -80.2)')
        self.test('AazimuthC', h.azimuthC, '-56.946607', prec=6)
        self.test('Agamma0', h.gamma0, '-39.985883', prec=6)
        f = h.forward('40 48N', -74)
        self.test('Aforward', _str(f), '(963436.092174, 4369142.81027, 4655443.688172, 1586767.317166, 40.8, -74.0, 1.030755, ...', known=startswith)
        r = h.reverse(f.x, f.y)
        x = '(963436.092174, 4369142.81027, 4655443.688172, 1586767.317166, 40.8, -74.0, 0.9996, ...'
        self.test('Areverse', _str(r), x, known=startswith)
        self.test('Aiteration', r.iteration, None)
        r = h.reverse(f.x, f.y, triPhi=9)
        self.test('Areverse', _str(r), x, known=startswith)
        self.test('Aiteration', r.iteration, r.iteration)

        h = HotineObliqueA((47.5, -122.3), (25.7, -80.2), earth=Es.WGS84, k0=0.9996, lat0=40, x0=4e6, y0=5e5)  # https://Neacsu.net/geodesy/snyder/numerics/ome/#sphere
        self.test('Ahotine2', h, '(47.5, -122.3), (25.7, -80.2), ', known=startswith, nl=1)
        self.test('AmapAB', h.ABmap.toRepr(), 'HotineObliqueA(A=6379223.515405, B=1.00116, D=1.304344, D2=1.701315, E=1.002162, G=1.088805, uC=0.0, zC=-0.993887, a0=0.698132, b0=-0.984984, g0=-0.697868, t0=0.46832, x0=4000000.0, y0=500000.0)')
        self.test('Alatlon0', h.latlon0, '(40.0, -56.435422)')
        self.test('Alatlon1', h.latlon1, '(47.5, -122.3)')
        self.test('Alatlon2', h.latlon2, '(25.7, -80.2)')
        self.test('AazimuthC', h.azimuthC, '-56.945553', prec=6)
        self.test('Agamma0', h.gamma0, '-39.984872', prec=6)
        f = h.forward(40.8, -74)
        self.test('Aforward', _str(f), '(963497.880351, 4369215.405042, 4655460.671537, 1586776.263063, 40.8, -74.0, 1.030756, ...', known=startswith)
        r = h.reverse(f.x, f.y)
        x = '(963497.880351, 4369215.405042, 4655460.671537, 1586776.263063, 40.8, -74.0, 0.9996, ...'
        self.test('Areverse', _str(r), x, known=startswith)
        self.test('Aiteration', r.iteration, None)
        r = h.reverse(f.x, f.y, triPhi=9)
        self.test('Areverse', _str(r), x, known=startswith)
        self.test('Aiteration', r.iteration, r.iteration)

        h = HotineObliqueB(36, -77.7610558, 14.3394883, earth=Es.Clarke1866)  # Snyder page 276, "Alternate B"
        self.test('Bhotine3', h, '36.0, -77.7611, 14.3395, ', known=startswith, nl=1)
        self.test('BmapAB', h.ABmap.toRepr(), 'HotineObliqueB(A=6380777.046609, B=1.001459, D=1.235119, D2=1.52552, E=1.001698, G=0.724928, uC=4092868.936638, zC=0.250271, a0=0.628319, b0=-1.505896, g0=0.201889, t0=0.511558, x0=0.0, y0=0.0)')
        self.test('Blatlon0', h.latlon0, '(0.0, -86.28148)')
        self.test('BlatlonC', h.latlonC, '(36.0, -77.761056)')
        self.test('BazimuthC', h.azimuthC, '14.33949', prec=5)
        self.test('Bgamma0', h.gamma0, '11.56740', prec=5)
        f = h.forward(38.8092128, -76.8797953)
        x = '(4414439.0, -2356.3, ...'
        self.test('Bforward', _str(f), x, known=True)  # x, y?
        r = h.reverse(f.x, f.y)
        self.test('Breverse', _str(r), x, known=True)
        self.test('Biteration', r.iteration, None)
        r = h.reverse(f.x, f.y, triPhi=9)
        self.test('Breverse', _str(r), x, known=True)
        self.test('Biteration', r.iteration, r.iteration)

        h = HotineObliqueA('47 30N, 122 18W', '25 42N, 80 12W', earth=Es.Sphere, lat0=40)
        self.test('Ahotine4', h, "(47.5, -122.3), (25.7, -80.2), earth=Ellipsoid(name='Sphere'", known=startswith, nl=1)
        self.test('AmapAB', h.ABmap.toRepr(), 'HotineObliqueA(A=6371008.771415, B=1.0, D=1.305407, D2=1.704088, E=1.0, G=1.091309, uC=0.0, zC=-0.992236, a0=0.698132, b0=-0.985706, g0=-0.696284, t0=0.466308, x0=0.0, y0=0.0)')

        h = HotineObliqueB(4, 115, '53 18 56.9537', earth=Es.Everest1967, k0=0.99984, x0=590476.87, y0=442857.65, name='Timbalai 1948 / R.S.O. Borneo (m)')  # IOGP page 63
        self.test('Bhotine5', h, "4.0, 115.0, 53.3158, earth=Ellipsoid(name='Everest1967', ", known=startswith, nl=1)
        self.test('Bxy0', h.xy0, '(590476.87, 442857.65)')
        self.test('BmapAB', h.ABmap.toRepr(), 'HotineObliqueB(A=6376278.686012, B=1.003303, D=1.002426, D2=1.004857, E=1.000003, G=0.069695, uC=738096.094116, zC=0.930537, a0=0.069813, b0=1.914373, g0=0.927295, t0=0.932947, x0=590476.87, y0=442857.65)')
        self.test('BlatlonC', h.latlonC, '(4.0, 115.0)')
        self.test('BazimuthC', h.azimuthC, '53.316', prec=3)  # 53 18 57...
        self.test('Bgamma0', h.gamma0, '53.130', prec=3)  # 53 07 48...
        h.xyuvIOGP = True
        self.test('BxyuvIOGP', h.xyuvIOGP, True)
        f = h.forward('5 23 14.1129N', '115 48 19.8196E')
        self.test('Bforward1', _str(f), '(679245.728142, 596562.777494, 163238.163018, -69702.78709, 5.387254, 115.805505, 1.009708, ...', known=startswith)  # IOGP E 679245.73, N 596562.78, ...
        r = h.reverse(f.x, f.y)  # triPhi=0
        x = '(679245.728142, 596562.777494, 901334.257135, -69702.78709, 5.387254, 115.805505, 0.99984, ...'
        self.test('Breverse1', _str(r), x, known=startswith)
        self.test('Biteration1', r.iteration, None)
        r = h.reverse(f.x, f.y, triPhi=9)
        self.test('Breverse1', _str(r), x, known=startswith)
        self.test('Biteration1', r.iteration, r.iteration)  # 5 macOS

        f = h.forward('4 39 20.783N', '114 28 10.539E')  # GeoTIFF.MapTools.org/proj_list/hotine_oblique_mercator.html
        self.test('Bforward2', _str(f), '(531404.80075, 515187.848228, -3859.536452, -93307.400133, 4.655773, 114.469594, 1.006658, ...', known=startswith)  # IOGP E 679245.73, N 596562.78, ...
        r = h.reverse(f.x, f.y)  # triPhi=0
        x = '(531404.80075, 515187.848228, 734236.557664, -93307.400133, 4.655773, 114.469594, 0.99984, ...'
        self.test('Breverse2', _str(r), x, known=startswith)
        self.test('Biteration2', r.iteration, None)
        r = h.reverse(f.x, f.y, triPhi=9)
        self.test('Breverse2', _str(r), x, known=startswith)
        self.test('Biteration2', r.iteration, r.iteration)  # 5 macOS

        h = HotineObliqueB(grades2degrees(-21), grades2degrees(51.5969213), grades2degrees(21), earth=Es.Intl1924, k0=0.9995, x0=4e5, y0=8e5, name='Tananarive (Paris)')  # IOGP page 65, Laborde
        self.test('Bhotine6', h.toRepr(), "HotineObliqueB(-18.9, 46.4372, 18.9, earth=Ellipsoid(name='Intl1924', a=6378388, f=0.003367, ", known=startswith, nl=1)
        self.test('BmapAB', h.ABmap.toRepr(), 'HotineObliqueB(A=6375433.454701, B=1.002708, D=1.056653, D2=1.116516, E=0.999703, G=-0.341344, uC=-2201599.156405, zC=0.329867, a0=-0.329867, b0=0.92034, g0=0.311567, t0=1.396319, x0=400000.0, y0=800000.0)')
        self.test('Blatlon0', h.latlon0, '(0.0, 52.731578)')
        f = h.forward('16 11 23.28S', '44 27 27.260E')
        x = '(188333.848, 1098841.091, ...'  # IOGP E 188333.848, N 1098841.091, ... Laborde!!
        self.test('BforwardL', _str(f), x, known=True)
        r = h.reverse(f.x, f.y)  # triPhi=0)
        self.test('BreverseL', _str(r), x, known=True)
        if isPython3:
            self.test('BreverseL', (latDMS(r.lat), lonDMS(r.lon)), "('16°11′23.28″S', '044°27′27.26″E')")

        f = h.forward('18 54S', '47 30E')  # IOGP page 66, Comparing the Hotine Oblique Mercator (variant B) ...
        self.test('BforwardL', _str(f), '(511921.054, 799665.520, ...', error=_err(f, 511921.054, 799665.520), known=True, nl=1)  # IOGP N/X? 799665.520, E/Y? 511921.054, ... XXX
        f = h.forward('16 12S', '44 24E')
        self.test('BforwardL', _str(f), '(182184.985, 1097651.426, ...', error=_err(f, 182184.985, 1097651.426), known=True)  # IOGP N/X? 1097651.426, E/Y? 182184.985, ... XXX
        f = h.forward('25 40S', '45 18E')
        self.test('BforwardL', _str(f), '(285294.788, 50636.850, ...', error=_err(f, 285294.788, 50636.850), known=True)  # IOGP N/X? 50636.850, E/Y? 285294.788, ... XXX
        f = h.forward('12 00S', '49 12E')
        self.test('BforwardL', _str(f), '(701352.935, 1561109.350, ...', error=_err(f, 701352.935, 1561109.350), known=True)  # IOGP N/X? 1561109.350, E/Y? 701352.935, ... XXX


if __name__ == '__main__':

    t = Tests(__file__, __version__, hom)
    t.testHotine123()
    t.testHotineAB()
    t.results()
    t.exit()
