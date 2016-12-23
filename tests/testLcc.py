
# -*- coding: utf-8 -*-

# Test LCC functions and methods.

__version__ = '16.12.23'

if __name__ == '__main__':

    from tests import Tests as _Tests

    from geodesy import Conic, Conics, Datums, F_D, F_DMS, Lcc, toLcc
    from geodesy.ellipsoidalNvector  import LatLon as _LatLon
    from geodesy.ellipsoidalVincenty import LatLon

    # Snyder, pp 297 <https://pubs.er.usgs.gov/djvu/PP/PP_1395.pdf>
    Snyder = Conic(LatLon(23, -96, datum=Datums.NAD27),
                   33, 45, E0=0, N0=0, name='Snyder')

    class Tests(_Tests):

        def testConic(self, LatLon, n=''):

            n = 'Snyder' + str(n)
            c = Conic(LatLon(23, -96, datum=Datums.NAD27), 33, 45, E0=0, N0=0, name=n)
            self.test(n, c, "lat0=23.0, lon0=-96.0, par1=33.0, par2=45.0, E0=0, N0=0, k0=1, SP=2, datum=(ellipsoid=Ellipsoids.Clarke1866, transform=Transforms.NAD27, name='NAD27'), name='%s'" % (n,))

        def testLcc(self):

            lb = Lcc(448251, 5411932.0001)
            self.test('lb1', lb.toStr(4), '448251.0 5411932.0001')
            self.test('lb1', lb.toStr(sep=', '), '448251, 5411932')
            self.test('lb1', lb.conic.name2, 'WRF_Lb.WGS84')

            ll = LatLon(46.5, 3)
            self.test('LatLon', ll, '46.5°N, 003.0°E')
            self.test('LatLon', ll.toStr(form=F_DMS), '46°30′00.0″N, 003°00′00.0″E')
            lb = toLcc(ll, conic=Conics.Fr93Lb)
            self.test('toLcc1', str(lb), '700000 6600000')
            self.test('toLcc1', lb.toLatLon(LatLon), '46.5°N, 003.0°E')

            lb = Lcc(1894410.9, 1564649.5, conic=Snyder)
            self.test('lb2', lb, '1894411 1564650')
            self.test('lb2', lb.conic.datum.ellipsoid.name, 'Clarke1866')
            ll = lb.toLatLon(LatLon)  # Clark1866
            self.test('toLatLon2', ll.toStr(prec=6, form=F_D), '35.0°N, 075.0°W')
            self.test('toLatLon2', ll.toStr(prec=4, form=F_DMS), '35°00′00.0007″N, 074°59′59.9997″W')
            self.test('toLatLon2', ll.datum.name, 'NAD27')
            lb = toLcc(ll, conic=Snyder)
            self.test('toLcc2', lb.toStr(prec=1), '1894410.9 1564649.5')
            self.test('toLcc2', lb.conic.name2, 'Snyder.NAD27')

            for n, c in sorted(Conics.items()):
                d = abs(c.par1 - c.par2)
                if d > 0:  # test corners of the conic
                    for ll in (LatLon(c.par1, c.lon0 - d, datum=c.datum),
                               LatLon(c.par1, c.lon0,     datum=c.datum),
                               LatLon(c.par1, c.lon0 + d, datum=c.datum),
                               LatLon(c.par2, c.lon0 - d, datum=c.datum),
                               LatLon(c.par2, c.lon0,     datum=c.datum),
                               LatLon(c.par2, c.lon0 + d, datum=c.datum)):
#                       self.test(n, ll, str(ll))  # PYCHOK expected
                        lb = toLcc(ll, conic=c)
#                       self.test(n, lb, '')
                        ll_ = lb.toLatLon(LatLon)
                        self.test(n, ll, str(ll_))
                        self.test(n, ll.datum.name, ll_.datum.name)

    t = Tests(__file__, __version__)
    t.testLcc()
    t.testConic( LatLon, 1)
    t.testConic(_LatLon, 2)
    t.results()
    t.exit()

    # Typical test results (on MacOS X):

    # testing testLcc.py version 16.12.15
    # test 1 lb1: 448251.0 5411932.0001
    # test 2 lb1: 448251, 5411932
    # test 3 lb1: WRF_Lb.WGS84
    # test 4 LatLon: 46.5°N, 003.0°E
    # test 5 LatLon: 46°30′00.0″N, 003°00′00.0″E
    # test 6 toLcc1: 700000 6600000
    # test 7 toLcc1: 46.5°N, 003.0°E
    # test 8 lb2: 1894411 1564650
    # test 9 lb2: Clarke1866
    # test 10 toLatLon2: 35.0°N, 075.0°W
    # test 11 toLatLon2: 35°00′00.0007″N, 074°59′59.9997″W
    # test 12 toLatLon2: NAD27
    # test 13 toLcc2: 1894410.9 1564649.5
    # test 14 toLcc2: Snyder.NAD27
    # test 15 Be72Lb: 49.833334°N, 003.034153°E
    # test 16 Be72Lb: NAD83
    # test 17 Be72Lb: 49.833334°N, 004.367487°E
    # test 18 Be72Lb: NAD83
    # test 19 Be72Lb: 49.833334°N, 005.70082°E
    # test 20 Be72Lb: NAD83
    # test 21 Be72Lb: 51.166667°N, 003.034153°E
    # test 22 Be72Lb: NAD83
    # test 23 Be72Lb: 51.166667°N, 004.367487°E
    # test 24 Be72Lb: NAD83
    # test 25 Be72Lb: 51.166667°N, 005.70082°E
    # test 26 Be72Lb: NAD83
    # test 27 Fr93Lb: 49.0°N, 002.0°W
    # test 28 Fr93Lb: WGS84
    # test 29 Fr93Lb: 49.0°N, 003.0°E
    # test 30 Fr93Lb: WGS84
    # test 31 Fr93Lb: 49.0°N, 008.0°E
    # test 32 Fr93Lb: WGS84
    # test 33 Fr93Lb: 44.0°N, 002.0°W
    # test 34 Fr93Lb: WGS84
    # test 35 Fr93Lb: 44.0°N, 003.0°E
    # test 36 Fr93Lb: WGS84
    # test 37 Fr93Lb: 44.0°N, 008.0°E
    # test 38 Fr93Lb: WGS84
    # test 39 MaNLb: 31.73°N, 008.54°W
    # test 40 MaNLb: NTF
    # test 41 MaNLb: 31.73°N, 005.4°W
    # test 42 MaNLb: NTF
    # test 43 MaNLb: 31.73°N, 002.26°W
    # test 44 MaNLb: NTF
    # test 45 MaNLb: 34.87°N, 008.54°W
    # test 46 MaNLb: NTF
    # test 47 MaNLb: 34.87°N, 005.4°W
    # test 48 MaNLb: NTF
    # test 49 MaNLb: 34.87°N, 002.26°W
    # test 50 MaNLb: NTF
    # test 51 MxLb: 17.5°N, 114.0°W
    # test 52 MxLb: WGS84
    # test 53 MxLb: 17.5°N, 102.0°W
    # test 54 MxLb: WGS84
    # test 55 MxLb: 17.5°N, 090.0°W
    # test 56 MxLb: WGS84
    # test 57 MxLb: 29.5°N, 114.0°W
    # test 58 MxLb: WGS84
    # test 59 MxLb: 29.5°N, 102.0°W
    # test 60 MxLb: WGS84
    # test 61 MxLb: 29.5°N, 090.0°W
    # test 62 MxLb: WGS84
    # test 63 PyT_Lb: 45.898939°N, 000.540154°E
    # test 64 PyT_Lb: NTF
    # test 65 PyT_Lb: 45.898939°N, 002.337229°E
    # test 66 PyT_Lb: NTF
    # test 67 PyT_Lb: 45.898939°N, 004.134305°E
    # test 68 PyT_Lb: NTF
    # test 69 PyT_Lb: 47.696014°N, 000.540154°E
    # test 70 PyT_Lb: NTF
    # test 71 PyT_Lb: 47.696014°N, 002.337229°E
    # test 72 PyT_Lb: NTF
    # test 73 PyT_Lb: 47.696014°N, 004.134305°E
    # test 74 PyT_Lb: NTF
    # test 75 Snyder: 33.0°N, 108.0°W
    # test 76 Snyder: NAD27
    # test 77 Snyder: 33.0°N, 096.0°W
    # test 78 Snyder: NAD27
    # test 79 Snyder: 33.0°N, 084.0°W
    # test 80 Snyder: NAD27
    # test 81 Snyder: 45.0°N, 108.0°W
    # test 82 Snyder: NAD27
    # test 83 Snyder: 45.0°N, 096.0°W
    # test 84 Snyder: NAD27
    # test 85 Snyder: 45.0°N, 084.0°W
    # test 86 Snyder: NAD27
    # test 87 USA_Lb: 33.0°N, 108.0°W
    # test 88 USA_Lb: WGS84
    # test 89 USA_Lb: 33.0°N, 096.0°W
    # test 90 USA_Lb: WGS84
    # test 91 USA_Lb: 33.0°N, 084.0°W
    # test 92 USA_Lb: WGS84
    # test 93 USA_Lb: 45.0°N, 108.0°W
    # test 94 USA_Lb: WGS84
    # test 95 USA_Lb: 45.0°N, 096.0°W
    # test 96 USA_Lb: WGS84
    # test 97 USA_Lb: 45.0°N, 084.0°W
    # test 98 USA_Lb: WGS84
    # test 99 WRF_Lb: 33.0°N, 109.0°W
    # test 100 WRF_Lb: WGS84
    # test 101 WRF_Lb: 33.0°N, 097.0°W
    # test 102 WRF_Lb: WGS84
    # test 103 WRF_Lb: 33.0°N, 085.0°W
    # test 104 WRF_Lb: WGS84
    # test 105 WRF_Lb: 45.0°N, 109.0°W
    # test 106 WRF_Lb: WGS84
    # test 107 WRF_Lb: 45.0°N, 097.0°W
    # test 108 WRF_Lb: WGS84
    # test 109 WRF_Lb: 45.0°N, 085.0°W
    # test 110 WRF_Lb: WGS84
    # test 111 Snyder1: lat0=23.0, lon0=-96.0, par1=33.0, par2=45.0, E0=0, N0=0, k0=1, SP=2, datum=(ellipsoid=Ellipsoids.Clarke1866, transform=Transforms.NAD27, name='NAD27'), name='Snyder1'
    # test 112 Snyder2: lat0=23.0, lon0=-96.0, par1=33.0, par2=45.0, E0=0, N0=0, k0=1, SP=2, datum=(ellipsoid=Ellipsoids.Clarke1866, transform=Transforms.NAD27, name='NAD27'), name='Snyder2'
    # all testLcc.py tests passed (Python 2.7.10)

    # testing testLcc.py version 16.12.15
    # test 1 lb1: 448251.0 5411932.0001
    # test 2 lb1: 448251, 5411932
    # test 3 lb1: WRF_Lb.WGS84
    # test 4 LatLon: 46.5°N, 003.0°E
    # test 5 LatLon: 46°30′00.0″N, 003°00′00.0″E
    # test 6 toLcc1: 700000 6600000
    # test 7 toLcc1: 46.5°N, 003.0°E
    # test 8 lb2: 1894411 1564650
    # test 9 lb2: Clarke1866
    # test 10 toLatLon2: 35.0°N, 075.0°W
    # test 11 toLatLon2: 35°00′00.0007″N, 074°59′59.9997″W
    # test 12 toLatLon2: NAD27
    # test 13 toLcc2: 1894410.9 1564649.5
    # test 14 toLcc2: Snyder.NAD27
    # test 15 Be72Lb: 49.833334°N, 003.034153°E
    # test 16 Be72Lb: NAD83
    # test 17 Be72Lb: 49.833334°N, 004.367487°E
    # test 18 Be72Lb: NAD83
    # test 19 Be72Lb: 49.833334°N, 005.70082°E
    # test 20 Be72Lb: NAD83
    # test 21 Be72Lb: 51.166667°N, 003.034153°E
    # test 22 Be72Lb: NAD83
    # test 23 Be72Lb: 51.166667°N, 004.367487°E
    # test 24 Be72Lb: NAD83
    # test 25 Be72Lb: 51.166667°N, 005.70082°E
    # test 26 Be72Lb: NAD83
    # test 27 Fr93Lb: 49.0°N, 002.0°W
    # test 28 Fr93Lb: WGS84
    # test 29 Fr93Lb: 49.0°N, 003.0°E
    # test 30 Fr93Lb: WGS84
    # test 31 Fr93Lb: 49.0°N, 008.0°E
    # test 32 Fr93Lb: WGS84
    # test 33 Fr93Lb: 44.0°N, 002.0°W
    # test 34 Fr93Lb: WGS84
    # test 35 Fr93Lb: 44.0°N, 003.0°E
    # test 36 Fr93Lb: WGS84
    # test 37 Fr93Lb: 44.0°N, 008.0°E
    # test 38 Fr93Lb: WGS84
    # test 39 MaNLb: 31.73°N, 008.54°W
    # test 40 MaNLb: NTF
    # test 41 MaNLb: 31.73°N, 005.4°W
    # test 42 MaNLb: NTF
    # test 43 MaNLb: 31.73°N, 002.26°W
    # test 44 MaNLb: NTF
    # test 45 MaNLb: 34.87°N, 008.54°W
    # test 46 MaNLb: NTF
    # test 47 MaNLb: 34.87°N, 005.4°W
    # test 48 MaNLb: NTF
    # test 49 MaNLb: 34.87°N, 002.26°W
    # test 50 MaNLb: NTF
    # test 51 MxLb: 17.5°N, 114.0°W
    # test 52 MxLb: WGS84
    # test 53 MxLb: 17.5°N, 102.0°W
    # test 54 MxLb: WGS84
    # test 55 MxLb: 17.5°N, 090.0°W
    # test 56 MxLb: WGS84
    # test 57 MxLb: 29.5°N, 114.0°W
    # test 58 MxLb: WGS84
    # test 59 MxLb: 29.5°N, 102.0°W
    # test 60 MxLb: WGS84
    # test 61 MxLb: 29.5°N, 090.0°W
    # test 62 MxLb: WGS84
    # test 63 PyT_Lb: 45.898939°N, 000.540154°E
    # test 64 PyT_Lb: NTF
    # test 65 PyT_Lb: 45.898939°N, 002.337229°E
    # test 66 PyT_Lb: NTF
    # test 67 PyT_Lb: 45.898939°N, 004.134305°E
    # test 68 PyT_Lb: NTF
    # test 69 PyT_Lb: 47.696014°N, 000.540154°E
    # test 70 PyT_Lb: NTF
    # test 71 PyT_Lb: 47.696014°N, 002.337229°E
    # test 72 PyT_Lb: NTF
    # test 73 PyT_Lb: 47.696014°N, 004.134305°E
    # test 74 PyT_Lb: NTF
    # test 75 Snyder: 33.0°N, 108.0°W
    # test 76 Snyder: NAD27
    # test 77 Snyder: 33.0°N, 096.0°W
    # test 78 Snyder: NAD27
    # test 79 Snyder: 33.0°N, 084.0°W
    # test 80 Snyder: NAD27
    # test 81 Snyder: 45.0°N, 108.0°W
    # test 82 Snyder: NAD27
    # test 83 Snyder: 45.0°N, 096.0°W
    # test 84 Snyder: NAD27
    # test 85 Snyder: 45.0°N, 084.0°W
    # test 86 Snyder: NAD27
    # test 87 USA_Lb: 33.0°N, 108.0°W
    # test 88 USA_Lb: WGS84
    # test 89 USA_Lb: 33.0°N, 096.0°W
    # test 90 USA_Lb: WGS84
    # test 91 USA_Lb: 33.0°N, 084.0°W
    # test 92 USA_Lb: WGS84
    # test 93 USA_Lb: 45.0°N, 108.0°W
    # test 94 USA_Lb: WGS84
    # test 95 USA_Lb: 45.0°N, 096.0°W
    # test 96 USA_Lb: WGS84
    # test 97 USA_Lb: 45.0°N, 084.0°W
    # test 98 USA_Lb: WGS84
    # test 99 WRF_Lb: 33.0°N, 109.0°W
    # test 100 WRF_Lb: WGS84
    # test 101 WRF_Lb: 33.0°N, 097.0°W
    # test 102 WRF_Lb: WGS84
    # test 103 WRF_Lb: 33.0°N, 085.0°W
    # test 104 WRF_Lb: WGS84
    # test 105 WRF_Lb: 45.0°N, 109.0°W
    # test 106 WRF_Lb: WGS84
    # test 107 WRF_Lb: 45.0°N, 097.0°W
    # test 108 WRF_Lb: WGS84
    # test 109 WRF_Lb: 45.0°N, 085.0°W
    # test 110 WRF_Lb: WGS84
    # test 111 Snyder1: lat0=23.0, lon0=-96.0, par1=33.0, par2=45.0, E0=0, N0=0, k0=1, SP=2, datum=(ellipsoid=Ellipsoids.Clarke1866, transform=Transforms.NAD27, name='NAD27'), name='Snyder1'
    # test 112 Snyder2: lat0=23.0, lon0=-96.0, par1=33.0, par2=45.0, E0=0, N0=0, k0=1, SP=2, datum=(ellipsoid=Ellipsoids.Clarke1866, transform=Transforms.NAD27, name='NAD27'), name='Snyder2'
    # all testLcc.py tests passed (Python 3.6.0)
