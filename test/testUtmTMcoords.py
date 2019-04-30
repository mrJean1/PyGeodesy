
# -*- coding: utf-8 -*-

u'''Test UTM function with the C(TMcoords.dat} from
U{C.F.F. Karney, "Test data for the transverse Mercator projection (2009)"
<http://GeographicLib.SourceForge.io/html/transversemercator.html>},
also available U{here<http://Zenodo.org/record/32470>}, file C{TMcoords.dat}.
'''

__all__ = ('Tests',)
__version__ = '19.04.26'

from base import TestsBase

from pygeodesy import Ellipsoids, RangeError, toUtm, utm


class Tests(TestsBase):

    def testUtmTMcoord(self, coord, line, fmt='%.4f', eps=1.5e-4):
        # format: lat lon easting northing convergence scale
        lat, lon, e1, n1, c1, s1 = map(float, coord.split())
        # skip tests with "out of range" lon
        if lon > 70.0:
            self.skip(line + repr(coord))
        else:
            try:
                _, e2, n2, _, c2, s2 = toUtm(lat, lon, cmoff=False)
                self.test(line + 'easting',     e2, e1, fmt=fmt, known=abs(e2 - e1) < eps)
                self.test(line + 'northing',    n2, n1, fmt=fmt, known=abs(e2 - e1) < eps)
                self.test(line + 'convergence', c2, c1, fmt=fmt)
                self.test(line + 'scale',       s2, s1, fmt=fmt)
            except RangeError as x:
                self.test(line + 'RangeError', x, None, known=True)


if __name__ == '__main__':

    import sys

    if len(sys.argv) > 1:  # get "TMcoords.dat" file
        coords = open(sys.argv[1], 'rb')
        v = False
    else:
        from testEpsg import _TMcoords as coords
        v = True

    t = Tests(__file__, __version__, utm, verbose=v)

    for n, coord in enumerate(coords.readlines()):
        t.testUtmTMcoord(coord.rstrip(), 'line %d ' % (n + 1,))

    # XXX Pythonista run_path doesn't reload modules
    E = Ellipsoids.WGS84
    t.test(E.name + '.KsOrder', E.KsOrder, 8)

    coords.close()

    t.results()
