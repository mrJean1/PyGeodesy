
# -*- coding: utf-8 -*-

u'''Test projection L{ExactTransverseMercator} with the C(TMcoords.dat} from
U{C.F.F. Karney, "Test data for the transverse Mercator projection (2009)"
<https://GeographicLib.SourceForge.io/html/transversemercator.html>},
also available U{here<https://Zenodo.org/record/32470>}, file C{TMcoords.dat}.
'''

__all__ = ('testExactTMcoords',)
__version__ = '19.05.21'

from testTMcoords import testTMcoords

from pygeodesy import etm, ExactTransverseMercator, LatLon_

_ETM = ExactTransverseMercator()  # default WGS84


class _ExactTM(object):
    '''Minimal C{Etm} class to test C{ExactTransverseMercator}.
    '''
    _datum  = _ETM.datum
    _scale0 = _ETM.k0

    def __init__(self, lat, lon, **unused):
        self.easting, self.northing, self.convergence, self.scale = _ETM.forward(lat, lon)

    def toLatLon(self, **unused):
        lat, lon, _, _ = _ETM.reverse(self.easting, self.northing)
        return LatLon_(lat, lon)


def testExactTMcoords(name):

    testTMcoords(etm, _ExactTM, name=name, eps1=4e-8, eps2=4e-8, lonE=360, Etm=_ExactTM)


if __name__ == '__main__':

    testExactTMcoords(__file__)
