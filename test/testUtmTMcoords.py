
# -*- coding: utf-8 -*-

# Test UTM functions with the C(TMcoords.dat} from U{C.F.F. Karney,
# "Test data for the transverse Mercator projection (2009)"
# <https://GeographicLib.SourceForge.io/html/transversemercator.html>},
# also available U{here<https://Zenodo.org/record/32470>}, file C{TMcoords.dat}.

__all__ = ('testUtmTMcoords',)
__version__ = '23.05.23'

from testTMcoords import testTMcoords


def testUtmTMcoords(name):

    from pygeodesy import toUtm8, Utm, utm

    testTMcoords(utm, toUtm8, name=name, Utm=Utm)


if __name__ == '__main__':

    testUtmTMcoords(__file__)
