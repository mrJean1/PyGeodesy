
# -*- coding: utf-8 -*-

# Test UTM/UPS function with the C(TMcoords.dat} from U{C.F.F. Karney,
# "Test data for the transverse Mercator projection (2009)"
# <https://GeographicLib.SourceForge.io/html/transversemercator.html>},
# also available U{here<https://Zenodo.org/record/32470>}, file C{TMcoords.dat}.

__all__ = ('testUtmUpsTMcoords',)
__version__ = '23.03.27'

from testTMcoords import testTMcoords


def testUtmUpsTMcoords(name):

    from pygeodesy import toUtmUps8, Ups, Utm, utmups

    testTMcoords(utmups, toUtmUps8, name=name, Utm=Utm, Ups=Ups)


if __name__ == '__main__':

    testUtmUpsTMcoords(__file__)
