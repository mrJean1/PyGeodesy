
# -*- coding: utf-8 -*-

u'''Test L{etm} module and L{Etm} class with the C(TMcoords.dat} from
U{C.F.F. Karney, "Test data for the transverse Mercator projection (2009)"
<https://GeographicLib.SourceForge.io/html/transversemercator.html>},
also available U{here<https://Zenodo.org/record/32470>}, file C{TMcoords.dat}.
'''

__all__ = ('testEtmTMcoords',)
__version__ = '20.12.06'

from testTMcoords import testTMcoords


def testEtmTMcoords(name):

    from pygeodesy import etm, Etm, toEtm8

    testTMcoords(etm, toEtm8, name=name, eps1=2e-7, eps2=1e-7, lonE=360, Etm=Etm)


if __name__ == '__main__':

    testEtmTMcoords(__file__)
