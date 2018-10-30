
# -*- coding: utf-8 -*-

u'''A pure Python implementation of geodesy tools for various ellipsoidal
and spherical earth models using precision trigonometric, vector-based
and approximate methods for geodetic (lat-/longitude) and geocentric
cartesian (x/y/z) coordinates.

Transcribed from U{JavaScript originals<http://GitHub.com/ChrisVeness/geodesy>}
by I{Chris Veness (C) 2005-2016} and published under the same U{MIT License
<http://OpenSource.org/licenses/MIT>}**.

There are three modules for ellipsoidal earth models, I{ellipsoidalKarney},
I{-Vincenty} and I{-Nvector} and two for spherical ones, I{sphericalTrigonometry}
and I{-Nvector}.  Each module provides a I{attributes-LatLon-html} class
with methods and functions to compute distance, initial and final bearing,
intermediate and nearest points, area, perimeter, conversions and unrolling,
among other things.  For more information and further details see the
U{documentation<http://mrJean1.GitHub.io/PyGeodesy>}, the descriptions of
U{Latitude/Longitude<http://www.Movable-Type.co.UK/scripts/latlong.html>},
U{Vincenty<http://www.Movable-Type.co.UK/scripts/latlong-vincenty.html>} and
U{Vector-based<http://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}
geodesy, the original U{JavaScript source<http://GitHub.com/ChrisVeness/geodesy>} or
U{docs<http://www.Movable-Type.co.UK/scripts/geodesy/docs>} and the Python
U{GeographicLib<http://PyPI.org/project/geographiclib>}.

Also included are modules for conversions to and from
U{UTM<http://www.Movable-Type.co.UK/scripts/latlong-utm-mgrs.html>}
(Universal Transverse Mercator) and U{Web Mercator
<http://WikiPedia.org/wiki/Web_Mercator>} (Pseudo-Mercator) coordinates,
U{MGRS<http://www.Movable-Type.co.UK/scripts/latlong-utm-mgrs.html>}
(NATO Military Grid Reference System) and
U{OSGR<http://www.Movable-Type.co.UK/scripts/latlong-os-gridref.html>}
(British Ordinance Survery Grid Reference) grid references and a module for
encoding and decoding U{Geohashes<http://www.Movable-Type.co.UK/scripts/geohash.html>}.

Two other modules provide Lambert conformal conic projections and positions
(from U{John P. Snyder, "Map Projections -- A Working Manual", 1987, pp 107-109
<http://pubs.er.USGS.gov/djvu/PP/PP_1395.pdf>}) and several functions to
U{simplify<http://Bost.Ocks.org/mike/simplify>} or linearize a path of
C{LatLon} points (or a U{NumPy array
<http://docs.SciPy.org/doc/numpy/reference/generated/numpy.array.html>}),
including implementations of the U{Ramer-Douglas-Peucker
<http://WikiPedia.org/wiki/Ramer-Douglas-Peucker_algorithm>}, the
U{Visvalingam-Whyatt<http://hydra.Hull.ac.UK/resources/hull:8338>} and the
U{Reumann-Witkam<http://psimpl.SourceForge.net/reumann-witkam.html>}
algorithms and modified versions of the former.

All Python source code has been statically U{checked
<http://GitHub.com/ActiveState/code/tree/master/recipes/Python/546532_PyChecker_postprocessor>}
with U{PyChecker<http://PyPI.org/project/pychecker>},
U{PyFlakes<http://PyPI.org/project/pyflakes>},
U{PyCodeStyle<http://PyPI.org/project/pycodestyle>} (formerly Pep8) and
U{McCabe<http://PyPI.org/project/mccabe>} using Python 2.7.15 and with
U{Flake8<http://PyPI.org/project/flake8>} using Python 3.7.1, both in
64-bit on macOS 10.13.6 High Sierra.

The tests have been run with Python 2.7.15 (with U{geographiclib
<http://PyPI.org/project/geographiclib>} 1.49 and U{numpy
<http://PyPI.org/project/numpy>} 1.15.2), with Python 3.7.1 (with
U{geographiclib<http://PyPI.org/project/geographiclib>} 1.49) and with
U{PyPy<http://PyPy.org>} 6.0.0 (Python 2.7.13 and 3.5.3) on macOS
10.13.6 High Sierra, with Python 2.6.9, 2.7.14, 3.5.6 and 3.6.3 (and
U{geographiclib<http://PyPI.org/project/geographiclib>} 1.49) on
U{Debian 8<http://Travis-CI.org/mrJean1/PyGeodesy>} and with Python
3.7.0 (and U{geographiclib<http://PyPI.org/project/geographiclib>} 1.49)
on U{Debian 9<http://Cirrus-CI.com/github/mrJean1/PyGeodesy/master>},
all in 64-bit only and with Python 2.7.15, 3.6.6 and 3.7.0 (all with
U{geographiclib<http://PyPI.org/project/geographiclib>} 1.49) on
U{Windows Server 2012R2<http://CI.AppVeyor.com/project/mrJean1/pygeodesy>}
in 32- and 64-bit.

Previously, the tests were run with Python 2.6.9 (and numpy 1.6.2), 2.7.10
(and numpy 1.8.0rc1), 2.7.13, 2.7.14 (and numpy 1.13.1 or 1.14.0), 3.5.3,
3.6.2, 3.6.3, 3.6.4, 3.6.5, 3.7.0 and U{Intel-Python
<http://software.Intel.com/en-us/distribution-for-python>} 3.5.3 (and
U{numpy<http://PyPI.org/project/numpy>} 1.11.3) on MacOS X 10.10 Yosemite,
MacOS X 10.11 El Capitan, macOS 10.12 Sierra, macOS 10.13.5 High Sierra and
macOS 10.14 Mojave, with U{Pythonista 3.1<http://OMZ-Software.com/pythonista>}
on iOS 10.3.3, 11.0.3, 11.1.2 and 11.3 on iPad4, with U{Pythonista 3.2
<http://OMZ-Software.com/pythonista>} on iOS 11.4.1 and 12.0 on iPad4,
iPhone7 and/or iPhone10, all in 64-bit only and with 32-bit Python 2.6.6
on Windows XP SP3 and with 32-bit Python 2.7.14 on Windows 10 Pro.

In addition to the U{PyGeodesy<http://PyPI.org/project/PyGeodesy>}
package, the distribution files contain the tests, the test results (on
macOS only) and the complete documentation (generated by U{Epydoc
<http://PyPI.org/project/epydoc>} using command line: C{epydoc --html
--no-private --no-source --name=PyGeodesy --url=... -v pygeodesy}).

To install PyGeodesy, type C{pip install PyGeodesy} or C{easy_install
PyGeodesy} in a terminal or command window.  Alternatively, download
C{PyGeodesy-yy.m.d.zip} from U{PyPI<http://PyPI.org/project/PyGeodesy>}
or U{GitHub<http://GitHub.com/mrJean1/PyGeodesy>}, C{unzip} the downloaded
file, C{cd} to directory C{Pygeodesy-yy.m.d} and type C{python setup.py
install}.  To run all PyGeodesy tests, type C{python setup.py test}
before installation.

Installation of U{NumPy<http://www.NumPy.org>} and U{GeographicLib
<http://PyPI.org/project/geographiclib>} is optional.  However, the
latter is required for module I{ellipsoidalKarney} classes C{LatLon} and
I{Cartesian} and functions I{areaOf} and I{perimeterOf}.

Some function and method names differ from the JavaScript version. In such
cases documentation tag B{JS name:} shows the original JavaScript name.

__

**) U{Copyright (C) 2016-2018 -- mrJean1 at Gmail dot com
<http://OpenSource.org/licenses/MIT>}

C{Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:}

C{The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.}

C{THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.}

@newfield example: Example, Examples
@newfield JSname: JS name, JS names

@var EPS:  System's M{epsilon} (C{float})
@var EPS1: M{1 - EPS} (C{float}), about 0.9999999999999998

@var F_D:   Format degrees as deg° (C{str}).
@var F_DM:  Format degrees as deg°min′ (C{str}).
@var F_DMS: Format degrees as deg°min′sec″ (C{str}).
@var F_DEG: Format degrees as [D]DD without symbol (C{str}).
@var F_MIN: Format degrees as [D]DDMM without symbols (C{str}).
@var F_SEC: Format degrees as [D]DDMMSS without symbols (C{str}).
@var F_RAD: Convert degrees to radians and format as RR (C{str}).

@var PI:   Constant M{math.pi} (C{float})
@var PI2:  Two PI, M{math.pi * 2} (C{float})
@var PI_2: Half PI, M{math.pi / 2} (C{float})

@var R_M:  Mean, spherical earth radius (C{meter}).
@var R_MA: Equatorial earth radius (C{meter}) WGS84, EPSG:3785.
@var R_MB: Polar earth radius (C{meter}) WGS84, EPSG:3785.
@var R_KM: Mean, spherical earth radius (C{km}, kilometer).
@var R_NM: Mean, spherical earth radius (C{NM}, nautical miles).
@var R_SM: Mean, spherical earth radius (C{SM}, statute miles).
@var R_FM: Former FAI Sphere earth radius (C{meter}).
@var R_VM: Aviation/Navigation earth radius (C{meter}).

@var S_DEG: Degrees symbol "°" (C{str}).
@var S_MIN: Minutes symbol "′" (C{str}).
@var S_SEC: Seconds symbol "″" (C{str}).
@var S_RAD: Radians symbol "" (C{str}).
@var S_SEP: Separator between deg°, min′ and sec″ "" (C{str}).

@var Conics:     Registered conics (C{enum-like}).
@var Datums:     Registered datums (C{enum-like}).
@var Ellipsoids: Registered ellipsoids (C{enum-like}).
@var Transforms: Registered transforms (C{enum-like}).

@var pygeodesy_abspath: Fully qualified C{pygeodesy} directory name (C{str}).

@var version: Normalized C{PyGeodesy} version (C{str}).

'''

from os.path import abspath, basename, dirname, splitext

_init_abspath     = abspath(__file__)
pygeodesy_abspath = dirname(_init_abspath)

# setting __path__ should make ...
__path__ = [pygeodesy_abspath]
try:  # ... this import work, ...
    import bases as _  # PYCHOK expected
    del _
except ImportError:  # ... if it doesn't, extend
    # sys.path to include this very directory such
    # that all public and private sub-modules can
    # be imported (and checked by PyChecker, etc.)
    import sys
    sys.path.insert(0, pygeodesy_abspath)  # XXX __path__[0]
    del sys

# keep ellipsoidal, spherical and vector modules as sub-modules
import ellipsoidalKarney  # PYCHOK false
import ellipsoidalNvector  # PYCHOK false
import ellipsoidalVincenty  # PYCHOK false
import geohash
import nvector  # PYCHOK false
import sphericalNvector  # PYCHOK false
import sphericalTrigonometry  # PYCHOK false
import vector3d  # PYCHOK false

CrossError    = vector3d.CrossError
crosserrors   = vector3d.crosserrors
Geohash       = geohash.Geohash
VincentyError = ellipsoidalVincenty.VincentyError

# all public sub-modules, contants, classes and functions
__all__ = ('bases', 'datum', 'dms', 'elevations',  # modules
           'ellipsoidalKarney', 'ellipsoidalNvector', 'ellipsoidalVincenty',
           'fmath', 'formy', 'geohash', 'lcc', 'mgrs',
           'nvector', 'osgr', 'points',
           'simplify', 'sphericalNvector', 'sphericalTrigonometry',
           'utily', 'utm', 'vector3d', 'webmercator',
           'CrossError', 'Geohash', 'VincentyError',  # classes
           'R_M',  # to avoid duplicates from .datum.py and .utily.py
           'pygeodesy_abspath',
           'version',
           'crosserrors')  # extended below
__version__ = '18.10.29'

# see setup.py for similar logic
version = '.'.join(map(str, map(int, __version__.split('.'))))

# lift all public classes, constants, functions, etc. but
# only from the following sub-modules ... (see also David
# Beazley's <http://DaBeaz.com/modulepackage/index.html>)
from bases       import *  # PYCHOK __all__
from datum       import *  # PYCHOK __all__
from dms         import *  # PYCHOK __all__
from elevations  import *  # PYCHOK __all__
from fmath       import *  # PYCHOK __all__
from formy       import *  # PYCHOK __all__
from lcc         import *  # PYCHOK __all__
from mgrs        import *  # PYCHOK __all__
from osgr        import *  # PYCHOK __all__
from points      import *  # PYCHOK __all__
from simplify    import *  # PYCHOK __all__
from utily       import *  # PYCHOK __all__
from utm         import *  # PYCHOK __all__
from webmercator import *  # PYCHOK __all__

import bases        # PYCHOK expected
import datum        # PYCHOK expected
import dms          # PYCHOK expected
import elevations   # PYCHOK expected
import fmath        # PYCHOK expected
import formy        # PYCHOK expected
import lcc          # PYCHOK expected
import mgrs         # PYCHOK expected
import osgr         # PYCHOK expected
import points       # PYCHOK expected
import simplify     # PYCHOK expected
import utily        # PYCHOK expected
import utm          # PYCHOK expected
import webmercator  # PYCHOK expected

# for backward compatibility with old, DEPRECATED names
areaof      = points.areaOf
perimeterof = points.perimeterOf

ellipsoidalVincenty.areaOf      = ellipsoidalKarney.areaOf
ellipsoidalVincenty.perimeterOf = ellipsoidalKarney.perimeterOf


def equirectangular3(lat1, lon1, lat2, lon2, **options):
    '''DEPRECATED, use function C{equirectangular_}.

       @return: 3-Tuple (distance2, delta_lat, delta_lon).
    '''
    return formy.equirectangular_(lat1, lon1, lat2, lon2, **options)[:3]


def _ismodule(m):
    p = abspath(m.__file__)  # PYCHOK undefined?
    if dirname(p) != pygeodesy_abspath:  # PYCHOK undefined?
        raise ImportError('foreign module %r from %r' % (m.__name__, p))


# check that all modules are from this very package, pygeodesy
for m in (ellipsoidalKarney, ellipsoidalNvector, ellipsoidalVincenty,
          geohash, nvector,
          sphericalNvector, sphericalTrigonometry,
          vector3d):
    _ismodule(m)

# concat __all__ with the public classes, constants,
# functions, etc. from the sub-modules mentioned above
for m in (bases, datum, dms, elevations, fmath, formy, lcc, mgrs, osgr,
          points, simplify, utily, utm, webmercator):
    __all__ += m.__all__
    _ismodule(m)

# remove any duplicates, only R_M?
__all__ = tuple(set(__all__))

if __name__ == '__main__':

    n = basename(pygeodesy_abspath)
    print('%s %s (%s)' % (n, __version__, pygeodesy_abspath))

# XXX del ellipsoidalBase, sphericalBase  # PYCHOK expected
del abspath, basename, dirname, _init_abspath, _ismodule, m, splitext

# **) MIT License
#
# Copyright (C) 2016-2018 -- mrJean1 at Gmail dot com
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
