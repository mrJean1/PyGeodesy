
# -*- coding: utf-8 -*-

u'''A pure Python implementation of geodesy tools for various ellipsoidal
and spherical earth models using precision trigonometric, vector-based,
elliptic and approximate methods for geodetic (lat-/longitude) and
geocentric (U{ECEF<https://WikiPedia.org/wiki/ECEF>} cartesian) coordinates.

Transcribed from U{JavaScript originals<https://GitHub.com/ChrisVeness/geodesy>}
by I{Chris Veness (C) 2005-2016} and several U{C++ classes
<https://GeographicLib.SourceForge.io/html/annotated.html>} by I{Charles Karney
(C) 2008-2019} and published under the same U{MIT License
<https://OpenSource.org/licenses/MIT>}**.

There are three modules for ellipsoidal earth models, C{ellipsoidalKarney},
C{-Vincenty} and C{-Nvector} and two for spherical ones, C{sphericalTrigonometry}
and C{-Nvector}.  Each module provides a geodetic B{C{LatLon}} and a geocentric
B{C{Cartesian}} class with methods and functions to compute distance, initial and
final bearing, intermediate and nearest points, area, perimeter, conversions and
unrolling, among other things.  For more information and further details see the
U{documentation<https://mrJean1.GitHub.io/PyGeodesy>}, the descriptions of
U{Latitude/Longitude<https://www.Movable-Type.co.UK/scripts/latlong.html>},
U{Vincenty<https://www.Movable-Type.co.UK/scripts/latlong-vincenty.html>} and
U{Vector-based<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>}
geodesy, the original U{JavaScript source<https://GitHub.com/ChrisVeness/geodesy>} or
U{docs<https://www.Movable-Type.co.UK/scripts/geodesy/docs>} and the Python
U{GeographicLib<https://PyPI.org/project/geographiclib>}.

Also included are modules for conversions to and from U{Cassini-Soldner
<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1CassiniSoldner.html>},
U{ECEF<https://WikiPedia.org/wiki/ECEF>} (Earth-Centered, Earth-Fixed cartesian),
U{UPS<https://WikiPedia.org/wiki/Universal_polar_stereographic_coordinate_system>}
(Universal Polar Stereographic), U{UTM
<https://www.Movable-Type.co.UK/scripts/latlong-utm-mgrs.html>} (U{Exact
<https://GeographicLib.SourceForge.io/html/classGeographicLib_1_1TransverseMercatorExact.html>}
and Universal Transverse Mercator) and U{Web Mercator
<https://WikiPedia.org/wiki/Web_Mercator>} (Pseudo-Mercator) coordinates,
U{MGRS<https://www.Movable-Type.co.UK/scripts/latlong-utm-mgrs.html>} (NATO
Military Grid Reference System) and U{OSGR
<https://www.Movable-Type.co.UK/scripts/latlong-os-gridref.html>}
(British Ordinance Survery Grid Reference) grid references, U{TRF
<http://ITRF.ENSG.IGN.FR>} (Terrestrial Reference Frames) and modules
to encode and decode U{EPSG<https://www.EPSG-Registry.org>}, U{Geohashes
<https://www.Movable-Type.co.UK/scripts/geohash.html>}, U{Georefs (WGRS)
<https://WikiPedia.org/wiki/World_Geographic_Reference_System>} and
U{Garefs (GARS)<https://WikiPedia.org/wiki/Global_Area_Reference_System>}.

Other modules provide azimuthal projections and Lambert conformal conic projections
and positions (from U{John P. Snyder, "Map Projections -- A Working Manual",
1987<https://pubs.er.USGS.gov/djvu/PP/PP_1395.pdf>}), functions to clip a path
or polygon of C{LatLon} points using the U{Cohen–Sutherland
<https://WikiPedia.org/wiki/Cohen-Sutherland_algorithm>} and the
U{Sutherland-Hodgman<https://WikiPedia.org/wiki/Sutherland-Hodgman_algorithm>}
methods, functions to U{simplify<https://Bost.Ocks.org/mike/simplify>} or
linearize a path of C{tLon} points (or a U{NumPy array
<https://docs.SciPy.org/doc/numpy/reference/generated/numpy.array.html>}),
including implementations of the U{Ramer-Douglas-Peucker
<https://WikiPedia.org/wiki/Ramer-Douglas-Peucker_algorithm>} the
U{Visvalingam-Whyatt<https://hydra.Hull.ac.UK/resources/hull:8338>} and
U{Reumann-Witkam<https://psimpl.SourceForge.net/reumann-witkam.html>}
the algorithms and modified versions of the former.  Other classes
U{interpolate<https://docs.SciPy.org/doc/scipy/reference/interpolate.html>}
the height of C{LatLon} points and several C{Geoid} models or compute
various U{Fréchet<https://WikiPedia.org/wiki/Frechet_distance>} or
U{Hausdorff<https://WikiPedia.org/wiki/Hausdorff_distance>} distances.

Installation
============

To install PyGeodesy, type C{pip install PyGeodesy} or C{easy_install
PyGeodesy} in a terminal or command window.

Alternatively, download C{PyGeodesy-yy.m.d.zip} from U{PyPI
<https://PyPI.org/project/PyGeodesy>} or U{GitHub
<https://GitHub.com/mrJean1/PyGeodesy>}, C{unzip} the downloaded file,
C{cd} to directory C{Pygeodesy-yy.m.d} and type C{python setup.py
install}.  To run all PyGeodesy tests, type C{python setup.py test}
before installation.

Installation of U{GeographicLib<https://PyPI.org/project/geographiclib>},
U{NumPy<https://www.NumPy.org>} and U{SciPy<https://SciPy.org>} is optional.
However, the former is required for classes L{CassiniSoldner} and L{Css} and
function L{toCss}, for module L{ellipsoidalKarney} classes C{LatLon} and
C{Cartesian} and functions C{areaOf} and C{perimeterOf} and for the
L{HeightIDWkarney} interpolator.  The latter are needed for the C{Geoid...}
and C{Height...} interpolator classes, except the L{GeoidKarney} and all
C{HeightIDW...} classes.

Documentation
=============

In addition to the C{pygeodesy} package, the U{PyGeodesy
<https://PyPI.org/project/PyGeodesy>} distribution files contain the tests,
the test results (on macOS only) and the complete U{documentation
<https://mrJean1.GitHub.io/PyGeodesy>} (generated by U{Epydoc
<https://PyPI.org/project/epydoc>} using command line: C{epydoc --html
--no-private --no-source --name=PyGeodesy --url=... -v pygeodesy}).

Tests
=====

The tests have been run with Python 3.8.5, 3.7.6 and 2.7.18 (all with
U{geographiclib<https://PyPI.org/project/geographiclib>} 1.50, U{numpy
<https://PyPI.org/project/numpy>} 1.19.0, 1.17.2 respectively 1.16.5 and
U{scipy<https://SciPy.org/scipylib/download.html>} 1.5.0, 1.3.1 respectively
1.2.2) and with Python 3.9.0b5 and macOS' Python 2.7.16 (both without
geographiclib, numpy and scipy), all on macOS 10.15.6 Catalina and all in
64-bit only.  The tests run with and without C{lazy import} for Python 3.
The results of those tests are included in the distribution files.

Test coverage has been measured with U{coverage
<https://PyPI.org/project/coverage>} 4.5.4 using Python 3.8.5 and 3.7.6
(both with U{geographiclib<https://PyPI.org/project/geographiclib>} 1.50,
U{numpy<https://PyPI.org/project/numpy>} 1.19.0 respectively 1.17.2 and
U{scipy<https://SciPy.org/scipylib/download.html>} 1.5.0 respectively
1.3.1) and macOS' Python 2.7.16 (without geographiclib, numpy and scipy).
The full HMTL report and a PDF summary are included in the distribution
files.

The tests also ran with Python 2.7.14, 3.5.6, 3.6.3, 3.7.1, 3.8.0 and
U{PyPy<https://PyPy.org>} 7.1.1 (Python 2.7.13 and 3.6.1) (and
U{geographiclib<https://PyPI.org/project/geographiclib>} 1.49 or 1.50) on
U{Ubuntu 14.04<https://Travis-CI.org/mrJean1/PyGeodesy>} and with Python
3.7.3 (and U{geographiclib<https://PyPI.org/project/geographiclib>} 1.49 or
1.50) on U{Debian 9<https://Cirrus-CI.com/github/mrJean1/PyGeodesy/master>}
all in 64-bit only and with Python 2.7.15, 3.6.8, 3.7.2, and 3.8.0 (all
with U{geographiclib<https://PyPI.org/project/geographiclib>} 1.49 or 1.50)
on U{Windows Server 2012R2<https://CI.AppVeyor.com/project/mrJean1/pygeodesy>}
in both 32- and 64-bit.

A single-File and single-Directory application with C{pygeodesy} has been
bundled using U{PyInstaller<https://www.PyInstaller.org>} 3.4 and 64-bit
Python 3.7.3 on macOS 10.13.6 High Sierra.

Previously, the tests were run with Python 2.6.9 (and numpy 1.6.2), 2.7.10
(and numpy 1.8.0rc1), 2.7.13 thru 2.7.17 (and numpy 1.13.1, 1.14.0, 1.15.2
or 1.16.2), 3.5.3, 3.6.2 thru 3.6.5, 3.7.0, 3.7.2 thru 3.7.5, 3.8 thru 3.8.3,
U{PyPy<https://PyPy.org>} 6.0.0 (Python 2.7.13 and 3.5.3), U{PyPy
<https://PyPy.org>} 7.3.0 (Python 2.7.13 and 3.6.9) and U{Intel-Python
<https://software.Intel.com/en-us/distribution-for-python>} 3.5.3 (and U{numpy
<https://PyPI.org/project/numpy>} 1.11.3) on MacOS X 10.10 Yosemite, MacOS X 10.11
El Capitan, macOS 10.12 Sierra, macOS 10.13.6 High Sierra, macOS 10.14 Mojave
and/or macOS 10.15.3 and 10.15.5 Catalina, with U{Pythonista 3.1
<https://OMZ-Software.com/pythonista>} on iOS 10.3.3, 11.0.3, 11.1.2 and 11.3
on iPad4, with U{Pythonista 3.2<https://OMZ-Software.com/pythonista>} (with
geographiclib 1.49 or 1.50 and numpy 1.8.0) on iOS 11.4.1, 12.0, 12.2 and 12.3
on iPad4, iPhone6 and/or iPhone10, all in 64-bit only and with 32-bit Python
2.6.6 on Windows XP SP3 and with 32-bit Python 2.7.14 on Windows 10 Pro.

Notes
=====

All Python source code has been statically U{checked
<https://GitHub.com/ActiveState/code/tree/master/recipes/Python/546532_PyChecker_postprocessor>}
with U{PyChecker<https://PyPI.org/project/pychecker>}, U{PyFlakes
<https://PyPI.org/project/pyflakes>}, U{PyCodeStyle
<https://PyPI.org/project/pycodestyle>} (formerly Pep8) and U{McCabe
<https://PyPI.org/project/mccabe>} using Python 2.7.18 and with U{Flake8
<https://PyPI.org/project/flake8>} using Python 3.8.3, both in 64-bit
on macOS 10.15.5 Catalina.

Some function and method names differ from the JavaScript version. In such
cases documentation tag B{JS name:} shows the original JavaScript name.

Classes with a name ending in C{-Karney} are transcribed from Karney's
U{C++ classes<https://GeographicLib.SourceForge.io/html/annotated.html>},
but there are more.  The complete list is in module L{karney}.

License
=======

**) U{Copyright (C) 2016-2020 -- mrJean1 at Gmail -- All Rights Reserved.
<https://OpenSource.org/licenses/MIT>}

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

@var EPS:    System's M{epsilon} ≈ 2.22e-16 (C{float}).
@var EPS_2:  Half system's M{epsilon} ≈ 1.11e-16 (C{float}).
@var EPS1:   M{1 - EPS} = 0.9999999999999998 (C{float}).
@var EPS1_2: M{1 - EPS_2} = 0.9999999999999999 (C{float}).
@var EPS2:   M{EPS * 2} = 4.440892098501e-16 (C{float}).

@var F_D:   Format degrees as unsigned "deg°" plus suffix (C{str}).
@var F_DM:  Format degrees as unsigned "deg°min′" plus suffix (C{str}).
@var F_DMS: Format degrees as unsigned "deg°min′sec″" plus suffix (C{str}).
@var F_DEG: Format degrees as unsigned "[D]DD" plus suffix without symbol (C{str}).
@var F_MIN: Format degrees as unsigned "[D]DDMM" plus suffix without symbols (C{str}).
@var F_SEC: Format degrees as unsigned "[D]DDMMSS" plus suffix without symbols (C{str}).
@var F__E:  Format degrees as unsigned "%E" plus suffix without symbol (C{str}).
@var F__F:  Format degrees as unsigned "%F" plus suffix without symbol (C{str}).
@var F__G:  Format degrees as unsigned "%G" plus suffix without symbol (C{str}).
@var F_RAD: Convert degrees to radians and format as unsigned "RR" plus suffix (C{str}).

@var F_D_:   Format degrees as signed "-/deg°" without suffix (C{str}).
@var F_DM_:  Format degrees as signed "-/deg°min′" without suffix (C{str}).
@var F_DMS_: Format degrees as signed "-/deg°min′sec″" without suffix (C{str}).
@var F_DEG_: Format degrees as signed "-/[D]DD" without suffix and symbol (C{str}).
@var F_MIN_: Format degrees as signed "-/[D]DDMM" without suffix and symbols (C{str}).
@var F_SEC_: Format degrees as signed "-/[D]DDMMSS" without suffix and symbols (C{str}).
@var F__E_:  Format degrees as signed "-%E" without suffix and symbol (C{str}).
@var F__F_:  Format degrees as signed "-%F" without suffix and symbol (C{str}).
@var F__G_:  Format degrees as signed "-%G" without suffix and symbol (C{str}).
@var F_RAD_: Convert degrees to radians and format as signed "-/RR" without suffix (C{str}).

@var F_D__:   Format degrees as signed "-/+deg°" without suffix (C{str}).
@var F_DM__:  Format degrees as signed "-/+deg°min′" without suffix (C{str}).
@var F_DMS__: Format degrees as signed "-/+deg°min′sec″" without suffix (C{str}).
@var F_DEG__: Format degrees as signed "-/+[D]DD" without suffix and symbol (C{str}).
@var F_MIN__: Format degrees as signed "-/+[D]DDMM" without suffix and symbols (C{str}).
@var F_SEC__: Format degrees as signed "-/+[D]DDMMSS" without suffix and symbols (C{str}).
@var F__E__:  Format degrees as signed "-/+%E" without suffix and symbol (C{str}).
@var F__F__:  Format degrees as signed "-/+%F" without suffix and symbol (C{str}).
@var F__G__:  Format degrees as signed "-/+%G" without suffix and symbol (C{str}).
@var F_RAD__: Convert degrees to radians and format as signed "-/+RR" without suffix (C{str}).

@var INF:    Infinity (C{float}), see function C{isinf}, C{isfinite}.
@var MANTIS: System's M{mantissa bits} ≈53 (C{int}).
@var MAX:    System's M{float max} ≈1.798e+308 (C{float}).
@var MIN:    System's M{float min} ≈2.225e-308 (C{float}).
@var NAN:    Not-A-Number (C{float}), see function C{isnan}.
@var NEG0:   Negative 0.0 (C{float}), see function C{isneg0}.
@var NN:     Empty (C{str}), U{Nomen Nescio<https://Wiktionary.org/wiki/N.N.>}.

@var PI:   Constant M{math.pi} (C{float}).
@var PI2:  Two PI, M{math.pi * 2} (C{float}).
@var PI_2: Half PI, M{math.pi / 2} (C{float}).
@var PI_4: Quarter PI, M{PI / 4} (C{float}).

@var R_M:  Mean (spherical) earth radius (C{meter}).
@var R_MA: Major (equatorial) earth radius (C{meter}) WGS84, EPSG:3785.
@var R_MB: Minor (polar) earth radius (C{meter}) WGS84, EPSG:3785.
@var R_KM: Mean (spherical) earth radius (C{km}, kilometer).
@var R_NM: Mean (spherical) earth radius (C{NM}, nautical miles).
@var R_SM: Mean (spherical) earth radius (C{SM}, statute miles).
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
@var RefFrames:  Registered reference frames (C{enum-like}).
@var Transforms: Registered transforms (C{enum-like}).

@var pygeodesy_abspath: Fully qualified C{pygeodesy} directory name (C{str}).

@var version: Normalized C{PyGeodesy} version (C{str}).

'''

from os.path import abspath, basename, dirname
import sys

# <https://PyInstaller.ReadTheDocs.io/en/stable/runtime-information.html>
_isfrozen         = getattr(sys, 'frozen', False)
pygeodesy_abspath = dirname(abspath(__file__))  # sys._MEIPASS + '/pygeodesy'
_pygeodesy        = __package__ or basename(pygeodesy_abspath)

__version__ = '20.08.04'
# see setup.py for similar logic
version = '.'.join(map(str, map(int, __version__.split('.'))))

if _isfrozen:  # avoid lazy import
    _lazy_import2 = None
else:
    # setting __path__ should ...
    __path__ = [pygeodesy_abspath]
    try:  # ... make this import work, ...
        import pygeodesy.lazily as _
    except ImportError:  # ... if it doesn't, extend
        # sys.path to include this very directory such
        # that all public and private sub-modules can
        # be imported (and checked by PyChecker, etc.)
        sys.path.insert(0, pygeodesy_abspath)  # XXX __path__[0]

    try:
        # lazily requires Python 3.7+, see lazily.__doc__
        from pygeodesy.lazily import LazyImportError, _lazy_import2  # PYCHOK expected
        _, __getattr__ = _lazy_import2(_pygeodesy)  # PYCHOK expected

    except (ImportError, LazyImportError, NotImplementedError):
        _lazy_import2 = None

if not _lazy_import2:  # import and set __all__

    # import all public modules and export as such
    import pygeodesy.azimuthal             as azimuthal              # PYCHOK exported
    import pygeodesy.bases                 as bases                  # PYCHOK DEPRECATED
    import pygeodesy.basics                as basics                 # PYCHOK exported
    import pygeodesy.clipy                 as clipy                  # PYCHOK exported
    import pygeodesy.css                   as css                    # PYCHOK exported
    import pygeodesy.datum                 as datum                  # PYCHOK exported
    import pygeodesy.deprecated            as deprecated             # PYCHOK exported
    import pygeodesy.dms                   as dms                    # PYCHOK exported
    import pygeodesy.ecef                  as ecef                   # PYCHOK exported
    import pygeodesy.elevations            as elevations             # PYCHOK exported
    import pygeodesy.ellipsoidalKarney     as ellipsoidalKarney      # PYCHOK exported
    import pygeodesy.ellipsoidalNvector    as ellipsoidalNvector     # PYCHOK exported
    import pygeodesy.ellipsoidalVincenty   as ellipsoidalVincenty    # PYCHOK exported
    import pygeodesy.elliptic              as elliptic               # PYCHOK exported
    import pygeodesy.epsg                  as epsg                   # PYCHOK exported
    import pygeodesy.etm                   as etm                    # PYCHOK exported
    import pygeodesy.errors                as errors                 # PYCHOK exported
    import pygeodesy.fmath                 as fmath                  # PYCHOK exported
    import pygeodesy.formy                 as formy                  # PYCHOK exported
    import pygeodesy.frechet               as frechet                # PYCHOK exported
    import pygeodesy.gars                  as gars                   # PYCHOK exported
    import pygeodesy.geohash               as geohash                # PYCHOK exported
    import pygeodesy.geoids                as geoids                 # PYCHOK exported
    import pygeodesy.hausdorff             as hausdorff              # PYCHOK exported
    import pygeodesy.heights               as heights                # PYCHOK exported
    import pygeodesy.interns               as interns                # PYCHOK exported
    import pygeodesy.karney                as karney                 # PYCHOK exported
    import pygeodesy.lazily                as lazily                 # PYCHOK exported
    import pygeodesy.lcc                   as lcc                    # PYCHOK exported
    import pygeodesy.mgrs                  as mgrs                   # PYCHOK exported
    import pygeodesy.named                 as named                  # PYCHOK exported
    import pygeodesy.nvector               as nvector                # PYCHOK DEPRECATED
    import pygeodesy.osgr                  as osgr                   # PYCHOK exported
    import pygeodesy.points                as points                 # PYCHOK exported
    import pygeodesy.simplify              as simplify               # PYCHOK exported
    import pygeodesy.sphericalNvector      as sphericalNvector       # PYCHOK exported
    import pygeodesy.sphericalTrigonometry as sphericalTrigonometry  # PYCHOK exported
    import pygeodesy.streprs               as streprs                # PYCHOK exported
    import pygeodesy.trf                   as trf                    # PYCHOK exported
    import pygeodesy.units                 as units                  # PYCHOK exported
    import pygeodesy.ups                   as ups                    # PYCHOK exported
    import pygeodesy.utily                 as utily                  # PYCHOK exported
    import pygeodesy.utm                   as utm                    # PYCHOK exported
    import pygeodesy.utmups                as utmups                 # PYCHOK exported
    import pygeodesy.vector3d              as vector3d               # PYCHOK exported
    import pygeodesy.webmercator           as webmercator            # PYCHOK exported
    import pygeodesy.wgrs                  as wgrs                   # PYCHOK exported

    # lift all public classes, constants, functions, etc. but ONLY
    # from the following modules ... (see also David Beazley's
    # talk <https://DaBeaz.com/modulepackage/index.html>) ... BUT
    # NOT modules ellipsoidal*, epsg, gars, geohash, spherical*,
    # vector and wgrs ... in order keep those as modules ONLY
    from pygeodesy.azimuthal             import *  # PYCHOK __all__
#   from pygeodesy.bases                 import *  # PYCHOK __all__
    from pygeodesy.basics                import *  # PYCHOK __all__
    from pygeodesy.clipy                 import *  # PYCHOK __all__
    from pygeodesy.css                   import *  # PYCHOK __all__
    from pygeodesy.datum                 import *  # PYCHOK __all__
    from pygeodesy.deprecated            import *  # PYCHOK __all__
    from pygeodesy.dms                   import *  # PYCHOK __all__
    from pygeodesy.ecef                  import *  # PYCHOK __all__
    from pygeodesy.elevations            import *  # PYCHOK __all__
#   from pygeodesy.ellipsoidalKarney     import -  # MODULE O_N_L_Y
#   from pygeodesy.ellipsoidalNvector    import -  # MODULE O_N_L_Y
    from pygeodesy.ellipsoidalVincenty   import VincentyError  # PYCHOK exported
    from pygeodesy.elliptic              import *  # PYCHOK __all__
    from pygeodesy.epsg                  import Epsg, EPSGError  # PYCHOK exported
    from pygeodesy.etm                   import *  # PYCHOK __all__
    from pygeodesy.errors                import *  # PYCHOK __all__
    from pygeodesy.fmath                 import *  # PYCHOK __all__
    from pygeodesy.formy                 import *  # PYCHOK __all__
    from pygeodesy.frechet               import *  # PYCHOK __all__
    from pygeodesy.gars                  import Garef, GARSError  # PYCHOK exported
    from pygeodesy.geohash               import Geohash, GeohashError  # PYCHOK exported
    from pygeodesy.geoids                import *  # PYCHOK __all__
    from pygeodesy.hausdorff             import *  # PYCHOK __all__
    from pygeodesy.heights               import *  # PYCHOK __all__
    from pygeodesy.interns               import *  # PYCHOK __all__
#   from pygeodesy.karney                import *  # MODULE O_N_L_Y
    from pygeodesy.lazily                import *  # PYCHOK __all__
    from pygeodesy.lcc                   import *  # PYCHOK __all__
    from pygeodesy.mgrs                  import *  # PYCHOK __all__
    from pygeodesy.named                 import *  # PYCHOK __all__
#   from pygeodesy.nvector               import -  # MODULE O_N_L_Y
    from pygeodesy.osgr                  import *  # PYCHOK __all__
    from pygeodesy.points                import *  # PYCHOK __all__
    from pygeodesy.simplify              import *  # PYCHOK __all__
#   from pygeodesy.sphericalNvector      import -  # MODULE O_N_L_Y
#   from pygeodesy.sphericalTrigonometry import -  # MODULE O_N_L_Y
    from pygeodesy.streprs               import *  # PYCHOK __all__
    from pygeodesy.trf                   import *  # PYCHOK __all__
    from pygeodesy.units                 import *  # PYCHOK __all__
    from pygeodesy.ups                   import *  # PYCHOK __all__
    from pygeodesy.utily                 import *  # PYCHOK __all__
    from pygeodesy.utm                   import *  # PYCHOK __all__
    from pygeodesy.utmups                import *  # PYCHOK __all__
    from pygeodesy.vector3d              import *  # PYCHOK __all__
    from pygeodesy.webmercator           import *  # PYCHOK __all__
    from pygeodesy.wgrs                  import Georef, WGRSError  # PYCHOK exported

    def _all(globalocals):
        from pygeodesy.interns import _dot_  # PYCHOK expected
        # collect all public module and attribute names and check
        # that modules are imported from this package, 'pygeodesy'
        # (but the latter only when not bundled with PyInstaller or
        # Py2Exe, since the file-layout is different.  Courtesy of
        # GilderGeek<https://GitHub.com/mrJean1/PyGeodesy/issues/31>)
        ns = list(lazily._ALL_INIT)
# XXX   ps = () if _isfrozen else set([_pygeodesy] + __name__.split('.'))
        for mod, attrs in lazily._ALL_LAZY.enums():
            if mod not in globalocals:
                raise ImportError('missing %s: %s' % ('module', _dot_(_pygeodesy, mod)))
            ns.append(mod)
            # check that all other public attributes do exist
            if attrs and isinstance(attrs, tuple):
                for a in attrs:
                    if a not in globalocals:
                        raise ImportError('missing %s: %s' % ('attribute', _dot_(mod, a)))
                ns.extend(attrs)
# XXX       if ps:  # check that mod is a _pygeodesy module
# XXX           m = globalocals[mod]  # assert(m.__name__ == mod)
# XXX           f = getattr(m, '__file__', NN)
# XXX           d = dirname(abspath(f)) if f else pygeodesy_abspath
# XXX           p = getattr(m, '__package__', NN) or _pygeodesy
# XXX           if p not in ps or d != pygeodesy_abspath:
# XXX               raise ImportError('foreign module: %s from %r' % (_dot_(p, mod), f or p))
        return tuple(set(ns))  # remove duplicates, only R_M?

    __all__ = _all(globals())  # or locals()

# XXX del ellipsoidalBase, sphericalBase  # PYCHOK expected
del abspath, basename, dirname, _lazy_import2, sys

# **) MIT License
#
# Copyright (C) 2016-2020 -- mrJean1 at Gmail -- All Rights Reserved.
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
