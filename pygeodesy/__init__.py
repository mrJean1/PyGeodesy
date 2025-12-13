
# -*- coding: utf-8 -*-

u'''A pure Python implementation of geodesy tools for various ellipsoidal and spherical earth models
using precision exact, elliptic, trigonometric, vector-based, iterative and approximate methods for
geodetic (lat-/longitude), geocentric (U{ECEF<https://WikiPedia.org/wiki/ECEF>} cartesian), local (U{LTP
<https://WikiPedia.org/wiki/Local_tangent_plane_coordinates>}) and certain U{triaxial ellipsoidal
<https://GeographicLib.SourceForge.io/1.44/triaxial.html>} coordinates.

Transcoded in part from U{JavaScript originals<https://GitHub.com/ChrisVeness/geodesy>} by I{Chris Veness (C)
2005-2024} and from several U{C++ classes<https://GeographicLib.SourceForge.io/C++/doc/annotated.html>} by I{Charles
F. F. Karney (C) 2008-2024} and published under the same U{MIT License<https://OpenSource.org/licenses/MIT>}**.

There are four modules for ellipsoidal earth models, C{ellipsoidalExact}, C{-Karney}, C{-Vincenty} and C{-Nvector}
and two for spherical ones, C{sphericalTrigonometry} and C{-Nvector}.  Each module provides a geodetic B{C{LatLon}}
and a geocentric B{C{Cartesian}} class with methods and functions to compute distance, surface area, perimeter,
forward and reverse azimuth, initial and final bearing, intermediate and nearest points, intersections of geodesic,
great circle and rhumb lines, circle intersections and secants, U{3-point resections
<https://WikiPedia.org/wiki/Position_resection_and_intersection>}, triangulation, trilateration (by intersection,
by overlap and in 3-D), among other things.

Also included are modules for conversions to and from U{Cassini-Soldner
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1CassiniSoldner.html>},
U{ECEF<https://WikiPedia.org/wiki/ECEF>} (Earth-Centered, Earth-Fixed cartesian), U{UTM
<https://www.Movable-Type.co.UK/scripts/latlong-utm-mgrs.html>} (Universal Transverse Mercator
and U{Exact<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1TransverseMercatorExact.html>}),
U{UPS<https://WikiPedia.org/wiki/Universal_polar_stereographic_coordinate_system>} (Universal Polar
Stereographic) and U{Web Mercator<https://WikiPedia.org/wiki/Web_Mercator>} (Pseudo-Mercator) coordinates,
U{MGRS<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1MGRS.html>} (Military Grid Reference
System, UTM I{and} UPS) and U{OSGR<https://www.Movable-Type.co.UK/scripts/latlong-os-gridref.html>} (British
Ordinance Survery Grid Reference) grid references, U{TRF<http://ITRF.ENSG.IGN.Fr>} (Terrestrial Reference
Frames) and modules to encode and decode U{EPSG<https://EPSG.org>}, U{Geohashes
<https://www.Movable-Type.co.UK/scripts/geohash.html>}, U{Georefs (WGRS)
<https://WikiPedia.org/wiki/World_Geographic_Reference_System>} and U{Garefs (GARS)
<https://WikiPedia.org/wiki/Global_Area_Reference_System>}.

Other modules provide U{Albers equal-area
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AlbersEqualArea.html>} projections, U{equidistant
<https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AzimuthalEquidistant.html>} and other I{azimuthal}
projections, Lambert I{conformal conic} projections and positions, functions to clip paths or polygons of C{LatLon}
points using the U{Cohen-Sutherland<https://WikiPedia.org/wiki/Cohen-Sutherland_algorithm>},
U{Forster-Hormann-Popa<https://www.ScienceDirect.com/science/article/pii/S259014861930007X>},
U{Greiner-Hormann<http://www.inf.USI.CH/hormann/papers/Greiner.1998.ECO.pdf>},
U{Liang-Barsky<https://www.CS.Helsinki.FI/group/goa/viewing/leikkaus/intro.html>} and
U{Sutherland-Hodgman<https://WikiPedia.org/wiki/Sutherland-Hodgman_algorithm>} methods,
functions to U{simplify<https://Bost.Ocks.org/mike/simplify>} or linearize a path of C{LatLon}
points (or a U{NumPy array<https://docs.SciPy.org/doc/numpy/reference/generated/numpy.array.html>}),
including implementations of the U{Ramer-Douglas-Peucker<https://WikiPedia.org/wiki/
Ramer-Douglas-Peucker_algorithm>}, the U{Visvalingam-Whyatt<https://hydra.Hull.ac.UK/
resources/hull:8338>} and the U{Reumann-Witkam<https://psimpl.SourceForge.net/reumann-witkam.html>}
algorithms and modified versions of the former.

Plus modules and classes to U{interpolate<https://docs.SciPy.org/doc/scipy/reference/interpolate.html>} the
L{height<pygeodesy.heights>} of C{LatLon} points and L{Geoid<pygeodesy.geoids>} models, to compute various U{Fréchet
<https://WikiPedia.org/wiki/Frechet_distance>} or U{Hausdorff<https://WikiPedia.org/wiki/Hausdorff_distance>}
distances or to perform I{boolean} operations between (composite) polygons of C{LatLon} points.

For further details see the U{documentation<https://mrJean1.GitHub.io/PyGeodesy>}, the descriptions of
U{Latitude/Longitude<https://www.Movable-Type.co.UK/scripts/latlong.html>}, U{Vincenty
<https://www.Movable-Type.co.UK/scripts/latlong-vincenty.html>} and U{Vector-based
<https://www.Movable-Type.co.UK/scripts/latlong-vectors.html>} geodesy, the original U{JavaScript source
<https://GitHub.com/ChrisVeness/geodesy>} or U{docs<https://www.Movable-Type.co.UK/scripts/geodesy/docs>}
and I{Karney}'s Python U{geographiclib<https://PyPI.org/project/geographiclib>} and C++ U{GeographicLib
<https://GeographicLib.SourceForge.io/C++/doc/index.html>}.

Installation
============

To install C{pygeodesy}, type C{python[3] -m pip install pygeodesy} or C{python[3] -m easy_install
pygeodesy} in a terminal or command window.

If the wheel C{pygeodesy-yy.m.d-py2.py3-none-any.whl} is missing in U{PyPI Download files<https://
PyPI.org/project/pygeodesy/#files>}, download the file from U{GitHub/dist<https://GitHub.com/mrJean1/
PyGeodesy/tree/master/dist>}.  Install that with C{python[3] -m pip install <path-to-downloaded-wheel>}
and verify with C{python[3] -m pygeodesy}.

Alternatively, download C{pygeodesy-yy.m.d.tar.gz} from U{PyPI<https://PyPI.org/project/pygeodesy>}
or U{GitHub<https://GitHub.com/mrJean1/PyGeodesy>}, C{unzip} the downloaded file, C{cd} to
directory C{pygeodesy-yy.m.d} and type C{python[3] setup.py install}.

To run all C{pygeodesy} tests, type C{python[3] test/run.py} or type C{python[3] test/unitTestSuite.py}
before or after installation.

Dependencies
============

Installation of I{Karney}'s Python package U{geographiclib<https://PyPI.org/project/geographiclib>}
is optional, but required for module L{ellipsoidalKarney}, L{azimuthal} classes L{EquidistantKarney}
and L{GnomonicKarney} and the L{HeightIDWkarney} interpolator.

Both U{numpy<https://PyPI.org/project/numpy>} and U{scipy<https://PyPI.org/project/scipy>} must be
installed for most L{Geoid...<pygeodesy.geoids>} and L{Height...<pygeodesy.heights>} interpolators,
except L{GeoidKarney} and the L{HeightIDW...<pygeodesy.heights>} ones.

Functions and C{LatLon} methods L{circin6}, L{circum3}, L{circum4_} and L{soddy4} and functions
L{triaxum5} and L{trilaterate3d2} require U{numpy<https://PyPI.org/project/numpy>} to be installed,
modules L{auxilats} and L{rhumb} may need U{numpy<https://PyPI.org/project/numpy>}.

Modules L{ellipsoidalGeodSolve} and L{geodsolve} and L{azimuthal} classes L{EquidistantGeodSolve}
and L{GnomonicGeodSolve} depend on I{Karney}'s C++ utility U{GeodSolve
<https://GeographicLib.SourceForge.io/C++/doc/GeodSolve.1.html>} to be executable and set with
env variable C{PYGEODESY_GEODSOLVE} or with property L{Ellipsoid.geodsolve}.

Class L{Intersectool} in module L{geodesici} needs I{Karney}'s C++ utility U{IntersectTool
<https://GeographicLib.SourceForge.io/C++/doc/IntersectTool.1.html>} to be executable and set with
env variable C{PYGEODESY_INTERSECTTOOL}.

To compare C{MGRS} results from modules L{mgrs} and C{testMgrs} with I{Karney}'s C++ utility
U{GeoConvert<https://GeographicLib.SourceForge.io/C++/doc/GeoConvert.1.html>}, the latter must
be executable and set with env variable C{PYGEODESY_GEOCONVERT}.

Module L{rhumb.solve} needs I{Karney}'s C++ utility U{RhumbSolve
<https://GeographicLib.SourceForge.io/C++/doc/RhumbSolve.1.html>} to be executable and set with
env variable C{PYGEODESY_RHUMBSOLVE} or with property L{Ellipsoid.rhumbsolve}.

Documentation
=============

In addition to the C{pygeodesy} package, the U{pygeodesy<https://PyPI.org/project/pygeodesy>}
U{distribution files<https://GitHub.com/mrJean1/PyGeodesy/tree/master/dist>} contain the tests,
the test results (on macOS only) and the complete U{documentation<https://mrJean1.GitHub.io/PyGeodesy>}
(generated by U{Epydoc<https://PyPI.org/project/epydoc>} using command line:
C{epydoc --html --no-private --no-source --name=pygeodesy --url=... -v pygeodesy}).

Tests
=====

The tests ran with Python 3.14.2 (with U{geographiclib<https://PyPI.org/project/geographiclib>} 2.1),
Python 3.13.9 (with U{geographiclib<https://PyPI.org/project/geographiclib>} 2.1),
U{numpy<https://PyPI.org/project/numpy>} 2.3.3, U{scipy<https://PyPI.org/project/scipy>} 1.16.2,
U{GeoConvert<https://GeographicLib.SourceForge.io/C++/doc/utilities.html>} 2.5 and
U{GeodSolve<https://GeographicLib.SourceForge.io/C++/doc/utilities.html>} 2.5),
Python 3.12.7 (with U{geographiclib<https://PyPI.org/project/geographiclib>} 2.0,
U{numpy<https://PyPI.org/project/numpy>} 2.1.0, U{scipy<https://PyPI.org/project/scipy>} 1.14.1,
U{GeodSolve<https://GeographicLib.SourceForge.io/C++/doc/utilities.html>} 2.5,
U{IntersectTool<https://GeographicLib.SourceForge.io/C++/doc/utilities.html>} 2.5 and
U{RhumbSolve<https://GeographicLib.SourceForge.io/C++/doc/utilities.html>} 2.5),
Python 3.11.5 (with U{geographiclib<https://PyPI.org/project/geographiclib>} 2.0,
U{numpy<https://PyPI.org/project/numpy>} 1.24.2 and U{scipy<https://PyPI.org/project/scipy>} 1.10.1),
and with Python 2.7.18 (with U{geographiclib<https://PyPI.org/project/geographiclib>} 1.50,
U{numpy<https://PyPI.org/project/numpy>} 1.16.6, U{scipy<https://PyPI.org/project/scipy>} 1.2.2,
U{GeoConvert<https://GeographicLib.SourceForge.io/C++/doc/utilities.html>} 2.5,
U{GeodSolve<https://GeographicLib.SourceForge.io/C++/doc/utilities.html>} 2.5,
U{IntersectTool<https://GeographicLib.SourceForge.io/C++/doc/utilities.html>} 2.5 and
U{RhumbSolve<https://GeographicLib.SourceForge.io/C++/doc/utilities.html>} 2.5), all in 64-bit on
macOS 26.1 Tahoe.

All tests ran with and without C{lazy import} for Python 3 and with command line option C{-W default} and
env variable C{PYGEODESY_WARNINGS=on} for all Python versions.  The results of those tests are included in
the distribution files.

Test coverage has been measured with U{coverage<https://PyPI.org/project/coverage>} 7.10.7 using Python
3.14.2, 3.13.9 and 3.12.7.  The complete coverage report in HTML and a PDF summary are included in the
distribution files.

Python 3.14.2, 3.13.9, 3.12.7 and 3.11.5 run on Apple M4 Si (C{arm64}), I{natively}.  Python 2.7.18 runs
on Intel (C{x86_64}) or Intel I{emulation} ("C{arm64_x86_64}", see function L{machine<pygeodesy.machine>}).

The tests also ran with Python 3.13.9 (and U{geographiclib<https://PyPI.org/project/geographiclib>} 2.1) on
U{Debian 12<https://Cirrus-CI.com/github/mrJean1/PyGeodesy/master>} in 64-bit only, with Python 3.12.8 (and
U{geographiclib<https://PyPI.org/project/geographiclib>} 2.0) on U{Windows 2019Server
<https://CI.AppVeyor.com/project/mrJean1/pygeodesy>} in 64-bit only and with Python 2.7.18 (and U{geographiclib
<https://PyPI.org/project/geographiclib>} 1.52) on U{Windows 10<https://CI.AppVeyor.com/project/mrJean1/pygeodesy>}
in 64- and 32-bit.

A single-File and single-Directory application with C{pygeodesy} has been bundled using U{PyInstaller
<https://PyPI.org/project/pyinstaller>} 3.4 and 64-bit Python 3.7.3 on macOS 10.13.6 High Sierra.

Previously, the tests were run with Python 3.13.0-7, 3.12.0-6, 3.11.2-4, 3.10.1-7, 3.9.6, 3.9.1, 3.8.7, 3.7.1, 2.7.15,
U{PyPy<https://PyPy.org>} 7.3.12 (Python 3.10.12), 7.3.1 (Python 3.6.9) and U{PyPy<https://PyPy.org>} 7.1.1 (Python
2.7.13) (and U{geographiclib <https://PyPI.org/project/geographiclib>} 1.52, U{numpy<https://PyPI.org/project/numpy>}
1.16.3, 1.16.4, 1.16.6, 1.19.0, 1.19.4, 1.19.5 or 1.22.4 and U{scipy<https://PyPI.org/project/scipy>} 1.2.1, 1.4.1,
1.5.2 or 1.8.1) on U{Ubuntu 16.04<https://Travis-CI.com/mrJean1/PyGeodesy>}, with Python 3.10.0-1, 3.9.0-5, 3.8.0-6,
3.7.2-6, 3.7.0, 3.6.2-5, 3.5.3, 2.7.13-17, 2.7.10 and 2.6.9 (and U{numpy<https://PyPI.org/project/numpy>} 1.19.0,
1.16.5, 1.16.2, 1.15.2, 1.14.0, 1.13.1, 1.8.0rc1 or 1.6.2 and U{scipy<https://PyPI.org/project/scipy>} 1.5.0), U{PyPy
<https://PyPy.org>} 7.3.0 (Python 2.7.13 and 3.6.9), U{PyPy<https://PyPy.org>} 6.0.0 (Python 2.7.13 and 3.5.3)
and U{Intel-Python<https://software.Intel.com/en-us/distribution-for-python>} 3.5.3 (and U{numpy
<https://PyPI.org/project/numpy>} 1.11.3) on macOS 15.0-6 Sequoia, 14.0-6.1 Sonoma, 13.0-5.2 Ventura, 12.1-6 Monterey,
11.0-5.2-6.1 Big Sur (aka 10.16), 10.15.3, 10.15.5-7 Catalina, 10.14 Mojave, 10.13.6 High Sierra and 10.12 Sierra,
MacOS X 10.11 El Capitan and/or MacOS X 10.10 Yosemite, with U{Pythonista<https://OMZ-Software.com/pythonista>}3.2
(with geographiclib 1.50 or 1.49 and numpy 1.8.0) on iOS 14.4.2, 11.4.1, 12.0-3 on iPad4, iPhone6, iPhone10 and/or
iPhone12, with U{Pythonista<https://OMZ-Software.com/pythonista>} 3.1 on iOS 10.3.3, 11.0.3, 11.1.2 and 11.3 on iPad4,
all in 64-bit only and with 32-bit Python 2.7.14 on Windows Server 2012R2, Windows 10 Pro and with 32-bit Python 2.6.6
on Windows XP SP3.

Notes
=====

All Python source code has been statically U{checked<https://GitHub.com/ActiveState/code/tree/master/recipes/Python/
546532_PyChecker_postprocessor>} with U{Ruff<https://GitHub.com/astral-sh/ruff>} using Python 3.13.9 and with
U{PyChecker<https://PyPI.org/project/pychecker>}, U{PyFlakes<https://PyPI.org/project/pyflakes>}, U{PyCodeStyle
<https://PyPI.org/project/pycodestyle>} (formerly Pep8) and U{McCabe<https://PyPI.org/project/mccabe>} using Python
2.7.18, both in 64-bit on macOS 26.1 Tahoe.

For a summary of all I{Karney}-based functionality in C{pygeodesy}, see module U{karney
<https://mrJean1.GitHub.io/PyGeodesy/docs/pygeodesy.karney-module.html>}.

In Python 2, symbols L{S_DEG}, L{S_MIN}, L{S_SEC}, L{S_RAD} and L{S_SEP} may be multi-byte, non-ascii
characters and if so, I{not} C{unicode}.

Env variables
=============

The following environment variables are observed by C{pygeodesy}:

 - C{PYGEODESY_EXCEPTION_CHAINING} - see module L{errors<pygeodesy.errors>}.
 - C{PYGEODESY_FMT_FORM} - see module L{dms<pygeodesy.dms>}.
 - C{PYGEODESY_FSUM_F2PRODUCT} - see module L{fsums<pygeodesy.fsums>} and method L{f2product<pygeodesy.Fsum.f2product>}.
 - C{PYGEODESY_FSUM_NONFINITES} - see module L{fsums<pygeodesy.fsums>} and method L{nonfinites<pygeodesy.Fsum.nonfinites>}.
 - C{PYGEODESY_FSUM_RESIDUAL} - see module L{fsums<pygeodesy.fsums>} and method L{RESIDUAL<pygeodesy.Fsum.RESIDUAL>}.
 - C{PYGEODESY_GEOCONVERT} - see module L{mgrs<pygeodesy.mgrs>}.
 - C{PYGEODESY_GEODSOLVE} - see module L{geodsolve<pygeodesy.geodsolve>}.
 - C{PYGEODESY_GEOD3SOLVE} - see module L{geodsolve<pygeodesy.geod3solve>}.
 - C{PYGEODESY_INTERSECTTOOL} - see module L{geodesici<pygeodesy.geodesici>}.
 - C{PYGEODESY_LAZY_IMPORT} - see module L{lazily<pygeodesy.lazily>} and variable L{isLazy<pygeodesy.isLazy>}.
 - C{PYGEODESY_NOTIMPLEMENTED} - C{__special__} methods return C{NotImplemented} if set to "std".
 - C{PYGEODESY_RHUMBSOLVE} - see module L{rhumb.solve<pygeodesy.rhumb.solve>}.
 - C{PYGEODESY_UPS_POLES} - see modules L{ups<pygeodesy.ups>} and L{mgrs<pygeodesy.mgrs>}.

and these to specify standard or I{named} C{repr}esentations:

 - C{PYGEODESY_AZIMUTH_STD_REPR} - see method L{Azimuth<pygeodesy.Azimuth>}C{.__repr__}.
 - C{PYGEODESY_BEARING_STD_REPR} - see method L{Bearing<pygeodesy.Bearing>}C{.__repr__}.
 - C{PYGEODESY_BOOL_STD_REPR} - see method L{Bool<pygeodesy.Bool>}C{.__repr__}.
 - C{PYGEODESY_DEGREES_STD_REPR} - see method L{Degrees<pygeodesy.Degrees>}C{.__repr__}.
 - C{PYGEODESY_EPOCH_STD_REPR} - see method L{Float<pygeodesy.Epoch>}C{.__repr__}.
 - C{PYGEODESY_FLOAT_STD_REPR} - see method L{Float<pygeodesy.Float>}C{.__repr__}.
 - C{PYGEODESY_INT_STD_REPR} - see method L{Int<pygeodesy.Int>}C{.__repr__}.
 - C{PYGEODESY_METER_STD_REPR} - see method L{Meter<pygeodesy.Meter>}C{.__repr__}.
 - C{PYGEODESY_RADIANS_STD_REPR} - see method L{Radians<pygeodesy.Radians>}C{.__repr__}.
 - C{PYGEODESY_STR_STD_REPR} - see method L{Str<pygeodesy.Str>}C{.__repr__}.

plus during development:

 - C{PYGEODESY_FOR_DOCS} - for extended documentation by C{epydoc}.
 - C{PYGEODESY_GEOGRAPHICLIB} - see module L{karney<pygeodesy.karney>}.
 - C{PYGEODESY_WARNINGS} - see module L{props<pygeodesy.props>} and function L{DeprecationWarnings<pygeodesy.DeprecationWarnings>}.
 - C{PYGEODESY_XPACKAGES} - see module L{basics<pygeodesy.basics>}.
 - C{PYTHONDEVMODE} - see modules L{errors<pygeodesy.errors>} and L{props<pygeodesy.props>}.

and:

 - C{PYGEODESY_INIT__ALL__} - Set env variable C{PYGEODESY_INIT__ALL__} to anything other than C{"__all__"} to avoid
   importing all C{pygeodesy} modules unnecessarily (in Python 2 or with C{PYGEODESY_LAZY_IMPORT} turned off in Python
   3).  However, to import a C{pygeodesy} item, the item name must be qualified with the C{module} name, for example
   C{ from pygeodesy.ellipsoidalExact import LatLon } or C{ from pygeodesy.deprecated import collins }

License
=======

**) U{Copyright (C) 2016-2026 -- mrJean1 at Gmail -- All Rights Reserved.<https://OpenSource.org/licenses/MIT>}

C{Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:}

C{The above copyright notice and this permission notice shall be included in all copies or substantial
portions of the Software.}

C{THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.}

@newfield example: Example, Examples

@var EPS:    System's M{epsilon} ≈ 2.22044604925e-16 (C{float}).
@var EPS0:   M{EPS**2}  ≈ 4.9e-32 for near-zero comparison
@var EPS02:  M{EPS**4}  ≈ 2.4e-63 for near-zero squared comparison
@var EPS1:   M{1 - EPS} ≈ 0.9999999999999998 (C{float}).
@var EPS2:   M{EPS * 2} ≈ 4.440892098501e-16 (C{float}).
@var EPS_2:  M{EPS / 2} ≈ 1.110223024625e-16 (C{float}).
@var EPS4:   M{EPS * 4} ≈ 8.881784197001e-16 (C{float}).
@var EPS8:   M{EPS * 8} ≈ 1.776356839400e-15 (C{float}).

@var F_D:   Format degrees as unsigned "deg°" with symbol, plus compass point suffix C{N, S, E} or C{W} (C{str}).
@var F_DM:  Format degrees as unsigned "deg°min′" with symbols, plus suffix (C{str}).
@var F_DMS: Format degrees as unsigned "deg°min′sec″" with symbols, plus suffix (C{str}).
@var F_DEG: Format degrees as unsigned "[D]DD" I{without} symbol, plus suffix (C{str}).
@var F_MIN: Format degrees as unsigned "[D]DDMM" I{without} symbols, plus suffix (C{str}).
@var F_SEC: Format degrees as unsigned "[D]DDMMSS" I{without} symbols, plus suffix (C{str}).
@var F_D60: Format degrees as unsigned "[D]DD.MMSS" C{sexagecimal} I{without} symbols, plus suffix (C{str}).
@var F__E:  Format degrees as unsigned "%E" I{without} symbols, plus suffix (C{str}).
@var F__F:  Format degrees as unsigned "%F" I{without} symbols, plus suffix (C{str}).
@var F__G:  Format degrees as unsigned "%G" I{without} symbols, plus suffix (C{str}).
@var F_RAD: Convert degrees to radians and format as unsigned "RR" with symbol, plus suffix (C{str}).

@var F_D_:   Format degrees as signed "-/deg°" with symbol, I{without} suffix (C{str}).
@var F_DM_:  Format degrees as signed "-/deg°min′" with symbols, I{without} suffix (C{str}).
@var F_DMS_: Format degrees as signed "-/deg°min′sec″" with symbols, I{without} suffix (C{str}).
@var F_DEG_: Format degrees as signed "-/[D]DD" I{without} symbol, I{without} suffix (C{str}).
@var F_MIN_: Format degrees as signed "-/[D]DDMM" I{without} symbols, I{without} suffix (C{str}).
@var F_SEC_: Format degrees as signed "-/[D]DDMMSS" I{without} symbols, I{without} suffix (C{str}).
@var F_D60_: Format degrees as signed "-/[D]DD.MMSS" C{sexagecimal} I{without} symbols, I{without} suffix (C{str}).
@var F__E_:  Format degrees as signed "-/%E" I{without} symbols, I{without} suffix (C{str}).
@var F__F_:  Format degrees as signed "-/%F" I{without} symbols, I{without} suffix (C{str}).
@var F__G_:  Format degrees as signed "-/%G" I{without} symbols, I{without} suffix (C{str}).
@var F_RAD_: Convert degrees to radians and format as signed "-/RR" I{without} symbol, I{without} suffix (C{str}).

@var F_D__:   Format degrees as signed "-/+deg°" with symbol, I{without} suffix (C{str}).
@var F_DM__:  Format degrees as signed "-/+deg°min′" with symbols, I{without} suffix (C{str}).
@var F_DMS__: Format degrees as signed "-/+deg°min′sec″" with symbols, I{without} suffix (C{str}).
@var F_DEG__: Format degrees as signed "-/+[D]DD" I{without} symbol, I{without} suffix (C{str}).
@var F_MIN__: Format degrees as signed "-/+[D]DDMM" I{without} symbols, without suffix (C{str}).
@var F_SEC__: Format degrees as signed "-/+[D]DDMMSS" I{without} symbols, I{without} suffix (C{str}).
@var F_D60__: Format degrees as signed "-/+[D]DD.MMSS" C{sexagecimal} I{without} symbols, I{without} suffix (C{str}).
@var F__E__:  Format degrees as signed "-/+%E" I{without} symbols, I{without} suffix (C{str}).
@var F__F__:  Format degrees as signed "-/+%F" I{without} symbols, I{without} suffix (C{str}).
@var F__G__:  Format degrees as signed "-/+%G" I{without} symbols, I{without} suffix (C{str}).
@var F_RAD__: Convert degrees to radians and format as signed "-/+RR" I{without} symbol, I{without} suffix (C{str}).

@var DIG:      System's M{float decimal digits} = 15 (C{int}).
@var INF:      Infinity (C{float}), see functions L{isinf<pygeodesy.isinf>} and L{isfinite<pygeodesy.isfinite>} and C{NINF}.
@var INT0:     C{int(0)}, missing Z-components, C{if B{z}=B{INT0}}, see functions L{isint0<pygeodesy.isint0>}, L{meeus2<pygeodesy.meeus2>}
@var MANT_DIG: System's M{float mantissa bits} = 53 (C{int}).
@var MAX:      System's M{float max} ≈ 1.798e+308 (C{float}).
@var MIN:      System's M{float min} ≈ 2.225e-308 (C{float}).
@var NAN:      Not-A-Number (C{float}), see function L{isnan<pygeodesy.isnan>}.
@var NEG0:     Negative 0.0 (C{float}), see function L{isneg0<pygeodesy.isneg0>}.
@var NINF:     Negative infinity (C{float}), see function L{isninf<pygeodesy.isninf>} and C{INF}.
@var NN:       Empty (C{str}), U{I{Nomen Nescio}<https://Wiktionary.org/wiki/N.N.>}.

@var OVERFLOW: Object representing C{overflow} (1 / L{EPS0<pygeodesy.constants.EPS0>} = 4.9e+32).

@var PI:    Constant M{math.pi} (C{float}).
@var PI2:   Two PI, M{PI * 2}, aka I{Tau} (C{float}).
@var PI_2:  Half PI, M{PI / 2} (C{float}).
@var PI3:   Three PI, M{PI * 3} (C{float}).
@var PI3_2: One and a half PI, M{PI * 3 / 2} (C{float}).
@var PI_3:  One third PI, M{PI / 3} (C{float}).
@var PI4:   Four PI, M{PI * 4} (C{float}).
@var PI_4:  Quarter PI, M{PI / 4} (C{float}).
@var PI_6:  One sixth PI, M{PI / 6} (C{float}).

@var R_MA: Equatorial earth radius (C{meter}), WGS84, EPSG:3785.
@var R_MB: Polar earth radius (C{meter}), WGS84, EPSG:3785.
@var R_M:  Mean (spherical) earth radius (C{meter}).
@var R_KM: Mean (spherical) earth radius (C{Km}, kilometer).
@var R_NM: Mean (spherical) earth radius (C{NM}, nautical miles).
@var R_SM: Mean (spherical) earth radius (C{SM}, statute miles).
@var R_FM: Former FAI-Sphere earth radius (C{meter}).
@var R_GM: Average earth radius, distance to geoid surface (C{meter})
@var R_QM: Earth' (triaxial) quadratic mean radius (C{meter})
@var R_VM: Aviation/Navigation earth radius (C{meter}).

@var S_DEG: Degrees symbol, default C{"°"}
@var S_MIN: Minutes symbol, default C{"′"} aka I{PRIME}
@var S_SEC: Seconds symbol, default C{"″"} aka I{DOUBLE_PRIME}
@var S_RAD: Radians symbol, default C{""} aka L{NN<pygeodesy.NN>}
@var S_DMS: If C{True}, include, otherwise cancel all DMS symbols, default C{True}.
@var S_SEP: Separator between C{deg°|min′|sec″|suffix}, default C{""} aka L{NN<pygeodesy.NN>}

@var Conics:     Registered, predefined conics (C{enum-like}).
@var Datums:     Registered, predefined datums (C{enum-like}).
@var Ellipsoids: Registered, predefined ellipsoids (C{enum-like}).
@var RefFrames:  Registered, predefined reference frames (C{enum-like}).
@var Transforms: Registered, predefined transforms (C{enum-like}).
@var Triaxials:  Registered, predefined triaxial ellipsoids (C{enum-like}).

@var isLazy: Lazy import setting (C{int} 0, 1, 2 or 3+) from C{env} variable C{PYGEODESY_LAZY_IMPORT}, or C{None} if C{lazy import} is not supported or not enabled, or C{False} if initializing C{lazy import} failed.

@var pygeodesy_abspath: Fully qualified C{pygeodesy} directory name (C{str}).

@var version: Normalized C{pygeodesy} version (C{str}).
'''

import os.path as _os_path
import sys as _sys

_init__all__      =  True
# <https://PyInstaller.ReadTheDocs.io/en/stable/runtime-information.html>
_isfrozen         =  getattr(_sys, 'frozen', False)
pygeodesy_abspath = _os_path.dirname(_os_path.abspath(__file__))  # _sys._MEIPASS + '/pygeodesy'
_pygeodesy_       = __package__ or   _os_path.basename(pygeodesy_abspath)

if _isfrozen:  # avoid lazy import
    _lazy_import2 = None
else:
    # setting __path__ should ...
    __path__ = [pygeodesy_abspath]
    try:  # ... make this import work, ...
        import pygeodesy.lazily as _
    except ImportError:  # ... if it doesn't, extend _sys.path to include
        # this very directory such that all public and private sub-modules
        # can be imported (by epydoc, checked by PyChecker, etc.)
        _sys.path.insert(0, pygeodesy_abspath)  # XXX __path__[0]

    try:  # lazily requires Python 3.7+, see lazily.__doc__
        from pygeodesy.lazily import _init__all__, _lazy_import2  # PYCHOK expected
        _, __getattr__ = _lazy_import2(_pygeodesy_)  # PYCHOK expected

    except (ImportError, NotImplementedError):  # LazyImportError
        _lazy_import2 = None

if _init__all__ and not _lazy_import2:  # import and set __all__

    # import all public modules and export as such
    import pygeodesy.albers                as albers                 # noqa: F401
    import pygeodesy.angles                as angles                 # noqa: F401
    import pygeodesy.auxilats              as auxilats               # noqa: F401
    import pygeodesy.azimuthal             as azimuthal              # noqa: F401
    import pygeodesy.basics                as basics                 # noqa: F401
    import pygeodesy.booleans              as booleans               # noqa: F401
    import pygeodesy.cartesianBase         as cartesianBase          # noqa: F401 INTERNAL
    import pygeodesy.clipy                 as clipy                  # noqa: F401
    import pygeodesy.constants             as constants              # noqa: F401
    import pygeodesy.css                   as css                    # noqa: F401
    import pygeodesy.datums                as datums                 # noqa: F401
    import pygeodesy.deprecated            as deprecated             # noqa: F401
    import pygeodesy.dms                   as dms                    # noqa: F401
    import pygeodesy.ecef                  as ecef                   # noqa: F401
    import pygeodesy.ecefLocals            as ecefLocals             # noqa: F401
    import pygeodesy.elevations            as elevations             # noqa: F401
    import pygeodesy.ellipsoidalBase       as ellipsoidalBase        # noqa: F401 INTERNAL
    import pygeodesy.ellipsoidalBaseDI     as ellipsoidalBaseDI      # noqa: F401 INTERNAL
    import pygeodesy.ellipsoidalExact      as ellipsoidalExact       # noqa: F401
    import pygeodesy.ellipsoidalGeodSolve  as ellipsoidalGeodSolve   # noqa: F401
    import pygeodesy.ellipsoidalKarney     as ellipsoidalKarney      # noqa: F401
    import pygeodesy.ellipsoidalNvector    as ellipsoidalNvector     # noqa: F401
    import pygeodesy.ellipsoidalVincenty   as ellipsoidalVincenty    # noqa: F401
    import pygeodesy.ellipsoids            as ellipsoids             # noqa: F401
    import pygeodesy.elliptic              as elliptic               # noqa: F401
    import pygeodesy.epsg                  as epsg                   # noqa: F401
    import pygeodesy.etm                   as etm                    # noqa: F401
    import pygeodesy.errors                as errors                 # noqa: F401
    import pygeodesy.fmath                 as fmath                  # noqa: F401
    import pygeodesy.formy                 as formy                  # noqa: F401
    import pygeodesy.frechet               as frechet                # noqa: F401
    import pygeodesy.fstats                as fstats                 # noqa: F401
    import pygeodesy.fsums                 as fsums                  # noqa: F401
    import pygeodesy.gars                  as gars                   # noqa: F401
    import pygeodesy.geodesici             as geodesici              # noqa: F401
    import pygeodesy.geodesicw             as geodesicw              # noqa: F401
    import pygeodesy.geodesicx             as geodesicx              # noqa: F401
    import pygeodesy.geodsolve             as geodsolve              # noqa: F401
    import pygeodesy.geod3solve            as geod3solve             # noqa: F401
    import pygeodesy.geohash               as geohash                # noqa: F401
    import pygeodesy.geoids                as geoids                 # noqa: F401
    import pygeodesy.hausdorff             as hausdorff              # noqa: F401
    import pygeodesy.heights               as heights                # noqa: F401
    import pygeodesy.internals             as internals              # noqa: F401
    import pygeodesy.interns               as interns                # noqa: F401
    import pygeodesy.iters                 as iters                  # noqa: F401
    import pygeodesy.karney                as karney                 # noqa: F401
    import pygeodesy.ktm                   as ktm                    # noqa: F401
    import pygeodesy.latlonBase            as latlonBase             # noqa: F401 INTERNAL
    import pygeodesy.lazily                as lazily                 # noqa: F401
    import pygeodesy.lcc                   as lcc                    # noqa: F401
    import pygeodesy.ltp                   as ltp                    # noqa: F401
    import pygeodesy.ltpTuples             as ltpTuples              # noqa: F401
    import pygeodesy.mgrs                  as mgrs                   # noqa: F401
    import pygeodesy.named                 as named                  # noqa: F401
    import pygeodesy.namedTuples           as namedTuples            # noqa: F401
    import pygeodesy.nvectorBase           as nvectorBase            # noqa: F401 INTERNAL
    import pygeodesy.osgr                  as osgr                   # noqa: F401
    import pygeodesy.points                as points                 # noqa: F401
    import pygeodesy.props                 as props                  # noqa: F401
    import pygeodesy.resections            as resections             # noqa: F401
    import pygeodesy.rhumb                 as rhumb                  # noqa: F401
    import pygeodesy.simplify              as simplify               # noqa: F401
    import pygeodesy.sphericalBase         as sphericalBase          # noqa: F401 INTERNAL
    import pygeodesy.sphericalNvector      as sphericalNvector       # noqa: F401
    import pygeodesy.sphericalTrigonometry as sphericalTrigonometry  # noqa: F401
    import pygeodesy.solveBase             as solveBase              # noqa: F401
    import pygeodesy.streprs               as streprs                # noqa: F401
    import pygeodesy.trf                   as trf                    # noqa: F401
    import pygeodesy.triaxials             as triaxials              # noqa: F401
    import pygeodesy.units                 as units                  # noqa: F401
    import pygeodesy.unitsBase             as unitsBase              # noqa: F401 INTERNAL
    import pygeodesy.ups                   as ups                    # noqa: F401
    import pygeodesy.utily                 as utily                  # noqa: F401
    import pygeodesy.utm                   as utm                    # noqa: F401
    import pygeodesy.utmups                as utmups                 # noqa: F401
    import pygeodesy.utmupsBase            as utmupsBase             # noqa: F401 INTERNAL
    import pygeodesy.vector2d              as vector2d               # noqa: F401
    import pygeodesy.vector3d              as vector3d               # noqa: F401
    import pygeodesy.vector3dBase          as vector3dBase           # noqa: F401 INTERNAL
    import pygeodesy.webmercator           as webmercator            # noqa: F401
    import pygeodesy.wgrs                  as wgrs                   # noqa: F401

    # lift all public classes, constants, functions, etc. but ONLY
    # from the following modules ... (see also David Beazley's
    # talk <https://DaBeaz.com/modulepackage/index.html>) ... BUT
    # NOT modules ellipsoidal*, epsg, gars, geohash, spherical*,
    # vector and wgrs ... in order keep those as modules ONLY
    from pygeodesy.albers                import *  # noqa: F403
    from pygeodesy.angles                import *  # noqa: F403
    from pygeodesy.azimuthal             import *  # noqa: F403
#   from pygeodesy.auxilats              import *  # noqa: F403
    from pygeodesy.basics                import *  # noqa: F403
    from pygeodesy.booleans              import *  # noqa: F403
    from pygeodesy.cartesianBase         import *  # noqa: F403 INTERNAL
    from pygeodesy.clipy                 import *  # noqa: F403
    from pygeodesy.constants             import *  # noqa: F403
    from pygeodesy.css                   import *  # noqa: F403
    from pygeodesy.datums                import *  # noqa: F403
    from pygeodesy.deprecated            import *  # noqa: F403 DEPRECATED
    from pygeodesy.dms                   import *  # noqa: F403
    from pygeodesy.ecef                  import *  # noqa: F403
#   from pygeodesy.ecefLocals            import *  # noqa: F403
    from pygeodesy.elevations            import *  # noqa: F403
#   from pygeodesy.ellipsoidalBase       import *  # noqa: F403 INTERNAL
#   from pygeodesy.ellipsoidalBaseDI     import *  # noqa: F403 INTERNAL
#   from pygeodesy.ellipsoidalExact      import *  # noqa: F403
#   from pygeodesy.ellipsoidalGeodSolve  import *  # noqa: F403
#   from pygeodesy.ellipsoidalKarney     import *  # noqa: F403
#   from pygeodesy.ellipsoidalNvector    import *  # noqa: F403
    from pygeodesy.ellipsoidalVincenty   import *  # noqa: F403
    from pygeodesy.ellipsoids            import *  # noqa: F403
    from pygeodesy.elliptic              import *  # noqa: F403
    from pygeodesy.epsg                  import *  # noqa: F403
    from pygeodesy.etm                   import *  # noqa: F403
    from pygeodesy.errors                import *  # noqa: F403
    from pygeodesy.fmath                 import *  # noqa: F403
    from pygeodesy.formy                 import *  # noqa: F403
    from pygeodesy.frechet               import *  # noqa: F403
    from pygeodesy.fstats                import *  # noqa: F403
    from pygeodesy.fsums                 import *  # noqa: F403
    from pygeodesy.gars                  import *  # noqa: F403
    from pygeodesy.geodesici             import *  # noqa: F403
    from pygeodesy.geodesicw             import *  # noqa: F403
    from pygeodesy.geodesicx             import *  # noqa: F403
    from pygeodesy.geodsolve             import *  # noqa: F403
    from pygeodesy.geod3solve            import *  # noqa: F403
    from pygeodesy.geohash               import *  # noqa: F403
    from pygeodesy.geoids                import *  # noqa: F403
    from pygeodesy.hausdorff             import *  # noqa: F403
    from pygeodesy.heights               import *  # noqa: F403
    from pygeodesy.internals             import *  # noqa: F403
    from pygeodesy.interns               import *  # noqa: F403
    from pygeodesy.iters                 import *  # noqa: F403
    from pygeodesy.karney                import *  # noqa: F403
    from pygeodesy.ktm                   import *  # noqa: F403
    from pygeodesy.latlonBase            import *  # noqa: F403 INTERNAL
    from pygeodesy.lazily                import *  # noqa: F403
    from pygeodesy.lcc                   import *  # noqa: F403
    from pygeodesy.ltp                   import *  # noqa: F403
    from pygeodesy.ltpTuples             import *  # noqa: F403
    from pygeodesy.mgrs                  import *  # noqa: F403
    from pygeodesy.named                 import *  # noqa: F403
    from pygeodesy.namedTuples           import *  # noqa: F403
    from pygeodesy.nvectorBase           import *  # noqa: F403 INTERNAL
    from pygeodesy.osgr                  import *  # noqa: F403
    from pygeodesy.points                import *  # noqa: F403
    from pygeodesy.props                 import *  # noqa: F403
    from pygeodesy.resections            import *  # noqa: F403
    from pygeodesy.rhumb                 import *  # noqa: F403
    from pygeodesy.simplify              import *  # noqa: F403
#   from pygeodesy.sphericalBase         import *  # noqa: F403 INTERNAL
#   from pygeodesy.sphericalNvector      import *  # noqa: F403
#   from pygeodesy.sphericalTrigonometry import *  # noqa: F403
#   from pygeodesy.solveBase             import *  # noqa: F403 INTERNAL
    from pygeodesy.streprs               import *  # noqa: F403
    from pygeodesy.trf                   import *  # noqa: F403
    from pygeodesy.triaxials             import *  # noqa: F403
#   from pygeodesy.triaxials.bases       import *  # noqa: F403
#   from pygeodesy.triaxials.conformals  import *  # noqa: F403
#   from pygeodesy.triaxials.spheres     import *  # noqa: F403
#   from pygeodesy.triaxials.triaxial3   import *  # noqa: F403
#   from pygeodesy.triaxials.triaxial5   import *  # noqa: F403
    from pygeodesy.units                 import *  # noqa: F403
    from pygeodesy.unitsBase             import *  # noqa: F403 Float, ...
    from pygeodesy.ups                   import *  # noqa: F403
    from pygeodesy.utily                 import *  # noqa: F403
    from pygeodesy.utm                   import *  # noqa: F403
    from pygeodesy.utmups                import *  # noqa: F403
#   from pygeodesy.utmupsBase            import *  # noqa: F403 INTERNAL
    from pygeodesy.vector2d              import *  # noqa: F403
    from pygeodesy.vector3d              import *  # noqa: F403
    from pygeodesy.vector3dBase          import *  # noqa: F403 INTERNAL
    from pygeodesy.webmercator           import *  # noqa: F403
    from pygeodesy.wgrs                  import *  # noqa: F403

    def _all(globalocals):
        from pygeodesy.internals import _headof,  _DOT_  # PYCHOK expected
        from pygeodesy.interns import NN as _NN, _attribute_, _COMMASPACE_, \
                                     _module_, _s_  # PYCHOK expected
        from pygeodesy.streprs import Fmt as _Fmt  # PYCHOK expected
        # collect all public module and attribute names and check
        # that modules are imported from this package, 'pygeodesy'
        # (but the latter only when not bundled with PyInstaller or
        # Py2Exe, since the file-layout is different.  Courtesy of
        # GilderGeek<https://GitHub.com/mrJean1/PyGeodesy/issues/31>)
        ns = list(lazily._ALL_INIT)
# XXX   ps = () if _isfrozen else set([_pygeodesy_] + __name__.split(_DOT_))
        for mod, attrs in lazily._all_enums():
            mod = _headof(_headof(mod))
            if mod not in globalocals:
                t = _DOT_(_pygeodesy_, mod)
                raise ImportError('missing %s%s: %s' % (_module_, _NN, t))
            ns.append(mod)
            # check that all other public attrs do exist
            if attrs and isinstance(attrs, tuple):
                attrs = lazily._ALL_ATTRS(attrs)
                t = tuple(a for a in attrs if a not in globalocals)
                if t:
                    s = _Fmt.SQUARE(_s_, len(t)) if len(t) > 1 else _NN
                    t = _COMMASPACE_.join(_DOT_(_pygeodesy_, mod, a) for a in t)
                    raise ImportError('missing %s%s: %s' % (_attribute_, s, t))
                ns.extend(attrs)
# XXX       if ps:  # check that mod is a _pygeodesy_ module
# XXX           m = globalocals[mod]  # assert(typename(m) == mod)
# XXX           f = getattr(m, _dunder_file_, _NN)
# XXX           d = _os_path.dirname(_os_path.abspath(f)) if f else pygeodesy_abspath
# XXX           p = getattr(m, _dunder_package_, _NN) or _pygeodesy_
# XXX           if p not in ps or d != pygeodesy_abspath:
# XXX               t = '%s from %r' % (_DOT_(p, mod), f or p)
# XXX               raise ImportError('foreign %s: %s' % (_module_, t))
        return tuple(set(ns))  # remove duplicates

    __all__ = _all(globals())  # or locals()
else:
    __all__ = ()
    _init__all__ = False

from pygeodesy.internals import _version2,  _DOT_  # noqa: E402
# from pygeodesy.interns import _DOT_  # from .internals
__version__ = '25.12.12'
# see setup.py for similar logic
version     = _DOT_(*_version2(__version__, n=3))

del _DOT_, _lazy_import2, _os_path, _sys, _version2

# **) MIT License
#
# Copyright (C) 2016-2026 -- mrJean1 at Gmail -- All Rights Reserved.
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
