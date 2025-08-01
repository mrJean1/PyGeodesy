
=========
PyGeodesy
=========

A pure Python implementation of ``geodesy`` tools for various ellipsoidal and spherical earth
models using precision exact, elliptic, trigonometric, vector-based, iterative and approximate
methods for geodetic (lat-/longitude), geocentric (ECEF_ cartesian), local (LTP_) and certain
`triaxial ellipsoidal`_ coordinates.

Transcoded in part from `JavaScript originals`_ by *Chris Veness (C) 2005-2024* and from several
`C++ classes`_ by *Charles F.F. Karney (C) 2008-2024* and published under the same `MIT License`_.

There are four modules for ellipsoidal earth models, *ellipsoidalExact*, *-Karney*, *-Vincenty*
and *-Nvector* and two for spherical ones, *sphericalTrigonometry* and *-Nvector*.  Each module
provides a geodetic LatLon_ and a geocentric Cartesian_ class with methods and functions to compute
distance, surface area, perimeter, forward and reverse azimuth, initial and final bearing, intermediate
and nearest points, intersections of geodesic, great circle and rhumb lines, circle intersections and
secants, `3-point resections`_, triangulation, trilateration (by intersection, by overlap and in
3d), conversions and unrolling, among other things.

Also included are modules for conversions to and from `Cassini-Soldner`_, ECEF_ (Earth-Centered,
Earth-Fixed cartesian), UTM_ (Universal Transverse Mercator and Exact_), UPS_ (Universal Polar
Stereographic) and `Web Mercator`_ (Pseudo-Mercator) coordinates, MGRS_ (Military Grid Reference
System, UTM *and* UPS) and OSGR_ (British Ordinance Survery Grid Reference) grid references, TRF_
(Terrestrial Reference Frames) and modules to encode and decode EPSG_, Geohashes_, `Georefs (WGRS)`_
and `Garefs (GARS)`_.

Other modules provide `Albers equal-area`_ projections, equidistant_ and other *azimuthal*
projections, Lambert *conformal conic* projections and positions, functions to clip paths or
polygons of *LatLon* points using the `Cohen-Sutherland`_, `Forster-Hormann-Popa`_,
`Greiner-Hormann`_, `Liang-Barsky`_ and `Sutherland-Hodgman`_ methods or to perform *boolean*
operations between (composite) polygons, functions to simplify_ or linearize a path of *LatLon*
points (or a `numpy array`_), including implementations of the `Ramer-Douglas-Peucker`_,
`Visvalingam-Whyatt`_ and `Reumann-Witkam`_ algorithms and modified versions of the former.

Plus modules and classes to interpolate_ the Height_ of *LatLon* points and Geoid_ models, compute
various Fréchet_ or Hausdorff_ distances or perform *boolean* operations between (composite) polygons.

For further details see the documentation_, the descriptions of `Latitude/Longitude`_, Vincenty_ and
`Vector-based`_ geodesy, the original `JavaScript source`_ or docs_ and *Karney*\'s Python geographiclib_
and C++ `GeographicLib.`_

Installation
============

To install ``pygeodesy``, type ``python[3] -m pip install pygeodesy`` or ``python[3] -m easy_install pygeodesy``
in a terminal or command window.

If the wheel ``pygeodesy-yy.m.d-py2.py3-none-any.whl`` is missing in `PyPI Download files`_, download
the file from `GitHub/dist`_.  Install that with ``python[3] -m pip install <path-to-downloaded-wheel>``
and verify with ``python[3] -m pygeodesy``.

Alternatively, download ``pygeodesy-yy.m.d.tar.gz`` from PyPI_ or GitHub_, ``unzip`` the downloaded file,
``cd`` to directory ``pygeodesy-yy.m.d`` and type ``python[3] setup.py install``.

To run all ``pygeodesy`` tests, type ``python[3] test/run.py`` or type ``python[3] test/unitTestSuite.py``
before or after installation.

Dependencies
============

Installation of *Karney*\'s Python package geographiclib_ is optional, but required for module
``ellipsoidalKarney``, ``azimuthal`` classes ``EquidistantKarney`` and ``GnomonicKarney`` and the
``HeightIDWkarney`` interpolator.

Both numpy_ and scipy_ must be installed for most Geoid_ and Height_ interpolators, except ``GeoidKarney``
and the ``HeigthIDW...`` ones.

Functions and ``LatLon`` methods ``circin6``, ``circum3``, ``circum4_`` and ``soddy4`` and functions
``triaxum5`` and ``trilaterate3d2`` require numpy_ to be installed, modules ``auxilats`` and ``rhumb``
may need numpy_.

Modules ``ellipsoidalGeodSolve`` and ``geodsolve`` and ``azimuthal`` classes ``EquidistantGeodSolve``
and ``GnomonicGeodSolve`` depend on *Karney*\'s C++ utility GeodSolve_ to be executable and set with
env variable ``PYGEODESY_GEODSOLVE`` or with property ``Ellipsoid.geodsolve``.

Class ``Intersectool`` in module ``geodesici`` needs *Karney*\'s C++ utility IntersectTool_ to be
executable and set with env variable ``PYGEODESY_INTERSECTTOOL``.

To compare ``MGRS`` results from modules ``mgrs`` and ``testMgrs`` with *Karney*\'s C++ utility
GeoConvert_, the latter must be executable and set with env variable ``PYGEODESY_GEOCONVERT``.

Module ``rhumbsolve`` needs *Karney*\'s C++ utility RhumbSolve_ to be executable and set with env
variable ``PYGEODESY_RHUMBSOLVE`` or with property ``Ellipsoid.rhumbsolve``.

Documentation
=============

In addition to the ``pygeodesy`` package, the pygeodesy_ `distribution files`_ contain the tests, the
test results (on macOS only) and the complete documentation_ generated by Epydoc_ using command line:
``epydoc --html --no-private --no-source --name=pygeodesy --url=... -v pygeodesy``.

Tests
=====

The tests ran with Python 3.13.5 (with geographiclib_ 2.0), 3.12.7 (with geographiclib_ 2.0, numpy_ 2.1.0,
scipy_ 1.14.1, GeodSolve_ 2.5, IntersectTool_ 2.5 and RhumbSolve_ 2.5), 3.11.5 (with geographiclib_ 2.0,
numpy_ 1.24.2 and scipy_ 1.10.1), Python 3.10.8 (with geographiclib_ 2.0, numpy_ 1.23.3, scipy_ 1.9.1,
GeoConvert_ 2.5, GeodSolve_ 2.5), Python 3.9.6 and Python 2.7.18 (with geographiclib_ 1.50, numpy_ 1.16.6,
scipy_ 1.2.2, GeoConvert_ 2.5, GeodSolve_ 2.5, IntersectTool_ 2.5 and RhumbSolve_ 2.5), all on macOS 15.5
Sequoia in 64-bit.

All tests ran with and without ``lazy import`` for Python 3 and with command line option ``-W default``
and env variable ``PYGEODESY_WARNINGS=on`` for all Python versions.  The results of those tests are
included in the distribution files.

Python 3.13.5, 3.12.7, 3.11.5 and 3.10.8 run on Apple M4 Si (``arm64``), *natively*.  Python 2.7.18 runs
on Intel (``x86_64``) or Intel *emulation* (\"``arm64_x86_64``\", see function `pygeodesy.machine`_).

Test coverage has been measured with coverage_ 7.6.1 using Python 3.13.4, 3.12.7, 3.11.5 and 3.10.8.  The
complete coverage report in HTML and a PDF summary are included in the distribution files.

The tests also ran with Python 3.13.5 (and geographiclib_ 2.0) on `Debian 12`_ in 64-bit only and with
Python 3.12.8 (and geographiclib_ 2.0) on `Windows 2019Server`_ in 64-bit only and with Python 2.7.18
(and with geographiclib_ 1.52) on `Windows 10`_ in 64- and 32-bit.

A single-File and single-Directory application with ``pygeodesy`` has been bundled using PyInstaller_ 3.4
and 64-bit Python 3.7.4 and 3.7.3 on macOS 10.13.6 High Sierra.

Previously, the tests were run with Python 3.13.0-4, 3.12.0-6, 3.11.2-4, 3.10.1-7, 3.9.1, 3.8.7, 3.7.1, 2.7.15,
PyPy_ 7.3.12 (Python 3.10.12), 7.3.1 (Python 3.6.9) and PyPy_ 7.1.1 (Python 2.7.13) (and geographiclib_ 1.52,
numpy_ 1.16.3, 1.16.4, 1.16.6, 1.19.0, 1.19.4, 1.19.5 or 1.22.4 and scipy_ 1.2.1, 1.4.1, 1.5.2 or 1.8.1) on
Ubuntu 16.04, with Python 3.10.0-1, 3.9.0-5, 3.8.0-6, 3.7.2-6, 3.7.0, 3.6.2-5, 3.5.3, 2.7.13-17, 2.7.10
and 2.6.9 (and numpy_ 1.19.0, 1.16.5, 1.16.2, 1.15.2, 1.14.0, 1.13.1, 1.8.0rc1 or 1.6.2 and scipy_ 1.5.0),
PyPy_ 7.3.0 (Python 2.7.13 and 3.6.9), PyPy_ 6.0.0 (Python 2.7.13 and 3.5.3) and `Intel-Python`_ 3.5.3 (and
numpy_ 1.11.3) on macOS 15.0-4 Sequoia, 14.0-6.1 Sonoma, 13.0-5.2 Ventura, 12.1-6 Monterey, 11.0-5.2-6.1 Big
Sur (aka 10.16), 10.15.3, 10.15.5-7 Catalina, 10.14 Mojave, 10.13.6 High Sierra and 10.12 Sierra, MacOS X
10.11 El Capitan and/or MacOS X 10.10 Yosemite, with Pythonista_ 3.2 (with geographiclib 1.50 or 1.49 and
numpy 1.8.0) on iOS 14.4.2, 11.4.1, 12.0-3 on iPad4, iPhone6, iPhone10 and/or iPhone12, with Pythonista_ 3.1
on iOS 10.3.3, 11.0.3, 11.1.2 and 11.3 on iPad4, all in 64-bit only and with 32-bit Python 2.7.14 on Windows
Server 2012R2, Windows 10 Pro and 32-bit Python 2.6.6 on Windows XP SP3.

Notes
=====

All Python source code has been statically checked_ with Ruff_ using Python 3.13.5 and with PyChecker_, PyFlakes_,
PyCodeStyle_ (formerly Pep8) and McCabe_ using Python 2.7.18, both in 64-bit on macOS 15.5 Sequoia only.

For a summary of all *Karney*-based functionality in ``pygeodesy``, see module karney_.

*Last updated: July 25, 2025.*

License
=======

**Copyright (C) 2016-2025 -\- mrJean1 at Gmail -\- All Rights Reserved.**

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and
to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions
of the Software.

THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.

.. image:: https://Img.Shields.io/pypi/pyversions/pygeodesy.svg?label=Python
  :target: https://PyPI.org/project/pygeodesy
.. image:: https://Img.Shields.io/pypi/v/pygeodesy.svg?label=PyPI
  :target: https://PyPI.org/project/pygeodesy
.. image:: https://Img.Shields.io/appveyor/ci/mrJean1/PyGeodesy.svg?branch=master&label=AppVeyor
  :target: https://CI.AppVeyor.com/project/mrJean1/PyGeodesy/branch/master
.. image:: https://Img.Shields.io/cirrus/github/mrJean1/PyGeodesy?branch=master&label=Cirrus
  :target: https://Cirrus-CI.com/github/mrJean1/pygeodesy
.. image:: https://Img.Shields.io/badge/coverage-95%25-brightgreen
  :target: https://GitHub.com/mrJean1/PyGeodesy/blob/master/testcoverage.pdf
.. image:: https://Img.Shields.io/pypi/wheel/pygeodesy.svg
  :target: https://PyPI.org/project/pygeodesy/#files
.. image:: https://Img.Shields.io/pypi/l/pygeodesy.svg
  :target: https://PyPI.org/project/pygeodesy
.. image:: https://img.shields.io/pypi/dm/pygeodesy
  :target: https://PyPI.org/project/pygeodesy

.. _Albers equal-area: https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AlbersEqualArea.html
.. _C++ classes: https://GeographicLib.SourceForge.io/C++/doc/annotated.html
.. _Cartesian: https://mrJean1.GitHub.io/PyGeodesy/docs/pygeodesy-Cartesian-attributes-table.html
.. _Cassini-Soldner: https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1CassiniSoldner.html
.. _checked: https://GitHub.com/ActiveState/code/tree/master/recipes/Python/546532_PyChecker_postprocessor
.. _Cohen-Sutherland: https://WikiPedia.org/wiki/Cohen-Sutherland_algorithm
.. _coverage: https://PyPI.org/project/coverage
.. _Debian 12: https://Cirrus-CI.com/github/mrJean1/pygeodesy/master
.. _distribution files: https://GitHub.com/mrJean1/PyGeodesy/tree/master/dist
.. _docs: https://www.Movable-Type.co.UK/scripts/geodesy/docs
.. _documentation: https://mrJean1.GitHub.io/PyGeodesy
.. _ECEF: https://WikiPedia.org/wiki/ECEF
.. _EPSG: https://EPSG.org
.. _Epydoc: https://PyPI.org/project/epydoc
.. _equidistant: https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1AzimuthalEquidistant.html
.. _Exact: https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1GeodesicExact.html
.. _Forster-Hormann-Popa: https://www.ScienceDirect.com/science/article/pii/S259014861930007X
.. _Fréchet: https://WikiPedia.org/wiki/Frechet_distance
.. _Garefs (GARS): https://WikiPedia.org/wiki/Global_Area_Reference_System
.. _GeoConvert: https://GeographicLib.SourceForge.io/C++/doc/utilities.html
.. _GeodSolve: https://GeographicLib.SourceForge.io/C++/doc/utilities.html
.. _geographiclib: https://PyPI.org/project/geographiclib
.. _GeographicLib.: https://GeographicLib.SourceForge.io/C++/doc/index.html
.. _Geohashes: https://www.Movable-Type.co.UK/scripts/geohash.html
.. _Geoid: https://mrJean1.GitHub.io/PyGeodesy/docs/pygeodesy.geoids-module.html
.. _Georefs (WGRS): https://WikiPedia.org/wiki/World_Geographic_Reference_System
.. _GitHub: https://GitHub.com/mrJean1/PyGeodesy
.. _GitHub/dist: https://GitHub.com/mrJean1/PyGeodesy/tree/master/dist
.. _Greiner-Hormann: http://www.inf.USI.CH/hormann/papers/Greiner.1998.ECO.pdf
.. _Hausdorff: https://WikiPedia.org/wiki/Hausdorff_distance
.. _Height: https://mrJean1.GitHub.io/PyGeodesy/docs/pygeodesy.heights-module.html
.. _Intel-Python: https://software.Intel.com/en-us/distribution-for-python
.. _interpolate: https://docs.SciPy.org/doc/scipy/reference/interpolate.html
.. _IntersectTool: https://GeographicLib.SourceForge.io/C++/doc/utilities.html
.. _JavaScript originals: https://GitHub.com/ChrisVeness/geodesy
.. _JavaScript source: https://GitHub.com/ChrisVeness/geodesy
.. _John P. Snyder: https://pubs.er.USGS.gov/djvu/PP/PP_1395.pdf
.. _karney: https://mrJean1.GitHub.io/PyGeodesy/docs/pygeodesy.karney-module.html
.. _Latitude/Longitude: https://www.Movable-Type.co.UK/scripts/latlong.html
.. _LatLon: https://mrJean1.GitHub.io/PyGeodesy/docs/pygeodesy-LatLon-attributes-table.html
.. _Liang-Barsky: https://www.CS.Helsinki.FI/group/goa/viewing/leikkaus/intro.html
.. _LTP: https://WikiPedia.org/wiki/Local_tangent_plane_coordinates
.. _McCabe: https://PyPI.org/project/mccabe
.. _MGRS: https://GeographicLib.SourceForge.io/C++/doc/classGeographicLib_1_1MGRS.html
.. _MIT License: https://OpenSource.org/licenses/MIT
.. _numpy: https://PyPI.org/project/numpy
.. _numpy array: https://docs.SciPy.org/doc/numpy/reference/generated/numpy.array.html
.. _OSGR: https://www.Movable-Type.co.UK/scripts/latlong-os-gridref.html
.. _3-point resections: https://WikiPedia.org/wiki/Position_resection_and_intersection
.. _PyChecker: https://PyPI.org/project/pychecker
.. _PyCodeStyle: https://PyPI.org/project/pycodestyle
.. _PyFlakes: https://PyPI.org/project/pyflakes
.. _pygeodesy: https://PyPI.org/project/pygeodesy
.. _pygeodesy.machine: https://mrJean1.GitHub.io/PyGeodesy/docs/pygeodesy.internals-module.html#machine
.. _PyInstaller: https://PyPI.org/project/pyinstaller
.. _PyPI: https://PyPI.org/project/pygeodesy
.. _PyPI Download files: https://PyPI.org/project/pygeodesy/#files
.. _PyPy: https://formulae.brew.sh/formula/pypy3.10
.. _Pythonista: https://OMZ-Software.com/pythonista
.. _Ramer-Douglas-Peucker: https://WikiPedia.org/wiki/Ramer-Douglas-Peucker_algorithm
.. _Reumann-Witkam: https://psimpl.SourceForge.net/reumann-witkam.html
.. _RhumbSolve: https://GeographicLib.SourceForge.io/C++/doc/utilities.html
.. _Ruff: https://GitHub.com/astral-sh/ruff
.. _scipy: https://PyPI.org/project/scipy
.. _simplify: https://Bost.Ocks.org/mike/simplify
.. _Sutherland-Hodgman: https://WikiPedia.org/wiki/Sutherland-Hodgman_algorithm
.. _TRF: http://ITRF.ENSG.IGN.FR
.. _triaxial ellipsoidal: https://GeographicLib.SourceForge.io/1.44/triaxial.html
.. _UPS: https://WikiPedia.org/wiki/Universal_polar_stereographic_coordinate_system
.. _UTM: https://www.Movable-Type.co.UK/scripts/latlong-utm-mgrs.html
.. _Vector-based: https://www.Movable-Type.co.UK/scripts/latlong-vectors.html
.. _Vincenty: https://www.Movable-Type.co.UK/scripts/latlong-vincenty.html
.. _Visvalingam-Whyatt: https://hydra.Hull.ac.UK/resources/hull:8338
.. _Web Mercator: https://WikiPedia.org/wiki/Web_Mercator
.. _Windows 10: https://CI.AppVeyor.com/project/mrJean1/pygeodesy
.. _Windows 2019Server: https://CI.AppVeyor.com/project/mrJean1/pygeodesy
