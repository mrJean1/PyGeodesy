
# -*- coding: utf-8 -*-

# The setuptools script to build, test and install a PyGeodesy distribution.

# Tested with 64-bit Python 3.13.0-1, 3.12.0-4, 3.11.0, 3.10.1-5, 3.9.6,
# 3.8.10, 3.7.0, 3.6.4-5, 3.5.3, 2.7,13 and 2.6.9 on macOS 10.12 Sierra,
# 10.1 High Sierra, 10.15.5-7 Catalina, 11.0-2 (10.16) Big Sur, 12.0-4
# (10.16) Monterey, 13.0.1 Ventura and 14.005 Sonoma on Intel (x86_64),
# Intel emulation ("arm64_x86_64") and/or Apple-Silicon M1 native (arm64)
# and with Pythonista 3.1 and 3.2 on iOS 10.3, 11.0, 11.1 and 11.3-4 and
# on iPhone 15.5.

# python setup.py sdist --formats=gztar,bztar,zip  # ztar,tar
# python setup.py bdist_wheel --universal  # XXX
# python setup.py test
# python setup.py install

# <https://packaging.Python.org/key_projects/#setuptools>
# <https://packaging.Python.org/distributing>
# <https://docs.Python.org/2/distutils/sourcedist.html>
# <https://docs.Python.org/3.6/distutils/sourcedist.html>
# <https://setuptools.ReadTheDocs.io/en/latest/setuptools.html#developer-s-guide>
# <https://setuptools.ReadTheDocs.io/en/latest/setuptools.html#test-build-package-and-run-a-unittest-suite>

from os import getenv
from setuptools import setup

__all__ = ()
__version__ = '25.04.24'

_PACKAGE = 'pygeodesy'  # 'PyGeodesy'


def _c(*names):
    return ' :: '.join(names)


def _long_description():
    with open('README.rst', 'rb') as f:
        t = f.read()
        if isinstance(t, bytes):
            t = t.decode('utf-8')
        return t


def _version():
    with open('pygeodesy/__init__.py') as f:
        for t in f.readlines():
            if t.startswith('__version__'):
                v = t.split('=')[-1].strip().strip('\'"')
                c = getenv('PYGEODESY_DIST_RC', '')
                return '.'.join(map(str, map(int, v.split('.')))) + c


_KeyWords = ('AER', 'Albers', 'altitude', 'Andoyer', 'annulus', 'antipode', 'area', 'attitude',
             'Authalic', 'auxiliary', 'azimuth', 'azimuthal', 'azimuth-elevation-range',
             'bearing', 'bank', 'Barsky', 'Barth', 'beta', 'bi-quadratic', 'boolean',
             'cached', 'Cagnoli', 'cartesian', 'Cassini', 'Cassini-Soldner', 'chord',
             'circle-intersections', 'circumcenter', 'circumcircle', 'circumradius',
             'clip', 'Cohen', 'Cohen-Sutherland', 'Collins', 'composite',
             'conformal', 'conic', 'constants', 'contact-triangle',
             'Cook', 'Correia', 'cosines-law', 'coverage', 'curvature', 'cylindrical',
             'datum', 'deprecation', 'deficit', 'development', 'discrete', 'distance', 'Douglas',
             'earth', 'east-north-up', 'eccentricity', 'ECEF', 'elevation', 'ellipsoid',
             'ellipsoidal-latitude-beta', 'ellipsoidal-longitude-omega', 'elliptic',
             'ENU', 'EPSG', 'equal-area', 'equidistant', 'equirectangular', 'ETM', 'ETRF',
             'Euclidean', 'even-odd-rule', 'ExactTM', 'excess', 'Farrell', 'Farrell-Barth',
             'Ferrari-solution', 'Field-Of-View', 'flattening', 'fma', 'fmath',
             'footpoint', 'footprint', 'Forster', 'Forster-Hormann-Popa', 'Forsythe', 'FOV',
             'fractional', 'Frechet', 'Fréchet', 'frustum', 'Fsum', 'fused-multiply-add',
             'GARS', 'geocentric', 'GeoConvert', 'GeodesicExact', 'geodesy', 'geodetic', 'GeodSolve', 'GeodTest',
             'geographiclib', 'Geohash', 'geoid', 'geoidHeight', 'GeoidHeights',
             'georef', 'Girard', 'gnomonic', 'gons', 'grades', 'gradians', 'Greiner', 'Greiner-Hormann',
             'Hartzell', 'Hausdorff', 'Haversine', 'heading', 'hectare', 'height', 'Heikkinen', 'Heron',
             'Hodgman', 'horizon', 'Hormann', 'Hubeny',
             'IDW', 'incenter', 'incirle', 'infix_@_operator', 'inradius', 'intermediate', 'interpolate',
             'intersect', 'intersection', 'intersection3d', 'intersections', 'IntersectTool',
             'Inverse-Distance-Weighting', 'Isometric', 'ITRF',
             'Jacobi', 'Jacobi-Conformal', 'Jarque-Bera', 'Jekel',
             'Karney', 'Krueger', 'Krüger', 'kurtosis',
             'Lambert', 'latitude', 'law-of-cosines', 'least-squares', 'Lesh',
             'L_Huilier', 'LHuilier', 'Liang', 'Liang-Barsky', 'linearize', 'Line-Of-Sight',
             'LocalCartesian', 'local-tangent-plane', 'local-x-y-z', 'longitude', 'LOS', 'loxodrome',
             'lstsq', 'LTP', 'lune', 'LV03', 'LV95',
             'mean', 'memoize', 'memoized', 'Mercator', 'Meeus', 'MGRS',
             'nearest', 'NED', 'Niemeyer', 'non-finite', 'normalize', 'Norrdine', 'north-east-down', 'numpy', 'n-vector', 'Nvector',
             'oblate', 'omega', 'orthographic', 'orthometric-height', 'OSGB', 'OSGR', 'overlap',
             'parallel', 'parallel-of-latitude', 'Parametric', 'path-intersection',
             'perimeter', 'Peucker', 'Pierlot', 'pitch', 'plumb', 'Point-Of-View', 'polar', 'Popa', 'POV',
             'precision-cubic-root', 'precision-hypotenuse', 'precision-powers',
             'precision-running-summation', 'precision-square-root', 'precision-summation',
             'prolate', 'Pseudo-Mercator', _PACKAGE, 'PyInstaller', 'PyPy', 'quartic',
             'radical', 'radii', 'radius', 'Ramer', 'Ramer-Douglas-Peucker', 'Rectifying',
             'Reduced', 'resect', 'resection', 'Rey-Jer', 'Reumann', 'Reumann-Witkam', 'rhumb', 'RhumbSolve',
             'running-linear-regression', 'running-statistics', 'running-stats', 'running-summation',
             'scipy', 'secant', 'semi-perimeter', 'sexagecimal', 'simplify', 'skewness',
             'Snellius', 'Snellius-Pothenot', 'Snyder', 'Soddy', 'Soddy-circles', 'Soldner',
             'sphere', 'sphere-intersections', 'spherical-deficit', 'spherical-excess', 'spherical-triangle',
             'squared-quartic', 'standard-deviation', 'stereographic', 'Sudano', 'surface-area', 'Sutherland',
             'Sutherland-Hodgman', 'tangent-circles', 'Terrestrial-Reference-Frame', 'Thomas', 'Tienstra',
             'tilt', 'TMcoords', 'TMExact', 'toise', 'transverse', 'TransverseMercatorExact', 'TRF',
             'triangle', 'triangulate', 'triaxial', 'triaxial-ellipsoid',
             'trigonometry', 'trilaterate', 'trilaterate-2d', 'trilaterate-3d', 'TwoProduct', 'TwoSum',
             'umbilic-point', 'unit', 'unroll', 'UPS', 'UTM', 'UTM/UPS',
             'variance', 'velocities', 'Veness', 'Vermeille', 'viewing-frustum', 'Vincenty',
             'Visvalingam', 'Visvalingam-Whyatt', 'volume', ' volumetric',
             'Web-Mercator', 'Welford', 'WGRS', 'WGS', 'Whyatt', 'Wildberger', 'Witkam', 'winding-number',
             'XYZ', 'yaw', 'You', 'zenzi-cubic', 'zenzi-quartic')

setup(name=_PACKAGE,
      packages=['pygeodesy', 'pygeodesy.auxilats', 'pygeodesy.deprecated', 'pygeodesy.geodesicx', 'pygeodesy.rhumb'],
      description='Pure Python geodesy tools',
      version=_version(),

      author='Jean M. Brouwers',
      author_email='mrJean1@Gmail.com',
      maintainer='Jean M. Brouwers',
      maintainer_email='mrJean1@Gmail.com',

      license='MIT',
      keywords=' '.join(_KeyWords),
      url='https://GitHub.com/mrJean1/PyGeodesy',

      long_description=_long_description(),

      package_data={_PACKAGE: ['LICENSE']},

      test_suite='test.TestSuite',

      zip_safe=False,

      # <https://PyPI.org/pypi?%3Aaction=list_classifiers>
      classifiers=[_c('Development Status', '5 - Production/Stable'),
                   _c('Environment', 'Console'),
                   _c('Intended Audience', 'Developers'),
                   _c('License', 'OSI Approved', 'MIT License'),
                   _c('Operating System', 'OS Independent'),
                   _c('Programming Language', 'Python'),
#                  _c('Programming Language', 'Python', '2.6'),
                   _c('Programming Language', 'Python', '2.7'),
#                  _c('Programming Language', 'Python', '3.5'),
#                  _c('Programming Language', 'Python', '3.6'),
#                  _c('Programming Language', 'Python', '3.7'),
#                  _c('Programming Language', 'Python', '3.8'),
#                  _c('Programming Language', 'Python', '3.9'),
                   _c('Programming Language', 'Python', '3.10'),
                   _c('Programming Language', 'Python', '3.11'),
                   _c('Programming Language', 'Python', '3.12'),
                   _c('Programming Language', 'Python', '3.13'),
                   _c('Topic', 'Software Development'),
                   _c('Topic', 'Scientific/Engineering', 'GIS'),
      ],  # PYCHOK indent

#     download_url='https://GitHub.com/mrJean1/PyGeodesy',
#     entry_points={},
#     include_package_data=False,
#     install_requires=[],
#     namespace_packages=[],
#     py_modules=[],
)  # PYCHOK indent
