
# -*- coding: utf-8 -*-

# The setuptools script to build, test and install a PyGeodesy
# distribution.

# Tested with 64-bit Python 3.11.0, 3.10.1-5, 3.9.6, 3.8.10,
# 3.7.0, 3.6.4-5, 3.5.3, 2.7,13 and 2.6.9 on macOS 10.12 Sierra,
# 10.1 High Sierra, 10.15.5-7 Catalina, 11.0-2 (10.16) Big Sur,
# 12.0-4 (10.16) Monterey and 13.0.1 Ventura on Intel (x86_64),
# Intel emulation ("arm64_x86_64") and/or Apple-Silicon M1 native
# (arm64) and with Pythonista 3.1 and 3.2 on iOS 10.3, 11.0, 11.1
# and 11.3-4 and on iPhone 15.5.

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

from setuptools import setup

__all__ = ()
__version__ = '23.09.15'


def _c2(*names):
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
                return '.'.join(map(str, map(int, v.split('.'))))


_KeyWords = ('AER', 'Albers', 'altitude', 'Andoyer', 'annulus', 'antipode', 'area', 'attitude',
             'Authalic', 'auxiliary', 'azimuth', 'azimuthal', 'azimuth-elevation-range',
             'bearing', 'bank', 'Barsky', 'Barth', 'beta', 'boolean',
             'cached', 'Cagnoli', 'cartesian', 'Cassini', 'Cassini-Soldner', 'chord',
             'circle-intersections', 'circumcenter', 'circumcircle', 'circumradius',
             'clip', 'Cohen', 'Cohen-Sutherland', 'Collins', 'composite',
             'conformal', 'conic', 'constants', 'contact-triangle',
             'Cook', 'Correia', 'cosines-law', 'coverage', 'curvature', 'cylindrical',
             'datum', 'deprecation', 'deficit', 'development', 'discrete', 'distance', 'Douglas',
             'earth', 'east-north-up', 'eccentricity', 'ECEF', 'elevation', 'ellipsoid',
             'ellipsoidal-latitude-beta', 'ellipsoidal-longitude-omega', 'elliptic',
             'ENU', 'EPSG', 'equal-area', 'equidistant', 'equirectangular', 'ETM', 'ETRF',
             'Euclidean', 'even-odd-rule', 'ExactTM', 'excess',
             'Farrell', 'Farrell-Barth', 'flattening', 'Field-Of-View', 'fmath',
             'footprint', 'Forster', 'Forster-Hormann-Popa', 'Forsythe', 'FOV',
             'fractional', 'Frechet', 'Fréchet', 'frustum', 'Fsum',
             'GARS', 'geocentric', 'GeoConvert', 'GeodesicExact', 'geodesy', 'geodetic', 'GeodSolve', 'GeodTest',
             'geographiclib', 'geohash', 'geoid', 'geoidHeight', 'GeoidHeights',
             'georef', 'Girard', 'gnomonic', 'gons', 'grades', 'gradians', 'Greiner', 'Greiner-Hormann',
             'Hartzell', 'Hausdorff', 'Haversine', 'heading', 'height', 'Heron',
             'Hodgman', 'horizon', 'Hormann', 'Hubeny',
             'IDW', 'incenter', 'incirle', 'infix_@_operator', 'inradius', 'intermediate', 'interpolate',
             'intersect', 'intersection', 'intersection3d', 'intersections',
             'Inverse-Distance-Weighting', 'Isometric', 'ITRF',
             'Jacobi', 'Jacobi-Conformal', 'Jarque-Bera', 'Jekel',
             'Karney', 'Krueger', 'Krüger', 'kurtosis',
             'Lambert', 'latitude', 'law-of-cosines', 'least-squares', 'Lesh',
             'L_Huilier', 'LHuilier', 'Liang', 'Liang-Barsky', 'linearize', 'Line-Of-Sight',
             'LocalCartesian', 'local-tangent-plane', 'local-x-y-z', 'longitude', 'LOS', 'loxodrome',
             'lstsq', 'LTP', 'lune', 'LV03', 'LV95',
             'mean', 'memoize', 'memoized', 'Mercator', 'Meeus', 'MGRS',
             'nearest', 'NED', 'normalize', 'Norrdine', 'north-east-down', 'numpy', 'n-vector', 'Nvector',
             'oblate', 'omega', 'orthographic', 'orthometric-height', 'OSGB', 'OSGR', 'overlap',
             'parallel', 'parallel-of-latitude', 'Parametric', 'path-intersection',
             'perimeter', 'Peucker', 'Pierlot', 'pitch', 'Point-Of-View', 'polar', 'Popa', 'POV',
             'precision-cubic-root', 'precision-hypotenuse', 'precision-powers',
             'precision-running-summation', 'precision-square-root', 'precision-summation',
             'prolate', 'Pseudo-Mercator', 'PyGeodesy', 'PyInstaller', 'PyPy',
             'radical', 'radii', 'radius', 'Ramer', 'Ramer-Douglas-Peucker', 'Rectifying',
             'Reduced', 'resect', 'resection', 'Rey-Jer', 'Reumann', 'Reumann-Witkam', 'rhumb', 'RhumbSolve',
             'running-linear-regression', 'running-statistics', 'running-stats', 'running-summation',
             'scipy', 'secant', 'semi-perimeter', 'sexagecimal', 'simplify', 'skewness',
             'Snellius', 'Snellius-Pothenot', 'Snyder', 'Soddy', 'Soddy-circles', 'Soldner',
             'sphere', 'sphere-intersections', 'spherical-deficit', 'spherical-excess', 'spherical-triangle',
             'standard-deviation', 'stereographic', 'Sudano', 'surface-area', 'Sutherland', 'Sutherland-Hodgman',
             'tangent-circles', 'Terrestrial-Reference-Frame', 'Thomas', 'Tienstra',
             'tilt', 'TMcoords', 'TMExact', 'toise', 'transverse', 'TransverseMercatorExact', 'TRF',
             'triangle', 'triangulate', 'triaxial', 'triaxial-ellipsoid',
             'trigonometry', 'trilaterate', 'trilaterate-2d', 'trilaterate-3d',
             'umbilic-point', 'unit', 'unroll', 'UPS', 'UTM', 'UTM/UPS',
             'variance', 'Veness', 'Vermeille', 'viewing-frustum', 'Vincenty', 'Visvalingam',
             'Visvalingam-Whyatt', 'volume', ' volumetric',
             'Web-Mercator', 'Welford', 'WGRS', 'WGS', 'Whyatt', 'Wildberger', 'Witkam', 'winding-number',
             'XYZ', 'yaw', 'You')

setup(name='PyGeodesy',
      packages=['pygeodesy', 'pygeodesy.auxilats', 'pygeodesy.deprecated', 'pygeodesy.geodesicx'],
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

      package_data={'pygeodesy': ['LICENSE']},

      test_suite='test.TestSuite',

      zip_safe=False,

      # <https://PyPI.org/pypi?%3Aaction=list_classifiers>
      classifiers=[_c2('Development Status', '5 - Production/Stable'),
                   _c2('Environment', 'Console'),
                   _c2('Intended Audience', 'Developers'),
                   _c2('License', 'OSI Approved', 'MIT License'),
                   _c2('Operating System', 'OS Independent'),
                   _c2('Programming Language', 'Python'),
#                  _c2('Programming Language', 'Python', '2.6'),
                   _c2('Programming Language', 'Python', '2.7'),
#                  _c2('Programming Language', 'Python', '3.5'),
#                  _c2('Programming Language', 'Python', '3.6'),
#                  _c2('Programming Language', 'Python', '3.7'),
#                  _c2('Programming Language', 'Python', '3.8'),
                   _c2('Programming Language', 'Python', '3.9'),
                   _c2('Programming Language', 'Python', '3.10'),
                   _c2('Programming Language', 'Python', '3.11'),
                   _c2('Topic', 'Software Development'),
                   _c2('Topic', 'Scientific/Engineering', 'GIS'),
      ],  # PYCHOK indent

#     download_url='https://GitHub.com/mrJean1/PyGeodesy',
#     entry_points={},
#     include_package_data=False,
#     install_requires=[],
#     namespace_packages=[],
#     py_modules=[],
)  # PYCHOK indent
