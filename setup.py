
# -*- coding: utf-8 -*-

# The setuptools script to build, install and test a PyGeodesy distribution.

# Tested with 64-bit Python 2.7.13-18, 3.6.1-2, 3.7.0-6 and 3.8.0-3
# (using setuptools 28.8.0), but only on macOS 10.12.3-6 Sierra,
# 10.13.0-6 High Sierra and 10.15.5 Catalina

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
# <https://ZetCode.com/articles/packageinpython>

from setuptools import setup

__all__ = ()
__version__ = '20.08.01'


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


_KeyWords = ('altitude', 'Andoyer', 'antipode', 'area', 'azimuth', 'azimuthal', 'bearing',
             'cartesian', 'Cassini', 'Cassini-Soldner', 'clip', 'Cohen',
             'Cohen-Sutherland', 'conic', 'cosines-law', 'coverage', 'curvature',
             'datum', 'development', 'discrete', 'distance', 'Douglas',
             'earth', 'ECEF', 'elevation', 'ellipsoid', 'elliptic', 'EPSG',
             'equal-area', 'equidistant', 'equirectangular', 'ETM', 'ETRF', 'Euclidean', 'ExactTM',
             'fmath', 'Forsythe', 'fractional', 'Fréchet',
             'GARS', 'geocentric', 'geodesy', 'geodetic', 'GeodTest', 'geographiclib',
             'geohash', 'geoid', 'geoidHeight', 'GeoidHeights', 'georef', 'gnomonic',
             'Hausdorff', 'Haversine', 'height', 'Hodgman', 'horizon', 'Hubeny',
             'IDW', 'interpolate', 'intersect', 'intersections', 'Inverse-Distance-Weighting', 'ITRF',
             'Karney', 'Krueger', 'Krüger',
             'Lambert', 'latitude', 'law-of-cosines', 'Lesh', 'linearize',
             'LocalCartesian', 'longitude',
             'Mercator', 'MGRS',
             'nearest', 'numpy', 'n-vector', 'Nvector',
             'orthographic', 'OSGR',
             'perimeter', 'Peucker', 'polar', 'Pseudo-Mercator',
             'PyGeodesy', 'PyInstaller', 'PyPy',
             'radii', 'radius', 'Ramer', 'Ramer-Douglas-Peucker',
             'Rey-Jer', 'Reumann', 'Reumann-Witkam', 'rhumb',
             'scipy', 'simplify', 'Snyder', 'Soldner', 'sphere', 'stereographic',
             'Sudano', 'Sutherland', 'Sutherland-Hodgman',
             'Terrestrial-Reference-Frame', 'Thomas', 'TMcoords', 'TMExact',
             'TransverseMercatorExact', 'TRF', 'trigonometry',
             'unroll', 'UPS', 'UTM', 'UTM/UPS',
             'Veness', 'Vermeille', 'Vincenty', 'Visvalingam', 'Visvalingam-Whyatt',
             'Web-Mercator', 'WGRS', 'WGS', 'Whyatt', 'Witkam', 'You')

setup(name='PyGeodesy',
      packages=['pygeodesy'],
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

#     data_files=[('docs',         ['docs/*.*']),
#                 ('images',       ['test/testRoute.jpg']),
#                 ('test',         ['test/test*.py']),
#                 ('testcoverage', ['testcoverage/*.*',
#                                   'testcoverage.pdf',
#                                   'testcoverage.rc']),
#                 ('testresults',  ['testresults/*.txt'])],
#     data_files fails somehow, see file MANIFEST.in

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
                   _c2('Programming Language', 'Python', '3.6'),
                   _c2('Programming Language', 'Python', '3.7'),
                   _c2('Programming Language', 'Python', '3.8'),
#                  _c2('Programming Language', 'Python', '3.9'),
                   _c2('Topic', 'Software Development'),
                   _c2('Topic', 'Scientific/Engineering', 'GIS'),
      ],

#     download_url='https://GitHub.com/mrJean1/PyGeodesy',
#     entry_points={},
#     include_package_data=False,
#     install_requires=[],
#     namespace_packages=[],
#     py_modules=[],
)
