
# -*- coding: utf-8 -*-

# The setuptools script to build, install and test a PyGeodesy distribution.

# Tested with 64-bit Python 2.7.13, 2.7.14, 3.6.1 and 3.6.2 (using setuptools
# 28.8.0) but only on macOS 10.12.3, 10.12.4 and 10.12.5, 10.12.6 Sierra.

# python setup.py sdist --formats=gztar,bztar,zip  # ztar,tar
# python setup.py bdist_wheel --universal  # XXX
# python setup.py test
# python setup.py install

# <https://packaging.python.org/key_projects/#setuptools>
# <https://packaging.python.org/distributing/>
# <https://docs.python.org/2/distutils/sourcedist.html>
# <https://docs.python.org/3.6/distutils/sourcedist.html>
# <https://setuptools.readthedocs.io/en/latest/setuptools.html#developer-s-guide>
# <https://setuptools.readthedocs.io/en/latest/setuptools.html#test-build-package-and-run-a-unittest-suite>
# <http://zetcode.com/articles/packageinpython/>

from setuptools import setup

__all__ = ()
__version__ = '17.08.24'


def _c2(*names):
    return ' :: '.join(names)


def _version():
    with open('pygeodesy/__init__.py') as f:
        for t in f.readlines():
            if t.startswith('__version__'):
                v = t.split('=')[-1].strip().strip("'").strip('"')
                return '.'.join(map(str, map(int, v.split('.'))))


def _long_description():
    with open('README.rst', 'rb') as f:
        t = f.read()
        if isinstance(t, bytes):
            t = t.decode('utf-8')
        return t


_KeyWords=('azimuth', 'bearing', 'cartesian', 'conic',
           'datum', 'development', 'distance', 'earth', 'ellipsoid',
           'geocentric', 'geodesy', 'geodetic', 'geohash', 'haversine',
           'Lambert', 'latitude', 'longitude', 'MGRS',
           'numpy', 'Nvector', 'OSGR', 'PyGeodesy',
           'Ramer-Douglas-Peucker', 'Reumann-Witkam', 'rhumb',
           'simplify', 'sphere', 'trigonometry', 'UTM',
           'Vincenty', 'Visvalingam-Whyatt', 'WGS')

setup(
    name='PyGeodesy',
    packages=['pygeodesy'],
    description='Pure Python geodesy tools',
    version=_version(),

    author='Jean M. Brouwers',
    author_email='mrJean1@Gmail.com',  # 'mrJean1 at Gmail dot com'
    maintainer='Jean M. Brouwers',
    maintainer_email='mrJean1@Gmail.com',  # 'mrJean1 at Gmail dot com'

    license='MIT',
    keywords=' '.join(_KeyWords),
    url='https://github.com/mrJean1/PyGeodesy',

    long_description=_long_description(),

    package_data={'pygeodesy': ['LICENSE']},

#   data_files=[('docs',        ['docs/*.*']),
#               ('images',      ['test/testRoute.jpg']),
#               ('test',        ['test/test*.py']),
#               ('testresults', ['testresults/*.txt'])],
#   data_files fails somehow, see file MANIFEST.in

    test_suite='test.TestSuite',

    zip_safe=False,
    classifiers=[
        _c2('Development Status', '5 - Production/Stable'),
        _c2('Environment', 'Console'),
        _c2('Intended Audience', 'Developers'),
        _c2('License', 'OSI Approved', 'MIT License'),
        _c2('Operating System', 'OS Independent'),
        _c2('Programming Language', 'Python'),
        _c2('Programming Language', 'Python', '2.6'),
        _c2('Programming Language', 'Python', '2.7'),
        _c2('Programming Language', 'Python', '3.5'),
        _c2('Programming Language', 'Python', '3.6'),
        _c2('Topic', 'Software Development'),
        _c2('Topic', 'Scientific/Engineering', 'GIS'),
    ],

#   download_url='https://github.com/mrJean1/PyGeodesy',
#   entry_points={},
#   include_package_data=False,
#   install_requires=[],
#   namespace_packages=[],
#   py_modules=[],
)
