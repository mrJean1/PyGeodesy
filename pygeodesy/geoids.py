
# -*- coding: utf-8 -*-

u'''Geoid models and geoid height interpolations.

Classes L{GeoidEGM96}, L{GeoidG2012B}, L{GeoidKarney} and L{GeoidPGM} to
interpolate the height of various U{geoid<https://WikiPedia.org/wiki/Geoid>}s
at C{LatLon} locations or separate lat-/longitudes using various interpolation
methods and C{geoid} model files.

L{GeoidKarney} is a transcoding of I{Charles Karney}'s C++ class U{Geoid
<https://GeographicLib.SourceForge.io/C++/doc/geoid.html>} to pure Python.

The L{GeoidEGM96}, L{GeoidG2012B} and L{GeoidPGM} interpolators both depend on
U{scipy<https://SciPy.org>} and U{numpy<https://PyPI.org/project/numpy>} and
require those packages to be installed.

In addition, each geoid interpolator needs C{grid knots} (down)loaded from
a C{geoid} model file, I{specific to the interpolator}.  More details below
and in the documentation of the interpolator class.  For each interpolator,
there are several interpolation choices, like I{linear}, I{cubic}, etc.

Typical usage
=============

1. Choose an interpolator class L{GeoidEGM96}, L{GeoidG2012B}, L{GeoidKarney}
or L{GeoidPGM} and download a C{geoid} model file, containing locations with
known heights also referred to as the C{grid knots}.  See the documentation of
the interpolator class for references to available C{grid} models.

C{>>> from pygeodesy import GeoidEGM96  # or -G2012B, -Karney or -PGM as GeoidXyz}

2. Instantiate an interpolator with the C{geoid} model file and use keyword
arguments to select different interpolation options

C{>>> ginterpolator = GeoidXyz(geoid_model_file, **options)}

3. Get the interpolated geoid height of C{LatLon} location(s) with

C{>>> ll = LatLon(1, 2, ...)}

C{>>> h = ginterpolator(ll)}

or

C{>>> h1, h2, h3, ... = ginterpolator(ll1, ll2, ll3, ...)}

or a list, tuple, generator, etc. of C{LatLon}s

C{>>> hs = ginterpolator(lls)}

4. For separate lat- and longitudes invoke the C{height} method as

C{>>> h = ginterpolator.height(lat, lon)}

or as 2 lists, 2 tuples, etc.

C{>>> hs = ginterpolator.height(lats, lons)}

or for several positionals use the C{height_} method

C{>>> h1, h2, ... = ginterpolator.height_(lat1, lon1, lat2, lon2, ...)}

5. An example is in U{issue #64<https://GitHub.com/mrJean1/PyGeodesy/issues/64>},
courtesy of SBFRF.

@note: Classes L{GeoidEGM96}, L{GeoidG2012B} and L{GeoidPGM} require both U{numpy
       <https://PyPI.org/project/numpy>} and U{scipy<https://PyPI.org/project/scipy>}
       to be installed.

@note: Errors from C{scipy} are raised as L{SciPyError}s.  Warnings issued by C{scipy} can
       be thrown as L{SciPyWarning} exceptions, provided Python C{warnings} are filtered
       accordingly, see L{SciPyWarning}.

@see: I{Karney}'s U{GeographicLib<https://GeographicLib.SourceForge.io/C++/doc/index.html>},
      U{Geoid height<https://GeographicLib.SourceForge.io/C++/doc/geoid.html>} and U{Installing
      the Geoid datasets<https://GeographicLib.SourceForge.io/C++/doc/geoid.html#geoidinst>},
      World Geodetic System 1984 (WG84) and U{Earth Gravitational Model 96 (EGM96) Data and
      Apps<https://earth-info.NGA.mil/index.php?dir=wgs84&action=wgs84>},
      U{SciPy<https://docs.SciPy.org/doc/scipy/reference/interpolate.html>} interpolation
      U{RectBivariateSpline<https://docs.SciPy.org/doc/scipy/reference/generated/scipy.interpolate.
      RectBivariateSpline.html>}, U{bisplrep/-ev<https://docs.scipy.org/doc/scipy/reference/generated/
      scipy.interpolate.bisplrep.html>} and U{interp2d<https://docs.SciPy.org/doc/scipy/reference/
      generated/scipy.interpolate.interp2d.html>}, functions L{elevations.elevation2} and
      L{elevations.geoidHeight2}, U{I{Ellispoid vs Orthometric Elevations}<https://www.YouTube.com/
      watch?v=dX6a6kCk3Po>} and U{6.22.1 Avoiding Pitfalls Related to Ellipsoid Height and Height
      Above Mean Sea Level<https://Wiki.ROS.org/mavros>}.
'''
# make sure int/int division yields float quotient, see .basics
from __future__ import division as _; del _  # noqa: E702 ;

from pygeodesy.basics import _isin, len2, min2, isodd, _splituple, \
                              ub2str as _ub2str
from pygeodesy.constants import EPS, _float as _F, _1_0, _N_90_0, _180_0, \
                                     _N_180_0, _360_0
from pygeodesy.datums import Datums, _ellipsoidal_datum, _WGS84
# from pygeodesy.dms import parseDMS2  # _MODS
from pygeodesy.errors import _incompatible, LenError, RangeError, _SciPyIssue, \
                             _xkwds_pop2
from pygeodesy.fmath import favg, Fdot, fdot, Fhorner, frange
# from pygoedesy.formy import heightOrthometric  # _MODS
from pygeodesy.heights import _as_llis2, _ascalar, _HeightBase, HeightError, \
                              _Wrap
# from pygeodesy.internals import typename, _version2  # _MODS
from pygeodesy.interns import NN, _COLONSPACE_, _COMMASPACE_, _DMAIN_, _E_, \
                             _height_, _in_, _kind_, _lat_, _lon_, _mean_, _N_, \
                             _n_a_, _numpy_, _on_, _outside_, _S_, _s_, _scipy_, \
                             _SPACE_, _stdev_, _tbd_, _W_, _width_, _4_
from pygeodesy.lazily import _ALL_DOCS, _ALL_LAZY, _ALL_MODS as _MODS, _FOR_DOCS
from pygeodesy.named import _name__, _Named, _NamedTuple
# from pygeodesy.namedTuples import LatLon3Tuple  # _MODS
from pygeodesy.props import Property_RO, property_RO, property_ROver
from pygeodesy.streprs import attrs, Fmt, fstr, pairs
from pygeodesy.units import Height, Int_, Lat, Lon
# from pygeodesy.utily import _Wrap  # from .heights

from math import floor as _floor
# from os import SEEK_CUR, SEEK_SET  # _MODS
# import os.path  # _MODS
from struct import calcsize as _calcsize, unpack as _unpack
try:
    from StringIO import StringIO as _BytesIO  # reads bytes
    _ub2str = str  # PYCHOK convert bytes to str for egm*.pgm text

except ImportError:  # Python 3+
    from io import BytesIO as _BytesIO  # PYCHOK expected

__all__ = _ALL_LAZY.geoids
__version__ = '25.05.21'

_assert_    = 'assert'
_bHASH_     = b'#'
_endian_    = 'endian'
_format_    = '%s %r'
_header_    = 'header'
_intCs      = {}  # cache int value, del below
_lli_       = 'lli'
_rb_        = 'rb'
_supported_ = 'supported'


class GeoidError(HeightError):
    '''Geoid interpolator C{Geoid...} or interpolation issue.
    '''
    pass


class _GeoidBase(_HeightBase):
    '''(INTERNAL) Base class for C{Geoid...}s.
    '''
#   _center   =  None   # (lat, lon, height)
    _cropped  =  None
#   _datum    = _WGS84  # from _HeightBase
    _egm      =  None   # open C{egm*.pgm} geoid file
    _endian   = _tbd_
    _Error    =  GeoidError  # in ._HeightBase._as_lls, ...
    _geoid    = _n_a_
    _hs_y_x   =  None   # numpy 2darray, row-major order
    _iscipy   =  True   # scipy or Karney's interpolation
    _kind     =  3      # order for interp2d, RectBivariateSpline
#   _kmin     =  2      # min number of knots
    _knots    =  0      # nlat * nlon
    _mean     =  None   # fixed in GeoidKarney
    _nBytes   =  0      # numpy size in bytes, float64
    _pgm      =  None   # PGM attributes, C{_PGM} or C{None}
    _sizeB    =  0      # geoid file size in bytes
    _smooth   =  0      # used only for RectBivariateSpline
    _stdev    =  None   # fixed in GeoidKarney
    _u2B      =  0      # np.itemsize or undefined
    _yx_hits  =  None   # cache hits, ala Karney's

#   _lat_d  = _0_0  # increment, +tive
#   _lat_lo = _0_0  # lower lat, south
#   _lat_hi = _0_0  # upper lat, noth
#   _lon_d  = _0_0  # increment, +tive
#   _lon_lo = _0_0  # left lon, west
#   _lon_hi = _0_0  # right lon, east
#   _lon_of = _0_0  # forward lon offset
#   _lon_og = _0_0  # reverse lon offset

    def __init__(self, hs, p):
        '''(INTERNAL) Set up the grid axes, the C{SciPy} interpolator and
           several internal geoid attributes.

           @arg hs: Grid knots with known height (C{numpy 2darray}).
           @arg p: The C{slat, wlon, nlat, nlon, dlat, dlon} and other
                   geoid parameters (C{INTERNAL}).
        '''
        spi = self.scipy_interpolate
        # for 2d scipy.interpolate.interp2d(xs, ys, hs, ...) and
        # scipy.interpolate.RectBivariateSpline(ys, xs, hs, ...)
        # require the shape of hs to be (len(ys), len(xs)), note
        # the different (xs, ys, ...) and (ys, xs, ...) orders
        if (p.nlat, p.nlon) != hs.shape:
            raise GeoidError(shape=hs.shape, txt=_incompatible((p.nlat, p.nlon)))

        # both axes and bounding box
        ys, self._lat_d = self._gaxis2(p.slat, p.dlat, p.nlat, _lat_ + _s_)
        xs, self._lon_d = self._gaxis2(p.wlon, p.dlon, p.nlon, _lon_ + _s_)

        bb = ys[0], ys[-1], xs[0], xs[-1] + p.dlon  # fudge lon_hi
        # geoid grids are typically stored in row-major order, some
        # with rows (90..-90) reversed and columns (0..360) wrapped
        # to Easten longitude, 0 <= east < 180 and 180 <= west < 360
        k = self.kind
        if k in self._k2interp2d:  # see _HeightBase
            self._interp2d(xs, ys, hs, kind=k)
        else:  # XXX order ys and xs, see HeightLSQBiSpline
            k = self._kxky(k)
            self._ev = spi.RectBivariateSpline(ys, xs, hs, bbox=bb, ky=k, kx=k,
                                                           s=self._smooth).ev
        self._hs_y_x = hs  # numpy 2darray, row-major
        self._nBytes = hs.nbytes  # numpy size in bytes
        self._knots  = p.knots  # grid knots == len(hs)
        self._lon_of =     float(p.flon)  # forward offset
        self._lon_og = g = float(p.glon)  # reverse offset
        # shrink the bounding box by 1 unit on every side:
        # +self._lat_d, -self._lat_d, +self._lon_d, -self._lon_d
        self._lat_lo, \
        self._lat_hi, \
        self._lon_lo, \
        self._lon_hi  = map(float, bb)
        self._lon_lo -= g
        self._lon_hi -= g

    def __call__(self, *llis, **wrap_H):
        '''Interpolate the geoid (or orthometric) height for one or more locations.

           @arg llis: One or several locations (each C{LatLon}), all positional.
           @kwarg wrap_H: Keyword arguments C{B{wrap}=False} (C{bool}) and
                       C{B{H}=False} (C{bool}).  Use C{B{wrap}=True} to wrap
                       or I{normalize} all B{C{llis}} locations.  If C{B{H}
                       is True}, return the I{orthometric} height instead of
                       the I{geoid} height at each location.

           @return: A single geoid (or orthometric) height (C{float}) or
                    a list or tuple of geoid (or orthometric) heights (each
                    C{float}).

           @raise GeoidError: Insufficient number of B{C{llis}}, an invalid
                              B{C{lli}} or the C{egm*.pgm} geoid file is closed.

           @raise RangeError: An B{C{lli}} is outside this geoid's lat- or
                              longitude range.

           @raise SciPyError: A C{scipy} issue.

           @raise SciPyWarning: A C{scipy} warning as exception.

           @note: To obtain I{orthometric} heights, each B{C{llis}} location
                  must have an ellipsoid C{height} or C{h} attribute, otherwise
                  C{height=0} is used.

           @see: Function L{pygeodesy.heightOrthometric}.
        '''
        return self._called(llis, **wrap_H)

    def __enter__(self):
        '''Open context.
        '''
        return self

    def __exit__(self, *unused):  # PYCHOK exc_type, exc_value, exc_traceback)
        '''Close context.
        '''
        self.close()
        # return None  # XXX False

    def __repr__(self):
        return self.toStr()

    def __str__(self):
        return Fmt.PAREN(self.classname, repr(self.name))

    def _called(self, llis, wrap=False, H=False):
        # handle __call__
        _H = self._heightOrthometric if H else None
        _as, llis = _as_llis2(llis, Error=GeoidError)
        _w, hs    = _Wrap._latlonop(wrap), []
        _h, _a    =  self._hGeoid, hs.append
        try:
            for i, lli in enumerate(llis):
                N = _h(*_w(lli.lat, lli.lon))
                # orthometric or geoid height
                _a(_H(lli, N) if _H else N)
            return _as(hs)
        except (GeoidError, RangeError) as x:
            # XXX avoid str(LatLon()) degree symbols
            n = _lli_ if _as is _ascalar else Fmt.INDEX(llis=i)
            t = fstr((lli.lat, lli.lon), strepr=repr)
            raise type(x)(n, t, wrap=wrap, H=H, cause=x)
        except Exception as x:
            if self._iscipy and self.scipy:
                raise _SciPyIssue(x, self._ev_name)
            else:
                raise

    @Property_RO
    def _center(self):
        '''(INTERNAL) Cache for method L{center}.
        '''
        return self._llh3(favg(self._lat_lo, self._lat_hi),
                          favg(self._lon_lo, self._lon_hi))

    def center(self, LatLon=None):
        '''Return the center location and height of this geoid.

           @kwarg LatLon: Optional class to return the location and height
                          (C{LatLon}) or C{None}.

           @return: If C{B{LatLon} is None}, a L{LatLon3Tuple}C{(lat, lon,
                    height)} otherwise a B{C{LatLon}} instance with the lat-,
                    longitude and geoid height of the center grid location.
        '''
        return self._llh3LL(self._center, LatLon)

    def close(self):
        '''Close the C{egm*.pgm} geoid file if open (and applicable).
        '''
        if not self.closed:
            self._egm.close()
            self._egm = None

    @property_RO
    def closed(self):
        '''Get the C{egm*.pgm} geoid file status.
        '''
        return self._egm is None

    @Property_RO
    def cropped(self):
        '''Is geoid cropped (C{bool} or C{None} if crop not supported).
        '''
        return self._cropped

    @Property_RO
    def dtype(self):
        '''Get the grid C{scipy} U{dtype<https://docs.SciPy.org/doc/numpy/
           reference/generated/numpy.ndarray.dtype.html>} (C{numpy.dtype}).
        '''
        return self._hs_y_x.dtype

    @Property_RO
    def endian(self):
        '''Get the geoid endianess and U{dtype<https://docs.SciPy.org/
           doc/numpy/reference/generated/numpy.dtype.html>} (C{str}).
        '''
        return self._endian

    def _ev(self, y, x):  # PYCHOK overwritten with .RectBivariateSpline.ev
        # see methods _HeightBase._ev and -._interp2d
        return self._ev2d(x, y)  # (y, x) flipped!

    def _gaxis2(self, lo, d, n, name):
        # build grid axis, hi = lo + (n - 1) * d
        m, a = len2(frange(lo, n, d))
        if m != n:
            raise LenError(type(self), grid=m, **{name: n})
        if d < 0:
            d, a = -d, list(reversed(a))
        a = self.numpy.array(a)
        m, i = min2(*map(float, a[1:] - a[:-1]))
        if m < EPS:  # non-increasing axis
            i = Fmt.INDEX(name, i + 1)
            raise GeoidError(i, m, txt_not_='increasing')
        return a, d

    def _g2ll2(self, lat, lon):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}.'''
        self._notOverloaded(lat, lon)

    def _gyx2g2(self, y, x):
        # convert grid (y, x) indices to grid (lat, lon)
        return ((self._lat_lo                + self._lat_d * y),
                (self._lon_lo + self._lon_of + self._lon_d * x))

    def height(self, lats, lons, **wrap):
        '''Interpolate the geoid height for one or several lat-/longitudes.

           @arg lats: Latitude or latitudes (each C{degrees}).
           @arg lons: Longitude or longitudes (each C{degrees}).
           @kwarg wrap: Use C{B{wrap}=True} to wrap or I{normalize} all
                        B{C{lats}} and B{C{lons}}.

           @return: A single geoid height (C{float}) or a list of geoid
                    heights (each C{float}).

           @raise GeoidError: Insufficient or unequal number of B{C{lats}}
                              and B{C{lons}}.

           @raise RangeError: A B{C{lat}} or B{C{lon}} is outside this geoid's
                              lat- or longitude range.

           @raise SciPyError: A C{scipy} issue.

           @raise SciPyWarning: A C{scipy} warning as exception.
        '''
        return _HeightBase.height(self, lats, lons, **wrap)

    def height_(self, *latlons, **wrap):
        '''Interpolate the geoid height for each M{(latlons[i], latlons[i+1])
           pair for i in range(0, len(latlons), B{2})}.

           @arg latlons: Alternating lat-/longitude pairs (each C{degrees}),
                         all positional.

           @see: Method L{height} for further details.

           @return: A tuple of geoid heights (each C{float}).
        '''
        lls = tuple(self._as_lls(latlons[0::2], *latlons[1::2]))
        return self._called(lls, **wrap)

    @property_ROver
    def _heightOrthometric(self):
        return _MODS.formy.heightOrthometric  # overwrite property_ROver

    def _hGeoid(self, lat, lon):
        out = self.outside(lat, lon)
        if out:  # XXX avoid str(LatLon()) degree symbols
            t = fstr((lat, lon), strepr=repr)
            raise RangeError(lli=t, txt=_SPACE_(_outside_, _on_, out))
        return float(self._ev(*self._ll2g2(lat, lon)))

    @Property_RO
    def _highest(self):
        '''(INTERNAL) Cache for C{.highest}.
        '''
        return self._LL3T(self._llh3minmax(True), name__=self.highest)

    def highest(self, LatLon=None, **unused):
        '''Return the location and largest height of this geoid.

           @kwarg LatLon: Optional class to return the location and height
                          (C{LatLon}) or C{None}.

           @return: If C{B{LatLon} is None}, a L{LatLon3Tuple}C{(lat, lon,
                    height)} otherwise a B{C{LatLon}} instance with the lat-,
                    longitude and geoid height of the highest grid location.
        '''
        return self._llh3LL(self._highest, LatLon)

    @Property_RO
    def hits(self):
        '''Get the number of cache hits (C{int} or C{None}).
        '''
        return self._yx_hits

    @Property_RO
    def kind(self):
        '''Get the interpolator kind and order (C{int}).
        '''
        return self._kind

    @Property_RO
    def knots(self):
        '''Get the number of grid knots (C{int}).
        '''
        return self._knots

    def _ll2g2(self, lat, lon):  # PYCHOK no cover
        '''(INTERNAL) I{Must be overloaded}.'''
        self._notOverloaded(lat, lon)

    @property_ROver
    def _LL3T(self):
        '''(INTERNAL) Get L{LatLon3Tuple}, I{once}.
        '''
        return _MODS.namedTuples.LatLon3Tuple  # overwrite property_ROver

    def _llh3(self, lat, lon):
        return self._LL3T(lat, lon, self._hGeoid(lat, lon), name=self.name)

    def _llh3LL(self, llh, LatLon):
        return llh if LatLon is None else self._xnamed(LatLon(*llh))

    def _llh3minmax(self, highest, *unused):
        hs, np = self._hs_y_x, self.numpy
        # <https://docs.SciPy.org/doc/numpy/reference/generated/
        #         numpy.argmin.html#numpy.argmin>
        arg  = np.argmax if highest else np.argmin
        y, x = np.unravel_index(arg(hs, axis=None), hs.shape)
        return self._g2ll2(*self._gyx2g2(y, x)) + (float(hs[y, x]),)

    def _load(self, g, dtype=float, n=-1, offset=0, **sep):  # sep=NN
        # numpy.fromfile, like .frombuffer
        g.seek(offset, _MODS.os.SEEK_SET)
        return self.numpy.fromfile(g, dtype, count=n, **sep)

    @Property_RO
    def _lowerleft(self):
        '''(INTERNAL) Cache for C{.lowerleft}.
        '''
        return self._llh3(self._lat_lo, self._lon_lo)

    def lowerleft(self, LatLon=None):
        '''Return the lower-left location and height of this geoid.

           @kwarg LatLon: Optional class to return the location
                          (C{LatLon}) and height or C{None}.

           @return: If C{B{LatLon} is None}, a L{LatLon3Tuple}C{(lat, lon, height)}
                    otherwise a B{C{LatLon}} instance with the lat-, longitude and
                    geoid height of the lower-left, SW grid corner.
        '''
        return self._llh3LL(self._lowerleft, LatLon)

    @Property_RO
    def _loweright(self):
        '''(INTERNAL) Cache for C{.loweright}.
        '''
        return self._llh3(self._lat_lo, self._lon_hi)

    def loweright(self, LatLon=None):
        '''Return the lower-right location and height of this geoid.

           @kwarg LatLon: Optional class to return the location and height
                          (C{LatLon}) or C{None}.

           @return: If C{B{LatLon} is None}, a L{LatLon3Tuple}C{(lat, lon, height)}
                    otherwise a B{C{LatLon}} instance with the lat-, longitude and
                    geoid height of the lower-right, SE grid corner.
        '''
        return self._llh3LL(self._loweright, LatLon)

    lowerright = loweright  # synonymous

    @Property_RO
    def _lowest(self):
        '''(INTERNAL) Cache for C{.lowest}.
        '''
        return self._LL3T(self._llh3minmax(False), name__=self.lowest)

    def lowest(self, LatLon=None, **unused):
        '''Return the location and lowest height of this geoid.

           @kwarg LatLon: Optional class to return the location and height
                          (C{LatLon}) or C{None}.

           @return: If C{B{LatLon} is None}, a L{LatLon3Tuple}C{(lat, lon,
                    height)} otherwise a B{C{LatLon}} instance with the lat-,
                    longitude and geoid height of the lowest grid location.
        '''
        return self._llh3LL(self._lowest, LatLon)

    @Property_RO
    def mean(self):
        '''Get the mean of this geoid's heights (C{float}).
        '''
        if self._mean is None:  # see GeoidKarney
            self._mean = float(self.numpy.mean(self._hs_y_x))
        return self._mean

    @property_RO
    def name(self):
        '''Get the name of this geoid (C{str}).
        '''
        return _HeightBase.name.fget(self) or self._geoid  # recursion

    @Property_RO
    def nBytes(self):
        '''Get the grid in-memory size in bytes (C{int}).
        '''
        return self._nBytes

    def _open(self, geoid, datum, kind, name, smooth):
        # open the geoid file
        try:
            self._geoid = _MODS.os.path.basename(geoid)
            self._sizeB = _MODS.os.path.getsize(geoid)
            g = open(geoid, _rb_)
        except (IOError, OSError) as x:
            raise GeoidError(geoid=geoid, cause=x)

        if not _isin(datum, None, self._datum):
            self._datum = _ellipsoidal_datum(datum, name=name)
        self._kind = int(kind)
        if name:
            _HeightBase.name.fset(self, name)  # rename
        if smooth:
            self._smooth = Int_(smooth=smooth, Error=GeoidError, low=0)

        return g

    def outside(self, lat, lon):
        '''Check whether a location is outside this geoid's lat-/longitude
           or crop range.

           @arg lat: The latitude (C{degrees}).
           @arg lon: The longitude (C{degrees}).

           @return: A 1- or 2-character C{str} if outside, an empty C{str} otherwise.
        '''
        lat = _S_ if lat < self._lat_lo else (_N_ if lat > self._lat_hi else NN)
        lon = _W_ if lon < self._lon_lo else (_E_ if lon > self._lon_hi else NN)
        return NN(lat, lon) if lat and lon else (lat or lon)

    @Property_RO
    def pgm(self):
        '''Get the PGM attributes (C{_PGM} or C{None} if not available/applicable).
        '''
        return self._pgm

    @Property_RO
    def sizeB(self):
        '''Get the geoid grid file size in bytes (C{int}).
        '''
        return self._sizeB

    @Property_RO
    def smooth(self):
        '''Get the C{RectBivariateSpline} smoothing (C{int}).
        '''
        return self._smooth

    @Property_RO
    def stdev(self):
        '''Get the standard deviation of this geoid's heights (C{float}) or C{None}.
        '''
        if self._stdev is None:  # see GeoidKarney
            self._stdev = float(self.numpy.std(self._hs_y_x))
        return self._stdev

    def _swne(self, crop):
        # crop box to 4-tuple (s, w, n, e)
        try:
            if len(crop) == 2:
                try:  # sw, ne LatLons
                    swne = (crop[0].lat, crop[0].lon,
                            crop[1].lat, crop[1].lon)
                except AttributeError:  # (s, w), (n, e)
                    swne = tuple(crop[0]) + tuple(crop[1])
            else:  # (s, w, n, e)
                swne = crop
            if len(swne) == 4:
                s, w, n, e = map(float, swne)
                if _N_90_0 <= s <= (n - _1_0) <=  89.0 and \
                  _N_180_0 <= w <= (e - _1_0) <= 179.0:
                    return s, w, n, e
        except (IndexError, TypeError, ValueError):
            pass
        raise GeoidError(crop=crop)

    def toStr(self, prec=3, sep=_COMMASPACE_):  # PYCHOK signature
        '''This geoid and all geoid attributes as a string.

           @kwarg prec: Number of decimal digits (0..9 or C{None} for
                        default).  Trailing zero decimals are stripped
                        for B{C{prec}} values of 1 and above, but kept
                        for negative B{C{prec}} values.
           @kwarg sep: Separator to join (C{str}).

           @return: Geoid name and attributes (C{str}).
        '''
        s =  1 if self.kind < 0 else 2
        t = _MODS.internals.typename
        t =  tuple(Fmt.PAREN(t(m), fstr(m(), prec=prec)) for m in
                                  (self.lowerleft, self.upperright,
                                   self.center,
                                   self.highest, self.lowest)) + \
             attrs( _mean_, _stdev_,           prec=prec, Nones=False) + \
             attrs((_kind_, 'smooth')[:s],     prec=prec, Nones=False) + \
             attrs( 'cropped', 'dtype', _endian_, 'hits', 'knots', 'nBytes',
                    'sizeB', _scipy_, _numpy_, prec=prec, Nones=False)
        return _COLONSPACE_(self, sep.join(t))

    @Property_RO
    def u2B(self):
        '''Get the PGM itemsize in bytes (C{int}).
        '''
        return self._u2B

    @Property_RO
    def _upperleft(self):
        '''(INTERNAL) Cache for C{.upperleft}.
        '''
        return self._llh3(self._lat_hi, self._lon_lo)

    def upperleft(self, LatLon=None):
        '''Return the upper-left location and height of this geoid.

           @kwarg LatLon: Optional class to return the location and height
                          (C{LatLon}) or C{None}.

           @return: If C{B{LatLon} is None}, a L{LatLon3Tuple}C{(lat, lon, height)}
                    otherwise a B{C{LatLon}} instance with the lat-, longitude and
                    geoid height of the upper-left, NW grid corner.
        '''
        return self._llh3LL(self._upperleft, LatLon)

    @Property_RO
    def _upperright(self):
        '''(INTERNAL) Cache for C{.upperright}.
        '''
        return self._llh3(self._lat_hi, self._lon_hi)

    def upperright(self, LatLon=None):
        '''Return the upper-right location and height of this geoid.

           @kwarg LatLon: Optional class to return the location and height
                          (C{LatLon}) or C{None}.

           @return: If C{B{LatLon} is None}, a L{LatLon3Tuple}C{(lat, lon, height)}
                    otherwise a B{C{LatLon}} instance with the lat-, longitude and
                    geoid height of the upper-right, NE grid corner.
        '''
        return self._llh3LL(self._upperright, LatLon)


class GeoidEGM96(_GeoidBase):
    '''Geoid height interpolator for the EGM96 U{15 Minute Interpolation Grid<https://earth-info.NGA.mil>}
       based on C{SciPy} interpolation U{RectBivariateSpline<https://docs.SciPy.org/doc/scipy/reference/
       generated/scipy.interpolate.RectBivariateSpline.html>}, U{interp2d<https://docs.SciPy.org/doc/scipy/
       reference/generated/scipy.interpolate.interp2d.html>} or U{bisplrep/-ev<https://docs.scipy.org/doc/
       scipy/reference/generated/scipy.interpolate.bisplrep.html>}.

       Use only the C{WW15MGH.GRD} file, unzipped from the EGM96 U{15 Minute Interpolation Grid
       <https://earth-info.NGA.mil/index.php?dir=wgs84&action=wgs84>} download.
    '''
    def __init__(self, EGM96_grd, datum=_WGS84, kind=3, smooth=0, **name_crop):
        '''New L{GeoidEGM96} interpolator.

           @arg EGM96_grd: An C{EGM96_grd} grid file name (C{.GRD}).
           @kwarg datum: Optional grid datum (L{Datum}, L{Ellipsoid}, L{Ellipsoid2} or L{a_f2Tuple}),
                         overriding C{WGS84}.
           @kwarg kind: C{scipy.interpolate} order (C{int}), use 1..5 for U{RectBivariateSpline
                        <https://docs.SciPy.org/doc/scipy/reference/generated/scipy.interpolate.
                        RectBivariateSpline.html>} or -1, -3 or -5 for U{bisplrep/-ev<https://
                        docs.SciPy.org/doc/scipy/reference/generated/scipy.interpolate.bisplrep.html>}
                        or U{interp2d<https://docs.SciPy.org/doc/scipy/reference/generated/scipy.
                        interpolate.interp2d.html>} C{linear}, C{cubic} respectively C{quintic},
                        see note for more details.
           @kwarg smooth: Smoothing factor for C{B{kind}=1..5} only (C{int}).
           @kwarg name_crop: Optional geoid C{B{name}=NN} (C{str}) and UNSUPPORTED keyword argument
                       C{B{crop}=None}.

           @raise GeoidError: Invalid B{C{crop}}, B{C{kind}} or B{C{smooth}} or a ECM96 grid file
                              B{C{ECM96_grd}} issue.

           @raise ImportError: Package C{numpy} or C{scipy} not found or not installed.

           @raise LenError: Grid file B{C{EGM96_grd}} axis mismatch.

           @raise SciPyError: A C{scipy} issue.

           @raise SciPyWarning: A C{scipy} warning as exception.

           @raise TypeError: Invalid B{C{datum}}.

           @note: Specify C{B{kind}=-1, -3 or -5} to use C{scipy.interpolate.interp2d}
                  before or C{scipy.interpolate.bisplrep/-ev} since C{Scipy} version 1.14.
        '''
        crop, name = _xkwds_pop2(name_crop, crop=None)
        if crop is not None:
            raise GeoidError(crop=crop, txt_not_=_supported_)

        g = self._open(EGM96_grd, datum, kind, _name__(**name), smooth)
        _ = self.numpy  # import numpy for .fromfile, .reshape

        try:
            p, hs = _Gpars(), self._load(g, sep=_SPACE_)  # text
            p.slat, n, p.wlon, e, p.dlat, p.dlon = hs[:6]  # n-s, 0-E
            p.nlat  =  int((n - p.slat) / p.dlat) + 1  # include S
            p.nlon  =  int((e - p.wlon) / p.dlon) + 1  # include W
            p.knots =  p.nlat * p.nlon  # inverted lats N downto S
            p.glon  = _180_0  # Eastern lons 0-360
            hs = hs[6:].reshape(p.nlat, p.nlon)
            _GeoidBase.__init__(self, hs, p)

        except Exception as x:
            raise _SciPyIssue(x, _in_, repr(EGM96_grd))
        finally:
            g.close()

    def _g2ll2(self, lat, lon):
        # convert grid (lat, lon) to earth (lat, lon)
        while lon > _180_0:  # Eastern
            lon -= _360_0
        return -lat, lon  # invert lat

    def _ll2g2(self, lat, lon):
        # convert earth (lat, lon) to grid (lat, lon)
        while lon < 0:  # Eastern
            lon += _360_0
        return -lat, lon  # invert lat

    if _FOR_DOCS:
        __call__ = _GeoidBase.__call__
        height   = _GeoidBase.height


class GeoidG2012B(_GeoidBase):
    '''Geoid height interpolator for U{GEOID12B Model
       <https://www.NGS.NOAA.gov/GEOID/GEOID12B/>} grids U{CONUS
       <https://www.NGS.NOAA.gov/GEOID/GEOID12B/GEOID12B_CONUS.shtml>},
       U{Alaska<https://www.NGS.NOAA.gov/GEOID/GEOID12B/GEOID12B_AK.shtml>},
       U{Hawaii<https://www.NGS.NOAA.gov/GEOID/GEOID12B/GEOID12B_HI.shtml>},
       U{Guam and Northern Mariana Islands
       <https://www.NGS.NOAA.gov/GEOID/GEOID12B/GEOID12B_GMNI.shtml>},
       U{Puerto Rico and U.S. Virgin Islands
       <https://www.NGS.NOAA.gov/GEOID/GEOID12B/GEOID12B_PRVI.shtml>} and
       U{American Samoa<https://www.NGS.NOAA.gov/GEOID/GEOID12B/GEOID12B_AS.shtml>}
       based on C{SciPy} interpolation U{RectBivariateSpline<https://docs.SciPy.org/doc/
       scipy/reference/generated/scipy.interpolate.RectBivariateSpline.html>}, U{interp2d
       <https://docs.SciPy.org/doc/scipy/reference/generated/scipy.interpolate.interp2d.html>}
       or U{bisplrep/-ev<https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.
       bisplrep.html>}.

       Use any of the C{le} (little endian) or C{be} (big endian) C{g2012b*.bin} binary grid files.
    '''
    _datum = Datums.NAD83

    def __init__(self, g2012b_bin, datum=Datums.NAD83, kind=3, smooth=0, **name_crop):
        '''New L{GeoidG2012B} interpolator.

           @arg g2012b_bin: A C{GEOID12B} grid file name (C{.bin}).
           @kwarg datum: Optional grid datum (L{Datum}, L{Ellipsoid}, L{Ellipsoid2} or
                         L{a_f2Tuple}), overriding C{NAD83}.
           @kwarg kind: C{scipy.interpolate} order (C{int}), use 1..5 for U{RectBivariateSpline
                        <https://docs.SciPy.org/doc/scipy/reference/generated/scipy.interpolate.
                        RectBivariateSpline.html>} or -1, -3 or -5 for U{bisplrep/-ev<https://
                        docs.SciPy.org/doc/scipy/reference/generated/scipy.interpolate.bisplrep.html>}
                        or U{interp2d<https://docs.SciPy.org/doc/scipy/reference/generated/scipy.
                        interpolate.interp2d.html>} C{linear}, C{cubic} respectively C{quintic},
                        see note for more details.
           @kwarg smooth: Smoothing factor for C{B{kind}=1..5} only (C{int}).
           @kwarg name_crop: Optional geoid C{B{name}=NN} (C{str}) and UNSUPPORTED keyword argument
                       C{B{crop}=None}.

           @raise GeoidError: Invalid B{C{crop}}, B{C{kind}} or B{C{smooth}} or a G2012B grid file
                              B{C{g2012b_bin}} issue.

           @raise ImportError: Package C{numpy} or C{scipy} not found or not installed.

           @raise LenError: Grid file B{C{g2012b_bin}} axis mismatch.

           @raise SciPyError: A C{scipy} issue.

           @raise SciPyWarning: A C{scipy} warning as exception.

           @raise TypeError: Invalid B{C{datum}}.

           @note: Specify C{B{kind}=-1, -3 or -5} to use C{scipy.interpolate.interp2d}
                  before or C{scipy.interpolate.bisplrep/-ev} since C{Scipy} version 1.14.
        '''
        crop, name = _xkwds_pop2(name_crop, crop=None)
        if crop is not None:
            raise GeoidError(crop=crop, txt_not_=_supported_)

        g = self._open(g2012b_bin, datum, kind, _name__(**name), smooth)
        _ = self.numpy  # import numpy for ._load and

        try:
            p = _Gpars()
            n = (self.sizeB // 4) - 11  # number of f4 heights
            # U{numpy dtype formats are different from Python struct formats
            # <https://docs.SciPy.org/doc/numpy-1.15.0/reference/arrays.dtypes.html>}
            for en_ in ('<', '>'):
                # skip 4xf8, get 3xi4
                p.nlat, p.nlon, ien = map(int, self._load(g, en_+'i4', 3, 32))
                if ien == 1:  # correct endian
                    p.knots = p.nlat * p.nlon
                    if p.knots == n and 1 < p.nlat < n \
                                    and 1 < p.nlon < n:
                        self._endian = en_+'f4'
                        break
            else:  # couldn't validate endian
                raise GeoidError(_endian_)

            # get the first 4xf8
            p.slat, p.wlon, p.dlat, p.dlon = map(float, self._load(g, en_+'f8', 4))
            # read all f4 heights, ignoring the first 4xf8 and 3xi4
            hs = self._load(g, self._endian, n, 44).reshape(p.nlat, p.nlon)
            p.wlon -= _360_0  # western-most East longitude to earth (..., lon)
            _GeoidBase.__init__(self, hs, p)

        except Exception as x:
            raise _SciPyIssue(x, _in_, repr(g2012b_bin))
        finally:
            g.close()

    def _g2ll2(self, lat, lon):
        # convert grid (lat, lon) to earth (lat, lon)
        return lat, lon

    def _ll2g2(self, lat, lon):
        # convert earth (lat, lon) to grid (lat, lon)
        return lat, lon

    if _FOR_DOCS:
        __call__ = _GeoidBase.__call__
        height   = _GeoidBase.height


class GeoidHeight5Tuple(_NamedTuple):  # .geoids.py
    '''5-Tuple C{(lat, lon, egm84, egm96, egm2008)} for U{GeoidHeights.dat
       <https://SourceForge.net/projects/geographiclib/files/testdata/>}
       tests with the heights for 3 different EGM grids at C{degrees90}
       and C{degrees180} degrees (after converting C{lon} from original
       C{0 <= EasterLon <= 360}).
    '''
    _Names_ = (_lat_, _lon_, 'egm84', 'egm96', 'egm2008')
    _Units_ = ( Lat,   Lon,   Height,  Height,  Height)


def _I(i):
    '''(INTERNAL) Cache a single C{int} constant.
    '''
    i = int(i)
    return _intCs.setdefault(i, i)  # noqa: F821 del


def _T(cs):
    '''(INTERNAL) Cache a tuple of single C{int} constants.
    '''
    return tuple(map(_I, _splituple(cs)))

_T0s12 = (_I(0),) * 12  # PYCHOK _T('0, 0, ..., 0')


class GeoidKarney(_GeoidBase):
    '''Geoid height interpolator for I{Karney}'s U{GeographicLib Earth
       Gravitational Model (EGM)<https://GeographicLib.SourceForge.io/C++/doc/
       geoid.html>} geoid U{egm*.pgm<https://GeographicLib.SourceForge.io/
       C++/doc/geoid.html#geoidinst>} datasets using bilinear or U{cubic
       <https://dl.ACM.org/citation.cfm?id=368443>} interpolation and U{caching
       <https://GeographicLib.SourceForge.io/C++/doc/geoid.html#geoidcache>}
       in pure Python, transcoded from I{Karney}'s U{C++ class Geoid
       <https://GeographicLib.SourceForge.io/C++/doc/geoid.html#geoidinterp>}.

       Use any of the geoid U{egm84-, egm96- or egm2008-*.pgm
       <https://GeographicLib.SourceForge.io/C++/doc/geoid.html#geoidinst>}
       datasets.
    '''
    _C0 = _F(372), _F(240), _F(372)  # n, _ and s common denominators
    # matrices c3n_, c3, c3s_, transposed from GeographicLib/Geoid.cpp
    _C3 = ((_T('0,    0,   62,  124,  124,  62,    0,   0,   0,   0,    0,   0'),
            _T0s12,
         _T('-131,    7,  -31,  -62,  -62, -31,   45, 216, 156, -45,  -55,  -7'),
            _T0s12,
          _T('138, -138,    0,    0,    0,   0, -183,  33, 153,  -3,   48, -48'),  # PYCHOK indent
          _T('144,   42,  -62, -124, -124, -62,   -9,  87,  99,   9,   42, -42'),
            _T0s12,
            _T('0,    0,    0,    0,    0,   0,   93, -93, -93,  93,    0,   0'),
         _T('-102,  102,    0,    0,    0,   0,   18,  12, -12, -18,  -84,  84'),
          _T('-31,  -31,   31,   62,   62,  31,    0, -93, -93,   0,   31,  31')),  # PYCHOK indent

           (_T('9,   -9,    9,  186,   54,  -9,   -9,  54, -54,   9,   -9,   9'),
          _T('-18,   18,  -88,  -42,  162, -32,    8, -78,  78,  -8,   18, -18'),
          _T('-88,    8,  -18,  -42,  -78,  18,   18, 162,  78, -18,  -32,  -8'),
            _T('0,    0,   90, -150,   30,  30,   30, -90,  90, -30,    0,   0'),
           _T('96,  -96,   96,  -96,  -24,  24,  -96, -24, 144, -24,   24, -24'),  # PYCHOK indent
           _T('90,   30,    0, -150,  -90,   0,    0,  30,  90,   0,   30, -30'),
            _T('0,    0,  -20,   60,  -60,  20,  -20,  60, -60,  20,    0,   0'),
            _T('0,    0,  -60,   60,   60, -60,   60, -60, -60,  60,    0,   0'),
          _T('-60,   60,    0,   60,  -60,   0,    0,  60, -60,   0,  -60,  60'),
          _T('-20,  -20,    0,   60,   60,   0,    0, -60, -60,   0,   20,  20')),

          (_T('18,  -18,   36,  210,  162, -36,    0,   0,   0,   0,  -18,  18'),  # PYCHOK indent
          _T('-36,   36, -165,   45,  141, -21,    0,   0,   0,   0,   36, -36'),
         _T('-122,   -2,  -27, -111,  -75,  27,   62, 124, 124,  62,  -64,   2'),
            _T('0,    0,   93,  -93,  -93,  93,    0,   0,   0,   0,    0,   0'),
          _T('120, -120,  147,  -57, -129,  39,    0,   0,   0,   0,   66, -66'),  # PYCHOK indent
          _T('135,   51,   -9, -192, -180,   9,   31,  62,  62,  31,   51, -51'),
            _T0s12,
            _T('0,    0,  -93,   93,   93, -93,    0,   0,   0,   0,    0,   0'),
          _T('-84,   84,   18,   12,  -12, -18,    0,   0,   0,   0, -102, 102'),
          _T('-31,  -31,    0,   93,   93,   0,  -31, -62, -62, -31,   31,  31')))

    _BT = (_T('0, 0'),  # bilinear 4-tuple [i, j] indices
           _T('1, 0'),
           _T('0, 1'),
           _T('1, 1'))

    _CM = (_T(' 0, -1'),  # 10x12 cubic matrix [i, j] indices
           _T(' 1, -1'),
           _T('-1,  0'),
           _T(' 0,  0'),
           _T(' 1,  0'),
           _T(' 2,  0'),
           _T('-1,  1'),
           _T(' 0,  1'),
           _T(' 1,  1'),
           _T(' 2,  1'),
           _T(' 0,  2'),
           _T(' 1,  2'))

#   _cropped =  None
    _endian  = '>H'   # struct.unpack 1 ushort (big endian, unsigned short)
    _4endian = '>4H'  # struct.unpack 4 ushorts
    _Rendian =  NN    # struct.unpack a row of ushorts
#   _highest = (-8.4,   147.367, 85.839) if egm2008-1.pgm else (
#              (-8.167, 147.25,  85.422) if egm96-5.pgm else
#              (-4.5,   148.75,  81.33))  # egm84-15.pgm
    _iscipy  =  False
#   _lowest  = (4.7,   78.767, -106.911) if egm2008-1.pgm else (
#              (4.667, 78.833, -107.043) if egm96-5.pgm else
#              (4.75,  79.25,  -107.34))  # egm84-15.pgm
    _mean    = _F(-1.317)  # from egm2008-1, -1.438 egm96-5, -0.855 egm84-15
    _nBytes  =  None  # not applicable
    _nterms  =  len(_C3[0])  # columns length, number of rows
    _smooth  =  None  # not applicable
    _stdev   = _F(29.244)  # from egm2008-1, 29.227 egm96-5, 29.183 egm84-15
    _u2B     = _calcsize(_endian)  # pixelsize_ in bytes
    _4u2B    = _calcsize(_4endian)  # 4 pixelsize_s in bytes
    _Ru2B    =  0  # row of pixelsize_s in bytes
    _yx_hits =  0  # cache hits
    _yx_i    = ()  # cached (y, x) indices
    _yx_t    = ()  # cached 4- or 10-tuple for _ev2k resp. _ev3k

    def __init__(self, egm_pgm, crop=None, datum=_WGS84, kind=3, **name_smooth):
        '''New L{GeoidKarney} interpolator.

           @arg egm_pgm: An U{EGM geoid dataset<https://GeographicLib.SourceForge.io/
                         C++/doc/geoid.html#geoidinst>} file name (C{egm*.pgm}), see
                         note below.
           @kwarg crop: Optional box to limit geoid locations, a 4-tuple (C{south,
                        west, north, east}), 2-tuple (C{(south, west), (north, east)})
                        with 2 C{degrees90} lat- and C{degrees180} longitudes or as
                        2-tuple (C{LatLonSW, LatLonNE}) of C{LatLon} instances.
           @kwarg datum: Optional grid datum (C{Datum}, L{Ellipsoid}, L{Ellipsoid2} or
                         L{a_f2Tuple}), overriding C{WGS84}.
           @kwarg kind: Interpolation order (C{int}), 2 for C{bilinear} or 3 for C{cubic}.
           @kwarg name_smooth: Optional geoid C{B{name}=NN} (C{str}) and UNSUPPORTED
                       keyword argument C{B{smooth}}, use C{B{smooth}=None} to ignore.

           @raise GeoidError: EGM dataset B{C{egm_pgm}} issue or invalid B{C{crop}},
                              B{C{kind}} or B{C{smooth}}.

           @raise TypeError: Invalid B{C{datum}}.

           @see: Class L{GeoidPGM} and function L{egmGeoidHeights}.

           @note: Geoid file B{C{egm_pgm}} remains open and I{must be closed} by calling
                  method C{close} or by using C{with B{GeoidKarney}(...) as ...:} context.
        '''
        smooth, name = _xkwds_pop2(name_smooth, smooth=None)
        if smooth is not None:
            raise GeoidError(smooth=smooth, txt_not_=_supported_)

        if _isin(kind, 2):
            self._ev2d = self._ev2k  # see ._ev_name
        elif not _isin(kind, 3):
            raise GeoidError(kind=kind)

        self._egm = g =  self._open(egm_pgm, datum, kind, _name__(**name), None)
        self._pgm = p = _PGM(g, pgm=egm_pgm, itemsize=self.u2B, sizeB=self.sizeB)

        self._Rendian =  self._4endian.replace(_4_, str(p.nlon))
        self._Ru2B    = _calcsize(self._Rendian)

        self._knots  = p.knots  # grid knots
        self._lon_of = float(p.flon)  # forward offset
        self._lon_og = float(p.glon)  # reverse offset
        # set earth (lat, lon) limits (s, w, n, e)
        self._lat_lo, self._lon_lo, \
        self._lat_hi, self._lon_hi = self._swne(crop if crop else p.crop4)
        self._cropped = bool(crop)

    def _c0c3v(self, y, x):
        # get the common denominator, the 10x12 cubic matrix and
        # the 12 cubic v-coefficients around geoid index (y, x)
        p = self._pgm
        if 0 < x < (p.nlon - 2) and 0 < y < (p.nlat - 2):
            # read 4x4 ushorts, drop the 4 corners
            S = _MODS.os.SEEK_SET
            e =  self._4endian
            g =  self._egm
            n =  self._4u2B
            R =  self._Ru2B
            b =  self._seek(y - 1, x - 1)
            v = _unpack(e, g.read(n))[1:3]
            b += R
            g.seek(b, S)
            v += _unpack(e, g.read(n))
            b += R
            g.seek(b, S)
            v += _unpack(e, g.read(n))
            b += R
            g.seek(b, S)
            v += _unpack(e, g.read(n))[1:3]
            j = 1

        else:  # likely some wrapped y and/or x's
            v = self._raws(y, x, GeoidKarney._CM)
            j = 0 if y < 1 else (1 if y < (p.nlat - 2) else 2)

        return GeoidKarney._C0[j], GeoidKarney._C3[j], v

    @Property_RO
    def dtype(self):
        '''Get the geoid's grid data type (C{str}).
        '''
        return 'ushort'

    def _ev(self, lat, lon):  # PYCHOK expected
        # interpolate the geoid height at grid (lat, lon)
        fy, fx = self._g2yx2(lat, lon)
        y, x = int(_floor(fy)), int(_floor(fx))
        fy -= y
        fx -= x
        H  = self._ev2d(fy, fx, y, x)  # PYCHOK ._ev2k or ._ev3k
        H *= self._pgm.Scale   # H.fmul(self._pgm.Scale)
        H += self._pgm.Offset  # H.fadd(self._pgm.Offset)
        return H.fsum()  # float(H)

    def _ev2k(self, fy, fx, *yx):
        # compute the bilinear 4-tuple and interpolate raw H
        if self._yx_i == yx:
            self._yx_hits += 1
        else:
            y, x = self._yx_i = yx
            self._yx_t = self._raws(y, x, GeoidKarney._BT)
        t  =  self._yx_t
        v  = _1_0, -fx, fx
        H  = Fdot(v, t[0], t[0], t[1]).fmul(_1_0 - fy)  # c = a * (1 - fy)
        H += Fdot(v, t[2], t[2], t[3]).fmul(fy)  # c += b * fy
        return H

    def _ev3k(self, fy, fx, *yx):
        # compute the cubic 10-tuple and interpolate raw H
        if self._yx_i == yx:
            self._yx_hits += 1
        else:
            c0, c3, v  = self._c0c3v(*yx)
            # assert len(c3) == self._nterms
            self._yx_t = tuple(fdot(v, *r3) / c0 for r3 in c3)
            self._yx_i = yx
        # GeographicLib/Geoid.cpp Geoid::height(lat, lon) ...
        # real h = t[0] + fx * (t[1] + fx * (t[3] + fx * t[6])) +
        #                 fy * (t[2] + fx * (t[4] + fx * t[7]) +
        #                 fy * (t[5] + fx *  t[8] + fy * t[9]));
        t  =  self._yx_t
        v  = _1_0, fx, fy
        H  = Fdot(v, t[5], t[8], t[9])
        H *= fy
        H += Fhorner(fx, t[2], t[4], t[7])
        H *= fy
        H += Fhorner(fx, t[0], t[1], t[3], t[6])
        return H

    _ev2d = _ev3k  # overriden for kind=2, see ._ev_name

    def _g2ll2(self, lat, lon):
        # convert grid (lat, lon) to earth (lat, lon), uncropped
        while lon > _180_0:
            lon -= _360_0
        return lat, lon

    def _g2yx2(self, lat, lon):
        # convert grid (lat, lon) to grid (y, x) indices
        p = self._pgm
        # note, slat = +90, rlat < 0 makes y >=0
        return ((lat - p.slat) * p.rlat), ((lon - p.wlon) * p.rlon)

    def _gyx2g2(self, y, x):
        # convert grid (y, x) indices to grid (lat, lon)
        p = self._pgm
        return (p.slat + p.dlat * y), (p.wlon + p.dlon * x)

    @Property_RO
    def _highest_ltd(self):
        '''(INTERNAL) Cache for C{.highest}.
        '''
        return self._LL3T(self._llh3minmax(True, -12, -4), name__=self.highest)

    def highest(self, LatLon=None, full=False):  # PYCHOK full
        '''Return the location and largest height of this geoid.

           @kwarg LatLon: Optional class to return the location and height
                          (C{LatLon}) or C{None}.
           @kwarg full: Search the full or limited latitude range (C{bool}).

           @return: If C{B{LatLon} is None}, a L{LatLon3Tuple}C{(lat, lon,
                    height)} otherwise a B{C{LatLon}} instance with the lat-,
                    longitude and geoid height of the highest grid location.
        '''
        llh = self._highest if full or self.cropped else self._highest_ltd
        return self._llh3LL(llh, LatLon)

    def _lat2y2(self, lat2):
        # convert earth lat(s) to min and max grid y indices
        ys, m = [], self._pgm.nlat - 1
        for lat in lat2:
            y, _ = self._g2yx2(*self._ll2g2(lat, 0))
            ys.append(max(min(int(y), m), 0))
        return min(ys), max(ys) + 1

    def _ll2g2(self, lat, lon):
        # convert earth (lat, lon) to grid (lat, lon), uncropped
        while lon < 0:
            lon += _360_0
        return lat, lon

    def _llh3minmax(self, highest, *lat2):
        # find highest or lowest, takes 10+ secs for egm2008-1.pgm geoid
        # (Python 2.7.16, macOS 10.13.6 High Sierra, iMac 3 GHz Core i3)
        if highest:
            def _mt(r, h):
                m = max(r)
                return m, (m > h)

        else:  # lowest
            def _mt(r, h):  # PYCHOK redef
                m = min(r)
                return m, (m < h)

        y = x = 0
        h = self._raw(y, x)
        for j, r in self._raw2(*lat2):
            m, t = _mt(r, h)
            if t:
                h, y, x = m, j, r.index(m)
        h *= self._pgm.Scale
        h += self._pgm.Offset
        return self._g2ll2(*self._gyx2g2(y, x)) + (h,)

    @Property_RO
    def _lowest_ltd(self):
        '''(INTERNAL) Cache for C{.lowest}.
        '''
        return self._LL3T(self._llh3minmax(False, 0, 8), name__=self.lowest)

    def lowest(self, LatLon=None, full=False):  # PYCHOK full
        '''Return the location and lowest height of this geoid.

           @kwarg LatLon: Optional class to return the location and height
                          (C{LatLon}) or C{None}.
           @kwarg full: Search the full or limited latitude range (C{bool}).

           @return: If C{B{LatLon} is None}, a L{LatLon3Tuple}C{(lat, lon,
                    height)} otherwise a B{C{LatLon}} instance with the lat-,
                    longitude and geoid height of the lowest grid location.
        '''
        llh = self._lowest if full or self.cropped else self._lowest_ltd
        return self._llh3LL(llh, LatLon)

    def _raw(self, y, x):
        # get the ushort geoid height at geoid index (y, x),
        # like GeographicLib/Geoid.hpp real rawval(is, iy)
        p = self._pgm
        if x < 0:
            x += p.nlon
        elif x >= p.nlon:
            x -= p.nlon
        h = p.nlon // 2
        if y < 0:
            y = -y
        elif y >= p.nlat:
            y = (p.nlat - 1) * 2 - y
        else:
            h = 0
        x += h if x < h else -h
        self._seek(y, x)
        h = _unpack(self._endian, self._egm.read(self._u2B))
        return h[0]

    def _raws(self, y, x, ijs):
        # get bilinear 4-tuple or 10x12 cubic matrix
        return tuple(self._raw(y + j, x + i) for i, j in ijs)

    def _raw2(self, *lat2):
        # yield a 2-tuple (y, ushorts) for each row or for
        # the rows between two (or more) earth lat values
        p = self._pgm
        g = self._egm
        e = self._Rendian
        n = self._Ru2B
        # min(lat2) <= lat <= max(lat2) or 0 <= y < p.nlat
        s, t = self._lat2y2(lat2) if lat2 else (0, p.nlat)
        self._seek(s, 0)  # to start of row s
        for y in range(s, t):
            yield y, _unpack(e, g.read(n))

    def _seek(self, y, x):
        # position geoid to grid index (y, x)
        p, g = self._pgm, self._egm
        if g:
            b = p.skip + (y * p.nlon + x) * self._u2B
            g.seek(b, _MODS.os.SEEK_SET)
            return b  # position
        raise GeoidError('closed file', txt=repr(p.egm))  # IOError


class GeoidPGM(_GeoidBase):
    '''Geoid height interpolator for I{Karney}'s U{GeographicLib Earth
       Gravitational Model (EGM)<https://GeographicLib.SourceForge.io/C++/doc/geoid.html>}
       geoid U{egm*.pgm<https://GeographicLib.SourceForge.io/C++/doc/geoid.html#geoidinst>}
       datasets but based on C{SciPy} U{RectBivariateSpline<https://docs.SciPy.org/doc/scipy/
       reference/generated/scipy.interpolate.RectBivariateSpline.html>}, U{bisplrep/-ev
       <https://docs.SciPy.org/doc/scipy/reference/generated/scipy.interpolate.bisplrep.html>}
       or U{interp2d<https://docs.SciPy.org/doc/scipy/reference/generated/scipy.interpolate.
       interp2d.html>} interpolation.

       Use any of the U{egm84-, egm96- or egm2008-*.pgm <https://GeographicLib.SourceForge.io/
       C++/doc/geoid.html#geoidinst>} datasets.  However, unless cropped, an entire C{egm*.pgm}
       dataset is loaded into the C{SciPy} interpolator and converted from 2-byte C{int} to
       8-byte C{dtype float64}.  Therefore, internal memory usage is 4x the U{egm*.pgm
       <https://GeographicLib.SourceForge.io/C++/doc/geoid.html#geoidinst>} file size and may
       exceed the available memory, especially with 32-bit Python, see properties C{.nBytes}
       and C{.sizeB}.
    '''
    _cropped =  False
    _endian  = '>u2'

    def __init__(self, egm_pgm, crop=None, datum=_WGS84, kind=3, smooth=0, **name):
        '''New L{GeoidPGM} interpolator.

           @arg egm_pgm: An U{EGM geoid dataset<https://GeographicLib.SourceForge.io/
                         C++/doc/geoid.html#geoidinst>} file name (C{egm*.pgm}).
           @kwarg crop: Optional box to crop B{C{egm_pgm}}, a 4-tuple (C{south, west,
                        north, east}) or 2-tuple (C{(south, west), (north, east)}),
                        in C{degrees90} lat- and C{degrees180} longitudes or a 2-tuple
                        (C{LatLonSW, LatLonNE}) of C{LatLon} instances.
           @kwarg datum: Optional grid datum (L{Datum}, L{Ellipsoid}, L{Ellipsoid2} or
                         L{a_f2Tuple}), overriding C{WGS84}.
           @kwarg kind: C{scipy.interpolate} order (C{int}), use 1..5 for U{RectBivariateSpline
                        <https://docs.SciPy.org/doc/scipy/reference/generated/scipy.interpolate.
                        RectBivariateSpline.html>} or -1, -3 or -5 for U{bisplrep/-ev<https://
                        docs.SciPy.org/doc/scipy/reference/generated/scipy.interpolate.bisplrep.html>}
                        or U{interp2d<https://docs.SciPy.org/doc/scipy/reference/generated/scipy.
                        interpolate.interp2d.html>} C{linear}, C{cubic} respectively C{quintic},
                        see note for more details.
           @kwarg smooth: Smoothing factor for C{B{kind}=1..5} only (C{int}).
           @kwarg name: Optional geoid C{B{name}=NN} (C{str}).

           @raise GeoidError: EGM dataset B{C{egm_pgm}} issue or invalid B{C{crop}}, B{C{kind}}
                              or B{C{smooth}}.

           @raise ImportError: Package C{numpy} or C{scipy} not found or not installed.

           @raise LenError: EGM dataset B{C{egm_pgm}} axis mismatch.

           @raise SciPyError: A C{scipy} issue.

           @raise SciPyWarning: A C{scipy} warning as exception.

           @raise TypeError: Invalid B{C{datum}} or unexpected argument.

           @note: Specify C{B{kind}=-1, -3 or -5} to use C{scipy.interpolate.interp2d}
                  before or C{scipy.interpolate.bisplrep/-ev} since C{Scipy} version 1.14.

           @note: The U{GeographicLib egm*.pgm<https://GeographicLib.SourceForge.io/C++/doc/
                  geoid.html#geoidinst>} file sizes are based on a 2-byte C{int} height
                  converted to 8-byte C{dtype float64} for C{scipy} interpolators.  Therefore,
                  internal memory usage is 4 times the C{egm*.pgm} file size and may exceed
                  the available memory, especially with 32-bit Python.  To reduce memory
                  usage, set keyword argument B{C{crop}} to the region of interest.  For
                  example C{B{crop}=(20, -125, 50, -65)} covers the U{conterminous US
                  <https://www.NGS.NOAA.gov/GEOID/GEOID12B/maps/GEOID12B_CONUS_grids.png>}
                  (CONUS), less than 3% of the entire C{egm2008-1.pgm} dataset.

           @see: Class L{GeoidKarney} and function L{egmGeoidHeights}.
        '''
        np = self.numpy
        self._u2B = np.dtype(self.endian).itemsize

        g = self._open(egm_pgm, datum, kind, _name__(**name), smooth)
        self._pgm = p = _PGM(g, pgm=egm_pgm, itemsize=self.u2B, sizeB=self.sizeB)
        if crop:
            g = p._cropped(g, abs(kind) + 1, *self._swne(crop))
            if _MODS.internals._version2(np.__version__) < (1, 9):
                g = open(g.name, _rb_)  # reopen tempfile for numpy 1.8.0-
            self._cropped = True
        try:
            # U{numpy dtype formats are different from Python struct formats
            # <https://docs.SciPy.org/doc/numpy-1.15.0/reference/arrays.dtypes.html>}
            # read all heights, skipping the PGM header lines, converted to float
            hs = self._load(g, self.endian, p.knots, p.skip).reshape(p.nlat, p.nlon) * p.Scale
            if p.Offset:  # offset
                hs = p.Offset + hs
            if p.dlat < 0:  # flip the rows
                hs = np.flipud(hs)
            _GeoidBase.__init__(self, hs, p)
        except Exception as x:
            raise _SciPyIssue(x, _in_, repr(egm_pgm))
        finally:
            g.close()

    def _g2ll2(self, lat, lon):
        # convert grid (lat, lon) to earth (lat, lon), un-/cropped
        if self._cropped:
            lon -= self._lon_of
        else:
            while lon > _180_0:
                lon -= _360_0
        return lat, lon

    def _ll2g2(self, lat, lon):
        # convert earth (lat, lon) to grid (lat, lon), un-/cropped
        if self._cropped:
            lon += self._lon_of
        else:
            while lon < 0:
                lon += _360_0
        return lat, lon

    if _FOR_DOCS:
        __call__ = _GeoidBase.__call__
        height   = _GeoidBase.height


class _Gpars(_Named):
    '''(INTERNAL) Basic geoid parameters.
    '''
    # interpolator parameters
    dlat = 0  # +/- latitude resolution in C{degrees}
    dlon = 0  # longitude resolution in C{degrees}
    nlat = 1  # number of latitude knots (C{int})
    nlon = 0  # number of longitude knots (C{int})
    rlat = 0  # +/- latitude resolution in C{float}, 1 / .dlat
    rlon = 0  # longitude resolution in C{float}, 1 / .dlon
    slat = 0  # nothern- or southern most latitude (C{degrees90})
    wlon = 0  # western-most longitude in Eastern lon (C{degrees360})

    flon = 0  # forward, earth to grid longitude offset
    glon = 0  # reverse, grid to earth longitude offset

    knots = 0  # number of knots, nlat * nlon (C{int})
    skip  = 0  # header bytes to skip (C{int})

    def __repr__(self):
        t = _COMMASPACE_.join(pairs((a, getattr(self, a)) for
                                     a in dir(self.__class__)
                                       if a[:1].isupper()))
        return _COLONSPACE_(self, t)

    def __str__(self):
        return Fmt.PAREN(self.classname, repr(self.name))


class _PGM(_Gpars):
    '''(INTERNAL) Parse an C{egm*.pgm} geoid dataset file.

       # Geoid file in PGM format for the GeographicLib::Geoid class
       # Description WGS84 EGM96, 5-minute grid
       # URL https://Earth-Info.NGA.mil/GandG/wgs84/gravitymod/egm96/egm96.html
       # DateTime 2009-08-29 18:45:03
       # MaxBilinearError 0.140
       # RMSBilinearError 0.005
       # MaxCubicError 0.003
       # RMSCubicError 0.001
       # Offset -108
       # Scale 0.003
       # Origin 90N 0E
       # AREA_OR_POINT Point
       # Vertical_Datum WGS84
       <width> <height>
       <pixel>
       ...
    '''
    crop4 = ()  # 4-tuple (C{south, west, north, east}).
    egm   = None
    glon  = 180  # reverse offset, uncropped
#   pgm   = NN   # name
    sizeB = 0
    u2B   = 2  # item size of grid height (C{int}).

    @staticmethod
    def _llstr2floats(latlon):
        # llstr to (lat, lon) floats
        lat, lon = latlon.split()
        return _MODS.dms.parseDMS2(lat, lon)

    # PGM file attributes, CamelCase but not .istitle()
    AREA_OR_POINT    =  str
    DateTime         =  str
    Description      =  str  # 'WGS84 EGM96, 5-minute grid'
    Geoid            =  str  # 'file in PGM format for the GeographicLib::Geoid class'
    MaxBilinearError =  float
    MaxCubicError    =  float
    Offset           =  float
    Origin           = _llstr2floats
    Pixel            =  0
    RMSBilinearError =  float
    RMSCubicError    =  float
    Scale            =  float
    URL              =  str  # 'https://Earth-Info.NGA.mil/GandG/wgs84/...'
    Vertical_Datum   =  str

    def __init__(self, g, pgm=NN, itemsize=0, sizeB=0):  # MCCABE 22
        '''(INTERNAL) New C{_PGM} parsed C{egm*.pgm} geoid dataset.
        '''
        self.name = pgm  # geoid file name
        if itemsize:
            self._u2B = itemsize
        if sizeB:
            self.sizeB = sizeB

        t = g.readline()  # make sure newline == '\n'
        if t != b'P5\n' and t.strip() != b'P5':
            raise self._Errorf(_format_, _header_, t)

        while True:  # read all # Attr ... lines,
            try:  # ignore empty ones or comments
                t = g.readline().strip()
                if t.startswith(_bHASH_):
                    t = t.lstrip(_bHASH_).lstrip()
                    a, v = map(_ub2str, t.split(None, 1))
                    f = getattr(_PGM, a, None)
                    if callable(f) and a[:1].isupper():
                        setattr(self, a, f(v))
                elif t:
                    break
            except (TypeError, ValueError):
                raise self._Errorf(_format_, 'Attr', t)
        else:  # should never get here
            raise self._Errorf(_format_, _header_, g.tell())

        try:  # must be (even) width and (odd) height
            nlon, nlat = map(int, t.split())
            if nlon < 2 or nlon > (360 * 60) or isodd(nlon) or \
               nlat < 2 or nlat > (181 * 60) or not isodd(nlat):
                raise ValueError
        except (TypeError, ValueError):
            raise self._Errorf(_format_, _SPACE_(_width_, _height_), t)

        try:  # must be 16 bit pixel height
            t = g.readline().strip()
            self.Pixel = int(t)
            if not 255 < self.Pixel < 65536:  # >u2 or >H only
                raise ValueError
        except (TypeError, ValueError):
            raise self._Errorf(_format_, 'pixel', t)

        for a in dir(_PGM):  # set undefined # Attr ... to None
            if a[:1].isupper() and callable(getattr(self, a)):
                setattr(self, a, None)

        if self.Origin is None:
            raise self._Errorf(_format_, 'Origin', self.Origin)
        if self.Offset is None or self.Offset > 0:
            raise self._Errorf(_format_, 'Offset', self.Offset)
        if self.Scale is None or self.Scale < EPS:
            raise self._Errorf(_format_, 'Scale', self.Scale)

        self.skip = g.tell()
        self.knots = nlat * nlon

        self.nlat, self.nlon = nlat, nlon
        self.slat, self.wlon = self.Origin
        # note, negative .dlat and .rlat since rows
        # are from .slat 90N down in decreasing lat
        self.dlat, self.dlon = _180_0 / (1 - nlat), _360_0 / nlon
        self.rlat, self.rlon = (1 - nlat) / _180_0, nlon / _360_0

        # grid corners in earth (lat, lon), .slat = 90, .dlat < 0
        n = float(self.slat)
        s = n + self.dlat * (nlat - 1)
        w = self.wlon - self.glon
        e = w + self.dlon * nlon
        self.crop4 = s, w, n, e

        n = self.sizeB - self.skip
        if n > 0 and n != (self.knots * self.u2B):
            raise self._Errorf('%s(%s x %s != %s)', _assert_, nlat, nlon, n)

    def _cropped(self, g, k1, south, west, north, east):  # MCCABE 15
        '''Crop the geoid to (south, west, north, east) box.
        '''
        # flon offset for both west and east
        f = 360 if west < 0 else 0
        # earth (lat, lon) to grid indices (y, x),
        # note y is decreasing, i.e. n < s
        s, w = self._lle2yx2(south, west, f)
        n, e = self._lle2yx2(north, east, f)
        s += 1  # s > n
        e += 1  # e > w

        hi, wi = self.nlat, self.nlon
        # handle special cases
        if (s - n) > hi:
            n, s = 0, hi  # entire lat range
        if (e - w) > wi:
            w, e, f = 0, wi, 180  # entire lon range
        if s == hi and w == n == 0 and e == wi:
            return g  # use entire geoid as-is

        if (e - w) < k1 or (s - n) < (k1 + 1):
            raise self._Errorf(_format_, 'swne', (north - south, east - west))

        if e > wi > w:  # wrap around
            # read w..wi and 0..e
            r, p = (wi - w), (e - wi)
        elif e > w:
            r, p = (e - w), 0
        else:
            raise self._Errorf('%s(%s < %s)', _assert_, w, e)

        # convert to bytes
        r *= self.u2B
        p *= self.u2B
        q = wi * self.u2B  # stride
        # number of rows and cols to skip from
        # the original (.slat, .wlon) origin
        z = self.skip + (n * wi + w) * self.u2B
        # sanity check
        if r < 2 or p < 0 or q < 2 or z < self.skip \
                                   or z > self.sizeB:
            raise self._Errorf(_format_, _assert_, (r, p, q, z))

        # can't use _BytesIO since numpy
        # needs .fileno attr in .fromfile
        t, c = 0, self._tmpfile()
        # reading (s - n) rows, forward
        for y in range(n, s):  # PYCHOK y unused
            g.seek(z, _MODS.os.SEEK_SET)
            # Python 2 tmpfile.write returns None
            t += c.write(g.read(r)) or r
            if p:  # wrap around to start of row
                g.seek(-q, _MODS.os.SEEK_CUR)
                # assert(g.tell() == (z - w * self.u2B))
                # Python 2 tmpfile.write returns None
                t += c.write(g.read(p)) or p
            z += q
        c.flush()
        g.close()

        s -= n  # nlat
        e -= w  # nlon
        k = s * e  # knots
        z = k * self.u2B
        if t != z:
            raise self._Errorf('%s(%s != %s) %s', _assert_, t, z, self)

        # update the _Gpars accordingly, note attributes
        # .dlat, .dlon, .rlat and .rlon remain unchanged
        self.slat += n * self.dlat
        self.wlon += w * self.dlon
        self.nlat = s
        self.nlon = e
        self.flon = self.glon = f

        self.crop4 = south, west, north, east
        self.knots = k
        self.skip  = 0  # no header lines in c

        c.seek(0, _MODS.os.SEEK_SET)
        # c = open(c.name, _rb_)  # reopen for numpy 1.8.0-
        return c

    def _Errorf(self, fmt, *args):  # PYCHOK no cover
        t = fmt % args
        e = self.pgm or NN
        if e:
            t = _SPACE_(t, _in_, repr(e))
        return PGMError(t)

    def _lle2yx2(self, lat, lon, flon):
        # earth (lat, lon) to grid indices (y, x)
        # with .dlat decreasing from 90N .slat
        lat -= self.slat
        lon += flon - self.wlon
        return (min(self.nlat - 1, max(0, int(lat * self.rlat))),
                                   max(0, int(lon * self.rlon)))

    def _tmpfile(self):
        # create a tmpfile to hold the cropped geoid grid
        try:
            from tempfile import NamedTemporaryFile as tmpfile
        except ImportError:  # Python 2.7.16-
            from _MODS.os import tmpfile
        t = _MODS.os.path.basename(self.pgm)
        t = _MODS.os.path.splitext(t)[0]
        f =  tmpfile(mode='w+b', prefix=t or 'egm')
        f.seek(0, _MODS.os.SEEK_SET)  # force overwrite
        return f

    @Property_RO
    def pgm(self):
        '''Get the geoid file name (C{str}).
        '''
        return self.name


class PGMError(GeoidError):
    '''An issue while parsing or cropping an C{egm*.pgm} geoid dataset.
    '''
    pass


def egmGeoidHeights(GeoidHeights_dat):
    '''Generate geoid U{egm*.pgm<https://GeographicLib.SourceForge.io/
       C++/doc/geoid.html#geoidinst>} height tests from U{GeoidHeights.dat
       <https://SourceForge.net/projects/geographiclib/files/testdata/>}
       U{Test data for Geoids<https://GeographicLib.SourceForge.io/C++/doc/
       geoid.html#testgeoid>}.

       @arg GeoidHeights_dat: The un-gz-ed C{GeoidHeights.dat} file
                              (C{str} or C{file} handle).

       @return: For each test, yield a L{GeoidHeight5Tuple}C{(lat, lon,
                egm84, egm96, egm2008)}.

       @raise GeoidError: Invalid B{C{GeoidHeights_dat}}.

       @note: Function L{egmGeoidHeights} is used to test the geoids
              L{GeoidKarney} and L{GeoidPGM}, see PyGeodesy module
              C{test/testGeoids.py}.
    '''
    dat = GeoidHeights_dat
    if isinstance(dat, bytes):
        dat = _BytesIO(dat)

    try:
        dat.seek(0, _MODS.os.SEEK_SET)  # reset
    except AttributeError as x:
        raise GeoidError(GeoidHeights_dat=type(dat), cause=x)

    for t in dat.readlines():
        t = t.strip()
        if t and not t.startswith(_bHASH_):
            lat, lon, egm84, egm96, egm2008 = map(float, t.split())
            while lon > _180_0:  # EasternLon to earth lon
                lon -= _360_0
            yield GeoidHeight5Tuple(lat, lon, egm84, egm96, egm2008)


__all__ += _ALL_DOCS(_GeoidBase)

if __name__ == _DMAIN_:  # MCCABE 14

    from pygeodesy.internals import printf, _secs2str, _versions,  _sys
    from time import time

    _crop     = {}
    _GeoidEGM = GeoidKarney
    _kind     = 3

    geoids = _sys.argv[1:]
    while geoids:
        G = geoids.pop(0)
        g = G.lower()

        if '-crop'.startswith(g):
            _crop = dict(crp=(20, -125, 50, -65))  # CONUS

        elif '-egm96'.startswith(g):
            _GeoidEGM = GeoidEGM96

        elif '-karney'.startswith(g):
            _GeoidEGM = GeoidKarney

        elif '-kind'.startswith(g):
            _kind = int(geoids.pop(0))

        elif '-pgm'.startswith(g):
            _GeoidEGM = GeoidPGM

        elif _isin(g[-4:], '.pgm', '.grd'):
            g = _GeoidEGM(G, kind=_kind, **_crop)
            t =  time()
            _ =  g.highest()
            t = _secs2str(time() - t)
            printf('%s: %s (%s)', g.toStr(), t, _versions(), nl=1, nt=1)
            t =  g.pgm
            if t:
                printf(repr(t), nt=1)
            # <https://GeographicLib.SourceForge.io/cgi-bin/GeoidEval>:
            # The height of the EGM96 geoid at Timbuktu
            #    echo 16:46:33N 3:00:34W | GeoidEval
            #    => 28.7068 -0.02e-6 -1.73e-6
            # The 1st number is the height of the geoid, the 2nd and
            # 3rd are its slopes in northerly and easterly direction
            t = 'Timbuktu %s' % (g,)
            k = {'egm84-15.pgm':  '31.2979',
                 'egm96-5.pgm':   '28.7067',
                 'egm2008-1.pgm': '28.7880'}.get(g.name.lower(), '28.7880')
            ll = _MODS.dms.parseDMS2('16:46:33N', '3:00:34W', sep=':')
            for ll in (ll, (16.776, -3.009),):
                try:
                    h, ll = g.height(*ll), fstr(ll, prec=6)
                    printf('%s.height(%s): %.4F vs %s', t, ll, h, k)
                except (GeoidError, RangeError) as x:
                    printf(_COLONSPACE_(t, str(x)))

        elif _isin(g[-4:], '.bin'):
            g = GeoidG2012B(G, kind=_kind)
            printf(g.toStr())

        else:
            raise GeoidError(grid=repr(G))

_I = int  # PYCHOK unused _I
del _intCs, _T, _T0s12  # trash ints cache and map


# <https://GeographicLib.SourceForge.io/cgi-bin/GeoidEval>
# _lowerleft = -90, -179, -30.1500  # egm2008-1.pgm
# _lowerleft = -90, -179, -29.5350  # egm96-5.pgm
# _lowerleft = -90, -179, -29.7120  # egm84-15.pgm

# _center = 0, 0, 17.2260  # egm2008-1.pgm
# _center = 0, 0, 17.1630  # egm96-5.pgm
# _center = 0, 0, 18.3296  # egm84-15.pgm

# _upperright = 90, 180, 14.8980  # egm2008-1.pgm
# _upperright = 90, 180, 13.6050  # egm96-5.pgm
# _upperright = 90, 180, 13.0980  # egm84-15.pgm


# % python3.12 -m pygeodesy.geoids -egm96 ../testGeoids/WW15MGH.GRD
#
# GeoidEGM96('WW15MGH.GRD'): lowerleft(-90.0, -180.0, -29.534), upperright(90.0, 180.25, 13.606), center(0.0, 0.125, 17.125), highest(-8.25, -32.75, 85.391), lowest(4.75, -101.25, -106.991): 1.267 ms (pygeodesy 24.12.24 Python 3.12.7 64bit arm64 macOS 14.6.1)
#
# Timbuktu GeoidEGM96('WW15MGH.GRD').height(16.775833, -3.009444): 28.7073 vs 28.7880
# Timbuktu GeoidEGM96('WW15MGH.GRD').height(16.776, -3.009): 28.7072 vs 28.7880


# % python3.12 -m pygeodesy.geoids -Karney ../testGeoids/egm*.pgm
#
# GeoidKarney('egm2008-1.pgm'): lowerleft(-90.0, -180.0, -30.15), upperright(90.0, 180.0, 14.898), center(0.0, 0.0, 17.226), highest(-8.4, 147.367, 85.839), lowest(4.7, 78.767, -106.911): 204.334 ms (pygeodesy 24.8.24 Python 3.12.5 64bit arm64 macOS 14.6.1)
#
# _PGM('../testGeoids/egm2008-1.pgm'): AREA_OR_POINT='Point', DateTime='2009-08-31 06:54:00', Description='WGS84 EGM2008, 1-minute grid', Geoid='file in PGM format for the GeographicLib::Geoid class', MaxBilinearError=0.025, MaxCubicError=0.003, Offset=-108.0, Origin=LatLon2Tuple(lat=90.0, lon=0.0), Pixel=65535, RMSBilinearError=0.001, RMSCubicError=0.001, Scale=0.003, URL='http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008', Vertical_Datum='WGS84'
#
# Timbuktu GeoidKarney('egm2008-1.pgm').height(16.775833, -3.009444): 28.7881 vs 28.7880
# Timbuktu GeoidKarney('egm2008-1.pgm').height(16.776, -3.009): 28.7880 vs 28.7880
#
# GeoidKarney('egm84-15.pgm'): lowerleft(-90.0, -180.0, -29.712), upperright(90.0, 180.0, 13.098), center(0.0, 0.0, 18.33), highest(-4.5, 148.75, 81.33), lowest(4.75, 79.25, -107.34): 1.007 ms (pygeodesy 24.8.24 Python 3.12.5 64bit arm64 macOS 14.6.1)
#
# _PGM('../testGeoids/egm84-15.pgm'): AREA_OR_POINT='Point', DateTime='2009-08-29 18:45:02', Description='WGS84 EGM84, 15-minute grid', Geoid='file in PGM format for the GeographicLib::Geoid class', MaxBilinearError=0.413, MaxCubicError=0.02, Offset=-108.0, Origin=LatLon2Tuple(lat=90.0, lon=0.0), Pixel=65535, RMSBilinearError=0.018, RMSCubicError=0.001, Scale=0.003, URL='http://earth-info.nga.mil/GandG/wgs84/gravitymod/wgs84_180/wgs84_180.html', Vertical_Datum='WGS84'
#
# Timbuktu GeoidKarney('egm84-15.pgm').height(16.775833, -3.009444): 31.2983 vs 31.2979
# Timbuktu GeoidKarney('egm84-15.pgm').height(16.776, -3.009): 31.2979 vs 31.2979
#
# GeoidKarney('egm96-5.pgm'): lowerleft(-90.0, -180.0, -29.535), upperright(90.0, 180.0, 13.605), center(0.0, 0.0, 17.163), highest(-8.167, 147.25, 85.422), lowest(4.667, 78.833, -107.043): 8.509 ms (pygeodesy 24.8.24 Python 3.12.5 64bit arm64 macOS 14.6.1)
#
# _PGM('../testGeoids/egm96-5.pgm'): AREA_OR_POINT='Point', DateTime='2009-08-29 18:45:03', Description='WGS84 EGM96, 5-minute grid', Geoid='file in PGM format for the GeographicLib::Geoid class', MaxBilinearError=0.14, MaxCubicError=0.003, Offset=-108.0, Origin=LatLon2Tuple(lat=90.0, lon=0.0), Pixel=65535, RMSBilinearError=0.005, RMSCubicError=0.001, Scale=0.003, URL='http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html', Vertical_Datum='WGS84'
#
# Timbuktu GeoidKarney('egm96-5.pgm').height(16.775833, -3.009444): 28.7068 vs 28.7067
# Timbuktu GeoidKarney('egm96-5.pgm').height(16.776, -3.009): 28.7067 vs 28.7067


# % python3.8 -m pygeodesy.geoids -egm96 ../testGeoids/WW15MGH.GRD
#
# GeoidEGM96('WW15MGH.GRD'): lowerleft(-90.0, -180.0, -29.534), upperright(90.0, 180.25, 13.606), center(0.0, 0.125, 17.125), highest(-8.25, -32.75, 85.391), lowest(4.75, -101.25, -106.991): 1.267 ms (pygeodesy 24.12.24 Python 3.8.10 64bit arm64_x86_64 macOS 10.16)
#
# Timbuktu GeoidEGM96('WW15MGH.GRD').height(16.775833, -3.009444): 28.7073 vs 28.7880
# Timbuktu GeoidEGM96('WW15MGH.GRD').height(16.776, -3.009): 28.7072 vs 28.7880


# % python3.8 -m pygeodesy.geoids -Karney ../testGeoids/egm*.pgm
#
# GeoidKarney('egm2008-1.pgm'): lowerleft(-90.0, -180.0, -30.15), upperright(90.0, 180.0, 14.898), center(0.0, 0.0, 17.226), highest(-8.4, 147.367, 85.839), lowest(4.7, 78.767, -106.911): 353.050 ms (pygeodesy 24.8.24 Python 3.8.10 64bit arm64_x86_64 macOS 10.16)
#
# _PGM('../testGeoids/egm2008-1.pgm'): AREA_OR_POINT='Point', DateTime='2009-08-31 06:54:00', Description='WGS84 EGM2008, 1-minute grid', Geoid='file in PGM format for the GeographicLib::Geoid class', MaxBilinearError=0.025, MaxCubicError=0.003, Offset=-108.0, Origin=LatLon2Tuple(lat=90.0, lon=0.0), Pixel=65535, RMSBilinearError=0.001, RMSCubicError=0.001, Scale=0.003, URL='http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008', Vertical_Datum='WGS84'
#
# Timbuktu GeoidKarney('egm2008-1.pgm').height(16.775833, -3.009444): 28.7881 vs 28.7880
# Timbuktu GeoidKarney('egm2008-1.pgm').height(16.776, -3.009): 28.7880 vs 28.7880
#
# GeoidKarney('egm84-15.pgm'): lowerleft(-90.0, -180.0, -29.712), upperright(90.0, 180.0, 13.098), center(0.0, 0.0, 18.33), highest(-4.5, 148.75, 81.33), lowest(4.75, 79.25, -107.34): 1.727 ms (pygeodesy 24.8.24 Python 3.8.10 64bit arm64_x86_64 macOS 10.16)
#
# _PGM('../testGeoids/egm84-15.pgm'): AREA_OR_POINT='Point', DateTime='2009-08-29 18:45:02', Description='WGS84 EGM84, 15-minute grid', Geoid='file in PGM format for the GeographicLib::Geoid class', MaxBilinearError=0.413, MaxCubicError=0.02, Offset=-108.0, Origin=LatLon2Tuple(lat=90.0, lon=0.0), Pixel=65535, RMSBilinearError=0.018, RMSCubicError=0.001, Scale=0.003, URL='http://earth-info.nga.mil/GandG/wgs84/gravitymod/wgs84_180/wgs84_180.html', Vertical_Datum='WGS84'
#
# Timbuktu GeoidKarney('egm84-15.pgm').height(16.775833, -3.009444): 31.2983 vs 31.2979
# Timbuktu GeoidKarney('egm84-15.pgm').height(16.776, -3.009): 31.2979 vs 31.2979
#
# GeoidKarney('egm96-5.pgm'): lowerleft(-90.0, -180.0, -29.535), upperright(90.0, 180.0, 13.605), center(0.0, 0.0, 17.163), highest(-8.167, 147.25, 85.422), lowest(4.667, 78.833, -107.043): 14.807 ms (pygeodesy 24.8.24 Python 3.8.10 64bit arm64_x86_64 macOS 10.16)
#
# _PGM('../testGeoids/egm96-5.pgm'): AREA_OR_POINT='Point', DateTime='2009-08-29 18:45:03', Description='WGS84 EGM96, 5-minute grid', Geoid='file in PGM format for the GeographicLib::Geoid class', MaxBilinearError=0.14, MaxCubicError=0.003, Offset=-108.0, Origin=LatLon2Tuple(lat=90.0, lon=0.0), Pixel=65535, RMSBilinearError=0.005, RMSCubicError=0.001, Scale=0.003, URL='http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html', Vertical_Datum='WGS84'
#
# Timbuktu GeoidKarney('egm96-5.pgm').height(16.775833, -3.009444): 28.7068 vs 28.7067
# Timbuktu GeoidKarney('egm96-5.pgm').height(16.776, -3.009): 28.7067 vs 28.7067


# % python2 -m pygeodesy.geoids -Karney ../testGeoids/egm*.pgm
#
# GeoidKarney('egm2008-1.pgm'): lowerleft(-90.0, -180.0, -30.15), upperright(90.0, 180.0, 14.898), center(0.0, 0.0, 17.226), highest(-8.4, 147.367, 85.839), lowest(4.7, 78.767, -106.911): 283.362 ms (pygeodesy 24.8.24 Python 2.7.18 64bit arm64_x86_64 macOS 10.16)
#
# _PGM('../testGeoids/egm2008-1.pgm'): AREA_OR_POINT='Point', DateTime='2009-08-31 06:54:00', Description='WGS84 EGM2008, 1-minute grid', Geoid='file in PGM format for the GeographicLib::Geoid class', MaxBilinearError=0.025, MaxCubicError=0.003, Offset=-108.0, Origin=LatLon2Tuple(lat=90.0, lon=0.0), Pixel=65535, RMSBilinearError=0.001, RMSCubicError=0.001, Scale=0.003, URL='http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008', Vertical_Datum='WGS84'
#
# Timbuktu GeoidKarney('egm2008-1.pgm').height(16.775833, -3.009444): 28.7881 vs 28.7880
# Timbuktu GeoidKarney('egm2008-1.pgm').height(16.776, -3.009): 28.7880 vs 28.7880
#
# GeoidKarney('egm84-15.pgm'): lowerleft(-90.0, -180.0, -29.712), upperright(90.0, 180.0, 13.098), center(0.0, 0.0, 18.33), highest(-4.5, 148.75, 81.33), lowest(4.75, 79.25, -107.34): 1.378 ms (pygeodesy 24.8.24 Python 2.7.18 64bit arm64_x86_64 macOS 10.16)
#
# _PGM('../testGeoids/egm84-15.pgm'): AREA_OR_POINT='Point', DateTime='2009-08-29 18:45:02', Description='WGS84 EGM84, 15-minute grid', Geoid='file in PGM format for the GeographicLib::Geoid class', MaxBilinearError=0.413, MaxCubicError=0.02, Offset=-108.0, Origin=LatLon2Tuple(lat=90.0, lon=0.0), Pixel=65535, RMSBilinearError=0.018, RMSCubicError=0.001, Scale=0.003, URL='http://earth-info.nga.mil/GandG/wgs84/gravitymod/wgs84_180/wgs84_180.html', Vertical_Datum='WGS84'
#
# Timbuktu GeoidKarney('egm84-15.pgm').height(16.775833, -3.009444): 31.2983 vs 31.2979
# Timbuktu GeoidKarney('egm84-15.pgm').height(16.776, -3.009): 31.2979 vs 31.2979
#
# GeoidKarney('egm96-5.pgm'): lowerleft(-90.0, -180.0, -29.535), upperright(90.0, 180.0, 13.605), center(0.0, 0.0, 17.163), highest(-8.167, 147.25, 85.422), lowest(4.667, 78.833, -107.043): 11.612 ms (pygeodesy 24.8.24 Python 2.7.18 64bit arm64_x86_64 macOS 10.16)
#
# _PGM('../testGeoids/egm96-5.pgm'): AREA_OR_POINT='Point', DateTime='2009-08-29 18:45:03', Description='WGS84 EGM96, 5-minute grid', Geoid='file in PGM format for the GeographicLib::Geoid class', MaxBilinearError=0.14, MaxCubicError=0.003, Offset=-108.0, Origin=LatLon2Tuple(lat=90.0, lon=0.0), Pixel=65535, RMSBilinearError=0.005, RMSCubicError=0.001, Scale=0.003, URL='http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html', Vertical_Datum='WGS84'
#
# Timbuktu GeoidKarney('egm96-5.pgm').height(16.775833, -3.009444): 28.7068 vs 28.7067
# Timbuktu GeoidKarney('egm96-5.pgm').height(16.776, -3.009): 28.7067 vs 28.7067


# % python3.12 -m pygeodesy.geoids -PGM ../testGeoids/egm*.pgm
#
# GeoidPGM('egm2008-1.pgm'): lowerleft(-90.0, -180.0, -30.15), upperright(90.0, 180.0, 14.898), center(0.0, 0.0, 17.226), highest(-8.4, -32.633, 85.839), lowest(4.683, -101.25, -106.911): 543.148 ms (pygeodesy 24.8.24 Python 3.12.5 64bit arm64 macOS 14.6.1)
#
# _PGM('../testGeoids/egm2008-1.pgm'): AREA_OR_POINT='Point', DateTime='2009-08-31 06:54:00', Description='WGS84 EGM2008, 1-minute grid', Geoid='file in PGM format for the GeographicLib::Geoid class', MaxBilinearError=0.025, MaxCubicError=0.003, Offset=-108.0, Origin=LatLon2Tuple(lat=90.0, lon=0.0), Pixel=65535, RMSBilinearError=0.001, RMSCubicError=0.001, Scale=0.003, URL='http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008', Vertical_Datum='WGS84'
#
# Timbuktu GeoidPGM('egm2008-1.pgm').height(16.775833, -3.009444): 28.7881 vs 28.7880
# Timbuktu GeoidPGM('egm2008-1.pgm').height(16.776, -3.009): 28.7880 vs 28.7880
#
# GeoidPGM('egm84-15.pgm'): lowerleft(-90.0, -180.0, -29.712), upperright(90.0, 180.0, 13.098), center(0.0, 0.0, 18.33), highest(-4.5, -31.25, 81.33), lowest(4.75, -100.75, -107.34): 1.762 ms (pygeodesy 24.8.24 Python 3.12.5 64bit arm64 macOS 14.6.1)
#
# _PGM('../testGeoids/egm84-15.pgm'): AREA_OR_POINT='Point', DateTime='2009-08-29 18:45:02', Description='WGS84 EGM84, 15-minute grid', Geoid='file in PGM format for the GeographicLib::Geoid class', MaxBilinearError=0.413, MaxCubicError=0.02, Offset=-108.0, Origin=LatLon2Tuple(lat=90.0, lon=0.0), Pixel=65535, RMSBilinearError=0.018, RMSCubicError=0.001, Scale=0.003, URL='http://earth-info.nga.mil/GandG/wgs84/gravitymod/wgs84_180/wgs84_180.html', Vertical_Datum='WGS84'
#
# Timbuktu GeoidPGM('egm84-15.pgm').height(16.775833, -3.009444): 31.2979 vs 31.2979
# Timbuktu GeoidPGM('egm84-15.pgm').height(16.776, -3.009): 31.2975 vs 31.2979
#
# GeoidPGM('egm96-5.pgm'): lowerleft(-90.0, -180.0, -29.535), upperright(90.0, 180.0, 13.605), center(0.0, -0.0, 17.179), highest(-8.167, -32.75, 85.422), lowest(4.667, -101.167, -107.043): 12.594 ms (pygeodesy 24.8.24 Python 3.12.5 64bit arm64 macOS 14.6.1)
#
# _PGM('../testGeoids/egm96-5.pgm'): AREA_OR_POINT='Point', DateTime='2009-08-29 18:45:03', Description='WGS84 EGM96, 5-minute grid', Geoid='file in PGM format for the GeographicLib::Geoid class', MaxBilinearError=0.14, MaxCubicError=0.003, Offset=-108.0, Origin=LatLon2Tuple(lat=90.0, lon=0.0), Pixel=65535, RMSBilinearError=0.005, RMSCubicError=0.001, Scale=0.003, URL='http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html', Vertical_Datum='WGS84'
#
# Timbuktu GeoidPGM('egm96-5.pgm').height(16.775833, -3.009444): 28.7065 vs 28.7067
# Timbuktu GeoidPGM('egm96-5.pgm').height(16.776, -3.009): 28.7064 vs 28.7067

# **) MIT License
#
# Copyright (C) 2016-2025 -- mrJean1 at Gmail -- All Rights Reserved.
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
