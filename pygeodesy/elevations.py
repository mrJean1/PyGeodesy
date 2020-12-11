
# -*- coding: utf-8 -*-

u'''Functions to obtain elevations and geoid heights thru web services,
for (lat, lon) locations, currently limited to the U{Conterminous US
(CONUS)<https://WikiPedia.org/wiki/Contiguous_United_States>}, see also
module L{geoids}, module L{heights} classes and U{USGS10mElev.py
<https://Gist.GitHub.com/pyRobShrk>}.

B{macOS}: If an C{SSLCertVerificationError} occurs, especially this
I{"[SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: self
"signed certificate in certificate chain ..."}, review U{this post
<https://StackOverflow.com/questions/27835619/urllib-and-ssl-certificate
-verify-failed-error>} for a remedy.  From a Terminal window run:
C{"/Applications/Python X.Y/Install Certificates.command"}

@newfield example: Example, Examples
'''

from pygeodesy.basics import clips, ub2str
from pygeodesy.errors import ParseError, _xkwds_get
from pygeodesy.interns import NN, _COLONSPACE_, _height_, \
                             _n_a_, _no_, _SPACE_
from pygeodesy.lazily import _ALL_LAZY
from pygeodesy.named import _NamedTuple
from pygeodesy.streprs import Fmt, fstr
from pygeodesy.units import Lat, Lon, Meter, Scalar, Str

__all__ = _ALL_LAZY.elevations
__version__ = '20.11.06'

try:
    from urllib2 import urlopen  # quote, urlcleanup
    from httplib import HTTPException as HTTPError

except (ImportError, NameError):  # Python 3+
    from urllib.request import urlopen  # urlcleanup
    # from urllib.parse import quote
    from urllib.error import HTTPError

_JSON_ = 'JSON'
_XML_  = 'XML'

try:
    from json import loads as _json
except ImportError:  # PYCHOK no cover

    from pygeodesy.interns import _COMMA_, _QUOTE2_

    def _json(ngs):
        '''(INTERNAL) Convert an NGS response in JSON to a C{dict}.
        '''
        # b'{"geoidModel": "GEOID12A",
        #    "station": "UserStation",
        #    "lat": 37.8816,
        #    "latDms": "N375253.76000",
        #    "lon": -121.9142,
        #    "lonDms": "W1215451.12000",
        #    "geoidHeight": -31.703,
        #    "error": 0.064
        #   }'
        #
        # or in case of errors:
        #
        # b'{"error": "No suitable Geoid model found for model 15"
        #   }'
        d = {}
        for t in ngs.strip().lstrip('{').rstrip('}').split(_COMMA_):
            t = t.strip()
            j = t.strip(_QUOTE2_).split('": ')
            if len(j) != 2:
                raise ParseError(json=t)
            k, v = j
            try:
                v = float(v)
            except (TypeError, ValueError):
                v = Str(ub2str(v.lstrip().lstrip(_QUOTE2_)), name=k)
            d[k] = v
        return d


def _error(fun, lat, lon, e):
    '''(INTERNAL) Format an error
    '''
    return _COLONSPACE_(Fmt.PAREN(fun.__name__, fstr((lat, lon))), e)


def _qURL(url, timeout=2, **params):
    '''(INTERNAL) Build B{C{url}} query, get and verify response.
    '''
    if params:  # build url query, don't map(quote, params)!
        p = '&'.join(Fmt.EQUAL(p, v) for p, v in params.items() if v)
        if p:
            url = NN(url, '?', p)
    u = urlopen(url, timeout=timeout)  # secs

    s = u.getcode()
    if s != 200:  # http.HTTPStatus.OK or http.client.OK
        raise IOError('code %d: %s' % (s, u.geturl()))

    r = u.read()
    u.close()
    # urlcleanup()
    return ub2str(r).strip()


def _xml(tag, xml):
    '''(INTERNAL) Get a <tag>value</tag> from XML.
    '''
    # b'<?xml version="1.0" encoding="utf-8" ?>
    #   <USGS_Elevation_Point_Query_Service>
    #    <Elevation_Query x="-121.914200" y="37.881600">
    #     <Data_Source>3DEP 1/3 arc-second</Data_Source>
    #     <Elevation>3851.03</Elevation>
    #     <Units>Feet</Units>
    #    </Elevation_Query>
    #   </USGS_Elevation_Point_Query_Service>'
    i = xml.find(Fmt.TAG(tag))
    if i > 0:
        i += len(tag) + 2
        j = xml.find(Fmt.TAGEND(tag), i)
        if j > i:
            return Str(xml[i:j].strip(), name=tag)
    return _no_(_XML_, Fmt.TAG(tag))


class Elevation2Tuple(_NamedTuple):  # .elevations.py
    '''2-Tuple C{(elevation, data_source)} in C{meter} and C{str}.
    '''
    _Names_ = ('elevation', 'data_source')
    _Units_ = ( Meter,       Str)


def elevation2(lat, lon, timeout=2.0):
    '''Get the geoid elevation at an C{NAD83} to C{NAVD88} location.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).
       @kwarg timeout: Optional, query timeout (seconds).

       @return: An L{Elevation2Tuple}C{(elevation, data_source)}
                or (C{None, "error"}) in case of errors.

       @raise ValueError: Invalid B{C{timeout}}.

       @note: The returned C{elevation} is C{None} if B{C{lat}} or B{C{lon}}
              is invalid or outside the C{Conterminous US (CONUS)},
              if conversion failed or if the query timed out.  The
              C{error} is the C{HTTP-, IO-, SSL-, Type-, URL-} or
              C{ValueError} as a string (C{str}).

       @see: U{USGS National Map<https://NationalMap.gov/epqs>},
             the U{FAQ<https://www.USGS.gov/faqs/what-are-projection-
             horizontal-and-vertical-datum-units-and-resolution-3dep-standard-dems>},
             U{geoid.py<https://Gist.GitHub.com/pyRobShrk>}, module
             L{geoids}, classes L{GeoidG2012B}, L{GeoidKarney} and
             L{GeoidPGM}.
    '''
    try:
        x = _qURL('https://NED.USGS.gov/epqs/pqs.php',  # 'https://NationalMap.gov/epqs/pqs.php'
                         x=Lon(lon).toStr(prec=6),
                         y=Lat(lat).toStr(prec=6),
                         units='Meters',  # 'Feet', capitalized
                         output=_XML_.lower(),  # _JSON_, lowercase only
                         timeout=Scalar(timeout=timeout))
        if x[:6] == '<?xml ':
            e = _xml('Elevation', x)
            try:
                e = float(e)
                if -1000000 < e < 1000000:
                    return Elevation2Tuple(e, _xml('Data_Source', x))
                e = 'non-CONUS %.2F' % (e,)
            except (TypeError, ValueError):
                pass
        else:
            e = _no_(_XML_, Fmt.QUOTE2(clips(x, limit=128, white=_SPACE_)))
    except (HTTPError, IOError, TypeError, ValueError) as x:
        e = repr(x)
    e = _error(elevation2, lat, lon, e)
    return Elevation2Tuple(None, e)


class GeoidHeight2Tuple(_NamedTuple):  # .elevations.py
    '''2-Tuple C{(height, model_name)}, geoid C{height} in C{meter}
       and C{model_name} as C{str}.
    '''
    _Names_ = (_height_, 'model_name')
    _Units_ = ( Meter,    Str)


def geoidHeight2(lat, lon, model=0, timeout=2.0):
    '''Get the C{NAVD88} geoid height at an C{NAD83} location.

       @arg lat: Latitude (C{degrees}).
       @arg lon: Longitude (C{degrees}).
       @kwarg model: Optional, geoid model ID (C{int}).
       @kwarg timeout: Optional, query timeout (seconds).

       @return: An L{GeoidHeight2Tuple}C{(height, model_name)}
                or C{(None, "error"}) in case of errors.

       @raise ValueError: Invalid B{C{timeout}}.

       @note: The returned C{height} is C{None} if B{C{lat}} or B{C{lon}} is
              invalid or outside the C{Conterminous US (CONUS)}, if the
              B{C{model}} was invalid, if conversion failed or if the
              query timed out.  The C{error} is the C{HTTP-, IO-, SSL-,
              Type-, URL-} or C{ValueError} as a string (C{str}).

       @see: U{NOAA National Geodetic Survey
             <https://www.NGS.NOAA.gov/INFO/geodesy.shtml>},
             U{Geoid<https://www.NGS.NOAA.gov/web_services/geoid.shtml>},
             U{USGS10mElev.py<https://Gist.GitHub.com/pyRobShrk>}, module
             L{geoids}, classes L{GeoidG2012B}, L{GeoidKarney} and
             L{GeoidPGM}.
    '''
    try:
        j = _qURL('https://Geodesy.NOAA.gov/api/geoid/ght',
                         lat=Lat(lat).toStr(prec=6),
                         lon=Lon(lon).toStr(prec=6),
                         model=(model if model else NN),
                         timeout=Scalar(timeout=timeout))  # PYCHOK indent
        if j[:1] == '{' and j[-1:] == '}' and j.find('"error":') > 0:
            d, e = _json(j), 'geoidHeight'
            if isinstance(_xkwds_get(d, error=_n_a_), float):
                h = d.get(e, None)
                if h is not None:
                    m = _xkwds_get(d, geoidModel=_n_a_)
                    return GeoidHeight2Tuple(h, m)
        else:
            e = _JSON_
        e = _no_(e, Fmt.QUOTE2(clips(j, limit=256, white=_SPACE_)))
    except (HTTPError, IOError, ParseError, TypeError, ValueError) as x:
        e = repr(x)
    e = _error(geoidHeight2, lat, lon, e)
    return GeoidHeight2Tuple(None, e)


if __name__ == '__main__':

    # <https://WikiPedia.org/wiki/Mount_Diablo>
    for f in (elevation2,     # (1173.79, '3DEP 1/3 arc-second')
              geoidHeight2):  # (-31.699, u'GEOID12B')
        t = f(37.8816, -121.9142)
        print(_COLONSPACE_(f.__name__, t))

# **) MIT License
#
# Copyright (C) 2016-2021 -- mrJean1 at Gmail -- All Rights Reserved.
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

# % python -m pygeodesy.elevations
# elevation2: (1173.79, '3DEP 1/3 arc-second')
# geoidHeight2: (-31.703, 'GEOID12B')
