
# -*- coding: utf-8 -*-

u'''Functions to obtain elevations and geoid heights thru web services,
for (lat, lon) locations, currently limited to the U{Conterminous US (CONUS)
<http://WikiPedia.org/wiki/Contiguous_United_States>}.

B{macOS}: If an C{SSLCertVerificationError} occurs, especially this
I{"[SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: self
"signed certificate in certificate chain ..."}, review U{this post
<http://StackOverflow.com/questions/27835619/urllib-and-ssl-certificate
-verify-failed-error>} for a remedy.  From a Terminal window run:
C{"/Applications/Python X.Y/Install Certificates.command"}

@newfield example: Example, Examples
'''

__all__ = 'elevation2', 'geoidHeight2'
__version__ = '18.10.04'

from fmath import fStr

try:
    from urllib2 import urlopen  # quote, urlcleanup
    _Bytes = unicode, bytearray  # PYCHOK expected

except (ImportError, NameError):  # Python 3+
    from urllib.request import urlopen  # urlcleanup
    # from urllib.parse import quote
    _Bytes = bytes, bytearray

try:
    from json import loads as _json
except ImportError:

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
        for t in ngs.strip().lstrip('{').rstrip('}').split(','):
            j = t.strip().strip('"').split('": ')
            if len(j) != 2:
                raise ValueError('json: %r' % (t,))
            k, v = j
            try:
                v = float(v)
            except ValueError:
                v = v.lstrip().lstrip('"')
            d[k] = v
        return d


def _error(fun, lat, lon, e):
    '''(INTERNAL) Format an error'
    '''
    return '%s(%s): %r' % (fun.__name__, fStr((lat, lon)), e)


def _qURL(url, *params, **timeout):
    '''(INTERNAL) Build I{url} query and get response.
    '''
    if params:  # build url query, do not map(quote, params)!
        url += '?' + '&'.join(_ for _ in map(str, params) if _)
    u = urlopen(url, **timeout)  # secs

    s = u.getcode()
    if s != 200:
        raise IOError('code %d: %s' % (s, u.geturl()))

    r = u.read()
    u.close()
    # urlcleanup()

    if isinstance(r, _Bytes):
        r = r.decode('utf-8')
    return r


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
    i = xml.find('<%s>' % (tag,))
    if i > 0:
        i += len(tag) + 2
        j = xml.find('</%s>' % (tag,), i)
        if j > i:
            return xml[i:j].strip()
    return 'no <%s>' % (tag,)


def elevation2(lat, lon, timeout=2.0):
    '''Get the geoid elevation at an C{NAD83} to C{NAVD88} location.

       @param lat: Latitude (C{degrees}).
       @param lon: Longitude (C{degrees}).
       @keyword timeout: Optional, query timeout (seconds).

       @return: 2-Tuple (elevation, data_source) in (C{meter}, C{str})
                or (C{None}, I{<error>}).

       @note: The returned elevation is C{None} if I{lat} or I{lon}
              is invalid or outside the C{Conterminous US (CONUS)},
              if conversion failed or if the query timed out.  The
              I{<error>} is the C{HTTP-, IO-, SSL-, Type-, URL-} or
              C{ValueError} as C{str}.

       @see: U{USGS National Map<http://NationalMap.gov/epqs>},
             the U{FAQ<http://www.USGS.gov/faqs/what-are-projection-
             horizontal-and-vertical-datum-units-and-resolution-3dep-standard-dems>}
             and U{USGS10mElev.py<http://gist.GitHub.com/pyRobShrk?page=2>}.
    '''
    try:
        x = _qURL('http://NED.USGS.gov/epqs/',
                         'x=%.6f' % (lon,),
                         'y=%.6f' % (lat,),
                         'units=Meters',
                         'output=xml',
                          timeout=float(timeout))

        e = _xml('Elevation', x)
        try:
            e = float(e)
            if -100000 < e < 1000000:
                return e, _xml('Data_Source', x)
            e = 'non-CONUS'
        except ValueError:
            pass

        raise ValueError(e)
    except (IOError, TypeError, ValueError) as x:
        return None, _error(elevation2, lat, lon, x)


def geoidHeight2(lat, lon, model=0, timeout=2.0):
    '''Get the C{NAVD88} geoid height at an C{NAD83} location.

       @param lat: Latitude (C{degrees}).
       @param lon: Longitude (C{degrees}).
       @keyword model: Optional, geoid model ID (C{int}).
       @keyword timeout: Optional, query timeout (seconds).

       @return: 2-Tuple (height, model_name) in (C{meter}, C{str}) or
                (C{None}, I{<error>}).

       @note: The returned height is C{None} if I{lat} or I{lon} is
              invalid or outside the C{Conterminous US (CONUS)}, if
              the I{model} was invalid, if conversion failed or if
              the query timed out.  The I{<error>} is the C{HTTP-,
              IO-, SSL-, Type-, URL-} or C{ValueError} as C{str}.

       @see: U{NOAA National Geodetic Survery
             <http://www.NGS.NOAA.gov/INFO/geodesy.shtml>} and
             U{Geoid<http://www.NGS.NOAA.gov/web_services/geoid.shtml>}
    '''
    try:
        d = _json(_qURL('http://Geodesy.NOAA.gov/api/geoid/ght',
                               'lat=%.6f' % (lat,),
                               'lon=%.6f' % (lon,),
                               'model=%d' % (model,) if model else '',
                                timeout=float(timeout)))

        e = d.get('error', 'N/A')
        if isinstance(e, float):
            h = d.get('geoidHeight', None)
            if h is not None:
                return h, d.get('geoidModel', 'N/A')

        raise ValueError(e)
    except (IOError, TypeError, ValueError) as x:
        return None, _error(geoidHeight2, lat, lon, x)


if __name__ == '__main__':

    # <http://WikiPedia.org/wiki/Mount_Diablo>
    for f in (elevation2,     # (1173.79, '3DEP 1/3 arc-second')
              geoidHeight2):  # (-31.703, 'GEOID12B')
        t = f(37.8816, -121.9142)
        print('%s: %s' % (f.__name__, t))
