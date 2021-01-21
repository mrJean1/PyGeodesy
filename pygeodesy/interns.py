# -*- coding: utf-8 -*-

u'''Single-instance floats and strings, C{intern}'ed across modules.
'''
from math import pi as PI, sqrt

__all__ = ('EPS', 'EPS_2', 'EPS0', 'EPS1', 'EPS1_2',
           'INF', 'MANTIS', 'MAX', 'MIN',  # not 'MISSING'!
           'NAN', 'NEG0', 'NN',
           'PI', 'PI2', 'PI4', 'PI_2', 'PI_4', 'PI3_2')  # imported by .lazily
__version__ = '21.01.19'


class _Join(str):
    '''(INTERNAL) Extended, callable C{str}.
    '''
    def join_(self, *args):
        '''Join all B{C{args}} like C{str.join(B{args})}.
        '''
        return _Join(str.join(self, map(str, args)))  # re-callable

    __call__ = join_


class _Prefix(_Join):
    '''(INTERNAL) Extended C{str} for prefix.
    '''
    def __call__(self, *args):
        '''Join C{self} plus all B{C{args}} like C{str.join((self,) + B{args})}.
        '''
        return _SPACE_.join_(self, *args)  # re-callable


class _Python_(str):  # overwritten below
    '''(INTERNAL) Extended C{str} for C{Python} and version.
    '''
    def __call__(self, sys):
        '''Return C{"Python <version>"}.
        '''
        return _SPACE_(self, sys.version.split()[0])


class _Slicer(str):
    '''(INTERNAL) String slicer C{.fromX} or C{.tillY}.
    '''
    def __getattr__(self, name):  # .fromX, .tillY
        if name.startswith(_till_):
            i = self.find(name[len(_till_):])
            if 0 < i < len(self):
                return _Slicer(self[:i + 1])
        elif name.startswith(_from_):
            i = self.find(name[len(_from_):])
            if 0 < (i + 1) < len(self):
                return _Slicer(self[i:])
        else:
            return getattr(str, name)
        return self


class MISSING(object):
    '''(INTERNAL) Singleton.
    '''
    def toRepr(self, **unused):
        return self.__class__.__name__

    __repr__ = toRepr
    __str__  = toRepr
    toStr    = toRepr

MISSING          = MISSING()  # PYCHOK singleton
MISSING.__name__ = str(MISSING)

NN = _Join('')  # Nomen Nescio <https://Wiktionary.org/wiki/N.N.>

# __DUNDER__-style names would get mangled in classes
_0_                   = '0'                  # PYCHOK expected
_0to9_                = '0123456789'         # PYCHOK expected
_1_                   = '1'                  # PYCHOK expected
_2_                   = '2'                  # PYCHOK expected
_3_                   = '3'                  # PYCHOK expected
_4_                   = '4'                  # PYCHOK expected
_a_                   = 'a'                  # PYCHOK expected
_A_                   = 'A'                  # PYCHOK expected
_Airy1830_            = 'Airy1830'           # PYCHOK expected
_AiryModified_        = 'AiryModified'       # PYCHOK expected
_an_                  = 'an'                 # PYCHOK expected
_angle_               = 'angle'              # PYCHOK expected
_areaOf_              = 'areaOf'             # PYCHOK expected
_ambiguous_           = 'ambiguous'          # PYCHOK expected
# _AMPERSAND_   = _Join('&')                 # PYCHOK expected
# _AND_               = _AMPERSAND_          # PYCHOK expected
_AT_            = _Join('@')                 # PYCHOK expected
_AtoZnoIO_    = _Slicer('ABCDEFGHJKLMNPQRSTUVWXYZ')  # PYCHOK in C{gars}, C{mgrs} and C{wgrs}
_attribute_           = 'attribute'          # PYCHOK expected
_azimuth_             = 'azimuth'            # PYCHOK expected
_band_                = 'band'               # PYCHOK expected
_BAR_           = _Join('|')                 # PYCHOK expected
_bearing_             = 'bearing'            # PYCHOK expected
_Bessel1841_          = 'Bessel1841'         # PYCHOK expected
_by_                  = 'by'                 # PYCHOK expected
_C_                   = 'C'                  # PYCHOK expected
_Cartesian_           = 'Cartesian'          # PYCHOK expected
_Clarke1866_          = 'Clarke1866'         # PYCHOK expected
_Clarke1880IGN_       = 'Clarke1880IGN'      # PYCHOK expected
_coincident_          = 'coincident'         # PYCHOK expected
_colinear_            = 'colinear'           # PYCHOK expected
_COLON_         = _Join(':')                 # PYCHOK expected
_COLONSPACE_    = _Join(': ')                # PYCHOK expected
_COMMA_         = _Join(',')                 # PYCHOK expected
_COMMASPACE_    = _Join(', ')                # PYCHOK expected
_convergence_ = _Prefix('convergence')       # PYCHOK expected
_conversion_          = 'conversion'         # PYCHOK expected
_convex_              = 'convex'             # PYCHOK expected
_cubic_               = 'cubic'              # PYCHOK expected
_DASH_          = _Join('-')                 # PYCHOK == _MINUS_
_datum_               = 'datum'              # PYCHOK expected
_decode3_             = 'decode3'            # PYCHOK expected
_deg_                 = 'deg'                # PYCHOK expected
_degrees_             = 'degrees'            # PYCHOK expected
_degrees2_            = 'degrees2'           # PYCHOK SQUARED
_DEQUALSPACED_  = _Join(' == ')              # PYCHOK expected
_distance_            = 'distance'           # PYCHOK expected
_distanceTo_          = 'distanceTo'         # PYCHOK expected
_distant_     = _Prefix('distant')           # PYCHOK expected
_doesn_t_exist_       = "doesn't exist"      # PYCHOK expected
_DOT_           = _Join('.')                 # PYCHOK expected
_e_                   = 'e'                  # PYCHOK expected
_E_                   = 'E'                  # PYCHOK expected
_easting_             = 'easting'            # PYCHOK expected
_ELLIPSIS_      = _Join('...')               # PYCHOK expected
# _ELLIPSISPACED_ = _Join(' ... ')           # PYCHOK <https://www.ThePunctuationGuide.com/ellipses.html>
_ellipsoid_           = 'ellipsoid'          # PYCHOK expected
_ellipsoidal_         = 'ellipsoidal'        # PYCHOK expected
_enabled_             = 'enabled'            # PYCHOK expected
_encode_              = 'encode'             # PYCHOK expected
_end_                 = 'end'                # PYCHOK expected
_epoch_               = 'epoch'              # PYCHOK expected
_EQUAL_         = _Join('=')                 # PYCHOK expected
_EQUALSPACED_   = _Join(' = ')               # PYCHOK expected
_exceed_PI_radians_   = 'exceed PI radians'  # PYCHOK expected
_exceeds_     = _Prefix('exceeds')           # PYCHOK expected
_f_                   = 'f'                  # PYCHOK expected
_feet_                = 'feet'               # PYCHOK expected
_few_                 = 'few'                # PYCHOK expected
_finite_              = 'finite'             # PYCHOK expected
_fraction_            = 'fraction'           # PYCHOK expected
_from_                = 'from'               # PYCHOK expected
_g_                   = 'g'                  # PYCHOK expected
_gamma_               = 'gamma'              # PYCHOK expected
_GRS80_               = 'GRS80'              # PYCHOK expected
_h_                   = 'h'                  # PYCHOK expected
_H_                   = 'H'                  # PYCHOK expected
_height_              = 'height'             # PYCHOK expected
_hemipole_            = 'hemipole'           # PYCHOK expected
_immutable_           = 'immutable'          # PYCHOK expected
_i_                   = 'i'                  # PYCHOK expected
_in_                  = 'in'                 # PYCHOK expected
_INF_                 = 'INF'                # PYCHOK expected
_initial_             = 'initial'            # PYCHOK expected
_inside_              = 'inside'             # PYCHOK expected
_intersection_        = 'intersection'       # PYCHOK expected
_Intl1924_            = 'Intl1924'           # PYCHOK expected
_invalid_             = 'invalid'            # PYCHOK expected
_isclockwise_         = 'isclockwise'        # PYCHOK expected
_ispolar_             = 'ispolar'            # PYCHOK expected
_j_                   = 'j'                  # PYCHOK expected
_k0_                  = 'k0'                 # PYCHOK expected
_knots_               = 'knots'              # PYCHOK expected
_Krassovski1940_      = 'Krassovski1940'     # PYCHOK expected
_Krassowsky1940_      = 'Krassowsky1940'     # PYCHOK expected
#_LANGLE_             = '<'                  # PYCHOK expected
_lam_                 = 'lam'                # PYCHOK expected
_lat_                 = 'lat'                # PYCHOK expected
_lat0_                = 'lat0'               # PYCHOK expected
_lat1_                = 'lat1'               # PYCHOK expected
_lat2_                = 'lat2'               # PYCHOK expected
_latlon_              = 'latlon'             # PYCHOK expected
_LatLon_              = 'LatLon'             # PYCHOK expected
#_LCURLY_             = '{'                  # PYCHOK LBRACE
_len_                 = 'len'                # PYCHOK expected
_linear_              = 'linear'             # PYCHOK expected
#_LPAREN_             = '('                  # PYCHOK expected
_lon_                 = 'lon'                # PYCHOK expected
_lon0_                = 'lon0'               # PYCHOK expected
#_LSQUARE_            = '['                  # PYCHOK LBRACK
_m_                   = 'm'                  # PYCHOK expected
_M_                   = 'M'                  # PYCHOK expected
_meanOf_              = 'meanOf'             # PYCHOK expected
_meridional_          = 'meridional'         # PYCHOK expected
_meter_               = 'meter'              # PYCHOK expected
_meter2_              = 'meter2'             # PYCHOK SQUARED
_MGRS_                = 'MGRS'               # PYCHOK expected
_MINUS_               = _DASH_
_module_              = 'module'             # PYCHOK expected
_n_                   = 'n'                  # PYCHOK expected
_N_                   = 'N'                  # PYCHOK expected
_n_a_                 = 'n/a'                # PYCHOK expected
_N_A_                 = 'N/A'                # PYCHOK expected
_NAD27_               = 'NAD27'              # PYCHOK expected
_NAD83_               = 'NAD83'              # PYCHOK expected
_name_                = 'name'               # PYCHOK expected
_NAN_                 = 'NAN'                # PYCHOK expected
_near_concentric_     = 'near-concentric'    # PYCHOK expected
_nearestOn2_          = 'nearestOn2'         # PYCHOK expected
_negative_            = 'negative'           # PYCHOK expected
_NL_            = _Join('\n')                # PYCHOK expected
_NL_hash_       = _Join(_NL_ + '# ')         # PYCHOK expected
_NL_var_        = _Join(_NL_ + '@var ')      # PYCHOK expected
_no_          = _Prefix('no')                # PYCHOK expected
_northing_            = 'northing'           # PYCHOK expected
_NorthPole_           = 'NorthPole'          # PYCHOK expected
_not_         = _Prefix('not')               # PYCHOK expected
_NTF_                 = 'NTF'                # PYCHOK expected
_number_              = 'number'             # PYCHOK expected
_numpy_               = 'numpy'              # PYCHOK expected
_OKd_                 = '._-'                # PYCHOK expected
_on_                  = 'on'                 # PYCHOK expected
_or_                  = 'or'                 # PYCHOK expected
_other_               = 'other'              # PYCHOK expected
_outside_             = 'outside'            # PYCHOK expected
_overlap_             = 'overlap'            # PYCHOK expected
_PERCENT_             = '%'                  # PYCHOK expected
_PERCENTDOTSTAR_      = '%.*'                # PYCHOK _DOT_(_PERCENT_, _STAR_)
_perimeterOf_         = 'perimeterOf'        # PYCHOK expected
_phi_                 = 'phi'                # PYCHOK expected
_PLUS_          = _Join('+')                 # PYCHOK expected
_PLUSMINUS_           = _PLUS_ + _MINUS_     # PYCHOK expected
_point_               = 'point'              # PYCHOK expected
_points_              = 'points'             # PYCHOK expected
_pole_                = 'pole'               # PYCHOK expected
_precision_           = 'precision'          # PYCHOK expected
_prime_vertical_      = 'prime_vertical'     # PYCHOK expected
_pygeodesy_abspath_   = 'pygeodesy_abspath'  # PYCHOK expected
_Python_     = _Python_('Python')            # PYCHOK singleton
# _QUOTE1_            = "'"                  # PYCHOK expected
_QUOTE2_              = '"'                  # PYCHOK expected
# _QUOTE3_            = "'''"                # PYCHOK expected
# _QUOTE6_            = '"""'                # PYCHOK expected
_radians_             = 'radians'            # PYCHOK expected
_radians2_            = 'radians2'           # PYCHOK SQUARED
_radius_              = 'radius'             # PYCHOK expected
_radius1_             = 'radius1'            # PYCHOK expected
_radius2_             = 'radius2'            # PYCHOK expected
# _range_        = _Range('range')           # moved down
#_RANGLE_             = '>'                  # PYCHOK expected
#_RCURLY_             = '}'                  # PYCHOK RBRACE
_reciprocal_          = 'reciprocal'         # PYCHOK expected
_resolution_          = 'resolution'         # PYCHOK expected
#_RPAREN_             = ')'                  # PYCHOK expected
#_RSQUARE_            = ']'                  # PYCHOK RBRACK
_s_                   = 's'                  # PYCHOK expected
_S_                   = 'S'                  # PYCHOK expected
_scalar_              = 'scalar'             # PYCHOK expected
_scale_               = 'scale'              # PYCHOK expected
_scipy_               = 'scipy'              # PYCHOK expected
_sep_                 = 'sep'                # PYCHOK expected
_singular_            = 'singular'           # PYCHOK expected
_small_               = 'small'              # PYCHOK expected
_Sphere_              = 'Sphere'             # PYCHOK expected
_spherical_           = 'spherical'          # PYCHOK expected
_SouthPole_           = 'SouthPole'          # PYCHOK expected
_SPACE_         = _Join(' ')                 # PYCHOK expected
_STAR_          = _Join('*')                 # PYCHOK expected
_start_               = 'start'              # PYCHOK expected
_std_                 = 'std'                # PYCHOK expected
_supported_           = 'supported'          # PYCHOK expected
_tbd_                 = 'tbd'                # PYCHOK expected
_till_                = 'till'               # PYCHOK expected
_to_                  = 'to'                 # PYCHOK expected
_too_         = _Prefix('too')               # PYCHOK expected
_transform_           = 'transform'          # PYCHOK expected
_UNDER_               = '_'                  # PYCHOK expected
_units_               = 'units'              # PYCHOK expected
_UPS_                 = 'UPS'                # PYCHOK expected
_utf_8_               = 'utf-8'              # PYCHOK expected
_UTM_                 = 'UTM'                # PYCHOK expected
_V_                   = 'V'                  # PYCHOK expected
_valid_               = 'valid'              # PYCHOK expected
_version_             = 'version'            # PYCHOK expected
_vs_                  = 'vs'                 # PYCHOK expected
_W_                   = 'W'                  # PYCHOK expected
_WGS72_               = 'WGS72'              # PYCHOK expected
_WGS84_               = 'WGS84'              # PYCHOK expected
_x_                   = 'x'                  # PYCHOK expected
_X_                   = 'X'                  # PYCHOK expected
_y_                   = 'y'                  # PYCHOK expected
_z_                   = 'z'                  # PYCHOK expected
_zone_                = 'zone'               # PYCHOK expected

_EW_                  = _E_  + _W_           # PYCHOK common cardinals
_NE_                  = _N_  + _E_           # PYCHOK expected
_NS_                  = _N_  + _S_           # PYCHOK expected
_NSEW_                = _NS_ + _EW_          # PYCHOK expected
_NW_                  = _N_  + _W_           # PYCHOK expected
_SE_                  = _S_  + _E_           # PYCHOK expected
_SW_                  = _S_  + _W_           # PYCHOK negative ones

_DDOT_          = _Join(_DOT_ * 2)           # PYCHOK expected
# _DEQUAL_      = _Join(_EQUAL_ * 2)         # PYCHOK expected
_DUNDER_        = _Join(_UNDER_ * 2)         # PYCHOK expected


class _Range(str):
    '''(INTERNAL) Extended C{str} for C{range} strings.
    '''
    def __call__(self, lo, hi, prec=0, lopen=False, ropen=False,
                               join=_COMMASPACE_):
        '''Return the range as C{"(lo, hi)"}, C{"(lo, hi]"},
           C{"[lo, hi)"} or C{"[lo, hi]"}.
        '''
        from pygeodesy.streprs import Fmt  # PYCHOK re-imported
        r = NN(Fmt.f(lo, prec=prec), join,
               Fmt.f(hi, prec=prec))
        if lopen:
            r = Fmt.PAREN(r) if ropen else Fmt.LOPEN(r)
        else:
            r = Fmt.ROPEN(r) if ropen else Fmt.SQUARE(r)
        return r

_range_ = _Range('range')  # PYCHOK expected


def _dunder_name(inst, *dflt):
    '''(INTERNAL) Get the double_underscore __name__ attr.
    '''
    try:
        return inst.__name__
    except AttributeError:
        pass
    return dflt[0] if dflt else inst.__class__.__name__


def _float(f):  # in .datums, .ellipsoids, .trf
    '''(INTERNAL) cache initial C{float}s.
    '''
    f = float(f)
    return _floats.setdefault(f, f)  # PYCHOK del _floats


def _floatuple(*fs):
    '''(INTERNAL) Cache a tuple of C{float}s.
    '''
    return tuple(map(_float, fs))


_floats = {}      # PYCHOK floats cache, in .__main__
# _float = float  # PYCHOK expected
# del _floats     # XXX zap floats cache never

_0_0    = _float(   0)      # PYCHOK expected
_0_001  = _float(   0.001)  # PYCHOK expected
_0_01   = _float(   0.01)   # PYCHOK expected
_0_1    = _float(   0.1)    # PYCHOK expected
_0_125  = _float(   0.125)  # PYCHOK expected
_0_25   = _float(   0.25)   # PYCHOK expected
_0_26   = _float(   0.26)   # PYCHOK expected
_0_5    = _float(   0.5)    # PYCHOK expected
_1_0    = _float(   1)      # PYCHOK expected
_1_5    = _float(   1.5)    # PYCHOK expected
_2_0    = _float(   2)      # PYCHOK expected
_3_0    = _float(   3)      # PYCHOK expected
_4_0    = _float(   4)      # PYCHOK expected
_5_0    = _float(   5)      # PYCHOK expected
_6_0    = _float(   6)      # PYCHOK expected
_8_0    = _float(   8)      # PYCHOK expected
_9_0    = _float(   9)      # PYCHOK expected
_16_0   = _float(  16)      # PYCHOK expected
_24_0   = _float(  24)      # PYCHOK expected
_32_0   = _float(  32)      # PYCHOK expected
_60_0   = _float(  60)      # PYCHOK expected
_90_0   = _float(  90)      # PYCHOK expected
_120_0  = _float( 120)      # PYCHOK expected
_180_0  = _float( 180)      # PYCHOK expected
_360_0  = _float( 360)      # PYCHOK expected
_720_0  = _float( 720)      # PYCHOK expected
_3600_0 = _float(3600)      # PYCHOK expected

try:
    from sys import float_info as _float_info

    EPS    = _float(_float_info.epsilon)   # PYCHOK system's EPSilon
    MANTIS =        _float_info.mant_dig   # PYCHOK system's mantissa bits
    MAX    = _float(_float_info.max)       # PYCHOK system's MAX float
    MIN    = _float(_float_info.min)       # PYCHOK system's MIN float
except (AttributeError, ImportError):  # PYCHOK no cover
    EPS    = _float(2.220446049250313e-16)  # PYCHOK EPSilon 2**-52, M{EPS +/- 1 != 1}
    MANTIS =  53  # PYCHOK mantissa bits ≈53 (C{int})
    MAX    = _float(pow(_2_0,  1023) * (_2_0 - EPS))  # PYCHOK ≈ 10**308
    MIN    = _float(pow(_2_0, -1022))  # PYCHOK ≈ 10**-308

EPS2     = _float(EPS * _2_0)      # PYCHOK ≈ 4.440892098501e-16
EPS_2    = _float(EPS / _2_0)      # PYCHOK ≈ 1.110223024625e-16
EPS1     = _float(_1_0 - EPS)      # PYCHOK ≈ 0.9999999999999998
EPS1_2   = _float(_1_0 - EPS_2)    # PYCHOK ≈ 0.9999999999999999
# _1EPS  = _float(_1_0 + EPS)      # PYCHOK ≈ 1.0000000000000002
_1_EPS   = _float(_1_0 / EPS)      # PYCHOK = 4503599627370496.0
# _2_EPS = _float(_2_0 / EPS)      # PYCHOK = 9007199254740992.0
_EPS4    = _float(EPS * _4_0)      # PYCHOK ≈ 8.881784197001e-16
_EPSqrt  = _float(sqrt(EPS))       # PYCHOK = 1.49011611938e5-08
_EPStol  = _float(_EPSqrt * _0_1)  # PYCHOK = 1.49011611938e5-09
_EPS__2  = _float(EPS**2)          # PYCHOK = 4.930380657631e-32
_EPS__4  = _float(EPS**4)          # PYCHOK = 2.430865342915e-63
_EPS0__2 = _EPS__4                 # PYCHOK near-zero square root

EPS0     = _EPS__2          # PYCHOK near-zero comparison, or EPS or EPS_2
INF      = _float( _INF_)   # PYCHOK INFinity, see function L{isinf}, L{isfinite}
NAN      = _float( _NAN_)   # PYCHOK Not-A-Number, see function L{isnan}
NEG0     =  float('-0.0')   # PYCHOK NEGative 0.0, see function L{isneg0}

PI2   = _float(PI * _2_0)  # PYCHOK Two PI, M{PI * 2} aka I{Tau}
PI4   = _float(PI * _4_0)  # PYCHOK Four PI, M{PI * 4}
PI_2  = _float(PI / _2_0)  # PYCHOK Half PI, M{PI / 2}
PI_4  = _float(PI / _4_0)  # PYCHOK Quarter PI, M{PI / 4}
PI3_2 = _float(PI * _1_5)  # PYCHOK PI and a half, M{PI * 3 / 2}

R_M   = _float(6371008.771415)  # PYCHOK mean, spherical earth radius (C{meter})

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
