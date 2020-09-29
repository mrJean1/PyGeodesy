# -*- coding: utf-8 -*-

u'''Single-instance floats and strings, C{intern}'ed across modules.
'''
from math import pi as PI

__all__ = ('EPS', 'EPS_2', 'EPS1', 'EPS1_2',
           'INF', 'MANTIS', 'MAX', 'MIN', 'NAN', 'NEG0',
           'NN',
           'PI', 'PI2', 'PI4', 'PI_2', 'PI_4',
           'R_M')  # import by .lazily
__version__ = '20.09.27'

NN = ''  # no name, empty str, Nomen Nescio <https://Wiktionary.org/wiki/N.N.>

# __DUNDER__ would get mangled in classes
_0_                  = '0'                    # PYCHOK expected
_1_                  = '1'                    # PYCHOK expected
_2_                  = '2'                    # PYCHOK expected
_3_                  = '3'                    # PYCHOK expected
_areaOf_             = 'areaOf'               # PYCHOK expected
_ambiguous_          = 'ambiguous'            # PYCHOK expected
_AT_                 = '@'                    # PYCHOK expected
_azimuth_            = 'azimuth'              # PYCHOK expected
_band_               = 'band'                 # PYCHOK expected
_bearing_            = 'bearing'              # PYCHOK expected
_C_                  = 'C'                    # PYCHOK expected
_Cartesian_          = 'Cartesian'            # PYCHOK expected
_coincident_         = 'coincident'           # PYCHOK expected
_colinear_           = 'colinear'             # PYCHOK expected
_COLON_              = ':'                    # PYCHOK expected
_COMMA_              = ','                    # PYCHOK expected
_convergence_        = 'convergence'          # PYCHOK expected
_cubic_              = 'cubic'                # PYCHOK expected
_CURLY_              = '{%s}'                 # PYCHOK expected
_datum_              = 'datum'                # PYCHOK expected
_decode3_            = 'decode3'              # PYCHOK expected
_deg_                = 'deg'                  # PYCHOK expected
_degrees_            = 'degrees'              # PYCHOK expected
_degrees2_           = 'degrees2'             # PYCHOK SQUARED
_distance_           = 'distance'             # PYCHOK expected
_distanceTo_         = 'distanceTo'           # PYCHOK expected
_doesn_t_exist_      = "doesn't exist"        # PYCHOK expected
_DOT_                = '.'                    # PYCHOK expected
_E_                  = 'E'                    # PYCHOK expected
_easting_            = 'easting'              # PYCHOK expected
_ellipsoid_          = 'ellipsoid'            # PYCHOK expected
_ellipsoidal_        = 'ellipsoidal'          # PYCHOK expected
_encode_             = 'encode'               # PYCHOK expected
_end_                = 'end'                  # PYCHOK expected
_epoch_              = 'epoch'                # PYCHOK expected
_EQUAL_              = '='                    # PYCHOK expected
_exceed_PI_radians_  = 'exceed PI radians'    # PYCHOK expected
_feet_               = 'feet'                 # PYCHOK expected
_fraction_           = 'fraction'             # PYCHOK expected
_gamma_              = 'gamma'                # PYCHOK expected
_h_                  = 'h'                    # PYCHOK expected
_height_             = 'height'               # PYCHOK expected
_hemipole_           = 'hemipole'             # PYCHOK expected
_intersection_       = 'intersection'         # PYCHOK expected
_inside_             = 'inside'               # PYCHOK expected
_invalid_            = 'invalid'              # PYCHOK expected
_isclockwise_        = 'isclockwise'          # PYCHOK expected
_ispolar_            = 'ispolar'              # PYCHOK expected
_k0_                 = 'k0'                   # PYCHOK expected
_knots_              = 'knots'                # PYCHOK expected
_lam_                = 'lam'                  # PYCHOK expected
_lat_                = 'lat'                  # PYCHOK expected
_lat0_               = 'lat0'                 # PYCHOK expected
_lat1_               = 'lat1'                 # PYCHOK expected
_lat2_               = 'lat2'                 # PYCHOK expected
_LatLon_             = 'LatLon'               # PYCHOK expected
_len_                = 'len'                  # PYCHOK expected
_linear_             = 'linear'               # PYCHOK expected
_lon_                = 'lon'                  # PYCHOK expected
_lon0_               = 'lon0'                 # PYCHOK expected
_m_                  = 'm'                    # PYCHOK expected
_M_                  = 'M'                    # PYCHOK expected
_meanOf_             = 'meanOf'               # PYCHOK expected
_meter_              = 'meter'                # PYCHOK expected
_MGRS_               = 'MGRS'                 # PYCHOK expected
_n_                  = 'n'                    # PYCHOK expected
_N_                  = 'N'                    # PYCHOK expected
_n_a_                = 'n/a'                  # PYCHOK expected
_name_               = 'name'                 # PYCHOK expected
_near_concentric_    = 'near-concentric'      # PYCHOK expected
_nearestOn2_         = 'nearestOn2'           # PYCHOK expected
_negative_           = 'negative'             # PYCHOK expected
_no_intersection_    = 'no intersection'      # PYCHOK expected
_no_convergence_     = 'no convergence'       # PYCHOK expected
_no_convergence_fmt_ = 'no convergence (%g)'  # PYCHOK expected
_no_conversion_      = 'no conversion'        # PYCHOK expected
_no_overlap_         = 'no overlap'           # PYCHOK expected
_northing_           = 'northing'             # PYCHOK expected
_NorthPole_          = 'NorthPole'            # PYCHOK expected
_not_convex_         = 'not convex'           # PYCHOK expected
_not_scalar_         = 'not scalar'           # PYCHOK expected
_number_             = 'number'               # PYCHOK expected
_Nvector_            = 'Nvector'              # PYCHOK expected
_other_              = 'other'                # PYCHOK expected
_outside_            = 'outside'              # PYCHOK expected
_PARENTH_            = '(%s)'                 # PYCHOK expected
_PERCENT_            = '%'                    # PYCHOK expected
_perimeterOf_        = 'perimeterOf'          # PYCHOK expected
_phi_                = 'phi'                  # PYCHOK expected
_PLUS_               = '+'                    # PYCHOK expected
_point_              = 'point'                # PYCHOK expected
_points_             = 'points'               # PYCHOK expected
_pole_               = 'pole'                 # PYCHOK expected
_precision_          = 'precision'            # PYCHOK expected
_radians_            = 'radians'              # PYCHOK expected
_radians2_           = 'radians2'             # PYCHOK SQUARED
_radius_             = 'radius'               # PYCHOK expected
_radius1_            = 'radius1'              # PYCHOK expected
_radius2_            = 'radius2'              # PYCHOK expected
_range_              = 'range'                # PYCHOK expected
_reciprocal_         = 'reciprocal'           # PYCHOK expected
_resolution_         = 'resolution'           # PYCHOK expected
_S_                  = 'S'                    # PYCHOK expected
_scalar_             = 'scalar'               # PYCHOK expected
_scale_              = 'scale'                # PYCHOK expected
_spherical_          = 'spherical'            # PYCHOK expected
_SouthPole_          = 'SouthPole'            # PYCHOK expected
_SPACE_              = ' '                    # PYCHOK expected
_SQUARE_             = '[%s]'                 # PYCHOK expected
_STAR_               = '*'                    # PYCHOK expected
_start_              = 'start'                # PYCHOK expected
_std_                = 'std'                  # PYCHOK expected
_sumOf_              = 'sumOf'                # PYCHOK expected
_too_distant_        = 'too distant'          # PYCHOK expected
_too_distant_fmt_    = 'too distant (%.3g)'   # PYCHOK expected
_too_few_            = 'too few'              # PYCHOK expected
_too_small_          = 'too small'            # PYCHOK expected
_transform_          = 'transform'            # PYCHOK expected
_UNDERSCORE_         = '_'                    # PYCHOK expected
_units_              = 'units'                # PYCHOK expected
_UPS_                = 'UPS'                  # PYCHOK expected
_utf_8_              = 'utf-8'                # PYCHOK expected
_UTM_                = 'UTM'                  # PYCHOK expected
_valid_              = 'valid'                # PYCHOK expected
_W_                  = 'W'                    # PYCHOK expected
_x_                  = 'x'                    # PYCHOK expected
_y_                  = 'y'                    # PYCHOK expected
_z_                  = 'z'                    # PYCHOK expected
_zone_               = 'zone'                 # PYCHOK expected

_COLON_SPACE_        = _COLON_ + _SPACE_      # PYCHOK expected
_COMMA_SPACE_        = _COMMA_ + _SPACE_      # PYCHOK expected
_DUNDER_             = _UNDERSCORE_ * 2       # PYCHOK expected
_EW_                 = _E_  + _W_             # PYCHOK common cardinals
_NE_                 = _N_  + _E_             # PYCHOK expected
_NS_                 = _N_  + _S_             # PYCHOK expected
_NSEW_               = _NS_ + _EW_            # PYCHOK expected
_NW_                 = _N_  + _W_             # PYCHOK expected
_SE_                 = _S_  + _E_             # PYCHOK expected
_SW_                 = _S_  + _W_             # PYCHOK negative ones


class _Missing(object):
    '''(INTERNAL) Singleton.
    '''
    def toRepr(self, **unused):
        return 'MISSING'  # self.__class__.__name__.strip(_UNDERSCORE_),lower()

    __repr__ = toRepr
    __str__  = toRepr
    toStr    = toRepr

_Missing = _Missing()  # PYCHOK singleton


def _dot_(*prefix_names):
    '''(INTERNAL) Period-join C{prefix} and C{name}.
    '''
    return _DOT_.join(prefix_names)


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

_0_0    = _float(   0)    # PYCHOK expected
_0_1    = _float(   0.1)  # PYCHOK expected
_0_5    = _float(   0.5)  # PYCHOK expected
_1_0    = _float(   1)    # PYCHOK expected
_2_0    = _float(   2)    # PYCHOK expected
_3_0    = _float(   3)    # PYCHOK expected
_4_0    = _float(   4)    # PYCHOK expected
_5_0    = _float(   5)    # PYCHOK expected
_6_0    = _float(   6)    # PYCHOK expected
_8_0    = _float(   8)    # PYCHOK expected
_16_0   = _float(  16)    # PYCHOK expected
_32_0   = _float(  32)    # PYCHOK expected
_60_0   = _float(  60)    # PYCHOK expected
_90_0   = _float(  90)    # PYCHOK expected
_180_0  = _float( 180)    # PYCHOK expected
_360_0  = _float( 360)    # PYCHOK expected
_3600_0 = _float(3600)    # PYCHOK expected

try:
    from sys import float_info as _float_info

    EPS    = _float(_float_info.epsilon)   # PYCHOK system's EPSilon
    MANTIS =        _float_info.mant_dig   # PYCHOK system's mantissa bits
    MAX    = _float(_float_info.max)       # PYCHOK system's MAX float
    MIN    = _float(_float_info.min)       # PYCHOK system's MIN float
except (AttributeError, ImportError):  # PYCHOK no cover
    EPS    = _float(2.220446049250313e-16)  # PYCHOK EPSilon 2**-52?, M{EPS +/- 1 != 1}
    MANTIS =  53  # PYCHOK mantissa bits ≈53 (C{int})
    MAX    = _float(pow(_2_0,  1023) * (_2_0 - EPS))  # PYCHOK ≈ 10**308, 2**1024?
    MIN    = _float(pow(_2_0, -1021))  # PYCHOK ≈ 10**-308, 2**-1021?

EPS2   = _float(EPS * _2_0)    # PYCHOK ≈ 4.440892098501e-16
EPS_2  = _float(EPS / _2_0)    # PYCHOK ≈ 1.110223024625e-16
EPS1   = _float(_1_0 - EPS)    # PYCHOK ≈ 0.9999999999999998
EPS1_2 = _float(_1_0 - EPS_2)  # PYCHOK ≈ 0.9999999999999999
# _1EPS  = _float(_1_0 + EPS)  # PYCHOK ≈ 1.0000000000000002

if not _0_0 < EPS < EPS1 < _1_0:  # for .frechet
    raise AssertionError('%s < %s: %s < %s < %.16f < %s' % ('EPS', 'EPS1', _0_0, EPS, EPS1, _1_0))

INF    = _float( 'INF')  # PYCHOK INFinity, see function L{isinf}, L{isfinite}
NAN    = _float( 'NAN')  # PYCHOK Not-A-Number, see function L{isnan}
NEG0   =  float('-0.0')  # PYCHOK NEGative 0.0, see function L{isneg0}

PI2    = _float(PI * _2_0)  # PYCHOK Two PI, M{PI * 2} aka I{Tau}
PI4    = _float(PI * _4_0)  # PYCHOK Four PI, M{PI * 4}
PI_2   = _float(PI / _2_0)  # PYCHOK Half PI, M{PI / 2}
PI_4   = _float(PI / _4_0)  # PYCHOK Quarter PI, M{PI / 4}

R_M    = _float(6371008.771415)  # PYCHOK mean, spherical earth radius


def _item_fmt(fmt, name_value_arg, name_value_kwd):
    '''(INTERNAL) Helper for C{_item_pr}, C{_item_ps} and C{_item_sq}.
    '''
    for n_v in name_value_kwd.items():
        break
    else:
        if len(name_value_arg) > 1:
            n_v = name_value_arg[:2]
        elif name_value_arg:
            n_v = name_value_arg, _Missing
        else:
            n_v = _Missing, _Missing
    return fmt % n_v


def _item_cs(*name_value_arg, **name_value_kwd):  # PYCHOK expected
    '''(INTERNAL) Return a named value string.
    '''
    return _item_fmt('%s: %s', name_value_arg, name_value_kwd)


def _item_pr(*name_value_arg, **name_value_kwd):  # PYCHOK expected
    '''(INTERNAL) Return a parenthesized name representation.
    '''
    return _item_fmt('%s(%r)', name_value_arg, name_value_kwd)


def _item_ps(*name_value_arg, **name_value_kwd):  # PYCHOK expected
    '''(INTERNAL) Return a parenthesized name string.
    '''
    return _item_fmt('%s(%s)', name_value_arg, name_value_kwd)


def _item_sq(*name_value_arg, **name_value_kwd):  # PYCHOK expected
    '''(INTERNAL) Return an indexed or keyed name.
    '''
    return _item_fmt('%s[%r]', name_value_arg, name_value_kwd)

# **) MIT License
#
# Copyright (C) 2016-2020 -- mrJean1 at Gmail -- All Rights Reserved.
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
