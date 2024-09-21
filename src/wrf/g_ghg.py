from __future__ import (absolute_import, division, print_function)
import numpy as np
from .extension import _xghg
from .metadecorators import copy_and_set_metadata
from .util import extract_vars

@copy_and_set_metadata(copy_varname='CH4_BCK', name='column_averaged_ghg',
                       description='column averaged ghg',
                       units='ppmv')
def get_xghg(wrfin, ghg="xch4", timeidx=0, method="cat", squeeze=True,
             cache=None, meta=True, _key=None):
    """Return column averaged greenhouse gas mixing ratio.

     This functions extracts the necessary variables from the NetCDF file
     object in order to perform the calculation.

     Args:

         wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
             iterable): WRF-ARW NetCDF
             data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile`
             or an iterable sequence of the aforementioned types.

         ghg (:obj:`str`, optional): The desired greenhouse gas. Must be
             one of 'xch4', 'xco2', or 'xco'. The default is 'xch4'. A
             non-supported choice raises a NotImplementedError. 

         timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The
             desired time index. This value can be a positive integer,
             negative integer, or
             :data:`wrf.ALL_TIMES` (an alias for None) to return
             all times in the file or sequence. The default is 0.

         method (:obj:`str`, optional): The aggregation method to use for
             sequences.  Must be either 'cat' or 'join'.
             'cat' combines the data along the Time dimension.
             'join' creates a new dimension for the file index.
             The default is 'cat'.

         squeeze (:obj:`bool`, optional): Set to False to prevent dimensions
             with a size of 1 from being automatically removed from the
             shape of the output. Default is True.

         cache (:obj:`dict`, optional): A dictionary of (varname, ndarray)
             that can be used to supply pre-extracted NetCDF variables to
             the computational routines.  It is primarily used for internal
             purposes, but can also be used to improve performance by
             eliminating the need to repeatedly extract the same variables
             used in multiple diagnostics calculations, particularly when
             using large sequences of files.
             Default is None.

         meta (:obj:`bool`, optional): Set to False to disable metadata and
             return :class:`numpy.ndarray` instead of
             :class:`xarray.DataArray`.  Default is True.

         _key (:obj:`int`, optional): A caching key. This is used for
             internal purposes only.  Default is None.

     Returns:
         :class:`xarray.DataArray` or :class:`numpy.ndarray`: Omega.
         If xarray is
         enabled and the *meta* parameter is True, then the result will be a
         :class:`xarray.DataArray` object.  Otherwise, the result will be a
         :class:`numpy.ndarray` object with no metadata.
    """
     
    match ghg:
        case 'xch4' | 'ch4' | 'methane':
            varnames = ('P', 'PSFC', 'CH4_ANT', 'CH4_BIO','CH4_BCK')
            ncvars = extract_vars(wrfin, timeidx, varnames, method, squeeze, cache,
                                  meta=False, _key=_key)
            pres = ncvars['P']
            sfc_p = ncvars['PSFC']
            ant = ncvars['CH4_ANT']
            bio = ncvars['CH4_BIO']
            bck = ncvars['CH4_BCK']
        case 'xco2' | 'co2' | 'carbon dioxide':
            varnames = ('P', 'PSFC', 'CO2_ANT', 'CO2_BIO','CO2_BCK')
            ncvars = extract_vars(wrfin, timeidx, varnames, method, squeeze, cache,
                                  meta=False, _key=_key)
            pres = ncvars['P']
            sfc_p = ncvars['PSFC']
            ant = ncvars['CO2_ANT']
            bio = ncvars['CO2_BIO']
            bck = ncvars['CO2_BCK']
        case 'xco' | 'co' | 'carbon monoxide':
            varnames = ('P', 'PSFC', 'CO_ANT', 'CO_BCK')
            ncvars = extract_vars(wrfin, timeidx, varnames, method, squeeze, cache,
                                  meta=False, _key=_key)
            pres = ncvars['P']
            sfc_p = ncvars['PSFC']
            ant = ncvars['CO_ANT']
            bio = np.zeros_like(ant)
            bck = ncvars['CO_BCK']
            
        case _:
            raise NotImplementedError('Chemistry variable not supported at this time.')
    
    xghg = _xghg(sfc_p, pres, ant, bio, bck)

    return xghg