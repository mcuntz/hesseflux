#!/usr/bin/env python
"""
nee2gpp : Estimates photosynthesis (GPP) and ecosystem respiration (RECO)
          from Eddy covariance CO2 flux data.

This module was written by Matthias Cuntz while at Department of
Computational Hydrosystems, Helmholtz Centre for Environmental
Research - UFZ, Leipzig, Germany, and continued while at Institut
National de Recherche pour l'Agriculture, l'Alimentation et
l'Environnement (INRAE), Nancy, France.

Copyright (c) 2012-2020 Matthias Cuntz - mc (at) macu (dot) de
Released under the MIT License; see LICENSE file for details.

* Written 2012 by Matthias Cuntz - mc (at) macu (dot) de
* Set default undef to NaN, Mar 2012, Arndt Piayda
* Add wrapper nee2gpp for individual routines, Nov 2012, Matthias Cuntz
* Ported to Python 3, Feb 2013, Matthias Cuntz
* Use generel cost function cost_abs from functions module, May 2013, Matthias Cuntz
* Use fmin_tnc to allow params < 0, Aug 2014, Arndt Piayda
* Keyword nogppnight, Aug 2014, Arndt Piayda
* Add wrapper nee2gpp for individual routines, Nov 2012, Matthias Cuntz
* Add wrapper nee2gpp for individual routines, Nov 2012, Matthias Cuntz
* Input can be pandas Dataframe or numpy array(s), Apr 2020, Matthias Cuntz
* Using numpy docstring format, May 2020, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz, Arndt Piayda

The following functions are provided

.. autosummary::
   nee2gpp
"""
from __future__ import division, absolute_import, print_function
import warnings
import numpy as np
import pandas as pd
import scipy.optimize as opt # curve_fit, fmin, fmin_tnc
try:    # import package
    from .mad import mad
    from .functions import cost_abs, lloyd_only_rref_p
    from .functions import cost_lloyd_fix, lloyd_fix, lloyd_fix_p
    from .functions import cost_lasslop, lasslop
except: # python nee2gpp.py
    from mad import mad
    from functions import cost_abs, lloyd_only_rref_p
    from functions import cost_lloyd_fix, lloyd_fix, lloyd_fix_p
    from functions import cost_lasslop, lasslop


__all__ = ['nee2gpp']


# ----------------------------------------------------------------------
def nee2gpp(dfin, flag=None, isday=None, date=None, timeformat='%Y-%m-%d %H:%M:%S', colhead=None,
            undef=-9999, method='reichstein', nogppnight=False, swthr=10.):
    """
    Calculate photosynthesis (GPP) and ecosystem respiration (RECO)
    from Eddy covariance CO2 flux data.

    It uses either

      1. a fit of Reco vs. temperature to all nighttime data (`method='falge'`), or

      2. several fits over the season of Reco vs. temperature as in Reichstein
         et al. (2005) (`method='reichstein'`), or

      3. the daytime method of Lasslop et al. (2010) (`method='lasslop'`),

    in order to calculate Reco and then GPP = Reco - NEE.

    Parameters
    ----------
    dfin : pandas.Dataframe or numpy.array
        time series of CO2 fluxes and air temperature, and possibly
        incoming shortwave radiation and air vapour pressure deficit.

        `dfin` can be a pandas.Dataframe with the columns
        'FC' or 'NEE' (or starting with 'FC\_' or 'NEE\_') for observed CO2 flux [umol(CO2) m-2 s-1]
        'TA'    (or starting with 'TA\_') for air temperature [K]

        `method='lasslop'` or `method='day'` needs also
        'SW_IN' (or starting with 'SW_IN') for incoming short-wave radiation [W m-2]
        'VPD'   (or starting with 'VPD') for air vapour deficit [Pa]
        The index is taken as date variable.

        `dfin` can also me a numpy array with the same columns. In this case
        `colhead`, `date`, and possibly `dateformat` must be given.
    flag : pandas.Dataframe or numpy.array, optional
        flag Dataframe or array has the same shape as `dfin`. Non-zero values in
        `flag` will be treated as missing values in `dfin`.

        `flag` must follow the same rules as `dfin` if pandas.Dataframe.

        If `flag` is numpy array, `df.columns.values` will be used as column heads
        and the index of `dfin` will be copied to `flag`.
    isday : array_like of bool, optional
        True when it is day, False when night. Must have the same length as `dfin.shape[0]`.

        If `isday` is not given, `dfin` must have a column with head 'SW_IN' or
        starting with 'SW_IN'. `isday` will then be `dfin['SW_IN'] > swthr`.
    date : array_like of string, optional
        1D-array_like of calendar dates in format given in `timeformat`.

        `date` must be given if `dfin` is numpy array.
    timeformat : str, optional
        Format of dates in `date`, if given (default: '%Y-%m-%d %H:%M:%S').
        See strftime documentation of Python's datetime module:
        https://docs.python.org/3/library/datetime.html#strftime-and-strptime-behavior
    colhed : array_like of str, optional
        column names if `dfin` is numpy array. See `dfin` for mandatory column names.
    undef : float, optional
        values having `undef` value are treated as missing values in `dfin` (default: -9999)
    method : str, optional
        method to use for partitioning. Possible values are:

        'global' or 'falge':     fit of Reco vs. temperature to all nighttime data

        'local' of 'reichstein': several fits over the season of Reco vs. temperature
                                 as in Reichstein et al. (2005) (default)

        'day' or 'lasslop':      method of Lasslop et al. (2010) fitting a light-response curve
    nogppnight : float, optional
        GPP will be set to zero at night. RECO will then equal NEE at night (default: False)
    swthr : float, optional
        Threshold to determine daytime from incoming shortwave radiation if `isday` not given (default: 10).

    Returns
    -------
    pandas.Dataframe or numpy arrays
        pandas.Dataframe with two columns 'GPP' and 'RECO' with estimated
        photosynthesis and ecosystem respiration, or two numpy arrays
        [GPP, RECO].

    Notes
    -----
    Negative respiration possible at night if GPP is forced to 0 with `nogppnight=True`.

    Literature
    ----------
    Falge et al. (2001)
        Gap filling strategies for defensible annual sums of net ecosystem exchange,
        Acricultural and Forest Meteorology 107, 43-69

    Lasslop et al. (2010)
        Separation of net ecosystem exchange into assimilation and respiration using
        a light response curve approach: critical issues and global evaluation,
        Global Change Biology 16, 187-208

    Reichstein et al. (2005)
        On the separation of net ecosystem exchange into assimilation and ecosystem
        respiration: review and improved algorithm,
        Global Change Biology 11, 1424-1439

    Examples
    --------
    >>> from fread import fread
    >>> from date2dec import date2dec
    >>> from dec2date import dec2date
    >>> ifile = 'test_nee2gpp.csv'
    >>> undef = -9999.
    >>> dat   = fread(ifile, skip=2, transpose=True)
    >>> ndat  = dat.shape[1]
    >>> head  = fread(ifile, skip=2, header=True)
    >>> head1 = head[0]
    >>> # date
    >>> jdate = date2dec(dy=dat[0,:], mo=dat[1,:], yr=dat[2,:], hr=dat[3,:], mi=dat[4,:])
    >>> adate = dec2date(jdate, eng=True)
    >>> # colhead
    >>> idx   = []
    >>> for i in head1:
    ...     if i in ['NEE', 'rg', 'Tair', 'VPD']: idx.append(head1.index(i))
    >>> colhead = ['FC', 'SW_IN', 'TA', 'VPD']
    >>> # data
    >>> dfin = dat[idx,:]
    >>> dfin[2,:] = np.where(dfin[2,:] == undef, undef, dfin[2,:]+273.15)
    >>> dfin[3,:] = np.where(dfin[3,:] == undef, undef, dfin[3,:]*100.)
    >>> # flag
    >>> flag = np.where(dfin == undef, 2, 0)
    >>> # partition
    >>> GPP, Reco = nee2gpp(dfin, flag=flag, date=adate, colhead=colhead, undef=undef, method='local')
    >>> print(GPP[1120:1128])
    [-9.99900000e+03 -9.99900000e+03 -9.99900000e+03  4.40606871e+00
      8.31942152e+00  1.06242542e+01  8.49245664e+00  1.12381973e+01]
    >>> print(Reco[1120:1128])
    [1.68311981 1.81012431 1.9874173  2.17108871 2.38759152 2.64372415
     2.90076664 3.18592735]

    >>> GPP, Reco = nee2gpp(dfin, flag=flag, date=adate, colhead=colhead, undef=undef, method='local')
    >>> print(GPP[1120:1128])
    [-9.99900000e+03 -9.99900000e+03 -9.99900000e+03  4.40606871e+00
      8.31942152e+00  1.06242542e+01  8.49245664e+00  1.12381973e+01]

    >>> GPP, Reco = nee2gpp(dfin, flag=flag, date=adate, colhead=colhead, undef=undef, method='global')
    >>> print(GPP[1120:1128])
    [-9.99900000e+03 -9.99900000e+03 -9.99900000e+03  4.33166157e+00
      8.18228013e+00  1.04092252e+01  8.19395317e+00  1.08427448e+01]

    >>> GPP, Reco = nee2gpp(dfin, flag=flag, date=adate, colhead=colhead, undef=undef, method='Reichstein')
    >>> print(GPP[1120:1128])
    [-9.99900000e+03 -9.99900000e+03 -9.99900000e+03  4.40606871e+00
      8.31942152e+00  1.06242542e+01  8.49245664e+00  1.12381973e+01]

    >>> GPP, Reco = nee2gpp(dfin, flag=flag, date=adate, colhead=colhead, undef=undef, method='reichstein')
    >>> print(GPP[1120:1128])
    [-9.99900000e+03 -9.99900000e+03 -9.99900000e+03  4.40606871e+00
      8.31942152e+00  1.06242542e+01  8.49245664e+00  1.12381973e+01]

    >>> GPP, Reco = nee2gpp(dfin, flag=flag, date=adate, colhead=colhead, undef=undef, method='day')
    >>> print(GPP[1120:1128])
    [-9.99900000e+03 -9.99900000e+03 -9.99900000e+03  2.78457540e+00
      6.63212545e+00  8.88902165e+00  6.74243873e+00  9.51364527e+00]
    >>> print(Reco[1120:1128])
    [0.28786696 0.34594516 0.43893276 0.5495954  0.70029545 0.90849165
     1.15074873 1.46137527]

    History
    -------
    Written  Matthias Cuntz, Mar 2012
    Modified Arndt Piayda,   Mar 2012 - undef=np.nan
             Matthias Cuntz, Nov 2012 - wrapper for individual routines nee2gpp_reichstein etc.
             Matthias Cuntz, Feb 2013 - ported to Python 3
             Matthias Cuntz, May 2013 - replaced cost functions by generel cost function cost_abs if possible
             Arndt Piayda,   Aug 2014 - replaced fmin with fmin_tnc to permit params<0,
                                        permit gpp<0 at any time if nogppnight=True 
    """
    # Check input
    # numpy or panda
    if isinstance(dfin, (np.ndarray, np.ma.MaskedArray)):
        isnumpy = True
        istrans = False
        assert colhead is not None, 'colhead must be given if input is numpy.ndarray.'
        if dfin.shape[0] == len(colhead):
            istrans = True
            df = pd.DataFrame(dfin.T, columns=colhead)
        elif dfin.shape[1] == len(colhead):
            df = pd.DataFrame(dfin, columns=colhead)
        else:
            raise ValueError('Length of colhead must be number of columns in input array. len(colhead)='+str(len(colhead))+' shape(input)=('+str(dfin.shape[0])+','+str(dfin.shape[1])+').')
        assert date is not None, 'Date must be given if input is numpy arrary.'
        df['Datetime'] = pd.to_datetime(date, format=timeformat)
        df.set_index('Datetime', drop=True, inplace=True)
    else:
        isnumpy = False
        istrans = False
        assert isinstance(dfin, pd.core.frame.DataFrame), 'Input must be either numpy.ndarray or pandas.DataFrame.'
        df = dfin.copy(deep=True)

    # Incoming flags
    if flag is not None:
        if isinstance(flag, (np.ndarray, np.ma.MaskedArray)):
            fisnumpy = True
            fistrans = False
            if flag.shape[0] == len(df):
                ff = pd.DataFrame(flag, columns=df.columns.values)
            elif flag.shape[1] == len(df):
                fistrans = True
                ff = pd.DataFrame(flag.T, columns=df.columns.values)
            else:
                raise ValueError('flag must have same shape as data array. data: ({:d},{:d}); flag: ({:d},{:d})'.format(dfin.shape[0], dfin.shape[1], flag.shape[0], flag.shape[1]))
            ff = ff.set_index(df.index)
        else:
            fisnumpy = False
            fistrans = False
            assert isinstance(flag, pd.core.frame.DataFrame), 'Flag must be either numpy.ndarray or pandas.DataFrame.'
            ff = flag.copy(deep=True)
    else:
        fisnumpy = isnumpy
        fistrans = istrans
        # flags: 0: good; 1: input flagged; 2: output flagged
        ff              = df.copy(deep=True).astype(int)
        ff[:]           = 0
        ff[df == undef] = 1
        ff[df.isna()]   = 1

    # day or night
    if isday is None:
        sw_id = ''
        for cc in df.columns:
            if cc.startswith('SW_IN'):
                sw_id = cc
                break
        assert sw_id, 'Global radiation with name SW or starting with SW_ must be in input if isday not given.'
        isday = df[sw_id] > swthr # Papale et al. (Biogeosciences, 2006): 20; REddyProc: 10
    if isinstance(isday, (pd.core.series.Series,pd.core.frame.DataFrame)):
        isday = isday.to_numpy()
    isday[isday == undef] = np.nan
    ff[np.isnan(isday)]   = 1

    # Global relationship in Reichstein et al. (2005)
    if ((method.lower() == 'global') | (method.lower() == 'falge')):
        dfout = _nee2gpp_falge(df, ff, isday, undef=undef)
    # Local relationship = Reichstein et al. (2005)
    elif ((method.lower() == 'local') | (method.lower() == 'reichstein')):
        dfout = _nee2gpp_reichstein(df, ff, isday, undef=undef, nogppnight=nogppnight)
    # Lasslop et al. (2010) daytime method
    elif ((method.lower() == 'day') | (method.lower() == 'lasslop')):
        dfout = _nee2gpp_lasslop(df, ff, isday, undef=undef, nogppnight=nogppnight)
    # Include new methods here
    else:
        raise ValueError('Error nee2gpp: method not implemented yet.')

    if isnumpy:
        return dfout['GPP'].to_numpy(), dfout['RECO'].to_numpy()
    else:
        return dfout


# ----------------------------------------------------------------------
def _nee2gpp_falge(df, ff, isday, undef=-9999):
    """
    Calculate photosynthesis (GPP) and ecosystem respiration (RECO) from original
    Eddy flux data, using a fit of Reco vs. temperature to all nighttime data,
    in order to calculate Reco and then GPP = Reco - NEE.

    Parameters
    ----------
    df : pandas.Dataframe
        time series of CO2 fluxes and air temperature.

        pandas.Dataframe with the columns
        'FC' or 'NEE' (or starting with 'FC\_' or 'NEE\_') for observed CO2 flux [umol(CO2) m-2 s-1]
        'TA'    (or starting with 'TA\_') for air temperature [K]
        The index is taken as date variable.
    ff : pandas.Dataframe
        flag Dataframe or array has the same shape as `df`. Non-zero values in
        `ff` will be treated as missing values in `df`.

        `ff` must follow the same rules as `df`.
    isday : array_like of bool
        True when it is day, False when night. Must have the same length as `df.shape[0].`
    undef : float, optional
        values having `undef` value are treated as missing values in `df` (default: -9999)

    Returns
    -------
    pandas.Dataframe or numpy arrays
        pandas.Dataframe with two columns 'GPP' and 'RECO' with estimated
        photosynthesis and ecosystem respiration.

    Literature
    ----------
    Falge et al. (2001)
        Gap filling strategies for defensible annual sums of net ecosystem exchange,
        Acricultural and Forest Meteorology 107, 43-69

    Examples
    --------
    >>> from fread import fread
    >>> from date2dec import date2dec
    >>> from dec2date import dec2date
    >>> ifile = 'test_nee2gpp.csv'
    >>> undef = -9999.
    >>> dat   = fread(ifile, skip=2, transpose=True)
    >>> ndat  = dat.shape[1]
    >>> head  = fread(ifile, skip=2, header=True)
    >>> head1 = head[0]
    >>> # date
    >>> jdate = date2dec(dy=dat[0,:], mo=dat[1,:], yr=dat[2,:], hr=dat[3,:], mi=dat[4,:])
    >>> adate = dec2date(jdate, eng=True)
    >>> # colhead
    >>> idx   = []
    >>> for i in head1:
    ...     if i in ['NEE', 'rg', 'Tair', 'VPD']: idx.append(head1.index(i))
    >>> colhead = ['FC', 'SW_IN', 'TA', 'VPD']
    >>> # data
    >>> dfin = dat[idx,:]
    >>> dfin[2,:] = np.where(dfin[2,:] == undef, undef, dfin[2,:]+273.15)
    >>> dfin[3,:] = np.where(dfin[3,:] == undef, undef, dfin[3,:]*100.)
    >>> # flag
    >>> flag = np.where(dfin == undef, 2, 0)
    >>> # partition
    >>> GPP, Reco = nee2gpp(dfin, flag=flag, date=adate, colhead=colhead, undef=undef, method='global')
    >>> print(GPP[1120:1128])
    [-9.99900000e+03 -9.99900000e+03 -9.99900000e+03  4.33166157e+00
      8.18228013e+00  1.04092252e+01  8.19395317e+00  1.08427448e+01]

    History
    -------
    Written  Matthias Cuntz, Mar 2012
    Modified Arndt Piayda,   Mar 2012 - undef=np.nan
             Matthias Cuntz, Nov 2012 - individual routine
             Matthias Cuntz, Feb 2013 - ported to Python 3
    """
    # Variables
    fc_id = ''
    for cc in df.columns:
        if cc.startswith('FC_') or (cc == 'FC') or cc.startswith('NEE_') or (cc == 'NEE'):
            fc_id = cc
            break
    ta_id = ''
    for cc in df.columns:
        if cc.startswith('TA_') or (cc == 'TA'):
            ta_id = cc
            break
    assert fc_id, 'Carbon net flux with name FC or NEE or starting with FC_ or NEE_ must be in input.'
    assert ta_id, 'Air temperature with name TA or starting with TA_ must be in input.'

    nee    = np.ma.array(df[fc_id], mask=(ff[fc_id] > 0))
    t      = np.ma.array(df[ta_id], mask=(ff[ta_id] > 0))
    misday = np.ma.array(isday, mask=((~np.isfinite(isday)) | (isday == undef)))

    # Partition - Global relationship as in Falge et al. (2001)

    # Select valid nighttime
    ndata = nee.size
    mask = misday | nee.mask | t.mask | misday.mask
    ii   = np.where(~mask)[0]
    tt   = np.ma.compressed(t[ii])
    net  = np.ma.compressed(nee[ii])
    p        = opt.fmin(cost_abs, [2., 200.],
                        args=(lloyd_fix_p, tt, net), disp=False)

    Reco     = np.ones(ndata)*undef
    ii       = np.where(~t.mask)[0]
    Reco[ii] = lloyd_fix(t[ii], p[0], p[1])

    # GPP
    GPP     = np.ones(ndata)*undef
    ii      = np.where(~(t.mask | nee.mask))[0]
    GPP[ii] = Reco[ii] - nee[ii]

    dfout = pd.DataFrame({'GPP':GPP, 'RECO':Reco}, index=df.index)

    return dfout


# ----------------------------------------------------------------------
def _nee2gpp_reichstein(df, ff, isday, undef=-9999, nogppnight=False):
    """
    Calculate photosynthesis (GPP) and ecosystem respiration (RECO) from original
    Eddy flux data, using several fits of Reco vs. temperature of nighttime data
    over the season, as in Reichstein et al. (2005), in order to calculate Reco
    and then GPP = Reco - NEE.

    Parameters
    ----------
    df : pandas.Dataframe
        time series of CO2 fluxes and air temperature.

        pandas.Dataframe with the columns
        'FC' or 'NEE' (or starting with 'FC\_' or 'NEE\_') for observed CO2 flux [umol(CO2) m-2 s-1]
        'TA'    (or starting with 'TA\_') for air temperature [K]
        The index is taken as date variable.
    ff : pandas.Dataframe
        flag Dataframe or array has the same shape as `df`. Non-zero values in
        `ff` will be treated as missing values in `df`.

        `ff` must follow the same rules as `df`.
    isday : array_like of bool
        True when it is day, False when night. Must have the same length as `df.shape[0].`
    undef : float, optional
        values having `undef` value are treated as missing values in `df` (default: -9999)
    nogppnight : float, optional
        GPP will be set to zero at night. RECO will then equal NEE at night (default: False)

    Returns
    -------
    pandas.Dataframe
        pandas.Dataframe with two columns 'GPP' and 'RECO' with estimated
        photosynthesis and ecosystem respiration.

    Literature
    ----------
    Reichstein et al. (2005)
        On the separation of net ecosystem exchange into assimilation and ecosystem
        respiration: review and improved algorithm,
        Global Change Biology 11, 1424-1439

    Examples
    --------
    >>> from fread import fread
    >>> from date2dec import date2dec
    >>> from dec2date import dec2date
    >>> ifile = 'test_nee2gpp.csv'
    >>> undef = -9999.
    >>> dat   = fread(ifile, skip=2, transpose=True)
    >>> ndat  = dat.shape[1]
    >>> head  = fread(ifile, skip=2, header=True)
    >>> head1 = head[0]
    >>> # date
    >>> jdate = date2dec(dy=dat[0,:], mo=dat[1,:], yr=dat[2,:], hr=dat[3,:], mi=dat[4,:])
    >>> adate = dec2date(jdate, eng=True)
    >>> # colhead
    >>> idx   = []
    >>> for i in head1:
    ...     if i in ['NEE', 'rg', 'Tair', 'VPD']: idx.append(head1.index(i))
    >>> colhead = ['FC', 'SW_IN', 'TA', 'VPD']
    >>> # data
    >>> dfin = dat[idx,:]
    >>> dfin[2,:] = np.where(dfin[2,:] == undef, undef, dfin[2,:]+273.15)
    >>> dfin[3,:] = np.where(dfin[3,:] == undef, undef, dfin[3,:]*100.)
    >>> # flag
    >>> flag = np.where(dfin == undef, 2, 0)
    >>> # partition
    >>> GPP, Reco = nee2gpp(dfin, flag=flag, date=adate, colhead=colhead, undef=undef, method='local')
    >>> print(GPP[1120:1128])
    [-9.99900000e+03 -9.99900000e+03 -9.99900000e+03  4.40606871e+00
      8.31942152e+00  1.06242542e+01  8.49245664e+00  1.12381973e+01]
    >>> print(Reco[1120:1128])
    [1.68311981 1.81012431 1.9874173  2.17108871 2.38759152 2.64372415
     2.90076664 3.18592735]

    >>> GPP, Reco = nee2gpp(dfin, flag=flag, date=adate, colhead=colhead, undef=undef, method='local')
    >>> print(GPP[1120:1128])
    [-9.99900000e+03 -9.99900000e+03 -9.99900000e+03  4.40606871e+00
      8.31942152e+00  1.06242542e+01  8.49245664e+00  1.12381973e+01]

    >>> GPP, Reco = nee2gpp(dfin, flag=flag, date=adate, colhead=colhead, undef=undef, method='Reichstein')
    >>> print(GPP[1120:1128])
    [-9.99900000e+03 -9.99900000e+03 -9.99900000e+03  4.40606871e+00
      8.31942152e+00  1.06242542e+01  8.49245664e+00  1.12381973e+01]

    >>> GPP, Reco = nee2gpp(dfin, flag=flag, date=adate, colhead=colhead, undef=undef, method='reichstein')
    >>> print(GPP[1120:1128])
    [-9.99900000e+03 -9.99900000e+03 -9.99900000e+03  4.40606871e+00
      8.31942152e+00  1.06242542e+01  8.49245664e+00  1.12381973e+01]

    History
    -------
    Written  Matthias Cuntz, Mar 2012
    Modified Arndt Piayda,   Mar 2012 - undef=np.nan
             Matthias Cuntz, Nov 2012 - individual routine
             Matthias Cuntz, Feb 2013 - ported to Python 3
    """
    # Variables
    fc_id = ''
    for cc in df.columns:
        if cc.startswith('FC_') or (cc == 'FC') or cc.startswith('NEE_') or (cc == 'NEE'):
            fc_id = cc
            break
    ta_id = ''
    for cc in df.columns:
        if cc.startswith('TA_') or (cc == 'TA'):
            ta_id = cc
            break
    assert fc_id, 'Carbon net flux with name FC or NEE or starting with FC_ or NEE_ must be in input.'
    assert ta_id, 'Air temperature with name TA or starting with TA_ must be in input.'

    nee    = np.ma.array(df[fc_id], mask=(ff[fc_id] > 0))
    t      = np.ma.array(df[ta_id], mask=(ff[ta_id] > 0))
    misday = np.ma.array(isday, mask=((~np.isfinite(isday)) | (isday == undef)))
    dates  = df.index.to_julian_date()

    # Partition - Local relationship = Reichstein et al. (2005)

    ndata = nee.size
    GPP   = np.ones(ndata)*undef
    Reco  = np.ones(ndata)*undef
    dfout = pd.DataFrame({'GPP':GPP, 'RECO':Reco}, index=df.index)

    # Select valid nighttime
    mask  = misday | nee.mask | t.mask | misday.mask
    ii    = np.where(~mask)[0]
    if (ii.size==0):
        # raise ValueError('Error _nee2gpp_reichstein: no valid nighttime data.')
        print('Warning _nee2gpp_reichstein: no valid nighttime data.')
        return dfout
    jul  = dates[ii]
    tt   = np.ma.compressed(t[ii])
    net  = np.ma.compressed(nee[ii])
    # 1. each 5 days, in 15 day period, fit if range of T > 5
    locp = [] # local param
    locs = [] # local err
    dmin = np.floor(np.amin(jul)).astype(np.int) # be aware that julian days starts at noon, i.e. 1.0 is 12h
    dmax = np.ceil(np.amax(jul)).astype(np.int)  # so the search will be from noon to noon and thus includes all nights
    for i in range(dmin,dmax,5):
        iii  = np.where((jul>=i) & (jul<(i+14)))[0]
        niii = iii.size
        if niii > 6:
            tt1  = tt[iii]
            net1 = net[iii]
            mm   = ~mad(net1, z=4.5) # make fit more robust by removing outliers
            if (np.ptp(tt[iii]) >= 5.) & (np.sum(mm) > 6):
                p, temp1, temp2 = opt.fmin_tnc(cost_lloyd_fix, [2., 200.],
                                               bounds=[[0.,None], [0.,None]],
                                               args=(tt1[mm], net1[mm]),
                                               approx_grad=True, disp=False)
                try:
                    p1, c = opt.curve_fit(lloyd_fix, tt1[mm], net1[mm], p0=p, maxfev=10000) # params, covariance
                    if np.all(np.isfinite(c)): # possible return of curvefit: c=inf
                        s = np.sqrt(np.diag(c))
                    else:
                        s = 10.*np.abs(p)
                except:
                    s = 10.*np.abs(p)
                locp += [p]
                locs += [s]
                # if ((s[1]/p[1])<0.5) & (p[1] > 0.): pdb.set_trace()
    if len(locp) == 0:
        # raise ValueError('Error _nee2gpp_reichstein: No local relationship found.')
        print('Warning _nee2gpp_reichstein: No local relationship found.')
        return dfout
    locp   = np.squeeze(np.array(locp).astype(np.float))
    locs   = np.squeeze(np.array(locs).astype(np.float))
    # 2. E0 = avg of best 3
    # Reichstein et al. (2005), p. 1430, 1st paragraph.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        iii  = np.where((locp[:,1] > 0.) & (locp[:,1] < 450.) & (np.abs(locs[:,1]/locp[:,1]) < 0.5))[0]
    niii = iii.size
    if niii==0:
        # raise ValueError('Error _nee2gpp_reichstein: No good local relationship found.')
        # loosen the criteria: take the best three estimates anyway
        iii   = np.where((locp[:,1] > 0.))[0]
        niii = iii.size
        if niii<1:
            # raise ValueError('Error _nee2gpp_reichstein: No E0>0 found.')
            print('Warning _nee2gpp_reichstein: No E0>0 found.')
            return dfout
        lp    = locp[iii,:]
        ls    = locs[iii,:]
        iis   = np.argsort(ls[:,1])
        bestp = np.mean(lp[iis[0:np.minimum(3,niii)],:],axis=0)
        bests = np.mean(ls[iis[0:np.minimum(3,niii)],:],axis=0)
    elif niii==1:
        bestp = np.squeeze(locp[iii,:])
        bests = np.squeeze(locs[iii,:])
    elif niii==2:
        bestp = np.mean(locp[iii,:],axis=0)
        bests = np.mean(locs[iii,:],axis=0)
        # ls    = locs[iii,:]
        # iis   = np.argsort(ls[:,1])
    else:
        lp    = locp[iii,:]
        ls    = locs[iii,:]
        iis   = np.argsort(ls[:,1])
        bestp = np.mean(lp[iis[0:3],:],axis=0)
        bests = np.mean(ls[iis[0:3],:],axis=0)

    # 3. Refit Rref with fixed E0, each 4 days
    refp  = [] # Rref param
    refii = [] # mean index of data points
    E0    = bestp[1]
    et    = lloyd_fix(tt, 1., E0)
    for i in range(dmin,dmax,4):
        iii  = np.where((jul>=i) & (jul<(i+4)))[0]
        niii = iii.size
        if niii > 3:
            # Calc directly minisation of (nee-p*et)**2
            p, temp1, temp2 = opt.fmin_tnc(cost_abs, [2.],
                                           bounds=[[0.,None]],
                                           args=(lloyd_only_rref_p, et[iii], net[iii]),
                                           approx_grad=True, disp=False)
            refp  += [p]
            refii += [np.int((iii[0]+iii[-1])//2)]
    if len(refp) == 0:
        # raise ValueError('Error _nee2gpp_reichstein: No ref relationship found.')
        print('Warning _nee2gpp_reichstein: No ref relationship found.')
        return dfout
    refp  = np.squeeze(np.array(refp))
    refii = np.squeeze(np.array(refii))

    # 4. Interpol Rref
    Rref = np.interp(dates, jul[refii], refp)

    # 5. Calc Reco
    Reco     = np.ones(ndata)*undef
    ii       = np.where(~t.mask)[0]
    Reco[ii] = lloyd_fix(t[ii], Rref[ii], E0)

    # 6. Calc GPP
    GPP     = np.ones(ndata)*undef
    ii      = np.where(~(t.mask | nee.mask))[0]
    GPP[ii] = Reco[ii] - nee[ii]

    # 7. Set GPP=0 at night, if wanted
    if nogppnight:
        mask = misday | nee.mask | t.mask | misday.mask # night
        ii   = np.where(~mask)[0]
        Reco[ii] = nee[ii]
        GPP[ii]  = 0.
        # and prohibit negative gpp at any time
        mask = nee.mask | t.mask | (GPP>0.)
        ii   = np.where(~mask)[0]
        Reco[ii] -= GPP[ii]
        GPP[ii]  = 0.

    dfout = pd.DataFrame({'GPP':GPP, 'RECO':Reco}, index=df.index)

    return dfout


# ----------------------------------------------------------------------
def _nee2gpp_lasslop(df, ff, isday, undef=-9999, nogppnight=False):
    """
    Calculate photosynthesis (GPP) and ecosystem respiration (RECO) from original
    Eddy flux data, using the daytime method of Lasslop et al. (2010),
    in order to calculate Reco and then GPP = Reco - NEE.

    Parameters
    ----------
    df : pandas.Dataframe or numpy.array
        time series of CO2 fluxes, air temperature,
        incoming shortwave radiation, and air vapour pressure deficit.

        `df` can be a pandas.Dataframe with the columns
        'FC' or 'NEE' (or starting with 'FC\_' or 'NEE\_') for observed CO2 flux [umol(CO2) m-2 s-1]
        'TA'    (or starting with 'TA\_') for air temperature [K]
        'SW_IN' (or starting with 'SW_IN') for incoming short-wave radiation [W m-2]
        'VPD'   (or starting with 'VPD') for air vapour deficit [Pa]
        The index is taken as date variable.
    ff : pandas.Dataframe or numpy.array, optional
        flag Dataframe or array has the same shape as `df`. Non-zero values in
        `ff` will be treated as missing values in `df`.

        `ff` must follow the same rules as `df` if pandas.Dataframe.
    isday : array_like of bool, optional
        True when it is day, False when night. Must have the same length as `df.shape[0]`.
    undef : float, optional
        values having `undef` value are treated as missing values in `df` (default: -9999)
    nogppnight : float, optional
        GPP will be set to zero at night. RECO will then equal NEE at night (default: False)

    Returns
    -------
    pandas.Dataframe
        pandas.Dataframe with two columns 'GPP' and 'RECO' with estimated
        photosynthesis and ecosystem respiration.

    Literature
    ----------
    Lasslop et al. (2010)
        Separation of net ecosystem exchange into assimilation and respiration using
        a light response curve approach: critical issues and global evaluation,
        Global Change Biology 16, 187-208

    Examples
    --------
    >>> from fread import fread
    >>> from date2dec import date2dec
    >>> from dec2date import dec2date
    >>> ifile = 'test_nee2gpp.csv'
    >>> undef = -9999.
    >>> dat   = fread(ifile, skip=2, transpose=True)
    >>> ndat  = dat.shape[1]
    >>> head  = fread(ifile, skip=2, header=True)
    >>> head1 = head[0]
    >>> # date
    >>> jdate = date2dec(dy=dat[0,:], mo=dat[1,:], yr=dat[2,:], hr=dat[3,:], mi=dat[4,:])
    >>> adate = dec2date(jdate, eng=True)
    >>> # colhead
    >>> idx   = []
    >>> for i in head1:
    ...     if i in ['NEE', 'rg', 'Tair', 'VPD']: idx.append(head1.index(i))
    >>> colhead = ['FC', 'SW_IN', 'TA', 'VPD']
    >>> # data
    >>> dfin = dat[idx,:]
    >>> dfin[2,:] = np.where(dfin[2,:] == undef, undef, dfin[2,:]+273.15)
    >>> dfin[3,:] = np.where(dfin[3,:] == undef, undef, dfin[3,:]*100.)
    >>> # flag
    >>> flag = np.where(dfin == undef, 2, 0)
    >>> # partition
    >>> GPP, Reco = nee2gpp(dfin, flag=flag, date=adate, colhead=colhead, undef=undef, method='day')
    >>> print(GPP[1120:1128])
    [-9.99900000e+03 -9.99900000e+03 -9.99900000e+03  2.78457540e+00
      6.63212545e+00  8.88902165e+00  6.74243873e+00  9.51364527e+00]
    >>> print(Reco[1120:1128])
    [0.28786696 0.34594516 0.43893276 0.5495954  0.70029545 0.90849165
     1.15074873 1.46137527]

    History
    -------
    Written  Matthias Cuntz, Mar 2012
    Modified Arndt Piayda,   Mar 2012 - undef=np.nan
             Matthias Cuntz, Nov 2012 - individual routine
             Matthias Cuntz, Feb 2013 - ported to Python 3
    """
    # Variables
    fc_id = ''
    for cc in df.columns:
        if cc.startswith('FC_') or (cc == 'FC') or cc.startswith('NEE_') or (cc == 'NEE'):
            fc_id = cc
            break
    ta_id = ''
    for cc in df.columns:
        if cc.startswith('TA_') or (cc == 'TA'):
            ta_id = cc
            break
    sw_id = ''
    for cc in df.columns:
        if cc.startswith('SW_IN_') or (cc == 'SW_IN'):
            sw_id = cc
            break
    vpd_id = ''
    for cc in df.columns:
        if cc.startswith('VPD_') or (cc == 'VPD'):
            vpd_id = cc
            break
    assert fc_id,  'Carbon net flux with name FC or NEE or starting with FC_ or NEE_ must be in input.'
    assert ta_id,  'Air temperature with name TA or starting with TA_ must be in input.'
    assert sw_id,  'Global radiation with name SW or starting with SW_ must be in input.'
    assert vpd_id, 'Vapour pressure deficit with name VPD or starting with VPD_ must be in input.'

    nee    = np.ma.array(df[fc_id],  mask=(ff[fc_id] > 0))
    t      = np.ma.array(df[ta_id],  mask=(ff[ta_id] > 0))
    sw     = np.ma.array(df[sw_id],  mask=(ff[sw_id] > 0))
    vpd    = np.ma.array(df[vpd_id], mask=(ff[vpd_id] > 0))
    misday = np.ma.array(isday, mask=((~np.isfinite(isday)) | (isday == undef)))
    dates  = df.index.to_julian_date()

    # Partition - Lasslop et al. (2010) method

    ndata = nee.size
    GPP   = np.ones(ndata)*undef
    Reco  = np.ones(ndata)*undef
    dfout = pd.DataFrame({'GPP':GPP, 'RECO':Reco}, index=df.index)

    do_lgpp = False
    mask  = nee.mask | t.mask | misday.mask | sw.mask | vpd.mask
    # night
    nmask = misday | mask
    nii   = np.ma.where(~nmask)[0]
    njul  = dates[nii]
    ntt   = np.ma.compressed(t[nii])
    nnet  = np.ma.compressed(nee[nii])
    aRref = np.mean(nnet)
    # day
    dmask = (~misday) | mask
    dii   = np.ma.where(~dmask)[0]
    djul  = dates[dii]
    dtt   = np.ma.compressed(t[dii])
    dnet  = np.ma.compressed(nee[dii])
    dsw   = np.ma.compressed(sw[dii])
    dvpd  = np.ma.compressed(vpd[dii])
    # starting values for optim
    aalpha = 0.01
    qnet   = np.sort(dnet)
    nqnet  = qnet.size
    abeta0 = np.abs(qnet[np.floor(0.97*nqnet).astype(np.int)]-qnet[np.ceil(0.03*nqnet).astype(np.int)])
    ak     = 0.
    # out
    lE0    = []
    lalpha = []
    if do_lgpp:
        lbeta0 = []
        lk     = []
    lRref  = []
    lii    = []
    dmin = np.floor(np.amin(dates)).astype(np.int)
    dmax = np.ceil(np.amax(dates)).astype(np.int)
    zaehl = -1
    for i in range(dmin,dmax,2):
        good = True
        # 1. Estimate E0 from nighttime data
        iii  = np.squeeze(np.where((njul>=i) & (njul<(i+12))))
        niii = iii.size
        if niii > 3:
            p, temp1, temp2 = opt.fmin_tnc(cost_abs, [aRref, 100.],
                                           bounds=[[0.,None], [0.,None]],
                                           args=(lloyd_fix_p, ntt[iii], nnet[iii]),
                                           approx_grad=True, disp=False)
            E0 = np.maximum(p[1], 50.)
        else:
            if zaehl >= 0:
                E0 = lE0[zaehl]
            else:
                # large gap at beginning of data set, i.e. skip the period
                good = False
                continue
        # 2. Estimate alpha, k, beta0, Rref from daytime data
        iii  = np.squeeze(np.where((djul>=i) & (djul<(i+4))))
        niii = iii.size
        if niii > 3:
            et     = lloyd_fix(dtt[iii], 1., E0)
            again  = True
            ialpha = aalpha
            ibeta0 = abeta0
            ik     = ak
            iRref  = aRref
            bounds = [[None,None], [None,None], [None,None], [None,None]]
            while again:
                again = False
                p, nfeval, rc  = opt.fmin_tnc(cost_lasslop, [ialpha, ibeta0, ik, iRref],
                                              bounds=bounds,
                                              args=(dsw[iii], et, dvpd[iii], dnet[iii]),
                                              approx_grad=True, disp=False)
                # if parameters beyond some bounds, set params and redo the optim or skip
                if ((p[0] < 0.) | (p[0] > 0.22)): # alpha
                    again = True
                    if zaehl >= 0:
                        bounds[0] = [lalpha[zaehl],lalpha[zaehl]]
                        ialpha    = lalpha[zaehl]
                    else:
                        bounds[0] = [0.,0.]
                        ialpha    = 0.
                if p[1] < 0.:                     # beta0
                    bounds[1] = [0.,0.]
                    ibeta0    = 0.
                    again = True
                if p[1] > 250.:
                    good = False
                    continue
                if p[2] < 0.:                     # k
                    bounds[2] = [0.,0.]
                    ik        = 0.
                    again = True
                if p[3] < 0:                      # Rref
                    good = False
                    continue
            if good:
                lalpha = lalpha + [p[0]]
                if do_lgpp:
                    lbeta0 = lbeta0 + [p[1]]
                    lk     = lk     + [p[2]]
                lRref  = lRref  + [p[3]]
                lii    = lii    + [np.int((iii[0]+iii[-1])/2)]
            else:
                continue
        else:
            continue
        lE0    = lE0 + [E0]
        zaehl += 1
    if len(lE0) == 0:
        # raise ValueError('Error _nee2gpp_lasslop: No day relationship found.')
        print('Warning _nee2gpp_lasslop: No day relationship found.')
        return dfout
    lE0 = np.squeeze(np.array(lE0))
    if do_lgpp:
        lalpha = np.squeeze(np.array(lalpha))
        lbeta0 = np.squeeze(np.array(lbeta0))
        lk     = np.squeeze(np.array(lk))
    lRref  = np.squeeze(np.array(lRref))
    lii    = np.squeeze(np.array(lii))

    # 3. Interpol E0 and Rref
    E0   = np.interp(dates, djul[lii], lE0)
    Rref = np.interp(dates, djul[lii], lRref)

    # 4. Calc Reco
    Reco     = np.ones(ndata)*undef
    ii       = np.squeeze(np.where(~t.mask))
    Reco[ii] = lloyd_fix(t[ii], Rref[ii], E0[ii])

    # 5. Calc GPP from light response for check
    if do_lgpp:
        alpha    = np.interp(dates, djul[lii], lE0)
        beta0    = np.interp(dates, djul[lii], lbeta0)
        k        = np.interp(dates, djul[lii], lk)
        et       = lloyd_fix(t, 1., E0)
        lmask    = t.mask | misday.mask | sw.mask | vpd.mask
        ii       = np.squeeze(np.where(~lmask))
        lgpp     = np.zeros(ndata)
        lgpp[ii] = lasslop(sw[ii], et[ii], vpd[ii], alpha[ii], beta0[ii], k[ii], Rref[ii]) - Reco[ii]

    # 6. GPP
    GPP     = np.ones(ndata)*undef
    ii      = np.squeeze(np.where(~(t.mask | nee.mask)))
    GPP[ii] = Reco[ii] - nee[ii]

    # 7. Set GPP=0 at night, if wanted
    if nogppnight:
        mask = misday | nee.mask | t.mask | misday.mask # night
        ii   = np.where(~mask)[0]
        Reco[ii] = nee[ii]
        GPP[ii]  = 0.
        # and prohibit negative gpp at any time
        mask = nee.mask | t.mask | (GPP>0.)
        ii   = np.where(~mask)[0]
        Reco[ii] -= GPP[ii]
        GPP[ii]   = 0.

    dfout = pd.DataFrame({'GPP':GPP, 'RECO':Reco}, index=df.index)

    return dfout


# -------------------------------------------------------------
if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # from fread import fread
    # from date2dec import date2dec
    # from dec2date import dec2date
    # ifile = 'test_nee2gpp.csv'
    # undef = -9999.
    # # Day,Month,Year,Hour,Minute,NEE,rg,Tair,VPD,GPP_f,Reco
    # # -,-,-,-,-,umolm-2s-1,Wm-2,degC,hPa,umol_m-2_s-1,umol_m-2_s-1
    # # 1,1,2006,0,15,-9999,0,5.235,0.14192,-0.42485,2.89716
    # dat   = fread(ifile, skip=2, transpose=True)
    # ndat  = dat.shape[1]
    # head  = fread(ifile, skip=2, header=True)
    # head1 = head[0]
    # # date
    # jdate = date2dec(dy=dat[0,:], mo=dat[1,:], yr=dat[2,:], hr=dat[3,:], mi=dat[4,:])
    # adate = dec2date(jdate, eng=True)
    # # colhead
    # idx   = []
    # for i in head1:
    #     if i in ['NEE', 'rg', 'Tair', 'VPD']: idx.append(head1.index(i))
    # colhead = ['FC', 'SW_IN', 'TA', 'VPD']
    # # data
    # dfin = dat[idx,:]
    # dfin[2,:] = np.where(dfin[2,:] == undef, undef, dfin[2,:]+273.15)
    # dfin[3,:] = np.where(dfin[3,:] == undef, undef, dfin[3,:]*100.)
    # # flag
    # flag = np.where(dfin == undef, 2, 0)
    # # isday = np.where(dfin[1,:] > 10., True, False)
    # GPP, Reco = nee2gpp(dfin, flag=flag, date=adate, colhead=colhead, undef=undef, method='local')
    # print(GPP[1120:1128])
    # # [-9.99900000e+03 -9.99900000e+03 -9.99900000e+03  4.40606871e+00
    # #   8.31942152e+00  1.06242542e+01  8.49245664e+00  1.12381973e+01]
    # print(Reco[1120:1128])
    # # [1.68311981 1.81012431 1.9874173  2.17108871 2.38759152 2.64372415
    # #  2.90076664 3.18592735]
    # GPP, Reco = nee2gpp(dfin, flag=flag, date=adate, colhead=colhead, undef=undef, method='local')
    # print(GPP[1120:1128])
    # # [-9.99900000e+03 -9.99900000e+03 -9.99900000e+03  4.40606871e+00
    # #   8.31942152e+00  1.06242542e+01  8.49245664e+00  1.12381973e+01]
    # GPP, Reco = nee2gpp(dfin, flag=flag, date=adate, colhead=colhead, undef=undef, method='global')
    # print(GPP[1120:1128])
    # # [-9.99900000e+03 -9.99900000e+03 -9.99900000e+03  4.33166157e+00
    # #   8.18228013e+00  1.04092252e+01  8.19395317e+00  1.08427448e+01]
    # GPP, Reco = nee2gpp(dfin, flag=flag, date=adate, colhead=colhead, undef=undef, method='Reichstein')
    # print(GPP[1120:1128])
    # # [-9.99900000e+03 -9.99900000e+03 -9.99900000e+03  4.406068706013192
    # #   8.319421516040766 10.624254150217764  8.492456637225963 11.238197347837367]
    # GPP, Reco = nee2gpp(dfin, flag=flag, date=adate, colhead=colhead, undef=undef, method='reichstein')
    # print(GPP[1120:1128])
    # # [-9.99900000e+03 -9.99900000e+03 -9.99900000e+03  4.40606871e+00
    # #   8.31942152e+00  1.06242542e+01  8.49245664e+00  1.12381973e+01]
    # GPP, Reco = nee2gpp(dfin, flag=flag, date=adate, colhead=colhead, undef=undef, method='day')
    # print(GPP[1120:1128])
    # # [-9.99900000e+03 -9.99900000e+03 -9.99900000e+03  2.78457540e+00
    # #   6.63212545e+00  8.88902165e+00  6.74243873e+00  9.51364527e+00]
    # print(Reco[1120:1128])
    # # [0.28786696 0.34594516 0.43893276 0.5495954  0.70029545 0.90849165
    # #  1.15074873 1.46137527]
