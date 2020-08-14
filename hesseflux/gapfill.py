#!/usr/bin/env python
"""
gapfill : Fills gaps of flux data from Eddy covariance measurements according
          to Reichstein et al. (Global Change Biology, 2005) or estimate flux
          uncertainties after Lasslop et al. (Biogeosciences, 2008).

This module was written by Matthias Cuntz while at Department of
Computational Hydrosystems, Helmholtz Centre for Environmental
Research - UFZ, Leipzig, Germany, and continued while at Institut
National de Recherche pour l'Agriculture, l'Alimentation et
l'Environnement (INRAE), Nancy, France.

Copyright (c) 2012-2020 Matthias Cuntz - mc (at) macu (dot) de
Released under the MIT License; see LICENSE file for details.

* Written Mar 2012 by Matthias Cuntz - mc (at) macu (dot) de
* Ported to Python 3, Feb 2013, Matthias Cuntz
* Input data can be ND-array, Apr 2014, Matthias Cuntz
* Bug in longestmarginalgap: was only working at time series edges, rename it
  to longgap, Apr 2014, Matthias Cuntz
* Keyword fullday, Apr 2014, Matthias Cuntz
* Input can be pandas Dataframe or numpy array(s), Apr 2020, Matthias Cuntz
* Using numpy docstring format, May 2020, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz

The following functions are provided

.. autosummary::
   gapfill
"""
from __future__ import division, absolute_import, print_function
import numpy as np
import pandas as pd


__all__ = ['gapfill']


def gapfill(dfin, flag=None, date=None, timeformat='%Y-%m-%d %H:%M:%S',
            colhead=None,
            sw_dev=50., ta_dev=2.5, vpd_dev=5.,
            longgap=60, fullday=False, undef=-9999, ddof=1,
            err=False, verbose=0):
    """
    Fills gaps in flux data from Eddy covariance measurements with
    Marginal Distribution Sampling (MDS) according to Reichstein et al.
    (Global Change Biology, 2005).

    This means, if there is a gap in the data, look for similar meteorological
    conditions (defined as maximum possible deviations) in a certain time
    window and fill with the average of these 'similar' values.

    The routine can also do the same search for similar meteorological
    conditions for every data point and calculate its standard deviation as a
    measure of uncertainty after Lasslop et al. (Biogeosciences, 2008).

    Parameters
    ----------
    dfin : pandas.Dataframe or numpy.array
        time series of fluxes to fill as well as
        meteorological variables incoming short-wave radiation,
        air temperature, air vapour pressure deficit.

        `dfin` can be a pandas.Dataframe with the columns
        'SW_IN' (or starting with 'SW_IN') for incoming short-wave radiation [W m-2]
        'TA'    (or starting with 'TA\_') for air temperature [deg C]
        'VPD'   (or starting with 'VPD') for air vapour deficit [hPa]
        and columns with ecosystem fluxes with possible missing values (gaps).
        The index is taken as date variable.

        `dfin` can also me a numpy array with the same columns. In this case
        `colhead`, `date`, and possibly `dateformat` must be given.
    flag : pandas.Dataframe or numpy.array, optional
        flag Dataframe or array has the same shape as dfin. Non-zero values in
        `flag` will be treated as missing values in `dfin`.

        `flag` must follow the same rules as `dfin` if pandas.Dataframe.

        If `flag` is numpy array, `df.columns.values` will be used as
        column heads and the index of `dfin` will be copied to `flag`.
    date : array_like of string, optional
        1D-array_like of calendar dates in format given in `timeformat`.

        `date` must be given if `dfin` is numpy array.
    timeformat : str, optional
        Format of dates in `date`, if given (default: '%Y-%m-%d %H:%M:%S').
        See strftime documentation of Python's datetime module:
        https://docs.python.org/3/library/datetime.html#strftime-and-strptime-behavior
    colhed : array_like of str, optional
        column names if `dfin` is numpy array. See `dfin` for mandatory
        column names.
    sw_dev : float, optional
        threshold for maximum deviation of global radiation (default: 50)
    ta_dev : float, optional
        threshold for maximum deviation of air Temperature (default: 2.5)
    vpd_dev : float, optional
        threshold for maximum deviation of vpd (default: 5.)
    longgap : int, optional
        avoid extraploation into a gap longer than `longgap` days (default: 60)
    fullday : bool, optional
        True: move beginning of large gap to start of next day and move end of
              large gap to end of last day (default: False)
    undef : float, optional
        values having `undef` value are treated as missing values in `dfin`
        (default: -9999)

        np.nan is not allowed (working).
    ddof : int, optional
        Delta Degrees of Freedom. The divisor used in calculation of standard
        deviation for error estimates (`err=True`) is ``N-ddof``, where ``N``
        represents the number of elements (default: 1).
    err : bool, optional
        True: fill every data point, i.e. used for error generation
        (default: False)
    shape : bool or tuple, optional
        True: output have the same shape as input data if `dfin` is
        numpy array; if a tuple is given, then this tuple is used to reshape.

        False: outputs are 1D arrays if `dfin` is numpy array (default: False).
    verbose : int, optional
        Verbosity level 0-3 (default: 0). 0 is no output; 3 is very verbose.

    Returns
    -------
    pandas.Dataframe(s) or numpy array(s)
        If `err==False`: filled_data, quality_class

        If `err==True`:  err_data

        pandas.Dataframe(s) will be returned if `dfin` was Dataframe.

        numpy array(s) will be returned if `dfin` was numpy array.

    Notes
    -----
    If `err==True`, there is no error estimate if there are no meteorological
    conditions in the vicinity of the data point.

    Routine Does not work with `undef=np.nan`.

    Reichstein et al. (2005)
        On the separation of net ecosystem exchange into assimilation and
        ecosystem respiration: review and improved algorithm.
        Global Change Biology 11, 1424-1439

    Examples
    --------
    >>> import numpy as np
    >>> from fread import fread
    >>> from date2dec import date2dec
    >>> from dec2date import dec2date
    >>> ifile = 'test_gapfill.csv' # Tharandt 1998 = Online tool example file
    >>> undef = -9999.
    >>> # data
    >>> dat   = fread(ifile, skip=2, transpose=True)
    >>> ndat  = dat.shape[1]
    >>> head  = fread(ifile, skip=2, header=True)
    >>> head1 = head[0]
    >>> # colhead
    >>> idx   = []
    >>> for i in head1:
    ...     if i in ['NEE', 'LE', 'H', 'Rg', 'Tair', 'VPD']:
    ...         idx.append(head1.index(i))
    >>> colhead = ['FC', 'LE', 'H', 'SW_IN', 'TA', 'VPD']
    >>> # data
    >>> dfin = dat[idx,:]
    >>> # flag
    >>> flag = np.where(dfin == undef, 2, 0)
    >>> flag[0,:] = dat[head1.index('qcNEE'),:].astype(np.int)
    >>> flag[1,:] = dat[head1.index('qcLE'),:].astype(np.int)
    >>> flag[2,:] = dat[head1.index('qcH'),:].astype(np.int)
    >>> flag[np.where(flag==1)] = 0
    >>> # date
    >>> day_id  = head1.index('Day')
    >>> hour_id = head1.index('Hour')
    >>> ntime   = dat.shape[1]
    >>> year  = np.ones(ntime, dtype=np.int) * 1998
    >>> hh    = dat[hour_id,:].astype(np.int)
    >>> mn    = np.rint((dat[hour_id,:]-hh)*60.).astype(np.int)
    >>> y0    = date2dec(yr=year[0], mo=1, dy=1, hr=hh, mi=mn)
    >>> jdate = y0 + dat[day_id,:]
    >>> adate = dec2date(jdate, eng=True)
    >>> # fill
    >>> dat_f, flag_f = gapfill(dfin, flag=flag, date=adate, colhead=colhead,
    ...                         undef=undef, verbose=0)
    >>> print('{:d} {:d} {:d} {:d} {:d} {:d}'.format(*flag_f[0,11006:11012]))
    1 1 1 2 2 2
    >>> print('{:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}'.format(*dat_f[0,11006:11012]))
    -18.68 -15.63 -19.61 -15.54 -12.40 -15.33

    >>> # 1D err
    >>> dat_std = gapfill(dfin, flag=flag, date=adate, colhead=colhead,
    ...                   undef=undef, verbose=0, err=True)
    >>> print('{:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f}'.format(*dat_std[0,11006:11012]))
    5.372 13.118 6.477 -9999.000 -9999.000 -9999.000

    >>> dat_err     = np.ones(ndat, dtype=np.int)*(-1)
    >>> kk          = np.where((dat_std[0,:] != undef) & (dat_f[0,:] != 0.))[0]
    >>> dat_err[kk] = np.abs(dat_std[0,kk]/dat_f[0,kk]*100.).astype(np.int)
    >>> print('{:d} {:d} {:d} {:d} {:d} {:d}'.format(*dat_err[11006:11012]))
    28 83 33 -1 -1 -1

    History
    -------
    Written,  Matthias Cuntz, Mar 2012 - modified gap_filling.py
    Modified, Matthias Cuntz, Feb 2013 - ported to Python 3
              Matthias Cuntz, Apr 2014 - assert
                                       - data ND-array
                                       - longestmarginalgap was only working at
                                         beginning and end of time series
                                         renamed to longgap
                                       - fullday
              Matthias Cuntz, Apr 2020 - Input can be pandas Dataframe or
                                         numpy array(s)
              Matthias Cuntz, May 2020 - numpy docstring format
    """
    # Check input
    # numpy or panda
    if isinstance(dfin, (np.ndarray, np.ma.MaskedArray)):
        isnumpy = True
        istrans = False
        astr = 'colhead must be given if input is numpy.ndarray.'
        assert colhead is not None, astr
        if dfin.shape[0] == len(colhead):
            istrans = True
            df = pd.DataFrame(dfin.T, columns=colhead)
        elif dfin.shape[1] == len(colhead):
            df = pd.DataFrame(dfin, columns=colhead)
        else:
            estr = 'Length of colhead must be number of columns in input'
            estr = estr + ' array. len(colhead)=' + str(len(colhead))
            estr = estr + ' shape(input)=(' + str(dfin.shape[0])
            estr = estr + ',' + str(dfin.shape[1]) + ').'
            raise ValueError(estr)
        assert date is not None, 'date must be given if input is numpy arrary.'
        df['Datetime'] = pd.to_datetime(date, format=timeformat)
        df.set_index('Datetime', drop=True, inplace=True)
    else:
        isnumpy = False
        istrans = False
        astr = 'Input must be either numpy.ndarray or pandas.DataFrame.'
        assert isinstance(dfin, pd.core.frame.DataFrame), astr
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
                estr = 'flag must have same shape as data array. data: ({:d},{:d}); flag: ({:d},{:d})'.format(dfin.shape[0], dfin.shape[1], flag.shape[0], flag.shape[1])
                raise ValueError(estr)
            ff = ff.set_index(df.index)
        else:
            fisnumpy = False
            fistrans = False
            astr = 'Flag must be either numpy.ndarray or pandas.DataFrame.'
            assert isinstance(flag, pd.core.frame.DataFrame), astr
            ff = flag.copy(deep=True)
    else:
        fisnumpy = isnumpy
        fistrans = istrans
        # flags: 0: good; 1: input flagged; 2: output flagged
        ff              = df.copy(deep=True).astype(int)
        ff[:]           = 0
        ff[df == undef] = 1
        ff[df.isna()]   = 1

    # Data and flags
    sw_id = ''
    for cc in df.columns:
        if cc.startswith('SW_IN_') or (cc == 'SW_IN'):
            sw_id = cc
            break
    ta_id = ''
    for cc in df.columns:
        if cc.startswith('TA_') or (cc == 'TA'):
            ta_id = cc
            break
    vpd_id = ''
    for cc in df.columns:
        if cc.startswith('VPD_') or (cc == 'VPD'):
            vpd_id = cc
            break
    astr = 'Global radiation with name SW or starting with SW_'
    astr = astr + ' must be in input.'
    assert sw_id,  astr
    astr = 'Air temperature with name TA or starting with TA_'
    astr = astr + ' must be in input.'
    assert ta_id,  astr
    astr = 'Vapour pressure deficit with name VPD or starting'
    astr = astr + ' with VPD_ must be in input.'
    assert vpd_id, astr

    sw      = df[sw_id].to_numpy()
    sw_flg  = ff[sw_id].to_numpy()
    ta      = df[ta_id].to_numpy()
    ta_flg  = ff[ta_id].to_numpy()
    vpd     = df[vpd_id].to_numpy()
    vpd_flg = ff[vpd_id].to_numpy()

    dfill    = df.copy(deep=True)
    ffill    = ff.copy(deep=True)
    ffill[:] = 0

    # Times
    # number of data points per week; basic factor of the time window
    week    = pd.Timedelta('1 W') / (df.index[1] - df.index[0])
    nperday = week // 7
    hour    = df.index.hour + df.index.minute/60.
    day     = (df.index.to_julian_date()-0.5).astype(np.int)

    # Filling variables
    ndata = len(df)
    for hcol in df.columns:

        if hcol.startswith('SW_IN_') or (hcol == 'SW_IN'):
            continue
        if hcol.startswith('TA_')    or (hcol == 'TA'):
            continue
        if hcol.startswith('VPD_')   or (hcol == 'VPD'):
            continue

        if verbose > 0:
            if err:
                print('  Error estimate ', str(hcol))
            else:
                print('  Filling ', str(hcol))

        data  = df[hcol].to_numpy()
        dflag = ff[hcol].to_numpy()

        data_f  = dfill[hcol].to_numpy()
        dflag_f = ffill[hcol].to_numpy()

        # Large margins

        # Check for large margins at beginning
        largegap   = np.zeros(ndata, dtype=np.bool)
        firstvalid = np.amin(np.where(dflag == 0)[0])
        lastvalid  = np.amax(np.where(dflag == 0)[0])
        nn         = int(nperday*longgap)
        if firstvalid > nn:
            if verbose > 1:
                print('    Large margin at beginning: ', firstvalid)
            largegap[0:(firstvalid-nn)] = True
        if lastvalid < (ndata-nn):
            if verbose > 1:
                print('    Large margin at end: ', lastvalid-nn)
            largegap[(lastvalid+nn):] = True

        # Large gaps

        # search largegap - code from maskgroup.py
        index  = []
        length = []
        count  = 0
        for i in range(ndata):
            if i == 0:
                if dflag[i] != 0:
                    index += [i]
                    count  = 1
            if i > 0:
                if (dflag[i] != 0) and (dflag[i-1] == 0):
                    index += [i]
                    count  = 1
                elif dflag[i] != 0:
                    count += 1
                elif (dflag[i] == 0) and (dflag[i-1] != 0):
                    length += [count]
                    count = 0
                else:
                    pass
        if count > 0:
            length += [count]

        # set largegap
        for i in range(len(index)):
            if length[i] > nn:
                if verbose > 1:
                    print('    Large gap: ', index[i], ':', index[i]+length[i])
                largegap[index[i]:index[i]+length[i]] = True

        # set or unset rest of days in large gaps
        if fullday:
            for i in range(ndata-1):
                # end of large margin
                if largegap[i] and not largegap[i+1]:
                    largegap[np.where(day == day[i])[0]] = False
                # beginning of large margin
                elif not largegap[i] and largegap[i+1]:
                    largegap[np.where(day == day[i])[0]] = False
                else:
                    continue

        # Gap filling

        # flag for all meteorological conditions
        meteo_flg = (ta_flg == 0) & (vpd_flg == 0) & (sw_flg == 0)
        # flag for all meteorological conditions and data
        total_flg = meteo_flg & (dflag == 0)

        # Fill loop over all data points
        for j in range(ndata):
            if not err:
                # no reason to go further, no gap -> continue
                if (dflag[j] == 0) | largegap[j]:
                    continue
            # 3 Methods
            #   1. ta, vpd and global radiation sw;
            #   2. just global radiation sw;
            #   3. no meteorolgical conditions: take the mean of +- hour

            # for better overview: dynamic calculation of radiation threshold
            # minimum 20; maximum 50 [Wm-2] according to private correspondence
            # with Markus Reichstein
            sw_devmax = np.maximum(20., np.minimum(sw[j], sw_dev))

            # Method 1: all met conditions
            if meteo_flg[j]:
                # search for values around the met-conditions
                # in a window of time
                # (one week in the first iteration and odd weeks in the next)
                j1  = j - np.arange(1, week+1, dtype=np.int) + 1
                j2  = j + np.arange(1, week, dtype=np.int)
                jj  = np.append(j1, j2)
                win = np.unique(np.sort(np.clip(jj, 0, ndata-1)))
                # get boolean array where meteo-conditions are in a given width
                conditions = ( (np.abs(sw[win]-sw[j])   < sw_devmax) &
                               (np.abs(ta[win]-ta[j])   < ta_dev) &
                               (np.abs(vpd[win]-vpd[j]) < vpd_dev) &
                               total_flg[win] )
                num4avg = np.sum(conditions)
                # we need at least two samples with similar conditions
                if num4avg >= 2:
                    dat = np.ma.array(data[win], mask=~conditions)
                    if verbose > 2:
                        print('    m1.1: ', j, win.size, dat.mean(),
                              dat.std(ddof=ddof))
                    if err:
                        data_f[j]  = dat.std(ddof=ddof)
                    else:
                        data_f[j]  = dat.mean()
                        # assign also quality category of gap filling
                        dflag_f[j] = 1
                    continue
                else:  # --> extend time window to two weeks
                    j1  = j - np.arange(1, 2*week+1, dtype=np.int) + 1
                    j2  = j + np.arange(1, 2*week, dtype=np.int)
                    jj  = np.append(j1, j2)
                    win = np.unique(np.sort(np.clip(jj, 0, ndata-1)))
                    conditions = ( (np.abs(sw[win]  - sw[j])  < sw_devmax) &
                                   (np.abs(ta[win]  - ta[j])  < ta_dev) &
                                   (np.abs(vpd[win] - vpd[j]) < vpd_dev) &
                                   total_flg[win] )
                    num4avg = np.sum(conditions)
                    if num4avg >= 2:
                        dat = np.ma.array(data[win], mask=~conditions)
                        if verbose > 2:
                            print('    m1.2: ', j, win.size, dat.mean(),
                                  dat.std(ddof=ddof))
                        if err:
                            data_f[j]  = dat.std(ddof=ddof)
                        else:
                            data_f[j]  = dat.mean()
                            # assign also quality category of gap filling
                            dflag_f[j] = 1
                        continue

            if err:
                continue
            # if you come here, gap-filling rather than error estimate

            # If nothing is found under similar meteo within two weeks,
            # look for global radiation within one week ->

            # Method 2: just global radiation available
            if sw_flg[j] == 0:
                j1  = j - np.arange(1, week+1, dtype=np.int) + 1
                j2  = j + np.arange(1, week, dtype=np.int)
                jj  = np.append(j1, j2)
                win = np.unique(np.sort(np.clip(jj, 0, ndata-1)))
                # get boolean array where meteo-conditions are in a given width
                conditions = ( (np.abs(sw[win]-sw[j]) < sw_devmax) &
                               total_flg[win] )
                num4avg = np.sum(conditions)
                # we need at least two samples with similar conditions
                if num4avg >= 2:
                    dat = np.ma.array(data[win], mask=~conditions)
                    if verbose > 2:
                        print('    m2: ', j, win.size, dat.mean(),
                              dat.std(ddof=ddof))
                    data_f[j]  = dat.mean()
                    dflag_f[j] = 1
                    continue

            # If still nothing is found under similar sw within one week,
            # take the same hour within 1-7 days

            # Method 3: same hour
            enough = False
            for i in range(2):
                t_win = (nperday * (2*i+1))//2
                j1  = j - np.arange(1, t_win+1, dtype=np.int) + 1
                j2  = j + np.arange(1, t_win, dtype=np.int)
                jj  = np.append(j1, j2)
                win = np.unique(np.sort(np.clip(jj, 0, ndata-1)))
                conditions = ( (np.abs(hour[win]-hour[j]) < 1.1)
                               & (dflag[win] == 0) )
                num4avg = np.sum(conditions)
                if num4avg >= 2:
                    dat = np.ma.array(data[win], mask=~conditions)
                    if verbose > 2:
                        print('    m3.{:d}: '.format(i), j, win.size,
                              dat.mean(), dat.std(ddof=ddof))
                    data_f[j] = dat.mean()
                    if i == 0:
                        dflag_f[j] = 1
                    else:
                        dflag_f[j] = 2
                    break

            # sanity check
            if dflag_f[j] > 0:
                continue

            # If still nothing is found, start a new cycle
            # with increased window size
            # Method 4: same as 1 but for 3-12 weeks
            if meteo_flg[j]:
                for multi in range(3, 12):
                    j1  = j - np.arange(1, multi*week+1, dtype=np.int) + 1
                    j2  = j + np.arange(1, multi*week, dtype=np.int)
                    jj  = np.append(j1, j2)
                    win = np.unique(np.sort(np.clip(jj, 0, ndata-1)))
                    conditions = ( (np.abs(sw[win]  - sw[j])  < sw_devmax) &
                                   (np.abs(ta[win]  - ta[j])  < ta_dev) &
                                   (np.abs(vpd[win] - vpd[j]) < vpd_dev) &
                                   total_flg[win] )
                    num4avg = np.sum(conditions)
                    # we need at least two samples with similar conditions
                    if num4avg >= 2:
                        dat = np.ma.array(data[win], mask=~conditions)
                        if verbose > 2:
                            print('    m4.{:d}: '.format(multi), j, win.size,
                                  dat.mean(), dat.std(ddof=ddof))
                        data_f[j] = dat.mean()
                        # assign also quality category of gap filling
                        if multi <= 2:
                            dflag_f[j] = 1
                        elif multi > 4:
                            dflag_f[j] = 3
                        else:
                            dflag_f[j] = 2
                        break

                # Check because continue does not support
                # to jump out of two loops
                if dflag_f[j] > 0:
                    continue

            # Method 5: same as 2 but for 2-12 weeks
            if sw_flg[j] == 0:
                for multi in range(2, 12):
                    j1  = j - np.arange(1, multi*week+1, dtype=np.int) + 1
                    j2  = j + np.arange(1, multi*week, dtype=np.int)
                    jj  = np.append(j1, j2)
                    win = np.unique(np.sort(np.clip(jj, 0, ndata-1)))
                    # get boolean array where meteo-conditions are
                    # in a given width
                    conditions = ( (np.abs(sw[win] - sw[j]) < sw_devmax) &
                                   total_flg[win] )
                    num4avg = np.sum(conditions)
                    # we need at least two samples with similar conditions
                    if num4avg >= 2:
                        dat = np.ma.array(data[win], mask=~conditions)
                        if verbose > 2:
                            print('    m5.{:d}: '.format(multi), j, win.size,
                                  dat.mean(), dat.std(ddof=ddof))
                        data_f[j] = dat.mean()
                        if multi == 0:
                            dflag_f[j] = 1
                        elif multi <= 2:
                            dflag_f[j] = 2
                        else:
                            dflag_f[j] = 3
                        break

                if dflag_f[j] > 0:
                    continue

            # Method 6: same as 3 but for 3-120 days
            for i in range(3, 120):
                t_win = nperday * (2*i+1)/2
                j1  = j - np.arange(1, t_win+1, dtype=np.int) + 1
                j2  = j + np.arange(1, t_win, dtype=np.int)
                jj  = np.append(j1, j2)
                win = np.unique(np.sort(np.clip(jj, 0, ndata-1)))
                conditions = ( (np.abs(hour[win]-hour[j]) < 1.1)
                               & (dflag[win] == 0) )
                num4avg = np.sum(conditions)
                if num4avg >= 2:
                    dat = np.ma.array(data[win], mask=~conditions)
                    if verbose > 2:
                        print('    m6.{:d}: '.format(i), j, win.size,
                              dat.mean(), dat.std(ddof=ddof))
                    data_f[j]  = dat.mean()
                    dflag_f[j] = 3
                    break

        dfill[hcol] = data_f
        ffill[hcol] = dflag_f

    # Finish

    if isnumpy:
        if istrans:
            dfout = dfill.to_numpy().T
        else:
            dfout = dfill.to_numpy()
    else:
        dfout = dfill

    if fisnumpy:
        if fistrans:
            ffout = ffill.to_numpy().T
        else:
            ffout = ffill.to_numpy()
    else:
        ffout = ffill

    if err:
        return dfout
    else:
        return dfout, ffout


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # import numpy as np
    # from fread import fread
    # from date2dec import date2dec
    # from dec2date import dec2date
    # from autostring import astr
    # ifile = 'test_gapfill.csv' # Tharandt 1998 = Online tool example file
    # undef = -9999.
    # # Day Hour NEE         qcNEE  LE    qcLE  H       qcH  Rg    Tair  Tsoil  rH     VPD  Ustar
    # # --  --   umolm-2s-1  --     Wm-2  --    Wm-2    --   Wm-2  degC  degC   %      hPa  ms-1
    # # 1   0.5  -1.21       1      1.49  1     -11.77  1    0     7.4   4.19   55.27  4.6  0.72
    # dat   = fread(ifile, skip=2, transpose=True)
    # # dat = dat[:,:1000]
    # ndat  = dat.shape[1]
    # head  = fread(ifile, skip=2, header=True)
    # head1 = head[0]
    # # colhead
    # idx   = []
    # for i in head1:
    #     if i in ['NEE', 'LE', 'H', 'Rg', 'Tair', 'VPD']: idx.append(head1.index(i))
    # colhead = ['FC', 'LE', 'H', 'SW_IN', 'TA', 'VPD']
    # # data
    # dfin = dat[idx,:]
    # # flag
    # flag = np.where(dfin == undef, 2, 0)
    # flag[0,:] = dat[head1.index('qcNEE'),:].astype(np.int)
    # flag[1,:] = dat[head1.index('qcLE'),:].astype(np.int)
    # flag[2,:] = dat[head1.index('qcH'),:].astype(np.int)
    # flag[np.where(flag==1)] = 0
    # # date
    # day_id  = head1.index('Day')
    # hour_id = head1.index('Hour')
    # ntime   = dat.shape[1]
    # year  = np.ones(ntime, dtype=np.int) * 1998
    # hh    = dat[hour_id,:].astype(np.int)
    # mn    = np.rint((dat[hour_id,:]-hh)*60.).astype(np.int)
    # y0    = date2dec(yr=year[0], mo=1, dy=1, hr=hh, mi=mn)
    # jdate = y0 + dat[day_id,:]
    # adate = dec2date(jdate, eng=True)
    # # fill
    # dat_f, flag_f = gapfill(dfin, flag=flag, date=adate, colhead=colhead, undef=undef, verbose=0)
    # print(astr(flag_f[0,11006:11012],0,pp=True))
    # # ['1' '1' '1' '2' '2' '2']
    # print(astr(dat_f[0,11006:11012],3,pp=True))
    # # ['-18.678' '-15.633' '-19.610' '-15.536' '-12.402' '-15.329']

    # # 1D err
    # dat_std = gapfill(dfin, flag=flag, date=adate, colhead=colhead, undef=undef, verbose=0, err=True)
    # print(astr(dat_std[0,11006:11012],3,pp=True))
    # # ['    5.372' '   13.118' '    6.477' '-9999.000' '-9999.000' '-9999.000']

    # dat_err     = np.ones(ndat, dtype=np.int)*(-1)
    # kk          = np.where((dat_std[0,:] != undef) & (dat_f[0,:] != 0.))[0]
    # dat_err[kk] = np.abs(dat_std[0,kk]/dat_f[0,kk]*100.).astype(np.int)
    # print(astr(dat_err[11006:11012],pp=True))
    # # [' 28' ' 83' ' 33' ' -1' ' -1' ' -1']
