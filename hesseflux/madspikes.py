#!/usr/bin/env python
"""
madspikes : Spike detection for using a moving median absolute difference filter.

This module was original written by Tino Rau and Matthias Cuntz, and
maintained by Arndt Piayda while at Department of Computational
Hydrosystems, Helmholtz Centre for Environmental Research - UFZ,
Leipzig, Germany, and continued by Matthias Cuntz while at Institut
National de Recherche pour l'Agriculture, l'Alimentation et
l'Environnement (INRAE), Nancy, France.

Copyright (c) 2008-2020 Matthias Cuntz - mc (at) macu (dot) de
Released under the MIT License; see LICENSE file for details.

* Written 2008 by Tino Rau and Matthias Cuntz - mc (at) macu (dot) de
* Maintained by Arndt Piayda since Aug 2014.
* Input can be pandas Dataframe or numpy array(s), Apr 2020, Matthias Cuntz
* Removed iteration, Apr 2020, Matthias Cuntz
* Using numpy docstring format, May 2020, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz, Arndt Piayda, Tino Rau

The following functions are provided

.. autosummary::
   madspikes
"""
from __future__ import division, absolute_import, print_function
import numpy as np
import pandas as pd
try:    # import package
    from .mad import mad
except: # python madspikes.py
    from mad import mad


__all__ = ['madspikes']


def madspikes(dfin, flag=None, isday=None,
              colhead=None, undef=-9999,
              nscan=15*48, nfill=1*48,
              z=7, deriv=2, swthr=10.,
              plot=False):
    """
    Spike detection for using a moving median absolute difference filter.
    Used with Eddy vovariance data in Papale et al. (Biogeosciences, 2006).

    Parameters
    ----------
    dfin : pandas.Dataframe or numpy.array
        time series of data where spike detection with MAD should be applied.

        `dfin` can be a pandas.Dataframe.

        `dfin` can also me a numpy array. In this case `colhead` must be given.
        MAD will be applied along axis=0, i.e. on each column of axis=1.
    flag : pandas.Dataframe or numpy.array, optional
        flag Dataframe or array has the same shape as dfin. Non-zero values in
        `flag` will be treated as missing values in `dfin`.

        If `flag` is numpy array, `df.columns.values` will be used as column heads.
    isday : array_like of bool, optional
        True when it is day, False when night. Must have the same length as dfin.shape[0].

        If `isday` is not given, `dfin` must have a column with head 'SW_IN' or
        starting with 'SW_IN'. `isday` will then be `dfin['SW_IN'] > swthr`.
    colhed : array_like of str, optional
        column names if `dfin` is numpy array.
    undef : float, optional
        values having `undef` value are treated as missing values in `dfin` (default: -9999)

        np.nan is not allowed (working).
    nscan : int, optional
        size of moving window to calculate mad in time steps (default: 15*48)
    nfill : int, optional
        step size of moving window to calculate mad in time steps (default: 1*48)

        mad will be calculated in `nscan` time window. Resulting mask will be applied
        only in `nfill` window in the middle of the `nscan` window. Then `nscan` window
        will be moved by `nfill` time steps.
    z : float, optional
        Input is allowed to deviate maximum `z` standard deviations from the median (default: 7)
    deriv : int, optional
        0: Act on raw input.

        1: Use first derivatives.

        2: Use 2nd derivatives (default).
    swthr : float, optional
        Threshold to determine daytime from incoming shortwave radiation if `isday` not given (default: 10).
    plot : bool, optional
        True: data and spikes are plotted into madspikes.pdf (default: False).

    Returns
    -------
    pandas.Dataframe or numpy array
        flags, 0 everywhere except detected spikes set to 2.

    History
    -------
    Written,    Matthias Cuntz & Tino Rau, 2008
    Maintained, Arndt Piayda,   Aug 2014
    Modified,   Matthias Cuntz, Apr 2020 - input can be pandas Dataframe or numpy array(s)
                                         - removed iteration
                Matthias Cuntz, May 2020 - numpy docstring format
    """
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

    # parameters
    nrow, ncol = df.shape
    half_scan_win = nscan//2
    half_fill_win = nfill//2

    # calculate dusk and dawn times and separate in day and night
    isdawn         = np.zeros(nrow, dtype=np.bool)
    isdusk         = np.zeros(nrow, dtype=np.bool)
    dis            = isday.astype(int) - np.roll(isday,-1).astype(int) # .astype(bool)
    isdawn[:-1]    = np.where(dis[:-1] == -1, True, False)
    isdusk[:-1]    = np.where(dis[:-1] ==  1, True, False)
    isddday        = isdawn
    tmp            = np.roll(isdusk,1)
    isddday[1:]   += tmp[1:] # start and end of day
    isddnight      = isdusk
    tmp            = np.roll(isdawn,1)
    isddnight[1:] += tmp[1:] # start and end of night

    # iterate over each column of data
    if plot:
        import matplotlib.pyplot as plt
        import matplotlib.backends.backend_pdf as pdf
        pd.plotting.register_matplotlib_converters()
        pp  = pdf.PdfPages('madspikes.pdf')

    cols = list(df.columns)
    for hcol in df.columns:

        if hcol.startswith == 'SW_IN': continue

        data  = df[hcol]
        dflag = ff[hcol]

        # get day and night data
        data_day = data.copy(deep=True)
        data_day[~(isday | isddday) | (dflag != 0) | (data == undef)] = np.nan
        data_night = data.copy(deep=True)
        data_night[~(~isday | isddnight) | (dflag != 0) | (data == undef)] = np.nan

        # iterate over fill window
        for j in range(half_fill_win, nrow-1, 2*half_fill_win):
            j1 = max(j - half_scan_win - 1, 0)
            j2 = min(j + half_scan_win + 1, nrow)
            fill_start = max(j - half_fill_win, 1)
            fill_end   = min(j + half_fill_win, nrow-1)

            dd = data_day[j1:j2].to_numpy()
            day_flag = mad(np.ma.masked_array(data=dd, mask=np.isnan(dd)), z=z, deriv=deriv)
            ff.iloc[fill_start:fill_end, cols.index(hcol)] += np.where(day_flag[fill_start-j1-1:fill_end-j1-1], 2, 0)

            nn = data_night[j1:j2]
            night_flag = mad(np.ma.masked_array(data=nn, mask=np.isnan(nn)), z=z, deriv=deriv)
            ff.iloc[fill_start:fill_end, cols.index(hcol)] += np.where(night_flag[fill_start-j1-1:fill_end-j1-1], 2, 0)

        if plot:
            fig = plt.figure(1)
            sub = fig.add_subplot(111)
            valid = ff[hcol]==0
            l1 = sub.plot(data[valid], 'ob')
            l3 = sub.plot(data[ff[hcol]==2], 'or')
            plt.title(hcol)
            pp.savefig(fig)
            plt.close(fig)

    # Finish

    if plot:
        pp.close()

    if fisnumpy:
        if fistrans:
            return ff.to_numpy().T
        else:
            return ff.to_numpy()
    else:
        return ff


if __name__ == '__main__':
    import doctest
    doctest.testmod()
