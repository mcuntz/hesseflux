#!/usr/bin/env python
"""
ustarfilter : Filter Eddy Covariance data with friction velocity
              following Papale et al. (Biogeosciences, 2006).

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
* Using numpy docstring format, May 2020, Matthias Cuntz
* No bootstrap by default, Jul 2020, Matthias Cuntz
* Optionally return thresholds and flags for each season, Jul 2020, Matthias Cuntz
* Bugfix if no threshold found, and for multi-year flags, Jul 2020, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz, Arndt Piayda, Tino Rau

The following functions are provided

.. autosummary::
   ustarfilter
"""
from __future__ import division, absolute_import, print_function
import numpy as np
import pandas as pd


__all__ = ['ustarfilter']


def ustarfilter(dfin, flag=None, isday=None, date=None, timeformat='%Y-%m-%d %H:%M:%S', colhead=None,
                ustarmin=0.01, nboot=1, undef=-9999, plot=False, seasonout=False,
                nmon=3, ntaclasses=7, corrcheck=0.5, nustarclasses=20, plateaucrit=0.95, swthr=10.):
    """
    Flag Eddy Covariance data using a threshold of friction velocity (u*)
    below which u* correlates with a reduction in CO2 flux. The algorithm
    follows the method presented in Papale et al. (Biogeosciences, 2006).

    Parameters
    ----------
    dfin : pandas.Dataframe or numpy.array
        time series of CO2 fluxes and friction velocity as well as air temperature.

        `dfin` can be a pandas.Dataframe with the columns
        'FC' or 'NEE' (or starting with 'FC\_' or 'NEE\_') for observed CO2 flux [umol(CO2) m-2 s-1]
        'USTAR' (or starting with 'USTAR') for friction velocity [m s-1]
        'TA'    (or starting with 'TA\_') for air temperature [deg C]
        The index is taken as date variable.

        `dfin` can also me a numpy array with the same columns. In this case
        `colhead`, `date`, and possibly `dateformat` must be given.
    flag : pandas.Dataframe or numpy.array, optional
        flag Dataframe or array has the same shape as dfin. Non-zero values in
        `flag` will be treated as missing values in `dfin`.

        `flag` must follow the same rules as `dfin` if pandas.Dataframe.

        If `flag` is numpy array, `df.columns.values` will be used as column heads
        and the index of `dfin` will be copied to `flag`.
    isday : array_like of bool, optional
        True when it is day, False when night. Must have the same length as dfin.shape[0].

        If `isday` is not given, `dfin` must have a column with head 'SW_IN' or
        starting with 'SW_IN'. `isday` will then be `dfin['SW_IN'] > swthr`.
    date : array_like of string, optional
        1D-array_like of calendar dates in format given in `timeformat`.

        `date` must be given if `dfin` is numpy array.
    timeformat : str, optional
        Format of dates in `date`, if given (default: '%Y-%m-%d %H:%M:%S').
        See strftime documentation of Python's datetime module:
        https://docs.python.org/3/library/datetime.html#strftime-and-strptime-behavior
    colhead : array_like of str, optional
        column names if `dfin` is numpy array. See `dfin` for mandatory column names.
    ustarmin : float, optional
        minimum ustar threshold (default: 0.01)

        Papale et al. (Biogeosciences, 2006) take 0.1 for forest ant 0.01 otherwise.
    nboot : int, optional
        number of boot straps for estimate of confidence interval of u* threshold (default: 1)
    undef : float, optional
        values having `undef` value are treated as missing values in `dfin` (default: -9999)
    plot : bool, optional
        True: data and u* thresholds are plotted into ustarfilter.pdf (default: False)
    seasonout : bool, optional
        True: return u* thresholds for each season (default: False)
    nmon : int, optional
        Number of month to combine for a season (default: 3).
    ntaclasses : int, optional
        Number of temperature classes per `nmon` months (default: 7).
    corrcheck : float, optional
        Skip temperature class if absolute of correlation coefficient between air temperature
        and ustar is greater equal `corrcheck` (default: 0.5).
    nustarclasses : int, optional
        Number of u* classes per temperature class (default: 20).
    plateaucrit : float, optional
        The u* threshold is the smallest u* class that has an average CO2 flux,
        which is higher than `plateaucrit` times the mean CO2 flux of all u* classes
        above this class (default: 0.95).
    swthr : float, optional
        Threshold to determine daytime from incoming shortwave radiation if `isday` not given (default: 10).

    Returns
    -------
    tuple of numpy arrays, pandas.Dataframe or numpy array
        numpy array with 5, 50 and 95 percentile of u* thresholds,
        flags: 0 everywhere except set to 2 where u* < u*-threshold.

        Either maximum threshold of all seasons or thresholds for each season, i.e.
        threshold array is `array(3,nseason)` if `seasonout` and `array(3)` otherwise.

    Notes
    -----
    Works ONLY for a data set of at least one full year.

    History
    -------
    Written,    Matthias Cuntz & Tino Rau, 2008
    Maintained, Arndt Piayda,   Aug 2014
    Modified,   Matthias Cuntz, Apr 2020 - input can be pandas Dataframe or numpy array(s)
                Matthias Cuntz, May 2020 - numpy docstring format
                Matthias Cuntz, Jul 2020 - default nboot=1
                Matthias Cuntz, Jul 2020 - seasonout
                Matthias Cuntz, Jul 2020 - bugfix if no threshold found flag_p -> flag_b
                                         - bugfix in multi-year flags
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
        isday = isday.to_numpy() # otherwise time indexes could not match anymore after shifting times
    isday[isday == undef] = np.nan
    ff[np.isnan(isday)]   = 1

    # data and flags
    fc_id = ''
    for cc in df.columns:
        if cc.startswith('FC_') or (cc == 'FC') or cc.startswith('NEE_') or (cc == 'NEE'):
            fc_id = cc
            break
    ustar_id = ''
    for cc in df.columns:
        if cc.startswith('USTAR_') or (cc == 'USTAR'):
            ustar_id = cc
            break
    ta_id = ''
    for cc in df.columns:
        if cc.startswith('TA_') or (cc == 'TA'):
            ta_id = cc
            break
    assert fc_id,    'Carbon net flux with name FC or NEE or starting with FC_ or NEE_ must be in input.'
    assert ustar_id, 'Friction velocity u* with name USTAR or starting with USTAR_ must be in input.'
    assert ta_id,    'Air temperature with name TA or starting with TA_ must be in input.'

    # Move time steps by half interval if end times
    ustar_in = df[ustar_id].copy(deep=True) # save for output with original DatetimeIndex
    times = df.index.minute.sort_values().unique()
    indx  = df.index.values
    nindx = indx.copy()
    if len(times) == 1:
        if 0 in times:
            # hourly time steps, shift by 30 min
            dt = pd.Timedelta('30 minute')
            if df.index[0].hour == 0:
                nindx += dt # beginning of time step
            else:
                nindx -= dt # end of time step
    elif len(times) == 2:
        if 0 in times:
            # half-hourly time steps, shift by 15 min
            dt = pd.Timedelta('15 minute')
            if df.index[0].minute == 0:
                nindx += dt # beginning of time step
            else:
                nindx -= dt # end of time step
    else:
        raise ValueError('Could not analyse time steps. Known time steps are hourly or half-hourly.')
    df['DateTime'] = nindx
    df = df.set_index('DateTime', drop=True)
    ff = ff.set_index(df.index)

    # Check time span
    yrmax  = df.index.max().year
    yrmin  = df.index.min().year
    nyears = yrmax - yrmin + 1
    ndays  = (df.index.max() - df.index.min()).days + 1
    assert ndays//nyears > 360, 'Full years must be given.'

    # calculate thresholds
    nperiod = 12//nmon  # number of nmon periods per year
    if seasonout:
        bustars = np.ones((nboot,nyears,nperiod)) * undef
    else:
        bustars = np.ones((nboot,nyears)) * undef
    for y in range(nyears):
        yy = yrmin + y
        iiyr    = df.index.year == yy
        iidf    = df.index[iiyr]
        df_f    = df.loc[iidf]
        ff_f    = ff.loc[iidf]
        isday_f = isday[iiyr]

        # bootstrap per year
        nstepsyear = len(df_f)
        for b in range(nboot):
            if b == 0:
                df_b    = df_f
                ff_b    = ff_f
                isday_b = isday_f
            else:
                iiboot  = np.random.randint(0, nstepsyear, size=nstepsyear)
                iidf    = df_f.index[iiboot]
                df_b    = df_f.loc[iidf]
                ff_b    = ff_f.loc[iidf]
                isday_b = isday_f[iiboot]

            # periods / seasons
            pustars = np.ones(nperiod) * undef
            for p in range(nperiod):
                flag_p   = ( (~isday_b) &
                            (ff_b[fc_id] == 0) & (ff_b[ustar_id] == 0) & (ff_b[ta_id] == 0) &
                            (df_b.index.month > p*nmon) & (df_b.index.month <= (p+1)*nmon) )
                fc_p    = df_b.loc[flag_p, fc_id]
                ustar_p = df_b.loc[flag_p, ustar_id]
                ta_p    = df_b.loc[flag_p, ta_id]

                # temperature classes
                custars = []
                ta_q     = np.quantile(ta_p, np.arange(ntaclasses+1,dtype=np.float)/np.float(ntaclasses))
                ta_q[0] -= 0.1 # 1st include min
                for t in range(ntaclasses):
                    iita    = (ta_p > ta_q[t]) & (ta_p <= ta_q[t+1])
                    fc_t    = fc_p[iita]
                    ustar_t = ustar_p[iita]
                    ta_t    = ta_p[iita]

                    # discard temperature class if correlation is strong
                    r_abs = np.abs(np.corrcoef(ta_t, ustar_t)[0,1])
                    if r_abs >= corrcheck:
                        continue

                    # ustar classes
                    ustar_q     = np.quantile(ustar_t, np.arange(nustarclasses+1,dtype=np.float)/np.float(nustarclasses))
                    ustar_q[0] -= 0.01 # 1st include min
                    for u in range(nustarclasses-1):
                        iiustar = (ustar_t > ustar_q[u]) & (ustar_t <= ustar_q[u+1])
                        fc_u    = fc_t[iiustar]
                        fc_a    = fc_t[ustar_t > ustar_q[u+1]]

                        if abs(fc_u.mean()) >= abs(plateaucrit*fc_a.mean()):
                            custars.append(ustar_q[u+1])
                            break

                # median of thresholds of all temperature classes = threshold of period
                if len(custars) > 0:
                    pustars[p] = np.median(custars)
                elif seasonout:
                    # Set threshold to 90% of data per season if no threshold found
                    pustars[p] = np.quantile(ustar_p, 0.9)
            # Take maximum of periods if any thresholds found,
            # otherwise set threshold to 90% of data
            ii = np.where(pustars != undef)[0]
            if ii.size > 0:
                if seasonout:
                    bustars[b,y,:] = pustars
                else:
                    bustars[b,y]   = pustars[ii].max()
            else:
                if seasonout:
                    raise ValueError('Should not be here.')
                flag_b       = ( (~isday_b) &
                                 (ff_b[fc_id] == 0) & (ff_b[ustar_id] == 0) & (ff_b[ta_id] == 0) )
                bustars[b,y] = np.quantile(df_b.loc[flag_b, ustar_id], 0.9)

    # set minimum ustar threshold
    bustars = np.maximum(bustars, ustarmin)

    # report 5, 50 and 95 percentile
    oustars = np.quantile(bustars, (0.05, 0.5, 0.95), axis=0)

    # flag out with original DatetimeIndex
    off    = ustar_in.astype(int)
    off[:] = 0
    ii     = np.zeros(len(off), dtype=np.bool)
    if seasonout:
        for y in range(nyears):
            yy = yrmin + y
            for p in range(nperiod):
                iiyr = ( (df.index.year == yy) &           # df DatetimeIndex
                         (df.index.month > p*nmon) &
                         (df.index.month <= (p+1)*nmon) )
                ii[iiyr] = ustar_in[iiyr] < oustars[1,y,p]
    else:
        for y in range(nyears):
            yy    = yrmin + y
            iiyr  = df.index.year == yy # df DatetimeIndex
            ii[iiyr] = ustar_in[iiyr] < oustars[1,y]
    off[ii] = 2 # original DatetimeIndex

    if plot:
        import matplotlib.pyplot as plt
        import matplotlib.backends.backend_pdf as pdf
        pd.plotting.register_matplotlib_converters()

        pp = pdf.PdfPages('ustarfilter.pdf')
        if seasonout:
            for y in range(nyears):
                yy = yrmin + y

                iiyr    = df.index.year == yy
                iidf    = df.index[iiyr]
                df_f    = df.loc[iidf]
                ff_f    = ff.loc[iidf]
                isday_f = isday[iiyr]

                for p in range(nperiod):
                    flag_p = ( (~isday_f) &
                               (ff_f[fc_id] == 0) & (ff_f[ustar_id] == 0) & (ff_f[ta_id] == 0) &
                               (df_f.index.month > p*nmon) & (df_f.index.month <= (p+1)*nmon) )
                    fc_p    = df_f.loc[flag_p, fc_id]
                    ustar_p = df_f.loc[flag_p, ustar_id]

                    fig  = plt.figure(1)
                    sub  = fig.add_subplot(111)
                    sub.plot(ustar_p, fc_p, 'bo')
                    sub.axvline(x=oustars[1,y,p], linewidth=0.75, color='r')
                    plt.ylabel('F CO2')
                    plt.xlabel('u_star')
                    plt.title('u_star thresh for season {:d} of year {:d}: {:5.3f}'.format(
                        p, yy, oustars[1,y,p]))

                    pp.savefig(fig)
                    plt.close(fig)
        else:
            for y in range(nyears):
                yy = yrmin + y

                iiyr    = (df.index.year == yy) & (~isday)
                fc_y    = df.loc[iiyr, fc_id]
                ffc_y   = ff.loc[iiyr, fc_id]
                ustar_y = df.loc[iiyr, ustar_id]
                ffu_y   = ff.loc[iiyr, ustar_id]
                # off_y   = off[iiyr & (isday == False)]

                fig  = plt.figure(1)
                sub  = fig.add_subplot(111)
                flag_p = (ffu_y == 0) & (ffc_y == 0)
                sub.plot(ustar_y[flag_p], fc_y[flag_p], 'bo')
                sub.axvline(x=oustars[1,y], linewidth=0.75, color='r')
                plt.ylabel('F CO2')
                plt.xlabel('u_star')
                plt.title('u_star thresh of year ${:d}: {:5.3f}'.format(yy, oustars[1,y]))

                pp.savefig(fig)
                plt.close(fig)
        pp.close()

    # out
    oustars = oustars.squeeze()

    if isnumpy:
        return oustars, off.to_numpy()
    else:
        return oustars, off


if __name__ == '__main__':
    import doctest
    doctest.testmod()
