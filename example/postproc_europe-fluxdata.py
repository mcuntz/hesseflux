#!/usr/bin/env python
"""
This scripts runs post-processing steps for Eddy covariance data coming
in one file in the format of europe-fluxdata.eu. This format is very similar
to the ICOS format (the only known difference is the unit of pressure,
which is hPa in europe-fluxdata.eu and kPa in ICOS).

The script covers the following steps:
- spike / outlier detection with mean absolute deviation filter
  after Papale et al. (Biogeosci, 2006)
- ustar filtering after Papale et al. (Biogeosci, 2006)
- carbon flux partitioning with the nighttime method
  of Reichstein et al. (Global Change Biolo, 2005) and
  the daytime method of Lasslop et al. (Global Change Biolo, 2010)
- gap filling with marginal distribution sampling (MDS)
  of Reichstein et al. (Global Change Biolo, 2005)
- flux error estimates using MDS after Lasslop et al. (Biogeosci, 2008)

The script is controlled by a config file in Python's standard configparser
format. The config file includes all possible parameters of used routines.
Default parameter values follow the package REddyProc where appropriate. See
comments in config file for details.

The script currently flags on input all NaN values and given *undefined*
values. Variables should be set to *undefined* in case of other existing flags
before calling the script. Otherwise it should be easy to set the appropriate
flags in the pandas DataFrame dff for the flags after its creation around line
160.

The output file can either have all flagged variables set to *undefined*
and/or can include flag columns for each variable (see config file).

Note, ustar filtering needs at least one full year.

Examples
--------
python postproc_europe-fluxdata.py hesseflux_example.cfg

History
-------
Written, Matthias Cuntz, April 2020
"""
from __future__ import division, absolute_import, print_function
import time as ptime
import sys
import configparser
import os.path
import datetime as dt
import numpy as np
import pandas as pd
import hesseflux as hf


#
# Find first elements in names that begin with elements of starts
def _findfirststart(starts, names):
    """
    Find first elements in names that begin with elements of starts

    Example
    -------
    >>> hout = _findfirststart(['TA', 'LE'],
                               ['TIMESTAMP', 'TAU_1_1_1',
                                'H_1_1_1', 'LE_1_1_1',
                                'TA_1_1_1', 'TA_1_2_1', 'LE_PI_1_1_1'])
    >>> print(hout)
    TAU_1_1_1 LE_1_1_1
    """
    hout = []
    for hh in starts:
        for cc in names:
            if cc.startswith(hh):
                hout.append(cc)
                break
    return hout


if __name__ == '__main__':
    t1 = ptime.time()

    # ToDo
    #   - Allow more flexibility in column names

    # Read config file
    if len(sys.argv) <= 1:
        raise IOError('Input configuration file must be given.')
    configfile = sys.argv[1]
    config = configparser.ConfigParser(interpolation=None)
    config.read(configfile)
    # file path
    outdir    = config['GENERAL'].get('outdir', ".")
    # program switches
    outlier   = config['POSTSWITCH'].getboolean('outlier',   True)
    ustar     = config['POSTSWITCH'].getboolean('ustar',     True)
    partition = config['POSTSWITCH'].getboolean('partition', True)
    fill      = config['POSTSWITCH'].getboolean('fill',      True)
    fluxerr   = config['POSTSWITCH'].getboolean('fluxerr',   True)
    # input file
    eufluxfile  = config['POSTIO'].get('inputfile',  '')
    timeformat  = config['POSTIO'].get('timeformat', '%Y%m%d%H%M')
    sep         = config['POSTIO'].get('sep',        ',')
    skiprows    = config['POSTIO'].get('skiprows',   '')
    undef       = config['POSTIO'].getfloat('undef', -9999.)
    swthr       = config['POSTIO'].getfloat('swthr', 10.)
    outputfile  = config['POSTIO'].get('outputfile'  '')
    outundef    = config['POSTSWITCH'].getboolean('outundef',    True)
    outflagcols = config['POSTSWITCH'].getboolean('outflagcols', False)
    # mad
    nscan = config['POSTMAD'].getint('nscan', 15)
    nfill = config['POSTMAD'].getint('nfill',  1)
    z     = config['POSTMAD'].getfloat('z',    7)
    deriv = config['POSTMAD'].getint('deriv',  2)
    # ustar
    ustarmin       = config['POSTUSTAR'].getfloat('ustarmin',    0.1)
    nboot          = config['POSTUSTAR'].getint('nboot',         1)
    plateaucrit    = config['POSTUSTAR'].getfloat('plateaucrit', 0.95)
    seasonout      = config['POSTUSTAR'].getboolean('seasonout', False)
    applyustarflag = config['POSTUSTAR'].getboolean('applyflag', False)
    # gap-filling
    sw_dev  = config['POSTGAP'].getfloat('sw_dev',  50.)
    ta_dev  = config['POSTGAP'].getfloat('ta_dev',  2.5)
    vpd_dev = config['POSTGAP'].getfloat('vpd_dev', 5.0)
    longgap = config['POSTGAP'].getint('longgap',   60)
    # partitioning
    nogppnight = config['POSTPARTITION'].getboolean('nogppnight', False)

    # ----------------------------------------------------------------
    # Check call

    # Assert iterable
    if ',' in eufluxfile:
        eufluxfile = eufluxfile.split(',')
        eufluxfile = [ ee.strip() for ee in eufluxfile ]
    else:
        if eufluxfile:
            eufluxfile = [eufluxfile]
        else:
            try:
                eufluxfile = hf.files_from_gui(
                    initialdir='.', title='europe-fluxdata.eu file(s)')
            except:
                raise IOError("GUI for europe-fluxdata.eu file(s) failed.")

    if skiprows == 'None':
        skiprows = ''
    if skiprows:
        import json  # to analyse int or list, tuple not working
        skiprows = json.loads(skiprows.replace('(', '[').replace(')', ']'))

    # ----------------------------------------------------------------
    # Read input files into Panda data frame and check variable availability

    print('Read data ', eufluxfile)
    t01 = ptime.time()
    # TIMESTAMP,TAU_1_1_1,H_1_1_1,LE_1_1_1,FC_1_1_1,...
    # 201901010030,0.0941,-11.0765,-9999.0000,-9999.0000,...
    # use lambda because of global var timeformat
    parser = lambda date: dt.datetime.strptime(date, timeformat)

    infile = eufluxfile[0]
    df = pd.read_csv(infile, sep, skiprows=skiprows, parse_dates=[0],
                     date_parser=parser, index_col=0, header=0)

    if len(eufluxfile) > 1:
        for infile in eufluxfile[1:]:
            df1 = pd.read_csv(infile, sep, skiprows=skiprows, parse_dates=[0],
                              date_parser=parser, index_col=0, header=0)
            df  = df.append(df1, sort=False)
    df.fillna(undef, inplace=True)
    # df.replace(-9999., np.nan, inplace=True)

    # Flag
    dff              = df.copy(deep=True).astype(int)
    dff[:]           = 0
    dff[df == undef] = 2
    # dff[df.isna()]   = 2

    # day / night
    isday = df['SW_IN_1_1_1'] > swthr

    # Check Ta in Kelvin
    hta = ['TA_']
    hout = _findfirststart(hta, df.columns)
    if df[hout[0]].max() < 100.:
        tkelvin = 273.15
    else:
        tkelvin = 0.
    # add tkelvin only where not flagged
    df.loc[dff[hout[0]] == 0, hout[0]] += tkelvin

    # add vpd if not given
    hvpd = ['VPD']
    hout = _findfirststart(hvpd, df.columns)
    if len(hout) == 0:
        hvpd = ['TA_', 'RH_']
        hout = _findfirststart(hvpd, df.columns)
        if len(hout) != 2:
            raise ValueError('Cannot calculate VPD.')
        ta_id = hout[0]
        rh_id = hout[1]
        if df[ta_id].max() < 100.:
            tk = df[ta_id] + 273.15
        else:
            tk = df[ta_id]
        if df[rh_id].max() > 10.:
            rh = df[rh_id] / 100.
        else:
            rh = df[rh_id]
        vpd = (1. - rh) * hf.esat(tk)
        vpd_id = 'VPD_PI_1_1_1'
        df[vpd_id] = vpd
        df[vpd_id].where((df[ta_id] != undef) | (df[rh_id] != undef),
                         other=undef, inplace=True)
        dff[vpd_id] = np.where((dff[ta_id] + dff[rh_id]) > 0, 2, 0)
        df.loc[dff[vpd_id] == 0, vpd_id] /= 100.

    # Check VPD in Pa
    hvpd = ['VPD']
    hout = _findfirststart(hvpd, df.columns)
    if df[hout[0]].max() < 10.:     # kPa
        vpdpa = 1000.
    elif df[hout[0]].max() < 100.:  # hPa
        vpdpa = 100.
    else:
        vpdpa = 1.                  # Pa
    df.loc[dff[hout[0]] == 0, hout[0]] *= vpdpa

    # time stepping
    dsec  = (df.index[1] - df.index[0]).seconds
    ntday = np.rint(86400 / dsec).astype(np.int)

    t02   = ptime.time()
    strin = ( '[m]: {:.1f}'.format((t02 - t01) / 60.)
              if (t02 - t01) > 60.
              else '[s]: {:d}'.format(int(t02 - t01)) )
    print('   in ', strin)

    # ----------------------------------------------------------------
    # Outlier detection

    if outlier:
        print('Spike detection')
        t11 = ptime.time()
        # assume *_PI variables after raw variables, e.g. LE before LE_PI,
        # if available
        houtlier = ['H_', 'LE', 'FC',
                    'H_PI', 'LE_PI', 'NEE']
        # houtlier = ['FC', 'NEE']
        hout = _findfirststart(houtlier, df.columns)
        print('  Using', hout)
        # ToDo
        #   - only one call to mad for all variables
        sflag = hf.madspikes(df[hout], flag=dff[hout], isday=isday,
                             undef=undef, nscan=nscan * ntday,
                             nfill=nfill * ntday, z=z, deriv=deriv, plot=False)
        for ii, hh in enumerate(hout):
            dff.loc[sflag[hh] == 2, hh] = 3

        t12   = ptime.time()
        strin = ( '[m]: {:.1f}'.format((t12 - t11) / 60.)
                  if (t12 - t11) > 60.
                  else '[s]: {:d}'.format(int(t12 - t11)) )
        print('   in ', strin)

    # ----------------------------------------------------------------
    # u* filtering

    if ustar:
        print('u* filtering')
        t21 = ptime.time()
        hfilt = ['NEE', 'USTAR', 'TA_']
        hout  = _findfirststart(hfilt, df.columns)
        if len(hout) == 2:  # take FC if NEE not in input file
            hfilt = ['FC', 'USTAR', 'TA_']
            hout  = _findfirststart(hfilt, df.columns)
        assert len(hout) == 3, 'Could not find CO2 flux (NEE or FC), USTAR or TA in input file.'
        hout = _findfirststart(hfilt, df.columns)
        print('  Using', hout)
        ffsave = dff[hout[0]].to_numpy()
        iic    = np.where((~isday) & (df[hout[0]] < 0.))[0]
        dff.iloc[iic, list(df.columns).index(hout[0])] = 4
        ustars, flag = hf.ustarfilter(df[hout], flag=dff[hout],
                                      isday=isday, undef=undef,
                                      ustarmin=ustarmin, nboot=nboot,
                                      plateaucrit=plateaucrit,
                                      seasonout=seasonout,
                                      plot=True)
        dff[hout[0]] = ffsave
        df  = df.assign(USTAR_TEST_1_1_1=flag)
        dff = dff.assign(USTAR_TEST_1_1_1=np.zeros(df.shape[0], dtype=np.int))
        if applyustarflag:
            # assume *_PI variables after raw variables, e.g. LE before LE_PI
            # if available
            hustar = ['H_', 'LE', 'FC',
                      'H_PI', 'LE_PI', 'NEE']
            hout = _findfirststart(hustar, df.columns)
            print('  Using', hout)
            for ii, hh in enumerate(hout):
                dff.loc[flag == 2, hh] = 5

        t22   = ptime.time()
        strin = ( '[m]: {:.1f}'.format((t22 - t21) / 60.)
                  if (t22 - t21) > 60.
                  else '[s]: {:d}'.format(int(t22 - t21)) )
        print('   in ', strin)

    # ----------------------------------------------------------------
    # Flux partitioning

    if partition:
        print('Flux partitioning')
        t41 = ptime.time()
        hpart = ['NEE', 'SW_IN', 'TA_', 'VPD']
        hout  = _findfirststart(hpart, df.columns)
        if len(hout) == 3:  # take FC if NEE not in input file
            hpart = ['FC', 'SW_IN', 'TA_', 'VPD']
            hout  = _findfirststart(hpart, df.columns)
        print('  Using', hout)
        astr = 'Could not find CO2 flux (NEE or FC), SW_IN, TA, or VPD in input file.'
        assert len(hout) == 4, astr
        # nighttime method
        print('  Nighttime partitioning')
        dfpartn = hf.nee2gpp(df[hout], flag=dff[hout], isday=isday,
                             undef=undef, method='reichstein',
                             nogppnight=nogppnight)
        if hout[0].startswith('NEE'):
            suff = hout[0][3:-1]
        else:
            suff = hout[0][2:-1]
        dfpartn.rename(columns=lambda c: c + suff + '1', inplace=True)
        # daytime method
        print('  Daytime partitioning')
        dfpartd = hf.nee2gpp(df[hout], flag=dff[hout], isday=isday,
                             undef=undef, method='lasslop',
                             nogppnight=nogppnight)
        dfpartd.rename(columns=lambda c: c + suff + '2', inplace=True)
        df = pd.concat([df, dfpartn, dfpartd],  axis=1)
        # take flags from NEE
        for dn in ['1', '2']:
            for gg in ['GPP', 'RECO']:
                dff[gg + suff + dn] = dff[hout[0]]
        # flag if partitioning was not possible
        for dn in ['1', '2']:
            for gg in ['GPP', 'RECO']:
                dff[df[gg + suff + dn] == undef] = 2

        t42   = ptime.time()
        strin = ( '[m]: {:.1f}'.format((t42 - t41) / 60.)
                  if (t42 - t41) > 60.
                  else '[s]: {:d}'.format(int(t42 - t41)) )
        print('   in ', strin)

    # ----------------------------------------------------------------
    # Gap-filling

    if fill:
        print('Gap-filling')
        t31 = ptime.time()
        hfill = ['SW_IN', 'TA_', 'VPD']
        hout  = _findfirststart(hfill, df.columns)
        assert len(hout) == 3, 'Could not find SW_IN, TA or VPD in input file.'
        # assume *_PI variables after raw variables, e.g. LE before LE_PI
        # if available
        hfill = ['H_', 'LE', 'FC',
                 'H_PI', 'LE_PI', 'NEE',
                 'GPP_1_1_1', 'RECO_1_1_1',
                 'GPP_1_1_2', 'RECO_1_1_2',
                 'GPP_PI_1_1_1', 'RECO_PI_1_1_1',
                 'GPP_PI_1_1_2', 'RECO_PI_1_1_2',
                 'SW_IN', 'TA_', 'VPD']
        # hfill = ['NEE', 'SW_IN', 'TA_', 'VPD']
        hout  = _findfirststart(hfill, df.columns)
        print('  Using', hout)
        df_f, dff_f = hf.gapfill(df[hout], flag=dff[hout],
                                 sw_dev=sw_dev, ta_dev=ta_dev, vpd_dev=vpd_dev,
                                 longgap=longgap, undef=undef, err=False,
                                 verbose=1)
        hdrop = ['SW_IN', 'TA_', 'VPD']
        hout = _findfirststart(hdrop, df.columns)
        df_f.drop(columns=hout,  inplace=True)
        dff_f.drop(columns=hout, inplace=True)

        def _add_f(c):
            return '_'.join(c.split('_')[:-3] + ['f'] + c.split('_')[-3:])
        df_f.rename(columns=_add_f,  inplace=True)
        dff_f.rename(columns=_add_f, inplace=True)
        df  = pd.concat([df,  df_f],  axis=1)
        dff = pd.concat([dff, dff_f], axis=1)

        t32   = ptime.time()
        strin = ( '[m]: {:.1f}'.format((t32 - t31) / 60.)
                  if (t32 - t31) > 60.
                  else '[s]: {:d}'.format(int(t32 - t31)) )
        print('   in ', strin)

    # ----------------------------------------------------------------
    # Error estimate

    if fluxerr:
        print('Flux error estimates')
        t51 = ptime.time()
        hfill = ['SW_IN', 'TA_', 'VPD']
        hout  = _findfirststart(hfill, df.columns)
        assert len(hout) == 3, 'Could not find SW_IN, TA or VPD in input file.'
        # assume *_PI variables after raw variables, e.g. LE before LE_PI
        # if available
        hfill = ['H_', 'LE', 'FC',
                 'H_PI', 'LE_PI', 'NEE',
                 'H_f', 'LE_f', 'FC_f',
                 'H_PI_f', 'LE_PI_f', 'NEE_f', 'NEE_PI_f',
                 'GPP_1_1_1', 'RECO_1_1_1',
                 'GPP_1_1_2', 'RECO_1_1_2',
                 'GPP_f_1_1_1', 'RECO_f_1_1_1',
                 'GPP_f_1_1_2', 'RECO_f_1_1_2',
                 'GPP_PI_1_1_1', 'RECO_PI_1_1_1',
                 'GPP_PI_1_1_2', 'RECO_PI_1_1_2',
                 'GPP_PI_f_1_1_1', 'RECO_PI_f_1_1_1',
                 'GPP_PI_f_1_1_2', 'RECO_PI_f_1_1_2',
                 'SW_IN', 'TA_', 'VPD']
        # hfill = ['NEE', 'GPP', 'SW_IN', 'TA_', 'VPD']
        hout = _findfirststart(hfill, df.columns)
        print('  Using', hout)
        df_f = hf.gapfill(df[hout], flag=dff[hout],
                          sw_dev=sw_dev, ta_dev=ta_dev, vpd_dev=vpd_dev,
                          longgap=longgap, undef=undef, err=True, verbose=1)
        hdrop = ['SW_IN', 'TA_', 'VPD']
        hout = _findfirststart(hdrop, df.columns)
        df_f.drop(columns=hout, inplace=True)
        colin = list(df_f.columns)
        # names such as: NEE_PI_err_1_1_1
        df_f.rename(columns=_add_f,  inplace=True)
        colout = list(df_f.columns)
        df = pd.concat([df, df_f], axis=1)
        # take flags of non-error columns
        for cc in range(len(colin)):
            dff[colout[cc]] = dff[colin[cc]]

        t52   = ptime.time()
        strin = ( '[m]: {:.1f}'.format((t52 - t51) / 60.)
                  if (t52 - t51) > 60.
                  else '[s]: {:d}'.format(int(t52 - t51)) )
        print('   in ', strin)

    # ----------------------------------------------------------------
    # Output

    if not outputfile:
        try:
            outputdir = hf.directory_from_gui(initialdir='.',
                                              title='Output directory')
        except:
            raise IOError("GUI for output directory failed.")
        outputfile = configfile[:configfile.rfind('.')]
        outputfile = outdir + '/' + os.path.basename(outputfile + '.csv')
    print('Write output ', outputfile)
    t61 = ptime.time()
    # Back to original units
    hta = ['TA_']
    hout = _findfirststart(hta, df.columns)
    df.loc[dff[hout[0]] == 0, hout[0]] -= tkelvin
    hvpd = ['VPD']
    hout = _findfirststart(hvpd, df.columns)
    df.loc[dff[hout[0]] == 0, hout[0]] /= vpdpa
    if outundef:
        print('   Set flags to undef.')
        for cc in df.columns:
            if cc.split('_')[-4] != 'f':  # exclude gap-filled columns
                df[cc].where(dff[cc] == 0, other=undef, inplace=True)
    if outflagcols:
        print('   Add flag columns.')

        def _add_flag(c):
            return 'flag_' + c
        dff.rename(columns=_add_flag, inplace=True)
        # no flag columns for flags
        dcol = []
        for hh in dff.columns:
            if '_TEST_' in hh:
                dcol.append(hh)
        if dcol:
            dff.drop(columns=dcol, inplace=True)
        df = pd.concat([df, dff], axis=1)
    else:
        print('   Add flag columns for gap-filled variables.')
        occ = []
        for cc in df.columns:
            if cc.split('_')[-4] == 'f':
                occ.append(cc)
        dff1 = dff[occ].copy(deep=True)
        dff1.rename(columns=lambda c: 'flag_' + c, inplace=True)
        df = pd.concat([df, dff1], axis=1)
    print('   Write.')
    df.to_csv(outputfile, sep=sep, na_rep=str(undef), index=True,
              date_format=timeformat)

    t62   = ptime.time()
    strin = ( '[m]: {:.1f}'.format((t62 - t61) / 60.)
              if (t62 - t61) > 60.
              else '[s]: {:d}'.format(int(t62 - t61)) )
    print('   in ', strin)

    # ----------------------------------------------------------------
    # Finish

    t2    = ptime.time()
    strin = ( '[m]: {:.1f}'.format((t2 - t1) / 60.)
              if (t2 - t1) > 60.
              else '[s]: {:d}'.format(int(t2 - t1)) )
    print('Time elapsed', strin)
