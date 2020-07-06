==========
User Guide
==========

``hesseflux`` collects functions used for processing Eddy covariance data of the `ICOS
<https://www.icos-cp.eu/>`_ ecosystem site `FR-Hes
<https://www.icos-france.fr/en/static3/the-network>`_.

The package uses several functions of the JAMS Python package

   https://github.com/mcuntz/jams_python

The JAMS package and hesseflux are synchronised irregularly.

``hesseflux`` includes a Python port of Logtools, the Logger Tools Software of
Olaf Kolle, MPI-BGC Jena, (c) 2012.

The post-processing functionality for Eddy flux data is similar to the R-package `REddyProc
<https://cran.r-project.org/web/packages/REddyProc/index.html>`_ and includes basically the steps
described in `Papale et al. (Biogeosciences, 2006) <https://doi.org/10.5194/bg-3-571-2006>`_ plus
some extensions such as the daytime method of flux partitioning (`Lasslop et al., Global Change
Biology 2010 <https://doi.org/10.1111/j.1365-2486.2009.02041.x>`_) and the estimation of
uncertainties on the fluxes as in `Lasslop et al. (Biogeosci, 2008)
<http://doi.org/10.5194/bg-5-1311-2008>`_.

Only the post-processing steps are described here. We are happy to discuss any processing or
post-processing directly. Contact us at mc (at) macu (dot) de.


europe-fluxdata.eu file format
==============================

The first processing steps at the ICOS ecosystem site FR-Hes (not shown) brings the data in a
format that can be submitted to the database `europe-fluxdata.eu <http://www.europe-fluxdata.eu>`_.
The database predates ICOS and is somewhat a precursor of the ICOS data processing.

The file format of europe-fluxdata.eu is hence very similar to the ICOS format. The only known
difference to us is the unit of atmospheric pressure, which is in hPa in `europe-fluxdata.eu
<http://www.europe-fluxdata.eu>`_ and in kPa in `ICOS ecosystems <http://www.icos-etc.eu>`_. The
file format has notably one header line with variable names. There are no units in the file.
``hesseflux`` provides a little helper script :any:`europe-fluxdata_units.py` in the `bin`
directory that adds a second header line with units. The script can be run on the output as:

.. code-block:: bash

    python europe-fluxdata_units.py output_file_of_postproc_europe-fluxdata.csv


Post-processing Eddy covariance data
====================================

The script :any:`postproc_europe-fluxdata.py` in the `example` directory provides a template for
post-processing data that is in the `europe-fluxdata.eu <http://www.europe-fluxdata.eu>`_ file
format. It basically makes all steps described in `Papale et al. (Biogeosciences, 2006)
<https://doi.org/10.5194/bg-3-571-2006>`_. The script is governed by a configuration file in
Python's standard :mod:`configparser` format. The example configuration file
:any:`hesseflux_example.cfg` in the `example` directory is highly commented and should be (almost)
self-explanatory. The script is called like:

.. code-block:: bash
		
    python postproc_europe-fluxdata.py hesseflux_example.cfg

This script should be taken as a template for one's own post-processing but includes most standard
post-processing steps.

Here we describe the main parts of the post-processing script.


Reading the configuration file
------------------------------

The script :any:`postproc_europe-fluxdata.py` starts by reading the configuration file :any:`hesseflux_example.cfg`:

.. code-block:: python

    import sys
    import configparser

    # Read config file
    if len(sys.argv) <= 1:
        raise IOError('Input configuration file must be given.')
    configfile = sys.argv[1]
    config = configparser.ConfigParser(interpolation=None)
    config.read(configfile)

It then analyses the configuration options. The first section in the configuration file are the
options controlling which steps shall be performed by the script. The section in the
:any:`hesseflux_example.cfg` looks like:

.. code-block:: python

    [POSTSWITCH]
    # spike detection (Papale et al., Biogeoci 2006)
    # bool
    outlier   = True
    # ustar filtering (Papale et al., Biogeoci 2006)
    # bool
    ustar     = True
    # flux partitioning (Reichstein et al., GCB 2005; Lasslop et al., GCB 2010)
    # bool
    partition = True
    # gap filling (Reichstein et al., GCB 2005)
    # bool
    fill      = True
    # error estimate of Eddy fluxes (Lasslop et al., Biogeosci 2008)
    # bool
    fluxerr   = False

And the code in :any:`postproc_europe-fluxdata.py` is:

.. code-block:: python

    # program switches
    outlier   = config['POSTSWITCH'].getboolean('outlier',   True)
    ustar     = config['POSTSWITCH'].getboolean('ustar',     True)
    partition = config['POSTSWITCH'].getboolean('partition', True)
    fill      = config['POSTSWITCH'].getboolean('fill',      True)
    fluxerr   = config['POSTSWITCH'].getboolean('fluxerr',   True)

All options are boolean and set to `True` by default if they are not given in the configuration
file. All post-processing steps except uncertainty estimation of flux data would be performed in
the given example.


Read the data
-------------

The script would then read in the data. The section in the configuration file is:

.. code-block:: python

    [POSTIO]
    # can be comma separated list or single file
    # str
    inputfile  = FR-Hes_europe-fluxdata_2016.txt
    # see strftime documentation of Python's datetime module
    # https://docs.python.org/3/library/datetime.html#strftime-and-strptime-behavior
    # str
    timeformat  = %Y%m%d%H%M
    # Delimter to use with pandas.read_csv. If None, Python’s builtin sniffer tool is used (slow)
    # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_csv.html
    # str
    sep         = ,
    # Line numbers to skip (0-indexed) or number of lines to skip (int) at the start of the file.
    # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_csv.html
    # list-like, int
    skiprows    = None
    # values being NaN and undef will be ignored
    # float
    undef       = -9999.
    # threshold of shortwave radiation for determining day/night. day is SW_IN > swthr
    # float
    swthr       = 10.

The analysis of the options in :any:`postproc_europe-fluxdata.py` is:

.. code-block:: python

    # input file
    eufluxfile  = config['POSTIO'].get('inputfile',  '')
    timeformat  = config['POSTIO'].get('timeformat', '%Y%m%d%H%M')
    sep         = config['POSTIO'].get('sep',        ',')
    skiprows    = config['POSTIO'].get('skiprows',   None)
    undef       = config['POSTIO'].getfloat('undef', -9999.)
    swthr       = config['POSTIO'].getfloat('swthr', 10.)

Note that strings are given without quotes in the configuration file.

`eufluxfile` can be a single filename or a comma-separated list of filenames. If it is missing or
empty, the script will try to open a GUI, where one can choose input files. The data will be
appended if several input files are given.

The (first) input file is read as:

.. code-block:: python

    import pandas as pd

    parser = lambda date: pd.datetime.strptime(date, timeformat)
    infile = eufluxfile[0]
    df     = pd.read_csv(infile, sep, skiprows=skiprows, parse_dates=[0], date_parser=parser, index_col=0, header=0)

:mod:`pandas` will use the first column as index (`index_col=0`), assuming that these are dates
(`parse_dates=[0]`) in the format `timeformat`, where columns are separated by `sep`. The defaults
follow the europe-fluxdata.eu format but similar formats may be used, and script and/or
configuration file can be adapted easily. Only variable names have to follow europe-fluxdata.eu,
ICOS or Ameriflux format at the moment. If the input file has a second header line with units, one
can skip it giving `skiprows=[1]` (not `skiprows=1`).

All input files are supposed to be in the same format, if `eufluxfile` is a comma-separated list of
filenames, and they will be read with the same command above. The :mod:`pandas` dataframes (`df`)
will simply be appended.


The flag dataframe
------------------

All Not-a-Number (NaN) values will be set to `undef` and will be ignored in the following.

This happens via a second dataframe (`dff`), having the same columns and index as the input
dataframe `df`, representing quality flags. All cells that have a value other than `0` in the flag
dataframe `dff` will be ignored in the dataframe `df`. This means all cells of `df` with `undef`
will be set to `2` in `dff` immediately:

.. code-block:: python

    # NaN -> undef
    df.fillna(undef, inplace=True)

    # Flag
    dff              = df.copy(deep=True).astype(int)
    dff[:]           = 0
    dff[df == undef] = 2


Day / night
-----------

Most post-processing routines differentiate between daytime and nighttime data. `Papale et al.
(Biogeosciences, 2006) <https://doi.org/10.5194/bg-3-571-2006>`_ use a threshold of 20 W m\
:sup:`-2` of global radiation to distinguish between day and night. `REddyProc
<https://cran.r-project.org/web/packages/REddyProc/index.html>`_ uses incoming shortwave radiation
greater than 10 W m\ :sup:`2` as daytime. The shortwave raditan threshold `swthr` (same name as in
ReddyProc) can be used to define the appropriate threshold. The default is 10 W m\ :sup:`2`. The
column `SW_IN_1_1_1` has to exist in the input data.

.. code-block:: python

    # day / night
    isday = df['SW_IN_1_1_1'] > swthr


Data check
----------

:any:`postproc_europe-fluxdata.py` checks the units of air temperature (i.e. the first column
starting with `TA_`).

.. code-block:: python

    # Check Ta in Kelvin
    hta = ['TA_']
    hout = _findfirststart(hta, df.columns)
    if df[hout[0]].max() < 100.:
        tkelvin = 273.15
    else:
        tkelvin = 0.
    df.loc[dff[hout[0]]==0, hout[0]] += tkelvin

:func:`_findfirststart(starts, names)` is a helper function that finds the first occurrence in
`names` that starts with the string `starts`. This is used for the moment until it was implemented
in `hesseflux`` that the user can give individual variable names.

The script calculates air vapour pressure deficit `VPD_PI_1_1_1` from air temperature and relative
humidity (i.e. the first column starting with `RH_`) if not given in input data using the function
:func:`esat` of ``hesseflux`` for saturation vapour pressure:

.. code-block:: python

    import numpy as np
    import hesseflux as hf

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
        vpd = (1.-rh) * hf.esat(tk)
        vpd_id = 'VPD_PI_1_1_1'
        df[vpd_id] = vpd
        df[vpd_id].where((df[ta_id]!=undef) | (df[rh_id]!=undef), other=undef, inplace=True)
        dff[vpd_id] = np.where((dff[ta_id] + dff[rh_id]) > 0, 2, 0)
        df.loc[dff[vpd_id]==0, vpd_id] /= 100. # hPa as in europe-fluxdata.eu

It further checks assures that VPD is in Pa for further calculations.

.. code-block:: python

   # Check VPD in Pa
    hvpd = ['VPD']
    hout = _findfirststart(hvpd, df.columns)
    if df[hout[0]].max() < 10.: # kPa
        vpdpa = 1000.
    elif df[hout[0]].max() < 100.: # hPa
        vpdpa = 100.
    else:
        vpdpa = 1.
    df.loc[dff[hout[0]]==0, hout[0]] *= vpdpa

And finally determines the time intervals of the input data `dtsec` (s) and the number of time
steps per day `ntday`.

.. code-block:: python

    # time stepping
    dsec  = (df.index[1] - df.index[0]).seconds
    ntday = np.rint(86400/dsec).astype(np.int)


Spike / outlier flagging
------------------------

If `outlier=True` is set in the configuration file, spikes will be detected with the method given
in `Papale et al. (Biogeosciences, 2006) <https://doi.org/10.5194/bg-3-571-2006>`_. A median
absolute deviation (MAD) filter will be used on the second derivatives of the time series in
two-week chunks. The section in :any:`hesseflux_example.cfg` looks like:

.. code-block:: python

    [POSTMAD]
    # spike / outlier detection, see help(hesseflux.madspikes)
    # scan window in days for spike detection
    # int
    nscan = 15
    # fill window in days for spike detection
    # int
    nfill = 1
    # spike will be set for values above z absolute deviations
    # float
    z     = 7.
    # 0: mad on raw values; 1, 2: mad on first or second derivatives
    # int
    deriv = 2

`nfill` is the number of days that are treated at once. `nfill=1` means that the time series will
be stepped through day by day. `nscan` are the days to be considered to calculate the mean absolute
deviations. `nscan=15` means that 7 days before the fill day, the fill day itself and 7 days after
the fill day will be used for the robust statistic. However, only spikes detected within the inner
`nfill` days will be flagged in the `nscan` days. Spikes will be detected if they deviate more than
`z` mean absolute deviations from the median. `deriv=2` applies the MAD filter to the second
derivatives. A spike has normally a strong curvature and hence a large second derivative. `deriv=1`
is currently not implemented. `deriv=0` applies the filter to the raw time series. This might be
useful to find outliers in smooth time series such as soil moisture. `deriv=0` is also used on the
20 Hz Eddy raw data in the quality and uncertainty strategy of `Mauder et al. (Agric Forest
Meteo, 2013) <http://doi.org/10.1016/j.agrformet.2012.09.006>`_.

The default values, if option are not given in the configuration file, are `nscan=15`, `nfill=1`,
`z=7`, and `deriv=2`.

:any:`postproc_europe-fluxdata.py` calls the spike detection like this:

.. code-block:: python

    houtlier = ['H_', 'LE', 'FC',       # assume *_PI variables after raw variables, e.g. LE before LE_PI
                'H_PI', 'LE_PI', 'NEE'] # if available
    hout = _findfirststart(houtlier, df.columns)
    sflag = hf.madspikes(df[hout], flag=dff[hout], isday=isday, undef=undef,
                         nscan=nscan*ntday, nfill=nfill*ntday, z=z, deriv=deriv,
                         plot=False)
    for ii, hh in enumerate(hout):
        dff.loc[sflag[hh]==2, hh] = 3

The function :func:`madspikes` returns flag columns for the input variables where spiked data is
flagged as 2. The scripts sets the corresponding columns in the flag dataframe `dff` to 3 (3 just
to keep track where the flag was set).


u* filtering
------------

If `ustar=True` is set in the configuration file, a u*-filter will be applied following `Papale et
al. (Biogeosciences, 2006) <https://doi.org/10.5194/bg-3-571-2006>`_.

The section in :any:`hesseflux_example.cfg` looks like:

.. code-block:: python

    [POSTUSTAR]
    # ustar filtering, see help(hesseflux.ustarfilter)
    # min ustar value. Papale et al. (Biogeosci 2006): 0.1 forest, 0.01 else
    # float
    ustarmin    = 0.1
    # number of boostraps for determining uncertainty of ustar threshold. 1 = no bootstrap
    # int
    nboot       = 1
    # significant difference between ustar class and mean of classes above
    # float
    plateaucrit = 0.95

A minimum threhold `ustarmin` is defined under which data is flagged by default. `Papale et al.
(Biogeosciences, 2006) <https://doi.org/10.5194/bg-3-571-2006>`_ suggest 0.1 for forests and 0.01
for other land cover types. :any:`postproc_europe-fluxdata.py` sets 0.01 as its default value.
Uncertainty of the u*-threshold is calculated via bootstraping in Papale et al. `nboot` gives the
number of bootstrepping for an the uncertainty estimate of the u*-threshold. The algorithm divides
the input data (per season) in 7 temperature classes and in 20 u*-classes within each temperature
class. It then determines the threshold as the average u* of the u*-class where the average CO2
flux is less than `plateaucrit` times the average of all CO2 fluxes with u* greater than the
u*-class. `Papale et al. (Biogeosciences, 2006) <https://doi.org/10.5194/bg-3-571-2006>`_ took
`plateaucrit=0.99`, while `REddyProc
<https://cran.r-project.org/web/packages/REddyProc/index.html>`_ takes `plateaucrit=0.95`, which
:any:`postproc_europe-fluxdata.py` also takes as its default.

The u*-filtering is then performed as:

.. code-block:: python

    # The algorithm uses NEE to determine u*-threshold
    hfilt = ['NEE', 'USTAR', 'TA_']
    hout = _findfirststart(hfilt, df.columns)
    # the algorithm does not work with carbon uptake at night -> flag it temporarily
    ffsave = dff[hout[0]].to_numpy()
    iic    = np.where((~isday) & (df[hout[0]] < 0.))[0]
    dff.iloc[iic, list(df.columns).index(hout[0])] = 4
    ustars, flag = hf.ustarfilter(df[hout], flag=dff[hout], isday=isday, undef=undef,
                                  ustarmin=ustarmin, nboot=nboot, plateaucrit=plateaucrit,
                                  plot=True)
    dff[hout[0]] = ffsave # set NEE<0 @ night flags back
    # The threshold is then applied to all eddy fluxes
    hustar = ['H_', 'LE', 'FC',       # assume *_PI variables after raw variables, e.g. LE before LE_PI
              'H_PI', 'LE_PI', 'NEE'] # if available
    hout = _findfirststart(hustar, df.columns)
    for ii, hh in enumerate(hout):
        dff.loc[flag==2, hh] = 5

The function :func:`ustarfilter` returns the ustar 5, 50 and 95 percentile of the bootstrapped
u*-thresholds and a flag columns, which is 0 except where u* is smaller than the median
u*-threshold. The scripts sets the columns of the Eddy fluxes in the flag dataframe `dff` to 5 (5
just to keep track where the flag was set).

One might not want to do u*-filtering, but use for example Integral Turbulence Characteristics
(ITC) that were calculated, for example, with EddyPro\ :sup:`(R)`. These should be set right at the
start after reading the input data into the dataframe `df` and producing the flag dataframe `dff`
like:

.. code-block:: python

    dff.loc[df['FC_SSITC_TEST_1_1_1']>0, 'FC_1_1_1'] = 2


Partitioning of Net Ecosystem Exchange
--------------------------------------

If `partition=True` is set in the configuration file, two estimates of Gross Primary Productivity
(GPP) and Ecosystem Respiration (RECO) are calculated: firstly with the method of `Reichstein et
al. (Glob Change Biolo, 2005) <http://doi.org/10.1111/j.1365-2486.2005.001002.x>`_ using nighttime
data only, and secondly with the method of `Lasslop et al. (Glob Change Biolo, 2010)
<http://doi.org/10.1111/j.1365-2486.2009.02041.x>`_ using a light-response curve on 'daytime' data.
The configuration :any:`hesseflux_example.cfg` gives only one option in this section:

.. code-block:: python

    [POSTPARTITION]
    # partitioning, see help(hesseflux.nee2gpp)
    # if True, set GPP=0 at night
    # bool
    nogppnight = False

Many people find it unaesthetic that the 'daytime' method gives negative GPP at night. We esteem
this the correct behaviour, reflecting the uncertainty in the gross flux estimates. However, one
can set `nogppnight=True` to set GPP=0 at night and RECO=NEE in this case, the latter having then
all variability of the net fluxes.

The partitioning is calculated as:

.. code-block:: python

    hpart = ['NEE', 'SW_IN', 'TA_', 'VPD']
    hout  = _findfirststart(hpart, df.columns)
    # nighttime method
    dfpartn = hf.nee2gpp(df[hout], flag=dff[hout], isday=isday,
                         undef=undef, method='reichstein', nogppnight=nogppnight)
    suff = hout[0][3:-1]
    dfpartn.rename(columns=lambda c: c+suff+'1', inplace=True)

    # daytime method
    dfpartd = hf.nee2gpp(df[hout], flag=dff[hout], isday=isday,
                         undef=undef, method='lasslop', nogppnight=nogppnight)
    dfpartd.rename(columns=lambda c: c+suff+'2', inplace=True)

    # add new columns to dataframe
    df = pd.concat([df, dfpartn, dfpartd],  axis=1)

    # take flags from NEE for the new columns
    for dn in ['1', '2']:
        for gg in ['GPP', 'RECO']:
            dff[gg+suff+dn] = dff[hout[0]]


Gap-filling / Imputation
------------------------

Marginal Distribution Sampling (MDS) of `Reichstein et al. (Glob Change Biolo, 2005)
<http://doi.org/10.1111/j.1365-2486.2005.001002.x>`_ is implemented as imputation or called
gap-filling algorithm. The algorithm looks for similar conditions in the vicinity of a missing data
point, if option `fill=True`. The configuration file is:

.. code-block:: python

    [POSTGAP]
    # gap-filling with MDS, see help(hesseflux.gapfill)
    # max deviation of SW_IN
    # float
    sw_dev  = 50.
    # max deviation of TA
    # float
    ta_dev  = 2.5
    # max deviation of VPD
    # float
    vpd_dev = 5.0
    # avoid extrapolation in gaps longer than longgap days
    longgap = 60

If a flux data point is missing, times with incoming shortwave radiation in the range of `sw_dev`
around the actual shortwave radiation will be looked for, as well as air temperatures within
`ta_dev` and air vapour pressure deficit within `vpd_dev`. The function does not fill long gaps
longer than `longgap` days. A good summary is given in Fig. A1 of `Reichstein et al. (Glob Change
Biolo, 2005) <http://doi.org/10.1111/j.1365-2486.2005.001002.x>`_.

The script invokes MDS as:

.. code-block:: python

    hfill = ['H_', 'LE', 'FC',       # assume *_PI variables after raw variables, e.g. LE before LE_PI
             'H_PI', 'LE_PI', 'NEE', # if available
             'GPP_1_1_1', 'RECO_1_1_1', 'GPP_1_1_2', 'RECO_1_1_2',
             'GPP_PI_1_1_1', 'RECO_PI_1_1_1', 'GPP_PI_1_1_2', 'RECO_PI_1_1_2',
             'SW_IN', 'TA_', 'VPD']
    hout  = _findfirststart(hfill, df.columns)
    df_f, dff_f = hf.gapfill(df[hout], flag=dff[hout],
                             sw_dev=sw_dev, ta_dev=ta_dev, vpd_dev=vpd_dev,
                             longgap=longgap, undef=undef, err=False, verbose=1)
    # remove meteorology columns
    hdrop = ['SW_IN', 'TA_', 'VPD']
    hout = _findfirststart(hdrop, df.columns)
    df_f.drop(columns=hout,  inplace=True)
    dff_f.drop(columns=hout, inplace=True)
    # we add _f to columns names of filled variables
    df_f.rename(columns=lambda c: '_'.join(c.split('_')[:-3]+['f']+c.split('_')[-3:]),  inplace=True)
    dff_f.rename(columns=lambda c: '_'.join(c.split('_')[:-3]+['f']+c.split('_')[-3:]), inplace=True)
    df  = pd.concat([df,  df_f],  axis=1)
    dff = pd.concat([dff, dff_f], axis=1)

The function :func:`gapfill` returns the filled columns `df_f` as well as flag
columns `dff_f` indicating fill quality. Fill quality A-C of Reichstein et al.
are translated to quality flags 1-3.


Uncertainty estimates of flux data
----------------------------------

`Lasslop et al. (Biogeosci, 2008) <http://doi.org/10.5194/bg-5-1311-2008>`_ presented an algorithm
to estimate uncertainties of Eddy covariance fluxes using Marginal Distribution Sampling (MDS). The
gap-filling function :func:`gapfill` can be used for uncertainty estimation giving the keyword
`err=True`. The same thresholds as for gap-filling are used.

The script :any:`postproc_europe-fluxdata.py` uses the function :func:`gapfill` to calculate flux
uncertainties like:

.. code-block:: python

    hfill = ['H_', 'LE', 'FC',       # assume *_PI variables after raw variables, e.g. LE before LE_PI
             'H_PI', 'LE_PI', 'NEE', # if available
             'H_f', 'LE_f', 'FC_f',
             'H_PI_f', 'LE_PI_f', 'NEE_f', 'NEE_PI_f',
             'GPP_1_1_1', 'RECO_1_1_1', 'GPP_1_1_2', 'RECO_1_1_2',
             'GPP_f_1_1_1', 'RECO_f_1_1_1', 'GPP_f_1_1_2', 'RECO_f_1_1_2',
             'GPP_PI_1_1_1', 'RECO_PI_1_1_1', 'GPP_PI_1_1_2', 'RECO_PI_1_1_2',
             'GPP_PI_f_1_1_1', 'RECO_PI_f_1_1_1', 'GPP_PI_f_1_1_2', 'RECO_PI_f_1_1_2',
             'SW_IN', 'TA_', 'VPD']
    hout = _findfirststart(hfill, df.columns)
    df_f = hf.gapfill(df[hout], flag=dff[hout],
                      sw_dev=sw_dev, ta_dev=ta_dev, vpd_dev=vpd_dev,
                      longgap=longgap, undef=undef, err=True, verbose=1)
    # remove meteorology columns
    hdrop = ['SW_IN', 'TA_', 'VPD']
    hout = _findfirststart(hdrop, df.columns)
    df_f.drop(columns=hout, inplace=True)
    colin = list(df_f.columns)
    # we add _err to columns names of uncertainty estimates, such as: NEE_PI_err_1_1_1
    df_f.rename(columns=lambda c: '_'.join(c.split('_')[:-3]+['err']+c.split('_')[-3:]),  inplace=True)
    colout = list(df_f.columns)
    df = pd.concat([df, df_f], axis=1)
    # take flags of non-error columns
    for cc in range(len(colin)):
        dff[colout[cc]] = dff[colin[cc]]

We recommend, however, to calculate flux uncertainties with the Eddy covariance raw data as
described in `Mauder et al. (Agric Forest Meteo, 2013)
<http://doi.org/10.1016/j.agrformet.2012.09.006>`_. This is, for example, implemented in the
processing software EddyPro\ :sup:`(R)`.


Writing the output file
-----------------------

The dataframe is written to the output file with :mod:`pandas` :func:`pandas.Dataframe.to_csv`:

.. code-block:: python

    df.to_csv(outputfile, sep=sep, na_rep=str(undef), index=True, date_format=timeformat)

using the same `sep` and `timeformat` as the input.

The configuration for output is:

.. code-block:: python

    [POSTIO]
    # if empty, write will ask for output path using the name of this config file with the suffix .csv
    outputfile  = FR-Hes_europe-fluxdata_2016-post.txt
    # if True, set variable to undef where flagged in output
    # bool
    outundef    = True
    # if True, add flag columns prepended with flag_ for each variable
    # bool
    outflagcols = False

If `outputfile` is missing or empty, the script will try to open a GUI, where one can choose an
output directory and the filename will then be name of the configuration file with the suffix
'.csv'.

If `outundef=True` then all values in `df` with a flag value in `dff` greater than zero will be set
to `undef`. The script can also add flag columns, prefixed with `flag_`, for each column in `df`,
if `outflagcols=True`. The script will always output the columns with the flags for fill quality,
if gap-filling was performed: option `fill=True`.

The code before :func:`pandas.Dataframe.to_csv` is then:

.. code-block:: python

    if outundef:
        for cc in df.columns:
            if cc.split('_')[-4] != 'f': # exclude gap-filled columns
                df[cc].where(dff[cc] == 0, other=undef, inplace=True)
    if outflagcols:
        dff.rename(columns=lambda c: 'flag_'+c, inplace=True)
        df = pd.concat([df, dff], axis=1)
    else:
        occ = []
        for cc in df.columns:
            if cc.split('_')[-4] == 'f': occ.append(cc)
        dff1 = dff[occ].copy(deep=True)
        dff1.rename(columns=lambda c: 'flag_'+c, inplace=True)
        df = pd.concat([df, dff1], axis=1)

That's all Folks!
