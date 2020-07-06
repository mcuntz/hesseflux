#!/usr/bin/env python
"""
date2dec : Converts calendar dates into decimal dates.

This module was written by Arndt Piayda and then enhanced and
maintained by Matthias Cuntz while at Department of Computational
Hydrosystems, Helmholtz Centre for Environmental Research - UFZ,
Leipzig, Germany, and continued by Matthias Cuntz while at Institut
National de Recherche pour l'Agriculture, l'Alimentation et
l'Environnement (INRAE), Nancy, France.

Copyright (c) 2010-2020 Matthias Cuntz - mc (at) macu (dot) de
Released under the MIT License; see LICENSE file for details.

* Written Jun 2010 by Arndt Piayda
* Input can be scalar, list, array, or mix of it, Feb 2012, Matthias Cuntz
* Changed checks, added calendars decimal and decimal360, Feb 2012, Matthias Cuntz
* Changed units of proleptic_gregorian calendar from days since 0001-01-01 00:00:00 to days since 0001-01-00 00:00:00, Dec 2010, Matthias Cuntz
* Deal with Excel leap year error, Feb 2013, Matthias Cuntz
* Ported to Python 3, Feb 2013, Matthias Cuntz
* ascii/eng without time defaults to 00:00:00, Jul 2013, Matthias Cuntz
* Excel starts at 1 not at 0 on 01 January 1900 or 1904, Oct 2013, Matthias Cuntz
* Bug: 01.01.0001 was substracted if Julian calendar even with units given, Oct 2013, Matthias Cuntz
* Removed remnant of time treatment before time check in eng keyword, Nov 2013, Matthias Cuntz
* Adapted to new netCDF4/netcdftime (>= v1.0) and datetime (>= Python v2.7.9), Jun 2015, Matthias Cuntz
* Call date2num with list instead of single netCDF4.datetime objects, Oct 2015, Matthias Cuntz
* mo for months always integer, Oct 2016, Matthias Cuntz
* 00, 01, etc. for integers not accepted by Python 3, removed from examples and code, Nov 2016, Matthias Cuntz
* Using numpy docstring format, May 2020, Matthias Cuntz
* Succeed eng by en keyword as in ascii2ascii and dec2date, Jul 2020, Matthias Cuntz
* proleptic_gregorian instead of gregorian calendar for Excel dates, Jul 2020, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz, Arndt Piayda

The following functions are provided

.. autosummary::
   date2dec
"""
from __future__ import division, absolute_import, print_function
import numpy as np


__all__ = ['date2dec']


def date2dec(calendar = 'standard', units=None,
             excelerr = True, yr=1,
             mo=1, dy=1, hr=0, mi=0, sc=0,
             ascii=None, en=None, eng=None):
    """
    Convert scalar and array_like with calendar dates into decimal
    dates. Supported calendar formats are standard, gregorian, julian,
    proleptic_gregorian, excel1900, excel1904, 365_day, noleap, 366_day,
    all_leap, 360_day, decimal, or decimal360.

    Input is year, month day, hour, minute, second or a combination of them.
    Input as date string is possible.

    Output is decimal date with day as unit.

    Parameters
    ----------
    yr : array_like, optional
        years (default: 1)
    mo : array_like, optional
        month (default: 1)
    dy : array_like, optional
        days (default: 1)
    hr : array_like, optional
        hours (default: 0)
    mi : array_like, optional
        minutes (default: 0)
    sc : array_like, optional
        seconds (default: 0)
    ascii : array_like, optional
        strings of the format 'dd.mm.yyyy hh:mm:ss'.
        Missing hour, minutes and/or seconds are set
        to their default values (0).

        `ascii` overwrites all other keywords.

        `ascii` and `eng` are mutually exclusive.
    en : array_like, optional
        strings of the format 'yyyy-mm-dd hh:mm:ss'.
        Missing hour, minutes and/or seconds are set
        to their default values (0).

        `en` overwrites all other keywords.

        `en` and `ascii` are mutually exclusive.
    eng : array_like, optional
        Same as en: obsolete.
    calendar : str, optional
        Calendar of output dates (default: 'standard').

        Possible values are:

        'standard', 'gregorian' = julian calendar from
        01.01.-4712 12:00:00 (BC) until 05.03.1583 00:00:00 and
        gregorian calendar from 15.03.1583 00:00:00 until now.
        Missing 10 days do not exsist.

        'julian' = julian calendar from 01.01.-4712 12:00:00 (BC)
         until now.

        'proleptic_gregorian' = gregorian calendar from
        01.01.0001 00:00:00 until now.

        'excel1900' = Excel dates with origin at
        01.01.1900 00:00:00.

        'excel1904' = Excel 1904 (Lotus) format.
        Same as excel1904 but with origin at
        01.01.1904 00:00:00.

        '365_day', 'noleap' = 365 days format,
        i.e. common years only (no leap years)
        with origin at 01.01.0001 00:00:00.

        '366_day', 'all_leap' = 366 days format,
        i.e. leap years only (no common years)
        with origin at 01.01.0001 00:00:00.

        '360_day' = 360 days format,
        i.e. years with only 360 days (30 days per month)
        with origin at 01.01.0001 00:00:00.

        'decimal' = decimal year instead of decimal days.

        'decimal360' = decimal year with a year of 360 days, i.e. 12 month with 30 days each.
    units : str, optional
        User set units of output dates. Must be a
        string in the format 'days since yyyy-mm-dd hh:mm:ss'.
        Default values are set automatically depending on `calendar`.
    excelerr : bool, optional
       In Excel, the year 1900 is normally considered a leap year,
       which it was not. By default, this error is taken into account
       if `calendar=='excel1900'` (default: True).

       1900 is not considered a leap year if `excelerr==False`.

    Returns
    -------
    array_like
       array_like with decimal dates. The type of output will be the same as the input type.

    Notes
    -----
    Most versions of `datetime` do not support negative years,
    i.e. Julian days < 1721423.5 = 01.01.0001 00:00:00.

    There is an issue in `netcdftime` version < 0.9.5 in proleptic_gregorian for dates before year 301:
    dec2date(date2dec(ascii='01.01.0300 00:00:00', calendar='proleptic_gregorian'), calendar='proleptic_gregorian')
    [300, 1, 2, 0, 0, 0]
    dec2date(date2dec(ascii='01.01.0301 00:00:00', calendar='proleptic_gregorian'), calendar='proleptic_gregorian')
    [301, 1, 1, 0, 0, 0]

    List input is only supported up to 2 dimensions.

    Requires `netcdftime.py` from module `netcdftime` available at:
    http://netcdf4-python.googlecode.com

    Examples
    --------
    # Some implementations of datetime do not support negative years
    >>> import datetime
    >>> if datetime.MINYEAR > 0:
    ...     print('The minimum year in your datetime implementation is ', datetime.MINYEAR)
    ...     print('i.e. it does not support negative years (BC).')
    The minimum year in your datetime implementation is  1
    i.e. it does not support negative years (BC).

    >>> if datetime.MINYEAR > 0:
    ...     year   = np.array([2000, 1810, 1630, 1510, 1271, 619, 2, 1])
    ... else:
    ...     year   = np.array([2000, 1810, 1630, 1510, 1271, 619, -1579, -4712])
    >>> month  = np.array([1, 4, 7, 9, 3, 8, 8, 1])
    >>> day    = np.array([5, 24, 15, 20, 18, 27, 23, 1])
    >>> hour   = np.array([12, 16, 10, 14, 19, 11, 20, 12])
    >>> minute = np.array([30, 15, 20, 35, 41, 8, 3, 0])
    >>> second = np.array([15, 10, 40, 50, 34, 37, 41, 0])
    >>> decimal = date2dec(calendar = 'standard', yr=year, mo=month, dy=day, hr=hour, mi=minute, sc=second)
    >>> nn = year.size
    >>> print('{:.14e} {:.14e} {:.14e} {:.14e}'.format(*decimal[:nn//2]))
    2.45154902100694e+06 2.38226217719907e+06 2.31660093101852e+06 2.27284810821759e+06
    >>> print('{:.14e} {:.14e}'.format(*decimal[nn//2:nn-2]))
    2.18536732053241e+06 1.94738596431713e+06
    >>> decimal = date2dec(calendar='standard', yr=year, mo=6, dy=15, hr=12, mi=minute, sc=second)
    >>> print('{:.14e} {:.14e} {:.14e} {:.14e}'.format(*decimal[:nn//2]))
    2.45171102100694e+06 2.38231401053241e+06 2.31657101435185e+06 2.27275102488426e+06
    >>> print('{:.14e} {:.14e}'.format(*decimal[nn//2:nn-2]))
    2.18545602886574e+06 1.94731300598380e+06

    # ascii input
    >>> if datetime.MINYEAR > 0:
    ...     a = np.array(['05.01.2000 12:30:15', '24.04.1810 16:15:10', '15.07.1630 10:20:40', '20.09.1510 14:35:50',
    ...                   '18.03.1271 19:41:34', '27.08. 619 11:08:37', '23.08.0002 20:03:41', '01.01.0001 12:00:00'])
    ... else:
    ...     a = np.array(['05.01.2000 12:30:15', '24.04.1810 16:15:10', '15.07.1630 10:20:40',  '20.09.1510 14:35:50',
    ...                   '18.03.1271 19:41:34', '27.08. 619 11:08:37', '23.08.-1579 20:03:41', '01.01.-4712 12:00:00'])
    >>> decimal = date2dec(calendar='standard', ascii=a)
    >>> nn = a.size
    >>> print('{:.14e} {:.14e} {:.14e} {:.14e}'.format(*decimal[:nn//2]))
    2.45154902100694e+06 2.38226217719907e+06 2.31660093101852e+06 2.27284810821759e+06
    >>> print('{:.14e} {:.14e}'.format(*decimal[nn//2:nn-2]))
    2.18536732053241e+06 1.94738596431713e+06

    # calendar = 'julian'
    >>> decimal = date2dec(calendar='julian', ascii=a)
    >>> print('{:.14e} {:.14e} {:.14e} {:.14e}'.format(*decimal[:nn//2]))
    2.45156202100694e+06 2.38227417719907e+06 2.31661093101852e+06 2.27284810821759e+06
    >>> print('{:.14e} {:.14e}'.format(*decimal[nn//2:nn-2]))
    2.18536732053241e+06 1.94738596431713e+06

    # calendar = 'proleptic_gregorian'
    >>> decimal = date2dec(calendar='proleptic_gregorian', ascii=a)
    >>> print('{:.7f} {:.7f} {:.7f} {:.7f}'.format(*decimal[:nn//2]))
    730123.5210069 660836.6771991 595175.4310185 551412.6082176
    >>> print('{:.7f} {:.7f}'.format(*decimal[nn//2:nn-2]))
    463934.8205324 225957.4643171

    # calendar = 'excel1900' WITH excelerr=True -> 1900 considered as leap year
    >>> d = np.array(['05.01.2000 12:30:15', '27.05.1950 16:25:10', '13.08.1910 10:40:55',
    ...               '01.03.1900 00:00:00', '29.02.1900 00:00:00', '28.02.1900 00:00:00',
    ...               '01.01.1900 00:00:00'])
    >>> decimal = date2dec(calendar='excel1900', ascii=d)
    >>> nn = d.size
    >>> print('{:.7f} {:.7f} {:.7f}'.format(*decimal[:nn//2]))
    36530.5210069 18410.6841435 3878.4450810
    >>> print('{:.1f} {:.1f} {:.1f} {:.1f}'.format(*decimal[nn//2:]))
    61.0 60.0 59.0 1.0

    # calendar = 'excel1900' WITH excelerr = False -> 1900 is NO leap year
    >>> decimal = date2dec(calendar='excel1900', ascii=d, excelerr=False)
    >>> print('{:.7f} {:.7f} {:.7f}'.format(*decimal[:nn//2]))
    36529.5210069 18409.6841435 3877.4450810
    >>> print('{:.1f} {:.1f} {:.1f} {:.1f}'.format(*decimal[nn//2:]))
    60.0 60.0 59.0 1.0

    # calendar = 'excel1904'
    >>> decimal = date2dec(calendar='excel1904', ascii=d[:nn//2])
    >>> print('{:.7f} {:.7f} {:.7f}'.format(*decimal[:nn//2]))
    35069.5210069 16949.6841435 2417.4450810

    # calendar = '365_day'
    >>> g = np.array(['18.08.1972 12:30:15', '25.10.0986 12:30:15', '28.11.0493 22:20:40', '01.01.0001 00:00:00'])
    >>> decimal = date2dec(calendar='365_day', ascii=g)
    >>> nn = g.size
    >>> print('{:.7f} {:.7f} {:.7f} {:.7f}'.format(*decimal[:nn]))
    719644.5210069 359822.5210069 179911.9310185 0.0000000

    # calendar = '366_day'
    >>> decimal = date2dec(calendar='366_day', ascii=g)
    >>> print('{:.7f} {:.7f} {:.7f} {:.7f}'.format(*decimal[:nn]))
    721616.5210069 360808.5210069 180404.9310185 0.0000000

    # 360_day does not work with netcdftime.py version equal or below 0.9.2
    # calendar = '360_day'
    >>> decimal = date2dec(calendar='360_day', ascii=g)
    >>> print('{:.7f} {:.7f} {:.7f} {:.7f}'.format(*decimal[:nn]))
    709787.5210069 354894.5210069 177447.9310185 0.0000000

    >>> print('{:.7f}'.format(date2dec(yr=1992, mo=1, dy=26, hr=2, mi=0, sc=0, calendar='decimal')))
    1992.0685337
    >>> print('{:.7f}'.format(date2dec(ascii='26.01.1992 02:00', calendar='decimal360')))
    1992.0696759
    >>> print('{:.7f} {:.7f}'.format(*date2dec(ascii=['26.01.1992 02:00','26.01.1992 02:00'], calendar='decimal360')))
    1992.0696759 1992.0696759
    >>> print('{:.7f} {:.7f}'.format(*date2dec(yr=[1992,1992], mo=1, dy=26, hr=2, mi=0, sc=0, calendar='decimal360')))
    1992.0696759 1992.0696759
    >>> print('{:.7f} {:.7f}'.format(*date2dec(yr=np.array([1992,1992]), mo=1, dy=26, hr=2, mi=0, sc=0, calendar='decimal360')))
    1992.0696759 1992.0696759
    >>> decimal = date2dec(ascii=[['26.01.1992 02:00','26.01.1992 02:00'],
    ...                           ['26.01.1992 02:00','26.01.1992 02:00'],
    ...                           ['26.01.1992 02:00','26.01.1992 02:00']],
    ...                    calendar='decimal360')
    >>> print('{:.7f} {:.7f}'.format(*decimal[0]))
    1992.0696759 1992.0696759
    >>> print('{:.7f} {:.7f}'.format(*decimal[2]))
    1992.0696759 1992.0696759
    >>> print((date2dec(ascii='01.03.2003 00:00:00') - date2dec(ascii='01.03.2003')) == 0.)
    True

    # en
    >>> decimal = date2dec(en='1992-01-26 02:00', calendar='decimal360')
    >>> print('{:.7f}'.format(decimal))
    1992.0696759
    >>> decimal = date2dec(eng='1992-01-26 02:00', calendar='decimal360')
    >>> print('{:.7f}'.format(decimal))
    1992.0696759

    History
    -------
    Written  Arndt Piayda, Jun 2010
    Modified Matthias Cuntz, Feb 2012 - All input can be scalar, list or array, also a mix
                                        - Changed checks for easier extension
                                        - decimal, decimal360
             Matthias Cuntz, Dec 2012 - change unit of proleptic_gregorian
                                        from 'days since 0001-01-01 00:00:00'
                                        to   'days since 0001-01-00 00:00:00'
             Matthias Cuntz, Feb 2013 - solved Excel leap year problem.
             Matthias Cuntz, Feb 2013 - ported to Python 3
             Matthias Cuntz, Jul 2013 - ascii/eng without time defaults to 00:00:00
             Matthias Cuntz, Oct 2013 - Excel starts at 1 not at 0
             Matthias Cuntz, Oct 2013 - units bugs, e.g. 01.01.0001 was substracted if Julian calendar even with units
             Matthias Cuntz, Nov 2013 - removed remnant of time treatment before time check in eng keyword
             Matthias Cuntz, Jun 2015 - adapted to new netCDF4/netcdftime (>= v1.0) and datetime (>= Python v2.7.9)
             Matthias Cuntz, Oct 2015 - call date2num with list instead of single netCDF4.datetime objects
             Matthias Cuntz, Oct 2016 - netcdftime provided even with netCDF4 > 1.0.0; make mo for months always integer
             Matthias Cuntz, Nov 2016 - 00, 01, etc. for integers not accepted by Python3
             Matthias Cuntz, May 2020 - numpy docstring format
             Matthias Cuntz, Jul 2020 - en for eng
             Matthias Cuntz, Jul 2020 - use proleptic_gregorian for Excel dates
    """
    #
    # Checks
    calendars = ['standard', 'gregorian', 'julian', 'proleptic_gregorian',
                 'excel1900', 'excel1904', '365_day', 'noleap', '366_day',
                 'all_leap', '360_day', 'decimal', 'decimal360']
    import netCDF4 as nt
    try:
        tst = nt.date2num
        tst = nt.datetime
    except:
        try:
            import netcdftime as nt
            if ((nt.__version__ <= '0.9.2') & (calendar == '360_day')):
                raise ValueError("date2dec error: Your version of netcdftime.py is equal"
                                 " or below 0.9.2. The 360_day calendar does not work with"
                                 " arrays here. Please download a newer one.")
        except:
            import cftime as nt
    #
    calendar = calendar.lower()
    if (calendar not in calendars):
        raise ValueError("date2dec error: Wrong calendar!"
                    " Choose: "+''.join([i+' ' for i in calendars]))
    # obsolete eng
    if (eng is not None):
        if (en is not None):
            raise ValueError("date2dec error: 'eng' was succeeded by 'en'. Only one can be given.")
        else:
            en = eng
    # if ascii input is given by user, other input will be neglected
    # calculation of input size and shape
    if (ascii is not None) and (en is not None):
        raise ValueError("date2dec error: 'ascii' and 'en' mutually exclusive")
    if (ascii is not None):
        islist = type(ascii) != type(np.array(ascii))
        isarr = np.ndim(ascii)
        if (islist & (isarr > 2)):
            raise ValueError("date2dec error: ascii input is list > 2D; Use array input")
        if isarr == 0:
            ascii = np.array([ascii])
        else:
            ascii = np.array(ascii)
        insize   = ascii.size
        outsize  = insize
        outshape = ascii.shape
        asciifl  = ascii.flatten()
        timeobj  = np.zeros(insize, dtype=object)
        # slicing of ascii strings to implement in datetime object. missing times
        # will be set to 0.
        yr = np.zeros(insize, dtype=np.int)
        mo = np.zeros(insize, dtype=np.int)
        dy = np.zeros(insize, dtype=np.int)
        hr = np.zeros(insize, dtype=np.int)
        mi = np.zeros(insize, dtype=np.int)
        sc = np.zeros(insize, dtype=np.int)
        for i in range(insize):
            aa      = asciifl[i].split('.')
            dy[i]   = int(aa[0])
            mo[i]   = int(aa[1])
            tail    = aa[2].split()
            yr[i]   = int(tail[0])
            if len(tail) > 1:
                tim     = tail[1].split(':')
                hr[i]   = int(tim[0])
                if len(tim) > 1:
                    mi[i] = int(tim[1])
                else:
                    mi[i] = 0
                if len(tim) > 2:
                    sc[i] = int(tim[2])
                else:
                    sc[i] = 0
            else:
                hr[i] = 0
                mi[i] = 0
                sc[i] = 0
            if ((yr[i]==1900) & (mo[i]==2) & (dy[i]==29)):
                timeobj[i] = nt.datetime(yr[i], 3, 1, hr[i], mi[i], sc[i])
            else:
                timeobj[i] = nt.datetime(yr[i], mo[i], dy[i], hr[i], mi[i], sc[i])
    if (en is not None):
        islist = type(en) != type(np.array(en))
        isarr  = np.ndim(en)
        if isarr == 0:
             en = np.array([en])
        else:
             en = np.array(en)
        if (islist & (isarr > 2)):
            raise ValueError("date2dec error: en input is list > 2D; Use array input")
        en = np.array(en)
        insize   = en.size
        outsize  = insize
        outshape = en.shape
        enfl     = en.flatten()
        timeobj  = np.zeros(insize, dtype=object)
        # slicing of en strings to implement in datetime object. missing times
        # will be set to 0.
        yr = np.zeros(insize, dtype=np.int)
        mo = np.zeros(insize, dtype=np.int)
        dy = np.zeros(insize, dtype=np.int)
        hr = np.zeros(insize, dtype=np.int)
        mi = np.zeros(insize, dtype=np.int)
        sc = np.zeros(insize, dtype=np.int)
        for i in range(insize):
            ee      = enfl[i].split('-')
            yr[i]   = int(ee[0])
            mo[i]   = int(ee[1])
            tail    = ee[2].split()
            dy[i]   = int(tail[0])
            if len(tail) > 1:
                tim     = tail[1].split(':')
                hr[i]   = int(tim[0])
                if len(tim) > 1:
                    mi[i] = int(tim[1])
                else:
                    mi[i] = 0
                if len(tim) > 2:
                    sc[i] = int(tim[2])
                else:
                    sc[i] = 0
            else:
                hr[i] = 0
                mi[i] = 0
                sc[i] = 0
            if ((yr[i]==1900) & (mo[i]==2) & (dy[i]==29)):
                timeobj[i] = nt.datetime(yr[i], 3, 1, hr[i], mi[i], sc[i])
            else:
                timeobj[i] = nt.datetime(yr[i], mo[i], dy[i], hr[i], mi[i], sc[i])
    # if no ascii input, other inputs will be concidered
    # calculation of input sizes, shapes and number of axis
    if (ascii is None) and (en is None):
        isnlist1 = type(yr) == type(np.array(yr))
        isarr1   = np.ndim(yr)
        if isarr1 == 0: yr = np.array([yr])
        else: yr = np.array(yr)
        isnlist2 = type(mo) == type(np.array(mo))
        isarr2   = np.ndim(mo)
        if isarr2 == 0: mo = np.array([mo], dtype=np.int)
        else: mo = np.array(mo, dtype=np.int)
        isnlist3 = type(dy) == type(np.array(dy))
        isarr3   = np.ndim(dy)
        if isarr3 == 0: dy = np.array([dy])
        else: dy = np.array(dy)
        isnlist4 = type(hr) == type(np.array(hr))
        isarr4   = np.ndim(hr)
        if isarr4 == 0: hr = np.array([hr])
        else: hr = np.array(hr)
        isnlist5 = type(mi) == type(np.array(mi))
        isarr5   = np.ndim(mi)
        if isarr5 == 0: mi = np.array([mi])
        else: mi = np.array(mi)
        isnlist6 = type(sc) == type(np.array(sc))
        isarr6   = np.ndim(sc)
        if isarr6 == 0: sc = np.array([sc])
        else: sc = np.array(sc)
        islist = not (isnlist1 | isnlist2 | isnlist3 | isnlist4 | isnlist5 | isnlist6)
        isarr  = isarr1 + isarr2 + isarr3 + isarr4 + isarr5 + isarr6
        shapes = [np.shape(yr), np.shape(mo), np.shape(dy), np.shape(hr), np.shape(mi), np.shape(sc)]
        nyr    = np.size(yr)
        nmo    = np.size(mo)
        ndy    = np.size(dy)
        nhr    = np.size(hr)
        nmi    = np.size(mi)
        nsc    = np.size(sc)
        sizes  = [nyr,nmo,ndy,nhr,nmi,nsc]
        nmax   = np.amax(sizes)
        ii     = np.argmax(sizes)
        outshape = shapes[ii]
        if (islist & (np.size(outshape) > 2)):
            raise ValueError("date2dec error: input is list > 2D; Use array input.")
        if nyr < nmax:
            if nyr == 1: yr  = np.ones(outshape,)*yr
            else: raise ValueError("date2dec error: size of yr != max input or 1.")
        if nmo < nmax:
            if nmo == 1: mo  = np.ones(outshape, dtype=np.int)*mo
            else: raise ValueError("date2dec error: size of mo != max input or 1.")
        if ndy < nmax:
            if ndy == 1: dy  = np.ones(outshape)*dy
            else: raise ValueError("date2dec error: size of dy != max input or 1.")
        if nhr < nmax:
            if nhr == 1: hr  = np.ones(outshape)*hr
            else: raise ValueError("date2dec error: size of hr != max input or 1.")
        if nmi < nmax:
            if nmi == 1: mi  = np.ones(outshape)*mi
            else: raise ValueError("date2dec error: size of mi != max input or 1.")
        if nsc < nmax:
            if nsc == 1: sc  = np.ones(outshape)*sc
            else: raise ValueError("date2dec error: size of sc != max input or 1.")
        indate  = [yr, mo, dy, hr, mi, sc]
        insize  = [np.size(i) for i in indate]
        inshape = [np.shape(i) for i in indate]
        dimnum  = [len(i) for i in inshape]
        # calculation of maximum input size and maximum number of axis for
        # reshaping the output
        indate  = [i.flatten() for i in indate]
        outsize = max(insize)
        timeobj = np.zeros(outsize, dtype=object)
        # datetime object is constructed
        for i in range(outsize):
            iyr = int(indate[0][i])
            imo = int(indate[1][i])
            idy = int(indate[2][i])
            ihr = int(indate[3][i])
            imi = int(indate[4][i])
            isc = int(indate[5][i])

            if ((iyr==1900) & (imo==2) & (idy==29)):
                timeobj[i] = nt.datetime(iyr, 3, 1, ihr, imi, isc)
            else:
                timeobj[i] = nt.datetime(iyr, imo, idy, ihr, imi, isc)
    # depending on chosen calendar and optional set of the time units
    # decimal date is calculated
    output = np.zeros(outsize)
    t0    = nt.datetime(1582, 10, 5, 0, 0, 0)
    t1    = nt.datetime(1582, 10, 15, 0, 0, 0)
    is121 = True if (min(timeobj)<t0) and (max(timeobj)>=t1) else False
    if (calendar == 'standard') or (calendar == 'gregorian'):
        if not units:
            units = 'days since 0001-01-01 12:00:00'
            dec0 = 1721424
        else:
            dec0 = 0
        if is121 and (nt.__version__ < '1.2.2'):
            for ii, tt in enumerate(timeobj): output[ii] = nt.date2num(tt, units, calendar='gregorian')+dec0
        else:
            output = nt.date2num(timeobj, units, calendar='gregorian')+dec0
    elif calendar == 'julian':
        if not units:
            units = 'days since 0001-01-01 12:00:00'
            dec0 = 1721424
        else:
            dec0 = 0
        if is121 and (nt.__version__ < '1.2.2'):
            for ii, tt in enumerate(timeobj): output[ii] = nt.date2num(tt, units, calendar='julian')+dec0
        else:
            output = nt.date2num(timeobj, units, calendar='julian')+dec0
    elif calendar == 'proleptic_gregorian':
        if not units: units = 'days since 0001-01-01 00:00:00'
        if is121 and (nt.__version__ < '1.2.2'):
            for ii, tt in enumerate(timeobj): output[ii] = nt.date2num(tt, units, calendar='proleptic_gregorian')
        else:
            output = nt.date2num(timeobj, units, calendar='proleptic_gregorian')
    elif calendar == 'excel1900':
        doerr = False
        if not units:
            units = 'days since 1899-12-31 00:00:00'
            if excelerr: doerr = True
        if is121 and (nt.__version__ < '1.2.2'):
            for ii, tt in enumerate(timeobj): output[ii] = nt.date2num(tt, units, calendar='proleptic_gregorian')
        else:
            output = nt.date2num(timeobj, units, calendar='proleptic_gregorian')
        if doerr:
            output = np.where(output >= 60., output+1., output)
            # date2num treats 29.02.1900 as 01.03.1990, i.e. is the same decimal number
            if np.any((output >= 61.) & (output < 62.)):
                for i in range(outsize):
                    # if (timeobj[i].year==1900) & (timeobj[i].month==2) & (timeobj[i].day==29):
                    #     output[i] -= 1.
                    if (yr[i]==1900) & (mo[i]==2) & (dy[i]==29):
                        output[i] -= 1.
    elif calendar == 'excel1904':
        if not units: units = 'days since 1903-12-31 00:00:00'
        if is121 and (nt.__version__ < '1.2.2'):
            for ii, tt in enumerate(timeobj): output[ii] = nt.date2num(tt, units, calendar='proleptic_gregorian')
        else:
            output = nt.date2num(timeobj, units, calendar='proleptic_gregorian')
    elif (calendar == '365_day') or (calendar == 'noleap'):
        if not units: units = 'days since 0001-01-01 00:00:00'
        if is121 and (nt.__version__ < '1.2.2'):
            for ii, tt in enumerate(timeobj): output[ii] = nt.date2num(tt, units, calendar='365_day')
        else:
            output = nt.date2num(timeobj, units, calendar='365_day')
    elif (calendar == '366_day') or (calendar == 'all_leap'):
        if not units: units = 'days since 0001-01-01 00:00:00'
        if is121 and (nt.__version__ < '1.2.2'):
            for ii, tt in enumerate(timeobj): output[ii] = nt.date2num(tt, units, calendar='366_day')
        else:
            output = nt.date2num(timeobj, units, calendar='366_day')
    elif calendar == '360_day':
        if not units: units = 'days since 0001-01-01 00:00:00'
        if is121 and (nt.__version__ < '1.2.2'):
            for ii, tt in enumerate(timeobj): output[ii] = nt.date2num(tt, units, calendar='360_day')
        else:
            output = nt.date2num(timeobj, units, calendar='360_day')
    elif calendar == 'decimal':
        ntime = np.size(yr)
        leap  = np.array((((yr%4)==0) & ((yr%100)!=0)) | ((yr%400)==0)).astype(np.int)
        tdy   = np.array(dy, dtype=np.float)
        diy   = np.array([ [-9,0, 31, 59, 90,120,151,181,212,243,273,304,334,365],
                           [-9,0, 31, 60, 91,121,152,182,213,244,274,305,335,366] ])
        for i in range(ntime):
            tdy[i] = tdy[i] + np.array(diy[leap[i],mo[i]], dtype=np.float)
        days_year = 365.
        output    = ( np.array(yr, dtype=np.float) +
                      ((tdy-1.)*24. + np.array(hr, dtype=np.float) +
                       np.array(mi, dtype=np.float)/60. +
                       np.array(sc, dtype=np.float)/3600.) /
                       ((days_year+np.array(leap, dtype=np.float))*24.) )
        # for numerical stability, i.e. back and forth transforms
        output += 1e-08 # 1/3 sec
    elif calendar == 'decimal360':
        ntime = np.size(yr)
        tdy   = np.array(dy, dtype=np.float)
        diy   = np.array([-9,  0, 30, 60, 90,120,150,180,210,240,270,300,330,360])
        for i in range(ntime):
            tdy[i] = tdy[i] + np.array(diy[mo[i]], dtype=np.float)
        days_year = 360.
        output    = ( np.array(yr, dtype=np.float) +
                      ((tdy-1.)*24. + np.array(hr, dtype=np.float) +
                       np.array(mi, dtype=np.float)/60. +
                       np.array(sc, dtype=np.float)/3600.) /
                       (days_year*24.) )
        # for numerical stability, i.e. back and forth transforms
        output += 1e-08 # 1/3 sec
    else:
        raise ValueError("date2dec error: calendar not implemented; should have been catched before.")


    # return of reshaped output
    output = np.reshape(output, outshape)
    if isarr == 0:
        output = np.float(output)
    else:
        if islist:
            ns = np.size(outshape)
            if ns == 1:
                output = [i for i in output]
            else:
                loutput = [ i for i in output[:,0]]
                for i in range(np.size(output[:,0])):
                    loutput[i] = list(np.squeeze(output[i,:]))
                output = loutput

    return output


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
