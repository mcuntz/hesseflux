#!/usr/bin/env python
"""
mad : Median absolute deviation test, either on raw values or on 1st or 2nd derivatives.

This module was written by Matthias Cuntz while at Department of
Computational Hydrosystems, Helmholtz Centre for Environmental
Research - UFZ, Leipzig, Germany, and continued while at Institut
National de Recherche pour l'Agriculture, l'Alimentation et
l'Environnement (INRAE), Nancy, France.

Copyright (c) 2011-2020 Matthias Cuntz - mc (at) macu (dot) de
Released under the MIT License; see LICENSE file for details.

* Written Nov 2011 by Matthias Cuntz - mc (at) macu (dot) de
* ND-array, act on axis=0, May 2012, Matthias Cuntz
* Removed bug in broadcasting, Jun 2012, Matthias Cuntz
* Better usage of numpy possibilities, e.g. using np.diff, Jun 2012, Matthias Cuntz
* Ported to Python 3, Feb 2013, Matthias Cuntz
* Use bottleneck for medians, otherwise loop over axis=1, Jul 2013, Matthias Cuntz and Juliane Mai
* Re-allow masked arrays and arrays with NaNs, Jul 2013, Matthias Cuntz
* Removed bug in NaN treatment, Oct 2013, Matthias Cuntz
* Keyword nonzero, Oct 2013, Matthias Cuntz
* Using numpy docstring format, May 2020, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz

The following functions are provided

.. autosummary::
   mad
"""
from __future__ import division, absolute_import, print_function
import numpy as np


__all__ = ['mad']


def mad(datin, z=7, deriv=0, nozero=False):
    """
    Median absolute deviation test, either on raw values,
    or on 1st or 2nd derivatives.

    Returns mask with False everywhere except where `<(median-MAD\*z/0.6745)`
    or `>(md+MAD\*z/0.6745)`.

    Parameters
    ----------
    datin : array or masked array
        `mad` acts on `axis=0`.
    z : float, optional
        Input is allowed to deviate maximum `z` standard deviations from the median (default: 7)
    deriv : int, optional
        0: Act on raw input (default).

        1: Use first derivatives.

        2: Use 2nd derivatives.
    nozero : bool, optional
        True: exclude zeros (0.) from input `datin`.

    Returns
    -------
    array of bool
        False everywhere except where input deviates more than `z` standard deviations from median

    Notes
    -----
    If input is an array then mad is checked along the zeroth axis for outlier.

    1st derivative is calculated as `d = datin[1:n]-datin[0:n-1]`
    because mean of left and right would give 0 for spikes.

    If `all(d.mask==True)` then return `d.mask`, which is all True.

    Examples
    --------
    >>> import numpy as np
    >>> y = np.array([-0.25,0.68,0.94,1.15,2.26,2.35,2.37,2.40,2.47,2.54,2.62,
    ...               2.64,2.90,2.92,2.92,2.93,3.21,3.26,3.30,3.59,3.68,4.30,
    ...               4.64,5.34,5.42,8.01],dtype=np.float)

    >>> # Normal MAD
    >>> print(mad(y))
    [False False False False False False False False False False False False
     False False False False False False False False False False False False
     False False]

    >>> print(mad(y,z=4))
    [False False False False False False False False False False False False
     False False False False False False False False False False False False
     False  True]

    >>> print(mad(y,z=3))
    [ True False False False False False False False False False False False
     False False False False False False False False False False False False
      True  True]

    >>> # MAD on 2nd derivatives
    >>> print(mad(y,z=4,deriv=2))
    [False False False False False False False False False False False False
     False False False False False False False False False False False  True]

    >>> # direct usage
    >>> my = np.ma.array(y, mask=mad(y,z=4))
    >>> print(my)
    [-0.25 0.68 0.94 1.15 2.26 2.35 2.37 2.4 2.47 2.54 2.62 2.64 2.9 2.92 2.92
     2.93 3.21 3.26 3.3 3.59 3.68 4.3 4.64 5.34 5.42 --]

    >>> # MAD on several dimensions
    >>> yy = np.transpose(np.array([y,y]))
    >>> print(np.transpose(mad(yy,z=4)))
    [[False False False False False False False False False False False False
      False False False False False False False False False False False False
      False  True]
     [False False False False False False False False False False False False
      False False False False False False False False False False False False
      False  True]]

    >>> yyy = np.transpose(np.array([y,y,y]))
    >>> print(np.transpose(mad(yyy,z=3)))
    [[ True False False False False False False False False False False False
      False False False False False False False False False False False False
       True  True]
     [ True False False False False False False False False False False False
      False False False False False False False False False False False False
       True  True]
     [ True False False False False False False False False False False False
      False False False False False False False False False False False False
       True  True]]

    >>> # Masked arrays
    >>> my = np.ma.array(y, mask=np.zeros(y.shape))
    >>> my.mask[-1] = True
    >>> print(mad(my,z=4))
    [True False False False False False False False False False False False
     False False False False False False False False False False False False
     False --]

    >>> print(mad(my,z=3))
    [True False False False False False False False False False False False
     False False False False False False False False False False False True
     True --]

    >>> # Arrays with NaNs
    >>> ny = y.copy()
    >>> ny[-1] = np.nan
    >>> print(mad(ny,z=4))
    [ True False False False False False False False False False False False
     False False False False False False False False False False False False
     False False]

    >>> print(mad(ny,z=3))
    [ True False False False False False False False False False False False
     False False False False False False False False False False False  True
      True False]

    >>> # Exclude zeros
    >>> zy = y.copy()
    >>> zy[1] = 0.
    >>> print(mad(zy,z=3))
    [ True  True False False False False False False False False False False
     False False False False False False False False False False False False
      True  True]

    >>> print(mad(zy,z=3,nozero=True))
    [ True False False False False False False False False False False False
     False False False False False False False False False False False False
      True  True]

    History
    -------
    Written,  Matthias Cuntz, Nov 2011
    Modified, Matthias Cuntz, May 2012 - act on axis=0 of array
              Matthias Cuntz, Jun 2012 - axis=0 did not always work: spread md and MAD to input dimensions
              Matthias Cuntz, Jun 2012 - use np.diff, remove spreads
              Matthias Cuntz, Feb 2013 - ported to Python 3
              Matthias Cuntz & Juliane Mai
                              Jul 2013 - loop over second dimension for medians, faster than array calculations :-(
                                         but use bottleneck for speed :-)
              Matthias Cuntz, Jul 2013 - (re-)allow masked arrays and NaNs in arrays
              Matthias Cuntz, Oct 2013 - nozero, bug in NaN treatment with dim=1
              Matthias Cuntz, May 2020 - numpy docstring format
    """
    if nozero:
        idatin = datin.copy()
        ii = np.where(idatin == 0.)[0]
        if ii.size > 0: idatin[ii] = np.nan
    else:
        idatin = datin
    sn = list(np.shape(idatin))
    n  = sn[0]
    if deriv == 0:
        m      = n
        d      = idatin
    elif deriv == 1:
        m      = n-1
        sm     = sn
        sm[0]  = m
        d      = np.diff(idatin, axis=0)
    elif deriv == 2:
        m      = n-2
        sm     = sn
        sm[0]  = m
        d      = np.diff(idatin, n=2, axis=0)
    else:
        raise ValueError('Unimplemented option.')


    # Shortcut if all masked
    ismasked = isinstance(d, np.ma.core.MaskedArray)
    if not ismasked:
        ii = np.where(~np.isfinite(d))[0]
        d  = np.ma.array(d)
        if ii.size > 0: d[ii] = np.ma.masked

    if np.all(d.mask == True):
        if ismasked:
            return d.mask
        else:
            return np.ones(d.shape, dtype=np.bool)

    # Median
    oldsettings = np.geterr()
    np.seterr(invalid='ignore')
    if d.ndim == 1:
        try:
            import bottleneck as bn
            dd = d.compressed()
            md = bn.median(dd)
            # Median absolute deviation
            MAD = bn.median(np.abs(dd-md))
            # Range around median
            thresh = MAD * (z/0.6745)
            # True where outside z-range
            res = (d<(md-thresh)) | (d>(md+thresh))
        except:
            dd = d.compressed()
            md = np.median(dd)
            # Median absolute deviation
            MAD = np.median(np.abs(dd-md))
            # Range around median
            thresh = MAD * (z/0.6745)
            # True where outside z-range
            res = (d<(md-thresh)) | (d>(md+thresh))
    elif d.ndim == 2:
        try:
            import bottleneck as bn
            res = np.empty(d.shape, dtype=np.bool)
            for i in range(d.shape[1]):
                di = d[:,i]
                dd = di.compressed()
                md = bn.median(dd)
                # Median absolute deviation
                MAD = bn.median(np.abs(dd-md))
                # Range around median
                thresh = MAD * (z/0.6745)
                # True where outside z-range
                res[:,i] = (d[:,i]<(md-thresh)) | (d[:,i]>(md+thresh))
        except:
            res = np.empty(d.shape, dtype=np.bool)
            for i in range(d.shape[1]):
                di = d[:,i]
                dd = di.compressed()
                md = np.median(dd)
                # Median absolute deviation
                MAD = np.median(np.abs(dd-md))
                # Range around median
                thresh = MAD * (z/0.6745)
                # True where outside z-range
                res[:,i] = (d[:,i]<(md-thresh)) | (d[:,i]>(md+thresh))
    else:
        np.seterr(**oldsettings)
        raise ValueError('datin.ndim must be <= 2')

    np.seterr(**oldsettings)
    if ismasked:
        return res
    else:
        if isinstance(res, np.ma.core.MaskedArray): # got masked because of NaNs
            return np.where(res.mask, False, res)
        else:
            return res


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # import numpy as np
    # y = np.array([-0.25,0.68,0.94,1.15,2.26,2.35,2.37,2.40,2.47,2.54,2.62,
    #                2.64,2.90,2.92,2.92,2.93,3.21,3.26,3.30,3.59,3.68,4.30,
    #                4.64,5.34,5.42,8.01],dtype=np.float)
    # print(mad(y))
    # #[False False False False False False False False False False False False
    # # False False False False False False False False False False False False
    # # False False]

    # print(mad(y,z=4))
    # #[False False False False False False False False False False False False
    # # False False False False False False False False False False False False
    # # False  True]

    # print(mad(y,z=3))
    # #[ True False False False False False False False False False False False
    # # False False False False False False False False False False False False
    # #  True  True]

    # print(mad(y,z=4,deriv=2))
    # #[False False False False False False False False False False False False
    # # False False False False False False False False False False False  True]

    # my = np.ma.array(y, mask=mad(y,z=4))
    # print(my)
    # #[-0.25 0.68 0.94 1.15 2.26 2.35 2.37 2.4 2.47 2.54 2.62 2.64 2.9 2.92 2.92
    # # 2.93 3.21 3.26 3.3 3.59 3.68 4.3 4.64 5.34 5.42 --]

    # yy = np.transpose(np.array([y,y]))
    # print(np.transpose(mad(yy,z=4)))
    # #[[False False False False False False False False False False False False
    # #  False False False False False False False False False False False False
    # #  False  True]
    # # [False False False False False False False False False False False False
    # #  False False False False False False False False False False False False
    # #  False  True]]

    # yyy = np.transpose(np.array([y,y,y]))
    # print(np.transpose(mad(yyy,z=3)))
    # #[[ True False False False False False False False False False False False
    # #  False False False False False False False False False False False False
    # #   True  True]
    # # [ True False False False False False False False False False False False
    # #  False False False False False False False False False False False False
    # #   True  True]
    # # [ True False False False False False False False False False False False
    # #  False False False False False False False False False False False False
    # #   True  True]]

    # my = np.ma.array(y, mask=np.zeros(y.shape))
    # my.mask[-1] = True
    # print(mad(my,z=4))
    # #[False False False False False False False False False False False False
    # # False False False False False False False False False False False False
    # # False --]

    # print(mad(my,z=3))
    # #[True False False False False False False False False False False False
    # # False False False False False False False False False False False True
    # # True --]

    # ny = y
    # ny[-1] = np.nan
    # print(mad(ny,z=4))
    # #[False False False False False False False False False False False False
    # # False False False False False False False False False False False False
    # # False False]

    # print(mad(ny,z=3))
    # #[True False False False False False False False False False False False
    # # False False False False False False False False False False False True
    # # True False]
