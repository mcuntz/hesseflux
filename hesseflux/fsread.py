#!/usr/bin/env python
"""
fsread : Read from a file numbers into 2D float array as well as characters into 2D string array.

This module was written by Matthias Cuntz while at Department of
Computational Hydrosystems, Helmholtz Centre for Environmental
Research - UFZ, Leipzig, Germany, and continued while at Institut
National de Recherche pour l'Agriculture, l'Alimentation et
l'Environnement (INRAE), Nancy, France.

Copyright (c) 2015-2020 Matthias Cuntz - mc (at) macu (dot) de
Released under the MIT License; see LICENSE file for details.

* Written Feb 2015 by Matthias Cuntz (mc (at) macu (dot) de)
* nc<=-1 rmoved in case of nc is list, Nov 2016, Matthias Cuntz
* range instead of np.arange, Nov 2017, Matthias Cuntz
* Keywords cname, sname, hstrip, rename file to infile, Nov 2017, Matthias Cuntz
* full_header=True returns vector of strings, Nov 2017, Matthias Cuntz
* Ignore unicode characters on read, Jun 2019, Matthias Cuntz
* Keywords encoding, errors with codecs module, Aug 2019, Matthias Cuntz
* return_list=False default, Jan 2020, Matthias Cuntz
* Using numpy docstring format, May 2020, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz

The following functions are provided

.. autosummary::
   fsread
"""
from __future__ import division, absolute_import, print_function
import numpy as np


__all__ = ['fsread']


def fsread(infile, nc=0, cname=None, snc=0, sname=None, skip=0, cskip=0, hskip=0, hstrip=True, separator=None,
           squeeze=False, reform=False, skip_blank=False, comment=None,
           fill=False, fill_value=0, sfill_value='', strip=None, encoding='ascii', errors='ignore',
           header=False, full_header=False,
           transpose=False, strarr=False):
    """
    Read from a file numbers into 2D float array
    as well as characters into 2D string array.

    Parameters
    ----------
    infile : str
        source file name
    nc : int or iterable, optional
        number of columns to be read as floats [default: all (`nc<=0`)].

        `nc` can be an int or a vector of column indexes, starting with 0;
        `cskip` will be ignored in the latter case.

        If `snc!=0`: `nc` must be iterable or -1 to read all other columns as floats.
    cname : iterable of str, optional
        float columns can be chosen by the values in the first header line;
        must be iterable with strings.
    snc : int or iterable, optional
        number of columns to be read as strings [default: none (`snc=0`)].

        `snc` can be an int or a vector of column indexes, starting with 0;
        `cskip` will be ignored in the latter case.

        If `nc!=0`: `snc` must be iterable or -1 to read all other columns as strings.
    sname : iterable of str, optional
        string columns can be chosen by the values in the first header line;
        must be iterable with strings.
    skip : int, optional
        number of lines to skip at the beginning of file (default: 0)
    cskip : int, optional
        number of columns to skip at the beginning of each line (default: 0)
    hskip : int, optional
        number of lines in skip that do not belong to header (default: 0)
    hstrip : bool, optional
        True: strip header cells to match with cname (default: True)
    separator : str, optional
        column separator. If not given, columns separator are (in order):
        comma (','), semicolon (';'), whitespace.
    comment : iterable, optional
         line gets excluded if first character of line is in comment sequence.
         Sequence must be iterable such as string, list and tuple.
    fill_value : float, optional
         value to fill in float array in empty cells or if not enough columns in line
         and `fill==True` (default: 0).
    sfill_value : str, optional
         value to fill in string array in empty cells or if not enough columns in line
         and `fill==True` (default: '').
    strip : str, optional
        Strip strings with str.strip(strip).

        None: strip quotes " and ' (default).

        False: no strip (~30% faster).

        str: strip character given by `strip`.
    encoding : str, optional
        Specifies the encoding which is to be used for the file (default: 'ascii').
        Any encoding that encodes to and decodes from bytes is allowed.
    errors : str, optional
        Errors may be given to define the error handling during encoding of the file (default: 'ignore').

        Possible values: 'strict', 'replace', 'ignore'.
    squeeze : bool, optional
        True:  2-dim array will be cleaned of degenerated dimension, i.e. results in a vector.

        False: array will be two-dimensional as read (default).
    reform : bool, optional
        Same as squeeze.
    skip_blank : bool, optional
        True:  continues reading after blank line.

        False: stops reading at first blank line (default).
    fill : bool, optional
        True:  fills in `fill_value` if not enough columns in line.

        False: stops execution and returns None if not enough columns in line (default).
    header : bool, optional
        True:  header strings will be returned.

        False: numbers in file will be returned (default).
    full_header : bool, optional
        True:  header is a string vector of the skipped rows.

        False: header will be split in columns, exactly as the data,
               and will hold only the selected columns (default).
    transpose : bool, optional
        True:  column-major format `output(0:ncolumns,0:nlines)`.

        False: row-major format `output(0:nlines,0:ncolumns)` (default).
    strarr : bool, optional
        True:  return header as numpy array of strings.

        False: return header as list.

    Returns
    -------
    array(s)
        1 output:  array of floats (`nc!=0` and `snc=0`)

        1 output:  array of strings (`nc=0` and `snc!=0`)

        2 outputs: array of floats, array of strings (`nc!=0` and `snc!=0`)

        1 output:  list/string array of header ((`nc=0` or `snc=0`) and `header=True`)

        2 outputs: list/string array of header for float array, list/string array of header for strarr ((`nc!=0` and `snc!=0`) and `header=True`)

        1 output:  String vector of full file header (`header=True` and `full_header=True`)

    Notes
    -----
    If `header==True` then skip is counterintuitive because it is
    actually the number of header rows to be read. This is to
    be able to have the exact same call of the function, once
    with `header=False` and once with `header=True`.
    
    If `fill==True`, blank lines are not filled but are taken as end of file.

    `transpose=True` has no effect on 1D output such as 1 header line.

    Passes file to :func:`fread` if `snc==0`.

    Passes file to :func:`sread` if `nc==0`.

    Examples
    --------
    >>> # Create some data
    >>> filename = 'test.dat'
    >>> ff = open(filename,'w')
    >>> ff.writelines('head1 head2 head3 head4\\n')
    >>> ff.writelines('1.1 1.2 1.3 1.4\\n')
    >>> ff.writelines('2.1 2.2 2.3 2.4\\n')
    >>> ff.close()

    >>> # Read sample with fread - see fread for more examples
    >>> print(fsread(filename, nc=[1,3], skip=1))
    [[1.2 1.4]
     [2.2 2.4]]
    >>> print(fsread(filename, nc=2, skip=1, header=True))
    ['head1', 'head2']

    >>> # Read sample with sread - see sread for more examples
    >>> print(fsread(filename, snc=[1,3], skip=1))
    [['1.2', '1.4'], ['2.2', '2.4']]

    >>> # Some mixed data
    >>> ff = open(filename,'w')
    >>> ff.writelines('head1 head2 head3 head4\\n')
    >>> ff.writelines('01.12.2012 1.2 name1 1.4\\n')
    >>> ff.writelines('01.01.2013 2.2 name2 2.4\\n')
    >>> ff.close()

    >>> # Read columns
    >>> print(fsread(filename, nc=[1,3], skip=1))
    [[1.2 1.4]
     [2.2 2.4]]
    >>> a, sa = fsread(filename, nc=[1,3], snc=[0,2], skip=1)
    >>> print(a)
    [[1.2 1.4]
     [2.2 2.4]]
    >>> print(sa[0][0])
    01.12.2012
    >>> print(sa[0][1])
    name1
    >>> print(sa[1][0])
    01.01.2013
    >>> print(sa[1][1])
    name2
    >>> a, sa = fsread(filename, nc=[1,3], snc=-1, skip=1)
    >>> print(a)
    [[1.2 1.4]
     [2.2 2.4]]
    >>> print(sa[0][0])
    01.12.2012
    >>> print(sa[0][1])
    name1
    >>> print(sa[1][0])
    01.01.2013
    >>> print(sa[1][1])
    name2
    >>> a, sa = fsread(filename, nc=-1, snc=[0,2], skip=1)
    >>> print(a)
    [[1.2 1.4]
     [2.2 2.4]]

    >>> # Read header
    >>> a, sa = fsread(filename, nc=[1,3], snc=[0,2], skip=1, header=True)
    >>> print(a)
    ['head2', 'head4']
    >>> print(sa)
    ['head1', 'head3']

    >>> # Some mixed data with missing values
    >>> ff = open(filename,'w')
    >>> ff.writelines('head1,head2,head3,head4\\n')
    >>> ff.writelines('01.12.2012,1.2,name1,1.4\\n')
    >>> ff.writelines('01.01.2013,,name2,2.4\\n')
    >>> ff.close()

    >>> print(fsread(filename, nc=[1,3], skip=1, fill=True, fill_value=-1))
    [[ 1.2  1.4]
     [-1.   2.4]]

    >>> # cname, sname
    >>> a, sa = fsread(filename, cname='head2', snc=[0,2], skip=1, fill=True, fill_value=-1, squeeze=True)
    >>> print(a)
    [ 1.2 -1. ]
    >>> print(sa)
    [['01.12.2012', 'name1'],
     ['01.01.2013', 'name2']]
    >>> a, sa = fsread(filename, cname=['head2','head4'], snc=-1, skip=1, fill=True, fill_value=-1)
    >>> print(a)
    [[ 1.2  1.4]
     [-1.   2.4]]
    >>> print(sa)
    [['01.12.2012', 'name1'],
     ['01.01.2013', 'name2']]
    >>> a, sa = fsread(filename, nc=[1,3], sname=['head1','head3'], skip=1, fill=True, fill_value=-1, strarr=True, header=True)
    >>> print(a)
    ['head2' 'head4']
    >>> print(sa)
    ['head1' 'head3']
    >>> print(fsread(filename, cname=['head2','head4'], snc=-1, skip=1, header=True, full_header=True))
    ['head1,head2,head3,head4']
    >>> print(fsread(filename, cname=['head2','head4'], snc=-1, skip=1, fill=True, fill_value=-1, header=True, full_header=True))
    ['head1,head2,head3,head4']
    >>> a, sa = fsread(filename, cname=['  head2','head4'], snc=-1, skip=1, fill=True, fill_value=-1, hstrip=False)
    >>> print(a)
    [[1.4]
     [2.4]]

    >>> # Clean up doctest
    >>> import os
    >>> os.remove(filename)

    History
    -------
    Written,  Matthias Cuntz, Feb 2015 - modified fread
    Modified, Matthias Cuntz, Nov 2016 - nc<=-1 does not work in Python 3 if nc is list
              Matthias Cuntz, Nov 2017 - use range instead of np.arange for producing indexes
              Matthias Cuntz, Nov 2017 - full_header=True returns vector of strings
                                       - cname, sname, file->infile, hstrip
              Matthias Cuntz, Jun 2019 - open(errors='ignore') to ignore unicode characters, for example, on read
                                         -> returns header in unicode in Python2
              Matthias Cuntz, Aug 2019 - use codecs module and allow user encoding and error handling
              Matthias Cuntz, May 2020 - numpy docstring format
    """

    # Input error
    if (nc == -1) and (snc == -1):
        raise ValueError('nc and snc must numbers or list of indices; -1 means read the rest. nc and snc cannot both be -1.')

    # wrap to fread/sread
    if (not isinstance(snc, (list, tuple, np.ndarray))) and (sname is None):
        if snc==0:
            try:    # import package
                from .fread import fread
            except: # python fsread.py
                from fread import fread
            return fread(infile, nc=nc, cname=cname, skip=skip, cskip=cskip, hskip=hskip, hstrip=hstrip,
                         separator=separator,
                         squeeze=squeeze, reform=reform, skip_blank=skip_blank, comment=comment,
                         fill=fill, fill_value=fill_value, strip=strip, encoding=encoding, errors=errors,
                         header=header, full_header=full_header,
                         transpose=transpose, strarr=strarr)
    # snc!=0
    if (not isinstance(nc, (list, tuple, np.ndarray))) and (cname is None):
        if nc==0:
            try:    # import package
                from .sread import sread
            except: # python fsread.py
                from sread import sread
            return sread(infile, nc=snc, cname=sname, skip=skip, cskip=cskip, hskip=hskip, hstrip=hstrip,
                         separator=separator,
                         squeeze=squeeze, reform=reform, skip_blank=skip_blank, comment=comment,
                         fill=fill, fill_value=sfill_value, strip=strip, encoding=encoding, errors=errors,
                         header=header, full_header=full_header,
                         transpose=transpose, strarr=strarr)

    # Open file
    import codecs
    f = codecs.open(infile, 'r', encoding=encoding, errors=errors)

    # Read header and Skip lines
    if hskip > 0:
        ihskip = 0
        while ihskip < hskip:
            tmp = f.readline().rstrip()
            ihskip += 1
    if skip > 0:
        head = ['']*(skip-hskip)
        iskip = 0
        while iskip < (skip-hskip):
            head[iskip] = str(f.readline().rstrip())
            iskip += 1

    # read first line to determine nc and separator (if not set)
    split = -1
    while True:
        s = str(f.readline().rstrip())
        if len(s) == 0:
            if skip_blank:
                continue
            else:
                break
        if comment is not None:
            if (s[0] in comment): continue
        break
    if separator is None:
        sep = ','
        res = s.split(sep)
        nres = len(res)
        if nres == 1:
            sep = ';'
            res = s.split(sep)
            nres = len(res)
            if nres == 1:
                sep = None
                res = s.split(sep)
                nres = len(res)
    else:
        sep = separator
        res = s.split(sep)
        nres = len(res)

    # Determine indices
    if nc != 0 and cname is not None:
        f.close()
        raise ValueError('nc and cname are mutually exclusive.')
    if snc != 0 and sname is not None:
        f.close()
        raise ValueError('snc and sname are mutually exclusive.')
    if (cname is not None) or (sname is not None):
        # from first header line
        if (skip-hskip) <= 0:
            f.close()
            raise IOError('No header line left for choosing columns by name.')
        hres = head[0].split(sep)
        if hstrip: hres = [ h.strip() for h in hres ]
    if cname is not None:
        if not isinstance(cname, (list, tuple, np.ndarray)): cname = [cname]
        if hstrip: cname = [ h.strip() for h in cname ]
        nc = []
        for k in range(len(hres)):
            if hres[k] in cname: nc.append(k)
    if sname is not None:
        if not isinstance(sname, (list, tuple, np.ndarray)): sname = [sname]
        if hstrip: sname = [ h.strip() for h in sname ]
        snc = []
        for k in range(len(hres)):
            if hres[k] in sname: snc.append(k)
    if isinstance(nc, (list, tuple, np.ndarray)) and isinstance(snc, (list, tuple, np.ndarray)):
        if np.in1d(nc, snc, assume_unique=True).any():
            raise ValueError('float and string indices overlap.')
        iinc  = nc
        iisnc = snc
    elif isinstance(nc, (list, tuple, np.ndarray)):
        iinc   = nc
        iirest = list(range(nres))
        for ii in iinc[::-1]: del iirest[ii]
        if snc <= -1:
            iisnc = iirest
        else:
            iisnc = iirest[:snc]
    elif isinstance(snc, (list, tuple, np.ndarray)):
        iisnc  = snc
        iirest = list(range(nres))
        for ii in iisnc[::-1]: del iirest[ii]
        if nc <= -1:
            iinc = iirest
        else:
            iinc = iirest[:nc]
    else:
        # cannot be nc=-1 and snc=-1
        if nc <= -1:
            iisnc = range(snc)
            iinc  = range(snc,nres)
        else:
            if snc <= -1:
                iinc  = range(nc)
                iisnc = range(nc,nres)
            else:
                # red snc first then nc
                iisnc = range(snc)
                iinc  = range(snc,snc+nc)
    nnc    = len(iinc)
    nsnc   = len(iisnc)
    miinc  = max(iinc)
    miisnc = max(iisnc)
    aiinc  = list()
    aiinc.extend(iinc)
    aiinc.extend(iisnc)
    miianc = max(aiinc)
    #
    # Header
    if header:
        # Split header
        var  = None
        svar = None
        if (skip-hskip) > 0:
            if full_header:
                var  = head
                if strarr:
                    var  = np.array(var, dtype=np.str)
                return var
            else:
                var  = list()
                svar = list()
                k = 0
                while k < (skip-hskip):
                    hres = head[k].split(sep)
                    nhres = len(hres)
                    if (miianc >= nhres) and (not fill):
                        f.close()
                        raise ValueError('Line has not enough columns to index: '+head[k])
                    null = line2var(hres, var, iinc, strip)
                    null = line2var(hres, svar, iisnc, False if strip is None else strip)
                    k += 1
                if (skip-hskip) == 1:
                    var  = var[0]
                    svar = svar[0]
        f.close()
        if strarr:
            var  = np.array(var, dtype=np.str)
            if transpose:
                var  = var.T
            if squeeze or reform:
                var  = var.squeeze()
            if fill:
                var  = np.where(var=='', fill_value, var)
            if svar is not None:
                svar = np.array(svar, dtype=np.str)
                if transpose:
                    svar = svar.T
                if squeeze or reform:
                    svar = svar.squeeze()
                if fill:
                    svar = np.where(svar=='', sfill_value, svar)
        else:
            if fill:
                var  = [ [ fill_value  if i=='' else i for i in row ] for row in var ]
            if squeeze or reform:
                maxi = max([ len(i) for i in var])
                if maxi==1: var = [ i[0] for i in var ]
            if transpose and isinstance(var[0], list):
                var = [list(i) for i in zip(*var)] # transpose
            if svar is not None:
                if fill:
                    svar = [ [ sfill_value if i=='' else i for i in row ] for row in svar ]
                if squeeze or reform:
                    maxi = max([ len(i) for i in svar])
                    if maxi==1: svar = [ i[0] for i in svar ]
                if transpose and isinstance(svar[0], list):
                    svar = [list(i) for i in zip(*svar)] # transpose
        if svar is not None:
            return var, svar
        else:
            return var
    #
    # Values - first line
    if (miianc >= nres) and (not fill):
        f.close()
        raise ValueError('Line has not enough columns to index: '+s)
    var  = list()
    svar = list()
    null = line2var(res, var, iinc, strip)
    null = line2var(res, svar, iisnc, False if strip is None else strip)
    #
    # Values - rest of file
    for line in f:
        s = str(line.rstrip())
        if len(s) == 0:
            if skip_blank:
                continue
            else:
                break
        if comment is not None:
            if (s[0] in comment): continue
        res = s.split(sep)
        nres = len(res)
        if (miianc >= nres) and (not fill):
            f.close()
            raise ValueError('Line has not enough columns to index: '+s)
        null = line2var(res, var, iinc, strip)
        null = line2var(res, svar, iisnc, False if strip is None else strip)
    f.close()
    # list -> array
    if fill:
        var = [ [ fill_value if i=='' else i for i in row ] for row in var ]
    # var  = np.array(var, dtype=np.float)
    # svar = np.array(svar, dtype=np.str)
    # if fill:
    #     svar = np.where(svar=='', sfill_value, svar)
    # if squeeze or reform:
    #     var  = var.squeeze()
    #     svar = svar.squeeze()
    # if transpose:
    #     var  = var.T
    #     svar = svar.T
    var  = np.array(var, dtype=np.float)
    if squeeze or reform: var  = var.squeeze()
    if transpose: var  = var.T
    if strarr:
        svar = np.array(svar, dtype=np.str)
        if transpose: svar = svar.T
        if squeeze or reform: svar = svar.squeeze()
        if fill: svar = np.where(svar=='', sfill_value, svar)
    else:
        if fill:
            svar = [ [ fill_value if i=='' else i for i in row ] for row in svar ]
        if squeeze or reform:
            maxi = max([ len(i) for i in svar])
            if maxi==1: svar = [ i[0] for i in svar ]
        if transpose and isinstance(svar[0], list):
            svar = [list(i) for i in zip(*svar)] # transpose

    return var, svar


# Helper for append var with current line already splitted into list
def line2var(res, var, iinc, strip):
    nres = len(res)
    if strip is None:
        tmp = [res[i].strip('"').strip("'") for i in iinc if i < nres]
    elif not strip:
        tmp = [res[i] for i in iinc if i < nres]
    else:
        tmp = [res[i].strip(strip) for i in iinc if i < nres]
    rest = len([ i for i in iinc if i >= nres ])
    if rest > 0: tmp.extend(['']*rest)
    var.append(tmp)
    return


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
