#!/usr/bin/env python
"""
sread : Read strings into array from a file.

This module was written by Matthias Cuntz while at Department of
Computational Hydrosystems, Helmholtz Centre for Environmental
Research - UFZ, Leipzig, Germany, and continued while at Institut
National de Recherche pour l'Agriculture, l'Alimentation et
l'Environnement (INRAE), Nancy, France.

Copyright (c) 2009-2020 Matthias Cuntz - mc (at) macu (dot) de
Released under the MIT License; see LICENSE file for details.

* Written Jul 2009 by Matthias Cuntz (mc (at) macu (dot) de)
* Keyword transpose, Feb 2012, Matthias Cuntz
* Ported to Python 3, Feb 2013, Matthias Cuntz
* Removed bug when nc is list and contains 0, Nov 2014, Matthias Cuntz
* Keyword hskip, Nov 2014, Matthias Cuntz
* Do not use function lif, Feb 2015, Matthias Cuntz
* nc can be tuple, Feb 2015, Matthias Cuntz
* Large rewrite of code to improve speed, Feb 2015, Matthias Cuntz
* range instead of np.arange, Nov 2017, Matthias Cuntz
* Keywords cname, sname, hstrip, rename file to infile, Nov 2017, Matthias Cuntz
* Ignore unicode characters on read, Jun 2019, Matthias Cuntz
* Make ignoring unicode characters campatible with Python 2 and Python 3, Jul 2019, Matthias Cuntz
* Keywords encoding, errors with codecs module, Aug 2019, Matthias Cuntz
* Using numpy docstring format, May 2020, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz

The following functions are provided

.. autosummary::
   sread
"""
from __future__ import division, absolute_import, print_function
import numpy as np


__all__ = ['sread']


def sread(infile, nc=0, cname=None, skip=0, cskip=0, hskip=0, hstrip=True, separator=None,
          squeeze=False, reform=False, skip_blank=False, comment=None,
          fill=False, fill_value='', strip=None, encoding='ascii', errors='ignore',
          header=False, full_header=False,
          transpose=False, strarr=False):
    """
    Read strings into string array from a file.

    Lines or columns can be skipped.
    Columns can be picked specifically.

    Blank (only whitespace) and comment lines can be excluded.

    The header of the file can be read separately.

    This routines is exactly the same as fread but reads
    everything as strings instead of floats.

    Parameters
    ----------
    infile : str
        source file name
    nc : int or iterable, optional
        number of columns to be read [default: all (`nc<=0`)].

        `nc` can be an int or a vector of column indexes, starting with 0;
        `cskip` will be ignored in the latter case.
    cname : iterable of str, optional
        columns can be chosen by the values in the first header line;
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
        column separator. If not given, columns separators are (in order):
        comma (','), semicolon (';'), whitespace.
    comment : iterable, optional
        line gets excluded if first character of line is in comment sequence.
        Sequence must be iterable such as string, list and tuple.
    fill_value : float, optional
        value to fill in array in empty cells or if not enough columns in line
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

        False: array will be two-dimensional as read (default)
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
        True:  return as numpy array of strings.

        False: return as list (default).

    Returns
    -------
    array of str
        Depending on options:

        List of strings if `header==False` and `strarr==False`.

        Array of strings if `header==False` and `strarr==True`.

        List with file header strings if `header==True` and `strarr==False`.

        String array of file header if `header==True` and `strarr==True`.

        List of lines of strings if `header=True` and `full_header=True`.

    Notes
    -----
    If `header=True` then skip is counterintuitive because it is
    actually the number of header rows to be read. This is to
    be able to have the exact same call of the function, once
    with `header=False` and once with `header=True`.

    If `fill==True`, blank lines are not filled but are takes as end of file.

    `transpose=True` has no effect on 1D output such as 1 header line.

    Examples
    --------
    >>> # Create some data
    >>> filename = 'test.dat'
    >>> ff = open(filename, 'w')
    >>> ff.writelines('head1 head2 head3 head4\\n')
    >>> ff.writelines('1.1 1.2 1.3 1.4\\n')
    >>> ff.writelines('2.1 2.2 2.3 2.4\\n')
    >>> ff.close()

    >>> # Read sample file in different ways
    >>> # header
    >>> print(sread(filename, nc=2, skip=1, header=True))
    ['head1', 'head2']
    >>> print(sread(filename, nc=2, skip=1, header=True, full_header=True))
    ['head1 head2 head3 head4']
    >>> print(sread(filename, nc=1, skip=2, header=True))
    [['head1'], ['1.1']]
    >>> print(sread(filename, nc=1, skip=2, header=True, squeeze=True))
    ['head1', '1.1']
    >>> print(sread(filename, nc=1, skip=2, header=True, squeeze=True, strarr=True))
    ['head1' '1.1']
    >>> print(sread(filename, nc=1, skip=2, header=True, squeeze=True, transpose=True))
    ['head1', '1.1']

    >>> # data
    >>> print(sread(filename, skip=1))
    [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4']]
    >>> print(sread(filename, skip=2))
    [['2.1', '2.2', '2.3', '2.4']]
    >>> print(sread(filename, skip=1, cskip=1))
    [['1.2', '1.3', '1.4'], ['2.2', '2.3', '2.4']]
    >>> print(sread(filename, nc=2, skip=1, cskip=1))
    [['1.2', '1.3'], ['2.2', '2.3']]
    >>> print(sread(filename, nc=[1,3], skip=1))
    [['1.2', '1.4'], ['2.2', '2.4']]
    >>> print(sread(filename, nc=1, skip=1))
    [['1.1'], ['2.1']]
    >>> print(sread(filename, nc=1, skip=1, reform=True))
    ['1.1', '2.1']

    >>> # skip blank lines
    >>> ff = open(filename, 'a')
    >>> ff.writelines('\\n')
    >>> ff.writelines('3.1 3.2 3.3 3.4\\n')
    >>> ff.close()
    >>> print(sread(filename, skip=1))
    [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4']]
    >>> print(sread(filename, skip=1, skip_blank=True))
    [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4'], ['3.1', '3.2', '3.3', '3.4']]
    >>> print(sread(filename, skip=1, strarr=True))
    [['1.1' '1.2' '1.3' '1.4']
     ['2.1' '2.2' '2.3' '2.4']]

    >>> print(sread(filename, skip=1, strarr=True, transpose=True))
    [['1.1' '2.1']
     ['1.2' '2.2']
     ['1.3' '2.3']
     ['1.4' '2.4']]
    >>> print(sread(filename, skip=1, transpose=True))
    [['1.1', '2.1'], ['1.2', '2.2'], ['1.3', '2.3'], ['1.4', '2.4']]

    >>> # skip comment lines
    >>> ff = open(filename, 'a')
    >>> ff.writelines('# First\\n')
    >>> ff.writelines('! Second second comment\\n')
    >>> ff.writelines('4.1 4.2 4.3 4.4\\n')
    >>> ff.close()
    >>> print(sread(filename, skip=1))
    [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4']]
    >>> print(sread(filename, skip=1, skip_blank=True, comment='#'))
    [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4'], ['3.1', '3.2', '3.3', '3.4'], ['!', 'Second', 'second', 'comment'], ['4.1', '4.2', '4.3', '4.4']]
    >>> print(sread(filename, skip=1, skip_blank=True, comment='#!'))
    [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4'], ['3.1', '3.2', '3.3', '3.4'], ['4.1', '4.2', '4.3', '4.4']]
    >>> print(sread(filename, skip=1, skip_blank=True, comment=('#','!')))
    [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4'], ['3.1', '3.2', '3.3', '3.4'], ['4.1', '4.2', '4.3', '4.4']]
    >>> print(sread(filename, skip=1, skip_blank=True, comment=['#','!']))
    [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4'], ['3.1', '3.2', '3.3', '3.4'], ['4.1', '4.2', '4.3', '4.4']]

    >>> # Create some more data with escaped numbers
    >>> filename2 = 'test2.dat'
    >>> ff = open(filename2, 'w')
    >>> ff.writelines('"head1" "head2" "head3" "head4"\\n')
    >>> ff.writelines('"1.1" "1.2" "1.3" "1.4"\\n')
    >>> ff.writelines('2.1 nan Inf "NaN"\\n')
    >>> ff.close()
    >>> print(sread(filename2, skip=1, strarr=True, transpose=True, strip='"'))
    [['1.1' '2.1']
     ['1.2' 'nan']
     ['1.3' 'Inf']
     ['1.4' 'NaN']]

    >>> # Create some more data with an extra (shorter) header line
    >>> filename3 = 'test3.dat'
    >>> ff = open(filename3, 'w')
    >>> ff.writelines('Extra header\\n')
    >>> ff.writelines('head1 head2 head3 head4\\n')
    >>> ff.writelines('1.1 1.2 1.3 1.4\\n')
    >>> ff.writelines('2.1 2.2 2.3 2.4\\n')
    >>> ff.close()

    >>> print(sread(filename3, skip=2))
    [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4']]
    >>> print(sread(filename3, skip=2, hskip=1))
    [['1.1', '1.2', '1.3', '1.4'], ['2.1', '2.2', '2.3', '2.4']]
    >>> print(sread(filename3, nc=2, skip=2, hskip=1, header=True))
    ['head1', 'head2']

    >>> # Create some more data with missing values
    >>> filename4 = 'test4.dat'
    >>> ff = open(filename4, 'w')
    >>> ff.writelines('Extra header\\n')
    >>> ff.writelines('head1,head2,head3,head4\\n')
    >>> ff.writelines('1.1,1.2,1.3,1.4\\n')
    >>> ff.writelines('2.1,,2.3,2.4\\n')
    >>> ff.close()

    >>> print(sread(filename4, skip=2))
    [['1.1', '1.2', '1.3', '1.4'], ['2.1', '', '2.3', '2.4']]
    >>> print(sread(filename4, skip=2, fill=True, fill_value='-1'))
    [['1.1', '1.2', '1.3', '1.4'], ['2.1', '-1', '2.3', '2.4']]

    >>> # cname
    >>> print(sread(filename, cname='head2', skip=1, skip_blank=True, comment='#!', squeeze=True))
    ['1.2', '2.2', '3.2', '4.2']
    >>> print(sread(filename, cname=['head1','head2'], skip=1, skip_blank=True, comment='#!'))
    [['1.1', '1.2'], ['2.1', '2.2'], ['3.1', '3.2'], ['4.1', '4.2']]
    >>> print(sread(filename, cname=['head1','head2'], skip=1, skip_blank=True, comment='#!', header=True))
    ['head1', 'head2']
    >>> print(sread(filename, cname=['head1','head2'], skip=1, skip_blank=True, comment='#!', header=True, full_header=True))
    ['head1 head2 head3 head4']
    >>> print(sread(filename, cname=['  head1','head2'], skip=1, skip_blank=True, comment='#!', hstrip=False))
    [['1.2'], ['2.2'], ['3.2'], ['4.2']]

    >>> # Clean up doctest
    >>> import os
    >>> os.remove(filename)
    >>> os.remove(filename2)
    >>> os.remove(filename3)
    >>> os.remove(filename4)

    History
    -------
    Written,  Matthias Cuntz, Jul 2009
    Modified, Matthias Cuntz, Feb 2012 - transpose
              Matthias Cuntz, Feb 2013 - ported to Python 3
              Matthias Cuntz, Nov 2014 - bug when nc is list and contains 0
              Matthias Cuntz, Nov 2014 - hskip
              Matthias Cuntz, Feb 2015 - no lif, nc can be tuple
                                       - large re-write
              Matthias Cuntz, Nov 2017 - use range instead of np.arange for producing indexes
              Matthias Cuntz, Nov 2017 - cname, sname, file->infile, hstrip
              Matthias Cuntz, Jun 2019 - open(errors='ignore') to ignore unicode characters, for example, on read
              Matthias Cuntz, Jul 2019 - errors='ignore' compatible with Python2 and Python3
                                         -> returns header in unicode in Python2
              Matthias Cuntz, Aug 2019 - use codecs module and allow user encoding and error handling
    """
    #
    # Open file
    import codecs
    f = codecs.open(infile, 'r', encoding=encoding, errors=errors)
    #
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
    #
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
    #
    # Determine indices
    if nc != 0 and cname is not None:
        f.close()
        raise ValueError('nc and cname are mutually exclusive.')
    if cname is not None:
        # from first header line
        if (skip-hskip) <= 0:
            f.close()
            raise IOError('No header line left for choosing columns by name.')
        if not isinstance(cname, (list, tuple, np.ndarray)): cname = [cname]
        if hstrip: cname = [ h.strip() for h in cname ]
        hres = head[0].split(sep)
        if hstrip: hres = [ h.strip() for h in hres ]
        iinc = []
        for k in range(len(hres)):
            if hres[k] in cname: iinc.append(k)
        nnc = len(iinc)
    else:
        # from nc keyword
        if isinstance(nc, (list, tuple, np.ndarray)):
            nnc  = len(nc)
            iinc = tuple(nc)
        else:
            if nc <= 0:
                iinc = range(cskip,nres)
                nnc  = nres-cskip
            else:
                iinc = range(cskip,cskip+nc)
                nnc = nc
    miinc = max(iinc)
    #
    # Header
    if header:
        # Split header
        var = None
        if (skip-hskip) > 0:
            if full_header:
                var = head
            else:
                var = list()
                k = 0
                while k < (skip-hskip):
                    hres = head[k].split(sep)
                    nhres = len(hres)
                    if (miinc >= nhres) and (not fill):
                        f.close()
                        raise ValueError('Line has not enough columns to index: '+head[k])
                    null = line2var(hres, var, iinc, strip)
                    k += 1
                if (skip-hskip) == 1: var = var[0]
        f.close()
        if strarr:
            var = np.array(var, dtype=np.str)
            if transpose: var = var.T
            if squeeze or reform: var = var.squeeze()
            if fill: var = np.where(var=='', fill_value, var)
        else:
            if fill:
                var = [ [ fill_value if i=='' else i for i in row ] for row in var ]
            if squeeze or reform:
                maxi = max([ len(i) for i in var])
                if maxi==1: var = [ i[0] for i in var ]
            if transpose and isinstance(var[0], list):
                var = [list(i) for i in zip(*var)] # transpose
        return var
    #
    # Values - first line
    if (miinc >= nres) and (not fill):
        f.close()
        raise ValueError('Line has not enough columns to index: '+s)
    var = list()
    null = line2var(res, var, iinc, strip)
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
        if (miinc >= nres) and (not fill):
            f.close()
            raise ValueError('Line has not enough columns to index: '+s)
        null = line2var(res, var, iinc, strip)
    f.close()
    if strarr:
        var = np.array(var, dtype=np.str)
        if transpose: var = var.T
        if squeeze or reform: var = var.squeeze()
        if fill: var = np.where(var=='', fill_value, var)
    else:
        if fill:
            var = [ [ fill_value if i=='' else i for i in row ] for row in var ]
        if squeeze or reform:
            maxi = max([ len(i) for i in var])
            if maxi==1: var = [ i[0] for i in var ]
        if transpose and isinstance(var[0], list):
            var = [list(i) for i in zip(*var)] # transpose

    return var


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
    if rest > 0:
        tmp.extend(['']*rest)
    var.append(tmp)
    return


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
