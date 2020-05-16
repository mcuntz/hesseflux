#!/usr/bin/env python
"""
fread : Read numbers into array from a file.

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
* Speed improvements: everything list until the very end, Feb 2015, Matthias Cuntz
* range instead of np.arange, Nov 2017, Matthias Cuntz
* Keywords cname, sname, hstrip, rename file to infile, Nov 2017, Matthias Cuntz
* Ignore unicode characters on read, Jun 2019, Matthias Cuntz
* Make ignoring unicode characters campatible with Python 2 and Python 3, Jul 2019, Matthias Cuntz
* Keywords encoding, errors with codecs module, Aug 2019, Matthias Cuntz
* Keyword return_list, Dec 2019, Stephan Thober
* return_list=False default, Jan 2020, Matthias Cuntz
* Using numpy docstring format, May 2020, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz

The following functions are provided

.. autosummary::
   fread
"""
from __future__ import division, absolute_import, print_function
import numpy as np


__all__ = ['fread']


def fread(infile, nc=0, cname=None, skip=0, cskip=0, hskip=0, hstrip=True, separator=None,
          squeeze=False, reform=False, skip_blank=False, comment=None,
          fill=False, fill_value=0, strip=None, encoding='ascii', errors='ignore',
          header=False, full_header=False,
          transpose=False, strarr=False, return_list=False):
    """
    Read numbers into array with floats from a file.

    Lines or columns can be skipped.
    Columns can be picked specifically.

    Blank (only whitespace) and comment lines can be excluded.

    The header of the file can be read separately.

    This routines is exactly the same as sread but transforms
    everything to floats, handling NaN and Inf.

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
         and `fill==True` (default: 0, and '' for header).
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
        True:  return header as numpy array of strings.

        False: return header as list (default).
    return_list : bool, optional
        True:  return file content as list.

        False: return file content as numpy array (default).

    Returns
    -------
    array of floats
        Depending on options:

        Array of floats if `header==False`.

        List of floats if `return_list==True`.

        List with file header strings if `header==True`.

        String array of file header if `header==True` and `strarr==True`.

        List of lines of strings if `header=True` and `full_header=True`.

    Notes
    -----
    If `header==True` then skip is counterintuitive because it is
    actually the number of header rows to be read. This is to
    be able to have the exact same call of the function, once
    with `header=False` and once with `header=True`.

    If `fill==True`, blank lines are not filled but are taken as end of file.

    `transpose=True` has no effect on 1D output such as 1 header line.

    Examples
    --------
    >>> # Create some data
    >>> filename = 'test.dat'
    >>> ff = open(filename,'w')
    >>> ff.writelines('head1 head2 head3 head4\\n')
    >>> ff.writelines('1.1 1.2 1.3 1.4\\n')
    >>> ff.writelines('2.1 2.2 2.3 2.4\\n')
    >>> ff.close()

    >>> # Read sample file in different ways
    >>> # header
    >>> print(fread(filename, nc=2, skip=1, header=True))
    ['head1', 'head2']
    >>> print(fread(filename, nc=2, skip=1, header=True, full_header=True))
    ['head1 head2 head3 head4']
    >>> print(fread(filename, nc=1, skip=2, header=True))
    [['head1'], ['1.1']]
    >>> print(fread(filename, nc=1, skip=2, header=True, squeeze=True))
    ['head1', '1.1']
    >>> print(fread(filename, nc=1, skip=2, header=True, strarr=True))
    [['head1']
     ['1.1']]

    >>> # data
    >>> print(fread(filename, skip=1))
    [[1.1 1.2 1.3 1.4]
     [2.1 2.2 2.3 2.4]]
    >>> print(fread(filename, skip=2))
    [[2.1 2.2 2.3 2.4]]
    >>> print(fread(filename, skip=1, cskip=1))
    [[1.2 1.3 1.4]
     [2.2 2.3 2.4]]
    >>> print(fread(filename, nc=2, skip=1, cskip=1))
    [[1.2 1.3]
     [2.2 2.3]]
    >>> print(fread(filename, nc=[1,3], skip=1))
    [[1.2 1.4]
     [2.2 2.4]]
    >>> print(fread(filename, nc=1, skip=1))
    [[1.1]
     [2.1]]
    >>> print(fread(filename, nc=1, skip=1, reform=True))
    [1.1 2.1]

    >>> # skip blank lines
    >>> ff = open(filename, 'a')
    >>> ff.writelines('\\n')
    >>> ff.writelines('3.1 3.2 3.3 3.4\\n')
    >>> ff.close()
    >>> print(fread(filename, skip=1))
    [[1.1 1.2 1.3 1.4]
     [2.1 2.2 2.3 2.4]]
    >>> print(fread(filename, skip=1, skip_blank=True, comment='#!'))
    [[1.1 1.2 1.3 1.4]
     [2.1 2.2 2.3 2.4]
     [3.1 3.2 3.3 3.4]]

    >>> # skip comment lines
    >>> ff = open(filename, 'a')
    >>> ff.writelines('# First comment\\n')
    >>> ff.writelines('! Second 2 comment\\n')
    >>> ff.writelines('4.1 4.2 4.3 4.4\\n')
    >>> ff.close()
    >>> print(fread(filename, skip=1))
    [[1.1 1.2 1.3 1.4]
     [2.1 2.2 2.3 2.4]]
    >>> print(fread(filename, skip=1, nc=[2], skip_blank=True, comment='#'))
    [[1.3]
     [2.3]
     [3.3]
     [2. ]
     [4.3]]
    >>> print(fread(filename, skip=1, skip_blank=True, comment='#!'))
    [[1.1 1.2 1.3 1.4]
     [2.1 2.2 2.3 2.4]
     [3.1 3.2 3.3 3.4]
     [4.1 4.2 4.3 4.4]]
    >>> print(fread(filename, skip=1, skip_blank=True, comment=('#','!')))
    [[1.1 1.2 1.3 1.4]
     [2.1 2.2 2.3 2.4]
     [3.1 3.2 3.3 3.4]
     [4.1 4.2 4.3 4.4]]
    >>> print(fread(filename, skip=1, skip_blank=True, comment=['#','!']))
    [[1.1 1.2 1.3 1.4]
     [2.1 2.2 2.3 2.4]
     [3.1 3.2 3.3 3.4]
     [4.1 4.2 4.3 4.4]]

    >>> # fill missing columns
    >>> ff = open(filename, 'a')
    >>> ff.writelines('5.1 5.2\\n')
    >>> ff.close()
    >>> print(fread(filename, skip=1))
    [[1.1 1.2 1.3 1.4]
     [2.1 2.2 2.3 2.4]]
    >>> print(fread(filename, skip=1, skip_blank=True, comment='#!', fill=True, fill_value=-1))
    [[ 1.1  1.2  1.3  1.4]
     [ 2.1  2.2  2.3  2.4]
     [ 3.1  3.2  3.3  3.4]
     [ 4.1  4.2  4.3  4.4]
     [ 5.1  5.2 -1.  -1. ]]

    >>> # transpose
    >>> print(fread(filename, skip=1))
    [[1.1 1.2 1.3 1.4]
     [2.1 2.2 2.3 2.4]]
    >>> print(fread(filename, skip=1, transpose=True))
    [[1.1 2.1]
     [1.2 2.2]
     [1.3 2.3]
     [1.4 2.4]]

    >>> # Create some more data with Nan and Inf
    >>> filename1 = 'test1.dat'
    >>> ff = open(filename1, 'w')
    >>> ff.writelines('head1 head2 head3 head4\\n')
    >>> ff.writelines('1.1 1.2 1.3 1.4\\n')
    >>> ff.writelines('2.1 nan Inf "NaN"\\n')
    >>> ff.close()

    >>> # Treat Nan and Inf with automatic strip of " and '
    >>> print(fread(filename1, skip=1, transpose=True))
    [[1.1 2.1]
     [1.2 nan]
     [1.3 inf]
     [1.4 nan]]

    >>> # Create some more data with escaped numbers
    >>> filename2 = 'test2.dat'
    >>> ff = open(filename2, 'w')
    >>> ff.writelines('head1 head2 head3 head4\\n')
    >>> ff.writelines('"1.1" "1.2" "1.3" "1.4"\\n')
    >>> ff.writelines('2.1 nan Inf "NaN"\\n')
    >>> ff.close()

    >>> # Strip
    >>> print(fread(filename2,  skip=1,  transpose=True,  strip='"'))
    [[1.1 2.1]
     [1.2 nan]
     [1.3 inf]
     [1.4 nan]]

    >>> # Create some more data with an extra (shorter) header line
    >>> filename3 = 'test3.dat'
    >>> ff = open(filename3, 'w')
    >>> ff.writelines('Extra header\\n')
    >>> ff.writelines('head1 head2 head3 head4\\n')
    >>> ff.writelines('1.1 1.2 1.3 1.4\\n')
    >>> ff.writelines('2.1 2.2 2.3 2.4\\n')
    >>> ff.close()

    >>> print(fread(filename3, skip=2, hskip=1))
    [[1.1 1.2 1.3 1.4]
     [2.1 2.2 2.3 2.4]]
    >>> print(fread(filename3, nc=2, skip=2, hskip=1, header=True))
    ['head1', 'head2']

    >>> # cname
    >>> print(fread(filename, cname='head2', skip=1, skip_blank=True, comment='#!', squeeze=True))
    [1.2 2.2 3.2 4.2 5.2]
    >>> print(fread(filename, cname=['head1','head2'], skip=1, skip_blank=True, comment='#!'))
    [[1.1 1.2]
     [2.1 2.2]
     [3.1 3.2]
     [4.1 4.2]
     [5.1 5.2]]
    >>> print(fread(filename, cname=['head1','head2'], skip=1, skip_blank=True, comment='#!', header=True))
    ['head1', 'head2']
    >>> print(fread(filename, cname=['head1','head2'], skip=1, skip_blank=True, comment='#!', header=True, full_header=True))
    ['head1 head2 head3 head4']
    >>> print(fread(filename, cname=['  head1','head2'], skip=1, skip_blank=True, comment='#!', hstrip=False))
    [[1.2]
     [2.2]
     [3.2]
     [4.2]
     [5.2]]

    >>> # Clean up doctest
    >>> import os
    >>> os.remove(filename)
    >>> os.remove(filename1)
    >>> os.remove(filename2)
    >>> os.remove(filename3)

    History
    -------
    Written,  Matthias Cuntz, Jul 2009
    Modified, Matthias Cuntz, Feb 2012 - transpose
              Matthias Cuntz, Feb 2013 - ported to Python 3
              Matthias Cuntz, Nov 2014 - bug when nc is list and contains 0
              Matthias Cuntz, Nov 2014 - hskip
              Matthias Cuntz, Feb 2015 - speed: everything list until very end
              Matthias Cuntz, Nov 2017 - use range instead of np.arange for producing indexes
              Matthias Cuntz, Nov 2017 - cname, sname, file->infile, hstrip
              Matthias Cuntz, Jun 2019 - open(errors='ignore') to ignore unicode characters, for example, on read
              Matthias Cuntz, Jul 2019 - errors='ignore' compatible with Python2 and Python3
                                         -> returns header in unicode in Python2
              Matthias Cuntz, Aug 2019 - use codecs module and allow user encoding and error handling
              Stephan Thober, Dec 2019 - added return_list flag
              Matthias Cuntz, Jan 2020 - default return_list=False
              Matthias Cuntz, May 2020 - numpy docstring format
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
        s = f.readline().rstrip()
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
                        raise ValueError('Line has not enough columns to index - 01: '+head[k])
                    null = line2var(hres, var, iinc, strip)
                    k += 1
                if (skip-hskip) == 1: var = var[0]
            f.close()
        else:
            var = None
            f.close()
            return var
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
                var = [ list(i) for i in zip(*var) ] # transpose
        return var
    #
    # Values - first line
    if (miinc >= nres) and (not fill):
        f.close()
        raise ValueError('Line has not enough columns to index - 02: '+s)
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
            raise ValueError('Line has not enough columns to index - 03: '+s)
        null = line2var(res, var, iinc, strip)

    f.close()
    # list -> array
    if fill:
        var = [ [ fill_value if i=='' else i for i in row ] for row in var ]
    var = np.array(var, dtype=np.float)
    if squeeze or reform: var = var.squeeze()
    if transpose: var = var.T
    if return_list:
        if var.ndim == 1:
            var = [ i for i in var ]
        else:
            var = [ [ var[i,j] for j in range(var.shape[1]) ] for i in range(var.shape[0]) ]

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
