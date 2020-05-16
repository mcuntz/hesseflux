#!/usr/bin/env python
"""
Takes file in convention of europe-fluxdata.eu,
which does not include unit information for the variables,
and adds a second header line with the units of the variables.

It uses the file 'europe-fluxdata_units.csv' with information on units,
supposed to be found in the run directory, where the script was started,
or the directory of this Python script.

Unknown variables will have no unit.

Output filename is the input filename with -units appended before the file suffix, e.g.
infile.csv -> infile-units.csv

If no input filename is given, a gui opens where you can chose an input file.


Example
-------
    python europe-fluxdata_units.py ../europe-fluxdata/FR-Hes_europe-fluxdata_????.csv


History
-------
Written, Matthias Cuntz, March 2020
"""
from __future__ import division, absolute_import, print_function
import sys
import os
import numpy as np
import hesseflux as hf

import time as ptime
t1 = ptime.time()

# units file
unitfile='europe-fluxdata_units.csv'

# -----------------------------------------------------------------------------
# Main

if __name__ == '__main__':

    curdir = os.path.dirname(sys.argv[0])
    if curdir == '': curdir='.'
    # no filename(s) given
    if len(sys.argv) == 1:
        try:
            infile = hf.file_from_gui(initialdir=curdir, title='Open europe-fluxdata.eu file.')
        except:
            raise IOError("File Gui for europe-fluxdata.eu file failed.")

        if not infile:
            raise IOError("No europe-fluxdata.eu file given.")

        if isinstance(infile, (list, tuple, np.ndarray)):
            infiles = infile
        else:
            infiles = [infile]
    # one filename given
    elif len(sys.argv) == 2:
        infiles = [sys.argv[1]]
    # several filenames given
    else:
        infiles = sys.argv[1:]

    # Read units file
    # Variable_Code;Description;Unit;Template_to_be_used
    # ALB;Albedo;%;Database
    uf = curdir+'/'+unitfile
    if not os.path.exists(uf): uf = './'+unitfile
    assert os.path.exists(uf), 'Unit file not in current folder nor in folder of Python script.'
    udat = hf.sread(uf, skip=1, nc=[0,2], separator=';', strarr=True)
    uvars  = list(udat[:,0])
    uunits = list(udat[:,1])

    # for all input files
    for ff in infiles:
        print(ff)
        # header
        hdat  = hf.sread(ff, skip=1, header=True)
        # Get units
        hunit = hdat.copy()
        for ii, hh in enumerate(hdat):
            ivar = hh.rsplit('_', maxsplit=3)[0]
            if ivar in uvars:
                hunit[ii] = uunits[uvars.index(ivar)]
            else:
                hunit[ii] = ''
        # output file name
        ofile = ff.rsplit('.', maxsplit=1)
        ofile = ofile[0]+'-units.'+ofile[1]
        # read / write
        iff = open(ff, 'r')
        off = open(ofile, 'w')
        # header line
        hh = iff.readline()
        print(hh, end='', file=off)
        # units line
        uu = ','.join(hunit)
        print(uu, file=off)
        # rest of file
        for line in iff:
            print(line, end='', file=off)
        iff.close()
        off.close()

    # -------------------------------------------------------------------------
    # Finish

    t2    = ptime.time()
    strin = '[m]: {:.1f}'.format((t2-t1)/60.) if (t2-t1)>60. else '[s]: {:d}'.format(int(t2-t1))
    print('Time elapsed', strin)
