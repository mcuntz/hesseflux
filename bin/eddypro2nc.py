#!/usr/bin/env python
"""
usage: eddypro2nc.py [-h] [-o outfile] [eddypro_file]

Make netCDF file from EddyPro output file.

positional arguments:
  eddypro_file          Input EddyPro file

optional arguments:
  -h, --help            show this help message and exit
  -o outfile, --outfile outfile
                        Name of netCDF output file (default: input file
                        with suffix replaced by .nc)

Example
-------
python eddypro2nc.py -o eddypro_2020.nc data/2020/eddypro_full_output_2020.csv

History
-------
Written,  Matthias Cuntz, Feb 2021 - read EddyPro with pandas from
    https://git.icare.univ-lille1.fr/demo/eddypro_nc
"""
from __future__ import division, absolute_import, print_function
import datetime as dt
import os
import pandas as pd


# -------------------------------------------------------------------------
# Functions
#

def _get_columns_names(ifile, delimiter=','):
    """
    Extract names from second line of EddyPro output file `ifile`.
    """
    # read second line
    with open(ifile, "r", encoding='utf8') as f_id:
        f_id.readline()
        header = f_id.readline().strip()
    # list of names
    col_names = header.split(delimiter)
    # make unique names
    for i in col_names:
        ii = col_names.count(i)
        if col_names.count(i) > 1:
            for j in range(ii):
                i1 = col_names.index(i)
                col_names[i1] = col_names[i1] + str(j+1)
    # clean names
    col_names = [ uu.replace('-', '_').replace('*', 'star').replace('%', 'percent').replace('(', '').replace(')', '').replace('/', '_')
                  for uu in col_names ]

    return col_names


def _get_columns_units(ifile, delimiter=','):
    """
    Extract units from third line of EddyPro output file `ifile`.
    """
    # read third line
    with open(ifile, "r", encoding='utf8') as f_id:
        f_id.readline()
        f_id.readline()
        header = f_id.readline().strip()
    # list of units
    col_units = header.split(delimiter)
    # remove brackets from unites and make empty units
    col_units = [ uu.replace('[', '').replace(']', '').replace('--', '')
                  for uu in col_units ]

    return col_units


# -------------------------------------------------------------------------
# EddyPro output to netCDF file
#

if __name__ == '__main__':

    import argparse

    outfile = ''
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Make netCDF file from EddyPro output file.""")
    hstr  = 'Name of netCDF output file (default: input file with'
    hstr += ' suffix replaced by .nc)'
    parser.add_argument('-o', '--outfile', action='store',
                        default=outfile, dest='outfile', metavar='outfile',
                        help=hstr)
    hstr = 'Input EddyPro file'
    parser.add_argument('ifile', nargs='?', default=None,
                        metavar='eddypro_file', help=hstr)

    args = parser.parse_args()
    outfile = args.outfile
    ifile   = args.ifile

    del parser, args

    print('Read ', ifile)
    # Cleaned variables names and units
    col_names = _get_columns_names(ifile, delimiter=',')
    col_units = _get_columns_units(ifile, delimiter=',')

    # Read EddyPro
    date_parser = lambda str_: dt.datetime.strptime(str_, "%Y-%m-%d %H:%M")
    df = pd.read_csv(ifile,
                     skiprows=3,
                     delimiter=',',
                     na_values='-9999',
                     names=col_names,
                     parse_dates=[[1, 2]],
                     date_parser=date_parser,
                     index_col='date_time')
    df.drop(columns='filename', inplace=True)
    col_names = col_names[3:]
    col_units = col_units[3:]

    # Transform to xarray
    ds = df.to_xarray()

    # Add units
    for inn, nn in enumerate(col_names):
        ds[nn].attrs = {'units': col_units[inn]}

    # Add file attributes
    ds.attrs = {'Conventions': 'CF-1.6',
                'title': os.path.basename(ifile)}

    # Write netCDF file
    if not outfile:
        outfile = ifile[:ifile.rfind('.')] + '.nc'
    print('Write', outfile)
    ds.to_netcdf(outfile)
