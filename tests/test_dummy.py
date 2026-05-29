#!/usr/bin/env python
"""
This is a dummy unittest.

python -m unittest -v test_dummy.py

"""
from pathlib import Path
import unittest

import pandas as pd

import sys
sys.path.insert(0, str((Path(__file__).parent.parent / "src").resolve()))
import hesseflux as hf

def load_data():
    fullout_file = Path(__file__).parent / "data" / "full_output.csv"
    fullout = pd.read_csv(fullout_file, skiprows=[0, 2], na_values=["-9999.00", "-9999"])
    # with open(fullout_file, "r") as f:
    #     f.readline()
    #     f.readline()
    #     units = f.readline().strip().split(",")
    #     units = {c:u.strip('[]') for c, u in zip(fullout.columns, units)}
    fullout["TIMESTAMP"] = pd.to_datetime(fullout["date"] + " " + fullout["time"])
    fullout = fullout.set_index("TIMESTAMP").sort_index()

    df = (
        fullout
        [["u*", "co2_flux", "air_temperature"]]
        .rename(columns={"u*":"USTAR", "co2_flux":"NEE", "air_temperature":"TA", "qc_co2_flux":"flag", "daytime":"isday"})
        .assign(TA = lambda df: df["TA"] - 273.15,)
    )

    df_flag = fullout[["qc_Tau", "qc_co2_flux"]].rename(columns={"qc_Tau":"USTAR", "qc_co2_flux":"NEE"})
    df_flag["TA"] = 0

    df_isday = fullout[["daytime"]].rename(columns={"daytime":"isday"}).astype(bool)

    return df, df_flag, df_isday

class TestDummy(unittest.TestCase):

    def test_dummy(self):
        import hesseflux as hf

    def test_ustarfilter(self):
        
        df, df_flag, df_isday = load_data()
        hf.ustarfilter(
            df, 
            flag=df_flag,
            isday=df_isday
        )




if __name__ == "__main__":
    df, df_flag, df_isday = load_data()
    hf.ustarfilter(df, flag=df_flag, isday=df_isday)
    # unittest.main()
