#!/usr/bin/env python
"""
This is a dummy unittest.

python -m unittest -v test_dummy.py

"""
from pathlib import Path
import unittest

import pandas as pd
import numpy as np

import sys
sys.path.insert(0, str((Path(__file__).parent.parent / "src").resolve()))
import hesseflux as hf

def load_data():
    fullout_file = Path(__file__).parent / "data" / "full_output.csv"
    fullout = pd.read_csv(fullout_file, skiprows=[0, 2], na_values=["-9999.00", "-9999"])
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
        out = hf.ustarfilter(df, flag=df_flag, isday=df_isday, minseasondata=20, mindaysperyear=25)
        self.assertEqual(out[1].sum(), 990, "ustarfilter did properly filter the data")
        self.assertTrue(np.all(np.isclose(out[0], 0.173353)), "ustarfilter did not return expected value")

        df, df_flag, df_isday = load_data()
        out = hf.ustarfilter(df, flag=df_flag, isday=df_isday, seasonout=True, nmon=1, nboot=10, ustardefault=0.5, randomstate=8472)
        print(out[1].sum())
        print(out[1].shape)
        self.assertTrue(
            np.allclose(
                out[0], 
                np.array([[0.5, 0.5, 0.2749906 , 0.163295  , 0.5,0.5, 0.5, 0.5, 0.5, 0.5,0.5, 0.5],
                        [0.5, 0.5, 0.5, 0.1807536 , 0.5,0.5, 0.5, 0.5, 0.5, 0.5,0.5, 0.5],
                        [0.5, 0.5, 0.5, 0.59664678, 0.5,0.5, 0.5, 0.5, 0.5, 0.5,0.5, 0.5]])
            ),
            "ustarfilter did not return expected ustar threshold values with seasonout=True"
        )
        self.assertEqual(out[1].sum(), 1304, "ustarfilter did properly filter the data with seasonout=True")
        



if __name__ == "__main__":
    unittest.main()
