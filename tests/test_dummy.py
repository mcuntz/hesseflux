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
        [["u*", "co2_flux", "air_temperature", "qc_co2_flux", "daytime"]]
        .dropna()
        .rename(columns={"u*":"USTAR", "co2_flux":"NEE", "air_temperature":"TA", "qc_co2_flux":"flag", "daytime":"isday"})
        .assign(
            TA = lambda df: df["TA"] - 273.15,
            flag = lambda df: df["flag"].astype(float),
            isday = lambda df: df["isday"].astype(float)
        )
    )

    return df

class TestDummy(unittest.TestCase):

    def test_dummy(self):
        import hesseflux as hf

    def test_ustarfilter(self):
        
        df = load_data()
        hf.ustarfilter(
            df[["USTAR", "NEE", "TA"]], 
            flag=df[["flag"]],
            isday=df[["isday"]]
        )




if __name__ == "__main__":
    unittest.main()
