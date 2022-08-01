#!/usr/bin/env python
"""
This is a dummy unittest.

python -m unittest -v test_dummy.py

"""
import unittest


class TestDummy(unittest.TestCase):

    def test_dummy(self):
        import hesseflux as hf


if __name__ == "__main__":
    unittest.main()
