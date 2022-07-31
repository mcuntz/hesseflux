#!/usr/bin/env python
"""
    This is the unittest for argsort module.

    python -m unittest -v test_argsort.py
"""
from __future__ import division, absolute_import, print_function
import unittest


# --------------------------------------------------------------------
# argsort.py
class TestArgsort(unittest.TestCase):

    def test_argmax(self):
        import os
        import numpy as np
        from hesseflux import argmax

        # One-dimensional array
        a = np.array([0,4,6,2,1,5,3,5])
        ii = argmax(a)
        self.assertEqual(ii, 2)
        self.assertEqual(a[ii], 6)

        # One-dimensional masked array
        a = np.ma.array([0,4,6,2,1,5,3,5], mask=[0,0,1,1,0,0,0,0])
        ii = argmax(a)
        self.assertEqual(ii, 5)
        self.assertEqual(a[ii], 5)
        ii = argmax(a, fill_value=6)
        self.assertEqual(ii, 2)

        # List
        a = [0,4,6,2,1,5,3,5]
        ii = argmax(a)
        self.assertEqual(ii, 2)
        self.assertEqual(a[ii], 6)

        # from numpy.argmax docstring
        a = np.arange(6).reshape(2,3) + 10
        ii = np.argmax(a)
        self.assertEqual(ii, 5)
        ii = np.argmax(a, axis=0)
        self.assertEqual(list(ii), [1, 1, 1])
        ii = np.argmax(a, axis=1)
        self.assertEqual(list(ii), [2, 2])

        # Indexes of the maximal elements of a N-dimensional array
        ii = np.unravel_index(np.argmax(a, axis=None), a.shape)
        self.assertEqual(list(ii), [1, 2])
        self.assertEqual(a[ii], 15)

        b = np.arange(6)
        b[1] = 5
        ii = np.argmax(b)  # Only the first occurrence is returned.
        self.assertEqual(ii, 1)


    def test_argmin(self):
        import os
        import numpy as np
        from hesseflux import argmin

        # One-dimensional array
        a = np.array([0,4,6,2,1,5,3,5])
        ii = argmin(a)
        self.assertEqual(ii, 0)
        self.assertEqual(a[ii], 0)

        # One-dimensional masked array
        a = np.ma.array([0,4,6,2,1,5,3,5], mask=[1,0,1,1,0,0,0,0])
        ii = argmin(a)
        self.assertEqual(ii, 4)
        self.assertEqual(a[ii], 1)
        ii = argmin(a, fill_value=1)
        self.assertEqual(ii, 0)

        # List
        a = [0,4,6,2,1,5,3,5]
        ii = argmin(a)
        self.assertEqual(ii, 0)
        self.assertEqual(a[ii], 0)

        # from numpy.argmin docstring
        a = np.arange(6).reshape(2,3) + 10
        ii = np.argmin(a)
        self.assertEqual(ii, 0)
        ii = np.argmin(a, axis=0)
        self.assertEqual(list(ii), [0, 0, 0])
        ii = np.argmin(a, axis=1)
        self.assertEqual(list(ii), [0, 0])

        # Indices of the minimum elements of a N-dimensional array:
        ii = np.unravel_index(np.argmin(a, axis=None), a.shape)
        self.assertEqual(list(ii), [0, 0])
        self.assertEqual(a[ii], 10)

        b = np.arange(6) + 10
        b[4] = 10
        ii = np.argmin(b)  # Only the first occurrence is returned.
        self.assertEqual(ii, 0)


    def test_argsort(self):
        import os
        import numpy as np
        from hesseflux import argsort

        # 1D array
        a = np.array([0,4,6,2,1,5,3,5])
        ii = argsort(a)
        self.assertEqual(list(a[ii]), [0, 1, 2, 3, 4, 5, 5, 6])
        ii = argsort(a, kind='quicksort')
        self.assertEqual(list(a[ii]), [0, 1, 2, 3, 4, 5, 5, 6])

        # 1D masked array
        a = np.ma.array([0,4,6,2,1,5,3,5], mask=[0,0,1,1,0,0,0,0])
        ii = argsort(a)
        self.assertEqual(list(a[ii]), [0, 1, 3, 4, 5, 5, np.ma.masked, np.ma.masked])
        ii = argsort(a, fill_value=1)
        self.assertEqual(list(a[ii]), [0, np.ma.masked, np.ma.masked, 1, 3, 4, 5, 5])

        # list
        a = [0,4,6,2,1,5,3,5]
        ii = argsort(a)
        b = [ a[i] for i in ii ]
        self.assertEqual(b, [0, 1, 2, 3, 4, 5, 5, 6])
        a = [0,4,6,2,1,5,3,5]
        ii = argsort(a, reverse=True)
        b = [ a[i] for i in ii ]
        self.assertEqual(b, [6, 5, 5, 4, 3, 2, 1, 0])
        self.assertRaises(KeyError, argsort, a, reverse=True, key='a')

        # from numpy.argsort docstring
        # One-dimensional array:
        x = np.array([3, 1, 2])
        ii = np.argsort(x)
        self.assertEqual(list(ii), [1, 2, 0])

        # Two-dimensional array:
        x = np.array([[0, 3], [2, 2]])
        ii = np.argsort(x, axis=0)  # sorts along first axis (down)
        self.assertEqual(list(ii.flatten()), [0, 1, 1, 0])
        # same as np.sort(x, axis=0)
        self.assertEqual(list(np.take_along_axis(x, ii, axis=0).flatten()), [0, 2, 2, 3])
        ii = np.argsort(x, axis=1)  # sorts along last axis (across)
        self.assertEqual(list(ii.flatten()), [0, 1, 0, 1])
        # same as np.sort(x, axis=1)
        self.assertEqual(list(np.take_along_axis(x, ii, axis=1).flatten()), [0, 3, 2, 2])

        # Indices of the sorted elements of a N-dimensional array:
        ii = np.unravel_index(np.argsort(x, axis=None), x.shape)
        self.assertEqual(list(ii[0]), [0, 1, 1, 0])
        self.assertEqual(list(ii[1]), [0, 0, 1, 1])
        self.assertEqual(list(x[ii]), [0, 2, 2, 3]) # same as np.sort(x, axis=None)

        # Sorting with keys:
        x = np.array([(1, 0), (0, 1)], dtype=[('x', '<i4'), ('y', '<i4')])
        ii = np.argsort(x, order=('x','y'))
        self.assertEqual(list(ii), [1, 0])
        ii = np.argsort(x, order=('y','x'))
        self.assertEqual(list(ii), [0, 1])


if __name__ == "__main__":
    unittest.main()
