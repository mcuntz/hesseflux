#!/usr/bin/env python
"""
    This is the unittest for the Sensitivity Analysis Test Functions module.

    python -m unittest -v test_sa_test_functions.py
    python -m pytest --cov-report term-missing -v tests/test_sa_test_functions.py
"""
from __future__ import division, absolute_import, print_function
import unittest

# --------------------------------------------------------------------
# sa_test_functions.py
class TestSATestFunctions(unittest.TestCase):

    def setUp(self):
        import numpy as np
        # seed for reproducible results
        seed = 1234
        np.random.seed(seed=seed)

    def test_sa_test_functions(self):
        import os
        import numpy as np
        from hesseflux.functions import B, g, G, Gstar, K, bratley, oakley_ohagan, ishigami_homma
        from hesseflux.functions import linear, product, ratio, ishigami_homma_easy, fmorris, morris

        # scalar
        self.assertEqual(B(np.arange(10)), 80)
        self.assertEqual(g(np.ones(5), np.zeros(5)), 32.0)
        self.assertEqual(G(np.ones(5), np.zeros(5)), 32.0)
        self.assertEqual(Gstar(np.ones(5), np.zeros(5), np.ones(5), np.zeros(5)), 1.0)
        self.assertEqual(Gstar(np.ones(5), [0.,0.,0.,0.,0.], np.ones(5), np.zeros(5)), 1.0)
        self.assertEqual(K(np.arange(5)+1.), -101.0)
        self.assertEqual(bratley(np.arange(5)+1.), -101.0)
        self.assertEqual(oakley_ohagan(np.zeros(15)), 15.75)
        self.assertEqual(ishigami_homma([np.pi/2.,np.pi/2.,1.], 1., 1.), 3.0)
        self.assertEqual(linear(np.ones(1), 1., 1.), 2.0)
        self.assertEqual(product(np.arange(2)+1.), 2.0)
        self.assertEqual(ratio(np.arange(2)+1.), 0.5)
        self.assertEqual(ishigami_homma_easy([np.pi/2.,1.]), 2.0)

        # vector
        self.assertEqual(list(B(np.arange(12).reshape(6,2))), [56, 89])
        self.assertEqual(list(g(np.ones((5,2)), np.zeros(5))), [32.0, 32.0])
        self.assertEqual(list(G(np.ones((5,2)), np.zeros(5))), [32.0, 32.0])
        self.assertEqual(list(Gstar(np.ones((5,2)), np.zeros(5), np.ones(5), np.zeros(5))), [1.0, 1.0])
        self.assertEqual(list(K(np.arange(8).reshape((4,2))+1.)), [92., 342.])
        self.assertEqual(list(bratley(np.arange(8).reshape((4,2))+1.)), [92., 342.])
        self.assertEqual(list(oakley_ohagan(np.zeros((15,2)))), [15.75, 15.75])
        self.assertEqual(list(ishigami_homma([[np.pi/2.,np.pi/2.],[np.pi/2.,np.pi/2.],[1.,1.]], 1., 1.)), [3.0, 3.0])
        self.assertEqual(list(linear(np.ones((1,2)), 1., 1.)), [2.0, 2.0])
        self.assertEqual(list(product(np.arange(4).reshape((2,2))+1.)), [3.0, 8.0])
        self.assertEqual(list(ratio(np.arange(2).repeat(2).reshape((2,2))+1.)), [0.5, 0.5])
        self.assertEqual(list(ishigami_homma_easy([[np.pi/2.,np.pi/2.],[1.,1.]])), [2.0, 2.0])

        # Morris
        npars = 20
        x0 = np.ones(npars)*0.5
        lb = np.zeros(npars)
        ub = np.ones(npars)
        beta0              = 0.
        beta1              = np.random.standard_normal(npars)
        beta1[:10]         = 20.
        beta2              = np.random.standard_normal((npars,npars))
        beta2[:6,:6]       = -15.
        beta3              = np.zeros((npars,npars,npars))
        beta3[:5,:5,:5]    = -10.
        beta4              = np.zeros((npars,npars,npars,npars))
        beta4[:4,:4,:4,:4] = 5.

        mm = fmorris(np.linspace(0,2*(npars-1),npars)/float(2*npars-1),
                     beta0, beta1, beta2, beta3, beta4)
        self.assertEqual(np.around(mm,3), -82.711)
        mm = fmorris(np.arange(2*npars,dtype=np.float).reshape((npars,2))/float(2*npars-1),
                     beta0, beta1, beta2, beta3, beta4)
        self.assertEqual(list(np.around(mm,3)), [-82.711, -60.589])

        mm = morris(np.linspace(0,2*(npars-1),npars)/float(2*npars-1),
                    beta0, beta1, beta2, beta3, beta4)
        self.assertEqual(np.around(mm,3), -82.711)
        mm = morris(np.arange(2*npars,dtype=np.float).reshape((npars,2))/float(2*npars-1),
                    beta0, beta1, beta2, beta3, beta4)
        self.assertEqual(list(np.around(mm,3)), [-82.711, -60.589])


if __name__ == "__main__":
    unittest.main()
