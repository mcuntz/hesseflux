#!/usr/bin/env python
"""
Module provides test functions for parameter sensitivity analysis from
    Ishigami and Homma (1990) An importance qualification technique in uncertainty analysis for computer models,
                              Proceedings of the isuma '90, First International Symposium on Uncertainty
                              Modelling and Analysis, University of Maryland, Dec. 03 - Dec 05 1990, 398-403

    Oakley and O'Hagan (2004) Probabilistic sensitivity analysis of complex models: a Bayesian approach
                              J. R. Statist. Soc. B 66, Part 3, 751-769.

    Morris (1991)             Factorial sampling plans for preliminary computational experiments,
                              Technometrics 33, 161-174.

    Saltelli et al. (2008)    Global Sensitivity Analysis. The Primer, John Wiley & Sons, pp. 292

    Saltelli et al. (2010)    Variance based sensitivity analysis of model output, Design and estimator
                              for the total sensitivity index, Comp. Phys. Comm. 181, 259-270.

    Sobol' (1990),            Sensitivity estimates for nonlinear mathematical models,
                              Matematicheskoe Modelirovanie 2, 112-118 (in Russian),
                              translated in English in Sobol' (1993).

    Sobol' (1993)             Sensitivity analysis for non-linear mathematical models,
                              Mathematical Modelling and Computational Experiment 1, 407-414,
                              English translation of Russian original paper Sobol' (1990).

Current functions are:
    B                     B of Saltelli et al. (2010)

    G / g                 G-function attributed to Sobol' (1990, 1993), given by Saltelli et al. (2008, 2010)

    Gstar                 G* of Saltelli et al. (2010)

    ishigami_homma        Ishigami and Homma (1990), given by Saltelli et al. (2008, page 179)

    K / bratley           K  of Saltelli et al. (2010)

    fmorris / morris      After Morris (1991)

    oakley_ohagan         Oakley and O'Hagan (2004), parameters given in Saltelli et al. (2008)
                          or on http://www.jeremy-oakley.staff.shef.ac.uk/psa_example.txt

This module was written by Matthias Cuntz & Juliane Mai while at
Department of Computational Hydrosystems, Helmholtz Centre for
Environmental Research - UFZ, Leipzig, Germany, and continued by
Matthias Cuntz while at Institut National de Recherche pour l'Agriculture,
l'Alimentation et l'Environnement (INRAE), Nancy, France.

Copyright (c) 2015-2020 Matthias Cuntz - mc (at) macu (dot) de
Released under the MIT License; see LICENSE file for details.

* Written Mar 2015 by Matthias Cuntz (mc (at) macu (dot) de) & Juliane Mai
* Added functions to properly test PAWN method: linear, product, ratio, and ishigami_homma_easy, Dec 2017, Juliane Mai
* Provide morris function under the name fmorris and the K function under the name bratley, Nov 2019, Matthias Cuntz
* Changed to Sphinx docstring and numpydoc, Dec 2019, Matthias Cuntz
* Distinguish iterable and array_like parameter types, Jan 2020, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz

The following functions are provided:

.. autosummary::
    B
    g
    G
    Gstar
    K
    bratley
    fmorris
    morris
    oakley_ohagan
    ishigami_homma
    linear
    product
    ratio
    ishigami_homma_easy
"""
from __future__ import division, absolute_import, print_function
import numpy as np


__all__ = ['B', 'g', 'G', 'Gstar', 'K', 'bratley', 'fmorris', 'morris', 'oakley_ohagan', 'ishigami_homma',
           'linear', 'product', 'ratio', 'ishigami_homma_easy']


# -----------------------------------------------------------

def B(X):
    """
    B function, Saltelli et al. (2010) Comp. Phys. Comm. 181, p. 259-270

    Parameters
    ----------
    X : array_like
        (nX,) or (nX,npoints) array of floats

    Returns
    -------
    B : float or ndarray
        float or (npoints,) floats of B function values at `X`
    """
    # Parameter sets are assumed to be in following ordering: (x_1, x_2, ..., X_m, w_1, w_2, ..., w_m)
    X = np.array(X)
    if X.ndim == 1:
        isone = True
        iX = X[:,np.newaxis]
    else:
        isone = False
        iX = X

    assert iX.shape[0] % 2 == 0, 'X.shape[0] must be even.'

    m  = iX.shape[0]//2
    y  = np.sum(iX[:m,:]*iX[m:,:], axis=0)

    if isone:
        return y[0]
    else:
        return y


# -----------------------------------------------------------


def g(X, a):
    """
    G-function

    Sobol' (1990) Matematicheskoe Modelirovanie 2, 112-118 (in Russian)

    Sobol' (1993) Mathematical Modelling and Computational Experiment 1, 407-414 (English translation)

    Parameters
    ----------
    X : array_like
        (nX,) or (nX,npoints) array of floats
    a : array_like
        (nX,) array of floats

    Returns
    -------
    G : float or ndarray
        float or (npoints,) floats of G function values at `X` with parameters `a`
    """
    return Gstar(X, np.ones(len(a)), np.zeros(len(a)), a)


# -----------------------------------------------------------


def G(X, a):
    """
    G-function

    Sobol' (1990) Matematicheskoe Modelirovanie 2, 112-118 (in Russian)

    Sobol' (1993) Mathematical Modelling and Computational Experiment 1, 407-414 (English translation)

    Parameters
    ----------
    X : array_like
        (nX,) or (nX,npoints) array of floats
    a : array_like
        (nX,) array of floats

    Returns
    -------
    g : float or ndarray
        float or (npoints,) floats of G function values at `X` with parameters `a`
    """
    return Gstar(X, np.ones(len(a)), np.zeros(len(a)), a)


# -----------------------------------------------------------


def Gstar(X, alpha, delta, a):
    """
    G* example, Saltelli et al. (2010) Comp. Phys. Comm., 181, p. 259-270

    Parameters
    ----------
    X : array_like
        (nX,) or (nX,npoints) array of floats
    alpha : array_like
        (nX,) array of floats
    delta : array_like
        (nX,) array of floats
    a : array_like
        (nX,) array of floats

    Returns
    -------
    G* : float or ndarray
        float or (npoints,) floats of G* function values at `X` with parameters `alpha`, `delta` and `a`
    """
    # Model output for given parameter set(s) is returned
    # X: dim1 = # of parameters = 10
    #       dim2 = # of parameter sets
    X = np.array(X)
    if X.ndim == 1:
        isone = True
        iX = X[:,np.newaxis]
    else:
        isone = False
        iX = X

    if isinstance(alpha, (list, tuple)):
        nalpha = len(alpha)
    else:
        nalpha = alpha.size

    assert iX.shape[0] == nalpha, 'X.shape[0] must len(alpha).'

    ialpha = np.array(alpha).reshape((nalpha,1))
    idelta = np.array(delta).reshape((nalpha,1))
    ia     = np.array(a).reshape((nalpha,1))
    yi = ((1.+ialpha)*np.abs(2.*(iX+idelta-np.trunc(iX+idelta))-1.)**ialpha+ia)/(1.+ia)
    y  = np.prod(yi, axis=0)

    if isone:
        return y[0]
    else:
        return y


# -----------------------------------------------------------


def linear(X, a, b):
    """
    Linear test function to test PAWN method:

        `Y = a*X + b`

    Parameters
    ----------
    X : array_like
        (1,) or (1,npoints) array of floats
    a : array_like
        float or (npoints,) array of floats
    b : array_like
        float or (npoints,) array of floats

    Returns
    -------
    linear : float or ndarray
        float or (npoints,) floats of linear function values at `X` with parameters `a` and `b`
    """
    # Model output for given parameter set(s) is returned
    # X:    dim1 = # of parameters = 1
    #       dim2 = # of parameter sets
    X = np.array(X)
    if X.ndim == 1:
        isone = True
        iX = X[:,np.newaxis]
    else:
        isone = False
        iX = X

    assert iX.shape[0] == 1, 'X.shape[0] must 1.'

    y = a * iX[0,:] + b

    if isone:
        return y[0]
    else:
        return y


# -----------------------------------------------------------


def product(X):
    """
    Product test function to test PAWN method:

        `Y = X[0] * X[1]`

    Parameters
    ----------
    X : array_like
        (2,) or (2,npoints) array of floats

    Returns
    -------
    product : float or ndarray
        float or (npoints,) floats of product function values at `X`
    """
    # Model output for given parameter set(s) is returned
    # X:    dim1 = # of parameters = 2
    #       dim2 = # of parameter sets
    X = np.array(X)
    if X.ndim == 1:
        isone = True
        iX = X[:,np.newaxis]
    else:
        isone = False
        iX = X

    assert iX.shape[0] == 2, 'X.shape[0] must 2.'

    y = iX[0,:] * iX[1,:]

    if isone:
        return y[0]
    else:
        return y


# -----------------------------------------------------------


def ratio(X):
    """
    Ratio test function:

        `Y = X[0] / X[1]`

    Simple nonlinear model proposed by Liu et al. (2006):

        Liu, H., Sudjianto, A., Chen, W., 2006.
        Relative entropy based method for probabilistic sensitivity analysis in engineering design.
        J. Mech. Des. 128, 326-336.

    Used by Pianosi & Wagener, Environmental Modelling & Software (2015)

        Pianosi, F. & Wagener T., 2015
        A simple and efficient method for global sensitivity analysis based on cumulative distribution functions.
        Environmental Modelling & Software 67, 1-11.

    Parameters
    ----------
    X : array_like
        (2,) or (2,npoints) array of floats

    Returns
    -------
    ratio : float or ndarray
        float or (npoints,) floats of ratio function values at `X`
    """
    # Model output for given parameter set(s) is returned
    # X:    dim1 = # of parameters = 2
    #       dim2 = # of parameter sets
    X = np.array(X)
    if X.ndim == 1:
        isone = True
        iX = X[:,np.newaxis]
    else:
        isone = False
        iX = X

    assert iX.shape[0] == 2, 'X.shape[0] must 2.'

    y = iX[0,:] / iX[1,:]

    if isone:
        return y[0]
    else:
        return y


# -----------------------------------------------------------


def ishigami_homma_easy(X):
    """
    Simplified Ishigami and Homma function to test PAWN method:

        `Y = sin(X[0]) + X[1]`

    with `X[0],X[1]` ~ `Uniform[-Pi,Pi]`

    Parameters
    ----------
    X : array_like
        (2,) or (2,npoints) array of floats

    Returns
    -------
    ishigami_homma_easy : float or ndarray
        float or (npoints,) floats of simplified Ishigami and Homma function values at `X`
    """
    # Model output for given parameter set(s) is returned
    # X:    dim1 = # of parameters = 2
    #       dim2 = # of parameter sets
    X = np.array(X)
    if X.ndim == 1:
        isone = True
        iX = X[:,np.newaxis]
    else:
        isone = False
        iX = X

    assert iX.shape[0] == 2, 'X.shape[0] must 2.'

    y = np.sin(iX[0,:])  + iX[1,:]

    if isone:
        return y[0]
    else:
        return y


# -----------------------------------------------------------


def ishigami_homma(X, a, b):
    """
    Ishigami and Homma (1990), given by Saltelli et al. (2008, page 179)

    Parameters
    ----------
    X : array_like
        (3,) or (3,npoints) array of floats
    a : array_like
        float or (npoints,) array of floats
    b : array_like
        float or (npoints,) array of floats

    Returns
    -------
    ishigami_homma : float or ndarray
        float or (npoints,) floats of Ishigami and Homma function values at `X` with parameters `a` and `b`
    """
    # Model output for given parameter set(s) is returned
    # X: dim1 = # of parameters = 3
    #       dim2 = # of parameter sets
    X = np.array(X)
    if X.ndim == 1:
        isone = True
        iX = X[:,np.newaxis]
    else:
        isone = False
        iX = X

    assert iX.shape[0] == 3, 'X.shape[0] must be 3.'

    y = np.sin(iX[0,:])  + a*(np.sin(iX[1,:]))**2  + b*iX[2,:]**4  * np.sin(iX[0,:])

    if isone:
        return y[0]
    else:
        return y


# -----------------------------------------------------------


def K(X):
    """
    K example, Saltelli et al. (2010) Comp. Phys. Comm., 181, p. 259-270

    Parameters
    ----------
    X : array_like
        (nX,) or (nX,npoints) array of floats

    Returns
    -------
    K : float or ndarray
        float or (npoints,) floats of K function values at `X`
    """
    # Model output for given parameter set(s) is returned
    # X: dim1 = # of parameters = 10
    #    dim2 = # of parameter sets
    X = np.array(X)
    if X.ndim == 1:
        isone = True
        iX = X[:,np.newaxis]
    else:
        isone = False
        iX = X

    assert iX.shape[0] > 1, 'X.shape[0] must be > 1.'

    nX = iX.shape[0]
    for ii in range(1,nX):
        iX[ii,:] = iX[ii-1,:] * iX[ii,:]
    for ii in range(nX):
        iX[ii,:] = (-1)**(ii+1) * iX[ii,:]
    y = iX.sum(axis=0)

    if isone:
        return y[0]
    else:
        return y


def bratley(*args):
    """
    K example, Saltelli et al. (2010) Comp. Phys. Comm., 181, p. 259-270

    Parameters
    ----------
    X : array_like
        (nX,) or (nX,npoints) array of floats

    Returns
    -------
    bratley : float or ndarray
        float or (npoints,) floats of K function values at `X`
    """
    return K(*args)


# -----------------------------------------------------------


def fmorris(X, beta0, beta1, beta2, beta3, beta4):
    """
    Morris-function, Morris (1991) Technometrics 33, 161-174

    Parameters
    ----------
    X : array_like
        (20,) or (20,npoints) array of floats
    beta0 : float
        float
    beta1 : array_like
        (20,) array of floats
    beta2 : array_like
        (20,20) array of floats
    beta3 : array_like
        (20,20,20) array of floats
    beta4 : array_like
        (20,20,20,20) array of floats

    Returns
    -------
    fmorris : float or ndarray
        float or (npoints,) floats of Morris function values at `X` with parameters `beta0-beta4`
    """
    X = np.array(X)
    if X.ndim == 1:
        isone = True
        iX = X[:,np.newaxis]
    else:
        isone = False
        iX = X

    om = 2.*(iX - 0.5)
    ii = np.array((2, 4, 6))
    om[ii,:] = 2.*(1.1*iX[ii,:]/(iX[ii,:]+0.1) - 0.5)

    nn  = 20

    assert iX.shape[0] == nn, 'X.shape[0] must '+str(nn)+'.'
    assert np.all(np.array(beta1.shape) == nn), 'beta1.shape must be all '+str(nn)+'.'
    assert np.all(np.array(beta2.shape) == nn), 'beta2.shape must be all '+str(nn)+'.'
    assert np.all(np.array(beta3.shape) == nn), 'beta3.shape must be all '+str(nn)+'.'
    assert np.all(np.array(beta4.shape) == nn), 'beta4.shape must be all '+str(nn)+'.'

    y    = np.empty(iX.shape[1])
    y[:] = beta0
    for i in range(0,nn):
        y[:] += beta1[i]*om[i,:]
        for j in range(i+1,nn):
            y[:] += beta2[i,j]*om[i,:]*om[j,:]
            for l in range(j+1,nn):
                y[:] += beta3[i,j,l]*om[i,:]*om[j,:]*om[l,:]
                for s in range(l+1,nn):
                    y[:] += beta4[i,j,l,s]*om[i,:]*om[j,:]*om[l,:]*om[s,:]

    if isone:
        return y[0]
    else:
        return y


def morris(*args):
    """
    Morris-function, Morris (1991) Technometrics 33, 161-174

    Parameters
    ----------
    X : array_like
        (20,) or (20,npoints) array of floats
    beta0 : float
        float
    beta1 : array_like
        (20,) array of floats
    beta2 : array_like
        (20,20) array of floats
    beta3 : array_like
        (20,20,20) array of floats
    beta4 : array_like
        (20,20,20,20) array of floats

    Returns
    -------
    morris : float or ndarray
        float or (npoints,) floats of Morris function values at `X` with parameters `beta0-beta4`
    """
    return fmorris(*args)


# -----------------------------------------------------------


def oakley_ohagan(X):
    """
    Oakley and O'Hagan (2004) J. R. Statist. Soc. B 66, Part 3, 751-769

    Parameters
    ----------
    X : array_like
        (15,) or (15,npoints) array of floats

    Returns
    -------
    oakley_ohagan : float or ndarray
        float or (npoints,) floats of Oakley and O'Hagan function values at `X`
    """
    X = np.array(X)
    if X.ndim == 1:
        isone = True
        iX = X[:,np.newaxis]
    else:
        isone = False
        iX = X

    nn = 15

    assert iX.shape[0] == nn, 'X.shape[0] must '+str(nn)+'.'

    a1 = np.array([0.01, 0.05, 0.23, 0.04, 0.12, 0.39, 0.39, 0.61, 0.62, 0.40, 1.07, 1.15, 0.79, 1.12, 1.20])
    a2 = np.array([0.43, 0.09, 0.05, 0.32, 0.15, 1.04, 0.99, 0.97, 0.90, 0.81, 1.84, 2.47, 2.39, 2.00, 2.26])
    a3 = np.array([0.10, 0.21, 0.08, 0.27, 0.13, 0.75, 0.86, 1.03, 0.84, 0.80, 2.21, 2.04, 2.40, 2.05, 1.98])
    M  = np.array([
        [-0.02, -0.19,   0.13,  0.37,  0.17,  0.14, -0.44,  -0.08,  0.71, -0.44,  0.5,  -0.02, -0.05,  0.22,  0.06],
        [ 0.26,  0.05,   0.26,  0.24, -0.59, -0.08, -0.29,   0.42,  0.5,   0.08, -0.11,  0.03, -0.14, -0.03, -0.22],
        [-0.06,  0.2,    0.1,  -0.29, -0.14,  0.22,  0.15,   0.29,  0.23, -0.32, -0.29, -0.21,  0.43,  0.02,  0.04],
        [ 0.66,  0.43,   0.3,  -0.16, -0.31, -0.39,  0.18,   0.06,  0.17,  0.13, -0.35,  0.25, -0.02,  0.36, -0.33],
        [-0.12,  0.12,   0.11,  0.05, -0.22,  0.19, -0.07,   0.02, -0.1,   0.19,  0.33,  0.31, -0.08, -0.25,  0.37],
        [-0.28, -0.33,  -0.1,  -0.22, -0.14, -0.14, -0.12,   0.22, -0.03, -0.52,  0.02,  0.04,  0.36,  0.31,  0.05],
        [-0.08,  0.004,  0.89, -0.27, -0.08, -0.04, -0.19,  -0.36, -0.17,  0.09,  0.4,  -0.06,  0.14,  0.21, -0.01],
        [-0.09,  0.59,   0.03, -0.03, -0.24, -0.1,   0.03,   0.1,  -0.34,  0.01, -0.61,  0.08,  0.89,  0.14,  0.15],
        [-0.13,  0.53,   0.13,  0.05,  0.58,  0.37,  0.11,  -0.29, -0.57,  0.46, -0.09,  0.14, -0.39, -0.45, -0.15],
        [ 0.06, -0.32,   0.09,  0.07, -0.57,  0.53,  0.24,  -0.01,  0.07,  0.08, -0.13,  0.23,  0.14, -0.45, -0.56],
        [ 0.66,  0.35,   0.14,  0.52, -0.28, -0.16, -0.07,  -0.2,   0.07,  0.23, -0.04, -0.16,  0.22,  0,    -0.09],
        [ 0.32, -0.03,   0.13,  0.13,  0.05, -0.17,  0.18,   0.06, -0.18, -0.31, -0.25,  0.03, -0.43, -0.62, -0.03],
        [-0.29,  0.03,   0.03, -0.12,  0.03, -0.34, -0.41,   0.05, -0.27, -0.03,  0.41,  0.27,  0.16, -0.19,  0.02],
        [-0.24, -0.44,   0.01,  0.25,  0.07,  0.25,  0.17,   0.01,  0.25, -0.15, -0.08,  0.37, -0.3,   0.11, -0.76],
        [ 0.04, -0.26,   0.46, -0.36, -0.95, -0.17,  0.003,  0.05,  0.23,  0.38,  0.46, -0.19,  0.01,  0.17,  0.16] ])

    y = np.dot(a1,iX) + np.dot(a2,np.sin(iX)) + np.dot(a3,np.cos(iX))
    for i in range(iX.shape[1]):
        y[i] += np.dot(iX[:,i].T,np.dot(M,iX[:,i]))

    if isone:
        return y[0]
    else:
        return y.squeeze()


# -----------------------------------------------------------


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # ishigami_homma([np.pi/2.,np.pi/2.,1.], 1., 1.)
    # # 3.0
    # B(np.arange(10))
    # # 80
    # B(np.arange(12).reshape(6,2))
    # # array([56, 89])
    # linear(np.ones(1), 1., 1.)
    # # 2.0
    # product(np.arange(2)+1.)
    # # 2.0
    # ratio(np.arange(2)+1.)
    # # 0.5
    # ishigami_homma_easy([np.pi/2.,1.])
    # # 2.0
    # bratley(np.arange(5)+1.)
    # K(np.arange(5)+1.)
    # # -101.0
    # K(np.arange(8).reshape((4,2))+1.)
    # # [ 92., 342.]
    # oakley_ohagan(np.zeros(15))
    # # 15.75
    # Gstar(np.ones(5), np.zeros(5), np.ones(5), np.zeros(5))
    # # 1.0
    # g(np.ones(5), np.zeros(5))
    # G(np.ones(5), np.zeros(5))
    # # 32.0

    # # Morris function
    # seed = 1234
    # np.random.seed(seed=seed)
    # npars = 20
    # x0 = np.ones(npars)*0.5
    # lb = np.zeros(npars)
    # ub = np.ones(npars)
    # beta0              = 0.
    # beta1              = np.random.standard_normal(npars)
    # beta1[:10]         = 20.
    # beta2              = np.random.standard_normal((npars,npars))
    # beta2[:6,:6]       = -15.
    # beta3              = np.zeros((npars,npars,npars))
    # beta3[:5,:5,:5]    = -10.
    # beta4              = np.zeros((npars,npars,npars,npars))
    # beta4[:4,:4,:4,:4] = 5.
    # mm = fmorris(np.linspace(0,2*(npars-1),npars)/float(2*npars-1), beta0, beta1, beta2, beta3, beta4)
    # print(np.around(mm,3))
    # # -82.711
    # mm = fmorris(np.arange(2*npars,dtype=np.float).reshape((npars,2))/float(2*npars-1),
    #              beta0, beta1, beta2, beta3, beta4)
    # print(np.around(mm,3))
    # # [-82.711 -60.589]
