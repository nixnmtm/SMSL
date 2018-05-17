# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# fluctmatch --- https://github.com/tclick/python-fluctmatch
# Copyright (c) 2013-2017 The fluctmatch Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the New BSD license.
#
# Please cite your use of fluctmatch in published work:
#
# Timothy H. Click, Nixon Raj, and Jhih-Wei Chu.
# Calculation of Enzyme Fluctuograms from All-Atom Molecular Dynamics
# Simulation. Meth Enzymology. 578 (2016), 327-342,
# doi:10.1016/bs.mie.2016.05.024.
#
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)
from future.builtins import range

import functools
import multiprocessing as mp
import os
from os import path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import linalg
from scipy.stats import (scoreatpercentile, t)
from sklearn.utils import extmath


def _rand(avg_kb, std_kb, x):
    return (x, np.random.normal(avg_kb, std_kb))


def randomize(table, ntrials=100):
    """

    Parameters
    ----------
    table : :class:`pandas.DataFrame`
        Table of coupling strengths
    ntrials : int, optional
        Number of trials for eigenvalues

    Returns
    -------
        Array of eigenvalues
    """
    _, Ntime = table.shape
    Lrand = []

    avg_kb = table.mean(axis=1)
    std_kb = table.std(axis=1)

    rand = functools.partial(_rand, avg_kb, std_kb)
    for _ in range(ntrials):
        pool = mp.Pool(maxtasksperchild=2)
        values = pool.map_async(rand, range(Ntime))
        pool.close()
        pool.join()
        kb_rand = pd.DataFrame.from_items(values.get())
        Lrand.append(linalg.svdvals(kb_rand))
    return np.array(Lrand)


def svd(kb):
    """Calculate the singular value decomposition with an appropriate sign flip.

    Parameters
    ----------
    a : (M, N) array_like
        Matrix to decompose.

    Returns
    -------
    U : ndarray
        Unitary matrix having left singular vectors as columns.
        Of shape ``(M, M)`` or ``(M, K)``, depending on `full_matrices`.
    s : ndarray
        The singular values, sorted in non-increasing order.
        Of shape (K,), with ``K = min(M, N)``.
    Vh : ndarray
        Unitary matrix having right singular vectors as rows.
        Of shape ``(N, N)`` or ``(K, N)`` depending on `full_matrices`.
    """
    U, W, Vt = linalg.svd(kb, full_matrices=False)
    U, Vt = extmath.svd_flip(U, Vt, u_based_decision=True)
    return U, W, Vt


def correlate(Usca, Lsca, kmax=6):
    """Calculate the correlation matrix of *Usca* with *Lsca* eigenvalues.

    Parameters
    ----------
    Usca : :class:`numpy.array`
        Eigenvector
    Lsca : :class:`numpy.array`
        Eigenvalue
    kmax : int, optional
        Number of eigenvectors/eigenvalues to use

    Returns
    -------
    Correlation matrix
    """
    S = np.power(Lsca, 2)
    Ucorr = [np.outer(Usca[:, _].dot(S[_]), Usca.T[_]) for _ in range(kmax)]
    return Ucorr


def chooseKpos(Lsca, Lrand, stddev=2):
    """Calculate the significant number of eigenvalues.

    Parameters
    ----------
    Lsca : :class:`numpy.array`
        Eigenvector from coupling strengths
    Lrand : :class:`numpy.array`
        Matrix of eigenvalues from randomized coupling strengths
    stddev : int, optional
        Number of standard deviations to use

    Returns
    -------
    Number of significant eigenvalues
    """
    value = Lrand[:, 1].mean() + ((stddev + 1) * Lrand[:, 1].std())
    return Lsca[Lsca > value].shape[0]


def figUnits(v1,
             v2,
             v3,
             units,
             filename,
             fig_path=os.getcwd(),
             marker='o',
             dotsize=9,
             notinunits=1):
    ''' 3d scatter plot specified by 'units', which must be a list of elements
    in the class Unit_. See figColors_ for the color code. Admissible color codes are in [0 1]
    (light/dark gray can also be obtained by using -1/+1).
    For instance: 0->red, 1/3->green, 2/3-> blue.

    .. _Unit: scaTools.html#scaTools.Unit
    .. _figColors: scaTools.html#scaTools.figColors

    **Arguments:**
       -  `v1` = xvals
       -  `v2` = yvals
       -  `units` = list of elements in units

    **Keyword Arguments**
       -  `marker` = plot marker symbol
       -  `dotsize` = specify marker/dotsize
       -  `notinunits` = if set to 1 : the elements not in a unit are represented in white, if set to 0
                         these elements are not represented, if set to [w1,w2] : elements with coordinates
                         w1,w2 are represented in white in the background.

    :Example:
      >>> figUnits(v1, v2, units, marker='o', gradcol=0, dotsize=9, notinunits=1)

     '''
    import colorsys

    # Plot all items in white:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.elev, ax.azim = 30., 60.
    ax.axes

    if notinunits == 1:
        ax.plot(
            v1,
            v2,
            v3,
            marker,
            markersize=dotsize,
            markerfacecolor='w',
            markeredgecolor='k')
    elif len(notinunits) == 3:
        ax.plot(
            notinunits[0],
            notinunits[1],
            notinunits[2],
            marker,
            markersize=dotsize,
            markerfacecolor='w',
            markeredgecolor='k')

    # Plot items in the units with colors:
    for u in units:
        items_list = list(u.items)
        if u.col >= 0 and u.col < 1:
            bgr = colorsys.hsv_to_rgb(u.col, 1, 1)
        if u.col == 1:
            bgr = [.3, .3, .3]
        if u.col < 0:
            bgr = [.7, .7, .7]
        ax.plot(
            v1[np.ix_(items_list)],
            v2[np.ix_(items_list)],
            v3[np.ix_(items_list)],
            marker,
            markersize=dotsize,
            markerfacecolor=bgr,
            markeredgecolor='k')

    ax.set_xlabel('IC{:d}'.format(1))
    ax.set_ylabel('IC{:d}'.format(2))
    ax.set_zlabel('IC{:d}'.format(3))
    fig.tight_layout()
    fig.savefig(path.join(fig_path, 'svd_ica', filename), dpi=600)


# From pySCA 6.0
class Unit(object):
    """ A class for units (sectors, sequence families, etc.)

        **Attributes:**
            -  `name`  = string describing the unit (ex: "firmicutes")
            -  `items` = set of member items (ex: indices for all firmicutes sequences in an alignment)
            -  `col`   = color code associated to the unit (for plotting)
            -  `vect`  = an additional vector describing the member items (ex: a list of sequence weights)

    """

    def __init__(self):
        self.name = ""
        self.items = set()
        self.col = 0
        self.vect = 0


def icList(Vpica, kpos, Csca, p_cut=0.95):
    """ Produces a list of positions contributing to each independent component (IC) above
    a defined statistical cutoff (p_cut, the cutoff on the CDF of the t-distribution
    fit to the histogram of each IC).  Any position above the cutoff on more than one IC
    are assigned to one IC based on which group of positions to which it shows a higher
    degree of coevolution. Additionally returns the numeric value of the cutoff for each IC, and the
    pdf fit, which can be used for plotting/evaluation.
    icList, icsize, sortedpos, cutoff, pd  = icList(Vsca,Lsca,Lrand) """
    # do the PDF/CDF fit, and assign cutoffs
    Npos = len(Vpica)
    cutoff = []
    scaled_pdf = []
    all_fits = []
    for k in range(kpos):
        pd = t.fit(Vpica[:, k])
        all_fits.append(pd)
        iqr = scoreatpercentile(Vpica[:, k], 75) - scoreatpercentile(
            Vpica[:, k], 25)
        binwidth = 2 * iqr * (len(Vpica[:, k])**(-0.33))
        nbins = round((max(Vpica[:, k]) - min(Vpica[:, k])) / binwidth)
        h_params = np.histogram(Vpica[:, k], nbins.astype(np.int))
        x_dist = np.linspace(min(h_params[1]), max(h_params[1]), num=100)
        area_hist = Npos * (h_params[1][2] - h_params[1][1])
        scaled_pdf.append(area_hist * (t.pdf(x_dist, pd[0], pd[1], pd[2])))
        cd = t.cdf(x_dist, pd[0], pd[1], pd[2])
        tmp = scaled_pdf[k].argmax()
        if abs(max(Vpica[:, k])) > abs(min(Vpica[:, k])):
            tail = cd[tmp:len(cd)]
        else:
            cd = 1 - cd
            tail = cd[0:tmp]
        diff = abs(tail - p_cut)
        x_pos = diff.argmin()
        cutoff.append(x_dist[x_pos + tmp])
    # select the positions with significant contributions to each IC
    ic_init = []
    for k in range(kpos):
        ic_init.append([i for i in range(Npos) if Vpica[i, k] > cutoff[k]])
    # construct the sorted, non-redundant iclist
    sortedpos = []
    icsize = []
    ics = []
    Csca_nodiag = Csca.copy()
    for i in range(Npos):
        Csca_nodiag[i, i] = 0
    for k in range(kpos):
        icpos_tmp = list(ic_init[k])
        for kprime in [kp for kp in range(kpos) if (kp != k)]:
            tmp = [v for v in icpos_tmp if v in ic_init[kprime]]
            for i in tmp:
                remsec = np.linalg.norm(Csca_nodiag[i, ic_init[k]]) \
                         < np.linalg.norm(Csca_nodiag[i, ic_init[kprime]])
                if remsec:
                    icpos_tmp.remove(i)
        sortedpos += sorted(icpos_tmp, key=lambda i: -Vpica[i, k])
        icsize.append(len(icpos_tmp))
        s = Unit()
        s.items = sorted(icpos_tmp, key=lambda i: -Vpica[i, k])
        s.col = k / kpos
        s.vect = -Vpica[s.items, k]
        ics.append(s)
    return ics, icsize, sortedpos, cutoff, scaled_pdf, all_fits


def basicICA(x, r, Niter):
    """ Basic ICA algorithm, based on work by Bell & Sejnowski (infomax). The input data should preferentially be sphered, i.e., x.T.dot(x) = 1

    **Arguments:**
      -  `x` = LxM input matrix where L = # features and M = # samples
      -  `r` = learning rate / relaxation parameter (e.g. r=.0001)
      -  `Niter` =  number of iterations (e.g. 1000)

    **Returns:**
      -  `w` = unmixing matrix
      -  `change` = record of incremental changes during the iterations.

    **Note:** r and Niter should be adjusted to achieve convergence, which
    should be assessed by visualizing "change" with plot(range(iter) ,change)

    **Example:**
      >>> [w, change] = basicICA(x, r, Niter)

    """
    [L, M] = x.shape
    w = np.eye(L)
    change = []
    for _ in range(Niter):
        w_old = np.copy(w)
        u = w.dot(x)
        w += r * (M * np.eye(L) + (1 - 2 * (1. /
                                            (1 + np.exp(-u)))).dot(u.T)).dot(w)
        delta = (w - w_old).ravel()
        change.append(delta.dot(delta.T))
    return [w, change]


def rotICA(V, kmax=6, learnrate=.0001, iterations=10000):
    """ ICA rotation (using basicICA) with default parameters and normalization of
    outputs.

    :Example:
       >>> Vica, W = rotICA(V, kmax=6, learnrate=.0001, iterations=10000)
    """
    V1 = V[:, :kmax].T
    [W, changes_s] = basicICA(V1, learnrate, iterations)
    Vica = (W.dot(V1)).T
    for n in range(kmax):
        imax = abs(Vica[:, n]).argmax()
        Vica[:, n] = np.sign(Vica[imax, n]) * Vica[:, n] / np.linalg.norm(
            Vica[:, n])
    return Vica, W
