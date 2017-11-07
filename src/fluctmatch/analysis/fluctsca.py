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

import functools
import multiprocessing as mp

import numpy as np
import pandas as pd
from future.builtins import range
from scipy import linalg
from scipy.stats import (scoreatpercentile, t)


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
    Ucorr = [
        np.outer(Usca[:, _].dot(S[_]), Usca.T[_])
        for _ in range(kmax)
    ]
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


def figUnits(v1, v2, v3, units, filename, marker='o', dotsize=9, notinunits=1):
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

    Ntot = len(v1)
    # Plot all items in white:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.elev, ax.azim = 30., 60.
    ax.axes

    if notinunits == 1:
        ax.plot(v1, v2, v3, marker, markersize=dotsize, markerfacecolor='w',
                markeredgecolor='k')
    elif len(notinunits) == 3:
        ax.plot(notinunits[0], notinunits[1], notinunits[2], marker, markersize=dotsize,
                markerfacecolor='w', markeredgecolor='k')

    # Plot items in the units with colors:
    for u in units:
        items_list = list(u.items)
        if u.col >= 0 and u.col < 1:
            bgr = colorsys.hsv_to_rgb(u.col, 1, 1)
        if u.col == 1:
            bgr = [.3, .3, .3]
        if u.col < 0:
            bgr = [.7, .7, .7]
        ax.plot(v1[np.ix_(items_list)], v2[np.ix_(items_list)], v3[np.ix_(items_list)],
                marker, markersize=dotsize, markerfacecolor=bgr, markeredgecolor='k')

    ax.set_xlabel('IC{:d}'.format(1))
    ax.set_ylabel('IC{:d}'.format(2))
    ax.set_zlabel('IC{:d}'.format(3))
    fig.tight_layout()
    fig.savefig(path.join(fig_path, 'svd_ica', filename), dpi=600)
