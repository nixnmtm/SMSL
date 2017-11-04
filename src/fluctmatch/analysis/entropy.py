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

import numpy as np
import pandas as pd
from scipy import stats

from fluctmatch.analysis.paramtable import ParamTable


class Entropy(object):
    """Calculate various entropic contributions from the coupling strengths.
    """
    def __init__(self, filename, ressep=3):
        """
        Parameters
        ----------
        filename : str
            Name of the file with a complete coupling strength time series
        ressep : int, optional
            Residue separation
        """
        self._table = ParamTable(ressep=ressep)
        self._table.from_file(filename=filename)

    def relative_entropy(self):
        """Calculate the relative entropy between windows.

        Calculate the relative entropy of a coarse-grain system by using the
        following equations:

        .. math::

            P^t_{IJ} = K^t_{IJ} / \sum {K^t_{IJ}}
            S = \sum {P^t_{IJ} ln(P^t_{IJ}}

        where :math:`t/2` is the previous window because of the overlap in
        trajectory times by :math:`t/2`. If :math:`Q^{t/2}_{IJ} = 0.0`, then
        it is replaced with a penalty value determined by a force constant value
        as determined by the highest probability of nonzero force constants
        within the overall time series.

        The required table must be a time series of an internal coordinates
        table containing bond force constants.

        Returns
        -------
        Timeseries of relative entropy per residue
        """
        header = ["segidI", "resI"]

        entropy = self._table._separate(self._table)
        entropy = entropy.groupby(level=header).apply(lambda x: x / x.sum())
        entropy = entropy.groupby(level=header).agg(stats.entropy)
        entropy.replace(-np.inf, 0., inplace=True)

        return entropy

    def coupling_entropy(self):
        """Calculate the entropic contributions between residues.

        Returns
        -------
        Entropy of residue-residue interactions over time
        """
        # Transpose matrix because stats.entropy does row-wise calculation.
        table = self._table.per_residue.T
        entropy = stats.entropy(table)
        return pd.DataFrame(entropy, index=table.index)

    def windiff_entropy(self, bins=100):
        """Calculate the relative entropy between windows.

        Calculate the relative entropy of a coarse-grain system by using the
        following equations:

        .. math::

            P^t_{IJ} = K^t_{IJ} / \leftangle K^t_{IJ}\rightangle
            Q^t/2_{IJ} = K^{t/2}_{IJ} / \leftangle K^{t/2}_{IJ}
            S = -0.5 * (\sum(P^t_{IJ}\ln(P^t_{IJ} / Q^{t/2}_{IJ})) +
                \sum(Q^{t/2}_{IJ})\ln(P^t_{IJ} / Q^{t/2}_{IJ})))

        where :math:`t/2` is the previous window because of the overlap in
        trajectory times by :math:`t/2`. If :math:`Q^{t/2}_{IJ} = 0.0`, then it
        is replaced with a penalty value determined by a force constant value as
        determined by the highest probability of nonzero force constants within
        the overall time series.

        The required table must be a time series of an internal coordinates table
        containing bond force constants.

        Parameters
        ----------
        bins : int, optional
            Number of bins for histogram determination

        Returns
        -------
        Entropy difference between two windows
        """
        # Calculate value of maximum probability and define penalty value.
        header = ["segidI", "resI"]

        table = self._table._separate(self._table)
        hist, edges = np.histogram(
            table, range=(1e-4, table.values.max()), bins=bins
        )
        hist = (hist / table.size).astype(dtype=np.float)
        xaxis = (edges[:-1] + edges[1:]) / 2
        try:
            penalty = xaxis[np.where(hist == hist.max())][0]
        except IndexError:
            penalty = xaxis[-1]

        # Calculate average coupling strength per residue.
        table[table == 0.] = penalty
        #     meanI = tmp.groupby(level=["resI"]).mean()
        table = table.groupby(level=header).transform(normalize)

        # Utilize the numpy arrays for calculations
        P = table[table.columns[1:]]
        Q = table[table.columns[:-1]]
        P.columns = Q.columns

        # Caclulate the relative entropy
        S_P = P * np.log(P / Q)
        S_Q = Q * np.log(P / Q)
        S_P.fillna(0., inplace=True)
        S_Q.fillna(0., inplace=True)
        entropy = -(
            S_P.groupby(level=header).sum() + S_Q.groupby(level=header).sum()
        )
        entropy[entropy == -0.0] = entropy[entropy == -0.0].abs()

        return entropy
