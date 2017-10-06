# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import glob
from os import path

from MDAnalysis.coordinates.core import reader
from MDAnalysis.lib.util import openany

from future.builtins import (
    dict,
)
from future.utils import native_str

import numpy as np
import pandas as pd


class ParamTable(object):
    """Create a parameter table time series for distance or coupling strength.

    """
    _header = ["I", "J"]
    _index = ["segidI", "resI", "I", "segidJ", "resJ", "J"]

    def __init__(self, topdir, prefix="fluctmatch", tbltype="Kb", ressep=3):
        """
        Parameters
        ----------
        topdir : str
            Parent directory for parameter files
        prefix : str, optional
            Filename prefix for files
        tbltype : {"Kb", "b0"}, optional
            Table to create (coupling strength or bond distance)
        ressep : int, optional
            Number of residues to exclude from interactions.
        """
        self._topdir = topdir
        self._prefix = prefix
        self._tbltype = tbltype
        self._ressep = ressep
        self.table = []
        self._filenames = dict(
            intcor="fluct.ic",
            param=".".join((self._prefix, "dist", "prm")),
        )

    def _separate(self, prm_table):
        index = prm_table.index.names
        table = prm_table.reset_index()
        tmp = table[table["segidI"] == table["segidJ"]]
        tmp = tmp[
            (tmp["resI"] >= tmp["resJ"] + self._ressep) |
            (tmp["resJ"] >= tmp["resI"] + self._ressep)
        ]
        diff = table[table["segidI"] != table["segidJ"]]
        table = pd.concat([tmp, diff], axis=0)
        table.set_index(index, inplace=True)
        return table

    def run(self, verbose=False):
        """Create the time series.

        Parameters
        ----------
        verbose : bool, optional
            Print each directory as it is being processed
        """
        revcol = ["segidJ", "resJ", "J", "segidI", "resI", "I"]

        directories = glob.iglob(path.join(self._topdir, "*"))
        for directory in directories:
            if path.isdir(directory):
                if verbose:
                    print("Reading directory {}".format(directory))
                with reader(path.join(directory, self._filenames["intcor"])) as ic_file:
                    if verbose:
                        print("    Processing {}...".format(path.join(directory, self._filenames["intcor"])))
                    ic_table = ic_file.read()
                    ic_table.set_index(self._header, inplace=True)
                with reader(path.join(directory, self._filenames["param"])) as prm_file:
                    if verbose:
                        print("    Processing {}...".format(path.join(directory, self._filenames["param"])))
                    prm_table = prm_file.read()["BONDS"].set_index(self._header)
                table = pd.concat([ic_table, prm_table], axis=1)
                table.reset_index(inplace=True)
                table = table.set_index(self._index)[self._tbltype].to_frame()
                table.columns = [path.basename(directory), ]
                self.table.append(table.copy(deep=True))

        self.table = pd.concat(self.table, axis=1)
        self.table.columns = self.table.columns.astype(np.int)
        self.table = self.table[np.sort(self.table.columns)]

        tmp = self.table.copy(deep=True)
        tmp.index.names = revcol
        self.table = pd.concat([self.table, tmp], axis=0)
        self.table.reset_index(inplace=True)
        self.table = self.table.drop_duplicates(subset=self._index)
        self.table.set_index(self._index, inplace=True)
        self.table.fillna(0., inplace=True)

    def from_file(self, filename):
        """Load a parameter table from a file.

        Parameters
        ----------
        filename : str or stream
            Filename of the parameter table.
        """
        with openany(filename, mode="r") as table:
            self.table = pd.read_csv(
                table,
                skipinitialspace=True,
                delim_whitespace=True,
                header=0,
                index_col=self._index,
            )

    def write(self, filename):
        """Write the parameter table to file.

        Parameters
        ----------
        filename : str or stream
            Location to write the parameter table.
        """
        with openany(filename, mode="w") as table:
            self.table.to_csv(
                table,
                sep=native_str(" "),
                header=True,
                index=True,
                float_format="%.6f",
            )

    @property
    def per_residue(self):
        """Create a single residue time series.

        Returns
        -------
        A table with values per residue.
        """
        # Separate by residue
        table = self._separate(self.table)
        table = 0.5 * table.groupby(level=["segidI", "resI"]).sum()
        table.sort_index(axis=1, inplace=True)
        table = table.reindex(index=table.index, columns=np.sort(table.columns))
        return table

    @property
    def interactions(self):
        """Create a time series for the residue-residue interactions.

        Returns
        -------
        A table with residue-residue values.
        """
        table = self._separate(self.table)
        table = table.groupby(level=["segidI", "resI", "segidJ", "resJ"]).sum()
        table.sort_index(axis=1, inplace=True)
        table = table.reindex(index=table.index, columns=np.sort(table.columns))
        return table
