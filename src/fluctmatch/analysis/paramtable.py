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
from future.builtins import (
    dict, )
from future.utils import native_str

import functools
import glob
import multiprocessing as mp
from os import path

import numpy as np
import pandas as pd
from MDAnalysis.coordinates.core import reader
from MDAnalysis.lib.util import openany

_header = ["I", "J"]
_index = dict(
    general=["segidI", "resI", "I", "segidJ", "resJ", "J"],
    complete=["segidI", "resI", "resnI", "I", "segidJ", "resJ", "resnJ", "J"],
)


def _create_table(directory,
                  intcor="average.ic",
                  parmfile="fluctmatch.dist.prm",
                  tbltype="Kb",
                  trimmed=False,
                  verbose=False):
    if path.isdir(directory):
        if verbose:
            print("Reading directory {}".format(directory))
        with reader(path.join(directory, intcor)) as ic_file:
            if verbose:
                print("    Processing {}...".format(
                    path.join(directory, intcor)))
            ic_table = ic_file.read()
            ic_table.set_index(_header, inplace=True)
        with reader(path.join(directory, parmfile)) as prm_file:
            if verbose:
                print("    Processing {}...".format(
                    path.join(directory, parmfile)))
            prm_table = prm_file.read()["BONDS"].set_index(_header)
        table = pd.concat([ic_table, prm_table], axis=1)
        table.reset_index(inplace=True)
        table = table.set_index(_index["general"])[tbltype].to_frame()
        if trimmed:
            table.columns = [
                path.basename(path.dirname(directory)),
            ]
        else:
            table.columns = [
                path.basename(directory),
            ]
        return table


class ParamTable(object):
    """Create a parameter table time series for distance or coupling strength.

    """

    def __init__(self,
                 prefix="fluctmatch",
                 tbltype="Kb",
                 ressep=3,
                 datadir=path.curdir,
                 start=None,
                 end=None,
                 trimmed=False,
                 trimcut=None):
        """
        Parameters
        ----------
        prefix : str, optional
            Filename prefix for files
        tbltype : {"Kb", "b0"}, optional
            Table to create (coupling strength or bond distance)
        ressep : int, optional
            Number of residues to exclude from interactions.
        datadir : str, optional
            Directory with data subdirectories
        start : str, optional
            Start data directory
        end : str, optional
            End data directory
        trimmed : bool, optional
            Look for trimmed folder
        trimcut : float, optional
            trimmed flag requires trim cutoff value
        """
        self._prefix = prefix
        self._tbltype = tbltype
        self._ressep = ressep
        self._datadir = datadir
        self.table = []
        self._filenames = dict(
            intcor="fluct.ic",
            param=".".join((self._prefix, "dist", "prm")),
        )
        self.start = start
        self.end = end
        self.trimmed = trimmed
        self.trimcut = trimcut

    def __add__(self, other):
        return self.table.add(other.table, fill_value=0.0)

    def __sub__(self, other):
        return self.table.subtract(other.table, fill_value=0.0)

    def __mul__(self, other):
        return self.table.multiply(other.table, fill_value=0.0)

    def __truediv__(self, other):
        return self.table.divide(other.table, fill_value=0.0)

    def _separate(self, prm_table):
        index = prm_table.index.names
        table = prm_table.reset_index()
        tmp = table[table["segidI"] == table["segidJ"]]
        tmp = tmp[(tmp["resI"] >= tmp["resJ"] + self._ressep)
                  | (tmp["resJ"] >= tmp["resI"] + self._ressep)]
        diff = table[table["segidI"] != table["segidJ"]]
        table = pd.concat([tmp, diff], axis=0)
        table.set_index(index, inplace=True)
        return table

    def _complete_table(self):
        """Create a full table by reversing I and J designations.

        Returns
        -------

        """
        revcol = ["segidJ", "resJ", "J", "segidI", "resI", "I"]

        columns = np.concatenate((revcol, self.table.columns[len(revcol):]))
        temp = self.table.copy(deep=True)
        same = temp[(temp["segidI"] == temp["segidJ"])
                    & (temp["resI"] != temp["resJ"])]

        diff = temp[temp["segidI"] != temp["segidJ"]]
        temp = pd.concat([same, diff], axis=0)
        temp.columns = columns
        self.table = pd.concat([self.table, temp], axis=0)

    def run(self, verbose=False):
        """Create the time series.

        Parameters
        ----------
        verbose : bool, optional
            Print each directory as it is being processed
        """

        if self.trimmed:
            if self.trimcut is None:
                raise ValueError('trimmed flag requires trimcut, please check ')
            directories = glob.iglob(path.join(self._datadir, f"*/trimmed_{self.trimcut}"))
        else:
            directories = glob.iglob(path.join(self._datadir, "*"))

        if self.start is not None and self.end is not None:
            directories = []
            for dirs in glob.iglob(path.join(self._datadir, "*")):
                basename = path.basename(dirs)
                try:
                    number = int(basename)
                except ValueError:
                    continue  # not numeric

                if self.start <= number <= self.end:
                    # process file
                    if self.trimmed:
                        pathbase = path.join(str(number), f"trimmed_{self.trimcut}")
                        directories.append(path.join(self._datadir, pathbase))
                    else:
                        directories.append(path.join(self._datadir, str(number)))

        create_table = functools.partial(
            _create_table,
            intcor=self._filenames["intcor"],
            parmfile=self._filenames["param"],
            tbltype=self._tbltype,
            trimmed=self.trimmed,
            verbose=verbose,
        )
        pool = mp.Pool()
        tables = pool.map_async(create_table, directories)
        pool.close()
        pool.join()
        tables.wait()

        self.table = pd.concat(tables.get(), axis=1)
        self.table.columns = self.table.columns.astype(np.int)
        self.table = self.table[np.sort(self.table.columns)]
        self.table.reset_index(inplace=True)

        self._complete_table()
        self.table.set_index(_index["general"], inplace=True)
        print(self.table.shape)
        # Nix, drop values with 0.0 in all columns of a row, Finalized to use hereafter
        self.table = self.table.where(self.table != 0.).dropna(how="all")
        print(self.table.shape)
        self.table.fillna(0., inplace=True)
        self.table.sort_index(kind="mergesort", inplace=True)

    def from_file(self, filename):
        """Load a parameter table from a file.

        Parameters
        ----------
        filename : str or stream
            Filename of the parameter table.
        """
        with open(filename, mode="rb") as table:
            self.table = pd.read_csv(
                table,
                skipinitialspace=True,
                delim_whitespace=True,
                header=0,
            )
            if "resnI" in self.table.columns:
                self.table.set_index(_index["complete"], inplace=True)
            else:
                self.table.set_index(_index["general"], inplace=True)
            self.table.columns = self.table.columns.astype(np.int)

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
                encoding="utf-8",
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
        table = table.reindex(
            index=table.index, columns=np.sort(table.columns))
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
        table = table.reindex(
            index=table.index, columns=np.sort(table.columns))
        return table

