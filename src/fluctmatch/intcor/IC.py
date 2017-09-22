# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import time
from os import environ

import numpy as np
import pandas as pd
from MDAnalysis.lib import util
from future.builtins import (
    dict,
    open,
)
from future.utils import (
    native_str,
)

from ..topology.base import (TopologyReaderBase, TopologyWriterBase)


class IntcorReader(TopologyReaderBase):
    format = "IC"
    units = dict(time=None, length="Angstrom")

    fmt = dict(
        # fortran_format = "(I5,1X,4(I3,1X,A4),F9.4,3F8.2,F9.4)"
        STD = "I5,1X,I3,1X,A4,I3,1X,A4,I3,1X,A4,I3,1X,A4,F9.4,3F8.2,F9.4",
        # fortran_format = "(I9,1X,4(I5,1X,A8),F9.4,3F8.2,F9.4)"
        STD_EXT = ("I9,1X,I5,1X,A8,I5,1X,A8,I5,1X,A8,I5,1X,A8,F9.4,"
                   "3F8.2,F9.4"),
        # fortran_format = "(I5,4(1X,A4,1X,A4,1X,A4,"":""),F12.6,3F12.4,F12.6)"
        RESID = ("I5,1X,A4,1X,A4,1X,A4,A1,1X,A4,1X,A4,1X,A4,A1,1X,A4,"
                 "1X,A4,1X,A4,A1,1X,A4,1X,A4,1X,A4,A1,F12.6,3F12.4,F12.6"),
        # fortran_format = "(I10,4(1X,A8,1X,A8,1X,A8,"":""),F12.6,3F12.4,F12.6)"
        RESID_EXT = ("I10,1X,A8,1X,A8,1X,A8,A1,1X,A8,1X,A8,1X,A8,A1,1X,"
                      "A8,1X,A8,1X,A8,A1,1X,A8,1X,A8,1X,A8,A1,F12.6,"
                      "3F12.4,F12.6"),
    )
    cols = np.asarray(["segidI", "resI", "I", "segidJ", "resJ", "J",
                       "segidK", "resK", "K", "segidL", "resL", "L",
                       "r_IJ", "T_IJK", "P_IJKL", "T_JKL", "r_KL"],)

    def __init__(self, filename):
        """
        Parameters
        ----------
        filename : str or :class:`~MDAnalysis.lib.util.NamedStream`
             name of the output file or a stream
        """

        self.filename = util.filename(filename, ext="ic")

    def read(self):
        """Read the internal coordinates file.

        :return: pandas.DataFrame
        """
        with open(self.filename, "rb") as icfile:
            for line in icfile:
                line = bytes_to_native_str(line).split("!")[0].strip()
                if line.startswith("*") or not line:
                    continue       # ignore TITLE and empty lines
                line = np.fromiter(line.strip().split(), dtype=np.int)
                key = "RESID" if line[1] == 2 else "STD"
                if line[0] == 30:
                    key += "_EXT"
                break

            TableParser = util.FORTRANReader(self.fmt[key])
            nlines, resid = np.array(icfile.readline().strip().split(), dtype=np.int)
            if resid != line[1]:
                raise IOError("A mismatch has occurred on determining the IC format.")

            table = pd.DataFrame([TableParser.read(bytes_to_native_str(line)) for line in icfile], dtype=np.unicode)
            table = table[table != ":"].dropna(axis=1).apply(pd.to_numeric, errors="ignore")
            table.set_index(0, inplace=True)
            if nlines != table.shape[0]:
                raise IOError(("A mismatch has occurred between the number "
                               "of lines expected and the number of lines "
                               "read. (%d != %d)") % (nlines, len(table)))

            if key == "STD":
                idx = np.where((self.cols != "segidI") & (self.cols != "segidJ") &
                               (self.cols != "segidK") & (self.cols != "segidL"))
                columns = self.cols[idx]
            else:
                columns = self.cols
            table.columns = columns
        return table


class IntcorWriter(TopologyWriterBase):
    format = "IC"
    units = {"time": "picosecond", "length": "Angstrom"}

    fmt = dict(
        # fortran_format = "(I5,1X,4(I3,1X,A4),F9.4,3F8.2,F9.4)"
        STD="%5d %3s %-4s%3s %-4%3s %-4%3s %-4%9.4f%8.2f%8.2f%8.2f%9.4f",
        # fortran_format = "(I9,1X,4(I5,1X,A8),F9.4,3F8.2,F9.4)"
        STD_EXT="%9d %5s %-8s%5s %-8s%5s %-8s%5s %-8s%9.4f%8.2f%8.2f%8.2f%9.4f",
        # fortran_format = "(I5,4(1X,A4,1X,A4,1X,A4,"":""),F12.6,3F12.4,F12.6)"
        RESID=("%5d %-4s %-4s %-4s: %-4s %-4s %-4s: %-4s %-4s %-4s:"
               " %-4s %-4s %-4s:%12.6f%12.4f%12.4f%12.4f%12.6f"),
        # fortran_format = "(I10,4(1X,A8,1X,A8,1X,A8,"":""),F12.6,3F12.4,F12.6)"
        RESID_EXT=("%10d %-8s %-8s %-8s: %-8s %-8s %-8s: "
                   "%-8s %-8s %-8s: %-8s %-8s %-8s:"
                   "%12.6f%12.4f%12.4f%12.4f%12.6f"),)

    def __init__(self, filename, extended=True, resid=True, title=None):
        self.filename = util.filename(filename, ext="ic")
        self.intcor = None
        self.extended = extended
        self.resid = resid
        self.key = "RESID" if self.resid else "STD"
        if self.extended:
            self.key += "_EXT"
        self.title = ("* Created by fluctmatch on {}".format(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())),
                      "* User: {}".format(environ["USER"]), "*") if title is None else title

    def write(self, table):
        # type: (object) -> object
        """Write an internal coordinates table.

        :param table: table of internal coordinates
        """
        ictable = table.copy(deep=True)
        rescol = ["resI", "resJ", "resK", "resL", ]
        ictable[rescol] = ictable[rescol].astype(np.unicode)

        with open(self.filename, "wb") as icfile:
            for _ in self.title:
                icfile.write((_ + "\n").encode())
            line = np.zeros(20, dtype=np.int)
            line[0] = 30 if self.extended else 20
            line[1] = 2 if self.resid else 1
            np.savetxt(icfile, line[np.newaxis, :], fmt=native_str("%4d"), delimiter=native_str(""))
            line = np.zeros(2, dtype=np.int)
            line[0] = ictable.shape[0]
            line[1] = 2 if self.resid else 1
            np.savetxt(icfile, line[np.newaxis, :], fmt=native_str("%4d"), delimiter=native_str(""))
            np.savetxt(icfile, ictable.reset_index(), fmt=native_str(self.fmt[self.key]))
