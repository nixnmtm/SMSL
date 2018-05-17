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
    dict,
    next,
    open,
)
from future.utils import (
    native_str, )

import logging
import time
from io import TextIOWrapper
from os import environ

import numpy as np
import pandas as pd
from MDAnalysis.lib import util
from fluctmatch.topology.base import (TopologyReaderBase, TopologyWriterBase)

logger = logging.getLogger(__name__)


class IntcorReader(TopologyReaderBase):
    """
    Parameters
    ----------
    filename : str or :class:`~MDAnalysis.lib.util.NamedStream`
         name of the output file or a stream
    """
    format = "IC"
    units = dict(time=None, length="Angstrom")

    fmt = dict(
        # fortran_format = "(I5,1X,4(I3,1X,A4),F9.4,3F8.2,F9.4)"
        STANDARD=("I5,1X,I3,1X,A4,I3,1X,A4,I3,1X,A4,I3,1X,A4,F9.4,3F8.2,F9.4"),
        # fortran_format = "(I9,1X,4(I5,1X,A8),F9.4,3F8.2,F9.4)"
        EXTENDED=("I9,1X,I5,1X,A8,I5,1X,A8,I5,1X,A8,I5,1X,A8,F9.4,3F8.2,F9.4"),
        # fortran_format = "(I5,4(1X,A4,1X,A4,1X,A4,"":""),F12.6,3F12.4,F12.6)"
        STANDARD_RESID=(
            "I5,1X,A4,1X,A4,1X,A4,A1,1X,A4,1X,A4,1X,A4,A1,1X,A4,1X,A4,1X,A4,A1,"
            "1X,A4,1X,A4,1X,A4,A1,F12.6,3F12.4,F12.6"),
        # fortran_format = "(I10,4(1X,A8,1X,A8,1X,A8,"":""),F12.6,3F12.4,F12.6)"
        EXTENDED_RESID=(
            "I10,1X,A8,1X,A8,1X,A8,A1,1X,A8,1X,A8,1X,A8,A1,1X,A8,1X,A8,1X,A8,"
            "A1,1X,A8,1X,A8,1X,A8,A1,F12.6,3F12.4,F12.6"),
    )
    cols = np.asarray([
        "segidI", "resI", "I", "segidJ", "resJ", "J", "segidK", "resK", "K",
        "segidL", "resL", "L", "r_IJ", "T_IJK", "P_IJKL", "T_JKL", "r_KL"
    ])

    def __init__(self, filename):
        self.filename = util.filename(filename, ext="ic")

    def read(self):
        """Read the internal coordinates file.

        Returns
        -------
        :class:`~pandas.DataFrame`
            An internal coordinates table.
        """
        table = pd.DataFrame()
        with open(self.filename, "rb") as icfile, TextIOWrapper(
                icfile, encoding="utf-8") as buf:
            logger.info("Writing to {}".format(self.filename))
            for line in buf:
                line = line.split("!")[0].strip()
                if line.startswith("*") or not line:
                    continue  # ignore TITLE and empty lines
                break
            line = np.fromiter(line.strip().split(), dtype=np.int)
            key = "EXTENDED" if line[0] == 30 else "STANDARD"
            key += "_RESID" if line[1] == 2 else ""
            resid_a = line[1]

            line = next(buf).strip().split()
            n_lines, resid_b = np.array(line, dtype=np.int)
            if resid_a != resid_b:
                raise IOError(
                    "A mismatch has occurred on determining the IC format.")

            TableParser = util.FORTRANReader(self.fmt[key])
            table = pd.DataFrame(
                [TableParser.read(line) for line in buf], dtype=np.object)
            table = table[table != ":"]
            table = table.dropna(axis=1).apply(pd.to_numeric, errors="ignore")
            table.set_index(0, inplace=True)
            if n_lines != table.shape[0]:
                raise IOError("A mismatch has occurred between the number "
                              "of lines expected and the number of lines "
                              "read. ({:d} != {:d})".format(
                                  n_lines, len(table)))

            if key == "STANDARD":
                idx = np.where(
                    (self.cols != "segidI") & (self.cols != "segidJ") &
                    (self.cols != "segidK") & (self.cols != "segidL"))
                columns = self.cols[idx]
            else:
                columns = self.cols
            table.columns = columns
            logger.info("Table read successfully.")
        return table


class IntcorWriter(TopologyWriterBase):
    """Write a CHARMM-formatted internal coordinate file.

    Parameters
    ----------
    filename : str
        Filename for output.
    n_atoms : int, optional
        The number of atoms in the output trajectory.
    extended : bool, optional
        Format with wider columns than the standard column width.
    resid : bool, optional
        Include segment names within each atom definition.
    title : str or list of str, optional
        A header section written at the beginning of the stream file.
        If no title is given, a default title will be written.
    """
    format = "IC"
    units = {"time": "picosecond", "length": "Angstrom"}

    fmt = dict(
        # fortran_format = "(I5,1X,4(I3,1X,A4),F9.4,3F8.2,F9.4)"
        STANDARD=(
            "%5d %3s %-4s%3s %-4%3s %-4%3s %-4%9.4f%8.2f%8.2f%8.2f%9.4f"),
        # fortran_format = "(I9,1X,4(I5,1X,A8),F9.4,3F8.2,F9.4)"
        EXTENDED=(
            "%10d %5s %-8s%5s %-8s%5s %-8s%5s %-8s%9.4f%8.2f%8.2f%8.2f%9.4f"),
        # fortran_format = "(I5,4(1X,A4,1X,A4,1X,A4,"":""),F12.6,3F12.4,F12.6)"
        STANDARD_RESID=("%5d %-4s %-4s %-4s: %-4s %-4s %-4s: %-4s %-4s %-4s: "
                        "%-4s %-4s %-4s:%12.6f%12.4f%12.4f%12.4f%12.6f"),
        # fortran_format = "(I10,4(1X,A8,1X,A8,1X,A8,"":""),F12.6,3F12.4,F12.6)"
        EXTENDED_RESID=("%10d %-8s %-8s %-8s: %-8s %-8s %-8s: %-8s %-8s %-8s: "
                        "%-8s %-8s %-8s:%12.6f%12.4f%12.4f%12.4f%12.6f"),
    )

    def __init__(self, filename, **kwargs):
        self.filename = util.filename(filename, ext="ic")
        self._intcor = None
        self._extended = kwargs.get("extended", True)
        self._resid = kwargs.get("resid", True)
        self.key = "EXTENDED" if self._extended else "STANDARD"
        self.key += "_RESID" if self._resid else ""

        date = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        user = environ["USER"]
        self._title = kwargs.get(
            "title", (
                "* Created by fluctmatch on {date}".format(date=date),
                "* User: {user}".format(user=user),
            ))
        if not util.iterable(self._title):
            self._title = util.asiterable(self._title)

    def write(self, table):
        """Write an internal coordinates table.

        Parameters
        ----------
        table : :class:`~pandas.DataFrame`
            A CHARMM-compliant internal coordinate table.
        """
        ictable = table.copy(deep=True)

        # Increment index.
        if ictable.index[0] == 0:
            ictable.index += 1
        rescol = [
            "resI",
            "resJ",
            "resK",
            "resL",
        ]
        ictable[rescol] = ictable[rescol].astype(np.unicode)

        with open(self.filename, "wb") as icfile:
            logger.info("Writing to {}".format(self.filename))
            for _ in self._title:
                icfile.write(_.encode())
                icfile.write("\n".encode())
            line = np.zeros(20, dtype=np.int)
            line[0] = 30 if self._extended else 20
            line[1] = 2 if self._resid else 1
            np.savetxt(
                icfile,
                line[np.newaxis, :],
                fmt=native_str("%4d"),
                delimiter=native_str(""))
            line = np.zeros(2, dtype=np.int)
            line[0] = ictable.shape[0]
            line[1] = 2 if self._resid else 1
            np.savetxt(
                icfile,
                line[np.newaxis, :],
                fmt=native_str("%5d"),
                delimiter=native_str(""))
            np.savetxt(
                icfile,
                ictable.reset_index(),
                fmt=native_str(self.fmt[self.key]))
            logger.info("Table successfully written.")
