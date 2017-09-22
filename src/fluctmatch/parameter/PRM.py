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
from collections import OrderedDict
from os import environ

import numpy as np
import pandas as pd
from MDAnalysis.lib import util
from future.builtins import (
    dict,
    open,
)
from future.utils import (
    bytes_to_native_str,
    native_str,
)

from ..topology.base import (TopologyReaderBase, TopologyWriterBase)


class ParamReader(TopologyReaderBase):
    format = "PRM"
    units = dict(time=None, length="Angstrom")

    parameters = OrderedDict(ATOMS=[], BONDS=[], ANGLES=[], DIHEDRALS=[], IMPROPER=[])
    _prmindex = dict(
        ATOMS=np.arange(1, 4),
        BONDS=np.arange(4),
        ANGLES=np.arange(5),
        DIHEDRALS=np.arange(6)
    )
    _prmcolumns = dict(
        ATOMS=["type", "atom", "mass"],
        BONDS=["I", "J", "Kb", "b0"],
        ANGLES=["I", "J", "K", "Ktheta", "theta0"],
        DIHEDRALS=["I", "J", "K", "L", "Kchi", "n", "delta"],
        IMPROPER=["I", "J", "K", "L", "Kchi", "n", "delta"]
    )

    def __init__(self, filename):
        """
        Parameters
        ----------
        filename : str or :class:`~MDAnalysis.lib.util.NamedStream`
             name of the output file or a stream
        """
        self.filename = util.filename(filename, ext="prm")

    def read(self):
        """Parse the parameter file.

        :return: parameters
        """
        headers = ("ATOMS", "BONDS", "ANGLES", "DIHEDRALS", "IMPROPER")
        with open(self.filename, "rb") as prmfile:
            for line in prmfile:
                line = bytes_to_native_str(line).split("!")[0].strip()
                if line.startswith("*") or not line:
                    continue       # ignore TITLE and empty lines

                # Parse sections
                if line in headers:
                    section = line
                    continue
                elif line.startswith("NONBONDED") or line.startswith("CMAP") or line.startswith("END"):
                    break

                # Removes any Urey-Bradley values
                field = line.split()
                if section == "ATOMS":
                    field = field[1:]
                if section == "ANGLES":
                    field = field[:5]
                if section == "DIHEDRALS" or section == "IMPROPER":
                    field = field[:7]
                self.parameters[section].append(field)

        for key, value in self.parameters.items():
            self.parameters[key] = pd.DataFrame(value, columns=self._prmcolumns[key])
            self.parameters[key] = self.parameters[key].apply(pd.to_numeric, errors="ignore")
        return self.parameters


class ParamWriter(TopologyWriterBase):
    format = "PRM"
    units = dict(time=None, length="Angstrom")

    _headers = ("ATOMS", "BONDS", "ANGLES", "DIHEDRALS", "IMPROPER")
    _fmt = dict(ATOMS="",
                ATOMS_C36="MASS %5d %-6s %9.5f",
                BONDS="%-4s %-4s %10.4f%10.4f",
                BONDS_C36="%-6s %-6s %10.4f%10.4f",
                ANGLES="%-4s %-4s %-4s %8.2f%10.2f",
                ANGLES_C36="%-6s %-6s %-6s %8.2f%10.2f",
                DIHEDRALS="%-4s %-4s %-4s %-4s %12.4f%3d%9.2f",
                DIHEDRALS_C36="%-6s %-6s %-6s %-6s %12.4f%3d%9.2f",
                IMPROPER="%-4s %-4s %-4s %-4s %12.4f%3d%9.2f",
                IMPROPER_C36="%-6s %-6s %-6s %-6s %12.4f%3d%9.2f",
                NONBONDED="%-4s %5.1f %13.4f %10.4f",
                NONBONDED_C36="%-6s %5.1f %13.4f %10.4f",)

    def __init__(self, filename, title=None, charmm36=True, nonbonded=True):
        """
        Parameters
        ----------
        filename : str or :class:`~MDAnalysis.lib.util.NamedStream`
             name of the output file or a stream
        title : str
             title lines at beginning of the file
        charmm36 : bool
             Use Charmm36 force field representation
        nonbonded : bool
             Add the nonbonded section
        """

        self.filename = util.filename(filename, ext="prm")
        self.charmm36 = charmm36
        self.nonbonded = nonbonded
        self.title = ("* Created by fluctmatch on {}".format(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())),
                      "* User: {}".format(environ["USER"])) if title is None else title

    def write(self, parameters, atomgroup=None):
        # type: (object, object) -> object
        """Write a CHARMM parameter file.

        :param parameters: CHARMM parameters
        :param atomgroup: optional Atomgroup for CHARMM36 atom list
        :return:
        """
        with open(self.filename, "wb") as prmfile:
            for title in self.title:
                prmfile.write((title + "\n").encode())
            prmfile.write("\n".encode())

            if self.charmm36 and parameters["ATOMS"].empty:
                if atomgroup:
                    if np.issubdtype(atomgroup.types.dtype, np.int):
                        atoms = [atomgroup.types, atomgroup.names, atomgroup.masses]
                    else:
                        atoms = [np.arange(atomgroup.n_atoms)+1, atomgroup.types, atomgroup.masses]
                    parameters["ATOMS"] = pd.concat([pd.Series(_) for _ in atoms], axis=1)
                    parameters["ATOMS"].columns = ["type", "atom", "mass"]
                else:
                    raise RuntimeError("Either define ATOMS parameter or provide a MDAnalsys.AtomGroup")

            for key, value in parameters.items():
                if value.empty or (not self.charmm36 and key == "ATOMS"):
                    continue
                fmt_key = key + "_C36" if self.charmm36 else key
                prmfile.write((key + "\n").encode())
                np.savetxt(prmfile, value, fmt=native_str(self._fmt[fmt_key]))
                prmfile.write("\n".encode())

            if self.nonbonded:
                key = "NONBONDED_C36" if self.charmm36 else "NONBONDED"
                n_atoms = parameters["ATOMS"].shape[0] if not parameters["ATOMS"].empty else atomgroup.n_atoms
                nblist = np.zeros((n_atoms, 3))
                if not parameters["ATOMS"].empty:
                    nblist = pd.concat([parameters["ATOMS"]["atom"], pd.DataFrame(nblist)], axis=1)
                else:
                    atomtypes = atomgroup.names if np.issubdtype(atomgroup.types.dtype, np.int) else atomgroup.types
                    nblist = pd.concat([pd.DataFrame(atomtypes), pd.DataFrame(nblist)], axis=1)
                prmfile.write("NONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch -\n".encode())
                prmfile.write("cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5\n\n".encode())
                np.savetxt(prmfile, nblist, fmt=native_str(self._fmt[key]))
            prmfile.write("\nEND\n".encode())
