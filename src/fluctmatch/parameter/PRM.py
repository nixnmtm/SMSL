# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import textwrap
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
    bytes_to_native_str,
    native_str,
)

from ..topology.base import (TopologyReaderBase, TopologyWriterBase)


class ParamReader(TopologyReaderBase):
    """Read a CHARMM-formated parameter file.

    Parameters
    ----------
    filename : str or :class:`~MDAnalysis.lib.util.NamedStream`
         name of the output file or a stream
    """
    format = "PRM"
    units = dict(time=None, length="Angstrom")

    parameters = dict(ATOMS=[], BONDS=[], ANGLES=[], DIHEDRALS=[], IMPROPER=[])
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
        self.filename = util.filename(filename, ext="prm")

    def read(self):
        """Parse the parameter file.

        Returns
        -------
        Dictionary with CHARMM parameters per key.
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
                if section == "BONDS":
                    field = field[:4]
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
    """Write a parameter dictionary to a CHARMM-formatted parameter file.

    Parameters
    ----------
    filename : str or :class:`~MDAnalysis.lib.util.NamedStream`
         name of the output file or a stream
    title : str
        Title lines at beginning of the file.
    charmm_version
        Version of CHARMM for formatting (default: 41)
    nonbonded
        Add the nonbonded section. (default: False)
    """
    format = "PRM"
    units = dict(time=None, length="Angstrom")

    _headers = ("ATOMS", "BONDS", "ANGLES", "DIHEDRALS", "IMPROPER")
    _fmt = dict(
        ATOMS="MASS %5d %-6s %9.5f",
        BONDS="%-6s %-6s %10.4f%10.4f",
        ANGLES="%-6s %-6s %-6s %8.2f%10.2f",
        DIHEDRALS="%-6s %-6s %-6s %-6s %12.4f%3d%9.2f",
        IMPROPER="%-6s %-6s %-6s %-6s %12.4f%3d%9.2f",
        NONBONDED="%-6s %5.1f %13.4f %10.4f",
    )

    def __init__(self, filename, **kwargs):
        self.filename = util.filename(filename, ext="prm")
        self._version = kwargs.get("charmm_version", 41)
        self._nonbonded = kwargs.get("nonbonded", False)

        date = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        user = environ["USER"]
        self._title = kwargs.get(
            "title",
            (
                "* Created by fluctmatch on {date}".format(date=date),
                "* User: {user}".format(user=user),
            )
        )
        if issubclass(type(self._title), str) or issubclass(type(self._title), np.unicode):
            self._title = (self._title,)

    def write(self, parameters, atomgroup=None):
        """Write a CHARMM-formatted parameter file.

        Parameters
        ----------
        parameters : dict
            Keys are the section names and the values are of class :class:`~pandas.DataFrame`,
            which contain the corresponding parameter data.
        atomgroup : :class:`~MDAnalysis.AtomGroup`, optional
            A collection of atoms in an AtomGroup to define the ATOMS section, if desired.
        """
        with open(self.filename, "wb") as prmfile:
            for title in self._title:
                prmfile.write((title + "\n").encode())
            prmfile.write("\n".encode())

            if self._version >= 36 and parameters["ATOMS"].empty:
                if atomgroup:
                    if np.issubdtype(atomgroup.types.dtype, np.int):
                        atom_types = atomgroup.types
                    else:
                        atom_types = np.arange(atomgroup.n_atoms)+1
                    atoms = [atom_types, atomgroup.types, atomgroup.masses]
                    parameters["ATOMS"] = pd.concat([pd.Series(_) for _ in atoms], axis=1)
                    parameters["ATOMS"].columns = ["type", "atom", "mass"]
                else:
                    raise RuntimeError("Either define ATOMS parameter or provide a MDAnalsys.AtomGroup")

            if self._version >= 39 and not parameters["ATOMS"].empty:
                parameters["ATOMS"]["type"] = -1

            for key in self._headers:
                value = parameters[key]
                if value.empty or (self._version < 36 and key == "ATOMS"):
                    continue
                prmfile.write((key + "\n").encode())
                np.savetxt(prmfile, value, fmt=native_str(self._fmt[key]), delimiter=native_str(""))
                prmfile.write("\n".encode())

            if self._nonbonded:
                nb_header = (
                    """
                    NONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch -
                    cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5

                    """
                )
                atom_list = np.concatenate(
                    (parameters["BONDS"]["I"].values, parameters["BONDS"]["J"].values),
                    axis=0,
                )
                atom_list = pd.DataFrame(np.unique(atom_list))
                nb_list = pd.DataFrame(np.zeros((atom_list.size, 3)))
                nb_list = pd.concat([atom_list, nb_list], axis=1)
                prmfile.write(textwrap.dedent(nb_header).encode())
                np.savetxt(prmfile, nb_list, fmt=native_str(self._fmt["NONBONDED"]), delimiter=native_str(""))
            prmfile.write("\nEND\n".encode())
