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
    open,
)
from future.utils import (
    native_str, )

import textwrap
import time
from io import (StringIO, TextIOWrapper)
from os import environ

import numpy as np
import pandas as pd
from MDAnalysis.lib import util
from fluctmatch.topology.base import (TopologyReaderBase, TopologyWriterBase)


class ParamReader(TopologyReaderBase):
    """Read a CHARMM-formated parameter file.

    Parameters
    ----------
    filename : str or :class:`~MDAnalysis.lib.util.NamedStream`
         name of the output file or a stream
    """
    format = "PRM"
    units = dict(time=None, length="Angstrom")

    _prmindex = dict(
        ATOMS=np.arange(1, 4),
        BONDS=np.arange(4),
        ANGLES=np.arange(7),
        DIHEDRALS=np.arange(6))
    _prmcolumns = dict(
        ATOMS=["hdr", "type", "atom", "mass"],
        BONDS=["I", "J", "Kb", "b0"],
        ANGLES=["I", "J", "K", "Ktheta", "theta0", "Kub", "S0"],
        DIHEDRALS=["I", "J", "K", "L", "Kchi", "n", "delta"],
        IMPROPER=["I", "J", "K", "L", "Kchi", "n", "delta"])
    _dtypes = dict(
        ATOMS=dict(hdr=np.str, type=np.int, atom=np.str, mass=np.float,),
        BONDS=dict(I=np.str, J=np.str, Kb=np.float, b0=np.float),
        ANGLES=dict(
            I=np.str, J=np.str, K=np.str, Ktheta=np.float, theta0=np.float,
            Kub=np.object, S0=np.object,
        ),
        DIHEDRALS=dict(
            I=np.str, J=np.str, K=np.str, L=np.str,
            Kchi=np.float, n=np.int, delta=np.float,
        ),
        IMPROPER=dict(
            I=np.str, J=np.str, K=np.str, L=np.str,
            Kchi=np.float, n=np.int, delta=np.float,
        ),
    )
    _na_values = dict(
        ATOMS=dict(type=-1, mass=0.0),
        BONDS=dict(Kb=0.0, b0=0.0),
        ANGLES=dict(Ktheta=0.0, theta0=0.0, Kub="", S0=""),
        DIHEDRALS=dict(Kchi=0.0, n=0, delta=0.0),
        IMPROPER=dict(Kchi=0.0, n=0, delta=0.0),
    )

    def __init__(self, filename):
        self.filename = util.filename(filename, ext="prm")
        self._prmbuffers = dict(
            ATOMS=StringIO(),
            BONDS=StringIO(),
            ANGLES=StringIO(),
            DIHEDRALS=StringIO(),
            IMPROPER=StringIO(),
        )

    def read(self):
        """Parse the parameter file.

        Returns
        -------
        Dictionary with CHARMM parameters per key.
        """
        parameters = dict.fromkeys(self._prmbuffers.keys())
        headers = ("ATOMS", "BONDS", "ANGLES", "DIHEDRALS", "IMPROPER")
        with open(self.filename, "rb") as prmfile, TextIOWrapper(
                prmfile, encoding="utf-8") as buf:
            for line in buf:
                line = line.strip()
                if line.startswith("*") or line.startswith("!") or not line:
                    continue  # ignore TITLE and empty lines

                # Parse sections
                if line in headers:
                    section = line
                    continue
                elif (line.startswith("NONBONDED") or
                      line.startswith("CMAP") or
                      line.startswith("END") or
                      line.startswith("end")):
                    break

                print(line, file=self._prmbuffers[section])

        for key, value in parameters.items():
            self._prmbuffers[key].seek(0)
            parameters[key] = pd.read_csv(
                self._prmbuffers[key], header=None, names=self._prmcolumns[key],
                skipinitialspace=True, delim_whitespace=True, comment="!",
                dtype=self._dtypes[key]
            )
            parameters[key].fillna(self._na_values[key], inplace=True)
        if not parameters["ATOMS"].empty:
            parameters["ATOMS"].drop("hdr", axis=1, inplace=True)
        return parameters


class PARReader(ParamReader):
    format = "PAR"

    def __init__(self, filename):
        super().__init__(filename)
        self.filename = util.filename(filename, ext="par")


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
        ANGLES="%-6s %-6s %-6s %8.2f%10.2f%10s%10s",
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
            "title", (
                "* Created by fluctmatch on {date}".format(date=date),
                "* User: {user}".format(user=user),
            ))
        if not util.iterable(self._title):
            self._title = util.asiterable(self._title)

    def write(self, parameters, atomgroup=None):
        """Write a CHARMM-formatted parameter file.

        Parameters
        ----------
        parameters : dict
            Keys are the section names and the values are of class
            :class:`~pandas.DataFrame`, which contain the corresponding parameter data.
        atomgroup : :class:`~MDAnalysis.AtomGroup`, optional
            A collection of atoms in an AtomGroup to define the ATOMS section,
            if desired.
        """
        with open(self.filename, "wb") as prmfile:
            for title in self._title:
                prmfile.write(title.encode())
                prmfile.write("\n".encode())
            prmfile.write("\n".encode())

            if self._version > 35 and parameters["ATOMS"].empty:
                if atomgroup:
                    if np.issubdtype(atomgroup.types.dtype, np.int):
                        atom_types = atomgroup.types
                    else:
                        atom_types = np.arange(atomgroup.n_atoms) + 1
                    atoms = [atom_types, atomgroup.types, atomgroup.masses]
                    parameters["ATOMS"] = pd.concat(
                        [pd.Series(_) for _ in atoms], axis=1)
                    parameters["ATOMS"].columns = ["type", "atom", "mass"]
                else:
                    raise RuntimeError(
                        "Either define ATOMS parameter or provide a "
                        "MDAnalsys.AtomGroup")

            if self._version >= 39 and not parameters["ATOMS"].empty:
                parameters["ATOMS"]["type"] = -1

            for key in self._headers:
                value = parameters[key]
                if self._version < 35 and key == "ATOMS":
                    continue
                prmfile.write(key.encode())
                prmfile.write("\n".encode())
                if value.empty:
                    prmfile.write("\n".encode())
                if not value.empty:
                    np.savetxt(prmfile, value, fmt=native_str(self._fmt[key]))
                    prmfile.write("\n".encode())

            nb_header = ("""
                NONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch -
                cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5
                """)
            prmfile.write(textwrap.dedent(nb_header[1:]).encode())

            if self._nonbonded:
                atom_list = np.concatenate(
                    (parameters["BONDS"]["I"].values,
                     parameters["BONDS"]["J"].values),
                    axis=0,
                )
                atom_list = pd.DataFrame(np.unique(atom_list))
                nb_list = pd.DataFrame(np.zeros((atom_list.size, 3)))
                nb_list = pd.concat([atom_list, nb_list], axis=1)
                np.savetxt(
                    prmfile,
                    nb_list,
                    fmt=native_str(self._fmt["NONBONDED"]),
                    delimiter=native_str(""))
            prmfile.write("\nEND\n".encode())


class PARWriter(ParamWriter):
    format = "PAR"
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)
        self.filename = util.filename(filename, ext="par")
