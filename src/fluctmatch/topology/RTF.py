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

import time
from io import TextIOWrapper
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

from fluctmatch.topology import base as topbase


class RTFWriter(topbase.TopologyWriterBase):
    """Write a CHARMM-formatted topology file.

    Parameters
    ----------
    filename : str
        Filename where to write the information.
    n_atoms : int, optional
        The number of atoms in the output trajectory.
    title
        A header section written at the beginning of the stream file.
        If no title is given, a default title will be written.
    charmm_version
        Version of CHARMM for formatting (default: 41)
    """
    format = "RTF"
    units = dict(time=None, length=None)
    fmt = dict(
        HEADER="{:>5d}{:>5d}\n",
        MASS="MASS %5d %-6s%12.5f",
        DECL="DECL +%s\nDECL -%s",
        RES="RESI {:<4s} {:>12.4f}\nGROUP\n",
        ATOM="ATOM %-6s %-6s %7.4f",
        IC="IC %-4s %-4s %-4s %-4s %7.4f %8.4f %9.4f %8.4f %7.4f",)
    bonds = (
        ("BOND", ("bonds", 8)),
        ("IMPH", ("impropers", 8)),
    )

    def __init__(self, filename, **kwargs):
        self.filename = util.filename(filename, ext="rtf")
        self._version = kwargs.get("charmm_version", 41)
        self._atoms = None

        date = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        user = environ["USER"]
        self._title = kwargs.get(
            "title",
            (
                "* Created by fluctmatch on {date}".format(date=date),
                "* User: {user}".format(user=user),
            )
        )
        if not util.iterable(self._title):
            self._title = util.asiterable(self._title)

    def _write_mass(self):
        _, idx = np.unique(self._atoms.names, return_index=True)
        try:
            atomtypes = self._atoms.types.astype(np.int)
        except ValueError:
            atomtypes = np.arange(self._atoms.n_atoms, dtype=np.int) + 1
        columns = (atomtypes, self._atoms.names, self._atoms.masses)
        columns = pd.concat([pd.DataFrame(_[idx]) for _ in columns], axis=1)
        columns.columns = ["itype", "stype", "mass"]

        if self._version >= 39:
            columns["itype"] = -1
        np.savetxt(
            self.rtffile,
            columns,
            fmt=native_str(self.fmt["MASS"]),
            delimiter=native_str("")
        )

    def _write_decl(self):
        names = np.unique(self._atoms.names)[:, np.newaxis]
        decl = np.concatenate((names, names), axis=1)
        np.savetxt(self.rtffile, decl, fmt=native_str(self.fmt["DECL"]))
        self.rtffile.write("\n".encode())

    def _write_residues(self, residue):
        self.rtffile.write(
            self.fmt["RES"].format(residue.resname, residue.charge).encode()
        )

        # Write the atom lines with site name, type, and charge.
        key = "ATOM"
        atoms = residue.atoms
        lines = (
            (atoms.names, atoms.types, atoms.charges)
            if not np.issubdtype(atoms.types.dtype, np.int)
            else (atoms.names, atoms.names, atoms.charges)
        )
        lines = pd.concat([pd.Series(_) for _ in lines], axis=1)
        np.savetxt(self.rtffile, lines, fmt=native_str(self.fmt[key]))

        # Write the bond, angle, dihedral, and improper dihedral lines.
        for key, value in self.bonds:
            attr, n_perline = value
            fmt = key + n_perline * "%6s"
            try:
                bonds = getattr(atoms, attr)
                if len(bonds) == 0:
                    continue

                # Create list of atom names and include "+" for atoms not
                # within the residue.
                names = np.concatenate([
                    _.atoms.names[np.newaxis, :]
                    for _ in bonds
                ], axis=0).astype(np.object)
                idx = np.any([
                    np.isin(_.atoms, atoms, invert=True)
                    for _ in bonds
                ], axis=1)
                pos_names = np.where(
                    np.isin(bonds[idx], atoms, invert=True),
                    "+", ""
                ).astype(np.object)
                names[idx] = pos_names + names[idx]
                names = names.astype(np.unicode)

                # Eliminate redundancies.
                # Code courtesy of Daniel F on
                # https://stackoverflow.com/questions/45005477/eliminating-redundant-numpy-rows/45006131?noredirect=1#comment76988894_45006131
                b = np.ascontiguousarray(
                    np.sort(names, -1)
                ).view(
                    np.dtype((np.void, names.dtype.itemsize * names.shape[1]))
                )
                _, idx = np.unique(b, return_index=True)
                names = names[idx]

                # Add padding for missing columns.
                n_rows, n_cols = names.shape
                n_values = n_perline // n_cols
                if n_rows % n_values > 0:
                    n_extra = n_values - (n_rows % n_values)
                    extras = np.full((n_extra, n_cols), native_str(""))
                    names = np.concatenate((names, extras), axis=0)
                names = names.reshape((names.shape[0] // n_values, n_perline))
                np.savetxt(self.rtffile, names, fmt=native_str(fmt))
            except (AttributeError,):
                continue
        self.rtffile.write("\n".encode())

    def write(self, universe):
        """Write a CHARMM-formatted RTF topology file.

        Parameters
        ----------
        universe : :class:`~MDAnalysis.Universe` or :class:`~MDAnalysis.AtomGroup`
            A collection of atoms in a universe or atomgroup with bond
            definitions.
        """
        self._atoms = universe.atoms
        with open(
            self.filename, "wb"
        ) as self.rtffile:
            # Write the title and header information.
            for _ in self._title:
                self.rtffile.write(_.encode())
                self.rtffile.write("\n".encode())
            self.rtffile.write(self.fmt["HEADER"].format(36, 1).encode())
            self.rtffile.write("\n".encode())

            # Write the atom mass and declaration sections
            self._write_mass()
            self.rtffile.write("\n".encode())
            self._write_decl()
            self.rtffile.write("DEFA FIRS NONE LAST NONE\n".encode())
            self.rtffile.write("AUTOGENERATE ANGLES DIHEDRAL\n\n".encode())

            # Write out the residue information
            _, idx = np.unique(self._atoms.residues.resnames, return_index=True)
            for residue in self._atoms.residues[idx]:
                self._write_residues(residue)
            self.rtffile.write("END\n".encode())
