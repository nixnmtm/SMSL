# -*- coding: utf-8 -*-
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import time
from os import (environ)

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

from . import base


class RTFWriter(base.TopologyWriterBase):
    format = "RTF"
    units = dict(time=None, length=None)
    fmt = dict(HEADER="{:<4d}{:>d}\n",
               MASS="MASS %5d %-4s %9.5f",
               DECL="DECL +%s\nDECL -%s",
               RES="RESI {:<4s} {:>12.2f}\nGROUP\n",
               ATOM="ATOM %-4s %-4s %7.2f",
               ATOM_C36="ATOM %-4s %-6s %7.2f",
               IC="IC %-4s %-4s %-4s %-4s %7.4f %8.4f %9.4f %8.4f %7.4f",)
    bonds = (("BOND", ("bonds", 8)), ("IMPH", ("impropers", 8)),)

    def __init__(self, filename, title=None, charmm39=False, **kwargs):
        self.filename = util.filename(filename, ext="rtf")
        self.title = ("* Created by fluctmatch on {}".format(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())),
                      "* User: {}".format(environ["USER"])) if title is None else title
        self.charmm39 = charmm39

    def _write_mass(self):
        _, idx = np.unique(self.atoms.names, return_index=True)
        try:
            atomtypes = self.atoms.types.astype(np.int)
        except ValueError:
            atomtypes = np.arange(self.atoms.n_atoms, dtype=np.int)+1
        columns = [atomtypes, self.atoms.names, self.atoms.masses]

        columns = np.concatenate([pd.DataFrame(_) for _ in columns], axis=1)
        if self.charmm39:
            columns[:, 0] = -1
        np.savetxt(self.rtffile, columns[idx], fmt=native_str(self.fmt["MASS"]), delimiter=native_str(""))

    def _write_decl(self):
        decl = np.concatenate((self.atoms.names[:, np.newaxis], self.atoms.names[:, np.newaxis]), axis=1)
        np.savetxt(self.rtffile, decl, fmt=native_str(self.fmt["DECL"]))
        self.rtffile.write("\n".encode())

    def _write_residues(self, residue):
        self.rtffile.write(self.fmt["RES"].format(residue.resname, residue.charge).encode())

        # Write the atom lines with site name, type, and charge.
        key = "ATOM"
        if not np.issubdtype(residue.atoms.types.dtype, np.int):
            key += "_C36" if np.any(np.where([len(_) for _ in residue.atoms.types.astype(np.unicode)])[0] > 4) else ""
        atoms = residue.atoms
        lines = ((atoms.names, atoms.types, atoms.charges) if np.issubdtype(atoms.types.dtype, np.int)
                 else (atoms.names, atoms.names, atoms.charges))
        lines = np.concatenate([pd.DataFrame(_) for _ in lines], axis=1)
        np.savetxt(self.rtffile, lines, fmt=native_str(self.fmt[key]))

        # Write the bond, angle, dihedral, and improper dihedral lines.
        for key, value in self.bonds:
            attr, n_perline = value
            fmt = key + n_perline * " %5s"
            try:
                bonds = getattr(atoms, attr)
                if len(bonds) == 0:
                    continue

                # Create list of atom names and include "+" for atoms not within the residue.
                names = np.concatenate([_.atoms.names[np.newaxis, :] for _ in bonds], axis=0)
                idx = np.any([np.isin(_.atoms, atoms, invert=True) for _ in bonds], axis=1)
                pos_names = np.where(np.isin(bonds[idx], atoms, invert=True), "+", "").astype(np.object)
                names[idx] = pos_names + names[idx]
                names = names.astype(np.unicode)

                # Eliminate redundances.
                # Code courtesy of Daneil F on
                # https://stackoverflow.com/questions/45005477/eliminating-redundant-numpy-rows/45006131?noredirect=1#comment76988894_45006131
                b = np.ascontiguousarray(np.sort(names, -1)).view(np.dtype((np.void, names.dtype.itemsize * names.shape[1])))
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
        # type: (object) -> object
        self.atoms = universe.atoms
        with open(self.filename, "wb") as self.rtffile:
            # Write the title and header information.
            for _ in self.title:
                self.rtffile.write((_ + "\n").encode())
            self.rtffile.write(self.fmt["HEADER"].format(36, 1).encode())

            # Write the atom mass and declaration sections
            self._write_mass()
            self.rtffile.write("\n".encode())
            self._write_decl()
            self.rtffile.write("DEFA FIRS NONE LAST NONE\n".encode())
            self.rtffile.write("AUTOGENERATE ANGLES DIHEDRAL\n\n".encode())

            # Write out the residue information
            _, idx = np.unique(self.atoms.residues.resnames, return_index=True)
            for residue in self.atoms.residues[idx]:
                self._write_residues(residue)
            self.rtffile.write("END\n".encode())
