# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import numpy as np
import pandas as pd
from future.utils import (
    raise_with_traceback,
)

_HEADER = ["segidI", "resI", "I", "segidJ", "resJ", "J",
           "segidK", "resK", "K", "segidL", "resL", "L",
           "r_IJ", "T_IJK", "P_IJKL", "T_JKL", "r_KL"]


def create_empty_table(universe):
    """Create an empty table of internal coordinates from an atomgroup

    Parameters
    ----------
    universe : :class:`~MDAnalysis.Universe` or :class:`~MDAnalysis.AtomGroup`
        A collection of atoms in a universe or atomgroup with bond definitions.

    Returns
    -------
    A :class:`~pandas.DataFrame` compliant with a CHARMM-formatted internal
    coordinates (IC) table. The table matches the 'resid' version of an IC table.
    """
    table = pd.DataFrame()
    atomgroup = universe.atoms
    try:
        dihedrals = atomgroup.dihedrals
        if len(dihedrals) == 0:
            raise AttributeError
    except AttributeError:
        try:
            angles = atomgroup.angles
            if len(angles) == 0:
                raise AttributeError
        except AttributeError:
            try:
                bonds = atomgroup.bonds
                if len(bonds) == 0:
                    raise AttributeError
            except AttributeError:
                raise_with_traceback(AttributeError("Bonds, angles, and torsions undefined"))
            else:
                n_bonds = len(bonds)
                atom1, atom2 = bonds.atom1, bonds.atom2
                zeros = pd.DataFrame(np.zeros((n_bonds, 5), dtype=np.float))
                cols = pd.DataFrame(
                    [atom1.segids, atom1.resnums, atom1.names,
                     atom2.segids, atom2.resnums, atom2.names,
                     ["??", ] * n_bonds, ["??", ] * n_bonds, ["??", ] * n_bonds,
                     ["??", ] * n_bonds, ["??", ] * n_bonds, ["??", ] * n_bonds]
                ).T
                table = pd.concat([table, cols, zeros], axis=1)
        else:
            n_angles = len(angles)
            atom1, atom2, atom3 = angles.atom1, angles.atom2, angles.atom3
            zeros = pd.DataFrame(np.zeros((n_angles, 5), dtype=np.float))
            cols = pd.DataFrame(
                [atom1.segids, atom1.resnums, atom1.names,
                 atom2.segids, atom2.resnums, atom2.names,
                 atom3.segids, atom3.resnums, atom3.names,
                 ["??", ] * n_angles, ["??", ] * n_angles, ["??", ] * n_angles]
            ).T
            table = pd.concat([table, cols, zeros], axis=1)
    else:
        n_dihedrals = len(dihedrals)
        atom1, atom2, atom3, atom4 = dihedrals.atom1, dihedrals.atom2, dihedrals.atom3, dihedrals.atom4
        zeros = pd.DataFrame(np.zeros((n_dihedrals, 5), dtype=np.float))
        cols = pd.DataFrame([
            atom1.segids, atom1.resnums, atom1.names,
            atom2.segids, atom2.resnums, atom2.names,
            atom3.segids, atom3.resnums, atom3.names,
            atom4.segids, atom4.resnums, atom4.names
        ]).T
        table = pd.concat([table, cols, zeros], axis=1)

    table.columns = _HEADER
    return table
