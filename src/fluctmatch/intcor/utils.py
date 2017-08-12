# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

from future.utils import (
    native_str,
    raise_from,
    reraise,
    with_metaclass,
)
from future.builtins import (
    ascii,
    bytes,
    chr,
    dict,
    filter,
    hex,
    input,
    map,
    next,
    oct,
    open,
    pow,
    range,
    round,
    str,
    super,
    zip,
)

import numpy as np
import pandas as pd

_HEADER = ["segidI", "resI", "I", "segidJ", "resJ", "J",
           "segidK", "resK", "K", "segidL", "resL", "L",
           "r_IJ", "T_IJK", "P_IJKL", "T_JKL", "r_KL"]

def create_empty_table(atomgroup):
    """Create an empty table of internal coordinates from an atomgroup

    :param atomgroup: An Atomgroup
    :return: pandas.Dataframe
    """
    table = pd.DataFrame()
    try:
        dihedrals = atomgroup.dihedrals
    except AttributeError:
        try:
            angles = atomgroup.angles
        except AttributeError:
            try:
                bonds = atomgroup.bonds
            except as exc:
                raise_from(AttributeError("Bonds, angles, and torsions undefined"), exc)
            else:
                n_bonds = len(bonds)
                atom1, atom2 = bonds.atom1, bonds.atom2
                zeros = np.zeros((n_bonds, 5), dtype=np.float)
                cols = pd.DataFrame(
                    [atom1.segids, atom1.resids, atom1.names,
                     atom2.segids, atom2.resids, atom2.names,
                     ["??", ] * n_bonds, ["??", ] * n_bonds, ["??", ] * n_bonds,
                     ["??", ] * n_bonds, ["??", ] * n_bonds, ["??", ] * n_bonds]
                ).T
                table = pd.concat([table, cols, zeros], axis=1)
        else:
            n_angles = len(angles)
            atom1, atom2, atom3 = angles.atom1, angles.atom2, angles.atom3
            zeros = np.zeros((n_angles, 5), dtype=np.float)
            cols = pd.DataFrame(
                [atom1.segids, atom1.resids, atom1.names,
                 atom2.segids, atom2.resids, atom2.names,
                 atom3.segids, atom3.resids, atom3.names,
                 ["??", ] * n_angles, ["??", ] * n_angles, ["??", ] * n_angles]
            ).T
            table = pd.concat([table, cols, zeros], axis=1)
    else:
        n_dihedrals = len(dihedrals)
        atom1, atom2, atom3, atom4 = dihedrals.atom1, dihedrals.atom2, dihedrals.atom3, dihedrals.atom4
        zeros = np.zeros((n_dihedrals, 5), dtype=np.float)
        cols = pd.DataFrame(
            [atom1.segids, atom1.resids, atom1.names,
             atom2.segids, atom2.resids, atom2.names,
             atom3.segids, atom3.resids, atom3.names,
             atom4.segids, atom4.resids, atom4.names]
        ).T
        table = pd.concat([table, cols, zeros], axis=1)

    table.columns = _HEADER
    return table
