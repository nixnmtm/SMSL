# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import MDAnalysis as mda
import numpy as np
import pandas as pd
from future.builtins import (
    dict,
)


def create_empty_parameters(universe):
    """Create an empty parameter set from a universe.

    :param universe: MDAnalysis.Universe
    :return: parameter dictionary
    """
    parameters = dict(
        ATOMS=pd.DataFrame(),
        BONDS=pd.DataFrame(),
        ANGLES=pd.DataFrame(),
        DIHEDRALS=pd.DataFrame(),
        IMPROPER=pd.DataFrame(),
    )
    prmcolumns = dict(
        ATOMS=["type", "atom", "mass"],
        BONDS=["I", "J", "Kb", "b0"],
        ANGLES=["I", "J", "K", "Ktheta", "theta0"],
        DIHEDRALS=["I", "J", "K", "L", "Kchi", "n", "delta"],
    )

    # Atoms
    print(universe.atoms.names)
    types = (
        universe.atoms.types
        if np.issubdtype(universe.atoms.types.dtype, np.int)
        else np.arange(universe.atoms.n_atoms) + 1
    )
    atoms = [types, universe.atoms.names, universe.atoms.masses]
    parameters["ATOMS"] = pd.concat([pd.DataFrame(_) for _ in atoms], axis=1)
    parameters["ATOMS"].columns = prmcolumns["ATOMS"]

    # Bonds
    try:
        bonds = [
            universe.bonds.atom1.names,
            universe.bonds.atom2.names,
            np.zeros((universe.bonds.atom1.names.size, 2), dtype=np.float),
        ]
        parameters["BONDS"] = pd.concat([pd.DataFrame(_) for _ in bonds], axis=1)
        parameters["BONDS"].columns = prmcolumns["BONDS"]
    except (mda.NoDataError, AttributeError):
        pass

    # Angles
    try:
        angles = [
            universe.angles.atom1.names,
            universe.angles.atom2.names,
            universe.angles.atom3.names,
            np.zeros((universe.angles.atom1.names.size, 2), dtype=np.float),
        ]
        parameters["ANGLES"] = pd.concat([pd.DataFrame(_) for _ in angles], axis=1)
        parameters["ANGLES"].columns = prmcolumns["ANGLES"]
    except (mda.NoDataError, AttributeError):
        pass

    # Dihedrals
    try:
        dihedrals = [
            universe.dihedrals.atom1.names,
            universe.dihedrals.atom2.names,
            universe.dihedrals.atom3.names,
            universe.dihedrals.atom4.names,
            np.zeros((universe.dihedrals.atom1.names.size, 1), dtype=np.float),
            np.zeros((universe.dihedrals.atom1.names.size, 1), dtype=np.int),
            np.zeros((universe.dihedrals.atom1.names.size, 1), dtype=np.float),
        ]
        parameters["DIHEDRALS"] = pd.concat([pd.DataFrame(_) for _ in dihedrals], axis=1)
        parameters["DIHEDRALS"].columns = prmcolumns["DIHEDRALS"]
    except (mda.NoDataError, AttributeError):
        pass

    # Impropers
    try:
        impropers = [
            universe.impropers.atom1.names,
            universe.impropers.atom2.names,
            universe.impropers.atom3.names,
            universe.impropers.atom4.names,
            np.zeros((universe.impropers.atom1.names.size, 1), dtype=np.float),
            np.zeros((universe.impropers.atom1.names.size, 1), dtype=np.int),
            np.zeros((universe.impropers.atom1.names.size, 1), dtype=np.float),
        ]
        parameters["IMPROPER"] = pd.concat([pd.DataFrame(_) for _ in impropers], axis=1)
        parameters["IMPROPER"].columns = prmcolumns["DIHEDRALS"]
    except (mda.NoDataError, AttributeError):
        pass

    return parameters
