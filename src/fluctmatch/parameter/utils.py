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

import MDAnalysis as mda
import numpy as np
import pandas as pd
from future.builtins import (
    dict, )


def create_empty_parameters(universe, **kwargs):
    """

    Parameters
    ----------
    universe : :class:`~MDAnalysis.Universe` or :class:`~MDAnalysis.AtomGroup`
        A collection of atoms in a universe or atomgroup with bond definitions.
    charmm40
        Allow for automatic atom typing (type set to -1).

    Returns
    -------
    An dictionary with keys defining the CHARMM parameter sections, and each
    associated value contianing a :class:`~pandas.DataFrame` with a table
    containing the necessary information. Tables for bonds, angles, dihedrals,
    and impropers will have values of 0 that can be filled by the user.
    """
    version = kwargs.get("charmm_version", 41)
    parameters = dict(
        ATOMS=pd.DataFrame(),
        BONDS=pd.DataFrame(),
        ANGLES=pd.DataFrame(),
        DIHEDRALS=pd.DataFrame(),
        IMPROPER=pd.DataFrame(),
    )
    param_columns = dict(
        ATOMS=["type", "atom", "mass"],
        BONDS=["I", "J", "Kb", "b0"],
        ANGLES=["I", "J", "K", "Ktheta", "theta0"],
        DIHEDRALS=["I", "J", "K", "L", "Kchi", "n", "delta"],
    )

    # Atoms
    types = (universe.atoms.types
             if np.issubdtype(universe.atoms.types.dtype, np.int) else
             np.arange(universe.atoms.n_atoms) + 1)
    atoms = [types, universe.atoms.names, universe.atoms.masses]
    parameters["ATOMS"] = pd.concat([pd.DataFrame(_) for _ in atoms], axis=1)
    parameters["ATOMS"].columns = param_columns["ATOMS"]
    if version > 39:
        parameters["ATOMS"]["type"] = -1

    # Bonds
    try:
        bonds = [
            universe.bonds.atom1.names,
            universe.bonds.atom2.names,
            np.zeros((universe.bonds.atom1.names.size, 2), dtype=np.float),
        ]
        parameters["BONDS"] = pd.concat(
            [pd.DataFrame(_) for _ in bonds], axis=1)
        parameters["BONDS"].columns = param_columns["BONDS"]
    except (mda.NoDataError, AttributeError, IndexError):
        pass

    # Angles
    try:
        angles = [
            universe.angles.atom1.names,
            universe.angles.atom2.names,
            universe.angles.atom3.names,
            np.zeros((universe.angles.atom1.names.size, 2), dtype=np.float),
        ]
        parameters["ANGLES"] = pd.concat(
            [pd.DataFrame(_) for _ in angles], axis=1)
        parameters["ANGLES"].columns = param_columns["ANGLES"]
    except (mda.NoDataError, AttributeError, IndexError):
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
        parameters["DIHEDRALS"] = pd.concat(
            [pd.DataFrame(_) for _ in dihedrals], axis=1)
        parameters["DIHEDRALS"].columns = param_columns["DIHEDRALS"]
    except (mda.NoDataError, AttributeError, IndexError):
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
        parameters["IMPROPER"] = pd.concat(
            [pd.DataFrame(_) for _ in impropers], axis=1)
        parameters["IMPROPER"].columns = param_columns["DIHEDRALS"]
    except (mda.NoDataError, AttributeError, IndexError):
        pass

    return parameters
