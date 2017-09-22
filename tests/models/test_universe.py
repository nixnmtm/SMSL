# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

from future.utils import native_str
from numpy import testing

from fluctmatch.models import (
    protein,
    solvent,
    ions,
    universe,
)
from fluctmatch.models.selection import *
from tests.datafiles import (
    PDB_prot,
    TIP3P,
    IONS
)


def test_universe():
    testing.assert_raises(TypeError, universe._Universe, PDB_prot)


def test_merge_creation():
    prot = protein.Ncsc(PDB_prot)
    water = solvent.Water(TIP3P)
    solvions = ions.SolventIons(IONS)

    cg_universe = universe.Merge(prot, water, solvions)
    cg_natoms = prot.atoms.n_atoms + water.atoms.n_atoms + solvions.atoms.n_atoms
    testing.assert_equal(
        cg_universe.atoms.n_atoms,
        cg_natoms,
        err_msg=native_str("Number of sites don't match."),
        verbose=True,
    )


def test_merge_positions():
    prot = protein.Ncsc(PDB_prot)
    water = solvent.Water(TIP3P)
    solvions = ions.SolventIons(IONS)

    cg_universe = universe.Merge(prot, water, solvions)
    positions = np.concatenate(
        (prot.atoms.positions, water.atoms.positions, solvions.atoms.positions),
        axis=0
    )
    testing.assert_allclose(
        cg_universe.atoms.positions,
        positions,
        err_msg=native_str("Coordinates don't match."),
    )


def test_rename_universe():
    cg_universe = protein.Ncsc(PDB_prot)
    universe.rename_universe(cg_universe)
    testing.assert_string_equal(
        native_str(cg_universe.atoms[0].name),
        native_str("A001"),
    )
    testing.assert_string_equal(
        native_str(cg_universe.residues[0].resname),
        native_str("A001"),
    )
