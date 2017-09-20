# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import pytest

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
    with pytest.raises(TypeError):
        u = universe._Universe(PDB_prot)


def test_merge_creation():
    prot = protein.Ncsc(PDB_prot)
    water = solvent.Water(TIP3P)
    solvions = ions.SolventIons(IONS)

    cg_universe = universe.Merge(prot, water, solvions)
    cg_natoms = prot.atoms.n_atoms + water.atoms.n_atoms + solvions.atoms.n_atoms
    assert cg_natoms == cg_universe.atoms.n_atoms


def test_merge_positions():
    prot = protein.Ncsc(PDB_prot)
    water = solvent.Water(TIP3P)
    solvions = ions.SolventIons(IONS)

    cg_universe = universe.Merge(prot, water, solvions)
    positions = np.concatenate(
        (prot.atoms.positions, water.atoms.positions, solvions.atoms.positions),
        axis=0
    )
    assert np.allclose(positions, cg_universe.atoms.positions)


def test_rename_universe():
    cg_universe = protein.Ncsc(PDB_prot)
    universe.rename_universe(cg_universe)
    assert (cg_universe.atoms[0].name == "A001") & (cg_universe.residues[0].resname == "A001")
