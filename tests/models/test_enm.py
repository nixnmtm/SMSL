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

from fluctmatch.models import enm
from fluctmatch.models.selection import *

from tests.datafiles import (
    ENM,
)


def test_enm_creation():
    aa_universe = mda.Universe(ENM)
    cg_universe = enm.Enm(ENM)
    cg_natoms = (
        aa_universe.select_atoms("all").n_atoms
    )
    assert cg_universe.atoms.n_atoms == cg_natoms
    assert cg_universe.atoms[0].name == "A001"
    assert cg_universe.residues[0].resname == "A001"


def test_enm_names():
    cg_universe = enm.Enm(ENM)
    assert (cg_universe.atoms[0].name == "A001") & (cg_universe.residues[0].resname == "A001")


def test_enm_positions():
    aa_universe = mda.Universe(ENM)
    cg_universe = enm.Enm(ENM)
    assert np.allclose(aa_universe.atoms.positions, cg_universe.atoms.positions)
