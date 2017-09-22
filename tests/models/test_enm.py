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
from future.utils import native_str
from numpy import testing

from fluctmatch.models import enm
from tests.datafiles import (
    NCSC,
)


def test_enm_creation():
    aa_universe = mda.Universe(NCSC)
    cg_universe = enm.Enm(NCSC)
    cg_natoms = (
        aa_universe.select_atoms("all").n_atoms
    )
    testing.assert_equal(
        cg_universe.atoms.n_atoms,
        cg_natoms,
        err_msg=native_str("The number of beads don't match."),
        verbose=True,
    )


def test_enm_names():
    cg_universe = enm.Enm(NCSC)
    testing.assert_string_equal(
        native_str(cg_universe.atoms[0].name),
        native_str("A001"),
    )
    testing.assert_string_equal(
        native_str(cg_universe.residues[0].resname),
        native_str("A001"),
    )


def test_enm_positions():
    aa_universe = mda.Universe(NCSC)
    cg_universe = enm.Enm(NCSC)
    testing.assert_allclose(
        cg_universe.atoms.positions,
        aa_universe.atoms.positions,
        err_msg=native_str("Coordinates don't match."),
    )
