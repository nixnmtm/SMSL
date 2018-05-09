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
from future.utils import native_str

import MDAnalysis as mda
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
