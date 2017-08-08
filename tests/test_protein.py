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
from MDAnalysis.tests.datafiles import (TPR, XTC)

from fluctmatch.models import protein


def test_calpha_creation():
    aa_universe = mda.Universe(TPR, XTC)
    ca_universe = protein.Calpha(TPR, XTC)
    assert ca_universe.atoms.n_atoms == aa_universe.select_atoms("protein and name CA").n_atoms


def test_calpha_trajectory():
    aa_universe = mda.Universe(TPR, XTC)
    ca_universe = protein.Calpha(TPR, XTC)
    assert aa_universe.trajectory.n_frames == ca_universe.trajectory.n_frames

