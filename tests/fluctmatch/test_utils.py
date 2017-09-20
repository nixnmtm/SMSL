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

from fluctmatch.fluctmatch import utils as fmutils
from ..datafiles import (
    TPR,
    XTC,
)


def test_average_structure():
    universe = mda.Universe(TPR, XTC)
    avg_positions = np.mean([
        universe.atoms.positions
        for _ in universe.trajectory
    ], axis=0)
    positions = fmutils.average_structure(TPR, XTC)
    assert np.allclose(positions, avg_positions)


def test_average_bonds():
    universe = mda.Universe(TPR, XTC)
    avg_bonds = np.mean([
        universe.bonds.bonds()
        for _ in universe.trajectory
    ], axis=0)
    bonds = fmutils.average_bonds(TPR, XTC)
    assert np.allclose(bonds, avg_bonds)
