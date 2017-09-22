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
import numpy.testing as testing
from future.utils import native_str

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
    testing.assert_allclose(
        positions,
        avg_positions,
        err_msg=native_str("Average coordinates don't match."),
    )


def test_average_bonds():
    universe = mda.Universe(TPR, XTC)
    avg_bonds = np.mean([
        universe.bonds.bonds()
        for _ in universe.trajectory
    ], axis=0)
    bonds = fmutils.average_bonds(TPR, XTC)
    testing.assert_allclose(
        bonds["r_IJ"],
        bond_fluct,
        err_msg=native_str("Average bond distances don't match."),
    )


def test_bond_fluctuation():
    universe = mda.Universe(TPR, XTC)
    bond_fluct = np.std([
        universe.bonds.bonds()
        for _ in universe.trajectory
    ], axis=0)
    bonds = fmutils.bond_fluctuation(TPR, XTC)
    testing.assert_allclose(
        bonds["r_IJ"],
        bond_fluct,
        err_msg=native_str("Bond fluctuations don't match."),
    )
