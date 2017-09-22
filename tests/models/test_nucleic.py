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

from fluctmatch.models import nucleic
from fluctmatch.models.selection import *
from tests.datafiles import (
    PDB_dna,
    TPR,
    XTC,
)


def test_nucleic3_creation():
    aa_universe = mda.Universe(PDB_dna)
    n_atoms = (
        aa_universe.select_atoms("nucleicphosphate").residues.n_residues +
        aa_universe.select_atoms("hnucleicsugar").residues.n_residues +
        aa_universe.select_atoms("hnucleicbase").residues.n_residues
    )
    cg_universe = nucleic.Nucleic3(PDB_dna)

    testing.assert_equal(
        cg_universe.atoms.n_atoms,
        n_atoms,
        err_msg=native_str("Number of sites do not match."),
        verbose=True,
    )


def test_nucleic3_positions():
    positions = []
    cg_universe = nucleic.Nucleic3(PDB_dna)
    for _ in mda.Universe(PDB_dna).select_atoms("nucleic or resname OXG").residues:
        positions.append(_.atoms.select_atoms("nucleicphosphate").center_of_mass())
        positions.append(_.atoms.select_atoms("hnucleicsugar").center_of_mass())
        positions.append(_.atoms.select_atoms("hnucleicbase").center_of_mass())
    testing.assert_allclose(
        np.array(positions),
        cg_universe.atoms.positions,
        err_msg=native_str("The coordinates do not match."),
    )


def test_nucleic3_trajectory():
    aa_universe = mda.Universe(TPR, XTC)
    cg_universe = nucleic.Nucleic3(TPR, XTC)
    testing.assert_equal(
        cg_universe.trajectory.n_frames,
        aa_universe.trajectory.n_frames,
        err_msg=native_str("All-atom and coarse-grain trajectories unequal."),
        verbose=True,
    )


def test_nucleic4_creation():
    aa_universe = mda.Universe(PDB_dna)
    n_atoms = (
        aa_universe.select_atoms("nucleicphosphate").residues.n_residues +
        aa_universe.select_atoms("sugarC4").residues.n_residues +
        aa_universe.select_atoms("sugarC3").residues.n_residues +
        aa_universe.select_atoms("nucleiccenter").residues.n_residues
    )
    cg_universe = nucleic.Nucleic4(PDB_dna)

    testing.assert_equal(
        cg_universe.atoms.n_atoms,
        n_atoms,
        err_msg=native_str("Number of sites do not match."),
        verbose=True,
    )


def test_nucleic4_positions():
    positions = []
    cg_universe = nucleic.Nucleic4(PDB_dna)
    for _ in mda.Universe(PDB_dna).select_atoms("nucleic or resname OXG").residues:
        positions.append(_.atoms.select_atoms("nucleicphosphate").center_of_mass())
        positions.append(_.atoms.select_atoms("sugarC4").center_of_mass())
        positions.append(_.atoms.select_atoms("sugarC3").center_of_mass())
        positions.append(_.atoms.select_atoms("nucleiccenter").center_of_mass())
    testing.assert_allclose(
        np.array(positions),
        cg_universe.atoms.positions,
        err_msg=native_str("The coordinates do not match."),
    )


def test_nucleic4_trajectory():
    aa_universe = mda.Universe(TPR, XTC)
    cg_universe = nucleic.Nucleic4(TPR, XTC)
    testing.assert_equal(
        cg_universe.trajectory.n_frames,
        aa_universe.trajectory.n_frames,
        err_msg=native_str("All-atom and coarse-grain trajectories unequal."),
        verbose=True,
    )
