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

    assert n_atoms == cg_universe.atoms.n_atoms


def test_nucleic3_positions():
    positions = []
    cg_universe = nucleic.Nucleic3(PDB_dna)
    for _ in mda.Universe(PDB_dna).select_atoms("nucleic or resname OXG").residues:
        positions.append(_.atoms.select_atoms("nucleicphosphate").center_of_mass())
        positions.append(_.atoms.select_atoms("hnucleicsugar").center_of_mass())
        positions.append(_.atoms.select_atoms("hnucleicbase").center_of_mass())
    assert np.allclose(np.array(positions), cg_universe.atoms.positions)


def test_nucleic3_trajectory():
    aa_universe = mda.Universe(TPR, XTC)
    cg_universe = nucleic.Nucleic3(TPR, XTC)
    assert aa_universe.trajectory.n_frames == cg_universe.trajectory.n_frames


def test_nucleic4_creation():
    aa_universe = mda.Universe(PDB_dna)
    n_atoms = (
        aa_universe.select_atoms("nucleicphosphate").residues.n_residues +
        aa_universe.select_atoms("sugarC4").residues.n_residues +
        aa_universe.select_atoms("sugarC3").residues.n_residues +
        aa_universe.select_atoms("nucleiccenter").residues.n_residues
    )
    cg_universe = nucleic.Nucleic4(PDB_dna)

    assert n_atoms == cg_universe.atoms.n_atoms


def test_nucleic4_positions():
    positions = []
    cg_universe = nucleic.Nucleic4(PDB_dna)
    for _ in mda.Universe(PDB_dna).select_atoms("nucleic or resname OXG").residues:
        positions.append(_.atoms.select_atoms("nucleicphosphate").center_of_mass())
        positions.append(_.atoms.select_atoms("sugarC4").center_of_mass())
        positions.append(_.atoms.select_atoms("sugarC3").center_of_mass())
        positions.append(_.atoms.select_atoms("nucleiccenter").center_of_mass())
    assert np.allclose(np.array(positions), cg_universe.atoms.positions)


def test_nucleic4_trajectory():
    aa_universe = mda.Universe(TPR, XTC)
    cg_universe = nucleic.Nucleic4(TPR, XTC)
    assert aa_universe.trajectory.n_frames == cg_universe.trajectory.n_frames
