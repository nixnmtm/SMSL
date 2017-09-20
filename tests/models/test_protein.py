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

from fluctmatch.models import protein
from fluctmatch.models.selection import *

from tests.datafiles import (
    PDB_prot,
    TPR,
    XTC,
)


def test_calpha_creation():
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Calpha(PDB_prot)
    cg_natoms = (
        aa_universe.select_atoms("calpha").n_atoms +
        aa_universe.select_atoms("bioion").n_atoms
    )
    assert cg_universe.atoms.n_atoms == cg_natoms


def test_calpha_positions():
    positions = []
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Calpha(PDB_prot)
    for _ in aa_universe.select_atoms("protein").residues:
        positions.append(_.atoms.select_atoms("calpha").center_of_mass())
    for _ in aa_universe.select_atoms("bioion").residues:
        positions.append(_.atoms.select_atoms("bioion").center_of_mass())
    assert np.allclose(np.array(positions), cg_universe.atoms.positions)


def test_calpha_trajectory():
    aa_universe = mda.Universe(TPR, XTC)
    cg_universe = protein.Calpha(TPR, XTC)
    assert aa_universe.trajectory.n_frames == cg_universe.trajectory.n_frames


def test_caside_creation():
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Caside(PDB_prot)
    cg_natoms = (
        aa_universe.select_atoms("calpha").n_atoms +
        aa_universe.select_atoms("cbeta").n_atoms +
        aa_universe.select_atoms("bioion").n_atoms
    )
    assert cg_universe.atoms.n_atoms == cg_natoms


def test_caside_positions():
    positions = []
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Caside(PDB_prot)
    for _ in aa_universe.select_atoms("protein").residues:
        positions.append(_.atoms.select_atoms("calpha").center_of_mass())
        if _.resname != "GLY":
            positions.append(_.atoms.select_atoms("hsidechain").center_of_mass())
    for _ in aa_universe.select_atoms("bioion").residues:
        positions.append(_.atoms.center_of_mass())
    assert np.allclose(np.array(positions), cg_universe.atoms.positions)


def test_caside_trajectory():
    aa_universe = mda.Universe(TPR, XTC)
    cg_universe = protein.Caside(TPR, XTC)
    assert aa_universe.trajectory.n_frames == cg_universe.trajectory.n_frames


def test_ncsc_creation():
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Ncsc(PDB_prot)
    cg_natoms = (
        aa_universe.select_atoms("protein and name N").n_atoms +
        aa_universe.select_atoms("protein and name C").n_atoms +
        aa_universe.select_atoms("cbeta").n_atoms +
        aa_universe.select_atoms("bioion").n_atoms
    )
    assert cg_universe.atoms.n_atoms == cg_natoms


def test_ncsc_positions():
    positions = []
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Ncsc(PDB_prot)
    for _ in aa_universe.select_atoms("protein").residues:
        positions.append(_.atoms.select_atoms("amine").center_of_mass())
        if _.resname != "GLY":
            positions.append(_.atoms.select_atoms("hsidechain").center_of_mass())
        positions.append(_.atoms.select_atoms("carboxyl").center_of_mass())
    for _ in aa_universe.select_atoms("bioion").residues:
        positions.append(_.atoms.center_of_mass())
    assert np.allclose(np.array(positions), cg_universe.atoms.positions)


def test_ncsc_trajectory():
    aa_universe = mda.Universe(TPR, XTC)
    cg_universe = protein.Ncsc(TPR, XTC)
    assert aa_universe.trajectory.n_frames == cg_universe.trajectory.n_frames
