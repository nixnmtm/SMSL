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

import MDAnalysis as mda
from future.utils import native_str
from numpy import testing

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
    testing.assert_equal(
        cg_universe.atoms.n_atoms,
        cg_natoms,
        err_msg=native_str("Number of sites not equal."),
        verbose=True,
    )


def test_calpha_positions():
    positions = []
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Calpha(PDB_prot)
    for _ in aa_universe.select_atoms("protein").residues:
        positions.append(_.atoms.select_atoms("calpha").center_of_mass())
    for _ in aa_universe.select_atoms("bioion").residues:
        positions.append(_.atoms.select_atoms("bioion").center_of_mass())
    testing.assert_allclose(
        np.array(positions),
        cg_universe.atoms.positions,
        err_msg=native_str("The coordinates do not match."),
    )


def test_calpha_trajectory():
    aa_universe = mda.Universe(TPR, XTC)
    cg_universe = protein.Calpha(TPR, XTC)
    testing.assert_equal(
        cg_universe.trajectory.n_frames,
        aa_universe.trajectory.n_frames,
        err_msg=native_str("All-atom and coarse-grain trajectories unequal."),
        verbose=True,
    )


def test_caside_creation():
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Caside(PDB_prot)
    cg_natoms = (
        aa_universe.select_atoms("calpha").n_atoms +
        aa_universe.select_atoms("cbeta").n_atoms +
        aa_universe.select_atoms("bioion").n_atoms
    )
    testing.assert_equal(
        cg_universe.atoms.n_atoms,
        cg_natoms,
        err_msg=native_str("Number of sites not equal."),
        verbose=True,
    )


def test_caside_positions():
    positions = []
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Caside(PDB_prot)
    for _ in aa_universe.select_atoms("protein").residues:
        positions.append(_.atoms.select_atoms("calpha").center_of_mass())
        if _.resname != "GLY":
            positions.append(
                _.atoms.select_atoms("hsidechain").center_of_mass()
            )
    for _ in aa_universe.select_atoms("bioion").residues:
        positions.append(_.atoms.center_of_mass())
    testing.assert_allclose(
        np.array(positions),
        cg_universe.atoms.positions,
        err_msg=native_str("The coordinates do not match."),
    )


def test_caside_trajectory():
    aa_universe = mda.Universe(TPR, XTC)
    cg_universe = protein.Caside(TPR, XTC)
    testing.assert_equal(
        cg_universe.trajectory.n_frames,
        aa_universe.trajectory.n_frames,
        err_msg=native_str("All-atom and coarse-grain trajectories unequal."),
        verbose=True,
    )


def test_ncsc_creation():
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Ncsc(PDB_prot)
    cg_natoms = (
        aa_universe.select_atoms("protein and name N").n_atoms +
        aa_universe.select_atoms("protein and name C").n_atoms +
        aa_universe.select_atoms("cbeta").n_atoms +
        aa_universe.select_atoms("bioion").n_atoms
    )
    testing.assert_equal(
        cg_universe.atoms.n_atoms,
        cg_natoms,
        err_msg=native_str("Number of sites not equal."),
        verbose=True,
    )


def test_ncsc_positions():
    positions = []
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Ncsc(PDB_prot)
    for _ in aa_universe.select_atoms("protein").residues:
        positions.append(_.atoms.select_atoms("amine").center_of_mass())
        if _.resname != "GLY":
            positions.append(
                _.atoms.select_atoms("hsidechain").center_of_mass()
            )
        positions.append(_.atoms.select_atoms("carboxyl").center_of_mass())
    for _ in aa_universe.select_atoms("bioion").residues:
        positions.append(_.atoms.center_of_mass())
    testing.assert_allclose(
        np.array(positions),
        cg_universe.atoms.positions,
        err_msg=native_str("The coordinates do not match."),
    )


def test_ncsc_trajectory():
    aa_universe = mda.Universe(TPR, XTC)
    cg_universe = protein.Ncsc(TPR, XTC)
    testing.assert_equal(
        cg_universe.trajectory.n_frames,
        aa_universe.trajectory.n_frames,
        err_msg=native_str("All-atom and coarse-grain trajectories unequal."),
        verbose=True,
    )
