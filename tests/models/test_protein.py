# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# SMSL -  https://github.com/nixnmtm/SMSL
# Copyright (c) 2013-2030 The fluctmatch Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the New BSD license.
#

import MDAnalysis as mda
from numpy import testing
from fluctmatch.models import protein
from fluctmatch.models.selection import *
from tests.datafiles import (
    PDB_prot,
    TPR,
    XTC,
)

def expected_counts_from_mapping(u, mapping):
    counts = {name: 0 for name in mapping}
    for res in u.residues:
        for name, sel in mapping.items():
            bead = res.atoms.select_atoms(sel)
            if bead.n_atoms > 0:
                counts[name] += 1
    return counts

def actual_counts_from_cg(cg):
    names = list(cg.atoms.names)
    return {name: names.count(name) for name in set(names)}


def test_calpha_creation():
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Calpha(PDB_prot)
    cg_natoms = (aa_universe.select_atoms("calpha").n_atoms +
                 aa_universe.select_atoms("bioion").n_atoms)
    testing.assert_equal(
        cg_universe.atoms.n_atoms,
        cg_natoms,
        err_msg="Number of sites not equal.",
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
        err_msg="The coordinates do not match.",
    )


def test_calpha_trajectory():
    aa_universe = mda.Universe(TPR, XTC)
    cg_universe = protein.Calpha(TPR, XTC)
    testing.assert_equal(
        cg_universe.trajectory.n_frames,
        aa_universe.trajectory.n_frames,
        err_msg="All-atom and coarse-grain trajectories unequal.",
        verbose=True,
    )

def test_caside_creation():
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Caside(PDB_prot)
    cg_natoms = (aa_universe.select_atoms("calpha").n_atoms +
                 aa_universe.select_atoms("cbeta").n_atoms +
                 aa_universe.select_atoms("bioion").n_atoms)
    testing.assert_equal(
        cg_universe.atoms.n_atoms,
        cg_natoms,
        err_msg="Number of sites not equal.",
        verbose=True,
    )


def test_caside_positions():
    positions = []
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Caside(PDB_prot)
    for _ in aa_universe.select_atoms("protein").residues:
        positions.append(_.atoms.select_atoms("calpha").center_of_mass())
        if _.resname != "GLY":
            cbeta = "hsidechain and not name H*"
            positions.append(_.atoms.select_atoms(cbeta).center_of_mass())
    for _ in aa_universe.select_atoms("bioion").residues:
        positions.append(_.atoms.center_of_mass())
    testing.assert_allclose(
        np.array(positions),
        cg_universe.atoms.positions,
        err_msg="The coordinates do not match.",
    )


def test_caside_trajectory():
    aa_universe = mda.Universe(TPR, XTC)
    cg_universe = protein.Caside(TPR, XTC)
    testing.assert_equal(
        cg_universe.trajectory.n_frames,
        aa_universe.trajectory.n_frames,
        err_msg="All-atom and coarse-grain trajectories unequal.",
        verbose=True,
    )

def _polar_cb_selection(cg_universe, res):
    """Return the Polar CB selection string for a residue, or None if absent."""
    return cg_universe._mapping["CB"].get(res.resname, None)


def _polar_has_cb(cg_universe, res):
    """Whether this residue should contribute a CB bead in Polar."""
    cb_sel = _polar_cb_selection(cg_universe, res)
    if cb_sel is None:
        return False
    return res.atoms.select_atoms(cb_sel).n_atoms > 0

def test_polar_creation():
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Polar(PDB_prot)

    expected = 0
    for res in aa_universe.residues:
        if res.atoms.select_atoms("protein and name N").n_atoms > 0:
            expected += 1

        if _polar_has_cb(cg_universe, res):
            expected += 1

        if res.atoms.select_atoms("protein and name O OT1 OT2 OXT").n_atoms > 0:
            expected += 1

    expected += aa_universe.select_atoms("bioion").n_atoms

    testing.assert_equal(
        cg_universe.atoms.n_atoms,
        expected,
        err_msg="Number of Polar sites not equal.",
        verbose=True,
    )


def test_polar_positions():
    positions = []
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Polar(PDB_prot)

    for res in aa_universe.select_atoms("protein").residues:
        positions.append(res.atoms.select_atoms("name N").center_of_mass())

        cb_sel = _polar_cb_selection(cg_universe, res)
        if cb_sel is not None:
            cb_atoms = res.atoms.select_atoms(cb_sel)
            if cb_atoms.n_atoms > 0:
                positions.append(cb_atoms.center_of_mass())

        positions.append(
            res.atoms.select_atoms("name O OT1 OT2 OXT").center_of_mass()
        )

    for ion in aa_universe.select_atoms("bioion").residues:
        positions.append(ion.atoms.center_of_mass())

    testing.assert_allclose(
        np.array(positions),
        cg_universe.atoms.positions,
        err_msg="The coordinates do not match.",
    )


def test_polar_trajectory():
    aa_universe = mda.Universe(TPR, XTC)
    cg_universe = protein.Polar(TPR, XTC)
    testing.assert_equal(
        cg_universe.trajectory.n_frames,
        aa_universe.trajectory.n_frames,
        err_msg="All-atom and coarse-grain trajectories unequal.",
        verbose=True,
    )