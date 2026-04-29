# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# SMSL -  https://github.com/nixnmtm/SMSL
# Copyright (c) 2013-2030 The fluctmatch Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the New BSD license.
#

import numpy as np
import MDAnalysis as mda
from numpy import testing

from fluctmatch.models import protein
from tests.datafiles import (
    PDB_prot,
    TPR,
    XTC,
)


# ---------- shared helpers for capped proteins ----------

def _trueprotein_residues(u):
    return list(u.select_atoms("trueprotein").residues)


def _cap_before(u, res):
    idx = res.resindex
    if idx > 0:
        prev_res = u.residues[idx - 1]
        if prev_res.resname == "ACE":
            return prev_res
    return None


def _cap_after(u, res):
    idx = res.resindex
    if idx < len(u.residues) - 1:
        next_res = u.residues[idx + 1]
        if next_res.resname == "NME":
            return next_res
    return None


def _calpha_group(u, res):
    group = res.atoms.select_atoms("calpha")
    cap_b = _cap_before(u, res)
    cap_a = _cap_after(u, res)

    if cap_b is not None:
        group = group + cap_b.atoms
    if cap_a is not None:
        group = group + cap_a.atoms

    return group.unique


def _caside_ca_group(u, res):
    group = res.atoms.select_atoms("calpha")
    cap_b = _cap_before(u, res)
    cap_a = _cap_after(u, res)

    if cap_b is not None:
        group = group + cap_b.atoms
    if cap_a is not None:
        group = group + cap_a.atoms

    return group.unique


def _caside_cb_group(res):
    return res.atoms.select_atoms("hsidechain and not name H*").unique


def _polar_cb_selection(cg_universe, res):
    return cg_universe._mapping["CB"].get(res.resname, None)


def _polar_has_cb(cg_universe, res):
    cb_sel = _polar_cb_selection(cg_universe, res)
    if cb_sel is None:
        return False
    return res.atoms.select_atoms(cb_sel).n_atoms > 0


def _polar_n_group(u, res):
    group = res.atoms.select_atoms("amine") + res.atoms.select_atoms("hcalpha")
    cap = _cap_before(u, res)
    if cap is not None:
        group = group + cap.atoms
    return group.unique


def _polar_o_group(u, res):
    group = res.atoms.select_atoms("carboxyl") + res.atoms.select_atoms("hcalpha")
    cap = _cap_after(u, res)
    if cap is not None:
        group = group + cap.atoms
    return group.unique


# ---------- Calpha ----------

def test_calpha_creation():
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Calpha(PDB_prot)

    expected = len(_trueprotein_residues(aa_universe))
    expected += aa_universe.select_atoms("bioion").n_atoms

    testing.assert_equal(
        cg_universe.atoms.n_atoms,
        expected,
        err_msg="Number of sites not equal.",
        verbose=True,
    )


def test_calpha_positions():
    positions = []
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Calpha(PDB_prot)

    for res in _trueprotein_residues(aa_universe):
        ca_group = _calpha_group(aa_universe, res)
        if len(ca_group) > 0:
            positions.append(ca_group.center_of_mass())

    for ion in aa_universe.select_atoms("bioion").residues:
        positions.append(ion.atoms.center_of_mass())

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


def test_calpha_terminal_caps_are_incorporated():
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Calpha(PDB_prot)

    true_res = _trueprotein_residues(aa_universe)
    first_res = true_res[0]
    last_res = true_res[-1]

    expected_first = _calpha_group(aa_universe, first_res).center_of_mass()
    expected_last = _calpha_group(aa_universe, last_res).center_of_mass()

    cg_ca = cg_universe.atoms.select_atoms("name CA")

    testing.assert_allclose(cg_ca.positions[0], expected_first)
    testing.assert_allclose(cg_ca.positions[-1], expected_last)


# ---------- Caside ----------

def test_caside_creation():
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Caside(PDB_prot)

    expected = 0
    for res in _trueprotein_residues(aa_universe):
        if len(_caside_ca_group(aa_universe, res)) > 0:
            expected += 1
        if len(_caside_cb_group(res)) > 0:
            expected += 1

    expected += aa_universe.select_atoms("bioion").n_atoms

    testing.assert_equal(
        cg_universe.atoms.n_atoms,
        expected,
        err_msg="Number of sites not equal.",
        verbose=True,
    )


def test_caside_positions():
    positions = []
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Caside(PDB_prot)

    for res in _trueprotein_residues(aa_universe):
        ca_group = _caside_ca_group(aa_universe, res)
        if len(ca_group) > 0:
            positions.append(ca_group.center_of_mass())

        cb_group = _caside_cb_group(res)
        if len(cb_group) > 0:
            positions.append(cb_group.center_of_mass())

    for ion in aa_universe.select_atoms("bioion").residues:
        positions.append(ion.atoms.center_of_mass())

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


def test_caside_terminal_caps_are_incorporated_into_ca_only():
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Caside(PDB_prot)

    true_res = _trueprotein_residues(aa_universe)
    first_res = true_res[0]
    last_res = true_res[-1]

    expected_first_ca = _caside_ca_group(aa_universe, first_res).center_of_mass()
    expected_last_ca = _caside_ca_group(aa_universe, last_res).center_of_mass()

    cg_ca = cg_universe.atoms.select_atoms("name CA")
    testing.assert_allclose(cg_ca.positions[0], expected_first_ca)
    testing.assert_allclose(cg_ca.positions[-1], expected_last_ca)

    # CB remains pure sidechain; no cap contribution
    first_cb = _caside_cb_group(first_res)
    if len(first_cb) > 0:
        cg_first_res = cg_universe.residues[0]
        cg_cb = cg_first_res.atoms.select_atoms("name CB")
        if len(cg_cb) > 0:
            testing.assert_allclose(cg_cb.positions[0], first_cb.center_of_mass())


# ---------- Polar ----------

def test_polar_creation():
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Polar(PDB_prot)

    expected = 0
    for res in _trueprotein_residues(aa_universe):
        if len(_polar_n_group(aa_universe, res)) > 0:
            expected += 1

        if _polar_has_cb(cg_universe, res):
            expected += 1

        if len(_polar_o_group(aa_universe, res)) > 0:
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

    for res in _trueprotein_residues(aa_universe):
        n_group = _polar_n_group(aa_universe, res)
        if len(n_group) > 0:
            positions.append(n_group.center_of_mass())

        cb_sel = _polar_cb_selection(cg_universe, res)
        if cb_sel is not None:
            cb_group = res.atoms.select_atoms(cb_sel)
            if len(cb_group) > 0:
                positions.append(cb_group.center_of_mass())

        o_group = _polar_o_group(aa_universe, res)
        if len(o_group) > 0:
            positions.append(o_group.center_of_mass())

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


def test_polar_terminal_caps_are_incorporated():
    aa_universe = mda.Universe(PDB_prot)
    cg_universe = protein.Polar(PDB_prot)

    true_res = _trueprotein_residues(aa_universe)
    first_res = true_res[0]
    last_res = true_res[-1]

    expected_first_n = _polar_n_group(aa_universe, first_res).center_of_mass()
    expected_last_o = _polar_o_group(aa_universe, last_res).center_of_mass()

    cg_n = cg_universe.atoms.select_atoms("name N")
    cg_o = cg_universe.atoms.select_atoms("name O")

    testing.assert_allclose(cg_n.positions[0], expected_first_n)
    testing.assert_allclose(cg_o.positions[-1], expected_last_o)