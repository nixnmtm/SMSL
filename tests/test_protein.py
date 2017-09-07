# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#

from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

import numpy as np
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import (
    PSF,
    CRD,
    TPR,
    XTC,
)


from fluctmatch.models import protein
from fluctmatch.models.selection import *


def test_calpha_creation():
    aa_universe = mda.Universe(PSF, CRD)
    ca_universe = protein.Calpha(PSF, CRD)
    ca = aa_universe.select_atoms("calpha")
    ca_natoms = ca.n_atoms
    assert ca_universe.atoms.n_atoms == ca_natoms

def test_calpha_positions():
    ca_pos = np.asarray([
        _.atoms.select_atoms("calpha").center_of_mass()
        for _ in mda.Universe(PSF, CRD).residues
    ])
    ca_universe = protein.Calpha(PSF, CRD).atoms.positions
    assert np.allclose(ca_pos, ca_universe)


def test_calpha_trajectory():
    aa_universe = mda.Universe(TPR, XTC)
    ca_universe = protein.Calpha(TPR, XTC)
    assert aa_universe.trajectory.n_frames == ca_universe.trajectory.n_frames


def test_caside_creation():
    aa_universe = mda.Universe(PSF, CRD)
    caside_universe = protein.Caside(PSF, CRD)
    caside_natoms = (
        aa_universe.select_atoms("calpha").n_atoms +
        aa_universe.select_atoms("cbeta").n_atoms
    )
    assert caside_universe.atoms.n_atoms == caside_natoms


def test_caside_positions():
    positions = list()
    caside_universe = protein.Caside(PSF, CRD)
    for _ in mda.Universe(PSF, CRD).residues:
        positions.append(_.atoms.select_atoms("calpha").center_of_mass())
        if _.resname != "GLY":
            positions.append(_.atoms.select_atoms("hsidechain").center_of_mass())
    assert np.allclose(np.asarray(positions), caside_universe.atoms.positions)


def test_caside_trajectory():
    aa_universe = mda.Universe(TPR, XTC)
    caside_universe = protein.Caside(TPR, XTC)
    assert aa_universe.trajectory.n_frames == caside_universe.trajectory.n_frames


def test_ncsc_creation():
    aa_universe = mda.Universe(PSF, CRD)
    ncsc_universe = protein.Ncsc(PSF, CRD)
    ncsc_natoms = (
        aa_universe.select_atoms("protein and name N").n_atoms +
        aa_universe.select_atoms("protein and name C").n_atoms +
        aa_universe.select_atoms("cbeta").n_atoms
    )
    assert ncsc_universe.atoms.n_atoms == ncsc_natoms


def test_ncsc_positions():
    positions = list()
    ncsc_universe = protein.Ncsc(PSF, CRD)
    for _ in mda.Universe(PSF, CRD).residues:
        positions.append(_.atoms.select_atoms("amine").center_of_mass())
        if _.resname != "GLY":
            positions.append(_.atoms.select_atoms("hsidechain").center_of_mass())
        positions.append(_.atoms.select_atoms("carboxyl").center_of_mass())
    assert np.allclose(np.asarray(positions), ncsc_universe.atoms.positions)


def test_ncsc_trajectory():
    aa_universe = mda.Universe(TPR, XTC)
    ncsc_universe = protein.Ncsc(TPR, XTC)
    assert aa_universe.trajectory.n_frames == ncsc_universe.trajectory.n_frames
