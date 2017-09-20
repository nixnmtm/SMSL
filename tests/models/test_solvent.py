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

from fluctmatch.models import (
    ions,
    solvent,
)
from fluctmatch.models.selection import *

from tests.datafiles import (
    TIP3P,
    TIP4P,
    IONS,
    DMA,
)


def test_water_from_tip3p_creation():
    aa_universe = mda.Universe(TIP3P)
    cg_universe = solvent.Water(TIP3P)
    cg_natoms = (
        aa_universe.select_atoms("name OW").n_atoms
    )
    assert cg_universe.atoms.n_atoms == cg_natoms


def test_water_from_tip3p_positions():
    positions = []
    aa_universe = mda.Universe(TIP3P)
    cg_universe = solvent.Water(TIP3P)
    for _ in aa_universe.select_atoms("name OW HW* MW").residues:
        positions.append(_.atoms.select_atoms("name OW HW* MW").center_of_mass())
    assert np.allclose(np.array(positions), cg_universe.atoms.positions)


def test_water_from_tip4p_creation():
    aa_universe = mda.Universe(TIP4P)
    cg_universe = solvent.Water(TIP4P)
    cg_natoms = (
        aa_universe.select_atoms("name OW").n_atoms
    )
    assert cg_universe.atoms.n_atoms == cg_natoms


def test_water_from_tip4p_positions():
    positions = []
    aa_universe = mda.Universe(TIP4P)
    cg_universe = solvent.Water(TIP4P)
    for _ in aa_universe.select_atoms("name OW HW* MW").residues:
        positions.append(_.atoms.select_atoms("name OW HW* MW").center_of_mass())
    assert np.allclose(np.array(positions), cg_universe.atoms.positions)


def test_tip3p_creation():
    aa_universe = mda.Universe(TIP3P)
    cg_universe = solvent.Tip3p(TIP3P)
    cg_natoms = (
        aa_universe.select_atoms("name OW").n_atoms +
        aa_universe.select_atoms("name HW1").n_atoms +
        aa_universe.select_atoms("name HW2").n_atoms
    )
    assert cg_universe.atoms.n_atoms == cg_natoms


def test_tip3p_positions():
    positions = []
    aa_universe = mda.Universe(TIP3P)
    cg_universe = solvent.Tip3p(TIP3P)
    for _ in aa_universe.select_atoms("name OW HW* MW").residues:
        positions.append(_.atoms.select_atoms("name OW MW").center_of_mass())
        positions.append(_.atoms.select_atoms("name HW1").center_of_mass())
        positions.append(_.atoms.select_atoms("name HW2").center_of_mass())
    assert np.allclose(np.array(positions), cg_universe.atoms.positions)


def test_ions_creation():
    aa_universe = mda.Universe(IONS)
    cg_universe = ions.SolventIons(IONS)
    cg_natoms = (
        aa_universe.select_atoms("name LI LIT K NA F CL BR I").n_atoms
    )
    assert cg_universe.atoms.n_atoms == cg_natoms


def test_ions_positions():
    positions = []
    aa_universe = mda.Universe(IONS)
    cg_universe = ions.SolventIons(IONS)
    for _ in aa_universe.select_atoms("name LI LIT K NA F CL BR I").residues:
        positions.append(_.atoms.select_atoms("name LI LIT K NA F CL BR I").center_of_mass())
    assert np.allclose(np.array(positions), cg_universe.atoms.positions)


def test_dma_creation():
    aa_universe = mda.Universe(DMA)
    cg_universe = solvent.Dma(DMA)
    cg_natoms = (
        aa_universe.select_atoms("name C1").n_atoms +
        aa_universe.select_atoms("name N").n_atoms +
        aa_universe.select_atoms("name C2").n_atoms +
        aa_universe.select_atoms("name C3").n_atoms
    )
    assert cg_universe.atoms.n_atoms == cg_natoms


def test_dma_positions():
    positions = []
    aa_universe = mda.Universe(DMA)
    cg_universe = solvent.Dma(DMA)
    for _ in aa_universe.select_atoms("resname DMA").residues:
        positions.append(_.atoms.select_atoms("resname DMA and name C1 H1*").center_of_mass())
        positions.append(_.atoms.select_atoms("resname DMA and name C N O").center_of_mass())
        positions.append(_.atoms.select_atoms("resname DMA and name C2 H2*").center_of_mass())
        positions.append(_.atoms.select_atoms("resname DMA and name C3 H3*").center_of_mass())
    assert np.allclose(np.array(positions), cg_universe.atoms.positions)
