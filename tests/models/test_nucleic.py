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
from future.utils import native_str

import MDAnalysis as mda
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
        aa_universe.select_atoms(
            "nucleicphosphate and not name H*"
        ).residues.n_residues +
        aa_universe.select_atoms(
            "hnucleicsugar and not name H*"
        ).residues.n_residues +
        aa_universe.select_atoms(
            "hnucleicbase and not name H*"
        ).residues.n_residues
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
    dna = mda.Universe(PDB_dna).select_atoms("nucleic or resname OXG")
    for _ in dna.residues:
        positions.append(
            _.atoms.select_atoms(
                "nucleicphosphate and not name H*"
            ).center_of_mass()
        )
        positions.append(_.atoms.select_atoms(
            "hnucleicsugar and not name H*"
        ).center_of_mass())
        positions.append(_.atoms.select_atoms(
            "hnucleicbase and not name H*"
        ).center_of_mass())
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
        aa_universe.select_atoms(
            "nucleicphosphate and not name H*"
        ).residues.n_residues +
        aa_universe.select_atoms(
            "sugarC4 and not name H*"
        ).residues.n_residues +
        aa_universe.select_atoms(
            "sugarC2 and not name H*"
        ).residues.n_residues +
        aa_universe.select_atoms(
            "nucleiccenter and not name H*"
        ).residues.n_residues
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
    dna = mda.Universe(PDB_dna).select_atoms("nucleic or resname OXG")
    for _ in dna.residues:
        positions.append(
            _.atoms.select_atoms("nucleicphosphate and not name H*").center_of_mass()
        )
        positions.append(
            _.atoms.select_atoms("name C4'").center_of_mass()
        )
        positions.append(
            _.atoms.select_atoms("name C2'").center_of_mass()
        )
        positions.append(
            _.atoms.select_atoms(
                "nucleiccenter and not name H*"
            ).center_of_mass()
        )
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

def test_nucleic6_creation():
    h1 = (
        "(resname ADE DA* RA* and name N6) or "
        "(resname OXG GUA DG* RG* and name O6) or "
        "(resname CYT DC* RC* and name N4) or "
        "(resname THY URA DT* RU* and name O4)"
    )
    h2 = (
        "(resname ADE DA* RA* OXG GUA DG* RG* and name N1) or "
        "(resname CYT DC* RC* THY URA DT* RU* and name N3)"
    )
    h3 = (
        "(resname ADE DA* RA* and name H2) or "
        "(resname OXG GUA DG* RG* and name N2) or "
        "(resname CYT DC* RC* THY URA DT* RU* and name O2)"
    )

    aa_universe = mda.Universe(PDB_dna)
    n_atoms = (
        aa_universe.select_atoms("name P or name H5T").residues.n_residues +
        aa_universe.select_atoms("name C4'").residues.n_residues +
        aa_universe.select_atoms("name C2'").residues.n_residues +
        aa_universe.select_atoms(h1).residues.n_residues +
        aa_universe.select_atoms(h2).residues.n_residues +
        aa_universe.select_atoms(h3).residues.n_residues
    )
    cg_universe = nucleic.Nucleic6(PDB_dna)

    testing.assert_equal(
        cg_universe.atoms.n_atoms,
        n_atoms,
        err_msg=native_str("Number of sites do not match."),
        verbose=True,
    )
