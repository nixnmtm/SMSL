import numpy as np
import MDAnalysis as mda
from numpy import testing

from fluctmatch.models import protein, ligand
from tests.datafiles import PDB  # replace with a ligand-containing structure if needed

from tests.datafiles import (
    PDB,
    TPR,
    XTC,
)

def _expected_count_from_mapping(universe, model):
    expected = 0
    for res in universe.residues:
        if model.resname is not None and res.resname != model.resname:
            continue
        for _, sel in model._mapping.items():
            if res.atoms.select_atoms(sel).n_atoms > 0:
                expected += 1
    return expected


def _expected_positions(universe, model):
    positions = []
    for _, bead in model._iter_beads(universe):
        positions.append(bead.center_of_mass())
    return np.asarray(positions)

def test_adp_creation():
    aa = mda.Universe(PDB)
    cg = ligand.ADP(PDB)
    expected = _expected_count_from_mapping(aa, cg)

    testing.assert_equal(
        cg.atoms.n_atoms,
        expected,
        err_msg="Number of ADP sites not equal.",
        verbose=True,
    )


def test_adp_positions():
    aa = mda.Universe(PDB)
    cg = ligand.ADP(PDB)

    testing.assert_allclose(
        _expected_positions(aa, cg),
        cg.atoms.positions,
        err_msg="ADP coordinates do not match.",
    )


def test_om_creation():
    aa = mda.Universe(PDB)
    cg = ligand.OM(PDB)
    expected = _expected_count_from_mapping(aa, cg)

    testing.assert_equal(
        cg.atoms.n_atoms,
        expected,
        err_msg="Number of OM sites not equal.",
        verbose=True,
    )


def test_om_positions():
    aa = mda.Universe(PDB)
    cg = ligand.OM(PDB)

    testing.assert_allclose(
        _expected_positions(aa, cg),
        cg.atoms.positions,
        err_msg="OM coordinates do not match.",
    )