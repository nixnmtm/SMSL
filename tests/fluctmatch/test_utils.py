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

import numpy as np
import MDAnalysis as mda
from numpy import testing
from fluctmatch.fluctmatch import utils as fmutils
from ..datafiles import (
    TPR,
    XTC,
)


def test_average_structure():
    universe = mda.Universe(TPR, XTC)
    avg_positions = np.mean(
        [universe.atoms.positions for _ in universe.trajectory], axis=0)
    positions = fmutils.AverageStructure(universe.atoms).run().result
    testing.assert_allclose(
        positions,
        avg_positions,
        err_msg="Average coordinates don't match.",
    )


def test_average_bonds():
    universe = mda.Universe(TPR, XTC)
    avg_bonds = np.mean(
        [universe.bonds.bonds() for _ in universe.trajectory], axis=0)
    bonds = fmutils.BondAverage(universe.atoms).run().result
    testing.assert_allclose(
        bonds["r_IJ"],
        avg_bonds,
        err_msg="Average bond distances don't match.",
    )


def test_bond_fluctuation():
    universe = mda.Universe(TPR, XTC)
    avg_bonds = np.mean(
        [universe.bonds.bonds() for _ in universe.trajectory], axis=0)
    bond_fluct = np.std(
        [universe.bonds.bonds() for _ in universe.trajectory], axis=0)
    bonds = fmutils.BondStd(universe.atoms, average=avg_bonds).run().result
    testing.assert_allclose(
        bonds["r_IJ"],
        bond_fluct,
        err_msg="Bond fluctuations don't match.",
    )


def test_bond_average_contains_atom_names_and_distances():
    universe = mda.Universe(TPR, XTC)
    bonds = fmutils.BondAverage(universe.atoms).run().result

    assert list(bonds.columns) == ["I", "J", "r_IJ"]
    assert len(bonds) == len(universe.bonds)
    testing.assert_equal(bonds["I"].values, universe.bonds.atom1.names)
    testing.assert_equal(bonds["J"].values, universe.bonds.atom2.names)


class DummyFrame(object):
    def __init__(self, frame, time):
        self.frame = frame
        self.time = time


class DummyTrajectory(object):
    def __init__(self, times):
        self._frames = [
            DummyFrame(frame, time) for frame, time in enumerate(times)
        ]

    def __iter__(self):
        return iter(self._frames)

    def __getitem__(self, key):
        return self._frames[key]


class DummyAtoms(object):
    n_atoms = 2


class DummyUniverse(object):
    def __init__(self, times):
        self.atoms = DummyAtoms()
        self.trajectory = DummyTrajectory(times)


class DummyWriter(object):
    def __init__(self, filename, n_atoms, precision):
        self.filename = filename
        self.n_atoms = n_atoms
        self.precision = precision
        self.writes = []

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return False

    def write(self, atoms):
        self.writes.append(atoms)


def _patch_mda_split_io(monkeypatch, times):
    writers = []

    def universe_factory(topology, trajectory):
        return DummyUniverse(times)

    def writer_factory(filename, n_atoms, precision):
        writer = DummyWriter(filename, n_atoms, precision)
        writers.append(writer)
        return writer

    monkeypatch.setattr(fmutils.mda, "Universe", universe_factory)
    monkeypatch.setattr(fmutils.mda, "Writer", writer_factory)
    return writers


def test_split_mda_selects_one_based_frame_window(tmp_path, monkeypatch):
    writers = _patch_mda_split_io(monkeypatch, times=[0, 10, 20, 30, 40])

    fmutils.split_mda(
        ("3", 2, 4),
        data_dir=str(tmp_path),
        topology="topology.tpr",
        trajectory="trajectory.xtc",
        outfile="window.xtc",
        logfile="window.log",
        precision=3,
    )

    assert len(writers) == 1
    assert writers[0].n_atoms == DummyAtoms.n_atoms
    assert writers[0].precision == 3
    assert len(writers[0].writes) == 3
    assert (tmp_path / "3" / "window.log").read_text().splitlines() == [
        "Frames written: 3",
        "Start: 2",
        "Stop: 4",
        "Frame based: True",
    ]


def test_split_mda_selects_time_window_when_frame_based_is_false(
        tmp_path, monkeypatch):
    writers = _patch_mda_split_io(monkeypatch, times=[0, 5, 10, 15, 20])

    fmutils.split_mda(
        ("time", 6, 16),
        data_dir=str(tmp_path),
        frame_based=False,
        outfile="window.xtc",
        logfile="window.log",
    )

    assert len(writers) == 1
    assert len(writers[0].writes) == 2
    assert (tmp_path / "time" / "window.log").read_text().splitlines() == [
        "Frames written: 2",
        "Start: 6",
        "Stop: 16",
        "Frame based: False",
    ]
