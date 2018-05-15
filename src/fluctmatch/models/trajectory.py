# -*- coding: utf-8 -*-
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
# The original code is from Richard J. Gowers.
# https://github.com/richardjgowers/MDAnalysis-coarsegraining
#
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)
from future.utils import (
    viewitems,
    viewvalues,
)

import itertools

import MDAnalysis
from MDAnalysis.coordinates import base
from MDAnalysis.core import groups


class _Trajectory(base.ReaderBase):
    """Fakes a coarse grained trajectory object

    Takes an atomistic trajectory and list of beads and manipulates both
    to recreate a reader for a coarse grained Universe.

    Would also probably work as a standalone thingy for writing out
    coarse grained trajectories.
    """

    def __init__(self, universe, mapping, n_atoms=1, com=True):
        """

        Parameters
        ----------
        universe : :class:`~MDAnalysis.Universe` / :class:`~MDAnalysis.AtomGroup`
            A collection of atoms in a universe or AtomGroup.
        mapping : dict
            Definitions of the beads.
        n_atoms : int, optional
            Number of atoms in the coarse-grain system.
        convert_units : bool, optional
            units are converted to the MDAnalysis base format; None selects the
            value of MDAnalysis.core.flags [‘convert_lengths’].
        com : bool, optional
            Calculate center of mass or center of geometry per bead definition.
        kwargs : dict, optional
            Additonal arguments for use within the MDAnalysis coordinate reader.
        """
        self._u = universe
        self._t = universe.trajectory
        self.__dict__.update(universe.trajectory.__dict__)
        self._mapping = mapping
        residue_selection = itertools.product(self._u.residues,
                                              viewitems(self._mapping))

        self._beads = []
        for res, (key, selection) in residue_selection:
            if key != "CB":
                self._beads.append(res.atoms.select_atoms(selection))
            elif key == "CB":
                if isinstance(selection, dict):
                    value = selection.get(res.resname,
                                          "hsidechain and not name H*")
                    self._beads.append(res.atoms.select_atoms(value))
                else:
                    self._beads.append(res.atoms.select_atoms(selection))
            self._beads = [_ for _ in self._beads if _]

        self.com = com
        self._auxs = self._t._auxs
        try:
            self._frame = self._t._frame
        except AttributeError:
            pass

        self.n_atoms = n_atoms
        self.format = self._t.format
        self.units.update(self._t.units)
        self.convert_units = MDAnalysis.core.flags["convert_lengths"]
        try:
            self.fixed = self._t.fixed
        except AttributeError:
            self.fixed = False
        try:
            self.periodic = self._t.periodic
        except AttributeError:
            self.periodic = True

        self.ts = self._Timestep(
            self.n_atoms,
            positions=self._t.ts.has_positions,
            velocities=self._t.ts.has_velocities,
            forces=self._t.ts.has_forces)
        self._fill_ts(self._t.ts)

    def __iter__(self):
        self._reopen()

        yield self.ts
        while True:
            try:
                yield self._read_next_timestep()
            except StopIteration:
                self._reopen()
                raise StopIteration

    def __len__(self):
        #         return self.n_frames
        return len(self._u.trajectory)

    def __repr__(self):
        return "<CG Trajectory doing {:d} beads >".format(self.n_atoms)

    def _fill_ts(self, other_ts):
        """Rip information from atomistic TS into our ts

        Parameters
        ----------
        other_ts : :class:`~MDAnalysis.coordinates.base.Timestep`
            Another timestep
        """
        self.ts.frame = other_ts.frame
        self.ts.order = self._t.ts.order
        self.ts._unitcell = other_ts._unitcell
        self.ts.time = other_ts.time
        try:
            self.ts.dimensions = other_ts.dimensions
        except ValueError:
            dim = other_ts.dimensions.size
            self.ts.dimensions[:dim] = other_ts.dimensions
        self.ts.dt = other_ts.dt

        residues = self._u.atoms.split("residue")
        if self.ts.has_positions:
            if self.com:
                self.ts._pos[:] = [
                    bead.center_of_mass() for bead in self._beads
                ]
            else:
                self.ts._pos[:] = [
                    bead.center_of_geometry() for bead in self._beads
                ]

        residue_selection = itertools.product(residues, viewitems(
            self._mapping))
        if self.ts.has_velocities:
            try:
                self.ts._velocities[:] = [
                    res.select_atoms(selection).velocities.sum()
                    for res, (_, selection) in residue_selection
                    if res.select_atoms(selection)
                ]
            except ValueError:
                pass

        if self.ts.has_forces:
            try:
                self.ts._forces[:] = [
                    res.select_atoms(selection).forces.sum()
                    for res, (_, selection) in residue_selection
                    if res.select_atoms(selection)
                ]
            except ValueError:
                pass

    def _read_next_timestep(self, ts=None):
        # Get the next TS from the atom trajectory
        at_ts = self._t.next()

        self._fill_ts(at_ts)
        return self.ts

    def _read_frame(self, frame):
        self._t._read_frame(frame)
        self._fill_ts(self._t.ts)

        return self.ts

    def _reopen(self):
        # Rewind my reference trajectory
        self._read_frame(0)

    def close(self):
        """Close the trajectory file.
        """
        self._t.close()

    def rewind(self):
        """Position at beginning of trajectory"""
        self._reopen()

    @property
    def dimensions(self):
        """unitcell dimensions (*A*, *B*, *C*, *alpha*, *beta*, *gamma*)
        """
        return self.ts.dimensions

    @property
    def dt(self):
        """timestep between frames"""
        return self.ts.dt

    @property
    def n_frames(self):
        return self._t.n_frames
