# -*- coding: utf-8 -*-
# coarse graining in MDAnalysis
# Copyright (c) 2015 Richard J Gowers
# Released under the GNU Lesser General Public License, version 2 or later.
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

from future.utils import (
    raise_with_traceback,
    viewitems,
    with_metaclass,
)
from future.builtins import (
    super,
)

import itertools

import MDAnalysis
from MDAnalysis.coordinates import base
from MDAnalysis.core import groups

from .selection import *


class _Trajectory(base.ReaderBase):
    """Fakes a coarse grained trajectory object

    Takes an atomistic trajectory and list of beads and manipulates both
    to recreate a reader for a coarse grained Universe.

    Would also probably work as a standalone thingy for writing out
    coarse grained trajectories.
    """

    def __init__(self, universe, mapping, n_atoms=1, convert_units=True, com=True, **kwargs):
        """
        Arguments:
            universe - the atomistic Universe you start with
            mapping - dictionary of selections
        """
        super().__init__(universe.trajectory.filename, convert_units=convert_units, **kwargs)
        self._u = universe
        self._t = universe.trajectory
        self._mapping = mapping
        self.com = com

        self.n_atoms = n_atoms
        self.n_frames = len(self._t)
        self.format = self._t.format
        self.units.update(self._t.units)
        self.convert_units = MDAnalysis.core.flags["convert_lengths"]

        self.fixed = False
        self.periodic = True
        self.skip = 1

        self.ts = self._Timestep(self.n_atoms, **self._t._ts_kwargs)
        self._frame = 0
        self.ts.dt = self._t.ts.dt
        self.ts.dimensions = self._t.ts.dimensions

        self._fill_ts(self._t.ts)

    def _read_next_timestep(self, ts=None):
        """Return the next timestep"""
        # Get the next TS from the atom trajectory
        at_ts = self._t.next()

        self._fill_ts(at_ts)

        return self.ts

    def _read_frame(self, frame):
        """Return a single frame"""
        at_ts = self._t[frame]

        self._fill_ts(at_ts)

        return self.ts

    def _fill_ts(self, other_ts):
        """Rip information from atomistic TS into our ts

        Make positions based on COM of our beads.
        """
        self.ts.frame = other_ts.frame
        self.ts._unitcell = other_ts._unitcell
        residues = self._u.atoms.split("residue")
        if self.com:
            self.ts._pos[:] = [
                res.select_atoms(selection).center_of_mass()
                for res, (_, selection) in itertools.product(residues, viewitems(self._mapping))
                if res.select_atoms(selection)
            ]
        else:
            self.ts._pos[:] = [
                res.select_atoms(selection).center_of_geometry()
                for res, (_, selection) in itertools.product(residues, viewitems(self._mapping))
                if res.select_atoms(selection)
            ]

        try:
            self.ts._velocities[:] = [
                res.select_atoms(selection).velocities.sum()
                for res, (_, selection) in itertools.product(residues, viewitems(self._mapping))
                if res.select_atoms(selection)
            ]
        except (AttributeError, MDAnalysis.NoDataError):
            pass

        try:
            self.ts._forces[:] = [
                res.select_atoms(selection).forces.sum()
                for res, (_, selection) in itertools.product(residues, viewitems(self._mapping))
                if res.select_atoms(selection)
            ]
        except (AttributeError, MDAnalysis.NoDataError):
            pass

    def _reopen(self):
        # Rewind my reference trajectory
        self._t.rewind()

    def __iter__(self):
        self._reopen()
        while True:
            try:
                yield self._read_next_timestep()
            except StopIteration:
                self.rewind()
                raise StopIteration

    def __len__(self):
        return self.n_frames

    def __repr__(self):
        return "<CG Trajectory doing {:d} beads >".format(self.n_atoms)

    def close(self):
        self._t.close()
