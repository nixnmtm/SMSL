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

from future.builtins import (
    super,
)

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

    def __init__(self, trajectory, beads, convert_units=True, com=True, **kwargs):
        """
        Arguments:
            universe - the atomistic Universe you start with
            mapping - list of list of indices
        """
        super().__init__(trajectory.filename, convert_units=convert_units, **kwargs)
        self._t = trajectory
        self.beads = beads
        self.com = com

        self.n_atoms = self.beads.n_atoms
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
        if isinstance(self.beads._beads[0], groups.AtomGroup):
            if self.com:
                self.ts._pos[:] = [bead.center_of_mass() for bead in self.beads._beads]
            else:
                self.ts._pos[:] = [bead.center_of_geometry() for bead in self.beads._beads]
        elif isinstance(self.beads._beads[0], groups.Atom):
            self.ts._pos[:] = [bead.position for bead in self.beads._beads]
        else:
            raise TypeError("The bead should be either an Atom or an AtomGroup.")
        if hasattr(self.beads, "_velocities"):
            self.ts._velocities[:] = [bead.velocities.sum() for bead in self.beads._beads]
        if hasattr(self.beads, "_forces"):
            self.ts._forces[:] = [bead.forces.sum() for bead in self.beads._beads]

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
