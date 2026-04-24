# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# SMSL - https://github.com/nixnmtm/SMSL
#

import itertools
from future.utils import viewitems
from MDAnalysis.coordinates import base


class _Trajectory(base.ReaderBase):
    """Fakes a coarse grained trajectory object."""

    def __init__(self, universe, mapping, n_atoms=1, com=True, convert_units=None):
        self._u = universe
        self._t = universe.atu.trajectory if hasattr(universe, "atu") else universe.trajectory
        self.__dict__.update(self._t.__dict__)
        self._mapping = mapping

        # Use model-defined bead iterator when available
        if hasattr(self._u, "_iter_beads"):
            self._beads = [bead for _, bead in self._u._iter_beads(self._u.atu)]
        else:
            residue_selection = itertools.product(
                self._u.residues, viewitems(self._mapping)
            )
            self._beads = []
            for res, (key, selection) in residue_selection:
                if key != "CB":
                    bead = res.atoms.select_atoms(selection)
                elif isinstance(selection, dict):
                    value = selection.get(res.resname, "hsidechain and not name H*")
                    bead = res.atoms.select_atoms(value)
                else:
                    bead = res.atoms.select_atoms(selection)

                if bead.n_atoms > 0:
                    self._beads.append(bead)

        self.com = com
        self._auxs = self._t._auxs
        try:
            self._frame = self._t._frame
        except AttributeError:
            pass

        self.n_atoms = n_atoms
        self.format = self._t.format
        self.units.update(self._t.units)

        if convert_units is None:
            convert_units = True
        self.convert_units = convert_units

        try:
            self.fixed = self._t.fixed
        except AttributeError:
            self.fixed = False
        try:
            self.periodic = self._t.periodic
        except AttributeError:
            self.periodic = True

        if len(self._beads) != self.n_atoms:
            raise ValueError(
                "_Trajectory bead mismatch: "
                f"{len(self._beads)} trajectory beads vs {self.n_atoms} topology atoms"
            )

        self.ts = self._Timestep(
            self.n_atoms,
            positions=self._t.ts.has_positions,
            velocities=self._t.ts.has_velocities,
            forces=self._t.ts.has_forces,
        )
        self._fill_ts(self._t.ts)

    def __iter__(self):
        self._reopen()
        yield self.ts
        while True:
            try:
                yield self._read_next_timestep()
            except StopIteration:
                self._reopen()
                return

    def __len__(self):
        return self._t.n_frames

    def __repr__(self):
        return "<CG Trajectory doing {:d} beads >".format(self.n_atoms)

    def _fill_ts(self, other_ts):
        self.ts.frame = other_ts.frame
        try:
            self.ts.order = self._t.ts.order
        except (AttributeError, TypeError):
            pass

        self.ts._unitcell = other_ts._unitcell
        self.ts.time = other_ts.time
        try:
            self.ts.dimensions = other_ts.dimensions
        except ValueError:
            dim = other_ts.dimensions.size
            self.ts.dimensions[:dim] = other_ts.dimensions
        self.ts.dt = other_ts.dt

        if self.ts.has_positions:
            if self.com:
                self.ts._pos[:] = [bead.center_of_mass() for bead in self._beads]
            else:
                self.ts._pos[:] = [bead.center_of_geometry() for bead in self._beads]

        if self.ts.has_velocities:
            try:
                self.ts._velocities[:] = [
                    bead.velocities.sum(axis=0) for bead in self._beads
                ]
            except (ValueError, AttributeError):
                pass

        if self.ts.has_forces:
            try:
                self.ts._forces[:] = [
                    bead.forces.sum(axis=0) for bead in self._beads
                ]
            except (ValueError, AttributeError):
                pass

    def _read_next_timestep(self, ts=None):
        at_ts = next(self._t)
        self._fill_ts(at_ts)
        return self.ts

    def _read_frame(self, frame):
        self._t._read_frame(frame)
        self._fill_ts(self._t.ts)
        return self.ts

    def _reopen(self):
        self._read_frame(0)

    def close(self):
        self._t.close()

    def rewind(self):
        self._reopen()

    @property
    def dimensions(self):
        return self.ts.dimensions

    @property
    def dt(self):
        return self.ts.dt

    @property
    def n_frames(self):
        return self._t.n_frames