# -*- coding: utf-8 -*-
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals,
)

from future.builtins import (
    super,
    zip,
)

from future.utils import (
    raise_with_traceback,
)

import numpy as np

import MDAnalysis as mda
from MDAnalysis.core import topologyattrs
from MDAnalysis.lib.distances import distance_array

from . import universe


class Enm(universe._Universe):
    _rmin = 0.
    _rmax = 10.

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._initialize(*args, **kwargs)

    def __repr__(self):
        message = "<CG Universe with {} beads".format(self.atoms.n_atoms)
        try:
            message += " and {:d} bonds".format(len(self._topology.bonds.values))
        except AttributeError as exc:
            pass
        finally:
            message += ">"
        return message

    def _initialize(self, *args, **kwargs):
        # Atomistic Universe
        try:
            self.atu = mda.Universe(*args, **kwargs)
        except (IOError, OSError, ValueError) as exc:
            raise_with_traceback(RuntimeError("Failed to create a universe."))
        self.__dict__.update(self.atu.__dict__)

        universe.rename_universe(self)
        self._generate_from_topology()

        self._add_bonds()
        if kwargs.get("guess_bonds", False):
            self._add_angles()
            self._add_dihedrals()
            self._add_impropers()

    def _add_bonds(self):
        dm = np.asarray([distance_array(self.atoms.positions, self.atoms.positions, backend="OpenMP") for _ in self.trajectory])
        if self._rmin > 0.:
            _, a0, a1 = np.where((dm >= self._rmin) & (dm <= self._rmax))
        else:
            _, a0, a1 = np.where((dm > self._rmin) & (dm <= self._rmax))
        bonds = topologyattrs.Bonds(set([(x, y) for x, y in zip(a0, a1) if y > x]))
        self._topology.add_TopologyAttr(bonds)
        self._generate_from_topology()

    @property
    def rmin(self):
        return self._rmin

    @rmin.setter
    def rmin(self, distance):
        self._rmin = distance

    @property
    def rmax(self):
        return self._rmax

    @rmax.setter
    def rmax(self, distance):
        self._rmax = distance
