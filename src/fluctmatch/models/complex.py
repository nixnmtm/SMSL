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
    raise_from,
)

import numpy as np

import MDAnalysis as mda
from MDAnalysis.core import topologyattrs
from MDAnalysis.lib.distances import distance_array

from . import universe


class Complex(universe._Universe):
    _rmin = 0.
    _rmax = 10.

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._initialize(*args, **kwargs)

        # Numeric types for CHARMM PSF
        np.copyto(self.atoms.numtypes, self.atoms.types)

    def _initialize(self, *args, **kwargs):
        # Atomistic Universe
        try:
            self.atu = mda.Universe(*args, **kwargs)
        except (IOError, OSError, ValueError) as exc:
            raise_from(RuntimeError("Failed to create a universe."), exc)
        self._topology = self.atu._topology

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
