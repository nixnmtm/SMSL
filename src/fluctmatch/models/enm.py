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

import numpy as np
import MDAnalysis as mda
from MDAnalysis.core import topologyattrs
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.lib.distances import distance_array

from fluctmatch.fluctmatch import utils as fmutils
from fluctmatch.models.base import (
    ModelBase,
    rename_universe,
)


class Enm(ModelBase):
    model = "ENM"
    describe = "Elastic network model"

    def __init__(self, *args, **kwargs):
        self._rmin = kwargs.pop("rmin", 0.0)
        self._rmax = kwargs.pop("rmax", 10.0)
        super().__init__(*args, **kwargs)
        self._initialize(*args, **kwargs)

    def __repr__(self):
        message = "<CG Universe with {} beads".format(self.atoms.n_atoms)
        try:
            message += " and {:d} bonds".format(len(self._topology.bonds.values))
        except AttributeError:
            pass
        finally:
            message += ">"
        return message

    def _attach_coords_from_atu(self):
        coords = self.atu.atoms.positions.copy()
        self.load_new(coords[np.newaxis, :, :], format=MemoryReader)

    def _initialize(self, *args, **kwargs):
        self.__dict__.update(self.atu.__dict__)

        rename_universe(self)

        charges = kwargs.get("charges", False)
        if not charges:
            self._topology.add_TopologyAttr(
                topologyattrs.Charges(np.zeros(self.atoms.n_atoms))
            )

        self._topology.add_TopologyAttr(
            topologyattrs.Atomtypes(np.arange(self.atoms.n_atoms) + 1)
        )
        self._topology.add_TopologyAttr(topologyattrs.Angles([]))
        self._topology.add_TopologyAttr(topologyattrs.Dihedrals([]))
        self._topology.add_TopologyAttr(topologyattrs.Impropers([]))

        # Rebuild topology view once after topology attrs are updated
        mda.Universe.__init__(self, self._topology)

        self._add_bonds()

        # FINAL step: attach coordinates after all topology rebuilds
        self._attach_coords_from_atu()

        if kwargs.get("guess_angles", False):
            self._add_angles()
            self._add_dihedrals()
            self._add_impropers()

    def _add_bonds(self):
        positions = fmutils.AverageStructure(self.atu.atoms).run().result
        distmat = distance_array(positions, positions, backend="OpenMP")

        if self._rmin > 0.0:
            a0, a1 = np.where((distmat >= self._rmin) & (distmat <= self._rmax))
        else:
            a0, a1 = np.where((distmat > self._rmin) & (distmat <= self._rmax))

        skeleton_bonds = set(
            (int(bond[0].ix), int(bond[1].ix)) for bond in self.atu.bonds
        )
        rmax_bonds = set(
            (int(x), int(y)) for x, y in zip(a0, a1) if y > x
        )
        all_bonds = skeleton_bonds.union(rmax_bonds)

        self._topology.add_TopologyAttr(topologyattrs.Bonds(all_bonds))

        # Rebuild topology after adding bonds
        mda.Universe.__init__(self, self._topology)

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